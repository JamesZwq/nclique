//
// R2 DCLP: Dual-CSR Lazy Peeling
//
// Key innovations over V8b:
//
// 1. Edge-Path Transpose Index (edgePathCSR):
//    For each edge, stores which CPI paths contain it + type (KK/PP/KP) + pivot info.
//    Replaces Phase 1 hash-based intersection (treeGraphV) with direct CSR lookup.
//    Complexity: O(paths_per_edge) sequential scan vs O(min(deg(u),deg(v))) hash probes.
//
// 2. CSR-based Case C:
//    Uses leafEdgeCSR (from V8b) to scan Case C edges directly,
//    avoiding O(|L|^2) edge pair enumeration + getEdgeCompressedId binary search.
//    For Case C alive (leaf survives with fewer pivots), scan CSR and apply
//    type-dependent deltas, checking removed pivots via getEdgeById (O(1)).
//
// 3. Case A / Case B: Same as V8b (CSR scan / BK respectively).
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <cstdlib>
#include <set>
#include <unordered_map>
#include <boost/heap/d_ary_heap.hpp>

#include "../BK/BronKerboschRmEdge.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

#ifndef R2_LAB_NAMESPACE
#define R2_LAB_NAMESPACE DCLPHybridLab
#endif

#ifndef R2_LAB_ENTRYPOINT
#define R2_LAB_ENTRYPOINT PlusNucleusEdgeCoreDecompositionSet_HybridLab
#endif

#ifndef R2_LAB_DCLP_LABEL
#define R2_LAB_DCLP_LABEL "DCLP"
#endif

#ifndef R2_LAB_HYBRID_LABEL
#define R2_LAB_HYBRID_LABEL "Hybrid R2 Lab"
#endif

#ifndef R2_LAB_ENABLE_CASEB_LEAF_REUSE
#define R2_LAB_ENABLE_CASEB_LEAF_REUSE 0
#endif

#ifndef R2_LAB_ENABLE_CASEB_ZEROEDGE_CSR
#define R2_LAB_ENABLE_CASEB_ZEROEDGE_CSR 0
#endif

#ifndef R2_LAB_ENABLE_DEFECT_D2_IMMEDIATE
#define R2_LAB_ENABLE_DEFECT_D2_IMMEDIATE 0
#endif

#ifndef R2_LAB_ENABLE_DEFECT_D2_ORBIT
#define R2_LAB_ENABLE_DEFECT_D2_ORBIT 0
#endif

#ifndef R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
#define R2_LAB_ENABLE_CASEB_DEFECT_GRAPH 0
#endif

#ifndef R2_LAB_DEFECT_MAX_CORE_VERTS
#define R2_LAB_DEFECT_MAX_CORE_VERTS 60
#endif

namespace R2_LAB_NAMESPACE {

namespace DCLP {

    enum EdgeType : uint8_t { KK = 0, PP = 1, KP = 2 };

    struct LeafEdgeEntry {
        daf::Size edgeId;
        EdgeType type;
        uint8_t kpPivotIsLow;
    };

    struct EdgeLeafEntry {
        daf::Size leafId;
        EdgeType type;
        uint8_t kpPivotIsLow;
    };

    struct EdgeLeafTranspose {
        std::vector<daf::Size> edgeLeafOff;
        std::vector<EdgeLeafEntry> edgeLeafData;
        std::vector<uint8_t> edgeIndexed;
        daf::Size indexedEdges = 0;
    };

    struct LeafRmInfo {
        bool removedKeepC;
        daf::StaticVector<daf::Size> removedPivots{0};
        daf::StaticVector<std::pair<daf::Size, daf::Size>> removedEdges{0};

        LeafRmInfo() : removedKeepC(false) {}

        bool empty() const {
            return !removedKeepC && removedPivots.empty() && removedEdges.empty();
        }

        void init(auto capacity = 400) {
            removedKeepC = false;
            removedPivots.reserve(capacity);
            removedEdges.reserve(capacity);
        }

        void clear() {
            removedKeepC = false;
            removedPivots.clear();
            removedEdges.clear();
        }
    };

    template<typename It1, typename It2, typename WeightT, typename UpdateFunc>
    inline void processEdgePairsImpl(It1 b1, It1 e1,
                                     It2 b2, It2 e2,
                                     WeightT weight,
                                     UpdateFunc &&upd) noexcept {
        if (weight < 0) return;
        if (b1 == b2 && e1 == e2 && b1 == e1 && b2 == e2) return;
        if (b1 == b2 && e1 == e2) {
            for (auto it = b1; it + 1 != e1; ++it) {
                auto u = *it;
                for (auto jt = it + 1; jt != e1; ++jt) {
                    upd(u, *jt, weight);
                }
            }
        } else {
            for (auto it = b1; it != e1; ++it) {
                auto u = *it;
                for (auto jt = b2; jt != e2; ++jt) {
                    upd(u, *jt, weight);
                }
            }
        }
    }

    template<typename Range1, typename Range2, typename WeightT, typename UpdateFunc>
    inline void processEdgePairs(const Range1 &r1, const Range2 &r2,
                                 WeightT weight, UpdateFunc &&upd) noexcept {
        processEdgePairsImpl(
            std::begin(r1), std::end(r1),
            std::begin(r2), std::end(r2),
            weight, std::forward<UpdateFunc>(upd));
    }

    template<typename Range, typename WeightT, typename UpdateFunc>
    inline void processEdgePairs(const Range &r, WeightT weight,
                                 UpdateFunc &&upd) noexcept {
        processEdgePairsImpl(
            std::begin(r), std::end(r),
            std::begin(r), std::end(r),
            weight, std::forward<UpdateFunc>(upd));
    }

#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
    struct DefectContributionWeights {
        bool supported = false;
        double cliqueCount = 0;
        double isolatedKeepPivotWeight = 0;
        double isolatedPairWeight = 0;
        std::vector<daf::Size> corePivotVerts;
        std::vector<daf::Size> isolatedPivotVerts;
        std::vector<double> coreKeepPivotWeights;
        std::vector<double> corePairWithIsolatedWeights;
        std::vector<std::tuple<int, int, double>> corePairWeights;
    };

    inline double chooseBounded(long long n, int r) {
        if (n < 0 || r < 0 || r > n) return 0;
        return nCr[n][r];
    }

    inline double evalDefectPolynomial(const std::vector<double> &poly,
                                       long long isolatedCount,
                                       int targetSize) {
        if (targetSize < 0) return 0;
        double total = 0;
        const int maxCoreSize = std::min<int>(targetSize, static_cast<int>(poly.size()) - 1);
        for (int coreTake = 0; coreTake <= maxCoreSize; ++coreTake) {
            if (poly[coreTake] == 0) continue;
            total += poly[coreTake] * chooseBounded(isolatedCount, targetSize - coreTake);
        }
        return total;
    }

    class DefectPolynomialCounter {
    public:
        DefectPolynomialCounter(std::vector<uint64_t> closedMasks, int maxNeed)
            : closedMasks_(std::move(closedMasks)), maxNeed_(std::max(0, maxNeed)) {}

        const std::vector<double> &solve(uint64_t mask) {
            auto it = memo_.find(mask);
            if (it != memo_.end()) return it->second;

            std::vector<double> result(maxNeed_ + 1, 0);
            if (mask == 0) {
                result[0] = 1;
                auto [insIt, _] = memo_.emplace(mask, std::move(result));
                return insIt->second;
            }

            const int v = __builtin_ctzll(mask);
            const auto &skipPoly = solve(mask & ~(uint64_t(1) << v));
            const auto &takePoly = solve(mask & ~closedMasks_[v]);

            for (int i = 0; i <= maxNeed_; ++i) {
                result[i] += skipPoly[i];
                if (i > 0) result[i] += takePoly[i - 1];
            }

            auto [insIt, _] = memo_.emplace(mask, std::move(result));
            return insIt->second;
        }

    private:
        std::vector<uint64_t> closedMasks_;
        int maxNeed_;
        std::unordered_map<uint64_t, std::vector<double>> memo_;
    };

    inline bool buildDefectContributionWeights(
        const daf::StaticVector<daf::Size> &keepVerts,
        const daf::StaticVector<daf::Size> &pivotVerts,
        const std::vector<std::pair<daf::Size, daf::Size>> &defectEdges,
        daf::CliqueSize k,
        DefectContributionWeights &out) {
        out = {};

        const int needPivot = static_cast<int>(k) - static_cast<int>(keepVerts.size());
        const int pivotCount = static_cast<int>(pivotVerts.size());
        if (needPivot < 0 || needPivot > pivotCount) {
            out.supported = true;
            return true;
        }

        std::unordered_map<daf::Size, int> pivotPos;
        pivotPos.reserve(static_cast<size_t>(pivotCount * 2 + 8));
        for (int i = 0; i < pivotCount; ++i)
            pivotPos.emplace(pivotVerts[i], i);

        std::vector<int> pivotIsCore(pivotCount, 0);
        std::vector<std::pair<int, int>> normalizedEdges;
        normalizedEdges.reserve(defectEdges.size());
        std::set<std::pair<int, int>> seenEdges;
        for (const auto &[u, v] : defectEdges) {
            auto itU = pivotPos.find(u);
            auto itV = pivotPos.find(v);
            if (itU == pivotPos.end() || itV == pivotPos.end()) continue;
            int a = itU->second;
            int b = itV->second;
            if (a == b) continue;
            if (a > b) std::swap(a, b);
            if (!seenEdges.emplace(a, b).second) continue;
            normalizedEdges.emplace_back(a, b);
            pivotIsCore[a] = 1;
            pivotIsCore[b] = 1;
        }

        std::vector<int> globalToCore(pivotCount, -1);
        for (int i = 0; i < pivotCount; ++i) {
            if (pivotIsCore[i]) {
                globalToCore[i] = static_cast<int>(out.corePivotVerts.size());
                out.corePivotVerts.push_back(pivotVerts[i]);
            } else {
                out.isolatedPivotVerts.push_back(pivotVerts[i]);
            }
        }

        const int coreCount = static_cast<int>(out.corePivotVerts.size());
#if R2_LAB_ENABLE_DEFECT_D2_ORBIT
        if (normalizedEdges.size() == 2) {
            std::vector<std::pair<int, int>> coreDefectEdges;
            coreDefectEdges.reserve(2);
            std::set<std::pair<int, int>> coreEdgeSet;
            for (const auto &[aGlobal, bGlobal] : normalizedEdges) {
                const int a = globalToCore[aGlobal];
                const int b = globalToCore[bGlobal];
                if (a < 0 || b < 0) continue;
                const auto edge = a < b ? std::pair<int, int>{a, b} : std::pair<int, int>{b, a};
                coreDefectEdges.push_back(edge);
                coreEdgeSet.insert(edge);
            }

            const long long isolatedCount = static_cast<long long>(out.isolatedPivotVerts.size());
            out.coreKeepPivotWeights.assign(coreCount, 0);
            out.corePairWithIsolatedWeights.assign(coreCount, 0);

            if (coreCount == 3) {
                std::array<int, 3> degree{0, 0, 0};
                for (const auto &[a, b] : coreDefectEdges) {
                    degree[a]++;
                    degree[b]++;
                }
                int center = -1;
                std::array<int, 2> leaves{-1, -1};
                int leafPos = 0;
                for (int i = 0; i < 3; ++i) {
                    if (degree[i] == 2) center = i;
                    else if (degree[i] == 1 && leafPos < 2) leaves[leafPos++] = i;
                }
                if (center >= 0 && leafPos == 2) {
                    out.cliqueCount =
                        chooseBounded(isolatedCount, needPivot) +
                        3.0 * chooseBounded(isolatedCount, needPivot - 1) +
                        chooseBounded(isolatedCount, needPivot - 2);
                    out.isolatedKeepPivotWeight =
                        chooseBounded(isolatedCount - 1, needPivot - 1) +
                        3.0 * chooseBounded(isolatedCount - 1, needPivot - 2) +
                        chooseBounded(isolatedCount - 1, needPivot - 3);
                    out.isolatedPairWeight =
                        chooseBounded(isolatedCount - 2, needPivot - 2) +
                        3.0 * chooseBounded(isolatedCount - 2, needPivot - 3) +
                        chooseBounded(isolatedCount - 2, needPivot - 4);

                    out.coreKeepPivotWeights[center] = chooseBounded(isolatedCount, needPivot - 1);
                    out.coreKeepPivotWeights[leaves[0]] =
                        chooseBounded(isolatedCount, needPivot - 1) +
                        chooseBounded(isolatedCount, needPivot - 2);
                    out.coreKeepPivotWeights[leaves[1]] = out.coreKeepPivotWeights[leaves[0]];

                    out.corePairWithIsolatedWeights[center] =
                        chooseBounded(isolatedCount - 1, needPivot - 2);
                    out.corePairWithIsolatedWeights[leaves[0]] =
                        chooseBounded(isolatedCount - 1, needPivot - 2) +
                        chooseBounded(isolatedCount - 1, needPivot - 3);
                    out.corePairWithIsolatedWeights[leaves[1]] =
                        out.corePairWithIsolatedWeights[leaves[0]];

                    const double leafLeafWeight = chooseBounded(isolatedCount, needPivot - 2);
                    if (leafLeafWeight > 0)
                        out.corePairWeights.emplace_back(leaves[0], leaves[1], leafLeafWeight);

                    out.supported = true;
                    return true;
                }
            } else if (coreCount == 4) {
                out.cliqueCount =
                    chooseBounded(isolatedCount, needPivot) +
                    4.0 * chooseBounded(isolatedCount, needPivot - 1) +
                    4.0 * chooseBounded(isolatedCount, needPivot - 2);
                out.isolatedKeepPivotWeight =
                    chooseBounded(isolatedCount - 1, needPivot - 1) +
                    4.0 * chooseBounded(isolatedCount - 1, needPivot - 2) +
                    4.0 * chooseBounded(isolatedCount - 1, needPivot - 3);
                out.isolatedPairWeight =
                    chooseBounded(isolatedCount - 2, needPivot - 2) +
                    4.0 * chooseBounded(isolatedCount - 2, needPivot - 3) +
                    4.0 * chooseBounded(isolatedCount - 2, needPivot - 4);

                const double coreKeepW =
                    chooseBounded(isolatedCount, needPivot - 1) +
                    2.0 * chooseBounded(isolatedCount, needPivot - 2);
                const double coreIsoW =
                    chooseBounded(isolatedCount - 1, needPivot - 2) +
                    2.0 * chooseBounded(isolatedCount - 1, needPivot - 3);
                for (int i = 0; i < 4; ++i) {
                    out.coreKeepPivotWeights[i] = coreKeepW;
                    out.corePairWithIsolatedWeights[i] = coreIsoW;
                }

                const double crossWeight = chooseBounded(isolatedCount, needPivot - 2);
                if (crossWeight > 0) {
                    for (int i = 0; i < 4; ++i) {
                        for (int j = i + 1; j < 4; ++j) {
                            if (coreEdgeSet.contains({i, j})) continue;
                            out.corePairWeights.emplace_back(i, j, crossWeight);
                        }
                    }
                }

                out.supported = true;
                return true;
            }
        }
#endif
        if (normalizedEdges.size() == 2) {
            std::vector<std::pair<int, int>> coreDefectEdges;
            coreDefectEdges.reserve(2);
            std::set<std::pair<int, int>> coreEdgeSet;
            for (const auto &[aGlobal, bGlobal] : normalizedEdges) {
                const int a = globalToCore[aGlobal];
                const int b = globalToCore[bGlobal];
                if (a < 0 || b < 0) continue;
                const auto edge = a < b ? std::pair<int, int>{a, b} : std::pair<int, int>{b, a};
                coreDefectEdges.push_back(edge);
                coreEdgeSet.insert(edge);
            }

            out.coreKeepPivotWeights.assign(coreCount, 0);
            out.corePairWithIsolatedWeights.assign(coreCount, 0);
            std::vector<double> accumCorePairWeights(static_cast<size_t>(coreCount * coreCount), 0);

            const long long isolatedCount = static_cast<long long>(out.isolatedPivotVerts.size());
            const uint32_t fullMask = coreCount == 0 ? 1u : (uint32_t(1) << coreCount);
            auto isValidMask = [&](uint32_t mask) {
                for (const auto &[a, b] : coreDefectEdges) {
                    if (((mask >> a) & 1u) && ((mask >> b) & 1u)) return false;
                }
                return true;
            };

            for (uint32_t mask = 0; mask < fullMask; ++mask) {
                if (!isValidMask(mask)) continue;
                const int coreTake = static_cast<int>(__builtin_popcount(mask));
                out.cliqueCount += chooseBounded(isolatedCount, needPivot - coreTake);
                out.isolatedKeepPivotWeight += chooseBounded(isolatedCount - 1, needPivot - coreTake - 1);
                out.isolatedPairWeight += chooseBounded(isolatedCount - 2, needPivot - coreTake - 2);

                for (int i = 0; i < coreCount; ++i) {
                    if (((mask >> i) & 1u) == 0) continue;
                    out.coreKeepPivotWeights[i] += chooseBounded(isolatedCount, needPivot - coreTake);
                    out.corePairWithIsolatedWeights[i] += chooseBounded(isolatedCount - 1, needPivot - coreTake - 1);
                }

                for (int i = 0; i < coreCount; ++i) {
                    if (((mask >> i) & 1u) == 0) continue;
                    for (int j = i + 1; j < coreCount; ++j) {
                        if (((mask >> j) & 1u) == 0) continue;
                        if (coreEdgeSet.contains({i, j})) continue;
                        const double w = chooseBounded(isolatedCount, needPivot - coreTake);
                        if (w > 0) accumCorePairWeights[static_cast<size_t>(i * coreCount + j)] += w;
                    }
                }
            }

            for (int i = 0; i < coreCount; ++i) {
                for (int j = i + 1; j < coreCount; ++j) {
                    const double w = accumCorePairWeights[static_cast<size_t>(i * coreCount + j)];
                    if (w > 0) out.corePairWeights.emplace_back(i, j, w);
                }
            }

            out.supported = true;
            return true;
        }
        if (coreCount > R2_LAB_DEFECT_MAX_CORE_VERTS) return false;

        std::vector<uint64_t> closedMasks(coreCount, 0);
        for (int i = 0; i < coreCount; ++i)
            closedMasks[i] = (uint64_t(1) << i);

        std::set<std::pair<int, int>> coreEdgeSet;
        for (const auto &[aGlobal, bGlobal] : normalizedEdges) {
            const int a = globalToCore[aGlobal];
            const int b = globalToCore[bGlobal];
            if (a < 0 || b < 0) continue;
            closedMasks[a] |= (uint64_t(1) << b);
            closedMasks[b] |= (uint64_t(1) << a);
            if (a < b) coreEdgeSet.emplace(a, b);
            else coreEdgeSet.emplace(b, a);
        }

        const long long isolatedCount = static_cast<long long>(out.isolatedPivotVerts.size());
        const uint64_t fullMask = coreCount == 64 ? ~uint64_t(0) : ((uint64_t(1) << coreCount) - 1);
        DefectPolynomialCounter counter(closedMasks, needPivot);
        const auto &basePoly = counter.solve(fullMask);
        out.cliqueCount = evalDefectPolynomial(basePoly, isolatedCount, needPivot);

        out.isolatedKeepPivotWeight = evalDefectPolynomial(basePoly, isolatedCount - 1, needPivot - 1);
        out.isolatedPairWeight = evalDefectPolynomial(basePoly, isolatedCount - 2, needPivot - 2);

        out.coreKeepPivotWeights.resize(coreCount, 0);
        out.corePairWithIsolatedWeights.resize(coreCount, 0);
        for (int i = 0; i < coreCount; ++i) {
            const auto &forcedPoly = counter.solve(fullMask & ~closedMasks[i]);
            out.coreKeepPivotWeights[i] = evalDefectPolynomial(forcedPoly, isolatedCount, needPivot - 1);
            out.corePairWithIsolatedWeights[i] = evalDefectPolynomial(forcedPoly, isolatedCount - 1, needPivot - 2);
        }

        for (int i = 0; i < coreCount; ++i) {
            for (int j = i + 1; j < coreCount; ++j) {
                if (coreEdgeSet.contains({i, j})) continue;
                const auto &forcedPairPoly = counter.solve(fullMask & ~(closedMasks[i] | closedMasks[j]));
                const double w = evalDefectPolynomial(forcedPairPoly, isolatedCount, needPivot - 2);
                if (w > 0)
                    out.corePairWeights.emplace_back(i, j, w);
            }
        }

        out.supported = true;
        return true;
    }
#endif

    struct InitResult {
        double *countingKE;
        // Path → Edge CSR
        std::vector<daf::Size> leafEdgeOff;
        std::vector<LeafEdgeEntry> leafEdgeData;
        std::vector<double> leafWKK, leafWPP, leafWKP;
    };

    struct DCLPOptions {
        const char *label;
        bool hybridPhase1;
        int hybridMinVertexLeafDegree;
        int defectMinLeafVerts;
        int defectMinPivots;
        int defectMinActiveEdges;
    };

    InitResult buildCSR(const DynamicGraph<TreeGraphNode> &treeGraph,
                            const Graph &edgeGraph,
                            const daf::CliqueSize k) {
        const daf::Size numEdges = edgeGraph.adj_list.size();
        const daf::Size numLeaves = treeGraph.adj_list.size();

        auto *countingE = new double[numEdges];
        memset(countingE, 0, numEdges * sizeof(double));
        std::vector<double> leafWKK(numLeaves, 0), leafWPP(numLeaves, 0), leafWKP(numLeaves, 0);

        // --- Pass 1: count entries per leaf ---
        std::vector<daf::Size> leafEdgeOff(numLeaves + 1, 0);
        daf::StaticVector<daf::Size> tPovit, tKeepC;
        for (daf::Size li = 0; li < numLeaves; ++li) {
            const auto &clique = treeGraph.adj_list[li];
            if (clique.size() < k) continue;
            tPovit.clear(); tKeepC.clear();
            for (const auto &node : clique) {
                if (node.isPivot) tPovit.push_back(node.v);
                else tKeepC.push_back(node.v);
            }
            int np = int(k) - int(tKeepC.size());
            daf::Size cnt = 0;
            if (np >= 0 && np <= int(tPovit.size()))
                cnt += tKeepC.size() * (tKeepC.size() - 1) / 2;
            if (np - 2 >= 0 && np - 2 <= int(tPovit.size()) - 2)
                cnt += tPovit.size() * (tPovit.size() - 1) / 2;
            if (np - 1 >= 0 && np - 1 <= int(tPovit.size()) - 1)
                cnt += tKeepC.size() * tPovit.size();
            leafEdgeOff[li + 1] = cnt;
        }
        for (daf::Size li = 0; li < numLeaves; ++li)
            leafEdgeOff[li + 1] += leafEdgeOff[li];
        size_t totalEntries = leafEdgeOff[numLeaves];
        bool buildCSRData = true;  // always build; caller decides to keep or free

        // --- Pass 2: fill leaf→edge CSR + compute countingE ---
        std::vector<LeafEdgeEntry> leafEdgeData(totalEntries);

        for (daf::Size li = 0; li < numLeaves; ++li) {
            const auto &clique = treeGraph.adj_list[li];
            if (clique.size() < k) continue;
            tPovit.clear(); tKeepC.clear();
            for (const auto &node : clique) {
                if (node.isPivot) tPovit.push_back(node.v);
                else tKeepC.push_back(node.v);
            }
            int needPivot = int(k) - int(tKeepC.size());
            daf::Size pos = leafEdgeOff[li];

            if (needPivot >= 0 && needPivot <= int(tPovit.size())) {
                double w = nCr[tPovit.size()][needPivot];
                leafWKK[li] = w;
                for (size_t i = 0; i < tKeepC.size(); ++i)
                    for (size_t j = i + 1; j < tKeepC.size(); ++j) {
                        auto eid = edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]);
                        countingE[eid] += w;
                        leafEdgeData[pos++] = {eid, KK, 0};
                    }
            }
            int needPP = needPivot - 2;
            if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
                double w = nCr[tPovit.size() - 2][needPP];
                leafWPP[li] = w;
                for (size_t i = 0; i < tPovit.size(); ++i)
                    for (size_t j = i + 1; j < tPovit.size(); ++j) {
                        auto eid = edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]);
                        countingE[eid] += w;
                        leafEdgeData[pos++] = {eid, PP, 0};
                    }
            }
            int needKP = needPivot - 1;
            if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
                double w = nCr[tPovit.size() - 1][needKP];
                leafWKP[li] = w;
                for (size_t i = 0; i < tKeepC.size(); ++i)
                    for (size_t j = 0; j < tPovit.size(); ++j) {
                        auto eid = edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]);
                        countingE[eid] += w;
                        leafEdgeData[pos++] = {eid, KP,
                                               static_cast<uint8_t>(tPovit[j] < tKeepC[i])};
                    }
            }
            assert(pos == leafEdgeOff[li + 1]);
        }
        tPovit.free(); tKeepC.free();


        return {countingE,
                std::move(leafEdgeOff), std::move(leafEdgeData),
                std::move(leafWKK), std::move(leafWPP), std::move(leafWKP)};
    }

    EdgeLeafTranspose buildEdgeLeafTranspose(
        const Graph &edgeGraph,
        const DynamicGraphSet<TreeGraphNode> &treeGraphV,
        int minVertexLeafDegree) {

        const daf::Size numEdges = edgeGraph.adj_list.size();
        const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;
        std::vector<uint8_t> edgeIndexed(numEdges, 0);
        std::vector<daf::Size> edgeLeafOff(numEdges + 1, 0);
        std::vector<daf::Size> vertexLeafDegree(numVertices, 0);
        std::vector<daf::Size> indexedEdgeIds;
        for (daf::Size v = 0; v < numVertices; ++v)
            vertexLeafDegree[v] = treeGraphV.adj_list[v].size();

        for (daf::Size edgeId = 0; edgeId < numEdges; ++edgeId) {
            auto [u, v] = edgeGraph.getEdgeById(edgeId);
            if (std::min(vertexLeafDegree[u], vertexLeafDegree[v]) < static_cast<daf::Size>(minVertexLeafDegree))
                continue;
            edgeIndexed[edgeId] = 1;
            indexedEdgeIds.push_back(edgeId);
        }
        const daf::Size indexedEdges = indexedEdgeIds.size();
        const long long numIndexed = static_cast<long long>(indexedEdgeIds.size());

#pragma omp parallel for if(numIndexed > 4096) schedule(dynamic, 1024)
        for (long long idx = 0; idx < numIndexed; ++idx) {
            const daf::Size eid = indexedEdgeIds[static_cast<size_t>(idx)];
            auto [u, v] = edgeGraph.getEdgeById(eid);
            auto &adjU = treeGraphV.adj_list[u];
            auto &adjV = treeGraphV.adj_list[v];
            daf::Size localCount = 0;
            daf::intersect_dense_sets(adjU, adjV,
                [&](const TreeGraphNode &, const TreeGraphNode &) {
                    localCount++;
                });
            edgeLeafOff[eid + 1] = localCount;
        }
        for (daf::Size eid = 0; eid < numEdges; ++eid) {
            edgeLeafOff[eid + 1] += edgeLeafOff[eid];
        }

        std::vector<EdgeLeafEntry> edgeLeafDataT(edgeLeafOff[numEdges]);
        if (indexedEdges == 0)
            return {std::move(edgeLeafOff), std::move(edgeLeafDataT), std::move(edgeIndexed), indexedEdges};

#pragma omp parallel for if(numIndexed > 4096) schedule(dynamic, 1024)
        for (long long idx = 0; idx < numIndexed; ++idx) {
            const daf::Size edgeId = indexedEdgeIds[static_cast<size_t>(idx)];
            auto [u, v] = edgeGraph.getEdgeById(edgeId);
            auto &adjU = treeGraphV.adj_list[u];
            auto &adjV = treeGraphV.adj_list[v];
            daf::Size cursor = edgeLeafOff[edgeId];
            daf::intersect_dense_sets(adjU, adjV,
                [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                    EdgeType type = DCLP::KP;
                    uint8_t kpPivotIsLow = 0;
                    if (!uClique.isPivot && !vClique.isPivot) {
                        type = DCLP::KK;
                    } else if (uClique.isPivot && vClique.isPivot) {
                        type = DCLP::PP;
                    } else {
                        const daf::Size pivotV = uClique.isPivot ? u : v;
                        const daf::Size keepV = uClique.isPivot ? v : u;
                        kpPivotIsLow = static_cast<uint8_t>(pivotV < keepV);
                    }
                    edgeLeafDataT[cursor++] = {uClique.v, type, kpPivotIsLow};
                });
        }

        return {std::move(edgeLeafOff), std::move(edgeLeafDataT), std::move(edgeIndexed), indexedEdges};
    }
}

static std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_DCLP_Impl(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k,
    const DCLP::DCLPOptions &options) {
    auto time_start = std::chrono::high_resolution_clock::now();

    // Build leaf→edge CSR + countingKE
    auto initResult = DCLP::buildCSR(tree, edgeGraph, k);
    auto *countingKE = initResult.countingKE;
    auto &leafEdgeOff = initResult.leafEdgeOff;
    auto &leafEdgeData = initResult.leafEdgeData;
    auto &leafWKK = initResult.leafWKK;
    auto &leafWPP = initResult.leafWPP;
    auto &leafWKP = initResult.leafWKP;
    DCLP::EdgeLeafTranspose edgeLeafTranspose;
    if (options.hybridPhase1)
        edgeLeafTranspose = DCLP::buildEdgeLeafTranspose(edgeGraph, treeGraphV, options.hybridMinVertexLeafDegree);

    auto time_init = std::chrono::high_resolution_clock::now();
    std::cout << options.label << " init took: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_init - time_start).count()
              << " ms, leafEdgeData=" << leafEdgeData.size();
    if (options.hybridPhase1) {
        std::cout << ", indexedEdges=" << edgeLeafTranspose.indexedEdges
                  << ", indexedPairs=" << edgeLeafTranspose.edgeLeafData.size()
                  << ", minVertexLeafDegree=" << options.hybridMinVertexLeafDegree;
    }
    std::cout << std::endl;

    const daf::Size initialNumLeaves = leafEdgeOff.size() > 0 ? leafEdgeOff.size() - 1 : 0;
    // leafEdgeAlive: whether CSR data is valid for this leaf
    std::vector<uint8_t> leafEdgeAlive(initialNumLeaves, 1);
    std::vector<uint8_t> leafStaticTransposeAlive;
    if (options.hybridPhase1)
        leafStaticTransposeAlive.assign(initialNumLeaves, 1);

    const daf::Size numEdgesInit = edgeGraph.adj_list.size();
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;

    auto *coreE = new double[numEdgesInit];
    memset(coreE, 0, numEdgesInit * sizeof(double));

    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keepC;
    daf::StaticVector<daf::Size> newPivot;
    daf::StaticVector<daf::Size> newKeepC;

    daf::StaticVector<daf::Size> currentRemoveEdgeIds(edgeGraph.adj_list.size());

    daf::StaticVector<uint8_t> edgeInHeap(edgeGraph.adj_list.size());
    edgeInHeap.c_size = edgeGraph.adj_list.size();
    memset(edgeInHeap.getData(), 1, edgeGraph.adj_list.size() * sizeof(uint8_t));

    daf::StaticVector<daf::Size> removedLeaf(tree.adj_list.size());
    daf::StaticVector<DCLP::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();
    std::vector<uint8_t> leafInDynamicPhase1(tree.adj_list.size(), 0);
    DynamicGraphSet<TreeGraphNode> dynamicTreeGraphV;
    if (options.hybridPhase1)
        dynamicTreeGraphV.adj_list.resize(numVertices);

    double currCore = 0;

    // --- Pure d-ary heap PQ (same as REF) ---
    const daf::Size numEdges = edgeGraph.adj_list.size();
    struct CmpEdge {
        const double *cnt;
        bool operator()(daf::Size a, daf::Size b) const { return cnt[a] > cnt[b]; }
    };
    using EdgeHeap = boost::heap::d_ary_heap<daf::Size, boost::heap::arity<8>,
        boost::heap::mutable_<true>, boost::heap::compare<CmpEdge>>;
    EdgeHeap pq{CmpEdge{countingKE}};
    pq.reserve(numEdges);
    std::vector<EdgeHeap::handle_type> pqHandles(numEdges);

    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numEdges; ++i) {
        if (countingKE[i] == 0) { edgeInHeap[i] = 0; continue; }
        pqHandles[i] = pq.push(i);
        remainingInHeap++;
    }

    auto bucketMove = [&](daf::Size id) {
        if (!edgeInHeap[id]) return;
        pq.update(pqHandles[id]);
    };

    // --- deltaAccum: dirty edge tracking ---
    std::vector<uint8_t> edgeDirty(numEdges, 0);
    std::vector<daf::Size> dirtyEdges;
    dirtyEdges.reserve(4096);

    auto markDirty = [&](daf::Size id) {
        if (!edgeDirty[id]) {
            edgeDirty[id] = 1;
            dirtyEdges.push_back(id);
        }
    };

    std::vector<uint8_t> leafAffected(leafRmInfo.size(), 0);
    std::vector<daf::Size> caseBLeafIds;
    caseBLeafIds.reserve(1024);
#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
    std::vector<daf::Size> defectLeafIds;
    defectLeafIds.reserve(1024);
    std::vector<std::vector<std::pair<daf::Size, daf::Size>>> leafDefectEdges(tree.adj_list.size());
    std::vector<uint8_t> leafDefectProbeCount(tree.adj_list.size(), 0);
#endif
#ifdef HARD_LEAF_LAB_PROFILE
    std::vector<uint8_t> profSeenCaseB(tree.adj_list.size(), 0);
    std::vector<uint8_t> profSeenCaseCAlive(tree.adj_list.size(), 0);
    std::vector<uint8_t> profSeenHard(tree.adj_list.size(), 0);
#endif

    auto ensureLeafTracking = [&](daf::Size leafId) {
        if (leafId >= leafRmInfo.size()) {
            const auto newCap = static_cast<size_t>(leafId * 3 / 2 + 8);
            removedLeaf.reserve(newCap);
            leafRmInfo.resize(newCap);
            leafAffected.resize(newCap, 0);
            leafInDynamicPhase1.resize(newCap, 0);
#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
            leafDefectEdges.resize(newCap);
            leafDefectProbeCount.resize(newCap, 0);
#endif
#ifdef HARD_LEAF_LAB_PROFILE
            profSeenCaseB.resize(newCap, 0);
            profSeenCaseCAlive.resize(newCap, 0);
            profSeenHard.resize(newCap, 0);
#endif
        } else if (leafId >= leafInDynamicPhase1.size()) {
            leafInDynamicPhase1.resize(leafId + 1, 0);
#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
            if (leafId >= leafDefectEdges.size())
                leafDefectEdges.resize(leafId + 1);
            if (leafId >= leafDefectProbeCount.size())
                leafDefectProbeCount.resize(leafId + 1, 0);
#endif
        }
    };

    auto removeDynamicLeafMembership = [&](daf::Size leafId) {
        if (!options.hybridPhase1) return;
        if (leafId >= leafInDynamicPhase1.size() || !leafInDynamicPhase1[leafId]) return;
        for (const auto &leafV : tree.adj_list[leafId])
            dynamicTreeGraphV.removeNbr(leafV.v, {leafId, leafV.isPivot});
        leafInDynamicPhase1[leafId] = 0;
    };

    auto addDynamicLeafMembership = [&](daf::Size leafId) {
        if (!options.hybridPhase1) return;
        ensureLeafTracking(leafId);
        for (const auto &leafV : tree.adj_list[leafId])
            dynamicTreeGraphV.addNbr(leafV.v, {leafId, leafV.isPivot});
        leafInDynamicPhase1[leafId] = 1;
        if (leafId < leafStaticTransposeAlive.size())
            leafStaticTransposeAlive[leafId] = 0;
    };

    auto invalidateStaticLeaf = [&](daf::Size leafId) {
        if (options.hybridPhase1 && leafId < leafStaticTransposeAlive.size())
            leafStaticTransposeAlive[leafId] = 0;
    };

    auto applyPlainLeafContribution = [&](const std::vector<TreeGraphNode> &leafNodes, auto &&upd) {
        newPivot.clear();
        newKeepC.clear();
        for (const auto &node : leafNodes) {
            if (node.isPivot) newPivot.push_back(node.v);
            else newKeepC.push_back(node.v);
        }
        const daf::Size need = k - newKeepC.size();
        if (need <= newPivot.size() && newKeepC.size() > 1) {
            const double w = nCr[newPivot.size()][need];
            DCLP::processEdgePairs(newKeepC, w, upd);
        }
        const int needPP = int(need) - 2;
        if (0 <= needPP && needPP <= int(newPivot.size()) - 2) {
            const double w = nCr[newPivot.size() - 2][needPP];
            DCLP::processEdgePairs(newPivot, w, upd);
        }
        const int needKP = int(need) - 1;
        if (0 <= needKP && needKP <= int(newPivot.size()) - 1) {
            const double w = nCr[newPivot.size() - 1][needKP];
            DCLP::processEdgePairs(newKeepC, newPivot, w, upd);
        }
        newPivot.clear();
        newKeepC.clear();
    };

#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
    auto applyDefectContribution = [&](const DCLP::DefectContributionWeights &weights,
                                       const daf::StaticVector<daf::Size> &keepVerts,
                                       auto &&upd) {
        if (weights.cliqueCount <= 0) return;
        if (weights.cliqueCount > 0 && keepVerts.size() > 1)
            DCLP::processEdgePairs(keepVerts, weights.cliqueCount, upd);

        if (weights.isolatedKeepPivotWeight > 0 && !weights.isolatedPivotVerts.empty() && !keepVerts.empty()) {
            for (const auto keepV : keepVerts)
                for (const auto isoV : weights.isolatedPivotVerts)
                    upd(keepV, isoV, weights.isolatedKeepPivotWeight);
        }

        for (size_t i = 0; i < weights.corePivotVerts.size(); ++i) {
            const double wKP = weights.coreKeepPivotWeights[i];
            if (wKP <= 0) continue;
            for (const auto keepV : keepVerts)
                upd(keepV, weights.corePivotVerts[i], wKP);
        }

        if (weights.isolatedPairWeight > 0 && weights.isolatedPivotVerts.size() > 1)
            DCLP::processEdgePairs(weights.isolatedPivotVerts, weights.isolatedPairWeight, upd);

        for (size_t i = 0; i < weights.corePivotVerts.size(); ++i) {
            const double wPI = weights.corePairWithIsolatedWeights[i];
            if (wPI <= 0) continue;
            for (const auto isoV : weights.isolatedPivotVerts)
                upd(weights.corePivotVerts[i], isoV, wPI);
        }

        for (const auto &[a, b, w] : weights.corePairWeights) {
            if (w <= 0) continue;
            upd(weights.corePivotVerts[a], weights.corePivotVerts[b], w);
        }
    };
#endif

    auto markAffectedLeaf = [&](daf::Size leafId, bool removedKeepC, bool removedBothPivots,
                                daf::Size edgeU, daf::Size edgeV, bool removedUPivot) {
        ensureLeafTracking(leafId);
        auto &lrm = leafRmInfo[leafId];
        if (lrm.empty()) {
            removedLeaf.push_back(leafId);
            leafAffected[leafId] = 1;
        }
        if (lrm.removedKeepC) return;
        if (removedKeepC) {
            lrm.removedKeepC = true;
            return;
        }
        if (removedBothPivots) {
            lrm.removedEdges.push_back({edgeU, edgeV});
            return;
        }
        lrm.removedPivots.push_back(removedUPivot ? edgeU : edgeV);
    };

    std::cout << "=========================begin=========================" << std::endl;
    double minCore = 0;
    long long cntA = 0, cntB = 0, cntC = 0, totalIters = 0;
    long long cntA_fast = 0, cntA_fallback = 0;
    long long cntC_csr = 0, cntC_fallback = 0;
    long long us_phase1 = 0, us_caseAC_delta = 0;
    long long us_caseB = 0, us_flush = 0;
    long long work_p1_alive = 0;
    long long work_caseA_csr = 0, work_caseA_fallback = 0;
    long long work_caseC_csr = 0, work_caseC_fallback = 0;
    long long cntB_closedForm = 0, cntB_bk = 0;
#if R2_LAB_ENABLE_CASEB_LEAF_REUSE
    long long cntB_reuseSingleLeaf = 0, cntB_reuseSplitLeaf = 0;
#endif
#if R2_LAB_ENABLE_CASEB_ZEROEDGE_CSR
    long long cntB_zeroEdgeCsr = 0;
    long long work_caseB_zeroEdgeCsr = 0;
#endif
#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
    long long cntDefect = 0, cntDefectFallbackBK = 0, cntDefectPersistent = 0;
    long long cntDefectCand2 = 0, cntDefectCand3 = 0, cntDefectCand4Plus = 0;
    long long cntDefectCand2Core3 = 0, cntDefectCand2Core4 = 0;
#endif
#ifdef HARD_LEAF_LAB_PROFILE
    long long prof_caseA_occ = 0, prof_caseB_occ = 0, prof_caseC_alive_occ = 0, prof_caseC_dead_occ = 0;
    long long prof_caseA_us = 0, prof_caseB_us = 0, prof_caseC_alive_us = 0, prof_caseC_dead_us = 0;
    long long prof_caseB_closedForm_us = 0, prof_caseB_bk_us = 0;
    long long prof_caseC_csr_us = 0, prof_caseC_fallback_us = 0;
    long long prof_caseB_activeRemovedEdges = 0, prof_caseC_removedPivots = 0;
    long long prof_caseB_leafVerts = 0, prof_caseC_alive_leafVerts = 0;
    daf::Size prof_caseB_maxLeafVerts = 0, prof_caseC_alive_maxLeafVerts = 0;
    long long prof_caseBUnique = 0, prof_caseCAliveUnique = 0, prof_hardUnique = 0;
#endif

    while (remainingInHeap > 0) {
        // --- Heap pop: batch all edges at current min level ---
        while (!pq.empty() && !edgeInHeap[pq.top()]) pq.pop();
        if (pq.empty()) break;

        minCore = std::max(countingKE[pq.top()], minCore);
        while (!pq.empty() && countingKE[pq.top()] <= minCore) {
            auto id = pq.top(); pq.pop();
            if (!edgeInHeap[id]) continue;
            edgeInHeap[id] = 0;
            currentRemoveEdgeIds.push_back(id);
            coreE[id] = minCore;
            remainingInHeap--;
        }

        currCore = minCore;

        if (remainingInHeap == 0) break;

        totalIters++;

        // ========== Phase 1: treeGraphV intersection (same as V8b) ==========
        auto _t0 = std::chrono::high_resolution_clock::now();
        for (int ei = 0; ei < (int)currentRemoveEdgeIds.size(); ++ei) {
            auto edgeId = currentRemoveEdgeIds[ei];
            auto [edgeU, edgeV] = edgeGraph.getEdgeById(edgeId);
            if (options.hybridPhase1 && edgeId < edgeLeafTranspose.edgeIndexed.size() &&
                edgeLeafTranspose.edgeIndexed[edgeId]) {
                const auto begin = edgeLeafTranspose.edgeLeafOff[edgeId];
                const auto end = edgeLeafTranspose.edgeLeafOff[edgeId + 1];
                for (daf::Size pos = begin; pos < end; ++pos) {
                    const auto &entry = edgeLeafTranspose.edgeLeafData[pos];
                    if (entry.leafId >= leafStaticTransposeAlive.size() ||
                        !leafStaticTransposeAlive[entry.leafId])
                        continue;
                    work_p1_alive++;
                    if (entry.type == DCLP::KK) {
                        markAffectedLeaf(entry.leafId, true, false, edgeU, edgeV, false);
                    } else if (entry.type == DCLP::PP) {
                        markAffectedLeaf(entry.leafId, false, true, edgeU, edgeV, false);
                    } else {
                        markAffectedLeaf(entry.leafId, false, false, edgeU, edgeV, entry.kpPivotIsLow != 0);
                    }
                }

                auto &adjU = dynamicTreeGraphV.getNbr(edgeU);
                auto &adjV = dynamicTreeGraphV.getNbr(edgeV);
                daf::intersect_dense_sets(adjU, adjV,
                    [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                        work_p1_alive++;
                        markAffectedLeaf(uClique.v,
                                         !uClique.isPivot && !vClique.isPivot,
                                         uClique.isPivot && vClique.isPivot,
                                         edgeU, edgeV, uClique.isPivot);
                    });
            }
            if (!options.hybridPhase1 || edgeId >= edgeLeafTranspose.edgeIndexed.size() ||
                !edgeLeafTranspose.edgeIndexed[edgeId]) {
                auto &adjU = treeGraphV.getNbr(edgeU);
                auto &adjV = treeGraphV.getNbr(edgeV);
                daf::intersect_dense_sets(adjU, adjV,
                    [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                        work_p1_alive++;
                        markAffectedLeaf(uClique.v,
                                         !uClique.isPivot && !vClique.isPivot,
                                         uClique.isPivot && vClique.isPivot,
                                         edgeU, edgeV, uClique.isPivot);
                    });
            }
        }

        // Pre-sort removedPivots
        auto _t1 = std::chrono::high_resolution_clock::now();
        us_phase1 += std::chrono::duration_cast<std::chrono::microseconds>(_t1 - _t0).count();
        for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
            auto &rp = leafRmInfo[removedLeaf[leafIdIdx]].removedPivots;
            if (rp.size() == 2) {
                if (rp[0] > rp[1]) std::swap(rp[0], rp[1]);
                if (rp[0] == rp[1]) rp.c_size = 1;
            } else if (rp.size() > 2) {
                std::ranges::sort(rp);
                rp.unique();
            }
        }

        // ========== Phase 2: MERGED delta + tree mutation for Case A & C ==========
        auto _t2 = std::chrono::high_resolution_clock::now();
        caseBLeafIds.clear();
        {
            daf::StaticVector<daf::Size> tPovit, tKeepC;

            auto directSub = [&](daf::Size idx, double w) {
                countingKE[idx] -= w;
                if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx])
                    pq.update(pqHandles[idx]);
            };

            for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                auto leafId = removedLeaf[leafIdIdx];
                DCLP::LeafRmInfo &leafRm = leafRmInfo[leafId];

                const auto& leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;
#ifdef HARD_LEAF_LAB_PROFILE
                auto leafProfileStart = std::chrono::high_resolution_clock::now();
#endif

                tPovit.clear(); tKeepC.clear();
                daf::Size numKeeps = 0;
                for (const auto& node : leaf) {
                    if (node.isPivot) tPovit.push_back(node.v);
                    else { tKeepC.push_back(node.v); numKeeps++; }
                }
                daf::Size needPivot = k - numKeeps;
                daf::Size numPivots = tPovit.size();

#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
                const bool hasPersistentDefect = leafId < leafDefectEdges.size() && !leafDefectEdges[leafId].empty();
                const bool survivesPivotRemovals = !leafRm.removedKeepC &&
                    needPivot <= numPivots - leafRm.removedPivots.size();
                const bool meetsDefectShapeThreshold =
                    survivesPivotRemovals &&
                    !leafRm.removedEdges.empty() &&
                    static_cast<int>(leaf.size()) >= options.defectMinLeafVerts &&
                    static_cast<int>(numPivots) >= options.defectMinPivots;
                int activeDefectEdges = 0;
                if (!hasPersistentDefect && meetsDefectShapeThreshold && options.defectMinActiveEdges > 0) {
                    for (const auto &[eu, ev] : leafRm.removedEdges) {
                        const bool euRemoved = std::binary_search(leafRm.removedPivots.begin(),
                                                                  leafRm.removedPivots.end(), eu);
                        const bool evRemoved = std::binary_search(leafRm.removedPivots.begin(),
                                                                  leafRm.removedPivots.end(), ev);
                        if (!euRemoved && !evRemoved) {
                            activeDefectEdges++;
                            if (activeDefectEdges >= options.defectMinActiveEdges) break;
                        }
                    }
                }
                const bool complexNewDefect =
                    meetsDefectShapeThreshold &&
                    activeDefectEdges >= options.defectMinActiveEdges;
                if (complexNewDefect) {
                    if (activeDefectEdges == 2) {
                        cntDefectCand2++;
                        std::array<daf::Size, 4> activeVerts{};
                        int activeVertCount = 0;
                        auto addActiveVert = [&](daf::Size v) {
                            for (int i = 0; i < activeVertCount; ++i) {
                                if (activeVerts[i] == v) return;
                            }
                            if (activeVertCount < 4) activeVerts[activeVertCount++] = v;
                        };
                        for (const auto &[eu, ev] : leafRm.removedEdges) {
                            const bool euRemoved = std::binary_search(leafRm.removedPivots.begin(),
                                                                      leafRm.removedPivots.end(), eu);
                            const bool evRemoved = std::binary_search(leafRm.removedPivots.begin(),
                                                                      leafRm.removedPivots.end(), ev);
                            if (!euRemoved && !evRemoved) {
                                addActiveVert(eu);
                                addActiveVert(ev);
                            }
                        }
                        if (activeVertCount == 3) cntDefectCand2Core3++;
                        else if (activeVertCount == 4) cntDefectCand2Core4++;
                    } else if (activeDefectEdges == 3) {
                        cntDefectCand3++;
                    } else {
                        cntDefectCand4Plus++;
                    }
                }
                bool eligibleNewDefect = false;
                if (complexNewDefect) {
                    if (R2_LAB_ENABLE_DEFECT_D2_IMMEDIATE && activeDefectEdges == 2) {
                        eligibleNewDefect = true;
                    } else {
                        auto &probeCount = leafDefectProbeCount[leafId];
                        if (probeCount > 0) {
                            eligibleNewDefect = true;
                        } else {
                            probeCount = 1;
                        }
                    }
                }
                if (hasPersistentDefect || eligibleNewDefect) {
                    defectLeafIds.push_back(leafId);
                    if (hasPersistentDefect) cntDefectPersistent++;
                    continue;
                }
#endif

                bool isDeadLeaf = leafRm.removedKeepC || needPivot > numPivots - leafRm.removedPivots.size();
                bool isCaseB = !leafRm.removedEdges.empty() && !isDeadLeaf;

                if (isCaseB) {
                    caseBLeafIds.push_back(leafId);
                    continue;
                }

                if (isDeadLeaf) {
                    // ---- Case A: CSR scan (same as V8b) ----
                    cntA++;

                    if (leafId < leafEdgeAlive.size() && leafEdgeAlive[leafId]) {
                        cntA_fast++;
                        double wKK = leafWKK[leafId];
                        double wPP = leafWPP[leafId];
                        double wKP = leafWKP[leafId];
                        daf::Size begin = leafEdgeOff[leafId];
                        daf::Size end = leafEdgeOff[leafId + 1];
                        work_caseA_csr += (end - begin);
                        for (daf::Size pos = begin; pos < end; ++pos) {
                            auto &entry = leafEdgeData[pos];
                            double w;
                            switch (entry.type) {
                                case DCLP::KK: w = wKK; break;
                                case DCLP::PP: w = wPP; break;
                                case DCLP::KP: w = wKP; break;
                            }
                            directSub(entry.edgeId, w);
                        }
                        leafEdgeAlive[leafId] = 0;
                    } else {
                        // Fallback for Case B overflow leaves
                        cntA_fallback++;
                        double KtoK = 0, KtoP = 0, PtoP = 0;
                        if (needPivot <= tPovit.size()) {
                            KtoK = nCr[tPovit.size()][needPivot];
                            for (daf::Size i = 0; i + 1 < tKeepC.size(); ++i)
                                for (daf::Size j = i + 1; j < tKeepC.size(); ++j) {
                                    work_caseA_fallback++;
                                    directSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]), KtoK);
                                }
                        }
                        int needPP = int(needPivot) - 2;
                        if (0 <= needPP && needPP <= int(tPovit.size()) - 2) {
                            PtoP = nCr[tPovit.size() - 2][needPP];
                            for (daf::Size i = 0; i + 1 < tPovit.size(); ++i)
                                for (daf::Size j = i + 1; j < tPovit.size(); ++j) {
                                    work_caseA_fallback++;
                                    directSub(edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]), PtoP);
                                }
                        }
                        int needKP = int(needPivot) - 1;
                        if (0 <= needKP && needKP <= int(tPovit.size()) - 1) {
                            KtoP = nCr[tPovit.size() - 1][needKP];
                            for (daf::Size i = 0; i < tKeepC.size(); ++i)
                                for (daf::Size j = 0; j < tPovit.size(); ++j) {
                                    work_caseA_fallback++;
                                    directSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]), KtoP);
                                }
                        }
                    }
                    // Tree mutation
                    if (options.hybridPhase1) {
                        removeDynamicLeafMembership(leafId);
                        invalidateStaticLeaf(leafId);
                    }
                    for (const auto& i : leaf)
                        treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                    tree.adj_list[leafId].clear();
#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
                    leafDefectProbeCount[leafId] = 0;
#endif
#ifdef HARD_LEAF_LAB_PROFILE
                    prof_caseA_occ++;
                    prof_caseA_us += std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - leafProfileStart).count();
#endif
                } else {
                    // ---- Case C: DCLP CSR-based delta ----
                    cntC++;
                    double KtoK = 0, KtoP = 0, PtoP = 0;
                    double RemovedKtoK = 0, RemovedKtoP = 0, RemovedPtoP = 0;

                    if (needPivot <= tPovit.size()) {
                        KtoK = nCr[tPovit.size()][needPivot] - nCr[tPovit.size() - leafRm.removedPivots.size()][needPivot];
                        RemovedKtoK = nCr[tPovit.size()][needPivot];
                    }
                    int needPP = int(needPivot) - 2;
                    if (0 <= needPP && needPP <= int(tPovit.size()) - 2) {
                        RemovedPtoP = nCr[tPovit.size() - 2][needPP];
                        PtoP = RemovedPtoP;
                        if (leafRm.removedPivots.size() + 2 <= tPovit.size())
                            PtoP -= nCr[tPovit.size() - 2 - leafRm.removedPivots.size()][needPP];
                    }
                    int needKP = int(needPivot) - 1;
                    if (0 <= needKP && needKP <= int(tPovit.size()) - 1) {
                        RemovedKtoP = nCr[tPovit.size() - 1][needKP];
                        KtoP = RemovedKtoP;
                        if (leafRm.removedPivots.size() + 1 <= tPovit.size())
                            KtoP -= nCr[tPovit.size() - 1 - leafRm.removedPivots.size()][needKP];
                    }

                    const bool isCaseCAlive = !leafRm.removedPivots.empty() &&
                                              needPivot <= tPovit.size() - leafRm.removedPivots.size();

                    if (isCaseCAlive) {
                        // Case C alive: leaf survives with fewer pivots
                        bool caseCCsrFastPath = false;
#ifdef HARD_LEAF_LAB_PROFILE
                        prof_caseC_alive_occ++;
                        prof_caseC_removedPivots += static_cast<long long>(leafRm.removedPivots.size());
                        prof_caseC_alive_leafVerts += static_cast<long long>(leaf.size());
                        prof_caseC_alive_maxLeafVerts = std::max<daf::Size>(prof_caseC_alive_maxLeafVerts, leaf.size());
                        if (!profSeenCaseCAlive[leafId]) {
                            profSeenCaseCAlive[leafId] = 1;
                            prof_caseCAliveUnique++;
                        }
                        if (!profSeenHard[leafId]) {
                            profSeenHard[leafId] = 1;
                            prof_hardUnique++;
                        }
#endif
                        if (leafId < leafEdgeAlive.size() && leafEdgeAlive[leafId]) {
                            // ---- DCLP: CSR scan for Case C ----
                            cntC_csr++;
                            caseCCsrFastPath = true;
                            daf::Size begin = leafEdgeOff[leafId];
                            daf::Size end = leafEdgeOff[leafId + 1];
                            work_caseC_csr += (end - begin);

                            for (daf::Size pos = begin; pos < end; ++pos) {
                                auto &entry = leafEdgeData[pos];
                                auto [eu, ev] = edgeGraph.getEdgeById(entry.edgeId);
                                bool eu_rm = std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), eu);
                                bool ev_rm = std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), ev);
                                double delta;

                                if (entry.type == DCLP::KK) {
                                    delta = KtoK;
                                } else if (entry.type == DCLP::PP) {
                                    if (eu_rm || ev_rm) delta = RemovedPtoP;
                                    else delta = PtoP;
                                } else { // KP
                                    // Determine which is the pivot
                                    bool pivotRemoved;
                                    // Check if eu is a pivot (in tPovit)
                                    if (std::find(tPovit.begin(), tPovit.end(), eu) != tPovit.end())
                                        pivotRemoved = eu_rm;
                                    else
                                        pivotRemoved = ev_rm;

                                    if (pivotRemoved) delta = RemovedKtoP;
                                    else delta = KtoP;
                                }
                                directSub(entry.edgeId, delta);
                            }
                            leafEdgeAlive[leafId] = 0;  // CSR weights now stale
                        } else {
                            // Fallback: same as V8b baseline Case C
                            cntC_fallback++;
                            daf::StaticVector<TreeGraphNode> newLeafF;
                            for (const auto& node : leaf) {
                                if (!std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), node.v))
                                    newLeafF.push_back(node);
                            }
                            daf::Size p = leafRm.removedPivots.size();
                            work_caseC_fallback += p*(p-1)/2 + newLeafF.size()*p + newLeafF.size()*(newLeafF.size()-1)/2;
                            for (daf::Size i = 0; i + 1 < leafRm.removedPivots.size(); ++i)
                                for (daf::Size j = i + 1; j < leafRm.removedPivots.size(); ++j)
                                    directSub(edgeGraph.getEdgeCompressedId(leafRm.removedPivots[i], leafRm.removedPivots[j]), RemovedPtoP);
                            for (daf::Size i = 0; i < newLeafF.size(); ++i)
                                for (daf::Size j = 0; j < leafRm.removedPivots.size(); ++j) {
                                    double d = newLeafF[i].isPivot ? RemovedPtoP : RemovedKtoP;
                                    directSub(edgeGraph.getEdgeCompressedId(newLeafF[i].v, leafRm.removedPivots[j]), d);
                                }
                            for (daf::Size i = 0; i + 1 < newLeafF.size(); ++i)
                                for (daf::Size j = i + 1; j < newLeafF.size(); ++j) {
                                    auto &u = newLeafF[i], &v = newLeafF[j];
                                    double d = (!u.isPivot && !v.isPivot) ? KtoK : (u.isPivot && v.isPivot) ? PtoP : KtoP;
                                    directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                                }
                            newLeafF.free();
                        }
                        // Tree mutation for Case C
                        if (options.hybridPhase1) {
                            removeDynamicLeafMembership(leafId);
                            invalidateStaticLeaf(leafId);
                        }
                        for (auto removedNbr : leafRm.removedPivots)
                            treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                        tree.removeNbrs(leafId, leafRm.removedPivots);
                        if (options.hybridPhase1)
                            addDynamicLeafMembership(leafId);
#ifdef HARD_LEAF_LAB_PROFILE
                        const auto caseCAliveUs = std::chrono::duration_cast<std::chrono::microseconds>(
                            std::chrono::high_resolution_clock::now() - leafProfileStart).count();
                        prof_caseC_alive_us += caseCAliveUs;
                        if (caseCCsrFastPath) prof_caseC_csr_us += caseCAliveUs;
                        else prof_caseC_fallback_us += caseCAliveUs;
#endif
                    } else {
                        // Case C dead: leaf dies (same as Case A fallback)
                        if (leafId < leafEdgeAlive.size() && leafEdgeAlive[leafId]) {
                            // Use CSR for dead Case C
                            double wKK = leafWKK[leafId];
                            double wPP = leafWPP[leafId];
                            double wKP = leafWKP[leafId];
                            daf::Size begin = leafEdgeOff[leafId];
                            daf::Size end = leafEdgeOff[leafId + 1];
                            for (daf::Size pos = begin; pos < end; ++pos) {
                                auto &entry = leafEdgeData[pos];
                                double w;
                                switch (entry.type) {
                                    case DCLP::KK: w = wKK; break;
                                    case DCLP::PP: w = wPP; break;
                                    case DCLP::KP: w = wKP; break;
                                }
                                directSub(entry.edgeId, w);
                            }
                            leafEdgeAlive[leafId] = 0;
                        } else {
                            for (daf::Size i = 0; i + 1 < leaf.size(); ++i)
                                for (daf::Size j = i + 1; j < leaf.size(); ++j) {
                                    auto &u = leaf[i], &v = leaf[j];
                                    double d = (!u.isPivot && !v.isPivot) ? RemovedKtoK : (u.isPivot && v.isPivot) ? RemovedPtoP : RemovedKtoP;
                                    directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                                }
                        }
                        if (options.hybridPhase1) {
                            removeDynamicLeafMembership(leafId);
                            invalidateStaticLeaf(leafId);
                        }
                        for (const auto& i : leaf)
                            treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                        tree.adj_list[leafId].clear();
#ifdef HARD_LEAF_LAB_PROFILE
                        prof_caseC_dead_occ++;
                        prof_caseC_dead_us += std::chrono::duration_cast<std::chrono::microseconds>(
                            std::chrono::high_resolution_clock::now() - leafProfileStart).count();
#endif
                    }
                }
                leafRmInfo[leafId].clear();
            }
            tPovit.free(); tKeepC.free();
        }

        // ========== Phase 2c: Case B — closed-form d=1, BK fallback d≥2 ==========
        auto _t3 = std::chrono::high_resolution_clock::now();
        us_caseAC_delta += std::chrono::duration_cast<std::chrono::microseconds>(_t3 - _t2).count();
#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
        cntDefect += defectLeafIds.size();
        for (int di = 0; di < (int)defectLeafIds.size(); ++di) {
            const auto leafId = defectLeafIds[di];
            DCLP::LeafRmInfo &leafRm = leafRmInfo[leafId];
            auto &leaf = tree.adj_list[leafId];
            auto &persistentDefects = leafDefectEdges[leafId];
            if (leaf.empty()) {
                persistentDefects.clear();
                leafDefectProbeCount[leafId] = 0;
                leafRm.clear();
                continue;
            }

            auto addW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] += w;
                if (edgeInHeap[idx])
                    pq.update(pqHandles[idx]);
            };
            auto removeW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] -= w;
                if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx])
                    pq.update(pqHandles[idx]);
            };

            auto applyDefectOrEnumeratedContribution = [&](const std::vector<TreeGraphNode> &leafNodes,
                                                           const std::vector<std::pair<daf::Size, daf::Size>> &defectEdges,
                                                           auto &&upd) -> bool {
                newKeepC.clear();
                newPivot.clear();
                for (const auto &node : leafNodes) {
                    if (node.isPivot) newPivot.push_back(node.v);
                    else newKeepC.push_back(node.v);
                }

                if (defectEdges.empty()) {
                    applyPlainLeafContribution(leafNodes, upd);
                    newKeepC.clear();
                    newPivot.clear();
                    return true;
                }

                DCLP::DefectContributionWeights weights;
                if (DCLP::buildDefectContributionWeights(newKeepC, newPivot, defectEdges, k, weights)) {
                    applyDefectContribution(weights, newKeepC, upd);
                    newKeepC.clear();
                    newPivot.clear();
                    return true;
                }

                daf::StaticVector<std::pair<daf::Size, daf::Size>> defectStatic;
                defectStatic.reserve(defectEdges.size());
                for (const auto &edge : defectEdges) defectStatic.push_back(edge);
                bkRmEdge::bronKerbosch(const_cast<std::vector<TreeGraphNode>&>(leafNodes), defectStatic, k,
                    [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                        auto subLeaf = bkRmEdge::coverToVertex(c, pivots, leafNodes);
                        applyPlainLeafContribution(subLeaf, upd);
                    });
                defectStatic.free();
                newKeepC.clear();
                newPivot.clear();
                return false;
            };

            applyDefectOrEnumeratedContribution(leaf, persistentDefects, removeW);

            if (options.hybridPhase1) {
                removeDynamicLeafMembership(leafId);
                invalidateStaticLeaf(leafId);
            }
            for (const auto &leafV : leaf)
                treeGraphV.removeNbr(leafV.v, {leafId, leafV.isPivot});
            if (!leafRm.removedPivots.empty())
                tree.removeNbrs(leafId, leafRm.removedPivots);
            if (leafId < leafEdgeAlive.size())
                leafEdgeAlive[leafId] = 0;

            auto &leafRef = tree.adj_list[leafId];
            auto isAliveVertex = [&](daf::Size v) {
                for (const auto &node : leafRef)
                    if (node.v == v) return true;
                return false;
            };

            std::vector<std::pair<daf::Size, daf::Size>> mergedDefects;
            mergedDefects.reserve(persistentDefects.size() + leafRm.removedEdges.size());
            auto pushActiveDefect = [&](daf::Size u, daf::Size v) {
                if (!isAliveVertex(u) || !isAliveVertex(v) || u == v) return;
                if (u > v) std::swap(u, v);
                mergedDefects.emplace_back(u, v);
            };
            for (const auto &[u, v] : persistentDefects) pushActiveDefect(u, v);
            for (const auto &[u, v] : leafRm.removedEdges) pushActiveDefect(u, v);
            std::sort(mergedDefects.begin(), mergedDefects.end());
            mergedDefects.erase(std::unique(mergedDefects.begin(), mergedDefects.end()), mergedDefects.end());

            newKeepC.clear();
            newPivot.clear();
            for (const auto &node : leafRef) {
                if (node.isPivot) newPivot.push_back(node.v);
                else newKeepC.push_back(node.v);
            }
            const daf::Size newNeedPivot = k - newKeepC.size();
            bool killLeaf = leafRm.removedKeepC || newNeedPivot > newPivot.size();

            if (!killLeaf) {
                DCLP::DefectContributionWeights weights;
                bool supported = DCLP::buildDefectContributionWeights(newKeepC, newPivot, mergedDefects, k, weights);
                if (supported) {
                    if (weights.cliqueCount <= 0) {
                        killLeaf = true;
                    } else {
                        for (const auto &leafV : leafRef)
                            treeGraphV.addNbr(leafV.v, {leafId, leafV.isPivot});
                        if (options.hybridPhase1) addDynamicLeafMembership(leafId);
                        applyDefectContribution(weights, newKeepC, addW);
                        persistentDefects = std::move(mergedDefects);
                        leafDefectProbeCount[leafId] = 0;
                    }
                } else {
                    cntDefectFallbackBK++;
                    daf::StaticVector<std::pair<daf::Size, daf::Size>> defectStatic;
                    defectStatic.reserve(mergedDefects.size());
                    for (const auto &edge : mergedDefects) defectStatic.push_back(edge);
                    bkRmEdge::bronKerbosch(leafRef, defectStatic, k,
                        [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                            auto subLeaf = bkRmEdge::coverToVertex(c, pivots, leafRef);
                            auto newId = tree.addNode(subLeaf);
                            ensureLeafTracking(newId);
                            for (const auto &node : tree.adj_list[newId])
                                treeGraphV.addNbr(node.v, {newId, node.isPivot});
                            if (options.hybridPhase1) addDynamicLeafMembership(newId);
                            applyPlainLeafContribution(tree.adj_list[newId], addW);
                        });
                    defectStatic.free();
                    persistentDefects.clear();
                    leafDefectProbeCount[leafId] = 0;
                    killLeaf = true;
                }
            }

            if (killLeaf) {
                leafRef.clear();
                persistentDefects.clear();
                leafDefectProbeCount[leafId] = 0;
            }
            newKeepC.clear();
            newPivot.clear();
            leafRm.clear();
        }
        defectLeafIds.clear();
#endif
        cntB += caseBLeafIds.size();

        for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
            auto leafId = caseBLeafIds[bi];
            DCLP::LeafRmInfo &leafRm = leafRmInfo[leafId];
            const auto& leaf = tree.adj_list[leafId];
#ifdef HARD_LEAF_LAB_PROFILE
            auto caseBLeafStart = std::chrono::high_resolution_clock::now();
            prof_caseB_occ++;
            prof_caseB_leafVerts += static_cast<long long>(leaf.size());
            prof_caseB_maxLeafVerts = std::max<daf::Size>(prof_caseB_maxLeafVerts, leaf.size());
            if (!profSeenCaseB[leafId]) {
                profSeenCaseB[leafId] = 1;
                prof_caseBUnique++;
            }
            if (!profSeenHard[leafId]) {
                profSeenHard[leafId] = 1;
                prof_hardUnique++;
            }
#endif

            povit.clear(); keepC.clear();
            for (const auto& node : leaf) {
                if (node.isPivot) povit.push_back(node.v);
                else keepC.push_back(node.v);
            }
            daf::Size needPivot = k - keepC.size();

            auto addW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] += w;
                if (edgeInHeap[idx])
                    pq.update(pqHandles[idx]);
            };
            auto directSubCaseB = [&](daf::Size idx, double w) {
                countingKE[idx] -= w;
                if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx])
                    pq.update(pqHandles[idx]);
            };

            // Remove old leaf from treeGraphV
            if (options.hybridPhase1) {
                removeDynamicLeafMembership(leafId);
                invalidateStaticLeaf(leafId);
            }
            for (const auto& leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }
            if (!leafRm.removedPivots.empty()) {
                tree.removeNbrs(leafId, leafRm.removedPivots);
            }

            // Invalidate CSR
            const bool hadLeafEdgeCsr = leafId < leafEdgeAlive.size() && leafEdgeAlive[leafId];
            if (leafId < leafEdgeAlive.size())
                leafEdgeAlive[leafId] = 0;

            auto &leafRef = tree.adj_list[leafId];
            bool keepLeafIdAlive = false;
            bool skipCaseBOldContributionRemoval = false;

            // Filter removedEdges: keep only those where both vertices still in leaf
            daf::StaticVector<std::pair<daf::Size,daf::Size>> activeRmEdges;
            activeRmEdges.reserve(leafRm.removedEdges.size());
            for (auto &[eu, ev] : leafRm.removedEdges) {
                bool foundU = false, foundV = false;
                for (const auto &n : leafRef) {
                    if (n.v == eu) foundU = true;
                    if (n.v == ev) foundV = true;
                }
                if (foundU && foundV) activeRmEdges.push_back({eu, ev});
            }
#ifdef HARD_LEAF_LAB_PROFILE
            const daf::Size activeRmEdgesCount = activeRmEdges.size();
            prof_caseB_activeRemovedEdges += static_cast<long long>(activeRmEdgesCount);
#endif

            auto activateLeafNode = [&](daf::Size targetId, bool addContribution = true) {
                for (const auto& i : tree.adj_list[targetId]) treeGraphV.addNbr(i.v, {targetId, i.isPivot});
                if (options.hybridPhase1) addDynamicLeafMembership(targetId);
                if (addContribution) {
                    newPivot.clear(); newKeepC.clear();
                    for (const auto& i : tree.adj_list[targetId]) {
                        if (i.isPivot) newPivot.push_back(i.v);
                        else newKeepC.push_back(i.v);
                    }
                    daf::Size np = k - newKeepC.size();
                    if (np <= newPivot.size() && newKeepC.size() > 1) {
                        double w = nCr[newPivot.size()][np];
                        DCLP::processEdgePairs(newKeepC, w, addW);
                    }
                    int nPP = int(np) - 2;
                    if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                        double w = nCr[newPivot.size() - 2][nPP];
                        DCLP::processEdgePairs(newPivot, w, addW);
                    }
                    int nKP = int(np) - 1;
                    if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                        double w = nCr[newPivot.size() - 1][nKP];
                        DCLP::processEdgePairs(newKeepC, newPivot, w, addW);
                    }
                    newPivot.clear(); newKeepC.clear();
                }
                ensureLeafTracking(targetId);
            };

            // Helper: create a sub-leaf by excluding one vertex
            auto createSubLeaf = [&](daf::Size excludeV) {
                std::vector<TreeGraphNode> sub;
                sub.reserve(leafRef.size() - 1);
                for (const auto &n : leafRef) {
                    if (n.v != excludeV) sub.push_back(n);
                }
                auto newId = tree.addNode(sub);
                activateLeafNode(newId);
            };

            if (activeRmEdges.size() == 0) {
                // All PP edges involved removed pivots — no split needed
                // The leaf (minus removed pivots) is still one clique
                // Re-add it as a "sub-leaf" = the current modified leaf
#if R2_LAB_ENABLE_CASEB_ZEROEDGE_CSR
                if (hadLeafEdgeCsr) {
                    double KtoK = 0, KtoP = 0, PtoP = 0;
                    double RemovedKtoP = 0, RemovedPtoP = 0;

                    if (needPivot <= povit.size()) {
                        KtoK = nCr[povit.size()][needPivot] -
                               nCr[povit.size() - leafRm.removedPivots.size()][needPivot];
                    }
                    int nPP = int(needPivot) - 2;
                    if (0 <= nPP && nPP <= int(povit.size()) - 2) {
                        RemovedPtoP = nCr[povit.size() - 2][nPP];
                        PtoP = RemovedPtoP;
                        if (leafRm.removedPivots.size() + 2 <= povit.size())
                            PtoP -= nCr[povit.size() - 2 - leafRm.removedPivots.size()][nPP];
                    }
                    int nKP = int(needPivot) - 1;
                    if (0 <= nKP && nKP <= int(povit.size()) - 1) {
                        RemovedKtoP = nCr[povit.size() - 1][nKP];
                        KtoP = RemovedKtoP;
                        if (leafRm.removedPivots.size() + 1 <= povit.size())
                            KtoP -= nCr[povit.size() - 1 - leafRm.removedPivots.size()][nKP];
                    }

                    daf::Size begin = leafEdgeOff[leafId];
                    daf::Size end = leafEdgeOff[leafId + 1];
                    work_caseB_zeroEdgeCsr += (end - begin);
                    for (daf::Size pos = begin; pos < end; ++pos) {
                        auto &entry = leafEdgeData[pos];
                        auto [eu, ev] = edgeGraph.getEdgeById(entry.edgeId);
                        bool eu_rm = std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), eu);
                        bool ev_rm = std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), ev);
                        double delta;

                        if (entry.type == DCLP::KK) {
                            delta = KtoK;
                        } else if (entry.type == DCLP::PP) {
                            delta = (eu_rm || ev_rm) ? RemovedPtoP : PtoP;
                        } else {
                            bool pivotRemoved;
                            if (std::find(povit.begin(), povit.end(), eu) != povit.end())
                                pivotRemoved = eu_rm;
                            else
                                pivotRemoved = ev_rm;
                            delta = pivotRemoved ? RemovedKtoP : KtoP;
                        }
                        directSubCaseB(entry.edgeId, delta);
                    }
#if R2_LAB_ENABLE_CASEB_LEAF_REUSE
                    activateLeafNode(leafId, false);
                    keepLeafIdAlive = true;
                    cntB_reuseSingleLeaf++;
#else
                    std::vector<TreeGraphNode> sub(leafRef.begin(), leafRef.end());
                    auto newId = tree.addNode(sub);
                    activateLeafNode(newId, false);
#endif
                    skipCaseBOldContributionRemoval = true;
                    cntB_zeroEdgeCsr++;
                    cntB_closedForm++;
                    activeRmEdges.free();
                    leafRmInfo[leafId].clear();
                    povit.clear(); keepC.clear();
#ifdef HARD_LEAF_LAB_PROFILE
                    const auto caseBLeafUs = std::chrono::duration_cast<std::chrono::microseconds>(
                        std::chrono::high_resolution_clock::now() - caseBLeafStart).count();
                    prof_caseB_us += caseBLeafUs;
                    prof_caseB_closedForm_us += caseBLeafUs;
#endif
                    continue;
                }
#endif
#if R2_LAB_ENABLE_CASEB_LEAF_REUSE
                activateLeafNode(leafId);
                keepLeafIdAlive = true;
                cntB_reuseSingleLeaf++;
#else
                std::vector<TreeGraphNode> sub(leafRef.begin(), leafRef.end());
                auto newId = tree.addNode(sub);
                activateLeafNode(newId);
#endif
                cntB_closedForm++;
            } else if (activeRmEdges.size() == 1) {
                // ---- CLOSED-FORM Case B: d=1 ----
                // Correct partition (no overlap in s-cliques):
                //   Sub-leaf 1 = L \ {rp2}  (rp1 stays as pivot)
                //   Sub-leaf 2 = L \ {rp1}, with rp2 changed to KEEP
                // Sub-leaf 1 covers cliques NOT using rp2
                // Sub-leaf 2 covers cliques USING rp2 but NOT rp1
                auto [rp1, rp2] = activeRmEdges[0];

                // Sub-leaf 1: remove rp2 entirely
#if R2_LAB_ENABLE_CASEB_LEAF_REUSE
                std::vector<TreeGraphNode> sub2;
                sub2.reserve(leafRef.size() - 1);
                for (const auto &n : leafRef) {
                    if (n.v == rp1) continue;
                    if (n.v == rp2) sub2.push_back({rp2, false}); // pivot → keep
                    else sub2.push_back(n);
                }
                {
                    std::vector<TreeGraphNode> sub1;
                    sub1.reserve(leafRef.size() - 1);
                    for (const auto &n : leafRef) {
                        if (n.v != rp2) sub1.push_back(n);
                    }
                    leafRef = std::move(sub1);
                }
                activateLeafNode(leafId);
                keepLeafIdAlive = true;
                cntB_reuseSplitLeaf++;
#else
                createSubLeaf(rp2);
                std::vector<TreeGraphNode> sub2;
                sub2.reserve(leafRef.size() - 1);
                for (const auto &n : leafRef) {
                    if (n.v == rp1) continue;
                    if (n.v == rp2) sub2.push_back({rp2, false}); // pivot → keep
                    else sub2.push_back(n);
                }
#endif

                // Sub-leaf 2: remove rp1, change rp2 from pivot to keep
                {
                    auto newId = tree.addNode(sub2);
                    activateLeafNode(newId);
                }
                cntB_closedForm++;
            } else {
                // ---- BK fallback: d≥2 ----
                bkRmEdge::bronKerbosch(leafRef, leafRm.removedEdges, k,
                    [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                        auto subLeaf = bkRmEdge::coverToVertex(c, pivots, leafRef);
                        auto newId = tree.addNode(subLeaf);
                        newPivot.clear(); newKeepC.clear();
                        for (const auto& i : tree.adj_list[newId]) {
                            if (i.isPivot) newPivot.push_back(i.v);
                            else newKeepC.push_back(i.v);
                        }
                        for (const auto& i : tree.adj_list[newId]) treeGraphV.addNbr(i.v, {newId, i.isPivot});
                        if (options.hybridPhase1) addDynamicLeafMembership(newId);
                        daf::Size np = k - newKeepC.size();
                        if (np <= newPivot.size() && newKeepC.size() > 1) {
                            double w = nCr[newPivot.size()][np];
                            DCLP::processEdgePairs(newKeepC, w, addW);
                        }
                        int nPP = int(np) - 2;
                        if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                            double w = nCr[newPivot.size() - 2][nPP];
                            DCLP::processEdgePairs(newPivot, w, addW);
                        }
                        int nKP = int(np) - 1;
                        if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                            double w = nCr[newPivot.size() - 1][nKP];
                            DCLP::processEdgePairs(newKeepC, newPivot, w, addW);
                        }
                        newPivot.clear(); newKeepC.clear();
                        ensureLeafTracking(newId);
                    });
                cntB_bk++;
            }
            activeRmEdges.free();

            // Remove old contribution
            if (!skipCaseBOldContributionRemoval) {
                auto removeW = [&](daf::Size u, daf::Size v, double w) {
                    auto idx = edgeGraph.getEdgeCompressedId(u, v);
                    countingKE[idx] -= w;
                    if (countingKE[idx] < 0) countingKE[idx] = 0;
                    if (edgeInHeap[idx])
                        pq.update(pqHandles[idx]);
                };
                if (needPivot <= povit.size()) {
                    double w = nCr[povit.size()][needPivot];
                    DCLP::processEdgePairs(keepC, w, removeW);
                }
                int needPP = int(needPivot) - 2;
                if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                    double w = nCr[povit.size() - 2][needPP];
                    DCLP::processEdgePairs(povit, w, removeW);
                }
                int needKP = int(needPivot) - 1;
                if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                    double w = nCr[povit.size() - 1][needKP];
                    DCLP::processEdgePairs(keepC, povit, w, removeW);
                }
            }

            if (!keepLeafIdAlive)
                tree.adj_list[leafId].clear();
#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
            if (!keepLeafIdAlive) leafDefectProbeCount[leafId] = 0;
#endif
            leafRmInfo[leafId].clear();
            povit.clear(); keepC.clear();
#ifdef HARD_LEAF_LAB_PROFILE
            const auto caseBLeafUs = std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - caseBLeafStart).count();
            prof_caseB_us += caseBLeafUs;
            if (activeRmEdgesCount > 1) prof_caseB_bk_us += caseBLeafUs;
            else prof_caseB_closedForm_us += caseBLeafUs;
#endif
        }

        // ========== Flush dirty edges ==========
        auto _t4 = std::chrono::high_resolution_clock::now();
        us_caseB += std::chrono::duration_cast<std::chrono::microseconds>(_t4 - _t3).count();

        for (auto eid : dirtyEdges) {
            edgeDirty[eid] = 0;
            if (edgeInHeap[eid])
                bucketMove(eid);
        }
        dirtyEdges.clear();

        auto _t5 = std::chrono::high_resolution_clock::now();
        us_flush += std::chrono::duration_cast<std::chrono::microseconds>(_t5 - _t4).count();

        currentRemoveEdgeIds.clear();
        for (auto leafId : removedLeaf) leafAffected[leafId] = 0;
        removedLeaf.clear();
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
    std::cout << "  Cases: A=" << cntA << " B=" << cntB << " C=" << cntC << " iters=" << totalIters << std::endl;
    std::cout << "  CaseA: csr=" << cntA_fast << " fallback=" << cntA_fallback << std::endl;
    std::cout << "  CaseC: csr=" << cntC_csr << " fallback=" << cntC_fallback << std::endl;
#if R2_LAB_ENABLE_CASEB_DEFECT_GRAPH
    std::cout << "  Defect: leaves=" << cntDefect
              << " persistentHits=" << cntDefectPersistent
              << " fallbackBK=" << cntDefectFallbackBK << std::endl;
    std::cout << "  DefectCand: d2=" << cntDefectCand2
              << " d3=" << cntDefectCand3
              << " d4p=" << cntDefectCand4Plus << std::endl;
    std::cout << "  DefectD2Shape: core3=" << cntDefectCand2Core3
              << " core4=" << cntDefectCand2Core4 << std::endl;
#endif
    std::cout << "  Phase1(edgePathCSR): " << us_phase1/1000 << " ms" << std::endl;
    std::cout << "  Phase2(A+C delta+mut): " << us_caseAC_delta/1000 << " ms" << std::endl;
    std::cout << "  CaseB: closedForm=" << cntB_closedForm << " bk=" << cntB_bk << std::endl;
#if R2_LAB_ENABLE_CASEB_LEAF_REUSE
    std::cout << "  CaseB reuse: singleLeaf=" << cntB_reuseSingleLeaf
              << " splitLeaf=" << cntB_reuseSplitLeaf << std::endl;
#endif
#if R2_LAB_ENABLE_CASEB_ZEROEDGE_CSR
    std::cout << "  CaseB zeroEdgeCSR: " << cntB_zeroEdgeCsr << std::endl;
#endif
    std::cout << "  CaseB(total): " << us_caseB/1000 << " ms" << std::endl;
    std::cout << "  Flush(bucketMove): " << us_flush/1000 << " ms" << std::endl;
    std::cout << "  P1: alive=" << work_p1_alive << std::endl;
    std::cout << "  CaseA_CSR=" << work_caseA_csr << " CaseA_fallback=" << work_caseA_fallback << std::endl;
    std::cout << "  CaseC_CSR=" << work_caseC_csr << " CaseC_fallback=" << work_caseC_fallback << std::endl;
#if R2_LAB_ENABLE_CASEB_ZEROEDGE_CSR
    std::cout << "  CaseB_zeroEdgeCSR_work=" << work_caseB_zeroEdgeCsr << std::endl;
#endif
#ifdef HARD_LEAF_LAB_PROFILE
    const long long prof_total_case_occ = prof_caseA_occ + prof_caseB_occ + prof_caseC_alive_occ + prof_caseC_dead_occ;
    const long long prof_hard_occ = prof_caseB_occ + prof_caseC_alive_occ;
    const long long prof_total_case_us = prof_caseA_us + prof_caseB_us + prof_caseC_alive_us + prof_caseC_dead_us;
    const long long prof_hard_us = prof_caseB_us + prof_caseC_alive_us;
    auto ratioPct = [](long long x, long long y) -> double {
        return y == 0 ? 0.0 : 100.0 * static_cast<double>(x) / static_cast<double>(y);
    };
    auto avgOrZero = [](long long x, long long y) -> double {
        return y == 0 ? 0.0 : static_cast<double>(x) / static_cast<double>(y);
    };
    std::cout << "  HardLeaf candidates: occ=" << prof_hard_occ << "/" << prof_total_case_occ
              << " (" << ratioPct(prof_hard_occ, prof_total_case_occ) << "%)"
              << ", unique=" << prof_hardUnique << std::endl;
    std::cout << "  HardLeaf time: " << prof_hard_us/1000 << " ms / " << prof_total_case_us/1000
              << " ms (" << ratioPct(prof_hard_us, prof_total_case_us) << "%)" << std::endl;
    std::cout << "  CaseB profile: occ=" << prof_caseB_occ
              << ", unique=" << prof_caseBUnique
              << ", avgLeafVerts=" << avgOrZero(prof_caseB_leafVerts, prof_caseB_occ)
              << ", maxLeafVerts=" << prof_caseB_maxLeafVerts
              << ", avgActiveRmEdges=" << avgOrZero(prof_caseB_activeRemovedEdges, prof_caseB_occ)
              << ", time=" << prof_caseB_us/1000 << " ms"
              << " [closed=" << prof_caseB_closedForm_us/1000
              << " ms, bk=" << prof_caseB_bk_us/1000 << " ms]" << std::endl;
    std::cout << "  CaseC-alive profile: occ=" << prof_caseC_alive_occ
              << ", unique=" << prof_caseCAliveUnique
              << ", avgLeafVerts=" << avgOrZero(prof_caseC_alive_leafVerts, prof_caseC_alive_occ)
              << ", maxLeafVerts=" << prof_caseC_alive_maxLeafVerts
              << ", avgRemovedPivots=" << avgOrZero(prof_caseC_removedPivots, prof_caseC_alive_occ)
              << ", time=" << prof_caseC_alive_us/1000 << " ms"
              << " [csr=" << prof_caseC_csr_us/1000
              << " ms, fallback=" << prof_caseC_fallback_us/1000 << " ms]" << std::endl;
    std::cout << "  Easy leaves: CaseA occ=" << prof_caseA_occ
              << ", time=" << prof_caseA_us/1000
              << " ms; CaseC-dead occ=" << prof_caseC_dead_occ
              << ", time=" << prof_caseC_dead_us/1000 << " ms" << std::endl;
#endif

    // Build output
    std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> sortedK;
    sortedK.reserve(edgeGraph.adj_list.size());
    const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;
    for (daf::Size u = 0; u < n; ++u) {
        const daf::Size start = edgeGraph.adj_list_offsets[u];
        const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
        for (daf::Size idx = start; idx < end; ++idx) {
            sortedK.emplace_back(
                std::make_pair(std::make_pair(u, edgeGraph.adj_list[idx]), (int)coreE[idx]));
        }
    }

    delete[] countingKE;
    delete[] coreE;
    povit.free(); keepC.free();
    newPivot.free(); newKeepC.free();
    return sortedK;
}

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DCLP(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    static const DCLP::DCLPOptions options{
        R2_LAB_DCLP_LABEL,
        false,
        0,
        0,
        0,
        0
    };
    return PlusNucleusEdgeCoreDecompositionSet_DCLP_Impl(tree, edgeGraph, treeGraphV, k, options);
}

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_Hybrid(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    int minVertexLeafDegree = 320;
    if (const char *env = std::getenv("PIVOTER_R2_HYBRID_LAB_MIN_VTX_LEAF")) {
        minVertexLeafDegree = std::max(0, std::atoi(env));
    } else if (const char *env = std::getenv("PIVOTER_R2_HYBRID_MIN_VTX_LEAF")) {
        minVertexLeafDegree = std::max(0, std::atoi(env));
    }
    const DCLP::DCLPOptions options{
        R2_LAB_HYBRID_LABEL,
        true,
        minVertexLeafDegree,
        0,
        0,
        0
    };
    return PlusNucleusEdgeCoreDecompositionSet_DCLP_Impl(tree, edgeGraph, treeGraphV, k, options);
}

} // namespace R2_LAB_NAMESPACE

namespace {

struct HybridLabAutoStats {
    int minVertexLeafDegree = 320;
    daf::Size estimatedIndexedEdges = 0;
    daf::Size sampledEdges = 0;
    double sampledAvgCommonLeaves = 0.0;
    double estimatedIndexedPairs = 0.0;
    double estimatedPairDensity = 0.0;
    bool shouldUseHybrid = false;
};

int getEnvIntOr(const char *name, int fallback) {
    if (const char *env = std::getenv(name))
        return std::max(0, std::atoi(env));
    return fallback;
}

double getEnvDoubleOr(const char *name, double fallback) {
    if (const char *env = std::getenv(name))
        return std::max(0.0, std::atof(env));
    return fallback;
}

HybridLabAutoStats estimateHybridLabAutoStats(
    const Graph &edgeGraph,
    const DynamicGraphSet<TreeGraphNode> &treeGraphV,
    int minVertexLeafDegree,
    int sampleBudget,
    int guardMinIndexedEdges,
    double minAvgCommonLeaves,
    double minPairDensity,
    double maxPairDensity) {
    HybridLabAutoStats stats;
    stats.minVertexLeafDegree = minVertexLeafDegree;

    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;
    const daf::Size numEdges = edgeGraph.adj_list.size();
    std::vector<daf::Size> vertexLeafDegree(numVertices, 0);
    for (daf::Size v = 0; v < numVertices; ++v)
        vertexLeafDegree[v] = treeGraphV.adj_list[v].size();

    for (daf::Size edgeId = 0; edgeId < numEdges; ++edgeId) {
        auto [u, v] = edgeGraph.getEdgeById(edgeId);
        if (std::min(vertexLeafDegree[u], vertexLeafDegree[v]) >=
            static_cast<daf::Size>(minVertexLeafDegree)) {
            stats.estimatedIndexedEdges++;
        }
    }

    if (stats.estimatedIndexedEdges == 0)
        return stats;
    if (stats.estimatedIndexedEdges < static_cast<daf::Size>(guardMinIndexedEdges))
        return stats;

    const daf::Size sampleTarget = std::min<daf::Size>(
        static_cast<daf::Size>(std::max(1, sampleBudget)),
        stats.estimatedIndexedEdges);
    const daf::Size stride = std::max<daf::Size>(1, stats.estimatedIndexedEdges / sampleTarget);

    daf::Size indexedSeen = 0;
    double sampledCommonLeaves = 0.0;
    for (daf::Size edgeId = 0; edgeId < numEdges; ++edgeId) {
        auto [u, v] = edgeGraph.getEdgeById(edgeId);
        if (std::min(vertexLeafDegree[u], vertexLeafDegree[v]) <
            static_cast<daf::Size>(minVertexLeafDegree)) {
            continue;
        }

        if ((indexedSeen % stride) != 0 && stats.sampledEdges + 1 < sampleTarget) {
            indexedSeen++;
            continue;
        }

        auto &adjU = treeGraphV.adj_list[u];
        auto &adjV = treeGraphV.adj_list[v];
        daf::Size localCount = 0;
        daf::intersect_dense_sets(adjU, adjV,
            [&](const TreeGraphNode &, const TreeGraphNode &) {
                localCount++;
            });
        sampledCommonLeaves += static_cast<double>(localCount);
        stats.sampledEdges++;
        indexedSeen++;

        if (stats.sampledEdges >= sampleTarget)
            break;
    }

    if (stats.sampledEdges == 0)
        return stats;

    stats.sampledAvgCommonLeaves = sampledCommonLeaves / static_cast<double>(stats.sampledEdges);
    stats.estimatedIndexedPairs =
        stats.sampledAvgCommonLeaves * static_cast<double>(stats.estimatedIndexedEdges);
    stats.estimatedPairDensity = numEdges == 0
        ? 0.0
        : stats.estimatedIndexedPairs / static_cast<double>(numEdges);
    stats.shouldUseHybrid =
        stats.estimatedIndexedEdges >= static_cast<daf::Size>(guardMinIndexedEdges) &&
        stats.sampledAvgCommonLeaves >= minAvgCommonLeaves &&
        stats.estimatedPairDensity >= minPairDensity &&
        stats.estimatedPairDensity <= maxPairDensity;
    return stats;
}

} // namespace

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > R2_LAB_ENTRYPOINT(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    int minVertexLeafDegree = 320;
    bool explicitThreshold = false;
    if (const char *env = std::getenv("PIVOTER_R2_HYBRID_LAB_MIN_VTX_LEAF")) {
        minVertexLeafDegree = std::max(0, std::atoi(env));
        explicitThreshold = true;
    } else if (const char *env = std::getenv("PIVOTER_R2_HYBRID_MIN_VTX_LEAF")) {
        minVertexLeafDegree = std::max(0, std::atoi(env));
        explicitThreshold = true;
    }

    const int guardMinIndexedEdges = getEnvIntOr(
        "PIVOTER_R2_HYBRID_LAB_GUARD_MIN_INDEXED_EDGES", 50000);
    const int sampleBudget = getEnvIntOr(
        "PIVOTER_R2_HYBRID_LAB_AUTO_SAMPLE_EDGES", 512);
    const double minAvgCommonLeaves = getEnvDoubleOr(
        "PIVOTER_R2_HYBRID_LAB_AUTO_MIN_AVG_COMMON", 128.0);
    const double minPairDensity = getEnvDoubleOr(
        "PIVOTER_R2_HYBRID_LAB_AUTO_MIN_PAIR_DENSITY", 3.0);
    const double maxPairDensity = getEnvDoubleOr(
        "PIVOTER_R2_HYBRID_LAB_AUTO_MAX_PAIR_DENSITY", 8.0);
    const bool forceHybrid = getEnvIntOr("PIVOTER_R2_HYBRID_LAB_FORCE", 0) != 0;

    const auto autoStats = estimateHybridLabAutoStats(
        edgeGraph,
        treeGraphV,
        minVertexLeafDegree,
        sampleBudget,
        guardMinIndexedEdges,
        minAvgCommonLeaves,
        minPairDensity,
        maxPairDensity);

    std::cout << "Hybrid R2 Lab auto: threshold=" << autoStats.minVertexLeafDegree
              << ", explicitThreshold=" << (explicitThreshold ? 1 : 0)
              << ", estimatedIndexedEdges=" << autoStats.estimatedIndexedEdges
              << ", sampledEdges=" << autoStats.sampledEdges
              << ", avgCommonLeaves=" << autoStats.sampledAvgCommonLeaves
              << ", estimatedIndexedPairs=" << autoStats.estimatedIndexedPairs
              << ", pairDensity=" << autoStats.estimatedPairDensity
              << ", densityRange=[" << minPairDensity << "," << maxPairDensity << "]"
              << std::endl;

    if (!forceHybrid && !autoStats.shouldUseHybrid) {
        std::cout << "Hybrid R2 Lab auto: fallback to DCLP baseline" << std::endl;
        return R2_LAB_NAMESPACE::PlusNucleusEdgeCoreDecompositionSet_DCLP(tree, edgeGraph, treeGraphV, k);
    }

    return R2_LAB_NAMESPACE::PlusNucleusEdgeCoreDecompositionSet_Hybrid(tree, edgeGraph, treeGraphV, k);
}
