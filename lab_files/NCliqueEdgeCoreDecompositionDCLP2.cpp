//
// R2 DCLP2: DCLP + Complement BK for Case B d≥2
//
// Key innovation over DCLP:
//   When Case B has d≥2 removed edges (both-pivot edges), the standard BK runs on
//   the full N-vertex near-complete leaf (O(3^{N/3})). DCLP2 instead identifies the
//   "affected" vertices A (endpoints of removed edges, |A|≤2d) and runs BK on only
//   |A| vertices (O(3^{|A|/3})). Since every maximal clique must include all
//   unaffected vertices, we combine BK results with unaffected vertices to get the
//   same sub-leaves. For web-it-2004 (omega=430), this turns 430-vertex BK into
//   ~10-40 vertex BK — a massive speedup.
//
// Inherited from DCLP:
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
#include <boost/heap/d_ary_heap.hpp>

#include "../BK/BronKerboschRmEdge.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace DCLP2 {

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
        for (daf::Size v = 0; v < numVertices; ++v)
            vertexLeafDegree[v] = treeGraphV.adj_list[v].size();

        for (daf::Size edgeId = 0; edgeId < numEdges; ++edgeId) {
            auto [u, v] = edgeGraph.getEdgeById(edgeId);
            if (std::min(vertexLeafDegree[u], vertexLeafDegree[v]) < static_cast<daf::Size>(minVertexLeafDegree))
                continue;
            edgeIndexed[edgeId] = 1;
        }
        daf::Size indexedEdges = 0;
        for (daf::Size eid = 0; eid < numEdges; ++eid) {
            if (edgeIndexed[eid]) indexedEdges++;
            if (!edgeIndexed[eid]) continue;
            auto [u, v] = edgeGraph.getEdgeById(eid);
            auto &adjU = treeGraphV.adj_list[u];
            auto &adjV = treeGraphV.adj_list[v];
            daf::intersect_dense_sets(adjU, adjV,
                [&](const TreeGraphNode &, const TreeGraphNode &) {
                    edgeLeafOff[eid + 1]++;
                });
        }
        for (daf::Size eid = 0; eid < numEdges; ++eid) {
            edgeLeafOff[eid + 1] += edgeLeafOff[eid];
        }

        std::vector<EdgeLeafEntry> edgeLeafDataT(edgeLeafOff[numEdges]);
        std::vector<daf::Size> cursor = edgeLeafOff;
        for (daf::Size edgeId = 0; edgeId < numEdges; ++edgeId) {
            if (!edgeIndexed[edgeId]) continue;
            auto [u, v] = edgeGraph.getEdgeById(edgeId);
            auto &adjU = treeGraphV.adj_list[u];
            auto &adjV = treeGraphV.adj_list[v];
            daf::intersect_dense_sets(adjU, adjV,
                [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                    EdgeType type = DCLP2::KP;
                    uint8_t kpPivotIsLow = 0;
                    if (!uClique.isPivot && !vClique.isPivot) {
                        type = DCLP2::KK;
                    } else if (uClique.isPivot && vClique.isPivot) {
                        type = DCLP2::PP;
                    } else {
                        const daf::Size pivotV = uClique.isPivot ? u : v;
                        const daf::Size keepV = uClique.isPivot ? v : u;
                        kpPivotIsLow = static_cast<uint8_t>(pivotV < keepV);
                    }
                    edgeLeafDataT[cursor[edgeId]++] = {uClique.v, type, kpPivotIsLow};
                });
        }

        return {std::move(edgeLeafOff), std::move(edgeLeafDataT), std::move(edgeIndexed), indexedEdges};
    }
}

static std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_DCLP_Impl(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k,
    const DCLP2::DCLPOptions &options) {
    auto time_start = std::chrono::high_resolution_clock::now();

    // Build leaf→edge CSR + countingKE
    auto initResult = DCLP2::buildCSR(tree, edgeGraph, k);
    auto *countingKE = initResult.countingKE;
    auto &leafEdgeOff = initResult.leafEdgeOff;
    auto &leafEdgeData = initResult.leafEdgeData;
    auto &leafWKK = initResult.leafWKK;
    auto &leafWPP = initResult.leafWPP;
    auto &leafWKP = initResult.leafWKP;
    DCLP2::EdgeLeafTranspose edgeLeafTranspose;
    if (options.hybridPhase1)
        edgeLeafTranspose = DCLP2::buildEdgeLeafTranspose(edgeGraph, treeGraphV, options.hybridMinVertexLeafDegree);

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
    daf::StaticVector<DCLP2::LeafRmInfo> leafRmInfo(tree.adj_list.size());
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

    auto ensureLeafTracking = [&](daf::Size leafId) {
        if (leafId >= leafRmInfo.size()) {
            const auto newCap = static_cast<size_t>(leafId * 3 / 2 + 8);
            removedLeaf.reserve(newCap);
            leafRmInfo.resize(newCap);
            leafAffected.resize(newCap, 0);
            leafInDynamicPhase1.resize(newCap, 0);
        } else if (leafId >= leafInDynamicPhase1.size()) {
            leafInDynamicPhase1.resize(leafId + 1, 0);
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
    long long cntB_closedForm = 0, cntB_bk = 0, cntB_complement = 0;
    long long complement_sumNA = 0, complement_maxNA = 0;
    long long complement_sumLeafSize = 0, complement_sumSubLeaves = 0;
    long long us_cbk_total = 0, us_cbk_bk = 0;

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
                    if (entry.type == DCLP2::KK) {
                        markAffectedLeaf(entry.leafId, true, false, edgeU, edgeV, false);
                    } else if (entry.type == DCLP2::PP) {
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
                DCLP2::LeafRmInfo &leafRm = leafRmInfo[leafId];

                const auto& leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;

                tPovit.clear(); tKeepC.clear();
                daf::Size numKeeps = 0;
                for (const auto& node : leaf) {
                    if (node.isPivot) tPovit.push_back(node.v);
                    else { tKeepC.push_back(node.v); numKeeps++; }
                }
                daf::Size needPivot = k - numKeeps;
                daf::Size numPivots = tPovit.size();

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
                                case DCLP2::KK: w = wKK; break;
                                case DCLP2::PP: w = wPP; break;
                                case DCLP2::KP: w = wKP; break;
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

                    if (!leafRm.removedPivots.empty() && needPivot <= tPovit.size() - leafRm.removedPivots.size()) {
                        // Case C alive: leaf survives with fewer pivots
                        if (leafId < leafEdgeAlive.size() && leafEdgeAlive[leafId]) {
                            // ---- DCLP: CSR scan for Case C ----
                            cntC_csr++;
                            daf::Size begin = leafEdgeOff[leafId];
                            daf::Size end = leafEdgeOff[leafId + 1];
                            work_caseC_csr += (end - begin);

                            for (daf::Size pos = begin; pos < end; ++pos) {
                                auto &entry = leafEdgeData[pos];
                                auto [eu, ev] = edgeGraph.getEdgeById(entry.edgeId);
                                bool eu_rm = std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), eu);
                                bool ev_rm = std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), ev);
                                double delta;

                                if (entry.type == DCLP2::KK) {
                                    delta = KtoK;
                                } else if (entry.type == DCLP2::PP) {
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
                                    case DCLP2::KK: w = wKK; break;
                                    case DCLP2::PP: w = wPP; break;
                                    case DCLP2::KP: w = wKP; break;
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
                    }
                }
                leafRmInfo[leafId].clear();
            }
            tPovit.free(); tKeepC.free();
        }

        // ========== Phase 2c: Case B — closed-form d=1, BK fallback d≥2 ==========
        auto _t3 = std::chrono::high_resolution_clock::now();
        us_caseAC_delta += std::chrono::duration_cast<std::chrono::microseconds>(_t3 - _t2).count();
        cntB += caseBLeafIds.size();

        for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
            auto leafId = caseBLeafIds[bi];
            DCLP2::LeafRmInfo &leafRm = leafRmInfo[leafId];
            const auto& leaf = tree.adj_list[leafId];

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
            if (leafId < leafEdgeAlive.size())
                leafEdgeAlive[leafId] = 0;

            auto &leafRef = tree.adj_list[leafId];

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

            // Helper: create a sub-leaf by excluding one vertex
            auto createSubLeaf = [&](daf::Size excludeV) {
                std::vector<TreeGraphNode> sub;
                sub.reserve(leafRef.size() - 1);
                for (const auto &n : leafRef) {
                    if (n.v != excludeV) sub.push_back(n);
                }
                auto newId = tree.addNode(sub);
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
                    DCLP2::processEdgePairs(newKeepC, w, addW);
                }
                int nPP = int(np) - 2;
                if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                    double w = nCr[newPivot.size() - 2][nPP];
                    DCLP2::processEdgePairs(newPivot, w, addW);
                }
                int nKP = int(np) - 1;
                if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                    double w = nCr[newPivot.size() - 1][nKP];
                    DCLP2::processEdgePairs(newKeepC, newPivot, w, addW);
                }
                newPivot.clear(); newKeepC.clear();
                ensureLeafTracking(newId);
            };

            if (activeRmEdges.size() == 0) {
                // All PP edges involved removed pivots — no split needed
                // The leaf (minus removed pivots) is still one clique
                // Re-add it as a "sub-leaf" = the current modified leaf
                std::vector<TreeGraphNode> sub(leafRef.begin(), leafRef.end());
                auto newId = tree.addNode(sub);
                newPivot.clear(); newKeepC.clear();
                for (const auto& i : tree.adj_list[newId]) {
                    if (i.isPivot) newPivot.push_back(i.v);
                    else newKeepC.push_back(i.v);
                }
                for (const auto& i : tree.adj_list[newId]) treeGraphV.addNbr(i.v, {newId, i.isPivot});
                if (options.hybridPhase1) addDynamicLeafMembership(newId);
                daf::Size np = k - newKeepC.size();
                if (np <= newPivot.size() && newKeepC.size() > 1)
                    DCLP2::processEdgePairs(newKeepC, nCr[newPivot.size()][np], addW);
                int nPP = int(np) - 2;
                if (0 <= nPP && nPP <= int(newPivot.size()) - 2)
                    DCLP2::processEdgePairs(newPivot, nCr[newPivot.size() - 2][nPP], addW);
                int nKP = int(np) - 1;
                if (0 <= nKP && nKP <= int(newPivot.size()) - 1)
                    DCLP2::processEdgePairs(newKeepC, newPivot, nCr[newPivot.size() - 1][nKP], addW);
                newPivot.clear(); newKeepC.clear();
                ensureLeafTracking(newId);
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
                createSubLeaf(rp2);

                // Sub-leaf 2: remove rp1, change rp2 from pivot to keep
                {
                    std::vector<TreeGraphNode> sub2;
                    sub2.reserve(leafRef.size() - 1);
                    for (const auto &n : leafRef) {
                        if (n.v == rp1) continue;
                        if (n.v == rp2) sub2.push_back({rp2, false}); // pivot → keep
                        else sub2.push_back(n);
                    }
                    auto newId = tree.addNode(sub2);
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
                        DCLP2::processEdgePairs(newKeepC, w, addW);
                    }
                    int nPP = int(np) - 2;
                    if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                        double w = nCr[newPivot.size() - 2][nPP];
                        DCLP2::processEdgePairs(newPivot, w, addW);
                    }
                    int nKP = int(np) - 1;
                    if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                        double w = nCr[newPivot.size() - 1][nKP];
                        DCLP2::processEdgePairs(newKeepC, newPivot, w, addW);
                    }
                    newPivot.clear(); newKeepC.clear();
                    ensureLeafTracking(newId);
                }
                cntB_closedForm++;
            } else {
                // ---- Complement BK: d≥2 ----
                auto _cbk_start = std::chrono::high_resolution_clock::now();
                int leafSize = (int)leafRef.size();

                // Collect affected vertex IDs (sorted, deduplicated)
                std::vector<daf::Size> affVerts;
                affVerts.reserve(2 * activeRmEdges.size());
                for (auto &[eu, ev] : activeRmEdges) {
                    affVerts.push_back(eu);
                    affVerts.push_back(ev);
                }
                std::sort(affVerts.begin(), affVerts.end());
                affVerts.erase(std::unique(affVerts.begin(), affVerts.end()), affVerts.end());
                int nA = (int)affVerts.size();

                complement_sumNA += nA;
                complement_sumLeafSize += leafSize;
                if (nA > complement_maxNA) complement_maxNA = nA;

                // Map: vertex ID -> small index via sorted position
                auto affIdx = [&](daf::Size v) -> int {
                    return (int)(std::lower_bound(affVerts.begin(), affVerts.end(), v) - affVerts.begin());
                };

                // Collect unaffected vertices (those NOT in affVerts)
                std::vector<TreeGraphNode> unaffected;
                unaffected.reserve(leafRef.size());
                for (const auto &nd : leafRef) {
                    if (!std::binary_search(affVerts.begin(), affVerts.end(), nd.v))
                        unaffected.push_back(nd);
                }

                // Build small adjacency (K_{nA} minus removed edges)
                using Bitset = bkRmEdge::Bitset;
                std::vector<Bitset> smallAdj(nA, Bitset(nA));
                for (int i = 0; i < nA; i++) { smallAdj[i].set(); smallAdj[i].reset(i); }
                for (auto &[eu, ev] : activeRmEdges) {
                    int ui = affIdx(eu), vi = affIdx(ev);
                    smallAdj[ui].reset(vi);
                    smallAdj[vi].reset(ui);
                }

                // Run BK on small graph
                Bitset Rinit(nA), Pinit(nA), pivInit(nA);
                Rinit.reset(); Pinit.set(); pivInit.reset();
                int minKSmall = std::max(1, (int)k - (int)unaffected.size());
                auto _cbk_bk_start = std::chrono::high_resolution_clock::now();
                int numSubLeaves_this = 0;

                bkRmEdge::bk_run(smallAdj, nA, minKSmall, Rinit, Pinit, pivInit,
                    [&](const Bitset &clique, const Bitset &pivBits) {
                        numSubLeaves_this++;
                        // Build sub-leaf: unaffected nodes + BK clique among affected
                        std::vector<TreeGraphNode> subLeaf;
                        subLeaf.reserve(unaffected.size() + clique.count());
                        for (const auto &nd : unaffected) subLeaf.push_back(nd);
                        for (size_t i = clique.find_first(); i != Bitset::npos && (int)i < nA;
                             i = clique.find_next(i))
                            subLeaf.push_back({affVerts[i], pivBits.test(i)});
                        if ((int)subLeaf.size() < (int)k) return;

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
                            DCLP2::processEdgePairs(newKeepC, w, addW);
                        }
                        int nPP = int(np) - 2;
                        if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                            double w = nCr[newPivot.size() - 2][nPP];
                            DCLP2::processEdgePairs(newPivot, w, addW);
                        }
                        int nKP = int(np) - 1;
                        if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                            double w = nCr[newPivot.size() - 1][nKP];
                            DCLP2::processEdgePairs(newKeepC, newPivot, w, addW);
                        }
                        newPivot.clear(); newKeepC.clear();
                        ensureLeafTracking(newId);
                    });
                auto _cbk_end = std::chrono::high_resolution_clock::now();
                complement_sumSubLeaves += numSubLeaves_this;
                us_cbk_total += std::chrono::duration_cast<std::chrono::microseconds>(_cbk_end - _cbk_start).count();
                us_cbk_bk += std::chrono::duration_cast<std::chrono::microseconds>(_cbk_bk_start - _cbk_start).count(); // setup
                cntB_complement++;
            }
            activeRmEdges.free();

            // Remove old contribution
            auto removeW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] -= w;
                if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx])
                    pq.update(pqHandles[idx]);
            };
            if (needPivot <= povit.size()) {
                double w = nCr[povit.size()][needPivot];
                DCLP2::processEdgePairs(keepC, w, removeW);
            }
            int needPP = int(needPivot) - 2;
            if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                double w = nCr[povit.size() - 2][needPP];
                DCLP2::processEdgePairs(povit, w, removeW);
            }
            int needKP = int(needPivot) - 1;
            if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                double w = nCr[povit.size() - 1][needKP];
                DCLP2::processEdgePairs(keepC, povit, w, removeW);
            }

            tree.adj_list[leafId].clear();
            leafRmInfo[leafId].clear();
            povit.clear(); keepC.clear();
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
    std::cout << "  Phase1(edgePathCSR): " << us_phase1/1000 << " ms" << std::endl;
    std::cout << "  Phase2(A+C delta+mut): " << us_caseAC_delta/1000 << " ms" << std::endl;
    std::cout << "  CaseB: closedForm=" << cntB_closedForm << " complement=" << cntB_complement << " bk=" << cntB_bk << std::endl;
    if (cntB_complement > 0) {
        std::cout << "  ComplementBK: avgNA=" << (complement_sumNA / cntB_complement)
                  << " maxNA=" << complement_maxNA
                  << " avgLeafSize=" << (complement_sumLeafSize / cntB_complement)
                  << " avgSubLeaves=" << (complement_sumSubLeaves / cntB_complement)
                  << " totalSubLeaves=" << complement_sumSubLeaves << std::endl;
        std::cout << "  ComplementBK timing: total=" << us_cbk_total/1000 << "ms setup=" << us_cbk_bk/1000 << "ms" << std::endl;
    }
    std::cout << "  CaseB(total): " << us_caseB/1000 << " ms" << std::endl;
    std::cout << "  Flush(bucketMove): " << us_flush/1000 << " ms" << std::endl;
    std::cout << "  P1: alive=" << work_p1_alive << std::endl;
    std::cout << "  CaseA_CSR=" << work_caseA_csr << " CaseA_fallback=" << work_caseA_fallback << std::endl;
    std::cout << "  CaseC_CSR=" << work_caseC_csr << " CaseC_fallback=" << work_caseC_fallback << std::endl;

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

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DCLP2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    static const DCLP2::DCLPOptions options{
        "DCLP2",
        false,
        0
    };
    return PlusNucleusEdgeCoreDecompositionSet_DCLP_Impl(tree, edgeGraph, treeGraphV, k, options);
}
