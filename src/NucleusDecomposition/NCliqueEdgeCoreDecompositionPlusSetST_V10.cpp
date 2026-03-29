//
// R2 ST V10: CSR fast-path for Case C + overflow CSR for Case B sub-leaves
//
// Two provably-beneficial optimizations over V8b:
//
// 1. Case C CSR fast-path:
//    When leafEdgeAlive[leafId], scan CSR entries with getEdgeById (O(1) array)
//    instead of O(|leaf|²) getEdgeCompressedId (binary search).
//    For each entry, check if endpoints are removed pivots to determine delta.
//    ~2.3x per-entry speedup (11 ns vs 25 ns).
//
// 2. Case B overflow CSR:
//    In BK callback's addW, capture edgeId and store in dynamic CSR buffer.
//    Later Case A of overflow leaf uses CSR fast-path (5 ns) instead of
//    fallback getEdgeCompressedId (20 ns). Zero additional getEdgeId calls.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>

#include "../BK/BronKerboschRmEdge.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace PlusECDSet_ST_V10 {

    enum EdgeType : uint8_t { KK = 0, PP = 1, KP = 2 };

    struct LeafEdgeEntry {
        daf::Size edgeId;
        EdgeType type;
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

    // CSR-based init result
    struct InitResult {
        double *countingKE;
        std::vector<daf::Size> leafEdgeOff;
        std::vector<LeafEdgeEntry> leafEdgeData;
        std::vector<double> leafWKK, leafWPP, leafWKP;
    };

    InitResult countingPerEdgeWithIndex(const DynamicGraph<TreeGraphNode> &treeGraph,
                                        const Graph &edgeGraph,
                                        const daf::CliqueSize k) {
        const daf::Size numEdges = edgeGraph.adj_list.size();
        const daf::Size numLeaves = treeGraph.adj_list.size();

        auto *countingE = new double[numEdges];
        memset(countingE, 0, numEdges * sizeof(double));

        std::vector<double> leafWKK(numLeaves, 0), leafWPP(numLeaves, 0), leafWKP(numLeaves, 0);

        // --- Pass 1: count edges per leaf ---
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

            int needPivot = int(k) - int(tKeepC.size());
            daf::Size cnt = 0;

            if (needPivot >= 0 && needPivot <= int(tPovit.size())) {
                cnt += tKeepC.size() * (tKeepC.size() - 1) / 2;
            }
            int needPP = needPivot - 2;
            if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
                cnt += tPovit.size() * (tPovit.size() - 1) / 2;
            }
            int needKP = needPivot - 1;
            if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
                cnt += tKeepC.size() * tPovit.size();
            }
            leafEdgeOff[li + 1] = cnt;
        }

        for (daf::Size li = 0; li < numLeaves; ++li) {
            leafEdgeOff[li + 1] += leafEdgeOff[li];
        }
        size_t totalEntries = leafEdgeOff[numLeaves];

        // --- Pass 2: fill data + compute countingE and weights ---
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
                double w = std::llround(nCr[tPovit.size()][needPivot]);
                leafWKK[li] = w;
                for (size_t i = 0; i < tKeepC.size(); ++i)
                    for (size_t j = i + 1; j < tKeepC.size(); ++j) {
                        auto eid = edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]);
                        countingE[eid] += w;
                        leafEdgeData[pos++] = {eid, KK};
                    }
            }

            int needPP = needPivot - 2;
            if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
                double w = std::llround(nCr[tPovit.size() - 2][needPP]);
                leafWPP[li] = w;
                for (size_t i = 0; i < tPovit.size(); ++i)
                    for (size_t j = i + 1; j < tPovit.size(); ++j) {
                        auto eid = edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]);
                        countingE[eid] += w;
                        leafEdgeData[pos++] = {eid, PP};
                    }
            }

            int needKP = needPivot - 1;
            if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
                double w = std::llround(nCr[tPovit.size() - 1][needKP]);
                leafWKP[li] = w;
                for (size_t i = 0; i < tKeepC.size(); ++i)
                    for (size_t j = 0; j < tPovit.size(); ++j) {
                        auto eid = edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]);
                        countingE[eid] += w;
                        leafEdgeData[pos++] = {eid, KP};
                    }
            }
            assert(pos == leafEdgeOff[li + 1]);
        }
        tPovit.free(); tKeepC.free();
        return {countingE, std::move(leafEdgeOff), std::move(leafEdgeData),
                std::move(leafWKK), std::move(leafWPP), std::move(leafWKP)};
    }
}

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V10(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    auto time_start = std::chrono::high_resolution_clock::now();

    // CSR init (same as V8b)
    auto initResult = PlusECDSet_ST_V10::countingPerEdgeWithIndex(tree, edgeGraph, k);
    auto *countingKE = initResult.countingKE;
    auto &leafEdgeOff = initResult.leafEdgeOff;
    auto &leafEdgeData = initResult.leafEdgeData;
    auto &leafWKK = initResult.leafWKK;
    auto &leafWPP = initResult.leafWPP;
    auto &leafWKP = initResult.leafWKP;

    auto time_init = std::chrono::high_resolution_clock::now();
    std::cout << "V10 init took: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_init - time_start).count()
              << " ms, leafEdgeData entries: " << leafEdgeData.size() << std::endl;

    // leafEdgeAlive[li]: whether leaf li's CSR data is still valid
    std::vector<uint8_t> leafEdgeAlive(leafEdgeOff.size() > 0 ? leafEdgeOff.size() - 1 : 0, 1);

    // --- V10 Optimization 2: Dynamic CSR for overflow leaves from Case B ---
    // Use vectors indexed by leafId. Original leaves use static CSR (leafEdge*).
    // Overflow leaves (from Case B) get dynamic per-leaf CSR stored here.
    daf::Size origNumLeaves = leafEdgeOff.size() > 0 ? leafEdgeOff.size() - 1 : 0;
    // Dynamic CSR for overflow leaves, indexed by leafId
    // Initially sized to tree capacity; resized when new IDs are allocated.
    daf::Size ovfCapacity = tree.adj_list.size();
    std::vector<std::vector<PlusECDSet_ST_V10::LeafEdgeEntry>> ovfCSR(ovfCapacity);
    std::vector<long long> ovfWKK(ovfCapacity, 0), ovfWPP(ovfCapacity, 0), ovfWKP(ovfCapacity, 0);
    std::vector<uint8_t> ovfLeafAlive(ovfCapacity, 0);

    const daf::Size numEdgesInit = edgeGraph.adj_list.size();

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
    daf::StaticVector<PlusECDSet_ST_V10::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();

    double currCore = 0;

    // --- Bucket array ---
    const daf::Size numEdges = edgeGraph.adj_list.size();
    int maxBucket = 0;
    for (daf::Size i = 0; i < numEdges; ++i) {
        if (countingKE[i] > 0)
            maxBucket = std::max(maxBucket, (int)countingKE[i]);
    }
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(numEdges);
    std::vector<daf::Size> pos_in_bucket(numEdges);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numEdges; ++i) {
        if (countingKE[i] == 0) {
            edgeInHeap[i] = 0;
            continue;
        }
        int b = (int)countingKE[i];
        bucket_of[i] = b;
        pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i);
        remainingInHeap++;
    }
    int curBucket = 0;

    auto bucketMove = [&](daf::Size id) {
        int newB = std::max(0, (int)countingKE[id]);
        int oldB = bucket_of[id];
        if (newB == oldB) return;
        auto& oldVec = buckets[oldB];
        daf::Size myPos = pos_in_bucket[id];
        if (myPos < oldVec.size() - 1) {
            daf::Size last = oldVec.back();
            oldVec[myPos] = last;
            pos_in_bucket[last] = myPos;
        }
        oldVec.pop_back();
        if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
        bucket_of[id] = newB;
        pos_in_bucket[id] = buckets[newB].size();
        buckets[newB].push_back(id);
        if (newB < curBucket) curBucket = newB;
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

    // Helper: ensure overflow vectors are large enough for leafId
    auto ensureOvfCapacity = [&](daf::Size leafId) {
        if (leafId >= ovfCapacity) {
            daf::Size newCap = leafId * 1.5 + 1;
            ovfCSR.resize(newCap);
            ovfWKK.resize(newCap, 0);
            ovfWPP.resize(newCap, 0);
            ovfWKP.resize(newCap, 0);
            ovfLeafAlive.resize(newCap, 0);
            ovfCapacity = newCap;
        }
    };

    // Helper: check if leaf has valid CSR (works for both original and overflow leaves)
    auto hasCSR = [&](daf::Size leafId) -> bool {
        if (leafId < origNumLeaves) {
            return leafEdgeAlive[leafId] != 0;
        } else {
            return leafId < ovfCapacity && ovfLeafAlive[leafId] != 0;
        }
    };

    // Helper: invalidate CSR for a leaf
    auto invalidateCSR = [&](daf::Size leafId) {
        if (leafId < origNumLeaves) {
            leafEdgeAlive[leafId] = 0;
        } else if (leafId < ovfCapacity) {
            ovfLeafAlive[leafId] = 0;
        }
    };

    std::cout << "=========================begin=========================" << std::endl;
    double minCore = 0;
    long long cntA = 0, cntB = 0, cntC = 0, totalIters = 0;
    long long cntA_fast = 0, cntA_fallback = 0;
    long long cntC_csr = 0, cntC_fallback = 0;
    long long us_phase1 = 0, us_caseA_delta = 0;
    long long us_caseB = 0, us_flush = 0;

    while (remainingInHeap > 0) {
        // --- Bucket pop ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((double)curBucket, minCore);

        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                edgeInHeap[id] = 0;
                currentRemoveEdgeIds.push_back(id);
                coreE[id] = minCore;
                remainingInHeap--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (curBucket + 1) <= (int)minCore) {
                curBucket++;
            } else break;
        }

        currCore = minCore;

        if (remainingInHeap == 0) break;

        totalIters++;

        // --- Phase 1: intersect edge adjacency lists ---
        auto _t0 = std::chrono::high_resolution_clock::now();
        for (int ei = 0; ei < (int)currentRemoveEdgeIds.size(); ++ei) {
            auto edgeId = currentRemoveEdgeIds[ei];
            auto [edgeU, edgeV] = edgeGraph.getEdgeById(edgeId);
            auto &adjEdgeU = treeGraphV.getNbr(edgeU);
            auto &adjEdgeV = treeGraphV.getNbr(edgeV);
            daf::intersect_dense_sets(adjEdgeU, adjEdgeV,
                                      [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                                          daf::Size leafId = uClique.v;
                                          auto &lrm = leafRmInfo[leafId];
                                          if (lrm.empty()) {
                                              removedLeaf.push_back(leafId);
                                              leafAffected[leafId] = 1;
                                          }
                                          if (!lrm.removedKeepC) {
                                              if (!uClique.isPivot && !vClique.isPivot) {
                                                  lrm.removedKeepC = true;
                                              } else if (uClique.isPivot && vClique.isPivot) {
                                                  daf::Size eu = edgeU, ev = edgeV;
                                                  lrm.removedEdges.push_back({eu, ev});
                                              } else if (uClique.isPivot) {
                                                  lrm.removedPivots.push_back(edgeU);
                                              } else {
                                                  lrm.removedPivots.push_back(edgeV);
                                              }
                                          }
                                      });
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
                    markDirty(idx);
            };

            for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                auto leafId = removedLeaf[leafIdIdx];
                PlusECDSet_ST_V10::LeafRmInfo &leafRm = leafRmInfo[leafId];

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
                    // ---- Case A ----
                    cntA++;

                    // Try CSR fast path (original leaf)
                    if (leafId < origNumLeaves && leafEdgeAlive[leafId]) {
                        cntA_fast++;
                        double wKK = leafWKK[leafId];
                        double wPP = leafWPP[leafId];
                        double wKP = leafWKP[leafId];
                        daf::Size begin = leafEdgeOff[leafId];
                        daf::Size end = leafEdgeOff[leafId + 1];
                        for (daf::Size pos = begin; pos < end; ++pos) {
                            auto &entry = leafEdgeData[pos];
                            double w;
                            switch (entry.type) {
                                case PlusECDSet_ST_V10::KK: w = wKK; break;
                                case PlusECDSet_ST_V10::PP: w = wPP; break;
                                case PlusECDSet_ST_V10::KP: w = wKP; break;
                            }
                            directSub(entry.edgeId, w);
                        }
                        leafEdgeAlive[leafId] = 0;
                    }
                    // V10: Try overflow CSR fast path
                    else if (leafId >= origNumLeaves && leafId < ovfCapacity && ovfLeafAlive[leafId]) {
                        cntA_fast++;
                        double wKK = ovfWKK[leafId];
                        double wPP = ovfWPP[leafId];
                        double wKP = ovfWKP[leafId];
                        auto &entries = ovfCSR[leafId];
                        for (auto &entry : entries) {
                            double w;
                            switch (entry.type) {
                                case PlusECDSet_ST_V10::KK: w = wKK; break;
                                case PlusECDSet_ST_V10::PP: w = wPP; break;
                                case PlusECDSet_ST_V10::KP: w = wKP; break;
                            }
                            directSub(entry.edgeId, w);
                        }
                        ovfLeafAlive[leafId] = 0;
                        entries.clear();
                        entries.shrink_to_fit();
                    }
                    else {
                        // Original leaf but CSR invalidated (by prior Case C)
                        cntA_fallback++;
                        double KtoK = 0, KtoP = 0, PtoP = 0;
                        if (needPivot <= tPovit.size()) {
                            KtoK = std::llround(nCr[tPovit.size()][needPivot]);
                            for (daf::Size i = 0; i + 1 < tKeepC.size(); ++i)
                                for (daf::Size j = i + 1; j < tKeepC.size(); ++j)
                                    directSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]), KtoK);
                        }
                        int needPP = int(needPivot) - 2;
                        if (0 <= needPP && needPP <= int(tPovit.size()) - 2) {
                            PtoP = std::llround(nCr[tPovit.size() - 2][needPP]);
                            for (daf::Size i = 0; i + 1 < tPovit.size(); ++i)
                                for (daf::Size j = i + 1; j < tPovit.size(); ++j)
                                    directSub(edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]), PtoP);
                        }
                        int needKP = int(needPivot) - 1;
                        if (0 <= needKP && needKP <= int(tPovit.size()) - 1) {
                            KtoP = std::llround(nCr[tPovit.size() - 1][needKP]);
                            for (daf::Size i = 0; i < tKeepC.size(); ++i)
                                for (daf::Size j = 0; j < tPovit.size(); ++j)
                                    directSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]), KtoP);
                        }
                    }
                    // Tree mutation
                    for (const auto& i : leaf)
                        treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                    tree.adj_list[leafId].clear();
                    tree.recycleNode(leafId);
                } else {
                    // ---- Case C ----
                    cntC++;
                    double KtoK = 0, KtoP = 0, PtoP = 0;
                    double RemovedKtoK = 0, RemovedKtoP = 0, RemovedPtoP = 0;

                    if (needPivot <= tPovit.size()) {
                        KtoK = std::llround(nCr[tPovit.size()][needPivot]) - std::llround(nCr[tPovit.size() - leafRm.removedPivots.size()][needPivot]);
                        RemovedKtoK = std::llround(nCr[tPovit.size()][needPivot]);
                    }
                    int needPP = int(needPivot) - 2;
                    if (0 <= needPP && needPP <= int(tPovit.size()) - 2) {
                        RemovedPtoP = std::llround(nCr[tPovit.size() - 2][needPP]);
                        PtoP = RemovedPtoP;
                        if (leafRm.removedPivots.size() + 2 <= tPovit.size())
                            PtoP -= std::llround(nCr[tPovit.size() - 2 - leafRm.removedPivots.size()][needPP]);
                    }
                    int needKP = int(needPivot) - 1;
                    if (0 <= needKP && needKP <= int(tPovit.size()) - 1) {
                        RemovedKtoP = std::llround(nCr[tPovit.size() - 1][needKP]);
                        KtoP = RemovedKtoP;
                        if (leafRm.removedPivots.size() + 1 <= tPovit.size())
                            KtoP -= std::llround(nCr[tPovit.size() - 1 - leafRm.removedPivots.size()][needKP]);
                    }

                    if (!leafRm.removedPivots.empty() && needPivot <= tPovit.size() - leafRm.removedPivots.size()) {
                        // Case C-alive: leaf survives with fewer pivots

                        // === V10 OPTIMIZATION 1: CSR fast-path for Case C ===
                        // Instead of O(|leaf|²) getEdgeCompressedId, scan CSR entries.
                        // For each entry, use getEdgeById (O(1)) to get endpoints,
                        // check if either endpoint is a removed pivot, and apply delta.
                        if (hasCSR(leafId)) {
                            cntC_csr++;

                            // Build a small lookup for removed pivots
                            // (typically 1-5 entries, linear scan is fastest)
                            auto &rmPivots = leafRm.removedPivots;

                            auto isRemovedPivot = [&](daf::Size v) -> bool {
                                for (daf::Size i = 0; i < rmPivots.size(); ++i)
                                    if (rmPivots[i] == v) return true;
                                return false;
                            };

                            // Scan CSR entries
                            if (leafId < origNumLeaves) {
                                daf::Size begin = leafEdgeOff[leafId];
                                daf::Size end = leafEdgeOff[leafId + 1];
                                for (daf::Size pos = begin; pos < end; ++pos) {
                                    auto &entry = leafEdgeData[pos];
                                    auto [eu, ev] = edgeGraph.getEdgeById(entry.edgeId);
                                    bool euRm = isRemovedPivot(eu);
                                    bool evRm = isRemovedPivot(ev);

                                    double delta;
                                    if (euRm || evRm) {
                                        // Edge involves at least one removed pivot: subtract full old weight
                                        switch (entry.type) {
                                            case PlusECDSet_ST_V10::KK: delta = RemovedKtoK; break;
                                            case PlusECDSet_ST_V10::PP: delta = RemovedPtoP; break;
                                            case PlusECDSet_ST_V10::KP: delta = RemovedKtoP; break;
                                        }
                                    } else {
                                        // Edge among remaining vertices: subtract (old - new)
                                        switch (entry.type) {
                                            case PlusECDSet_ST_V10::KK: delta = KtoK; break;
                                            case PlusECDSet_ST_V10::PP: delta = PtoP; break;
                                            case PlusECDSet_ST_V10::KP: delta = KtoP; break;
                                        }
                                    }
                                    if (delta > 0)
                                        directSub(entry.edgeId, delta);
                                }
                            } else {
                                // Overflow leaf CSR
                                auto &entries = ovfCSR[leafId];
                                for (auto &entry : entries) {
                                    auto [eu, ev] = edgeGraph.getEdgeById(entry.edgeId);
                                    bool euRm = isRemovedPivot(eu);
                                    bool evRm = isRemovedPivot(ev);

                                    double delta;
                                    if (euRm || evRm) {
                                        switch (entry.type) {
                                            case PlusECDSet_ST_V10::KK: delta = RemovedKtoK; break;
                                            case PlusECDSet_ST_V10::PP: delta = RemovedPtoP; break;
                                            case PlusECDSet_ST_V10::KP: delta = RemovedKtoP; break;
                                        }
                                    } else {
                                        switch (entry.type) {
                                            case PlusECDSet_ST_V10::KK: delta = KtoK; break;
                                            case PlusECDSet_ST_V10::PP: delta = PtoP; break;
                                            case PlusECDSet_ST_V10::KP: delta = KtoP; break;
                                        }
                                    }
                                    if (delta > 0)
                                        directSub(entry.edgeId, delta);
                                }
                                // Free overflow CSR
                                entries.clear();
                                entries.shrink_to_fit();
                            }
                            // Invalidate CSR after Case C (pivot composition changed)
                            invalidateCSR(leafId);
                        } else {
                            // No CSR — fallback to original O(|leaf|²) path
                            cntC_fallback++;
                            daf::StaticVector<TreeGraphNode> newLeafF;
                            for (const auto& node : leaf) {
                                if (!std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), node.v))
                                    newLeafF.push_back(node);
                            }
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
                        // Tree mutation for Case C-alive
                        for (auto removedNbr : leafRm.removedPivots)
                            treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                        tree.removeNbrs(leafId, leafRm.removedPivots);
                        invalidateCSR(leafId);
                    } else {
                        // Case C-dead: leaf dies (all pivots removed or removedPivots empty)
                        // Same as Case A but counted as C
                        // Try CSR
                        if (hasCSR(leafId)) {
                            cntC_csr++;
                            if (leafId < origNumLeaves) {
                                daf::Size begin = leafEdgeOff[leafId];
                                daf::Size end = leafEdgeOff[leafId + 1];
                                for (daf::Size pos = begin; pos < end; ++pos) {
                                    auto &entry = leafEdgeData[pos];
                                    double w;
                                    switch (entry.type) {
                                        case PlusECDSet_ST_V10::KK: w = RemovedKtoK; break;
                                        case PlusECDSet_ST_V10::PP: w = RemovedPtoP; break;
                                        case PlusECDSet_ST_V10::KP: w = RemovedKtoP; break;
                                    }
                                    directSub(entry.edgeId, w);
                                }
                            } else {
                                auto &entries = ovfCSR[leafId];
                                for (auto &entry : entries) {
                                    double w;
                                    switch (entry.type) {
                                        case PlusECDSet_ST_V10::KK: w = RemovedKtoK; break;
                                        case PlusECDSet_ST_V10::PP: w = RemovedPtoP; break;
                                        case PlusECDSet_ST_V10::KP: w = RemovedKtoP; break;
                                    }
                                    directSub(entry.edgeId, w);
                                }
                                entries.clear();
                                entries.shrink_to_fit();
                            }
                            invalidateCSR(leafId);
                        } else {
                            cntC_fallback++;
                            for (daf::Size i = 0; i + 1 < leaf.size(); ++i)
                                for (daf::Size j = i + 1; j < leaf.size(); ++j) {
                                    auto &u = leaf[i], &v = leaf[j];
                                    double d = (!u.isPivot && !v.isPivot) ? RemovedKtoK : (u.isPivot && v.isPivot) ? RemovedPtoP : RemovedKtoP;
                                    directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                                }
                        }
                        for (const auto& i : leaf)
                            treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                        tree.adj_list[leafId].clear();
                        tree.recycleNode(leafId);
                        invalidateCSR(leafId);
                    }
                }
                leafRmInfo[leafId].clear();
            }
            tPovit.free(); tKeepC.free();
        }

        // ========== Case B: BK + apply ==========
        auto _t3 = std::chrono::high_resolution_clock::now();
        us_caseA_delta += std::chrono::duration_cast<std::chrono::microseconds>(_t3 - _t2).count();
        cntB += caseBLeafIds.size();

        for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
            auto leafId = caseBLeafIds[bi];
            PlusECDSet_ST_V10::LeafRmInfo &leafRm = leafRmInfo[leafId];
            const auto& leaf = tree.adj_list[leafId];

            povit.clear(); keepC.clear();
            for (const auto& node : leaf) {
                if (node.isPivot) povit.push_back(node.v);
                else keepC.push_back(node.v);
            }
            daf::Size needPivot = k - keepC.size();
            double KtoK = 0, KtoP = 0, PtoP = 0;

            // Remove old leaf from treeGraphV
            for (const auto& leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }
            if (!leafRm.removedPivots.empty()) {
                tree.removeNbrs(leafId, leafRm.removedPivots);
            }

            // Invalidate CSR for old leaf
            invalidateCSR(leafId);

            // === V10 OPTIMIZATION 2: Build CSR for overflow sub-leaves in BK callback ===
            auto &leafRef = tree.adj_list[leafId];
            bkRmEdge::bronKerbosch(leafRef, leafRm.removedEdges, k,
                [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                    auto subLeaf = bkRmEdge::coverToVertex(c, pivots, leafRef);
                    auto newId = tree.addNode(subLeaf);
                    newPivot.clear(); newKeepC.clear();
                    for (const auto& i : tree.adj_list[newId]) {
                        if (i.isPivot) { newPivot.push_back(i.v); treeGraphV.addNbr(i.v, {newId, true}); }
                        else { newKeepC.push_back(i.v); treeGraphV.addNbr(i.v, {newId, false}); }
                    }
                    daf::Size np = k - newKeepC.size();

                    // Build CSR for this overflow leaf while computing addW
                    std::vector<PlusECDSet_ST_V10::LeafEdgeEntry> newCSREntries;
                    {
                        daf::Size nk = newKeepC.size(), npv = newPivot.size();
                        newCSREntries.reserve(nk*(nk-1)/2 + npv*(npv-1)/2 + nk*npv);
                    }
                    long long newWKK = 0, newWPP = 0, newWKP = 0;

                    // addW that also captures edgeId into CSR
                    auto addW_csr = [&](daf::Size u, daf::Size v, long long w,
                                        PlusECDSet_ST_V10::EdgeType type) {
                        auto idx = edgeGraph.getEdgeCompressedId(u, v);
                        countingKE[idx] += w;
                        if (edgeInHeap[idx])
                            markDirty(idx);
                        newCSREntries.push_back({idx, type});
                    };

                    double wKK = 0, wPP = 0, wKP = 0;
                    if (np <= newPivot.size() && newKeepC.size() > 1) {
                        wKK = std::llround(nCr[newPivot.size()][np]);
                        newWKK = wKK;
                        for (daf::Size i = 0; i < newKeepC.size(); ++i)
                            for (daf::Size j = i + 1; j < newKeepC.size(); ++j)
                                addW_csr(newKeepC[i], newKeepC[j], wKK, PlusECDSet_ST_V10::KK);
                    }
                    int nPP = int(np) - 2;
                    if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                        wPP = std::llround(nCr[newPivot.size() - 2][nPP]);
                        newWPP = wPP;
                        for (daf::Size i = 0; i < newPivot.size(); ++i)
                            for (daf::Size j = i + 1; j < newPivot.size(); ++j)
                                addW_csr(newPivot[i], newPivot[j], wPP, PlusECDSet_ST_V10::PP);
                    }
                    int nKP = int(np) - 1;
                    if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                        wKP = std::llround(nCr[newPivot.size() - 1][nKP]);
                        newWKP = wKP;
                        for (daf::Size i = 0; i < newKeepC.size(); ++i)
                            for (daf::Size j = 0; j < newPivot.size(); ++j)
                                addW_csr(newKeepC[i], newPivot[j], wKP, PlusECDSet_ST_V10::KP);
                    }

                    // Store overflow CSR
                    if (!newCSREntries.empty()) {
                        ensureOvfCapacity(newId);
                        ovfCSR[newId] = std::move(newCSREntries);
                        ovfWKK[newId] = newWKK;
                        ovfWPP[newId] = newWPP;
                        ovfWKP[newId] = newWKP;
                        ovfLeafAlive[newId] = 1;
                    }

                    newPivot.clear(); newKeepC.clear();
                    if (newId >= leafRmInfo.size()) {
                        removedLeaf.reserve(newId * 1.5);
                        leafRmInfo.resize(newId * 1.5);
                        leafAffected.resize(newId * 1.5, 0);
                    }
                });

            // Remove old contribution
            auto removeW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] -= w;
                if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx])
                    markDirty(idx);
            };
            if (needPivot <= povit.size()) {
                KtoK = std::llround(nCr[povit.size()][needPivot]);
                PlusECDSet_ST_V10::processEdgePairs(keepC, KtoK, removeW);
            }
            int needPP = int(needPivot) - 2;
            if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                PtoP = std::llround(nCr[povit.size() - 2][needPP]);
                PlusECDSet_ST_V10::processEdgePairs(povit, PtoP, removeW);
            }
            int needKP = int(needPivot) - 1;
            if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                KtoP = std::llround(nCr[povit.size() - 1][needKP]);
                PlusECDSet_ST_V10::processEdgePairs(keepC, povit, KtoP, removeW);
            }

            tree.removeNode(leafId);
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
    std::cout << "  CaseA: fast=" << cntA_fast << " fallback=" << cntA_fallback << std::endl;
    std::cout << "  CaseC: csr=" << cntC_csr << " fallback=" << cntC_fallback << std::endl;
    std::cout << "  Phase1(intersect): " << us_phase1/1000 << " ms" << std::endl;
    std::cout << "  Phase2(A+C delta+mut): " << us_caseA_delta/1000 << " ms" << std::endl;
    std::cout << "  CaseB(BK+apply): " << us_caseB/1000 << " ms" << std::endl;
    std::cout << "  Flush(bucketMove): " << us_flush/1000 << " ms" << std::endl;

    daf::Size numCounting = 0;
    std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> sortedK;
    sortedK.reserve(edgeGraph.adj_list.size());

    const daf::Size m = edgeGraph.adj_list.size();
    const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;
    for (daf::Size u = 0; u < n; ++u) {
        const daf::Size start = edgeGraph.adj_list_offsets[u];
        const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
        for (daf::Size idx = start; idx < end; ++idx) {
            sortedK.emplace_back(
                std::make_pair(std::make_pair(u, edgeGraph.adj_list[idx]), (int) coreE[idx]));
        }
    }

    assert(numCounting == 0);

    delete[] countingKE;
    delete[] coreE;
    povit.free();
    keepC.free();
    newPivot.free();
    newKeepC.free();
    return sortedK;
}
