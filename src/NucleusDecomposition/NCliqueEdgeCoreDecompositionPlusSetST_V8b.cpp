//
// R2 ST V8b: CSR init + deltaAccum batch bucketMove
//
// Two optimizations layered on V8's Case A fast-path:
//
// 1. CSR init: Replace vector<vector<LeafEdgeEntry>> with flat CSR
//    (leafEdgeOff[numLeaves+1] + leafEdgeData[totalEntries]).
//    Single allocation instead of millions of small vectors.
//
// 2. deltaAccum: In directSub, accumulate deltas without calling bucketMove.
//    Flush all dirty edges once at end of Phase 2 (after Case A/C/B).
//    Eliminates redundant bucket moves for edges hit multiple times.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>

#include "../BK/BronKerboschRmEdge.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace PlusECDSet_ST_V8b {

    enum EdgeType : uint8_t { KK = 0, PP = 1, KP = 2 };

    struct LeafEdgeEntry {
        daf::Size edgeId;
        EdgeType type;
    };

    struct LeafRmInfo {
        bool removedKeepC;
        daf::StaticVector<daf::Size> removedPivots{0};
        daf::StaticVector<std::pair<daf::Size, daf::Size> > removedEdges{0};

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
        // Flat CSR: leafEdgeData[leafEdgeOff[li] .. leafEdgeOff[li+1])
        std::vector<daf::Size> leafEdgeOff;   // size = numLeaves + 1
        std::vector<LeafEdgeEntry> leafEdgeData; // size = totalEntries
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
                cnt += tKeepC.size() * (tKeepC.size() - 1) / 2; // KK pairs
            }
            int needPP = needPivot - 2;
            if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
                cnt += tPovit.size() * (tPovit.size() - 1) / 2; // PP pairs
            }
            int needKP = needPivot - 1;
            if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
                cnt += tKeepC.size() * tPovit.size(); // KP pairs
            }
            leafEdgeOff[li + 1] = cnt;
        }

        // Prefix sum
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

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_ST_V8b(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    auto time_start = std::chrono::high_resolution_clock::now();

    // V8b: CSR init
    auto initResult = PlusECDSet_ST_V8b::countingPerEdgeWithIndex(tree, edgeGraph, k);
    auto *countingKE = initResult.countingKE;
    auto &leafEdgeOff = initResult.leafEdgeOff;
    auto &leafEdgeData = initResult.leafEdgeData;
    auto &leafWKK = initResult.leafWKK;
    auto &leafWPP = initResult.leafWPP;
    auto &leafWKP = initResult.leafWKP;

    auto time_init = std::chrono::high_resolution_clock::now();
    std::cout << "V8b init took: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_init - time_start).count()
              << " ms, leafEdgeData entries: " << leafEdgeData.size() << std::endl;

    // leafEdgeAlive[li]: whether leaf li's CSR data is still valid
    // Set to 0 when leaf dies (Case A), or when leaf is modified (Case C surviving / Case B)
    std::vector<uint8_t> leafEdgeAlive(leafEdgeOff.size() > 0 ? leafEdgeOff.size() - 1 : 0, 1);

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
    daf::StaticVector<PlusECDSet_ST_V8b::LeafRmInfo> leafRmInfo(tree.adj_list.size());
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

    // Bucket move helper
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

    std::cout << "=========================begin=========================" << std::endl;
    double minCore = 0;
    long long cntA = 0, cntB = 0, cntC = 0, totalIters = 0;
    long long cntA_fast = 0, cntA_fallback = 0;
    long long us_phase1 = 0, us_caseA_delta = 0;
    long long us_caseB = 0, us_flush = 0;
    // Detailed work counters
    long long work_p1_iterate = 0, work_p1_found = 0;  // Phase 1: iterations & matches
    long long work_caseA_csr = 0;          // CSR entries scanned
    long long work_caseA_fallback_edges = 0; // getEdgeCompressedId calls (fallback)
    long long work_caseA_removeNbr = 0;    // removeNbr calls (Case A)
    long long work_caseC_edges = 0;        // getEdgeCompressedId calls (Case C)
    long long work_caseC_removeNbr = 0;    // removeNbr calls (Case C)
    long long work_caseC_alive = 0, work_caseC_dead = 0; // C sub-cases
    long long work_caseC_alive_hasCSR = 0, work_caseC_alive_noCSR = 0;
    long long work_caseC_alive_hasCSR_edges = 0, work_caseC_alive_noCSR_edges = 0;
    long long work_flush_dirty = 0;        // dirty edges flushed
    long long us_caseA_csr = 0, us_caseA_fallback_t = 0, us_caseA_mut = 0;
    long long us_caseC_delta_t = 0, us_caseC_mut_t = 0;

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

        // --- Phase 1: intersect edge adjacency lists (identical to ST baseline) ---
        auto _t0 = std::chrono::high_resolution_clock::now();
        for (int ei = 0; ei < (int)currentRemoveEdgeIds.size(); ++ei) {
            auto edgeId = currentRemoveEdgeIds[ei];
            auto [edgeU, edgeV] = edgeGraph.getEdgeById(edgeId);
            auto &adjEdgeU = treeGraphV.getNbr(edgeU);
            auto &adjEdgeV = treeGraphV.getNbr(edgeV);
            work_p1_iterate += std::min(adjEdgeU.size(), adjEdgeV.size());
            daf::intersect_dense_sets(adjEdgeU, adjEdgeV,
                                      [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                                          work_p1_found++;
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
                                                  lrm.removedEdges.push_back({edgeU, edgeV});
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
        // V8b: directSub accumulates deltas; bucketMove is deferred to flush
        auto _t2 = std::chrono::high_resolution_clock::now();
        caseBLeafIds.clear();
        {
            daf::StaticVector<daf::Size> tPovit, tKeepC;

            // V8b: directSub does NOT call bucketMove — just accumulates
            auto directSub = [&](daf::Size idx, double w) {
                countingKE[idx] -= w;
                if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx])
                    markDirty(idx);
            };

            for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                auto leafId = removedLeaf[leafIdIdx];
                PlusECDSet_ST_V8b::LeafRmInfo &leafRm = leafRmInfo[leafId];

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
                    // ---- Case A: V8b FAST PATH — CSR scan ----
                    cntA++;

                    if (leafId < leafEdgeAlive.size() && leafEdgeAlive[leafId]) {
                        // V8b fast path: linear scan of CSR entries
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
                                case PlusECDSet_ST_V8b::KK: w = wKK; break;
                                case PlusECDSet_ST_V8b::PP: w = wPP; break;
                                case PlusECDSet_ST_V8b::KP: w = wKP; break;
                            }
                            directSub(entry.edgeId, w);
                        }
                        leafEdgeAlive[leafId] = 0;
                    } else {
                        // Fallback: overflow leaves from Case B (no CSR data)
                        cntA_fallback++;
                        double KtoK = 0, KtoP = 0, PtoP = 0;
                        if (needPivot <= tPovit.size()) {
                            KtoK = std::llround(nCr[tPovit.size()][needPivot]);
                            for (daf::Size i = 0; i + 1 < tKeepC.size(); ++i)
                                for (daf::Size j = i + 1; j < tKeepC.size(); ++j) {
                                    work_caseA_fallback_edges++;
                                    directSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]), KtoK);
                                }
                        }
                        int needPP = int(needPivot) - 2;
                        if (0 <= needPP && needPP <= int(tPovit.size()) - 2) {
                            PtoP = std::llround(nCr[tPovit.size() - 2][needPP]);
                            for (daf::Size i = 0; i + 1 < tPovit.size(); ++i)
                                for (daf::Size j = i + 1; j < tPovit.size(); ++j) {
                                    work_caseA_fallback_edges++;
                                    directSub(edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]), PtoP);
                                }
                        }
                        int needKP = int(needPivot) - 1;
                        if (0 <= needKP && needKP <= int(tPovit.size()) - 1) {
                            KtoP = std::llround(nCr[tPovit.size() - 1][needKP]);
                            for (daf::Size i = 0; i < tKeepC.size(); ++i)
                                for (daf::Size j = 0; j < tPovit.size(); ++j) {
                                    work_caseA_fallback_edges++;
                                    directSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]), KtoP);
                                }
                        }
                    }
                    // Tree mutation (identical to baseline)
                    work_caseA_removeNbr += leaf.size();
                    for (const auto& i : leaf)
                        treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                    tree.adj_list[leafId].clear();
                    tree.recycleNode(leafId);
                } else {
                    // ---- Case C: UNCHANGED from ST baseline ----
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
                        work_caseC_alive++;
                        if (leafId < leafEdgeAlive.size() && leafEdgeAlive[leafId]) {
                            work_caseC_alive_hasCSR++;
                            work_caseC_alive_hasCSR_edges += (leafEdgeOff[leafId+1] - leafEdgeOff[leafId]);
                        } else {
                            work_caseC_alive_noCSR++;
                        }
                        daf::StaticVector<TreeGraphNode> newLeafF;
                        for (const auto& node : leaf) {
                            if (!std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), node.v))
                                newLeafF.push_back(node);
                        }
                        daf::Size p = leafRm.removedPivots.size();
                        work_caseC_edges += p*(p-1)/2 + newLeafF.size()*p + newLeafF.size()*(newLeafF.size()-1)/2;
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
                        // Tree mutation for Case C
                        work_caseC_removeNbr += leafRm.removedPivots.size();
                        for (auto removedNbr : leafRm.removedPivots)
                            treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                        tree.removeNbrs(leafId, leafRm.removedPivots);
                        // Invalidate CSR data for this leaf (pivots changed)
                        if (leafId < leafEdgeAlive.size())
                            leafEdgeAlive[leafId] = 0;
                    } else {
                        work_caseC_dead++;
                        work_caseC_edges += leaf.size()*(leaf.size()-1)/2;
                        for (daf::Size i = 0; i + 1 < leaf.size(); ++i)
                            for (daf::Size j = i + 1; j < leaf.size(); ++j) {
                                auto &u = leaf[i], &v = leaf[j];
                                double d = (!u.isPivot && !v.isPivot) ? RemovedKtoK : (u.isPivot && v.isPivot) ? RemovedPtoP : RemovedKtoP;
                                directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                            }
                        work_caseC_removeNbr += leaf.size();
                        for (const auto& i : leaf)
                            treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                        tree.adj_list[leafId].clear();
                        tree.recycleNode(leafId);
                        if (leafId < leafEdgeAlive.size())
                            leafEdgeAlive[leafId] = 0;
                    }
                }
                leafRmInfo[leafId].clear();
            }
            tPovit.free(); tKeepC.free();
        }

        // ========== Phase 2c: BK + apply for Case B (identical to ST baseline) ==========
        auto _t3 = std::chrono::high_resolution_clock::now();
        us_caseA_delta += std::chrono::duration_cast<std::chrono::microseconds>(_t3 - _t2).count();
        cntB += caseBLeafIds.size();

        for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
            auto leafId = caseBLeafIds[bi];
            PlusECDSet_ST_V8b::LeafRmInfo &leafRm = leafRmInfo[leafId];
            const auto& leaf = tree.adj_list[leafId];

            povit.clear(); keepC.clear();
            for (const auto& node : leaf) {
                if (node.isPivot) povit.push_back(node.v);
                else keepC.push_back(node.v);
            }
            daf::Size needPivot = k - keepC.size();
            double KtoK = 0, KtoP = 0, PtoP = 0;

            // Case B uses immediate addW (needs correct bucket for new sub-leaves)
            auto addW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] += w;
                if (edgeInHeap[idx])
                    markDirty(idx);
            };

            // Remove old leaf from treeGraphV
            for (const auto& leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }
            if (!leafRm.removedPivots.empty()) {
                tree.removeNbrs(leafId, leafRm.removedPivots);
            }

            // Invalidate CSR data for old leaf
            if (leafId < leafEdgeAlive.size())
                leafEdgeAlive[leafId] = 0;

            // Inline BK + apply new sub-leaves
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
                    double KtoK = 0, KtoP = 0, PtoP = 0;
                    if (np <= newPivot.size() && newKeepC.size() > 1) {
                        KtoK = std::llround(nCr[newPivot.size()][np]);
                        PlusECDSet_ST_V8b::processEdgePairs(newKeepC, KtoK, addW);
                    }
                    int nPP = int(np) - 2;
                    if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                        PtoP = std::llround(nCr[newPivot.size() - 2][nPP]);
                        PlusECDSet_ST_V8b::processEdgePairs(newPivot, PtoP, addW);
                    }
                    int nKP = int(np) - 1;
                    if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                        KtoP = std::llround(nCr[newPivot.size() - 1][nKP]);
                        PlusECDSet_ST_V8b::processEdgePairs(newKeepC, newPivot, KtoP, addW);
                    }
                    newPivot.clear(); newKeepC.clear();
                    if (newId >= leafRmInfo.size()) {
                        removedLeaf.reserve(newId * 1.5);
                        leafRmInfo.resize(newId * 1.5);
                        leafAffected.resize(newId * 1.5, 0);
                        // New sub-leaves from Case B don't have CSR data
                        // leafEdgeAlive doesn't need resize — newId >= original numLeaves
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
                PlusECDSet_ST_V8b::processEdgePairs(keepC, KtoK, removeW);
            }
            int needPP = int(needPivot) - 2;
            if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                PtoP = std::llround(nCr[povit.size() - 2][needPP]);
                PlusECDSet_ST_V8b::processEdgePairs(povit, PtoP, removeW);
            }
            int needKP = int(needPivot) - 1;
            if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                KtoP = std::llround(nCr[povit.size() - 1][needKP]);
                PlusECDSet_ST_V8b::processEdgePairs(keepC, povit, KtoP, removeW);
            }

            tree.removeNode(leafId);
            leafRmInfo[leafId].clear();
            povit.clear(); keepC.clear();
        }

        // ========== V8b: Flush dirty edges — single batch of bucketMoves ==========
        auto _t4 = std::chrono::high_resolution_clock::now();
        us_caseB += std::chrono::duration_cast<std::chrono::microseconds>(_t4 - _t3).count();

        work_flush_dirty += dirtyEdges.size();
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
    std::cout << "  Phase1(intersect): " << us_phase1/1000 << " ms" << std::endl;
    std::cout << "  Phase2(A+C delta+mut): " << us_caseA_delta/1000 << " ms" << std::endl;
    std::cout << "  CaseB(BK+apply): " << us_caseB/1000 << " ms" << std::endl;
    std::cout << "  Flush(bucketMove): " << us_flush/1000 << " ms" << std::endl;
    std::cout << "  --- Work Counters ---" << std::endl;
    std::cout << "  P1: iterate=" << work_p1_iterate << " found=" << work_p1_found
              << " hitRate=" << (work_p1_iterate > 0 ? 100.0*work_p1_found/work_p1_iterate : 0) << "%" << std::endl;
    std::cout << "  CaseA_CSR: entries=" << work_caseA_csr << std::endl;
    std::cout << "  CaseA_fallback: getEdgeId=" << work_caseA_fallback_edges << std::endl;
    std::cout << "  CaseA_removeNbr: " << work_caseA_removeNbr << std::endl;
    std::cout << "  CaseC: alive=" << work_caseC_alive << " dead=" << work_caseC_dead
              << " getEdgeId=" << work_caseC_edges << " removeNbr=" << work_caseC_removeNbr << std::endl;
    std::cout << "  CaseC_alive: hasCSR=" << work_caseC_alive_hasCSR
              << " (CSR_entries=" << work_caseC_alive_hasCSR_edges << ")"
              << " noCSR=" << work_caseC_alive_noCSR << std::endl;
    std::cout << "  Flush: dirtyEdges=" << work_flush_dirty << std::endl;

    daf::Size numCounting = 0;
    std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > sortedK;
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
