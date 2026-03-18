//
// Single-thread optimized r=2 nucleus decomposition.
// Stripped of all OMP infrastructure: no locks, no per-thread vectors,
// no deferred dirty tracking, no atomic directives.
// Direct countingKE updates + immediate bucketMove.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>

#include "../BK/BronKerboschRmEdge.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace PlusECDSet_ST {
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

    std::pair<double *, daf::Size *> countingPerEdge(const DynamicGraph<TreeGraphNode> &treeGraph,
                                                     const Graph &edgeGraph,
                                                     const daf::CliqueSize k) {
        double *countingE = new double[edgeGraph.adj_list.size()];
        daf::Size *degreeE = new daf::Size[edgeGraph.adj_list.size()];
        memset(countingE, 0, edgeGraph.adj_list.size() * sizeof(double));
        memset(degreeE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));

        const int numLeaves = (int)treeGraph.adj_list.size();
        daf::StaticVector<daf::Size> tPovit, tKeepC;

        for (int li = 0; li < numLeaves; ++li) {
            const auto &clique = treeGraph.adj_list[li];
            if (clique.size() < k) continue;

            tPovit.clear(); tKeepC.clear();
            for (const auto &node : clique) {
                if (node.isPivot) tPovit.push_back(node.v);
                else tKeepC.push_back(node.v);
            }

            int needPivot = int(k) - int(tKeepC.size());

            // 1) keep-keep
            if (needPivot >= 0 && needPivot <= int(tPovit.size())) {
                double totalKcliques = nCr[tPovit.size()][needPivot];
                for (size_t i = 0; i < tKeepC.size(); ++i)
                    for (size_t j = i + 1; j < tKeepC.size(); ++j) {
                        auto index = edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]);
                        countingE[index] += totalKcliques;
                        degreeE[index]++;
                    }
            }

            // 2) pivot-pivot
            int needPP = needPivot - 2;
            if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
                double eachPP = nCr[tPovit.size() - 2][needPP];
                for (size_t i = 0; i < tPovit.size(); ++i)
                    for (size_t j = i + 1; j < tPovit.size(); ++j) {
                        auto index = edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]);
                        countingE[index] += eachPP;
                        degreeE[index]++;
                    }
            }

            // 3) keep-pivot cross
            int needKP = needPivot - 1;
            if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
                double eachKP = nCr[tPovit.size() - 1][needKP];
                for (size_t i = 0; i < tKeepC.size(); ++i)
                    for (size_t j = 0; j < tPovit.size(); ++j) {
                        auto index = edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]);
                        countingE[index] += eachKP;
                        degreeE[index]++;
                    }
            }
        }
        tPovit.free(); tKeepC.free();
        return {countingE, degreeE};
    }
}

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_ST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    auto time_start = std::chrono::high_resolution_clock::now();
    auto [countingKE_dbl, degreeERemove] = PlusECDSet_ST::countingPerEdge(tree, edgeGraph, k);

    // Convert double support counts to exact integer arithmetic
    const daf::Size numEdgesInit = edgeGraph.adj_list.size();
    auto *countingKE = new long long[numEdgesInit];
    for (daf::Size i = 0; i < numEdgesInit; ++i)
        countingKE[i] = std::llround(countingKE_dbl[i]);
    delete[] countingKE_dbl;

    auto *coreE = new long long[numEdgesInit];
    memset(coreE, 0, numEdgesInit * sizeof(long long));

    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keepC;
    daf::StaticVector<daf::Size> newPivot;
    daf::StaticVector<daf::Size> newKeepC;

    daf::StaticVector<daf::Size> currentRemoveEdgeIds(edgeGraph.adj_list.size());

    daf::StaticVector<bool> edgeInHeap(edgeGraph.adj_list.size());
    edgeInHeap.c_size = edgeGraph.adj_list.size();
    memset(edgeInHeap.getData(), true, edgeGraph.adj_list.size() * sizeof(bool));

    daf::StaticVector<daf::Size> removedLeaf(tree.adj_list.size());
    daf::StaticVector<PlusECDSet_ST::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();

    long long currCore = 0;

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
            edgeInHeap[i] = false;
            continue;
        }
        int b = (int)countingKE[i];
        bucket_of[i] = b;
        pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i);
        remainingInHeap++;
    }
    int curBucket = 0;

    // Bucket move helper — immediate, no deferred tracking
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

    std::vector<bool> leafAffected(leafRmInfo.size(), false);
    std::vector<daf::Size> caseBLeafIds;
    caseBLeafIds.reserve(1024);

    std::cout << "=========================begin=========================" << std::endl;
    long long minCore = 0;
    long long cntA = 0, cntB = 0, cntC = 0, totalIters = 0;

    while (remainingInHeap > 0) {
        // --- Bucket pop ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((long long)curBucket, minCore);

        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                edgeInHeap[id] = false;
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

        // --- Phase 1: intersect edge adjacency lists (serial, no locks) ---
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
                                              leafAffected[leafId] = true;
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

        // Pre-sort removedPivots for all leaves
        for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
            auto leafId = removedLeaf[leafIdIdx];
            std::ranges::sort(leafRmInfo[leafId].removedPivots);
            leafRmInfo[leafId].removedPivots.unique();
        }

        // ========== Phase 2a: delta computation for Case A & C (direct updates) ==========
        {
            daf::StaticVector<daf::Size> tPovit, tKeepC;

            auto directSub = [&](daf::Size idx, long long w) {
                countingKE[idx] -= w;
                if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx])
                    bucketMove(idx);
            };

            for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                auto leafId = removedLeaf[leafIdIdx];
                PlusECDSet_ST::LeafRmInfo &leafRm = leafRmInfo[leafId];

                auto leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;

                tPovit.clear(); tKeepC.clear();
                for (auto node : leaf) {
                    if (node.isPivot) tPovit.push_back(node.v);
                    else tKeepC.push_back(node.v);
                }
                daf::Size needPivot = k - tKeepC.size();

                bool isDeadLeaf = leafRm.removedKeepC || needPivot > tPovit.size() - leafRm.removedPivots.size();
                bool isCaseB = !leafRm.removedEdges.empty() && !isDeadLeaf;
                if (isCaseB) continue; // Case B handled in Phase 2c

                if (isDeadLeaf) {
                    // ---- Case A: leaf dies — subtract full contribution ----
                    long long KtoK = 0, KtoP = 0, PtoP = 0;
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
                } else {
                    // ---- Case C: only pivots removed ----
                    long long KtoK = 0, KtoP = 0, PtoP = 0;
                    long long RemovedKtoK = 0, RemovedKtoP = 0, RemovedPtoP = 0;

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
                        daf::StaticVector<TreeGraphNode> newLeafF;
                        for (auto node : leaf) {
                            if (!std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), node.v))
                                newLeafF.push_back(node);
                        }
                        // removed x removed
                        for (daf::Size i = 0; i + 1 < leafRm.removedPivots.size(); ++i)
                            for (daf::Size j = i + 1; j < leafRm.removedPivots.size(); ++j)
                                directSub(edgeGraph.getEdgeCompressedId(leafRm.removedPivots[i], leafRm.removedPivots[j]), RemovedPtoP);
                        // newLeaf x removed
                        for (daf::Size i = 0; i < newLeafF.size(); ++i)
                            for (daf::Size j = 0; j < leafRm.removedPivots.size(); ++j) {
                                long long d = newLeafF[i].isPivot ? RemovedPtoP : RemovedKtoP;
                                directSub(edgeGraph.getEdgeCompressedId(newLeafF[i].v, leafRm.removedPivots[j]), d);
                            }
                        // newLeaf x newLeaf
                        for (daf::Size i = 0; i + 1 < newLeafF.size(); ++i)
                            for (daf::Size j = i + 1; j < newLeafF.size(); ++j) {
                                auto &u = newLeafF[i], &v = newLeafF[j];
                                long long d = (!u.isPivot && !v.isPivot) ? KtoK : (u.isPivot && v.isPivot) ? PtoP : KtoP;
                                directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                            }
                        newLeafF.free();
                    } else {
                        // Full removal (leaf dies)
                        for (daf::Size i = 0; i + 1 < leaf.size(); ++i)
                            for (daf::Size j = i + 1; j < leaf.size(); ++j) {
                                auto &u = leaf[i], &v = leaf[j];
                                long long d = (!u.isPivot && !v.isPivot) ? RemovedKtoK : (u.isPivot && v.isPivot) ? RemovedPtoP : RemovedKtoP;
                                directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                            }
                    }
                }
            }
            tPovit.free(); tKeepC.free();
        }

        // ========== Phase 2b: tree mutations for Case A & C ==========
        {
            for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                auto leafId = removedLeaf[leafIdIdx];
                PlusECDSet_ST::LeafRmInfo &leafRm = leafRmInfo[leafId];
                auto leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;

                daf::Size numKeeps = 0;
                for (auto node : leaf) if (!node.isPivot) numKeeps++;
                daf::Size needPivot = k - numKeeps;
                daf::Size numPivots = leaf.size() - numKeeps;

                bool isDeadLeaf = leafRm.removedKeepC || needPivot > numPivots - leafRm.removedPivots.size();
                bool isCaseB = !leafRm.removedEdges.empty() && !isDeadLeaf;
                if (isCaseB) continue;

                if (isDeadLeaf) {
                    cntA++;
                    for (auto i : leaf) {
                        treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                    }
                    tree.adj_list[leafId].clear();
                    tree.recycleNode(leafId);
                } else {
                    cntC++;
                    if (!leafRm.removedPivots.empty() && needPivot <= numPivots - leafRm.removedPivots.size()) {
                        for (auto removedNbr : leafRm.removedPivots) {
                            treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                        }
                        tree.removeNbrs(leafId, leafRm.removedPivots);
                    } else {
                        for (auto i : leaf) {
                            treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                        }
                        tree.adj_list[leafId].clear();
                        tree.recycleNode(leafId);
                    }
                }
                leafRmInfo[leafId].clear();
            }
        }

        // ========== Phase 2c: BK + apply for Case B (always serial) ==========
        caseBLeafIds.clear();
        for (daf::Size leafIdIdx = 0; leafIdIdx < removedLeaf.size(); ++leafIdIdx) {
            auto leafId = removedLeaf[leafIdIdx];
            PlusECDSet_ST::LeafRmInfo &leafRm = leafRmInfo[leafId];
            if (leafRm.removedEdges.empty()) continue;
            auto &leaf = tree.adj_list[leafId];
            if (leaf.empty()) continue;
            daf::Size numKeepsB = 0;
            for (auto node : leaf) if (!node.isPivot) numKeepsB++;
            daf::Size needPivotB = k - numKeepsB;
            daf::Size numPivotsB = leaf.size() - numKeepsB;
            if (leafRm.removedKeepC || needPivotB > numPivotsB - leafRm.removedPivots.size()) continue;
            caseBLeafIds.push_back(leafId);
        }
        cntB += caseBLeafIds.size();

        // Serial BK path
        for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
            auto leafId = caseBLeafIds[bi];
            PlusECDSet_ST::LeafRmInfo &leafRm = leafRmInfo[leafId];
            auto leaf = tree.adj_list[leafId];

            povit.clear(); keepC.clear();
            for (auto node : leaf) {
                if (node.isPivot) povit.push_back(node.v);
                else keepC.push_back(node.v);
            }
            daf::Size needPivot = k - keepC.size();
            long long KtoK = 0, KtoP = 0, PtoP = 0;

            auto addW = [&](daf::Size u, daf::Size v, long long w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] += w;
                if (edgeInHeap[idx])
                    bucketMove(idx);
            };

            // Remove old leaf from treeGraphV
            for (auto leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }
            if (!leafRm.removedPivots.empty()) {
                tree.removeNbrs(leafId, leafRm.removedPivots);
            }

            // Inline BK + apply new sub-leaves
            auto &leafRef = tree.adj_list[leafId];
            bkRmEdge::bronKerbosch(leafRef, leafRm.removedEdges, k,
                [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                    auto subLeaf = bkRmEdge::coverToVertex(c, pivots, leafRef);
                    auto newId = tree.addNode(subLeaf);
                    newPivot.clear(); newKeepC.clear();
                    for (auto i : tree.adj_list[newId]) {
                        if (i.isPivot) { newPivot.push_back(i.v); treeGraphV.addNbr(i.v, {newId, true}); }
                        else { newKeepC.push_back(i.v); treeGraphV.addNbr(i.v, {newId, false}); }
                    }
                    daf::Size np = k - newKeepC.size();
                    long long KtoK = 0, KtoP = 0, PtoP = 0;
                    if (np <= newPivot.size() && newKeepC.size() > 1) {
                        KtoK = std::llround(nCr[newPivot.size()][np]);
                        PlusECDSet_ST::processEdgePairs(newKeepC, KtoK, addW);
                    }
                    int nPP = int(np) - 2;
                    if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                        PtoP = std::llround(nCr[newPivot.size() - 2][nPP]);
                        PlusECDSet_ST::processEdgePairs(newPivot, PtoP, addW);
                    }
                    int nKP = int(np) - 1;
                    if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                        KtoP = std::llround(nCr[newPivot.size() - 1][nKP]);
                        PlusECDSet_ST::processEdgePairs(newKeepC, newPivot, KtoP, addW);
                    }
                    newPivot.clear(); newKeepC.clear();
                    if (newId >= leafRmInfo.size()) {
                        removedLeaf.reserve(newId * 1.5);
                        leafRmInfo.resize(newId * 1.5);
                        leafAffected.resize(newId * 1.5, false);
                    }
                });

            // Remove old contribution
            auto removeW = [&](daf::Size u, daf::Size v, long long w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] -= w;
                if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx])
                    bucketMove(idx);
            };
            if (needPivot <= povit.size()) {
                KtoK = std::llround(nCr[povit.size()][needPivot]);
                PlusECDSet_ST::processEdgePairs(keepC, KtoK, removeW);
            }
            int needPP = int(needPivot) - 2;
            if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                PtoP = std::llround(nCr[povit.size() - 2][needPP]);
                PlusECDSet_ST::processEdgePairs(povit, PtoP, removeW);
            }
            int needKP = int(needPivot) - 1;
            if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                KtoP = std::llround(nCr[povit.size() - 1][needKP]);
                PlusECDSet_ST::processEdgePairs(keepC, povit, KtoP, removeW);
            }

            tree.removeNode(leafId);
            leafRmInfo[leafId].clear();
            povit.clear(); keepC.clear();
        }

        currentRemoveEdgeIds.clear();
        for (auto leafId : removedLeaf) leafAffected[leafId] = false;
        removedLeaf.clear();
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
    std::cout << "  Cases: A=" << cntA << " B=" << cntB << " C=" << cntC << " iters=" << totalIters << std::endl;

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
    delete[] degreeERemove;
    povit.free();
    keepC.free();
    newPivot.free();
    newKeepC.free();
    return sortedK;
}
