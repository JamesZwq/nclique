//
// R=2 Lazy Verification Edge Peeling.
//
// Same batch-pop + path classification as SIGMOD baseline.
// Key difference: In Case A/C, instead of O(|P|²) edge-pair updates
// with hash lookups, we do O(|P|) vertex updates + mark affected edges
// dirty. Dirty edges get their support recomputed on-demand via
// treeGraphV intersection when they're next needed for bucket placement.
//
// countingKE[] is maintained lazily: only recomputed for dirty edges.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>

#include "../BK/BronKerboschRmEdge.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace PlusECDSet_Lazy {

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

    // Initial per-edge support computation (same as SIGMOD baseline)
    double * countingPerEdge(const DynamicGraph<TreeGraphNode> &treeGraph,
                                const Graph &edgeGraph,
                                const daf::CliqueSize k) {
        auto *countingE = new double[edgeGraph.adj_list.size()];
        memset(countingE, 0, edgeGraph.adj_list.size() * sizeof(double));

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

            if (needPivot >= 0 && needPivot <= int(tPovit.size())) {
                long long totalKcliques = std::llround(nCr[tPovit.size()][needPivot]);
                for (size_t i = 0; i < tKeepC.size(); ++i)
                    for (size_t j = i + 1; j < tKeepC.size(); ++j)
                        countingE[edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j])] += totalKcliques;
            }

            int needPP = needPivot - 2;
            if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
                long long eachPP = std::llround(nCr[tPovit.size() - 2][needPP]);
                for (size_t i = 0; i < tPovit.size(); ++i)
                    for (size_t j = i + 1; j < tPovit.size(); ++j)
                        countingE[edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j])] += eachPP;
            }

            int needKP = needPivot - 1;
            if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
                long long eachKP = std::llround(nCr[tPovit.size() - 1][needKP]);
                for (size_t i = 0; i < tKeepC.size(); ++i)
                    for (size_t j = 0; j < tPovit.size(); ++j)
                        countingE[edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j])] += eachKP;
            }
        }
        tPovit.free(); tKeepC.free();
        return countingE;
    }

    // On-demand edge support recomputation via treeGraphV intersection
    long long computeEdgeSupport(daf::Size u, daf::Size v,
                                  DynamicGraphSet<TreeGraphNode> &treeGraphV,
                                  const std::vector<int> &leafPivotCount,
                                  const std::vector<int> &leafNeedPivot,
                                  const std::vector<uint8_t> &leafAlive) {
        long long support = 0;
        auto &adjU = treeGraphV.getNbr(u);
        auto &adjV = treeGraphV.getNbr(v);
        daf::intersect_dense_sets(adjU, adjV,
            [&](const TreeGraphNode &uEntry, const TreeGraphNode &vEntry) {
                daf::Size pathId = uEntry.v;
                if (!leafAlive[pathId]) return;
                int p = leafPivotCount[pathId];
                int np = leafNeedPivot[pathId];
                if (uEntry.isPivot && vEntry.isPivot) {
                    int npp = np - 2;
                    if (npp >= 0 && npp <= p - 2)
                        support += std::llround(nCr[p - 2][npp]);
                } else if (!uEntry.isPivot && !vEntry.isPivot) {
                    if (np >= 0 && np <= p)
                        support += std::llround(nCr[p][np]);
                } else {
                    int nkp = np - 1;
                    if (nkp >= 0 && nkp <= p - 1)
                        support += std::llround(nCr[p - 1][nkp]);
                }
            });
        return support;
    }
}

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_Lazy(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    auto time_start = std::chrono::high_resolution_clock::now();
    auto *countingKE = PlusECDSet_Lazy::countingPerEdge(tree, edgeGraph, k);
    const daf::Size numEdgesInit = edgeGraph.adj_list.size();

    // Per-path counters
    const daf::Size numLeavesInit = tree.adj_list.size();
    std::vector<int> leafPivotCount(numLeavesInit);
    std::vector<int> leafNeedPivot(numLeavesInit);
    std::vector<uint8_t> leafAlive(numLeavesInit);

    for (daf::Size L = 0; L < numLeavesInit; ++L) {
        int keeps = 0, pivots = 0;
        for (auto &node : tree.adj_list[L]) {
            if (node.isPivot) pivots++;
            else keeps++;
        }
        leafPivotCount[L] = pivots;
        leafNeedPivot[L] = (int)k - keeps;
        leafAlive[L] = (leafNeedPivot[L] >= 0 && leafNeedPivot[L] <= pivots) ? 1 : 0;
    }

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
    daf::StaticVector<PlusECDSet_Lazy::LeafRmInfo> leafRmInfo(tree.adj_list.size());
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

    // --- Dirty edge tracking ---
    // When a path dies (Case A) or shrinks (Case C), we mark all edges on
    // that path as dirty. Dirty edges get their countingKE recomputed via
    // on-demand intersection before the next bucketMove.
    std::vector<uint8_t> edgeDirty(numEdges, 0);
    std::vector<daf::Size> dirtyEdges;
    dirtyEdges.reserve(65536);

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

    // Mark an edge dirty (needs support recomputation)
    auto markEdgeDirty = [&](daf::Size eid) {
        if (!edgeDirty[eid] && edgeInHeap[eid]) {
            edgeDirty[eid] = 1;
            dirtyEdges.push_back(eid);
        }
    };

    // Mark all edges on a path as dirty using O(|P|) vertex scan
    // For each vertex on the path, scan its CSR edges and mark each
    // that also has its other endpoint on the path.
    // This is O(|P| * avg_graph_degree) but avoids O(|P|²) hash lookups.
    auto markPathEdgesDirty = [&](const auto &leaf) {
        // Use a temporary set of path vertices for fast membership check
        // For small paths (<= 64), use a flat sorted array
        daf::StaticVector<daf::Size> pathVerts;
        for (auto &node : leaf) pathVerts.push_back(node.v);
        std::ranges::sort(pathVerts);

        for (auto &node : leaf) {
            daf::Size w = node.v;
            daf::Size start = edgeGraph.adj_list_offsets[w];
            daf::Size end = edgeGraph.adj_list_offsets[w + 1];
            for (daf::Size idx = start; idx < end; ++idx) {
                if (!edgeInHeap[idx]) continue;
                daf::Size nbr = edgeGraph.adj_list[idx];
                // Check if nbr is also on this path
                if (std::binary_search(pathVerts.begin(), pathVerts.end(), nbr)) {
                    markEdgeDirty(idx);
                }
            }
        }
        pathVerts.free();
    };

    std::vector<uint8_t> leafAffected(leafRmInfo.size(), 0);
    std::vector<daf::Size> caseBLeafIds;
    caseBLeafIds.reserve(1024);

    std::cout << "=========================begin=========================" << std::endl;
    double minCore = 0;
    long long cntA = 0, cntB = 0, cntC = 0, totalIters = 0;
    long long cntDirtyRecompute = 0;

    while (remainingInHeap > 0) {
        // --- Bucket pop (same as SIGMOD) ---
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

        // --- Phase 1: intersect edge adjacency lists (same as SIGMOD) ---
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

        // ========== Phase 2: LAZY delta + tree mutation for Case A & C ==========
        caseBLeafIds.clear();
        {
            for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                auto leafId = removedLeaf[leafIdIdx];
                PlusECDSet_Lazy::LeafRmInfo &leafRm = leafRmInfo[leafId];

                const auto& leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;

                daf::Size numKeeps = 0;
                daf::Size numPivots = 0;
                for (const auto& node : leaf) {
                    if (node.isPivot) numPivots++;
                    else numKeeps++;
                }
                daf::Size needPivot = k - numKeeps;

                bool isDeadLeaf = leafRm.removedKeepC || needPivot > numPivots - leafRm.removedPivots.size();
                bool isCaseB = !leafRm.removedEdges.empty() && !isDeadLeaf;

                if (isCaseB) {
                    caseBLeafIds.push_back(leafId);
                    continue;
                }

                if (isDeadLeaf) {
                    // ---- Case A: leaf dies ----
                    // Instead of O(|P|²) edge-pair directSub, mark all path edges dirty
                    cntA++;

                    markPathEdgesDirty(leaf);

                    // Update per-path state
                    leafAlive[leafId] = 0;

                    // Tree mutation (same as SIGMOD)
                    for (const auto& i : leaf)
                        treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                    tree.adj_list[leafId].clear();
                    tree.recycleNode(leafId);
                } else {
                    // ---- Case C: only pivots removed ----
                    cntC++;

                    int old_p = leafPivotCount[leafId];
                    int d = (int)leafRm.removedPivots.size();
                    int new_p = old_p - d;
                    int np = leafNeedPivot[leafId];

                    if (!leafRm.removedPivots.empty() && np <= new_p) {
                        // Mark all path edges dirty
                        markPathEdgesDirty(leaf);

                        // Update per-path counter
                        leafPivotCount[leafId] = new_p;

                        // Tree mutation for Case C
                        for (auto removedNbr : leafRm.removedPivots)
                            treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                        tree.removeNbrs(leafId, leafRm.removedPivots);
                    } else {
                        // Path dies — treat as Case A
                        markPathEdgesDirty(leaf);
                        leafAlive[leafId] = 0;
                        for (const auto& i : leaf)
                            treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                        tree.adj_list[leafId].clear();
                        tree.recycleNode(leafId);
                    }
                }
                leafRmInfo[leafId].clear();
            }
        }

        // ========== Phase 2c: BK + apply for Case B (same as SIGMOD) ==========
        cntB += caseBLeafIds.size();

        for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
            auto leafId = caseBLeafIds[bi];
            PlusECDSet_Lazy::LeafRmInfo &leafRm = leafRmInfo[leafId];
            const auto& leaf = tree.adj_list[leafId];

            povit.clear(); keepC.clear();
            for (const auto& node : leaf) {
                if (node.isPivot) povit.push_back(node.v);
                else keepC.push_back(node.v);
            }
            daf::Size needPivot = k - keepC.size();

            // Mark all path edges dirty (old contribution will be removed)
            markPathEdgesDirty(leaf);

            leafAlive[leafId] = 0;

            // Remove old leaf from treeGraphV
            for (const auto& leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }
            if (!leafRm.removedPivots.empty()) {
                tree.removeNbrs(leafId, leafRm.removedPivots);
            }

            // BK + apply new sub-leaves
            auto &leafRef = tree.adj_list[leafId];
            bkRmEdge::bronKerbosch(leafRef, leafRm.removedEdges, k,
                [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                    auto subLeaf = bkRmEdge::coverToVertex(c, pivots, leafRef);
                    auto newId = tree.addNode(subLeaf);
                    newPivot.clear(); newKeepC.clear();
                    int newPivotC = 0, newKeepC_count = 0;
                    for (const auto& i : tree.adj_list[newId]) {
                        if (i.isPivot) { newPivotC++; newPivot.push_back(i.v); treeGraphV.addNbr(i.v, {newId, true}); }
                        else { newKeepC_count++; newKeepC.push_back(i.v); treeGraphV.addNbr(i.v, {newId, false}); }
                    }
                    int newNP = (int)k - newKeepC_count;

                    // Grow per-path counter arrays if needed
                    if (newId >= leafPivotCount.size()) {
                        size_t newSize = (size_t)(newId * 1.5) + 1;
                        leafPivotCount.resize(newSize, 0);
                        leafNeedPivot.resize(newSize, 0);
                        leafAlive.resize(newSize, 0);
                    }
                    leafPivotCount[newId] = newPivotC;
                    leafNeedPivot[newId] = newNP;
                    leafAlive[newId] = (newNP >= 0 && newNP <= newPivotC) ? 1 : 0;

                    // Mark new sub-path edges dirty too (they get new contributions)
                    markPathEdgesDirty(tree.adj_list[newId]);

                    newPivot.clear(); newKeepC.clear();
                    if (newId >= leafRmInfo.size()) {
                        removedLeaf.reserve(newId * 1.5);
                        leafRmInfo.resize(newId * 1.5);
                        leafAffected.resize(newId * 1.5, 0);
                    }
                });

            tree.removeNode(leafId);
            leafRmInfo[leafId].clear();
            povit.clear(); keepC.clear();
        }

        // ========== Phase 3: Recompute dirty edges + bucketMove ==========
        for (auto eid : dirtyEdges) {
            if (!edgeInHeap[eid]) { edgeDirty[eid] = 0; continue; }
            auto [u, v] = edgeGraph.getEdgeById(eid);
            long long newSupport = PlusECDSet_Lazy::computeEdgeSupport(
                u, v, treeGraphV, leafPivotCount, leafNeedPivot, leafAlive);
            countingKE[eid] = newSupport;
            if (countingKE[eid] < 0) countingKE[eid] = 0;
            bucketMove(eid);
            edgeDirty[eid] = 0;
            cntDirtyRecompute++;
        }
        dirtyEdges.clear();

        currentRemoveEdgeIds.clear();
        for (auto leafId : removedLeaf) leafAffected[leafId] = 0;
        removedLeaf.clear();
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
    std::cout << "  Cases: A=" << cntA << " B=" << cntB << " C=" << cntC << " iters=" << totalIters << std::endl;
    std::cout << "  DirtyRecompute=" << cntDirtyRecompute << std::endl;

    daf::Size numCounting = 0;
    std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > sortedK;
    sortedK.reserve(edgeGraph.adj_list.size());

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
