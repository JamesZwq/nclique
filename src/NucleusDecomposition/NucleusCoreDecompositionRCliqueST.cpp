//
// Single-thread optimized r≥3 nucleus decomposition.
// Optimizations over parallel NucleusCoreDecompositionRClique:
//   - No OMP: no thread_pairs, no parallel merge/sort, no atomic directives
//   - Integer arithmetic (long long countingRClique)
//   - Merged phases: intersect → BK → tree mutation → support update in single pass
//   - Immediate bucketMove (no deferred dirty tracking)
//   - Leaf-death fast path: skip BK when leaf can't possibly survive
//     (uses BK's own forced-removal criterion: vertex in ALL its r-cliques)
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <cstring>
#include <set>

#include "../BK/BronKerboschRmRClique.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace RCliqueST {

// Serial counting with integer arithmetic
std::vector<double> countingPerRClique(
    const DynamicGraph<TreeGraphNode> &treeGraph,
    const StaticCliqueIndex &cliqueIndex,
    const daf::CliqueSize r,
    const daf::CliqueSize s) {
    const daf::Size nClique = cliqueIndex.size();
    std::vector<double> counting(nClique, 0);
    for (const auto &leaf : treeGraph.adj_list) {
        if (leaf.size() < r) continue;
        daf::CliqueSize pivotC = 0, keepC = 0;
        for (const auto &i : leaf) {
            if (i.isPivot) pivotC++;
            else keepC++;
        }
        int needPivot = s - static_cast<int>(keepC);
        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            daf::CliqueSize subNumPovit = 0;
            for (const auto &node : rClique) if (node.isPivot) subNumPovit++;
            double ncrValue = nCr[pivotC - subNumPovit][needPivot - subNumPovit];
            auto id = cliqueIndex.byClique(rClique);
            if (id < nClique) counting[id] += ncrValue;
            return true;
        });
    }
    return counting;
}

} // namespace RCliqueST


std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    // ============ TIMERS ============
    long long duration_init = 0;
    long long duration_pop = 0;
    long long duration_intersect = 0;
    long long duration_bk = 0;
    long long duration_support = 0;
    long long duration_heap = 0;
    long long cntLeafDeath = 0, cntBK = 0;
    // ================================

    auto time_start = std::chrono::high_resolution_clock::now();

    StaticCliqueIndex cliqueIndex(r);
    const daf::Size numLeaves = tree.adj_list.size();
    const daf::Size maxV = edgeGraph.adj_list.size();

    // Pre-compute per-leaf pivot/keep counts
    std::vector<daf::CliqueSize> leafPivotC(numLeaves, 0);
    std::vector<int> leafNeedPivot(numLeaves, 0);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        const auto &leaf = tree.adj_list[L];
        daf::CliqueSize pC = 0, kC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pC++; else kC++;
        }
        leafPivotC[L] = pC;
        leafNeedPivot[L] = s - static_cast<int>(kC);
    }

    std::vector<double> countingRClique;
    daf::timeCount("fused build+counting (ST)", [&]() {
        cliqueIndex.buildWithFullEnum(tree, maxV,
            [&](daf::Size leafIdx, StaticCliqueIndex::Id cliqueId, daf::CliqueSize subNumPivot,
                const uint8_t* positions) {
                if (cliqueId >= countingRClique.size())
                    countingRClique.resize(cliqueId + 1, 0);
                int np = leafNeedPivot[leafIdx];
                int remPivot = (int)leafPivotC[leafIdx] - (int)subNumPivot;
                int remNeed = np - (int)subNumPivot;
                if (remPivot < 0 || remNeed < 0 || remPivot < remNeed) return;
                countingRClique[cliqueId] += nCr[remPivot][remNeed];
            });
    });

    daf::log_memory("r-clique index");

    std::vector<double> coreRClique(countingRClique.size(), 0);
    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<std::vector<daf::Size>> removedRCliqueIdForLeaf;
    std::vector<daf::Size> changedLeaf;
    std::vector<daf::Size> currentRemoveRcliqueIds;

    removedRCliqueIdForLeaf.reserve(tree.adj_list.size());
    changedLeaf.reserve(tree.adj_list.size());
    currentRemoveRcliqueIds.reserve(cliqueIndex.size());

    daf::StaticVector<bool> rCliqueInHeap(cliqueIndex.size());
    rCliqueInHeap.resize(cliqueIndex.size());
    memset(rCliqueInHeap.getData(), true, cliqueIndex.size() * sizeof(bool));

    const daf::Size graphN = edgeGraph.n;

    // --- Hybrid Bucket + Set priority queue ---
    constexpr double BUCKET_THRESHOLD = 5000000;
    double rawMaxBucket = 0;
    for (daf::Size i = 0; i < cliqueIndex.size(); ++i)
        rawMaxBucket = std::max(rawMaxBucket, countingRClique[i]);
    int maxBucket = (int)std::min(rawMaxBucket, BUCKET_THRESHOLD);

    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::set<std::pair<double, daf::Size>> overflowSet;
    std::vector<int> bucket_of(cliqueIndex.size(), -1);
    std::vector<daf::Size> pos_in_bucket(cliqueIndex.size());
    std::vector<double> overflowStoredVal(cliqueIndex.size(), -1);

    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
        if (countingRClique[i] <= BUCKET_THRESHOLD) {
            int b = (int)countingRClique[i];
            bucket_of[i] = b;
            pos_in_bucket[i] = buckets[b].size();
            buckets[b].push_back(i);
        } else {
            overflowSet.insert({countingRClique[i], i});
            overflowStoredVal[i] = countingRClique[i];
        }
    }
    int curBucket = 0;
    daf::Size remainingInHeap = cliqueIndex.size();

    auto bucketMove = [&](daf::Size id) {
        if (!rCliqueInHeap[id]) return;
        double val = std::max(0.0, countingRClique[id]);
        int oldB = bucket_of[id];

        if (oldB == -1) {
            overflowSet.erase({overflowStoredVal[id], id});
        }

        if (val <= BUCKET_THRESHOLD) {
            int newB = (int)val;
            if (oldB >= 0 && newB == oldB) return;
            if (oldB >= 0) {
                auto &oldVec = buckets[oldB];
                daf::Size myPos = pos_in_bucket[id];
                if (myPos < oldVec.size() - 1) {
                    daf::Size last = oldVec.back();
                    oldVec[myPos] = last;
                    pos_in_bucket[last] = myPos;
                }
                oldVec.pop_back();
            }
            bucket_of[id] = newB;
            pos_in_bucket[id] = buckets[newB].size();
            buckets[newB].push_back(id);
            if (newB < curBucket) curBucket = newB;
        } else {
            if (oldB >= 0) {
                auto &oldVec = buckets[oldB];
                daf::Size myPos = pos_in_bucket[id];
                if (myPos < oldVec.size()) {
                    if (myPos < oldVec.size() - 1) {
                        daf::Size last = oldVec.back();
                        oldVec[myPos] = last;
                        pos_in_bucket[last] = myPos;
                    }
                    oldVec.pop_back();
                }
            }
            overflowSet.insert({val, id});
            overflowStoredVal[id] = val;
            bucket_of[id] = -1;
        }
    };

    auto drainOverflowToBucket = [&]() {
        while (!overflowSet.empty()) {
            daf::Size id = overflowSet.begin()->second;
            if (!rCliqueInHeap[id]) { overflowSet.erase(overflowSet.begin()); continue; }
            double val = countingRClique[id];
            if (val <= BUCKET_THRESHOLD) {
                overflowSet.erase(overflowSet.begin());
                int b = (int)val;
                bucket_of[id] = b;
                pos_in_bucket[id] = buckets[b].size();
                buckets[b].push_back(id);
            } else {
                break;
            }
        }
    };

    daf::log_memory("Other index (include bucket)");

    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    // Per-vertex conflict degree for leaf-death fast path
    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);

    std::cout << "=========================begin (r≥3 ST)=========================" << std::endl;
    double minCore = 0;
    long long totalIters = 0;

    while (remainingInHeap > 0) {
        auto t_loop_start = std::chrono::high_resolution_clock::now();

        for (auto &leafId : changedLeaf)
            changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        changedLeaf.clear();
        removedRCliqueIdForLeaf.clear();
        currentRemoveRcliqueIds.clear();

        // --- Drain overflow → buckets ---
        drainOverflowToBucket();

        // --- Bucket pop ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) {
            // All buckets empty — drain overflow set
            if (!overflowSet.empty()) {
                while (!overflowSet.empty()) {
                    daf::Size id = overflowSet.begin()->second;
                    overflowSet.erase(overflowSet.begin());
                    if (!rCliqueInHeap[id]) continue;
                    double val = countingRClique[id];
                    minCore = std::max(val, minCore);
                    rCliqueInHeap[id] = false;
                    currentRemoveRcliqueIds.push_back(id);
                    coreRClique[id] = minCore;
                    remainingInHeap--;
                    while (!overflowSet.empty()) {
                        daf::Size next = overflowSet.begin()->second;
                        if (!rCliqueInHeap[next]) { overflowSet.erase(overflowSet.begin()); continue; }
                        if (countingRClique[next] <= minCore) {
                            overflowSet.erase(overflowSet.begin());
                            rCliqueInHeap[next] = false;
                            currentRemoveRcliqueIds.push_back(next);
                            coreRClique[next] = minCore;
                            remainingInHeap--;
                        } else break;
                    }
                    break;
                }
                if (currentRemoveRcliqueIds.empty()) break;
                goto phase1_st;
            }
            break;
        }

        minCore = std::max((double)curBucket, minCore);

        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                rCliqueInHeap[id] = false;
                currentRemoveRcliqueIds.push_back(id);
                coreRClique[id] = minCore;
                remainingInHeap--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (curBucket + 1) <= (int)minCore) {
                curBucket++;
            } else break;
        }
        phase1_st:

        duration_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_loop_start).count();

        if (remainingInHeap == 0) break;
        totalIters++;

        // --- Phase 1: Serial intersect (no parallel overhead) ---
        auto t_intersect = std::chrono::high_resolution_clock::now();
        for (auto rmRCliqueId : currentRemoveRcliqueIds) {
            auto rClique = cliqueIndex.byId(rmRCliqueId);
            daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                [&](const TreeGraphNode &uClique) {
                    auto &leafId = changedLeafIndex[uClique.v];
                    if (leafId == std::numeric_limits<daf::Size>::max()) {
                        leafId = removedRCliqueIdForLeaf.size();
                        removedRCliqueIdForLeaf.emplace_back();
                        changedLeaf.push_back(uClique.v);
                        removedRCliqueIdForLeaf.back().reserve(10);
                    }
                    removedRCliqueIdForLeaf[leafId].emplace_back(rmRCliqueId);
                });
        }
        duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_intersect).count();

        // --- Phase 2: Merged BK + tree mutation + support update per leaf ---
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            auto leafId = changedLeaf[idx];
            const auto& leaf = tree.adj_list[leafId];
            if (leaf.empty()) continue;
            auto leafIndex = changedLeafIndex[leafId];
            auto &removedR = removedRCliqueIdForLeaf[leafIndex];

            // Classify leaf: count keeps/pivots
            daf::CliqueSize keepC = 0, pivotC = 0;
            for (const auto &node : leaf) {
                if (node.isPivot) pivotC++;
                else keepC++;
            }
            int needPivot = s - static_cast<int>(keepC);
            int n = (int)leaf.size();

            // === Leaf-death fast path ===
            // Use BK's own criterion: a vertex is forced-removed if it appears in
            // ALL possible r-cliques containing it within this leaf.
            // maxRClique(v) = C(n-1, r-1) for each vertex in leaf of size n.
            // If conflictDeg(v) >= maxRClique(v), vertex v is forced out.
            // If a forced-out vertex is a keep → leaf dies immediately.
            // If too many pivots are forced out → leaf dies.
            auto t_bk = std::chrono::high_resolution_clock::now();

            daf::Size maxRCliquePerVertex = (daf::Size)(nCr[n - 1][r - 1] + 0.5);

            // Build vertex index and count conflict degree
            vertexConflictDeg.assign(n, 0);
            daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
            for (int i = 0; i < n; ++i)
                mapRef[leaf[i].v] = (daf::Size)i;

            for (auto rmId : removedR) {
                auto rClique = cliqueIndex.byId(rmId);
                for (auto v : rClique) {
                    daf::Size pos = mapRef[v];
                    if (pos < (daf::Size)n)
                        vertexConflictDeg[pos]++;
                }
            }

            // Check for forced removals
            bool leafDead = false;
            int forcedPivotRemove = 0;
            for (int i = 0; i < n; ++i) {
                if (vertexConflictDeg[i] >= maxRCliquePerVertex) {
                    if (!leaf[i].isPivot) {
                        // Keep vertex forced out → leaf dies
                        leafDead = true;
                        break;
                    }
                    forcedPivotRemove++;
                }
            }
            if (!leafDead) {
                // Check if remaining pivots can still satisfy needPivot
                int remainingPivots = (int)pivotC - forcedPivotRemove;
                int remainingTotal = n - forcedPivotRemove;
                if (remainingTotal < (int)s || remainingPivots < needPivot)
                    leafDead = true;
            }

            if (leafDead) {
                // === Leaf death: subtract full contribution, skip BK ===
                cntLeafDeath++;
                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                auto t_supp = std::chrono::high_resolution_clock::now();
                daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &clique) {
                    auto cliqueIndexId = cliqueIndex.byClique(clique);
                    if (!rCliqueInHeap[cliqueIndexId]) return true;
                    daf::CliqueSize subP = 0;
                    for (const auto &node : clique) if (node.isPivot) subP++;
                    double ncrValue = nCr[pivotC - subP][needPivot - subP];
                    countingRClique[cliqueIndexId] -= ncrValue;
                    if (countingRClique[cliqueIndexId] < 0) countingRClique[cliqueIndexId] = 0;
                    bucketMove(cliqueIndexId);
                    return true;
                });
                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

                // Tree mutation: remove leaf entirely
                auto t_struct = std::chrono::high_resolution_clock::now();
                for (const auto& i : leaf)
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                tree.adj_list[leafId].clear();
                tree.recycleNode(leafId);
                duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_struct).count();
            } else {
                // === General case: BK required ===
                cntBK++;

                // Remove old leaf from treeGraphV
                for (const auto& leafV : leaf) {
                    if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                    else treeGraphV.removeNbr(leafV.v, {leafId, false});
                }

                // Run BK to find surviving sub-leaves
                auto mapped = removedR | std::views::transform(
                    [&](const daf::Size id) { return cliqueIndex.byId(id); });

                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                        auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                        auto newId = tree.addNode(newLeaf);
                        const auto &stored = tree.adj_list[newId];
                        daf::CliqueSize newPivotC = 0, newKeepC = 0;
                        for (const auto &i : stored) {
                            if (i.isPivot) { treeGraphV.addNbr(i.v, {newId, true}); newPivotC++; }
                            else { treeGraphV.addNbr(i.v, {newId, false}); newKeepC++; }
                        }
                        daf::Size np = s - newKeepC;
                        daf::enumerateCombinations(stored, r, [&](const daf::StaticVector<TreeGraphNode> &rclique) {
                            daf::CliqueSize subP = 0;
                            for (const auto &node : rclique) if (node.isPivot) subP++;
                            if (subP <= np && newPivotC - subP < 1001 && np - subP < 401) {
                                double ncrVal = nCr[newPivotC - subP][np - subP];
                                countingRClique[cliqueIndex.byClique(rclique)] += ncrVal;
                            }
                            return true;
                        });
                        if (newId >= changedLeafIndex.size())
                            changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                    });

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                // Subtract old leaf contribution
                auto t_supp = std::chrono::high_resolution_clock::now();
                daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &clique) {
                    auto cliqueIndexId = cliqueIndex.byClique(clique);
                    if (!rCliqueInHeap[cliqueIndexId]) return true;
                    daf::CliqueSize subP = 0;
                    for (const auto &node : clique) if (node.isPivot) subP++;
                    double ncrVal = nCr[pivotC - subP][needPivot - subP];
                    countingRClique[cliqueIndexId] -= ncrVal;
                    if (countingRClique[cliqueIndexId] < 0) countingRClique[cliqueIndexId] = 0;
                    bucketMove(cliqueIndexId);
                    return true;
                });
                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

                // Remove old leaf node
                auto t_struct = std::chrono::high_resolution_clock::now();
                tree.removeNode(leafId);
                duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_struct).count();
            }
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::cout << "time: " << elapsed << " ms" << std::endl;
    std::cout << "Time Breakdown (ms):" << std::endl;
    std::cout << "  Init:      " << duration_init / 1000000.0 << std::endl;
    std::cout << "  Pop:       " << duration_pop / 1000000.0 << std::endl;
    std::cout << "  Intersect: " << duration_intersect / 1000000.0 << std::endl;
    std::cout << "  BK:        " << duration_bk / 1000000.0 << std::endl;
    std::cout << "  Support:   " << duration_support / 1000000.0 << std::endl;
    std::cout << "  Heap:      " << duration_heap / 1000000.0 << std::endl;
    std::cout << "  Cases: LeafDeath=" << cntLeafDeath
              << " BK=" << cntBK
              << " iters=" << totalIters << std::endl;

    // Build output
    std::vector<std::pair<std::vector<daf::Size>, double>> sortedK;
    sortedK.reserve(countingRClique.size());
    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueCopy(clique.begin(), clique.end());
        sortedK.emplace_back(cliqueCopy, coreRClique[i]);
    }
    std::sort(sortedK.begin(), sortedK.end(),
              [](const auto &a, const auto &b) { return a.second < b.second; });
    return sortedK;
}
