//
// ST V12: On-Demand r-Clique Peeling (Memory-Optimized) + Per-Leaf Cache
//
// Eliminates leafCliqueInfo entirely. During peeling, r-cliques are
// re-enumerated on-demand via enumerateCombinationsWithIdx + lookupRaw().
// Trades modestly more compute for dramatically less memory.
//
// Optimization: Per-leaf clique cache — enumerate + lookupRaw once per leaf,
// cache (cliqueId, subP) in a transient vector. BK children and subtraction
// reuse the cache via deterministic enumeration order, avoiding redundant
// lookupRaw calls. Reduces hash probes from (K+1)*C(n,r) to C(n,r) per leaf.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <set>
#include <cstring>
#include <cassert>

#include "../BK/BronKerboschRmRClique.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V12(
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

    // Build clique index — NO leafCliqueInfo
    StaticCliqueIndex cliqueIndex(r);
    const daf::Size numLeaves = tree.adj_list.size();
    const daf::Size maxV = edgeGraph.adj_list.size();

    std::vector<std::vector<daf::Size>> cliqueLeafIds;
    std::vector<double> countingRClique;

    // Pre-compute per-leaf pivotC/keepC for nCr calculation
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

    daf::timeCount("fused build (ST_V12, no leafCliqueInfo)", [&]() {
        cliqueIndex.buildWithFullEnum(tree, maxV,
            [&](daf::Size leafIdx, StaticCliqueIndex::Id cliqueId, daf::CliqueSize subNumPivot,
                const uint8_t* positions) {
                // Grow vectors if needed
                if (cliqueId >= countingRClique.size()) {
                    daf::Size newSz = cliqueId + 1;
                    countingRClique.resize(newSz, 0);
                    cliqueLeafIds.resize(newSz);
                }
                int np = leafNeedPivot[leafIdx];
                int remPivot = (int)leafPivotC[leafIdx] - (int)subNumPivot;
                int remNeed = np - (int)subNumPivot;
                if (remPivot < 0 || remNeed < 0 || remPivot < remNeed) return;
                double ncrVal = nCr[remPivot][remNeed];
                countingRClique[cliqueId] += ncrVal;

                // NO leafCliqueInfo — this is the memory saving
                cliqueLeafIds[cliqueId].push_back(leafIdx);
            });
    });

    // V12: do NOT call cliqueIndex.freeMapList() — need lookupRaw() during peeling
    daf::log_memory("r-clique index (V12, mapList_ kept)");

    const daf::Size nClique = cliqueIndex.size();
    std::vector<double> coreRClique(nClique, 0);
    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<std::vector<daf::Size>> removedRCliqueIdForLeaf;
    std::vector<daf::Size> changedLeaf;
    std::vector<daf::Size> currentRemoveRcliqueIds;

    removedRCliqueIdForLeaf.reserve(tree.adj_list.size());
    changedLeaf.reserve(tree.adj_list.size());
    currentRemoveRcliqueIds.reserve(nClique);

    daf::StaticVector<bool> rCliqueInHeap(nClique);
    rCliqueInHeap.resize(nClique);
    memset(rCliqueInHeap.getData(), true, nClique * sizeof(bool));

    // --- Bucket array ---
    int maxBucket = 0;
    for (daf::Size i = 0; i < nClique; ++i)
        if (countingRClique[i] <= 5e6) maxBucket = std::max(maxBucket, (int)countingRClique[i]);
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(nClique);
    std::vector<daf::Size> pos_in_bucket(nClique);
    std::set<std::pair<double, daf::Size>> overflowSet;
    std::vector<double> overflowStored(nClique, -1);
    for (daf::Size i = 0; i < nClique; ++i) {
        if (countingRClique[i] > 5e6) {
            overflowSet.insert({countingRClique[i], i});
            bucket_of[i] = -1;
        } else {
            int b = (int)countingRClique[i];
            bucket_of[i] = b;
            pos_in_bucket[i] = buckets[b].size();
            buckets[b].push_back(i);
        }
    }
    int curBucket = 0;
    daf::Size remainingInHeap = nClique;

    auto bucketMove = [&](daf::Size id) {
        if (!rCliqueInHeap[id]) return;
        double val = std::max(0.0, countingRClique[id]);
        int oldB = bucket_of[id];
        if (oldB == -1) overflowSet.erase({overflowStored[id], id});
        if (val <= 5e6) {
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
                    if (myPos < oldVec.size() - 1) { auto last = oldVec.back(); oldVec[myPos] = last; pos_in_bucket[last] = myPos; }
                    oldVec.pop_back();
                }
            }
            overflowSet.insert({val, id});
            overflowStored[id] = val;
            bucket_of[id] = -1;
        }
    };

    std::vector<TreeGraphNode> reusableLeaf;
    reusableLeaf.reserve(400);

    // Per-leaf clique cache: built once via lookupRaw, reused for BK children + subtraction
    // Stores posMask so BK children can do containment via linear scan (no re-enumerate)
    struct CachedEntry {
        daf::Size cliqueId;      // 4B
        daf::CliqueSize subP;    // 1B — number of pivots in this r-subset
        uint8_t pad[3];          // 3B
        uint64_t posMask;        // 8B — bitmask of leaf-local positions
    };  // 16B total
    std::vector<CachedEntry> leafCache;
    leafCache.reserve(4096);

    daf::log_memory("Other index (include bucket, V12)");

    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);

    // Reusable buffer for lookupRaw vertex array
    daf::Size vertBuf[8];

    std::cout << "=========================begin (r>=3 ST_V12)===========================" << std::endl;
    double minCore = 0;
    long long totalIters = 0;

    while (remainingInHeap > 0) {
        auto t_loop_start = std::chrono::high_resolution_clock::now();

        for (auto &leafId : changedLeaf)
            changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        changedLeaf.clear();
        removedRCliqueIdForLeaf.clear();
        currentRemoveRcliqueIds.clear();

        // --- Bucket pop ---

        // Drain overflow items that dropped below threshold
        while (!overflowSet.empty()) {
            daf::Size oid = overflowSet.begin()->second;
            if (!rCliqueInHeap[oid]) { overflowSet.erase(overflowSet.begin()); continue; }
            if (countingRClique[oid] <= 5e6) {
                overflowSet.erase(overflowSet.begin());
                int b = (int)countingRClique[oid];
                bucket_of[oid] = b;
                pos_in_bucket[oid] = buckets[b].size();
                buckets[b].push_back(oid);
            } else break;
        }
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) {
            // Process overflow set
            if (overflowSet.empty()) break;
            while (!overflowSet.empty()) {
                daf::Size oid = overflowSet.begin()->second;
                overflowSet.erase(overflowSet.begin());
                if (!rCliqueInHeap[oid]) continue;
                minCore = std::max(countingRClique[oid], minCore);
                rCliqueInHeap[oid] = false;
                currentRemoveRcliqueIds.push_back(oid);
                coreRClique[oid] = minCore;
                remainingInHeap--;
                while (!overflowSet.empty()) {
                    daf::Size nid = overflowSet.begin()->second;
                    if (!rCliqueInHeap[nid]) { overflowSet.erase(overflowSet.begin()); continue; }
                    if (countingRClique[nid] <= minCore) {
                        overflowSet.erase(overflowSet.begin());
                        rCliqueInHeap[nid] = false;
                        currentRemoveRcliqueIds.push_back(nid);
                        coreRClique[nid] = minCore;
                        remainingInHeap--;
                    } else break;
                }
                break;
            }
            goto pop_done_label;
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
        pop_done_label:
        duration_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_loop_start).count();
        if (remainingInHeap == 0) break;
        totalIters++;

        // --- Phase 1: Intersect ---
        auto t_intersect = std::chrono::high_resolution_clock::now();
        for (size_t ri = 0; ri < currentRemoveRcliqueIds.size(); ++ri) {
            auto rmRCliqueId = currentRemoveRcliqueIds[ri];
            if (ri + 1 < currentRemoveRcliqueIds.size()) {
                auto nextId = currentRemoveRcliqueIds[ri + 1];
                if (nextId < cliqueLeafIds.size())
                    __builtin_prefetch(&cliqueLeafIds[nextId], 0, 1);
            }
            if (rmRCliqueId < cliqueLeafIds.size()) {
                for (auto leafId : cliqueLeafIds[rmRCliqueId]) {
                    if (tree.adj_list[leafId].empty()) continue;
                    auto &leafIdx = changedLeafIndex[leafId];
                    if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                        leafIdx = removedRCliqueIdForLeaf.size();
                        removedRCliqueIdForLeaf.emplace_back();
                        changedLeaf.push_back(leafId);
                        removedRCliqueIdForLeaf.back().reserve(10);
                    }
                    removedRCliqueIdForLeaf[leafIdx].emplace_back(rmRCliqueId);
                }
            }
        }
        duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_intersect).count();

        // --- Phase 2: BK + support ---
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            auto leafId = changedLeaf[idx];
            const auto& leaf = tree.adj_list[leafId];
            if (leaf.empty()) continue;
            auto leafIndex = changedLeafIndex[leafId];
            auto &removedR = removedRCliqueIdForLeaf[leafIndex];

            daf::CliqueSize keepC = 0, pivotC = 0;
            for (const auto &node : leaf) {
                if (node.isPivot) pivotC++;
                else keepC++;
            }
            int needPivot = s - static_cast<int>(keepC);
            int n = (int)leaf.size();

            auto t_bk = std::chrono::high_resolution_clock::now();

            daf::Size maxRCliquePerVertex = (daf::Size)((double)(nCr[n - 1][r - 1] + 0.5));

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

            bool leafDead = false;
            int forcedPivotRemove = 0;
            for (int i = 0; i < n; ++i) {
                if (vertexConflictDeg[i] >= maxRCliquePerVertex) {
                    if (!leaf[i].isPivot) { leafDead = true; break; }
                    forcedPivotRemove++;
                }
            }
            if (!leafDead) {
                int remainingPivots = (int)pivotC - forcedPivotRemove;
                int remainingTotal = n - forcedPivotRemove;
                if (remainingTotal < (int)s || remainingPivots < needPivot)
                    leafDead = true;
            }

            if (leafDead) {
                cntLeafDeath++;
                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                // === Leaf-death: direct enumerate + lookupRaw + subtract (no cache) ===
                auto t_supp = std::chrono::high_resolution_clock::now();

                daf::enumerateCombinationsWithIdx(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &rClique, const size_t* idxArr) {
                        daf::CliqueSize subP = 0;
                        for (daf::Size j = 0; j < r; ++j) {
                            vertBuf[j] = rClique[j].v;
                            if (rClique[j].isPivot) subP++;
                        }
                        auto id = cliqueIndex.lookupRaw(vertBuf);
                        if (id >= nClique || !rCliqueInHeap[id]) return true;
                        int remP = (int)pivotC - (int)subP;
                        int remN = needPivot - (int)subP;
                        if (remP < 0 || remN < 0 || remP < remN) return true;
                        double ncrVal = nCr[remP][remN];
                        countingRClique[id] -= ncrVal;
                        if (countingRClique[id] < 0) countingRClique[id] = 0;
                        bucketMove(id);
                        return true;
                    });

                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

                auto t_struct = std::chrono::high_resolution_clock::now();
                tree.adj_list[leafId].clear();
                duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_struct).count();
            } else {
                // === BK case: cache with posMask, linear scan for children ===
                cntBK++;

                // Step 1: Build cache — enumerate once, lookupRaw + store posMask
                leafCache.clear();
                daf::enumerateCombinationsWithIdx(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &rClique, const size_t* idxArr) {
                        daf::CliqueSize subP = 0;
                        uint64_t posMask = 0;
                        for (daf::Size j = 0; j < r; ++j) {
                            vertBuf[j] = rClique[j].v;
                            if (rClique[j].isPivot) subP++;
                            posMask |= (uint64_t(1) << idxArr[j]);
                        }
                        auto id = cliqueIndex.lookupRaw(vertBuf);
                        leafCache.push_back({id, subP, {}, posMask});
                        return true;
                    });

                // Step 2: BK children — linear scan over cache (NO re-enumerate!)
                auto mapped = removedR | std::views::transform(
                    [&](const daf::Size id) { return cliqueIndex.byId(id); });

                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &bkPivots) {
                        reusableLeaf.clear();
                        daf::CliqueSize newPivotC = 0, newKeepC = 0;

                        uint64_t childMask = c.data[0];
                        uint64_t pivotMask = bkPivots.data[0];

                        c.for_each_bit([&](size_t parentPos) {
                            bool isP = bkPivots.test(parentPos);
                            reusableLeaf.push_back({leaf[parentPos].v, isP});
                            if (isP) newPivotC++; else newKeepC++;
                        });

                        auto newId = tree.addNodePresorted(reusableLeaf);
                        daf::Size np = s - newKeepC;

                        // Linear scan over cache — no enumerate, no lookupRaw
                        for (const auto &entry : leafCache) {
                            // Bitmask containment
                            if ((entry.posMask & childMask) != entry.posMask) continue;

                            // Branchless pivot count
                            daf::CliqueSize subP = (daf::CliqueSize)__builtin_popcountll(entry.posMask & pivotMask);

                            if (subP > np || newPivotC - subP >= 1001 || np - subP >= 401) continue;

                            // Add new contribution
                            double ncrVal = nCr[newPivotC - subP][np - subP];
                            countingRClique[entry.cliqueId] += ncrVal;

                            if (entry.cliqueId < cliqueLeafIds.size()) {
                                cliqueLeafIds[entry.cliqueId].push_back(newId);
                            }
                        }

                        // Update leafPivotC/leafNeedPivot for new leaf
                        if (newId >= leafPivotC.size()) {
                            leafPivotC.resize(newId + 1, 0);
                            leafNeedPivot.resize(newId + 1, 0);
                        }
                        leafPivotC[newId] = newPivotC;
                        leafNeedPivot[newId] = np;

                        if (newId >= changedLeafIndex.size())
                            changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                    });

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                // Step 3: Subtract old leaf — linear scan over cache (no re-enumeration!)
                auto t_supp = std::chrono::high_resolution_clock::now();

                for (size_t ei = 0; ei < leafCache.size(); ++ei) {
                    if (ei + 1 < leafCache.size()) {
                        __builtin_prefetch(&countingRClique[leafCache[ei + 1].cliqueId], 1, 1);
                    }
                    const auto &entry = leafCache[ei];
                    if (entry.cliqueId >= nClique || !rCliqueInHeap[entry.cliqueId]) continue;
                    int remP = (int)pivotC - (int)entry.subP;
                    int remN = needPivot - (int)entry.subP;
                    if (remP < 0 || remN < 0 || remP < remN) continue;
                    double ncrVal = nCr[remP][remN];
                    countingRClique[entry.cliqueId] -= ncrVal;
                    if (countingRClique[entry.cliqueId] < 0) countingRClique[entry.cliqueId] = 0;
                    bucketMove(entry.cliqueId);
                }

                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

                auto t_struct = std::chrono::high_resolution_clock::now();
                tree.adj_list[leafId].clear();
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
