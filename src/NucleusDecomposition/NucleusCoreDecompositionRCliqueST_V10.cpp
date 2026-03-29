//
// ST V10: BK Callback Micro-Optimizations
// Opt E: Reuse newEntries vector in BK callback (eliminates per-callback heap alloc)
// Opt F: addNodePresorted in BK callback (skip sort on already-sorted data)
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

namespace RCliqueSTv10 {

struct LeafCliqueEntry {
    daf::Size cliqueId;
    double ncrValue;
    uint8_t positions[8]; // leaf-local positions of r-clique vertices
};

struct DualIndex {
    std::vector<std::vector<LeafCliqueEntry>> leafCliqueInfo;
    std::vector<std::vector<daf::Size>> cliqueLeafIds;
    std::vector<double> counting;
};

DualIndex buildDualIndex(
    const DynamicGraph<TreeGraphNode> &treeGraph,
    const StaticCliqueIndex &cliqueIndex,
    const daf::CliqueSize r,
    const daf::CliqueSize s) {

    const daf::Size nClique = cliqueIndex.size();
    const daf::Size numLeaves = treeGraph.adj_list.size();

    DualIndex result;
    result.leafCliqueInfo.resize(numLeaves);
    result.cliqueLeafIds.resize(nClique);
    result.counting.assign(nClique, 0);

    for (daf::Size leafId = 0; leafId < numLeaves; ++leafId) {
        const auto &leaf = treeGraph.adj_list[leafId];
        if (leaf.size() < r) continue;
        assert(leaf.size() <= 255 && "leaf size must fit in uint8_t positions");

        daf::CliqueSize pivotC = 0, keepC = 0;
        for (const auto &i : leaf) {
            if (i.isPivot) pivotC++;
            else keepC++;
        }
        int needPivot = s - static_cast<int>(keepC);

        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            daf::CliqueSize subNumPivot = 0;
            for (const auto &node : rClique) if (node.isPivot) subNumPivot++;
            double ncrValue = nCr[pivotC - subNumPivot][needPivot - subNumPivot];
            auto id = cliqueIndex.byClique(rClique);
            if (id < nClique) {
                result.counting[id] += ncrValue;

                LeafCliqueEntry entry;
                entry.cliqueId = id;
                entry.ncrValue = ncrValue;
                for (int j = 0, li = 0; j < (int)r; ++j) {
                    while (leaf[li].v != rClique[j].v) li++;
                    entry.positions[j] = (uint8_t)li;
                    li++;
                }

                result.leafCliqueInfo[leafId].push_back(entry);
                result.cliqueLeafIds[id].push_back(leafId);
            }
            return true;
        });
    }

    return result;
}

} // namespace RCliqueSTv10


std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V10(
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

    StaticCliqueIndex localIndex(r);
    StaticCliqueIndex &cliqueIndex = prebuiltIndex ? *prebuiltIndex : localIndex;
    daf::timeCount("clique Index build", [&]() {
        if (!prebuiltIndex) localIndex.build(tree, edgeGraph.adj_list.size());
    });

    daf::log_memory("r-clique index");

    auto dualIndex = daf::timeCount("buildDualIndex (ST_V10)", [&]() {
        return RCliqueSTv10::buildDualIndex(tree, cliqueIndex, r, s);
    });

    auto &leafCliqueInfo = dualIndex.leafCliqueInfo;
    auto &cliqueLeafIds = dualIndex.cliqueLeafIds;
    auto countingRClique = std::move(dualIndex.counting);

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

    // --- Bucket array (integer) ---
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

    // Opt C: Pre-allocate reusable vector for BK callback
    std::vector<TreeGraphNode> reusableLeaf;
    reusableLeaf.reserve(400);

    // Opt E: Pre-allocate reusable vector for newEntries (eliminates per-callback heap alloc)
    std::vector<RCliqueSTv10::LeafCliqueEntry> reusableEntries;
    reusableEntries.reserve(512);

    daf::log_memory("Other index (include bucket)");

    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);

    std::cout << "=========================begin (r>=3 ST_V10)=========================" << std::endl;
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

        // --- Phase 1: Intersect via cliqueLeafIds lookup ---
        auto t_intersect = std::chrono::high_resolution_clock::now();
        for (auto rmRCliqueId : currentRemoveRcliqueIds) {
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

        // --- Phase 2: BK + tree mutation + support per leaf ---
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            auto leafId = changedLeaf[idx];
            const auto &leaf = tree.adj_list[leafId];
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

            // === Leaf-death fast path ===
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
                    if (!leaf[i].isPivot) {
                        leafDead = true;
                        break;
                    }
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

                auto t_supp = std::chrono::high_resolution_clock::now();
                if (leafId < leafCliqueInfo.size()) {
                    for (const auto &entry : leafCliqueInfo[leafId]) {
                        if (!rCliqueInHeap[entry.cliqueId]) continue;
                        countingRClique[entry.cliqueId] -= entry.ncrValue;
                        if (countingRClique[entry.cliqueId] < 0) countingRClique[entry.cliqueId] = 0;
                        bucketMove(entry.cliqueId);
                    }
                }
                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

                auto t_struct = std::chrono::high_resolution_clock::now();
                tree.adj_list[leafId].clear();
                if (leafId < leafCliqueInfo.size())
                    leafCliqueInfo[leafId].clear();
                duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_struct).count();
            } else {
                // === General case: BK required ===
                cntBK++;

                auto mapped = removedR | std::views::transform(
                    [&](const daf::Size id) { return cliqueIndex.byId(id); });

                // removeRClique sets mapRef[leaf[i].v] = i (parent-leaf positions)
                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                        // Opt C: Reuse pre-allocated vector instead of coverToVertex heap alloc
                        reusableLeaf.clear();
                        daf::CliqueSize newPivotC = 0, newKeepC = 0;
                        uint8_t parentToSubPos[400];
                        uint8_t subIdx = 0;
                        c.for_each_bit([&](size_t parentPos) {
                            bool isP = pivots.test(parentPos);
                            reusableLeaf.push_back({leaf[parentPos].v, isP});
                            parentToSubPos[parentPos] = subIdx++;
                            if (isP) newPivotC++; else newKeepC++;
                        });

                        // Opt F: addNodePresorted — for_each_bit iterates ascending, so output is sorted
                        auto newId = tree.addNodePresorted(reusableLeaf);
                        daf::Size np = s - newKeepC;

                        // Opt E: Reuse pre-allocated newEntries vector
                        reusableEntries.clear();
                        if (leafId < leafCliqueInfo.size()) {
                            for (const auto &entry : leafCliqueInfo[leafId]) {
                                bool allIn = true;
                                daf::CliqueSize subP = 0;
                                for (int j = 0; j < (int)r; ++j) {
                                    if (!c.test(entry.positions[j])) {
                                        allIn = false;
                                        break;
                                    }
                                    if (pivots.test(entry.positions[j])) subP++;
                                }
                                if (allIn && subP <= np && newPivotC - subP < 1001 && np - subP < 401) {
                                    double ncrVal = nCr[newPivotC - subP][np - subP];
                                    countingRClique[entry.cliqueId] += ncrVal;

                                    // Build new entry with remapped positions
                                    RCliqueSTv10::LeafCliqueEntry newEntry;
                                    newEntry.cliqueId = entry.cliqueId;
                                    newEntry.ncrValue = ncrVal;
                                    for (int j = 0; j < (int)r; ++j) {
                                        newEntry.positions[j] = parentToSubPos[entry.positions[j]];
                                    }
                                    reusableEntries.push_back(newEntry);

                                    if (entry.cliqueId < cliqueLeafIds.size()) {
                                        cliqueLeafIds[entry.cliqueId].push_back(newId);
                                    }
                                }
                            }
                        }

                        if (newId >= leafCliqueInfo.size())
                            leafCliqueInfo.resize(newId + 1);
                        // Swap: gives reusableEntries the old (empty) storage back
                        leafCliqueInfo[newId].swap(reusableEntries);

                        if (newId >= changedLeafIndex.size())
                            changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                    });

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                // Subtract old leaf contribution
                auto t_supp = std::chrono::high_resolution_clock::now();
                if (leafId < leafCliqueInfo.size()) {
                    for (const auto &entry : leafCliqueInfo[leafId]) {
                        if (!rCliqueInHeap[entry.cliqueId]) continue;
                        countingRClique[entry.cliqueId] -= entry.ncrValue;
                        if (countingRClique[entry.cliqueId] < 0) countingRClique[entry.cliqueId] = 0;
                        bucketMove(entry.cliqueId);
                    }
                }
                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

                auto t_struct = std::chrono::high_resolution_clock::now();
                tree.adj_list[leafId].clear();
                if (leafId < leafCliqueInfo.size())
                    leafCliqueInfo[leafId].clear();
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
