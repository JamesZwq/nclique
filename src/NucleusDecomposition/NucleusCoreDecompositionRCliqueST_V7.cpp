//
// ST V7: Relaxed Case C — Pivot-Only Conflict
//
// Key innovation: When all vertices touched by conflict r-cliques are pivots
// (regardless of conflict degree), we can avoid BK entirely. Instead:
//   1. Identify which pivots are "fully conflicted" (conflictDeg == max) → remove them
//   2. For remaining partial-conflict pivots, the leaf survives unchanged
//      (BK would produce the same leaf minus fully-conflicted pivots)
//   3. Compute support delta using leafCliqueInfo + subPivotCount
//
// This is the r≥3 generalization of r=2's Case C (removedPivots only, no removedEdges).
//
// Additionally: for Case C, the containment check uses cliqueIndex.byId() to get
// vertex IDs and checks against the removed pivot set — avoiding BK + hash map entirely.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <cstring>

#include "../BK/BronKerboschRmRClique.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace RCliqueSTv7 {

struct LeafCliqueEntry {
    daf::Size cliqueId;
    long long ncrValue;
    uint8_t subPivotCount;  // pivots in this r-clique within this leaf
};

struct DualIndex {
    std::vector<std::vector<LeafCliqueEntry>> leafCliqueInfo;
    std::vector<std::vector<daf::Size>> cliqueLeafIds;
    std::vector<long long> counting;
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

        daf::CliqueSize pivotC = 0, keepC = 0;
        for (const auto &i : leaf) {
            if (i.isPivot) pivotC++;
            else keepC++;
        }
        int needPivot = s - static_cast<int>(keepC);

        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            daf::CliqueSize subNumPivot = 0;
            for (const auto &node : rClique) if (node.isPivot) subNumPivot++;
            if ((int)subNumPivot > needPivot) return true;
            if (pivotC - subNumPivot >= 1001 || (daf::Size)needPivot - subNumPivot >= 401) return true;
            long long ncrValue = (long long)(nCr[pivotC - subNumPivot][needPivot - subNumPivot] + 0.5);
            auto id = cliqueIndex.byClique(rClique);
            if (id < nClique) {
                result.counting[id] += ncrValue;
                result.leafCliqueInfo[leafId].push_back({id, ncrValue, (uint8_t)subNumPivot});
                result.cliqueLeafIds[id].push_back(leafId);
            }
            return true;
        });
    }

    return result;
}

} // namespace RCliqueSTv7


std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V7(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {

    // ============ TIMERS ============
    long long duration_init = 0;
    long long duration_pop = 0;
    long long duration_intersect = 0;
    long long duration_bk = 0;
    long long duration_support = 0;
    long long duration_heap = 0;
    long long cntLeafDeath = 0, cntBK = 0, cntCaseC = 0;
    // ================================

    auto time_start = std::chrono::high_resolution_clock::now();

    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });

    daf::log_memory("r-clique index");

    auto dualIndex = daf::timeCount("buildDualIndex (ST_V7)", [&]() {
        return RCliqueSTv7::buildDualIndex(tree, cliqueIndex, r, s);
    });

    auto &leafCliqueInfo = dualIndex.leafCliqueInfo;
    auto &cliqueLeafIds = dualIndex.cliqueLeafIds;
    auto countingRClique = std::move(dualIndex.counting);

    const daf::Size nClique = cliqueIndex.size();
    std::vector<daf::Size> coreRClique(nClique);
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
        maxBucket = std::max(maxBucket, (int)countingRClique[i]);
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(nClique);
    std::vector<daf::Size> pos_in_bucket(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        int b = (int)countingRClique[i];
        bucket_of[i] = b;
        pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i);
    }
    int curBucket = 0;
    daf::Size remainingInHeap = nClique;

    auto bucketMove = [&](daf::Size id) {
        int newB = std::max(0, (int)countingRClique[id]);
        int oldB = bucket_of[id];
        if (newB == oldB) return;
        auto &oldVec = buckets[oldB];
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

    daf::log_memory("Other index (include bucket)");

    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);

    // Reusable buffer: set of unique conflict vertex IDs (global) for Case C
    std::vector<daf::Size> conflictPivotVertices;
    conflictPivotVertices.reserve(64);

    std::cout << "=========================begin (r≥3 ST_V7)=========================" << std::endl;
    long long minCore = 0;
    long long totalIters = 0;

    while (remainingInHeap > 0) {
        auto t_loop_start = std::chrono::high_resolution_clock::now();

        for (auto &leafId : changedLeaf)
            changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        changedLeaf.clear();
        removedRCliqueIdForLeaf.clear();
        currentRemoveRcliqueIds.clear();

        // --- Bucket pop ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((long long)curBucket, minCore);

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

        // --- Phase 2: Case A / Case C / Case B per leaf ---
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            auto leafId = changedLeaf[idx];
            auto leaf = tree.adj_list[leafId];
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

            // Build vertex→position map
            daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
            for (int i = 0; i < n; ++i)
                mapRef[leaf[i].v] = (daf::Size)i;

            // Collect unique conflict vertices and check if all are pivots
            // Use vertexConflictDeg to track which positions have conflicts
            vertexConflictDeg.assign(n, 0);
            for (auto rmId : removedR) {
                auto rClique = cliqueIndex.byId(rmId);
                for (auto v : rClique) {
                    daf::Size pos = mapRef[v];
                    if (pos < (daf::Size)n)
                        vertexConflictDeg[pos]++;
                }
            }

            // Check leaf-death and classify
            daf::Size maxRCliquePerVertex = (daf::Size)((long long)(nCr[n - 1][r - 1] + 0.5));
            bool leafDead = false;
            int forcedPivotRemove = 0;
            bool allConflictArePivots = true;
            bool allConflictFullOrZero = true; // no partial conflicts

            for (int i = 0; i < n; ++i) {
                if (vertexConflictDeg[i] > 0) {
                    if (!leaf[i].isPivot) {
                        allConflictArePivots = false;
                    }
                    if (vertexConflictDeg[i] < maxRCliquePerVertex) {
                        allConflictFullOrZero = false;
                    }
                }
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

            // === Case A: Leaf dies ===
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
                for (auto i : leaf)
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                tree.adj_list[leafId].clear();
                if (leafId < leafCliqueInfo.size())
                    leafCliqueInfo[leafId].clear();
                duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_struct).count();
                continue;
            }

            // === Case C: All conflict vertices are pivots ===
            // This means the removed r-cliques only involve pivot vertices in this leaf.
            // No keep vertex is affected → no structural split needed.
            // We can just remove the fully-conflicted pivots and update support via delta.
            //
            // Note: even if forcedPivotRemove == 0, if allConflictArePivots is true,
            // no vertices are fully forced out, meaning the conflict r-cliques don't
            // force any vertex removal — the leaf stays as-is minus the removed r-cliques.
            // In that case, we just subtract the removed r-cliques' contributions.
            if (allConflictArePivots && allConflictFullOrZero && forcedPivotRemove > 0) {
                cntCaseC++;
                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                auto t_supp = std::chrono::high_resolution_clock::now();

                int d = forcedPivotRemove;
                int newPivotC = (int)pivotC - d;

                if (leafId < leafCliqueInfo.size()) {
                    // Mark which positions are removed pivots
                    // (vertexConflictDeg[i] >= maxRCliquePerVertex && leaf[i].isPivot)

                    // Precompute new nCr values per pivot class
                    long long newNcrByClass[9] = {};
                    bool validClass[9] = {};
                    for (int p = 0; p <= (int)r && p <= 8; ++p) {
                        if (p <= newPivotC && p <= needPivot &&
                            newPivotC - p >= 0 && newPivotC - p < 1001 &&
                            needPivot - p >= 0 && needPivot - p < 401) {
                            newNcrByClass[p] = (long long)(nCr[newPivotC - p][needPivot - p] + 0.5);
                            validClass[p] = true;
                        }
                    }

                    std::vector<RCliqueSTv7::LeafCliqueEntry> newEntries;
                    newEntries.reserve(leafCliqueInfo[leafId].size());

                    for (const auto &entry : leafCliqueInfo[leafId]) {
                        // Check if this r-clique contains any removed pivot
                        auto cliqueSpan = cliqueIndex.byId(entry.cliqueId);
                        bool containsRemovedPivot = false;
                        for (auto v : cliqueSpan) {
                            daf::Size pos = mapRef[v];
                            if (pos < (daf::Size)n && vertexConflictDeg[pos] >= maxRCliquePerVertex) {
                                containsRemovedPivot = true;
                                break;
                            }
                        }

                        if (containsRemovedPivot) {
                            // LOST: subtract old contribution
                            if (rCliqueInHeap[entry.cliqueId]) {
                                countingRClique[entry.cliqueId] -= entry.ncrValue;
                                if (countingRClique[entry.cliqueId] < 0)
                                    countingRClique[entry.cliqueId] = 0;
                                bucketMove(entry.cliqueId);
                            }
                        } else {
                            // SURVIVES: delta = nCr[newPivotC - p][needPivot - p] - old
                            int p = entry.subPivotCount;
                            if (p <= 8 && validClass[p]) {
                                long long newNcr = newNcrByClass[p];
                                long long delta = newNcr - entry.ncrValue;
                                if (delta != 0 && rCliqueInHeap[entry.cliqueId]) {
                                    countingRClique[entry.cliqueId] += delta;
                                    if (countingRClique[entry.cliqueId] < 0)
                                        countingRClique[entry.cliqueId] = 0;
                                    bucketMove(entry.cliqueId);
                                }
                                newEntries.push_back({entry.cliqueId, newNcr, entry.subPivotCount});
                            } else {
                                // nCr out of range → contribution becomes 0
                                if (rCliqueInHeap[entry.cliqueId] && entry.ncrValue > 0) {
                                    countingRClique[entry.cliqueId] -= entry.ncrValue;
                                    if (countingRClique[entry.cliqueId] < 0)
                                        countingRClique[entry.cliqueId] = 0;
                                    bucketMove(entry.cliqueId);
                                }
                            }
                        }
                    }

                    leafCliqueInfo[leafId] = std::move(newEntries);
                }

                // Tree mutation: remove only the dead pivot vertices
                for (int i = 0; i < n; ++i) {
                    if (vertexConflictDeg[i] >= maxRCliquePerVertex) {
                        treeGraphV.removeNbr(leaf[i].v, static_cast<TreeGraphNode>(leafId));
                    }
                }
                // Shrink the leaf (remove dead pivots)
                {
                    auto &adjList = tree.adj_list[leafId];
                    adjList.erase(
                        std::remove_if(adjList.begin(), adjList.end(),
                            [&](const TreeGraphNode &node) {
                                daf::Size pos = mapRef[node.v];
                                return pos < (daf::Size)n && vertexConflictDeg[pos] >= maxRCliquePerVertex;
                            }),
                        adjList.end());
                }

                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();
                continue;
            }

            // === Case B: General case — BK required ===
            cntBK++;

            for (auto leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }

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

                    // Build leafCliqueInfo for new sub-leaf
                    std::vector<RCliqueSTv7::LeafCliqueEntry> newEntries;
                    if (leafId < leafCliqueInfo.size()) {
                        robin_hood::unordered_flat_map<daf::Size, bool> vertPivotMap;
                        vertPivotMap.reserve(stored.size());
                        for (const auto &v : stored) vertPivotMap[v.v] = v.isPivot;

                        for (const auto &entry : leafCliqueInfo[leafId]) {
                            auto cliqueSpan = cliqueIndex.byId(entry.cliqueId);
                            bool allIn = true;
                            daf::CliqueSize subP = 0;
                            for (auto v : cliqueSpan) {
                                auto it = vertPivotMap.find(v);
                                if (it == vertPivotMap.end()) {
                                    allIn = false;
                                    break;
                                }
                                if (it->second) subP++;
                            }
                            if (allIn && subP <= np && newPivotC - subP < 1001 && np - subP < 401) {
                                long long ncrVal = (long long)(nCr[newPivotC - subP][np - subP] + 0.5);
                                countingRClique[entry.cliqueId] += ncrVal;
                                newEntries.push_back({entry.cliqueId, ncrVal, (uint8_t)subP});
                                if (entry.cliqueId < cliqueLeafIds.size())
                                    cliqueLeafIds[entry.cliqueId].push_back(newId);
                            }
                        }
                    }

                    if (newId >= leafCliqueInfo.size())
                        leafCliqueInfo.resize(newId + 1);
                    leafCliqueInfo[newId] = std::move(newEntries);

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
              << " CaseC=" << cntCaseC
              << " BK=" << cntBK
              << " iters=" << totalIters << std::endl;

    // Build output
    std::vector<std::pair<std::vector<daf::Size>, int>> sortedK;
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
