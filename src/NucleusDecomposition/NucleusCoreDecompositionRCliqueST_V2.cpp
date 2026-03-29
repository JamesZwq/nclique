//
// ST V2: Opt 1 — Leaf-Clique Reverse Index
// Eliminates enumerateCombinations + byClique hash lookups in Support phase.
// Builds leafCliqueInfo[leafId] = vector<pair<cliqueId, ncrValue>> during Init.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <set>
#include <cstring>

#include "../BK/BronKerboschRmRClique.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace RCliqueSTv2 {

// Leaf-Clique reverse index entry: (cliqueId, nCr contribution)
struct LeafCliqueEntry {
    daf::Size cliqueId;
    double ncrValue;
};

// Build leafCliqueInfo and countingPerRClique in one pass
std::pair<std::vector<std::vector<LeafCliqueEntry>>, std::vector<double>>
buildLeafCliqueIndex(
    const DynamicGraph<TreeGraphNode> &treeGraph,
    const StaticCliqueIndex &cliqueIndex,
    const daf::CliqueSize r,
    const daf::CliqueSize s) {

    const daf::Size nClique = cliqueIndex.size();
    const daf::Size numLeaves = treeGraph.adj_list.size();

    std::vector<std::vector<LeafCliqueEntry>> leafCliqueInfo(numLeaves);
    std::vector<double> counting(nClique, 0);

    for (daf::Size leafId = 0; leafId < numLeaves; ++leafId) {
        const auto &leaf = treeGraph.adj_list[leafId];
        if (leaf.size() < r) continue;

        daf::CliqueSize pivotC = 0, keepC = 0;
        for (const auto &i : leaf) {
            if (i.isPivot) pivotC++;
            else keepC++;
        }
        int needPivot = s - static_cast<int>(keepC);

        // Reserve estimate for leafCliqueInfo
        // C(n, r) combinations per leaf
        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            daf::CliqueSize subNumPivot = 0;
            for (const auto &node : rClique) if (node.isPivot) subNumPivot++;
            double ncrValue = nCr[pivotC - subNumPivot][needPivot - subNumPivot];
            auto id = cliqueIndex.byClique(rClique);
            if (id < nClique) {
                counting[id] += ncrValue;
                leafCliqueInfo[leafId].push_back({id, ncrValue});
            }
            return true;
        });
    }

    return {std::move(leafCliqueInfo), std::move(counting)};
}

} // namespace RCliqueSTv2


std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V2(
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

    // Build leaf-clique reverse index and counting in one pass
    auto [leafCliqueInfo, countingRClique] = daf::timeCount("buildLeafCliqueIndex (ST_V2)", [&]() {
        return RCliqueSTv2::buildLeafCliqueIndex(tree, cliqueIndex, r, s);
    });

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

    // --- Bucket array (integer) ---
    int maxBucket = 0;
    for (daf::Size i = 0; i < cliqueIndex.size(); ++i)
        if (countingRClique[i] <= 5e6) maxBucket = std::max(maxBucket, (int)countingRClique[i]);
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(cliqueIndex.size());
    std::vector<daf::Size> pos_in_bucket(cliqueIndex.size());
    std::set<std::pair<double, daf::Size>> overflowSet;
    std::vector<double> overflowStored(cliqueIndex.size(), -1);
    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
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
    daf::Size remainingInHeap = cliqueIndex.size();

    // Immediate bucket move helper
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

    daf::log_memory("Other index (include bucket)");

    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    // Per-vertex conflict degree for leaf-death fast path
    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);

    std::cout << "=========================begin (r≥3 ST_V2)=========================" << std::endl;
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

        // --- Phase 1: Serial intersect (same as baseline ST) ---
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
            const auto &leaf = tree.adj_list[leafId];
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
                // === Leaf death: use leafCliqueInfo instead of enumerateCombinations ===
                cntLeafDeath++;
                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                auto t_supp = std::chrono::high_resolution_clock::now();
                // V2 optimization: direct lookup from leafCliqueInfo
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

                // Tree mutation: remove leaf entirely
                auto t_struct = std::chrono::high_resolution_clock::now();
                for (const auto &i : leaf)
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                tree.adj_list[leafId].clear();
                tree.recycleNode(leafId);
                // Clear leafCliqueInfo for this dead leaf
                if (leafId < leafCliqueInfo.size())
                    leafCliqueInfo[leafId].clear();
                duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_struct).count();
            } else {
                // === General case: BK required ===
                cntBK++;

                // Remove old leaf from treeGraphV
                for (const auto &leafV : leaf) {
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

                        // Build leafCliqueInfo for new sub-leaf
                        // We need to check which cliques from parent are fully contained
                        // Build vertex set for fast containment check
                        // Use mapRef (vListMap) - set positions for new leaf vertices
                        std::vector<RCliqueSTv2::LeafCliqueEntry> newEntries;
                        if (leafId < leafCliqueInfo.size()) {
                            newEntries.reserve(leafCliqueInfo[leafId].size());
                            // Build a set of vertices in the new leaf for fast membership test
                            robin_hood::unordered_flat_set<daf::Size> newLeafVerts;
                            newLeafVerts.reserve(stored.size());
                            for (const auto &v : stored) newLeafVerts.insert(v.v);

                            for (const auto &entry : leafCliqueInfo[leafId]) {
                                // Check if all vertices of this r-clique are in the new leaf
                                auto clique = cliqueIndex.byId(entry.cliqueId);
                                bool allIn = true;
                                for (auto v : clique) {
                                    if (newLeafVerts.find(v) == newLeafVerts.end()) {
                                        allIn = false;
                                        break;
                                    }
                                }
                                if (allIn) {
                                    // Recompute nCr for this sub-leaf context
                                    auto cliqueSpan = cliqueIndex.byId(entry.cliqueId);
                                    daf::CliqueSize subP = 0;
                                    for (auto v : cliqueSpan) {
                                        // Find if this vertex is a pivot in the new leaf
                                        for (const auto &sv : stored) {
                                            if (sv.v == v) {
                                                if (sv.isPivot) subP++;
                                                break;
                                            }
                                        }
                                    }
                                    if (subP <= np && newPivotC - subP < 1001 && np - subP < 401) {
                                        double ncrVal = nCr[newPivotC - subP][np - subP];
                                        countingRClique[entry.cliqueId] += ncrVal;
                                        newEntries.push_back({entry.cliqueId, ncrVal});
                                    }
                                }
                            }
                        }

                        // Ensure leafCliqueInfo is large enough
                        if (newId >= leafCliqueInfo.size())
                            leafCliqueInfo.resize(newId + 1);
                        leafCliqueInfo[newId] = std::move(newEntries);

                        if (newId >= changedLeafIndex.size())
                            changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                    });

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                // Subtract old leaf contribution using leafCliqueInfo
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

                // Remove old leaf node
                auto t_struct = std::chrono::high_resolution_clock::now();
                tree.removeNode(leafId);
                // Clear leafCliqueInfo for this old leaf
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
