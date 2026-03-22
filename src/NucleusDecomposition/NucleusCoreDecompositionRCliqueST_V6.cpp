//
// ST V6: Case C Extraction — r=2-style delta formula for r≥3
//
// Key innovation over V4:
//   When all conflict vertices in a leaf are pivots with full conflict degree
//   (conflictDeg == maxRCliquePerVertex), the leaf doesn't need BK.
//   Instead, remove those pivot vertices and compute support delta using
//   nCr[old] - nCr[new] formula, just like r=2's Case C.
//
//   For the remaining "true BK" cases, use V4's approach with positional
//   containment from V5 in the BK callback.
//
// Data structure: LeafCliqueEntry stores subPivotCount to enable delta formula.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <cstring>

#include "../BK/BronKerboschRmRClique.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace RCliqueSTv6 {

struct LeafCliqueEntry {
    daf::Size cliqueId;
    long long ncrValue;
    uint8_t subPivotCount;  // how many vertices in this r-clique are pivots in this leaf
};

struct DualIndex {
    std::vector<std::vector<LeafCliqueEntry>> leafCliqueInfo;
    std::vector<std::vector<daf::Size>> cliqueLeafIds;
    std::vector<long long> counting;
    // Per-clique: how many of its r vertices are pivots in the *global* vertex set?
    // Not needed — subPivotCount is per (clique, leaf) pair, stored in LeafCliqueEntry.
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

} // namespace RCliqueSTv6


std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V6(
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

    // Build both indices and counting in one pass
    auto dualIndex = daf::timeCount("buildDualIndex (ST_V6)", [&]() {
        return RCliqueSTv6::buildDualIndex(tree, cliqueIndex, r, s);
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

    // --- Bucket array (integer) ---
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

    // Per-leaf metadata for immutable-tree Case C
    const daf::Size initNumLeaves = tree.adj_list.size();
    // leafPivotC[L] = current remaining pivots in leaf L (decremented on Case C)
    std::vector<int> leafPivotC(initNumLeaves);
    std::vector<int> leafKeepC(initNumLeaves);
    for (daf::Size L = 0; L < initNumLeaves; ++L) {
        int p = 0, k = 0;
        for (const auto &node : tree.adj_list[L]) {
            if (node.isPivot) p++;
            else k++;
        }
        leafPivotC[L] = p;
        leafKeepC[L] = k;
    }

    daf::log_memory("Other index (include bucket)");

    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);

    std::cout << "=========================begin (r≥3 ST_V6)=========================" << std::endl;
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

            int pivotC = leafPivotC[leafId];
            int keepC = leafKeepC[leafId];
            int needPivot = s - keepC;
            int n = (int)leaf.size();

            auto t_bk = std::chrono::high_resolution_clock::now();

            // Compute per-vertex conflict degree
            daf::Size maxRCliquePerVertex = (daf::Size)((long long)(nCr[n - 1][r - 1] + 0.5));

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

            // Classify: count forced pivot removals, check for keep removal
            bool hasKeepConflict = false;
            int forcedPivotRemove = 0;
            for (int i = 0; i < n; ++i) {
                if (vertexConflictDeg[i] >= maxRCliquePerVertex) {
                    if (!leaf[i].isPivot) {
                        hasKeepConflict = true;
                        break;
                    }
                    forcedPivotRemove++;
                }
            }

            // === Case A: Leaf dies ===
            bool leafDead = hasKeepConflict;
            if (!leafDead) {
                int remainingPivots = pivotC - forcedPivotRemove;
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
                for (auto i : leaf)
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                tree.adj_list[leafId].clear();
                if (leafId < leafCliqueInfo.size())
                    leafCliqueInfo[leafId].clear();
                duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_struct).count();
                continue;
            }

            // === Case C check: only pivots with full conflict degree removed ===
            // Condition: forcedPivotRemove > 0, AND all conflict vertices are pivots
            // with conflictDeg == maxRCliquePerVertex (no partial pivot conflicts),
            // AND no non-forced vertex has any conflict (they can stay peacefully).
            bool isCaseC = false;
            if (forcedPivotRemove > 0) {
                // Check that NO vertex has a conflict degree between 1 and maxRCliquePerVertex-1
                // (i.e., every conflicted vertex is fully conflicted)
                bool allFullOrZero = true;
                for (int i = 0; i < n; ++i) {
                    if (vertexConflictDeg[i] > 0 && vertexConflictDeg[i] < maxRCliquePerVertex) {
                        allFullOrZero = false;
                        break;
                    }
                }
                isCaseC = allFullOrZero;
            }

            if (isCaseC) {
                // === Case C: Only pivots fully removed — delta formula, no BK ===
                cntCaseC++;
                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                auto t_supp = std::chrono::high_resolution_clock::now();

                int d = forcedPivotRemove;  // number of pivots being removed
                int oldPivotC = pivotC;
                int newPivotC = pivotC - d;

                // For each r-clique in this leaf, compute delta based on subPivotCount.
                //
                // An r-clique with subPivotCount = p has:
                //   old ncrValue = nCr[oldPivotC - p][needPivot - p]
                //   After removing d pivots:
                //     - If any of the d removed pivots is in this r-clique → r-clique LOST
                //       delta = -ncrValue (subtract it entirely)
                //     - If none of the d removed pivots is in this r-clique → r-clique SURVIVES
                //       new ncrValue = nCr[newPivotC - p][needPivot - p]
                //       delta = newNcrValue - oldNcrValue
                //
                // How do we know if a removed pivot is in the r-clique?
                // We DON'T need to check per-clique! We use combinatorics:
                //
                // For pivot class p: out of C(oldPivotC, p) * C(keepC, r-p) total r-cliques,
                //   - SURVIVING: C(newPivotC, p) * C(keepC, r-p) (choose p pivots from survivors)
                //   - LOST: the rest
                //
                // But we don't have classes — we have individual r-clique entries.
                // So we need to check per-entry whether the r-clique contains any removed pivot.
                //
                // FAST CHECK: use the cliqueIndex.byId() to get vertex IDs, then check
                // if any vertex is a removed pivot. Since removed pivots have
                // vertexConflictDeg == maxRCliquePerVertex, we can mark them.

                // Mark removed pivot positions
                // (vertexConflictDeg[i] >= maxRCliquePerVertex && leaf[i].isPivot already confirmed)

                if (leafId < leafCliqueInfo.size()) {
                    // Precompute new nCr values per pivot class
                    // newNcr[p] = nCr[newPivotC - p][needPivot - p] for surviving r-cliques
                    long long newNcrByClass[9] = {}; // r <= 8
                    bool validClass[9] = {};
                    for (int p = 0; p <= (int)r; ++p) {
                        if (p <= newPivotC && p <= needPivot &&
                            newPivotC - p < 1001 && needPivot - p >= 0 && needPivot - p < 401) {
                            newNcrByClass[p] = (long long)(nCr[newPivotC - p][needPivot - p] + 0.5);
                            validClass[p] = true;
                        }
                    }

                    // Build new leafCliqueInfo (only surviving entries, with updated ncrValue)
                    std::vector<RCliqueSTv6::LeafCliqueEntry> newEntries;
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
                            // Don't add to newEntries — r-clique is gone from this leaf
                        } else {
                            // SURVIVES: compute delta
                            int p = entry.subPivotCount;
                            if (validClass[p]) {
                                long long newNcr = newNcrByClass[p];
                                long long delta = newNcr - entry.ncrValue;
                                if (delta != 0 && rCliqueInHeap[entry.cliqueId]) {
                                    countingRClique[entry.cliqueId] += delta;
                                    if (countingRClique[entry.cliqueId] < 0)
                                        countingRClique[entry.cliqueId] = 0;
                                    bucketMove(entry.cliqueId);
                                }
                                // Update entry with new ncrValue
                                newEntries.push_back({entry.cliqueId, newNcr, entry.subPivotCount});
                            } else {
                                // nCr out of range — r-clique contribution becomes 0
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

                // Update leaf counters (immutable tree — don't modify tree/treeGraphV)
                leafPivotC[leafId] = newPivotC;

                // Tree mutation: remove only the dead pivot vertices
                for (int i = 0; i < n; ++i) {
                    if (vertexConflictDeg[i] >= maxRCliquePerVertex) {
                        treeGraphV.removeNbr(leaf[i].v, static_cast<TreeGraphNode>(leafId));
                    }
                }
                // Remove pivot vertices from the actual tree leaf
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

                    // Build leafCliqueInfo for new sub-leaf from parent's entries
                    std::vector<RCliqueSTv6::LeafCliqueEntry> newEntries;
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
                                if (entry.cliqueId < cliqueLeafIds.size()) {
                                    cliqueLeafIds[entry.cliqueId].push_back(newId);
                                }
                            }
                        }
                    }

                    if (newId >= leafCliqueInfo.size())
                        leafCliqueInfo.resize(newId + 1);
                    leafCliqueInfo[newId] = std::move(newEntries);

                    // Extend per-leaf metadata
                    if (newId >= leafPivotC.size()) {
                        leafPivotC.resize(newId + 1);
                        leafKeepC.resize(newId + 1);
                    }
                    leafPivotC[newId] = newPivotC;
                    leafKeepC[newId] = newKeepC;

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
