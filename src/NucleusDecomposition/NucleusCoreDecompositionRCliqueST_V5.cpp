//
// ST V5: Δ-Support with Positional Containment
// Key innovations over V4:
//   1. LeafCliqueEntry stores vertex positions within leaf (uint16_t[8])
//   2. BK path uses bitset positional containment (O(r) bit tests) instead of hash map lookup
//   3. Δ-Support: compute net support change in one pass over leafCliqueInfo
//   4. Union pre-filter: skip LOST r-cliques with O(r) bit tests
//   5. Merged delta computation + leafCliqueInfo construction in single traversal
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

namespace RCliqueSTv5 {

struct LeafCliqueEntry {
    daf::Size cliqueId;      // 4 bytes
    double ncrValue;       // 8 bytes
    uint16_t positions[8];    // 16 bytes: vertex positions within leaf (max r=8)
    // No pivotMask needed — sub-leaf pivot assignment differs from parent
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

        daf::CliqueSize pivotC = 0, keepC = 0;
        for (const auto &i : leaf) {
            if (i.isPivot) pivotC++;
            else keepC++;
        }
        int needPivot = s - static_cast<int>(keepC);
        int n = (int)leaf.size();

        // Build vertex→position map for this leaf
        daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
        for (int i = 0; i < n; ++i)
            mapRef[leaf[i].v] = (daf::Size)i;

        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            daf::CliqueSize subNumPivot = 0;
            for (const auto &node : rClique) if (node.isPivot) subNumPivot++;
            if (subNumPivot > (daf::CliqueSize)needPivot) return true;
            if (pivotC - subNumPivot >= 1001 || (daf::Size)needPivot - subNumPivot >= 401) return true;
            double ncrValue = nCr[pivotC - subNumPivot][needPivot - subNumPivot];
            auto id = cliqueIndex.byClique(rClique);
            if (id < nClique) {
                result.counting[id] += ncrValue;

                LeafCliqueEntry entry;
                entry.cliqueId = id;
                entry.ncrValue = ncrValue;
                for (int i = 0; i < (int)r; ++i)
                    entry.positions[i] = (uint16_t)mapRef[rClique[i].v];
                // Zero remaining positions
                for (int i = (int)r; i < 8; ++i)
                    entry.positions[i] = 0;

                result.leafCliqueInfo[leafId].push_back(entry);
                result.cliqueLeafIds[id].push_back(leafId);
            }
            return true;
        });
    }

    return result;
}

// Sub-leaf info collected from BK callback
struct SubLeafInfo {
    DynBitset cover;                  // sub-leaf vertex bitset (positions 0..n-1 in parent)
    DynBitset pivots;                 // which positions are pivots in sub-leaf
    daf::Size newId;                  // tree node ID for this sub-leaf
    daf::CliqueSize newPivotC;
    daf::CliqueSize newKeepC;
};

} // namespace RCliqueSTv5


std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V5(
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

    // Build both indices and counting in one pass
    auto dualIndex = daf::timeCount("buildDualIndex (ST_V5)", [&]() {
        return RCliqueSTv5::buildDualIndex(tree, cliqueIndex, r, s);
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

    daf::log_memory("Other index (include bucket)");

    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);

    // Reusable buffer for sub-leaf info
    std::vector<RCliqueSTv5::SubLeafInfo> subLeaves;
    subLeaves.reserve(64);

    // Reusable oldPos→newPos mapping buffer
    uint16_t oldToNew[400];

    std::cout << "=========================begin (r≥3 ST_V5)=========================" << std::endl;
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

        // --- Phase 2: BK + Δ-Support per leaf ---
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

                // Subtract old leaf contribution for all r-cliques
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

                // Tree mutation
                auto t_struct = std::chrono::high_resolution_clock::now();
                for (const auto& i : leaf)
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                tree.adj_list[leafId].clear();
                if (leafId < leafCliqueInfo.size())
                    leafCliqueInfo[leafId].clear();
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

                auto mapped = removedR | std::views::transform(
                    [&](const daf::Size id) { return cliqueIndex.byId(id); });

                // Collect all sub-leaves from BK
                subLeaves.clear();

                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivotsBs) {
                        auto newLeaf = bkRmClique::coverToVertex(c, pivotsBs, leaf);
                        auto newId = tree.addNode(newLeaf);
                        const auto &stored = tree.adj_list[newId];
                        daf::CliqueSize newPivotC = 0, newKeepC = 0;
                        for (const auto &i : stored) {
                            if (i.isPivot) { treeGraphV.addNbr(i.v, {newId, true}); newPivotC++; }
                            else { treeGraphV.addNbr(i.v, {newId, false}); newKeepC++; }
                        }

                        RCliqueSTv5::SubLeafInfo sl;
                        sl.cover = c;
                        sl.pivots = pivotsBs;
                        sl.newId = newId;
                        sl.newPivotC = newPivotC;
                        sl.newKeepC = newKeepC;
                        subLeaves.push_back(std::move(sl));

                        // Ensure changedLeafIndex capacity
                        if (newId >= changedLeafIndex.size())
                            changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                    });

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                // === Δ-Support with Positional Containment ===
                auto t_supp = std::chrono::high_resolution_clock::now();

                if (leafId < leafCliqueInfo.size() && !subLeaves.empty()) {
                    const int numSL = (int)subLeaves.size();

                    // 1. Compute survivor union (bitwise OR of all sub-leaf covers)
                    DynBitset survivorUnion;
                    survivorUnion.setSize(n);
                    survivorUnion.reset();
                    for (const auto &sl : subLeaves)
                        survivorUnion |= sl.cover;

                    // 2. Pre-compute oldPos→newPos mapping for each sub-leaf
                    // Store as vector of arrays
                    std::vector<std::vector<uint16_t>> oldToNewMaps(numSL);
                    for (int j = 0; j < numSL; ++j) {
                        oldToNewMaps[j].resize(n, 0);
                        uint16_t newPos = 0;
                        subLeaves[j].cover.for_each_bit([&](size_t bit) {
                            oldToNewMaps[j][bit] = newPos++;
                        });
                    }

                    // 3. Prepare new leafCliqueInfo vectors for sub-leaves
                    std::vector<std::vector<RCliqueSTv5::LeafCliqueEntry>> newEntriesVec(numSL);

                    // Ensure leafCliqueInfo has capacity for new sub-leaf IDs
                    daf::Size maxNewId = 0;
                    for (const auto &sl : subLeaves)
                        maxNewId = std::max(maxNewId, sl.newId);
                    if (maxNewId >= leafCliqueInfo.size())
                        leafCliqueInfo.resize(maxNewId + 1);

                    // 4. Single pass over parent's leafCliqueInfo:
                    //    - compute delta support
                    //    - build new leafCliqueInfo entries
                    //    - update cliqueLeafIds
                    for (const auto &entry : leafCliqueInfo[leafId]) {
                        const bool inHeap = rCliqueInHeap[entry.cliqueId];

                        // Union pre-filter: check if r-clique survives in ANY sub-leaf
                        bool inUnion = true;
                        for (int i = 0; i < (int)r; ++i) {
                            if (!survivorUnion.test(entry.positions[i])) {
                                inUnion = false;
                                break;
                            }
                        }

                        if (!inUnion) {
                            // LOST r-clique: subtract old contribution
                            if (inHeap) {
                                countingRClique[entry.cliqueId] -= entry.ncrValue;
                                if (countingRClique[entry.cliqueId] < 0)
                                    countingRClique[entry.cliqueId] = 0;
                                bucketMove(entry.cliqueId);
                            }
                            continue;
                        }

                        // Surviving r-clique: check each sub-leaf
                        long long totalNew = 0;
                        for (int j = 0; j < numSL; ++j) {
                            const auto &sl = subLeaves[j];

                            // Positional containment check: O(r) bit tests
                            bool contained = true;
                            daf::CliqueSize subP = 0;
                            for (int i = 0; i < (int)r; ++i) {
                                uint16_t pos = entry.positions[i];
                                if (!sl.cover.test(pos)) {
                                    contained = false;
                                    break;
                                }
                                if (sl.pivots.test(pos)) subP++;
                            }
                            if (!contained) continue;

                            daf::Size np = s - sl.newKeepC;
                            if (subP > np || sl.newPivotC - subP >= 1001 || np - subP >= 401)
                                continue;

                            double ncrVal = nCr[sl.newPivotC - subP][np - subP];
                            totalNew += ncrVal;

                            // Build leafCliqueInfo entry for this sub-leaf
                            RCliqueSTv5::LeafCliqueEntry ne;
                            ne.cliqueId = entry.cliqueId;
                            ne.ncrValue = ncrVal;
                            for (int i = 0; i < (int)r; ++i)
                                ne.positions[i] = oldToNewMaps[j][entry.positions[i]];
                            for (int i = (int)r; i < 8; ++i)
                                ne.positions[i] = 0;
                            newEntriesVec[j].push_back(ne);

                            // Update cliqueLeafIds
                            if (entry.cliqueId < cliqueLeafIds.size())
                                cliqueLeafIds[entry.cliqueId].push_back(sl.newId);
                        }

                        // Apply delta: new contributions - old contribution
                        if (inHeap) {
                            countingRClique[entry.cliqueId] += (totalNew - entry.ncrValue);
                            if (countingRClique[entry.cliqueId] < 0)
                                countingRClique[entry.cliqueId] = 0;
                            bucketMove(entry.cliqueId);
                        }
                    }

                    // 5. Install new leafCliqueInfo for each sub-leaf
                    for (int j = 0; j < numSL; ++j) {
                        leafCliqueInfo[subLeaves[j].newId] = std::move(newEntriesVec[j]);
                    }
                } else if (leafId < leafCliqueInfo.size()) {
                    // No sub-leaves produced — subtract all old contributions
                    for (const auto &entry : leafCliqueInfo[leafId]) {
                        if (!rCliqueInHeap[entry.cliqueId]) continue;
                        countingRClique[entry.cliqueId] -= entry.ncrValue;
                        if (countingRClique[entry.cliqueId] < 0)
                            countingRClique[entry.cliqueId] = 0;
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
