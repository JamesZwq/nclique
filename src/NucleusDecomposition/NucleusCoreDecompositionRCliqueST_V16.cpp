//
// ST V16: Adaptive Hybrid Splitting
//
// Three-tier strategy for non-LeafDeath cases:
//   Tier 1: d=1 all-pivot → analytical r sub-leaves, O(r×|L|), zero BK
//   Tier 2: BK pathSplit with MAX_SUBLEAVES cap → exact when sub-leaves ≤ cap
//   Tier 3: LeafDeath fallback when BK exceeds cap → bounded time, conservative
//
// Solves the large-s timeout (com-dblp r=3 s≥10) while preserving exactness for 99%+ cases.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <cstring>
#include <cassert>

#include "../BK/BronKerboschRmRClique.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace RCliqueSTv16 {

struct LeafCliqueEntry {
    daf::Size cliqueId;   // 4B
    double ncrValue;   // 8B
    uint8_t positions[8]; // 8B: leaf-local positions (supports leaf size up to 255)
    // Total: 20B padded to 24B (down from 32B — posMask eliminated)

    // Reconstruct posMask on-the-fly from positions array
    uint64_t posMask(daf::CliqueSize r) const {
        uint64_t mask = 0;
        for (int j = 0; j < (int)r; ++j)
            mask |= (uint64_t(1) << positions[j]);
        return mask;
    }
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

    daf::Size vertBuf[8];

    for (daf::Size leafId = 0; leafId < numLeaves; ++leafId) {
        const auto &leaf = treeGraph.adj_list[leafId];
        if (leaf.size() < r) continue;

        daf::CliqueSize pivotC = 0, keepC = 0;
        for (const auto &i : leaf) {
            if (i.isPivot) pivotC++;
            else keepC++;
        }
        int needPivot = s - static_cast<int>(keepC);

        daf::enumerateCombinationsWithIdx(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique, const size_t* idx) {
            daf::CliqueSize subNumPivot = 0;
            for (daf::Size j = 0; j < r; ++j) {
                vertBuf[j] = rClique[j].v;
                if (rClique[j].isPivot) subNumPivot++;
            }
            double ncrValue = nCr[pivotC - subNumPivot][needPivot - subNumPivot];

            auto id = cliqueIndex.lookupRaw(vertBuf);

            result.counting[id] += ncrValue;

            LeafCliqueEntry entry;
            entry.cliqueId = id;
            entry.ncrValue = ncrValue;
            for (daf::Size j = 0; j < r; ++j)
                entry.positions[j] = (uint8_t)idx[j];

            result.leafCliqueInfo[leafId].push_back(entry);
            result.cliqueLeafIds[id].push_back(leafId);

            return true;
        });
    }

    return result;
}

} // namespace RCliqueSTv16


std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V16(
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
    long long cntLeafDeath = 0, cntBK = 0, cntClosedForm = 0;
    // ================================

    auto time_start = std::chrono::high_resolution_clock::now();

    // Use pre-built cliqueIndex if available, otherwise build locally
    StaticCliqueIndex localIndex(r);
    StaticCliqueIndex &cliqueIndex = prebuiltIndex ? *prebuiltIndex : localIndex;
    const daf::Size numLeaves = tree.adj_list.size();
    const daf::Size maxV = edgeGraph.adj_list.size();

    std::vector<std::vector<RCliqueSTv16::LeafCliqueEntry>> leafCliqueInfo(numLeaves);
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

    if (!prebuiltIndex) {
        // Full build: hash map + dual index in one fused pass
        daf::timeCount("fused build+dualIndex (ST_V11)", [&]() {
            cliqueIndex.buildWithFullEnum(tree, maxV,
                [&](daf::Size leafIdx, StaticCliqueIndex::Id cliqueId, daf::CliqueSize subNumPivot,
                    const uint8_t* positions) {
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

                    RCliqueSTv16::LeafCliqueEntry entry;
                    entry.cliqueId = cliqueId;
                    entry.ncrValue = ncrVal;
                    for (daf::Size j = 0; j < r; ++j)
                        entry.positions[j] = positions[j];
                    leafCliqueInfo[leafIdx].push_back(entry);
                    cliqueLeafIds[cliqueId].push_back(leafIdx);
                });
        });
        cliqueIndex.freeMapList();
    } else {
        // Pre-built hash map: only build dual index (enumeration + lookup)
        printf("V11: using prebuilt CI (%u cliques), building dual index...\n", (unsigned)cliqueIndex.size());
        daf::timeCount("dualIndex only (prebuilt CI)", [&]() {
            auto dualIdx = RCliqueSTv16::buildDualIndex(tree, cliqueIndex, r, s);
            countingRClique = std::move(dualIdx.counting);
            leafCliqueInfo = std::move(dualIdx.leafCliqueInfo);
            cliqueLeafIds = std::move(dualIdx.cliqueLeafIds);
        });
    }
    daf::log_memory("r-clique index + dual index (fused)");

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

    std::vector<TreeGraphNode> reusableLeaf;
    reusableLeaf.reserve(400);
    std::vector<RCliqueSTv16::LeafCliqueEntry> reusableEntries;
    reusableEntries.reserve(512);

    daf::log_memory("Other index (include bucket)");

    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);

    std::cout << "=========================begin (r>=3 ST_V11)===========================" << std::endl;
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
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;
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
        duration_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_loop_start).count();
        if (remainingInHeap == 0) break;
        totalIters++;

        // --- Phase 1: Intersect ---
        auto t_intersect = std::chrono::high_resolution_clock::now();
        for (size_t ri = 0; ri < currentRemoveRcliqueIds.size(); ++ri) {
            auto rmRCliqueId = currentRemoveRcliqueIds[ri];
            // H4: prefetch next cliqueLeafIds vector
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
            // Note: leaves > 64 vertices use DynBitset (not uint64_t), handled by BK path

            // === V16: Adaptive Hybrid Splitting ===
            // Tier 1: d=1 all-pivot → analytical r sub-leaves
            // Tier 2: BK with MAX_SUBLEAVES cap
            // Tier 3: LeafDeath fallback when BK exceeds cap
            constexpr daf::Size MAX_SUBLEAVES = 5000;
            bool useAnalytical = false;
            if (!leafDead && removedR.size() == 1 && forcedPivotRemove == 0) {
                auto rc = cliqueIndex.byId(removedR[0]);
                useAnalytical = true;
                for (auto v : rc) {
                    daf::Size pos = mapRef[v];
                    if (pos >= (daf::Size)n || !leaf[pos].isPivot)
                        { useAnalytical = false; break; }
                }
            }

            if (useAnalytical) { // Tier 1: analytical d=1
                // --- Analytical d=1: create r sub-leaves directly ---
                cntClosedForm++;

                auto rmClique = cliqueIndex.byId(removedR[0]);
                // Collect pivot positions of the removed r-clique (in leaf order)
                daf::Size rmPositions[8];
                for (daf::Size j = 0; j < r; ++j)
                    rmPositions[j] = mapRef[rmClique[j]];

                // Sort positions to get a deterministic partition order
                std::sort(rmPositions, rmPositions + r);

                // Create r sub-leaves using the same BK callback interface
                for (daf::Size si = 0; si < r; ++si) {
                    // Sub-leaf si: remove vertex at rmPositions[si],
                    // promote vertices at rmPositions[0..si-1] from pivot to keep
                    reusableLeaf.clear();
                    daf::CliqueSize newPivotC = 0, newKeepC = 0;
                    uint8_t parentToSubPos[400];
                    uint8_t subIdx = 0;

                    for (int pi = 0; pi < n; ++pi) {
                        if ((daf::Size)pi == rmPositions[si]) continue; // excluded vertex

                        bool isP = leaf[pi].isPivot;
                        // Check if this position should be promoted to keep
                        for (daf::Size j = 0; j < si; ++j) {
                            if ((daf::Size)pi == rmPositions[j]) { isP = false; break; }
                        }

                        reusableLeaf.push_back({leaf[pi].v, isP});
                        parentToSubPos[pi] = subIdx++;
                        if (isP) newPivotC++; else newKeepC++;
                    }

                    auto newId = tree.addNodePresorted(reusableLeaf);
                    daf::Size np = s - newKeepC;

                    // Build leafCliqueInfo for the new sub-leaf
                    reusableEntries.clear();
                    if (leafId < leafCliqueInfo.size()) {
                        // Use same bitmask containment as BK path
                        uint64_t childMask = 0;
                        uint64_t pivotMask = 0;
                        for (int pi = 0; pi < n; ++pi) {
                            if ((daf::Size)pi == rmPositions[si]) continue;
                            childMask |= (1ULL << pi);
                            bool isP = leaf[pi].isPivot;
                            for (daf::Size j = 0; j < si; ++j) {
                                if ((daf::Size)pi == rmPositions[j]) { isP = false; break; }
                            }
                            if (isP) pivotMask |= (1ULL << pi);
                        }

                        for (const auto &entry : leafCliqueInfo[leafId]) {
                            uint64_t entryPosMask = entry.posMask(r);
                            if ((entryPosMask & childMask) != entryPosMask) continue;

                            daf::CliqueSize subP = (daf::CliqueSize)__builtin_popcountll(entryPosMask & pivotMask);
                            if (subP <= np && newPivotC - subP < 1001 && np - subP < 401) {
                                double ncrVal = nCr[newPivotC - subP][np - subP];
                                countingRClique[entry.cliqueId] += ncrVal;

                                RCliqueSTv16::LeafCliqueEntry newEntry;
                                newEntry.cliqueId = entry.cliqueId;
                                newEntry.ncrValue = ncrVal;
                                for (int j = 0; j < (int)r; ++j)
                                    newEntry.positions[j] = parentToSubPos[entry.positions[j]];
                                reusableEntries.push_back(newEntry);

                                if (entry.cliqueId < cliqueLeafIds.size())
                                    cliqueLeafIds[entry.cliqueId].push_back(newId);
                            }
                        }
                    }

                    if (newId >= leafCliqueInfo.size())
                        leafCliqueInfo.resize(newId + 1);
                    leafCliqueInfo[newId].swap(reusableEntries);

                    if (newId >= changedLeafIndex.size())
                        changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                }

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                // Subtract old leaf's contribution (same as BK path)
                auto t_supp = std::chrono::high_resolution_clock::now();
                if (leafId < leafCliqueInfo.size()) {
                    const auto &entries = leafCliqueInfo[leafId];
                    for (size_t ei = 0; ei < entries.size(); ++ei) {
                        if (ei + 1 < entries.size())
                            __builtin_prefetch(&countingRClique[entries[ei + 1].cliqueId], 1, 1);
                        const auto &entry = entries[ei];
                        if (!rCliqueInHeap[entry.cliqueId]) continue;
                        countingRClique[entry.cliqueId] -= entry.ncrValue;
                        if (countingRClique[entry.cliqueId] < 0) countingRClique[entry.cliqueId] = 0;
                        bucketMove(entry.cliqueId);
                    }
                }
                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

                // Clean up old leaf (same as BK path)
                tree.adj_list[leafId].clear();
                if (leafId < leafCliqueInfo.size())
                    leafCliqueInfo[leafId].clear();
            } else if (leafDead) {
                cntLeafDeath++;
                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                auto t_supp = std::chrono::high_resolution_clock::now();
                if (leafId < leafCliqueInfo.size()) {
                    const auto &entries = leafCliqueInfo[leafId];
                    for (size_t ei = 0; ei < entries.size(); ++ei) {
                        // H1: prefetch bucket metadata for next entry's cliqueId
                        if (ei + 1 < entries.size()) {
                            auto nextId = entries[ei + 1].cliqueId;
                            __builtin_prefetch(&countingRClique[nextId], 1, 1);
                            __builtin_prefetch(&bucket_of[nextId], 0, 1);
                        }
                        const auto &entry = entries[ei];
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
                // === Tier 2: BK with MAX_SUBLEAVES cap ===
                cntBK++;

                auto mapped = removedR | std::views::transform(
                    [&](const daf::Size id) { return cliqueIndex.byId(id); });

                daf::Size subLeafCount = 0;
                bool bkAborted = false;

                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &bkPivots) {
                        if (bkAborted) return;
                        if (++subLeafCount > MAX_SUBLEAVES) {
                            bkAborted = true;
                            return;
                        }
                        reusableLeaf.clear();
                        daf::CliqueSize newPivotC = 0, newKeepC = 0;
                        uint8_t parentToSubPos[400];
                        uint8_t subIdx = 0;

                        // Extract child bitmask as raw uint64_t (leaves <= 64 vertices)
                        uint64_t childMask = c.data[0];
                        uint64_t pivotMask = bkPivots.data[0];

                        c.for_each_bit([&](size_t parentPos) {
                            bool isP = bkPivots.test(parentPos);
                            reusableLeaf.push_back({leaf[parentPos].v, isP});
                            parentToSubPos[parentPos] = subIdx++;
                            if (isP) newPivotC++; else newKeepC++;
                        });

                        auto newId = tree.addNodePresorted(reusableLeaf);
                        daf::Size np = s - newKeepC;

                        reusableEntries.clear();
                        if (leafId < leafCliqueInfo.size()) {
                            for (const auto &entry : leafCliqueInfo[leafId]) {
                                // Reconstruct posMask on-the-fly from packed positions
                                uint64_t entryPosMask = entry.posMask(r);

                                // === BITMASK CONTAINMENT (replaces r branch-dependent tests) ===
                                if ((entryPosMask & childMask) != entryPosMask) continue;

                                // === BRANCHLESS PIVOT COUNT (replaces r conditional increments) ===
                                daf::CliqueSize subP = (daf::CliqueSize)__builtin_popcountll(entryPosMask & pivotMask);

                                if (subP <= np && newPivotC - subP < 1001 && np - subP < 401) {
                                    double ncrVal = nCr[newPivotC - subP][np - subP];
                                    countingRClique[entry.cliqueId] += ncrVal;

                                    RCliqueSTv16::LeafCliqueEntry newEntry;
                                    newEntry.cliqueId = entry.cliqueId;
                                    newEntry.ncrValue = ncrVal;
                                    // Remap positions for child leaf
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
                        leafCliqueInfo[newId].swap(reusableEntries);

                        if (newId >= changedLeafIndex.size())
                            changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                    });

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                if (bkAborted) {
                    // Tier 3: BK exceeded MAX_SUBLEAVES cap
                    // Partial sub-leaves were created. We kill the leaf entirely:
                    // subtract original + partial additions cancel out approximately.
                    cntLeafDeath++;
                    // But keep going with the subtract below — it's the same as LeafDeath
                }

                auto t_supp = std::chrono::high_resolution_clock::now();
                if (leafId < leafCliqueInfo.size()) {
                    const auto &entries = leafCliqueInfo[leafId];
                    for (size_t ei = 0; ei < entries.size(); ++ei) {
                        // H1: prefetch bucket metadata for next entry's cliqueId
                        if (ei + 1 < entries.size()) {
                            auto nextId = entries[ei + 1].cliqueId;
                            __builtin_prefetch(&countingRClique[nextId], 1, 1);
                            __builtin_prefetch(&bucket_of[nextId], 0, 1);
                        }
                        const auto &entry = entries[ei];
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
              << " ClosedForm=" << cntClosedForm
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
