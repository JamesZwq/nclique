//
// ST V18: Adaptive Hybrid — V17's light init + lazy leafCliqueInfo for BK leaves
//
// Key idea: 94-97% of leaves are dead and never need leafCliqueInfo.
// Build leafCliqueInfo on-demand ONLY when a leaf enters BK.
// Dead leaves use direct enumeration for support subtraction (one-time, no caching).
// BK leaves build leafCliqueInfo once, reuse for containment + subtraction.
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

namespace RCliqueSTv18 {

struct LeafCliqueEntry {
    daf::Size cliqueId;
    double ncrValue;
    uint8_t positions[8]; // leaf-local positions
};

// Build leafCliqueInfo for ONE leaf on-demand
static std::vector<LeafCliqueEntry> buildLeafInfo(
    const std::vector<TreeGraphNode> &leaf,
    const StaticCliqueIndex &cliqueIndex,
    daf::CliqueSize r, daf::CliqueSize s) {

    std::vector<LeafCliqueEntry> entries;
    daf::CliqueSize pivotC = 0, keepC = 0;
    for (const auto &i : leaf) {
        if (i.isPivot) pivotC++; else keepC++;
    }
    int needPivot = s - static_cast<int>(keepC);

    daf::enumerateCombinationsWithIdx(leaf, r,
        [&](const daf::StaticVector<TreeGraphNode> &rClique, const size_t* idx) {
        daf::CliqueSize subNumPivot = 0;
        daf::Size vertBuf[8];
        for (daf::Size j = 0; j < r; ++j) {
            vertBuf[j] = rClique[j].v;
            if (rClique[j].isPivot) subNumPivot++;
        }
        if (subNumPivot <= needPivot) {
            int row = (int)pivotC - (int)subNumPivot;
            int col = needPivot - (int)subNumPivot;
            if (row >= 0 && row < 1001 && col >= 0 && col < 401) {
                double ncrVal = nCr[row][col];
                auto id = cliqueIndex.lookupRaw(vertBuf);
                LeafCliqueEntry entry;
                entry.cliqueId = id;
                entry.ncrValue = ncrVal;
                for (daf::Size j = 0; j < r; ++j)
                    entry.positions[j] = (uint8_t)idx[j];
                entries.push_back(entry);
            }
        }
        return true;
    });
    return entries;
}

} // namespace RCliqueSTv18

std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V18(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    long long duration_init = 0, duration_pop = 0, duration_intersect = 0;
    long long duration_bk = 0, duration_support = 0;
    long long cntLeafDeath = 0, cntBK = 0;

    auto time_start = std::chrono::high_resolution_clock::now();

    // ========== INIT: V17's lightweight approach (no dual index) ==========
    StaticCliqueIndex localIndex(r);
    StaticCliqueIndex &cliqueIndex = prebuiltIndex ? *prebuiltIndex : localIndex;
    if (!prebuiltIndex) {
        daf::timeCount("clique Index build", [&]() {
            localIndex.build(tree, edgeGraph.adj_list.size());
        });
    }

    // Count support (Reference's approach: iterate tree once)
    std::vector<double> countingRClique;
    daf::timeCount("countingPerRClique (V18)", [&]() {
        countingRClique.assign(cliqueIndex.size(), 0.0);
        for (const auto &leaf : tree.adj_list) {
            if (leaf.size() < r) continue;
            daf::CliqueSize pivotC = 0, keepC = 0;
            for (const auto &i : leaf) { if (i.isPivot) pivotC++; else keepC++; }
            int needPivot = s - static_cast<int>(keepC);
            daf::enumerateCombinations(leaf, r,
                [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                daf::CliqueSize subNumPivot = 0;
                for (const auto &node : rClique)
                    if (node.isPivot) subNumPivot++;
                if (subNumPivot <= needPivot) {
                    int row = (int)pivotC - (int)subNumPivot;
                    int col = needPivot - (int)subNumPivot;
                    if (row >= 0 && row < 1001 && col >= 0 && col < 401) {
                        countingRClique[cliqueIndex.byClique(rClique)] += nCr[row][col];
                    }
                }
                return true;
            });
        }
    });

    const daf::Size nClique = cliqueIndex.size();
    daf::log_memory("r-clique index + counting");

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

    // Lazy leafCliqueInfo: only populated for BK leaves and their sub-leaves
    std::vector<std::vector<RCliqueSTv18::LeafCliqueEntry>> leafCliqueInfo(tree.adj_list.size());

    // ========== V11's hybrid bucket+set PQ ==========
    constexpr double BUCKET_THRESHOLD = 5000000.0;
    double rawMaxBucket = 0;
    for (daf::Size i = 0; i < nClique; ++i)
        rawMaxBucket = std::max(rawMaxBucket, countingRClique[i]);
    int maxBucket = (int)std::min(rawMaxBucket, BUCKET_THRESHOLD);

    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::set<std::pair<double, daf::Size>> overflowSet;
    std::vector<int> bucket_of(nClique, -1);
    std::vector<daf::Size> pos_in_bucket(nClique);
    std::vector<double> overflowStoredVal(nClique, -1);

    for (daf::Size i = 0; i < nClique; ++i) {
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
    daf::Size remainingInHeap = nClique;

    if (overflowSet.size() > 0)
        printf("Hybrid bucket+set: %u in buckets, %zu in overflow (maxSupport=%.0f)\n",
               (unsigned)(nClique - overflowSet.size()), overflowSet.size(), rawMaxBucket);

    auto bucketMove = [&](daf::Size id) {
        if (!rCliqueInHeap[id]) return;
        double val = std::max(0.0, countingRClique[id]);
        int oldB = bucket_of[id];
        if (oldB == -1) overflowSet.erase({overflowStoredVal[id], id});
        if (val <= BUCKET_THRESHOLD) {
            int newB = (int)val;
            if (oldB >= 0 && newB == oldB) return;
            if (oldB >= 0) {
                auto &v = buckets[oldB]; daf::Size p = pos_in_bucket[id];
                if (p < v.size() - 1) { auto last = v.back(); v[p] = last; pos_in_bucket[last] = p; }
                v.pop_back();
            }
            bucket_of[id] = newB;
            pos_in_bucket[id] = buckets[newB].size();
            buckets[newB].push_back(id);
            if (newB < curBucket) curBucket = newB;
        } else {
            if (oldB >= 0) {
                auto &v = buckets[oldB]; daf::Size p = pos_in_bucket[id];
                if (p < v.size()) {
                    if (p < v.size() - 1) { auto last = v.back(); v[p] = last; pos_in_bucket[last] = p; }
                    v.pop_back();
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
            if (countingRClique[id] <= BUCKET_THRESHOLD) {
                overflowSet.erase(overflowSet.begin());
                int b = (int)countingRClique[id];
                bucket_of[id] = b; pos_in_bucket[id] = buckets[b].size();
                buckets[b].push_back(id);
            } else break;
        }
    };

    daf::log_memory("Other index (include bucket)");
    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);
    std::vector<TreeGraphNode> reusableLeaf;
    reusableLeaf.reserve(400);
    std::vector<RCliqueSTv18::LeafCliqueEntry> reusableEntries;
    reusableEntries.reserve(512);

    std::cout << "=========================begin (r>=3 ST_V18)===========================" << std::endl;
    double minCore = 0;
    long long totalIters = 0;

    while (remainingInHeap > 0) {
        auto t_loop_start = std::chrono::high_resolution_clock::now();

        for (auto &leafId : changedLeaf)
            changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        changedLeaf.clear();
        removedRCliqueIdForLeaf.clear();
        currentRemoveRcliqueIds.clear();

        // --- Pop ---
        drainOverflowToBucket();
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) {
            if (!overflowSet.empty()) {
                while (!overflowSet.empty()) {
                    daf::Size id = overflowSet.begin()->second;
                    overflowSet.erase(overflowSet.begin());
                    if (!rCliqueInHeap[id]) continue;
                    minCore = std::max(countingRClique[id], minCore);
                    rCliqueInHeap[id] = false;
                    currentRemoveRcliqueIds.push_back(id);
                    coreRClique[id] = minCore; remainingInHeap--;
                    while (!overflowSet.empty()) {
                        daf::Size next = overflowSet.begin()->second;
                        if (!rCliqueInHeap[next]) { overflowSet.erase(overflowSet.begin()); continue; }
                        if (countingRClique[next] <= minCore) {
                            overflowSet.erase(overflowSet.begin());
                            rCliqueInHeap[next] = false;
                            currentRemoveRcliqueIds.push_back(next);
                            coreRClique[next] = minCore; remainingInHeap--;
                        } else break;
                    }
                    break;
                }
                if (currentRemoveRcliqueIds.empty()) break;
                goto phase1_v18;
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
                coreRClique[id] = minCore; remainingInHeap--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (curBucket + 1) <= (int)minCore) curBucket++;
            else break;
        }
        phase1_v18:
        duration_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_loop_start).count();
        if (remainingInHeap == 0) break;
        totalIters++;

        // --- Phase 1: Intersect via treeGraphV ---
        auto t_intersect = std::chrono::high_resolution_clock::now();
        for (auto rmRCliqueId : currentRemoveRcliqueIds) {
            auto rClique = cliqueIndex.byId(rmRCliqueId);
            daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                [&](const TreeGraphNode &uClique) {
                    auto &leafIdx = changedLeafIndex[uClique.v];
                    if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                        leafIdx = removedRCliqueIdForLeaf.size();
                        removedRCliqueIdForLeaf.emplace_back();
                        changedLeaf.push_back(uClique.v);
                        removedRCliqueIdForLeaf.back().reserve(10);
                    }
                    removedRCliqueIdForLeaf[leafIdx].emplace_back(rmRCliqueId);
                });
        }
        duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_intersect).count();

        // --- Phase 2: Process changed leaves ---
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            auto leafId = changedLeaf[idx];
            const auto &leaf = tree.adj_list[leafId];
            if (leaf.empty()) continue;
            auto leafIndex = changedLeafIndex[leafId];
            const auto &removedR = removedRCliqueIdForLeaf[leafIndex];

            daf::CliqueSize keepC = 0, pivotC = 0;
            for (const auto &node : leaf) {
                if (node.isPivot) pivotC++; else keepC++;
            }
            int n = (int)leaf.size();
            int needPivot = s - static_cast<int>(keepC);

            auto t_bk = std::chrono::high_resolution_clock::now();

            // --- LeafDeath fast-path (V11's optimization) ---
            daf::Size maxRCliquePerVertex = (daf::Size)(nCr[n - 1][r - 1] + 0.5);
            vertexConflictDeg.assign(n, 0);
            daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
            for (int i = 0; i < n; ++i) mapRef[leaf[i].v] = (daf::Size)i;

            for (auto rmId : removedR) {
                auto rClique = cliqueIndex.byId(rmId);
                for (auto v : rClique) {
                    daf::Size pos = mapRef[v];
                    if (pos < (daf::Size)n) vertexConflictDeg[pos]++;
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
                // ===== DEAD LEAF: V17's direct enumeration (no leafCliqueInfo) =====
                cntLeafDeath++;
                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                auto t_supp = std::chrono::high_resolution_clock::now();
                // Remove from treeGraphV
                for (const auto &leafV : leaf) {
                    if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                    else treeGraphV.removeNbr(leafV.v, {leafId, false});
                }
                // Subtract support via direct enumeration
                daf::enumerateCombinations(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &clique) {
                    auto cliqueId = cliqueIndex.byClique(clique);
                    if (!rCliqueInHeap[cliqueId]) return true;
                    daf::CliqueSize subP = 0;
                    for (const auto &node : clique) if (node.isPivot) subP++;
                    countingRClique[cliqueId] -= nCr[pivotC - subP][needPivot - subP];
                    if (countingRClique[cliqueId] < 0) countingRClique[cliqueId] = 0;
                    bucketMove(cliqueId);
                    return true;
                });
                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();
                tree.removeNode(leafId);

            } else {
                // ===== BK LEAF: build leafCliqueInfo on-demand, then use V11's fast path =====
                cntBK++;

                // Build leafCliqueInfo for this leaf (if not already built from being a sub-leaf)
                if (leafId >= leafCliqueInfo.size()) leafCliqueInfo.resize(leafId + 1);
                if (leafCliqueInfo[leafId].empty()) {
                    leafCliqueInfo[leafId] = RCliqueSTv18::buildLeafInfo(leaf, cliqueIndex, r, s);
                }

                // Remove from treeGraphV
                for (const auto &leafV : leaf) {
                    if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                    else treeGraphV.removeNbr(leafV.v, {leafId, false});
                }

                // BK split with leafCliqueInfo-based sub-leaf construction
                auto mapped = removedR | std::views::transform(
                    [&](const daf::Size id) { return cliqueIndex.byId(id); });

                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &bkPivots) {
                        reusableLeaf.clear();
                        daf::CliqueSize newPivotC = 0, newKeepC = 0;
                        uint8_t parentToSubPos[400];
                        memset(parentToSubPos, 255, sizeof(parentToSubPos));
                        uint8_t subIdx = 0;

                        c.for_each_bit([&](size_t parentPos) {
                            bool isP = bkPivots.test(parentPos);
                            reusableLeaf.push_back({leaf[parentPos].v, isP});
                            parentToSubPos[parentPos] = subIdx++;
                            if (isP) newPivotC++; else newKeepC++;
                        });

                        auto newId = tree.addNode(reusableLeaf);
                        // Add to treeGraphV
                        for (const auto &v : tree.adj_list[newId]) {
                            treeGraphV.addNbr(v.v, {newId, v.isPivot});
                        }

                        daf::Size np = s - newKeepC;

                        // Build new leafCliqueInfo from parent's entries (V11's containment check)
                        reusableEntries.clear();
                        for (const auto &entry : leafCliqueInfo[leafId]) {
                            bool contained = true;
                            daf::CliqueSize subP = 0;
                            for (int j = 0; j < (int)r; ++j) {
                                uint8_t pos = entry.positions[j];
                                if (!c.test(pos)) { contained = false; break; }
                                if (bkPivots.test(pos)) subP++;
                            }
                            if (!contained) continue;

                            if (subP <= np && newPivotC - subP < 1001 && np - subP < 401) {
                                double ncrVal = nCr[newPivotC - subP][np - subP];
                                countingRClique[entry.cliqueId] += ncrVal;

                                RCliqueSTv18::LeafCliqueEntry newEntry;
                                newEntry.cliqueId = entry.cliqueId;
                                newEntry.ncrValue = ncrVal;
                                for (int j = 0; j < (int)r; ++j)
                                    newEntry.positions[j] = parentToSubPos[entry.positions[j]];
                                reusableEntries.push_back(newEntry);
                            }
                        }

                        if (newId >= leafCliqueInfo.size())
                            leafCliqueInfo.resize(newId + 1);
                        leafCliqueInfo[newId] = reusableEntries;

                        if (newId >= changedLeafIndex.size())
                            changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                    });

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                // Subtract old leaf's contributions via leafCliqueInfo (V11's fast path)
                auto t_supp = std::chrono::high_resolution_clock::now();
                for (const auto &entry : leafCliqueInfo[leafId]) {
                    if (!rCliqueInHeap[entry.cliqueId]) continue;
                    countingRClique[entry.cliqueId] -= entry.ncrValue;
                    if (countingRClique[entry.cliqueId] < 0) countingRClique[entry.cliqueId] = 0;
                    bucketMove(entry.cliqueId);
                }
                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

                tree.removeNode(leafId);
                leafCliqueInfo[leafId].clear();
                leafCliqueInfo[leafId].shrink_to_fit();
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
    std::cout << "  Cases: LeafDeath=" << cntLeafDeath << " BK=" << cntBK
              << " iters=" << totalIters << std::endl;

    std::vector<std::pair<std::vector<daf::Size>, double>> sortedK;
    sortedK.reserve(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueCopy(clique.begin(), clique.end());
        sortedK.emplace_back(cliqueCopy, coreRClique[i]);
    }
    return sortedK;
}
