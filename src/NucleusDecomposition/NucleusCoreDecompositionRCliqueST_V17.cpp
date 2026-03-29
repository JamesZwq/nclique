//
// ST V17: Hybrid — Reference's lightweight init + V11's bucket PQ
//
// Key idea: Reference's init (no dual index) eliminates the 65-75% Init overhead.
// V11's hybrid bucket+set PQ gives O(1) peeling ops.
// BK pruning (singleton, subsumption, component LB) comes from shared BronKerboschRmRClique.hpp.
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

std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V17(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    long long duration_init = 0, duration_pop = 0, duration_intersect = 0;
    long long duration_bk = 0, duration_support = 0;
    long long cntLeafDeath = 0, cntBK = 0;

    auto time_start = std::chrono::high_resolution_clock::now();

    // ========== INIT: Reference's lightweight approach ==========
    // Use prebuilt index or build locally — NO dual index, NO leafCliqueInfo
    StaticCliqueIndex localIndex(r);
    StaticCliqueIndex &cliqueIndex = prebuiltIndex ? *prebuiltIndex : localIndex;
    if (!prebuiltIndex) {
        daf::timeCount("clique Index build", [&]() {
            localIndex.build(tree, edgeGraph.adj_list.size());
        });
    }

    // Count support: iterate tree, enumerate r-cliques per leaf, accumulate nCr
    // (Reference's approach — no dual index needed)
    std::vector<double> countingRClique;
    daf::timeCount("countingPerRClique (V17)", [&]() {
        countingRClique.assign(cliqueIndex.size(), 0.0);
        for (const auto &leaf : tree.adj_list) {
            if (leaf.size() < r) continue;
            daf::CliqueSize pivotC = 0, keepC = 0;
            for (const auto &i : leaf) {
                if (i.isPivot) pivotC++; else keepC++;
            }
            int needPivot = s - static_cast<int>(keepC);
            daf::enumerateCombinations(leaf, r,
                [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                daf::CliqueSize subNumPivot = 0;
                for (const auto &node : rClique)
                    if (node.isPivot) subNumPivot++;
                if (subNumPivot <= needPivot) {
                    auto ncrValue = nCr[pivotC - subNumPivot][needPivot - subNumPivot];
                    auto cliqueId = cliqueIndex.byClique(rClique);
                    countingRClique[cliqueId] += ncrValue;
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
        if (oldB == -1)
            overflowSet.erase({overflowStoredVal[id], id});
        if (val <= BUCKET_THRESHOLD) {
            int newB = (int)val;
            if (oldB >= 0 && newB == oldB) return;
            if (oldB >= 0) {
                auto &v = buckets[oldB];
                daf::Size p = pos_in_bucket[id];
                if (p < v.size() - 1) { auto last = v.back(); v[p] = last; pos_in_bucket[last] = p; }
                v.pop_back();
            }
            bucket_of[id] = newB;
            pos_in_bucket[id] = buckets[newB].size();
            buckets[newB].push_back(id);
            if (newB < curBucket) curBucket = newB;
        } else {
            if (oldB >= 0) {
                auto &v = buckets[oldB];
                daf::Size p = pos_in_bucket[id];
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
                bucket_of[id] = b;
                pos_in_bucket[id] = buckets[b].size();
                buckets[b].push_back(id);
            } else break;
        }
    };

    daf::log_memory("Other index (include bucket)");

    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::cout << "=========================begin (r>=3 ST_V17)===========================" << std::endl;
    double minCore = 0;
    long long totalIters = 0;

    while (remainingInHeap > 0) {
        auto t_loop_start = std::chrono::high_resolution_clock::now();

        for (auto &leafId : changedLeaf)
            changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        changedLeaf.clear();
        removedRCliqueIdForLeaf.clear();
        currentRemoveRcliqueIds.clear();

        // --- Drain + Pop ---
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
                goto phase1_v17;
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
                && (curBucket + 1) <= (int)minCore)
                curBucket++;
            else break;
        }
        phase1_v17:
        duration_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_loop_start).count();
        if (remainingInHeap == 0) break;
        totalIters++;

        // --- Phase 1: Intersect via treeGraphV (Reference's approach) ---
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

        // --- Phase 2: BK split + support update (Reference's approach + bucket PQ) ---
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

            auto t_bk = std::chrono::high_resolution_clock::now();

            // --- initCore callback: Reference's approach ---
            // Called by BK for each new sub-leaf
            auto initCore = [&](const std::vector<TreeGraphNode> &newLeaf, const daf::Size &newLeafId) {
                daf::CliqueSize newPivotC = 0, newKeepC = 0;
                for (const auto &i : newLeaf) {
                    if (i.isPivot) {
                        treeGraphV.addNbr(i.v, {newLeafId, true});
                        newPivotC++;
                    } else {
                        treeGraphV.addNbr(i.v, {newLeafId, false});
                        newKeepC++;
                    }
                }
                daf::Size needPivot = s - newKeepC;

                daf::enumerateCombinations(newLeaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &rclique) {
                    daf::CliqueSize subNumPivot = 0;
                    for (const auto &node : rclique)
                        if (node.isPivot) subNumPivot++;
                    if (subNumPivot <= needPivot) {
                        int row = (int)newPivotC - (int)subNumPivot;
                        int col = (int)needPivot - (int)subNumPivot;
                        if (row >= 0 && row < 1001 && col >= 0 && col < 401) {
                            auto ncrValue = nCr[row][col];
                            auto cliqueId = cliqueIndex.byClique(rclique);
                            countingRClique[cliqueId] += ncrValue;
                        }
                    }
                    return true;
                });
            };

            // Remove old leaf from treeGraphV
            for (const auto &leafV : leaf) {
                if (leafV.isPivot)
                    treeGraphV.removeNbr(leafV.v, {leafId, true});
                else
                    treeGraphV.removeNbr(leafV.v, {leafId, false});
            }

            // BK split (shared code with all pruning)
            auto mapped = removedR | std::views::transform(
                [&](const daf::Size id) { return cliqueIndex.byId(id); });

            bkRmClique::removeRClique(leaf, mapped, r, s,
                [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &bkPivots) {
                    cntBK++;
                    auto newLeaf = bkRmClique::coverToVertex(c, bkPivots, leaf);
                    auto newId = tree.addNode(newLeaf);
                    initCore(tree.adj_list[newId], newId);
                    if (newId >= changedLeafIndex.size())
                        changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                });

            duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - t_bk).count();

            // Subtract old leaf's contributions + bucketMove
            auto t_supp = std::chrono::high_resolution_clock::now();
            daf::enumerateCombinations(leaf, r,
                [&](const daf::StaticVector<TreeGraphNode> &clique) {
                auto cliqueId = cliqueIndex.byClique(clique);
                if (!rCliqueInHeap[cliqueId]) return true;
                daf::CliqueSize subNumPivot = 0;
                for (const auto &node : clique)
                    if (node.isPivot) subNumPivot++;
                auto ncrValue = nCr[pivotC - subNumPivot][s - keepC - subNumPivot];
                countingRClique[cliqueId] -= ncrValue;
                if (countingRClique[cliqueId] < 0) countingRClique[cliqueId] = 0;
                bucketMove(cliqueId);
                return true;
            });
            duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - t_supp).count();

            tree.removeNode(leafId);
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
    std::cout << "  Cases: BK=" << cntBK << " iters=" << totalIters << std::endl;

    std::vector<std::pair<std::vector<daf::Size>, double>> sortedK;
    sortedK.reserve(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueCopy(clique.begin(), clique.end());
        sortedK.emplace_back(cliqueCopy, coreRClique[i]);
    }
    return sortedK;
}
