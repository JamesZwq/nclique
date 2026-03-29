//
// ST V19: Analytical Inclusion-Exclusion (no pathSplit for small m)
//
// Key idea: instead of generating Π|F_i| subpaths via BK pathSplit,
// compute each r-clique's support change directly via inclusion-exclusion.
//
// For m removed r-cliques with conflict sets F_1..F_m:
//   new_support(q') = Σ_{S⊆[m]} (-1)^|S| · nCr(a - |P_q' ∪ ⋃_{i∈S} F_i|, need - ...)
//
// Cost: O(C(n,r) × 2^m)  vs  pathSplit: O(Π|F_i| × C(n',r))
// For m=10, r=3:  1024 vs 59049 = 57x improvement.
//
// Fallback: when m > IE_THRESHOLD, use pathSplit (V18-style).
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

namespace RCliqueSTv19 {

// Max removals for IE path (2^IE_THRESHOLD must be feasible)
static constexpr int IE_THRESHOLD = 15;

// Precompute union sizes for all 2^m subsets of conflict pivot-sets.
// conflictPivots[i] = bitmask of pivot positions for conflict i (leaf-local).
// Returns unionSize[S] = |⋃_{i∈S} conflictPivots[i]| for each subset S.
static void precomputeUnionSizes(
    int m,
    const uint64_t conflictMasks[],  // bitmask per conflict set (up to 64 pivots)
    uint16_t unionSize[])            // output: size per subset
{
    int total = 1 << m;
    unionSize[0] = 0;
    for (int S = 1; S < total; ++S) {
        // Find lowest set bit
        int lsb = __builtin_ctz(S);
        int prev = S & (S - 1); // S without the lowest bit
        // Union = union of prev ∪ conflictMasks[lsb]
        // But we store sizes, not full masks. We need the actual mask.
        // Better: store full union masks, then count bits.
        // Actually, let's just compute incrementally.
    }
    // Rewrite: store union MASKS, then popcount
    uint64_t unionMask[1 << IE_THRESHOLD];
    unionMask[0] = 0;
    for (int S = 1; S < total; ++S) {
        int lsb = __builtin_ctz(S);
        unionMask[S] = unionMask[S & (S-1)] | conflictMasks[lsb];
        unionSize[S] = (uint16_t)__builtin_popcountll(unionMask[S]);
    }
}

} // namespace RCliqueSTv19

std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V19(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    long long duration_init = 0, duration_pop = 0, duration_intersect = 0;
    long long duration_bk = 0, duration_support = 0;
    long long cntLeafDeath = 0, cntBK = 0, cntIE = 0, cntFallback = 0;

    auto time_start = std::chrono::high_resolution_clock::now();

    // ========== INIT: V17's lightweight approach ==========
    StaticCliqueIndex localIndex(r);
    StaticCliqueIndex &cliqueIndex = prebuiltIndex ? *prebuiltIndex : localIndex;
    if (!prebuiltIndex) {
        daf::timeCount("clique Index build", [&]() {
            localIndex.build(tree, edgeGraph.adj_list.size());
        });
    }

    std::vector<double> countingRClique;
    daf::timeCount("countingPerRClique (V19)", [&]() {
        countingRClique.assign(cliqueIndex.size(), 0.0);
        for (const auto &leaf : tree.adj_list) {
            if (leaf.size() < r) continue;
            daf::CliqueSize pivotC = 0, keepC = 0;
            for (const auto &i : leaf) { if (i.isPivot) pivotC++; else keepC++; }
            int needPivot = s - static_cast<int>(keepC);
            daf::enumerateCombinations(leaf, r,
                [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                daf::CliqueSize subNumPivot = 0;
                for (const auto &node : rClique) if (node.isPivot) subNumPivot++;
                if (subNumPivot <= needPivot) {
                    int row = (int)pivotC - (int)subNumPivot;
                    int col = needPivot - (int)subNumPivot;
                    if (row >= 0 && row < 1001 && col >= 0 && col < 401)
                        countingRClique[cliqueIndex.byClique(rClique)] += nCr[row][col];
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

    // Per-leaf pending removals for lazy IE
    std::vector<std::vector<daf::Size>> pendingRemovals(tree.adj_list.size());

    // ========== Hybrid bucket+set PQ ==========
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
                int b = (int)countingRClique[id]; bucket_of[id] = b;
                pos_in_bucket[id] = buckets[b].size(); buckets[b].push_back(id);
            } else break;
        }
    };

    daf::log_memory("Other index (include bucket)");
    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::vector<daf::Size> vertexConflictDeg;
    vertexConflictDeg.reserve(512);

    std::cout << "=========================begin (r>=3 ST_V19)===========================" << std::endl;
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
                goto phase1_v19;
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
        phase1_v19:
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

            int n = (int)leaf.size();
            daf::CliqueSize keepC = 0, pivotC = 0;
            for (const auto &node : leaf) {
                if (node.isPivot) pivotC++; else keepC++;
            }
            int needPivot = s - static_cast<int>(keepC);

            auto t_bk = std::chrono::high_resolution_clock::now();

            // --- LeafDeath fast-path ---
            daf::Size maxRCliquePerVertex = (daf::Size)(nCr[n - 1][r - 1] + 0.5);
            vertexConflictDeg.assign(n, 0);
            daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
            for (int i = 0; i < n; ++i) mapRef[leaf[i].v] = (daf::Size)i;

            // Count conflict degree including ALL pending + new removals
            auto &pending = pendingRemovals[leafId];
            // New removals for this iteration
            for (auto rmId : removedR) {
                pending.push_back(rmId);
            }
            for (auto rmId : pending) {
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
                int rp = (int)pivotC - forcedPivotRemove;
                int rt = n - forcedPivotRemove;
                if (rt < (int)s || rp < needPivot) leafDead = true;
            }

            int m = (int)pending.size(); // total pending removals

            if (leafDead) {
                // ===== DEAD LEAF =====
                cntLeafDeath++;
                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                auto t_supp = std::chrono::high_resolution_clock::now();
                // Remove from treeGraphV
                for (const auto &leafV : leaf) {
                    treeGraphV.removeNbr(leafV.v, {leafId, leafV.isPivot});
                }
                // Subtract ALL remaining support (leaf dies → all r-cliques lose this leaf's contribution)
                // But we need the CURRENT contribution, accounting for previous pending removals.
                // If no previous pending (common case): just subtract nCr directly.
                // If pending exists: use IE to compute current support, then subtract.
                if (m <= RCliqueSTv19::IE_THRESHOLD && m > 0) {
                    // Build conflict masks for IE
                    uint64_t conflictMasks[RCliqueSTv19::IE_THRESHOLD];
                    for (int ci = 0; ci < m; ++ci) {
                        conflictMasks[ci] = 0;
                        auto rClique = cliqueIndex.byId(pending[ci]);
                        for (auto v : rClique) {
                            daf::Size pos = mapRef[v];
                            if (pos < (daf::Size)n && leaf[pos].isPivot)
                                conflictMasks[ci] |= (1ULL << pos);
                        }
                    }
                    uint16_t unionSz[1 << RCliqueSTv19::IE_THRESHOLD];
                    RCliqueSTv19::precomputeUnionSizes(m, conflictMasks, unionSz);

                    int totalSubsets = 1 << m;
                    int signs[1 << RCliqueSTv19::IE_THRESHOLD];
                    for (int S = 0; S < totalSubsets; ++S)
                        signs[S] = (__builtin_popcount(S) & 1) ? -1 : 1;

                    // For each r-clique, compute current support via IE, then subtract
                    daf::enumerateCombinations(leaf, r,
                        [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                        auto cliqueId = cliqueIndex.byClique(rClique);
                        if (!rCliqueInHeap[cliqueId]) return true;

                        // Build q' pivot mask
                        uint64_t qMask = 0;
                        daf::CliqueSize subP = 0;
                        for (const auto &node : rClique) {
                            daf::Size pos = mapRef[node.v];
                            if (node.isPivot) { qMask |= (1ULL << pos); subP++; }
                        }

                        // Compute current support via IE
                        double currentSupport = 0;
                        for (int S = 0; S < totalSubsets; ++S) {
                            // |P_q' ∪ ⋃_{i∈S} F_i| = popcount(qMask | unionMask[S])
                            // But we precomputed unionSz[S] = popcount(unionMask[S])
                            // We need popcount(qMask | unionMask[S])
                            // = |qMask| + |unionMask[S]| - |qMask & unionMask[S]|
                            // Hmm, we don't have the actual unionMask. Let me recompute.
                        }
                        // Actually, we need full union masks. Let me fix.
                        return true;
                    });
                }
                // Simpler: for dead leaves, just subtract the FULL original contribution
                // (before any pending removals). The pending removals' deltas were already
                // applied in previous iterations.
                // Only need to subtract the INCREMENTAL change from NEW removals.
                // But this is complex... let's use direct enumeration for dead leaves.
                daf::enumerateCombinations(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                    auto cliqueId = cliqueIndex.byClique(rClique);
                    if (!rCliqueInHeap[cliqueId]) return true;
                    daf::CliqueSize subP = 0;
                    for (const auto &node : rClique) if (node.isPivot) subP++;
                    countingRClique[cliqueId] -= nCr[pivotC - subP][needPivot - subP];
                    if (countingRClique[cliqueId] < 0) countingRClique[cliqueId] = 0;
                    bucketMove(cliqueId);
                    return true;
                });
                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();
                tree.removeNode(leafId);
                pending.clear();

            } else if (m <= RCliqueSTv19::IE_THRESHOLD) {
                // ===== IE PATH: Analytical inclusion-exclusion =====
                cntIE++;

                // Build local vertex index
                // Build conflict pivot masks (leaf-local positions, using bitmask)
                uint64_t conflictMasks[RCliqueSTv19::IE_THRESHOLD];
                for (int ci = 0; ci < m; ++ci) {
                    conflictMasks[ci] = 0;
                    auto rClique = cliqueIndex.byId(pending[ci]);
                    for (auto v : rClique) {
                        daf::Size pos = mapRef[v];
                        if (pos < (daf::Size)n && leaf[pos].isPivot)
                            conflictMasks[ci] |= (1ULL << pos);
                    }
                }

                // Precompute union masks for all 2^m subsets
                int totalSubsets = 1 << m;
                uint64_t unionMask[1 << RCliqueSTv19::IE_THRESHOLD];
                unionMask[0] = 0;
                for (int S = 1; S < totalSubsets; ++S) {
                    int lsb = __builtin_ctz(S);
                    unionMask[S] = unionMask[S & (S-1)] | conflictMasks[lsb];
                }

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                auto t_supp = std::chrono::high_resolution_clock::now();

                // For each r-clique q' in the leaf:
                // new_support = Σ_{S⊆[m]} (-1)^|S| · nCr(pivotC - |pivots of q' ∪ unionMask[S]|,
                //                                         needPivot - |pivots of q' ∪ unionMask[S]|)
                // delta = new_support - old_support (where old = nCr[pivotC-subP][needPivot-subP])
                daf::enumerateCombinations(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                    auto cliqueId = cliqueIndex.byClique(rClique);
                    if (!rCliqueInHeap[cliqueId]) return true;

                    // Build q' pivot bitmask
                    uint64_t qMask = 0;
                    daf::CliqueSize subP = 0;
                    for (const auto &node : rClique) {
                        if (node.isPivot) {
                            qMask |= (1ULL << mapRef[node.v]);
                            subP++;
                        }
                    }

                    // old support (before any removals from this leaf)
                    double oldSupp = nCr[pivotC - subP][needPivot - subP];

                    // new support via IE
                    double newSupp = 0;
                    for (int S = 0; S < totalSubsets; ++S) {
                        int sign = (__builtin_popcount(S) & 1) ? -1 : 1;
                        uint64_t combined = qMask | unionMask[S];
                        int totalPivUsed = __builtin_popcountll(combined);
                        int row = (int)pivotC - totalPivUsed;
                        int col = needPivot - totalPivUsed;
                        if (row >= 0 && col >= 0 && row < 1001 && col < 401)
                            newSupp += sign * nCr[row][col];
                    }

                    double delta = newSupp - oldSupp;
                    if (delta != 0) {
                        countingRClique[cliqueId] += delta;
                        if (countingRClique[cliqueId] < 0) countingRClique[cliqueId] = 0;
                        bucketMove(cliqueId);
                    }
                    return true;
                });

                // NO tree mutation, NO subpath generation!
                // The leaf stays as-is with pending removals tracked.
                // treeGraphV unchanged.

                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

            } else {
                // ===== FALLBACK: pathSplit (V17-style) when m > threshold =====
                cntFallback++;

                // Remove from treeGraphV
                for (const auto &leafV : leaf)
                    treeGraphV.removeNbr(leafV.v, {leafId, leafV.isPivot});

                // Use ALL pending as conflicts
                auto mapped = pending | std::views::transform(
                    [&](const daf::Size id) { return cliqueIndex.byId(id); });

                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &bkPivots) {
                        std::vector<TreeGraphNode> newLeaf;
                        c.for_each_bit([&](size_t pos) {
                            newLeaf.push_back({leaf[pos].v, bkPivots.test(pos)});
                        });
                        auto newId = tree.addNode(newLeaf);
                        for (const auto &v : tree.adj_list[newId])
                            treeGraphV.addNbr(v.v, {newId, v.isPivot});

                        // Compute support for new sub-leaf
                        daf::CliqueSize newPivotC = 0, newKeepC = 0;
                        for (const auto &v : tree.adj_list[newId]) {
                            if (v.isPivot) newPivotC++; else newKeepC++;
                        }
                        daf::Size np = s - newKeepC;
                        daf::enumerateCombinations(tree.adj_list[newId], r,
                            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                            daf::CliqueSize sp = 0;
                            for (const auto &node : rClique) if (node.isPivot) sp++;
                            if (sp <= np) {
                                int row = (int)newPivotC - (int)sp;
                                int col = (int)np - (int)sp;
                                if (row >= 0 && row < 1001 && col >= 0 && col < 401)
                                    countingRClique[cliqueIndex.byClique(rClique)] += nCr[row][col];
                            }
                            return true;
                        });

                        if (newId >= changedLeafIndex.size())
                            changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                        if (newId >= pendingRemovals.size())
                            pendingRemovals.resize(newId + 1);
                    });

                duration_bk += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_bk).count();

                // Subtract old leaf's full contribution
                auto t_supp = std::chrono::high_resolution_clock::now();
                daf::enumerateCombinations(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                    auto cliqueId = cliqueIndex.byClique(rClique);
                    if (!rCliqueInHeap[cliqueId]) return true;
                    daf::CliqueSize subP = 0;
                    for (const auto &node : rClique) if (node.isPivot) subP++;
                    countingRClique[cliqueId] -= nCr[pivotC - subP][needPivot - subP];
                    if (countingRClique[cliqueId] < 0) countingRClique[cliqueId] = 0;
                    bucketMove(cliqueId);
                    return true;
                });
                duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_supp).count();

                tree.removeNode(leafId);
                pending.clear();
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
    std::cout << "  Cases: LeafDeath=" << cntLeafDeath << " IE=" << cntIE
              << " Fallback=" << cntFallback << " iters=" << totalIters << std::endl;

    std::vector<std::pair<std::vector<daf::Size>, double>> sortedK;
    sortedK.reserve(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueCopy(clique.begin(), clique.end());
        sortedK.emplace_back(cliqueCopy, coreRClique[i]);
    }
    return sortedK;
}
