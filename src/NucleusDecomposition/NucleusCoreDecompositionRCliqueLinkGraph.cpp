//
// Framework 1: Link-Graph Peeling
// THIS IS A NEW FRAMEWORK — NOT an SDCT/peeling optimization.
//
// Core idea: Instead of peeling r-cliques from an SDCT (which requires BK
// re-splitting), we compute pairwise co-occurrence weights on-the-fly using
// the IMMUTABLE SDCT and do standard weighted bucket peeling.
//
// When r-clique C is removed, for each co-leaf r-clique D:
//   Δ(D,C,leaf) = nCr[pivotC - unionPivots(C,D)][needPivot - unionPivots(C,D)]
// This is the exact number of s-cliques containing BOTH C and D in that leaf.
// No BK, no tree mutation, no sub-leaf enumeration.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <cstring>
#include <vector>

#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace LinkGraphPeel {

struct LeafCliqueEntry {
    daf::Size cliqueId;
    double ncrValue; // single-clique nCr contribution (for init support)
};

struct DualIndex {
    std::vector<std::vector<LeafCliqueEntry>> leafCliqueInfo;
    std::vector<std::vector<daf::Size>> cliqueLeafIds;
    std::vector<double> counting;
};

DualIndex buildDualIndex(
    const DynamicGraph<TreeGraphNode> &treeGraph,
    const StaticCliqueIndex &cliqueIndex,
    daf::CliqueSize r, daf::CliqueSize s) {

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
            daf::CliqueSize subP = 0;
            for (const auto &node : rClique) if (node.isPivot) subP++;
            double ncrValue = nCr[pivotC - subP][needPivot - subP];
            auto id = cliqueIndex.byClique(rClique);
            if (id < nClique) {
                result.counting[id] += ncrValue;
                result.leafCliqueInfo[leafId].push_back({id, ncrValue});
                result.cliqueLeafIds[id].push_back(leafId);
            }
            return true;
        });
    }
    return result;
}

// Compute pairwise co-occurrence weight: number of s-cliques containing
// BOTH r-cliques C and D in the same leaf.
// Uses the SDCT leaf's pivot/keep structure.
inline long long pairwiseWeight(
    const std::span<const daf::Size> &cVerts,
    const std::span<const daf::Size> &dVerts,
    const std::vector<TreeGraphNode> &leaf,
    daf::StaticVector<daf::Size> &mapRef, // vertex -> leaf position
    int pivotC, int needPivot)
{
    // Count unique pivot vertices in C ∪ D
    int unionPivots = 0;
    int unionSize = 0;

    // C's contribution
    for (auto v : cVerts) {
        daf::Size pos = mapRef[v];
        if (leaf[pos].isPivot) unionPivots++;
        unionSize++;
    }
    // D's contribution (skip vertices already in C)
    for (auto v : dVerts) {
        daf::Size pos = mapRef[v];
        // Check if v is already in C (positions would have been set)
        bool inC = false;
        for (auto cv : cVerts) {
            if (cv == v) { inC = true; break; }
        }
        if (!inC) {
            if (leaf[pos].isPivot) unionPivots++;
            unionSize++;
        }
    }

    // Number of s-cliques containing both C and D:
    // All keeps are always included. We need (needPivot - unionPivots) more
    // pivots from the remaining (pivotC - unionPivots) pivots.
    int remainPivots = pivotC - unionPivots;
    int remainNeed = needPivot - unionPivots;

    if (remainNeed < 0 || remainPivots < 0 || remainPivots < remainNeed) return 0;
    if (remainPivots >= 1001 || remainNeed >= 401) return 0;

    return (double)(nCr[remainPivots][remainNeed] + 0.5);
}

} // namespace LinkGraphPeel


std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLinkPeel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    auto time_start = std::chrono::high_resolution_clock::now();

    // === Init: same as V4 ===
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build (LinkPeel)", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });

    auto dualIndex = daf::timeCount("buildDualIndex (LinkPeel)", [&]() {
        return LinkGraphPeel::buildDualIndex(tree, cliqueIndex, r, s);
    });

    auto &leafCliqueInfo = dualIndex.leafCliqueInfo;
    auto &cliqueLeafIds = dualIndex.cliqueLeafIds;
    auto countingRClique = std::move(dualIndex.counting);

    const daf::Size nClique = cliqueIndex.size();
    std::vector<double> coreRClique(nClique, 0);

    // Precompute per-leaf metadata (immutable — tree never changes!)
    const daf::Size numLeaves = tree.adj_list.size();
    std::vector<int> leafPivotC(numLeaves), leafNeedPivot(numLeaves);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        int keeps = 0, pivots = 0;
        for (const auto &node : tree.adj_list[L]) {
            if (node.isPivot) pivots++;
            else keeps++;
        }
        leafPivotC[L] = pivots;
        leafNeedPivot[L] = (int)s - keeps;
    }

    // Build position map for each leaf's vertices (reusable)
    daf::StaticVector<daf::Size> &mapRef = daf::vListMap;

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

    auto duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::cout << "=========================begin (Link-Graph Peeling r>=" << (int)r << ")=========================" << std::endl;

    double minCore = 0;
    long long totalIters = 0;
    long long duration_pop = 0, duration_update = 0;

    while (remainingInHeap > 0) {
        auto t_pop = std::chrono::high_resolution_clock::now();

        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((double)curBucket, minCore);

        // Pop ONE r-clique at a time to avoid double-counting
        auto id = buckets[curBucket].back();
        buckets[curBucket].pop_back();
        rCliqueInHeap[id] = false;
        coreRClique[id] = minCore;
        remainingInHeap--;

        duration_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_pop).count();

        if (remainingInHeap == 0) break;
        totalIters++;

        // Update support of co-leaf neighbors
        auto t_update = std::chrono::high_resolution_clock::now();

        auto cVerts = cliqueIndex.byId(id);

        for (auto leafId : cliqueLeafIds[id]) {
            const auto &leaf = tree.adj_list[leafId];
            if (leaf.empty()) continue;

            int pC = leafPivotC[leafId];
            int nP = leafNeedPivot[leafId];
            int n = (int)leaf.size();

            // Set up mapRef for this leaf
            for (int i = 0; i < n; ++i)
                mapRef[leaf[i].v] = (daf::Size)i;

            // Count C's pivots in this leaf
            int cPivots = 0;
            for (auto v : cVerts) {
                if (leaf[mapRef[v]].isPivot) cPivots++;
            }

            // For each OTHER r-clique D in this leaf:
            for (const auto &entry : leafCliqueInfo[leafId]) {
                if (entry.cliqueId == id) continue;
                if (!rCliqueInHeap[entry.cliqueId]) continue;

                auto dVerts = cliqueIndex.byId(entry.cliqueId);

                // Count union pivots of C ∪ D
                int unionPivots = cPivots;
                for (auto v : dVerts) {
                    bool inC = false;
                    for (auto cv : cVerts) {
                        if (cv == v) { inC = true; break; }
                    }
                    if (!inC && leaf[mapRef[v]].isPivot) {
                        unionPivots++;
                    }
                }

                int remainPivots = pC - unionPivots;
                int remainNeed = nP - unionPivots;
                if (remainNeed < 0 || remainPivots < 0 || remainPivots < remainNeed) continue;
                if (remainPivots >= 1001 || remainNeed >= 401) continue;

                long long w = (double)(nCr[remainPivots][remainNeed] + 0.5);
                if (w <= 0) continue;

                countingRClique[entry.cliqueId] -= w;
                if (countingRClique[entry.cliqueId] < 0) countingRClique[entry.cliqueId] = 0;
                bucketMove(entry.cliqueId);
            }
        }

        duration_update += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_update).count();
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::cout << "time: " << elapsed << " ms" << std::endl;
    std::cout << "Time Breakdown (ms):" << std::endl;
    std::cout << "  Init:      " << duration_init / 1000000.0 << std::endl;
    std::cout << "  Pop:       " << duration_pop / 1000000.0 << std::endl;
    std::cout << "  Update:    " << duration_update / 1000000.0 << std::endl;
    std::cout << "  iters=" << totalIters << std::endl;

    // Build output
    std::vector<std::pair<std::vector<daf::Size>, double>> sortedK;
    sortedK.reserve(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> copy(clique.begin(), clique.end());
        sortedK.emplace_back(std::move(copy), (int)coreRClique[i]);
    }
    std::sort(sortedK.begin(), sortedK.end(),
              [](const auto &a, const auto &b) { return a.second < b.second; });
    return sortedK;
}
