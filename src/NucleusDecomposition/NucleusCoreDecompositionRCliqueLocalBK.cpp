//
// Local H-index for r≥3 with exact s-clique enumeration (Online BK).
//
// Tree is immutable. Each r-clique independently computes its core value
// via H-index. Instead of vertex-proxy, we enumerate ALL s-cliques
// containing the target r-clique within each leaf and compute exact
// support at each threshold level.
//
// For each s-clique S containing target r-clique C:
//   level = min(coreC[C'] for all other r-cliques C' in S)
//   bucket[level] += 1
// Then scan buckets for H-index.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <vector>
#include <cstring>
#include <queue>
#include <cmath>

#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace RCliqueLocalBK {

// Initial support counting (same as RCliqueST)
std::vector<long long> countingPerRClique(
    const DynamicGraph<TreeGraphNode> &treeGraph,
    const StaticCliqueIndex &cliqueIndex,
    const daf::CliqueSize r,
    const daf::CliqueSize s) {
    const daf::Size nClique = cliqueIndex.size();
    std::vector<long long> counting(nClique, 0);
    for (const auto &leaf : treeGraph.adj_list) {
        if (leaf.size() < r) continue;
        daf::CliqueSize pivotC = 0, keepC = 0;
        for (const auto &i : leaf) {
            if (i.isPivot) pivotC++;
            else keepC++;
        }
        int needPivot = s - static_cast<int>(keepC);
        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            daf::CliqueSize subNumPovit = 0;
            for (const auto &node : rClique) if (node.isPivot) subNumPovit++;
            long long ncrValue = (long long)(nCr[pivotC - subNumPovit][needPivot - subNumPovit] + 0.5);
            auto id = cliqueIndex.byClique(rClique);
            if (id < nClique) counting[id] += ncrValue;
            return true;
        });
    }
    return counting;
}

} // namespace RCliqueLocalBK


// ============================================================
// Compute H-index for r-clique C using exact s-clique enumeration.
// For each leaf containing all vertices of C:
//   - others = leaf \ C, classified as keeps/pivots
//   - enumerate all (s-r)-subsets from others (respecting keep/pivot)
//   - for each such s-clique, find min coreC among other r-cliques
//   - bucket[level] += 1
// ============================================================
static long long computeRCliqueHIndexBK(
    daf::Size cliqueId,
    const StaticCliqueIndex &cliqueIndex,
    const DynamicGraph<TreeGraphNode> &tree,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    const long long *coreC,
    const daf::CliqueSize r,
    const daf::CliqueSize s,
    const std::vector<int> &leafNeedPivot,
    const std::vector<uint8_t> &leafValid,
    long long currentCore,
    long long &totalSupportOut,
    std::vector<long long> &buckets,
    std::vector<TreeGraphNode> &otherKeeps,
    std::vector<TreeGraphNode> &otherPivots,
    std::vector<TreeGraphNode> &sCliqueOthers)
{
    totalSupportOut = 0;
    if (currentCore <= 0) return 0;

    int bucketSize = (int)currentCore;
    if (bucketSize < 1) bucketSize = 1;
    if ((int)buckets.size() < bucketSize + 1)
        buckets.resize(bucketSize + 1, 0);
    for (int i = 0; i <= bucketSize; ++i) buckets[i] = 0;

    long long rawTotalSupport = 0;

    auto cliqueVerts = cliqueIndex.byId(cliqueId);
    const daf::Size nClique = cliqueIndex.size();

    // Find leaves containing ALL vertices of C
    daf::intersect_dense_sets_multi(cliqueVerts, treeGraphV.adj_list,
        [&](const TreeGraphNode &leafNode) {
            daf::Size leafId = leafNode.v;
            if (!leafValid[leafId]) return;

            const auto &leaf = tree.adj_list[leafId];
            int needPivot = leafNeedPivot[leafId];

            // Build set of clique vertex ids for fast lookup
            // and count how many of C's vertices are pivots in this leaf
            int cPivotsInLeaf = 0;
            robin_hood::unordered_flat_set<daf::Size> cliqueSet;
            for (auto v : cliqueVerts) cliqueSet.insert(v);

            // Find pivot status of C's vertices in this leaf
            for (const auto &node : leaf) {
                if (cliqueSet.count(node.v) && node.isPivot) cPivotsInLeaf++;
            }

            // Classify other vertices
            otherKeeps.clear();
            otherPivots.clear();
            for (const auto &node : leaf) {
                if (cliqueSet.count(node.v)) continue;
                if (node.isPivot) otherPivots.push_back(node);
                else otherKeeps.push_back(node);
            }

            // effectiveNeedPivot for the (s-r) additional vertices
            int effectiveNeedPivot = needPivot - cPivotsInLeaf;
            int numOtherKeeps = (int)otherKeeps.size();
            int numOtherPivots = (int)otherPivots.size();
            int need = (int)s - (int)r; // how many more vertices we need

            // We need all keeps + effectiveNeedPivot pivots from others
            // Total others needed = s - r
            // keeps contribute: numOtherKeeps (all must be included)
            // pivots needed: need - numOtherKeeps = effectiveNeedPivot
            if (effectiveNeedPivot < 0 || effectiveNeedPivot > numOtherPivots) return;
            if (numOtherKeeps + numOtherPivots < need) return;
            if (numOtherKeeps > need) return; // can't fit all keeps

            // Helper: given chosen pivots, build s-clique, find min other r-clique core
            auto processSClique = [&](const TreeGraphNode *chosenPivotData, int chosenPivotCount) {
                sCliqueOthers.clear();
                for (const auto &k : otherKeeps) sCliqueOthers.push_back(k);
                for (int i = 0; i < chosenPivotCount; ++i) sCliqueOthers.push_back(chosenPivotData[i]);

                daf::StaticVector<TreeGraphNode> fullSClique;
                for (const auto &node : leaf) {
                    if (cliqueSet.count(node.v)) fullSClique.push_back(node);
                }
                for (const auto &node : sCliqueOthers) fullSClique.push_back(node);
                // Must be sorted by vertex ID for tryByClique to work
                std::sort(fullSClique.begin(), fullSClique.end());

                long long minOtherCore = std::numeric_limits<long long>::max();
                daf::enumerateCombinations(fullSClique, r,
                    [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                        auto id = cliqueIndex.tryByClique(rClique);
                        if (id >= nClique) return true;
                        if (id == cliqueId) return true; // skip self
                        long long c = coreC[id];
                        if (c < minOtherCore) minOtherCore = c;
                        return true;
                    });
                fullSClique.free();

                if (minOtherCore == std::numeric_limits<long long>::max()) return;
                if (minOtherCore < 1) return;

                int level = (minOtherCore > bucketSize) ? bucketSize : (int)minOtherCore;
                buckets[level] += 1;
                rawTotalSupport += 1;
            };

            if (effectiveNeedPivot == 0) {
                processSClique(nullptr, 0);
            } else {
                daf::enumerateCombinations(otherPivots, effectiveNeedPivot,
                    [&](const daf::StaticVector<TreeGraphNode> &chosenPivots) {
                        processSClique(chosenPivots.getData(), (int)chosenPivots.size());
                        return true;
                    });
            }
        });

    totalSupportOut = rawTotalSupport;
    if (rawTotalSupport < 1) return 0;

    // Bucket H-index scan
    long long accumulated = 0;
    for (int c = bucketSize; c >= 1; --c) {
        accumulated += buckets[c];
        if (accumulated >= c) return c;
    }
    return 0;
}


// ============================================================
// r≥3 Local H-index with exact enumeration (Online BK)
// ============================================================
std::vector<std::pair<std::vector<daf::Size>, int>>
NucleusCoreDecompositionRCliqueLocal_BK(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numLeaves = tree.adj_list.size();

    // Build clique index
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build (LocalBK)", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });
    const daf::Size numCliques = cliqueIndex.size();

    // Initial support
    auto countingVec = RCliqueLocalBK::countingPerRClique(tree, cliqueIndex, r, s);
    auto *coreC = new long long[numCliques];
    for (daf::Size i = 0; i < numCliques; ++i)
        coreC[i] = countingVec[i];

    // Precompute per-leaf metadata
    std::vector<int> leafNeedPivot(numLeaves);
    std::vector<uint8_t> leafValid(numLeaves);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        int keeps = 0, pivots = 0;
        for (const auto &node : tree.adj_list[L]) {
            if (node.isPivot) pivots++;
            else keeps++;
        }
        leafNeedPivot[L] = (int)s - keeps;
        leafValid[L] = (leafNeedPivot[L] >= 0 && leafNeedPivot[L] <= pivots) ? 1 : 0;
    }

    // Scratch buffers
    std::vector<long long> buckets;
    std::vector<TreeGraphNode> otherKeeps, otherPivots, sCliqueOthers;
    buckets.reserve(4096);
    otherKeeps.reserve(64);
    otherPivots.reserve(64);
    sCliqueOthers.reserve(64);

    // Dirty queue
    std::vector<uint8_t> inQueue(numCliques, 0);
    std::queue<daf::Size> dirtyQueue;

    std::cout << "=========================begin (Local H-index r>=" << (int)r << " BK)=========================" << std::endl;
    std::cout << "r-cliques: " << numCliques << ", leaves: " << numLeaves << std::endl;

    auto initEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Init time: " << std::chrono::duration_cast<std::chrono::milliseconds>(initEnd - time_start).count() << " ms" << std::endl;

    // Phase 1: initial full scan
    for (daf::Size cid = 0; cid < numCliques; ++cid) {
        if (coreC[cid] <= 0) continue;

        long long totalSup = 0;
        long long newCore = computeRCliqueHIndexBK(
            cid, cliqueIndex, tree, treeGraphV, coreC, r, s,
            leafNeedPivot, leafValid, coreC[cid], totalSup,
            buckets, otherKeeps, otherPivots, sCliqueOthers);
        newCore = std::min(newCore, coreC[cid]);

        if (newCore != coreC[cid]) {
            long long oldCore = coreC[cid];
            coreC[cid] = newCore;
            // Enqueue other r-cliques in affected leaves
            auto cliqueVerts = cliqueIndex.byId(cid);
            daf::intersect_dense_sets_multi(cliqueVerts, treeGraphV.adj_list,
                [&](const TreeGraphNode &leafNode) {
                    daf::Size leafId = leafNode.v;
                    if (!leafValid[leafId]) return;
                    const auto &leaf = tree.adj_list[leafId];
                    daf::enumerateCombinations(leaf, r,
                        [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                            auto id = cliqueIndex.tryByClique(rClique);
                            if (id >= numCliques || id == cid) return true;
                            if (coreC[id] > 0 && oldCore >= coreC[id] && !inQueue[id]) {
                                inQueue[id] = 1;
                                dirtyQueue.push(id);
                            }
                            return true;
                        });
                });
        }
    }

    // Phase 2: dirty queue propagation
    long long recomputeCount = 0;
    int iteration = 1;

    while (!dirtyQueue.empty()) {
        iteration++;
        size_t roundSize = dirtyQueue.size();
        for (size_t qi = 0; qi < roundSize; ++qi) {
            daf::Size cid = dirtyQueue.front();
            dirtyQueue.pop();
            inQueue[cid] = 0;
            recomputeCount++;

            if (coreC[cid] <= 0) continue;

            long long totalSup = 0;
            long long newCore = computeRCliqueHIndexBK(
                cid, cliqueIndex, tree, treeGraphV, coreC, r, s,
                leafNeedPivot, leafValid, coreC[cid], totalSup,
                buckets, otherKeeps, otherPivots, sCliqueOthers);
            newCore = std::min(newCore, coreC[cid]);

            if (newCore != coreC[cid]) {
                long long oldCore = coreC[cid];
                coreC[cid] = newCore;
                auto cliqueVerts = cliqueIndex.byId(cid);
                daf::intersect_dense_sets_multi(cliqueVerts, treeGraphV.adj_list,
                    [&](const TreeGraphNode &leafNode) {
                        daf::Size leafId = leafNode.v;
                        if (!leafValid[leafId]) return;
                        const auto &leaf = tree.adj_list[leafId];
                        daf::enumerateCombinations(leaf, r,
                            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                                auto id = cliqueIndex.tryByClique(rClique);
                                if (id >= numCliques || id == cid) return true;
                                if (coreC[id] > 0 && oldCore >= coreC[id] && !inQueue[id]) {
                                    inQueue[id] = 1;
                                    dirtyQueue.push(id);
                                }
                                return true;
                            });
                    });
            }
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "Local H-index r>=" << (int)r << " BK converged in " << iteration
              << " iterations, " << recomputeCount << " r-clique re-evaluations" << std::endl;
    std::cout << "time: " << elapsed << " ms" << std::endl;

    // Build output
    std::vector<std::pair<std::vector<daf::Size>, int>> result;
    result.reserve(numCliques);
    for (daf::Size i = 0; i < numCliques; ++i) {
        auto verts = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueCopy(verts.begin(), verts.end());
        result.emplace_back(std::move(cliqueCopy), (int)coreC[i]);
    }

    delete[] coreC;
    return result;
}
