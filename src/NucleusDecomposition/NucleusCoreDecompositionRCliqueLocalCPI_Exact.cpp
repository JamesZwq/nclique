//
// Local H-index for r≥3 — Variant B: Exact CPI
//
// Tree is IMMUTABLE. No peeling, no BK, no tree mutation.
// Each r-clique computes its EXACT H-index by enumerating
// (s-r)-subsets of "other" vertices per leaf (CPI-guided),
// then finding min core among other r-subcliques in each s-clique.
//
// This is the same algorithm as RCliqueLocal but with two CPI optimizations:
//   1. Uses cliqueLeafIds[] instead of intersect_dense_sets_multi
//   2. Uses leafCliqueInfo[] for dirty-queue propagation instead of
//      enumerating all C(n,r) per leaf
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

namespace RCliqueLocalExact {

struct LeafCliqueEntry {
    daf::Size cliqueId;
    double ncrValue;
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

} // namespace RCliqueLocalExact


// ============================================================
// Compute EXACT H-index for r-clique C.
// For each leaf containing C, enumerate s-cliques via CPI
// (C's vertices + all otherKeeps + choose effectiveNeedPivot from otherPivots).
// For each s-clique, enumerate all C(s,r) r-subcliques, find min other core.
// ============================================================
static long long computeHIndexExact(
    daf::Size cliqueId,
    const StaticCliqueIndex &cliqueIndex,
    const DynamicGraph<TreeGraphNode> &tree,
    const long long *coreC,
    const std::vector<std::vector<daf::Size>> &cliqueLeafIds,
    daf::CliqueSize r, daf::CliqueSize s,
    const std::vector<int> &leafNeedPivot,
    const std::vector<uint8_t> &leafValid,
    long long currentCore,
    std::vector<long long> &buckets,
    std::vector<TreeGraphNode> &otherKeeps,
    std::vector<TreeGraphNode> &otherPivots)
{
    if (currentCore <= 0) return 0;

    int bucketSize = (int)currentCore;
    if (bucketSize < 1) bucketSize = 1;
    if ((int)buckets.size() < bucketSize + 1)
        buckets.resize(bucketSize + 1, 0);
    for (int i = 0; i <= bucketSize; ++i) buckets[i] = 0;

    auto cliqueVerts = cliqueIndex.byId(cliqueId);
    const daf::Size numCliques = cliqueIndex.size();

    // Build fast lookup for C's vertices
    robin_hood::unordered_flat_set<daf::Size> cliqueSet;
    for (auto v : cliqueVerts) cliqueSet.insert(v);

    long long rawTotalSupport = 0;

    if (cliqueId >= cliqueLeafIds.size()) return 0;

    for (auto leafId : cliqueLeafIds[cliqueId]) {
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.empty()) continue;
        if (!leafValid[leafId]) continue;

        int needPivot = leafNeedPivot[leafId];

        // Count C's pivots and classify other vertices
        int cPivotsInLeaf = 0;
        for (const auto &node : leaf) {
            if (cliqueSet.count(node.v) && node.isPivot) cPivotsInLeaf++;
        }

        otherKeeps.clear();
        otherPivots.clear();
        for (const auto &node : leaf) {
            if (cliqueSet.count(node.v)) continue;
            if (node.isPivot) otherPivots.push_back(node);
            else otherKeeps.push_back(node);
        }

        int effectiveNeedPivot = needPivot - cPivotsInLeaf;
        int numOtherKeeps = (int)otherKeeps.size();
        int numOtherPivots = (int)otherPivots.size();
        int need = (int)s - (int)r; // how many "other" vertices needed

        if (effectiveNeedPivot < 0 || effectiveNeedPivot > numOtherPivots) continue;
        if (numOtherKeeps + numOtherPivots < need) continue;
        if (numOtherKeeps > need) continue; // more keeps than slots

        // Process one s-clique given chosen pivots
        auto processSClique = [&](const TreeGraphNode *chosenPivotData, int chosenPivotCount) {
            // Build full s-clique
            daf::StaticVector<TreeGraphNode> fullSClique;
            for (const auto &node : leaf) {
                if (cliqueSet.count(node.v)) fullSClique.push_back(node);
            }
            for (const auto &kv : otherKeeps) fullSClique.push_back(kv);
            for (int pi = 0; pi < chosenPivotCount; ++pi) fullSClique.push_back(chosenPivotData[pi]);
            std::sort(fullSClique.begin(), fullSClique.end());

            long long minOtherCore = std::numeric_limits<long long>::max();
            daf::enumerateCombinations(fullSClique, r,
                [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                    auto id = cliqueIndex.tryByClique(rClique);
                    if (id >= numCliques) return true;
                    if (id == cliqueId) return true;
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
    }

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
// Main entry: Exact CPI Local H-index for r≥3
// ============================================================
std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLocalCPI_Exact(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numLeaves = tree.adj_list.size();

    // Build clique index
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build (LocalCPI_Exact)", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });
    const daf::Size numCliques = cliqueIndex.size();

    // Build dual index
    auto dualIndex = daf::timeCount("buildDualIndex (LocalCPI_Exact)", [&]() {
        return RCliqueLocalExact::buildDualIndex(tree, cliqueIndex, r, s);
    });
    auto &leafCliqueInfo = dualIndex.leafCliqueInfo;
    auto &cliqueLeafIds = dualIndex.cliqueLeafIds;

    // Initial support → core estimates
    auto *coreC = new long long[numCliques];
    for (daf::Size i = 0; i < numCliques; ++i)
        coreC[i] = dualIndex.counting[i];

    // Precompute per-leaf metadata (immutable)
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
    std::vector<TreeGraphNode> otherKeeps, otherPivots;
    buckets.reserve(4096);
    otherKeeps.reserve(64);
    otherPivots.reserve(64);

    // Dirty queue
    std::vector<uint8_t> inQueue(numCliques, 0);
    std::queue<daf::Size> dirtyQueue;

    std::cout << "=========================begin (Local CPI Exact r>=" << (int)r << ")=========================" << std::endl;
    std::cout << "r-cliques: " << numCliques << ", leaves: " << numLeaves << std::endl;

    auto initEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Init time: " << std::chrono::duration_cast<std::chrono::milliseconds>(initEnd - time_start).count() << " ms" << std::endl;

    // Phase 1: initial full scan
    for (daf::Size cid = 0; cid < numCliques; ++cid) {
        if (coreC[cid] <= 0) continue;

        long long newCore = computeHIndexExact(
            cid, cliqueIndex, tree, coreC,
            cliqueLeafIds, r, s, leafNeedPivot, leafValid,
            coreC[cid], buckets, otherKeeps, otherPivots);
        newCore = std::min(newCore, coreC[cid]);

        if (newCore != coreC[cid]) {
            long long oldCore = coreC[cid];
            coreC[cid] = newCore;
            // Enqueue co-leaf r-cliques via leafCliqueInfo (CPI optimization)
            for (auto leafId : cliqueLeafIds[cid]) {
                if (!leafValid[leafId]) continue;
                for (const auto &entry : leafCliqueInfo[leafId]) {
                    if (entry.cliqueId == cid) continue;
                    if (coreC[entry.cliqueId] > 0 && oldCore >= coreC[entry.cliqueId] && !inQueue[entry.cliqueId]) {
                        inQueue[entry.cliqueId] = 1;
                        dirtyQueue.push(entry.cliqueId);
                    }
                }
            }
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

            long long newCore = computeHIndexExact(
                cid, cliqueIndex, tree, coreC,
                cliqueLeafIds, r, s, leafNeedPivot, leafValid,
                coreC[cid], buckets, otherKeeps, otherPivots);
            newCore = std::min(newCore, coreC[cid]);

            if (newCore != coreC[cid]) {
                long long oldCore = coreC[cid];
                coreC[cid] = newCore;
                for (auto leafId : cliqueLeafIds[cid]) {
                    if (!leafValid[leafId]) continue;
                    for (const auto &entry : leafCliqueInfo[leafId]) {
                        if (entry.cliqueId == cid) continue;
                        if (coreC[entry.cliqueId] > 0 && oldCore >= coreC[entry.cliqueId] && !inQueue[entry.cliqueId]) {
                            inQueue[entry.cliqueId] = 1;
                            dirtyQueue.push(entry.cliqueId);
                        }
                    }
                }
            }
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "Local CPI Exact r>=" << (int)r << " converged in " << iteration
              << " iterations, " << recomputeCount << " r-clique re-evaluations" << std::endl;
    std::cout << "time: " << elapsed << " ms" << std::endl;

    // Build output
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    result.reserve(numCliques);
    for (daf::Size i = 0; i < numCliques; ++i) {
        auto verts = cliqueIndex.byId(i);
        std::vector<daf::Size> copy(verts.begin(), verts.end());
        result.emplace_back(std::move(copy), (int)coreC[i]);
    }

    delete[] coreC;
    return result;
}
