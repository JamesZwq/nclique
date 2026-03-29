//
// Local H-index for r≥3 — Variant A: Vertex-Proxy CPI
//
// Tree is IMMUTABLE. No peeling, no BK, no tree mutation.
// Each r-clique computes its H-index using CPI breakpoint counting
// with vertex-level proxy cores (same trick as r=1 Local H-index).
//
// vertexProxy[v] = max(coreC[id] for all r-cliques containing v)
// For each leaf containing r-clique C, classify "other" vertices
// and use nCr breakpoints over their proxy cores.
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

namespace RCliqueLocalVP {

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

} // namespace RCliqueLocalVP


// ============================================================
// Compute H-index for r-clique C using vertex-proxy CPI breakpoints.
// Same breakpoint trick as r=1, but applied to "other" vertices' proxy cores.
// ============================================================
static long long computeHIndexVP(
    daf::Size cliqueId,
    const StaticCliqueIndex &cliqueIndex,
    const DynamicGraph<TreeGraphNode> &tree,
    const long long *coreC,
    const long long *vertexProxy,
    const std::vector<std::vector<daf::Size>> &cliqueLeafIds,
    daf::CliqueSize r, daf::CliqueSize s,
    const std::vector<int> &leafNeedPivot,
    const std::vector<uint8_t> &leafValid,
    long long currentCore,
    std::vector<long long> &pivotProxies,   // scratch
    std::vector<long long> &buckets)        // scratch
{
    if (currentCore <= 0) return 0;

    int bucketSize = (int)currentCore;
    if (bucketSize < 1) bucketSize = 1;
    if ((int)buckets.size() < bucketSize + 1)
        buckets.resize(bucketSize + 1, 0);
    for (int i = 0; i <= bucketSize; ++i) buckets[i] = 0;

    auto cliqueVerts = cliqueIndex.byId(cliqueId);

    // Build fast lookup for C's vertices
    // Use mapRef for O(1) membership test
    daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
    static constexpr daf::Size SENTINEL = std::numeric_limits<daf::Size>::max();
    // We mark C's vertices with a special tag
    for (auto v : cliqueVerts) mapRef[v] = 0; // 0 = "in C"

    long long rawTotalSupport = 0;

    if (cliqueId >= cliqueLeafIds.size()) {
        for (auto v : cliqueVerts) mapRef[v] = SENTINEL;
        return 0;
    }

    for (auto leafId : cliqueLeafIds[cliqueId]) {
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.empty()) continue;
        if (!leafValid[leafId]) continue;

        int needPivot = leafNeedPivot[leafId];

        // Count C's pivots and classify other vertices
        int cPivotsInLeaf = 0;
        long long minOtherKeepProxy = std::numeric_limits<long long>::max();
        bool anyKeepDead = false;
        pivotProxies.clear();

        for (const auto &node : leaf) {
            if (mapRef[node.v] == 0) {
                // This vertex is in C
                if (node.isPivot) cPivotsInLeaf++;
            } else {
                // Other vertex
                if (node.isPivot) {
                    pivotProxies.push_back(vertexProxy[node.v]);
                } else {
                    long long proxy = vertexProxy[node.v];
                    if (proxy < 1) { anyKeepDead = true; break; }
                    if (proxy < minOtherKeepProxy) minOtherKeepProxy = proxy;
                }
            }
        }
        if (anyKeepDead) continue;

        int effectiveNeedPivot = needPivot - cPivotsInLeaf;
        int numOtherPivots = (int)pivotProxies.size();
        if (effectiveNeedPivot < 0 || effectiveNeedPivot > numOtherPivots) continue;

        if (effectiveNeedPivot == 0) {
            // Single breakpoint at level = minOtherKeepProxy
            long long level = minOtherKeepProxy;
            if (level > bucketSize) level = bucketSize;
            if (level >= 1) {
                buckets[level] += 1;
                rawTotalSupport += 1;
            }
            continue;
        }

        // Sort descending for breakpoint generation
        std::sort(pivotProxies.begin(), pivotProxies.end(), std::greater<long long>());

        long long prevSupport = 0;
        int idx = 0;
        while (idx < numOtherPivots) {
            long long threshold = pivotProxies[idx];
            if (threshold < 1) break;
            if (threshold > minOtherKeepProxy) threshold = minOtherKeepProxy;

            int countAtThreshold = idx + 1;
            while (idx + 1 < numOtherPivots && pivotProxies[idx + 1] >= threshold) {
                idx++;
                countAtThreshold = idx + 1;
            }

            if (countAtThreshold >= effectiveNeedPivot) {
                long long support = (double)(nCr[countAtThreshold][effectiveNeedPivot] + 0.5);
                long long delta = support - prevSupport;
                if (delta > 0) {
                    int level = (threshold > bucketSize) ? bucketSize : (int)threshold;
                    if (level >= 1) {
                        buckets[level] += delta;
                        rawTotalSupport += delta;
                    }
                }
                prevSupport = support;
            }
            idx++;
        }
    }

    // Restore mapRef
    for (auto v : cliqueVerts) mapRef[v] = SENTINEL;

    if (rawTotalSupport < 1) return 0;

    // Bucket H-index scan (high to low)
    long long accumulated = 0;
    for (int c = bucketSize; c >= 1; --c) {
        accumulated += buckets[c];
        if (accumulated >= c) return c;
    }
    return 0;
}


// ============================================================
// Main entry: Local H-index with vertex-proxy CPI for r≥3
// ============================================================
std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLocalCPI_VP(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numLeaves = tree.adj_list.size();
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;

    // Build clique index
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build (LocalCPI_VP)", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });
    const daf::Size numCliques = cliqueIndex.size();

    // Build dual index
    auto dualIndex = daf::timeCount("buildDualIndex (LocalCPI_VP)", [&]() {
        return RCliqueLocalVP::buildDualIndex(tree, cliqueIndex, r, s);
    });
    auto &leafCliqueInfo = dualIndex.leafCliqueInfo;
    auto &cliqueLeafIds = dualIndex.cliqueLeafIds;

    // Initial support → core estimates
    auto *coreC = new long long[numCliques];
    for (daf::Size i = 0; i < numCliques; ++i)
        coreC[i] = dualIndex.counting[i];

    // Build vertex proxy: vertexProxy[v] = max(coreC[id] for all r-cliques containing v)
    auto *vertexProxy = new long long[numVertices + 1];
    memset(vertexProxy, 0, (numVertices + 1) * sizeof(long long));
    for (daf::Size cid = 0; cid < numCliques; ++cid) {
        auto verts = cliqueIndex.byId(cid);
        for (auto v : verts) {
            if (coreC[cid] > vertexProxy[v])
                vertexProxy[v] = coreC[cid];
        }
    }

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
    std::vector<long long> pivotProxies, buckets;
    pivotProxies.reserve(512);
    buckets.reserve(4096);

    // Dirty queue
    std::vector<uint8_t> inQueue(numCliques, 0);
    std::queue<daf::Size> dirtyQueue;

    std::cout << "=========================begin (Local CPI VP r>=" << (int)r << ")=========================" << std::endl;
    std::cout << "r-cliques: " << numCliques << ", leaves: " << numLeaves << ", vertices: " << numVertices << std::endl;

    auto initEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Init time: " << std::chrono::duration_cast<std::chrono::milliseconds>(initEnd - time_start).count() << " ms" << std::endl;

    // Phase 1: initial full scan
    for (daf::Size cid = 0; cid < numCliques; ++cid) {
        if (coreC[cid] <= 0) continue;

        long long newCore = computeHIndexVP(
            cid, cliqueIndex, tree, coreC, vertexProxy,
            cliqueLeafIds, r, s, leafNeedPivot, leafValid,
            coreC[cid], pivotProxies, buckets);
        newCore = std::min(newCore, coreC[cid]);

        if (newCore != coreC[cid]) {
            long long oldCore = coreC[cid];
            coreC[cid] = newCore;
            // Update vertexProxy for C's vertices
            auto verts = cliqueIndex.byId(cid);
            for (auto v : verts) {
                // Recompute proxy (conservative: just cap at newCore if it was the max)
                if (vertexProxy[v] == oldCore) {
                    // Need to recompute — find max among all r-cliques containing v
                    // For efficiency, just set to newCore (may be too low, but monotone decreasing is safe)
                    vertexProxy[v] = newCore;
                }
            }
            // Enqueue co-leaf r-cliques
            for (auto leafId : cliqueLeafIds[cid]) {
                if (!leafValid[leafId]) continue;
                for (const auto &entry : leafCliqueInfo[leafId]) {
                    if (entry.cliqueId == cid) continue;
                    if (coreC[entry.cliqueId] > 0 && !inQueue[entry.cliqueId]) {
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

            long long newCore = computeHIndexVP(
                cid, cliqueIndex, tree, coreC, vertexProxy,
                cliqueLeafIds, r, s, leafNeedPivot, leafValid,
                coreC[cid], pivotProxies, buckets);
            newCore = std::min(newCore, coreC[cid]);

            if (newCore != coreC[cid]) {
                long long oldCore = coreC[cid];
                coreC[cid] = newCore;
                // Update vertexProxy
                auto verts = cliqueIndex.byId(cid);
                for (auto v : verts) {
                    if (vertexProxy[v] == oldCore)
                        vertexProxy[v] = newCore;
                }
                // Enqueue co-leaf r-cliques
                for (auto leafId : cliqueLeafIds[cid]) {
                    if (!leafValid[leafId]) continue;
                    for (const auto &entry : leafCliqueInfo[leafId]) {
                        if (entry.cliqueId == cid) continue;
                        if (coreC[entry.cliqueId] > 0 && !inQueue[entry.cliqueId]) {
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
    std::cout << "Local CPI VP r>=" << (int)r << " converged in " << iteration
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
    delete[] vertexProxy;
    return result;
}
