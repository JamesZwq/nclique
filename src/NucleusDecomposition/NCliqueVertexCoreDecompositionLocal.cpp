//
// Local H-index iterative convergence for r=1 nucleus decomposition.
//
// Tree is immutable. Each vertex independently computes its core value
// via H-index over its neighboring leaves, iterating until convergence.
// Core values can only decrease (monotone), so convergence is guaranteed.
//
// Optimizations:
//   1. Merged single-pass leaf scanning (minKeepCore + pivotCores together)
//   2. Integer bucket H-index (replaces O(B log B) sort with O(B + maxLevel))
//   3. Early termination when H-index >= currentCore (core can't drop)
//   4. supportV-based smart enqueuing (only enqueue if support < core)
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <vector>
#include <cstring>
#include <queue>
#include <cmath>

#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace VCD_Local {

// Initial support counting (same as ST version)
double * countingPerVertex(const DynamicGraph<TreeGraphNode> &treeGraph,
                           const Graph &edgeGraph,
                           const daf::CliqueSize k) {
    const daf::Size n = edgeGraph.adj_list_offsets.size();
    auto *countingV = new double[n];
    memset(countingV, 0, n * sizeof(double));
    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keepC;
    for (const auto &clique : treeGraph.adj_list) {
        povit.clear();
        keepC.clear();
        if (clique.size() < k) continue;
        for (auto &i : clique) {
            if (i.isPivot) povit.push_back(i.v);
            else keepC.push_back(i.v);
        }
        int needPivot = int(k) - int(keepC.size());
        double keepVal = nCr[povit.size()][needPivot];
        for (const auto &v : keepC)
            countingV[v] += keepVal;
        double pivotVal = 0;
        const int npv = needPivot - 1;
        if (npv >= 0 && npv <= static_cast<int>(povit.size()) - 1)
            pivotVal = nCr[povit.size() - 1][npv];
        for (const auto &v : povit)
            countingV[v] += pivotVal;
    }
    povit.free();
    keepC.free();
    return countingV;
}

} // namespace VCD_Local


// ============================================================
// Optimized H-index computation with:
//   - Merged single-pass leaf scanning
//   - Integer bucket accumulation (no sort)
//   - Early termination via currentCore
// Returns the H-index (floor'd integer). Also writes totalSupport out.
// ============================================================
static double computeHIndex(
    daf::Size v,
    const DynamicGraph<TreeGraphNode> &tree,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    const double *coreV,
    const std::vector<int> &leafNeedPivot,
    const std::vector<uint8_t> &leafValid,
    double currentCore,                       // Opt 3: early termination
    double &totalSupportOut,                  // Opt 4: output total support
    std::vector<double> &pivotCores,          // scratch
    std::vector<double> &buckets)             // scratch: bucket[level] = support delta
{
    const auto &adjCliques = treeGraphV.getNbr(v);
    if (adjCliques.empty()) { totalSupportOut = 0.0; return 0.0; }

    // Determine max possible level for bucket sizing
    int maxLevel = 0;

    // We'll collect breakpoints into buckets. First pass: generate breakpoints,
    // track maxLevel. Use a two-phase approach: first collect into a temp vector,
    // then bucket them (since we don't know maxLevel upfront).
    // Actually, currentCore is an upper bound on the H-index, so maxLevel <= currentCore.
    int bucketSize = (int)currentCore + 1;
    if (bucketSize < 1) bucketSize = 1;

    // Ensure buckets vector is large enough and zeroed
    if ((int)buckets.size() < bucketSize + 1) {
        buckets.resize(bucketSize + 1, 0.0);
    }
    // Zero the range we'll use
    for (int i = 0; i <= bucketSize; ++i) buckets[i] = 0.0;

    double rawTotalSupport = 0.0;

    for (const auto &clique : adjCliques) {
        daf::Size leafId = clique.v;
        if (!leafValid[leafId]) continue;

        int needPivot = leafNeedPivot[leafId];

        // Opt 1: MERGED single pass — compute minKeepCore AND collect pivotCores
        double minKeepCore = 1e18;
        bool anyKeepDead = false;
        pivotCores.clear();

        for (auto &node : tree.adj_list[leafId]) {
            if (node.v == v) continue;
            if (node.isPivot) {
                pivotCores.push_back(coreV[node.v]);
            } else {
                double c = coreV[node.v];
                if (c < 1) { anyKeepDead = true; break; }
                if (c < minKeepCore) minKeepCore = c;
            }
        }
        if (anyKeepDead) continue;

        int effectiveNeedPivot = clique.isPivot ? needPivot - 1 : needPivot;
        int numOtherPivots = (int)pivotCores.size();
        if (effectiveNeedPivot < 0 || effectiveNeedPivot > numOtherPivots) continue;

        if (effectiveNeedPivot == 0) {
            // Single breakpoint at level = floor(minKeepCore), delta = 1
            int level = (int)std::floor(minKeepCore);
            if (level > bucketSize) level = bucketSize;
            if (level >= 1) {
                buckets[level] += 1.0;
            }
            rawTotalSupport += 1.0;
            continue;
        }

        // Sort descending for breakpoint generation
        std::sort(pivotCores.begin(), pivotCores.end(), std::greater<double>());

        double prevSupport = 0.0;
        int idx = 0;
        while (idx < numOtherPivots) {
            double threshold = pivotCores[idx];
            if (threshold < 1) break;
            if (threshold > minKeepCore) threshold = minKeepCore;

            int countAtThreshold = idx + 1;
            while (idx + 1 < numOtherPivots && pivotCores[idx + 1] >= threshold) {
                idx++;
                countAtThreshold = idx + 1;
            }

            double support = nCr[countAtThreshold][effectiveNeedPivot];
            double delta = support - prevSupport;
            if (delta > 0) {
                int level = (int)std::floor(threshold);
                if (level > bucketSize) level = bucketSize;
                if (level >= 1) {
                    buckets[level] += delta;
                }
                rawTotalSupport += delta;
            }
            prevSupport = support;
            idx++;
        }
    }

    totalSupportOut = rawTotalSupport;

    if (rawTotalSupport < 1.0) return 0.0;

    // Opt 2: Integer bucket H-index scan (from high to low)
    // H-index = max c such that accumulated support from levels >= c is >= c
    double accumulated = 0.0;
    // Clamp scan to bucketSize
    for (int c = bucketSize; c >= 1; --c) {
        accumulated += buckets[c];
        if (accumulated >= (double)c) {
            return (double)c;
        }
    }

    return 0.0;
}


// ============================================================
// Naive version: full scan every iteration (no dirty queue).
// Kept for experiments / baseline comparison.
// ============================================================
double * NCliqueVertexCoreDecomposition_LocalNaive(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numLeaves = tree.adj_list.size();
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;

    std::vector<int> leafNeedPivot(numLeaves);
    std::vector<uint8_t> leafValid(numLeaves);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        int keeps = 0, pivots = 0;
        for (auto &node : tree.adj_list[L]) {
            if (node.isPivot) pivots++;
            else keeps++;
        }
        leafNeedPivot[L] = (int)k - keeps;
        leafValid[L] = (leafNeedPivot[L] >= 0 && leafNeedPivot[L] <= pivots) ? 1 : 0;
    }

    auto initSupport = VCD_Local::countingPerVertex(tree, edgeGraph, k);
    auto *coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = 0.0;
    for (daf::Size i = 0; i < numVertices; ++i)
        coreV[i] = initSupport[i];
    delete[] initSupport;

    std::vector<double> pivotCores;
    std::vector<double> buckets;
    pivotCores.reserve(512);
    buckets.reserve(4096);

    std::cout << "=========================begin (Local H-index Naive)=========================" << std::endl;
    std::cout << "vertices: " << numVertices << ", leaves: " << numLeaves << std::endl;

    int iteration = 0;
    bool changed = true;

    while (changed) {
        changed = false;
        iteration++;

        for (daf::Size v = 0; v < numVertices; ++v) {
            if (coreV[v] <= 0) continue;

            double totalSup = 0.0;
            double newCore = computeHIndex(v, tree, treeGraphV, coreV,
                                            leafNeedPivot, leafValid,
                                            coreV[v], totalSup,
                                            pivotCores, buckets);
            newCore = std::min(newCore, coreV[v]);

            if (newCore != coreV[v]) {
                coreV[v] = newCore;
                changed = true;
            }
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "Local H-index Naive converged in " << iteration << " iterations" << std::endl;
    std::cout << "time: " << elapsed << " ms" << std::endl;

    for (daf::Size i = 0; i < numVertices; ++i) {
        if (coreV[i] <= 0) coreV[i] = -1.0;
    }

    return coreV;
}

// ============================================================
// Optimized version: dirty-queue propagation.
// Opts 1-3: merged leaf pass, bucket H-index, early termination.
// ============================================================
double * NCliqueVertexCoreDecomposition_Local(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numLeaves = tree.adj_list.size();
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;

    // --- Precompute per-leaf metadata (immutable) ---
    std::vector<int> leafNeedPivot(numLeaves);
    std::vector<uint8_t> leafValid(numLeaves);

    for (daf::Size L = 0; L < numLeaves; ++L) {
        int keeps = 0, pivots = 0;
        for (auto &node : tree.adj_list[L]) {
            if (node.isPivot) pivots++;
            else keeps++;
        }
        leafNeedPivot[L] = (int)k - keeps;
        leafValid[L] = (leafNeedPivot[L] >= 0 && leafNeedPivot[L] <= pivots) ? 1 : 0;
    }

    // --- Initial per-vertex core = initial support ---
    auto initSupport = VCD_Local::countingPerVertex(tree, edgeGraph, k);
    auto *coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = 0.0;
    for (daf::Size i = 0; i < numVertices; ++i)
        coreV[i] = initSupport[i];
    delete[] initSupport;

    // Scratch buffers
    std::vector<double> pivotCores;
    std::vector<double> buckets;
    pivotCores.reserve(512);
    buckets.reserve(4096);

    // --- Dirty queue: vertices that need re-evaluation ---
    std::vector<uint8_t> inQueue(numVertices, 0);
    std::queue<daf::Size> dirtyQueue;

    std::cout << "=========================begin (Local H-index)=========================" << std::endl;
    std::cout << "vertices: " << numVertices << ", leaves: " << numLeaves << std::endl;

    // Phase 1: initial full scan to get first H-index values
    for (daf::Size v = 0; v < numVertices; ++v) {
        if (coreV[v] <= 0) continue;

        double totalSup = 0.0;
        double newCore = computeHIndex(v, tree, treeGraphV, coreV,
                                        leafNeedPivot, leafValid,
                                        coreV[v], totalSup,
                                        pivotCores, buckets);
        newCore = std::min(newCore, coreV[v]); // monotone
        if (newCore != coreV[v]) {
            coreV[v] = newCore;
            // Enqueue neighbors that share a leaf with v
            auto &adjCliques = treeGraphV.getNbr(v);
            for (const auto &clique : adjCliques) {
                daf::Size leafId = clique.v;
                if (!leafValid[leafId]) continue;
                for (auto &node : tree.adj_list[leafId]) {
                    daf::Size u = node.v;
                    if (u == v) continue;
                    if (coreV[u] > 0 && !inQueue[u]) {
                        inQueue[u] = 1;
                        dirtyQueue.push(u);
                    }
                }
            }
        }
    }

    // Phase 2: process dirty queue until empty
    long long recomputeCount = 0;
    int iteration = 1;

    while (!dirtyQueue.empty()) {
        iteration++;
        size_t roundSize = dirtyQueue.size();
        for (size_t qi = 0; qi < roundSize; ++qi) {
            daf::Size v = dirtyQueue.front();
            dirtyQueue.pop();
            inQueue[v] = 0;
            recomputeCount++;

            if (coreV[v] <= 0) continue;

            double totalSup = 0.0;
            double newCore = computeHIndex(v, tree, treeGraphV, coreV,
                                            leafNeedPivot, leafValid,
                                            coreV[v], totalSup,
                                            pivotCores, buckets);
            newCore = std::min(newCore, coreV[v]); // monotone

            if (newCore != coreV[v]) {
                coreV[v] = newCore;
                // Propagate: enqueue co-leaf neighbors
                auto &adjCliques = treeGraphV.getNbr(v);
                for (const auto &clique : adjCliques) {
                    daf::Size leafId = clique.v;
                    if (!leafValid[leafId]) continue;
                    for (auto &node : tree.adj_list[leafId]) {
                        daf::Size u = node.v;
                        if (u == v) continue;
                        if (coreV[u] > 0 && !inQueue[u]) {
                            inQueue[u] = 1;
                            dirtyQueue.push(u);
                        }
                    }
                }
            }
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "Local H-index converged in " << iteration << " iterations, "
              << recomputeCount << " vertex re-evaluations" << std::endl;
    std::cout << "time: " << elapsed << " ms" << std::endl;

    // Set vertices with core=0 to -1 for consistency with ST output
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (coreV[i] <= 0) coreV[i] = -1.0;
    }

    return coreV;
}
