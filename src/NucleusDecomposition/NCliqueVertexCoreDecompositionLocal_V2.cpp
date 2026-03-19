//
// Local H-index V2 — optimized convergence for r=1 nucleus decomposition.
//
// Improvements over V1:
//   1. Core-level enqueue filter: when v's core drops from oldCore to newCore,
//      only enqueue neighbor u if oldCore >= coreV[u]. If v was already below
//      u's core, v's drop can't affect u's H-index.
//   2. Leaf-level timestamp skip: track when each leaf last had a member change.
//      Skip re-evaluation if no leaf of v has changed since v's last eval.
//   3. FIFO queue + merged pass + bucket H-index (proven from V1).
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <vector>
#include <cstring>
#include <queue>

#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace VCD_LocalV2 {

static double * countingPerVertex(const DynamicGraph<TreeGraphNode> &treeGraph,
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

} // namespace VCD_LocalV2


// H-index: merged pass + bucket accumulation.
static double computeHIndexV2(
    daf::Size v,
    const DynamicGraph<TreeGraphNode> &tree,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    const double *coreV,
    const std::vector<int> &leafNeedPivot,
    const std::vector<uint8_t> &leafValid,
    double currentCore,
    std::vector<double> &pivotCores,
    std::vector<double> &buckets)
{
    const auto &adjCliques = treeGraphV.getNbr(v);
    if (adjCliques.empty()) return 0.0;

    int bucketSize = (int)currentCore + 1;
    if (bucketSize < 1) bucketSize = 1;

    if ((int)buckets.size() < bucketSize + 1)
        buckets.resize(bucketSize + 1, 0.0);
    for (int i = 0; i <= bucketSize; ++i) buckets[i] = 0.0;

    double rawTotalSupport = 0.0;

    for (const auto &clique : adjCliques) {
        daf::Size leafId = clique.v;
        if (!leafValid[leafId]) continue;

        int needPivot = leafNeedPivot[leafId];

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
            int level = (int)minKeepCore;
            if (level > bucketSize) level = bucketSize;
            if (level >= 1) buckets[level] += 1.0;
            rawTotalSupport += 1.0;
            continue;
        }

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
                int level = (int)threshold;
                if (level > bucketSize) level = bucketSize;
                if (level >= 1) buckets[level] += delta;
                rawTotalSupport += delta;
            }
            prevSupport = support;
            idx++;
        }
    }

    if (rawTotalSupport < 1.0) return 0.0;

    double accumulated = 0.0;
    for (int c = bucketSize; c >= 1; --c) {
        accumulated += buckets[c];
        if (accumulated >= (double)c)
            return (double)c;
    }
    return 0.0;
}


// ============================================================
// Local H-index V2: core-level enqueue filter + timestamp skip.
// ============================================================
double * NCliqueVertexCoreDecomposition_LocalV2(
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

    auto initSupport = VCD_LocalV2::countingPerVertex(tree, edgeGraph, k);
    auto *coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = 0.0;
    for (daf::Size i = 0; i < numVertices; ++i)
        coreV[i] = initSupport[i];
    delete[] initSupport;

    std::vector<double> pivotCores;
    std::vector<double> buckets;
    pivotCores.reserve(512);
    buckets.reserve(4096);

    // Leaf-level timestamp for skip optimization
    std::vector<uint32_t> leafTimestamp(numLeaves, 0);
    std::vector<uint32_t> vertexTimestamp(numVertices, 0);
    uint32_t globalClock = 0;

    std::vector<uint8_t> inQueue(numVertices, 0);
    std::queue<daf::Size> dirtyQueue;

    std::cout << "=========================begin (Local H-index V2)=========================" << std::endl;
    std::cout << "vertices: " << numVertices << ", leaves: " << numLeaves << std::endl;

    // Helper: when v's core drops from oldCore, enqueue neighbors with core-level filter
    auto propagate = [&](daf::Size v, double oldCore) {
        globalClock++;
        auto &adjCliques = treeGraphV.getNbr(v);
        for (const auto &clique : adjCliques) {
            daf::Size leafId = clique.v;
            if (!leafValid[leafId]) continue;
            leafTimestamp[leafId] = globalClock;
            for (auto &node : tree.adj_list[leafId]) {
                daf::Size u = node.v;
                if (u == v) continue;
                if (coreV[u] <= 0 || inQueue[u]) continue;
                // Core-level filter: only enqueue if v's OLD core was >= u's core.
                // If oldCore < coreV[u], then v was already below u's H-index threshold,
                // so v's further drop can't affect u.
                if (oldCore >= coreV[u]) {
                    inQueue[u] = 1;
                    dirtyQueue.push(u);
                }
            }
        }
    };

    // Phase 1: initial full scan
    for (daf::Size v = 0; v < numVertices; ++v) {
        if (coreV[v] <= 0) continue;

        double newCore = computeHIndexV2(v, tree, treeGraphV, coreV,
                                          leafNeedPivot, leafValid,
                                          coreV[v], pivotCores, buckets);
        newCore = std::min(newCore, coreV[v]);
        vertexTimestamp[v] = globalClock;
        if (newCore != coreV[v]) {
            double oldCore = coreV[v];
            coreV[v] = newCore;
            propagate(v, oldCore);
        }
    }

    // Phase 2: dirty queue
    long long recomputeCount = 0;
    long long skipCount = 0;
    long long enqueueFilteredCount = 0;

    while (!dirtyQueue.empty()) {
        daf::Size v = dirtyQueue.front();
        dirtyQueue.pop();
        inQueue[v] = 0;

        if (coreV[v] <= 0) continue;

        // Timestamp skip: if no leaf changed since v's last eval, skip
        bool anyLeafChanged = false;
        auto &adjCliques = treeGraphV.getNbr(v);
        for (const auto &clique : adjCliques) {
            if (leafValid[clique.v] && leafTimestamp[clique.v] > vertexTimestamp[v]) {
                anyLeafChanged = true;
                break;
            }
        }
        if (!anyLeafChanged) {
            skipCount++;
            continue;
        }

        recomputeCount++;

        double newCore = computeHIndexV2(v, tree, treeGraphV, coreV,
                                          leafNeedPivot, leafValid,
                                          coreV[v], pivotCores, buckets);
        newCore = std::min(newCore, coreV[v]);
        vertexTimestamp[v] = globalClock;

        if (newCore != coreV[v]) {
            double oldCore = coreV[v];
            coreV[v] = newCore;
            propagate(v, oldCore);
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "Local H-index V2 converged: "
              << recomputeCount << " re-evals, "
              << skipCount << " skipped" << std::endl;
    std::cout << "time: " << elapsed << " ms" << std::endl;

    for (daf::Size i = 0; i < numVertices; ++i) {
        if (coreV[i] <= 0) coreV[i] = -1.0;
    }

    return coreV;
}
