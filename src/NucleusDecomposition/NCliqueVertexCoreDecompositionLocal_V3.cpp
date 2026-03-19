//
// Local H-index V3 — parallel convergence for r=1 nucleus decomposition.
//
// Parallelization strategy: round-based synchronous.
//   - Each round: all queued vertices compute H-index in parallel (read frozen coreV)
//   - After round: apply updates sequentially, collect dirty neighbors for next round
//   - Phase 1 (initial scan) fully parallelized over all vertices
//   - Per-thread scratch buffers (pivotCores, buckets) — no contention
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <vector>
#include <cstring>
#include <omp.h>
#include <atomic>

#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace VCD_LocalV3 {

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

} // namespace VCD_LocalV3


// H-index: merged pass + bucket accumulation. Thread-safe (uses caller's scratch).
static double computeHIndexV3(
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


// Collect all co-leaf neighbors of changed vertices into nextWork.
// Uses atomic CAS on inQueue for dedup.
static void collectDirtyNeighbors(
    const std::vector<daf::Size> &changedVertices,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    const DynamicGraph<TreeGraphNode> &tree,
    const std::vector<uint8_t> &leafValid,
    const double *coreV,
    std::atomic<uint8_t> *inQueue,
    std::vector<daf::Size> &nextWork)
{
    for (auto v : changedVertices) {
        auto &adjCliques = treeGraphV.getNbr(v);
        for (const auto &clique : adjCliques) {
            daf::Size leafId = clique.v;
            if (!leafValid[leafId]) continue;
            for (auto &node : tree.adj_list[leafId]) {
                daf::Size u = node.v;
                if (u == v) continue;
                if (coreV[u] > 0) {
                    uint8_t expected = 0;
                    if (inQueue[u].compare_exchange_strong(expected, 1,
                            std::memory_order_relaxed)) {
                        nextWork.push_back(u);
                    }
                }
            }
        }
    }
}


// ============================================================
// Local H-index V3: OpenMP parallel round-based convergence.
// ============================================================
double * NCliqueVertexCoreDecomposition_LocalV3(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const int nThreads = omp_get_max_threads();
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

    auto initSupport = VCD_LocalV3::countingPerVertex(tree, edgeGraph, k);
    auto *coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = 0.0;
    for (daf::Size i = 0; i < numVertices; ++i)
        coreV[i] = initSupport[i];
    delete[] initSupport;

    // Per-thread scratch buffers
    std::vector<std::vector<double>> threadPivotCores(nThreads);
    std::vector<std::vector<double>> threadBuckets(nThreads);
    for (int t = 0; t < nThreads; ++t) {
        threadPivotCores[t].reserve(512);
        threadBuckets[t].reserve(4096);
    }

    // newCoreV: parallel round writes here, applied to coreV after barrier
    auto *newCoreV = new double[numVertices];

    auto *inQueue = new std::atomic<uint8_t>[numVertices];
    for (daf::Size i = 0; i < numVertices; ++i)
        inQueue[i].store(0, std::memory_order_relaxed);

    std::vector<daf::Size> currentWork;
    std::vector<daf::Size> nextWork;
    std::vector<daf::Size> changedVertices;
    currentWork.reserve(numVertices);
    nextWork.reserve(numVertices);
    changedVertices.reserve(numVertices);

    std::cout << "=========================begin (Local H-index V3 parallel)=========================" << std::endl;
    std::cout << "vertices: " << numVertices << ", leaves: " << numLeaves
              << ", threads: " << nThreads << std::endl;

    // ======== Phase 1: parallel initial full scan ========
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto &pc = threadPivotCores[tid];
        auto &bk = threadBuckets[tid];

        #pragma omp for schedule(dynamic, 256)
        for (daf::Size v = 0; v < numVertices; ++v) {
            if (coreV[v] <= 0) {
                newCoreV[v] = coreV[v];
                continue;
            }
            double h = computeHIndexV3(v, tree, treeGraphV, coreV,
                                        leafNeedPivot, leafValid,
                                        coreV[v], pc, bk);
            newCoreV[v] = std::min(h, coreV[v]);
        }
    }

    // Apply changes, track who changed
    changedVertices.clear();
    for (daf::Size v = 0; v < numVertices; ++v) {
        if (newCoreV[v] != coreV[v]) {
            coreV[v] = newCoreV[v];
            changedVertices.push_back(v);
        }
    }

    // Collect dirty neighbors of changed vertices
    collectDirtyNeighbors(changedVertices, treeGraphV, tree, leafValid,
                          coreV, inQueue, currentWork);

    // ======== Phase 2: parallel iterative rounds ========
    long long totalRecomputeCount = 0;
    int iteration = 1;

    while (!currentWork.empty()) {
        iteration++;
        totalRecomputeCount += (long long)currentWork.size();

        // Parallel H-index computation — all read frozen coreV[]
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto &pc = threadPivotCores[tid];
            auto &bk = threadBuckets[tid];

            #pragma omp for schedule(dynamic, 64)
            for (int wi = 0; wi < (int)currentWork.size(); ++wi) {
                daf::Size v = currentWork[wi];
                if (coreV[v] <= 0) {
                    newCoreV[v] = coreV[v];
                    continue;
                }
                double h = computeHIndexV3(v, tree, treeGraphV, coreV,
                                            leafNeedPivot, leafValid,
                                            coreV[v], pc, bk);
                newCoreV[v] = std::min(h, coreV[v]);
            }
        }

        // Clear inQueue for current batch
        for (auto v : currentWork)
            inQueue[v].store(0, std::memory_order_relaxed);

        // Apply updates, collect changed vertices
        changedVertices.clear();
        for (auto v : currentWork) {
            if (newCoreV[v] != coreV[v]) {
                coreV[v] = newCoreV[v];
                changedVertices.push_back(v);
            }
        }

        // Collect dirty neighbors for next round
        nextWork.clear();
        collectDirtyNeighbors(changedVertices, treeGraphV, tree, leafValid,
                              coreV, inQueue, nextWork);

        std::swap(currentWork, nextWork);
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "Local H-index V3 converged in " << iteration << " iterations, "
              << totalRecomputeCount << " re-evaluations" << std::endl;
    std::cout << "time: " << elapsed << " ms" << std::endl;

    for (daf::Size i = 0; i < numVertices; ++i) {
        if (coreV[i] <= 0) coreV[i] = -1.0;
    }

    delete[] newCoreV;
    delete[] inQueue;
    return coreV;
}
