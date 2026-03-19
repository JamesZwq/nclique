//
// Local H-index V4 — asynchronous parallel for r=1 nucleus decomposition.
//
// Key insight: coreV[] is monotone decreasing. Any thread reading a stale
// value sees an UPPER bound, which is safe — it may delay convergence by
// one round but never produces incorrect results.
//
// Design:
//   1. Async in-place updates: threads write coreV[] directly (no double buffer)
//      → eliminates the doubled re-evaluation cost of V3's synchronous rounds
//   2. Parallel neighbor enqueuing: per-thread local work lists, merged at round end
//   3. Core-level enqueue filter: only enqueue neighbor u if oldCore >= coreV[u]
//   4. Phase 1: parallel initial scan with immediate coreV[] update + parallel enqueue
//   5. Phase 2: parallel rounds processing work list, each thread updates coreV in place
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

namespace VCD_LocalV4 {

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

} // namespace VCD_LocalV4


// H-index: merged pass + bucket accumulation. Thread-safe (uses caller's scratch).
// Reads coreV[] which may be concurrently updated — safe because monotone decreasing.
static double computeHIndexV4(
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
// Local H-index V4: asynchronous parallel with in-place updates.
// ============================================================
double * NCliqueVertexCoreDecomposition_LocalV4(
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

    auto initSupport = VCD_LocalV4::countingPerVertex(tree, edgeGraph, k);
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

    // Per-thread local work lists for parallel enqueuing
    std::vector<std::vector<daf::Size>> threadLocalWork(nThreads);
    for (int t = 0; t < nThreads; ++t)
        threadLocalWork[t].reserve(numVertices / nThreads + 1024);

    // Atomic inQueue flag
    auto *inQueue = new std::atomic<uint8_t>[numVertices];
    for (daf::Size i = 0; i < numVertices; ++i)
        inQueue[i].store(0, std::memory_order_relaxed);

    std::vector<daf::Size> currentWork;
    currentWork.reserve(numVertices);

    std::cout << "=========================begin (Local H-index V4 async)=========================" << std::endl;
    std::cout << "vertices: " << numVertices << ", leaves: " << numLeaves
              << ", threads: " << nThreads << std::endl;

    // ======== Phase 1: parallel initial scan with in-place update ========
    // Each thread computes H-index and updates coreV[] immediately.
    // Other threads may read stale values — safe because monotone decreasing.
    // Per-thread local enqueue lists merged after.
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto &pc = threadPivotCores[tid];
        auto &bk = threadBuckets[tid];
        auto &localWork = threadLocalWork[tid];
        localWork.clear();

        #pragma omp for schedule(dynamic, 256)
        for (daf::Size v = 0; v < numVertices; ++v) {
            if (coreV[v] <= 0) continue;

            double oldCore = coreV[v];
            double h = computeHIndexV4(v, tree, treeGraphV, coreV,
                                        leafNeedPivot, leafValid,
                                        oldCore, pc, bk);
            double newCore = std::min(h, oldCore);

            if (newCore != oldCore) {
                coreV[v] = newCore;  // in-place update (monotone ↓)

                // Enqueue co-leaf neighbors (core-level filter)
                auto &adjCliques = treeGraphV.getNbr(v);
                for (const auto &clique : adjCliques) {
                    daf::Size leafId = clique.v;
                    if (!leafValid[leafId]) continue;
                    for (auto &node : tree.adj_list[leafId]) {
                        daf::Size u = node.v;
                        if (u == v) continue;
                        if (coreV[u] <= 0) continue;
                        if (oldCore >= coreV[u]) {
                            uint8_t expected = 0;
                            if (inQueue[u].compare_exchange_strong(expected, 1,
                                    std::memory_order_relaxed)) {
                                localWork.push_back(u);
                            }
                        }
                    }
                }
            }
        }
    }

    // Merge per-thread work lists
    for (int t = 0; t < nThreads; ++t) {
        currentWork.insert(currentWork.end(),
                           threadLocalWork[t].begin(), threadLocalWork[t].end());
    }

    // ======== Phase 2: async parallel iterative rounds ========
    long long totalRecomputeCount = 0;
    int iteration = 1;

    while (!currentWork.empty()) {
        iteration++;
        totalRecomputeCount += (long long)currentWork.size();

        // Parallel: compute H-index AND update coreV in place AND enqueue neighbors
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto &pc = threadPivotCores[tid];
            auto &bk = threadBuckets[tid];
            auto &localWork = threadLocalWork[tid];
            localWork.clear();

            #pragma omp for schedule(dynamic, 64)
            for (int wi = 0; wi < (int)currentWork.size(); ++wi) {
                daf::Size v = currentWork[wi];
                inQueue[v].store(0, std::memory_order_relaxed);

                double oldCore = coreV[v];
                if (oldCore <= 0) continue;

                double h = computeHIndexV4(v, tree, treeGraphV, coreV,
                                            leafNeedPivot, leafValid,
                                            oldCore, pc, bk);
                double newCore = std::min(h, oldCore);

                if (newCore != oldCore) {
                    coreV[v] = newCore;  // in-place update

                    // Enqueue neighbors with core-level filter
                    auto &adjCliques = treeGraphV.getNbr(v);
                    for (const auto &clique : adjCliques) {
                        daf::Size leafId = clique.v;
                        if (!leafValid[leafId]) continue;
                        for (auto &node : tree.adj_list[leafId]) {
                            daf::Size u = node.v;
                            if (u == v) continue;
                            if (coreV[u] <= 0) continue;
                            if (oldCore >= coreV[u]) {
                                uint8_t expected = 0;
                                if (inQueue[u].compare_exchange_strong(expected, 1,
                                        std::memory_order_relaxed)) {
                                    localWork.push_back(u);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Merge per-thread work lists for next round
        currentWork.clear();
        for (int t = 0; t < nThreads; ++t) {
            currentWork.insert(currentWork.end(),
                               threadLocalWork[t].begin(), threadLocalWork[t].end());
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "Local H-index V4 converged in " << iteration << " iterations, "
              << totalRecomputeCount << " re-evaluations" << std::endl;
    std::cout << "time: " << elapsed << " ms" << std::endl;

    for (daf::Size i = 0; i < numVertices; ++i) {
        if (coreV[i] <= 0) coreV[i] = -1.0;
    }

    delete[] inQueue;
    return coreV;
}
