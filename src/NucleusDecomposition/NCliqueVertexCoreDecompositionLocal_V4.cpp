//
// Local H-index V4 — asynchronous parallel for r=1 nucleus decomposition.
//
// Reads paths from the same flat dual CSR (ST_V2_Data) used by SPIN★, so
// the SPIN vs SPIN★ memory comparison reflects algorithmic differences
// (per-event bookkeeping, peel-order state) rather than data-structure
// layout overhead. Earlier versions used the legacy DynamicGraph /
// DynamicGraphSet vector-of-vectors / set-per-vertex structures inherited
// from the reference baseline; on dense graphs those carried tens of
// GB of STL container overhead unrelated to the algorithm.
//
// Key insight: coreV[] is monotone decreasing. Any thread reading a stale
// value sees an UPPER bound, which is safe — it may delay convergence by
// one round but never produces incorrect results.
//
// Design:
//   1. Async in-place updates: threads write coreV[] directly (no double buffer)
//   2. Parallel neighbor enqueuing: per-thread local work lists, merged at round end
//   3. Core-level enqueue filter: only enqueue neighbor u if oldCore >= coreV[u]
//   4. Phase 1: parallel initial scan with immediate coreV[] update + parallel enqueue
//   5. Phase 2: parallel rounds processing work list, each thread updates coreV in place
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include "../PhaseLogger.h"
#include <chrono>
#include <algorithm>
#include <vector>
#include <map>
#include <cstring>
#include <omp.h>
#include <atomic>

extern double nCr[1001][401];


// H-index: merged pass + bucket accumulation. Thread-safe (uses caller's scratch).
// Reads coreV[] which may be concurrently updated — safe because monotone decreasing.
//
// int64 LEVELS: legacy code cast `(int)currentCore + 1` for bucketSize, which
// silently wrapped to negative when currentCore = sdeg(v) exceeded INT_MAX
// on dense graphs at moderate s (com-dblp s=8 → sdeg ~ 1e14). The wrap
// poisoned the iterative h-index, causing infinite-loop non-convergence at
// high s. We use int64_t throughout for level / bucketSize so the cap is
// effectively unbounded.
//
// SPARSE buckets (std::map<int64_t,double>): distinct contribution levels
// per call are bounded by O(d_Σ(v)) (paths through v), typically <<
// bucketSize. The sparse map caps memory at the actual number of touched
// levels — this avoids the multi-GB allocations of the legacy dense vector.
//
// CSR PATH ACCESS: paths are read directly from ST_V2_Data's flat dual CSR
// (vtxLeafIds / leafVtxIds + bit-packed isPivot). No vector-of-vectors,
// no per-vertex std::set — same compact layout SPIN★ uses.
static double computeHIndexV4(
    daf::Size v,
    const ST_V2_Data &data,
    const double *coreV,
    const std::vector<int> &leafNeedPivot,
    const std::vector<uint8_t> &leafValid,
    double currentCore,
    std::vector<double> &pivotCores,
    std::map<int64_t, double> &buckets)
{
    const daf::Size vBeg = data.vtxLeafOff[v];
    const daf::Size vEnd = data.vtxLeafOff[v + 1];
    if (vBeg == vEnd) return 0.0;

    // bucketSize cap on level (= h[v] + 1). int64 avoids the legacy overflow.
    // Keep a finite ceiling so casts stay safe; INT64_MAX/2 is plenty.
    constexpr int64_t kLevelCap = (int64_t)1 << 60;
    int64_t bucketSize = (currentCore > (double)kLevelCap)
                       ? kLevelCap : (int64_t)currentCore + 1;
    if (bucketSize < 1) bucketSize = 1;

    buckets.clear();

    double rawTotalSupport = 0.0;

    for (daf::Size ei = vBeg; ei < vEnd; ++ei) {
        daf::Size leafId = data.vtxLeafIds[ei];
        if (!leafValid[leafId]) continue;
        const bool v_isPivot = STV3_getBit(data.vtxLeafIsPivot, ei);

        int needPivot = leafNeedPivot[leafId];

        double minKeepCore = 1e18;
        bool anyKeepDead = false;
        pivotCores.clear();

        const daf::Size lBeg = data.leafVtxOff[leafId];
        const daf::Size lEnd = data.leafVtxOff[leafId + 1];
        for (daf::Size lj = lBeg; lj < lEnd; ++lj) {
            daf::Size u = data.leafVtxIds[lj];
            if (u == v) continue;
            if (STV3_getBit(data.leafVtxIsPivot, lj)) {
                pivotCores.push_back(coreV[u]);
            } else {
                double c = coreV[u];
                if (c < 1) { anyKeepDead = true; break; }
                if (c < minKeepCore) minKeepCore = c;
            }
        }
        if (anyKeepDead) continue;

        int effectiveNeedPivot = v_isPivot ? needPivot - 1 : needPivot;
        int numOtherPivots = (int)pivotCores.size();
        if (effectiveNeedPivot < 0 || effectiveNeedPivot > numOtherPivots) continue;

        if (effectiveNeedPivot == 0) {
            int64_t level = (int64_t)minKeepCore;
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
                int64_t level = (int64_t)threshold;
                if (level > bucketSize) level = bucketSize;
                if (level >= 1) buckets[level] += delta;
                rawTotalSupport += delta;
            }
            prevSupport = support;
            idx++;
        }
    }

    if (rawTotalSupport < 1.0) return 0.0;

    // SPARSE descending accumulation. Walk non-empty levels high→low; between
    // two adjacent non-empty levels, the cumulative is constant — so an
    // intermediate (empty) level k satisfies prev_acc ≥ k iff
    // k ≤ floor(prev_acc). The "gap check" before processing each level
    // catches answers in those gaps.
    double prev_acc = 0.0;
    int64_t prev_level = bucketSize + 1;
    for (auto it = buckets.rbegin(); it != buckets.rend(); ++it) {
        int64_t c = it->first;
        // Gap (c, prev_level - 1]: cumulative was prev_acc.
        if (prev_acc >= 1.0) {
            int64_t gap_ans = (prev_acc >= (double)(prev_level - 1))
                            ? (prev_level - 1) : (int64_t)prev_acc;
            if (gap_ans > c && gap_ans >= 1) return (double)gap_ans;
        }
        double new_acc = prev_acc + it->second;
        if (new_acc >= (double)c) return (double)c;
        prev_acc = new_acc;
        prev_level = c;
    }
    // Tail gap [1, prev_level - 1].
    if (prev_acc >= 1.0) {
        int64_t gap_ans = (prev_acc >= (double)(prev_level - 1))
                        ? (prev_level - 1) : (int64_t)prev_acc;
        if (gap_ans >= 1) return (double)gap_ans;
    }
    return 0.0;
}


// ============================================================
// Local H-index V4: asynchronous parallel with in-place updates.
// Reads paths from ST_V2_Data dual CSR (same layout as SPIN★).
// ============================================================
double * NCliqueVertexCoreDecomposition_LocalV4_Peel(
    ST_V2_Data &data, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const int nThreads = omp_get_max_threads();
    const daf::Size numLeaves = data.numLeaves;
    const daf::Size numVertices = data.numVertices;

    // Compute per-leaf needPivot and validity from CSR (single pass over
    // leaf vertex sets to count pivots vs keeps).
    std::vector<int> leafNeedPivot(numLeaves);
    std::vector<uint8_t> leafValid(numLeaves);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        int keeps = 0, pivots = 0;
        const daf::Size lBeg = data.leafVtxOff[L];
        const daf::Size lEnd = data.leafVtxOff[L + 1];
        for (daf::Size lj = lBeg; lj < lEnd; ++lj) {
            if (STV3_getBit(data.leafVtxIsPivot, lj)) pivots++;
            else keeps++;
        }
        leafNeedPivot[L] = (int)k - keeps;
        leafValid[L] = (leafNeedPivot[L] >= 0 && leafNeedPivot[L] <= pivots) ? 1 : 0;
    }

    // Initial support comes from data.countingV (already populated by Build).
    auto *coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = 0.0;
    if (data.countingV != nullptr) {
        for (daf::Size i = 0; i < numVertices; ++i)
            coreV[i] = data.countingV[i];
        delete[] data.countingV;
        data.countingV = nullptr;
    }

    // Per-thread scratch buffers. SPARSE buckets: std::map<int64_t,double>
    // instead of dense std::vector<double>(sdeg+1). int64 levels avoid the
    // legacy int-overflow that poisoned the iterative h-index at high s.
    std::vector<std::vector<double>> threadPivotCores(nThreads);
    std::vector<std::map<int64_t, double>> threadBuckets(nThreads);
    for (int t = 0; t < nThreads; ++t) {
        threadPivotCores[t].reserve(512);
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
            double h = computeHIndexV4(v, data, coreV,
                                        leafNeedPivot, leafValid,
                                        oldCore, pc, bk);
            double newCore = std::min(h, oldCore);

            if (newCore != oldCore) {
                coreV[v] = newCore;  // in-place update (monotone ↓)

                // Enqueue co-leaf neighbors (core-level filter) via dual CSR.
                const daf::Size vBeg = data.vtxLeafOff[v];
                const daf::Size vEnd = data.vtxLeafOff[v + 1];
                for (daf::Size ei = vBeg; ei < vEnd; ++ei) {
                    daf::Size leafId = data.vtxLeafIds[ei];
                    if (!leafValid[leafId]) continue;
                    const daf::Size lBeg = data.leafVtxOff[leafId];
                    const daf::Size lEnd = data.leafVtxOff[leafId + 1];
                    for (daf::Size lj = lBeg; lj < lEnd; ++lj) {
                        daf::Size u = data.leafVtxIds[lj];
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

                double h = computeHIndexV4(v, data, coreV,
                                            leafNeedPivot, leafValid,
                                            oldCore, pc, bk);
                double newCore = std::min(h, oldCore);

                if (newCore != oldCore) {
                    coreV[v] = newCore;  // in-place update

                    // Enqueue neighbors with core-level filter via dual CSR.
                    const daf::Size vBeg = data.vtxLeafOff[v];
                    const daf::Size vEnd = data.vtxLeafOff[v + 1];
                    for (daf::Size ei = vBeg; ei < vEnd; ++ei) {
                        daf::Size leafId = data.vtxLeafIds[ei];
                        if (!leafValid[leafId]) continue;
                        const daf::Size lBeg = data.leafVtxOff[leafId];
                        const daf::Size lEnd = data.leafVtxOff[leafId + 1];
                        for (daf::Size lj = lBeg; lj < lEnd; ++lj) {
                            daf::Size u = data.leafVtxIds[lj];
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

// Combined entry: build dual CSR + run peel.
double * NCliqueVertexCoreDecomposition_LocalV4(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto data = NCliqueVertexCoreDecomposition_ST_V3_Build(edgeGraph, k);
    return NCliqueVertexCoreDecomposition_LocalV4_Peel(data, k);
}
