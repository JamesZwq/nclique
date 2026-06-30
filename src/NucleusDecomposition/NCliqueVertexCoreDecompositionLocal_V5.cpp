//
// Local H-index V5 — V4 + leafMinCore saturated-leaf filter.
//
// V4 reminder: per-vertex h-index recompute scans v's leaves; per leaf
// it walks the L→V CSR collecting cores for keep-min + pivot-cores
// threshold sweep. Per-leaf cost = O(|L|).
//
// V5 insight: when computing h(v) ≤ currentCore, a leaf L is "saturated"
// iff every other vertex u in L has coreV[u] ≥ currentCore. A saturated
// leaf contributes its full clique count C(pivotC, effectiveNeedPivot)
// at level bucketSize in one shot — there is no need to iterate L.
//
// Cache: leafMinCore[L] = current min coreV over all vertices in L.
//   - Monotone non-increasing: when any v's core drops, leafMinCore[L]
//     can only stay or shrink.
//   - Race-safe with atomic CAS-min: stale reads may over-saturate
//     (delay convergence one round), but the monotone-decreasing h-index
//     property guarantees eventual correctness.
//
// Per-leaf maintenance cost on a vertex update: O(|v's leaves|) = O(Σ_v)
// — paid only once when v's core actually drops, same order as V4's
// propagate-to-nbr pass.
//
// Stats reported: saturated_leaf_visits, non_saturated_leaf_visits.
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
#include <iostream>

extern double nCr[1001][401];


// Atomic min-update on a double (monotone non-increasing).
// Returns true if the slot was updated.
static inline bool atomic_min_double(std::atomic<double>& slot, double val) {
    double cur = slot.load(std::memory_order_relaxed);
    while (val < cur) {
        if (slot.compare_exchange_weak(cur, val,
                std::memory_order_relaxed, std::memory_order_relaxed)) {
            return true;
        }
    }
    return false;
}


// Per-leaf saturated counts: precomputed at init.
//   leafKeepFull[L]  = nCr[pivotC][needPivot]   (clique count if v is keep)
//   leafPivotFull[L] = nCr[pivotC-1][needPivot-1] (clique count if v is pivot)
// Both clamp to 0 when the binomial is invalid.
struct LeafFullCount {
    double keep;
    double pivot;
};


// V5 compute h-index. Identical to V4 except: per-leaf saturation
// fast-path that skips L's vertex iteration entirely.
static double computeHIndexV5(
    daf::Size v,
    const ST_V2_Data &data,
    const double *coreV,
    const std::vector<int>          &leafNeedPivot,
    const std::vector<uint8_t>      &leafValid,
    const std::vector<int>          &leafPivotCount,
    const std::vector<LeafFullCount>&leafFull,
    const std::vector<std::atomic<double>>& leafMinCore,
    double currentCore,
    std::vector<double> &pivotCores,
    std::map<int64_t, double> &buckets,
    long long &satVisits, long long &nonSatVisits)
{
    const daf::Size vBeg = data.vtxLeafOff[v];
    const daf::Size vEnd = data.vtxLeafOff[v + 1];
    if (vBeg == vEnd) return 0.0;

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

        // ---------- saturated leaf fast path ----------
        // Note: we compare leafMin to currentCore. v itself is in the
        // leaf, but currentCore == coreV[v] (just read by caller); v's
        // own core in leafMin reflects either v's current value or a
        // recent stale read — either way ≥ currentCore is fine because
        // we treat v as a participant in the clique, not a constraint.
        // The constraint is "all OTHER vertices in L have core ≥ currentCore".
        // leafMinCore includes v, so leafMin ≥ currentCore implies "all
        // others ≥ currentCore" since v's contribution to the min is
        // exactly currentCore (or higher if stale).
        double leafMin = leafMinCore[leafId].load(std::memory_order_relaxed);
        if (leafMin >= currentCore) {
            // All vertices in L (including v) have core ≥ currentCore.
            // Place the full clique-count contribution at level bucketSize.
            double contrib = v_isPivot ? leafFull[leafId].pivot
                                        : leafFull[leafId].keep;
            if (contrib > 0) {
                buckets[bucketSize] += contrib;
                rawTotalSupport += contrib;
            }
            ++satVisits;
            continue;
        }
        ++nonSatVisits;

        // ---------- V4 original non-saturated path ----------
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

    // h-index extraction (identical to V4)
    double prev_acc = 0.0;
    int64_t prev_level = bucketSize + 1;
    for (auto it = buckets.rbegin(); it != buckets.rend(); ++it) {
        int64_t c = it->first;
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
    if (prev_acc >= 1.0) {
        int64_t gap_ans = (prev_acc >= (double)(prev_level - 1))
                        ? (prev_level - 1) : (int64_t)prev_acc;
        if (gap_ans >= 1) return (double)gap_ans;
    }
    return 0.0;
}


double * NCliqueVertexCoreDecomposition_LocalV5_Peel(
    ST_V2_Data &data, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    const int nThreads = omp_get_max_threads();
    const daf::Size numLeaves = data.numLeaves;
    const daf::Size numVertices = data.numVertices;

    // Per-leaf static metadata (recomputed from CSR if not cached).
    std::vector<int>     leafNeedPivot(numLeaves);
    std::vector<int>     leafPivotCount(numLeaves);
    std::vector<uint8_t> leafValid(numLeaves);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        int keeps = 0, pivots = 0;
        const daf::Size lBeg = data.leafVtxOff[L];
        const daf::Size lEnd = data.leafVtxOff[L + 1];
        for (daf::Size lj = lBeg; lj < lEnd; ++lj) {
            if (STV3_getBit(data.leafVtxIsPivot, lj)) pivots++;
            else                                       keeps++;
        }
        leafNeedPivot[L]  = (int)k - keeps;
        leafPivotCount[L] = pivots;
        leafValid[L] = (leafNeedPivot[L] >= 0 && leafNeedPivot[L] <= pivots) ? 1 : 0;
    }

    // Precompute the saturated full-counts so the fast path is O(1).
    std::vector<LeafFullCount> leafFull(numLeaves);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        if (!leafValid[L]) { leafFull[L] = {0, 0}; continue; }
        int p  = leafPivotCount[L];
        int np = leafNeedPivot[L];
        double wKeep  = (np >= 0 && np <= p) ? nCr[p][np] : 0.0;
        double wPivot = (np >= 1 && p >= 1) ? nCr[p - 1][np - 1] : 0.0;
        leafFull[L] = {wKeep, wPivot};
    }

    // Initialize coreV from countingV.
    auto *coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = 0.0;
    if (data.countingV != nullptr) {
        for (daf::Size i = 0; i < numVertices; ++i)
            coreV[i] = data.countingV[i];
        delete[] data.countingV;
        data.countingV = nullptr;
    }

    // Initialize leafMinCore = min over all vertices in L of coreV[u].
    // Parallel over leaves.
    std::vector<std::atomic<double>> leafMinCore(numLeaves);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (daf::Size L = 0; L < numLeaves; ++L) {
        if (!leafValid[L]) {
            leafMinCore[L].store(0.0, std::memory_order_relaxed);
            continue;
        }
        double mn = 1e18;
        const daf::Size lBeg = data.leafVtxOff[L];
        const daf::Size lEnd = data.leafVtxOff[L + 1];
        for (daf::Size lj = lBeg; lj < lEnd; ++lj) {
            daf::Size u = data.leafVtxIds[lj];
            if (coreV[u] < mn) mn = coreV[u];
        }
        leafMinCore[L].store(mn, std::memory_order_relaxed);
    }

    auto t_init = std::chrono::high_resolution_clock::now();
    long long initMs = std::chrono::duration_cast<std::chrono::milliseconds>(t_init - time_start).count();
    std::cout << "LocalV5: per-leaf init (leafMin + leafFull) took " << initMs << " ms" << std::endl;

    // Per-thread scratch buffers + atomic stat counters.
    std::vector<std::vector<double>>           threadPivotCores(nThreads);
    std::vector<std::map<int64_t, double>>     threadBuckets(nThreads);
    std::vector<std::vector<daf::Size>>        threadLocalWork(nThreads);
    std::vector<long long>                     threadSatVisits(nThreads, 0);
    std::vector<long long>                     threadNonSatVisits(nThreads, 0);
    for (int t = 0; t < nThreads; ++t) {
        threadPivotCores[t].reserve(512);
        threadLocalWork[t].reserve(numVertices / nThreads + 1024);
    }

    auto *inQueue = new std::atomic<uint8_t>[numVertices];
    for (daf::Size i = 0; i < numVertices; ++i)
        inQueue[i].store(0, std::memory_order_relaxed);

    std::vector<daf::Size> currentWork;
    currentWork.reserve(numVertices);

    std::cout << "=========================begin (Local H-index V5)=========================" << std::endl;
    std::cout << "vertices: " << numVertices << ", leaves: " << numLeaves
              << ", threads: " << nThreads << std::endl;

    auto updateLeafMinsForVertex = [&](daf::Size v, double newCore) {
        const daf::Size vBeg = data.vtxLeafOff[v];
        const daf::Size vEnd = data.vtxLeafOff[v + 1];
        for (daf::Size ei = vBeg; ei < vEnd; ++ei) {
            daf::Size leafId = data.vtxLeafIds[ei];
            atomic_min_double(leafMinCore[leafId], newCore);
        }
    };

    // ======== Phase 1: parallel initial scan ========
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto &pc = threadPivotCores[tid];
        auto &bk = threadBuckets[tid];
        auto &localWork = threadLocalWork[tid];
        long long &sat = threadSatVisits[tid];
        long long &nsat = threadNonSatVisits[tid];
        localWork.clear();

        #pragma omp for schedule(dynamic, 256)
        for (daf::Size v = 0; v < numVertices; ++v) {
            if (coreV[v] <= 0) continue;

            double oldCore = coreV[v];
            double h = computeHIndexV5(v, data, coreV,
                                        leafNeedPivot, leafValid,
                                        leafPivotCount, leafFull, leafMinCore,
                                        oldCore, pc, bk, sat, nsat);
            double newCore = std::min(h, oldCore);

            if (newCore != oldCore) {
                coreV[v] = newCore;
                updateLeafMinsForVertex(v, newCore);

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

    for (int t = 0; t < nThreads; ++t) {
        currentWork.insert(currentWork.end(),
                           threadLocalWork[t].begin(), threadLocalWork[t].end());
    }

    // ======== Phase 2: async rounds ========
    long long totalRecomputeCount = 0;
    int iteration = 1;

    while (!currentWork.empty()) {
        iteration++;
        totalRecomputeCount += (long long)currentWork.size();

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto &pc = threadPivotCores[tid];
            auto &bk = threadBuckets[tid];
            auto &localWork = threadLocalWork[tid];
            long long &sat = threadSatVisits[tid];
            long long &nsat = threadNonSatVisits[tid];
            localWork.clear();

            #pragma omp for schedule(dynamic, 64)
            for (int wi = 0; wi < (int)currentWork.size(); ++wi) {
                daf::Size v = currentWork[wi];
                inQueue[v].store(0, std::memory_order_relaxed);

                double oldCore = coreV[v];
                if (oldCore <= 0) continue;

                double h = computeHIndexV5(v, data, coreV,
                                            leafNeedPivot, leafValid,
                                            leafPivotCount, leafFull, leafMinCore,
                                            oldCore, pc, bk, sat, nsat);
                double newCore = std::min(h, oldCore);

                if (newCore != oldCore) {
                    coreV[v] = newCore;
                    updateLeafMinsForVertex(v, newCore);

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

        currentWork.clear();
        for (int t = 0; t < nThreads; ++t) {
            currentWork.insert(currentWork.end(),
                               threadLocalWork[t].begin(), threadLocalWork[t].end());
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    long long elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - time_start).count();
    long long peelMs    = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_init).count();
    long long elapsedUs = std::chrono::duration_cast<std::chrono::microseconds>(t_end - time_start).count();

    long long totalSat = 0, totalNonSat = 0;
    for (int t = 0; t < nThreads; ++t) {
        totalSat += threadSatVisits[t];
        totalNonSat += threadNonSatVisits[t];
    }
    long long totalLeafVisits = totalSat + totalNonSat;
    double satRate = totalLeafVisits > 0 ? (double)totalSat / totalLeafVisits : 0.0;

    std::cout << "Local H-index V5 converged in " << iteration << " iterations, "
              << totalRecomputeCount << " re-evaluations" << std::endl;
    std::cout << "  saturated leaf visits = " << totalSat
              << ", non-saturated = " << totalNonSat
              << ", saturation rate = " << satRate << std::endl;
    std::cout << "time: " << elapsedMs << " ms (init " << initMs
              << " ms, peel " << peelMs << " ms)" << std::endl;
    std::cout << "LOCALV5_TOTAL_US: " << elapsedUs << std::endl;
    daf::phaseMark("LocalV5_peel");

    for (daf::Size i = 0; i < numVertices; ++i) {
        if (coreV[i] <= 0) coreV[i] = -1.0;
    }

    delete[] inQueue;
    return coreV;
}

double * NCliqueVertexCoreDecomposition_LocalV5(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto data = NCliqueVertexCoreDecomposition_ST_V3_Build(edgeGraph, k);
    return NCliqueVertexCoreDecomposition_LocalV5_Peel(data, k);
}
