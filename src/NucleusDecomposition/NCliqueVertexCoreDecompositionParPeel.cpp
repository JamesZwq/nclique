//
// NCliqueVertexCoreDecompositionParPeel.cpp
//
// Parallel peel for SPIN★ (r=1 (1,s)-nucleus) on the preserved compact
// index (V3 dual CSR). Within each peeling round (one minCore level),
// the O(Σ)-work Phase 1 (affected-leaf marking) and Phase 2 (closed-form
// binomial delta propagation) are parallelized. Rounds remain sequential
// because minCore is monotone non-decreasing.
//
// Distinction vs prior parallel nucleus decomposition (Shi, Dhulipala,
// Shun, VLDB 2021): that work re-enumerates s-cliques during each peel
// update, giving O(mα^{s-2} + ρ log n) work = Θ(#s-cliques). We peel on
// the pivot-compressed compact index with closed-form per-leaf deltas,
// giving O(Σ + ρ log n) work where Σ ≤ s·#s-cliques and Σ << #s-cliques
// on clique-dense graphs. No clique is ever re-enumerated during peel.
//
// Correctness (bit-identical to serial V3):
//   - Drain (serial): all vertices with key ≤ minCore get coreV = minCore.
//   - Phase 1 (parallel over batch): leafRemovedPivots is an atomic SUM,
//     leafDies an atomic OR — both commutative, so the final per-leaf
//     state equals the serial result regardless of thread interleaving.
//   - Phase 2 (parallel over affected leaves, after a barrier): each
//     vertex counter update is max(0, c - delta) applied atomically;
//     since max(0, max(0, c-a)-b) = max(0, c-a-b), the composition is
//     order-independent and equals max(0, c - Σdeltas) — identical to V3.
//   - Bucket moves (serial): re-bucket every vertex whose counter changed.
//
// Contention-reduction techniques borrowed from the literature:
//   - Per-thread local buffers for affected-leaf and dirty-vertex lists
//     (list-buffer, Shi et al. VLDB'21) — avoids contention on shared
//     output vectors during enumeration.
//   - "First-to-mark" CAS on leafAffected[]/dirtyMark[] so each leaf /
//     vertex is appended to exactly one thread-local list.
//

#include "NCliqueCoreDecomposition.h"
#include "../PhaseLogger.h"
#include <chrono>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <atomic>
#include <iostream>
#include <omp.h>

extern double nCr[1001][401];

// Atomic double "subtract then clamp at 0", relaxed ordering, via a CAS
// loop on the 64-bit pattern. Uses __atomic builtins (always available on
// clang/gcc, no libc++ atomic_ref dependency).
static inline void atomicSubClamp(double &slot, double delta) {
    uint64_t *bits = reinterpret_cast<uint64_t *>(&slot);
    uint64_t curBits = __atomic_load_n(bits, __ATOMIC_RELAXED);
    for (;;) {
        double cur;
        std::memcpy(&cur, &curBits, sizeof(double));
        double next = cur - delta;
        if (next < 0.0) next = 0.0;
        uint64_t nextBits;
        std::memcpy(&nextBits, &next, sizeof(double));
        if (__atomic_compare_exchange_n(bits, &curBits, nextBits, true,
                __ATOMIC_RELAXED, __ATOMIC_RELAXED))
            break;
        // curBits updated with the current value on failure; retry.
    }
}

// First-to-mark CAS on a uint8_t flag: returns true iff this call flipped
// it from 0 to 1 (i.e. this thread is the unique "owner" of the mark).
static inline bool markFirst(uint8_t &flag) {
    uint8_t expected = 0;
    return __atomic_compare_exchange_n(&flag, &expected, (uint8_t)1, false,
                __ATOMIC_RELAXED, __ATOMIC_RELAXED);
}

double * NCliqueVertexCoreDecomposition_ParPeel(ST_V2_Data &d, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numVertices = d.numVertices;
    const size_t numLeaves = d.numLeaves;
    double *countingV = d.countingV;
    d.countingV = nullptr;

    const int nThreads = omp_get_max_threads();

    std::vector<uint8_t> leafAlive(numLeaves);
    std::vector<int>     leafRemainPivots(numLeaves);
    for (size_t L = 0; L < numLeaves; ++L) {
        leafRemainPivots[L] = d.leafPivotCount[L];
        int np = d.leafNeedPivot[L];
        leafAlive[L] = (np >= 0 && np <= d.leafPivotCount[L]) ? 1 : 0;
    }

    auto coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;

    static constexpr int64_t kBucketKeyClamp = (int64_t)1e18;
    auto supportToKey = [](double s) -> int64_t {
        if (s <= 0.0) return 0;
        if (s >= (double)kBucketKeyClamp) return kBucketKeyClamp;
        return (int64_t)s;
    };

    // Tombstone sentinel: an invalid vertex id placed in a bucket slot when
    // its vertex is re-bucketed (decrease-key). Removal is thus a single slot
    // write (parallelizable: each dirty vertex owns a distinct slot), and the
    // drain skips sentinels. bucket_vec[v] points to the std::vector holding v
    // (the mapped value is node-stable in std::map, so the pointer survives
    // push_back reallocations and only dies when the key is erased on drain —
    // by which point v has been peeled and is never touched again).
    static constexpr daf::Size SENTINEL = (daf::Size)-1;

    std::map<int64_t, std::vector<daf::Size>> buckets;
    std::vector<int64_t> bucket_of(numVertices, -1);
    std::vector<daf::Size> pos_in_bucket(numVertices, 0);
    std::vector<std::vector<daf::Size>*> bucket_vec(numVertices, nullptr);
    std::vector<uint8_t> vertexInHeap(numVertices, 0);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] <= 0) continue;
        int64_t b = supportToKey(countingV[i]);
        auto &lst = buckets[b];
        bucket_of[i] = b;
        pos_in_bucket[i] = lst.size();
        bucket_vec[i] = &lst;
        lst.push_back(i);
        vertexInHeap[i] = 1;
        remainingInHeap++;
    }

    // Serial insert of a dirty vertex into its new (decreased-key) bucket.
    // The matching removal was already done in Phase 2 via tombstoning.
    auto bucketInsert = [&](daf::Size id) {
        int64_t newB = supportToKey(countingV[id]);
        auto &newVec = buckets[newB];
        bucket_of[id] = newB;
        pos_in_bucket[id] = newVec.size();
        bucket_vec[id] = &newVec;
        newVec.push_back(id);
    };

    // Per-leaf affected flag + per-vertex dirty flag (CAS-deduped).
    std::vector<uint8_t> leafAffected(numLeaves, 0);
    std::vector<uint8_t> dirtyMark(numVertices, 0);

    // Per-thread persistent buffers (cleared, not freed, each round).
    std::vector<std::vector<daf::Size>> tAffected(nThreads);
    std::vector<std::vector<daf::Size>> tDirty(nThreads);
    for (int t = 0; t < nThreads; ++t) {
        tAffected[t].reserve(4096);
        tDirty[t].reserve(8192);
    }

    // Global concatenated lists for the serial bucket-move + cleanup.
    std::vector<daf::Size> affectedLeaves;
    std::vector<daf::Size> dirtyVertices;
    affectedLeaves.reserve(1 << 16);
    dirtyVertices.reserve(1 << 16);

    std::vector<daf::Size> currentRemoveVertexIds;
    currentRemoveVertexIds.reserve(1 << 16);

    std::cout << "=========================begin (ParPeel, "
              << nThreads << " threads)=========================" << std::endl;
    std::cout << "vertices in heap: " << remainingInHeap << std::endl;

    // Phase timing accumulators (microseconds).
    long long t_phase1 = 0, t_phase2 = 0, t_bucket = 0, t_drain = 0;
    long long t_p1work = 0, t_p2work = 0;  // actual parallel-for wall time
    long long t_cleanup = 0;               // flag-reset cleanup wall time
    using clk = std::chrono::high_resolution_clock;
    clk::time_point ts_p1start, ts_p2start;

    double minCore = 0;
    long long roundCount = 0;

    // Shared loop state (written in single regions, read after barriers).
    bool   finished = false;
    int    batchSize = 0;
    int    numAffected = 0;

    // Persistent thread pool: one parallel region for the whole peel. Phase
    // boundaries use cheap omp barriers (hundreds of cycles) instead of
    // fork/join (~15k cycles each). Serial sections (drain, concatenation,
    // bucket moves) run in `single` so threads stay alive across all rounds.
    #pragma omp parallel num_threads(nThreads)
    {
        const int tid = omp_get_thread_num();
        while (true) {
            // ---- Drain (serial): collect this round's batch ----
            #pragma omp single
            {
                auto td0 = clk::now();
                // Find the smallest bucket that still holds a live (non-
                // sentinel) vertex. All-sentinel buckets (every vertex
                // re-bucketed away) are erased without advancing minCore.
                bool foundLive = false;
                int64_t curKey = 0;
                while (!buckets.empty()) {
                    auto &frontVec = buckets.begin()->second;
                    bool hasLive = false;
                    for (daf::Size id : frontVec)
                        if (id != SENTINEL) { hasLive = true; break; }
                    if (!hasLive) { buckets.erase(buckets.begin()); continue; }
                    curKey = buckets.begin()->first;
                    foundLive = true;
                    break;
                }
                if (!foundLive || remainingInHeap == 0) {
                    finished = true;
                } else {
                    minCore = std::max((double)curKey, minCore);
                    currentRemoveVertexIds.clear();
                    while (!buckets.empty()
                           && (double)buckets.begin()->first <= minCore) {
                        auto &lst = buckets.begin()->second;
                        for (daf::Size id : lst) {
                            if (id == SENTINEL) continue;
                            vertexInHeap[id] = 0;
                            bucket_of[id] = -1;
                            currentRemoveVertexIds.push_back(id);
                            coreV[id] = minCore;
                            remainingInHeap--;
                        }
                        buckets.erase(buckets.begin());
                    }
                    if (remainingInHeap == 0) finished = true;
                    batchSize = (int)currentRemoveVertexIds.size();
                    ++roundCount;
                }
                t_drain += std::chrono::duration_cast<std::chrono::microseconds>(clk::now() - td0).count();
                ts_p1start = clk::now();
            }
            // implicit barrier after single: finished/batchSize visible to all.
            if (finished) break;

            // ---- Phase 1 (parallel): mark affected leaves ONLY ----
            // No atomic accumulation here. Each (v, leaf) incidence does at
            // most one CAS (markFirst); deduped so each affected leaf is
            // appended to exactly one thread-local list. The per-leaf counts
            // are recomputed in Phase 2 by scanning the leaf (vertexInHeap is
            // read-only there) — this removes the light-work / heavy-atomic
            // pass that does not scale.
            tAffected[tid].clear();
            {
                auto &myAff = tAffected[tid];
                #pragma omp for schedule(dynamic, 64)
                for (int vi = 0; vi < batchSize; ++vi) {
                    daf::Size v = currentRemoveVertexIds[vi];
                    const daf::Size begin = d.vtxLeafOff[v];
                    const daf::Size end   = d.vtxLeafOff[v + 1];
                    for (daf::Size ei = begin; ei < end; ++ei) {
                        daf::Size leafId = d.vtxLeafIds[ei];
                        if (!leafAlive[leafId]) continue;
                        if (markFirst(leafAffected[leafId]))
                            myAff.push_back(leafId);
                    }
                }
                // implicit barrier after omp for
            }

            // ---- Concatenate affected leaves (serial) ----
            #pragma omp single
            {
                t_p1work += std::chrono::duration_cast<std::chrono::microseconds>(clk::now() - ts_p1start).count();
                auto tc = clk::now();
                affectedLeaves.clear();
                for (int t = 0; t < nThreads; ++t)
                    affectedLeaves.insert(affectedLeaves.end(),
                                          tAffected[t].begin(), tAffected[t].end());
                numAffected = (int)affectedLeaves.size();
                t_phase1 += std::chrono::duration_cast<std::chrono::microseconds>(clk::now() - tc).count();
                ts_p2start = clk::now();
            }

            // ---- Phase 2 (parallel): propagate closed-form deltas ----
            tDirty[tid].clear();
            {
                auto &myDirty = tDirty[tid];
                #pragma omp for schedule(dynamic, 32)
                for (int ai = 0; ai < numAffected; ++ai) {
                    daf::Size leafId = affectedLeaves[ai];
                    const int old_rp = leafRemainPivots[leafId];
                    const int np     = d.leafNeedPivot[leafId];
                    const daf::Size lBegin = d.leafVtxOff[leafId];
                    const daf::Size lEnd   = d.leafVtxOff[leafId + 1];

                    // Pass 1: recompute alive-pivot count and keep-death by
                    // scanning the leaf. vertexInHeap is read-only here
                    // (set during the serial drain), so no contention.
                    int new_rp = 0;
                    bool keepDied = false;
                    for (daf::Size li = lBegin; li < lEnd; ++li) {
                        daf::Size vtx = d.leafVtxIds[li];
                        bool alive = vertexInHeap[vtx];
                        if (STV3_getBit(d.leafVtxIsPivot, li)) {
                            if (alive) ++new_rp;
                        } else {
                            if (!alive) keepDied = true;
                        }
                    }
                    bool dies = keepDied || np > new_rp || new_rp < 0;

                    double deltaKeep, deltaPivot;
                    if (dies) {
                        deltaKeep  = (np >= 0 && np <= old_rp) ? nCr[old_rp][np] : 0.0;
                        deltaPivot = (np >= 1 && old_rp >= 1) ? nCr[old_rp - 1][np - 1] : 0.0;
                    } else {
                        deltaKeep  = nCr[old_rp][np] - nCr[new_rp][np];
                        deltaPivot = (np >= 1) ? nCr[old_rp - 1][np - 1] - nCr[new_rp - 1][np - 1] : 0.0;
                    }

                    // Pass 2: scatter the closed-form delta to alive vertices.
                    for (daf::Size li = lBegin; li < lEnd; ++li) {
                        daf::Size vtx = d.leafVtxIds[li];
                        if (!vertexInHeap[vtx]) continue;
                        double delta = STV3_getBit(d.leafVtxIsPivot, li) ? deltaPivot : deltaKeep;
                        if (delta > 0) {
                            atomicSubClamp(countingV[vtx], delta);
                            if (markFirst(dirtyMark[vtx])) {
                                // Parallel tombstone removal: write SENTINEL to
                                // vtx's own (unique) slot in its current bucket.
                                // The map is read-only during Phase 2, and each
                                // dirty vertex owns a distinct slot, so this is
                                // race-free. Re-insertion happens serially after.
                                (*bucket_vec[vtx])[pos_in_bucket[vtx]] = SENTINEL;
                                myDirty.push_back(vtx);
                            }
                        }
                    }
                    leafRemainPivots[leafId] = new_rp;   // owner-only write
                    if (dies) leafAlive[leafId] = 0;      // owner-only write
                }
                // implicit barrier after omp for
            }

            // ---- Concatenate dirty + bucket moves + cleanup (serial) ----
            #pragma omp single
            {
                t_p2work += std::chrono::duration_cast<std::chrono::microseconds>(clk::now() - ts_p2start).count();
                auto tp2c = clk::now();
                dirtyVertices.clear();
                for (int t = 0; t < nThreads; ++t)
                    dirtyVertices.insert(dirtyVertices.end(),
                                         tDirty[t].begin(), tDirty[t].end());
                t_phase2 += std::chrono::duration_cast<std::chrono::microseconds>(clk::now() - tp2c).count();

                auto tb0 = clk::now();
                // Serial re-insert: tombstone removal was already done in
                // Phase 2; here we only append each dirty vertex to its new
                // (decreased-key) bucket.
                for (auto v : dirtyVertices) {
                    if (vertexInHeap[v]) bucketInsert(v);
                }
                t_bucket += std::chrono::duration_cast<std::chrono::microseconds>(clk::now() - tb0).count();
                auto tcl = clk::now();
                for (auto v : dirtyVertices) dirtyMark[v] = 0;
                for (auto leafId : affectedLeaves) leafAffected[leafId] = 0;
                t_cleanup += std::chrono::duration_cast<std::chrono::microseconds>(clk::now() - tcl).count();
            }
            // implicit barrier after single
        }
    }

    auto _peel_us = std::chrono::duration_cast<std::chrono::microseconds>(
        clk::now() - time_start).count();
    std::cout << "ParPeel peeling time: " << (_peel_us / 1000) << " ms" << std::endl;
    std::cout << "PARPEEL_PEEL_US: " << _peel_us << std::endl;
    std::cout << "  [breakdown] drain=" << (t_drain/1000) << "ms"
              << " p1work=" << (t_p1work/1000) << "ms"
              << " p1concat=" << (t_phase1/1000) << "ms"
              << " p2work=" << (t_p2work/1000) << "ms"
              << " p2concat=" << (t_phase2/1000) << "ms"
              << " bucket=" << (t_bucket/1000) << "ms"
              << " cleanup=" << (t_cleanup/1000) << "ms"
              << " rounds=" << roundCount << std::endl;
    daf::phaseMark("ParPeel_peel_loop");

    delete[] countingV;
    return coreV;
}
