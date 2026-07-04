//
// NCliqueVertexCoreDecompositionST_V3.cpp
//
// V3 = V2 + segfault fix (Phase 1a of the 6-step optimisation plan).
//
// Difference vs V2 (bug fix):
//   - Bucket queue switched from dense std::vector<std::vector<vid>> to
//     sparse std::map<int64_t, std::vector<vid>>.  All bucket indices
//     are int64_t.  Fixes the segfault that happened on web-Stanford at
//     s>=8 when (int)countingV[v] overflowed for support values above
//     INT_MAX, and removes the pathological dense-vector allocation
//     (was sized maxBucket+2 = up to 10^10 vector headers).
//
// Theoretical guarantee (Phase 1a non-regression):
//   - Per-relocation cost: O(log L) where L = #distinct active support
//     levels simultaneously, <= n always, typically L << n.  Heap (V2's
//     fallback) is O(log_8 n) per decrease-key.  Both reduce to O(log)
//     in worst case but bucket has lower constant and admits Phase 4's
//     batched-update trick (already in V2: dirtyMark deferral).
//   - Same number of relocations as V2 (no algorithmic change).
//   - Memory: O(L + n) vs V2's would-be O(kappa_max + n).  Strictly
//     better.
//
// (Phase 3-6 added in subsequent revisions of this file.)
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include "../PhaseLogger.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <map>

extern double nCr[1001][401];

// =====================================================================
//  Phase 0: Build all data structures via Augmented SDCT callback.
//  MUST be called BEFORE edgeGraph.beSingleEdge() mutates the graph.
// =====================================================================

ST_V2_Data NCliqueVertexCoreDecomposition_ST_V3_Build(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    ST_V3_Data d;
    d.numVertices = edgeGraph.getGraphNodeSize();

    // COO buffer
    struct COOEntry {
        daf::Size vertex;
        daf::Size leafId;
        uint8_t isPivot;
    };
    std::vector<COOEntry> cooBuf;
    cooBuf.reserve(1 << 20);

    d.countingV = new double[d.numVertices];
    memset(d.countingV, 0, d.numVertices * sizeof(double));

    // Run tree-free SDCT with inline callback
    // max_k=s(=k), min_k=1 (r=1 for vertex decomposition)
    d.numLeaves = SDCT_Augmented_NoTree(edgeGraph, k, /*min_k=*/1,
        [&](daf::Size leafId, const daf::StaticVector<int> &keepV,
            const daf::StaticVector<int> &dropV)
        {
            int pivotC = (int)dropV.size();
            int keepC = (int)keepV.size();
            int needP = (int)k - keepC;

            // Ensure metadata vectors are large enough
            if (leafId >= d.leafPivotCount.size()) {
                size_t newSz = std::max<size_t>(leafId + 1, d.leafPivotCount.size() * 2);
                d.leafPivotCount.resize(newSz, 0);
                d.leafNeedPivot.resize(newSz, 0);
            }

            d.leafPivotCount[leafId] = pivotC;
            d.leafNeedPivot[leafId] = needP;

            double wKeep = (needP >= 0 && needP <= pivotC) ? nCr[pivotC][needP] : 0.0;
            double wPivot = (needP >= 1 && pivotC >= 1) ? nCr[pivotC - 1][needP - 1] : 0.0;

            for (int i = 0; i < keepC; ++i) {
                daf::Size v = keepV[i];
                d.countingV[v] += wKeep;
                cooBuf.push_back({v, leafId, 0});
            }
            for (int i = 0; i < pivotC; ++i) {
                daf::Size v = dropV[i];
                d.countingV[v] += wPivot;
                cooBuf.push_back({v, leafId, 1});
            }
        });

    d.leafPivotCount.resize(d.numLeaves);
    d.leafNeedPivot.resize(d.numLeaves);

    auto time_sdct = std::chrono::high_resolution_clock::now();
    std::cout << "ST_V3: SDCT+callback took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_sdct - time_start).count()
              << " ms, leaves=" << d.numLeaves << ", COO entries=" << cooBuf.size() << std::endl;
    // Mark "Sigma" = total vertex-leaf incidences for this run.
    daf::phaseMark("STV3_SDCT_walk", (long)(cooBuf.capacity() * sizeof(COOEntry)));

    // --- Build dual CSR from COO (packed) ---
    // Offsets / per-vertex counts must be size_t: cumulative offset = total Σ
    // exceeds uint32 max on billion-edge graphs at moderate s.
    //
    // Memory-conscious construction order:
    //   (1) build vtxLeafOff via prefix-sum over cooBuf
    //   (2) fill packed vtxLeafIds + vtxLeafIsPivot while tallying per-leaf counts
    //   (3) FREE cooBuf  (saves ~12 bytes per incidence at billion-edge s)
    //   (4) build packed leafVtxIds + leafVtxIsPivot by transposing vtxLeafIds
    //
    // Packed layout: {ids: uint32[Σ], isPivot: bit[Σ]} = 4.125 bytes/incidence,
    // half of the 8-byte legacy {id, role-byte} struct.
    d.vtxLeafOff.assign(d.numVertices + 2, 0);
    for (auto &e : cooBuf)
        if (e.vertex < d.numVertices) d.vtxLeafOff[e.vertex + 1]++;
    for (daf::Size i = 1; i <= d.numVertices; ++i)
        d.vtxLeafOff[i] += d.vtxLeafOff[i - 1];

    const size_t totalIncidence = d.vtxLeafOff[d.numVertices];
    d.vtxLeafIds.resize(totalIncidence);
    d.vtxLeafIsPivot.assign((totalIncidence + 63) >> 6, 0);
    d.leafVtxOff.assign(d.numLeaves + 1, 0);
    {
        std::vector<size_t> pos(d.numVertices, 0);
        for (auto &e : cooBuf) {
            daf::Size v = e.vertex;
            if (v < d.numVertices) {
                size_t p = d.vtxLeafOff[v] + pos[v]++;
                d.vtxLeafIds[p] = e.leafId;
                if (e.isPivot) STV3_setBit(d.vtxLeafIsPivot, p);
                if (e.leafId < d.numLeaves) d.leafVtxOff[e.leafId + 1]++;
            }
        }
    }

    // Free cooBuf before allocating leafVtx* arrays. On billion-edge graphs at
    // moderate s, this drops the build-phase peak by ~Σ * sizeof(COOEntry).
    std::vector<COOEntry>().swap(cooBuf);

    for (size_t i = 1; i <= d.numLeaves; ++i)
        d.leafVtxOff[i] += d.leafVtxOff[i - 1];
    const size_t totalLeafIncidence = d.leafVtxOff[d.numLeaves];
    d.leafVtxIds.resize(totalLeafIncidence);
    d.leafVtxIsPivot.assign((totalLeafIncidence + 63) >> 6, 0);
    {
        std::vector<size_t> pos(d.numLeaves, 0);
        for (daf::Size v = 0; v < d.numVertices; ++v) {
            const size_t row_begin = d.vtxLeafOff[v];
            const size_t row_end = d.vtxLeafOff[v + 1];
            for (size_t k = row_begin; k < row_end; ++k) {
                daf::Size L = d.vtxLeafIds[k];
                if (L < d.numLeaves) {
                    size_t p = d.leafVtxOff[L] + pos[L]++;
                    d.leafVtxIds[p] = v;
                    if (STV3_getBit(d.vtxLeafIsPivot, k)) STV3_setBit(d.leafVtxIsPivot, p);
                }
            }
        }
    }

    auto time_csr = std::chrono::high_resolution_clock::now();
    std::cout << "ST_V3: dual CSR built in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_csr - time_sdct).count()
              << " ms" << std::endl;
    daf::log_memory("ST_V2 after CSR build");

    // Component-byte attribution for the dual CSR + initial-support array.
    const long bytesVtxLeaf = (long)(d.vtxLeafOff.capacity() * sizeof(d.vtxLeafOff[0])
                            + d.vtxLeafIds.capacity() * sizeof(d.vtxLeafIds[0])
                            + d.vtxLeafIsPivot.capacity() * sizeof(d.vtxLeafIsPivot[0]));
    const long bytesLeafVtx = (long)(d.leafVtxOff.capacity() * sizeof(d.leafVtxOff[0])
                            + d.leafVtxIds.capacity() * sizeof(d.leafVtxIds[0])
                            + d.leafVtxIsPivot.capacity() * sizeof(d.leafVtxIsPivot[0]));
    const long bytesSupport = (long)(d.numVertices * sizeof(double));
    const long bytesLeafMeta = (long)(d.leafPivotCount.capacity() * sizeof(int)
                            + d.leafNeedPivot.capacity() * sizeof(int));
    daf::phaseMark("STV3_CSR_build", bytesVtxLeaf + bytesLeafVtx + bytesSupport + bytesLeafMeta);

    return d;
}

// =====================================================================
//  Phase 1: Peeling.  Graph-independent — uses only dual CSR + metadata.
// =====================================================================

double * NCliqueVertexCoreDecomposition_ST_V3_Peel(ST_V2_Data &d, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numVertices = d.numVertices;
    const size_t numLeaves = d.numLeaves;
    double *countingV = d.countingV;
    d.countingV = nullptr; // transfer ownership

    std::vector<uint8_t> leafAlive(numLeaves);
    std::vector<int> leafRemainPivots(numLeaves);
    for (size_t L = 0; L < numLeaves; ++L) {
        leafRemainPivots[L] = d.leafPivotCount[L];
        int np = d.leafNeedPivot[L];
        leafAlive[L] = (np >= 0 && np <= d.leafPivotCount[L]) ? 1 : 0;
    }

    auto coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;
    // Optional peel-order (pop-rank) capture for the v5 order-certificate
    // experiment. Rank = position in the pop sequence (0 = first popped).
    // No effect on the algorithm; only recorded when PIVOTER_DUMP_ORDER set.
    const bool _dumpOrder = std::getenv("PIVOTER_DUMP_ORDER") != nullptr;
    std::vector<daf::Size> _popRank;
    daf::Size _popCtr = 0;
    if (_dumpOrder) _popRank.assign(numVertices, (daf::Size)-1);

    // ---- Sparse bucket queue (Phase 1a fix vs V2) ----
    // std::map<int64_t, vector<vid>>; per-relocation O(log L).
    // Theoretical non-regression vs heap: heap is O(log_8 n) per
    // decrease-key.  Map is O(log L) where L <= n is #active support
    // levels.  Asymptotically same; constants close.  Wins big when
    // (a) Phase 4's batched updates reduce #relocations and (b) kappa_max
    // exceeds INT_MAX (no overflow / no giant dense allocation).
    static constexpr int64_t kBucketKeyClamp = (int64_t)1e18;
    auto supportToKey = [](double s) -> int64_t {
        if (s <= 0.0) return 0;
        if (s >= (double)kBucketKeyClamp) return kBucketKeyClamp;
        return (int64_t)s;
    };

    std::map<int64_t, std::vector<daf::Size>> buckets;
    std::vector<int64_t> bucket_of(numVertices, -1);
    std::vector<daf::Size> pos_in_bucket(numVertices, 0);
    std::vector<uint8_t> vertexInHeap(numVertices, 0);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] <= 0) continue;
        int64_t b = supportToKey(countingV[i]);
        auto &lst = buckets[b];
        bucket_of[i] = b;
        pos_in_bucket[i] = lst.size();
        lst.push_back(i);
        vertexInHeap[i] = 1;
        remainingInHeap++;
    }

    auto bucketMove = [&](daf::Size id) {
        int64_t newB = supportToKey(countingV[id]);
        int64_t oldB = bucket_of[id];
        if (newB == oldB) return;
        auto oldIt = buckets.find(oldB);
        auto &oldVec = oldIt->second;
        daf::Size myPos = pos_in_bucket[id];
        if (myPos + 1 < oldVec.size()) {
            daf::Size last = oldVec.back();
            oldVec[myPos] = last;
            pos_in_bucket[last] = myPos;
        }
        oldVec.pop_back();
        if (oldVec.empty()) buckets.erase(oldIt);
        bucket_of[id] = newB;
        auto &newVec = buckets[newB];
        pos_in_bucket[id] = newVec.size();
        newVec.push_back(id);
    };

    // Per-leaf batch tracking
    std::vector<int> leafRemovedPivots(numLeaves, 0);
    std::vector<uint8_t> leafDies(numLeaves, 0);
    std::vector<uint8_t> leafAffected(numLeaves, 0);
    std::vector<daf::Size> affectedLeaves;
    affectedLeaves.reserve(4096);

    // Component-byte attribution for the bucket array + per-vertex book-keeping.
    {
        long bytesBucket = 0;
        for (auto &kv : buckets) bytesBucket += (long)(kv.second.capacity() * sizeof(daf::Size));
        bytesBucket += (long)buckets.size() * 64;  // red-black tree node overhead
        const long bytesVertexAux = (long)(bucket_of.capacity() * sizeof(int64_t)
                                  + pos_in_bucket.capacity() * sizeof(daf::Size)
                                  + vertexInHeap.capacity() * sizeof(uint8_t));
        const long bytesLeafAux = (long)(leafAlive.capacity() * sizeof(uint8_t)
                                + leafRemainPivots.capacity() * sizeof(int)
                                + leafRemovedPivots.capacity() * sizeof(int)
                                + leafDies.capacity() * sizeof(uint8_t)
                                + leafAffected.capacity() * sizeof(uint8_t));
        daf::phaseMark("STV3_peel_init", bytesBucket + bytesVertexAux + bytesLeafAux);
    }

    std::vector<uint8_t> dirtyMark(numVertices, 0);
    std::vector<daf::Size> dirtyVertices;
    dirtyVertices.reserve(8192);

    daf::StaticVector<daf::Size> currentRemoveVertexIds(numVertices);

    std::cout << "=========================begin=========================" << std::endl;
    std::cout << "vertices in heap: " << remainingInHeap << std::endl;

    double minCore = 0;

    while (remainingInHeap > 0) {
        while (!buckets.empty() && buckets.begin()->second.empty()) {
            buckets.erase(buckets.begin());
        }
        if (buckets.empty()) break;

        int64_t curKey = buckets.begin()->first;
        minCore = std::max((double)curKey, minCore);

        while (!buckets.empty()
               && (double)buckets.begin()->first <= minCore) {
            auto &lst = buckets.begin()->second;
            while (!lst.empty()) {
                auto id = lst.back();
                lst.pop_back();
                vertexInHeap[id] = 0;
                currentRemoveVertexIds.push_back(id);
                coreV[id] = minCore;
                if (_dumpOrder) _popRank[id] = _popCtr++;
                remainingInHeap--;
            }
            buckets.erase(buckets.begin());
        }

        if (remainingInHeap == 0) break;

        // Phase 1: find affected leaves via Vertex→Leaf CSR
        for (int vi = 0; vi < (int)currentRemoveVertexIds.size(); ++vi) {
            auto v = currentRemoveVertexIds[vi];
            const daf::Size begin = d.vtxLeafOff[v];
            const daf::Size end = d.vtxLeafOff[v + 1];
            if (vi + 1 < (int)currentRemoveVertexIds.size()) {
                auto nextV = currentRemoveVertexIds[vi + 1];
                __builtin_prefetch(&d.vtxLeafIds[d.vtxLeafOff[nextV]], 0, 1);
            }
            for (daf::Size ei = begin; ei < end; ++ei) {
                daf::Size leafId = d.vtxLeafIds[ei];
                if (!leafAlive[leafId]) continue;
                if (!leafAffected[leafId]) {
                    leafAffected[leafId] = 1;
                    affectedLeaves.push_back(leafId);
                }
                if (!STV3_getBit(d.vtxLeafIsPivot, ei)) {
                    leafDies[leafId] = 1;
                } else {
                    leafRemovedPivots[leafId]++;
                }
            }
        }

        // Phase 2: process each affected leaf via Leaf→Vertex CSR
        for (auto leafId : affectedLeaves) {
            int old_rp = leafRemainPivots[leafId];
            int np = d.leafNeedPivot[leafId];
            int new_rp = old_rp - leafRemovedPivots[leafId];
            bool dies = leafDies[leafId] || np > new_rp || new_rp < 0;

            double deltaKeep, deltaPivot;
            if (dies) {
                deltaKeep = (np >= 0 && np <= old_rp) ? nCr[old_rp][np] : 0.0;
                deltaPivot = (np >= 1 && old_rp >= 1) ? nCr[old_rp - 1][np - 1] : 0.0;
            } else {
                deltaKeep = nCr[old_rp][np] - nCr[new_rp][np];
                deltaPivot = (np >= 1) ? nCr[old_rp - 1][np - 1] - nCr[new_rp - 1][np - 1] : 0.0;
            }

            const daf::Size lBegin = d.leafVtxOff[leafId];
            const daf::Size lEnd = d.leafVtxOff[leafId + 1];
            for (daf::Size li = lBegin; li < lEnd; ++li) {
                daf::Size vtx = d.leafVtxIds[li];
                if (!vertexInHeap[vtx]) continue;
                double delta = STV3_getBit(d.leafVtxIsPivot, li) ? deltaPivot : deltaKeep;
                if (delta > 0) {
                    countingV[vtx] -= delta;
                    if (countingV[vtx] < 0) countingV[vtx] = 0;
                    if (!dirtyMark[vtx]) {
                        dirtyMark[vtx] = 1;
                        dirtyVertices.push_back(vtx);
                    }
                }
            }

            leafRemainPivots[leafId] = new_rp;
            if (dies) leafAlive[leafId] = 0;
        }

        // Batch bucket moves
        for (auto v : dirtyVertices) {
            if (vertexInHeap[v]) bucketMove(v);
            dirtyMark[v] = 0;
        }
        dirtyVertices.clear();

        // Cleanup
        for (auto leafId : affectedLeaves) {
            leafRemovedPivots[leafId] = 0;
            leafDies[leafId] = 0;
            leafAffected[leafId] = 0;
        }
        affectedLeaves.clear();
        currentRemoveVertexIds.clear();
    }

    {
        auto _peel_us = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - time_start).count();
        std::cout << "ST_V2 peeling time: " << (_peel_us / 1000) << " ms" << std::endl;
        // Microsecond-precision line for the bench harness (small-peel cells).
        std::cout << "STV3_PEEL_US: " << _peel_us << std::endl;
    }

    daf::phaseMark("STV3_peel_loop");

    if (_dumpOrder) {
        const char *op = std::getenv("PIVOTER_DUMP_ORDER");
        if (FILE *of = std::fopen(op, "w")) {
            std::fprintf(of, "# internal_id\tpop_rank\n");
            for (daf::Size i = 0; i < numVertices; ++i)
                if (coreV[i] >= 0)
                    std::fprintf(of, "%llu\t%llu\n", (unsigned long long)i,
                                 (unsigned long long)_popRank[i]);
            std::fclose(of);
            std::cerr << "PIVOTER_DUMP_ORDER: wrote " << _popCtr
                      << " pop ranks to " << op << std::endl;
        }
    }

    delete[] countingV;
    currentRemoveVertexIds.free();
    return coreV;
}

// =====================================================================
//  Combined entry point (for use when graph is NOT yet mutated)
// =====================================================================
double * NCliqueVertexCoreDecomposition_ST_V3(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto data = NCliqueVertexCoreDecomposition_ST_V3_Build(edgeGraph, k);
    return NCliqueVertexCoreDecomposition_ST_V3_Peel(data, k);
}
