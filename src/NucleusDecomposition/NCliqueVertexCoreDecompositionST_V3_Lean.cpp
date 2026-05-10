//
// NCliqueVertexCoreDecompositionST_V3_Lean.cpp
//
// SPIN★ Lean — memory-tighter variant of ST_V3.
//
// Drops per-leaf persistent state arrays:
//   - leafPivotCount, leafNeedPivot   (build-time, 8 bytes/leaf)
//   - leafAlive, leafRemainPivots     (peel-time persistent, 5 bytes/leaf)
// Replaced by 1-bit-per-leaf dead bitmap (~0.125 bytes/leaf) and on-the-fly
// recomputation of np and old_rp from a per-event scan of each affected
// leaf's vertex set.
//
// Net saving: ~12.875 bytes per CPI leaf (permanent state).
// Cost: each affected leaf gets one extra pass over its vertex list per
// event to count alive pivots, ~1.5-2x peel CPU vs ST_V3.
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

ST_V2_Data NCliqueVertexCoreDecomposition_ST_V3_Lean_Build(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    ST_V2_Data d;
    d.numVertices = edgeGraph.getGraphNodeSize();

    struct COOEntry { daf::Size vertex; daf::Size leafId; uint8_t isPivot; };
    std::vector<COOEntry> cooBuf;
    cooBuf.reserve(1 << 20);

    d.countingV = new double[d.numVertices];
    memset(d.countingV, 0, d.numVertices * sizeof(double));

    d.numLeaves = SDCT_Augmented_NoTree(edgeGraph, k, /*min_k=*/1,
        [&](daf::Size leafId, const daf::StaticVector<int> &keepV,
            const daf::StaticVector<int> &dropV)
        {
            int pivotC = (int)dropV.size();
            int keepC = (int)keepV.size();
            int needP = (int)k - keepC;

            // LEAN: skip leafPivotCount/leafNeedPivot writes — derived at peel time.

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

    auto time_sdct = std::chrono::high_resolution_clock::now();
    std::cout << "ST_V3_Lean: SDCT+callback took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_sdct - time_start).count()
              << " ms, leaves=" << d.numLeaves << ", COO entries=" << cooBuf.size() << std::endl;
    daf::phaseMark("STV3_LEAN_SDCT_walk", (long)(cooBuf.capacity() * sizeof(COOEntry)));

    // Build dual CSR — identical to ST_V3.
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
    std::cout << "ST_V3_Lean: dual CSR built in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_csr - time_sdct).count()
              << " ms" << std::endl;
    daf::log_memory("ST_V3_Lean after CSR build");

    return d;
}

double * NCliqueVertexCoreDecomposition_ST_V3_Lean_Peel(ST_V2_Data &d, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numVertices = d.numVertices;
    const size_t numLeaves = d.numLeaves;
    double *countingV = d.countingV;
    d.countingV = nullptr;

    auto coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;

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

    // LEAN: per-event scratch. Cleared at the end of each event via affectedLeaves walk.
    std::vector<int> leafRemovedPivots(numLeaves, 0);
    std::vector<uint8_t> leafKeepRemoved(numLeaves, 0);
    std::vector<uint8_t> leafAffected(numLeaves, 0);
    std::vector<daf::Size> affectedLeaves;
    affectedLeaves.reserve(4096);

    // LEAN: persistent dead-leaf bitmap (1 bit/leaf vs ST_V3's 1 byte/leaf).
    std::vector<uint64_t> leafDead((numLeaves + 63) >> 6, 0);
    auto leafIsDead = [&](size_t L) -> bool {
        return (leafDead[L >> 6] >> (L & 63)) & 1;
    };
    auto markLeafDead = [&](size_t L) {
        leafDead[L >> 6] |= (uint64_t(1) << (L & 63));
    };

    {
        long bytesBucket = 0;
        for (auto &kv : buckets) bytesBucket += (long)(kv.second.capacity() * sizeof(daf::Size));
        bytesBucket += (long)buckets.size() * 64;
        const long bytesVertexAux = (long)(bucket_of.capacity() * sizeof(int64_t)
                                  + pos_in_bucket.capacity() * sizeof(daf::Size)
                                  + vertexInHeap.capacity() * sizeof(uint8_t));
        const long bytesLeafAux = (long)(leafDead.capacity() * sizeof(uint64_t)
                                + leafRemovedPivots.capacity() * sizeof(int)
                                + leafKeepRemoved.capacity() * sizeof(uint8_t)
                                + leafAffected.capacity() * sizeof(uint8_t));
        daf::phaseMark("STV3_LEAN_peel_init", bytesBucket + bytesVertexAux + bytesLeafAux);
    }

    std::vector<uint8_t> dirtyMark(numVertices, 0);
    std::vector<daf::Size> dirtyVertices;
    dirtyVertices.reserve(8192);

    daf::StaticVector<daf::Size> currentRemoveVertexIds(numVertices);

    std::cout << "=========================begin (Lean)=========================" << std::endl;
    std::cout << "vertices in heap: " << remainingInHeap << std::endl;

    double minCore = 0;

    while (remainingInHeap > 0) {
        while (!buckets.empty() && buckets.begin()->second.empty()) {
            buckets.erase(buckets.begin());
        }
        if (buckets.empty()) break;

        int64_t curKey = buckets.begin()->first;
        minCore = std::max((double)curKey, minCore);

        while (!buckets.empty() && (double)buckets.begin()->first <= minCore) {
            auto &lst = buckets.begin()->second;
            while (!lst.empty()) {
                auto id = lst.back();
                lst.pop_back();
                vertexInHeap[id] = 0;
                currentRemoveVertexIds.push_back(id);
                coreV[id] = minCore;
                remainingInHeap--;
            }
            buckets.erase(buckets.begin());
        }

        if (remainingInHeap == 0) break;

        // Phase 1: walk vertex incidences, mark affected leaves.
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
                if (leafIsDead(leafId)) continue;
                if (!leafAffected[leafId]) {
                    leafAffected[leafId] = 1;
                    affectedLeaves.push_back(leafId);
                }
                if (STV3_getBit(d.vtxLeafIsPivot, ei)) {
                    leafRemovedPivots[leafId]++;
                } else {
                    leafKeepRemoved[leafId] = 1;
                }
            }
        }

        // Phase 2 (LEAN): for each affected leaf, recompute alive_pivot/total_keep
        // by scanning the leaf's vertex set, then apply delta in a second pass.
        for (auto leafId : affectedLeaves) {
            const daf::Size lBegin = d.leafVtxOff[leafId];
            const daf::Size lEnd = d.leafVtxOff[leafId + 1];

            // Pass A: count alive_pivot and total_keep.
            int alivePivot = 0;
            int totalKeep = 0;
            for (daf::Size li = lBegin; li < lEnd; ++li) {
                bool isPiv = STV3_getBit(d.leafVtxIsPivot, li);
                bool inHeap = vertexInHeap[d.leafVtxIds[li]];
                if (isPiv) {
                    if (inHeap) alivePivot++;
                } else {
                    totalKeep++;
                }
            }

            int new_rp = alivePivot;
            int old_rp = alivePivot + leafRemovedPivots[leafId];
            int np = (int)k - totalKeep;
            bool dies = leafKeepRemoved[leafId] || (np > new_rp) || (new_rp < 0);

            double deltaKeep, deltaPivot;
            if (dies) {
                deltaKeep = (np >= 0 && np <= old_rp) ? nCr[old_rp][np] : 0.0;
                deltaPivot = (np >= 1 && old_rp >= 1) ? nCr[old_rp - 1][np - 1] : 0.0;
            } else {
                deltaKeep = nCr[old_rp][np] - nCr[new_rp][np];
                deltaPivot = (np >= 1) ? nCr[old_rp - 1][np - 1] - nCr[new_rp - 1][np - 1] : 0.0;
            }

            // Pass B: apply delta to alive vertices.
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

            if (dies) markLeafDead(leafId);
        }

        for (auto v : dirtyVertices) {
            if (vertexInHeap[v]) bucketMove(v);
            dirtyMark[v] = 0;
        }
        dirtyVertices.clear();

        for (auto leafId : affectedLeaves) {
            leafRemovedPivots[leafId] = 0;
            leafKeepRemoved[leafId] = 0;
            leafAffected[leafId] = 0;
        }
        affectedLeaves.clear();
        currentRemoveVertexIds.clear();
    }

    std::cout << "ST_V3_Lean peeling time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    daf::phaseMark("STV3_LEAN_peel_loop");

    delete[] countingV;
    currentRemoveVertexIds.free();
    return coreV;
}

double * NCliqueVertexCoreDecomposition_ST_V3_Lean(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto data = NCliqueVertexCoreDecomposition_ST_V3_Lean_Build(edgeGraph, k);
    return NCliqueVertexCoreDecomposition_ST_V3_Lean_Peel(data, k);
}
