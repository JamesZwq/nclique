//
// NCliqueVertexCoreDecompositionOnDemand.cpp
//
// Minimal-storage R=1 nucleus decomposition.
//
// Same CSR approach as ST_V2 but eliminates per-leaf metadata vectors
// (leafPivotCount, leafNeedPivot, leafAlive, leafRemainPivots).
// Instead, each leaf's original keepCount/pivotCount is stored compactly,
// and pivot removals are tracked cumulatively via a single counter per leaf.
// Leaf alive status is computed on-the-fly.
//
// Storage savings vs ST_V2:
//   - No leafAlive vector (O(numLeaves))
//   - No leafNeedPivot vector (O(numLeaves))
//   - leafPivotCount and cumulative removed pivots merged into one vector
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstring>

extern double nCr[1001][401];

// =====================================================================
//  Main entry point
//  Must be called BEFORE edgeGraph.beSingleEdge() (needs original graph).
// =====================================================================

double* NCliqueVertexCoreDecomposition_OnDemand(
    Graph& edgeGraph, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numVertices = edgeGraph.getGraphNodeSize();

    // =====================================================================
    //  Phase 0: Build dual CSR via SDCT callback + accumulate countingV.
    //  Store per-leaf: original keepCount and pivotCount (2 ints per leaf).
    // =====================================================================
    struct COOEntry {
        daf::Size vertex;
        daf::Size leafId;
        uint8_t isPivot;
    };
    std::vector<COOEntry> cooBuf;
    cooBuf.reserve(1 << 20);

    auto* countingV = new double[numVertices];
    memset(countingV, 0, numVertices * sizeof(double));

    // Per-leaf: store origKeepCount and origPivotCount
    std::vector<int> leafOrigKeepCount;
    std::vector<int> leafOrigPivotCount;
    leafOrigKeepCount.reserve(1 << 16);
    leafOrigPivotCount.reserve(1 << 16);

    size_t numLeaves = SDCT_Augmented_NoTree(edgeGraph, k, /*min_k=*/1,
        [&](daf::Size leafId, const daf::StaticVector<int>& keepV,
            const daf::StaticVector<int>& dropV) {
            int pivotC = (int)dropV.size();
            int keepC = (int)keepV.size();
            int needP = (int)k - keepC;

            // Ensure metadata vectors are large enough
            if (leafId >= leafOrigKeepCount.size()) {
                size_t newSz = std::max<size_t>(leafId + 1, leafOrigKeepCount.size() * 2);
                leafOrigKeepCount.resize(newSz, 0);
                leafOrigPivotCount.resize(newSz, 0);
            }
            leafOrigKeepCount[leafId] = keepC;
            leafOrigPivotCount[leafId] = pivotC;

            double wKeep = (needP >= 0 && needP <= pivotC) ? nCr[pivotC][needP] : 0.0;
            double wPivot = (needP >= 1 && pivotC >= 1) ? nCr[pivotC - 1][needP - 1] : 0.0;

            for (int i = 0; i < keepC; ++i) {
                daf::Size v = keepV[i];
                countingV[v] += wKeep;
                cooBuf.push_back({v, leafId, 0});
            }
            for (int i = 0; i < pivotC; ++i) {
                daf::Size v = dropV[i];
                countingV[v] += wPivot;
                cooBuf.push_back({v, leafId, 1});
            }
        });

    leafOrigKeepCount.resize(numLeaves);
    leafOrigPivotCount.resize(numLeaves);

    auto time_sdct = std::chrono::high_resolution_clock::now();

    // Build Vertex→Leaf CSR
    struct VLeafEntry { daf::Size leafId; uint8_t isPivot; };
    std::vector<daf::Size> vtxLeafOff(numVertices + 2, 0);
    for (auto& e : cooBuf)
        if (e.vertex < numVertices) vtxLeafOff[e.vertex + 1]++;
    for (daf::Size i = 1; i <= numVertices; ++i)
        vtxLeafOff[i] += vtxLeafOff[i - 1];
    std::vector<VLeafEntry> vtxLeafData(vtxLeafOff[numVertices]);
    {
        std::vector<daf::Size> pos(numVertices, 0);
        for (auto& e : cooBuf) {
            daf::Size v = e.vertex;
            if (v < numVertices) {
                daf::Size p = vtxLeafOff[v] + pos[v]++;
                vtxLeafData[p] = {e.leafId, e.isPivot};
            }
        }
    }

    // Build Leaf→Vertex CSR
    struct LeafVtxEntry { daf::Size vertex; uint8_t isPivot; };
    std::vector<daf::Size> leafVtxOff(numLeaves + 1, 0);
    for (auto& e : cooBuf)
        if (e.leafId < numLeaves) leafVtxOff[e.leafId + 1]++;
    for (size_t i = 1; i <= numLeaves; ++i)
        leafVtxOff[i] += leafVtxOff[i - 1];
    std::vector<LeafVtxEntry> leafVtxData(leafVtxOff[numLeaves]);
    {
        std::vector<daf::Size> pos(numLeaves, 0);
        for (auto& e : cooBuf) {
            daf::Size L = e.leafId;
            if (L < numLeaves) {
                daf::Size p = leafVtxOff[L] + pos[L]++;
                leafVtxData[p] = {e.vertex, e.isPivot};
            }
        }
    }

    cooBuf.clear();
    cooBuf.shrink_to_fit();

    auto time_csr = std::chrono::high_resolution_clock::now();
    std::cout << "OnDemand R1: SDCT+callback took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_sdct - time_start).count()
              << " ms, CSR built in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_csr - time_sdct).count()
              << " ms, leaves=" << numLeaves << std::endl;
    daf::log_memory("OnDemand R1 after CSR build");

    // =====================================================================
    //  Phase 1: Peeling — same proven logic as ST_V2 but with on-demand
    //  leaf state computation.
    //
    //  Per-leaf tracking:
    //    leafRemainPivots[L] — current number of alive pivots
    //    leafAlive[L] — 0 if any keep was removed or needPivot > remainPivots
    //  These replace ST_V2's identical vectors, computed from origKeep/origPivot.
    // =====================================================================

    auto coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;

    // Initialize per-leaf state from origKeep/origPivot
    std::vector<uint8_t> leafAlive(numLeaves);
    std::vector<int> leafRemainPivots(numLeaves);
    for (size_t L = 0; L < numLeaves; ++L) {
        leafRemainPivots[L] = leafOrigPivotCount[L];
        int np = (int)k - leafOrigKeepCount[L];
        leafAlive[L] = (np >= 0 && np <= leafOrigPivotCount[L]) ? 1 : 0;
    }

    // Compute needPivot per leaf (constant — doesn't change)
    std::vector<int> leafNeedPivot(numLeaves);
    for (size_t L = 0; L < numLeaves; ++L)
        leafNeedPivot[L] = (int)k - leafOrigKeepCount[L];

    // Free origKeep/origPivot — no longer needed
    leafOrigKeepCount.clear(); leafOrigKeepCount.shrink_to_fit();
    leafOrigPivotCount.clear(); leafOrigPivotCount.shrink_to_fit();

    // Bucket array
    int maxBucket = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] > 0)
            maxBucket = std::max(maxBucket, (int)countingV[i]);
    }
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(numVertices);
    std::vector<daf::Size> pos_in_bucket(numVertices);
    std::vector<uint8_t> vertexInHeap(numVertices, 0);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] <= 0) continue;
        int b = (int)countingV[i];
        bucket_of[i] = b;
        pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i);
        vertexInHeap[i] = 1;
        remainingInHeap++;
    }
    int curBucket = 0;

    auto bucketMove = [&](daf::Size id) {
        int newB = std::max(0, (int)countingV[id]);
        int oldB = bucket_of[id];
        if (newB == oldB) return;
        auto& oldVec = buckets[oldB];
        daf::Size myPos = pos_in_bucket[id];
        if (myPos + 1 < oldVec.size()) {
            daf::Size last = oldVec.back();
            oldVec[myPos] = last;
            pos_in_bucket[last] = myPos;
        }
        oldVec.pop_back();
        if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
        bucket_of[id] = newB;
        pos_in_bucket[id] = buckets[newB].size();
        buckets[newB].push_back(id);
        if (newB < curBucket) curBucket = newB;
    };

    // Per-leaf batch tracking (reused each iteration)
    std::vector<int> leafRemovedPivots(numLeaves, 0);
    std::vector<uint8_t> leafDies(numLeaves, 0);
    std::vector<uint8_t> leafAffected(numLeaves, 0);
    std::vector<daf::Size> affectedLeaves;
    affectedLeaves.reserve(4096);

    std::vector<uint8_t> dirtyMark(numVertices, 0);
    std::vector<daf::Size> dirtyVertices;
    dirtyVertices.reserve(8192);

    daf::StaticVector<daf::Size> currentRemoveVertexIds(numVertices);

    std::cout << "=========================begin=========================" << std::endl;
    std::cout << "vertices in heap: " << remainingInHeap << std::endl;

    double minCore = 0;

    while (remainingInHeap > 0) {
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((double)curBucket, minCore);

        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && (double)curBucket <= minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                vertexInHeap[id] = 0;
                currentRemoveVertexIds.push_back(id);
                coreV[id] = minCore;
                remainingInHeap--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (double)(curBucket + 1) <= minCore) {
                curBucket++;
            } else break;
        }

        if (remainingInHeap == 0) break;

        // Phase 1: find affected leaves via Vertex→Leaf CSR
        for (int vi = 0; vi < (int)currentRemoveVertexIds.size(); ++vi) {
            auto v = currentRemoveVertexIds[vi];
            const daf::Size begin = vtxLeafOff[v];
            const daf::Size end = vtxLeafOff[v + 1];
            if (vi + 1 < (int)currentRemoveVertexIds.size()) {
                auto nextV = currentRemoveVertexIds[vi + 1];
                __builtin_prefetch(&vtxLeafData[vtxLeafOff[nextV]], 0, 1);
            }
            for (daf::Size ei = begin; ei < end; ++ei) {
                const auto& entry = vtxLeafData[ei];
                daf::Size leafId = entry.leafId;
                if (!leafAlive[leafId]) continue;
                if (!leafAffected[leafId]) {
                    leafAffected[leafId] = 1;
                    affectedLeaves.push_back(leafId);
                }
                if (!entry.isPivot) {
                    leafDies[leafId] = 1;
                } else {
                    leafRemovedPivots[leafId]++;
                }
            }
        }

        // Phase 2: process each affected leaf via Leaf→Vertex CSR
        for (auto leafId : affectedLeaves) {
            int old_rp = leafRemainPivots[leafId];
            int np = leafNeedPivot[leafId];
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

            const daf::Size lBegin = leafVtxOff[leafId];
            const daf::Size lEnd = leafVtxOff[leafId + 1];
            for (daf::Size li = lBegin; li < lEnd; ++li) {
                auto& nd = leafVtxData[li];
                if (!vertexInHeap[nd.vertex]) continue;
                double delta = nd.isPivot ? deltaPivot : deltaKeep;
                if (delta > 0) {
                    countingV[nd.vertex] -= delta;
                    if (countingV[nd.vertex] < 0) countingV[nd.vertex] = 0;
                    if (!dirtyMark[nd.vertex]) {
                        dirtyMark[nd.vertex] = 1;
                        dirtyVertices.push_back(nd.vertex);
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

    auto time_end = std::chrono::high_resolution_clock::now();
    std::cout << "time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count()
              << " ms" << std::endl;
    daf::log_memory("OnDemand R1 final");

    delete[] countingV;
    currentRemoveVertexIds.free();
    return coreV;
}
