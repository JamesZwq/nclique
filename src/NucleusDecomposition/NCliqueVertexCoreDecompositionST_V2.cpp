//
// NCliqueVertexCoreDecompositionST_V2.cpp
//
// Tree-Free R=1 nucleus decomposition via Augmented SDCT.
//
// Key innovation: SDCT_Augmented_NoTree() builds dual CSR indices
// (Vertex→Leaf + Leaf→Vertex) inline during BK enumeration. The tree
// (DynamicGraph<TreeGraphNode>.adj_list) is NEVER materialised, saving
// ~53 → 16 bytes per leaf-vertex entry (3.3× memory reduction).
//
// Three redundant post-SDCT passes are eliminated:
//   1. treeGraphV construction  (vertex→leaf hash index)
//   2. countingPerVertex        (initial support scan)
//   3. CSR build pass           (flat vertex→leaf index)
// All three are fused into the single SDCT callback.
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstring>

extern double nCr[1001][401];

// =====================================================================
//  Phase 0: Build all data structures via Augmented SDCT callback.
//  MUST be called BEFORE edgeGraph.beSingleEdge() mutates the graph.
// =====================================================================

ST_V2_Data NCliqueVertexCoreDecomposition_ST_V2_Build(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    ST_V2_Data d;
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
    std::cout << "ST_V2: SDCT+callback took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_sdct - time_start).count()
              << " ms, leaves=" << d.numLeaves << ", COO entries=" << cooBuf.size() << std::endl;

    // --- Build dual CSR from COO ---
    d.vtxLeafOff.assign(d.numVertices + 2, 0);
    for (auto &e : cooBuf)
        if (e.vertex < d.numVertices) d.vtxLeafOff[e.vertex + 1]++;
    for (daf::Size i = 1; i <= d.numVertices; ++i)
        d.vtxLeafOff[i] += d.vtxLeafOff[i - 1];
    d.vtxLeafData.resize(d.vtxLeafOff[d.numVertices]);
    {
        std::vector<daf::Size> pos(d.numVertices, 0);
        for (auto &e : cooBuf) {
            daf::Size v = e.vertex;
            if (v < d.numVertices) {
                daf::Size p = d.vtxLeafOff[v] + pos[v]++;
                d.vtxLeafData[p] = {e.leafId, e.isPivot};
            }
        }
    }

    d.leafVtxOff.assign(d.numLeaves + 1, 0);
    for (auto &e : cooBuf)
        if (e.leafId < d.numLeaves) d.leafVtxOff[e.leafId + 1]++;
    for (size_t i = 1; i <= d.numLeaves; ++i)
        d.leafVtxOff[i] += d.leafVtxOff[i - 1];
    d.leafVtxData.resize(d.leafVtxOff[d.numLeaves]);
    {
        std::vector<daf::Size> pos(d.numLeaves, 0);
        for (auto &e : cooBuf) {
            daf::Size L = e.leafId;
            if (L < d.numLeaves) {
                daf::Size p = d.leafVtxOff[L] + pos[L]++;
                d.leafVtxData[p] = {e.vertex, e.isPivot};
            }
        }
    }

    auto time_csr = std::chrono::high_resolution_clock::now();
    std::cout << "ST_V2: dual CSR built in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_csr - time_sdct).count()
              << " ms" << std::endl;
    daf::log_memory("ST_V2 after CSR build");

    return d;
}

// =====================================================================
//  Phase 1: Peeling.  Graph-independent — uses only dual CSR + metadata.
// =====================================================================

double * NCliqueVertexCoreDecomposition_ST_V2_Peel(ST_V2_Data &d, daf::CliqueSize k)
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
        auto &oldVec = buckets[oldB];
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

    // Per-leaf batch tracking
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
            const daf::Size begin = d.vtxLeafOff[v];
            const daf::Size end = d.vtxLeafOff[v + 1];
            if (vi + 1 < (int)currentRemoveVertexIds.size()) {
                auto nextV = currentRemoveVertexIds[vi + 1];
                __builtin_prefetch(&d.vtxLeafData[d.vtxLeafOff[nextV]], 0, 1);
            }
            for (daf::Size ei = begin; ei < end; ++ei) {
                const auto &entry = d.vtxLeafData[ei];
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
                auto &nd = d.leafVtxData[li];
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

    std::cout << "ST_V2 peeling time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    delete[] countingV;
    currentRemoveVertexIds.free();
    return coreV;
}

// =====================================================================
//  Combined entry point (for use when graph is NOT yet mutated)
// =====================================================================
double * NCliqueVertexCoreDecomposition_ST_V2(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto data = NCliqueVertexCoreDecomposition_ST_V2_Build(edgeGraph, k);
    return NCliqueVertexCoreDecomposition_ST_V2_Peel(data, k);
}
