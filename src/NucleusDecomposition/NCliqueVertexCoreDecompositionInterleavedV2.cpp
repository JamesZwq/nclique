//
// NCliqueVertexCoreDecompositionInterleavedV2.cpp
//
// Optimized Interleaved Construction-Decomposition for R=1 nucleus decomposition.
//
// Improvements over Interleaved.cpp:
//   Opt 1: Fix drainHeap — advance minCore to curBucket (was breaking immediately)
//   Opt 2: Guard drain calls — only call drainHeap when new vertex might trigger peeling
//   Opt 3: Pre-reserve vtxLeafAdj based on vertex degree (reduces scattered realloc)
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstring>

extern double nCr[1001][401];

double * NCliqueVertexCoreDecomposition_InterleavedV2(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numVertices = edgeGraph.getGraphNodeSize();

    // --- Per-vertex support tracking ---
    auto *countingV = new double[numVertices];
    memset(countingV, 0, numVertices * sizeof(double));

    auto *pendingDelta = new double[numVertices];
    memset(pendingDelta, 0, numVertices * sizeof(double));

    auto *coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;

    // --- Dynamic adjacency lists ---
    struct VLeafEntry { daf::Size leafId; uint8_t isPivot; };
    struct LeafVtxEntry { daf::Size vertex; uint8_t isPivot; };

    std::vector<std::vector<VLeafEntry>> vtxLeafAdj(numVertices);
    std::vector<std::vector<LeafVtxEntry>> leafVtxAdj;
    leafVtxAdj.reserve(1 << 16);

    // Opt 3: Pre-reserve vtxLeafAdj based on vertex degree
    for (daf::Size v = 0; v < numVertices; ++v) {
        int deg = (int)edgeGraph.getNbrCount(v);
        vtxLeafAdj[v].reserve(std::max(1, deg));
    }

    // Per-leaf metadata
    std::vector<int> leafPivotCount;
    std::vector<int> leafNeedPivot;
    std::vector<int> leafRemainPivots;
    std::vector<uint8_t> leafAlive;

    leafPivotCount.reserve(1 << 16);
    leafNeedPivot.reserve(1 << 16);
    leafRemainPivots.reserve(1 << 16);
    leafAlive.reserve(1 << 16);

    // --- Peeling state ---
    std::vector<uint8_t> finalized(numVertices, 0);
    std::vector<uint8_t> vertexInHeap(numVertices, 0);
    std::vector<uint8_t> peeled(numVertices, 0);

    // Bucket array
    std::vector<std::vector<daf::Size>> buckets;
    std::vector<int> bucket_of(numVertices, 0);
    std::vector<daf::Size> pos_in_bucket(numVertices, 0);
    daf::Size remainingInHeap = 0;
    int curBucket = 0;

    auto ensureBucketSize = [&](int b) {
        if (b >= (int)buckets.size()) buckets.resize(b + 1);
    };

    auto bucketMove = [&](daf::Size id) {
        double eff = countingV[id] - pendingDelta[id];
        int newB = std::max(0, (int)eff);
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
        ensureBucketSize(newB);
        bucket_of[id] = newB;
        pos_in_bucket[id] = buckets[newB].size();
        buckets[newB].push_back(id);
        if (newB < curBucket) curBucket = newB;
    };

    // Per-leaf batch tracking
    std::vector<int> leafRemovedPivots;
    std::vector<uint8_t> leafDies;
    std::vector<uint8_t> leafAffected;
    std::vector<daf::Size> affectedLeaves;
    affectedLeaves.reserve(4096);

    std::vector<uint8_t> dirtyMark(numVertices, 0);
    std::vector<daf::Size> dirtyVertices;
    dirtyVertices.reserve(4096);

    double minCore = 0;

    // --- Drain peelable vertices from heap ---
    // Fixed Opt 1: During construction, only peel bucket 0 (zero effective support).
    // We cannot advance minCore during construction because not all vertices are
    // finalized — later leaves may increase a vertex's support. Only bucket-0
    // vertices are guaranteed to have core=0 regardless of future leaves.
    auto drainHeap = [&]() {
        while (remainingInHeap > 0) {
            // Only drain bucket 0 during construction
            if (curBucket > 0) {
                curBucket = 0;
            }
            if (buckets.empty() || buckets[0].empty()) break;

            // Collect all vertices at bucket 0
            std::vector<daf::Size> currentBatch;
            currentBatch.reserve(256);
            while (!buckets[0].empty()) {
                auto id = buckets[0].back();
                buckets[0].pop_back();
                vertexInHeap[id] = 0;
                peeled[id] = 1;
                currentBatch.push_back(id);
                coreV[id] = 0;  // zero support = core 0
                remainingInHeap--;
            }

            if (currentBatch.empty()) break;
            if (remainingInHeap == 0) break;

            // Phase 1: find affected leaves via vtxLeafAdj
            for (auto v : currentBatch) {
                for (auto &entry : vtxLeafAdj[v]) {
                    daf::Size leafId = entry.leafId;
                    if (leafId >= leafAlive.size() || !leafAlive[leafId]) continue;
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

            // Phase 2: process affected leaves
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

                for (auto &nd : leafVtxAdj[leafId]) {
                    daf::Size u = nd.vertex;
                    double delta = nd.isPivot ? deltaPivot : deltaKeep;
                    if (delta <= 0) continue;

                    if (finalized[u] && vertexInHeap[u]) {
                        countingV[u] -= delta;
                        if (countingV[u] - pendingDelta[u] < 0) {
                            countingV[u] = pendingDelta[u];
                        }
                        if (!dirtyMark[u]) {
                            dirtyMark[u] = 1;
                            dirtyVertices.push_back(u);
                        }
                    } else if (!finalized[u] && !peeled[u]) {
                        pendingDelta[u] += delta;
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
        }
    };

    // Track statistics
    size_t numPeeledDuring = 0;
    size_t numPeeledAfter = 0;

    using hrc = std::chrono::high_resolution_clock;
    long long t_onLeaf_ns = 0;
    long long t_onVertexDone_ns = 0;
    long long t_drain_ns = 0;
    size_t drainCallCount = 0;
    size_t drainPeelEvents = 0;
    size_t drainSkipCount = 0;  // Opt 2: how many drain calls were skipped
    size_t totalCOOEntries = 0;

    // --- SDCT with interleaved peeling ---
    size_t numLeaves = SDCT_Augmented_NoTree_Interleaved(edgeGraph, k, /*min_k=*/1,
        // onLeaf callback
        [&](daf::Size leafId, const daf::StaticVector<int> &keepV,
            const daf::StaticVector<int> &dropV)
        {
            auto t0 = hrc::now();
            int pivotC = (int)dropV.size();
            int keepC = (int)keepV.size();
            int needP = (int)k - keepC;

            // Grow per-leaf arrays if needed
            if (leafId >= leafVtxAdj.size()) {
                size_t newSz = std::max<size_t>(leafId + 1, leafVtxAdj.size() * 2);
                leafVtxAdj.resize(newSz);
                leafPivotCount.resize(newSz, 0);
                leafNeedPivot.resize(newSz, 0);
                leafRemainPivots.resize(newSz, 0);
                leafAlive.resize(newSz, 0);
                leafRemovedPivots.resize(newSz, 0);
                leafDies.resize(newSz, 0);
                leafAffected.resize(newSz, 0);
            }

            leafPivotCount[leafId] = pivotC;
            leafNeedPivot[leafId] = needP;
            leafRemainPivots[leafId] = pivotC;
            leafAlive[leafId] = (needP >= 0 && needP <= pivotC) ? 1 : 0;

            double wKeep = (needP >= 0 && needP <= pivotC) ? nCr[pivotC][needP] : 0.0;
            double wPivot = (needP >= 1 && pivotC >= 1) ? nCr[pivotC - 1][needP - 1] : 0.0;

            // Build dynamic adjacency
            leafVtxAdj[leafId].reserve(keepC + pivotC);
            for (int i = 0; i < keepC; ++i) {
                daf::Size v = keepV[i];
                countingV[v] += wKeep;
                vtxLeafAdj[v].push_back({leafId, 0});
                leafVtxAdj[leafId].push_back({v, 0});
            }
            for (int i = 0; i < pivotC; ++i) {
                daf::Size v = dropV[i];
                countingV[v] += wPivot;
                vtxLeafAdj[v].push_back({leafId, 1});
                leafVtxAdj[leafId].push_back({v, 1});
            }
            totalCOOEntries += keepC + pivotC;
            t_onLeaf_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - t0).count();
        },
        // onVertexDone callback
        [&](int vertex)
        {
            auto t0 = hrc::now();
            daf::Size v = (daf::Size)vertex;
            finalized[v] = 1;

            double effSupport = countingV[v] - pendingDelta[v];
            if (effSupport <= 0) {
                coreV[v] = minCore;
                peeled[v] = 1;
                numPeeledDuring++;
                return;
            }

            // Insert into bucket array
            int b = (int)effSupport;
            ensureBucketSize(b);
            bucket_of[v] = b;
            pos_in_bucket[v] = buckets[b].size();
            buckets[b].push_back(v);
            vertexInHeap[v] = 1;
            remainingInHeap++;

            // Opt 2: Only call drain when new vertex is at bucket 0
            // (drain only peels bucket 0 during construction)
            if (b == 0) {
                auto td0 = hrc::now();
                size_t heapBefore = remainingInHeap;
                drainHeap();
                drainPeelEvents += (heapBefore - remainingInHeap);
                t_drain_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - td0).count();
                drainCallCount++;
            } else {
                drainSkipCount++;
            }
            t_onVertexDone_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - t0).count();
        });

    // Trim per-leaf arrays
    leafVtxAdj.resize(numLeaves);
    leafPivotCount.resize(numLeaves);
    leafNeedPivot.resize(numLeaves);
    leafRemainPivots.resize(numLeaves);
    leafAlive.resize(numLeaves);
    leafRemovedPivots.resize(numLeaves);
    leafDies.resize(numLeaves);
    leafAffected.resize(numLeaves);

    auto time_sdct = std::chrono::high_resolution_clock::now();

    // --- Final peeling: drain remaining vertices in heap ---
    if (remainingInHeap > 0) {
        curBucket = 0;
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;

        while (remainingInHeap > 0) {
            while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
            if (curBucket >= (int)buckets.size()) break;

            minCore = std::max(minCore, (double)curBucket);

            std::vector<daf::Size> currentBatch;
            currentBatch.reserve(256);
            while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
                   && (double)curBucket <= minCore) {
                while (!buckets[curBucket].empty()) {
                    auto id = buckets[curBucket].back();
                    buckets[curBucket].pop_back();
                    vertexInHeap[id] = 0;
                    peeled[id] = 1;
                    currentBatch.push_back(id);
                    coreV[id] = minCore;
                    remainingInHeap--;
                    numPeeledAfter++;
                }
                if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                    && (double)(curBucket + 1) <= minCore) {
                    curBucket++;
                } else break;
            }

            if (currentBatch.empty() || remainingInHeap == 0) {
                if (remainingInHeap == 0) break;
                continue;
            }

            // Phase 1: find affected leaves
            for (auto v : currentBatch) {
                for (auto &entry : vtxLeafAdj[v]) {
                    daf::Size leafId = entry.leafId;
                    if (leafId >= leafAlive.size() || !leafAlive[leafId]) continue;
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

            // Phase 2: process affected leaves
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

                for (auto &nd : leafVtxAdj[leafId]) {
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

            for (auto leafId : affectedLeaves) {
                leafRemovedPivots[leafId] = 0;
                leafDies[leafId] = 0;
                leafAffected[leafId] = 0;
            }
            affectedLeaves.clear();
        }
    }

    auto time_end = std::chrono::high_resolution_clock::now();
    auto sdct_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_sdct - time_start).count();
    auto peel_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_sdct).count();
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();

    std::cout << "InterleavedV2 R=1: SDCT+inline_peel=" << sdct_ms << "ms, final_peel=" << peel_ms
              << "ms, total=" << total_ms << "ms" << std::endl;
    std::cout << "  leaves=" << numLeaves
              << ", peeled_during=" << numPeeledDuring
              << ", peeled_after=" << numPeeledAfter
              << " (" << (numPeeledDuring * 100.0 / std::max<size_t>(1, numPeeledDuring + numPeeledAfter))
              << "% during construction)" << std::endl;
    std::cout << "  [Breakdown] onLeaf=" << (t_onLeaf_ns / 1000000) << "ms"
              << ", onVertexDone=" << (t_onVertexDone_ns / 1000000) << "ms"
              << " (drain=" << (t_drain_ns / 1000000) << "ms)"
              << ", BK_pure=" << (sdct_ms - t_onLeaf_ns/1000000 - t_onVertexDone_ns/1000000) << "ms" << std::endl;
    std::cout << "  [Stats] COO_entries=" << totalCOOEntries
              << ", drainCalls=" << drainCallCount
              << ", drainSkipped=" << drainSkipCount
              << ", drainPeelEvents=" << drainPeelEvents << std::endl;
    std::cout << "time: " << total_ms << " ms" << std::endl;

    delete[] countingV;
    delete[] pendingDelta;
    return coreV;
}
