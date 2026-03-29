//
// NCliqueVertexCoreDecompositionInterleaved.cpp
//
// Interleaved Construction-Decomposition for R=1 nucleus decomposition.
//
// Key innovation: vertices are peeled DURING SDCT construction, not after.
// After each top-level vertex v's subtree completes in the degeneracy-ordered
// BK, countingV[v] is finalized (v will never appear in any future leaf).
// If effectiveSupport[v] = countingV[v] - pendingDelta[v] <= currentMinCore,
// v can be peeled immediately.
//
// Two-tier support tracking:
//   countingV[u]     — accumulated from SDCT callback (only increases)
//   pendingDelta[u]  — accumulated from peeling of finalized neighbors (only increases)
//   effectiveSupport = countingV[u] - pendingDelta[u]
//
// Only peel u when: (1) u is finalized AND (2) effectiveSupport[u] <= minCore
//
// Uses vector-of-vectors (dynamic adjacency lists) built incrementally:
//   vtxLeafAdj[v]  — vertex v's list of (leafId, isPivot) entries
//   leafVtxAdj[L]  — leaf L's list of (vertex, isPivot) entries
//
// Must be called BEFORE edgeGraph.beSingleEdge() (needs original graph for SDCT).
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstring>
#include <numeric>

extern double nCr[1001][401];

double * NCliqueVertexCoreDecomposition_Interleaved(
    Graph &edgeGraph, daf::CliqueSize k)
{
    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numVertices = edgeGraph.getGraphNodeSize();

    // --- Per-vertex support tracking ---
    auto *countingV = new double[numVertices];
    memset(countingV, 0, numVertices * sizeof(double));

    // pendingDelta[u] = support reduction from peeling of already-finalized neighbors
    auto *pendingDelta = new double[numVertices];
    memset(pendingDelta, 0, numVertices * sizeof(double));

    auto *coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;

    // --- Dynamic adjacency lists (built incrementally during onLeaf) ---
    struct VLeafEntry { daf::Size leafId; uint8_t isPivot; };
    struct LeafVtxEntry { daf::Size vertex; uint8_t isPivot; };

    std::vector<std::vector<VLeafEntry>> vtxLeafAdj(numVertices);
    std::vector<std::vector<LeafVtxEntry>> leafVtxAdj;
    leafVtxAdj.reserve(1 << 16);

    // Per-leaf metadata (grown dynamically)
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
    std::vector<uint8_t> peeled(numVertices, 0);  // already assigned core value

    // Bucket array — dynamically resized
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

    // Per-leaf batch tracking (reused across drain iterations)
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
    auto drainHeap = [&]() {
        while (remainingInHeap > 0) {
            while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
            if (curBucket >= (int)buckets.size()) break;
            if ((double)curBucket > minCore) break;  // can't peel above minCore

            // Collect vertices at current bucket level
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
                }
                if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                    && (double)(curBucket + 1) <= minCore) {
                    curBucket++;
                } else break;
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
                        // u is in heap — subtract and mark dirty
                        countingV[u] -= delta;
                        if (countingV[u] - pendingDelta[u] < 0) {
                            countingV[u] = pendingDelta[u]; // floor at 0 effective
                        }
                        if (!dirtyMark[u]) {
                            dirtyMark[u] = 1;
                            dirtyVertices.push_back(u);
                        }
                    } else if (!finalized[u] && !peeled[u]) {
                        // u not finalized yet — accumulate pending delta
                        pendingDelta[u] += delta;
                    }
                    // if peeled[u], ignore (already has core value)
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

    // Detailed timing accumulators (nanoseconds)
    using hrc = std::chrono::high_resolution_clock;
    long long t_onLeaf_ns = 0;        // time in onLeaf callback (support + adj build)
    long long t_onVertexDone_ns = 0;   // time in onVertexDone callback (bucket insert + drain)
    long long t_drain_ns = 0;          // time in drainHeap (subset of onVertexDone)
    size_t drainCallCount = 0;
    size_t drainPeelEvents = 0;        // how many vertices were actually peeled in drain
    size_t totalCOOEntries = 0;        // total vertex-leaf entries created

    // --- SDCT with interleaved peeling ---
    size_t numLeaves = SDCT_Augmented_NoTree_Interleaved(edgeGraph, k, /*min_k=*/1,
        // onLeaf callback: accumulate support + build dynamic adjacency lists
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
        // onVertexDone callback: vertex finalized, insert into heap and drain
        [&](int vertex)
        {
            auto t0 = hrc::now();
            daf::Size v = (daf::Size)vertex;
            finalized[v] = 1;

            double effSupport = countingV[v] - pendingDelta[v];
            if (effSupport <= 0) {
                coreV[v] = minCore;  // assign current minCore (may be 0)
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

            // Try to drain peelable vertices
            auto td0 = hrc::now();
            size_t heapBefore = remainingInHeap;
            drainHeap();
            drainPeelEvents += (heapBefore - remainingInHeap);
            t_drain_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - td0).count();
            drainCallCount++;
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
    // All vertices are finalized now, so we do standard peeling
    // Update minCore to allow full peeling
    if (remainingInHeap > 0) {
        // Find the minimum non-empty bucket
        curBucket = 0;
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;

        while (remainingInHeap > 0) {
            while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
            if (curBucket >= (int)buckets.size()) break;

            minCore = std::max(minCore, (double)curBucket);

            // Pop from bucket
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

    std::cout << "Interleaved R=1: SDCT+inline_peel=" << sdct_ms << "ms, final_peel=" << peel_ms
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
              << ", drainPeelEvents=" << drainPeelEvents << std::endl;
    std::cout << "time: " << total_ms << " ms" << std::endl;

    delete[] countingV;
    delete[] pendingDelta;
    return coreV;
}
