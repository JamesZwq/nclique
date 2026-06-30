//
// NCliqueVertexCoreDecompositionLazyPop.cpp
//
// R=1 (1,s)-nucleus: V3 build + LAZY refresh on pop. EXACT.
//
// Strategy:
//   - Build: V3 dual CSR + per-leaf state.
//   - Peel: update per-leaf state on each pop (Phase 1 + 2). NO vertex
//     propagation, NO Δub precompute.
//   - On pop: if v's support could be stale (any batch happened since its
//     last refresh), recompute support[v] from current per-leaf state.
//     If true < curKey, bucket-move down and retry (bounce). Else pop.
//
// Correctness: refresh on pop guarantees exact bucket key at pop time;
// minCore evolution matches V3.
//
// Must be called BEFORE edgeGraph.beSingleEdge() (V3 Build needs original).
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include "../PhaseLogger.h"
#include <chrono>
#include <vector>
#include <map>
#include <algorithm>
#include <cstring>
#include <iostream>

extern double nCr[1001][401];

// =====================================================================
// LazyPop Lean Build: ONLY V→L CSR (vtxLeafOff/Ids/IsPivot) +
// per-leaf metadata. Skips L→V CSR (leafVtx*) entirely; LazyPop peel
// never reads it. Saves ~half the dual-CSR memory at build peak.
// =====================================================================
static ST_V2_Data lazypop_leanBuild(Graph& edgeGraph, daf::CliqueSize k)
{
    auto t_start = std::chrono::high_resolution_clock::now();

    ST_V2_Data d;
    d.numVertices = edgeGraph.getGraphNodeSize();

    struct COOEntry {
        daf::Size vertex;
        daf::Size leafId;
        uint8_t   isPivot;
    };
    std::vector<COOEntry> cooBuf;
    cooBuf.reserve(1 << 20);

    d.countingV = new double[d.numVertices];
    std::memset(d.countingV, 0, d.numVertices * sizeof(double));

    d.numLeaves = SDCT_Augmented_NoTree(edgeGraph, k, /*min_k=*/1,
        [&](daf::Size leafId,
            const daf::StaticVector<int>& keepV,
            const daf::StaticVector<int>& dropV)
        {
            int pivotC = (int)dropV.size();
            int keepC  = (int)keepV.size();
            int needP  = (int)k - keepC;

            if (leafId >= d.leafPivotCount.size()) {
                size_t newSz = std::max<size_t>(leafId + 1,
                                                d.leafPivotCount.size() * 2);
                d.leafPivotCount.resize(newSz, 0);
                d.leafNeedPivot.resize(newSz, 0);
            }
            d.leafPivotCount[leafId] = pivotC;
            d.leafNeedPivot[leafId]  = needP;

            double wKeep  = (needP >= 0 && needP <= pivotC) ? nCr[pivotC][needP] : 0.0;
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

    auto t_sdct = std::chrono::high_resolution_clock::now();
    std::cout << "LazyPop-Lean: SDCT+callback took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t_sdct - t_start).count()
              << " ms, leaves=" << d.numLeaves
              << ", COO entries=" << cooBuf.size() << std::endl;
    daf::phaseMark("LazyPop_Lean_SDCT",
                   (long)(cooBuf.capacity() * sizeof(COOEntry)));

    // ---- Build ONLY V→L CSR (vtxLeafOff/Ids/IsPivot). L→V skipped. ----
    d.vtxLeafOff.assign(d.numVertices + 2, 0);
    for (auto& e : cooBuf)
        if (e.vertex < d.numVertices) d.vtxLeafOff[e.vertex + 1]++;
    for (daf::Size i = 1; i <= d.numVertices; ++i)
        d.vtxLeafOff[i] += d.vtxLeafOff[i - 1];

    const size_t totalIncidence = d.vtxLeafOff[d.numVertices];
    d.vtxLeafIds.resize(totalIncidence);
    d.vtxLeafIsPivot.assign((totalIncidence + 63) >> 6, 0);
    {
        std::vector<size_t> pos(d.numVertices, 0);
        for (auto& e : cooBuf) {
            daf::Size v = e.vertex;
            if (v < d.numVertices) {
                size_t p = d.vtxLeafOff[v] + pos[v]++;
                d.vtxLeafIds[p] = e.leafId;
                if (e.isPivot) STV3_setBit(d.vtxLeafIsPivot, p);
            }
        }
    }
    std::vector<COOEntry>().swap(cooBuf);

    auto t_csr = std::chrono::high_resolution_clock::now();
    std::cout << "LazyPop-Lean: V→L CSR built in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t_csr - t_sdct).count()
              << " ms (L→V skipped: saves ~50% dual-CSR)" << std::endl;
    daf::log_memory("LazyPop-Lean after V→L CSR");

    const long bytesVtxLeaf = (long)(d.vtxLeafOff.capacity() * sizeof(d.vtxLeafOff[0])
                            + d.vtxLeafIds.capacity() * sizeof(d.vtxLeafIds[0])
                            + d.vtxLeafIsPivot.capacity() * sizeof(d.vtxLeafIsPivot[0]));
    const long bytesSupport = (long)(d.numVertices * sizeof(double));
    const long bytesLeafMeta = (long)(d.leafPivotCount.capacity() * sizeof(int)
                            + d.leafNeedPivot.capacity() * sizeof(int));
    daf::phaseMark("LazyPop_Lean_CSR", bytesVtxLeaf + bytesSupport + bytesLeafMeta);

    return d;
}


double* NCliqueVertexCoreDecomposition_LazyPop(Graph& g, daf::CliqueSize s) {
    auto t0 = std::chrono::high_resolution_clock::now();

    auto d = lazypop_leanBuild(g, s);
    auto t1 = std::chrono::high_resolution_clock::now();

    const daf::Size n = d.numVertices;
    const size_t Lcnt = d.numLeaves;

    double* support = d.countingV;
    d.countingV = nullptr;

    std::vector<uint8_t> peeled(n, 0);
    std::vector<uint8_t> leafAlive(Lcnt);
    std::vector<int>     leafRemainPivots(Lcnt);
    for (size_t i = 0; i < Lcnt; ++i) {
        leafRemainPivots[i] = d.leafPivotCount[i];
        int np = d.leafNeedPivot[i];
        leafAlive[i] = (np >= 0 && np <= d.leafPivotCount[i]) ? 1 : 0;
    }

    auto refresh = [&](daf::Size u) -> double {
        double acc = 0.0;
        size_t b = d.vtxLeafOff[u];
        size_t e = d.vtxLeafOff[u + 1];
        for (size_t k = b; k < e; ++k) {
            daf::Size L = d.vtxLeafIds[k];
            if (!leafAlive[L]) continue;
            int rp = leafRemainPivots[L];
            int np = d.leafNeedPivot[L];
            bool isPivot = STV3_getBit(d.vtxLeafIsPivot, k);
            if (!isPivot) {
                if (np >= 0 && np <= rp) acc += nCr[rp][np];
            } else {
                if (np >= 1 && rp >= 1) acc += nCr[rp - 1][np - 1];
            }
        }
        return acc;
    };

    auto* coreV = new double[n + 1];
    for (daf::Size i = 0; i <= n; ++i) coreV[i] = -1.0;

    auto toKey = [](double x) -> int64_t {
        if (x <= 0.0) return 0;
        if (x >= 1e18) return (int64_t)1e18;
        return (int64_t)x;
    };

    std::map<int64_t, std::vector<daf::Size>> buckets;
    std::vector<int64_t>   bucket_of(n, -1);
    std::vector<daf::Size> pos_in_bucket(n, 0);
    daf::Size remaining = 0;
    for (daf::Size v = 0; v < n; ++v) {
        if (support[v] <= 0.0) { coreV[v] = 0.0; peeled[v] = 1; continue; }
        int64_t b = toKey(support[v]);
        auto& lst = buckets[b];
        bucket_of[v] = b;
        pos_in_bucket[v] = lst.size();
        lst.push_back(v);
        ++remaining;
    }

    auto bucketMove = [&](daf::Size v) {
        int64_t newB = toKey(support[v]);
        int64_t oldB = bucket_of[v];
        if (newB == oldB) return;
        auto it = buckets.find(oldB);
        auto& ov = it->second;
        daf::Size mp = pos_in_bucket[v];
        if (mp + 1 < ov.size()) {
            daf::Size last = ov.back();
            ov[mp] = last;
            pos_in_bucket[last] = mp;
        }
        ov.pop_back();
        if (ov.empty()) buckets.erase(it);
        bucket_of[v] = newB;
        auto& nv = buckets[newB];
        pos_in_bucket[v] = nv.size();
        nv.push_back(v);
    };

    // Staleness: vertex v is stale iff lastRefreshBatch[v] < currentBatch.
    // Initially support[v] is exact at batch 0, so set lastRefresh=0; bump
    // currentBatch after each pop so subsequent pops trigger refresh check.
    std::vector<int> lastRefreshBatch(n, 0);
    int currentBatch = 0;

    double minCore = 0.0;
    long long popRefreshCnt = 0, popBounceCnt = 0;

    std::vector<int>       leafRemovedPivots(Lcnt, 0);
    std::vector<uint8_t>   leafDies(Lcnt, 0);
    std::vector<uint8_t>   leafAffected(Lcnt, 0);
    std::vector<daf::Size> affectedLeaves;
    affectedLeaves.reserve(4096);

    auto tPeel0 = std::chrono::high_resolution_clock::now();

    while (remaining > 0) {
        while (!buckets.empty() && buckets.begin()->second.empty()) {
            buckets.erase(buckets.begin());
        }
        if (buckets.empty()) break;

        int64_t curKey = buckets.begin()->first;
        daf::Size v = buckets.begin()->second.back();

        if (lastRefreshBatch[v] < currentBatch) {
            double newSup = refresh(v);
            ++popRefreshCnt;
            support[v] = newSup;
            lastRefreshBatch[v] = currentBatch;
            int64_t newKey = toKey(newSup);
            if (newKey < curKey) {
                bucketMove(v);
                ++popBounceCnt;
                continue;
            }
        }

        // Commit pop.
        buckets.begin()->second.pop_back();
        if (buckets.begin()->second.empty()) buckets.erase(buckets.begin());
        bucket_of[v] = -1;
        --remaining;

        minCore = std::max(minCore, (double)curKey);
        coreV[v] = minCore;
        peeled[v] = 1;

        if (remaining == 0) break;

        // Phase 1: mark v's affected leaves.
        {
            size_t lb = d.vtxLeafOff[v];
            size_t le = d.vtxLeafOff[v + 1];
            for (size_t k = lb; k < le; ++k) {
                daf::Size leafId = d.vtxLeafIds[k];
                if (!leafAlive[leafId]) continue;
                if (!leafAffected[leafId]) {
                    leafAffected[leafId] = 1;
                    affectedLeaves.push_back(leafId);
                }
                if (!STV3_getBit(d.vtxLeafIsPivot, k)) leafDies[leafId] = 1;
                else                                    leafRemovedPivots[leafId]++;
            }
        }
        // Phase 2: commit per-leaf state changes.
        for (daf::Size leafId : affectedLeaves) {
            int new_rp = leafRemainPivots[leafId] - leafRemovedPivots[leafId];
            int np = d.leafNeedPivot[leafId];
            bool dies = leafDies[leafId] || (np > new_rp) || (new_rp < 0);
            leafRemainPivots[leafId] = new_rp;
            if (dies) leafAlive[leafId] = 0;
        }
        // Cleanup.
        for (daf::Size leafId : affectedLeaves) {
            leafRemovedPivots[leafId] = 0;
            leafDies[leafId] = 0;
            leafAffected[leafId] = 0;
        }
        affectedLeaves.clear();

        ++currentBatch;
    }

    auto tEnd = std::chrono::high_resolution_clock::now();
    auto build_us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    auto peel_us  = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tPeel0).count();
    std::cout << "LazyPop: Lean Build "
              << (build_us / 1000) << " ms"
              << ", peel "
              << (peel_us / 1000) << " ms"
              << " (pop_refresh=" << popRefreshCnt
              << ", pop_bounce=" << popBounceCnt << ")" << std::endl;
    // Microsecond-precision lines for the bench harness (small-peel cells).
    std::cout << "LAZYPOP_BUILD_US: " << build_us << std::endl;
    std::cout << "LAZYPOP_PEEL_US: "  << peel_us  << std::endl;
    daf::phaseMark("LazyPop_peel");

    delete[] support;
    return coreV;
}
