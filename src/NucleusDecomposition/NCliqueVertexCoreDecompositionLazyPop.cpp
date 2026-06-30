//
// NCliqueVertexCoreDecompositionLazyPop.cpp
//
// R=1 (1,s)-nucleus decomposition with LB-key (lower-bound) lazy peel.
//
// Algorithm:
//   Invariant: support[u] <= true_support[u] at all times (LB).
//     - Initially support[u] = exact (from SDCT build).
//     - On peel of v: for each edge (v, u), support[u] -= Δub[v,u]
//       (Δub = INITIAL #s-cliques through (v,u), >= Δ_actual). Bucket
//       key tracks support[u], so the bucket holds LBs.
//   Because min(LB) <= min(true), the front bucket is guaranteed to
//   contain a candidate for the global true min.
//   On pop u:
//     - Refresh u (walk u's leaves to compute current true support).
//     - If true == curKey: pop. coreV[u] = max(minCore, curKey).
//     - If true >  curKey: bounce UP to bucket toKey(true), retry.
//   This is provably correct: each pop assigns true min support at that
//   moment, matching V3 exactly.
//
// Compared to V3 push:
//   - peel update is per EDGE (deg(v)) instead of per (affected leaf,
//     vertex in leaf), avoiding the leaf-size-squared work in V3 Phase 2;
//   - refresh on pop is per VERTEX (walk u's V→L row), so each pop costs
//     O(|u's incidences|);
//   - bounce-up adds wasted refreshes when Δub is loose late in peel.
//
// Memory: Lean Build skips L→V CSR (only V→L kept). Plus O(m) for Δub.
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include "../PhaseLogger.h"
#include <chrono>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <cstring>
#include <iostream>

extern double nCr[1001][401];

// =====================================================================
// Lean build with inline per-edge Δub accumulation.
// Builds:
//   d.vtxLeafOff / vtxLeafIds / vtxLeafIsPivot  (V→L CSR)
//   d.leafPivotCount / leafNeedPivot            (per-leaf static)
//   d.countingV                                 (initial true support)
// Returns: ST_V2_Data (with L→V fields empty), plus edgeMap by reference.
// =====================================================================
static ST_V2_Data lazypop_leanBuild_withDub(
    Graph& edgeGraph, daf::CliqueSize k,
    std::unordered_map<uint64_t, double>& edgeMap)
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

    auto encodeEdge = [](uint32_t a, uint32_t b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return ((uint64_t)a << 32) | (uint64_t)b;
    };

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

            // -------- inline Δub accumulation per leaf pair --------
            if (needP < 0) return;
            double wKK = (needP <= pivotC) ? nCr[pivotC][needP] : 0.0;
            double wKP = (needP >= 1 && pivotC >= 1) ? nCr[pivotC - 1][needP - 1] : 0.0;
            double wPP = (needP >= 2 && pivotC >= 2) ? nCr[pivotC - 2][needP - 2] : 0.0;
            if (wKK > 0.0) {
                for (int i = 0; i < keepC; ++i)
                    for (int j = i + 1; j < keepC; ++j)
                        edgeMap[encodeEdge(keepV[i], keepV[j])] += wKK;
            }
            if (wKP > 0.0) {
                for (int i = 0; i < keepC; ++i)
                    for (int j = 0; j < pivotC; ++j)
                        edgeMap[encodeEdge(keepV[i], dropV[j])] += wKP;
            }
            if (wPP > 0.0) {
                for (int i = 0; i < pivotC; ++i)
                    for (int j = i + 1; j < pivotC; ++j)
                        edgeMap[encodeEdge(dropV[i], dropV[j])] += wPP;
            }
        });

    d.leafPivotCount.resize(d.numLeaves);
    d.leafNeedPivot.resize(d.numLeaves);

    auto t_sdct = std::chrono::high_resolution_clock::now();
    std::cout << "LazyPop-Lean: SDCT+callback took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t_sdct - t_start).count()
              << " ms, leaves=" << d.numLeaves
              << ", COO entries=" << cooBuf.size()
              << ", Δub edges=" << edgeMap.size() << std::endl;

    // ---- Build ONLY V→L CSR ----
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
              << " ms (L→V skipped)" << std::endl;
    daf::log_memory("LazyPop-Lean after V→L CSR");

    return d;
}


double* NCliqueVertexCoreDecomposition_LazyPop(Graph& g, daf::CliqueSize s) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // ============ Phase A: build (Lean + inline Δub) ============
    std::unordered_map<uint64_t, double> edgeMap;
    edgeMap.reserve(1 << 16);
    auto d = lazypop_leanBuild_withDub(g, s, edgeMap);
    auto t1 = std::chrono::high_resolution_clock::now();

    // ============ Phase B: convert edgeMap -> per-vertex CSR ============
    std::vector<size_t>   edgeOff(d.numVertices + 1, 0);
    for (auto& kv : edgeMap) {
        uint32_t a = (uint32_t)(kv.first >> 32);
        uint32_t b = (uint32_t)(kv.first & 0xffffffffULL);
        edgeOff[a + 1]++;
        edgeOff[b + 1]++;
    }
    for (daf::Size i = 1; i <= d.numVertices; ++i) edgeOff[i] += edgeOff[i - 1];
    size_t totalE = edgeOff[d.numVertices];
    std::vector<uint32_t> edgeNbr(totalE);
    std::vector<double>   edgeDub(totalE);
    {
        std::vector<size_t> pos(d.numVertices, 0);
        for (auto& kv : edgeMap) {
            uint32_t a = (uint32_t)(kv.first >> 32);
            uint32_t b = (uint32_t)(kv.first & 0xffffffffULL);
            size_t pa = edgeOff[a] + pos[a]++;
            size_t pb = edgeOff[b] + pos[b]++;
            edgeNbr[pa] = b; edgeDub[pa] = kv.second;
            edgeNbr[pb] = a; edgeDub[pb] = kv.second;
        }
    }
    std::unordered_map<uint64_t, double>().swap(edgeMap);

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "LazyPop: edge CSR built in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
              << " ms (unique edges = " << (totalE / 2) << ")" << std::endl;

    // ============ Phase C: LB-key lazy peel ============
    const daf::Size n = d.numVertices;
    const size_t Lcnt = d.numLeaves;

    double* support = d.countingV;
    d.countingV = nullptr;
    // Invariant: support[u] <= true_support[u]. Starts EXACT (LB tight).

    std::vector<uint8_t> peeled(n, 0);
    std::vector<uint8_t> leafAlive(Lcnt);
    std::vector<int>     leafRemainPivots(Lcnt);
    for (size_t i = 0; i < Lcnt; ++i) {
        leafRemainPivots[i] = d.leafPivotCount[i];
        int np = d.leafNeedPivot[i];
        leafAlive[i] = (np >= 0 && np <= d.leafPivotCount[i]) ? 1 : 0;
    }

    // refresh(u): true_support[u] given current per-leaf state.
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

    // Per-vertex flag: "support[u] is up-to-date with all Δub applied since
    // last refresh, but may still be a STRICT LB on true". Used to skip
    // redundant refresh when u has not been touched since the previous one.
    std::vector<int> lastRefreshBatch(n, 0);
    std::vector<int> lastTouchBatch(n, 0);
    int currentBatch = 0;

    double minCore = 0.0;
    long long popRefreshCnt = 0, bounceUpCnt = 0, peelEventCnt = 0;
    long long batchCount = 0;

    std::vector<int>       leafRemovedPivots(Lcnt, 0);
    std::vector<uint8_t>   leafDies(Lcnt, 0);
    std::vector<uint8_t>   leafAffected(Lcnt, 0);
    std::vector<daf::Size> affectedLeaves;
    affectedLeaves.reserve(4096);

    // Per-batch buffers.
    std::vector<daf::Size> batch;
    batch.reserve(1024);
    std::vector<uint8_t>   touched(n, 0);
    std::vector<daf::Size> touchedList;
    touchedList.reserve(4096);

    auto tPeel0 = std::chrono::high_resolution_clock::now();

    while (remaining > 0) {
        while (!buckets.empty() && buckets.begin()->second.empty()) {
            buckets.erase(buckets.begin());
        }
        if (buckets.empty()) break;

        int64_t curKey = buckets.begin()->first;
        double  minCore_now = std::max(minCore, (double)curKey);

        batch.clear();

        // ---- Drain entire front bucket (key == curKey) ----
        // Stale entries either pop with their true sup or bounce UP.
        while (remaining > 0 && !buckets.empty()
               && buckets.begin()->first == curKey) {
            auto& lst = buckets.begin()->second;
            if (lst.empty()) {
                buckets.erase(buckets.begin());
                continue;
            }
            daf::Size v = lst.back();

            if (lastTouchBatch[v] > lastRefreshBatch[v]) {
                double trueSup = refresh(v);
                ++popRefreshCnt;
                support[v] = trueSup;
                lastRefreshBatch[v] = currentBatch;
                int64_t newKey = toKey(trueSup);
                if (newKey > curKey) {
                    bucketMove(v);
                    ++bounceUpCnt;
                    continue;
                }
                // newKey == curKey: pop normally (LB was tight at int).
            }

            // Commit pop.
            lst.pop_back();
            if (lst.empty()) buckets.erase(buckets.begin());
            bucket_of[v] = -1;
            --remaining;
            coreV[v] = minCore_now;
            peeled[v] = 1;
            batch.push_back(v);
        }

        minCore = minCore_now;
        if (batch.empty()) continue;
        if (remaining == 0) break;

        // ---- Phase 1: per-leaf state update (deduped across batch) ----
        for (daf::Size v : batch) {
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
        for (daf::Size leafId : affectedLeaves) {
            int new_rp = leafRemainPivots[leafId] - leafRemovedPivots[leafId];
            int np = d.leafNeedPivot[leafId];
            bool dies = leafDies[leafId] || (np > new_rp) || (new_rp < 0);
            leafRemainPivots[leafId] = new_rp;
            if (dies) leafAlive[leafId] = 0;
        }
        for (daf::Size leafId : affectedLeaves) {
            leafRemovedPivots[leafId] = 0;
            leafDies[leafId] = 0;
            leafAffected[leafId] = 0;
        }
        affectedLeaves.clear();

        ++currentBatch;
        ++batchCount;

        // ---- Phase 2: eager Δub apply, dedup u via touched[] ----
        // Accumulate ALL drops to support[u] first, then ONE bucketMove
        // per touched u at the end. Avoids repeated O(log L) churn when u
        // is hit by multiple batch members.
        for (daf::Size v : batch) {
            size_t eb = edgeOff[v], ee = edgeOff[v + 1];
            for (size_t k = eb; k < ee; ++k) {
                daf::Size u = edgeNbr[k];
                if (peeled[u]) continue;
                double dub = edgeDub[k];
                if (dub <= 0.0) continue;
                support[u] -= dub;
                if (support[u] < 0.0) support[u] = 0.0;
                if (!touched[u]) {
                    touched[u] = 1;
                    touchedList.push_back(u);
                }
                lastTouchBatch[u] = currentBatch;
                ++peelEventCnt;
            }
        }
        for (daf::Size u : touchedList) {
            bucketMove(u);
            touched[u] = 0;
        }
        touchedList.clear();
    }

    auto tEnd = std::chrono::high_resolution_clock::now();
    auto build_us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    auto dub_us   = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    auto peel_us  = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tPeel0).count();
    std::cout << "LazyPop: Lean Build "
              << (build_us / 1000) << " ms"
              << ", Δub CSR " << (dub_us / 1000) << " ms"
              << ", peel " << (peel_us / 1000) << " ms"
              << " (batches=" << batchCount
              << ", peel_events=" << peelEventCnt
              << ", pop_refresh=" << popRefreshCnt
              << ", bounce_up=" << bounceUpCnt << ")"
              << std::endl;
    std::cout << "LAZYPOP_BUILD_US: " << build_us << std::endl;
    std::cout << "LAZYPOP_DUB_US: "   << dub_us   << std::endl;
    std::cout << "LAZYPOP_PEEL_US: "  << peel_us  << std::endl;
    daf::phaseMark("LazyPop_peel");

    delete[] support;
    return coreV;
}
