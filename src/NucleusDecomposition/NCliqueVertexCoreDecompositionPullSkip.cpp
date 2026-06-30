//
// NCliqueVertexCoreDecompositionPullSkip.cpp
//
// R=1 (1,s)-nucleus decomposition with PULL + SKIP peel.
//
// V3 baseline: peel v -> CSR lookup affected leaves -> push delta to each
// vertex in those leaves (Phase 2 in V3).
//
// Here we keep V3 build (dual CSR + per-leaf state). Peel changes:
//   - precompute per-edge Δub = initial #s-cliques through (v,u) edge.
//   - per peel v: scan N(v) via our edge list; for each live u:
//       LB = support[u] - pendingUB[u] - Δub[v,u]
//       if LB > minCore + eps -> SKIP (pendingUB[u] += Δub[v,u])
//       else -> mark u for refresh
//   - refresh u: walk u's leaves via V3 dual CSR + current per-leaf state.
//   - on pop: if pendingUB[v] > 0, refresh first. If true < curKey, move
//     v to bucket toKey(true) and retry. Otherwise pop normally.
//
// Correctness invariant: true_support[u] in [support[u]-pendingUB, support[u]].
// Skip is safe because LB > minCore => true > minCore.
// On pop, refresh recovers true; bounce maintains correct peel order.
//
// Must be called BEFORE edgeGraph.beSingleEdge() (delegates to V3 Build
// which needs the original graph).
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

double* NCliqueVertexCoreDecomposition_PullSkip(Graph& g, daf::CliqueSize s) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // ====================================================================
    // Phase A: V3 build (SDCT + dual CSR + per-leaf metadata + countingV).
    // ====================================================================
    auto d = NCliqueVertexCoreDecomposition_ST_V3_Build(g, s);
    auto t1 = std::chrono::high_resolution_clock::now();

    // ====================================================================
    // Phase B: precompute per-edge Δub from leaves.
    //   leaf (V_h, V_p), |V_p|=p, needP=np:
    //     keep-keep contribution per pair: C(p, np)
    //     keep-pivot per pair:             C(p-1, np-1)
    //     pivot-pivot per pair:            C(p-2, np-2)
    // ====================================================================
    std::unordered_map<uint64_t, double> edgeMap;
    edgeMap.reserve(d.numLeaves * 4);

    auto encodeEdge = [](uint32_t a, uint32_t b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return ((uint64_t)a << 32) | (uint64_t)b;
    };

    for (size_t L = 0; L < d.numLeaves; ++L) {
        size_t lb = d.leafVtxOff[L];
        size_t le = d.leafVtxOff[L+1];
        if (le - lb < 2) continue;
        int p = d.leafPivotCount[L];
        int np = d.leafNeedPivot[L];
        if (np < 0) continue;
        double wKK = (np <= p) ? nCr[p][np] : 0.0;
        double wKP = (np >= 1 && p >= 1) ? nCr[p-1][np-1] : 0.0;
        double wPP = (np >= 2 && p >= 2) ? nCr[p-2][np-2] : 0.0;
        if (wKK == 0 && wKP == 0 && wPP == 0) continue;
        for (size_t i = lb; i < le; ++i) {
            uint32_t u = d.leafVtxIds[i];
            bool u_p = STV3_getBit(d.leafVtxIsPivot, i);
            for (size_t j = i + 1; j < le; ++j) {
                uint32_t v = d.leafVtxIds[j];
                bool v_p = STV3_getBit(d.leafVtxIsPivot, j);
                double w;
                if (!u_p && !v_p) w = wKK;
                else if (u_p && v_p) w = wPP;
                else                  w = wKP;
                if (w > 0.0) edgeMap[encodeEdge(u, v)] += w;
            }
        }
    }

    // Convert edgeMap -> per-vertex CSR (nbr+dub), each edge stored twice.
    std::vector<size_t> edgeOff(d.numVertices + 1, 0);
    for (auto& kv : edgeMap) {
        uint32_t a = (uint32_t)(kv.first >> 32);
        uint32_t b = (uint32_t)(kv.first & 0xffffffffULL);
        edgeOff[a + 1]++;
        edgeOff[b + 1]++;
    }
    for (daf::Size i = 1; i <= d.numVertices; ++i) edgeOff[i] += edgeOff[i-1];
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
    std::unordered_map<uint64_t, double>().swap(edgeMap);  // free

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "PullSkip: V3 Build "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms"
              << ", Δub precompute "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms"
              << ", unique edges = " << (totalE / 2) << std::endl;
    daf::phaseMark("PullSkip_dub",
                   (long)(edgeOff.capacity() * sizeof(size_t)
                        + edgeNbr.capacity() * sizeof(uint32_t)
                        + edgeDub.capacity() * sizeof(double)));

    // ====================================================================
    // Phase C: pull-skip peel.
    // ====================================================================
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

    std::vector<double> pendingUB(n, 0.0);

    // Refresh: recompute support[u] from current per-leaf state.
    auto refresh = [&](daf::Size u) -> double {
        double acc = 0.0;
        size_t b = d.vtxLeafOff[u];
        size_t e = d.vtxLeafOff[u+1];
        for (size_t k = b; k < e; ++k) {
            daf::Size L = d.vtxLeafIds[k];
            if (!leafAlive[L]) continue;
            int rp = leafRemainPivots[L];
            int np = d.leafNeedPivot[L];
            bool isPivot = STV3_getBit(d.vtxLeafIsPivot, k);
            if (!isPivot) {
                if (np >= 0 && np <= rp) acc += nCr[rp][np];
            } else {
                if (np >= 1 && rp >= 1) acc += nCr[rp-1][np-1];
            }
        }
        return acc;
    };

    auto* coreV = new double[n + 1];
    for (daf::Size i = 0; i <= n; ++i) coreV[i] = -1.0;

    // Sparse bucket queue (same pattern as V3).
    auto toKey = [](double x) -> int64_t {
        if (x <= 0.0) return 0;
        if (x >= 1e18) return (int64_t)1e18;
        return (int64_t)x;
    };

    std::map<int64_t, std::vector<daf::Size>> buckets;
    std::vector<int64_t>    bucket_of(n, -1);
    std::vector<daf::Size>  pos_in_bucket(n, 0);
    std::vector<uint8_t>    in_bucket(n, 0);
    daf::Size remaining = 0;
    for (daf::Size v = 0; v < n; ++v) {
        if (support[v] <= 0.0) { coreV[v] = 0.0; peeled[v] = 1; continue; }
        int64_t b = toKey(support[v]);
        auto& lst = buckets[b];
        bucket_of[v] = b;
        pos_in_bucket[v] = lst.size();
        lst.push_back(v);
        in_bucket[v] = 1;
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

    double minCore = 0.0;
    long long skipCnt = 0, refreshCnt = 0, popRefreshCnt = 0, popBounceCnt = 0;

    std::vector<int>       leafRemovedPivots(Lcnt, 0);
    std::vector<uint8_t>   leafDies(Lcnt, 0);
    std::vector<uint8_t>   leafAffected(Lcnt, 0);
    std::vector<daf::Size> affectedLeaves;
    affectedLeaves.reserve(4096);
    std::vector<uint8_t>   mustRefreshMark(n, 0);
    std::vector<daf::Size> mustRefreshSet;
    mustRefreshSet.reserve(4096);

    auto tPeel0 = std::chrono::high_resolution_clock::now();

    while (remaining > 0) {
        while (!buckets.empty() && buckets.begin()->second.empty()) {
            buckets.erase(buckets.begin());
        }
        if (buckets.empty()) break;

        int64_t curKey = buckets.begin()->first;
        daf::Size v = buckets.begin()->second.back();

        // Lazy refresh on pop: if v is stale, refresh and bounce if needed.
        if (pendingUB[v] > 0.0) {
            double newSup = refresh(v);
            ++popRefreshCnt;
            int64_t newKey = toKey(newSup);
            support[v] = newSup;
            pendingUB[v] = 0.0;
            if (newKey < curKey) {
                bucketMove(v);
                ++popBounceCnt;
                continue;  // retry: this bucket / vertex changed
            }
            // newKey == curKey (true equals stale): proceed to pop.
        }

        // Commit pop.
        buckets.begin()->second.pop_back();
        if (buckets.begin()->second.empty()) buckets.erase(buckets.begin());
        in_bucket[v] = 0;
        bucket_of[v] = -1;
        --remaining;

        minCore = std::max(minCore, (double)curKey);
        coreV[v] = minCore;
        peeled[v] = 1;

        if (remaining == 0) break;

        // Phase 1: mark v's leaves as affected, update kill/removedPivots.
        {
            size_t lb = d.vtxLeafOff[v];
            size_t le = d.vtxLeafOff[v+1];
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
        // Phase 2: commit per-leaf state changes (NO vertex propagation).
        for (daf::Size leafId : affectedLeaves) {
            int new_rp = leafRemainPivots[leafId] - leafRemovedPivots[leafId];
            int np = d.leafNeedPivot[leafId];
            bool dies = leafDies[leafId] || (np > new_rp) || (new_rp < 0);
            leafRemainPivots[leafId] = new_rp;
            if (dies) leafAlive[leafId] = 0;
        }

        // Phase 3: PULL via edge list. Skip if provably safe, else refresh.
        {
            size_t eb = edgeOff[v], ee = edgeOff[v+1];
            for (size_t k = eb; k < ee; ++k) {
                daf::Size u = edgeNbr[k];
                if (peeled[u]) continue;
                double dub = edgeDub[k];
                double lb_new = support[u] - pendingUB[u] - dub;
                if (lb_new > minCore + 1e-9) {
                    pendingUB[u] += dub;
                    ++skipCnt;
                } else {
                    if (!mustRefreshMark[u]) {
                        mustRefreshMark[u] = 1;
                        mustRefreshSet.push_back(u);
                    }
                }
            }
        }
        // Refresh must-refreshes (uses CURRENT per-leaf state from Phase 2).
        for (daf::Size u : mustRefreshSet) {
            double newSup = refresh(u);
            ++refreshCnt;
            support[u] = newSup;
            pendingUB[u] = 0.0;
            if (in_bucket[u]) bucketMove(u);
            mustRefreshMark[u] = 0;
        }
        mustRefreshSet.clear();

        // Cleanup affected leaves bookkeeping.
        for (daf::Size leafId : affectedLeaves) {
            leafRemovedPivots[leafId] = 0;
            leafDies[leafId] = 0;
            leafAffected[leafId] = 0;
        }
        affectedLeaves.clear();
    }

    auto tEnd = std::chrono::high_resolution_clock::now();
    long long total_events = skipCnt + refreshCnt;
    double skipRate = total_events > 0 ? (double)skipCnt / total_events : 0.0;
    std::cout << "PullSkip: peel "
              << std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tPeel0).count() << " ms"
              << " (skips=" << skipCnt
              << ", refreshes=" << refreshCnt
              << ", pop_refresh=" << popRefreshCnt
              << ", pop_bounce=" << popBounceCnt
              << ", skip_rate=" << skipRate << ")"
              << std::endl;
    daf::phaseMark("PullSkip_peel");

    delete[] support;
    return coreV;
}
