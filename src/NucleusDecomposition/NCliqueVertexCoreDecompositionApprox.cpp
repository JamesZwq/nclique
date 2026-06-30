//
// NCliqueVertexCoreDecompositionApprox.cpp
//
// R=1 (1,s)-nucleus: APPROXIMATE peel via Δub-blind apply.
//
// Strategy:
//   - Build: V3 dual CSR + per-leaf state + per-edge Δub precompute.
//   - Peel: per peel of v, for each u in N(v) (via edge list):
//       support[u] -= Δub[v,u]    // apply UB on Δ as if it were exact
//       bucketMove(u)             // update key with possibly-too-low value
//     NO refresh ever.
//   - coreV[u] = max(minCore, support[u]_at_pop)
//
// Correctness:
//   support[u] becomes a LOWER bound on true (Δub >= Δ_actual).
//   Cores are SYSTEMATICALLY UNDER-ASSIGNED vs V3.
//   Use only if approximate cores are acceptable.
//
// Must be called BEFORE edgeGraph.beSingleEdge() (V3 Build needs original).
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

double* NCliqueVertexCoreDecomposition_Approx(Graph& g, daf::CliqueSize s) {
    auto t0 = std::chrono::high_resolution_clock::now();

    auto d = NCliqueVertexCoreDecomposition_ST_V3_Build(g, s);
    auto t1 = std::chrono::high_resolution_clock::now();

    // ---- Per-edge Δub precompute (same as PullSkip) ----
    std::unordered_map<uint64_t, double> edgeMap;
    edgeMap.reserve(d.numLeaves * 4);
    auto encodeEdge = [](uint32_t a, uint32_t b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return ((uint64_t)a << 32) | (uint64_t)b;
    };
    for (size_t L = 0; L < d.numLeaves; ++L) {
        size_t lb = d.leafVtxOff[L];
        size_t le = d.leafVtxOff[L + 1];
        if (le - lb < 2) continue;
        int p = d.leafPivotCount[L];
        int np = d.leafNeedPivot[L];
        if (np < 0) continue;
        double wKK = (np <= p) ? nCr[p][np] : 0.0;
        double wKP = (np >= 1 && p >= 1) ? nCr[p - 1][np - 1] : 0.0;
        double wPP = (np >= 2 && p >= 2) ? nCr[p - 2][np - 2] : 0.0;
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

    std::vector<size_t> edgeOff(d.numVertices + 1, 0);
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
    std::cout << "Approx: V3 Build "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms"
              << ", Δub precompute "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms"
              << ", unique edges = " << (totalE / 2) << std::endl;

    // ---- Approximate peel ----
    const daf::Size n = d.numVertices;
    double* support = d.countingV;
    d.countingV = nullptr;
    std::vector<uint8_t> peeled(n, 0);

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

    double minCore = 0.0;
    long long propagateCnt = 0;

    auto tPeel0 = std::chrono::high_resolution_clock::now();

    while (remaining > 0) {
        while (!buckets.empty() && buckets.begin()->second.empty()) {
            buckets.erase(buckets.begin());
        }
        if (buckets.empty()) break;
        int64_t curKey = buckets.begin()->first;
        daf::Size v = buckets.begin()->second.back();
        buckets.begin()->second.pop_back();
        if (buckets.begin()->second.empty()) buckets.erase(buckets.begin());
        bucket_of[v] = -1;
        --remaining;

        minCore = std::max(minCore, (double)curKey);
        coreV[v] = minCore;
        peeled[v] = 1;

        if (remaining == 0) break;

        // Approximate propagation: apply Δub blindly to each live nbr.
        size_t eb = edgeOff[v], ee = edgeOff[v + 1];
        for (size_t k = eb; k < ee; ++k) {
            daf::Size u = edgeNbr[k];
            if (peeled[u]) continue;
            double dub = edgeDub[k];
            support[u] -= dub;
            if (support[u] < 0.0) support[u] = 0.0;
            bucketMove(u);
            ++propagateCnt;
        }
    }

    auto tEnd = std::chrono::high_resolution_clock::now();
    std::cout << "Approx: peel "
              << std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tPeel0).count() << " ms"
              << " (propagate_events=" << propagateCnt << ")" << std::endl;
    daf::phaseMark("Approx_peel");

    delete[] support;
    return coreV;
}
