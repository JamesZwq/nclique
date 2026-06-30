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

double* NCliqueVertexCoreDecomposition_LazyPop(Graph& g, daf::CliqueSize s) {
    auto t0 = std::chrono::high_resolution_clock::now();

    auto d = NCliqueVertexCoreDecomposition_ST_V3_Build(g, s);
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
    std::cout << "LazyPop: V3 Build "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << " ms"
              << ", peel "
              << std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tPeel0).count() << " ms"
              << " (pop_refresh=" << popRefreshCnt
              << ", pop_bounce=" << popBounceCnt << ")" << std::endl;
    daf::phaseMark("LazyPop_peel");

    delete[] support;
    return coreV;
}
