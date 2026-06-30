//
// NCliqueVertexCoreDecompositionLean.cpp
//
// SPIN★-Lean: V→L CSR only (no L→V CSR). Vertex-pull peel.
//
// Memory: ~50% of SPIN★ (drops leafVtxOff/Ids/IsPivot).
// Time:   Θ(d̄ · Σ) vs SPIN★'s Θ(s · Σ). Asymptotically equivalent
//         in Σ; the leading constant is the average degree d̄
//         (Lean) vs the clique size s (SPIN★).
// Crossover: d̄ < s favors Lean, d̄ > s favors SPIN★.
//
// Peel structure:
//   Phase 1: for each peeled v, walk vtxLeafIds[v] (V→L CSR) to
//            mark affected leaves and update per-leaf counters
//            (leafRemainPivots, leafAlive). Compute deltaKeep/Pivot
//            for each affected leaf.
//   Phase 2: vertex-pull. For each peeled v in batch, walk graph
//            adjacency N(v). For each alive nbr u, mark u touched.
//            Per u (deduped), walk u's V→L CSR; for each leaf in
//            u that is affected, apply the leaf's delta (keep or
//            pivot variant) to support[u].
//
// Must be called BEFORE edgeGraph.beSingleEdge() (needs original
// graph adjacency during peel).
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

// ----------------------------------------------------------------------
// Lean build: V→L CSR + per-leaf static metadata + countingV ONLY.
// Does NOT build leafVtxOff/leafVtxIds/leafVtxIsPivot (L→V CSR).
// ----------------------------------------------------------------------
static ST_V2_Data lean_buildVL(Graph& edgeGraph, daf::CliqueSize k) {
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
    std::cout << "Lean: SDCT+callback took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t_sdct - t_start).count()
              << " ms, leaves=" << d.numLeaves
              << ", COO entries=" << cooBuf.size() << std::endl;
    daf::phaseMark("Lean_SDCT", (long)(cooBuf.capacity() * sizeof(COOEntry)));

    // Build V→L CSR ONLY.
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
    std::cout << "Lean: V→L CSR built in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t_csr - t_sdct).count()
              << " ms (L→V CSR skipped: ~50% memory saved on dual CSR)"
              << std::endl;
    daf::log_memory("Lean after V→L CSR");

    const long bytesVL = (long)(d.vtxLeafOff.capacity() * sizeof(d.vtxLeafOff[0])
                              + d.vtxLeafIds.capacity() * sizeof(d.vtxLeafIds[0])
                              + d.vtxLeafIsPivot.capacity() * sizeof(d.vtxLeafIsPivot[0]));
    daf::phaseMark("Lean_CSR", bytesVL);

    return d;
}


double* NCliqueVertexCoreDecomposition_Lean(Graph& g, daf::CliqueSize s) {
    auto t0 = std::chrono::high_resolution_clock::now();

    auto d = lean_buildVL(g, s);
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

    // Phase 1/2 scratch.
    std::vector<int>       leafRemovedPivots(Lcnt, 0);
    std::vector<uint8_t>   leafDies(Lcnt, 0);
    std::vector<uint8_t>   leafAffected(Lcnt, 0);
    std::vector<double>    leafDeltaKeep(Lcnt, 0.0);
    std::vector<double>    leafDeltaPivot(Lcnt, 0.0);
    std::vector<daf::Size> affectedLeaves;
    affectedLeaves.reserve(4096);

    std::vector<daf::Size> batch;
    batch.reserve(1024);
    std::vector<uint8_t>   touched(n, 0);
    std::vector<daf::Size> touchedList;
    touchedList.reserve(4096);
    std::vector<uint8_t>   dirtyMark(n, 0);
    std::vector<daf::Size> dirtyVertices;
    dirtyVertices.reserve(4096);

    double minCore = 0.0;
    long long peelEvents = 0;
    long long checkOps   = 0;

    auto tPeel0 = std::chrono::high_resolution_clock::now();

    while (remaining > 0) {
        while (!buckets.empty() && buckets.begin()->second.empty()) {
            buckets.erase(buckets.begin());
        }
        if (buckets.empty()) break;

        int64_t curKey = buckets.begin()->first;
        minCore = std::max(minCore, (double)curKey);

        // Drain ALL buckets with key ≤ minCore into batch (V3-style cascade).
        batch.clear();
        while (!buckets.empty() && (double)buckets.begin()->first <= minCore + 1e-9) {
            auto& lst = buckets.begin()->second;
            while (!lst.empty()) {
                daf::Size v = lst.back();
                lst.pop_back();
                bucket_of[v] = -1;
                --remaining;
                coreV[v] = minCore;
                peeled[v] = 1;
                batch.push_back(v);
            }
            if (buckets.begin()->second.empty()) buckets.erase(buckets.begin());
        }
        if (remaining == 0) break;
        if (batch.empty()) continue;

        // Phase 1: mark affected leaves + collect counter changes.
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
        // Compute per-affected-leaf delta + commit state change.
        for (daf::Size leafId : affectedLeaves) {
            int old_rp = leafRemainPivots[leafId];
            int np     = d.leafNeedPivot[leafId];
            int new_rp = old_rp - leafRemovedPivots[leafId];
            bool dies  = leafDies[leafId] || (np > new_rp) || (new_rp < 0);

            double deltaKeep, deltaPivot;
            if (dies) {
                deltaKeep  = (np >= 0 && np <= old_rp) ? nCr[old_rp][np] : 0.0;
                deltaPivot = (np >= 1 && old_rp >= 1) ? nCr[old_rp - 1][np - 1] : 0.0;
            } else {
                deltaKeep  = nCr[old_rp][np] - nCr[new_rp][np];
                deltaPivot = (np >= 1) ? nCr[old_rp - 1][np - 1] - nCr[new_rp - 1][np - 1] : 0.0;
            }
            leafDeltaKeep[leafId]  = deltaKeep;
            leafDeltaPivot[leafId] = deltaPivot;
            leafRemainPivots[leafId] = new_rp;
            if (dies) leafAlive[leafId] = 0;
        }

        // Phase 2: vertex-pull. For each batch v, scan graph nbrs.
        // For each unique alive u, walk u's V→L; apply delta for
        // those of u's leaves that are affected.
        for (daf::Size v : batch) {
            const daf::Size adjStart = g.adj_list_offsets[v];
            const daf::Size adjEnd   = g.adj_list_offsets[v + 1];
            for (daf::Size i = adjStart; i < adjEnd; ++i) {
                daf::Size u = g.adj_list[i];
                if (peeled[u]) continue;
                if (touched[u]) continue;
                touched[u] = 1;
                touchedList.push_back(u);

                size_t uLb = d.vtxLeafOff[u];
                size_t uLe = d.vtxLeafOff[u + 1];
                for (size_t k = uLb; k < uLe; ++k) {
                    ++checkOps;
                    daf::Size leafId = d.vtxLeafIds[k];
                    if (!leafAffected[leafId]) continue;

                    double delta = STV3_getBit(d.vtxLeafIsPivot, k)
                                 ? leafDeltaPivot[leafId]
                                 : leafDeltaKeep[leafId];
                    if (delta > 0.0) {
                        support[u] -= delta;
                        if (support[u] < 0.0) support[u] = 0.0;
                        if (!dirtyMark[u]) {
                            dirtyMark[u] = 1;
                            dirtyVertices.push_back(u);
                        }
                        ++peelEvents;
                    }
                }
            }
        }
        // Batch bucket moves.
        for (daf::Size u : dirtyVertices) {
            if (bucket_of[u] != -1) bucketMove(u);
            dirtyMark[u] = 0;
        }
        dirtyVertices.clear();

        // Cleanup.
        for (daf::Size u : touchedList) touched[u] = 0;
        touchedList.clear();
        for (daf::Size leafId : affectedLeaves) {
            leafRemovedPivots[leafId] = 0;
            leafDies[leafId] = 0;
            leafAffected[leafId] = 0;
            leafDeltaKeep[leafId] = 0.0;
            leafDeltaPivot[leafId] = 0.0;
        }
        affectedLeaves.clear();
    }

    auto tEnd = std::chrono::high_resolution_clock::now();
    auto build_us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    auto peel_us  = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tPeel0).count();
    std::cout << "Lean: V→L Build " << (build_us / 1000) << " ms"
              << ", peel " << (peel_us / 1000) << " ms"
              << " (peel_events=" << peelEvents
              << ", check_ops=" << checkOps << ")" << std::endl;
    std::cout << "LEAN_BUILD_US: " << build_us << std::endl;
    std::cout << "LEAN_PEEL_US: "  << peel_us  << std::endl;
    daf::phaseMark("Lean_peel");

    delete[] support;
    return coreV;
}
