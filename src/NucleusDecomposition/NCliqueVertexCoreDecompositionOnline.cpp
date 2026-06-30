//
// NCliqueVertexCoreDecompositionOnline.cpp
//
// Online R=1 (1,s)-nucleus decomposition.
//
// Zero CSR storage: no vtxLeafIds / leafVtxIds, no per-leaf metadata.
// Memory: O(n) -- support[], peeled[], coreV[].
//
// Algorithm (levelwise peel via repeated SDCT walks):
//   sup[v] = #s-cliques containing v over the live subgraph
//   Repeat:
//     k = min sup[v] over live v
//     minCore = max(minCore, k)
//     pop all v with sup[v] <= minCore, set coreV[v] = minCore
//     if no live remains: done
//     reset sup[live] = 0; re-walk SDCT to recompute sup[live]
//
// Walk filter (correctness for the live subgraph):
//   skip a leaf if ANY V_h vertex is peeled (the leaf's s-cliques
//   all contain the dead held vertex and no longer exist).
//   count only live pivots in V_p; if needP > p_live skip.
//
// Time: O(K * Sigma_walk) where K = #distinct popped minCore values.
// Trades a large memory win for K SDCT walks (vs V3's single walk).
//
// Must be called BEFORE edgeGraph.beSingleEdge() (needs original graph).
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include "../PhaseLogger.h"
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <cstdlib>

extern double nCr[1001][401];

double* NCliqueVertexCoreDecomposition_Online(
    Graph& edgeGraph, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numVertices = edgeGraph.getGraphNodeSize();

    // O(n) state -- only data structures we ever allocate.
    std::vector<double> sup(numVertices, 0.0);
    std::vector<uint8_t> peeled(numVertices, 0);
    auto* coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;

    long long total_leaves_skipped = 0;
    long long total_leaves_kept = 0;

    auto walk = [&](std::vector<double>& s_out) {
        long long skipped = 0, kept = 0;
        SDCT_Augmented_NoTree(edgeGraph, k, /*min_k=*/1,
            [&](daf::Size /*leafId*/,
                const daf::StaticVector<int>& keepV,
                const daf::StaticVector<int>& dropV) {
                const int keepC = (int)keepV.size();
                const int pivotC = (int)dropV.size();
                const int needP = (int)k - keepC;

                // Skip if any held vertex is peeled: all s-cliques from
                // this leaf contain it and are gone.
                for (int i = 0; i < keepC; ++i) {
                    if (peeled[keepV[i]]) { ++skipped; return; }
                }

                // Count live pivots.
                int p_live = 0;
                for (int i = 0; i < pivotC; ++i) {
                    if (!peeled[dropV[i]]) ++p_live;
                }

                if (needP < 0 || needP > p_live) { ++skipped; return; }

                ++kept;
                const double wKeep  = nCr[p_live][needP];
                const double wPivot = (needP >= 1) ? nCr[p_live - 1][needP - 1] : 0.0;

                for (int i = 0; i < keepC; ++i) {
                    s_out[keepV[i]] += wKeep;
                }
                if (wPivot > 0.0) {
                    for (int i = 0; i < pivotC; ++i) {
                        if (!peeled[dropV[i]]) s_out[dropV[i]] += wPivot;
                    }
                }
            });
        total_leaves_skipped += skipped;
        total_leaves_kept    += kept;
    };

    // Phase 0: initial walk.
    walk(sup);

    auto time_first = std::chrono::high_resolution_clock::now();
    std::cout << "ONLINE: initial SDCT walk took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_first - time_start).count()
              << " ms" << std::endl;
    daf::phaseMark("ONLINE_initial_walk",
                   (long)(sup.capacity() * sizeof(double)
                        + peeled.capacity() * sizeof(uint8_t)
                        + (numVertices + 1) * sizeof(double)));

    // Phase 1: levelwise peel.
    double minCore = 0.0;
    daf::Size alive = numVertices;
    int walks = 1;
    long long total_pops = 0;

    // Initial alive count: only vertices with positive sup participate.
    // Vertices with sup = 0 get coreV = 0 immediately.
    {
        for (daf::Size v = 0; v < numVertices; ++v) {
            if (sup[v] <= 0.0) {
                peeled[v] = 1;
                coreV[v] = 0.0;
                --alive;
            }
        }
    }

    auto t_peel0 = std::chrono::high_resolution_clock::now();

    while (alive > 0) {
        // Find min sup among live.
        double mn = 1e300;
        for (daf::Size v = 0; v < numVertices; ++v) {
            if (!peeled[v] && sup[v] < mn) mn = sup[v];
        }
        if (mn >= 1e299) break;

        minCore = std::max(minCore, mn);

        // Pop all live vertices with sup <= minCore.
        daf::Size popped = 0;
        for (daf::Size v = 0; v < numVertices; ++v) {
            if (!peeled[v] && sup[v] <= minCore + 1e-9) {
                peeled[v] = 1;
                coreV[v] = minCore;
                --alive;
                ++popped;
            }
        }
        total_pops += popped;
        if (alive == 0) break;

        // Refresh sup for live vertices via a fresh SDCT walk.
        for (daf::Size v = 0; v < numVertices; ++v) {
            if (!peeled[v]) sup[v] = 0.0;
        }
        walk(sup);
        ++walks;

        if (walks % 20 == 0) {
            auto t_now = std::chrono::high_resolution_clock::now();
            long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_now - t_peel0).count();
            std::cout << "  ONLINE: walks=" << walks
                      << " minCore=" << minCore
                      << " alive=" << alive
                      << " peel_ms=" << ms << std::endl;
        }
    }

    auto time_end = std::chrono::high_resolution_clock::now();
    long long total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    long long peel_ms  = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - t_peel0).count();
    std::cout << "ONLINE: walks=" << walks
              << " total_pops=" << total_pops
              << " kept_leaves(sum)=" << total_leaves_kept
              << " skipped_leaves(sum)=" << total_leaves_skipped
              << " total=" << total_ms << " ms"
              << " peel=" << peel_ms << " ms" << std::endl;
    daf::phaseMark("ONLINE_peel_loop");

    return coreV;
}
