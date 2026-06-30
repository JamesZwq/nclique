//
// MaxCliqueLeafProfile.cpp
//
// Quick profiler: count per-vertex membership in "Type A" (true max
// clique) SDCT leaves vs "Type B" (depth-bounded at |keepV|==s).
//
// From SDCT_Augmented.inl bkRecurse_NoTree:
//   if (|keepV| == max_k) emit with empty dropV   ← Type B (NOT max clique)
//   else                  emit with actual dropV  ← Type A (max clique:
//                                                   BK natural terminate,
//                                                   H ∪ P is maximal)
//
// Gated by PIVOTER_PROFILE_MAX_CLIQUE=1.
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include <chrono>
#include <vector>
#include <map>
#include <iostream>

void profileMaxCliqueLeaves(Graph& g, daf::CliqueSize s) {
    auto t0 = std::chrono::high_resolution_clock::now();

    const daf::Size n = g.getGraphNodeSize();
    std::vector<int> numAllLeaves(n, 0);
    std::vector<int> numMaxCliqueLeaves(n, 0);

    long long totalLeaves = 0, typeALeaves = 0, typeBLeaves = 0;
    long long typeATotalIncidence = 0, typeBTotalIncidence = 0;

    SDCT_Augmented_NoTree(g, s, /*min_k=*/1,
        [&](daf::Size /*leafId*/,
            const daf::StaticVector<int>& keepV,
            const daf::StaticVector<int>& dropV)
        {
            ++totalLeaves;
            const int keepC  = (int)keepV.size();
            const int pivotC = (int)dropV.size();
            const bool isMaxClique = (keepC < (int)s);

            if (isMaxClique) {
                ++typeALeaves;
                typeATotalIncidence += keepC + pivotC;
            } else {
                ++typeBLeaves;
                typeBTotalIncidence += keepC + pivotC;
            }

            for (int i = 0; i < keepC; ++i) {
                ++numAllLeaves[keepV[i]];
                if (isMaxClique) ++numMaxCliqueLeaves[keepV[i]];
            }
            for (int i = 0; i < pivotC; ++i) {
                ++numAllLeaves[dropV[i]];
                if (isMaxClique) ++numMaxCliqueLeaves[dropV[i]];
            }
        });

    auto t1 = std::chrono::high_resolution_clock::now();
    long long walkMs = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    std::map<int, int> histAll, histMax;
    int aliveVertices = 0;
    for (daf::Size v = 0; v < n; ++v) {
        if (numAllLeaves[v] > 0) ++aliveVertices;
        histAll[numAllLeaves[v]]++;
        histMax[numMaxCliqueLeaves[v]]++;
    }

    std::cout << "================ MAX_CLIQUE_LEAF_PROFILE (s=" << (int)s << ") ================" << std::endl;
    std::cout << "  walk time: " << walkMs << " ms" << std::endl;
    std::cout << "  total leaves: " << totalLeaves << std::endl;
    std::cout << "    Type A (max clique, |keepV|<s): " << typeALeaves
              << "  (" << (100.0 * typeALeaves / std::max<long long>(1, totalLeaves)) << "%)"
              << "  incidence=" << typeATotalIncidence << std::endl;
    std::cout << "    Type B (depth-bound, |keepV|==s): " << typeBLeaves
              << "  (" << (100.0 * typeBLeaves / std::max<long long>(1, totalLeaves)) << "%)"
              << "  incidence=" << typeBTotalIncidence << std::endl;
    std::cout << "  alive vertices (in any leaf): " << aliveVertices << " / " << n << std::endl;
    std::cout << std::endl;

    std::cout << "  per-vertex #Type-A-leaves (max-clique-leaf) histogram:" << std::endl;
    std::cout << "    in 0 Type-A leaves: " << histMax[0]
              << "  (" << (100.0 * histMax[0] / n) << "% of all vertices)" << std::endl;
    long long inAtLeastOne = 0;
    for (auto& kv : histMax) if (kv.first > 0) inAtLeastOne += kv.second;
    std::cout << "    in >=1 Type-A leaves: " << inAtLeastOne
              << "  (" << (100.0 * inAtLeastOne / n) << "%)" << std::endl;
    std::cout << "    in exactly 1 (= singleton in max-clique sense): " << histMax[1]
              << "  (" << (100.0 * histMax[1] / n) << "% of all,  "
              << (100.0 * histMax[1] / std::max<long long>(1, inAtLeastOne)) << "% of alive)" << std::endl;
    std::cout << "    in exactly 2: " << histMax[2] << std::endl;
    std::cout << "    in exactly 3: " << histMax[3] << std::endl;
    std::cout << "    in exactly 4: " << histMax[4] << std::endl;
    std::cout << "    in exactly 5: " << histMax[5] << std::endl;
    std::cout << "    in [6, 10]: ";
    long long bucket = 0;
    for (auto& kv : histMax) if (kv.first >= 6 && kv.first <= 10) bucket += kv.second;
    std::cout << bucket << std::endl;
    std::cout << "    in [11, 100]: ";
    bucket = 0;
    for (auto& kv : histMax) if (kv.first >= 11 && kv.first <= 100) bucket += kv.second;
    std::cout << bucket << std::endl;
    std::cout << "    in >100: ";
    bucket = 0;
    for (auto& kv : histMax) if (kv.first > 100) bucket += kv.second;
    std::cout << bucket << std::endl;

    std::cout << std::endl;
    std::cout << "  per-vertex #ALL-leaves (= existing vtxLeafOff) histogram:" << std::endl;
    std::cout << "    in 0 leaves: " << histAll[0] << std::endl;
    long long inAllAtLeastOne = 0;
    for (auto& kv : histAll) if (kv.first > 0) inAllAtLeastOne += kv.second;
    std::cout << "    in >=1 leaves: " << inAllAtLeastOne << std::endl;
    std::cout << "    in exactly 1 leaf (singleton, ANY-leaf sense): " << histAll[1] << std::endl;
    std::cout << "==========================================================================" << std::endl;
}
