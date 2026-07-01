//
// CompressionProfile.cpp
//
// Measures the per-leaf compression of the SDCT compact index to estimate
// how much an adaptive enumerate-vs-compress build could save.
//
// For each SDCT leaf (hold set H, pivot set P, η = s - |H|):
//   cliques(L)   = C(|P|, η)       # s-cliques the leaf encodes
//   cpi(L)       = |H| + |P|       # incidences stored by the compact index
//   enum(L)      = s * cliques(L)  # flat storage if we enumerated instead
//
// A leaf "compresses" iff cliques(L) > cpi(L). We report:
//   - global compression ratio  = Σ cliques / Σ cpi
//   - fraction of Σ (index size) in low-compression leaves (cliques ≤ cpi)
//   - Σ_leaves min(cpi, enum-as-storage) vs Σ cpi  (memory-side adaptive gain)
//
// Gated by PIVOTER_PROFILE_COMPRESSION=1.
//

#include "NCliqueCoreDecomposition.h"
#include "SDCT_Augmented.h"
#include <chrono>
#include <vector>
#include <cstdint>
#include <iostream>

extern double nCr[1001][401];

void profileCompression(Graph& g, daf::CliqueSize s) {
    auto t0 = std::chrono::high_resolution_clock::now();

    long long numLeaves = 0;
    long double totalCliques = 0.0L;   // Σ C(p, η)
    long long   totalCpi = 0;          // Σ (h + p)
    long long   cpiInLowCompress = 0;  // Σ (h+p) over leaves with cliques ≤ cpi
    long long   lowCompressLeaves = 0;
    long double cliquesInLowCompress = 0.0L;
    // adaptive memory: Σ min(cpi_incidences, enum_incidences)
    // enum stores each clique's s vertices => s * cliques; but per leaf we
    // compare cpi=(h+p) vs a flat clique list = s * C(p,η).
    long double adaptiveMem = 0.0L;
    long double pureCpiMem = 0.0L;

    SDCT_Augmented_NoTree(g, s, /*min_k=*/1,
        [&](daf::Size /*leafId*/,
            const daf::StaticVector<int>& keepV,
            const daf::StaticVector<int>& dropV)
        {
            int h = (int)keepV.size();
            int p = (int)dropV.size();
            int eta = (int)s - h;
            if (eta < 0 || eta > p) return;
            double cl = nCr[p][eta];              // cliques encoded
            long long cpi = h + p;                // compact incidences
            double enumMem = (double)s * cl;      // flat clique storage

            ++numLeaves;
            totalCliques += (long double)cl;
            totalCpi += cpi;
            pureCpiMem += (long double)cpi;
            adaptiveMem += (long double)std::min((double)cpi, enumMem);
            if (cl <= (double)cpi) {
                ++lowCompressLeaves;
                cpiInLowCompress += cpi;
                cliquesInLowCompress += (long double)cl;
            }
        });

    auto t1 = std::chrono::high_resolution_clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    double compRatio = totalCpi > 0 ? (double)(totalCliques / totalCpi) : 0.0;
    std::cout << "================ COMPRESSION_PROFILE (s=" << (int)s << ") ================" << std::endl;
    std::cout << "  walk " << ms << " ms, leaves=" << numLeaves << std::endl;
    std::cout << "  #s-cliques (Σ C(p,η)) = " << (double)totalCliques << std::endl;
    std::cout << "  Σ (compact incidences) = " << totalCpi << std::endl;
    std::cout << "  GLOBAL compression ratio (#cliques/Σ) = " << compRatio << "x" << std::endl;
    std::cout << "  low-compression leaves (cliques ≤ incidences): "
              << lowCompressLeaves << " / " << numLeaves
              << "  (" << (100.0 * lowCompressLeaves / std::max<long long>(1, numLeaves)) << "%)" << std::endl;
    std::cout << "    their Σ share = " << cpiInLowCompress << " / " << totalCpi
              << "  (" << (100.0 * cpiInLowCompress / std::max<long long>(1, totalCpi)) << "% of index)" << std::endl;
    std::cout << "    their clique share = " << (double)cliquesInLowCompress
              << "  (" << (100.0 * (double)(cliquesInLowCompress / (totalCliques>0?totalCliques:1))) << "% of cliques)" << std::endl;
    std::cout << "  --- memory-side adaptive estimate ---" << std::endl;
    std::cout << "  pure-CPI incidences   = " << (double)pureCpiMem << std::endl;
    std::cout << "  adaptive Σmin(cpi,enum) = " << (double)adaptiveMem
              << "  (save " << (100.0 * (1.0 - (double)(adaptiveMem / (pureCpiMem>0?pureCpiMem:1)))) << "%)" << std::endl;
    std::cout << "==========================================================================" << std::endl;
}
