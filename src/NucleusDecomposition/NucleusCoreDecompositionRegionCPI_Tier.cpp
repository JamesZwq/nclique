// ============================================================
// PIVOTER_TIER dispatcher: 4-tier RegND family ablation entry.
//
//   T1 (\regnd)        : tuple peel + direct s-clique support, NO CPI,
//                        NO dead-box IE, NO V-safe.
//   T2 (\regndplus)    : tuple peel + CPI AggrCount, full re-compute
//                        per peel (NO dead-box IE), NO V-safe.
//   T3 (\regndplusplus): tuple peel + CPI + dead-box IE update,
//                        NO V-safe / private cloud.
//   T4 (\regndstar)    : tuple peel + CPI + dead-box IE + V-safe
//                        direct-binning (headline RegNDC).
//
// Selection via:  PIVOTER_RUN_REGION_TIER=1 ./binary  (registered in
//                 degeneracy_cliques.cpp dispatch table)
//                 PIVOTER_TIER={1,2,3,4}   (default 4 for back-compat)
//
// T3 and T4 reuse the existing V3LM implementation by toggling the
// pre-existing env vars `PIVOTER_V3_NO_PRIVATE` and `PIVOTER_VSAFE_CLOUD`
// before invoking it.
//
// T1 and T2 dispatch into a recompute-peel variant guarded by
// PIVOTER_RECOMPUTE_PEEL (added to V3LM and V3LM_NoCPI as a Step-5/6
// substitute that re-runs Step-4 with rPeeled[] gating each batch).
// ============================================================

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "NCliqueCoreDecomposition.h"

namespace {

// setenv wrapper that does not stomp values the user explicitly set.
// We want PIVOTER_TIER to be the SOLE knob in ablation runs, so we set
// `overwrite=1` here -- a user mixing PIVOTER_TIER with PIVOTER_V3_NO_PRIVATE
// or PIVOTER_VSAFE_CLOUD by mistake would otherwise produce nonsense rows.
inline void forceSetenv(const char *name, const char *value) {
    setenv(name, value, 1);
}

inline void forceUnsetenv(const char *name) {
    unsetenv(name);
}

inline int decodeTier() {
    if (const char *t = std::getenv("PIVOTER_TIER")) {
        int v = std::atoi(t);
        if (v >= 1 && v <= 4) return v;
        std::cerr << "[Tier] WARN: PIVOTER_TIER='" << t
                  << "' out of range, defaulting to 4." << std::endl;
    }
    return 4;
}

inline const char *tierName(int tier) {
    switch (tier) {
        case 1: return "T1 / RegND      (tuple peel, direct s-enum, recompute)";
        case 2: return "T2 / RegND+     (tuple peel, CPI AggrCount, recompute)";
        case 3: return "T3 / RegND++    (tuple peel, CPI, dead-box IE)";
        case 4: return "T4 / RegND*     (T3 + V-safe direct-binning, headline)";
        default: return "T? / unknown";
    }
}

}  // namespace

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionCPI_Tier(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    const int tier = decodeTier();
    std::cout << "[Tier] selected: " << tierName(tier) << std::endl;

    // Reset all tier-related env vars to a known state.  We then set the
    // ones this tier requires.  V3LM and V3LM_NoCPI read these env vars
    // at the top of their functions, so the ordering "set env then call"
    // is correct.
    forceUnsetenv("PIVOTER_V3_NO_PRIVATE");
    forceUnsetenv("PIVOTER_VSAFE_CLOUD");
    forceUnsetenv("PIVOTER_RECOMPUTE_PEEL");

    switch (tier) {
        case 4: {
            // T4 = V3LM with V-safe direct-binning (= current RegNDC).
            // private cloud also enabled (default ON inside V3LM when
            // PIVOTER_V3_NO_PRIVATE is unset).
            forceSetenv("PIVOTER_VSAFE_CLOUD", "1");
            return NucleusCoreDecompositionRClique_RegionCPI_LowMem(
                tree, edgeGraph, treeGraphV, r, s, prebuiltIndex);
        }
        case 3: {
            // T3 = V3LM with NO private-cloud and NO V-safe.  Dead-box IE
            // peel still active.
            forceSetenv("PIVOTER_V3_NO_PRIVATE", "1");
            return NucleusCoreDecompositionRClique_RegionCPI_LowMem(
                tree, edgeGraph, treeGraphV, r, s, prebuiltIndex);
        }
        case 2: {
            // T2 = V3LM with CPI initial-support + RECOMPUTE peel
            // (skip Step 5+6 dead-box machinery, re-run Step 4 with
            // rPeeled[] gating after every batch).  No private cloud,
            // no V-safe.
            forceSetenv("PIVOTER_V3_NO_PRIVATE", "1");
            forceSetenv("PIVOTER_RECOMPUTE_PEEL", "1");
            return NucleusCoreDecompositionRClique_RegionCPI_LowMem(
                tree, edgeGraph, treeGraphV, r, s, prebuiltIndex);
        }
        case 1:
        default: {
            // T1 = V3LM_NoCPI (direct s-clique enumeration for support)
            // + RECOMPUTE peel.  No private cloud, no V-safe.
            forceSetenv("PIVOTER_V3_NO_PRIVATE", "1");
            forceSetenv("PIVOTER_RECOMPUTE_PEEL", "1");
            return NucleusCoreDecompositionRClique_RegionCPI_LowMem_NoCPI(
                tree, edgeGraph, treeGraphV, r, s, prebuiltIndex);
        }
    }
}
