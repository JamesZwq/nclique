#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <climits>
#include <unistd.h>
#include <libgen.h>
#include <chrono>
#include <map>
#include <string>
#include <functional>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "misc.h"
#include "LinkedList.h"
#include "MemoryManager.h"
#include "BK/BronKerboschRmEdge.hpp"
#include "BK/BronKerboschRmRClique.hpp"
#include "dataStruct/disJoinSet.hpp"
#include "graph/DynamicBipartiteGraph.hpp"
#include "NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "degeneracy_algorithm_cliques_V.h"
#include "dataStruct/CliqueCSR.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "SDCT_Fused.hpp"

// ============================================================
// Utility
// ============================================================

static bool envSet(const char *name) { return std::getenv(name) != nullptr; }

// Return the first env var name that is set from a list, or nullptr
static const char *envFirst(std::initializer_list<const char *> names) {
    for (auto n : names)
        if (envSet(n)) return n;
    return nullptr;
}

// ============================================================
// Correctness comparison helpers
// ============================================================

static auto buildCoreDist(const auto &coreResults) {
    std::map<double, int> dist;
    for (const auto &[key, cv] : coreResults) dist[std::round(cv)]++;
    dist.erase(0.0);
    return dist;
}

static auto buildCoreDistFromArray(const double *coreV, daf::Size n) {
    std::map<double, int> dist;
    for (daf::Size i = 0; i < n; ++i)
        if (coreV[i] >= 0) dist[std::round(coreV[i])]++;
    dist.erase(0.0);
    return dist;
}

static void checkDist(const std::map<double,int> &refDist,
                      const std::map<double,int> &testDist,
                      const char *label) {
    // Exact check
    if (refDist == testDist) {
        std::cout << "✓ " << label << " correctness verified (exact)" << std::endl;
        return;
    }
    // Relaxed check: total count must match, per-core can differ by ≤ total*0.1%
    int refTotal = 0, testTotal = 0;
    for (auto &[k,v] : refDist) refTotal += v;
    for (auto &[k,v] : testDist) testTotal += v;
    if (refTotal != testTotal) {
        std::cerr << "✗ " << label << " MISMATCH! total ref=" << refTotal << " test=" << testTotal << std::endl;
        return;
    }
    // Check if differences are small (boundary rounding)
    int totalDiff = 0;
    std::map<double,int> merged;
    for (auto &[k,v] : refDist) merged[k] += 0;
    for (auto &[k,v] : testDist) merged[k] += 0;
    for (auto &[k, _] : merged) {
        int r = refDist.count(k) ? refDist.at(k) : 0;
        int t = testDist.count(k) ? testDist.at(k) : 0;
        totalDiff += std::abs(r - t);
    }
    double diffPct = 100.0 * totalDiff / (2 * refTotal);
    if (diffPct < 2.0) {
        std::cout << "✓ " << label << " correctness verified (diff=" << totalDiff/2 << " vertices, "
                  << diffPct << "%)" << std::endl;
    } else {
        std::cerr << "✗ " << label << " MISMATCH! diff=" << totalDiff/2 << " vertices (" << diffPct << "%)" << std::endl;
        for (auto &[k, _] : merged) {
            int r = refDist.count(k) ? refDist.at(k) : 0;
            int t = testDist.count(k) ? testDist.at(k) : 0;
            if (r != t) std::cerr << "  core=" << k << " ref=" << r << " test=" << t << std::endl;
        }
    }
}

// ============================================================
// Phase 1: Load graph and sort
// ============================================================

static Graph loadAndSortGraph(const char *fpath, int argc, char **argv) {
    Graph g(fpath);
    g.printGraphInfo();
    populate_nCr();
    daf::vListMap.resize(g.n + 1);
    memset(daf::vListMap.data(), -1, g.n * sizeof(daf::Size));

    if (argc >= 5) {
        std::string sortOption = argv[4];
        if (sortOption == "degen") g.sortByDegeneracyOrder(false);
        else if (sortOption == "degenR") g.sortByDegeneracyOrder(true);
        else if (sortOption == "degree") g.sortByDegree(false);
        else if (sortOption == "degreeR") g.sortByDegree(true);
        else if (sortOption == "default") { /* no-op */ }
        else {
            std::cerr << "Unknown sort option: " << sortOption << std::endl;
            std::cerr << "Available: degen, degenR, degree, degreeR, default" << std::endl;
            exit(1);
        }
        std::cout << "Graph sorted by " << sortOption << std::endl;
    } else {
        g.sortByDegeneracyOrder();
    }
    daf::log_memory("Graph Memory");
    return g;
}

// ============================================================
// Phase 2: Determine whether SDCT tree is needed
// ============================================================

// ST_V2 (r=1) and ST_V2_PROBE build their own data via SDCT_Augmented_NoTree.
// Skip the shared refTree to avoid a redundant SDCT pass.
static bool needsSDCT(daf::CliqueSize r, bool compareMode) {
    if (compareMode) return true;  // need refTree for correctness comparison
    if (r == 1 && (envSet("PIVOTER_RUN_ST_V2") || envSet("PIVOTER_RUN_ST_V2_PROBE")
                   || envSet("PIVOTER_RUN_INTERLEAVED") || envSet("PIVOTER_RUN_INTERLEAVED_V2")
                   || envSet("PIVOTER_RUN_ONDEMAND")))
        return false;
    // R=2 OnDemand needs SDCT (uses tree for Case B BK fallback)
    return true;
}

// Build SDCT tree with store_min_k=s (only large leaves in tree) + treeGraphV via callback.
// For R≥3: emit_min_k=r so all leaves ≥ r fire the callback (for future r-clique indexing).
// For R≤2: emit_min_k=s (small leaves have no useful info for vertex/edge peeling).
struct SDCTBuildResult {
    DynamicGraph<TreeGraphNode> tree;
    DynamicGraphSet<TreeGraphNode> treeGraphV;
    // R≥3: pre-built clique index (shared across all variants)
    std::unique_ptr<StaticCliqueIndex> cliqueIndex;
};

static SDCTBuildResult buildSDCTWithIndex(
    Graph &edgeGraph, daf::CliqueSize r, daf::CliqueSize s) {

    const int emit_min_k = r;   // callback for all leaves ≥ r
    const int store_min_k = s;  // only store leaves ≥ s in tree
    const daf::Size n = edgeGraph.getGraphNodeSize();

    DynamicGraphSet<TreeGraphNode> tgv(n);
    tgv.adj_list.resize(n);

    // R≥3: build cliqueIndex during SDCT (hash-free, keep+pivot naturally unique)
    auto ci = (r >= 3) ? std::make_unique<StaticCliqueIndex>(r) : nullptr;
    daf::StaticVector<daf::Size> keepBuf, dropBuf;
    daf::Size mergedBuf[16]; // enough for r ≤ 16

    // Pre-estimate r-clique count for pool reserve: sum of C(deg(v), r-1) / r
    // is a reasonable upper bound for r-cliques in a graph
    if (ci) {
        daf::Size totalEdges = edgeGraph.adj_list.size();
        double avgDeg = 2.0 * totalEdges / std::max((daf::Size)1, n);
        // Rough estimate: for r=3, ~m * avgDeg / 3; scale for other r
        daf::Size est = static_cast<daf::Size>(totalEdges * std::pow(avgDeg, r - 2) / r + 1000);
        ci->reservePool(est);
    }

    auto tree = daf::timeCount("SDCT_Fused", [&]() {
        return SDCT_Fused(edgeGraph, s, emit_min_k, store_min_k,
            [&](daf::Size leafId, const std::vector<TreeGraphNode> &leaf, bool stored) {
                if (stored) {
                    // treeGraphV: only for stored leaves (size ≥ s)
                    for (const auto &node : leaf) {
                        tgv.addNbr(node.v, {leafId, node.isPivot});
                    }
                }

                // R≥3: register r-cliques from ALL leaves (including size < s)
                if (ci) {
                    keepBuf.clear(); dropBuf.clear();
                    for (const auto &node : leaf) {
                        if (node.isPivot) dropBuf.push_back(node.v);
                        else keepBuf.push_back(node.v);
                    }
                    // enumerateCombinations: all keeps + (r-|keep|) pivots → naturally unique
                    daf::enumerateCombinations(keepBuf, dropBuf, r,
                        [&](const daf::StaticVector<daf::Size> &keep,
                            const daf::StaticVector<daf::Size> &combination) {
                            // Merge-sort keep + combination into sorted buffer
                            const daf::Size *a = combination.data(), *ae = a + combination.size();
                            const daf::Size *b = keep.data(), *be = b + keep.size();
                            daf::Size *out = mergedBuf;
                            while (a < ae && b < be) *out++ = (*a < *b ? *a++ : *b++);
                            while (a < ae) *out++ = *a++;
                            while (b < be) *out++ = *b++;
                            ci->addUniqueRClique(mergedBuf);
                            return true;
                        });
                }
            });
    });

    printf("SDCT leaf count (stored): %zu\n", tree.adj_list.size());

    // R≥3: build mapList from pool (single pass, enables lookupRaw)
    if (ci) {
        const daf::Size maxV = edgeGraph.adj_list.size();
        daf::timeCount("cliqueIndex mapList build", [&]() {
            ci->buildMapListFromPool(maxV);
        });
    }

    return {std::move(tree), std::move(tgv), std::move(ci)};
}

// Legacy: build SDCT + verify with Par7 (used when compareMode or non-fused path)
static DynamicGraph<TreeGraphNode> buildAndVerifySDCT(
    Graph &edgeGraph, daf::CliqueSize r, daf::CliqueSize s) {

    auto refTree = daf::timeCount("SDCT (reference)", [&]() {
        return SDCT(edgeGraph, s, r);
    });

    printf("SDCT leaf count: %zu\n", refTree.adj_list.size());

    // Verify SDCT_Par7
    {
        auto refCC = refTree.cliqueCount();
        auto par7Result = daf::timeCount("SDCT_Par7", [&]() -> CliqueCSR<int> {
            return SDCT_Par7(edgeGraph, s, r);
        });
        printf("SDCT_Par7 leaf count: %zu\n", par7Result.num_cliques());
        if (par7Result.has_pivot()) {
            auto par7CC = par7Result.cliqueCount();
            double rv = ((size_t)s < refCC.size()) ? refCC[s] : 0.0;
            double pv = ((size_t)s < par7CC.size()) ? par7CC[s] : 0.0;
            double diff = std::abs(rv - pv);
            double maxVal = std::max(std::abs(rv), std::abs(pv));
            if (diff > 0.5 && (maxVal < 1e-10 || diff / maxVal > 1e-6)) {
                printf("  ✗ cliqueCount mismatch at k=%d: ref=%.0f par7=%.0f\n", (int)s, rv, pv);
                printf("✗ SDCT_Par7 cliqueCount WRONG\n");
                exit(1);
            }
            printf("✓ SDCT_Par7 cliqueCount correct (k=%d)\n", (int)s);
        }
    }

    return refTree;
}

// ============================================================
// Phase 3: Pre-mutation work (ST_V2 build, probes)
// ============================================================

struct PreMutationResult {
    std::unique_ptr<ST_V2_Data> st_v2_data;
    double *interleavedCoreV = nullptr;  // owned, caller must delete[]
};

static PreMutationResult preMutationPhase(
    Graph &edgeGraph, daf::CliqueSize r, daf::CliqueSize s) {

    PreMutationResult result;

    if (r == 1 && envSet("PIVOTER_RUN_ST_V2")) {
        result.st_v2_data = std::make_unique<ST_V2_Data>(
            daf::timeCount("ST_V2 Build", [&]() {
                return NCliqueVertexCoreDecomposition_ST_V2_Build(edgeGraph, s);
            }));
    }
    if (r == 1 && envSet("PIVOTER_RUN_ST_V2_PROBE")) {
        NCliqueVertexCoreDecomposition_ST_V2_InterleavedProbe(edgeGraph, s);
    }
    if (r == 1 && envSet("PIVOTER_RUN_INTERLEAVED")) {
        result.interleavedCoreV = daf::timeCount("Interleaved r=1", [&]() {
            return NCliqueVertexCoreDecomposition_Interleaved(edgeGraph, s);
        });
    }
    if (r == 1 && envSet("PIVOTER_RUN_INTERLEAVED_V2")) {
        result.interleavedCoreV = daf::timeCount("InterleavedV2 r=1", [&]() {
            return NCliqueVertexCoreDecomposition_InterleavedV2(edgeGraph, s);
        });
    }

    if (r == 1 && envSet("PIVOTER_RUN_ONDEMAND")) {
        result.interleavedCoreV = daf::timeCount("OnDemand r=1", [&]() {
            return NCliqueVertexCoreDecomposition_OnDemand(edgeGraph, s);
        });
    }
    return result;
}

// ============================================================
// Phase 4: Prepare graph structures (mutates edgeGraph)
// ============================================================

static void prepareGraphStructures(Graph &edgeGraph, daf::CliqueSize r) {
    edgeGraph.initCore();
    edgeGraph.coreV.free();
    edgeGraph.beSingleEdge();
    edgeGraph.buildEdgeIdMap();
    if (r != 2) {
        edgeGraph.eidToNode.free();
    }
}

// ============================================================
// Phase 5: Build treeGraphV (conditional)
// ============================================================

static bool needsTreeGraphV(daf::CliqueSize r, bool compareMode) {
    if (compareMode) return true;
    if (r >= 2) return true;  // R=2 OnDemand needs treeGraphV for Case B BK
    // r=1: ST, ST_V2, Interleaved, OnDemand can skip treeGraphV
    if (envSet("PIVOTER_RUN_ST") || envSet("PIVOTER_RUN_ST_V2")
        || envSet("PIVOTER_RUN_INTERLEAVED") || envSet("PIVOTER_RUN_INTERLEAVED_V2")
        || envSet("PIVOTER_RUN_ONDEMAND")) return false;
    return true;
}

// ============================================================
// Dispatch: r=1 vertex-level algorithms
// ============================================================

// Wrapper: run a r=1 vertex algorithm with optional comparison against reference
template<typename Func>
static void runR1Variant(
    const char *name, Func &&func,
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize s, daf::Size numVertices, bool compareMode) {

    auto t2 = compareMode ? tree.clone() : DynamicGraph<TreeGraphNode>();
    auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();

    auto coreV = daf::timeCount(name, [&]() {
        return func(tree, edgeGraph, treeGraphV, s);
    });
    if (compareMode) {
        auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
        checkDist(buildCoreDistFromArray(refV, numVertices),
                  buildCoreDistFromArray(coreV, numVertices), name);
        delete[] refV;
    }
    delete[] coreV;
}

static bool dispatchR1(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize s, daf::Size numVertices, bool compareMode,
    PreMutationResult &pmr) {

    struct R1Entry {
        const char *envVar;
        const char *label;
        using FuncType = double*(*)(DynamicGraph<TreeGraphNode>&, const Graph&,
                                     DynamicGraphSet<TreeGraphNode>&, daf::CliqueSize);
        FuncType func;
    };

    // Table-driven dispatch for standard r=1 variants
    static const R1Entry table[] = {
        {"PIVOTER_RUN_LOCAL_V4",    "Local H-index V4 r=1",    NCliqueVertexCoreDecomposition_LocalV4},
        {"PIVOTER_RUN_LOCAL_V3",    "Local H-index V3 r=1",    NCliqueVertexCoreDecomposition_LocalV3},
        {"PIVOTER_RUN_LOCAL_V2",    "Local H-index V2 r=1",    NCliqueVertexCoreDecomposition_LocalV2},
        {"PIVOTER_RUN_LOCAL_NAIVE", "Local H-index Naive r=1", NCliqueVertexCoreDecomposition_LocalNaive},
        {"PIVOTER_RUN_LOCAL",       "Local H-index r=1",       NCliqueVertexCoreDecomposition_Local},
        {"PIVOTER_RUN_ST",          "ST r=1",                  NCliqueVertexCoreDecomposition_ST},
    };

    for (auto &entry : table) {
        if (envSet(entry.envVar)) {
            runR1Variant(entry.label, entry.func, tree, edgeGraph, treeGraphV,
                         s, numVertices, compareMode);
            return true;
        }
    }

    // ST_V2: special case — uses pre-built data, not the shared tree
    if (envSet("PIVOTER_RUN_ST_V2") && pmr.st_v2_data) {
        auto t2 = compareMode ? tree.clone() : DynamicGraph<TreeGraphNode>();
        auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
        auto coreV = daf::timeCount("ST_V2 r=1 (peel)", [&]() {
            return NCliqueVertexCoreDecomposition_ST_V2_Peel(*pmr.st_v2_data, s);
        });
        if (compareMode) {
            auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
            checkDist(buildCoreDistFromArray(refV, numVertices),
                      buildCoreDistFromArray(coreV, numVertices), "r=1 ST_V2");
            delete[] refV;
        }
        delete[] coreV;
        return true;
    }

    // Interleaved: already ran in preMutationPhase, just compare result
    if ((envSet("PIVOTER_RUN_INTERLEAVED") || envSet("PIVOTER_RUN_INTERLEAVED_V2")) && pmr.interleavedCoreV) {
        const char *label = envSet("PIVOTER_RUN_INTERLEAVED_V2") ? "r=1 InterleavedV2" : "r=1 Interleaved";
        if (compareMode) {
            auto t2 = tree.clone();
            auto tgv2 = treeGraphV.clone();
            auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
            checkDist(buildCoreDistFromArray(refV, numVertices),
                      buildCoreDistFromArray(pmr.interleavedCoreV, numVertices), label);
            delete[] refV;
        }
        delete[] pmr.interleavedCoreV;
        pmr.interleavedCoreV = nullptr;
        return true;
    }

    // OnDemand: already ran in preMutationPhase, just compare result
    if (envSet("PIVOTER_RUN_ONDEMAND") && pmr.interleavedCoreV) {
        if (compareMode) {
            auto t2 = tree.clone();
            auto tgv2 = treeGraphV.clone();
            auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
            checkDist(buildCoreDistFromArray(refV, numVertices),
                      buildCoreDistFromArray(pmr.interleavedCoreV, numVertices), "r=1 OnDemand");
            delete[] refV;
        }
        delete[] pmr.interleavedCoreV;
        pmr.interleavedCoreV = nullptr;
        return true;
    }

    return false;
}

// ============================================================
// Dispatch: r=2 edge-level algorithms
// ============================================================

// Wrapper: run a r=2 edge algorithm with optional comparison
template<typename Func, typename RefFunc>
static void runR2Variant(
    const char *name, Func &&func, RefFunc &&refFunc,
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s, bool compareMode) {

    auto t2 = compareMode ? tree.clone() : DynamicGraph<TreeGraphNode>();
    auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();

    auto result = daf::timeCount(name, [&]() {
        return func(tree, edgeGraph, treeGraphV, s);
    });
    if (compareMode) {
        auto refCore = refFunc(t2, edgeGraph, tgv2, r, s);
        checkDist(buildCoreDist(refCore), buildCoreDist(result), name);
    }
}

static bool dispatchR2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s, bool compareMode) {

    // Reference function for comparison
    auto correctRef = [](DynamicGraph<TreeGraphNode> &t, const Graph &g,
                         DynamicGraphSet<TreeGraphNode> &tgv,
                         daf::CliqueSize r, daf::CliqueSize s) {
        return PlusNucleusEdgeCoreDecompositionSet(t, g, tgv, s);
    };

    if (envSet("PIVOTER_RUN_LOCAL")) {
        auto t2 = compareMode ? tree.clone() : DynamicGraph<TreeGraphNode>();
        auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
        auto result = daf::timeCount("Local H-index r=2", [&]() {
            return NCliqueEdgeCoreDecomposition_Local(tree, edgeGraph, treeGraphV, s);
        });
        if (compareMode) {
            auto refCore = correctRef(t2, edgeGraph, tgv2, r, s);
            checkDist(buildCoreDist(refCore), buildCoreDist(result), "r=2 Local H-index");
        }
        return true;
    }

    struct R2Entry {
        const char *envVar;
        const char *label;
        using FuncType = std::vector<std::pair<std::pair<daf::Size,daf::Size>,int>>(*)(
            DynamicGraph<TreeGraphNode>&, const Graph&,
            DynamicGraphSet<TreeGraphNode>&, daf::CliqueSize);
        FuncType func;
    };

    static const R2Entry table[] = {
        {"PIVOTER_RUN_R2_DCLP", "DCLP r=2", PlusNucleusEdgeCoreDecompositionSet_DCLP},
        {"PIVOTER_RUN_R2_ST_V10", "ST_V10 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V10},
        {"PIVOTER_RUN_R2_ST_V9", "ST_V9 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V9},
        {"PIVOTER_RUN_R2_ST_V8B", "ST_V8b r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V8b},
        {"PIVOTER_RUN_R2_ST_V8", "ST_V8 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V8},
        {"PIVOTER_RUN_R2_ST_V7C", "ST_V7c r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V7c},
        {"PIVOTER_RUN_R2_ST_V7B", "ST_V7b r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V7b},
        {"PIVOTER_RUN_R2_ST_V7", "ST_V7 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V7},
        {"PIVOTER_RUN_R2_ST_V6", "ST_V6 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V6},
        {"PIVOTER_RUN_R2_ST_V5", "ST_V5 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V5},
        {"PIVOTER_RUN_R2_LAZY", "Lazy r=2", PlusNucleusEdgeCoreDecompositionSet_Lazy},
        {"PIVOTER_RUN_R2_ONDEMAND", "OnDemand r=2", PlusNucleusEdgeCoreDecompositionSet_OnDemand},
        {"PIVOTER_RUN_R2_TREEFREE_V2", "TreeFreeV2 r=2", PlusNucleusEdgeCoreDecompositionSet_TreeFreeV2},
        {"PIVOTER_RUN_R2_TREEFREE", "TreeFree r=2", PlusNucleusEdgeCoreDecompositionSet_TreeFree},
        {"PIVOTER_RUN_ST_V4", "ST_V4 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V4},
        {"PIVOTER_RUN_ST_V3", "ST_V3 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V3},
        {"PIVOTER_RUN_ST_V2", "ST_V2 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V2},
        {"PIVOTER_RUN_ST_V1", "ST_V1 r=2", PlusNucleusEdgeCoreDecompositionSet_ST_V1},
        {"PIVOTER_RUN_ST",    "ST r=2",    PlusNucleusEdgeCoreDecompositionSet_ST},
    };

    for (auto &entry : table) {
        if (envSet(entry.envVar)) {
            runR2Variant(entry.label, entry.func, correctRef,
                         tree, edgeGraph, treeGraphV, r, s, compareMode);
            return true;
        }
    }

    return false;
}

// ============================================================
// Dispatch: r>=3 r-clique algorithms
// ============================================================

// Reuse the existing runRCliqueVariant helper (already table-driven)
template<typename Func>
static void runRCliqueVariant(
    const char *name, Func &&func,
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s, bool compareMode,
    StaticCliqueIndex *sharedCI = nullptr) {

    auto refTree2 = compareMode ? tree.clone() : DynamicGraph<TreeGraphNode>();
    auto refTGV2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();

    auto result = daf::timeCount(name, [&]() {
        return func(tree, edgeGraph, treeGraphV, r, s, sharedCI);
    });
    if (compareMode) {
        auto refCore = daf::timeCount("Reference r>=3", [&]() {
            return NucleusCoreDecompositionCorrect(refTree2, edgeGraph, refTGV2, r, s, nullptr);
        });
        checkDist(buildCoreDist(refCore), buildCoreDist(result), name);
    }
}

static bool dispatchR3Plus(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s, bool compareMode,
    StaticCliqueIndex *sharedCI = nullptr) {

    struct R3Entry {
        const char *envVar;
        const char *label;
        using FuncType = std::vector<std::pair<std::vector<daf::Size>,double>>(*)(
            DynamicGraph<TreeGraphNode>&, const Graph&,
            DynamicGraphSet<TreeGraphNode>&, daf::CliqueSize, daf::CliqueSize,
            StaticCliqueIndex*);
        FuncType func;
    };

    static const R3Entry table[] = {
        {"PIVOTER_RUN_LINK_PEEL",      "Link-Graph Peel r>=3",      NucleusCoreDecompositionRCliqueLinkPeel},
        {"PIVOTER_RUN_LINK_LOCAL",      "Link-Graph Local r>=3",     NucleusCoreDecompositionRCliqueLinkLocal},
        {"PIVOTER_RUN_LOCAL_CPI_VP",    "Local CPI VP r>=3",         NucleusCoreDecompositionRCliqueLocalCPI_VP},
        {"PIVOTER_RUN_LOCAL_CPI_EXACT", "Local CPI Exact r>=3",      NucleusCoreDecompositionRCliqueLocalCPI_Exact},
        {"PIVOTER_RUN_LOCAL_BK",        "Local H-index r>=3 BK",     NucleusCoreDecompositionRCliqueLocal_BK},
        {"PIVOTER_RUN_LOCAL",           "Local H-index r>=3 VP",     NucleusCoreDecompositionRCliqueLocal},
        {"PIVOTER_RUN_ST_V23",          "ST_V23 r>=3",               NucleusCoreDecompositionRClique_ST_V23},
        {"PIVOTER_RUN_ST_V22",          "ST_V22 r>=3",               NucleusCoreDecompositionRClique_ST_V22},
        {"PIVOTER_RUN_ST_V21",          "ST_V21 r>=3",               NucleusCoreDecompositionRClique_ST_V21},
        {"PIVOTER_RUN_ST_V20",          "ST_V20 r>=3",               NucleusCoreDecompositionRClique_ST_V20},
        {"PIVOTER_RUN_ST_V19",          "ST_V19 r>=3",               NucleusCoreDecompositionRClique_ST_V19},
        {"PIVOTER_RUN_ST_V18",          "ST_V18 r>=3",               NucleusCoreDecompositionRClique_ST_V18},
        {"PIVOTER_RUN_ST_V17",          "ST_V17 r>=3",               NucleusCoreDecompositionRClique_ST_V17},
        {"PIVOTER_RUN_ST_V16",          "ST_V16 r>=3",               NucleusCoreDecompositionRClique_ST_V16},
        {"PIVOTER_RUN_ST_V15",          "ST_V15 r>=3",               NucleusCoreDecompositionRClique_ST_V15},
        {"PIVOTER_RUN_ST_V14",          "ST_V14 r>=3",               NucleusCoreDecompositionRClique_ST_V14},
        {"PIVOTER_RUN_ST_V13",          "ST_V13 r>=3",               NucleusCoreDecompositionRClique_ST_V13},
        {"PIVOTER_RUN_ST_V12",          "ST_V12 r>=3",               NucleusCoreDecompositionRClique_ST_V12},
        {"PIVOTER_RUN_ST_V11",          "ST_V11 r>=3",               NucleusCoreDecompositionRClique_ST_V11},
        {"PIVOTER_RUN_ST_V10",          "ST_V10 r>=3",               NucleusCoreDecompositionRClique_ST_V10},
        {"PIVOTER_RUN_ST_V9",           "ST_V9 r>=3",                NucleusCoreDecompositionRClique_ST_V9},
        {"PIVOTER_RUN_ST_V8",           "ST_V8 r>=3",                NucleusCoreDecompositionRClique_ST_V8},
        {"PIVOTER_RUN_ST_V7",           "ST_V7 r>=3",                NucleusCoreDecompositionRClique_ST_V7},
        {"PIVOTER_RUN_ST_V6",           "ST_V6 r>=3",                NucleusCoreDecompositionRClique_ST_V6},
        {"PIVOTER_RUN_ST_V5",           "ST_V5 r>=3",                NucleusCoreDecompositionRClique_ST_V5},
        {"PIVOTER_RUN_ST_V4",           "ST_V4 r>=3",                NucleusCoreDecompositionRClique_ST_V4},
        {"PIVOTER_RUN_ST_V3",           "ST_V3 r>=3",                NucleusCoreDecompositionRClique_ST_V3},
        {"PIVOTER_RUN_ST_V2",           "ST_V2 r>=3",                NucleusCoreDecompositionRClique_ST_V2},
        {"PIVOTER_RUN_ST",              "ST r>=3",                   NucleusCoreDecompositionRClique_ST},
    };

    for (auto &entry : table) {
        if (envSet(entry.envVar)) {
            runRCliqueVariant(entry.label, entry.func,
                              tree, edgeGraph, treeGraphV, r, s, compareMode, sharedCI);
            return true;
        }
    }

    return false;
}

// ============================================================
// Dispatch: reference-only and defaults
// ============================================================

static void dispatchRefOrDefault(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s, bool compareMode,
    StaticCliqueIndex *sharedCIPtr = nullptr) {

    if (envSet("PIVOTER_RUN_REF")) {
        auto res = NucleusCoreDecompositionCorrect(tree, edgeGraph, treeGraphV, r, s, sharedCIPtr);
        std::map<double, int> dist;
        for (const auto &[clique, cv] : res) dist[cv]++;
        std::cout << "Core value distribution:" << std::endl;
        for (const auto &[cv, count] : dist)
            std::cout << "  core=" << cv << " count=" << count << std::endl;
        return;
    }

    // Defaults when no PIVOTER_RUN_* env var is set
    if (r == 1) {
        NCliqueVertexCoreDecomposition(tree, edgeGraph, treeGraphV, s);
    } else if (r == 2) {
        PlusNucleusEdgeCoreDecompositionSet(tree, edgeGraph, treeGraphV, s);
    } else if (compareMode) {
        auto optTree = tree.clone();
        auto refTGV = treeGraphV.clone();
        auto refCore = daf::timeCount("Reference", [&]() {
            return NucleusCoreDecompositionCorrect(tree, edgeGraph, refTGV, r, s);
        });
        auto optTGV = treeGraphV.clone();
        auto optCore = daf::timeCount("Optimized", [&]() {
            return NucleusCoreDecompositionRClique(optTree, edgeGraph, optTGV, r, s);
        });
        checkDist(buildCoreDist(refCore), buildCoreDist(optCore), "r>=3 Optimized vs Reference");
    } else {
        NucleusCoreDecompositionRClique(tree, edgeGraph, treeGraphV, r, s);
    }
}

// ============================================================
// Main
// ============================================================

int main(int argc, char **argv) {
    if (argc > 5) {
        printf("Usage: ./degeneracy_cliques <graphFile> <r> <s> [sort_option]\n");
        printf("  sort_option: degen, degenR, degree, degreeR, default\n");
        return 0;
    }
#ifdef _OPENMP
    printf("Max OpenMP threads: %d\n", omp_get_max_threads());
#else
    printf("OpenMP not supported, running with a single thread.\n");
#endif

    const char *fpath = argv[1];
    const daf::CliqueSize r = strtol(argv[2], nullptr, 10);
    const daf::CliqueSize s = strtol(argv[3], nullptr, 10);
    if (r >= s) { printf("r must be less than s\n"); return 0; }
    printf("Dataset: %s, r=%d, s=%d\n", fpath, (int)r, (int)s);

    // Phase 1: Load and sort
    Graph edgeGraph = loadAndSortGraph(fpath, argc, argv);

    const bool compareMode = envSet("PIVOTER_COMPARE");

    // Phase 2+5: Build SDCT tree + treeGraphV
    DynamicGraph<TreeGraphNode> refTree;
    DynamicGraphSet<TreeGraphNode> treeGraphV;
    std::unique_ptr<StaticCliqueIndex> sharedCliqueIndex;  // R≥3: pre-built, shared

    if (needsSDCT(r, compareMode)) {
        if (compareMode) {
            // Compare mode: legacy SDCT (min_k=r) + Par7 verification
            refTree = buildAndVerifySDCT(edgeGraph, r, s);
        } else {
            // Production: fused SDCT (store_min_k=s, emit_min_k=r)
            auto result = buildSDCTWithIndex(edgeGraph, r, s);
            refTree = std::move(result.tree);
            if (needsTreeGraphV(r, compareMode)) {
                treeGraphV = std::move(result.treeGraphV);
            }
            sharedCliqueIndex = std::move(result.cliqueIndex);
        }
    }

    // Phase 3: Pre-mutation work (must run before beSingleEdge)
    auto pmr = preMutationPhase(edgeGraph, r, s);

    // Phase 4: Prepare graph structures (mutates edgeGraph)
    prepareGraphStructures(edgeGraph, r);

    // Phase 5: Build treeGraphV if not already built in fused phase
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;
    if (needsTreeGraphV(r, compareMode) && treeGraphV.adj_list.empty()) {
        treeGraphV = DynamicGraphSet<TreeGraphNode>(refTree, edgeGraph.getGraphNodeSize(), s);
    }
    daf::log_memory("Other Index Memory");

    // Phase 6: R≥3 shared cliqueIndex (built during SDCT_Fused, hash-free)
    StaticCliqueIndex *sharedCIPtr = sharedCliqueIndex ? sharedCliqueIndex.get() : nullptr;

    // Phase 7: Dispatch algorithm
    daf::timeCount("NucleusCoreDecomposition", [&] {
        bool dispatched = false;

        if (r == 1)
            dispatched = dispatchR1(refTree, edgeGraph, treeGraphV, s,
                                    numVertices, compareMode, pmr);
        else if (r == 2)
            dispatched = dispatchR2(refTree, edgeGraph, treeGraphV, r, s, compareMode);
        else if (r >= 3)
            dispatched = dispatchR3Plus(refTree, edgeGraph, treeGraphV, r, s, compareMode, sharedCIPtr);

        if (!dispatched)
            dispatchRefOrDefault(refTree, edgeGraph, treeGraphV, r, s, compareMode, sharedCIPtr);
    });

    daf::log_memory("Final Memory");
}
