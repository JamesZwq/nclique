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
    std::map<int, int> dist;
    for (const auto &[key, cv] : coreResults) dist[cv]++;
    dist.erase(0);
    return dist;
}

static auto buildCoreDistFromArray(const double *coreV, daf::Size n) {
    std::map<int, int> dist;
    for (daf::Size i = 0; i < n; ++i)
        if (coreV[i] >= 0) dist[(int)coreV[i]]++;
    dist.erase(0);
    return dist;
}

static void checkDist(const std::map<int,int> &refDist,
                      const std::map<int,int> &testDist,
                      const char *label) {
    if (refDist == testDist) {
        std::cout << "✓ " << label << " correctness verified" << std::endl;
    } else {
        std::cerr << "✗ " << label << " MISMATCH!" << std::endl;
        for (auto &[k, v] : refDist)
            if (!testDist.count(k) || testDist.at(k) != v)
                std::cerr << "  core=" << k << " ref=" << v
                          << " test=" << (testDist.count(k) ? testDist.at(k) : 0) << std::endl;
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
    if (r == 1 && (envSet("PIVOTER_RUN_ST_V2") || envSet("PIVOTER_RUN_ST_V2_PROBE")))
        return false;
    return true;
}

static DynamicGraph<TreeGraphNode> buildAndVerifySDCT(
    Graph &edgeGraph, daf::CliqueSize r, daf::CliqueSize s) {

    auto refTree = daf::timeCount("SDCT (reference)", [&]() {
        return SDCT(edgeGraph, s, r);
    });

    auto refCC = refTree.cliqueCount();
    printf("SDCT leaf count: %zu, maxK: %zu\n", refTree.adj_list.size(), (size_t)refCC.size()-1);
    printf("cliqueCount:");
    for (size_t k = 1; k < refCC.size(); k++)
        if (refCC[k] > 0) printf(" [%zu]=%.0f", k, refCC[k]);
    printf("\n");

    // Verify SDCT_Par7
    {
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

static std::unique_ptr<ST_V2_Data> preMutationPhase(
    Graph &edgeGraph, daf::CliqueSize r, daf::CliqueSize s) {

    std::unique_ptr<ST_V2_Data> st_v2_data;

    if (r == 1 && envSet("PIVOTER_RUN_ST_V2")) {
        st_v2_data = std::make_unique<ST_V2_Data>(
            daf::timeCount("ST_V2 Build", [&]() {
                return NCliqueVertexCoreDecomposition_ST_V2_Build(edgeGraph, s);
            }));
    }
    if (r == 1 && envSet("PIVOTER_RUN_ST_V2_PROBE")) {
        NCliqueVertexCoreDecomposition_ST_V2_InterleavedProbe(edgeGraph, s);
    }

    return st_v2_data;
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
    if (r >= 2) return true;
    // r=1: only ST and ST_V2 can skip treeGraphV
    if (envSet("PIVOTER_RUN_ST") || envSet("PIVOTER_RUN_ST_V2")) return false;
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
    ST_V2_Data *st_v2_data) {

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
    if (envSet("PIVOTER_RUN_ST_V2") && st_v2_data) {
        auto t2 = compareMode ? tree.clone() : DynamicGraph<TreeGraphNode>();
        auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
        auto coreV = daf::timeCount("ST_V2 r=1 (peel)", [&]() {
            return NCliqueVertexCoreDecomposition_ST_V2_Peel(*st_v2_data, s);
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
        return NucleusCoreDecompositionCorrect(t, g, tgv, r, s);
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
    daf::CliqueSize r, daf::CliqueSize s, bool compareMode) {

    auto refTree2 = compareMode ? tree.clone() : DynamicGraph<TreeGraphNode>();
    auto refTGV2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();

    auto result = daf::timeCount(name, [&]() {
        return func(tree, edgeGraph, treeGraphV, r, s);
    });
    if (compareMode) {
        auto refCore = daf::timeCount("Reference r>=3", [&]() {
            return NucleusCoreDecompositionCorrect(refTree2, edgeGraph, refTGV2, r, s);
        });
        checkDist(buildCoreDist(refCore), buildCoreDist(result), name);
    }
}

static bool dispatchR3Plus(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s, bool compareMode) {

    struct R3Entry {
        const char *envVar;
        const char *label;
        using FuncType = std::vector<std::pair<std::vector<daf::Size>,int>>(*)(
            DynamicGraph<TreeGraphNode>&, const Graph&,
            DynamicGraphSet<TreeGraphNode>&, daf::CliqueSize, daf::CliqueSize);
        FuncType func;
    };

    static const R3Entry table[] = {
        {"PIVOTER_RUN_LINK_PEEL",      "Link-Graph Peel r>=3",      NucleusCoreDecompositionRCliqueLinkPeel},
        {"PIVOTER_RUN_LINK_LOCAL",      "Link-Graph Local r>=3",     NucleusCoreDecompositionRCliqueLinkLocal},
        {"PIVOTER_RUN_LOCAL_CPI_VP",    "Local CPI VP r>=3",         NucleusCoreDecompositionRCliqueLocalCPI_VP},
        {"PIVOTER_RUN_LOCAL_CPI_EXACT", "Local CPI Exact r>=3",      NucleusCoreDecompositionRCliqueLocalCPI_Exact},
        {"PIVOTER_RUN_LOCAL_BK",        "Local H-index r>=3 BK",     NucleusCoreDecompositionRCliqueLocal_BK},
        {"PIVOTER_RUN_LOCAL",           "Local H-index r>=3 VP",     NucleusCoreDecompositionRCliqueLocal},
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
                              tree, edgeGraph, treeGraphV, r, s, compareMode);
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
    daf::CliqueSize r, daf::CliqueSize s, bool compareMode) {

    if (envSet("PIVOTER_RUN_REF")) {
        auto res = NucleusCoreDecompositionCorrect(tree, edgeGraph, treeGraphV, r, s);
        std::map<daf::Size, int> dist;
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

    // Phase 2: Build SDCT (skip if ST_V2 only — it builds its own data)
    DynamicGraph<TreeGraphNode> refTree;
    if (needsSDCT(r, compareMode)) {
        refTree = buildAndVerifySDCT(edgeGraph, r, s);
    }

    // Phase 3: Pre-mutation work (must run before beSingleEdge)
    auto st_v2_data = preMutationPhase(edgeGraph, r, s);

    // Phase 4: Prepare graph structures (mutates edgeGraph)
    prepareGraphStructures(edgeGraph, r);

    // Phase 5: Build treeGraphV (conditional)
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;

    DynamicGraphSet<TreeGraphNode> treeGraphV;
    if (needsTreeGraphV(r, compareMode)) {
        treeGraphV = DynamicGraphSet<TreeGraphNode>(refTree, edgeGraph.getGraphNodeSize(), s);
    }
    daf::log_memory("Other Index Memory");

    // Phase 6: Dispatch algorithm
    daf::timeCount("NucleusCoreDecomposition", [&] {
        bool dispatched = false;

        if (r == 1)
            dispatched = dispatchR1(refTree, edgeGraph, treeGraphV, s,
                                    numVertices, compareMode, st_v2_data.get());
        else if (r == 2)
            dispatched = dispatchR2(refTree, edgeGraph, treeGraphV, r, s, compareMode);
        else if (r >= 3)
            dispatched = dispatchR3Plus(refTree, edgeGraph, treeGraphV, r, s, compareMode);

        if (!dispatched)
            dispatchRefOrDefault(refTree, edgeGraph, treeGraphV, r, s, compareMode);
    });

    daf::log_memory("Final Memory");
}