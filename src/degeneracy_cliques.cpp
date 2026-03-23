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
// Helpers: compare core-value distributions
// ============================================================

// Compare r-clique results (vector<pair<vector<Size>, int>>)
static void compareRCliqueDist(
    const std::vector<std::pair<std::vector<daf::Size>, int>> &refCore,
    const std::vector<std::pair<std::vector<daf::Size>, int>> &testCore,
    const char *label) {
    std::map<int, int> refDist, testDist;
    for (const auto &[clique, cv] : refCore) refDist[cv]++;
    for (const auto &[clique, cv] : testCore) testDist[cv]++;
    refDist.erase(0);
    testDist.erase(0);
    if (refDist == testDist) {
        std::cout << "✓ " << label << " correctness verified" << std::endl;
    } else {
        std::cerr << "✗ " << label << " MISMATCH!" << std::endl;
        for (auto &[k, v] : refDist)
            if (testDist.count(k) == 0 || testDist[k] != v)
                std::cerr << "  core=" << k << " ref=" << v
                          << " test=" << (testDist.count(k) ? testDist[k] : 0) << std::endl;
    }
}

// Compare edge results (vector<pair<pair<Size,Size>, int>>)
static void compareEdgeDist(
    const std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> &refCore,
    const std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> &testCore,
    const char *label) {
    std::map<int, int> refDist, testDist;
    for (const auto &[edge, cv] : refCore) refDist[cv]++;
    for (const auto &[edge, cv] : testCore) testDist[cv]++;
    refDist.erase(0);
    testDist.erase(0);
    if (refDist == testDist) {
        std::cout << "✓ " << label << " correctness verified" << std::endl;
    } else {
        std::cerr << "✗ " << label << " MISMATCH!" << std::endl;
        for (auto &[k, v] : refDist)
            if (testDist.count(k) == 0 || testDist[k] != v)
                std::cerr << "  core=" << k << " ref=" << v
                          << " test=" << (testDist.count(k) ? testDist[k] : 0) << std::endl;
    }
}

// Compare vertex results (double* vs reference r-clique result)
static void compareVertexDist(
    const double *testCoreV, daf::Size numVertices,
    const std::vector<std::pair<std::vector<daf::Size>, int>> &refCore,
    const char *label) {
    std::map<int, int> refDist, testDist;
    for (const auto &[clique, cv] : refCore) refDist[cv]++;
    for (daf::Size i = 0; i < numVertices; ++i)
        if (testCoreV[i] >= 0) testDist[(int)testCoreV[i]]++;
    refDist.erase(0);
    testDist.erase(0);
    if (refDist == testDist) {
        std::cout << "✓ " << label << " correctness verified" << std::endl;
    } else {
        std::cerr << "✗ " << label << " MISMATCH!" << std::endl;
        std::cerr << "refDist: " << refDist << std::endl;
        std::cerr << "testDist: " << testDist << std::endl;
    }
}

// Compare two double* vertex results
static void compareVertexDistBoth(
    const double *refV, const double *testV, daf::Size numVertices,
    const char *label) {
    std::map<int, int> refDist, testDist;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (refV[i] >= 0) refDist[(int)refV[i]]++;
        if (testV[i] >= 0) testDist[(int)testV[i]]++;
    }
    refDist.erase(0);
    testDist.erase(0);
    if (refDist == testDist) {
        std::cout << "✓ " << label << " correctness verified" << std::endl;
    } else {
        std::cerr << "✗ " << label << " MISMATCH!" << std::endl;
        std::cerr << "refDist: " << refDist << std::endl;
        std::cerr << "testDist: " << testDist << std::endl;
    }
}

static bool envSet(const char *name) { return std::getenv(name) != nullptr; }

// ============================================================
// Helpers: run r>=3 variant with optional comparison
// ============================================================

// Run r>=3 variant with optional comparison against reference
template<typename Func>
static void runRCliqueVariant(
    const char *name,
    Func &&func,
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
        compareRCliqueDist(refCore, result, name);
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

    // --- Load graph ---
    Graph edgeGraph(fpath);
    edgeGraph.printGraphInfo();
    populate_nCr();
    daf::vListMap.resize(edgeGraph.n + 1);
    memset(daf::vListMap.data(), -1, edgeGraph.n * sizeof(daf::Size));

    // --- Sort ---
    if (argc >= 5) {
        std::string sortOption = argv[4];
        if (sortOption == "degen") edgeGraph.sortByDegeneracyOrder(false);
        else if (sortOption == "degenR") edgeGraph.sortByDegeneracyOrder(true);
        else if (sortOption == "degree") edgeGraph.sortByDegree(false);
        else if (sortOption == "degreeR") edgeGraph.sortByDegree(true);
        else if (sortOption == "default") { /* no-op */ }
        else {
            std::cout << "Unknown sort option: " << sortOption << std::endl;
            std::cout << "Available: degen, degenR, degree, degreeR, default" << std::endl;
            return 1;
        }
        std::cout << "Graph sorted by " << sortOption << std::endl;
    } else {
        edgeGraph.sortByDegeneracyOrder();
    }
    daf::log_memory("Graph Memory");

    // --- Build SDCT ---
    auto refTree = daf::timeCount("SDCT (reference)", [&]() {
        return SDCT(edgeGraph, s, r);
    });
    auto refCC = refTree.cliqueCount();
    printf("SDCT leaf count: %zu, maxK: %zu\n", refTree.adj_list.size(), (size_t)refCC.size()-1);
    printf("cliqueCount:");
    for (size_t k = 1; k < refCC.size(); k++)
        if (refCC[k] > 0) printf(" [%zu]=%.0f", k, refCC[k]);
    printf("\n");

    // --- Verify SDCT_Par7 ---
    {
        auto par7Result = daf::timeCount("SDCT_Par7", [&]() -> CliqueCSR<int> {
            return SDCT_Par7(edgeGraph, s, r);
        });
        printf("SDCT_Par7 leaf count: %zu\n", par7Result.num_cliques());
        if (par7Result.has_pivot()) {
            auto par7CC = par7Result.cliqueCount();
            bool match = true;
            size_t maxCCSize = refCC.size() > par7CC.size() ? refCC.size() : par7CC.size();
            for (size_t k = (size_t)s; k <= (size_t)s && k < maxCCSize; k++) {
                double rv = (k < refCC.size()) ? refCC[k] : 0.0;
                double pv = (k < par7CC.size()) ? par7CC[k] : 0.0;
                double diff = std::abs(rv - pv);
                double maxVal = std::max(std::abs(rv), std::abs(pv));
                if (diff > 0.5 && (maxVal < 1e-10 || diff / maxVal > 1e-6)) {
                    printf("  ✗ cliqueCount mismatch at k=%zu: ref=%.0f par7=%.0f\n", k, rv, pv);
                    match = false;
                }
            }
            if (match) printf("✓ SDCT_Par7 cliqueCount correct (k=%d)\n", (int)s);
            else { printf("✗ SDCT_Par7 cliqueCount WRONG\n"); return 1; }
        }
    } // par7Result freed here

    // --- ST_V2: Build phase MUST run before beSingleEdge() mutates the graph ---
    const bool st_v2_mode = (r == 1 && envSet("PIVOTER_RUN_ST_V2"));
    const bool st_v2_probe = (r == 1 && envSet("PIVOTER_RUN_ST_V2_PROBE"));
    std::unique_ptr<ST_V2_Data> st_v2_data;
    if (st_v2_mode) {
        st_v2_data = std::make_unique<ST_V2_Data>(
            daf::timeCount("ST_V2 Build", [&]() {
                return NCliqueVertexCoreDecomposition_ST_V2_Build(edgeGraph, s);
            }));
    }
    if (st_v2_probe) {
        NCliqueVertexCoreDecomposition_ST_V2_InterleavedProbe(edgeGraph, s);
    }

    // --- Prepare graph structures ---
    edgeGraph.initCore();
    edgeGraph.coreV.free(); // Opt 2: coreV only needed by sortByDegeneracyOrder, already done
    edgeGraph.beSingleEdge();
    edgeGraph.buildEdgeIdMap();
    // Opt 3: eidToNode only needed by getEdgeById(), which r=1 and r>=3 never call
    if (r != 2) {
        edgeGraph.eidToNode.free();
    }

    // --- Read mode flags (moved before treeGraphV to enable conditional construction) ---
    const bool compareMode = envSet("PIVOTER_COMPARE");
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;

    // Opt 6: R1 ST/ST_V2 never uses treeGraphV — skip construction to save n × ~60B hash set overhead.
    // Build treeGraphV only when needed: r>=2, or r=1 non-ST variants, or compare mode.
    const bool r1_st_only = (r == 1 && (envSet("PIVOTER_RUN_ST") || envSet("PIVOTER_RUN_ST_V2")) && !compareMode);
    DynamicGraphSet<TreeGraphNode> treeGraphV;
    if (!r1_st_only) {
        treeGraphV = DynamicGraphSet<TreeGraphNode>(refTree, edgeGraph.getGraphNodeSize(), s);
    }
    daf::log_memory("Other Index Memory");

    // ============================================================
    // Dispatch: use PIVOTER_RUN_* env vars to select algorithm
    // ============================================================
    daf::timeCount("NucleusCoreDecomposition", [&] {

        // ---- r=1 vertex-level algorithms ----
        if (r == 1 && envSet("PIVOTER_RUN_LOCAL_V4")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto coreV = daf::timeCount("Local H-index V4 r=1", [&]() {
                return NCliqueVertexCoreDecomposition_LocalV4(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
                compareVertexDistBoth(refV, coreV, numVertices, "r=1 Local H-index V4");
                delete[] refV;
            }
            delete[] coreV;
        } else if (r == 1 && envSet("PIVOTER_RUN_LOCAL_V3")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto coreV = daf::timeCount("Local H-index V3 r=1", [&]() {
                return NCliqueVertexCoreDecomposition_LocalV3(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
                compareVertexDistBoth(refV, coreV, numVertices, "r=1 Local H-index V3");
                delete[] refV;
            }
            delete[] coreV;
        } else if (r == 1 && envSet("PIVOTER_RUN_LOCAL_V2")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto coreV = daf::timeCount("Local H-index V2 r=1", [&]() {
                return NCliqueVertexCoreDecomposition_LocalV2(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
                compareVertexDistBoth(refV, coreV, numVertices, "r=1 Local H-index V2");
                delete[] refV;
            }
            delete[] coreV;
        } else if (r == 1 && envSet("PIVOTER_RUN_LOCAL_NAIVE")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto coreV = daf::timeCount("Local H-index Naive r=1", [&]() {
                return NCliqueVertexCoreDecomposition_LocalNaive(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
                compareVertexDistBoth(refV, coreV, numVertices, "r=1 Local H-index Naive");
                delete[] refV;
            }
            delete[] coreV;
        } else if (r == 1 && envSet("PIVOTER_RUN_LOCAL")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto coreV = daf::timeCount("Local H-index r=1", [&]() {
                return NCliqueVertexCoreDecomposition_Local(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
                compareVertexDistBoth(refV, coreV, numVertices, "r=1 Local H-index");
                delete[] refV;
            }
            delete[] coreV;
        } else if (r == 1 && envSet("PIVOTER_RUN_ST_V2")) {
            // ST V2: Tree-free R1 — SDCT already built before beSingleEdge()
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto coreV = daf::timeCount("ST_V2 r=1 (peel)", [&]() {
                return NCliqueVertexCoreDecomposition_ST_V2_Peel(*st_v2_data, s);
            });
            if (compareMode) {
                auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
                compareVertexDistBoth(refV, coreV, numVertices, "r=1 ST_V2");
                delete[] refV;
            }
            delete[] coreV;
        } else if (r == 1 && envSet("PIVOTER_RUN_ST")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto coreV = daf::timeCount("ST r=1", [&]() {
                return NCliqueVertexCoreDecomposition_ST(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refV = NCliqueVertexCoreDecomposition(t2, edgeGraph, tgv2, s);
                compareVertexDistBoth(refV, coreV, numVertices, "r=1 ST");
                delete[] refV;
            }
            delete[] coreV;

        // ---- r=2 edge-level algorithms ----
        } else if (r == 2 && envSet("PIVOTER_RUN_LOCAL")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto result = daf::timeCount("Local H-index r=2", [&]() {
                return NCliqueEdgeCoreDecomposition_Local(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refCore = NucleusCoreDecompositionCorrect(t2, edgeGraph, tgv2, r, s);
                // Convert r-clique ref to edge dist for comparison
                std::map<int, int> refDist, testDist;
                for (const auto &[c, cv] : refCore) refDist[cv]++;
                for (const auto &[e, cv] : result) testDist[cv]++;
                refDist.erase(0); testDist.erase(0);
                if (refDist == testDist) std::cout << "✓ r=2 Local H-index correctness verified" << std::endl;
                else std::cerr << "✗ r=2 Local H-index MISMATCH!" << std::endl;
            }
        } else if (r == 2 && envSet("PIVOTER_RUN_ST_V4")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto result = daf::timeCount("ST_V4 r=2", [&]() {
                return PlusNucleusEdgeCoreDecompositionSet_ST_V4(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refCore = daf::timeCount("ref r=2", [&]() {
                    return PlusNucleusEdgeCoreDecompositionSet(t2, edgeGraph, tgv2, s);
                });
                std::map<int, int> refDist, testDist;
                for (const auto &[c, cv] : refCore) refDist[cv]++;
                for (const auto &[e, cv] : result) testDist[cv]++;
                if (refDist == testDist) std::cout << "✓ r=2 ST_V4 correctness verified" << std::endl;
                else std::cerr << "✗ r=2 ST_V4 MISMATCH!" << std::endl;
            }
        } else if (r == 2 && envSet("PIVOTER_RUN_ST_V3")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto result = daf::timeCount("ST_V3 r=2", [&]() {
                return PlusNucleusEdgeCoreDecompositionSet_ST_V3(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refCore = NucleusCoreDecompositionCorrect(t2, edgeGraph, tgv2, r, s);
                std::map<int, int> refDist, testDist;
                for (const auto &[c, cv] : refCore) refDist[cv]++;
                for (const auto &[e, cv] : result) testDist[cv]++;
                if (refDist == testDist) std::cout << "✓ r=2 ST_V3 correctness verified" << std::endl;
                else std::cerr << "✗ r=2 ST_V3 MISMATCH!" << std::endl;
            }
        } else if (r == 2 && envSet("PIVOTER_RUN_ST_V2")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto result = daf::timeCount("ST_V2 r=2", [&]() {
                return PlusNucleusEdgeCoreDecompositionSet_ST_V2(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refCore = NucleusCoreDecompositionCorrect(t2, edgeGraph, tgv2, r, s);
                std::map<int, int> refDist, testDist;
                for (const auto &[c, cv] : refCore) refDist[cv]++;
                for (const auto &[e, cv] : result) testDist[cv]++;
                if (refDist == testDist) std::cout << "✓ r=2 ST_V2 correctness verified" << std::endl;
                else std::cerr << "✗ r=2 ST_V2 MISMATCH!" << std::endl;
            }
        } else if (r == 2 && envSet("PIVOTER_RUN_ST_V1")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto result = daf::timeCount("ST_V1 r=2", [&]() {
                return PlusNucleusEdgeCoreDecompositionSet_ST_V1(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refCore = NucleusCoreDecompositionCorrect(t2, edgeGraph, tgv2, r, s);
                std::map<int, int> refDist, testDist;
                for (const auto &[c, cv] : refCore) refDist[cv]++;
                for (const auto &[e, cv] : result) testDist[cv]++;
                if (refDist == testDist) std::cout << "✓ r=2 ST_V1 correctness verified" << std::endl;
                else std::cerr << "✗ r=2 ST_V1 MISMATCH!" << std::endl;
            }
        } else if (r == 2 && envSet("PIVOTER_RUN_ST")) {
            auto t2 = compareMode ? refTree.clone() : DynamicGraph<TreeGraphNode>();
            auto tgv2 = compareMode ? treeGraphV.clone() : DynamicGraphSet<TreeGraphNode>();
            auto result = daf::timeCount("ST r=2", [&]() {
                return PlusNucleusEdgeCoreDecompositionSet_ST(refTree, edgeGraph, treeGraphV, s);
            });
            if (compareMode) {
                auto refCore = daf::timeCount("ref r=2", [&]() {
                    return PlusNucleusEdgeCoreDecompositionSet(t2, edgeGraph, tgv2, s);
                });
                std::map<int, int> refDist, testDist;
                for (const auto &[c, cv] : refCore) refDist[cv]++;
                for (const auto &[e, cv] : result) testDist[cv]++;
                if (refDist == testDist) std::cout << "✓ r=2 ST correctness verified" << std::endl;
                else std::cerr << "✗ r=2 ST MISMATCH!" << std::endl;
            }

        // ---- r>=3 algorithms (table-driven via runRCliqueVariant) ----
        } else if (r >= 3 && envSet("PIVOTER_RUN_LINK_PEEL")) {
            runRCliqueVariant("Link-Graph Peel r>=3",
                NucleusCoreDecompositionRCliqueLinkPeel,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_LINK_LOCAL")) {
            runRCliqueVariant("Link-Graph Local r>=3",
                NucleusCoreDecompositionRCliqueLinkLocal,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_LOCAL_CPI_VP")) {
            runRCliqueVariant("Local CPI VP r>=3",
                NucleusCoreDecompositionRCliqueLocalCPI_VP,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_LOCAL_CPI_EXACT")) {
            runRCliqueVariant("Local CPI Exact r>=3",
                NucleusCoreDecompositionRCliqueLocalCPI_Exact,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_LOCAL_BK")) {
            runRCliqueVariant("Local H-index r>=3 BK",
                NucleusCoreDecompositionRCliqueLocal_BK,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_LOCAL")) {
            runRCliqueVariant("Local H-index r>=3 vertex-proxy",
                NucleusCoreDecompositionRCliqueLocal,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V11")) {
            runRCliqueVariant("ST_V11 r>=3",
                NucleusCoreDecompositionRClique_ST_V11,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V10")) {
            runRCliqueVariant("ST_V10 r>=3",
                NucleusCoreDecompositionRClique_ST_V10,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V9")) {
            runRCliqueVariant("ST_V9 r>=3",
                NucleusCoreDecompositionRClique_ST_V9,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V8")) {
            runRCliqueVariant("ST_V8 r>=3",
                NucleusCoreDecompositionRClique_ST_V8,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V7")) {
            runRCliqueVariant("ST_V7 r>=3",
                NucleusCoreDecompositionRClique_ST_V7,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V6")) {
            runRCliqueVariant("ST_V6 r>=3",
                NucleusCoreDecompositionRClique_ST_V6,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V5")) {
            runRCliqueVariant("ST_V5 r>=3",
                NucleusCoreDecompositionRClique_ST_V5,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V4")) {
            runRCliqueVariant("ST_V4 r>=3",
                NucleusCoreDecompositionRClique_ST_V4,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V3")) {
            runRCliqueVariant("ST_V3 r>=3",
                NucleusCoreDecompositionRClique_ST_V3,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST_V2")) {
            runRCliqueVariant("ST_V2 r>=3",
                NucleusCoreDecompositionRClique_ST_V2,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);
        } else if (r >= 3 && envSet("PIVOTER_RUN_ST")) {
            runRCliqueVariant("ST r>=3",
                NucleusCoreDecompositionRClique_ST,
                refTree, edgeGraph, treeGraphV, r, s, compareMode);

        // ---- Reference-only mode ----
        } else if (envSet("PIVOTER_RUN_REF")) {
            auto res = NucleusCoreDecompositionCorrect(refTree, edgeGraph, treeGraphV, r, s);
            std::map<daf::Size, int> dist;
            for (const auto &[clique, cv] : res) dist[cv]++;
            std::cout << "Core value distribution:" << std::endl;
            for (const auto &[cv, count] : dist)
                std::cout << "  core=" << cv << " count=" << count << std::endl;

        // ---- Defaults (no env var) ----
        } else if (r == 1) {
            NCliqueVertexCoreDecomposition(refTree, edgeGraph, treeGraphV, s);
        } else if (r == 2) {
            PlusNucleusEdgeCoreDecompositionSet(refTree, edgeGraph, treeGraphV, s);
        } else if (compareMode) {
            // r>=3 default compare: optimized vs reference
            auto optTree = refTree.clone();
            auto refTGV = treeGraphV.clone();
            auto refCore = daf::timeCount("Reference", [&]() {
                return NucleusCoreDecompositionCorrect(refTree, edgeGraph, refTGV, r, s);
            });
            auto optTGV = treeGraphV.clone();
            auto optCore = daf::timeCount("Optimized", [&]() {
                return NucleusCoreDecompositionRClique(optTree, edgeGraph, optTGV, r, s);
            });
            compareRCliqueDist(refCore, optCore, "r>=3 Optimized vs Reference");
        } else {
            NucleusCoreDecompositionRClique(refTree, edgeGraph, treeGraphV, r, s);
        }
    });

    daf::log_memory("Final Memory");
}