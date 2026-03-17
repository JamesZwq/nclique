#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<cstdlib>
#include<limits.h>
#include<unistd.h>
#include<libgen.h>
#include<chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
// #include <boost/version.hpp>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include "BK/BronKerboschRmEdge.hpp"
#include "BK/BronKerboschRmRClique.hpp"
#include "dataStruct/disJoinSet.hpp"
#include "graph/DynamicBipartiteGraph.hpp"
#include "NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "degeneracy_algorithm_cliques_V.h"
#include "dataStruct/CliqueCSR.hpp"



int main(int argc, char **argv) {


    // std::cout << "Boost version: " << BOOST_LIB_VERSION << std::endl;
    if (argc > 5 ) {
        printf("Incorrect number of arguments.\n");
        printf("./main <graphFile> <r> <s>\n");
        printf("graphFile: path to graph\n");
        printf("r: r\n");
        printf("s: s\n");

        return 0;
    }

    // char *opt = NULL;
    const char *fpath = argv[1];
    const daf::CliqueSize r = strtol(argv[2], nullptr, 10);
    const daf::CliqueSize s = strtol(argv[3], nullptr, 10);
    if (r >= s) {
        printf("r must be less than s\n");
        return 0;
    }
    printf("about to call runAndPrint for dataset %s\n", fpath);

    Graph edgeGraph(fpath);
    edgeGraph.printGraphInfo();
    populate_nCr();
    daf::vListMap.resize(edgeGraph.n + 1);
    // std::numeric_limits<daf::Size>::max();
    memset(daf::vListMap.data(), -1, edgeGraph.n * sizeof(daf::Size));
    std::cout << "Han Nbr C: " << edgeGraph.getNbrCount(154) << std::endl;
    // edgeGraph.sortByDegeneracyOrder(false);
    // edgeGraph.sortByDegeneracyOrder(true);
    // edgeGraph.sortByDegree(false);
    // edgeGraph.sortByDegree(true);

    if (argc >= 5) {
        auto sortOption = std::string(argv[4]);
        if (sortOption == "degen") {
            edgeGraph.sortByDegeneracyOrder(false);
        } else if (sortOption == "degenR") {
            edgeGraph.sortByDegeneracyOrder(true);
        } else if (sortOption == "degree") {
            edgeGraph.sortByDegree(false);
        } else if (sortOption == "degreeR") {
            edgeGraph.sortByDegree(true);
        } else if (sortOption == "default") {
            // do nothing
        } else {
            std::cout << "Unknown sort option: " << sortOption << std::endl;
            std::cout << "Available options: degen, degenR, degree, degreeR, default" << std::endl;
            return 1;
        }

        std::cout << "Graph sorted by " << sortOption << std::endl;
    } else {
        edgeGraph.sortByDegeneracyOrder();
    }
    // edgeGraph.sortByDegree();

    daf::log_memory("Graph Memory");

    // Get reference cliqueCount from original SDCT
    auto refTree = daf::timeCount("SDCT (reference)", [&]() {
        return SDCT(edgeGraph, 1000000, 0);
    });
    auto refCC = refTree.cliqueCount();
    printf("Reference SDCT leaf count: %zu, maxK: %zu\n", refTree.adj_list.size(), (size_t)refCC.size()-1);
    printf("Reference cliqueCount:");
    for (size_t k = 1; k < refCC.size(); k++) {
        if (refCC[k] > 0) printf(" [%zu]=%.0f", k, refCC[k]);
    }
    printf("\n");

    // Test SDCT_Par7 (stateless independent parallel BK)
    // Compare cliqueCount vectors (pivot-invariant) with floating-point tolerance
    auto par7Result = daf::timeCount("SDCT_Par7", [&]() -> CliqueCSR<int> {
        return SDCT_Par7(edgeGraph, 1000000, 0);
    });
    printf("SDCT_Par7 leaf count: %zu\n", par7Result.num_cliques());

    if (par7Result.has_pivot()) {
        auto par7CC = par7Result.cliqueCount();
        printf("Par7 cliqueCount:");
        for (size_t k = 1; k < par7CC.size(); k++) {
            if (par7CC[k] > 0) printf(" [%zu]=%.0f", k, par7CC[k]);
        }
        printf("\n");

        // Compare with tolerance: relative error < 1e-6 or absolute error < 0.5
        bool match = true;
        size_t maxK = std::max((size_t)refCC.size(), par7CC.size());
        for (size_t k = 0; k < maxK; k++) {
            double rv = (k < refCC.size()) ? refCC[k] : 0.0;
            double pv = (k < par7CC.size()) ? par7CC[k] : 0.0;
            double diff = std::abs(rv - pv);
            double maxVal = std::max(std::abs(rv), std::abs(pv));
            if (diff > 0.5 && (maxVal < 1e-10 || diff / maxVal > 1e-6)) {
                printf("  ✗ cliqueCount mismatch at k=%zu: ref=%.0f par7=%.0f diff=%.0f\n",
                       k, rv, pv, rv - pv);
                match = false;
            }
        }
        if (match)
            printf("✓ SDCT_Par7 cliqueCount correct\n");
        else {
            printf("✗ SDCT_Par7 cliqueCount WRONG\n");
            return 1;
        }
    } else {
        printf("  (no pivot flags — cannot verify cliqueCount)\n");
    }

    edgeGraph.initCore();
    // treeGraph.printGraphPerV();

    edgeGraph.beSingleEdge();
    edgeGraph.buildEdgeIdMap();
    // DynamicBipartiteGraph BGraph(treeGraph, edgeGraph);

    // DynamicGraphSet<TreeGraphNode> treeGraphV(treeGraph, edgeGraph.getGraphNodeSize(), s);
    DynamicGraphSet<TreeGraphNode> treeGraphV(refTree, edgeGraph.getGraphNodeSize(), s);


    // StaticCliqueIndex cliqueIndex(r);
    // daf::timeCount("clique Index build",
    //                [&]() {
    //                    cliqueIndex.build(treeGraph, edgeGraph.adj_list.size());
    //                });
    //
    // cliqueIndex.verify();
    // print all leaf with node 154
    for (daf::Size leafId = 0; leafId < refTree.adj_list.size(); ++leafId) {
        const auto &leaf = refTree.adj_list[leafId];
        for (const auto &node: leaf) {
            if (node.v == 154) {
                std::cout << "leafId: " << leafId << " leaf Size: " << leaf.size() << std::endl;
            }
        }
    }
    daf::log_memory("Other Index Memory");
    const bool compareMode = std::getenv("PIVOTER_COMPARE") != nullptr;
    const bool referenceOnlyMode = std::getenv("PIVOTER_RUN_REF") != nullptr;
    daf::timeCount("NucleusCoreDecomposition", [&] {
        if (r == 1 && compareMode) {
            // Compare r=1: optimized vs reference
            auto optTree = refTree.clone();
            auto optTreeGraphV = treeGraphV.clone();
            auto refCoreV = daf::timeCount("Reference r=1", [&]() {
                return NucleusCoreDecompositionCorrect(refTree, edgeGraph, treeGraphV, r, s);
            });
            auto optCoreV = daf::timeCount("Optimized r=1", [&]() {
                return NCliqueVertexCoreDecomposition(optTree, edgeGraph, optTreeGraphV, s);
            });
            // Compare: reference returns vector<pair<vector<Size>, int>>, optimized returns double*
            std::map<int, int> refDist, optDist;
            for (const auto &[clique, coreValue]: refCoreV) {
                refDist[coreValue]++;
            }
            for (daf::Size i = 0; i < edgeGraph.adj_list_offsets.size() - 1; ++i) {
                if (optCoreV[i] >= 0) optDist[(int)optCoreV[i]]++;
            }
            // Exclude core=0 from comparison (vertices with 0 nCr contribution)
            refDist.erase(0);
            optDist.erase(0);
            if (refDist == optDist) {
                std::cout << "✓ r=1 correctness verified: distributions match" << std::endl;
            } else {
                std::cerr << "✗ r=1 MISMATCH!" << std::endl;
                std::cerr << "Reference dist size: " << refDist.size() << " Optimized dist size: " << optDist.size() << std::endl;
                for (auto &[k, v] : refDist) {
                    if (optDist.count(k) == 0 || optDist[k] != v)
                        std::cerr << "  core=" << k << " ref=" << v << " opt=" << (optDist.count(k) ? optDist[k] : 0) << std::endl;
                }
                for (auto &[k, v] : optDist) {
                    if (refDist.count(k) == 0)
                        std::cerr << "  core=" << k << " ref=0 opt=" << v << " (extra in opt)" << std::endl;
                }
            }
            delete[] optCoreV;
        } else if (r == 2 && compareMode) {
            // Compare r=2: optimized vs reference
            auto optTree = refTree.clone();
            auto optTreeGraphV = treeGraphV.clone();
            auto refCore = daf::timeCount("Reference r=2", [&]() {
                return NucleusCoreDecompositionCorrect(refTree, edgeGraph, treeGraphV, r, s);
            });
            auto optCore = daf::timeCount("Optimized r=2", [&]() {
                return PlusNucleusEdgeCoreDecompositionSet(optTree, edgeGraph, optTreeGraphV, s);
            });
            // Compare distributions
            std::map<int, int> refDist, optDist;
            for (const auto &[clique, coreValue]: refCore) {
                refDist[coreValue]++;
            }
            for (const auto &[edge, coreValue]: optCore) {
                optDist[coreValue]++;
            }
            if (refDist == optDist) {
                std::cout << "✓ r=2 correctness verified: distributions match" << std::endl;
            } else {
                std::cerr << "✗ r=2 MISMATCH!" << std::endl;
                std::cerr << "Reference dist size: " << refDist.size() << " Optimized dist size: " << optDist.size() << std::endl;
                for (auto &[k, v] : refDist) {
                    if (optDist.count(k) == 0 || optDist[k] != v)
                        std::cerr << "  core=" << k << " ref=" << v << " opt=" << (optDist.count(k) ? optDist[k] : 0) << std::endl;
                }
            }
        } else if (r == 2) {
            PlusNucleusEdgeCoreDecompositionSet(refTree, edgeGraph, treeGraphV, s);
        } else if (r == 1) {
            NCliqueVertexCoreDecomposition(refTree, edgeGraph, treeGraphV, s);
        } else if (referenceOnlyMode) {
            auto res = NucleusCoreDecompositionCorrect(refTree, edgeGraph, treeGraphV, r, s);
            std::map<daf::Size, int> coreValueCount;
            for (const auto &[clique, coreValue]: res) {
                coreValueCount[coreValue]++;
            }
            std::cout << "Core value distribution (vertex degree in tree):" << std::endl;
            for (const auto &[coreValue, count]: coreValueCount) {
                std::cout << "Core value: " << coreValue << " Count: " << count << std::endl;
            }
        } else if (compareMode) {
            auto optTree = refTree.clone();
            // auto refTree = refTree.clone();
            auto refTreeGraphV = treeGraphV.clone();
            auto correctCore = daf::timeCount("Reference NucleusCoreDecomposition", [&]() {
                return NucleusCoreDecompositionCorrect(refTree, edgeGraph, refTreeGraphV, r, s);
            });

            auto optTreeGraphV = treeGraphV.clone();
            // auto optCore = NucleusCoreDecompositionRClique(optTree, edgeGraph, optTreeGraphV, r, s);
            auto optCore = daf::timeCount("Optimized NucleusCoreDecomposition", [&]() {
                return NucleusCoreDecompositionRClique(optTree, edgeGraph, optTreeGraphV, r, s);
            });
            std::map<daf::Size, int> correctCoreValueCount;
            for (const auto &[clique, coreValue]: correctCore) {
                correctCoreValueCount[coreValue]++;
            }
            std::map<daf::Size, int> optCoreValueCount;
            for (const auto &[clique, coreValue]: optCore) {
                optCoreValueCount[coreValue]++;
            }
            std::vector<int> corrDist, optDist;
            for (const auto &[coreValue, count]: correctCoreValueCount) {
                corrDist.push_back(count);
            }
            for (const auto &[coreValue, count]: optCoreValueCount) {
                optDist.push_back(count);
            }

            std::sort(corrDist.begin(), corrDist.end());
            std::sort(optDist.begin(), optDist.end());

            if (corrDist != optDist) {
                std::cerr << "Comparison failed: optimized core decomposition differs from reference." << std::endl;
                const auto mismatchCount = std::min(correctCore.size(), optCore.size());
                for (size_t i = 0; i < mismatchCount; ++i) {
                    if (correctCore[i] != optCore[i]) {
                        std::cerr << "First mismatch at index " << i << std::endl;
                        std::cerr << "Reference clique:";
                        for (auto v: correctCore[i].first) std::cerr << ' ' << v;
                        std::cerr << std::endl;
                        std::cerr << "Optimized clique:";
                        for (auto v: optCore[i].first) std::cerr << ' ' << v;
                        std::cerr << std::endl;
                        std::cerr << "Reference core: " << correctCore[i].second << " Optimized core: " << optCore[i].second << std::endl;
                        break;
                    }
                }
                std::abort();
            }
            std::cout << "Comparison passed: optimized result matches reference." << std::endl;
        } else {
            // r>=3 default: run optimized NucleusCoreDecompositionRClique
            NucleusCoreDecompositionRClique(refTree, edgeGraph, treeGraphV, r, s);
        }
    });

    daf::log_memory("Final Memory");

    // Print thread count and single/multi-threaded timing for r>=3


    // auto corePlus = daf::timeCount("NucleusCoreDecomposition", [&] {
    //     return PlusNucleusEdgeCoreDecompositionSet(treeGraph, edgeGraph, treeGraphV, s);
    // });

    // std::cout << "corePlus: " << corePlus << std::endl;
     //
     // auto coreBase = daf::timeCount("PlusNucleusEdgeCoreDecomposition", [&] {
     //     return baseNucleusCoreDecompositionLeaf(treeGraph, edgeGraph, treeGraphVSize, s);
     // });
//
//     std::ranges::sort(corePlus,
//                       [vertexMap](const auto &a, const auto &b) {
//                           if (a.second != b.second) return a.second < b.second;
//                           auto a_from = vertexMap[a.first.first];
//                           auto a_to = vertexMap[a.first.second];
//                           auto b_from = vertexMap[b.first.first];
//                           auto b_to = vertexMap[b.first.second];
//                           if (a_from != b_from) return a_from < b_from;
//                           if (a_to != b_to) return a_to < b_to;
//                           return a.second < b.second;
//     );
//     std::ranges::sort(coreBase,
//                       [vertexMap](const auto &a, const auto &b) {
//                           if (a.second != b.second) return a.second < b.second;
//                           auto a_from = vertexMap[a.first.first]//                       };
    //                           auto a_to = vertexMap[a.first.second];
//                           auto b_from = vertexMap[b.first.first];
//                           auto b_to = vertexMap[b.first.second];
//                           if (a_from != b_from) return a_from < b_from;
//                           if (a_to != b_to) return a_to < b_to;
//                           return a.second < b.second;
//                       }
//     );
//
// #ifndef NDEBUG
//     std::cout << "corePlus: " << corePlus << std::endl;
//     std::cout << "coreBase: " << coreBase << std::endl;
// #endif
//     for (auto i = 0; i < corePlus.size(); ++i) {
//         if (corePlus[i] != coreBase[i]) {
//             std::cout << "corePlus: " << corePlus[i].first.first << " " << corePlus[i].first.second << " "
//                     << corePlus[i].second << std::endl;
//             std::cout << "coreBase: " << coreBase[i].first.first << " " << coreBase[i].first.second << " "
//                     << coreBase[i].second << std::endl;
//             return 1;
//         }
//         coreV[corePlus[i].first.first] = std::max(coreV[corePlus[i].first.first], (double) corePlus[i].second);
//         coreV[corePlus[i].first.second] = std::max(coreV[corePlus[i].first.second], (double) corePlus[i].second);
//     }
    //  std::vector<double> coreV(edgeGraph.n);
    // for (daf::Size i = 0; i < edgeGraph.n; ++i) {
    //     coreV[i] = corePlus[i];
    //     // std::cout << "coreV[" << i << "]: " << coreV[i] << std::endl;
    // }
    //  std::ranges::sort(coreV);
    //  // ~/_/pivoter/a
    //  auto fileOutput = fopen("~/_/pivoter/a", "w");
    //  // for (const auto &i: corePlus) {
    //  //     fprintf(fileOutput, "%d %d %d\n", i.first.first, i.first.second, i.second);
    //  // }
    //  for (const auto &i: coreV) {
    //      fprintf(fileOutput, "%lf\n", i);
    //  }
    //
    // return 0;
}