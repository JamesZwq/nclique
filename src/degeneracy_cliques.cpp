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
    DynamicGraph<TreeGraphNode> treeGraph = daf::timeCount("Tree Build", [&]() -> DynamicGraph<TreeGraphNode> {
        return SDCT_Parallel(edgeGraph, 1000000, 0);  // 并行版本（已验证正确性和性能）
    });
    std::cout << "TreeGraph Clique Count: \n" << treeGraph.cliqueCount() << std::endl;
    // std::cout << "TreeGraphPerV Clique Count: \n" << treeGraphPar.cliqueCount() << std::endl;

    // daf::log_memory("Tree Memory");
    // std::cout << s << "-Clique count: "<< treeGraph.cliqueCount(s) << std::endl;
    // std::cout << "max clique: " << treeGraph.maxDegree() << std::endl;
    // // if (s >
    // treeGraph.printGraphPerV();
    // for (auto leaf: treeGraph.adj_list) {
    //     if (leaf[0].isPivot) {
    //         std::cout << "leaf: " << leaf[0].v << " is pivot" << std::endl;
    //     }
    // }

    // return 0;

    edgeGraph.initCore();
    // treeGraph.printGraphPerV();

    edgeGraph.beSingleEdge();
    edgeGraph.buildEdgeIdMap();
    // DynamicBipartiteGraph BGraph(treeGraph, edgeGraph);


    std::cout << "nun Leaf: " << treeGraph.adj_list.size() << std::endl;
    // DynamicGraphSet<TreeGraphNode> treeGraphV(treeGraph, edgeGraph.getGraphNodeSize(), s);
    DynamicGraphSet<TreeGraphNode> treeGraphV(treeGraph, edgeGraph.getGraphNodeSize(), s);


    // StaticCliqueIndex cliqueIndex(r);
    // daf::timeCount("clique Index build",
    //                [&]() {
    //                    cliqueIndex.build(treeGraph, edgeGraph.adj_list.size());
    //                });
    //
    // cliqueIndex.verify();
    // print all leaf with node 154
    for (daf::Size leafId = 0; leafId < treeGraph.adj_list.size(); ++leafId) {
        const auto &leaf = treeGraph.adj_list[leafId];
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
        if (r == 2) {
            PlusNucleusEdgeCoreDecompositionSet(treeGraph, edgeGraph, treeGraphV, s);
        } else if (r == 1) {
            NCliqueVertexCoreDecomposition(treeGraph, edgeGraph, treeGraphV, s);
        } else if (referenceOnlyMode) {
            NucleusCoreDecompositionRCliqueRef(treeGraph, edgeGraph, treeGraphV, r, s);
            std::map<int, int> coreValueCount;
            for (const auto &leaf: treeGraph.adj_list) {
                for (const auto &node: leaf) {
                    if (!node.isPivot) {
                        coreValueCount[node.v]++;
                    }
                }
            }
            std::cout << "Core value distribution (vertex degree in tree):" << std::endl;
            for (const auto &[coreValue, count]: coreValueCount) {
                std::cout << "Core value: " << coreValue << " Count: " << count << std::endl;
            }
        } else if (compareMode) {
            auto refTree = treeGraph.clone();
            auto refTreeGraphV = treeGraphV.clone();
            // auto refCore = NucleusCoreDecompositionRCliqueRef(refTree, edgeGraph, refTreeGraphV, r, s);
            auto refCore = daf::timeCount("Reference NucleusCoreDecomposition", [&]() {
                return NucleusCoreDecompositionRCliqueRef(refTree, edgeGraph, refTreeGraphV, r, s);
            });

            auto optTree = treeGraph.clone();
            auto optTreeGraphV = treeGraphV.clone();
            // auto optCore = NucleusCoreDecompositionRClique(optTree, edgeGraph, optTreeGraphV, r, s);
            auto optCore = daf::timeCount("Optimized NucleusCoreDecomposition", [&]() {
                return NucleusCoreDecompositionRClique(optTree, edgeGraph, optTreeGraphV, r, s);
            });

            auto canonicalLess = [](const auto &a, const auto &b) {
                if (a.second != b.second) {
                    return a.second < b.second;
                }
                return a.first < b.first;
            };
            std::sort(refCore.begin(), refCore.end(), canonicalLess);
            std::sort(optCore.begin(), optCore.end(), canonicalLess);

            if (refCore != optCore) {
                std::cerr << "Comparison failed: optimized core decomposition differs from reference." << std::endl;
                const auto mismatchCount = std::min(refCore.size(), optCore.size());
                for (size_t i = 0; i < mismatchCount; ++i) {
                    if (refCore[i] != optCore[i]) {
                        std::cerr << "First mismatch at index " << i << std::endl;
                        std::cerr << "Reference clique:";
                        for (auto v: refCore[i].first) std::cerr << ' ' << v;
                        std::cerr << std::endl;
                        std::cerr << "Optimized clique:";
                        for (auto v: optCore[i].first) std::cerr << ' ' << v;
                        std::cerr << std::endl;
                        std::cerr << "Reference core: " << refCore[i].second << " Optimized core: " << optCore[i].second << std::endl;
                        break;
                    }
                }
                std::abort();
            }
            std::cout << "Comparison passed: optimized result matches reference." << std::endl;
        } else {
            // 使用优化的 Ref 版本
            std::cout << "Using optimized Ref version" << std::endl;
            NucleusCoreDecompositionRCliqueRef(treeGraph, edgeGraph, treeGraphV, r, s);
        }
    });

    daf::log_memory("Final Memory");

    // Print thread count and single/multi-threaded timing for r>=3
    if (r >= 3) {
#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
        
        // 测试不同线程数的性能
        std::vector<int> thread_counts = {1, 8, 16, 32, 64};
        std::vector<long long> times;
        
        std::cout << "\n=== Multi-threading Performance Test ===" << std::endl;
        std::cout << "Testing with different thread counts..." << std::endl;
        
        for (int threads : thread_counts) {
            auto treeCopy = treeGraph.clone();
            auto treeGraphVCopy = treeGraphV.clone();
            omp_set_num_threads(threads);
            
            auto t1 = std::chrono::high_resolution_clock::now();
            NucleusCoreDecompositionRCliqueRef(treeCopy, edgeGraph, treeGraphVCopy, r, s);
            auto t2 = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
            times.push_back(elapsed);
            
            std::cout << "Threads: " << threads << " -> Time: " << elapsed << " ms";
            if (threads > 1) {
                double speedup = (double)times[0] / elapsed;
                std::cout << " (speedup: " << speedup << "x)";
            }
            std::cout << std::endl;
        }
        
        // 恢复原来的线程数
        omp_set_num_threads(nthreads);
#else
        std::cout << "--- Timing Summary ---" << std::endl;
        std::cout << "Thread count: 1 (OpenMP not enabled)" << std::endl;
        std::cout << "Single-threaded time: N/A" << std::endl;
        std::cout << "Multi-threaded time: N/A" << std::endl;
#endif
    }

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