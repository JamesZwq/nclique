#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#include<unistd.h>
#include<libgen.h>
// #include <boost/version.hpp>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include "BK/BronKerboschRmEdge.hpp"
#include "BK/BronKerboschRmRClique.hpp"
#include "dataStruct/disJoinSet.hpp"
#include "graph/DynamicBipartiteGraph.hpp"
#include "NucleusDecomposition/NCliqueCoreDecomposition.h"



int main(int argc, char **argv) {
    // populate_nCr();
    // daf::vListMap.resize(100);
    // // bkRmEdge::bronKerboschFromFile("/Users/zhangwenqian/UNSW/pivoter/a.edge", 1,
    // //                                [](const bkRmEdge::Bitset &clique, const bkRmEdge::Bitset &pivots) {
    // //                                    std::cout << "Find clique: " << clique << " pivots: " << pivots << std::endl;
    // //                                    return true;
    // //                                });
    //
    // bkRmClique::testBronKerbosch();
    // return 0;


    // std::cout << "Boost version: " << BOOST_LIB_VERSION << std::endl;
    if (argc != 4) {
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

    printf("about to call runAndPrint for dataset %s\n", fpath);

    Graph edgeGraph(fpath);
    edgeGraph.printGraphInfo();
    populate_nCr();
    daf::vListMap.resize(edgeGraph.n + 1);
    // std::numeric_limits<daf::Size>::max();
    memset(daf::vListMap.data(), -1, edgeGraph.n * sizeof(daf::Size));
    // auto vertexMap = edgeGraph.sortByDegeneracyOrder();

    DynamicGraph<TreeGraphNode> treeGraph = daf::timeCount("Tree Build", [&] {
        return SDCT(edgeGraph, s, s);
    });
    std::cout << s << "-Clique count: "<< treeGraph.cliqueCount(2) << std::endl;
    std::cout << "max clique: " << treeGraph.maxDegree() << std::endl;
    // if (s >
    // treeGraph.printGraphPerV();
    for (auto leaf: treeGraph.adj_list) {
        if (leaf[0].isPivot) {
            std::cout << "leaf: " << leaf[0].v << " is pivot" << std::endl;
        }
    }

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

    daf::timeCount("NucleusCoreDecomposition", [&] {
        if (r == 2) {
            PlusNucleusEdgeCoreDecompositionSet(treeGraph, edgeGraph, treeGraphV, s);
        } else if (r == 1) {
            NCliqueVertexCoreDecomposition(treeGraph, edgeGraph, treeGraphV, s);
        } else {
            // NucleusCoreDecomposition(treeGraph, edgeGraph, treeGraphV, r, s);
            NucleusCoreDecompositionRClique(treeGraph, edgeGraph, treeGraphV, r, s);
            // NucleusCoreDecompositionHierarchy(treeGraph, edgeGraph, treeGraphV, r, s);
        }
    });
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
    //  // /Users/zhangwenqian/UNSW/pivoter/a
    //  auto fileOutput = fopen("/Users/zhangwenqian/UNSW/pivoter/a", "w");
    //  // for (const auto &i: corePlus) {
    //  //     fprintf(fileOutput, "%d %d %d\n", i.first.first, i.first.second, i.second);
    //  // }
    //  for (const auto &i: coreV) {
    //      fprintf(fileOutput, "%lf\n", i);
    //  }
    //
    // return 0;
}