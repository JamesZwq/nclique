#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#include<unistd.h>
#include<libgen.h>
#include <boost/version.hpp>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include "graph/DynamicBipartiteGraph.hpp"
#include "tree/NCliqueCoreDecomposition.h"

int main(int argc, char **argv) {
    std::cout << "Boost version: " << BOOST_LIB_VERSION << std::endl;
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
    memset(daf::vListMap.data, std::numeric_limits<daf::Size>::max(), edgeGraph.n * sizeof(daf::Size));
    auto vertexMap = edgeGraph.sortByDegeneracyOrder();

    DynamicGraph<TreeGraphNode> treeGraph = daf::timeCount("Tree Build", [&] {
        return listAllCliquesDegeneracy_VedgeGraph(edgeGraph, s, s);
    });
    edgeGraph.beSingleEdge();
    edgeGraph.buildEdgeIdMap();
    // DynamicBipartiteGraph BGraph(treeGraph, edgeGraph);


    std::cout << "nun Leaf: " << treeGraph.adj_list.size() << std::endl;
    DynamicGraph<TreeGraphNode> treeGraphV(treeGraph, edgeGraph.getGraphNodeSize(), s);

    // auto core = baseNucleusEdgeCoreDecomposition(treeGraph, edgeGraph, treeGraphV, s);


    auto treeGraphClone = treeGraph.clone();
    DynamicGraph<daf::Size> treeGraphVSize(treeGraph, edgeGraph.getGraphNodeSize(), s);


    auto corePlus = daf::timeCount("Plus Core Decomposition", [&] {
        return PlusNucleusEdgeCoreDecomposition(treeGraphClone, edgeGraph, treeGraphV, s);
    });


    auto coreBase = daf::timeCount("Base Core Decomposition", [&] {
        return baseNucleusEdgeCoreDecomposition(treeGraph, edgeGraph, treeGraphVSize, s);
    });

    std::ranges::sort(corePlus,
                      [vertexMap](const auto &a, const auto &b) {
                          if (a.second != b.second) return a.second < b.second;
                          auto a_from = vertexMap[a.first.first];
                          auto a_to = vertexMap[a.first.second];
                          auto b_from = vertexMap[b.first.first];
                          auto b_to = vertexMap[b.first.second];
                          if (a_from != b_from) return a_from < b_from;
                          if (a_to != b_to) return a_to < b_to;
                          return a.second < b.second;
                      }
    );
    std::ranges::sort(coreBase,
                      [vertexMap](const auto &a, const auto &b) {
                          if (a.second != b.second) return a.second < b.second;
                          auto a_from = vertexMap[a.first.first];
                          auto a_to = vertexMap[a.first.second];
                          auto b_from = vertexMap[b.first.first];
                          auto b_to = vertexMap[b.first.second];
                          if (a_from != b_from) return a_from < b_from;
                          if (a_to != b_to) return a_to < b_to;
                          return a.second < b.second;
                      }
    );

#ifndef NDEBUG
    std::cout << "corePlus: " << corePlus << std::endl;
    std::cout << "coreBase: " << coreBase << std::endl;
#endif
    for (auto i = 0; i < corePlus.size(); ++i) {
        if (corePlus[i] != coreBase[i]) {
            std::cout << "corePlus: " << corePlus[i].first.first << " " << corePlus[i].first.second << " "
                    << corePlus[i].second << std::endl;
            std::cout << "coreBase: " << coreBase[i].first.first << " " << coreBase[i].first.second << " "
                    << coreBase[i].second << std::endl;
            return 1;
        }
    }


    return 0;
}
