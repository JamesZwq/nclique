//
// Created by 张文谦 on 25-3-18.
//e

#include<assert.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<limits.h>
#include<unistd.h>
#include<libgen.h>
#include <graph/Graph.h>
#include <tree/NCliqueCoreDecomposition.h>
#include <tree/NucleusCoreDecomposition.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include "graph/DynamicGraph.h"
#include "tree/BronKerbosch.h"

// 使用模板和可变参数

int main(int argc, char **argv) {
    // usage: ./main treeFile r s
    if (argc < 4) {
        printf("Incorrect number of arguments.\n");
        printf("./main <treeFile> <r> <s>\n");
        printf("treeFile: path to file\n");
        printf("r: r\n");
        printf("s: s\n");
        return 0;
    }
    const char *fpath = argv[1];
    const daf::CliqueSize r = strtol(argv[2], nullptr, 10);
    const char *graphFile;
    if (r == 2) {
        if (argc < 5) {
            printf("Incorrect number of arguments.\n");
            printf("./main <treeFile> <r> <s> <graphFile>\n");
            printf("treeFile: path to file\n");
            printf("r: r\n");
            printf("s: s\n");
            printf("graphFile: path to graph file if r == 2\n");
            return 0;
        }
        graphFile = argv[4];
    }
    const daf::CliqueSize s = strtol(argv[3], nullptr, 10);
    std::cout << fpath << " " << r << " " << s << std::endl;
    populate_nCr();



    const daf::CliqueSize minK = s;
    MultiBranchTree tree = MultiBranchTree::deserialize(fpath);
    const daf::Size numLeaf = tree.initLeafsParentAndId(minK);
    std::cout << "numLeaf: " << numLeaf << std::endl;
    std::cout << "nun of nodes: " << tree.getRoot()->children.size() << std::endl;

    daf::StaticVector<TreeNode *> leafList = tree.getLeafsList(numLeaf, minK);

// if debug
#ifndef NDEBUG
    tree.printTree();
    MultiBranchTree::printCliques(leafList, minK);
#endif

    if (r == 2) {
        DynamicGraph<TreeGraphNode> treeGraph = DynamicGraph<TreeGraphNode>(leafList, minK);
        Graph edgeGraph = Graph(graphFile, true);
        DynamicGraph<daf::Size> treeGraphV(leafList, edgeGraph.getGraphNodeSize(), minK);
        daf::vListMap.resize(edgeGraph.n + 1);

#ifndef NDEBUG
        treeGraphV.printGraphPerV();
#endif
        PlusNucleusEdgeCoreDecomposition(treeGraph, edgeGraph, treeGraphV, minK);
        leafList.free();
        daf::vListMap.free();
        return 0;
    }



    // std::cout << "avg degree: " << treeGraphV.getAvgDegree() << std::endl;
    // daf::StaticVector<double> list = tree.cliqueCount();
    //
    //
    // double uCount = tree.unionCliqueCount(treeGraphV, minK);
    // std::cout << "list: " << list << std::endl;
    // std::cout << "uCount: " << uCount << std::endl;
    // std::cout << "difference: " << list[minK] - uCount << " " << (list[minK] - uCount) / list[minK] << "%" << std::endl;
    // // print the number of vertices that degree == 1
    // tree.unionCliqueList(treeGraphV, minK);
    // treeGraphV.printGraphPerV();
    // //
    // MultiBranchTree::printCliques(leafList, minK);
    // // tree中每一个leaf的nbr
    // const Graph leafGraph = daf::timeCount("cliqueGraph", [&] {
    //      return Graph(tree, treeGraphV, leafList, r, s);
    // });
    // leafGraph.printGraphPerV();
    // std::cout << "avg degree cliqueGraph: " << treeGraphV.getAvgDegree() << std::endl;

    // std::vector<std::vector<daf::Size>> maxCliques = daf::timeCount("BronKerboschPivot", [&] {
    //     return cliqueGraph.BronKerboschPivot();
    // });
    // MultiBranchTree::printCliques(leafList, minK);
    // // sort maxCliques
    // for (auto &clique: maxCliques) {
    //     std::ranges::sort(clique);
    // }
    // std::cout << "maxCliques: " << maxCliques << std::endl;
    // Graph tmpG = Graph("/Users/zhangwenqian/UNSW/KClique/newSmallGraph.txt");
    // tmpG.printGraphPerV();
    // daf::timeCount("BronKerboschPivot", [&] {
    //     tmpG.BronKerboschPivot();
    // });


    // NucleusCoreDecomposition(tree, leafList, treeGraphV, leafGraph, r, s);

    leafList.free();
    // list.free();
}