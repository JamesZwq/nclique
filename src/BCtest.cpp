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
#include <NucleusDecomposition/NCliqueCoreDecomposition.h>
#include <NucleusDecomposition/NucleusCoreDecomposition.h>

#include"misc.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include "graph/DynamicGraph.h"
#include "BK/BronKerboschRmEdge.hpp"
// #include "nCr.h"

// 使用模板和可变参数
extern double nCr[1001][401];

int main() {
    populate_nCr();
    using Node = TreeGraphNode;
    // -------- ❶ 手动拼一个小图 ------------
    // /Users/zhangwenqian/UNSW/pivoter/delEdges
    // std::vector<std::pair<int,int>> edges;
    daf::StaticVector<std::pair<daf::Size, daf::Size> > delEdges;
    // 删掉 {0,1} 这一条边
    auto file = fopen("/Users/zhangwenqian/UNSW/pivoter/b", "r");
    int u, v;
    int cliqueSize;
    //frist line is size
    if (fscanf(file, "%d", &cliqueSize) != 1) {
        std::cerr << "Error reading clique size" << std::endl;
        return 1;
    }
    while (fscanf(file, "%d %d", &u, &v) == 2) {
        delEdges.emplace_back(u, v);
    }
    fclose(file);
    // std::cout << "edges: " << delEdges << std::endl;

    std::vector<Node> vList;
    daf::vListMap.resize(cliqueSize + 1);
    vList.reserve(cliqueSize);
    for (int i = 0; i < cliqueSize; ++i) {
        vList.emplace_back(i, true);
    }

    // -------- ❷ 跑 bronKerbosch ----------
    std::vector<double> cliqueCount;
    cliqueCount.resize(cliqueSize + 1);
    memset(cliqueCount.data(), 0, (cliqueSize + 1) * sizeof(double));
    std::vector<std::vector<daf::Size>> cliques;
    bkRmEdge::bronKerbosch(
        vList, delEdges, /*minK=*/1,
        [&](const bkRmEdge::Bitset &R, const bkRmEdge::Bitset &piv) {
            auto clique = bkRmEdge::coverToVertex(R, piv, vList);
            std::cout << "clique: " << clique << " Size: " << clique.size() << std::endl;

            int numPivots = 0;
            int numKeep = 0;
            for (int i = 0; i < clique.size(); i++) {
                if (clique[i].isPivot) {
                    numPivots++;
                } else {
                    numKeep++;
                }
            }
            int rSize = numPivots + numKeep;
            for (int i = numPivots; i >= 0; i--) {
                int k = rSize - i;
                cliqueCount[k] += nCr[numPivots][i];
                // std::cout << "cliqueCount[" << k << "] = " << cliqueCount[k] << std::endl;
            }
            // show all k clique, add to edges
            int k = 6;
            daf::StaticVector<daf::Size> keepIdx, pivIdx;
            for (int i = 0; i < clique.size(); ++i) {
                if (clique[i].isPivot) {
                    // pivot
                    pivIdx.push_back(i);
                } else {
                    // keep
                    keepIdx.push_back(i);
                }
            }
            daf::enumerateCombinations(keepIdx,
                                       pivIdx,
                                       k,
                                       [&](daf::StaticVector<daf::Size> &keep, daf::StaticVector<daf::Size> &combination) {
                                            std::vector<daf::Size> curr_clique;
                                           for (auto &i: keep) {
                                               curr_clique.push_back(clique[i]);
                                           }
                                             for (auto &i: combination) {
                                                  curr_clique.push_back(clique[i]);
                                             }
                                           if (curr_clique.size() == k) {
                                               // std::cout << "clique: " << clique << std::endl;
                                               std::ranges::sort(curr_clique);
                                               cliques.emplace_back(curr_clique);
                                           }
                                           return true;
                                       });
        });
    std::cout << "num cliques: " << cliques.size() << std::endl;
    // sort
    std::ranges::sort(cliques, [](const std::vector<daf::Size> &a, const std::vector<daf::Size> &b) {
        for (int i = 0; i < a.size() && i < b.size(); i++) {
            if (a[i] != b[i]) {
                return a[i] < b[i];
            }
        }
        return a.size() < b.size();
    });
    // check duplicates
    for (size_t i = 1; i < cliques.size(); i++) {
        if (cliques[i] == cliques[i - 1]) {
            std::cout << "duplicate clique: ";
            for (auto v : cliques[i]) std::cout << v << ' ';
            std::cout << std::endl;
        }
    }
    // daf::vListMap.print();
    //remove ending 0
    for (int i = cliqueCount.size() - 1; i >= 0; i--) {
        if (cliqueCount[i] == 0) {
            cliqueCount.pop_back();
        } else {
            break;
        }
    }
    // print as lf
    cliqueCount[0] = 0;
    auto file1 = fopen("/Users/zhangwenqian/UNSW/pivoter/outA.txt", "w");
    for (auto &i: cliqueCount) {
        // printf("%lf\n", i);
        fprintf(file1, "%lf\n", i);
    }
    fclose(file1);

    // sort edges
    // std::sort(edges.begin(), edges.end());
    // for (auto & i: edges) {
    //     std::cout << i.first << " " << i.second << std::endl;
    // }
    return 0;
}
