//
// Created by 张文谦 on 25-3-4.
//

#include "../tree/NCliqueCoreDecomposition.h"

#include "dataStruct/StaticVector.h"


extern double nCr[1001][401];
namespace baseCD{
void countingPerVertexHelp(const TreeNode &node,
                           const daf::CliqueSize k,
                           daf::Size *core,
                           StaticVector<daf::Size> &povit,
                           StaticVector<daf::Size> &keepC
) {
    daf::Size cliqueSize = povit.size() + keepC.size();
    if (node.children.empty() && cliqueSize >= k && keepC.size() <= k) {
        const int needPivot = k - keepC.size(); // 还需从 pivot 中选的顶点数
        daf::Size totalKcliques = 0;
        if (needPivot >= 0 && needPivot <= povit.size()) {
            totalKcliques = nCr[povit.size()][needPivot];
        }
        for (const auto v: keepC) {
            // if (v == 0) {
                // std::cout << "keep: " << v << " totalKcliques: " << totalKcliques << std::endl;
            // }
            core[v] += totalKcliques;
        }

        daf::Size eachPivotKcliques = 0;
        const int needPivotWithV = needPivot - 1;
        if (needPivotWithV >= 0 && needPivotWithV <= povit.size() - 1) {
            eachPivotKcliques = nCr[povit.size() - 1][needPivotWithV];
        }

        for (auto v: povit) {
            // std::cout << "povit: " << v << " eachPivotKcliques: " << eachPivotKcliques << std::endl;
            core[v] += eachPivotKcliques;
        }

        return;
    }

    for (const auto &child: node.children) {
        if (child->MaxDeep < k) {
            continue;
        }
        if (child->isPivot) {
            povit.push_back(child->v);
            countingPerVertexHelp(*child, k, core, povit, keepC);
            povit.pop_back();
        } else {
            keepC.push_back(child->v);
            countingPerVertexHelp(*child, k, core, povit, keepC);
            keepC.pop_back();
        }
    }
}

daf::Size *countingPerVertex(const MultiBranchTree &tree, const daf::CliqueSize k) {
    auto *core = new daf::Size[tree.getRoot()->children.size()];
    //init 0
    std::memset(core, 0, tree.getRoot()->children.size() * sizeof(daf::Size));
    StaticVector<daf::Size> povitC;
    StaticVector<daf::Size> keepC;
    for (auto node: tree.getRoot()->children) {
        if (node->MaxDeep < k) {
            continue;
        }
        keepC.push_back(node->v);
        countingPerVertexHelp(*node, k, core, povitC, keepC);
        keepC.pop_back();
    }
    return core;
}



void computeSupHelp(const TreeNode &node,
                           const daf::CliqueSize k,
                           const daf::Size *core,
                           StaticVector<daf::Size> &povit,
                           StaticVector<daf::Size> &keepC,
                           std::vector<std::map<daf::Size, daf::Size, std::greater<>>> &sup
) {
    daf::Size cliqueSize = povit.size() + keepC.size();
    if (node.children.empty() && cliqueSize >= k && keepC.size() <= k) {
        assert(!keepC.empty());

        StaticVector<daf::Size> orderPovit = povit;
        // sort from large to small
        std::ranges::sort(orderPovit, [core](daf::Size a, daf::Size b) {
            return core[a] > core[b];
        });
        const int needPivot = k - keepC.size(); // 还需从 pivot 中选的顶点数


        daf::Size maxK = core[*std::ranges::min_element(keepC, [core](daf::Size a, daf::Size b) {
                                                return core[a] < core[b];
                                            })];
        if (needPivot == 0) {
            for (const auto v: keepC) {
                sup[v][maxK]++;
            }
            return;
        }

        maxK = std::min(maxK, core[orderPovit[needPivot - 1]]);

        int numGraterThanMaxKinKeepC = needPivot;
        for (daf::Size i = needPivot; i < orderPovit.size(); i++) {
            auto v = orderPovit[i];
            if (core[v] >= maxK) {
                numGraterThanMaxKinKeepC++;
            } else {
                break;
            }
        }

        daf::Size prveCore = maxK;
        daf::Size prvePovit = 0;
        daf::Size prveKeepCliqueCount = 0;
        daf::Size prvepovitCliqueCount = 0;
        for (daf::Size i = numGraterThanMaxKinKeepC; i <= orderPovit.size(); i++) {
            if (i == orderPovit.size() || core[orderPovit[i]] < prveCore) {
                const daf::Size keepCliqueCount = nCr[i][needPivot];
                const daf::Size povitCliqueCount = nCr[i - 1][needPivot - 1];
                for (const auto v: keepC) {
                    sup[v][prveCore] += keepCliqueCount - prveKeepCliqueCount;
                }
                for (daf::Size v = 0; v < prvePovit; v++) {
                    sup[orderPovit[v]][prveCore] += povitCliqueCount - prvepovitCliqueCount;
                }
                for (daf::Size v = prvePovit; v < i; v++) {
                    sup[orderPovit[v]][prveCore] += povitCliqueCount;
                }
                if (i == orderPovit.size()) break;
                prveCore = core[orderPovit[i]];
                prvepovitCliqueCount = povitCliqueCount;
                prveKeepCliqueCount = keepCliqueCount;
                prvePovit = i;
            }
        }

        return;
    }

    for (const auto &child: node.children) {
        if (child->MaxDeep < k) {
            continue;
        }
        if (child->isPivot) {
            povit.push_back(child->v);
            computeSupHelp(*child, k, core, povit, keepC, sup);
            povit.pop_back();
        } else {
            keepC.push_back(child->v);
            computeSupHelp(*child, k, core, povit, keepC, sup);
            keepC.pop_back();
        }
    }
}

std::vector<std::map<daf::Size, daf::Size, std::greater<>>> computeSup(const MultiBranchTree &tree,
                const daf::CliqueSize k,
                const daf::Size *core) {
    StaticVector<daf::Size> povitC;
    StaticVector<daf::Size> keepC;
    std::vector<std::map<daf::Size, daf::Size, std::greater<>>> sup(tree.getRoot()->children.size());
    for (const auto node: tree.getRoot()->children) {
        if (node->MaxDeep < k) {
            continue;
        }
        keepC.push_back(node->v);
        computeSupHelp(*node, k, core, povitC, keepC, sup);
        keepC.pop_back();
    }
    return sup;
}
}

void baseNucleusCoreDecomposition(const MultiBranchTree &tree, daf::CliqueSize k) {
    // daf::Size numNodes = tree.getRoot()->children.size();
    // tree.printTree();
    auto time_start = std::chrono::high_resolution_clock::now();
    auto *core = baseCD::countingPerVertex(tree, k);
    // std::cout << "init core: ";
    // daf::printArray(core, tree.getRoot()->children.size());
    bool update = true;
    daf::Size iter = 0;
    while (update) {
        update = false;
        std::cout << ++iter << " " << std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
        time_start = std::chrono::high_resolution_clock::now();
        auto sup = baseCD::computeSup(tree, k, core);
        for (daf::Size v = 0; v < tree.getRoot()->children.size(); v++) {
            auto &supV = sup[v];
            std::cout << v << " supV: " << supV << std::endl;
            daf::Size prveCount = 0;
            for (auto &[c, count]: supV) {
                if (prveCount + count >= c) {
                    daf::Size newCore = std::max(c,prveCount);
                    if (core[v] > newCore) {
                        core[v] = newCore;
                        update = true;
                    }
                    break;
                }
                prveCount += count;
                // prveCount = count;
            }
        }
        // std::cout << "core: ";
        // daf::printArray(core, tree.getRoot()->children.size());
    }

    auto file = fopen("/Users/zhangwenqian/UNSW/pivoter/a", "w");
    std::sort(core, core + tree.getRoot()->children.size());
    for (daf::Size i = 0; i < tree.getRoot()->children.size(); i++) {
        fprintf(file, "%d\n", core[i]);
    }
    fclose(file);
    delete[] core;
}
