//
// Created by _ on 25-3-4.
//

#include <graph/Graph.h>

#include "NCliqueCoreDecomposition.h"


extern double nCr[1001][401];

namespace baseCDLeaf {
    double *countingPerVertex(const MultiBranchTree &tree, const daf::CliqueSize k,
                              const daf::StaticVector<TreeNode *> &leafs) {
        auto *core = new double[tree.getRoot()->children.size()];
        //init 0
        std::memset(core, 0, tree.getRoot()->children.size() * sizeof(daf::Size));
        daf::StaticVector<daf::Size> povitC{};
        daf::StaticVector<daf::Size> keepC{};
        int count = 0;
        for (const auto leaf: leafs) {
            // std::cout << "leaf: " << count++ << std::endl;
            if (leaf->MaxDeep < k) continue;
            povitC.clear();
            keepC.clear();
            auto node = leaf;
            while (node != nullptr) {
                if (node->v == ROOTID) {
                    if (keepC.size() > k) {
                        break;
                    }
                    const int needPivot = k - keepC.size(); //  pivot 
                    double totalKcliques = 0;
                    if (needPivot >= 0 && needPivot <= povitC.size()) {
                        totalKcliques = nCr[povitC.size()][needPivot];
                    }
                    for (const auto v: keepC) {
                        core[v] += totalKcliques;
                    }

                    double eachPivotKcliques = 0;
                    const int needPivotWithV = needPivot - 1;
                    if (needPivotWithV >= 0 && needPivotWithV <= povitC.size() - 1) {
                        eachPivotKcliques = nCr[povitC.size() - 1][needPivotWithV];
                    }

                    for (auto v: povitC) {
                        core[v] += eachPivotKcliques;
                    }
                    break;
                }
                if (node->isPivot) povitC.push_back(node->v);
                else keepC.push_back(node->v);
                node = node->parent;
            }
        }
        povitC.free();
        keepC.free();
        return core;
    }

    std::vector<std::map<double, double, std::greater<> > > computeSupLeaf(const MultiBranchTree &tree,
                                                                           const daf::CliqueSize k,
                                                                           const double *core,
                                                                           const daf::StaticVector<TreeNode *> &leafs) {
        daf::StaticVector<daf::Size> povitC{};
        daf::StaticVector<daf::Size> keepC{};
        std::vector<std::map<double, double, std::greater<> > > sup(tree.getRoot()->children.size());
        for (const auto leaf: leafs) {
            if (leaf->MaxDeep < k) continue;
            povitC.clear();
            keepC.clear();
            auto node = leaf;
            while (node != nullptr) {
                if (node->v == ROOTID) {
                    assert(!keepC.empty());

                    // sort from large to small
                    std::ranges::sort(povitC, [core](daf::Size a, daf::Size b) {
                        return core[a] > core[b];
                    });
                    const daf::Size needPivot = k - keepC.size(); //  pivot 


                    double maxK = core[*std::ranges::min_element(keepC, [core](daf::Size a, daf::Size b) {
                        return core[a] < core[b];
                    })];
                    if (needPivot == 0) {
                        for (const auto v: keepC) {
                            sup[v][maxK]++;
                        }
                        break;
                    }

                    maxK = std::min(maxK, core[povitC[needPivot - 1]]);

                    daf::Size numGraterThanMaxKinKeepC = needPivot;
                    for (daf::Size i = needPivot; i < povitC.size(); i++) {
                        auto v = povitC[i];
                        if (core[v] >= maxK) {
                            numGraterThanMaxKinKeepC++;
                        } else {
                            break;
                        }
                    }

                    double prveCore = maxK;
                    daf::Size prvePovit = 0;
                    double prveKeepCliqueCount = 0;
                    double prvepovitCliqueCount = 0;
                    for (daf::Size i = numGraterThanMaxKinKeepC; i <= povitC.size(); i++) {
                        if (i == povitC.size() || core[povitC[i]] < prveCore) {
                            const double keepCliqueCount = nCr[i][needPivot];
                            const double povitCliqueCount = nCr[i - 1][needPivot - 1];
                            for (const auto v: keepC) {
                                sup[v][prveCore] += keepCliqueCount - prveKeepCliqueCount;
                            }
                            for (daf::Size v = 0; v < prvePovit; v++) {
                                sup[povitC[v]][prveCore] += povitCliqueCount - prvepovitCliqueCount;
                            }
                            for (daf::Size v = prvePovit; v < i; v++) {
                                sup[povitC[v]][prveCore] += povitCliqueCount;
                            }
                            if (i == povitC.size()) break;
                            prveCore = core[povitC[i]];
                            prvepovitCliqueCount = povitCliqueCount;
                            prveKeepCliqueCount = keepCliqueCount;
                            prvePovit = i;
                        }
                    }
                }
                if (node->isPivot) povitC.push_back(node->v);
                else keepC.push_back(node->v);
                node = node->parent;
            }
        }
        povitC.free();
        keepC.free();
        return sup;
    }
}

double * baseNucleusCoreDecompositionLeaf(const MultiBranchTree &tree, daf::CliqueSize k) {
    // daf::Size numNodes = tree.getRoot()->children.size();
    // tree.printTree();
    auto time_start = std::chrono::high_resolution_clock::now();

    const auto numLeaf = tree.initLeafsParentAndId();
    std::cout << "numLeaf: " << numLeaf << std::endl;
    daf::StaticVector<TreeNode *> leafList = tree.getLeafsList(numLeaf);
    auto *core = baseCDLeaf::countingPerVertex(tree, k, leafList);
    // daf::printArray(core, tree.getRoot()->children.size());
    bool update = true;
    daf::Size iter = 0;
    while (update) {
        update = false;
        std::cout << ++iter << " " << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
        time_start = std::chrono::high_resolution_clock::now();
        auto sup = baseCDLeaf::computeSupLeaf(tree, k, core, leafList);
        for (daf::Size v = 0; v < tree.getRoot()->children.size(); v++) {
            auto &supV = sup[v];
            // std::cout << v << " supV: " << supV << std::endl;
            double prveCount = 0;
            for (auto &[c, count]: supV) {
                if (prveCount + count >= c) {
                    double newCore = std::max(c, prveCount);
                    if (core[v] > newCore) {
                        core[v] = newCore;
                        update = true;
                    }
                    break;
                }
                prveCount += count;
            }
        }
    }

    // auto file = fopen("~/_/pivoter/b", "w");
    // std::sort(core, core + tree.getRoot()->children.size());
    // for (daf::Size i = 0; i < tree.getRoot()->children.size(); i++) {
    //     fprintf(file, "%f\n", core[i]);
    // }
    // fclose(file);
    delete[] core;
    leafList.free();

    return core;
}
