//
// Created by 张文谦 on 25-3-25.
//

#include "NucleusCoreDecomposition.h"
extern double nCr[1001][401];

namespace NucleusCD {
    void countingPerVertexHelp(const TreeNode &node,
                               const daf::CliqueSize k,
                               double *core,
                               daf::StaticVector<daf::Size> &povit,
                               daf::StaticVector<daf::Size> &keepC
    ) {
        daf::Size cliqueSize = povit.size() + keepC.size();
        if (node.children.empty() && cliqueSize >= k && keepC.size() <= k) {
            const int needPivot = k - keepC.size(); // 还需从 pivot 中选的顶点数
            double totalKcliques = 0;
            if (needPivot >= 0 && needPivot <= povit.size()) {
                totalKcliques = nCr[povit.size()][needPivot];
            }
            for (const auto v: keepC) {
                core[v] += totalKcliques;
            }

            double eachPivotKcliques = 0;
            const int needPivotWithV = needPivot - 1;
            if (needPivotWithV >= 0 && needPivotWithV <= povit.size() - 1) {
                eachPivotKcliques = nCr[povit.size() - 1][needPivotWithV];
            }

            for (auto v: povit) {
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

    double *countingPerVertex(const MultiBranchTree &tree, const daf::CliqueSize k) {
        auto *core = new double[tree.getRoot()->children.size()];
        //init 0
        std::memset(core, 0, tree.getRoot()->children.size() * sizeof(daf::Size));
        daf::StaticVector<daf::Size> povitC;
        daf::StaticVector<daf::Size> keepC;
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
                        const daf::CliqueSize r,
                        const daf::CliqueSize s,
                        const double *core,
                        daf::StaticVector<daf::Size> &povit,
                        daf::StaticVector<daf::Size> &keepC,
                        std::vector<std::map<double, double, std::greater<> > > &sup
    ) {
        daf::Size cliqueSize = povit.size() + keepC.size();
        if (node.children.empty() && cliqueSize >= s && keepC.size() <= s) {
            assert(!keepC.empty());

            daf::StaticVector<daf::Size> orderPovit = povit.deepCopy();
            // sort from large to small
            std::ranges::sort(orderPovit, [core](daf::Size a, daf::Size b) {
                return core[a] > core[b];
            });
            const daf::Size needPivot = s - keepC.size(); // 还需从 pivot 中选的顶点数


            double maxK = core[*std::ranges::min_element(keepC, [core](daf::Size a, daf::Size b) {
                return core[a] < core[b];
            })];

            // if s == keepC, then for each r clique, only in one s clique
            if (needPivot == 0) {
                for (const auto v: keepC) {
                    sup[v][maxK]++;
                }
                return;
            }

            maxK = std::min(maxK, core[orderPovit[needPivot - 1]]);

            daf::Size numGraterThanMaxKinKeepC = needPivot;
            for (daf::Size i = needPivot; i < orderPovit.size(); i++) {
                auto v = orderPovit[i];
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
            for (daf::Size i = numGraterThanMaxKinKeepC; i <= orderPovit.size(); i++) {
                if (i == orderPovit.size() || core[orderPovit[i]] < prveCore) {
                    const double keepCliqueCount = nCr[i - r + 1][needPivot - r + 1];
                    const double povitCliqueCount = nCr[i - 1 - r + 1][needPivot - 1 - r + 1];
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
            orderPovit.free();
            return;
        }

        for (const auto &child: node.children) {
            if (child->MaxDeep < s) {
                continue;
            }
            if (child->isPivot) {
                povit.push_back(child->v);
                // std::cout << "Apovit: " << povit << std::endl;
                computeSupHelp(*child, r, s, core, povit, keepC, sup);
                // std::cout << "Bpovit: " << povit << std::endl;
                povit.pop_back();
            } else {
                keepC.push_back(child->v);
                computeSupHelp(*child, r, s, core, povit, keepC, sup);
                keepC.pop_back();
            }
        }
    }

    std::vector<std::map<double, double, std::greater<> > > computeSupOld(const MultiBranchTree &tree,
                                                                          const daf::CliqueSize r,
                                                                          const daf::CliqueSize s,
                                                                          const double *core) {
        daf::StaticVector<daf::Size> povitC;
        daf::StaticVector<daf::Size> keepC;
        std::vector<std::map<double, double, std::greater<> > > sup(tree.getRoot()->children.size());
        for (const auto node: tree.getRoot()->children) {
            if (node->MaxDeep < s) {
                continue;
            }
            keepC.push_back(node->v);
            computeSupHelp(*node, r, s, core, povitC, keepC, sup);
            keepC.pop_back();
        }
        return sup;
    }

    double coreEvaluation(
        const daf::StaticVector<daf::Size> &intersectIds,
        const daf::StaticVector<TreeNode *> &leafList,
        const daf::CliqueSize r,
        const daf::CliqueSize s,
        const daf::Size currV,
        const double *core) {
        daf::StaticVector<daf::Size> povitC;
        daf::StaticVector<daf::Size> keepC;
        std::map<double, double, std::greater<> > sup;
        for (const auto leafId: intersectIds) {
            auto leaf = leafList[leafId];
            if (leaf->MaxDeep < s) continue;
            povitC.clear();
            keepC.clear();
            auto node = leaf;
            double minKeepCore = std::numeric_limits<double>::max();
            bool isPivotVertex = false;
            while (!node->isRoot()) {
                if (node->isPivot) {
                    povitC.push_back(node->v);
                } else {
                    keepC.push_back(node->v);
                    minKeepCore = std::min(minKeepCore, core[node->v]);
                }
                if (node->v == currV) {
                    isPivotVertex = node->isPivot;
                }
                node = node->parent;
            }
            const daf::Size needPivot = s - keepC.size(); // 还需从 pivot 中选的顶点数
            if (needPivot == 0) {
                if (!isPivotVertex) {
                    sup[minKeepCore]++;
                }
                continue;
            }

            std::ranges::sort(povitC, [core](daf::Size a, daf::Size b) {
                return core[a] > core[b];
            });

            daf::Size numBiggerInPovit = 0;
            for (const auto v: povitC) {
                if (core[v] >= core[currV] && core[v] >= minKeepCore) {
                    numBiggerInPovit++;
                } else {
                    break;
                }
            }

            double prveCore = minKeepCore;
            double prveCliqueCount = 0;
            // povitC.print();
            // keepC.print();
            for (daf::Size i = numBiggerInPovit; i <= povitC.size(); i++) {
                if (i == povitC.size() || core[povitC[i]] < prveCore) {
                    // const double povitCliqueCount = nCr[i - 1 - r + 1][needPivot - 1 - r + 1];
                    if (!isPivotVertex) {
                        const double keepCliqueCount = nCr[i - r + 1][needPivot - r + 1];
                        sup[prveCore] += keepCliqueCount - prveCliqueCount;
                        prveCliqueCount = keepCliqueCount;
                    } else {
                        const double povitCliqueCount = nCr[i - 1 - r + 1][needPivot - 1 - r + 1];
                        sup[prveCore] += povitCliqueCount - prveCliqueCount;
                        prveCliqueCount = povitCliqueCount;
                    }
                    if (i == povitC.size()) break;
                    prveCore = core[povitC[i]];
                }
            }
            // std::cout << "vertex: " << currV
            // << " leaf: " << leaf->leafId
            // << " sup: " << sup << std::endl;
        }
        povitC.free();
        keepC.free();
        double prveCount = 0;
        for (auto &[c, count]: sup) {
            if (prveCount + count >= c) {
                double newCore = std::max(c, prveCount);
                if (core[currV] > newCore) {
                    return newCore;
                }
            }
            prveCount += count;
        }
        if (prveCount < core[currV]) {
            return prveCount;
        }
        return core[currV];
    }
}


void NucleusCoreDecomposition(MultiBranchTree &tree,
                              daf::StaticVector<TreeNode *> &leafList,
                              Graph &treeGraphV,
                              Graph &leafGraph,
                              daf::CliqueSize r, daf::CliqueSize s) {
    auto time_start = std::chrono::high_resolution_clock::now();
    auto *core = new double[treeGraphV.getGraphNodeSize()];
    for (daf::Size i = 0; i < treeGraphV.getGraphNodeSize(); i++) {
        const auto [nbr_start, nbr_end] = treeGraphV.getNbr(i);
        if (nbr_end == nbr_start) {
            core[i] = 0;
        } else if (nbr_end - nbr_start == 1) {
            const auto cliqueSize = leafList[treeGraphV.adj_list[nbr_start]]->MaxDeep;
            core[i] = nCr[cliqueSize - r][s - r];
        } else {
            core[i] = std::numeric_limits<double>::max();
        }
    }

    // daf::printArray(core, treeGraphV.getGraphNodeSize());
    bool update = true;
    daf::Size iter = 0;
    while (update) {
        update = false;
        std::cout << ++iter << " " << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
        time_start = std::chrono::high_resolution_clock::now();
        for (daf::Size v = 0; v < treeGraphV.getGraphNodeSize(); v++) {
            if (treeGraphV.getNbrCount(v) <= 1) {
                continue;
            }
            auto [leafStart, leafEnd] = treeGraphV.getNbr(v);
            double maxCoreForV = std::numeric_limits<double>::max();
            for (daf::Size i = leafStart; i < leafEnd; i++) {
                const auto leafId = treeGraphV.adj_list[i];
                const auto [leafNbrStart, leafNbrEnd] = leafGraph.getNbr(leafId);
                daf::StaticVector<daf::Size> intersect(std::min(leafNbrEnd - leafNbrStart, leafEnd - leafStart));
                std::ranges::set_intersection(
                    leafGraph.getNbrVec(leafId),
                    treeGraphV.getNbrVec(v),
                    std::back_inserter(intersect)
                );

                auto newCore = NucleusCD::coreEvaluation(intersect, leafList, r, s, v, core);

                // intersect.print();

                maxCoreForV = std::min(maxCoreForV, newCore);
                intersect.free();
                // 需要计算这个 leafId 的 nbr和节点的leaf nbr
            }
            if (core[v] > maxCoreForV) {
                core[v] = maxCoreForV;
                update = true;
            }
            //
        }
        // std::cout << "core: ";
        daf::printArray(core, treeGraphV.getGraphNodeSize());
    }

    auto file = fopen("/Users/zhangwenqian/UNSW/pivoter/a", "w");
    std::sort(core, core + tree.getRoot()->children.size());
    for (daf::Size i = 0; i < tree.getRoot()->children.size(); i++) {
        fprintf(file, "%f\n", core[i]);
    }
    fclose(file);
    delete[] core;
}
