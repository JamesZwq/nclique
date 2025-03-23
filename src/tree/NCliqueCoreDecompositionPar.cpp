//
// Created by 张文谦 on 25-3-4.
//

#include "../tree/NCliqueCoreDecomposition.h"
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <mutex>

#include "dataStruct/StaticVector.h"

extern double nCr[1001][401];
namespace basePar{
void countingPerVertexHelp(const TreeNode &node,
                           const daf::CliqueSize k,
                           MutexStaticVector<double>& core,
                           StaticVector<daf::Size> &povit,
                           StaticVector<daf::Size> &keepC
) {
    daf::Size cliqueSize = povit.size() + keepC.size();
    if (node.children.empty() && cliqueSize >= k && keepC.size() <= k) {
        const int needPivot = k - keepC.size(); // 还需从 pivot 中选的顶点数
        double totalKcliques = 0;
        if (needPivot >= 0 && needPivot <= povit.size()) {
            totalKcliques = nCr[povit.size()][needPivot];
        }
        for (const auto v: keepC) {
            // if (v == 0) {
                // std::cout << "keep: " << v << " totalKcliques: " << totalKcliques << std::endl;
            // }
            // core[v] += totalKcliques;
            core.add(v, totalKcliques);
        }

        double eachPivotKcliques = 0;
        const int needPivotWithV = needPivot - 1;
        if (needPivotWithV >= 0 && needPivotWithV <= povit.size() - 1) {
            eachPivotKcliques = nCr[povit.size() - 1][needPivotWithV];
        }

        for (auto v: povit) {
            // std::cout << "povit: " << v << " eachPivotKcliques: " << eachPivotKcliques << std::endl;
            // core[v] += eachPivotKcliques;
            core.add(v, eachPivotKcliques);
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
    // auto *core = new daf::Size[tree.getRoot()->children.size()];
    auto core = new MutexStaticVector<double>(tree.getRoot()->children.size());
    //init 0
    // std::memset(core, 0, tree.getRoot()->children.size() * sizeof(daf::Size));
    // for (auto node: tree.getRoot()->children) {
    //     if (node->MaxDeep < k) {
    //         continue;
    //     }
    //     keepC.push_back(node->v);
    //     countingPerVertexHelp(*node, k, core, povitC, keepC);
    //     keepC.pop_back();
    // }

    tbb::parallel_for(tbb::blocked_range<size_t>(0, tree.getRoot()->children.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto node = tree.getRoot()->children[i];
                if (node->MaxDeep < k) {
                    continue;
                }
                StaticVector<daf::Size> povitC;
                StaticVector<daf::Size> keepC;
                keepC.push_back(node->v);
                countingPerVertexHelp(*node, k, *core, povitC, keepC);
            }
        });
    return core->data;
}



void computeSupHelp(const TreeNode &node,
                    const daf::CliqueSize k,
                    const double *core,
                    StaticVector<daf::Size> &povit,
                    StaticVector<daf::Size> &keepC,
                    std::vector<ThreadSafeMap<double, double>> &sup
) {
    daf::Size cliqueSize = povit.size() + keepC.size();
    if (node.children.empty() && cliqueSize >= k && keepC.size() <= k) {
        assert(!keepC.empty());

        StaticVector<daf::Size> orderPovit = povit;
        // sort from large to small
        std::ranges::sort(orderPovit, [core](daf::Size a, daf::Size b) {
            return core[a] > core[b];
        });
        const daf::Size needPivot = k - keepC.size(); // 还需从 pivot 中选的顶点数


        double maxK = core[*std::ranges::min_element(keepC, [core](daf::Size a, daf::Size b) {
                                                return core[a] < core[b];
                                            })];
        if (needPivot == 0) {
            for (const auto v: keepC) {
                // sup[v][maxK]++;
                sup[v].add(maxK, 1);
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
                const double keepCliqueCount = nCr[i][needPivot];
                const double povitCliqueCount = nCr[i - 1][needPivot - 1];
                for (const auto v: keepC) {
                    // sup[v][prveCore] += keepCliqueCount - prveKeepCliqueCount;
                    sup[v].add(prveCore, keepCliqueCount - prveKeepCliqueCount);
                }
                for (daf::Size v = 0; v < prvePovit; v++) {
                    // sup[orderPovit[v]][prveCore] += povitCliqueCount - prvepovitCliqueCount;
                    sup[orderPovit[v]].add(prveCore, povitCliqueCount - prvepovitCliqueCount);
                }
                for (daf::Size v = prvePovit; v < i; v++) {
                    // sup[orderPovit[v]][prveCore] += povitCliqueCount;
                    sup[orderPovit[v]].add(prveCore, povitCliqueCount);
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

std::vector<ThreadSafeMap<double,double>> computeSup(const MultiBranchTree &tree,
                const daf::CliqueSize k,
                const double *core) {
    auto &children = tree.getRoot()->children;
    std::vector<ThreadSafeMap<double, double>> sup(children.size());

    // 并行遍历每个子节点
    tbb::parallel_for(tbb::blocked_range<size_t>(0, children.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                auto node = children[i];
                if (node->MaxDeep < k) {
                    continue;
                }
                // 为每个迭代创建独立的局部变量，防止多线程竞争
                StaticVector<daf::Size> povitC;
                StaticVector<daf::Size> keepC;
                keepC.push_back(node->v);
                computeSupHelp(*node, k, core, povitC, keepC, sup);
            }
        });
    return sup;
}

}
void baseNucleusCoreDecompositionPar(const MultiBranchTree &tree, daf::CliqueSize k){
    daf::Size numThreads = 1;
    tbb::global_control gc(tbb::global_control::max_allowed_parallelism, numThreads);
    std::cout << "numThreads: " << tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism) << std::endl;

    auto time_start = std::chrono::high_resolution_clock::now();
    auto *core = basePar::countingPerVertex(tree, k);
    // std::cout << "init core: ";
    // daf::printArray(core, tree.getRoot()->children.size());
    bool update = true;
    daf::Size iter = 0;
    while (update) {
        update = false;
        std::cout << ++iter << " " << std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
        time_start = std::chrono::high_resolution_clock::now();
        auto sup = basePar::computeSup(tree, k, core);
        // std::cout << sup << std::endl;
        for (daf::Size v = 0; v < tree.getRoot()->children.size(); v++) {
            auto &supV = sup[v];
            double prveCount = 0;
            for (auto &[c, count]: supV.map) {
                if (prveCount + count >= c) {
                    auto newCore = std::max(c,prveCount);
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
    }

    // auto file = fopen("/Users/zhangwenqian/UNSW/pivoter/b", "w");
    // std::sort(core, core + tree.getRoot()->children.size());
    // for (daf::Size i = 0; i < tree.getRoot()->children.size(); i++) {
    //     fprintf(file, "%d\n", core[i]);
    // }
    // fclose(file);
    delete[] core;
}
