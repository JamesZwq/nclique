//
// Created by 张文谦 on 25-3-4.
//

#include "../tree/NCliqueCoreDecomposition.h"
#include <boost/heap/d_ary_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraph.h"

extern double nCr[1001][401];
// 放在你的函数外（比如文件顶部），保证编译时可见并内联
template<typename It1, typename It2, typename UpdateFunc>
inline void processEdgePairs(It1 b1, It1 e1, bool isPovit1,
                             It2 b2, It2 e2, bool isPovit2,
                             double weight,
                             UpdateFunc &&upd) noexcept {
    if (weight < 0.0) return;

    // 判断两个区间迭代器是否相同
    // if all same, do nothing
    if (b1 == b2 && e1 == e2 && b1 == e1 && b2 == e2) {
        return;
    }
    if (b1 == b2 && e1 == e2) {
        // 同一范围：i < j
        for (auto it = b1; it + 1 != e1; ++it) {
            auto u = *it;
            for (auto jt = it + 1; jt != e1; ++jt) {
                upd(u, isPovit1, *jt, isPovit2, weight);
            }
        }
    } else {
        // 不同范围：笛卡尔积
        for (auto it = b1; it != e1; ++it) {
            auto u = *it;
            for (auto jt = b2; jt != e2; ++jt) {
                upd(u, isPovit1, *jt, isPovit2, weight);
            }
        }
    }
}


namespace baseECD {
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

    void countingPerEdgeHelp(const TreeNode &node,
                             const daf::CliqueSize k,
                             const Graph &edgeGraph,
                             double *coreE,
                             daf::Size *degreeE,
                             // EdgeHashMap<double> coreE,
                             daf::StaticVector<daf::Size> &povit,
                             daf::StaticVector<daf::Size> &keepC
    ) {
        daf::Size cliqueSize = povit.size() + keepC.size();
        if (node.children.empty() && cliqueSize >= k && keepC.size() <= k) {
            // 还需从 pivot 里选的点数
            int needPivot = int(k) - int(keepC.size());
            // 1) keep-keep 边：两端都在 keepC 中
            double totalKcliques = -1;
            // 从 pivot.size() 个点里选 needPivot 个，再分配给每条 keep‑keep 边
            if (needPivot >= 1 && needPivot <= int(povit.size())) {
                totalKcliques = nCr[povit.size()][needPivot];
                for (size_t i = 0; i < keepC.size(); ++i) {
                    for (size_t j = i + 1; j < keepC.size(); ++j) {
                        daf::Size u = keepC[i], v = keepC[j];
                        auto index = edgeGraph.getEdgeIndex(u, v);
                        coreE[index] += totalKcliques;
                        degreeE[index]++;
                    }
                }
            }

            // 2) pivot‑pivot 边：两端都在 povit 中
            double eachPivotKcliques = -1;
            int needPivotWithV = needPivot - 2;
            if (needPivotWithV >= 0 && needPivotWithV <= static_cast<int>(povit.size()) - 2) {
                eachPivotKcliques = nCr[povit.size() - 2][needPivotWithV];
                for (size_t i = 0; i < povit.size(); ++i) {
                    for (size_t j = i + 1; j < povit.size(); ++j) {
                        daf::Size u = povit[i], v = povit[j];
                        auto index = edgeGraph.getEdgeIndex(u, v);
                        coreE[index] += eachPivotKcliques;
                        degreeE[index]++;
                    }
                }
            }


            // 3) cross 边：一端在 keepC，一端在 povit
            double eachKeepPivotKcliques = -1;
            int needKeepPivotWithV = needPivot - 1;
            if (needKeepPivotWithV >= 0 && needKeepPivotWithV <= static_cast<int>(povit.size()) - 1) {
                eachKeepPivotKcliques = nCr[povit.size() - 1][needKeepPivotWithV];
                for (size_t i = 0; i < keepC.size(); ++i) {
                    for (size_t j = 0; j < povit.size(); ++j) {
                        daf::Size u = keepC[i], v = povit[j];
                        auto index = edgeGraph.getEdgeIndex(u, v);
                        coreE[index] += eachKeepPivotKcliques;
                        degreeE[index]++;
                    }
                }
            }
            // std::cout << "keep: " << keepC << " povit: " << povit << " k: " << k
            // << " needPivot: " << needPivot
            //           << " totalKcliques: " << totalKcliques
            //           << " eachPivotKcliques: " << eachPivotKcliques
            //           << " eachKeepPivotKcliques: " << eachKeepPivotKcliques
            //           << std::endl;

            return;
        }
        for (const auto &child: node.children) {
            if (child->MaxDeep < k) {
                continue;
            }
            if (child->isPivot) {
                povit.push_back(child->v);
                countingPerEdgeHelp(*child, k, edgeGraph, coreE, degreeE, povit, keepC);
                povit.pop_back();
            } else {
                keepC.push_back(child->v);
                countingPerEdgeHelp(*child, k, edgeGraph, coreE, degreeE, povit, keepC);
                keepC.pop_back();
            }
        }
    }

    std::pair<double *, daf::Size *> countingPerEdge(const MultiBranchTree &tree, const Graph &edgeGraph,
                                                     const daf::CliqueSize k) {
        // EdgeHashMap<double> coreE(edgeGraph.adj_list.size());
        double *coreE = new double[edgeGraph.adj_list.size()];
        daf::Size *degreeE = new daf::Size[edgeGraph.adj_list.size()];
        memset(coreE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));
        memset(degreeE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));
        daf::StaticVector<daf::Size> povitC;
        daf::StaticVector<daf::Size> keepC;
        daf::Size count = 0;
        for (auto node: tree.getRoot()->children) {
            if (node->MaxDeep < k) {
                continue;
            }
            keepC.push_back(node->v);
            countingPerEdgeHelp(*node, k, edgeGraph, coreE, degreeE, povitC, keepC);
            keepC.pop_back();
        }
        return {coreE, degreeE};
    }


    std::pair<double *, daf::Size *> countingPerEdge(const DynamicGraph<TreeGraphNode> &treeGraph, const Graph &edgeGraph,
                                                     const daf::CliqueSize k) {
        // EdgeHashMap<double> coreE(edgeGraph.adj_list.size());
        double *countingE = new double[edgeGraph.adj_list.size()];
        daf::Size *degreeE = new daf::Size[edgeGraph.adj_list.size()];
        memset(countingE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));
        memset(degreeE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));
        daf::StaticVector<daf::Size> povit;
        daf::StaticVector<daf::Size> keepC;
        daf::Size count = 0;
        for (const auto& clique: treeGraph.adj_list) {
            povit.clear();
            keepC.clear();
            if (clique.size() < k) {
                continue;
            }
            for (auto & i : clique) {
                if (i.isPivot) {
                    povit.push_back(i.v);
                } else {
                    keepC.push_back(i.v);
                }
            }

            int needPivot = int(k) - int(keepC.size());
            // 1) keep-keep 边：两端都在 keepC 中
            double totalKcliques = -1;
            // 从 pivot.size() 个点里选 needPivot 个，再分配给每条 keep‑keep 边
            if (needPivot >= 1 && needPivot <= int(povit.size())) {
                totalKcliques = nCr[povit.size()][needPivot];
                for (size_t i = 0; i < keepC.size(); ++i) {
                    for (size_t j = i + 1; j < keepC.size(); ++j) {
                        daf::Size u = keepC[i], v = keepC[j];
                        auto index = edgeGraph.getEdgeIndex(u, v);
                        countingE[index] += totalKcliques;
                        degreeE[index]++;
                        if (std::min(u,v) == 0 && std::max(u,v) == 1) {
                            std::cout << "1111 leaf: " << clique << " degreeE: " << degreeE[index] << std::endl;
                        }
                    }
                }
            }

            // 2) pivot‑pivot 边：两端都在 povit 中
            double eachPivotKcliques = -1;
            int needPivotWithV = needPivot - 2;
            if (needPivotWithV >= 0 && needPivotWithV <= static_cast<int>(povit.size()) - 2) {
                eachPivotKcliques = nCr[povit.size() - 2][needPivotWithV];
                for (size_t i = 0; i < povit.size(); ++i) {
                    for (size_t j = i + 1; j < povit.size(); ++j) {
                        daf::Size u = povit[i], v = povit[j];
                        auto index = edgeGraph.getEdgeIndex(u, v);
                        countingE[index] += eachPivotKcliques;
                        degreeE[index]++;
                        if (std::min(u,v) == 0 && std::max(u,v) == 1) {
                            std::cout << "1111 leaf: " << clique << " degreeE: " << degreeE[index] << std::endl;
                        }
                    }
                }
            }


            // 3) cross 边：一端在 keepC，一端在 povit
            double eachKeepPivotKcliques = -1;
            int needKeepPivotWithV = needPivot - 1;
            if (needKeepPivotWithV >= 0 && needKeepPivotWithV <= static_cast<int>(povit.size()) - 1) {
                eachKeepPivotKcliques = nCr[povit.size() - 1][needKeepPivotWithV];
                for (size_t i = 0; i < keepC.size(); ++i) {
                    for (size_t j = 0; j < povit.size(); ++j) {
                        daf::Size u = keepC[i], v = povit[j];
                        auto index = edgeGraph.getEdgeIndex(u, v);
                        countingE[index] += eachKeepPivotKcliques;
                        degreeE[index]++;
                        if (std::min(u,v) == 0 && std::max(u,v) == 1) {
                            std::cout << "1111 leaf: " << clique << " degreeE: " << degreeE[index] << std::endl;
                        }
                    }
                }
            }

        }
        povit.free();
        keepC.free();
        return {countingE, degreeE};
    }

    double *initEdgeCore(const Graph &edgeGraph, const double *coreV) {
        const daf::Size m = edgeGraph.adj_list.size();
        const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;
        double *coreE = new double[m];

        // Number of vertices in edgeGraph

        // For each vertex u, go over its outgoing edge‑entries [start..end)
        // and set coreE[idx] = min(coreV[u], coreV[v])
        for (daf::Size u = 0; u < n; ++u) {
            const daf::Size start = edgeGraph.adj_list_offsets[u];
            const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
            for (daf::Size idx = start; idx < end; ++idx) {
                const daf::Size v = edgeGraph.adj_list[idx];
                coreE[idx] = std::min(coreV[u], coreV[v]);
            }
        }
        return coreE;
    }


    double *initVertexCore(const Graph &edgeGraph, const double *coreE) {
        const daf::Size m = edgeGraph.adj_list.size();
        const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;
        double *coreV = new double[n];
        memset(coreV, 0, n * sizeof(daf::Size));
        // Number of vertices in edgeGraph

        // For each vertex u, go over its outgoing edge‑entries [start..end)
        // and set coreE[idx] = min(coreV[u], coreV[v])
        for (daf::Size u = 0; u < n; ++u) {
            const daf::Size start = edgeGraph.adj_list_offsets[u];
            const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
            for (daf::Size idx = start; idx < end; ++idx) {
                // coreV[idx] = std::max(coreV[u], coreE[idx]);
                coreV[edgeGraph.adj_list[idx]] = std::max(coreV[edgeGraph.adj_list[idx]], coreE[idx]);
                coreV[u] = std::max(coreV[u], coreE[idx]);
            }
        }
        return coreV;
    }

    template<typename T>
    void printEdgeCore(const Graph &edgeGraph, const T *coreE) {
        const daf::Size m = edgeGraph.adj_list.size();
        const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;
        for (daf::Size u = 0; u < n; ++u) {
            const daf::Size start = edgeGraph.adj_list_offsets[u];
            const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
            for (daf::Size idx = start; idx < end; ++idx) {
                std::cout << "[" << u << ", " << edgeGraph.adj_list[idx] << "] " << coreE[idx] << " ";
            }
            std::cout << std::endl;
        }
    }

    double *initLeafCore(const DynamicGraph<TreeGraphNode> &tree, double * &coreE, daf::Size k,
                         const Graph &edgeGraph) {
        // init as the min one in the edge
        // daf::em
        // daf::enumerateCombinations()
        auto *coreLeaf = new double[tree.adj_list.size()];
        // memset(coreLeaf, std::numeric_limits<double>::max(), sizeof(double) * leafList.size());
        const daf::Size numLeaf = tree.adj_list.size();
        for (daf::Size i = 0; i < numLeaf; ++i) {
            auto leaf = tree.adj_list[i];
            // TODO: add lowerBound
            // double lowerBound = nCr[povit.size() - 2][k - 2];



            // std::cout << "leaf: " << leaf->leafId << " keepC: " << keepC << " povit: " << povit
            //           << " k: " << k << std::endl;
            double minCore = std::numeric_limits<double>::max();
            // daf::enumerateCombinations(keepC, povit, 2, k,
            //                            [&](const daf::StaticVector<daf::Size> &keepC,
            //                                const daf::StaticVector<daf::Size> &povit) {
            //                                if (keepC.size() + povit.size() != 2) {
            //                                    std::cerr << "Error: [" << k << ", " << povit << "], ["
            //                                            << keepC << "], the size is not 1" << std::endl;
            //                                }
            //                                // auto edgeCore = coreE[edgeGraph.getEdgeIndex(keepC[0], povit[0])];
            //                                double edgeCore;
            //                                edgeCore = coreE[edgeGraph.getEdgeIndex(povit[0], povit[1])];
            //                                // std::cout << "keepC: " << keepC << " povit: " << povit << " edgeCore: "
            //                                //           << edgeCore << " lowerBound: " << lowerBound << std::endl;
            //                                if (edgeCore == lowerBound) {
            //                                    minCore = edgeCore;
            //                                    return false;
            //                                }
            //                                minCore = std::min(minCore, edgeCore);
            //                                return true;
            //                            });

            for (daf::Size j = 0; j < leaf.size(); ++j) {
                for (daf::Size k = j + 1; k < leaf.size(); ++k) {
                    auto u = leaf[j];
                    auto v = leaf[k];
                    auto index = edgeGraph.getEdgeIndex(u.v, v.v);
                    double edgeCore = coreE[index];
                    // if (edgeCore == lowerBound) {
                    //     minCore = edgeCore;
                    //     break;
                    // }
                    minCore = std::min(minCore, edgeCore);
                }
            }
            coreLeaf[i] = minCore;
        }

        // std::cout << "numLeaf: " << numLeaf << std::endl;

        return coreLeaf;
    }


    struct CompareLeaf {
        const double *coreLeaf; // 指向外部数组
        CompareLeaf(const double *cl) : coreLeaf(cl) {
        }

        // 注意：这里要返回 “a 排在前面” 的条件，为最小堆写成 coreLeaf[a] > coreLeaf[b]
        bool operator()(daf::Size const &a, daf::Size const &b) const {
            return coreLeaf[a] > coreLeaf[b];
        }
    };
}

void baseNucleusEdgeCoreDecomposition(DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
                                      DynamicGraph<daf::Size> &treeGraphV, daf::CliqueSize k) {

    tree.printGraphPerV();
    // daf::Size numNodes = tree.getRoot()->children.size();
    // tree.printTree();
    auto time_start = std::chrono::high_resolution_clock::now();
    auto [countingKE, degreeE] = baseECD::countingPerEdge(tree, edgeGraph, k);

    double *leafCore = baseECD::initLeafCore(tree, countingKE, k, edgeGraph);

#ifndef NDEBUG
    daf::printArray(countingKE, edgeGraph.adj_list.size());
    daf::printArray(degreeE, tree.adj_list.size());
    baseECD::printEdgeCore(edgeGraph, countingKE);
    baseECD::printEdgeCore(edgeGraph, degreeE);
    daf::printArray(leafCore, tree.adj_list.size());
#endif

    using DHeap = boost::heap::d_ary_heap<
        daf::Size,
        boost::heap::arity<8>,
        boost::heap::mutable_<true>,
        boost::heap::compare<baseECD::CompareLeaf>
    >;

    DHeap heap{baseECD::CompareLeaf(leafCore)};
    daf::StaticVector<DHeap::handle_type> heapHandles(tree.adj_list.size());
    heapHandles.c_size = tree.adj_list.size();

    for (daf::Size i = 0; i < tree.adj_list.size(); ++i) {
        auto leaf = tree.adj_list[i];
        if (leaf.size() < k) {
            std::cerr << "Error: leaf id is not equal to index" << std::endl;
        }
        heapHandles[i] = heap.push(i);
    }
    daf::StaticVector<bool> removedLeaf(tree.adj_list.size());
    removedLeaf.c_size = tree.adj_list.size();
    memset(removedLeaf.getData(), false, tree.adj_list.size() * sizeof(bool));

    auto *coreE = new double[edgeGraph.adj_list.size()];
    memset(coreE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));

    auto *degreeV = new daf::Size[treeGraphV.adj_list.size()];
    // std::cout << "adj_list_offsets: " << treeGraphV.adj_list_offsets << std::endl;
    for (daf::Size i = 0; i < treeGraphV.adj_list.size() - 1; ++i) {
        degreeV[i] = treeGraphV.adj_list[i].size();
    }

    // std::cout << "degreeV: " << std::endl;
    // daf::printArray(degreeV, treeGraphV.adj_list_offsets.size() - 1);
    // std::cout << std::endl;

    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keepC;
    daf::StaticVector<daf::Size> removedPovit;
    daf::StaticVector<bool> isRemoved(treeGraphV.adj_list.size());
    isRemoved.c_size = treeGraphV.adj_list.size();
    memset(isRemoved.getData(), false, treeGraphV.adj_list.size() * sizeof(bool));

    daf::StaticVector<daf::Size> currentRemoveLeafIds(tree.adj_list.size());
    bool removedKeepC = false;
    daf::StaticVector<std::pair<daf::Size, double> > updateLeaf(tree.adj_list.size() * 10);

    double currCore = 0;
    std::cout << "=========================begin=========================" << std::endl;
    while (!heap.empty()) {
        double minCore = leafCore[heap.top()];
        // 一次循环把所有 core==minCore 的 leaf 全部 pop 出来
        while (!heap.empty() && leafCore[heap.top()] == minCore) {
            auto id = heap.top();
            removedLeaf[id] = true;
            heap.pop();
            currentRemoveLeafIds.push_back(id);
        }

        for (auto currLeafId: currentRemoveLeafIds) {
            auto leaf = tree.adj_list[currLeafId];
            daf::StaticVector<daf::Size> povit;
            daf::StaticVector<daf::Size> keepC;
            for (auto i = 0; i < leaf.size(); ++i) {
                --degreeV[leaf[i].v];
                if (leaf[i].isPivot) {
                    povit.push_back(leaf[i].v);
                } else {
                    keepC.push_back(leaf[i].v);
                }
            }
            // for (auto j = i + 1; j < leaf.size(); ++j) {
            //     auto u = leaf[i].v, v = leaf[j].v;
            //     auto idx = edgeGraph.getEdgeIndex(u, v);
            //     --degreeE[idx];
            // }
        }
        std::cout << "currentRemoveLeafIds: " << currentRemoveLeafIds << std::endl;
        for (auto leafId: currentRemoveLeafIds) {
            auto leaf = tree.adj_list[leafId];
            currCore = std::max(currCore, leafCore[leafId]);
            for (auto node : leaf) {
                if (node.isPivot) {
                    povit.push_back(node.v);
                    if (degreeV[node.v] == 0) {
                        removedPovit.push_back(node.v);
                        isRemoved[node.v] = true;
                    }
                } else {
                    keepC.push_back(node.v);
                    if (degreeV[node.v] == 0) {
                        removedKeepC = true;
                    }
                }
            }


            removedPovit.print("removedPovit1: ");
            // double totalKcliques = -1;totalKcliques = nCr[povit.size()][needPivot];

            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            for (auto i = 0; i < keepC.size() && !removedKeepC; ++i) {
                for (auto j = i + 1; j < keepC.size() && !removedKeepC; ++j) {
                    auto u = keepC[i], v = keepC[j];
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    if (degreeE[idx] == 0) {
                        removedKeepC = true;
                    }
                }
            }

            for (daf::Size i = 0; i < povit.size(); ++i) {
                auto u = povit[i];
                if (isRemoved[u]) continue;
                for (auto j = 0; j < keepC.size(); ++j) {
                    auto v = keepC[j];
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    if (degreeE[idx] == 0 && degreeV[u] != 0 && degreeV[v] != 0) {
                        removedPovit.push_back(u);
                        isRemoved[u] = true;
                        break;
                    }
                }
            }

            for (daf::Size i = 0; i < povit.size(); ++i) {
                for (daf::Size j = i + 1; j < povit.size(); ++j) {
                    auto u = povit[i], v = povit[j];
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    if (degreeE[idx] == 0 && degreeV[u] != 0 && degreeV[v] != 0 && !isRemoved[u] && !isRemoved[v]) {
                        // TODO
                        std::cout << "###Error: degreeE is equal to 0 : " << degreeE[idx] << " u: " << u
                                << " v: " << v << std::endl;
                    }
                }
            }
            removedPovit.print("removedPovit2: ");
            isRemoved.print("isRemoved: ");
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            double KtoK = 0;
            double KtoP = 0;
            double PtoP = 0;
            daf::Size needPivot = k - keepC.size();
            if (needPivot >= 1 && needPivot <= povit.size()) {
                if (removedKeepC) {
                    KtoK = nCr[povit.size()][needPivot];
                } else {
                    KtoK = nCr[povit.size()][needPivot] - nCr[povit.size() - removedPovit.size()][needPivot];
                }
            }
            int needPivotWithV = needPivot - 2;
            if (needPivotWithV >= 0 && needPivotWithV <= povit.size() - 2) {
                // PtoP = nCr[povit.size() - 2][needPivotWithV];
                if (removedKeepC) {
                    PtoP = nCr[povit.size() - 2][needPivotWithV];
                } else {
                    PtoP = nCr[povit.size() - 2][needPivotWithV] - nCr[povit.size() - 2 - removedPovit.size()][
                               needPivotWithV];
                }
            }
            int needKeepPivotWithV = needPivot - 1;
            if (needKeepPivotWithV >= 0 && needKeepPivotWithV <= povit.size() - 1) {
                if (removedKeepC) {
                    KtoP = nCr[povit.size() - 1][needKeepPivotWithV];
                } else {
                    KtoP = nCr[povit.size() - 1][needKeepPivotWithV] - nCr[povit.size() - 1 - removedPovit.size()][
                               needKeepPivotWithV];
                }
            }


            auto minCounting = std::numeric_limits<double>::max();
            auto updateEdge = [&](daf::Size u, bool uP, daf::Size v, bool vP, double w) {
                auto idx = edgeGraph.getEdgeIndex(u, v);
                countingKE[idx] -= w;
                // if (u == 0 && v == 1) {
                //     std::cout << "Error: u or v is 0 degree: " << degreeE[idx] << " " << countingKE[idx] << std::endl;
                // }
                if (degreeE[idx] == 0 || countingKE[idx] <= currCore) {
                    coreE[idx] = currCore;
                } else {
                    minCounting = std::min(minCounting, countingKE[idx]);
                }
                auto uIter = treeGraphV.getNbr(u);
                auto vIter = treeGraphV.getNbr(v);
                daf::intersect_with_callback(uIter->begin(), uIter->end(),
                                             vIter->begin(), vIter->end(),
                                             [&](const daf::Size &x) {
                                                 if (leafCore[x] > countingKE[idx] && !removedLeaf[x]) {
                                                     updateLeaf.push_back(std::make_pair(x, countingKE[idx]));
                                                 }
                                             });
            };

            // 1) keep–keep
            processEdgePairs(
                keepC.begin(), keepC.end(), false,
                keepC.begin(), keepC.end(), false,
                KtoK,
                updateEdge
            );

            // 2) pivot–pivot
            processEdgePairs(
                povit.begin(), povit.end(), true,
                povit.begin(), povit.end(), true,
                PtoP,
                updateEdge
            );

            // 3) keep–pivot
            processEdgePairs(
                keepC.begin(), keepC.end(), false,
                povit.begin(), povit.end(), true,
                KtoP,
                updateEdge
            );

            std::ranges::sort(updateLeaf,
                              [](auto &a, auto &b) {
                                  if (a.first != b.first) return a.first < b.first;
                                  return a.second < b.second;
                              }
            );

            updateLeaf.unique(
                [](auto &a, auto &b) {
                    return a.first == b.first;
                }
            );

            for (auto &u: updateLeaf) {
                daf::Size leafId = u.first;
                double newCore = u.second;
                if (leafCore[leafId] > newCore && !removedLeaf[leafId]) {
                    leafCore[leafId] = newCore;
                    heap.update(heapHandles[leafId]);
                }
            }

            if (removedKeepC) {
                leafCore[leafId] = 0;
                tree.removeNode(leafId);
            } else if (!removedPovit.empty() && needPivot <= povit.size() - removedPovit.size()) {
                const auto newLeaf = tree.removeNbrs(leafId, removedPovit);
                if (!newLeaf.empty()) {
                    leafCore[leafId] = minCounting;
                    heapHandles[leafId] = heap.push(leafId);
                    removedLeaf[leafId] = false;
                    for (auto i = 0; i < newLeaf.size(); ++i) {
                        degreeV[newLeaf[i].v]++;
                        for (auto j = i + 1; j < newLeaf.size(); ++j) {
                            auto u = newLeaf[i].v, v = newLeaf[j].v;
                            auto idx = edgeGraph.getEdgeIndex(u, v);
                            degreeE[idx]++;
                        }
                    }
                } else {
                    leafCore[leafId] = 0;
                }
            } else {
                leafCore[leafId] = 0;
                tree.removeNode(leafId);
            }


#ifndef NDEBUG
            std::cout << "\n\n leafId: " << leafId << " leafCore: " << leafCore[leafId] << " currCore: " << currCore
                    << " KtoK: " << KtoK << " KtoP: " << KtoP << " PtoP: " << PtoP
                    << std::endl;
            updateLeaf.print("updateLeaf");
            removedLeaf.print("removedLeaf");
            removedPovit.print("removedPovit");
            std::cout << removedKeepC << std::endl;
            std::cout << "countingKE: ";
            baseECD::printEdgeCore(edgeGraph, countingKE);

            std::cout << "degreeE: ";
            baseECD::printEdgeCore(edgeGraph, degreeE);
            std::cout << "degreeV: ";
            daf::printArray(degreeV, treeGraphV.adj_list.size() - 1);

            std::cout << "coreE: ";
            baseECD::printEdgeCore(edgeGraph, coreE);

            std::cout << "leafCore: ";
            daf::printArray(leafCore, tree.adj_list.size());

            removedLeaf.print("removedLeaf");

            tree.printGraphPerV();

            std::cout << "============================================================" << std::endl;
#endif


            updateLeaf.clear();
            povit.clear();
            keepC.clear();
            removedPovit.clear();
            // leafIds.clear();
            removedKeepC = false;
            for (auto v : removedPovit) {
                isRemoved[v] = false;
            }
        }
        currentRemoveLeafIds.clear();
        // sleep
        // std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    // coreE
    // daf::printArray(coreE, edgeGraph.adj_list.size());

#ifndef NDEBUG
    baseECD::printEdgeCore(edgeGraph, coreE);
#endif

    // /Users/zhangwenqian/UNSW/pivoter/a
    auto file = fopen("/Users/zhangwenqian/UNSW/pivoter/a", "w");
    std::sort(coreE, coreE + edgeGraph.adj_list.size());
    for (daf::Size i = 0; i < edgeGraph.adj_list.size(); i++) {
        // cover to int
        fprintf(file, "%d\n", (int) coreE[i]);
    }
    fclose(file);


    delete[] countingKE;
    delete[] leafCore;
    delete[] degreeV;
    delete[] coreE;
    updateLeaf.free();
    povit.free();
    keepC.free();
    removedPovit.free();
    currentRemoveLeafIds.free();
    // leafIds.free();
}
