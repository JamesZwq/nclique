//
// Created by 张文谦 on 25-3-4.
//

#include "../tree/NCliqueCoreDecomposition.h"
#include <boost/heap/d_ary_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <set>

#include "BronKerbosch.h"
#include "dataStruct/CliqueHashMap.h"
#include "debug/EdgeSet.h"
#include "graph/DynamicGraph.h"

extern double nCr[1001][401];
// 放在你的函数外（比如文件顶部），保证编译时可见并内联
// #ifndef NDEBUG
// set NOEBUG as trus
namespace PlusECD {
    template<typename It1, typename It2, typename UpdateFunc>
    inline void processEdgePairsImpl(It1 b1, It1 e1,
                                     It2 b2, It2 e2,
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
                    upd(u, *jt, weight);
                }
            }
        } else {
            // 不同范围：笛卡尔积
            for (auto it = b1; it != e1; ++it) {
                auto u = *it;
                for (auto jt = b2; jt != e2; ++jt) {
                    upd(u, *jt, weight);
                }
            }
        }
    }

    template<
        typename Range1, typename Range2,
        typename UpdateFunc
    >
    inline void processEdgePairs(const Range1 &r1,
                                 const Range2 &r2,
                                 double weight,
                                 UpdateFunc &&upd) noexcept {
        processEdgePairsImpl(
            std::begin(r1), std::end(r1),
            std::begin(r2), std::end(r2),
            weight,
            std::forward<UpdateFunc>(upd)
        );
    }

    template<
        typename Range,
        typename UpdateFunc
    >
    inline void processEdgePairs(const Range &r,
                                 double weight,
                                 UpdateFunc &&upd) noexcept {
        processEdgePairsImpl(
            std::begin(r), std::end(r),
            std::begin(r), std::end(r),
            weight,
            std::forward<UpdateFunc>(upd)
        );
    }


    struct CompareLeaf {
        const std::vector<double> &coreLeaf; // 指向外部数组
        explicit CompareLeaf(const std::vector<double> &coreLeaf) : coreLeaf(coreLeaf) {
        }

        // 注意：这里要返回 “a 排在前面” 的条件，为最小堆写成 coreLeaf[a] > coreLeaf[b]
        bool operator()(daf::Size const &a, daf::Size const &b) const {
            return coreLeaf[a] > coreLeaf[b];
        }
    };

    using DHeap = boost::heap::d_ary_heap<
        daf::Size,
        boost::heap::arity<8>,
        boost::heap::mutable_<true>,
        boost::heap::compare<CompareLeaf>
    >;


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


    std::pair<double *, daf::Size *> countingPerEdge(const DynamicGraph<TreeGraphNode> &treeGraph,
                                                     const Graph &edgeGraph,
                                                     const daf::CliqueSize k) {
        // EdgeHashMap<double> coreE(edgeGraph.adj_list.size());
        double *countingE = new double[edgeGraph.adj_list.size()];
        daf::Size *degreeE = new daf::Size[edgeGraph.adj_list.size()];
        memset(countingE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));
        memset(degreeE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));
        daf::StaticVector<daf::Size> povit;
        daf::StaticVector<daf::Size> keepC;
        daf::Size count = 0;
        for (const auto &clique: treeGraph.adj_list) {
            povit.clear();
            keepC.clear();
            if (clique.size() < k) {
                continue;
            }
            for (auto &i: clique) {
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
            if (needPivot >= 0 && needPivot <= int(povit.size())) {
                totalKcliques = nCr[povit.size()][needPivot];
                for (size_t i = 0; i < keepC.size(); ++i) {
                    for (size_t j = i + 1; j < keepC.size(); ++j) {
                        daf::Size u = keepC[i], v = keepC[j];
                        auto index = edgeGraph.getEdgeIndex(u, v);
                        countingE[index] += totalKcliques;
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
                        countingE[index] += eachPivotKcliques;
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
                        countingE[index] += eachKeepPivotKcliques;
                        degreeE[index]++;
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

    std::vector<double> initLeafCore(const DynamicGraph<TreeGraphNode> &tree, double * &coreE, daf::Size k,
                                     const Graph &edgeGraph) {
        // init as the min one in the edge
        // daf::em
        // daf::enumerateCombinations()
        // memset(coreLeaf, std::numeric_limits<double>::max(), sizeof(double) * leafList.size());
        std::vector<double> leafCore(tree.adj_list.size());
        const daf::Size numLeaf = tree.adj_list.size();
        for (daf::Size i = 0; i < numLeaf; ++i) {
            auto leaf = tree.adj_list[i];
            // TODO: add lowerBound
            // double lowerBound = nCr[povit.size() - 2][k - 2];


            // std::cout << "leaf: " << leaf->leafId << " keepC: " << keepC << " povit: " << povit
            //           << " k: " << k << std::endl;
            double minCore = std::numeric_limits<double>::max();

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
            leafCore[i] = minCore;
        }
        return leafCore;
    }
}

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> >  PlusNucleusEdgeCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraph<daf::Size> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();
    auto [countingKE, degreeERemove] = PlusECD::countingPerEdge(tree, edgeGraph, k);

    std::vector<double> leafCore = PlusECD::initLeafCore(tree, countingKE, k, edgeGraph);

#ifndef NDEBUG
    tree.printGraphPerV();
    daf::printArray(countingKE, edgeGraph.adj_list.size());
    PlusECD::printEdgeCore(edgeGraph, countingKE);
    std::cout << "leafCore: " << leafCore << std::endl;
    // PlusECD::printEdgeCore(edgeGraph, degreeE);
#endif


    std::vector<bool> removedLeaf(tree.adj_list.size());
    // memset(removedLeaf.getData(), false, tree.adj_list.size() * sizeof(bool));

    auto *coreE = new double[edgeGraph.adj_list.size()];
    memset(coreE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));

    // auto *degreeV = new daf::Size[treeGraphV.adj_list.size()];
    // // std::cout << "adj_list_offsets: " << treeGraphV.adj_list_offsets << std::endl;
    // for (daf::Size i = 0; i < treeGraphV.adj_list.size(); ++i) {
    //     degreeV[i] = treeGraphV.adj_list[i].size();
    // }


    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keepC;
    daf::StaticVector<daf::Size> newPivot;
    daf::StaticVector<daf::Size> newKeepC;
    daf::StaticVector<daf::Size> removedPovit;
    daf::StaticVector<bool> isRemovedV(treeGraphV.adj_list.size());
    isRemovedV.c_size = treeGraphV.adj_list.size();
    memset(isRemovedV.getData(), false, treeGraphV.adj_list.size() * sizeof(bool));

    daf::StaticVector<daf::Size> currentRemoveLeafIds(tree.adj_list.size());
    bool removedKeepC = false;


    // daf::StaticVector<std::pair<daf::Size, double> > updateLeaf(tree.adj_list.size() * 10);
    std::map<daf::Size, double> updateLeaf;


    double currCore = 0;
    PlusECD::DHeap heap{PlusECD::CompareLeaf(leafCore)};
    std::vector<PlusECD::DHeap::handle_type> heapHandles(tree.adj_list.size());

    for (daf::Size i = 0; i < tree.adj_list.size(); ++i) {
        auto leaf = tree.adj_list[i];
        if (leaf.size() < k) {
            std::cerr << "Error: leaf id is not equal to index" << std::endl;
        }
        heapHandles[i] = heap.push(i);
    }
    std::cout << "=========================begin=========================" << std::endl;
    daf::StaticVector<std::pair<daf::Size, daf::Size> > removedEdges(1000);
    while (!heap.empty()) {
        double minCore = leafCore[heap.top()];
        // 一次循环把所有 core==minCore 的 leaf 全部 pop 出来
        while (!heap.empty() && leafCore[heap.top()] == minCore) {
            auto id = heap.top();
            removedLeaf[id] = true;
            heap.pop();
            currentRemoveLeafIds.push_back(id);
        }

        // for (auto currLeafId: currentRemoveLeafIds) {
        //     auto leaf = tree.adj_list[currLeafId];
        //     for (auto i = 0; i < leaf.size(); ++i) {
        //         --degreeV[leaf[i].v];
        //     }
        // }
        std::cout << "currentCore: " << currCore << std::endl;
#ifndef NDEBUG
        std::cout << "currentRemoveLeafIds: " << currentRemoveLeafIds << std::endl;
        std::cout << "currentRemoveLeafIds: " << currentRemoveLeafIds << std::endl;
#endif
        for (auto leafId: currentRemoveLeafIds) {
        //     // if (leafId == 2677) {
        //     //     std::cerr << "Error: leaf id is not found" << std::endl;
        //     // }
            auto leaf = tree.adj_list[leafId];
        // std::cout << leafId << " leaf: " << leaf << std::endl;
            currCore = std::max(currCore, leafCore[leafId]);
            for (auto node: leaf) {
                if (node.isPivot) {
                    povit.push_back(node.v);
                } else {
                    keepC.push_back(node.v);
                }
            }


            // removedPovit.print("removedPovit1: ");
            // double totalKcliques = -1;totalKcliques = nCr[povit.size()][needPivot];


            daf::Size needPivot = k - keepC.size();
            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            for (auto i = 0; i < keepC.size() && !removedKeepC; ++i) {
                for (auto j = i + 1; j < keepC.size() && !removedKeepC; ++j) {
                    auto u = keepC[i], v = keepC[j];
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    if (countingKE[idx] <= currCore) {
                        removedKeepC = true;
                        coreE[idx] = currCore;
                    }
                }
            }

            for (daf::Size i = 0; i < povit.size(); ++i) {
                auto u = povit[i];
                if (isRemovedV[u]) continue;
                for (auto j = 0; j < keepC.size(); ++j) {
                    auto v = keepC[j];
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    if (countingKE[idx] <= currCore) {
                        removedPovit.push_back(u);
                        isRemovedV[u] = true;
                        break;
                    }
                    coreE[idx] = currCore;
                }
            }

            int removePPEdgeCount = 0;
            for (daf::Size i = 0; i < povit.size(); ++i) {
                for (daf::Size j = i + 1; j < povit.size(); ++j) {
                    auto u = povit[i], v = povit[j];
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    if (countingKE[idx] <= currCore && !isRemovedV[u] && !
                        isRemovedV
                        [v]) {
                        removedEdges.push_back(std::make_pair(u, v));
                        ++removePPEdgeCount;
                    }
                    coreE[idx] = currCore;
                }
            }
            // std::cout << "coreE: " << std::endl;
            // daf::printArray(coreE, edgeGraph.adj_list.size());
            if (removedKeepC || needPivot > povit.size() - removedPovit.size()) {
                auto removeW = [&](daf::Size u, daf::Size v, double w) {
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    countingKE[idx] -= w;
                    if (countingKE[idx] <= currCore) {
                        coreE[idx] = currCore;
                    }
                    auto uIter = treeGraphV.getNbr(u);
                    auto vIter = treeGraphV.getNbr(v);
                    daf::intersect_with_callback(uIter->begin(), uIter->end(),
                                                 vIter->begin(), vIter->end(),
                                                 [&](const daf::Size &x) {
                                                     if (leafCore[x] > countingKE[idx] && !removedLeaf[x]) {
                                                         auto pvre = updateLeaf.find(x);
                                                         if (pvre == updateLeaf.end()) {
                                                             updateLeaf[x] = countingKE[idx];
                                                         } else {
                                                             pvre->second = std::min(pvre->second, countingKE[idx]);
                                                         }
                                                     }
                                                 });
                };
                double KtoK = 0;
                double KtoP = 0;
                double PtoP = 0;
                if (needPivot <= povit.size()) {
                    KtoK = nCr[povit.size()][needPivot];
                    PlusECD::processEdgePairs(keepC, KtoK, removeW);
                }
                int needPP = int(needPivot) - 2;
                if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                    PtoP = nCr[povit.size() - 2][needPP];
                    PlusECD::processEdgePairs(povit, PtoP, removeW);
                }
                int needKP = int(needPivot) - 1;
                if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                    KtoP = nCr[povit.size() - 1][needKP];
                    PlusECD::processEdgePairs(keepC, povit, KtoP, removeW);
                }
                tree.removeNode(leafId);
                removedLeaf[leafId] = true;
                for (auto i: leaf) {
                    treeGraphV.removeNbr(i.v, leafId);
                }
            } else if (!removedEdges.empty()) {
                // std::cout << "============================================================" << std::endl;
                //     std::cout << "removedEdges: " << removedEdges << std::endl;

                double KtoK = 0;
                double KtoP = 0;
                double PtoP = 0;


                auto addW = [&](daf::Size u, daf::Size v, double w) {
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    countingKE[idx] += w;
                };
                auto initCore = [&](const std::vector<TreeGraphNode> &leaf) {
                    for (auto i: leaf) {
                        if (i.isPivot) {
                            newPivot.push_back(i.v);
                        } else {
                            newKeepC.push_back(i.v);
                        }
                    }

                    daf::Size needPivot = k - newKeepC.size();
                    if (needPivot <= newPivot.size() && newKeepC.size() > 1) {
                        KtoK = nCr[newPivot.size()][needPivot];
                        PlusECD::processEdgePairs(newKeepC, KtoK, addW);
                    }
                    int needPP = int(needPivot) - 2;
                    if (0 <= needPP && needPP <= int(newPivot.size()) - 2) {
                        PtoP = nCr[newPivot.size() - 2][needPP];
                        PlusECD::processEdgePairs(newPivot, PtoP, addW);
                    }
                    int needKP = int(needPivot) - 1;
                    if (0 <= needKP && needKP <= int(newPivot.size()) - 1) {
                        KtoP = nCr[newPivot.size() - 1][needKP];
                        PlusECD::processEdgePairs(newKeepC, newPivot, KtoP, addW);
                    }
                    newPivot.clear();
                    newKeepC.clear();
                };

                // if (!removedPovit.empty() && needPivot <= povit.size() - removedPovit.size())
                daf::StaticVector<daf::Size> newLeafIds;
                auto &newLeaf = leaf;
                if (!removedPovit.empty()) {
                    newLeaf = tree.removeNbrs(leafId, removedPovit);
                }
                bk::bronKerbosch(newLeaf, removedEdges, k,
                                 [&](const bk::Bitset &c, const bk::Bitset &pivots) {
                                     std::vector<TreeGraphNode> newSubLeaf = bk::coverToVertex(c, pivots, leaf);
                                     auto newId = tree.addNode(newSubLeaf);
                                     // std::cout << "newSubLeaf: " << tree.adj_list[newId] << " newId: " << newId << std::endl;
                                     // for (auto i: newLeaf) {
                                     //     ++degreeV[i.v];
                                     // }
                                     initCore(tree.adj_list[newId]);
                                     newLeafIds.push_back(newId);
                                     if (newId >= leafCore.size()) {
                                         leafCore.resize(newId + 10);
                                         removedLeaf.resize(newId + 10);
                                         heapHandles.resize(newId + 10);
                                     }
                                 }
                );


                auto removeW = [&](daf::Size u, daf::Size v, double w) {
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    countingKE[idx] -= w;

                    if (countingKE[idx] <= currCore) {
                        coreE[idx] = currCore;
                    }
                    auto index = edgeGraph.getEdgeIndex(u, v);
                    auto nbrU = treeGraphV.getNbr(u);
                    auto nbrV = treeGraphV.getNbr(v);
                    daf::intersect_with_callback(nbrU->begin(), nbrU->end(),
                                                 nbrV->begin(), nbrV->end(),
                                                 [&](const daf::Size &x) {
                                                     if (leafCore[x] > countingKE[index] && !removedLeaf[x]) {
                                                         auto prve = updateLeaf.find(x);
                                                         if (prve == updateLeaf.end()) {
                                                             updateLeaf[x] = countingKE[index];
                                                         } else {
                                                             prve->second = std::min(
                                                                 prve->second, countingKE[index]);
                                                         }
                                                     }
                                                 });
                };
                // daf::Size needPivot = k - keepC.size();
                if (needPivot <= povit.size()) {
                    KtoK = nCr[povit.size()][needPivot];
                    PlusECD::processEdgePairs(keepC, KtoK, removeW);
                }
                int needPP = int(needPivot) - 2;
                if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                    PtoP = nCr[povit.size() - 2][needPP];
                    PlusECD::processEdgePairs(povit, PtoP, removeW);
                }
                int needKP = int(needPivot) - 1;
                if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                    KtoP = nCr[povit.size() - 1][needKP];
                    PlusECD::processEdgePairs(keepC, povit, KtoP, removeW);
                }

                for (auto i: newLeafIds) {
                    auto newLeaf = tree.adj_list[i];
                    double minCore = std::numeric_limits<double>::max();
                    for (auto j = 0; j < newLeaf.size(); ++j) {
                        for (auto k = j + 1; k < newLeaf.size(); ++k) {
                            auto u = newLeaf[j], v = newLeaf[k];
                            auto index = edgeGraph.getEdgeIndex(u.v, v.v);
                            minCore = std::min(minCore, countingKE[index]);
                        }
                        treeGraphV.addNbr(newLeaf[j].v, i);
                    }
                    leafCore[i] = minCore;
                    heapHandles[i] = heap.push(i);
                    removedLeaf[i] = false;
                }


                if (removedPovit.empty()) {
                    for (auto v: leaf) {
                        treeGraphV.removeNbr(v.v, leafId);
                    }
                    tree.removeNode(leafId);
                }
                newLeafIds.free();
            } else {
                // only povit removed

                ////////////////////////////////////////////////////////////////////////////////////////////////////////


                double KtoK = 0;
                double KtoP = 0;
                double PtoP = 0;


                double RemovedKtoK = 0;
                double RemovedKtoP = 0;
                double RemovedPtoP = 0;

                if (needPivot <= povit.size()) {
                    KtoK = removedKeepC
                               ? nCr[povit.size()][needPivot]
                               : nCr[povit.size()][needPivot] - nCr[povit.size() - removedPovit.size()][needPivot];
                    RemovedKtoK = nCr[povit.size()][needPivot];
                }
                int needPP = int(needPivot) - 2;
                if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                    PtoP = removedKeepC
                               ? nCr[povit.size() - 2][needPP]
                               : nCr[povit.size() - 2][needPP] - nCr[povit.size() - 2 - removedPovit.size()][needPP];
                    RemovedPtoP = nCr[povit.size() - 2][needPP];
                }
                int needKP = int(needPivot) - 1;
                if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                    KtoP = removedKeepC
                               ? nCr[povit.size() - 1][needKP]
                               : nCr[povit.size() - 1][needKP] - nCr[povit.size() - 1 - removedPovit.size()][needKP];
                    RemovedKtoP = nCr[povit.size() - 1][needKP];
                }


                auto minCounting = std::numeric_limits<double>::max();
                auto updateEdgeTemp = [&](daf::Size u, daf::Size v, double w, double rW) {
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    // if (debugSet.contains(u,v) && coreE[idx] <= 188) {
                    //     std::cout << "debugSet: " << u << " " << v << std::endl;
                    // }

                    if (isRemovedV[u] || isRemovedV[v]) {
                        countingKE[idx] -= rW;
                    } else {
                        countingKE[idx] -= w;
                        if (countingKE[idx] > currCore) {
                            minCounting = std::min(minCounting, countingKE[idx]);
                        }
                    }
                    if (countingKE[idx] <= currCore) {
                        coreE[idx] = currCore;
                    }
                    auto uIter = treeGraphV.getNbr(u);
                    auto vIter = treeGraphV.getNbr(v);
                    daf::intersect_with_callback(uIter->begin(), uIter->end(),
                                                 vIter->begin(), vIter->end(),
                                                 [&](const daf::Size &x) {
                                                     if (leafCore[x] > countingKE[idx] && !removedLeaf[x]) {
                                                         // updateLeaf.push_back(std::make_pair(x, countingKE[idx]));
                                                         // updateLeaf[x] = countingKE[idx];
                                                         auto prve = updateLeaf.find(x);
                                                         if (prve == updateLeaf.end()) {
                                                             updateLeaf[x] = countingKE[idx];
                                                         } else {
                                                             prve->second = std::min(prve->second, countingKE[idx]);
                                                         }
                                                     }
                                                 });
                };

                // 1) keep–keep
                PlusECD::processEdgePairs(keepC, KtoK, [&](daf::Size u, daf::Size v, double w) {
                    updateEdgeTemp(u, v, w, RemovedKtoK);
                });
                PlusECD::processEdgePairs(povit, PtoP, [&](daf::Size u, daf::Size v, double w) {
                    updateEdgeTemp(u, v, w, RemovedPtoP);
                });
                PlusECD::processEdgePairs(keepC, povit, KtoP, [&](daf::Size u, daf::Size v, double w) {
                    updateEdgeTemp(u, v, w, RemovedKtoP);
                });

                std::ranges::sort(removedPovit);
                if (!removedPovit.empty() && needPivot <= povit.size() - removedPovit.size()) {
                    const auto newLeaf = tree.removeNbrs(leafId, removedPovit);
                    for (auto removedNbr: removedPovit) {
                        treeGraphV.removeNbr(removedNbr, leafId);
                    }
                    if (!newLeaf.empty()) {
                        leafCore[leafId] = minCounting;
                        heapHandles[leafId] = heap.push(leafId);
                        removedLeaf[leafId] = false;
                        // daf::StaticVector<daf::Size> newPovit;
                        // for (auto i: newLeaf) {
                        //     ++degreeV[i.v];
                        // }
                        if (!removedEdges.empty()) {
                            removedEdges.clear();
                        }
                    } else {
                        leafCore[leafId] = 0;
                    }
                }
            }


            // std::cout << "updateLeaf: " << std::endl;
            for (auto &u: updateLeaf) {
                daf::Size updateLeafId = u.first;
                double newCore = u.second;
                if (leafCore[updateLeafId] > newCore && !removedLeaf[updateLeafId]) {
                    leafCore[updateLeafId] = newCore;
                    heap.update(heapHandles[updateLeafId]);
                }
            }


#ifndef NDEBUG
            std::cout << "\n\n leafId: " << leafId << " leafCore: " << leafCore[leafId] << " currCore: " << currCore
            << std::endl;
            std::cout << "countingKE: ";
            PlusECD::printEdgeCore(edgeGraph, countingKE);
                // updateLeaf.print("updateLeaf");
                // removedLeaf.print("removedLeaf");
            std::cout << "removedLeaf: " << removedLeaf << std::endl;

            removedPovit.print("removedPovit");
            std::cout << removedKeepC << std::endl;
            // //
            // // std::cout << "degreeE: ";
            // // PlusECD::printEdgeCore(edgeGraph, degreeE);
            // std::cout << "degreeV: ";
            // daf::printArray(degreeV, treeGraphV.adj_list.size());
            //
            std::cout << "coreE: ";
            PlusECD::printEdgeCore(edgeGraph, coreE);
            //
            std::cout << "leafCore: " << leafCore << std::endl;

            // removedLeaf.print("removedLeaf");

            tree.printGraphPerV();
            treeGraphV.printGraphPerV();
            std::cout << "--------------------------------------------------------------------------------------------------------------------"<< std::endl;
#endif


            updateLeaf.clear();
            povit.clear();
            keepC.clear();
            removedEdges.clear();
            removedKeepC = false;
            for (auto v: removedPovit) {
                isRemovedV[v] = false;
            }
            removedPovit.clear();
        }
        currentRemoveLeafIds.clear();
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    // coreE
    // daf::printArray(coreE, edgeGraph.adj_list.size());

#ifndef NDEBUG
    PlusECD::printEdgeCore(edgeGraph, coreE);
#endif

    // /Users/zhangwenqian/UNSW/pivoter/a
    daf::Size numCounting = 0;
    // std::sort(coreE, coreE + edgeGraph.adj_list.size());
    std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > sortedK;
    sortedK.reserve(edgeGraph.adj_list.size());
    // for (daf::Size i = 0; i < edgeGraph.adj_list.size(); i++) {
    //     // cover to int
    //     fprintf(file, "%d\n", (int) coreE[i]);
    //     numCounting += countingKE[i];
    // }

    const daf::Size m = edgeGraph.adj_list.size();
    const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;
    for (daf::Size u = 0; u < n; ++u) {
        const daf::Size start = edgeGraph.adj_list_offsets[u];
        const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
        for (daf::Size idx = start; idx < end; ++idx) {
            // std::cout << "[" << u << ", " << edgeGraph.adj_list[idx] << "] " << coreE[idx] << " ";
            sortedK.emplace_back(
                std::make_pair(std::make_pair(u, edgeGraph.adj_list[idx]), (int) coreE[idx]));
        }
    }

    // for (auto i: sortedK) {
    //     fprintf(file, "%d %d %d\n", i.first.first, i.first.second, i.second);
    // }
    // fclose(file);

    // if (numCounting != 0) {
    //     // exit 1
    //     std::cerr << "Error: numCounting != 0" << std::endl;
    //     std::cerr << "numCounting: " << numCounting << std::endl;
    // }
    assert(numCounting == 0);

    delete[] countingKE;
    // delete[] degreeV;
    delete[] coreE;
    delete[] degreeERemove;
    povit.free();
    keepC.free();
    removedPovit.free();
    newPivot.free();
    newKeepC.free();
    currentRemoveLeafIds.free();
    removedEdges.free();
    // leafIds.free();
    return sortedK;
}


template<class Bitset>
void print_clique(const Bitset &bs) {
    std::cout << '[';
    bool first = true;
    bk::for_each_bit(bs, (int) bs.size(), [&](int v) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << v;
        return true;
    });
    std::cout << "]\n";
}
