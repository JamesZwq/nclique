//
// Created by 张文谦 on 25-3-4.
//

#include "NCliqueCoreDecomposition.h"
#include <boost/heap/d_ary_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <set>

#include "../BK/BronKerboschRmEdge.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "debug/EdgeSet.h"
#include "graph/DynamicBipartiteGraph.hpp"
// #include "graph/DynamicGraph.h"
#include "dataStruct/coreDisJoin.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];
// 放在你的函数外（比如文件顶部），保证编译时可见并内联
// #ifndef NDEBUG
// set NOEBUG as trus


namespace VCD {
    struct LeafRmInfo {
        bool removedKeepC;
        daf::StaticVector<daf::Size> removedPivots{0};

        LeafRmInfo() : removedKeepC(false) {}
        bool empty() const {
            return !removedKeepC && removedPivots.empty();
        }

        void init(auto capacity = 400) {
            removedKeepC = false;
            removedPivots.reserve(capacity);
        }

        void clear() {
            removedKeepC = false;
            removedPivots.clear();
        }

        friend std::ostream &operator<<(std::ostream &os, const LeafRmInfo &info) {
            os << "removedKeepC: " << info.removedKeepC << "\n removedPivots: " << info.removedPivots;
            return os;
        }
    };
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


    struct CompareVertex {
        const double *vertexCounting; // 指向外部数组
        explicit CompareVertex(const double *coreLeaf) : vertexCounting(coreLeaf) {
        }

        // 注意：这里要返回 “a 排在前面” 的条件，为最小堆写成 coreLeaf[a] > coreLeaf[b]
        bool operator()(daf::Size const &a, daf::Size const &b) const {
            return vertexCounting[a] > vertexCounting[b];
        }
    };

    using DHeap = boost::heap::d_ary_heap<
        daf::Size,
        boost::heap::arity<8>,
        boost::heap::mutable_<true>,
        boost::heap::compare<CompareVertex>
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
        std::memset(core, 0, tree.getRoot()->children.size() * sizeof(double));
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
                        auto index = edgeGraph.getEdgeCompressedId(u, v);
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
                        auto index = edgeGraph.getEdgeCompressedId(u, v);
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
                        auto index = edgeGraph.getEdgeCompressedId(u, v);
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



    // std::pair<double *, daf::Size *> countingPerEdge(const MultiBranchTree &tree, const Graph &edgeGraph,
    //                                                  const daf::CliqueSize k) {
    //     // EdgeHashMap<double> coreE(edgeGraph.adj_list.size());
    //     double *coreE = new double[edgeGraph.adj_list.size()];
    //     daf::Size *degreeE = new daf::Size[edgeGraph.adj_list.size()];
    //     memset(coreE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));
    //     memset(degreeE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));
    //     daf::StaticVector<daf::Size> povitC;
    //     daf::StaticVector<daf::Size> keepC;
    //     daf::Size count = 0;
    //     for (auto node: tree.getRoot()->children) {
    //         if (node->MaxDeep < k) {
    //             continue;
    //         }
    //         keepC.push_back(node->v);
    //         countingPerEdgeHelp(*node, k, edgeGraph, coreE, degreeE, povitC, keepC);
    //         keepC.pop_back();
    //     }
    //     return {coreE, degreeE};
    // }

    double * countingPerVertex(const DynamicGraph<TreeGraphNode> &treeGraph,
                                                     const Graph &edgeGraph,
                                                     const daf::CliqueSize k) {
        double *countingV = new double[edgeGraph.adj_list_offsets.size()];
        memset(countingV, 0, edgeGraph.adj_list_offsets.size() * sizeof(double));
        daf::StaticVector<daf::Size> povit;
        daf::StaticVector<daf::Size> keepC;
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
            for (const auto &v: keepC) {
                countingV[v] += nCr[povit.size()][needPivot];
            }
            double eachPivotKcliques = 0;
            const int needPivotWithV = needPivot - 1;
            if (needPivotWithV >= 0 && needPivotWithV <= static_cast<int>(povit.size()) - 1) {
                eachPivotKcliques = nCr[povit.size() - 1][needPivotWithV];
            }
            for (const auto &v: povit) {
                countingV[v] += eachPivotKcliques;
            }
        }
        povit.free();
        keepC.free();
        return countingV;
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
        memset(coreV, 0, n * sizeof(double));
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

    daf::StaticVector<daf::Size> initLeafCore(const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph) {
        // init as the min one in the edge
        // daf::em
        // daf::enumerateCombinations()
        // memset(coreLeaf, std::numeric_limits<double>::max(), sizeof(double) * leafList.size());
        daf::StaticVector<daf::Size> leafCore(tree.adj_list.size() * 1.5);
        leafCore.c_size = tree.adj_list.size() * 1.5;
        const daf::Size numLeaf = tree.adj_list.size();
        for (daf::Size i = 0; i < numLeaf; ++i) {
            auto leaf = tree.adj_list[i];
            for (auto node : leaf) {
                if (!node.isPivot) {
                    leafCore[i] = edgeGraph.coreV[node.v];
                    break;
                }
            }
        }
        return leafCore;
    }
}

double *  NCliqueVertexCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    auto countingV = VCD::countingPerVertex(tree, edgeGraph, k);
    auto coreV = new double[edgeGraph.adj_list_offsets.size()];
    memset(coreV, 0, edgeGraph.adj_list_offsets.size() * sizeof(double));
    // std::vector<double> leafCore = VCD::initLeafCore(tree, countingKE, k, edgeGraph);

#ifndef NDEBUG
    tree.printGraphPerV();
    daf::printArray(countingV, edgeGraph.adj_list.size());
    // VCD::printEdgeCore(edgeGraph, degreeE);
#endif

    // memset(removedLeaf.getData(), false, tree.adj_list.size() * sizeof(bool));

    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keepC;
    daf::StaticVector<daf::Size> newPivot;
    daf::StaticVector<daf::Size> newKeepC;

    daf::StaticVector<daf::Size> currentRemoveVertexIds(edgeGraph.adj_list_offsets.size());

    daf::StaticVector<bool> vertexInHeap(edgeGraph.adj_list_offsets.size());
    vertexInHeap.c_size = edgeGraph.adj_list_offsets.size();
    memset(vertexInHeap.getData(), true, edgeGraph.adj_list_offsets.size() * sizeof(bool));

    // daf::StaticVector<std::pair<daf::Size, double> > updateLeaf(tree.adj_list.size() * 10);
    daf::StaticVector<daf::Size> removedLeaf(tree.adj_list.size());
    daf::StaticVector<VCD::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();

    VCD::DHeap heap{VCD::CompareVertex(countingV)};
    std::vector<VCD::DHeap::handle_type> heapHandles(edgeGraph.adj_list_offsets.size() - 1);
    for (daf::Size i = 0; i < edgeGraph.adj_list_offsets.size() - 1; ++i) {
        heapHandles[i] = heap.push(i);
    }
    // std::cout << "tree: ";
    // tree.printGraphPerV();
    //
    // std::cout << "treeGraphV: ";
    // treeGraphV.printGraphPerV();
#ifndef NDEBUG

    std::cout << "countingV: ";
    daf::printArray(countingV, edgeGraph.adj_list_offsets.size() - 1);
    std::cout << "tree: ";
    tree.printGraphPerV();

    std::cout << "treeGraphV: ";
    treeGraphV.printGraphPerV();
#endif
    std::cout << "=========================begin=========================" << std::endl;
    // daf::StaticVector<std::pair<daf::Size, daf::Size> > removedEdges(1000);
    double minCore = 0;
    int numProgress = 0;
    CoreDisJoin hierarchyBuilder(edgeGraph.adj_list_offsets.size() - 1, 20);
    while (!heap.empty()) {

        minCore = std::max( countingV[heap.top()], minCore );
        // 一次循环把所有 core==minCore 的 leaf 全部 pop 出来
        // printf("minCore: %.2f, heap size: %zu\n", minCore, heap.size());

        std::cout << "minCore: " << minCore
        << " heap size: " << heap.size()
        << " num Leaf: " << tree.size() << " "
        << k << "-Clique count: " << tree.cliqueCount(k)
        << std::endl;
        while (!heap.empty() && countingV[heap.top()] <= minCore) {
            auto id = heap.top();
            vertexInHeap[id] = false;
            heap.pop();
            currentRemoveVertexIds.push_back(id);
            coreV[id] = minCore;
            // daf::printProgress(numProgress++, edgeGraph.adj_list.size());
        }

        // std::cout << "currentCore: " << currCore << std::endl;
#ifndef NDEBUG
        std::cout << "currentRemoveVertexIds: " << currentRemoveVertexIds << std::endl;
#endif
        for (auto v: currentRemoveVertexIds) {
            auto &adjClique = treeGraphV.getNbr(v);
            for (const auto &clique: adjClique) {
                if (leafRmInfo[clique.v].empty()) {
                    removedLeaf.push_back(clique.v);
                    leafRmInfo[clique.v].init(tree.adj_list[clique.v].size());
                }
                if (leafRmInfo[clique.v].removedKeepC) {
                    continue;
                }
                if (!clique.isPivot) {
                    leafRmInfo[clique.v].removedKeepC = true;
                } else {
                    leafRmInfo[clique.v].removedPivots.push_back(v);
                }
            }
        }
        // removedLeaf.print("removedLeaf");
        // for (auto leafId : removedLeaf) {
        for (auto leafIdIdx = 0; leafIdIdx < removedLeaf.size(); ++leafIdIdx) {
            auto leafId = removedLeaf[leafIdIdx];
            auto leaf = tree.adj_list[leafId];
            // std::cout << "leafId: " << leafId << " leaf: " << leaf << std::endl;
            VCD::LeafRmInfo &leafRm = leafRmInfo[leafId];
            if (leaf.empty()) {
                std::cout << "Error: newLeaf is empty" << std::endl;
                std::cout << "leafId: " << leafId << std::endl;
                std::cout << "leafRm.removedPivots: " << leafRm.removedPivots << std::endl;
                std::cout << "leaf: " << leaf << std::endl;
                std::cerr << "Error: leaf is empty" << std::endl;
                std::exit(1);
            }

            std::ranges::sort(leafRm.removedPivots);
            leafRm.removedPivots.unique();
#ifndef NDEBUG
            std::cout << leafId << " leaf: " << leaf << "\n leafRm: " << leafRm << std::endl;
#endif
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
            // std::cout << "coreE: " << std::endl;
            // daf::printArray(coreE, edgeGraph.adj_list.size());
            if (leafRm.removedKeepC || needPivot > povit.size() - leafRm.removedPivots.size()) {
                // auto removeW = [&](daf::Size u, daf::Size v, double w) {
                //     auto idx = edgeGraph.getEdgeCompressedId(u, v);
                //     countingV[idx] -= w;
                //     if (edgeInHeap[idx]) {
                //         heap.update(heapHandles[idx]);
                //     }
                // };
                // -----------------------------------------------------------------

                auto keepValue = nCr[povit.size()][needPivot];
                for (int i = 0; i < keepC.size(); ++i) {
                    countingV[keepC[i]] -= keepValue;
                    countingV[keepC[i]] = std::max(countingV[keepC[i]], 0.0);
                    if (vertexInHeap[keepC[i]]) {
                        heap.update(heapHandles[keepC[i]]);
                    }
                }

                if (needPivot > 0) {
                    auto pivotValue = nCr[povit.size() - 1][needPivot - 1];
                    for (int i = 0; i < povit.size(); ++i) {
                        countingV[povit[i]] -= pivotValue;
                        countingV[povit[i]] = std::max(countingV[povit[i]], 0.0);
                        if (vertexInHeap[povit[i]]) {
                            heap.update(heapHandles[povit[i]]);
                        }
                    }
                }
                for (auto i: leaf) {
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                }
                tree.removeNode(leafId);
            } else {


                auto KeepCount = nCr[povit.size()][needPivot] - nCr[povit.size() - leafRm.removedPivots.size()][needPivot];
                if (!leafRm.removedPivots.empty() && needPivot <= povit.size() - leafRm.removedPivots.size()) {
                    for (auto removedNbr: leafRm.removedPivots) {
                        treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                    }
                    const auto newLeaf = tree.removeNbrs(leafId, leafRm.removedPivots);


                    auto PivotCount = nCr[povit.size() - 1][needPivot - 1] - nCr[povit.size() - leafRm.removedPivots.size() - 1][needPivot - 1];
                    auto RemovedPivot = nCr[povit.size() - 1][needPivot - 1];
                    for (auto removedPivot : leafRm.removedPivots) {
                        // check(countingV[removedPivot], RemovedPivot);
                        countingV[removedPivot] -= RemovedPivot;
                        countingV[removedPivot] = std::max(countingV[removedPivot], 0.0);
                        if (vertexInHeap[removedPivot]) {
                            heap.update(heapHandles[removedPivot]);
                        }
                        // std::cout << "removedPivot: " << removedPivot << " RemovedPivot: " << RemovedPivot << std::endl;
                    }

                    for (auto v : newLeaf) {
                        if (v.isPivot) {
                            // check(countingV[v.v], PivotCount);
                            countingV[v.v] -= PivotCount;
                            countingV[v.v] = std::max(countingV[v.v], 0.0);
                        } else {
                            // check(countingV[v.v], KeepCount);
                            countingV[v.v] -= KeepCount;
                            countingV[v.v] = std::max(countingV[v.v], 0.0);
                        }
                        if (vertexInHeap[v.v]) {
                            heap.update(heapHandles[v.v]);
                        }
                    }

                }
            }


// #ifndef NDEBUG
//             std::cout << "countingKE: ";
//             VCD::printEdgeCore(edgeGraph, countingKE);
//                 // updateLeaf.print("updateLeaf");
//                 // removedLeaf.print("removedLeaf");
//             std::cout << "removedLeaf: " << removedLeaf << std::endl;
//
//             // //
//             // // std::cout << "degreeE: ";
//             // // VCD::printEdgeCore(edgeGraph, degreeE);
//             // std::cout << "degreeV: ";
//             // daf::printArray(degreeV, treeGraphV.adj_list.size());
//             //
//             std::cout << "coreE: ";
//             VCD::printEdgeCore(edgeGraph, coreE);
//             //
//
//             // removedLeaf.print("removedLeaf");
//
//             std::cout << "--------------------------------------------------------------------------------------------------------------------"<< std::endl;
            // #endif
#ifndef NDEBUG
            std::cout << "coreV: " << std::endl;
            daf::printArray(coreV, edgeGraph.adj_list_offsets.size() - 1);
            std::cout << "countingV: ";
            daf::printArray(countingV, edgeGraph.adj_list_offsets.size() - 1);
            std::cout << "tree: ";
            tree.printGraphPerV();

            std::cout << "treeGraphV: ";
            treeGraphV.printGraphPerV();
#endif

            leafRmInfo[leafId].clear();
            povit.clear();
            keepC.clear();
            // removedEdges.clear();
            // removedKeepC = false;
            // for (auto v: removedPovit) {
            //     isRemovedV[v] = false;
            // }
            // removedPovit.clear();
        }
        currentRemoveVertexIds.clear();
        removedLeaf.clear();
        // currentRemoveLeafIds.clear();
#ifndef NDEBUG
        std::cout << "coreV: " << std::endl;
        daf::printArray(coreV, edgeGraph.adj_list_offsets.size() - 1);
        std::cout << "countingV: ";
        daf::printArray(countingV, edgeGraph.adj_list_offsets.size() - 1);
        std::cout << "tree: ";
        tree.printGraphPerV();

        std::cout << "treeGraphV: ";
        treeGraphV.printGraphPerV();
#endif

    }
    //
    // for (auto i = 0;  i < edgeGraph.adj_list_offsets.size(); ++i) {
    //     if (countingV[i] != 0) {
    //         std::cerr << "Error: countingKE != 0" << std::endl;
    //         std::cerr << "countingKE: " << countingV[i] << " at index " << i << std::endl;
    //         std::exit(1);
    //     }
    // }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;


    delete[] countingV;
    // delete[] degreeV;
    povit.free();
    keepC.free();
    newPivot.free();
    newKeepC.free();
    currentRemoveVertexIds.free();
    removedLeaf.free();
    leafRmInfo.free();
    return coreV;
}


template<class Bitset>
void print_clique(const Bitset &bs) {
    std::cout << '[';
    bool first = true;
    bkRmEdge::for_each_bit(bs, (int) bs.size(), [&](int v) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << v;
        return true;
    });
    std::cout << "]\n";
}