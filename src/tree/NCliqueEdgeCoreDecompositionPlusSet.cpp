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
#include "graph/DynamicBipartiteGraph.hpp"
#include "graph/DynamicGraph.h"

class DynamicBipartiteGraph;
extern double nCr[1001][401];
// 放在你的函数外（比如文件顶部），保证编译时可见并内联
// #ifndef NDEBUG
// set NOEBUG as trus


namespace PlusECDSet {
    struct LeafRmInfo {
        bool removedKeepC;
        daf::StaticVector<daf::Size> removedPivots{0};
        daf::StaticVector<std::pair<daf::Size, daf::Size>> removedEdges{0};

        LeafRmInfo() : removedKeepC(false) {}
        bool empty() const {
            return !removedKeepC && removedPivots.empty() && removedEdges.empty();
        }

        void init() {
            removedKeepC = false;
            removedPivots.reserve(400);
            removedEdges.reserve(400);
        }

        void clear() {
            removedKeepC = false;
            removedPivots.clear();
            removedEdges.clear();
        }

        ~LeafRmInfo() {
            removedEdges.free();
            removedPivots.free();
        }

        friend std::ostream &operator<<(std::ostream &os, const LeafRmInfo &info) {
            os << "removedKeepC: " << info.removedKeepC << "\n removedPivots: " << info.removedPivots
               << "\n removedEdges: " << info.removedEdges;
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


    struct CompareEdge {
        const double *edgeCounting; // 指向外部数组
        explicit CompareEdge(const double *coreLeaf) : edgeCounting(coreLeaf) {
        }

        // 注意：这里要返回 “a 排在前面” 的条件，为最小堆写成 coreLeaf[a] > coreLeaf[b]
        bool operator()(daf::Size const &a, daf::Size const &b) const {
            return edgeCounting[a] > edgeCounting[b];
        }
    };

    using DHeap = boost::heap::d_ary_heap<
        daf::Size,
        boost::heap::arity<8>,
        boost::heap::mutable_<true>,
        boost::heap::compare<CompareEdge>
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

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet(
    DynamicGraphSet<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();
    auto [countingKE, degreeERemove] = PlusECDSet::countingPerEdge(tree, edgeGraph, k);

    // std::vector<double> leafCore = PlusECDSet::initLeafCore(tree, countingKE, k, edgeGraph);

#ifndef NDEBUG
    tree.printGraphPerV();
    daf::printArray(countingKE, edgeGraph.adj_list.size());
    PlusECDSet::printEdgeCore(edgeGraph, countingKE);
    // PlusECDSet::printEdgeCore(edgeGraph, degreeE);
#endif

    // memset(removedLeaf.getData(), false, tree.adj_list.size() * sizeof(bool));

    auto *coreE = new double[edgeGraph.adj_list.size()];
    memset(coreE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));

    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keepC;
    daf::StaticVector<daf::Size> newPivot;
    daf::StaticVector<daf::Size> newKeepC;
    // daf::StaticVector<daf::Size> removedPovit;
    // daf::StaticVector<bool> isRemovedV(edgeGraph.adj_list_offsets.size());
    // isRemovedV.c_size = edgeGraph.adj_list_offsets.size();
    // memset(isRemovedV.getData(), false, edgeGraph.adj_list_offsets.size() * sizeof(bool));

    daf::StaticVector<daf::Size> currentRemoveEdgeIds(edgeGraph.adj_list.size());

    daf::StaticVector<bool> edgeInHeap(edgeGraph.adj_list.size());
    edgeInHeap.c_size = edgeGraph.adj_list.size();
    memset(edgeInHeap.getData(), true, edgeGraph.adj_list.size() * sizeof(bool));

    // daf::StaticVector<std::pair<daf::Size, double> > updateLeaf(tree.adj_list.size() * 10);
    std::map<daf::Size, double> updateLeaf;

    daf::StaticVector<daf::Size> removedLeaf(tree.adj_list.size());
    daf::StaticVector<PlusECDSet::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();


    double currCore = 0;
    PlusECDSet::DHeap heap{PlusECDSet::CompareEdge(countingKE)};
    std::vector<PlusECDSet::DHeap::handle_type> heapHandles(edgeGraph.adj_list.size());

    for (daf::Size i = 0; i < edgeGraph.adj_list.size(); ++i) {
        heapHandles[i] = heap.push(i);
    }
    // std::cout << "tree: ";
    // tree.printGraphPerV();
    //
    // std::cout << "treeGraphV: ";
    // treeGraphV.printGraphPerV();
#ifndef NDEBUG
    std::cout << "coreE: ";
    PlusECDSet::printEdgeCore(edgeGraph, coreE);
    std::cout << "countingKE: ";
    PlusECDSet::printEdgeCore(edgeGraph, countingKE);
    std::cout << "tree: ";
    tree.printGraphPerV();

    std::cout << "treeGraphV: ";
    treeGraphV.printGraphPerV();
#endif
    std::cout << "=========================begin=========================" << std::endl;
    // daf::StaticVector<std::pair<daf::Size, daf::Size> > removedEdges(1000);
    double minCore = 0;
    int numProgress = 0;
    while (!heap.empty()) {

        minCore = std::max( countingKE[heap.top()], minCore );
        // 一次循环把所有 core==minCore 的 leaf 全部 pop 出来
        while (!heap.empty() && countingKE[heap.top()] <= minCore) {
            auto id = heap.top();
            edgeInHeap[id] = false;
            heap.pop();
            currentRemoveEdgeIds.push_back(id);
            coreE[id] = minCore;
            // std::cout << "progress: " << numProgress++ << "/" << edgeGraph.adj_list.size() << std::flush;
            daf::printProgress(numProgress++, edgeGraph.adj_list.size());
        }

        currCore = minCore;

        // std::cout << "currentCore: " << currCore << std::endl;
#ifndef NDEBUG
        std::cout << "currentRemoveEdge: " << std::endl;
#endif
        for (auto edgeId: currentRemoveEdgeIds) {
            //     // if (leafId == 2677) {
            //     //     std::cerr << "Error: leaf id is not found" << std::endl;
            //     // }
            auto [edgeU, edgeV] = edgeGraph.getEdgeById(edgeId);
            #ifndef NDEBUG
                        std::cout << edgeU << " " << edgeV << std::endl;
            #endif
            std::vector<TreeGraphNode> &adjEdgeU = treeGraphV.getNbr(edgeU);
            std::vector<TreeGraphNode> &adjEdgeV = treeGraphV.getNbr(edgeV);

            daf::intersect_with_callback_pair(adjEdgeU, adjEdgeV, [&](const TreeGraphNode &uClique,const TreeGraphNode &vClique) {
                if (leafRmInfo[uClique.v].empty()) {
                    removedLeaf.push_back(uClique.v);
                    leafRmInfo[uClique.v].init();
                }
                if (leafRmInfo[uClique.v].removedKeepC) return;
                if (!uClique.isPivot && !vClique.isPivot) {
                    leafRmInfo[uClique.v].removedKeepC = true;
                } else if (uClique.isPivot && vClique.isPivot) {
                    leafRmInfo[uClique.v].removedEdges.push_back({edgeU, edgeV});
                } else if (uClique.isPivot && !vClique.isPivot) {
                    leafRmInfo[uClique.v].removedPivots.push_back(edgeU);
                } else if (!uClique.isPivot && vClique.isPivot) {
                    leafRmInfo[uClique.v].removedPivots.push_back(edgeV);
                }
            });
        }
        // removedLeaf.print("removedLeaf");
        for (auto leafId : removedLeaf) {
            auto leaf = tree.adj_list[leafId];


            if (leaf.empty()) {
                // std::cerr << "Error: newLeaf is empty" << std::endl;
                // std::cout << "leafId: " << leafId << std::endl;
                // std::cout << "leafRm.removedPivots: " << leafRm.removedPivots << std::endl;
                // std::cout << "leafRm.removedEdges: " << leafRm.removedEdges << std::endl;
                // std::cout << "leaf: " << leaf << std::endl;
                std::cerr << "Error: leaf is empty" << std::endl;
                std::exit(1);
            }
            PlusECDSet::LeafRmInfo &leafRm = leafRmInfo[leafId];
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


            std::ranges::sort(leafRm.removedPivots);
            leafRm.removedPivots.unique();

            daf::Size needPivot = k - keepC.size();
            // std::cout << "coreE: " << std::endl;
            // daf::printArray(coreE, edgeGraph.adj_list.size());
            if (leafRm.removedKeepC || needPivot > povit.size() - leafRm.removedPivots.size()) {
                auto removeW = [&](daf::Size u, daf::Size v, double w) {
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    countingKE[idx] -= w;
                    if (edgeInHeap[idx]) {
                        heap.update(heapHandles[idx]);
                    }

                    if (u == 6 && v == 7) {
                        std::cout << "asdfasdf: u: " << u << " v: " << v << " idx: " << idx
                                  << " w: " << w << " countingKE[idx]: " << countingKE[idx]
                        << " leafId: " << leafId << " leafRm.removedPivots: " << leafRm.removedPivots
                                  << std::endl;

                    }
                };
                double KtoK = 0;
                double KtoP = 0;
                double PtoP = 0;

                if (needPivot <= povit.size()) {
                    KtoK = nCr[povit.size()][needPivot];
                    PlusECDSet::processEdgePairs(keepC, KtoK, removeW);
                }
                int needPP = int(needPivot) - 2;
                if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                    PtoP = nCr[povit.size() - 2][needPP];
                    PlusECDSet::processEdgePairs(povit, PtoP, removeW);
                }
                int needKP = int(needPivot) - 1;
                if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                    KtoP = nCr[povit.size() - 1][needKP];
                    PlusECDSet::processEdgePairs(keepC, povit, KtoP, removeW);
                }
                tree.removeNode(leafId);
                for (auto i: leaf) {
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                }
            } else if (!leafRm.removedEdges.empty()) {
                // std::cout << "============================================================" << std::endl;
                //     std::cout << "removedEdges: " << removedEdges << std::endl;

                double KtoK = 0;
                double KtoP = 0;
                double PtoP = 0;


                auto addW = [&](daf::Size u, daf::Size v, double w) {
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    countingKE[idx] += w;
                };

                auto initCore = [&](const std::vector<TreeGraphNode> &leaf, const daf::Size &leafId) {
                    for (auto i: leaf) {
                        if (i.isPivot) {
                            newPivot.push_back(i.v);
                            treeGraphV.addNbr(i.v, {leafId, true});
                        } else {
                            newKeepC.push_back(i.v);
                            treeGraphV.addNbr(i.v, {leafId, false});
                        }
                    }
                    daf::Size needPivot = k - newKeepC.size();
                    if (needPivot <= newPivot.size() && newKeepC.size() > 1) {
                        KtoK = nCr[newPivot.size()][needPivot];
                        PlusECDSet::processEdgePairs(newKeepC, KtoK, addW);
                    }
                    int needPP = int(needPivot) - 2;
                    if (0 <= needPP && needPP <= int(newPivot.size()) - 2) {
                        PtoP = nCr[newPivot.size() - 2][needPP];
                        PlusECDSet::processEdgePairs(newPivot, PtoP, addW);
                    }
                    int needKP = int(needPivot) - 1;
                    if (0 <= needKP && needKP <= int(newPivot.size()) - 1) {
                        KtoP = nCr[newPivot.size() - 1][needKP];
                        PlusECDSet::processEdgePairs(newKeepC, newPivot, KtoP, addW);
                    }
                    newPivot.clear();
                    newKeepC.clear();
                };

                // if (!removedPovit.empty() && needPivot <= povit.size() - removedPovit.size())
                // daf::StaticVector<daf::Size> newLeafIds;
                auto &newLeaf = leaf;
                for (auto leafV : leaf) {
                    if (leafV.isPivot) {
                        treeGraphV.removeNbr(leafV.v, {leafId, true});
                    } else {
                        treeGraphV.removeNbr(leafV.v, {leafId, false});
                    }
                }
                std::ranges::sort(leafRm.removedPivots);
                leafRm.removedPivots.unique();
                if (!leafRm.removedPivots.empty()) {
                    // 91528
                    // if (leafRm.removedPivots[0] == 91528) {
                    //     std::cerr << "Error: leafRm.removedPivots[0] == 91528" << std::endl;
                    //     std::cerr << "Error: newLeaf is empty" << std::endl;
                    //     std::cout << "leafId: " << leafId << std::endl;
                    //     std::cout << "leafRm.removedPivots: " << leafRm.removedPivots << std::endl;
                    //     std::cout << "leafRm.removedEdges: " << leafRm.removedEdges << std::endl;
                    //     std::cout << "leaf: " << leaf << std::endl;
                    // }
                    newLeaf = tree.removeNbrs(leafId, leafRm.removedPivots);
                }
                bk::bronKerbosch(newLeaf, leafRm.removedEdges, k,
                                 [&](const bk::Bitset &c, const bk::Bitset &pivots) {
                                     std::vector<TreeGraphNode> newSubLeaf = bk::coverToVertex(c, pivots, leaf);
                                     auto newId = tree.addNode(newSubLeaf);
                                     initCore(tree.adj_list[newId], newId);
                                     if (newId >= leafRmInfo.size()) {
                                         removedLeaf.resize(newId * 1.5);
                                         leafRmInfo.resize(newId * 1.5);
                                         heapHandles.resize(newId * 1.5);
                                     }
                                 }
                );


                auto removeW = [&](daf::Size u, daf::Size v, double w) {
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    countingKE[idx] -= w;
                    if (edgeInHeap[idx]) {
                        heap.update(heapHandles[idx]);
                    }
                };


                // daf::Size needPivot = k - keepC.size();
                if (needPivot <= povit.size()) {
                    KtoK = nCr[povit.size()][needPivot];
                    PlusECDSet::processEdgePairs(keepC, KtoK, removeW);
                }
                int needPP = int(needPivot) - 2;
                if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                    PtoP = nCr[povit.size() - 2][needPP];
                    PlusECDSet::processEdgePairs(povit, PtoP, removeW);
                }
                int needKP = int(needPivot) - 1;
                if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                    KtoP = nCr[povit.size() - 1][needKP];
                    PlusECDSet::processEdgePairs(keepC, povit, KtoP, removeW);
                }
                tree.removeNode(leafId);
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
                    KtoK = nCr[povit.size()][needPivot] - nCr[povit.size() - leafRm.removedPivots.size()][needPivot];
                    RemovedKtoK = nCr[povit.size()][needPivot];
                }
                int needPP = int(needPivot) - 2;
                if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                    RemovedPtoP = nCr[povit.size() - 2][needPP];
                    PtoP = RemovedPtoP;
                    if (leafRm.removedPivots.size() + 2 <= povit.size()) {
                        PtoP -= nCr[povit.size() - 2 - leafRm.removedPivots.size()][needPP];
                    }
                }
                int needKP = int(needPivot) - 1;
                if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                    RemovedKtoP = nCr[povit.size() - 1][needKP];
                    KtoP = RemovedKtoP;
                    if (leafRm.removedPivots.size() + 1 <= povit.size()) {
                        KtoP -= nCr[povit.size() - 1 - leafRm.removedPivots.size()][needKP];
                    }
                }

                auto updateEdgeTemp = [&](daf::Size u, daf::Size v, double w, double rW) {
                    auto idx = edgeGraph.getEdgeIndex(u, v);
                    countingKE[idx] -= w;
                    if (edgeInHeap[idx]) {
                        heap.update(heapHandles[idx]);
                    }
                };
                // λ: 根据 (u,v) 的 pivot 状态与给定增量，统一更新 countingKE & 堆
                auto adjustKE = [&](const TreeGraphNode& u,
                                    const TreeGraphNode& v,
                                    daf::Size idx,
                                    double deltaKtoK,
                                    double deltaPtoP,
                                    double deltaKtoP) {
                    if (!u.isPivot && !v.isPivot) {                 // keep–keep
                        countingKE[idx] -= deltaKtoK;
                    } else if (u.isPivot && v.isPivot) {            // pivot–pivot
                        countingKE[idx] -= deltaPtoP;
                    } else {                                        // mixed
                        countingKE[idx] -= deltaKtoP;
                    }

                    if ((u.v == 6 && v.v == 7) ) {
                        std::cout << "asdfasdf: u: " << u << " v: " << v << " idx: " << idx
                                  << " deltaKtoK: " << deltaKtoK
                                  << " deltaPtoP: " << deltaPtoP
                                  << " deltaKtoP: " << deltaKtoP
                        << "leafId: " << leafId << " leafRm.removedPivots: " << leafRm.removedPivots
                                  << std::endl;

                    }
                };


                if (!leafRm.removedPivots.empty() && needPivot <= povit.size() - leafRm.removedPivots.size()) {
                    for (auto removedNbr: leafRm.removedPivots) {
                        treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                        // 1) keep–keep
                    }


                    const auto newLeaf = tree.removeNbrs(leafId, leafRm.removedPivots);
                    for (daf::Size i = 0; i < leafRm.removedPivots.size(); ++i) {
                        for (daf::Size j = i + 1; j < leafRm.removedPivots.size(); ++j) {
                            auto u = leafRm.removedPivots[i];
                            auto v = leafRm.removedPivots[j];
                            auto index = edgeGraph.getEdgeIndex(u, v);
                            countingKE[index] -= RemovedPtoP;
                            if (edgeInHeap[index]) {
                                heap.update(heapHandles[index]);
                            }
                        }
                    }

                    for (daf::Size i = 0; i < newLeaf.size(); ++i) {
                        for (daf::Size j = 0; j < leafRm.removedPivots.size(); ++j) {
                            auto u = newLeaf[i];
                            auto v = leafRm.removedPivots[j];
                            auto index = edgeGraph.getEdgeIndex(u.v, v);
                            if (u.isPivot) {            // pivot–pivot
                                countingKE[index] -= RemovedPtoP;
                            } else {                                        // mixed
                                countingKE[index] -= RemovedKtoP;
                            }
                            if (edgeInHeap[index]) {
                                heap.update(heapHandles[index]);
                            }
                        }
                    }


                    for (daf::Size i = 0; i < newLeaf.size(); ++i) {
                        for (daf::Size j = i + 1; j < newLeaf.size(); ++j) {
                            auto u = newLeaf[i];
                            auto v = newLeaf[j];
                            auto index = edgeGraph.getEdgeIndex(u.v, v.v);
                            adjustKE(u, v, index,
                                KtoK, PtoP, KtoP);
                            if (edgeInHeap[index]) {
                                heap.update(heapHandles[index]);
                            }
                        }
                    }
                } else {
                    for (auto removedNbr: leaf) {
                        treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                    }
                    for (daf::Size i = 0; i < leaf.size(); ++i) {
                        for (daf::Size j = i + 1; j < leaf.size(); ++j) {
                            auto u = leaf[i];
                            auto v = leaf[j];
                            auto index = edgeGraph.getEdgeIndex(u.v, v.v);
                            adjustKE(u, v, index,
                                RemovedKtoK, RemovedPtoP, RemovedKtoP);
                            if (edgeInHeap[index]) {
                                heap.update(heapHandles[index]);
                            }
                        }
                    }
                    tree.removeNode(leafId);
                }
            }


// #ifndef NDEBUG
//             std::cout << "countingKE: ";
//             PlusECDSet::printEdgeCore(edgeGraph, countingKE);
//                 // updateLeaf.print("updateLeaf");
//                 // removedLeaf.print("removedLeaf");
//             std::cout << "removedLeaf: " << removedLeaf << std::endl;
//
//             // //
//             // // std::cout << "degreeE: ";
//             // // PlusECDSet::printEdgeCore(edgeGraph, degreeE);
//             // std::cout << "degreeV: ";
//             // daf::printArray(degreeV, treeGraphV.adj_list.size());
//             //
//             std::cout << "coreE: ";
//             PlusECDSet::printEdgeCore(edgeGraph, coreE);
//             //
//
//             // removedLeaf.print("removedLeaf");
//
//             std::cout << "--------------------------------------------------------------------------------------------------------------------"<< std::endl;
// #endif

            leafRmInfo[leafId].clear();
            updateLeaf.clear();
            povit.clear();
            keepC.clear();
            // removedEdges.clear();
            // removedKeepC = false;
            // for (auto v: removedPovit) {
            //     isRemovedV[v] = false;
            // }
            // removedPovit.clear();
        }
        currentRemoveEdgeIds.clear();
        removedLeaf.clear();
        // currentRemoveLeafIds.clear();
#ifndef NDEBUG
        std::cout << "coreE: ";
        PlusECDSet::printEdgeCore(edgeGraph, coreE);
        std::cout << "countingKE: ";
        PlusECDSet::printEdgeCore(edgeGraph, countingKE);
        std::cout << "tree: ";
        tree.printGraphPerV();

        std::cout << "treeGraphV: ";
        treeGraphV.printGraphPerV();
#endif

    }

    for (auto i = 0;  i < edgeGraph.adj_list.size(); ++i) {
        auto counting = countingKE[i];
        if (counting != 0) {
            std::cerr << "Error: countingKE != 0" << std::endl;
            std::cerr << "countingKE: " << counting << std::endl;
            std::exit(1);
        }
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    // coreE
    // daf::printArray(coreE, edgeGraph.adj_list.size());


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
    newPivot.free();
    newKeepC.free();
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
