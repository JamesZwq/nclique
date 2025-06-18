//
// Created by 张文谦 on 25-3-4.
//

#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include <boost/heap/d_ary_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <set>

#include "../BK/BronKerboschRmEdge.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "debug/EdgeSet.h"
#include "graph/DynamicBipartiteGraph.hpp"
// #include "graph/DynamicGraph.h"
#include "dataStruct/disJoinSet.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];
// 放在你的函数外（比如文件顶部），保证编译时可见并内联
// #ifndef NDEBUG
// set NOEBUG as trus


namespace CDSet {
    struct LeafRmInfo {
        bool removedKeepC;
        daf::StaticVector<daf::Size> removedPivots{0};
        daf::StaticVector<std::pair<daf::Size, daf::Size> > removedEdges{0};

        LeafRmInfo() : removedKeepC(false) {
        }

        bool empty() const {
            return !removedKeepC && removedPivots.empty() && removedEdges.empty();
        }

        void init(daf::CliqueSize k) {
            removedKeepC = false;
            removedPivots.reserve(k);
            removedEdges.reserve(k * (k - 1) / 2); // k choose 2
        }


        void clear() {
            removedKeepC = false;
            removedPivots.clear();
            removedEdges.clear();
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


    struct CompareRClique {
        const double *RCliqueCounting; // 指向外部数组
        explicit CompareRClique(const double *coreLeaf) : RCliqueCounting(coreLeaf) {
        }

        // 注意：这里要返回 “a 排在前面” 的条件，为最小堆写成 coreLeaf[a] > coreLeaf[b]
        bool operator()(daf::Size const &a, daf::Size const &b) const {
            return RCliqueCounting[a] > RCliqueCounting[b];
        }
    };

    using DHeap = boost::heap::d_ary_heap<
        daf::Size,
        boost::heap::arity<8>,
        boost::heap::mutable_<true>,
        boost::heap::compare<CompareRClique>
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

    std::pair<std::vector<double>, std::vector<double> > countingPerEdgeAndRClique(
        const DynamicGraph<TreeGraphNode> &treeGraph,
        StaticCliqueIndex &cliqueHashmap,
        const Graph &edgeGraph,
        const daf::CliqueSize r,
        const daf::CliqueSize s) {
        // EdgeHashMap<double> coreE(edgeGraph.adj_list.size());
        // double *countingE = new double[edgeGraph.adj_list.size()];
        // memset(countingE, 0, edgeGraph.adj_list.size() * sizeof(double));
        std::vector<double> countingE(edgeGraph.adj_list.size(), 0.0);


        // double *rCliqueSCounting = new double[cliqueHashmap.size()];
        // memset(rCliqueSCounting, 0, cliqueHashmap.size() * sizeof(double));
        std::vector<double> rCliqueSCounting(cliqueHashmap.size() * 1.1, 0.0);
        daf::StaticVector<daf::Size> povit;
        daf::StaticVector<daf::Size> keepC;
        daf::Size count = 0;
        for (const auto &clique: treeGraph.adj_list) {
            povit.clear();
            keepC.clear();
            if (clique.size() < r) {
                continue;
            }
            for (auto &i: clique) {
                if (i.isPivot) {
                    povit.push_back(i.v);
                } else {
                    keepC.push_back(i.v);
                }
            }

            int needPivot = int(r) - int(keepC.size());

            // 1) keep-keep 边：两端都在 keepC 中
            double totalKcliques = -1;
            // 从 pivot.size() 个点里选 needPivot 个，再分配给每条 keep‑keep 边
            if (needPivot >= 0 && needPivot <= int(povit.size())) {
                totalKcliques = nCr[povit.size()][needPivot];
                for (size_t i = 0; i < keepC.size(); ++i) {
                    for (size_t j = i + 1; j < keepC.size(); ++j) {
                        daf::Size u = keepC[i], v = keepC[j];
                        auto index = edgeGraph.getEdgeCompressedId(u, v);
                        countingE[index] += totalKcliques;
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
                        countingE[index] += eachPivotKcliques;
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
                        countingE[index] += eachKeepPivotKcliques;
                    }
                }
            }
            // clique
            needPivot = int(s) - int(keepC.size());
            daf::enumerateCombinations(clique, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                daf::CliqueSize subNumKeepC = 0;
                daf::CliqueSize subNumPovit = 0;
                for (const auto &node: rClique) {
                    if (node.isPivot) {
                        subNumPovit++;
                    } else {
                        subNumKeepC++;
                    }
                }

                auto ncrValue = nCr[povit.size() - subNumPovit][needPivot - subNumPovit];
                // rCliqueSCounting[cliqueHashmap.byNewClique(rClique)] += ncrValue;
                auto [id, isNew] = cliqueHashmap.byNewClique(rClique);
                if (isNew) {
                    for (size_t i = 0; i < rClique.size(); ++i) {
                        for (size_t j = i + 1; j < rClique.size(); ++j) {
                            daf::Size u = rClique[i].v, v = rClique[j].v;
                            auto index = edgeGraph.getEdgeCompressedId(u, v);
                            countingE[index] += 1;
                        }
                    }
                    if (rCliqueSCounting.size() <= id) {
                        rCliqueSCounting.resize(id * 1.1, 0.0);
                    }
                }
                rCliqueSCounting[id] += ncrValue;


                return true;
            });
        }
        rCliqueSCounting.shrink_to_fit();
        povit.free();
        keepC.free();
        return {countingE, rCliqueSCounting};
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

    template<typename T>
    void printEdgeCore(const Graph &edgeGraph, const std::vector<T> coreE) {
        printEdgeCore(edgeGraph, coreE.data());
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
                    auto index = edgeGraph.getEdgeCompressedId(u.v, v.v);
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


std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > NucleusCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    auto time_start = std::chrono::high_resolution_clock::now();
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build",
                   [&]() {
                       cliqueIndex.build(tree, edgeGraph.adj_list.size());
                   });
    // tree.printGraphPerV();
    // cliqueIndex.print();
    // cliqueIndex.verify();


    auto [countingKE, countingRClique] = daf::timeCount("countingPerEdgeAndRClique",
                                                        [&]() {
                                                            return CDSet::countingPerEdgeAndRClique(
                                                                tree, cliqueIndex, edgeGraph, r, s);
                                                        });



    auto *coreE = new double[edgeGraph.adj_list.size()];
    memset(coreE, 0, edgeGraph.adj_list.size() * sizeof(double));
    // std::vector<double> leafCore = CDSet::initLeafCore(tree, countingKE, k, edgeGraph);
#ifndef NDEBUG
    tree.printGraphPerV();
    // daf::printArray(countingKE, edgeGraph.adj_list.size());
    CDSet::printEdgeCore(edgeGraph, countingKE);
    // CDSet::printEdgeCore(edgeGraph, degreeE);
#endif

    daf::StaticVector<daf::Size> povit;
    daf::StaticVector<daf::Size> keepC;
    daf::StaticVector<daf::Size> newPivot;
    daf::StaticVector<daf::Size> currentRemoveRcliqueIds(cliqueIndex.size());
    daf::StaticVector<daf::Size> currentRemoveEdgeIds(edgeGraph.adj_list.size());

    daf::StaticVector<bool> rCliqueInHeap(cliqueIndex.size());
    rCliqueInHeap.c_size = cliqueIndex.size();
    memset(rCliqueInHeap.getData(), true, cliqueIndex.size() * sizeof(bool));
    // daf::StaticVector<std::pair<daf::Size, double> > updateLeaf(tree.adj_list.size() * 10);
    std::map<daf::Size, double> updateLeaf;

    daf::StaticVector<daf::Size> removedLeaf(tree.adj_list.size());
    daf::StaticVector<CDSet::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();

    CDSet::DHeap heap{CDSet::CompareRClique(countingRClique.data())};
    std::vector<CDSet::DHeap::handle_type> heapHandles(cliqueIndex.size());

    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
        heapHandles[i] = heap.push(i);
    }
#ifndef NDEBUG
    std::cout << "coreE: ";
    CDSet::printEdgeCore(edgeGraph, coreE);
    std::cout << "countingKE: ";
    CDSet::printEdgeCore(edgeGraph, countingKE);

    std::cout << "countingRClique: " << countingRClique << std::endl;

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
        minCore = std::max(countingRClique[heap.top()], minCore);
        // 一次循环把所有 core==minCore 的 leaf 全部 pop 出来
        std::cout << "minCore: " << minCore << std::endl;
        while (!heap.empty() && countingRClique[heap.top()] <= minCore) {
            auto id = heap.top();
            rCliqueInHeap[id] = false;
            heap.pop();
            currentRemoveRcliqueIds.push_back(id);
            // std::cout << "progress: " << numProgress++ << "/" << edgeGraph.adj_list.size() << std::flush;
            // daf::printProgress(numProgress++, edgeGraph.adj_list.size());
        }

        for (auto rmRCliqueId: currentRemoveRcliqueIds) {
            auto rClique = cliqueIndex.byId(rmRCliqueId);
            // std::cout << "rClique: " << rClique << std::endl;
            for (int i = 0; i < rClique.size(); ++i) {
                for (int j = i + 1; j < rClique.size(); ++j) {
                    auto u = rClique[i];
                    auto v = rClique[j];
                    auto index = edgeGraph.getEdgeCompressedId(u, v);
                    countingKE[index] -= 1;
                    if (countingKE[index] == 0) {
                        currentRemoveEdgeIds.push_back(index);
                        coreE[index] = minCore;
                    }
                }
            }
        }

        // std::cout << "currentCore: " << currCore << std::endl;
#ifndef NDEBUG
        std::cout << "currentRemoveEdge: " << std::endl;
#endif
        for (auto edgeId: currentRemoveEdgeIds) {
            auto [edgeU, edgeV] = edgeGraph.getEdgeById(edgeId);
#ifndef NDEBUG
            std::cout << edgeU << " " << edgeV << std::endl;
#endif
            auto &adjEdgeU = treeGraphV.getNbr(edgeU);
            auto &adjEdgeV = treeGraphV.getNbr(edgeV);

            daf::intersect_dense_sets(adjEdgeU, adjEdgeV,
                                      [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                                          // std::cout << uClique << " " << vClique << std::endl;
                                          if (leafRmInfo[uClique.v].empty()) {
                                              removedLeaf.push_back(uClique.v);
                                              leafRmInfo[uClique.v].init(tree.adj_list[uClique.v].size());
                                          }
                                          if (leafRmInfo[uClique.v].removedKeepC) return;
                                          if (!uClique.isPivot && !vClique.isPivot) {
                                              leafRmInfo[uClique.v].removedKeepC = true;
                                          } else if (uClique.isPivot && vClique.isPivot) {
                                              leafRmInfo[uClique.v].removedEdges.push_back_with_check({edgeU, edgeV});
                                          } else if (uClique.isPivot && !vClique.isPivot) {
                                              leafRmInfo[uClique.v].removedPivots.push_back_with_check(edgeU);
                                          } else if (!uClique.isPivot && vClique.isPivot) {
                                              leafRmInfo[uClique.v].removedPivots.push_back_with_check(edgeV);
                                          }
                                      });
        }
        // removedLeaf.print("removedLeaf");
        // for (auto leafId : removedLeaf) {
        for (auto leafIdIdx = 0; leafIdIdx < removedLeaf.size(); ++leafIdIdx) {
            auto leafId = removedLeaf[leafIdIdx];
            auto leaf = tree.adj_list[leafId];

            // std::cout << "leafId: " << leafId << " leaf: " << leaf << std::endl;
            CDSet::LeafRmInfo &leafRm = leafRmInfo[leafId];
            if (leaf.empty()) {
                std::cout << "Error: newLeaf is empty" << std::endl;
                std::cout << "leafId: " << leafId << std::endl;
                std::cout << "leafRm.removedPivots: " << leafRm.removedPivots << std::endl;
                std::cout << "leafRm.removedEdges: " << leafRm.removedEdges << std::endl;
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


            daf::Size needPivot = s - keepC.size();
            // std::cout << "coreE: " << std::endl;
            // daf::printArray(coreE, edgeGraph.adj_list.size());
            if (leafRm.removedKeepC || needPivot > povit.size() - leafRm.removedPivots.size()) {
                daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &newLeaf) {
                    daf::CliqueSize subNumKeepC = 0;
                    daf::CliqueSize subNumPovit = 0;
                    for (const auto &node: newLeaf) {
                        if (node.isPivot) {
                            subNumPovit++;
                        } else {
                            subNumKeepC++;
                        }
                    }

                    auto ncrValue = nCr[povit.size() - subNumPovit][needPivot - subNumPovit];
                    auto cliqueIndexId = cliqueIndex.byClique(newLeaf);

                    countingRClique[cliqueIndexId] -= ncrValue;
                    if (rCliqueInHeap[cliqueIndexId]) {
                        heap.update(heapHandles[cliqueIndexId]);
                    }

                    return true;
                });

                for (auto i: leaf) {
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                }
                tree.removeNode(leafId);
            } else if (!leafRm.removedEdges.empty()) {
                // std::cout << "============================================================" << std::endl;

                auto initCore = [&](const std::vector<TreeGraphNode> &leaf, const daf::Size &leafId) {
                    daf::CliqueSize newPivotC = 0, newKeepC = 0;
                    for (auto i: leaf) {
                        if (i.isPivot) {
                            treeGraphV.addNbr(i.v, {leafId, true});
                            newPivotC++;
                        } else {
                            treeGraphV.addNbr(i.v, {leafId, false});
                            newKeepC++;
                        }
                    }
                    daf::Size needPivot = s - newKeepC;

                    daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &newLeaf) {
                        daf::CliqueSize subNumKeepC = 0;
                        daf::CliqueSize subNumPovit = 0;
                        for (const auto &node: newLeaf) {
                            if (node.isPivot) {
                                subNumPovit++;
                            } else {
                                subNumKeepC++;
                            }
                        }
                        auto ncrValue = nCr[newPivotC - subNumPovit][needPivot - subNumPovit];
                        auto cliqueIndexId = cliqueIndex.byClique(newLeaf);
                        countingRClique[cliqueIndexId] += ncrValue;

                        return true;
                    });
                };

                // if (!removedPovit.empty() && needPivot <= povit.size() - removedPovit.size())
                // daf::StaticVector<daf::Size> newLeafIds;
                auto &newLeaf = leaf;
                auto oldLeaf = leaf;
                for (auto leafV: leaf) {
                    if (leafV.isPivot) {
                        treeGraphV.removeNbr(leafV.v, {leafId, true});
                    } else {
                        treeGraphV.removeNbr(leafV.v, {leafId, false});
                    }
                }
                if (!leafRm.removedPivots.empty()) {
                    newLeaf = tree.removeNbrs(leafId, leafRm.removedPivots);
                }
                bkRmEdge::bronKerbosch(newLeaf, leafRm.removedEdges, s,
                                       [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                                           std::vector<TreeGraphNode> newSubLeaf = bkRmEdge::coverToVertex(
                                               c, pivots, leaf);
                                           auto newId = tree.addNode(newSubLeaf);
                                           initCore(tree.adj_list[newId], newId);
                                           if (newId >= leafRmInfo.size()) {
                                               removedLeaf.reserve(newId * 1.5);
                                               leafRmInfo.resize(newId * 1.5);
                                               heapHandles.resize(newId * 1.5);
                                           }
                                       }
                );

                daf::enumerateCombinations(oldLeaf, r, [&](const daf::StaticVector<TreeGraphNode> &newLeaf) {
                    daf::CliqueSize subNumKeepC = 0;
                    daf::CliqueSize subNumPovit = 0;
                    for (const auto &node: newLeaf) {
                        if (node.isPivot) {
                            subNumPovit++;
                        } else {
                            subNumKeepC++;
                        }
                    }
                    auto ncrValue = nCr[povit.size() - subNumPovit][needPivot - subNumPovit];
                    auto cliqueIndexId = cliqueIndex.byClique(newLeaf);
                    countingRClique[cliqueIndexId] -= ncrValue;
                    if (rCliqueInHeap[cliqueIndexId]) {
                        heap.update(heapHandles[cliqueIndexId]);
                    }

                    return true;
                });
                tree.removeNode(leafId);
            } else {
                // only povit removed


                if (!leafRm.removedPivots.empty() && needPivot <= povit.size() - leafRm.removedPivots.size()) {
                    for (auto removedNbr: leafRm.removedPivots) {
                        treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                        // 1) keep–keep
                    }

                    tree.removeNbrs(leafId, leafRm.removedPivots);
                    newPivot.clear();
                    for (const auto &node: tree.adj_list[leafId]) {
                        if (node.isPivot) {
                            newPivot.push_back(node.v);
                        }
                    }

                    // for removed r cliques
                    daf::enumAtLeastOneFromTwo(leafRm.removedPivots, tree.adj_list[leafId], r,
                                               [&](daf::Size *pivots, size_t pCount, TreeGraphNode *cliqueList,
                                                   daf::Size cliqueSize) {
                                                   daf::CliqueSize subNumKeepC = 0;
                                                   daf::CliqueSize subNumPovit = pCount;
                                                   for (daf::Size i = 0; i < cliqueSize; ++i) {
                                                       if (cliqueList[i].isPivot) {
                                                           subNumPovit++;
                                                       } else {
                                                           subNumKeepC++;
                                                       }
                                                   }

                                                   const auto cliqueId = cliqueIndex.byClique(
                                                       pivots, pCount, cliqueList, cliqueSize);
                                                   const auto ncrValue = nCr[povit.size() - subNumPovit][
                                                       needPivot - subNumPovit];
                                                   countingRClique[cliqueId] -= ncrValue;
                                                   if (rCliqueInHeap[cliqueId]) {
                                                       heap.update(heapHandles[cliqueId]);
                                                   }
                                                   return true;
                                               });
                    // for existing r cliques
                    daf::enumerateCombinations(keepC, newPivot, r,
                                               [&](const daf::StaticVector<daf::Size> &newKeep,
                                                   const daf::StaticVector<daf::Size> &newDrop) {
                                                   daf::CliqueSize subNumPovit = newDrop.size();
                                                   auto ncrValue =
                                                           nCr[povit.size() - subNumPovit][needPivot - subNumPovit] -
                                                           nCr[povit.size() - leafRm.removedPivots.size() - subNumPovit]
                                                           [needPivot - subNumPovit];

                                                   auto cliqueIndexId = cliqueIndex.byClique(newKeep, newDrop);
                                                   countingRClique[cliqueIndexId] -= ncrValue;
                                                   if (rCliqueInHeap[cliqueIndexId]) {
                                                       heap.update(heapHandles[cliqueIndexId]);
                                                   }
                                                   return true;
                                               });
                } else {
                    throw std::runtime_error(
                        "Error: leafRm.removedPivots.size() + needPivot > povit.size() - leafRm.removedPivots.size()");
                }
            }


#ifndef NDEBUG
            std::cout << "coreE: ";
            CDSet::printEdgeCore(edgeGraph, coreE);
            std::cout << "countingKE: ";
            CDSet::printEdgeCore(edgeGraph, countingKE);
            std::cout << "tree: ";
            tree.printGraphPerV();

            std::cout << "treeGraphV: ";
            treeGraphV.printGraphPerV();
#endif

            leafRmInfo[leafId].clear();
            updateLeaf.clear();
            povit.clear();
            keepC.clear();
        }
        currentRemoveEdgeIds.clear();
        currentRemoveRcliqueIds.clear();
        removedLeaf.clear();
        // currentRemoveLeafIds.clear();
#ifndef NDEBUG
        std::cout << "coreE: ";
        CDSet::printEdgeCore(edgeGraph, coreE);
        std::cout << "countingKE: ";
        CDSet::printEdgeCore(edgeGraph, countingKE);
        std::cout << "tree: ";
        tree.printGraphPerV();

        std::cout << "treeGraphV: ";
        treeGraphV.printGraphPerV();
#endif
    }

    for (auto i = 0; i < edgeGraph.adj_list.size(); ++i) {
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
    for (daf::Size i = 0; i < edgeGraph.adj_list.size(); i++) {
        // cover to int
        numCounting += countingKE[i];
    }

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

    if (numCounting != 0) {
        // exit 1
        std::cerr << "Error: numCounting != 0" << std::endl;
        std::cerr << "numCounting: " << numCounting << std::endl;
        throw std::runtime_error("numCounting != 0");
    }
    assert(numCounting == 0);

    // delete[] degreeV;
    delete[] coreE;
    povit.free();
    keepC.free();
    newPivot.free();
    // leafIds.free();
    return sortedK;
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
