//
// Created by _ on 25-3-4.
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
#include "graph/DynamicGraphSet.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];
// （），
// #ifndef NDEBUG
// set NOEBUG as trus


namespace PlusECDSet {
    struct LeafRmInfo {
        bool removedKeepC;
        daf::StaticVector<daf::Size> removedPivots{0};
        daf::StaticVector<std::pair<daf::Size, daf::Size> > removedEdges{0};

        LeafRmInfo() : removedKeepC(false) {
        }

        bool empty() const {
            return !removedKeepC && removedPivots.empty() && removedEdges.empty();
        }

        void init(auto capacity = 400) {
            removedKeepC = false;
            removedPivots.reserve(capacity);
            removedEdges.reserve(capacity);
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

    template<typename It1, typename It2, typename WeightT, typename UpdateFunc>
    inline void processEdgePairsImpl(It1 b1, It1 e1,
                                     It2 b2, It2 e2,
                                     WeightT weight,
                                     UpdateFunc &&upd) noexcept {
        if (weight < 0) return;

        // 
        // if all same, do nothing
        if (b1 == b2 && e1 == e2 && b1 == e1 && b2 == e2) {
            return;
        }
        if (b1 == b2 && e1 == e2) {
            // ：i < j
            for (auto it = b1; it + 1 != e1; ++it) {
                auto u = *it;
                for (auto jt = it + 1; jt != e1; ++jt) {
                    upd(u, *jt, weight);
                }
            }
        } else {
            // ：
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
        typename WeightT, typename UpdateFunc
    >
    inline void processEdgePairs(const Range1 &r1,
                                 const Range2 &r2,
                                 WeightT weight,
                                 UpdateFunc &&upd) noexcept {
        processEdgePairsImpl(
            std::begin(r1), std::end(r1),
            std::begin(r2), std::end(r2),
            weight,
            std::forward<UpdateFunc>(upd)
        );
    }

    template<
        typename Range, typename WeightT,
        typename UpdateFunc
    >
    inline void processEdgePairs(const Range &r,
                                 WeightT weight,
                                 UpdateFunc &&upd) noexcept {
        processEdgePairsImpl(
            std::begin(r), std::end(r),
            std::begin(r), std::end(r),
            weight,
            std::forward<UpdateFunc>(upd)
        );
    }


    struct CompareEdge {
        const double *edgeCounting; // 
        explicit CompareEdge(const double *coreLeaf) : edgeCounting(coreLeaf) {
        }

        // ： “a ” ， coreLeaf[a] > coreLeaf[b]
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
            const int needPivot = k - keepC.size(); //  pivot 
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
            //  pivot 
            int needPivot = int(k) - int(keepC.size());
            // 1) keep-keep ： keepC 
            double totalKcliques = -1;
            //  pivot.size()  needPivot ， keep‑keep 
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

            // 2) pivot‑pivot ： povit 
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


            // 3) cross ： keepC， povit
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

    std::pair<double *, daf::Size *> countingPerEdge(const DynamicGraph<TreeGraphNode> &treeGraph,
                                                     const Graph &edgeGraph,
                                                     const daf::CliqueSize k) {
        double *countingE = new double[edgeGraph.adj_list.size()];
        daf::Size *degreeE = new daf::Size[edgeGraph.adj_list.size()];
        memset(countingE, 0, edgeGraph.adj_list.size() * sizeof(double));
        memset(degreeE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));

        const int numLeaves = (int)treeGraph.adj_list.size();
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            daf::StaticVector<daf::Size> tPovit, tKeepC;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 256)
#endif
            for (int li = 0; li < numLeaves; ++li) {
                const auto &clique = treeGraph.adj_list[li];
                if (clique.size() < k) continue;

                tPovit.clear(); tKeepC.clear();
                for (const auto &node : clique) {
                    if (node.isPivot) tPovit.push_back(node.v);
                    else tKeepC.push_back(node.v);
                }

                int needPivot = int(k) - int(tKeepC.size());

                // 1) keep-keep
                if (needPivot >= 0 && needPivot <= int(tPovit.size())) {
                    double totalKcliques = nCr[tPovit.size()][needPivot];
                    for (size_t i = 0; i < tKeepC.size(); ++i)
                        for (size_t j = i + 1; j < tKeepC.size(); ++j) {
                            auto index = edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]);
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            countingE[index] += totalKcliques;
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            degreeE[index]++;
                        }
                }

                // 2) pivot-pivot
                int needPP = needPivot - 2;
                if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
                    double eachPP = nCr[tPovit.size() - 2][needPP];
                    for (size_t i = 0; i < tPovit.size(); ++i)
                        for (size_t j = i + 1; j < tPovit.size(); ++j) {
                            auto index = edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]);
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            countingE[index] += eachPP;
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            degreeE[index]++;
                        }
                }

                // 3) keep-pivot cross
                int needKP = needPivot - 1;
                if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
                    double eachKP = nCr[tPovit.size() - 1][needKP];
                    for (size_t i = 0; i < tKeepC.size(); ++i)
                        for (size_t j = 0; j < tPovit.size(); ++j) {
                            auto index = edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]);
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            countingE[index] += eachKP;
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            degreeE[index]++;
                        }
                }
            }
            tPovit.free(); tKeepC.free();
        }
        return {countingE, degreeE};
    }

    std::pair<double *, daf::Size *> countingPerEdge(const DynamicGraphSet<TreeGraphNode> &treeGraph,
                                                     const Graph &edgeGraph,
                                                     const daf::CliqueSize k) {
        double *countingE = new double[edgeGraph.adj_list.size()];
        daf::Size *degreeE = new daf::Size[edgeGraph.adj_list.size()];
        memset(countingE, 0, edgeGraph.adj_list.size() * sizeof(double));
        memset(degreeE, 0, edgeGraph.adj_list.size() * sizeof(daf::Size));

        const int numLeaves = (int)treeGraph.adj_list.size();
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            daf::StaticVector<daf::Size> tPovit, tKeepC;
#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 256)
#endif
            for (int li = 0; li < numLeaves; ++li) {
                const auto &clique = treeGraph.adj_list[li];
                if (clique.size() < k) continue;

                tPovit.clear(); tKeepC.clear();
                for (const auto &node : clique) {
                    if (node.isPivot) tPovit.push_back(node.v);
                    else tKeepC.push_back(node.v);
                }

                int needPivot = int(k) - int(tKeepC.size());

                // 1) keep-keep
                if (needPivot >= 0 && needPivot <= int(tPovit.size())) {
                    double totalKcliques = nCr[tPovit.size()][needPivot];
                    for (size_t i = 0; i < tKeepC.size(); ++i)
                        for (size_t j = i + 1; j < tKeepC.size(); ++j) {
                            auto index = edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]);
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            countingE[index] += totalKcliques;
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            degreeE[index]++;
                        }
                }

                // 2) pivot-pivot
                int needPP = needPivot - 2;
                if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
                    double eachPP = nCr[tPovit.size() - 2][needPP];
                    for (size_t i = 0; i < tPovit.size(); ++i)
                        for (size_t j = i + 1; j < tPovit.size(); ++j) {
                            auto index = edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]);
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            countingE[index] += eachPP;
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            degreeE[index]++;
                        }
                }

                // 3) keep-pivot cross
                int needKP = needPivot - 1;
                if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
                    double eachKP = nCr[tPovit.size() - 1][needKP];
                    for (size_t i = 0; i < tKeepC.size(); ++i)
                        for (size_t j = 0; j < tPovit.size(); ++j) {
                            auto index = edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]);
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            countingE[index] += eachKP;
#ifdef _OPENMP
                            #pragma omp atomic
#endif
                            degreeE[index]++;
                        }
                }
            }
            tPovit.free(); tKeepC.free();
        }
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

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {
    auto time_start = std::chrono::high_resolution_clock::now();
    auto [countingKE_dbl, degreeERemove] = PlusECDSet::countingPerEdge(tree, edgeGraph, k);

    // Convert double support counts to exact integer arithmetic to eliminate float drift
    const daf::Size numEdgesInit = edgeGraph.adj_list.size();
    auto *countingKE = new long long[numEdgesInit];
    for (daf::Size i = 0; i < numEdgesInit; ++i)
        countingKE[i] = std::llround(countingKE_dbl[i]);
    delete[] countingKE_dbl;

#ifndef NDEBUG
    tree.printGraphPerV();
    // PlusECDSet::printEdgeCore(edgeGraph, countingKE);
#endif

    // memset(removedLeaf.getData(), false, tree.adj_list.size() * sizeof(bool));

    auto *coreE = new long long[numEdgesInit];
    memset(coreE, 0, numEdgesInit * sizeof(long long));

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
    // std::map<daf::Size, double> updateLeaf;

    daf::StaticVector<daf::Size> removedLeaf(tree.adj_list.size());
    daf::StaticVector<PlusECDSet::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();


    long long currCore = 0;

    // --- Bucket array replacing heap ---
    const daf::Size numEdges = edgeGraph.adj_list.size();
    int maxBucket = 0;
    for (daf::Size i = 0; i < numEdges; ++i) {
        if (countingKE[i] > 0)
            maxBucket = std::max(maxBucket, (int)countingKE[i]);
    }
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(numEdges);
    std::vector<daf::Size> pos_in_bucket(numEdges);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numEdges; ++i) {
        if (countingKE[i] == 0) {
            edgeInHeap[i] = false;
            continue;
        }
        int b = (int)countingKE[i];
        bucket_of[i] = b;
        pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i);
        remainingInHeap++;
    }
    int curBucket = 0;

    // Bucket move helper
    auto bucketMove = [&](daf::Size id) {
        int newB = std::max(0, (int)countingKE[id]);
        int oldB = bucket_of[id];
        if (newB == oldB) return;
        auto& oldVec = buckets[oldB];
        daf::Size myPos = pos_in_bucket[id];
        if (myPos < oldVec.size() - 1) {
            daf::Size last = oldVec.back();
            oldVec[myPos] = last;
            pos_in_bucket[last] = myPos;
        }
        oldVec.pop_back();
        if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
        bucket_of[id] = newB;
        pos_in_bucket[id] = buckets[newB].size();
        buckets[newB].push_back(id);
        if (newB < curBucket) curBucket = newB;
    };

    // Striped vertex locks for parallel Phase 2b
    static constexpr int NUM_VERTEX_LOCKS = 1024; // must be power of 2
#ifdef _OPENMP
    omp_lock_t vertexLocks[NUM_VERTEX_LOCKS];
    for (int i = 0; i < NUM_VERTEX_LOCKS; i++) omp_init_lock(&vertexLocks[i]);
#endif

    // Deferred bucket move tracking
    std::vector<bool> dirtyMark(numEdges, false);
    std::vector<daf::Size> dirtyEdges;
    dirtyEdges.reserve(8192);
    auto deferMove = [&](daf::Size idx) {
        if (edgeInHeap[idx] && !dirtyMark[idx]) {
            dirtyMark[idx] = true;
            dirtyEdges.push_back(idx);
        }
    };

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
    long long minCore = 0;

    long long cntA = 0, cntB = 0, cntC = 0, totalIters = 0;

    // Per-leaf striped locks for lock-free merge in Phase 1
    static constexpr int NUM_LEAF_LOCKS = 4096; // must be power of 2
#ifdef _OPENMP
    omp_lock_t leafLocks[NUM_LEAF_LOCKS];
    for (int i = 0; i < NUM_LEAF_LOCKS; i++) omp_init_lock(&leafLocks[i]);
    const int numThreadsE = omp_get_max_threads();
#else
    const int numThreadsE = 1;
#endif
    // Track which leaves are newly affected (per-thread, for removedLeaf collection)
    std::vector<std::vector<daf::Size>> threadNewLeaves(numThreadsE);
    for (auto& v : threadNewLeaves) v.reserve(4096);
    // Per-leaf "is affected this iteration" flag
    std::vector<bool> leafAffected(leafRmInfo.size(), false);
    // Pre-allocated for Phase 2c (reused across iterations)
    std::vector<daf::Size> caseBLeafIds;
    caseBLeafIds.reserve(1024);

    while (remainingInHeap > 0) {
        // --- Bucket pop ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((long long)curBucket, minCore);

        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                edgeInHeap[id] = false;
                currentRemoveEdgeIds.push_back(id);
                coreE[id] = minCore;
                remainingInHeap--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (curBucket + 1) <= (int)minCore) {
                curBucket++;
            } else break;
        }

        currCore = minCore;

        if (remainingInHeap == 0) break;

        totalIters++;

        // --- Phase 1 (parallel): intersect edge adjacency lists with direct-write ---
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 8) if(currentRemoveEdgeIds.size() > 16)
#endif
        for (int ei = 0; ei < (int)currentRemoveEdgeIds.size(); ++ei) {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            auto edgeId = currentRemoveEdgeIds[ei];
            auto [edgeU, edgeV] = edgeGraph.getEdgeById(edgeId);
            auto &adjEdgeU = treeGraphV.getNbr(edgeU);
            auto &adjEdgeV = treeGraphV.getNbr(edgeV);
            daf::intersect_dense_sets(adjEdgeU, adjEdgeV,
                                      [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                                          daf::Size leafId = uClique.v;
                                          auto lockIdx = leafId & (NUM_LEAF_LOCKS - 1);
#ifdef _OPENMP
                                          omp_set_lock(&leafLocks[lockIdx]);
#endif
                                          auto &lrm = leafRmInfo[leafId];
                                          bool wasEmpty = lrm.empty();
                                          if (!lrm.removedKeepC) {
                                              if (!uClique.isPivot && !vClique.isPivot) {
                                                  lrm.removedKeepC = true;
                                              } else if (uClique.isPivot && vClique.isPivot) {
                                                  lrm.removedEdges.push_back({edgeU, edgeV});
                                              } else if (uClique.isPivot) {
                                                  lrm.removedPivots.push_back(edgeU);
                                              } else {
                                                  lrm.removedPivots.push_back(edgeV);
                                              }
                                          }
#ifdef _OPENMP
                                          omp_unset_lock(&leafLocks[lockIdx]);
#endif
                                          if (wasEmpty) threadNewLeaves[tid].push_back(leafId);
                                      });
        }
        // Fast merge: only collect leaf IDs (updates already applied)
        for (auto& newLeaves : threadNewLeaves) {
            for (auto leafId : newLeaves) {
                if (!leafAffected[leafId]) {
                    leafAffected[leafId] = true;
                    removedLeaf.push_back(leafId);
                }
            }
            newLeaves.clear();
        }

        // Pre-sort removedPivots for all leaves (can parallelize)
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 64) if(removedLeaf.size() > 64)
#endif
        for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
            auto leafId = removedLeaf[leafIdIdx];
            std::ranges::sort(leafRmInfo[leafId].removedPivots);
            leafRmInfo[leafId].removedPivots.unique();
        }


        // ========== Phase 2a: PARALLEL delta computation for Case A & C ==========
        // Each leaf's delta is independent — countingKE[] updates use omp atomic.
        // Tree mutations are deferred to Phase 2b (serial).
#ifdef _OPENMP
        #pragma omp parallel if(removedLeaf.size() > 32)
#endif
        {
            daf::StaticVector<daf::Size> tPovit, tKeepC;
#ifdef _OPENMP
            int tid2a = omp_get_thread_num();
#else
            int tid2a = 0;
#endif
            auto &myDirty = threadNewLeaves[tid2a]; // reuse threadNewLeaves as dirty list

            auto atomicSub = [&](daf::Size idx, long long w) {
#ifdef _OPENMP
                #pragma omp atomic
#endif
                countingKE[idx] -= w;
                if (edgeInHeap[idx]) myDirty.push_back(idx);
            };

#ifdef _OPENMP
            #pragma omp for schedule(dynamic, 64)
#endif
            for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                auto leafId = removedLeaf[leafIdIdx];
                PlusECDSet::LeafRmInfo &leafRm = leafRmInfo[leafId];

                auto leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;

                tPovit.clear(); tKeepC.clear();
                for (auto node : leaf) {
                    if (node.isPivot) tPovit.push_back(node.v);
                    else tKeepC.push_back(node.v);
                }
                daf::Size needPivot = k - tKeepC.size();

                // Determine true case
                bool isDeadLeaf = leafRm.removedKeepC || needPivot > tPovit.size() - leafRm.removedPivots.size();
                bool isCaseB = !leafRm.removedEdges.empty() && !isDeadLeaf;
                if (isCaseB) continue; // Case B handled in Phase 2c

                if (isDeadLeaf) {
                    // ---- Case A: leaf dies — subtract full contribution ----
                    long long KtoK = 0, KtoP = 0, PtoP = 0;
                    if (needPivot <= tPovit.size()) {
                        KtoK = std::llround(nCr[tPovit.size()][needPivot]);
                        for (daf::Size i = 0; i + 1 < tKeepC.size(); ++i)
                            for (daf::Size j = i + 1; j < tKeepC.size(); ++j)
                                atomicSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]), KtoK);
                    }
                    int needPP = int(needPivot) - 2;
                    if (0 <= needPP && needPP <= int(tPovit.size()) - 2) {
                        PtoP = std::llround(nCr[tPovit.size() - 2][needPP]);
                        for (daf::Size i = 0; i + 1 < tPovit.size(); ++i)
                            for (daf::Size j = i + 1; j < tPovit.size(); ++j)
                                atomicSub(edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]), PtoP);
                    }
                    int needKP = int(needPivot) - 1;
                    if (0 <= needKP && needKP <= int(tPovit.size()) - 1) {
                        KtoP = std::llround(nCr[tPovit.size() - 1][needKP]);
                        for (daf::Size i = 0; i < tKeepC.size(); ++i)
                            for (daf::Size j = 0; j < tPovit.size(); ++j)
                                atomicSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]), KtoP);
                    }
                } else {
                    // ---- Case C: only pivots removed ----
                    long long KtoK = 0, KtoP = 0, PtoP = 0;
                    long long RemovedKtoK = 0, RemovedKtoP = 0, RemovedPtoP = 0;

                    if (needPivot <= tPovit.size()) {
                        KtoK = std::llround(nCr[tPovit.size()][needPivot]) - std::llround(nCr[tPovit.size() - leafRm.removedPivots.size()][needPivot]);
                        RemovedKtoK = std::llround(nCr[tPovit.size()][needPivot]);
                    }
                    int needPP = int(needPivot) - 2;
                    if (0 <= needPP && needPP <= int(tPovit.size()) - 2) {
                        RemovedPtoP = std::llround(nCr[tPovit.size() - 2][needPP]);
                        PtoP = RemovedPtoP;
                        if (leafRm.removedPivots.size() + 2 <= tPovit.size())
                            PtoP -= std::llround(nCr[tPovit.size() - 2 - leafRm.removedPivots.size()][needPP]);
                    }
                    int needKP = int(needPivot) - 1;
                    if (0 <= needKP && needKP <= int(tPovit.size()) - 1) {
                        RemovedKtoP = std::llround(nCr[tPovit.size() - 1][needKP]);
                        KtoP = RemovedKtoP;
                        if (leafRm.removedPivots.size() + 1 <= tPovit.size())
                            KtoP -= std::llround(nCr[tPovit.size() - 1 - leafRm.removedPivots.size()][needKP]);
                    }

                    if (!leafRm.removedPivots.empty() && needPivot <= tPovit.size() - leafRm.removedPivots.size()) {
                        // Build filtered newLeaf (read-only on tree)
                        daf::StaticVector<TreeGraphNode> newLeafF;
                        for (auto node : leaf) {
                            if (!std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), node.v))
                                newLeafF.push_back(node);
                        }
                        // removed×removed
                        for (daf::Size i = 0; i + 1 < leafRm.removedPivots.size(); ++i)
                            for (daf::Size j = i + 1; j < leafRm.removedPivots.size(); ++j)
                                atomicSub(edgeGraph.getEdgeCompressedId(leafRm.removedPivots[i], leafRm.removedPivots[j]), RemovedPtoP);
                        // newLeaf×removed
                        for (daf::Size i = 0; i < newLeafF.size(); ++i)
                            for (daf::Size j = 0; j < leafRm.removedPivots.size(); ++j) {
                                long long d = newLeafF[i].isPivot ? RemovedPtoP : RemovedKtoP;
                                atomicSub(edgeGraph.getEdgeCompressedId(newLeafF[i].v, leafRm.removedPivots[j]), d);
                            }
                        // newLeaf×newLeaf
                        for (daf::Size i = 0; i + 1 < newLeafF.size(); ++i)
                            for (daf::Size j = i + 1; j < newLeafF.size(); ++j) {
                                auto &u = newLeafF[i], &v = newLeafF[j];
                                long long d = (!u.isPivot && !v.isPivot) ? KtoK : (u.isPivot && v.isPivot) ? PtoP : KtoP;
                                atomicSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                            }
                        newLeafF.free();
                    } else {
                        // Full removal (leaf dies)
                        for (daf::Size i = 0; i + 1 < leaf.size(); ++i)
                            for (daf::Size j = i + 1; j < leaf.size(); ++j) {
                                auto &u = leaf[i], &v = leaf[j];
                                long long d = (!u.isPivot && !v.isPivot) ? RemovedKtoK : (u.isPivot && v.isPivot) ? RemovedPtoP : RemovedKtoP;
                                atomicSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                            }
                    }
                }
            } // end omp for
            tPovit.free(); tKeepC.free();
        } // end omp parallel

        // Merge per-thread dirty edges from Phase 2a into dirtyEdges
        for (int t = 0; t < numThreadsE; ++t) {
            for (auto idx : threadNewLeaves[t]) {
                if (!dirtyMark[idx]) {
                    dirtyMark[idx] = true;
                    dirtyEdges.push_back(idx);
                }
            }
            threadNewLeaves[t].clear();
        }

        // ========== Phase 2b: PARALLEL tree mutations for Case A & C ==========
        {
            long long localCntA = 0, localCntC = 0;
#ifdef _OPENMP
            #pragma omp parallel reduction(+:localCntA,localCntC) if(removedLeaf.size() > 32)
#endif
            {
                std::vector<daf::Size> myDeadNodes;
#ifdef _OPENMP
                #pragma omp for schedule(dynamic, 64)
#endif
                for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                    auto leafId = removedLeaf[leafIdIdx];
                    PlusECDSet::LeafRmInfo &leafRm = leafRmInfo[leafId];
                    auto leaf = tree.adj_list[leafId];
                    if (leaf.empty()) continue;

                    daf::Size numKeeps = 0;
                    for (auto node : leaf) if (!node.isPivot) numKeeps++;
                    daf::Size needPivot = k - numKeeps;
                    daf::Size numPivots = leaf.size() - numKeeps;

                    bool isDeadLeaf = leafRm.removedKeepC || needPivot > numPivots - leafRm.removedPivots.size();
                    bool isCaseB = !leafRm.removedEdges.empty() && !isDeadLeaf;
                    if (isCaseB) continue;

                    if (isDeadLeaf) {
                        localCntA++;
                        for (auto i : leaf) {
                            auto lockIdx = i.v & (NUM_VERTEX_LOCKS - 1);
#ifdef _OPENMP
                            omp_set_lock(&vertexLocks[lockIdx]);
#endif
                            treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
#ifdef _OPENMP
                            omp_unset_lock(&vertexLocks[lockIdx]);
#endif
                        }
                        tree.adj_list[leafId].clear();
                        myDeadNodes.push_back(leafId);
                    } else {
                        localCntC++;
                        if (!leafRm.removedPivots.empty() && needPivot <= numPivots - leafRm.removedPivots.size()) {
                            for (auto removedNbr : leafRm.removedPivots) {
                                auto lockIdx = removedNbr & (NUM_VERTEX_LOCKS - 1);
#ifdef _OPENMP
                                omp_set_lock(&vertexLocks[lockIdx]);
#endif
                                treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
#ifdef _OPENMP
                                omp_unset_lock(&vertexLocks[lockIdx]);
#endif
                            }
                            tree.removeNbrs(leafId, leafRm.removedPivots);
                        } else {
                            for (auto i : leaf) {
                                auto lockIdx = i.v & (NUM_VERTEX_LOCKS - 1);
#ifdef _OPENMP
                                omp_set_lock(&vertexLocks[lockIdx]);
#endif
                                treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
#ifdef _OPENMP
                                omp_unset_lock(&vertexLocks[lockIdx]);
#endif
                            }
                            tree.adj_list[leafId].clear();
                            myDeadNodes.push_back(leafId);
                        }
                    }
                    leafRmInfo[leafId].clear();
                }
                // Merge dead nodes into tree's free list
#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    for (auto id : myDeadNodes)
                        tree.recycleNode(id);
                }
            }
            cntA += localCntA;
            cntC += localCntC;
        }


        // ========== Phase 2c: BK + apply for Case B ==========
        caseBLeafIds.clear();
        for (daf::Size leafIdIdx = 0; leafIdIdx < removedLeaf.size(); ++leafIdIdx) {
            auto leafId = removedLeaf[leafIdIdx];
            PlusECDSet::LeafRmInfo &leafRm = leafRmInfo[leafId];
            if (leafRm.removedEdges.empty()) continue;
            auto &leaf = tree.adj_list[leafId];
            if (leaf.empty()) continue;
            daf::Size numKeepsB = 0;
            for (auto node : leaf) if (!node.isPivot) numKeepsB++;
            daf::Size needPivotB = k - numKeepsB;
            daf::Size numPivotsB = leaf.size() - numKeepsB;
            if (leafRm.removedKeepC || needPivotB > numPivotsB - leafRm.removedPivots.size()) continue;
            caseBLeafIds.push_back(leafId);
        }
        cntB += caseBLeafIds.size();

        if (numThreadsE > 1 && caseBLeafIds.size() > 4) {
            // Parallel BK computation
            std::vector<std::vector<std::vector<TreeGraphNode>>> bkResults(caseBLeafIds.size());
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 1)
#endif
            for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
                auto leafId = caseBLeafIds[bi];
                auto &leafRm = leafRmInfo[leafId];
                auto leafCopy = tree.adj_list[leafId];
                if (!leafRm.removedPivots.empty()) {
                    std::vector<TreeGraphNode> filtered;
                    filtered.reserve(leafCopy.size());
                    for (auto &node : leafCopy)
                        if (!std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), node.v))
                            filtered.push_back(node);
                    leafCopy = std::move(filtered);
                }
                bkRmEdge::bronKerbosch(leafCopy, leafRm.removedEdges, k,
                    [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                        bkResults[bi].push_back(bkRmEdge::coverToVertex(c, pivots, leafCopy));
                    });
            }

            // Serial apply — old contribution already subtracted in Phase 2a
            for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
                auto leafId = caseBLeafIds[bi];
                PlusECDSet::LeafRmInfo &leafRm = leafRmInfo[leafId];
                auto leaf = tree.adj_list[leafId];

                povit.clear(); keepC.clear();
                for (auto node : leaf) {
                    if (node.isPivot) povit.push_back(node.v);
                    else keepC.push_back(node.v);
                }
                daf::Size needPivot = k - keepC.size();
                long long KtoK = 0, KtoP = 0, PtoP = 0;

                auto addW = [&](daf::Size u, daf::Size v, long long w) {
                    auto idx = edgeGraph.getEdgeCompressedId(u, v);
                    countingKE[idx] += w;
                };

                // Add new sub-leaves + their contributions
                for (auto &subLeaf : bkResults[bi]) {
                    auto newId = tree.addNode(subLeaf);
                    newPivot.clear(); newKeepC.clear();
                    for (auto i : tree.adj_list[newId]) {
                        if (i.isPivot) { newPivot.push_back(i.v); treeGraphV.addNbr(i.v, {newId, true}); }
                        else { newKeepC.push_back(i.v); treeGraphV.addNbr(i.v, {newId, false}); }
                    }
                    daf::Size np = k - newKeepC.size();
                    long long KtoK = 0, KtoP = 0, PtoP = 0;
                    if (np <= newPivot.size() && newKeepC.size() > 1) {
                        KtoK = std::llround(nCr[newPivot.size()][np]);
                        PlusECDSet::processEdgePairs(newKeepC, KtoK, addW);
                    }
                    int nPP = int(np) - 2;
                    if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                        PtoP = std::llround(nCr[newPivot.size() - 2][nPP]);
                        PlusECDSet::processEdgePairs(newPivot, PtoP, addW);
                    }
                    int nKP = int(np) - 1;
                    if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                        KtoP = std::llround(nCr[newPivot.size() - 1][nKP]);
                        PlusECDSet::processEdgePairs(newKeepC, newPivot, KtoP, addW);
                    }
                    newPivot.clear(); newKeepC.clear();
                    if (newId >= leafRmInfo.size()) {
                        removedLeaf.reserve(newId * 1.5);
                        leafRmInfo.resize(newId * 1.5);
                        leafAffected.resize(newId * 1.5, false);
                    }
                }

                // Remove old leaf from treeGraphV
                for (auto leafV : leaf) {
                    if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                    else treeGraphV.removeNbr(leafV.v, {leafId, false});
                }
                if (!leafRm.removedPivots.empty()) {
                    tree.removeNbrs(leafId, leafRm.removedPivots);
                }

                // Remove old contribution
                auto removeW = [&](daf::Size u, daf::Size v, long long w) {
                    auto idx = edgeGraph.getEdgeCompressedId(u, v);
                    countingKE[idx] -= w;
                    deferMove(idx);
                };
                if (needPivot <= povit.size()) {
                    KtoK = std::llround(nCr[povit.size()][needPivot]);
                    PlusECDSet::processEdgePairs(keepC, KtoK, removeW);
                }
                int needPP = int(needPivot) - 2;
                if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                    PtoP = std::llround(nCr[povit.size() - 2][needPP]);
                    PlusECDSet::processEdgePairs(povit, PtoP, removeW);
                }
                int needKP = int(needPivot) - 1;
                if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                    KtoP = std::llround(nCr[povit.size() - 1][needKP]);
                    PlusECDSet::processEdgePairs(keepC, povit, KtoP, removeW);
                }

                tree.removeNode(leafId);
                leafRmInfo[leafId].clear();
                povit.clear(); keepC.clear();
            }
        } else {
            // Serial BK path (no bkResults allocation, BK + apply inline)
            for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
                auto leafId = caseBLeafIds[bi];
                PlusECDSet::LeafRmInfo &leafRm = leafRmInfo[leafId];
                auto leaf = tree.adj_list[leafId];

                povit.clear(); keepC.clear();
                for (auto node : leaf) {
                    if (node.isPivot) povit.push_back(node.v);
                    else keepC.push_back(node.v);
                }
                daf::Size needPivot = k - keepC.size();
                long long KtoK = 0, KtoP = 0, PtoP = 0;

                auto addW = [&](daf::Size u, daf::Size v, long long w) {
                    auto idx = edgeGraph.getEdgeCompressedId(u, v);
                    countingKE[idx] += w;
                };

                // Remove old leaf from treeGraphV
                for (auto leafV : leaf) {
                    if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                    else treeGraphV.removeNbr(leafV.v, {leafId, false});
                }
                if (!leafRm.removedPivots.empty()) {
                    tree.removeNbrs(leafId, leafRm.removedPivots);
                }

                // Inline BK + apply new sub-leaves
                auto &leafRef = tree.adj_list[leafId];
                bkRmEdge::bronKerbosch(leafRef, leafRm.removedEdges, k,
                    [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                        auto subLeaf = bkRmEdge::coverToVertex(c, pivots, leafRef);
                        auto newId = tree.addNode(subLeaf);
                        newPivot.clear(); newKeepC.clear();
                        for (auto i : tree.adj_list[newId]) {
                            if (i.isPivot) { newPivot.push_back(i.v); treeGraphV.addNbr(i.v, {newId, true}); }
                            else { newKeepC.push_back(i.v); treeGraphV.addNbr(i.v, {newId, false}); }
                        }
                        daf::Size np = k - newKeepC.size();
                        long long KtoK = 0, KtoP = 0, PtoP = 0;
                        if (np <= newPivot.size() && newKeepC.size() > 1) {
                            KtoK = std::llround(nCr[newPivot.size()][np]);
                            PlusECDSet::processEdgePairs(newKeepC, KtoK, addW);
                        }
                        int nPP = int(np) - 2;
                        if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                            PtoP = std::llround(nCr[newPivot.size() - 2][nPP]);
                            PlusECDSet::processEdgePairs(newPivot, PtoP, addW);
                        }
                        int nKP = int(np) - 1;
                        if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                            KtoP = std::llround(nCr[newPivot.size() - 1][nKP]);
                            PlusECDSet::processEdgePairs(newKeepC, newPivot, KtoP, addW);
                        }
                        newPivot.clear(); newKeepC.clear();
                        if (newId >= leafRmInfo.size()) {
                            removedLeaf.reserve(newId * 1.5);
                            leafRmInfo.resize(newId * 1.5);
                            leafAffected.resize(newId * 1.5, false);
                        }
                    });

                // Remove old contribution
                auto removeW = [&](daf::Size u, daf::Size v, long long w) {
                    auto idx = edgeGraph.getEdgeCompressedId(u, v);
                    countingKE[idx] -= w;
                    deferMove(idx);
                };
                if (needPivot <= povit.size()) {
                    KtoK = std::llround(nCr[povit.size()][needPivot]);
                    PlusECDSet::processEdgePairs(keepC, KtoK, removeW);
                }
                int needPP = int(needPivot) - 2;
                if (0 <= needPP && needPP <= int(povit.size()) - 2) {
                    PtoP = std::llround(nCr[povit.size() - 2][needPP]);
                    PlusECDSet::processEdgePairs(povit, PtoP, removeW);
                }
                int needKP = int(needPivot) - 1;
                if (0 <= needKP && needKP <= int(povit.size()) - 1) {
                    KtoP = std::llround(nCr[povit.size() - 1][needKP]);
                    PlusECDSet::processEdgePairs(keepC, povit, KtoP, removeW);
                }

                tree.removeNode(leafId);
                leafRmInfo[leafId].clear();
                povit.clear(); keepC.clear();
            }
        }
        // Deferred bucket move sweep (integer arithmetic — no rounding needed)
        for (auto idx : dirtyEdges) {
            if (countingKE[idx] < 0) countingKE[idx] = 0;
            if (edgeInHeap[idx]) bucketMove(idx);
            dirtyMark[idx] = false;
        }
        dirtyEdges.clear();
        currentRemoveEdgeIds.clear();
        // Reset leafAffected flags for next iteration
        for (auto leafId : removedLeaf) leafAffected[leafId] = false;
        removedLeaf.clear();
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
    std::cout << "  Cases: A=" << cntA << " B=" << cntB << " C=" << cntC << " iters=" << totalIters << std::endl;

    // for (auto i = 0;  i < edgeGraph.adj_list.size(); ++i) {
    //     auto counting = countingKE[i];
    //     if (counting != 0) {
    //         std::cerr << "Error: countingKE != 0" << std::endl;
    //         std::cerr << "countingKE: " << counting << std::endl;
    //         std::exit(1);
    //     }
    // }

    // coreE
    // daf::printArray(coreE, edgeGraph.adj_list.size());


    // ~/_/pivoter/a
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

#ifdef _OPENMP
    for (int i = 0; i < NUM_VERTEX_LOCKS; i++) omp_destroy_lock(&vertexLocks[i]);
    for (int i = 0; i < NUM_LEAF_LOCKS; i++) omp_destroy_lock(&leafLocks[i]);
#endif
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
    bkRmEdge::for_each_bit(bs, (int) bs.size(), [&](int v) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << v;
        return true;
    });
    std::cout << "]\n";
}