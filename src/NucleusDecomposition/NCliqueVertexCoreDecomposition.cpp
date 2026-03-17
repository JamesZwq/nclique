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
#include "dataStruct/coreDisJoin.hpp"
#include "graph/DynamicGraphSet.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];
// （），
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
        const double *vertexCounting; // 
        explicit CompareVertex(const double *coreLeaf) : vertexCounting(coreLeaf) {
        }

        // ： “a ” ， coreLeaf[a] > coreLeaf[b]
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

std::chrono::steady_clock::time_point get_value(std::chrono::steady_clock::time_point time_start) {
    return time_start;
}

double *  NCliqueVertexCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    auto countingV = VCD::countingPerVertex(tree, edgeGraph, k);
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;
    auto coreV = new double[numVertices + 1];
    // Initialize to -1 to distinguish non-participating vertices from core=0
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;

#ifdef _OPENMP
    const int numThreads = omp_get_max_threads();
#else
    const int numThreads = 1;
#endif

    // --- Bucket array replacing O(log n) heap ---
    int maxBucket = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] > 0)
            maxBucket = std::max(maxBucket, (int)countingV[i]);
    }
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(numVertices);
    std::vector<daf::Size> pos_in_bucket(numVertices);
    std::vector<bool> vertexInHeap(numVertices, false);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] <= 0) continue;
        int b = (int)countingV[i];
        bucket_of[i] = b;
        pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i);
        vertexInHeap[i] = true;
        remainingInHeap++;
    }
    int curBucket = 0;

    // Bucket move helper (serial only)
    auto bucketMove = [&](daf::Size id) {
        int newB = std::max(0, (int)countingV[id]);
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

    // Deferred bucket move tracking
    std::vector<bool> dirtyMark(numVertices, false);
    std::vector<daf::Size> dirtyVertices;
    dirtyVertices.reserve(8192);

    // Per-thread dirty lists for parallel phases
    std::vector<std::vector<daf::Size>> threadDirty(numThreads);
    for (auto& v : threadDirty) v.reserve(4096);

    // Per-leaf locks for Phase 1
    static constexpr int NUM_LEAF_LOCKS = 4096;
#ifdef _OPENMP
    omp_lock_t leafLocks[NUM_LEAF_LOCKS];
    for (int i = 0; i < NUM_LEAF_LOCKS; i++) omp_init_lock(&leafLocks[i]);
#endif

    // Leaf tracking
    daf::StaticVector<daf::Size> removedLeaf(tree.adj_list.size());
    daf::StaticVector<VCD::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();
    std::vector<bool> leafAffected(tree.adj_list.size(), false);

    // Per-thread new leaf lists for Phase 1
    std::vector<std::vector<daf::Size>> threadNewLeaves(numThreads);
    for (auto& v : threadNewLeaves) v.reserve(4096);

    daf::StaticVector<daf::Size> currentRemoveVertexIds(numVertices);

    std::cout << "=========================begin=========================" << std::endl;
    std::cout << "vertices in heap: " << remainingInHeap << std::endl;

    double minCore = 0;
    // Pre-allocate per-thread recycle lists outside the loop
    std::vector<std::vector<daf::Size>> threadRecycle(numThreads);

    while (remainingInHeap > 0) {
        // --- Bucket pop ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((double)curBucket, minCore);

        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && (double)curBucket <= minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                vertexInHeap[id] = false;
                currentRemoveVertexIds.push_back(id);
                coreV[id] = minCore;
                remainingInHeap--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (double)(curBucket + 1) <= minCore) {
                curBucket++;
            } else break;
        }

        if (remainingInHeap == 0) break;

        // --- Phase 1: identify affected leaves ---
        if (numThreads > 1 && currentRemoveVertexIds.size() > 16) {
            // Parallel path with per-leaf locks
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 8)
#endif
            for (int vi = 0; vi < (int)currentRemoveVertexIds.size(); ++vi) {
#ifdef _OPENMP
                int tid = omp_get_thread_num();
#else
                int tid = 0;
#endif
                auto v = currentRemoveVertexIds[vi];
                auto &adjClique = treeGraphV.getNbr(v);
                for (const auto &clique : adjClique) {
                    daf::Size leafId = clique.v;
                    if (tree.adj_list[leafId].empty()) continue;
                    auto lockIdx = leafId & (NUM_LEAF_LOCKS - 1);
#ifdef _OPENMP
                    omp_set_lock(&leafLocks[lockIdx]);
#endif
                    auto &lrm = leafRmInfo[leafId];
                    bool wasEmpty = lrm.empty();
                    if (!lrm.removedKeepC) {
                        if (!clique.isPivot) {
                            lrm.removedKeepC = true;
                        } else {
                            lrm.removedPivots.push_back(v);
                        }
                    }
#ifdef _OPENMP
                    omp_unset_lock(&leafLocks[lockIdx]);
#endif
                    if (wasEmpty) threadNewLeaves[tid].push_back(leafId);
                }
            }
            // Merge per-thread new leaves
            for (int t = 0; t < numThreads; ++t) {
                for (auto leafId : threadNewLeaves[t]) {
                    if (!leafAffected[leafId]) {
                        leafAffected[leafId] = true;
                        removedLeaf.push_back(leafId);
                    }
                }
                threadNewLeaves[t].clear();
            }
        } else {
            // Serial path (no locks, no OMP overhead)
            for (int vi = 0; vi < (int)currentRemoveVertexIds.size(); ++vi) {
                auto v = currentRemoveVertexIds[vi];
                auto &adjClique = treeGraphV.getNbr(v);
                for (const auto &clique : adjClique) {
                    daf::Size leafId = clique.v;
                    if (tree.adj_list[leafId].empty()) continue;
                    auto &lrm = leafRmInfo[leafId];
                    if (lrm.empty()) {
                        removedLeaf.push_back(leafId);
                        leafAffected[leafId] = true;
                    }
                    if (!lrm.removedKeepC) {
                        if (!clique.isPivot) {
                            lrm.removedKeepC = true;
                        } else {
                            lrm.removedPivots.push_back(v);
                        }
                    }
                }
            }
        }

        // --- Phase 2 (parallel): sort + countingV deltas + tree mutations ---
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 64) if(removedLeaf.size() > 32)
#endif
        for (int li = 0; li < (int)removedLeaf.size(); ++li) {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            auto leafId = removedLeaf[li];
            auto &leaf = tree.adj_list[leafId];
            VCD::LeafRmInfo &leafRm = leafRmInfo[leafId];
            if (leaf.empty()) continue;

            // Sort removedPivots inline
            std::ranges::sort(leafRm.removedPivots);
            leafRm.removedPivots.unique();

            int numPivots = 0, numKeeps = 0;
            for (auto &node : leaf) {
                if (node.isPivot) numPivots++;
                else numKeeps++;
            }

            int needPivot = (int)k - numKeeps;
            int remainPivots = numPivots - (int)leafRm.removedPivots.size();
            bool leafDies = leafRm.removedKeepC || needPivot < 0 || needPivot > remainPivots;

            if (leafDies) {
                double keepValue = (needPivot >= 0 && needPivot <= numPivots)
                                   ? nCr[numPivots][needPivot] : 0.0;
                double pivotValue = (needPivot >= 1 && needPivot - 1 <= numPivots - 1)
                                    ? nCr[numPivots - 1][needPivot - 1] : 0.0;

                for (auto &node : leaf) {
                    double delta = node.isPivot ? pivotValue : keepValue;
                    if (delta > 0) {
#ifdef _OPENMP
                        #pragma omp atomic
#endif
                        countingV[node.v] -= delta;
                        if (vertexInHeap[node.v])
                            threadDirty[tid].push_back(node.v);
                    }
                }
                leaf.clear();
                threadRecycle[tid].push_back(leafId);
            } else if (!leafRm.removedPivots.empty() && needPivot <= remainPivots) {
                double KeepDelta = nCr[numPivots][needPivot] - nCr[remainPivots][needPivot];
                double RemovedPivotFull = (needPivot >= 1) ? nCr[numPivots - 1][needPivot - 1] : 0.0;
                double PivotDelta = (needPivot >= 1)
                    ? nCr[numPivots - 1][needPivot - 1] - nCr[remainPivots - 1][needPivot - 1]
                    : 0.0;

                for (auto rp : leafRm.removedPivots) {
                    if (RemovedPivotFull > 0) {
#ifdef _OPENMP
                        #pragma omp atomic
#endif
                        countingV[rp] -= RemovedPivotFull;
                        if (vertexInHeap[rp])
                            threadDirty[tid].push_back(rp);
                    }
                }

                for (auto &node : leaf) {
                    if (node.isPivot && std::binary_search(
                            leafRm.removedPivots.begin(),
                            leafRm.removedPivots.end(), node.v))
                        continue;
                    double delta = node.isPivot ? PivotDelta : KeepDelta;
                    if (delta > 0) {
#ifdef _OPENMP
                        #pragma omp atomic
#endif
                        countingV[node.v] -= delta;
                        if (vertexInHeap[node.v])
                            threadDirty[tid].push_back(node.v);
                    }
                }
                tree.removeNbrs(leafId, leafRm.removedPivots);
            }
        }

        // --- Serial: merge recycle + dirty + bucket moves + cleanup ---
        for (int t = 0; t < numThreads; ++t) {
            for (auto id : threadRecycle[t])
                tree.recycleNode(id);
            threadRecycle[t].clear();
        }
        for (int t = 0; t < numThreads; ++t) {
            for (auto v : threadDirty[t]) {
                if (!dirtyMark[v]) {
                    dirtyMark[v] = true;
                    dirtyVertices.push_back(v);
                }
            }
            threadDirty[t].clear();
        }
        for (auto v : dirtyVertices) {
            countingV[v] = std::max(countingV[v], 0.0);
            if (vertexInHeap[v]) bucketMove(v);
            dirtyMark[v] = false;
        }
        dirtyVertices.clear();

        for (int li = 0; li < (int)removedLeaf.size(); ++li) {
            leafRmInfo[removedLeaf[li]].clear();
            leafAffected[removedLeaf[li]] = false;
        }
        currentRemoveVertexIds.clear();
        removedLeaf.clear();
    }

    // Cleanup locks
#ifdef _OPENMP
    for (int i = 0; i < NUM_LEAF_LOCKS; i++) omp_destroy_lock(&leafLocks[i]);
#endif

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    delete[] countingV;
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