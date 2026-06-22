//
// Created by _ on 25-3-4.
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
#include "BK/BronKerboschRmRClique.hpp"
#include "dataStruct/disJoinSet.hpp"
#include "graph/DynamicGraphSet.h"
// timing
#include <chrono>
#include <unordered_map>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];
// §105 M1: per-vertex class id (set in degeneracy_cliques.cpp under
// PIVOTER_M1_TUPLE_PROBE). Empty unless the probe is active.
extern std::vector<int> g_m1ClassOf;

namespace CDSetRS {
    template<typename It1, typename It2, typename UpdateFunc>
    inline void processEdgePairsImpl(It1 b1, It1 e1,
                                     It2 b2, It2 e2,
                                     double weight,
                                     UpdateFunc &&upd) noexcept {
        if (weight < 0.0) return;
        if (b1 == b2 && e1 == e2 && b1 == e1 && b2 == e2) {
            return;
        }
        if (b1 == b2 && e1 == e2) {
            for (auto it = b1; it + 1 != e1; ++it) {
                auto u = *it;
                for (auto jt = it + 1; jt != e1; ++jt) {
                    upd(u, *jt, weight);
                }
            }
        } else {
            for (auto it = b1; it != e1; ++it) {
                auto u = *it;
                for (auto jt = b2; jt != e2; ++jt) {
                    upd(u, *jt, weight);
                }
            }
        }
    }

    template<typename Range1, typename Range2, typename UpdateFunc>
    inline void processEdgePairs(const Range1 &r1, const Range2 &r2, double weight, UpdateFunc &&upd) noexcept {
        processEdgePairsImpl(std::begin(r1), std::end(r1), std::begin(r2), std::end(r2), weight, std::forward<UpdateFunc>(upd));
    }

    template<typename Range, typename UpdateFunc>
    inline void processEdgePairs(const Range &r, double weight, UpdateFunc &&upd) noexcept {
        processEdgePairsImpl(std::begin(r), std::end(r), std::begin(r), std::end(r), weight, std::forward<UpdateFunc>(upd));
    }


    struct CompareRClique {
        const double *RCliqueCounting;
        explicit CompareRClique(const double *coreLeaf) : RCliqueCounting(coreLeaf) {}
        bool operator()(daf::Size const &a, daf::Size const &b) const {
            if (RCliqueCounting[a] != RCliqueCounting[b])
                return RCliqueCounting[a] > RCliqueCounting[b];
            return a < b; // deterministic tie-break: smaller id has higher priority (popped first)
        }
    };

    using DHeap = boost::heap::d_ary_heap<
        daf::Size,
        boost::heap::arity<8>,
        boost::heap::mutable_<true>,
        boost::heap::compare<CompareRClique>
    >;

    std::vector<double> countingPerRClique(
        const DynamicGraph<TreeGraphNode> &treeGraph,
        StaticCliqueIndex &cliqueHashmap,
        const daf::CliqueSize r,
        const daf::CliqueSize s) {
        std::vector<double> rCliqueSCounting(cliqueHashmap.size(), 0.0);
        for (const auto &leaf: treeGraph.adj_list) {
            if (leaf.size() < r) continue;

            daf::CliqueSize pivotC = 0, keepC = 0;
            for (auto &i: leaf) {
                if (i.isPivot) pivotC++;
                else keepC++;
            }

            int needPivot = s - static_cast<int>(keepC);
            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                daf::CliqueSize subNumPovit = 0;
                for (const auto &node: rClique) {
                    if (node.isPivot) subNumPovit++;
                }

                auto ncrValue = nCr[pivotC - subNumPovit][needPivot - subNumPovit];
                auto [id, isNew] = cliqueHashmap.byNewClique(rClique);
                if (isNew) {
                    if (rCliqueSCounting.size() <= id) rCliqueSCounting.push_back(0.0);
                    if (rCliqueSCounting.capacity() <= id) rCliqueSCounting.reserve(std::max(id + 2, id * 2));
                }
                rCliqueSCounting[id] += ncrValue;
                return true;
            });
        }
        rCliqueSCounting.shrink_to_fit();
        return rCliqueSCounting;
    }

    // Parallel counting after index is built: each thread sums into local array, then merge (same result as countingPerRClique).
    std::vector<double> countingPerRCliqueParallel(
        const DynamicGraph<TreeGraphNode> &treeGraph,
        const StaticCliqueIndex &cliqueIndex,
        const daf::CliqueSize r,
        const daf::CliqueSize s) {
        const daf::Size nClique = cliqueIndex.size();
        std::vector<double> rCliqueSCounting(nClique, 0.0);
#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
        std::vector<std::vector<double>> thread_locals(nthreads, std::vector<double>(nClique, 0.0));
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto &local = thread_locals[tid];
#pragma omp for schedule(dynamic, 64)
            for (daf::Size leafIdx = 0; leafIdx < treeGraph.adj_list.size(); ++leafIdx) {
                const auto &leaf = treeGraph.adj_list[leafIdx];
                if (leaf.size() < r) continue;
                daf::CliqueSize pivotC = 0, keepC = 0;
                for (const auto &i : leaf) {
                    if (i.isPivot) pivotC++;
                    else keepC++;
                }
                int needPivot = s - static_cast<int>(keepC);
                daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                    daf::CliqueSize subNumPovit = 0;
                    for (const auto &node : rClique) if (node.isPivot) subNumPovit++;
                    double ncrValue = nCr[pivotC - subNumPovit][needPivot - subNumPovit];
                    auto id = cliqueIndex.byClique(rClique);
                    if (id < nClique) local[id] += ncrValue;
                    return true;
                });
            }
        }
        // 并行合并
#pragma omp parallel for schedule(static)
        for (daf::Size i = 0; i < nClique; ++i) {
            for (int t = 0; t < nthreads; ++t) {
                rCliqueSCounting[i] += thread_locals[t][i];
            }
        }
#else
        for (const auto &leaf : treeGraph.adj_list) {
            if (leaf.size() < r) continue;
            daf::CliqueSize pivotC = 0, keepC = 0;
            for (const auto &i : leaf) { if (i.isPivot) pivotC++; else keepC++; }
            int needPivot = s - static_cast<int>(keepC);
            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                daf::CliqueSize subNumPovit = 0;
                for (const auto &node : rClique) if (node.isPivot) subNumPovit++;
                rCliqueSCounting[cliqueIndex.byClique(rClique)] += nCr[pivotC - subNumPovit][needPivot - subNumPovit];
                return true;
            });
        }
#endif
        return rCliqueSCounting;
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
    void printEdgeCore(const Graph &edgeGraph, const std::vector<T> &coreE) {
        printEdgeCore(edgeGraph, coreE.data());
    }
}


std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionRClique(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    // ============ TIMERS ============
    long long duration_init = 0;
    long long duration_pop = 0;
    long long duration_structure = 0;
    long long duration_support = 0;
    long long duration_heap = 0;
    // ================================

    auto time_start = std::chrono::high_resolution_clock::now();

    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });

    daf::log_memory("r-clique index");

    auto countingRClique = daf::timeCount("countingPerRClique", [&]() {
        return CDSetRS::countingPerRCliqueParallel(tree, cliqueIndex, r, s);
    });

    // §105 M1: tag each r-clique with its TUPLE = sorted class-multiset, then
    // count distinct tuples (must == region_native #patterns) and total
    // support-bearing r-cliques (must == region_native #r-cliques). r-cliques
    // with zero s-clique support (in no maximal clique >= s) are NOT region_native
    // patterns, so they are skipped -- mirroring region_native's region-derived
    // enumeration. Read-only, env-gated, never touches coreRClique -> corehash unchanged.
    if (std::getenv("PIVOTER_M1_TUPLE_PROBE") && !g_m1ClassOf.empty()) {
        std::unordered_map<std::string, long long> tupleMap;
        long long totalRC = 0, skipNoSupport = 0, skipNoClass = 0;
        std::vector<int> cls;
        std::string key;
        for (daf::Size id = 0; id < cliqueIndex.size(); ++id) {
            if (id >= countingRClique.size() || countingRClique[id] <= 0) { skipNoSupport++; continue; }
            auto span = cliqueIndex.byId(id);
            cls.clear();
            bool ok = true;
            for (daf::Size v : span) {
                int c = (v < g_m1ClassOf.size()) ? g_m1ClassOf[v] : -1;
                if (c < 0) { ok = false; break; }
                cls.push_back(c);
            }
            if (!ok) { skipNoClass++; continue; }
            std::sort(cls.begin(), cls.end());
            key.clear();
            key.reserve(cls.size() * sizeof(int));
            for (int x : cls) key.append((const char *)&x, sizeof(int));
            tupleMap[key]++;
            totalRC++;
        }
        printf("[m1-tuple] patterns=%zu  r-cliques=%lld  (indexed=%lld skipNoSupport=%lld skipNoClass=%lld)\n",
               tupleMap.size(), totalRC, (long long)cliqueIndex.size(), skipNoSupport, skipNoClass);
        fflush(stdout);
    }

    std::vector<double> coreRClique(countingRClique.size(), 0);
    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<std::vector<daf::Size> > removedRCliqueIdForLeaf;
    std::vector<daf::Size> changedLeaf;
    std::vector<daf::Size> currentRemoveRcliqueIds;

    removedRCliqueIdForLeaf.reserve(tree.adj_list.size());
    changedLeaf.reserve(tree.adj_list.size());
    currentRemoveRcliqueIds.reserve(cliqueIndex.size());

    daf::StaticVector<bool> rCliqueInHeap(cliqueIndex.size());
    rCliqueInHeap.resize(cliqueIndex.size());
    memset(rCliqueInHeap.getData(), true, cliqueIndex.size() * sizeof(bool));

    int maxBucket = 0;
    for (daf::Size i = 0; i < cliqueIndex.size(); ++i)
        maxBucket = std::max(maxBucket, (int)countingRClique[i]);
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(cliqueIndex.size());
    std::vector<daf::Size> pos_in_bucket(cliqueIndex.size());
    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
        int b = (int)countingRClique[i];
        bucket_of[i] = b;
        pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i);
    }
    int curBucket = 0;
    daf::Size remainingInHeap = cliqueIndex.size();

    daf::log_memory("Other index(incloud head)");

    // Measure Init Time
    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - time_start).count();

    // O5: Pre-allocate LeafResult struct outside loop
    struct LeafResult {
        std::vector<std::vector<TreeGraphNode>> newLeaves;
        std::vector<std::pair<daf::Size, double>> incr;
        std::vector<std::pair<daf::Size, double>> decr;
        void clear() { newLeaves.clear(); incr.clear(); decr.clear(); }
    };

    const daf::Size graphN = edgeGraph.n;

    double minCore = 0;
    while (remainingInHeap > 0) {
        // [Timer] Start Pop
        auto t_loop_start = std::chrono::high_resolution_clock::now();

        for (auto &leafId: changedLeaf) {
            changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        }
        changedLeaf.clear();
        removedRCliqueIdForLeaf.clear();
        currentRemoveRcliqueIds.clear();

        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((double)curBucket, minCore);

        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && (double)curBucket <= minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                rCliqueInHeap[id] = false;
                currentRemoveRcliqueIds.push_back(id);
                coreRClique[id] = minCore;
                remainingInHeap--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (double)(curBucket + 1) <= minCore) {
                curBucket++;
            } else break;
        }

        // [Timer] End Pop
        duration_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_loop_start).count();

        if (remainingInHeap == 0) break;

        // [Timer] Start Structure (Part A - Intersect)
        auto t_struct_A = std::chrono::high_resolution_clock::now();

#ifdef _OPENMP
        using PairOV = std::pair<daf::Size, daf::Size>;
        int nthreads = omp_get_max_threads();
        std::vector<std::vector<PairOV>> thread_pairs(nthreads);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            thread_pairs[tid].reserve(2048);
#pragma omp for schedule(dynamic, 64)
            for (daf::Size origIdx = 0; origIdx < currentRemoveRcliqueIds.size(); ++origIdx) {
                auto rmRCliqueId = currentRemoveRcliqueIds[origIdx];
                auto rClique = cliqueIndex.byId(rmRCliqueId);
                daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                    [&](const TreeGraphNode &uClique) {
                        thread_pairs[tid].emplace_back(origIdx, uClique.v);
                    });
            }
        }
        size_t total_size = 0;
        for (const auto &tp : thread_pairs) total_size += tp.size();
        std::vector<PairOV> allPairs;
        allPairs.reserve(total_size);
        for (auto &tp : thread_pairs) {
            allPairs.insert(allPairs.end(), tp.begin(), tp.end());
        }
        std::sort(allPairs.begin(), allPairs.end());
        for (const auto &p : allPairs) {
            daf::Size origIdx = p.first;
            daf::Size leafV = p.second;
            auto rmRCliqueId = currentRemoveRcliqueIds[origIdx];
            auto &leafId = changedLeafIndex[leafV];
            if (leafId == std::numeric_limits<daf::Size>::max()) {
                leafId = removedRCliqueIdForLeaf.size();
                removedRCliqueIdForLeaf.emplace_back();
                changedLeaf.push_back(leafV);
                removedRCliqueIdForLeaf.back().reserve(10);
            }
            removedRCliqueIdForLeaf[leafId].emplace_back(rmRCliqueId);
        }
#else
        for (auto rmRCliqueId: currentRemoveRcliqueIds) {
            auto rClique = cliqueIndex.byId(rmRCliqueId);
            daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                                            [&](const TreeGraphNode &uClique) {
                                                auto &leafId = changedLeafIndex[uClique.v];
                                                if (leafId == std::numeric_limits<daf::Size>::max()) {
                                                    leafId = removedRCliqueIdForLeaf.size();
                                                    removedRCliqueIdForLeaf.emplace_back();
                                                    changedLeaf.push_back(uClique.v);
                                                    removedRCliqueIdForLeaf.back().reserve(10);
                                                }
                                                removedRCliqueIdForLeaf[leafId].emplace_back(rmRCliqueId);
                                            });
        }
#endif

        // [Timer] End Structure A
        duration_structure += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_struct_A).count();

        // O5: Resize leafResults to match, clearing old data
        std::vector<LeafResult> leafResults(changedLeaf.size());

#ifdef _OPENMP
        // ============ Phase B: Parallel BK + enumeration ============
        auto t_par = std::chrono::high_resolution_clock::now();
#pragma omp parallel
        {
            static thread_local daf::StaticVector<daf::Size> threadMap;
            if (threadMap.maxSize < graphN + 1) threadMap.resize(graphN + 1);
#pragma omp for schedule(dynamic, 64)
            for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
                auto leafId = changedLeaf[idx];
                const auto &leaf = tree.adj_list[leafId];
                auto leafIndex = changedLeafIndex[leafId];
                const auto &removedR = removedRCliqueIdForLeaf[leafIndex];
                auto mapped = removedR | std::views::transform([&](const daf::Size id) { return cliqueIndex.byId(id); });

                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                        auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                        leafResults[idx].newLeaves.push_back(std::move(newLeaf));
                    },
                    &threadMap);

                for (const auto &newLeaf : leafResults[idx].newLeaves) {
                    daf::CliqueSize newPivotC = 0, newKeepC = 0;
                    for (const auto &j : newLeaf) { if (j.isPivot) newPivotC++; else newKeepC++; }
                    daf::Size needPivot = s - newKeepC;
                    daf::enumerateCombinations(newLeaf, r, [&](const daf::StaticVector<TreeGraphNode> &rclique) {
                        daf::CliqueSize subNumPovit = 0;
                        for (const auto &node : rclique) if (node.isPivot) subNumPovit++;
                        if (subNumPovit <= needPivot && newPivotC - subNumPovit >= 0 && newPivotC - subNumPovit < 1001 &&
                            needPivot - subNumPovit >= 0 && needPivot - subNumPovit < 401) {
                            double ncrValue = nCr[newPivotC - subNumPovit][needPivot - subNumPovit];
                            leafResults[idx].incr.emplace_back(cliqueIndex.byClique(rclique), ncrValue);
                        }
                        return true;
                    });
                }
                daf::CliqueSize keepC = 0, pivotC = 0;
                for (const auto &node : leaf) { if (node.isPivot) pivotC++; else keepC++; }
                daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &clique) {
                    auto cliqueIndexId = cliqueIndex.byClique(clique);
                    if (!rCliqueInHeap[cliqueIndexId]) return true;
                    daf::CliqueSize subNumPovit = 0;
                    for (const auto &node : clique) if (node.isPivot) subNumPovit++;
                    leafResults[idx].decr.emplace_back(cliqueIndexId, nCr[pivotC - subNumPovit][s - keepC - subNumPovit]);
                    return true;
                });
            }
        }
        duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_par).count();

        // ============ Phase C: Per-leaf serial mutations + support + heap ============
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            auto leafId = changedLeaf[idx];
            const auto &leaf = tree.adj_list[leafId];
            LeafResult &res = leafResults[idx];

            auto t_struct_B = std::chrono::high_resolution_clock::now();
            for (const auto &leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }
            for (auto &newLeaf : res.newLeaves) {
                auto newId = tree.addNode(newLeaf);
                const auto &stored = tree.adj_list[newId];
                for (const auto &i : stored) {
                    if (i.isPivot) treeGraphV.addNbr(i.v, {newId, true});
                    else treeGraphV.addNbr(i.v, {newId, false});
                }
                if (newId >= changedLeafIndex.size())
                    changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
            }
            duration_structure += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_struct_B).count();

            auto t_supp = std::chrono::high_resolution_clock::now();
            for (const auto &p : res.incr) countingRClique[p.first] += p.second;
            duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_supp).count();

            auto t_heap = std::chrono::high_resolution_clock::now();
            for (const auto &p : res.decr) {
                countingRClique[p.first] -= p.second;
                int newB = std::max(0, (int)countingRClique[p.first]);
                int oldB = bucket_of[p.first];
                if (newB != oldB) {
                    auto& oldVec = buckets[oldB];
                    daf::Size myPos = pos_in_bucket[p.first];
                    if (myPos < oldVec.size() - 1) {
                        daf::Size last = oldVec.back();
                        oldVec[myPos] = last;
                        pos_in_bucket[last] = myPos;
                    }
                    oldVec.pop_back();
                    if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
                    bucket_of[p.first] = newB;
                    pos_in_bucket[p.first] = buckets[newB].size();
                    buckets[newB].push_back(p.first);
                    if (newB < curBucket) curBucket = newB;
                }
            }
            duration_heap += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_heap).count();

            auto t_struct_C = std::chrono::high_resolution_clock::now();
            tree.removeNode(leafId);
            duration_structure += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_struct_C).count();
        }

#else
        for (auto leafId: changedLeaf) {
            auto t_struct_B = std::chrono::high_resolution_clock::now();
            const auto &leaf = tree.adj_list[leafId];
            auto leafIndex = changedLeafIndex[leafId];
            auto initCore = [&](const std::vector<TreeGraphNode> &newLeaf, const daf::Size &newLeafId) {
                daf::CliqueSize newPivotC = 0, newKeepC = 0;
                for (const auto &i: newLeaf) {
                    if (i.isPivot) { treeGraphV.addNbr(i.v, {newLeafId, true}); newPivotC++; }
                    else { treeGraphV.addNbr(i.v, {newLeafId, false}); newKeepC++; }
                }
                daf::Size needPivot = s - newKeepC;
                daf::enumerateCombinations(newLeaf, r, [&](const daf::StaticVector<TreeGraphNode> &rclique) {
                    daf::CliqueSize subNumPovit = 0;
                    for (const auto &node: rclique) if (node.isPivot) ++subNumPovit;
                    if (subNumPovit <= needPivot) {
                        countingRClique[cliqueIndex.byClique(rclique)] += nCr[newPivotC - subNumPovit][needPivot - subNumPovit];
                    }
                    return true;
                });
            };
            for (const auto &leafV: leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }
            const auto &removedR = removedRCliqueIdForLeaf[leafIndex];
            auto mapped = removedR | std::views::transform([&](const daf::Size id) { return cliqueIndex.byId(id); });
            bkRmClique::removeRClique(leaf, mapped, r, s, [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                auto newId = tree.addNode(newLeaf);
                initCore(tree.adj_list[newId], newId);
                if (newId >= changedLeafIndex.size()) changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
            });
            duration_structure += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_struct_B).count();
            daf::CliqueSize keepC = 0, pivotC = 0;
            for (const auto &node: leaf) { if (node.isPivot) ++pivotC; else ++keepC; }
            auto t_supp_dec = std::chrono::high_resolution_clock::now();
            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &clique) {
                auto cliqueIndexId = cliqueIndex.byClique(clique);
                if (!rCliqueInHeap[cliqueIndexId]) return true;
                daf::CliqueSize subNumPovit = 0;
                for (const auto &node: clique) if (node.isPivot) ++subNumPovit;
                countingRClique[cliqueIndexId] -= nCr[pivotC - subNumPovit][s - keepC - subNumPovit];
                int newB = std::max(0, (int)countingRClique[cliqueIndexId]);
                int oldB = bucket_of[cliqueIndexId];
                if (newB != oldB) {
                    auto& oldVec = buckets[oldB];
                    daf::Size myPos = pos_in_bucket[cliqueIndexId];
                    if (myPos < oldVec.size() - 1) {
                        daf::Size last = oldVec.back();
                        oldVec[myPos] = last;
                        pos_in_bucket[last] = myPos;
                    }
                    oldVec.pop_back();
                    if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
                    bucket_of[cliqueIndexId] = newB;
                    pos_in_bucket[cliqueIndexId] = buckets[newB].size();
                    buckets[newB].push_back(cliqueIndexId);
                    if (newB < curBucket) curBucket = newB;
                }
                return true;
            });
            duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_supp_dec).count();
            auto t_struct_C = std::chrono::high_resolution_clock::now();
            tree.removeNode(leafId);
            duration_structure += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t_struct_C).count();
        }
#endif
    }


    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    // Output Detailed Breakdown
    std::cout << "Time Breakdown (ms):" << std::endl;
    std::cout << "  Init:      " << duration_init / 1000000.0 << std::endl;
    std::cout << "  Pop:       " << duration_pop / 1000000.0 << std::endl;
    std::cout << "  Structure: " << duration_structure / 1000000.0 << std::endl;
    std::cout << "  Support:   " << duration_support / 1000000.0 << std::endl;
    std::cout << "  Heap:      " << duration_heap / 1000000.0 << std::endl;

    std::vector<std::pair<std::vector<daf::Size>, double> > sortedK;
    sortedK.reserve(countingRClique.size());

    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueCopy(clique.begin(), clique.end());
        sortedK.emplace_back(cliqueCopy, coreRClique[i]);
    }
    std::sort(sortedK.begin(), sortedK.end(),
              [](const auto &a, const auto &b) {
                  return a.second < b.second;
              });

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