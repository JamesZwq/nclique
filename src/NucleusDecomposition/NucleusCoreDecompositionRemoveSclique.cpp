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
#include <unordered_set>
#include <string>
#include <functional>
#include <array>
#include <span>
#include <algorithm>
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

    // §105 M2.6 (validate the M3 premise BEFORE the rewrite): does CND's per-r-clique
    // peel assign the SAME core to every r-clique of a tuple? If yes, whole-tuple
    // peeling is provably correct (rawSupport_T = sum support_c = mult_T*support_T,
    // bucket by support_T). Snapshot the support>0 set now (region_native scope);
    // check uniformity post-peel. Read-only, env-gated.
    const bool m3probe = std::getenv("PIVOTER_M3_INVARIANT_PROBE") && !g_m1ClassOf.empty();
    std::vector<char> m3HadSupport;
    if (m3probe) {
        m3HadSupport.assign(cliqueIndex.size(), 0);
        for (daf::Size id = 0; id < cliqueIndex.size() && id < countingRClique.size(); ++id)
            m3HadSupport[id] = countingRClique[id] > 0 ? 1 : 0;
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

    // §105 M2 (measure-before-build): per-round reprocess accounting. Counts the
    // REAL #(r-clique,leaf) reprocess events CND pays (incr + decr enumerations)
    // vs the #(distinct-tuple,leaf) a tuple-batched peel would pay. The ratio is
    // the TRUE per-leaf reprocess saving; the global #rclq/#pat is only an upper
    // bound. OMP-path only (uses res.incr/res.decr); read-only on peel state;
    // env-gated -> corehash unchanged.
    const bool m2on = std::getenv("PIVOTER_M2_REPROCESS_PROBE") && !g_m1ClassOf.empty();
    long long m2_rcliqueReproc = 0, m2_tupleReproc = 0, m2_leafEvents = 0, m2_rounds = 0;
    auto m2KeyOf = [&](daf::Size cid, std::vector<int> &cls, std::string &key) {
        auto span = cliqueIndex.byId(cid);
        cls.clear();
        for (daf::Size v : span) cls.push_back(v < g_m1ClassOf.size() ? g_m1ClassOf[v] : -1);
        std::sort(cls.begin(), cls.end());
        key.clear();
        key.reserve(cls.size() * sizeof(int));
        for (int x : cls) key.append((const char *)&x, sizeof(int));
    };

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
        if (m2on) m2_rounds++;

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

            if (m2on) {
                // incr (new sub-leaves) and decr (old leaf) are separate batched
                // groups; distinct tuples within each are what a tuple-batched
                // peel would touch.
                std::unordered_set<std::string> incrTup, decrTup;
                std::vector<int> cls;
                std::string key;
                for (const auto &p : res.incr) { m2KeyOf(p.first, cls, key); incrTup.insert(key); }
                for (const auto &p : res.decr) { m2KeyOf(p.first, cls, key); decrTup.insert(key); }
                m2_rcliqueReproc += (long long)res.incr.size() + (long long)res.decr.size();
                m2_tupleReproc   += (long long)incrTup.size() + (long long)decrTup.size();
                m2_leafEvents++;
            }

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

    if (m3probe) {
        // group final cores by tuple; a tuple is uniform iff all its r-cliques
        // share one core (the §105 whole-tuple-peeling premise).
        std::unordered_map<std::string, double> tupMin, tupMax;
        std::vector<int> cls;
        std::string key;
        long long checked = 0;
        for (daf::Size id = 0; id < cliqueIndex.size(); ++id) {
            if (id >= m3HadSupport.size() || !m3HadSupport[id]) continue;
            m2KeyOf(id, cls, key);
            double cv = coreRClique[id];
            auto it = tupMin.find(key);
            if (it == tupMin.end()) { tupMin.emplace(key, cv); tupMax.emplace(key, cv); }
            else { it->second = std::min(it->second, cv); auto &mx = tupMax[key]; mx = std::max(mx, cv); }
            checked++;
        }
        long long nonUniform = 0;
        double maxSpread = 0;
        for (auto &kv : tupMin) {
            double spread = tupMax[kv.first] - kv.second;
            if (spread > 0) { nonUniform++; maxSpread = std::max(maxSpread, spread); }
        }
        printf("[m3-invariant] tuples=%zu  checked_rcliques=%lld  nonUniformTuples=%lld  maxCoreSpread=%.0f  %s\n",
               tupMin.size(), checked, nonUniform, maxSpread,
               nonUniform == 0 ? "=> PREMISE HOLDS (whole-tuple peel valid)"
                               : "=> PREMISE VIOLATED (per-tuple peel would be WRONG)");
        fflush(stdout);
    }

    if (m2on) {
        printf("[m2-reproc] rclique_reproc=%lld  tuple_reproc=%lld  ratio=%.2fx"
               "  (leafEvents=%lld rounds=%lld)\n",
               m2_rcliqueReproc, m2_tupleReproc,
               m2_tupleReproc ? (double)m2_rcliqueReproc / (double)m2_tupleReproc : 0.0,
               m2_leafEvents, m2_rounds);
        if (m2_leafEvents == 0)
            printf("[m2-reproc] NOTE: 0 leaf events -- serial (#else) path not instrumented; run with OpenMP.\n");
        fflush(stdout);
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


// ============================================================
// §105 M3: tuple-batched r-clique nucleus decomposition (PIVOTER_RUN_TUPLE_BATCH).
// Peels whole TUPLES (class-multisets) instead of individual r-cliques. By the
// M2.6-validated symmetry, all r-cliques of a tuple share one support at ALL
// times, so CND would peel them in the same round. We maintain ONE counter per
// tuple: rawSupport_T = sum_{c in T} support_c = mult_T * support_T, bucket by
// support_T = rawSupport_T / mult_T. Removed r-cliques drive the SAME vertex-level
// clean-split as CND -> the tree evolves identically -> cores are bit-identical.
// support-0 r-cliques (in no s-clique) get core 0 and are never peeled: their
// removal is net-zero on every support, so skipping them leaves support>0 cores
// unchanged.
// v1 (this) computes the per-leaf reprocess deltas by ENUMERATING r-cliques and
// grouping by tuple (correct, no speedup yet) -- it validates the peel/bucket/tree
// machinery. v2 replaces the enumeration with a combinatorial pivot-convolution.
// Requires g_m1ClassOf (PIVOTER_M1 class computation). Serial.
// ============================================================
std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionRCliqueTupleBatch(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    if (g_m1ClassOf.empty()) {
        std::cerr << "[tuple-batch] ERROR: g_m1ClassOf empty; set PIVOTER_M1_TUPLE_PROBE "
                     "(class computation) when running the tuple-batch peel.\n";
        return {};
    }

    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build", [&]() { cliqueIndex.build(tree, edgeGraph.adj_list.size()); });
    auto counting = daf::timeCount("countingPerRClique", [&]() {
        return CDSetRS::countingPerRCliqueParallel(tree, cliqueIndex, r, s);
    });
    const daf::Size nClique = cliqueIndex.size();

    // ---- build tuples over support>0 r-cliques (region scope) ----
    auto tupleKeyOf = [&](daf::Size id, std::vector<int> &cls, std::string &key) {
        auto span = cliqueIndex.byId(id);
        cls.clear();
        for (daf::Size v : span) cls.push_back(v < g_m1ClassOf.size() ? g_m1ClassOf[v] : -1);
        std::sort(cls.begin(), cls.end());
        key.clear();
        key.reserve(cls.size() * sizeof(int));
        for (int x : cls) key.append((const char *)&x, sizeof(int));
    };
    std::vector<double> coreRClique(nClique, 0.0);
    std::vector<int> tupleOf(nClique, -1);
    std::unordered_map<std::string, int> tupleOfKey;
    std::vector<long long> tupleMult;
    std::vector<double> rawSupport;
    // pass 1: assign each support>0 r-clique to its tuple, count mult, sum rawSupport
    {
        std::vector<int> cls;
        std::string key;
        for (daf::Size id = 0; id < nClique; ++id) {
            if (id >= counting.size() || counting[id] <= 0) continue; // support 0 -> core 0, never peeled
            tupleKeyOf(id, cls, key);
            auto it = tupleOfKey.find(key);
            int t;
            if (it == tupleOfKey.end()) {
                t = (int)tupleMult.size();
                tupleOfKey.emplace(std::move(key), t);
                tupleMult.push_back(0);
                rawSupport.push_back(0.0);
            } else t = it->second;
            tupleOf[id] = t;
            tupleMult[t]++;
            rawSupport[t] += counting[id];
        }
    }
    const int nTuples = (int)tupleMult.size();
    // tuple -> member r-cliques as a flat CSR (replaces vector<vector>, which cost
    // +12-30% RSS over CND from per-tuple heap chunks). M4-step-1.
    std::vector<daf::Size> memberOff(nTuples + 1, 0);
    for (int t = 0; t < nTuples; ++t) memberOff[t + 1] = memberOff[t] + (daf::Size)tupleMult[t];
    std::vector<daf::Size> memberFlat(memberOff[nTuples]);
    {
        std::vector<daf::Size> cursor(memberOff.begin(), memberOff.end() - 1);
        for (daf::Size id = 0; id < nClique; ++id) {
            int t = tupleOf[id];
            if (t >= 0) memberFlat[cursor[t]++] = id;
        }
    }

    std::vector<double> supportT(nTuples);
    int maxBucket = 0;
    long long sumMult = 0;
    for (int t = 0; t < nTuples; ++t) {
        supportT[t] = rawSupport[t] / (double)tupleMult[t];
        maxBucket = std::max(maxBucket, (int)supportT[t]);
        sumMult += tupleMult[t];
    }
    // Adaptive reprocess: the combinatorial per-leaf convolution pays off only when
    // r-cliques compress into far fewer tuples (high RS / symmetry); with little
    // compression it is pure overhead (measured 3x slower on dense ca-AstroPh).
    // avgMult = #support>0 r-cliques / #tuples = the batching benefit. Use the
    // combinatorial path iff avgMult >= threshold; else enumerate (which is already
    // <= CND). PIVOTER_TB_V1 forces enumeration; PIVOTER_TB_THRESHOLD overrides.
    const double avgMult = nTuples > 0 ? (double)sumMult / (double)nTuples : 1.0;
    double comboThreshold = 2.0;
    if (const char *th = std::getenv("PIVOTER_TB_THRESHOLD")) comboThreshold = atof(th);
    const bool forceEnum = std::getenv("PIVOTER_TB_V1") != nullptr;
    const bool forceCombo = std::getenv("PIVOTER_TB_COMBO") != nullptr;
    const bool useEnum = forceCombo ? false : (forceEnum || avgMult < comboThreshold);

    // ---- tuple buckets ----
    std::vector<std::vector<int> > buckets(maxBucket + 2);
    std::vector<int> bucket_of(nTuples), pos_in_bucket(nTuples);
    std::vector<char> tupleAlive(nTuples, 1);
    for (int t = 0; t < nTuples; ++t) {
        int b = (int)supportT[t];
        bucket_of[t] = b;
        pos_in_bucket[t] = (int)buckets[b].size();
        buckets[b].push_back(t);
    }
    auto bucketMove = [&](int t, int newB) {
        int oldB = bucket_of[t];
        if (newB == oldB) return;
        auto &oldVec = buckets[oldB];
        int myPos = pos_in_bucket[t];
        if (myPos < (int)oldVec.size() - 1) {
            int last = oldVec.back();
            oldVec[myPos] = last;
            pos_in_bucket[last] = myPos;
        }
        oldVec.pop_back();
        if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
        bucket_of[t] = newB;
        pos_in_bucket[t] = (int)buckets[newB].size();
        buckets[newB].push_back(t);
    };

    std::unordered_map<int, double> tupleDelta; // affected tuple -> raw support delta this round

    // ---- combinatorial per-leaf reprocess (v2; the speedup) ----
    // For a (sub)leaf L, A_L(T) = sum_{p<=needPivot} W_p(T) * nCr[pivotC-p][needPivot-p]
    // where W_p = [x^p] prod_{class c in T} ( sum_j C(kp_c,j) C(kk_c, m_c - j) x^j ),
    // needPivot = s - keepC_L. This sums CND's per-r-clique nCr over T's r-cliques in
    // L without enumerating each one -- enumerate the DISTINCT class-multisets (the
    // M2-reduced count) instead. accumulate sign*A into tupleDelta. Over-pivot terms
    // (p>needPivot) contribute 0 (guarded), matching v1/REF.
    constexpr int RMAX = 16;
    std::vector<int> rl_classId, rl_kp, rl_kk, rl_mult, rl_keyExp;
    std::string rl_key;
    std::function<void(int, int, int, int)> rl_rec; // (idx, remaining, pivotC, needPivot)
    auto reprocessLeaf = [&](const auto &leafVerts, double sign) {
        rl_classId.clear(); rl_kp.clear(); rl_kk.clear();
        static thread_local std::unordered_map<int, int> classIdx;
        classIdx.clear();
        int pivotC = 0, keepC = 0;
        for (const auto &node : leafVerts) {
            if (node.isPivot) pivotC++; else keepC++;
            int c = (node.v < g_m1ClassOf.size()) ? g_m1ClassOf[node.v] : -1;
            if (c < 0) continue; // support-0 vertex: counts in pivotC/keepC, not a tuple class
            auto it = classIdx.find(c);
            int ci;
            if (it == classIdx.end()) {
                ci = (int)rl_classId.size();
                classIdx.emplace(c, ci);
                rl_classId.push_back(c); rl_kp.push_back(0); rl_kk.push_back(0);
            } else ci = it->second;
            if (node.isPivot) rl_kp[ci]++; else rl_kk[ci]++;
        }
        int needPivot = (int)s - keepC;
        int K = (int)rl_classId.size();
        rl_mult.assign(K, 0);
        rl_rec = [&](int idx, int remaining, int pC, int nP) {
            if (remaining == 0) {
                rl_keyExp.clear();
                for (int i = 0; i < K; ++i) for (int m = 0; m < rl_mult[i]; ++m) rl_keyExp.push_back(rl_classId[i]);
                std::sort(rl_keyExp.begin(), rl_keyExp.end());
                rl_key.clear(); rl_key.reserve(rl_keyExp.size() * sizeof(int));
                for (int x : rl_keyExp) rl_key.append((const char *)&x, sizeof(int));
                auto it = tupleOfKey.find(rl_key);
                if (it == tupleOfKey.end()) return;        // all members support-0
                int t = it->second;
                if (!tupleAlive[t]) return;
                // W_p convolution over classes with mult>0
                std::array<double, RMAX> W{}; W[0] = 1.0; int deg = 0;
                for (int i = 0; i < K; ++i) {
                    int m = rl_mult[i]; if (m == 0) continue;
                    std::array<double, RMAX> P{};
                    for (int j = 0; j <= m; ++j) P[j] = nCr[rl_kp[i]][j] * nCr[rl_kk[i]][m - j];
                    std::array<double, RMAX> NW{};
                    for (int a = 0; a <= deg; ++a) if (W[a] != 0.0)
                        for (int j = 0; j <= m; ++j) NW[a + j] += W[a] * P[j];
                    W = NW; deg += m;
                }
                double A = 0.0;
                int pmax = std::min(nP, deg);
                for (int p = 0; p <= pmax; ++p) {
                    if (W[p] == 0.0) continue;
                    A += W[p] * nCr[pC - p][nP - p];
                }
                if (A != 0.0) tupleDelta[t] += sign * A;
                return;
            }
            if (idx == K) return;
            int cap = std::min(rl_kp[idx] + rl_kk[idx], remaining);
            for (int m = 0; m <= cap; ++m) { rl_mult[idx] = m; rl_rec(idx + 1, remaining - m, pC, nP); }
            rl_mult[idx] = 0;
        };
        rl_rec(0, (int)r, pivotC, needPivot);
    };

    // ---- peel ----
    auto tb_peel_start = std::chrono::high_resolution_clock::now();
    long long remaining = nTuples;
    int curBucket = 0;
    double minCore = 0;
    const daf::Size graphN = edgeGraph.n;
    static thread_local daf::StaticVector<daf::Size> threadMap;
    if (threadMap.maxSize < graphN + 1) threadMap.resize(graphN + 1);

    std::vector<daf::Size> changedLeaf;
    std::vector<std::vector<daf::Size> > removedRForLeaf;
    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<daf::Size> removedRCliques;

    while (remaining > 0) {
        for (auto leafId : changedLeaf) changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        changedLeaf.clear();
        removedRForLeaf.clear();
        removedRCliques.clear();
        tupleDelta.clear();

        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;
        minCore = std::max((double)curBucket, minCore);

        // pop all tuples at buckets [curBucket .. minCore]; assign core to members
        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty() && (double)curBucket <= minCore) {
            while (!buckets[curBucket].empty()) {
                int t = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                tupleAlive[t] = 0;
                remaining--;
                for (daf::Size i = memberOff[t]; i < memberOff[t + 1]; ++i) {
                    daf::Size c = memberFlat[i];
                    coreRClique[c] = minCore;
                    removedRCliques.push_back(c);
                }
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (double)(curBucket + 1) <= minCore)
                curBucket++;
            else break;
        }

        if (removedRCliques.empty()) {
            if (remaining == 0) break;
            continue;
        }

        // find affected leaves (same as CND: intersect removed members' vertices via treeGraphV)
        for (daf::Size rmId : removedRCliques) {
            auto rClique = cliqueIndex.byId(rmId);
            daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                [&](const TreeGraphNode &uClique) {
                    auto &leafSlot = changedLeafIndex[uClique.v];
                    if (leafSlot == std::numeric_limits<daf::Size>::max()) {
                        leafSlot = removedRForLeaf.size();
                        removedRForLeaf.emplace_back();
                        changedLeaf.push_back(uClique.v);
                    }
                    removedRForLeaf[leafSlot].push_back(rmId);
                });
        }

        // per affected leaf: clean-split (vertex-level, identical to CND) + tuple deltas
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            daf::Size leafId = changedLeaf[idx];
            const auto &leaf = tree.adj_list[leafId];
            daf::Size leafSlot = changedLeafIndex[leafId];
            const auto &removedR = removedRForLeaf[leafSlot];

            daf::CliqueSize pivotC = 0, keepC = 0;
            for (const auto &node : leaf) { if (node.isPivot) pivotC++; else keepC++; }

            for (const auto &leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }

            auto mapped = removedR | std::views::transform([&](const daf::Size id) { return cliqueIndex.byId(id); });
            bkRmClique::removeRClique(leaf, mapped, r, s,
                [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                    auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                    auto newId = tree.addNode(newLeaf);
                    const auto &stored = tree.adj_list[newId];
                    daf::CliqueSize newPivotC = 0, newKeepC = 0;
                    for (const auto &i : stored) {
                        if (i.isPivot) { treeGraphV.addNbr(i.v, {newId, true}); newPivotC++; }
                        else { treeGraphV.addNbr(i.v, {newId, false}); newKeepC++; }
                    }
                    if (useEnum) {
                        daf::Size needPivot = s - newKeepC;
                        daf::enumerateCombinations(stored, r, [&](const daf::StaticVector<TreeGraphNode> &rclique) {
                            daf::CliqueSize sub = 0;
                            for (const auto &node : rclique) if (node.isPivot) sub++;
                            if (sub <= needPivot) {
                                auto id = cliqueIndex.byClique(rclique);
                                int t = (id < (daf::Size)tupleOf.size()) ? tupleOf[id] : -1;
                                if (t >= 0 && tupleAlive[t]) tupleDelta[t] += nCr[newPivotC - sub][needPivot - sub];
                            }
                            return true;
                        });
                    } else {
                        reprocessLeaf(stored, +1.0);
                    }
                    if (newId >= changedLeafIndex.size())
                        changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                },
                &threadMap);

            if (useEnum) {
                daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &clique) {
                    auto id = cliqueIndex.byClique(clique);
                    int t = (id < (daf::Size)tupleOf.size()) ? tupleOf[id] : -1;
                    if (t < 0 || !tupleAlive[t]) return true;
                    daf::CliqueSize sub = 0;
                    for (const auto &node : clique) if (node.isPivot) sub++;
                    tupleDelta[t] -= nCr[pivotC - sub][s - keepC - sub];
                    return true;
                });
            } else {
                reprocessLeaf(leaf, -1.0);
            }

            tree.removeNode(leafId);
        }

        // apply tuple deltas, move buckets, pull curBucket back
        for (auto &kv : tupleDelta) {
            int t = kv.first;
            if (!tupleAlive[t]) continue;
            rawSupport[t] += kv.second;
            int newB = std::max(0, (int)(rawSupport[t] / (double)tupleMult[t]));
            if (newB != bucket_of[t]) {
                bucketMove(t, newB);
                if (newB < curBucket) curBucket = newB;
            }
        }
    }

    {
        auto tb_peel_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - tb_peel_start).count();
        printf("[tuple-batch] peel: %lld ms  (tuples=%d, r-cliques=%lld, mode=%s)\n",
               (long long)tb_peel_ms, nTuples, (long long)nClique, useEnum ? "v1-enum" : "v2-combinatorial");
        fflush(stdout);
    }

    // ---- assemble per-r-clique result (same format as CND) ----
    std::vector<std::pair<std::vector<daf::Size>, double> > sortedK;
    sortedK.reserve(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueCopy(clique.begin(), clique.end());
        sortedK.emplace_back(cliqueCopy, coreRClique[i]);
    }
    std::sort(sortedK.begin(), sortedK.end(), [](const auto &a, const auto &b) { return a.second < b.second; });
    return sortedK;
}


// ============================================================
// §105 M4: tuple-NATIVE r-clique nucleus decomposition (PIVOTER_RUN_TUPLE_NATIVE).
// The cliqueIndex-FREE engine: never materialises the per-r-clique index, so peak
// memory is O(#tuples + tuple-leaf incidences + tree) -- below CND (whose #r-clique
// CPI explodes to 94GB on com-dblp 5,6). It also batches the REMOVAL-side leaf
// finding via a maintained tuple->leaves index (L1), not per-r-clique intersect.
//   counting: rawSupport_T = sum_leaves A_L(T) (the validated convolution).
//   mult_T   = prod_c C(classSize[c], m_c) (region_native's formula; same-class
//              vertices are mutually adjacent, comp classes share a host region).
//   support_T = rawSupport_T / mult_T  (the bucket key).
//   removal:  affected leaves from tupleLeaves[T]; enumerate the tuple's concrete
//             members IN each leaf on the fly to feed bkRmClique::removeRClique.
//   reprocess: combinatorial A_L(T) deltas; maintain tupleLeaves on split.
//   output:   per-tuple core -> distribution (Sum mult), validated vs region_native.
// Requires g_m1ClassOf. Serial. Combinatorial regime (high RS) is the target.
// ============================================================
std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionRCliqueTupleNative(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    if (g_m1ClassOf.empty()) {
        std::cerr << "[tuple-native] ERROR: g_m1ClassOf empty; set the class computation gate.\n";
        return {};
    }
    const std::vector<int> &classOf = g_m1ClassOf;
    int nClasses = 0;
    for (int c : classOf) nClasses = std::max(nClasses, c + 1);
    std::vector<long long> classSize(nClasses, 0);
    for (int c : classOf) if (c >= 0) classSize[c]++;

    auto binomD = [](long long n, long long k) -> double {
        if (k < 0 || k > n) return 0.0;
        if (k > n - k) k = n - k;
        double res = 1.0;
        for (long long i = 0; i < k; ++i) res = res * (double)(n - i) / (double)(i + 1);
        return res;
    };

    const daf::Size graphN = edgeGraph.n;
    static thread_local daf::StaticVector<daf::Size> threadMap;
    if (threadMap.maxSize < graphN + 1) threadMap.resize(graphN + 1);

    constexpr int RMAX = 16;
    // fixed-size tuple key (sorted class-multiset, exactly r entries + padding) -- no
    // string allocation, no hashing of a heap buffer.
    struct ArrKey {
        std::array<int, RMAX> a;
        bool operator==(const ArrKey &o) const { return a == o.a; }
    };
    struct ArrHash {
        size_t operator()(const ArrKey &k) const {
            size_t h = 1469598103934665603ULL;
            for (int x : k.a) { h ^= (size_t)(unsigned)x; h *= 1099511628211ULL; }
            return h;
        }
    };
    // ---- tuple registry ----
    std::unordered_map<ArrKey, int, ArrHash> tupleOfKey;
    std::vector<double> rawSupport;
    std::vector<long long> tupleMult;                       // prod C(classSize,m)
    std::vector<std::vector<std::pair<int, int> > > tupleComp; // (classId, mult), sorted by classId
    // robin_hood flat set: ~13B/incidence (open addressing) vs std::unordered_set's
    // ~40B (node + pointers) -> ~3x leaner on the high-RS inverted index (181M
    // incidences = 11.4GB -> ~3.8GB on com-dblp 5,6). Same insert/erase/iterate API.
    std::vector<robin_hood::unordered_flat_set<daf::Size> > tupleLeaves;
    std::vector<int> classToLocal(nClasses, -1);            // flat class->local-index (reset per leaf via pl_cid)

    // process one (sub)leaf: enumerate its class-multisets, compute A_L(T), and either
    // (counting) create/update tuples + register tupleLeaves, or (reprocess) update an
    // existing alive tuple's delta + maintain tupleLeaves. Only A_L(T)>0 is registered.
    std::vector<int> pl_cid, pl_kp, pl_kk, pl_mult, pl_keyExp;
    auto processLeaf = [&](const auto &leafVerts, daf::Size leafId, double sign, bool counting,
                           std::vector<char> *alivePtr, std::unordered_map<int, double> *deltaPtr) {
        pl_cid.clear(); pl_kp.clear(); pl_kk.clear();
        int pivotC = 0, keepC = 0;
        for (const auto &node : leafVerts) {
            if (node.isPivot) pivotC++; else keepC++;
            int c = (node.v < classOf.size()) ? classOf[node.v] : -1;
            if (c < 0) continue;
            int ci = classToLocal[c];
            if (ci < 0) { ci = (int)pl_cid.size(); classToLocal[c] = ci; pl_cid.push_back(c); pl_kp.push_back(0); pl_kk.push_back(0); }
            if (node.isPivot) pl_kp[ci]++; else pl_kk[ci]++;
        }
        int needPivot = (int)s - keepC;
        int K = (int)pl_cid.size();
        pl_mult.assign(K, 0);
        auto pl_rec = [&](auto &&self, int idx, int remaining, int pC, int nP) -> void {
            if (remaining == 0) {
                // A_L(T) convolution
                std::array<double, RMAX> W{}; W[0] = 1.0; int deg = 0;
                for (int i = 0; i < K; ++i) {
                    int m = pl_mult[i]; if (m == 0) continue;
                    std::array<double, RMAX> P{};
                    for (int j = 0; j <= m; ++j) P[j] = nCr[pl_kp[i]][j] * nCr[pl_kk[i]][m - j];
                    std::array<double, RMAX> NW{};
                    for (int a = 0; a <= deg; ++a) if (W[a] != 0.0)
                        for (int j = 0; j <= m; ++j) NW[a + j] += W[a] * P[j];
                    W = NW; deg += m;
                }
                double A = 0.0;
                int pmax = std::min(nP, deg);
                for (int p = 0; p <= pmax; ++p) if (W[p] != 0.0) A += W[p] * nCr[pC - p][nP - p];
                if (A == 0.0) return;                       // only A>0 leaves are registered
                // build canonical key (expanded sorted class list)
                pl_keyExp.clear();
                for (int i = 0; i < K; ++i) for (int m = 0; m < pl_mult[i]; ++m) pl_keyExp.push_back(pl_cid[i]);
                std::sort(pl_keyExp.begin(), pl_keyExp.end());
                ArrKey pl_key;
                pl_key.a.fill(-1);
                for (int i = 0; i < (int)pl_keyExp.size(); ++i) pl_key.a[i] = pl_keyExp[i];
                if (counting) {
                    auto it = tupleOfKey.find(pl_key);
                    int t;
                    if (it == tupleOfKey.end()) {
                        t = (int)tupleMult.size();
                        tupleOfKey.emplace(pl_key, t);
                        rawSupport.push_back(0.0);
                        tupleLeaves.emplace_back();
                        // comp = run-length of pl_keyExp; mult = prod C(classSize,m)
                        std::vector<std::pair<int, int> > comp;
                        long long mu = 1;
                        for (int a = 0; a < (int)pl_keyExp.size();) {
                            int c = pl_keyExp[a], b = a;
                            while (b < (int)pl_keyExp.size() && pl_keyExp[b] == c) b++;
                            int m = b - a;
                            comp.emplace_back(c, m);
                            mu *= (long long)llround(binomD(classSize[c], m));
                            a = b;
                        }
                        tupleComp.push_back(std::move(comp));
                        tupleMult.push_back(mu);
                    } else t = it->second;
                    rawSupport[t] += A;
                    tupleLeaves[t].insert(leafId);
                } else {
                    auto it = tupleOfKey.find(pl_key);
                    if (it == tupleOfKey.end()) return;     // not a support>0 tuple
                    int t = it->second;
                    if (!(*alivePtr)[t]) return;
                    (*deltaPtr)[t] += sign * A;
                    if (sign > 0) tupleLeaves[t].insert(leafId);
                    else tupleLeaves[t].erase(leafId);
                }
                return;
            }
            if (idx == K) return;
            int cap = std::min(pl_kp[idx] + pl_kk[idx], remaining);
            for (int m = 0; m <= cap; ++m) { pl_mult[idx] = m; self(self, idx + 1, remaining - m, pC, nP); }
            pl_mult[idx] = 0;
        };
        pl_rec(pl_rec, 0, (int)r, pivotC, needPivot);
        for (int c : pl_cid) classToLocal[c] = -1;   // reset only touched entries
    };

    // ---- counting: build tuples + rawSupport + tupleLeaves from the initial tree ----
    daf::timeCount("tuple-native count", [&]() {
        for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
            const auto &leaf = tree.adj_list[leafId];
            if (leaf.size() < r) continue;
            processLeaf(leaf, leafId, +1.0, true, nullptr, nullptr);
        }
    });
    const int nTuples = (int)tupleMult.size();
    if (std::getenv("PIVOTER_TN_PROFILE")) {
        long long inc = 0;
        for (const auto &ls : tupleLeaves) inc += (long long)ls.size();
        printf("[tn-mem] after-count: tuples=%d tupleLeaves_incidences=%lld\n", nTuples, inc);
        daf::log_memory("tn after-count");
        fflush(stdout);
    }

    std::vector<double> coreTuple(nTuples, 0.0);
    std::vector<double> supportT(nTuples);
    std::vector<char> tupleAlive(nTuples, 1);
    int maxBucket = 0;
    for (int t = 0; t < nTuples; ++t) {
        supportT[t] = rawSupport[t] / (double)tupleMult[t];
        maxBucket = std::max(maxBucket, (int)supportT[t]);
    }

    // ---- tuple buckets ----
    std::vector<std::vector<int> > buckets(maxBucket + 2);
    std::vector<int> bucket_of(nTuples), pos_in_bucket(nTuples);
    for (int t = 0; t < nTuples; ++t) {
        int b = (int)supportT[t];
        bucket_of[t] = b; pos_in_bucket[t] = (int)buckets[b].size();
        buckets[b].push_back(t);
    }
    auto bucketMove = [&](int t, int newB) {
        int oldB = bucket_of[t]; if (newB == oldB) return;
        auto &ov = buckets[oldB]; int mp = pos_in_bucket[t];
        if (mp < (int)ov.size() - 1) { int last = ov.back(); ov[mp] = last; pos_in_bucket[last] = mp; }
        ov.pop_back();
        if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
        bucket_of[t] = newB; pos_in_bucket[t] = (int)buckets[newB].size(); buckets[newB].push_back(t);
    };

    // Allocation-free member enumeration (was the 80%-of-peel hot spot): reused
    // scratch instead of a per-leaf unordered_map + a per-member heap vector.
    // Members are appended as a flat CSR (cflat = concatenated sorted vertex sets,
    // coff = offsets). compVerts groups the leaf's vertices by the tuple's classes.
    std::vector<daf::Size> cflat;                       // all members' vertices, concatenated
    std::vector<size_t> coff;                           // CSR offsets; coff.size()-1 = #members
    std::vector<std::vector<daf::Size> > compVerts(std::max<int>(1, (int)r)); // per comp-class leaf verts
    std::vector<daf::Size> emCur;                       // current partial member
    auto enumMembers = [&](const std::vector<TreeGraphNode> &leaf, const std::vector<std::pair<int, int> > &comp) {
        const int K = (int)comp.size();
        for (int i = 0; i < K; ++i) compVerts[i].clear();
        for (const auto &node : leaf) {
            int c = classOf[node.v];
            if (c < 0) continue;
            for (int i = 0; i < K; ++i) if (comp[i].first == c) { compVerts[i].push_back(node.v); break; }
        }
        emCur.clear();
        auto em_rec = [&](auto &&self_em, int ci) -> void {
            if (ci == K) {
                size_t base = cflat.size();
                cflat.insert(cflat.end(), emCur.begin(), emCur.end());
                std::sort(cflat.begin() + base, cflat.end());
                coff.push_back(cflat.size());
                return;
            }
            const auto &verts = compVerts[ci];
            int m = comp[ci].second, nv = (int)verts.size();
            auto choose = [&](auto &&self_ch, int start, int k) -> void {
                if (k == 0) { self_em(self_em, ci + 1); return; }
                for (int i = start; i + k <= nv; ++i) { emCur.push_back(verts[i]); self_ch(self_ch, i + 1, k - 1); emCur.pop_back(); }
            };
            choose(choose, 0, m);
        };
        em_rec(em_rec, 0);
    };
    // lightweight indexable view over (cflat, coff) for bkRmClique::removeRClique
    std::vector<std::span<const daf::Size> > conflictSpans; // reused range view for removeRClique

    // ---- peel ----
    auto tn_peel_start = std::chrono::high_resolution_clock::now();
    const bool tnprof = std::getenv("PIVOTER_TN_PROFILE") != nullptr;
    long long pT_find = 0, pT_enum = 0, pT_decr = 0, pT_split = 0, pT_apply = 0;
    auto pNow = [] { return std::chrono::high_resolution_clock::now(); };
    auto pNs = [](auto a, auto b) { return std::chrono::duration_cast<std::chrono::nanoseconds>(b - a).count(); };
    long long remaining = nTuples;
    int curBucket = 0;
    double minCore = 0;

    std::vector<daf::Size> changedLeaf;
    std::vector<std::vector<int> > perLeafPopped;
    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::unordered_map<int, double> tupleDelta;

    while (remaining > 0) {
        for (auto leafId : changedLeaf) changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        changedLeaf.clear();
        perLeafPopped.clear();
        tupleDelta.clear();

        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;
        minCore = std::max((double)curBucket, minCore);

        std::vector<int> popped;
        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty() && (double)curBucket <= minCore) {
            while (!buckets[curBucket].empty()) {
                int t = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                tupleAlive[t] = 0;
                coreTuple[t] = minCore;
                remaining--;
                popped.push_back(t);
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (double)(curBucket + 1) <= minCore)
                curBucket++;
            else break;
        }
        if (popped.empty()) { if (remaining == 0) break; continue; }

        // collect affected leaves per popped tuple (from tupleLeaves), then free
        auto _tf = pNow();
        for (int t : popped) {
            for (daf::Size L : tupleLeaves[t]) {
                auto &slot = changedLeafIndex[L];
                if (slot == std::numeric_limits<daf::Size>::max()) {
                    slot = perLeafPopped.size();
                    perLeafPopped.emplace_back();
                    changedLeaf.push_back(L);
                }
                perLeafPopped[slot].push_back(t);
            }
            tupleLeaves[t] = robin_hood::unordered_flat_set<daf::Size>{}; // free
        }
        pT_find += pNs(_tf, pNow());

        // per affected leaf: enumerate popped members, decr, clean-split (incr), mutate tree
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            daf::Size leafId = changedLeaf[idx];
            const auto &leaf = tree.adj_list[leafId];
            daf::Size slot = changedLeafIndex[leafId];

            auto _te = pNow();
            cflat.clear();
            coff.assign(1, 0);
            for (int t : perLeafPopped[slot]) enumMembers(leaf, tupleComp[t]);
            conflictSpans.clear();
            for (size_t i = 0; i + 1 < coff.size(); ++i)
                conflictSpans.emplace_back(cflat.data() + coff[i], coff[i + 1] - coff[i]);
            auto &conflictSets = conflictSpans;
            pT_enum += pNs(_te, pNow());

            // decr on the old leaf (erase L from alive tuples' tupleLeaves, accumulate -A)
            auto _td = pNow();
            processLeaf(leaf, leafId, -1.0, false, &tupleAlive, &tupleDelta);
            pT_decr += pNs(_td, pNow());

            auto _ts = pNow();
            // remove old leaf's vertices from treeGraphV
            for (const auto &lv : leaf) {
                if (lv.isPivot) treeGraphV.removeNbr(lv.v, {leafId, true});
                else treeGraphV.removeNbr(lv.v, {leafId, false});
            }

            // clean-split (vertex-level); incr on each new sub-leaf
            bkRmClique::removeRClique(leaf, conflictSets, (int)r, (int)s,
                [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                    auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                    auto newId = tree.addNode(newLeaf);
                    const auto &stored = tree.adj_list[newId];
                    for (const auto &i : stored) {
                        if (i.isPivot) treeGraphV.addNbr(i.v, {newId, true});
                        else treeGraphV.addNbr(i.v, {newId, false});
                    }
                    if (stored.size() >= r) processLeaf(stored, newId, +1.0, false, &tupleAlive, &tupleDelta);
                    if (newId >= changedLeafIndex.size())
                        changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                },
                &threadMap);

            tree.removeNode(leafId);
            pT_split += pNs(_ts, pNow());
        }

        // apply tuple deltas, move buckets, pull curBucket back
        auto _ta = pNow();
        for (auto &kv : tupleDelta) {
            int t = kv.first;
            if (!tupleAlive[t]) continue;
            rawSupport[t] += kv.second;
            int newB = std::max(0, (int)(rawSupport[t] / (double)tupleMult[t]));
            if (newB != bucket_of[t]) { bucketMove(t, newB); if (newB < curBucket) curBucket = newB; }
        }
        pT_apply += pNs(_ta, pNow());
    }

    {
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - tn_peel_start).count();
        printf("[tuple-native] peel: %lld ms  (tuples=%d)\n", (long long)ms, nTuples);
        if (tnprof)
            printf("[tn-profile] find=%.0fms enumMembers=%.0fms decr=%.0fms split+incr=%.0fms apply=%.0fms\n",
                   pT_find / 1e6, pT_enum / 1e6, pT_decr / 1e6, pT_split / 1e6, pT_apply / 1e6);
    }

    // ---- output: core distribution weighted by mult (matches region_native) ----
    std::map<double, double> dist;
    for (int t = 0; t < nTuples; ++t) dist[coreTuple[t]] += (double)tupleMult[t];
    for (auto &kv : dist) printf("[tuple-native] core=%.0f count=%.0f\n", kv.first, kv.second);
    fflush(stdout);

    // return empty (this engine validates via the printed distribution, not per-r-clique)
    return {};
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