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

extern double nCr[1001][401];
// （），
// #ifndef NDEBUG
// set NOEBUG as trus


namespace CDSetRS {
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


    struct CompareRClique {
        const double *RCliqueCounting; // 
        explicit CompareRClique(const double *coreLeaf) : RCliqueCounting(coreLeaf) {
        }

        // ： “a ” ， coreLeaf[a] > coreLeaf[b]
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

    std::vector<double> countingPerRClique(
        const DynamicGraph<TreeGraphNode> &treeGraph,
        StaticCliqueIndex &cliqueHashmap,
        const daf::CliqueSize r,
        const daf::CliqueSize s) {
        // double *rCliqueSCounting = new double[cliqueHashmap.size()];
        // memset(rCliqueSCounting, 0, cliqueHashmap.size() * sizeof(double));
        std::vector<double> rCliqueSCounting(cliqueHashmap.size(), 0.0);
        daf::Size count = 0;
        for (const auto &leaf: treeGraph.adj_list) {
            if (leaf.size() < r) {
                continue;
            }
            daf::CliqueSize pivotC = 0, keepC = 0;
            for (auto &i: leaf) {
                if (i.isPivot) pivotC++;
                else keepC++;
            }
            // std::cout << "leaf: " << leaf << " pivotC: " << pivotC
            //           << " keepC: " << keepC << " size: " << leaf.size() << std::endl;
            // DEBUG_BREAK_IF(leaf.size() == 4);
            // clique
            int needPivot = s - static_cast<int>(keepC);
            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            // daf::enumerateCombinationsIdx(leaf, r, [&](const std::size_t* idx) {
                daf::CliqueSize subNumKeepC = 0;
                daf::CliqueSize subNumPovit = 0;
                for (const auto &node: rClique) {
                    if (node.isPivot) {
                        subNumPovit++;
                    } else {
                        subNumKeepC++;
                    }
                }

                auto ncrValue = nCr[pivotC - subNumPovit][needPivot - subNumPovit];
                // rCliqueSCounting[cliqueHashmap.byNewClique(rClique)] += ncrValue;
                auto [id, isNew] = cliqueHashmap.byNewClique(rClique);
                if (isNew) {
                    if (rCliqueSCounting.size() <= id) {
                        rCliqueSCounting.push_back(0.0);
                    }
                    if (rCliqueSCounting.capacity() <= id) {
                        rCliqueSCounting.reserve(std::max(id + 2, id * 2));
                    }
                }
                rCliqueSCounting[id] += ncrValue;
                return true;
            });
        }
        rCliqueSCounting.shrink_to_fit();
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
    void printEdgeCore(const Graph &edgeGraph, const std::vector<T> coreE) {
        printEdgeCore(edgeGraph, coreE.data());
    }
}


std::vector<std::pair<std::vector<daf::Size>, int> > NucleusCoreDecompositionRClique(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    auto time_start = std::chrono::high_resolution_clock::now();
    // Accumulate total time (nanoseconds) spent updating countingRClique in this function
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build",
                   [&]() {
                       cliqueIndex.build(tree, edgeGraph.adj_list.size());
                   });

    daf::log_memory("r-clique index");
    // tree.printGraphPerV();
    // cliqueIndex.print();
    // cliqueIndex.verify();


    auto countingRClique = daf::timeCount("countingPerRClique",
                                          [&]() {
                                              return CDSetRS::countingPerRClique(
                                                  tree, cliqueIndex, r, s);
                                          });

    std::vector<daf::Size> coreRClique(countingRClique.size());
#ifndef NDEBUG
    tree.printGraphPerV();
    // daf::printArray(countingKE, edgeGraph.adj_list.size());
    // CDSetRS::printEdgeCore(edgeGraph, countingKE);
    // std::cout << "countingRClique" << countingRClique << std::endl;
    for (daf::Size i = 0; i < countingRClique.size(); ++i) {
        std::cout << "Clique: " << cliqueIndex.byId(i) << " id: " << i
                  << " count: " << countingRClique[i] << std::endl;
    }
    std::cout << "cliqueIndex Size : " << cliqueIndex.size() << std::endl;
    // CDSetRS::printEdgeCore(edgeGraph, degreeE);
#endif

    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<std::vector<daf::Size> > removedRCliqueIdForLeaf;
    std::vector<daf::Size> changedLeaf;
    std::vector<daf::Size> currentRemoveRcliqueIds;

    removedRCliqueIdForLeaf.reserve(tree.adj_list.size());
    changedLeaf.reserve(tree.adj_list.size());
    currentRemoveRcliqueIds.reserve(cliqueIndex.size());


    daf::StaticVector<bool> rCliqueInHeap(cliqueIndex.size());
    // rCliqueInHeap.fill(true);
    rCliqueInHeap.resize(cliqueIndex.size());
    memset(rCliqueInHeap.getData(), true, cliqueIndex.size() * sizeof(bool));

    CDSetRS::DHeap heap{CDSetRS::CompareRClique(countingRClique.data())};
    heap.reserve(cliqueIndex.size());

    std::vector<CDSetRS::DHeap::handle_type> heapHandles(cliqueIndex.size());

    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
        heapHandles[i] = heap.push(i);
    }
#ifndef NDEBUG
    std::cout << "countingKE: ";
    // CDSetRS::printEdgeCore(edgeGraph, countingKE);

    std::cout << "countingRClique: " << countingRClique << std::endl;

    std::cout << "tree: ";
    tree.printGraphPerV();

    std::cout << "treeGraphV: ";
    treeGraphV.printGraphPerV();
#endif

    daf::log_memory("Other index(incloud head)");
    std::cout << "=========================begin=========================" << std::endl;
    double minCore = 0;
    while (!heap.empty()) {
        for (auto &leafId: changedLeaf) {
            changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        }
        changedLeaf.clear();
        removedRCliqueIdForLeaf.clear();
        currentRemoveRcliqueIds.clear();

        minCore = std::max(countingRClique[heap.top()], minCore);
        //  core==minCore  leaf  pop 
        std::cout << "minCore: " << minCore
        << " heap size: " << heap.size()
        << " num Leaf: " << tree.size() << " "
        << s << "-Clique count: " << tree.cliqueCount(s)
        << std::endl;

        daf::log_memory("while loop top");
        // if (minCore == 99) break;
        while (!heap.empty() && countingRClique[heap.top()] <= minCore) {
            auto id = heap.top();
            rCliqueInHeap[id] = false;
            heap.pop();
            currentRemoveRcliqueIds.push_back(id);
            coreRClique[id] = minCore;
#ifndef NDEBUG
            std::cout << "removed Clique: " << cliqueIndex.byId(id) << " id: " << id
                      << " core: " << countingRClique[id] << std::endl;
#endif
        }

        if (heap.empty()) {
            break;
        }

        for (auto rmRCliqueId: currentRemoveRcliqueIds) {
            auto rClique = cliqueIndex.byId(rmRCliqueId);
            // std::cout << "rClique: " << rClique << std::endl;
            daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                                            [&](const TreeGraphNode &uClique) {
                                                // std::cout << "uClique: " << uClique << std::endl;
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
        // std::cout << "changedLeaf: " << changedLeaf << std::endl;
        // std::cout << "changedLeafIndex: " << changedLeafIndex << std::endl;
        for (auto leafId: changedLeaf) {
            auto leaf = tree.adj_list[leafId];
            auto leafIndex = changedLeafIndex[leafId];
            // std::cout << "============================================================" << std::endl;
            // std::cout << "changed leafId: " << leafId << " leaf index: " << leafIndex
            //           << " leaf: " << leaf << std::endl;
            // std::cout << "removedRCliqueIdForLeaf: ";
            // for (const auto &id: removedRCliqueIdForLeaf[leafIndex]) {
            //     std::cout << id << " (" << cliqueIndex.byId(id) << ") ";
            // }

            // std::cout << std::endl;

            auto initCore = [&](const std::vector<TreeGraphNode> &newLeaf, const daf::Size &newLeafId) {
                daf::CliqueSize newPivotC = 0, newKeepC = 0;
                for (auto i: newLeaf) {
                    if (i.isPivot) {
                        treeGraphV.addNbr(i.v, {newLeafId, true});
                        newPivotC++;
                    } else {
                        treeGraphV.addNbr(i.v, {newLeafId, false});
                        newKeepC++;
                    }
                }
                daf::Size needPivot = s - newKeepC;


                daf::enumerateCombinations(newLeaf, r, [&](const daf::StaticVector<TreeGraphNode> &rclique) {
                    daf::CliqueSize subNumKeepC = 0;
                    daf::CliqueSize subNumPovit = 0;
                    for (const auto &node: rclique) {
                        if (node.isPivot) {
                            subNumPovit++;
                        } else {
                            subNumKeepC++;
                        }
                    }

                    if (subNumPovit <= needPivot) {
                        if (newPivotC - subNumPovit < 0 || newPivotC - subNumPovit >= 1001 ||
                            needPivot - subNumPovit < 0 || needPivot - subNumPovit >= 401) {
                            // std::cerr << "nCr index out of range: row
                            std::cerr << "nCr index out of range: row=" << newPivotC - subNumPovit
                                    << " col=" << needPivot - subNumPovit
                                    << " newPivotC=" << newPivotC
                                    << " subNumPovit=" << subNumPovit
                                    << " needPivot=" << needPivot
                                    << " subNumKeepC=" << subNumKeepC
                                    << " newLeaf: " << rclique
                                    << std::endl;
                            for (auto &node: newLeaf) {
                                std::cout << "node: " << node.v << " isPivot: " << node.isPivot << std::endl;
                            }
                            // << "leaf: " << newLeaf

                            std::abort();
                        }
                        // timed update to countingRClique (measure whole enumerateCombinations call by outer wrapper)
                        auto ncrValue = nCr[newPivotC - subNumPovit][needPivot - subNumPovit];
                        auto cliqueIndexId = cliqueIndex.byClique(rclique);
                        countingRClique[cliqueIndexId] += ncrValue;
                    }

                    return true;
                });
                // end of enumerateCombinations in initCore (timed by outer scope)
            };

            // if (!removedPovit.empty() && needPivot <= povit.size() - removedPovit.size())
            // daf::StaticVector<daf::Size> newLeafIds;
            for (auto leafV: leaf) {
                if (leafV.isPivot) {
                    treeGraphV.removeNbr(leafV.v, {leafId, true});
                } else {
                    treeGraphV.removeNbr(leafV.v, {leafId, false});
                }
            }

            auto removedR = removedRCliqueIdForLeaf[leafIndex];
            auto mapped = removedR | std::views::transform(
                              [&](const daf::Size id) {
                                  return cliqueIndex.byId(id);
                              }
                          );

            // DEBUG_BREAK_IF(leafId == 5);
            bkRmClique::removeRClique(leaf, mapped, r, s, [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                DEBUG_BREAK_IF(newLeaf.size() < s);
                // std::cout << " newLeaf: " << newLeaf << std::endl;
                auto newId = tree.addNode(newLeaf);
                // std::cout << " newId: " << newId << std::endl;
                initCore(tree.adj_list[newId], newId);
                if (newId >= changedLeafIndex.size()) {
                    changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
                }
            });

            daf::CliqueSize keepC = 0, pivotC = 0;
            for (const auto &node: leaf) {
                if (node.isPivot) pivotC++;
                else keepC++;
            }
            // time the enumeration that decrements countingRClique
            {
                auto __t1 = std::chrono::high_resolution_clock::now();
                daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &clique) {
                     auto cliqueIndexId = cliqueIndex.byClique(clique);
                     // std::cout << "cliqueIndexId: " << cliqueIndexId
                     //           << " clique: " << clique << std::endl;
                     // rCliqueInHeap.print("rCliqueInHeap: ");
                     if (!rCliqueInHeap[cliqueIndexId]) return true;
                     daf::CliqueSize subNumKeepC = 0;
                     daf::CliqueSize subNumPovit = 0;
                     for (const auto &node: clique) {
                         if (node.isPivot) subNumPovit++;
                         else subNumKeepC++;
                     }
                     auto ncrValue = nCr[pivotC - subNumPovit][s - keepC - subNumPovit];
                     countingRClique[cliqueIndexId] -= ncrValue;
                     // std::cout << "ncrValue: " << ncrValue
                     //           << " pivotC: " << pivotC
                     //           << " subNumPovit: " << subNumPovit
                     //           << " keepC: " << keepC
                     //           << " subNumKeepC: " << subNumKeepC
                     //           << std::endl;
                     heap.update(heapHandles[cliqueIndexId]);
                     return true;
                 });
            }
            // end of timed enumerateCombinations that decrements countingRClique

             tree.removeNode(leafId);
        }


#ifndef NDEBUG
            std::cout << "tree: ";
            tree.printGraphPerV();

            std::cout << "treeGraphV: ";
            treeGraphV.printGraphPerV();


            std::cout << "clique countingRClique: ";
            for (daf::Size i = 0; i < countingRClique.size(); ++i) {
                std::cout << i << ": " << countingRClique[i] << " " << cliqueIndex.byId(i) << std::endl;
            }
#endif
    }
    // currentRemoveLeafIds.clear();
#ifndef NDEBUG
        std::cout << "tree: ";
        tree.printGraphPerV();

        std::cout << "treeGraphV: ";
        treeGraphV.printGraphPerV();
#endif


    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;


    // coreE
    // daf::printArray(coreE, edgeGraph.adj_list.size());


    // ~/_/pivoter/a
    // std::sort(coreE, coreE + edgeGraph.adj_list.size());
    std::vector<std::pair<std::vector<daf::Size>, int> > sortedK;
    sortedK.reserve(countingRClique.size());

    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueCopy(clique.begin(), clique.end());
        sortedK.emplace_back(cliqueCopy, coreRClique[i]);
    }
    std::sort(sortedK.begin(), sortedK.end(),
              [](const auto &a, const auto &b) {
                  return a.second < b.second; //  core 
              });
    // auto file = fopen("~/_/pivoter/a.out", "w");
    // for (auto i: sortedK) {
    //     printf("[");
    //     for (std::size_t j = 0; j < i.first.size(); ++j) {
    //         if (j != 0) printf(",");
    //         printf("%zu", i.first[j]);
    //     }
    //     printf("] %d\n", (int) i.second);
    //     // fprintf(file, "%d\n", (int) i.second);
    // }
    // fclose(file);
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