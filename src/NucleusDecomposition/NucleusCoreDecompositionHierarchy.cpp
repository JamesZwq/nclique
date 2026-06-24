//
// Created by _ on 25-3-4.
//

#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include <boost/heap/d_ary_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <set>
#include <fstream>
#include <filesystem>
#include <vector>
#include <map>
#include <limits>
#include <cstdlib>

#include "../BK/BronKerboschRmEdge.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "debug/EdgeSet.h"
#include "graph/DynamicBipartiteGraph.hpp"
// #include "graph/DynamicGraph.h"
#include "BK/BronKerboschRmRClique.hpp"
#include "dataStruct/coreDisJoin.hpp"
#include "dataStruct/disJoinSet.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace CDSetRSH {
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
        const double *RCliqueCounting;
        explicit CompareRClique(const double *coreLeaf) : RCliqueCounting(coreLeaf) {
        }

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
        std::vector<double> rCliqueSCounting(cliqueHashmap.size(), 0.0);
        daf::Size count = 0;
        for (const auto &leaf: treeGraph.adj_list) {
            if (leaf.size() < r) {
                continue;
            }
            daf::CliqueSize pivotC = 0, keepC = 0;
            for (const auto &i: leaf) {
                if (i.isPivot) pivotC++;
                else keepC++;
            }
            int needPivot = s - static_cast<int>(keepC);
            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                daf::CliqueSize subNumKeepC = 0;
                daf::CliqueSize subNumPovit = 0;
                for (const auto &node: rClique) {
                    if (node.isPivot) {
                        subNumPovit++;
                    } else {
                        subNumKeepC++;
                    }
                }
                if (subNumPovit > needPivot) {
                    return true;
                }
                auto ncrValue = nCr[pivotC - subNumPovit][needPivot - subNumPovit];
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
    void printEdgeCore(const Graph &edgeGraph, const std::vector<T> &coreE) {
        printEdgeCore(edgeGraph, coreE.data());
    }
}


std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionHierarchy(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {
    auto time_start = std::chrono::high_resolution_clock::now();
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build",
                   [&]() {
                       cliqueIndex.build(tree, edgeGraph.adj_list.size());
                   });

    auto countingRClique = daf::timeCount("countingPerEdgeAndRClique",
                                          [&]() {
                                              return CDSetRSH::countingPerRClique(
                                                  tree, cliqueIndex, r, s);
                                          });

    std::vector<double> coreRClique(countingRClique.size());
#ifndef NDEBUG
    tree.printGraphPerV();
    for (daf::Size i = 0; i < countingRClique.size(); ++i) {
        std::cout << "Clique: " << cliqueIndex.byId(i) << " id: " << i
                << " count: " << countingRClique[i] << std::endl;
    }
    std::cout << "cliqueIndex Size : " << cliqueIndex.size() << std::endl;
#endif

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

    CDSetRSH::DHeap heap{CDSetRSH::CompareRClique(countingRClique.data())};
    heap.reserve(cliqueIndex.size());

    std::vector<CDSetRSH::DHeap::handle_type> heapHandles(cliqueIndex.size());

    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
        heapHandles[i] = heap.push(i);
    }

    CoreDisJoin hierarchyBuilder(cliqueIndex.size(), 10);
#ifndef NDEBUG
    std::cout << "countingKE: ";
    std::cout << "countingRClique: " << countingRClique << std::endl;
    std::cout << "tree: ";
    tree.printGraphPerV();
    std::cout << "treeGraphV: ";
    treeGraphV.printGraphPerV();
#endif

    std::cout << "=========================begin=========================" << std::endl;
    double minCore = 0;
    std::vector<std::pair<double, daf::Size>> core_stats_map;
    std::map<double, size_t> node_stats_map;
    // edgeGraph.getNbrCount(154);
    // std::cout << "154 nbr count: " << edgeGraph.getNbrCount(154) << std::endl;
    while (!heap.empty()) {
        for (auto &leafId: changedLeaf) {
            changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
        }
        changedLeaf.clear();
        removedRCliqueIdForLeaf.clear();
        currentRemoveRcliqueIds.clear();

        auto currMinCore = countingRClique[heap.top()];
        if (currMinCore > minCore) {
                std::set<daf::Size> unique_nodes;
                for (const auto& adj : tree.adj_list) {
                    for (const auto& node : adj) {
                        unique_nodes.insert(node.v);
                    }
                }
                size_t current_living_nodes = unique_nodes.size();
                node_stats_map[currMinCore] = current_living_nodes;
                if (unique_nodes.contains(154)) {
                    std::cout << "154 is alive at core " << currMinCore << std::endl;
                } else {
                    std::cout << "154 is dead at core " << currMinCore << std::endl;
                }

                core_stats_map.push_back({currMinCore, hierarchyBuilder.currKIndex});
                hierarchyBuilder.addK(currMinCore);
            minCore = currMinCore;
        }

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
        for (auto leafId: changedLeaf) {
            const auto &leaf = tree.adj_list[leafId];
            auto leafIndex = changedLeafIndex[leafId];

            auto initCore = [&](const std::vector<TreeGraphNode> &newLeaf, const daf::Size &newLeafId) {
                daf::CliqueSize newPivotC = 0, newKeepC = 0;
                for (const auto &i: newLeaf) {
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

                            std::abort();
                        }
                        auto ncrValue = nCr[newPivotC - subNumPovit][needPivot - subNumPovit];
                        auto cliqueIndexId = cliqueIndex.byClique(rclique);
                        countingRClique[cliqueIndexId] += ncrValue;
                    }

                    return true;
                });
            };

            for (const auto &leafV: leaf) {
                if (leafV.isPivot) {
                    treeGraphV.removeNbr(leafV.v, {leafId, true});
                } else {
                    treeGraphV.removeNbr(leafV.v, {leafId, false});
                }
            }

            const auto &removedR = removedRCliqueIdForLeaf[leafIndex];
            auto mapped = removedR | std::views::transform(
                              [&](const daf::Size id) {
                                  return cliqueIndex.byId(id);
                              }
                          );

            bkRmClique::removeRClique(leaf, mapped, r, s,
                                      [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                                          auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                                          DEBUG_BREAK_IF(newLeaf.size() < s);
                                          auto newId = tree.addNode(newLeaf);
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

            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &clique) {
                auto cliqueIndexId = cliqueIndex.byClique(clique);
                hierarchyBuilder.unite(cliqueIndexId, removedR[0]);

                if (!rCliqueInHeap[cliqueIndexId]) return true;
                daf::CliqueSize subNumKeepC = 0;
                daf::CliqueSize subNumPovit = 0;
                for (const auto &node: clique) {
                    if (node.isPivot) subNumPovit++;
                    else subNumKeepC++;
                }

                if (subNumPovit > s - keepC) {
                    return true;
                }

                auto ncrValue = nCr[pivotC - subNumPovit][s - keepC - subNumPovit];
                countingRClique[cliqueIndexId] -= ncrValue;
                heap.update(heapHandles[cliqueIndexId]);
                return true;
            });

            tree.removeNode(leafId);
        }
        std::cout << "Current Leaves: " << tree.size() << std::endl;

#ifndef NDEBUG
        std::cout << "tree: ";
        tree.printGraphPerV();
        std::cout << "treeGraphV: ";
        treeGraphV.printGraphPerV();
#endif
    }
#ifndef NDEBUG
    std::cout << "tree: ";
    tree.printGraphPerV();

    std::cout << "treeGraphV: ";
    treeGraphV.printGraphPerV();
#endif


    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;


    std::cout << "=========================Core Component Stats=========================" << std::endl;
    for (const auto& [core, level] : core_stats_map) {
        auto &ds = hierarchyBuilder.codeDisjointSets[level];
        size_t l_nodes = node_stats_map[core];
        std::cout << "Core: " << core << " Components: " << ds.live_count << " LivingNodes: " << l_nodes << std::endl;
    }
    std::cout << "======================================================================" << std::endl;

    std::cout << "Largest Core Value: " << coreRClique[std::size(coreRClique) - 1] << std::endl;

    // Parseable per-r-clique core distribution for correctness verification:
    // must be identical to the default CND engine (NucleusCoreDecompositionRClique).
    {
        std::map<long long, long long> cndh_dist;
        long long cndh_max = 0;
        for (double cv : coreRClique) {
            long long k = (long long)(cv + 0.5);
            cndh_dist[k]++;
            if (k > cndh_max) cndh_max = k;
        }
        for (const auto &kv : cndh_dist)
            std::cout << "core=" << kv.first << " count=" << kv.second << std::endl;
        std::cout << "[cnd-hier] Max core: " << cndh_max << std::endl;
    }


    std::vector<std::pair<std::vector<daf::Size>, double> > sortedK;

    // ================= START OF LOGGING BLOCK (Data Export for Python Vis) =================
    // 目标：导出 Nodes.csv (Id, MaxCore, ClusterId) 和 Edges.csv (Source, Target)
    // Gated: the in-memory hierarchy is fully built above; this block only
    // serialises it to CSV. Skipped by default so benchmarks measure the build
    // (not disk I/O) and never blow up disk on large no-timeout cells.
    if (std::getenv("PIVOTER_CND_HIER_DUMP")) {

    std::string outDir = "case_study_output";
    if (!std::filesystem::exists(outDir)) {
        std::filesystem::create_directory(outDir);
    }

    std::string fileSuffix = "_r" + std::to_string(r) + "_s" + std::to_string(s);
    std::cout << "[Logging] Processing and Exporting data for Visualization..." << std::endl;

    // --- 1. Map Core Values to Hierarchy Levels ---
    // CoreDisJoin uses levels to store disjoint sets. We need to map core values to these levels.
    std::map<double, size_t> coreToLevel;
    for(const auto& p : core_stats_map) {
        coreToLevel[p.first] = p.second;
    }

    // --- 2. Compute Node Properties (Max Core & Cluster ID) ---
    // 每个节点可能属于多个 r-clique (nucleus)。我们需要找到它所属的 Max Core。
    // 并在这个 Max Core 级别上，找到它所属的 Cluster ID (Connected Component Root)。

    struct NodeInfo {
        double maxCore = 0.0;
        daf::Size bestCliqueId = std::numeric_limits<daf::Size>::max(); // Representative clique
    };

    // 使用 vector 存储节点信息，index 即为 node_id (0 to N)
    const daf::Size numNodes = edgeGraph.adj_list_offsets.size() - 1;
    std::vector<NodeInfo> nodeStats(numNodes);

    // 遍历所有 r-clique，更新包含的节点的 Max Core
    for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
        double core = coreRClique[i];
        if (core == 0) continue; // 忽略 core 为 0 的

        auto nodes = cliqueIndex.byId(i);
        for (const auto& u : nodes) {
            if (u < numNodes) { // 安全检查
                if (core >= nodeStats[u].maxCore) {
                    nodeStats[u].maxCore = core;
                    nodeStats[u].bestCliqueId = i;
                    // print info for node 154
                }
            }
        }
    }

    // --- 3. Export Nodes CSV ---
    // Format: node_id, core, cluster_id
    std::string nodePath = outDir + "/nodes" + fileSuffix + ".csv";
    std::ofstream nodeOfs(nodePath);
    nodeOfs << "node_id,core,cluster_id\n";

    for (daf::Size nodeId = 0; nodeId < numNodes; ++nodeId) {
        double core = nodeStats[nodeId].maxCore;
        long long clusterId = -1;

        // 只有当 core > 0 时才有 cluster
        if (core > 0 && coreToLevel.count(core)) {
            size_t level = coreToLevel[core];
            daf::Size cliqueId = nodeStats[nodeId].bestCliqueId;

            // Access DSU at this level
            // codeDisjointSets 是 public 的，可以直接访问
            auto& dsu = hierarchyBuilder.codeDisjointSets[level];

            // 使用 find() 获取 Cluster ID (Component Root)
            // 假设 disJoinSet 实现了 standard find()。
            // 如果 find 是 private, 可能需要手动 traverse parent_。
            // 根据之前代码的注释 `ds.parent_` 是可访问的。
            // 这里我们尝试使用 dsu.find(cliqueId)。
            clusterId = dsu.find(cliqueId);
        }

        // 仅导出 core > 0 的节点 (或者全部导出让 python 过滤)
        // 为了作图完整性，全部导出，Python 侧过滤噪音。
        nodeOfs << nodeId << "," << core << "," << clusterId << "\n";
    }
    nodeOfs.close();
    std::cout << "[Logging] Node attributes saved: " << nodePath << std::endl;

    // --- 4. Export Edges CSV ---
    // Format: source, target
    std::string edgePath = outDir + "/edges" + fileSuffix + ".csv";
    std::ofstream edgeOfs(edgePath);
    edgeOfs << "source,target\n";

    for (daf::Size u = 0; u < numNodes; ++u) {
        const daf::Size start = edgeGraph.adj_list_offsets[u];
        const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
        for (daf::Size idx = start; idx < end; ++idx) {
            daf::Size v = edgeGraph.adj_list[idx];
            if (u < v) { // Undirected edges, write once
                edgeOfs << u << "," << v << "\n";
            }
        }
    }
    edgeOfs.close();
    std::cout << "[Logging] Graph edges saved: " << edgePath << std::endl;
    } // end if PIVOTER_CND_HIER_DUMP

    // --- 5. Export Hierarchy Metadata (Optional but useful for stats) ---
    // 仍然保留 Metadata json 以备不时之需，或用于 debugging
    // std::string metaPath = outDir + "/metadata" + fileSuffix + ".json";
    // std::ofstream metaOfs(metaPath);
    // metaOfs << "[\n";
    // bool first = true;
    // for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
    //     if (coreRClique[i] == 0) continue; // 只导出有意义的
    //
    //     if (!first) metaOfs << ",\n";
    //     first = false;
    //
    //     auto nodes = cliqueIndex.byId(i);
    //     metaOfs << "  {\"id\": " << i
    //             << ", \"core\": " << coreRClique[i]
    //             << ", \"nodes\": [";
    //     for (size_t j = 0; j < nodes.size(); ++j) {
    //         metaOfs << nodes[j] << (j < nodes.size() - 1 ? "," : "");
    //     }
    //     metaOfs << "]}";
    // }
    // metaOfs << "\n]\n";
    // metaOfs.close();

    // ================= END OF LOGGING BLOCK =================

    return sortedK;
}


namespace CDSetRSH {
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
}