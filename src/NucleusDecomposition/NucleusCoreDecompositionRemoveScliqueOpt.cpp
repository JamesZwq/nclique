// //
// // Created by 张文谦 on 25-3-4.
// //
//
// #include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
// #include <boost/heap/d_ary_heap.hpp>
// #include <boost/heap/pairing_heap.hpp>
// #include <boost/heap/fibonacci_heap.hpp>
// #include <set>
//
// #include "../BK/BronKerboschRmEdge.hpp"
// #include "dataStruct/CliqueHashMap.h"
// #include "debug/EdgeSet.h"
// #include "graph/DynamicBipartiteGraph.hpp"
// // #include "graph/DynamicGraph.h"
// #include "BK/BronKerboschRmRClique.hpp"
// #include "dataStruct/disJoinSet.hpp"
// #include "graph/DynamicGraphSet.h"
//
// extern double nCr[1001][401];
// // 放在你的函数外（比如文件顶部），保证编译时可见并内联
// // #ifndef NDEBUG
// // set NOEBUG as trus
//
//
//
//
// namespace CDSetRSO {
//     struct LeafRmInfo {
//         bool removedKeepC;
//         daf::StaticVector<daf::Size> removedPivots{0};
//         daf::StaticVector<std::pair<daf::Size, daf::Size>> removedEdges{0};
//         daf::StaticVector<daf::Size> removedRClique{0};
//         daf::StaticVector<daf::Size> vertexRCliqueRemining{0};
//         //存储edge的地方, 对于
//         robin_hood::unordered_flat_map<uint64_t, daf::Size> adj; // 邻接表，存储每个点的邻居
//         LeafRmInfo() : removedKeepC(false) {
//         }
//
//         bool empty() const {
//             return !removedKeepC && removedPivots.empty() && removedEdges.empty();
//         }
//
//         void init(const daf::CliqueSize r, const std::vector<daf::Size> &leaf) {
//             removedKeepC = false;
//             auto leafSize = leaf.size();
//             removedPivots.reserve(leafSize);
//             vertexRCliqueRemining.resize(leafSize);
//             for (daf::Size i = 0; i < leafSize; ++i) {
//                 vertexRCliqueRemining[i] = nCr[leafSize-1][r-1]; // 每个点在 r 大小的 clique 中的数量
//             }
//             removedEdges.reserve(leafSize * (leafSize - 1) / 2); // k choose 2
//             removedRClique.reserve(leafSize * leafSize); // k choose 2
//
//             adj.reserve(leafSize * (leafSize - 1) / 2);
//         }
//
//         void reduceEdgeCount(daf::Size u, daf::Size v, daf::Size leafSize, daf::Size r,  daf::Size reduceCount = 1) {
//             if (u > v) std::swap(u, v);
//             auto key = (static_cast<uint64_t>(u) << 32) | v; // 使用 u 和 v 的组合作为键
//             auto it = adj.try_emplace(key, nCr[leafSize-2][r-2]); // 尝试插入，如果不存在则初始化为 nCr[leafSize-2][r-2]
//
//             if (reduceCount > it.first->second) {
//                 std::cerr << "Error: reduceCount exceeds current count for edge (" << u << ", " << v << ")." << std::endl;
//                 exit(1);
//             }
//             it.first->second -= reduceCount; // 如果已经存在，减少计数
//
//             if (it.first->second == 0) {
//                 removedEdges.emplace_back(u, v); // 如果计数为 0，添加到 removedEdges
//             }
//         }
//
//         void clear() {
//             removedKeepC = false;
//             removedPivots.clear();
//             removedEdges.clear();
//             vertexRCliqueRemining.clear();
//             removedRClique.clear();
//             adj.clear();
//         }
//
//         ~LeafRmInfo() {
//             removedPivots.free();
//             removedEdges.free();
//             vertexRCliqueRemining.free();
//             removedRClique.free();
//         }
//
//         friend std::ostream &operator<<(std::ostream &os, const LeafRmInfo &info) {
//             os << "removedKeepC: " << info.removedKeepC << "\n removedPivots: " << info.removedPivots
//                     << "\n removedEdges: " << info.removedEdges << "\n vertexRCliqueRemining: " << info.
//                     vertexRCliqueRemining
//                     << "\n removedRClique: " << info.removedRClique;
//             return os;
//         }
//     };
//
//     template<typename It1, typename It2, typename UpdateFunc>
//     inline void processEdgePairsImpl(It1 b1, It1 e1,
//                                      It2 b2, It2 e2,
//                                      double weight,
//                                      UpdateFunc &&upd) noexcept {
//         if (weight < 0.0) return;
//
//         // 判断两个区间迭代器是否相同
//         // if all same, do nothing
//         if (b1 == b2 && e1 == e2 && b1 == e1 && b2 == e2) {
//             return;
//         }
//         if (b1 == b2 && e1 == e2) {
//             // 同一范围：i < j
//             for (auto it = b1; it + 1 != e1; ++it) {
//                 auto u = *it;
//                 for (auto jt = it + 1; jt != e1; ++jt) {
//                     upd(u, *jt, weight);
//                 }
//             }
//         } else {
//             // 不同范围：笛卡尔积
//             for (auto it = b1; it != e1; ++it) {
//                 auto u = *it;
//                 for (auto jt = b2; jt != e2; ++jt) {
//                     upd(u, *jt, weight);
//                 }
//             }
//         }
//     }
//
//     template<
//         typename Range1, typename Range2,
//         typename UpdateFunc
//     >
//     inline void processEdgePairs(const Range1 &r1,
//                                  const Range2 &r2,
//                                  double weight,
//                                  UpdateFunc &&upd) noexcept {
//         processEdgePairsImpl(
//             std::begin(r1), std::end(r1),
//             std::begin(r2), std::end(r2),
//             weight,
//             std::forward<UpdateFunc>(upd)
//         );
//     }
//
//     template<
//         typename Range,
//         typename UpdateFunc
//     >
//     inline void processEdgePairs(const Range &r,
//                                  double weight,
//                                  UpdateFunc &&upd) noexcept {
//         processEdgePairsImpl(
//             std::begin(r), std::end(r),
//             std::begin(r), std::end(r),
//             weight,
//             std::forward<UpdateFunc>(upd)
//         );
//     }
//
//
//     struct CompareRClique {
//         const double *RCliqueCounting; // 指向外部数组
//         explicit CompareRClique(const double *coreLeaf) : RCliqueCounting(coreLeaf) {
//         }
//
//         // 注意：这里要返回 “a 排在前面” 的条件，为最小堆写成 coreLeaf[a] > coreLeaf[b]
//         bool operator()(daf::Size const &a, daf::Size const &b) const {
//             return RCliqueCounting[a] > RCliqueCounting[b];
//         }
//     };
//
//     using DHeap = boost::heap::d_ary_heap<
//         daf::Size,
//         boost::heap::arity<8>,
//         boost::heap::mutable_<true>,
//         boost::heap::compare<CompareRClique>
//     >;
//
//     std::vector<double> countingPerRClique(
//         const DynamicGraph<TreeGraphNode> &treeGraph,
//         StaticCliqueIndex &cliqueHashmap,
//         const daf::CliqueSize r,
//         const daf::CliqueSize s) {
//         // double *rCliqueSCounting = new double[cliqueHashmap.size()];
//         // memset(rCliqueSCounting, 0, cliqueHashmap.size() * sizeof(double));
//         std::vector<double> rCliqueSCounting(cliqueHashmap.size(), 0.0);
//         daf::Size count = 0;
//         for (const auto &leaf: treeGraph.adj_list) {
//             if (leaf.size() < r) {
//                 continue;
//             }
//             daf::CliqueSize pivotC = 0, keepC = 0;
//             for (auto &i: leaf) {
//                 if (i.isPivot) pivotC++;
//                 else keepC++;
//             }
//
//             // clique
//             int needPivot = s - static_cast<int>(keepC);
//             daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
//                 daf::CliqueSize subNumKeepC = 0;
//                 daf::CliqueSize subNumPovit = 0;
//                 for (const auto &node: rClique) {
//                     if (node.isPivot) {
//                         subNumPovit++;
//                     } else {
//                         subNumKeepC++;
//                     }
//                 }
//
//                 auto ncrValue = nCr[pivotC - subNumPovit][needPivot - subNumPovit];
//                 // rCliqueSCounting[cliqueHashmap.byNewClique(rClique)] += ncrValue;
//                 auto [id, isNew] = cliqueHashmap.byNewClique(rClique);
//                 if (isNew) {
//                     if (rCliqueSCounting.size() <= id) {
//                         rCliqueSCounting.resize(std::min(id + 1.1, id * 1.3), 0.0);
//                     }
//                 }
//                 rCliqueSCounting[id] += ncrValue;
//
//
//                 return true;
//             });
//         }
//         rCliqueSCounting.shrink_to_fit();
//         return rCliqueSCounting;
//     }
//
//
//     template<typename T>
//     void printEdgeCore(const Graph &edgeGraph, const T *coreE) {
//         const daf::Size m = edgeGraph.adj_list.size();
//         const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;
//         for (daf::Size u = 0; u < n; ++u) {
//             const daf::Size start = edgeGraph.adj_list_offsets[u];
//             const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
//             for (daf::Size idx = start; idx < end; ++idx) {
//                 std::cout << "[" << u << ", " << edgeGraph.adj_list[idx] << "] " << coreE[idx] << " ";
//             }
//             std::cout << std::endl;
//         }
//     }
//
//     template<typename T>
//     void printEdgeCore(const Graph &edgeGraph, const std::vector<T> coreE) {
//         printEdgeCore(edgeGraph, coreE.data());
//     }
// }
//
//
// std::vector<std::pair<std::vector<daf::Size>, int> > NucleusCoreDecompositionRCliqueOpt(
//     DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
//     DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
//     auto time_start = std::chrono::high_resolution_clock::now();
//     StaticCliqueIndex cliqueIndex(r);
//     daf::timeCount("clique Index build",
//                    [&]() {
//                        cliqueIndex.build(tree, edgeGraph.adj_list.size());
//                    });
//     // tree.printGraphPerV();
//     // cliqueIndex.print();
//     // cliqueIndex.verify();
//
//
//     auto countingRClique = daf::timeCount("countingPerEdgeAndRClique",
//                                           [&]() {
//                                               return CDSetRSO::countingPerRClique(
//                                                   tree, cliqueIndex, r, s);
//                                           });
//
//     std::vector<daf::Size> coreRClique(countingRClique.size());
// #ifndef NDEBUG
//     tree.printGraphPerV();
//     // daf::printArray(countingKE, edgeGraph.adj_list.size());
//     // CDSetRSO::printEdgeCore(edgeGraph, countingKE);
//     // std::cout << "countingRClique" << countingRClique << std::endl;
//     for (daf::Size i = 0; i < countingRClique.size(); ++i) {
//         std::cout << "Clique: " << cliqueIndex.byId(i) << " id: " << i
//                   << " count: " << countingRClique[i] << std::endl;
//     }
//     std::cout << "cliqueIndex Size : " << cliqueIndex.size() << std::endl;
//     // CDSetRSO::printEdgeCore(edgeGraph, degreeE);
// #endif
//
//     std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
//     std::vector<CDSetRSO::LeafRmInfo> leafRmInfos;
//     std::vector<daf::Size> changedLeaf;
//     std::vector<daf::Size> currentRemoveRcliqueIds;
//
//     leafRmInfos.reserve(tree.adj_list.size());
//     changedLeaf.reserve(tree.adj_list.size());
//     currentRemoveRcliqueIds.reserve(cliqueIndex.size());
//
//
//     daf::StaticVector<bool> rCliqueInHeap(cliqueIndex.size());
//     // rCliqueInHeap.fill(true);
//     rCliqueInHeap.resize(cliqueIndex.size());
//     memset(rCliqueInHeap.getData(), true, cliqueIndex.size() * sizeof(bool));
//
//     CDSetRSO::DHeap heap{CDSetRSO::CompareRClique(countingRClique.data())};
//     heap.reserve(cliqueIndex.size());
//
//     std::vector<CDSetRSO::DHeap::handle_type> heapHandles(cliqueIndex.size());
//
//     for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
//         heapHandles[i] = heap.push(i);
//     }
// #ifndef NDEBUG
//     std::cout << "countingKE: ";
//     // CDSetRSO::printEdgeCore(edgeGraph, countingKE);
//
//     std::cout << "countingRClique: " << countingRClique << std::endl;
//
//     std::cout << "tree: ";
//     tree.printGraphPerV();
//
//     std::cout << "treeGraphV: ";
//     treeGraphV.printGraphPerV();
// #endif
//     std::cout << "=========================begin=========================" << std::endl;
//     double minCore = 0;
//
//
//     while (!heap.empty()) {
//         for (auto &leafId: changedLeaf) {
//             changedLeafIndex[leafId] = std::numeric_limits<daf::Size>::max();
//         }
//         changedLeaf.clear();
//         leafRmInfos.clear();
//         currentRemoveRcliqueIds.clear();
//
//         minCore = std::max(countingRClique[heap.top()], minCore);
//         // 一次循环把所有 core==minCore 的 leaf 全部 pop 出来
//         std::cout << "minCore: " << minCore << std::endl;
//         while (!heap.empty() && countingRClique[heap.top()] <= minCore) {
//             auto id = heap.top();
//             rCliqueInHeap[id] = false;
//             heap.pop();
//             currentRemoveRcliqueIds.push_back(id);
//             coreRClique[id] = minCore;
// #ifndef NDEBUG
//             std::cout << "removed Clique: " << cliqueIndex.byId(id) << " id: " << id
//                       << " core: " << countingRClique[id] << std::endl;
// #endif
//         }
//
//         for (auto rmRCliqueId: currentRemoveRcliqueIds) {
//             auto rClique = cliqueIndex.byId(rmRCliqueId);
//             // std::cout << "rClique: " << rClique << std::endl;
//             daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
//                                             [&](const TreeGraphNode &uClique) {
//                                                 // std::cout << "uClique: " << uClique << std::endl;
//                                                 auto &leafId = changedLeafIndex[uClique.v];
//                                                 if (leafId == std::numeric_limits<daf::Size>::max()) {
//                                                     leafId = leafRmInfos.size();
//                                                     leafRmInfos.emplace_back();
//                                                     changedLeaf.push_back(uClique.v);
//                                                 }
//                                                 auto &currRmInfo = leafRmInfos[leafId];
//
//                                             });
//         }
//         // std::cout << "changedLeaf: " << changedLeaf << std::endl;
//         // std::cout << "changedLeafIndex: " << changedLeafIndex << std::endl;
//         for (auto leafId: changedLeaf) {
//             auto leaf = tree.adj_list[leafId];
//             auto leafIndex = changedLeafIndex[leafId];
//             // std::cout << "============================================================" << std::endl;
//             // std::cout << "changed leafId: " << leafId << " leaf index: " << leafIndex
//             //           << " leaf: " << leaf << std::endl;
//             // std::cout << "removedRCliqueIdForLeaf: ";
//             // for (const auto &id: removedRCliqueIdForLeaf[leafIndex]) {
//             //     std::cout << id << " (" << cliqueIndex.byId(id) << ") ";
//             // }
//
//             // std::cout << std::endl;
//
//             auto initCore = [&](const std::vector<TreeGraphNode> &leaf, const daf::Size &leafId) {
//                 daf::CliqueSize newPivotC = 0, newKeepC = 0;
//                 for (auto i: leaf) {
//                     if (i.isPivot) {
//                         treeGraphV.addNbr(i.v, {leafId, true});
//                         newPivotC++;
//                     } else {
//                         treeGraphV.addNbr(i.v, {leafId, false});
//                         newKeepC++;
//                     }
//                 }
//                 daf::Size needPivot = s - newKeepC;
//
//                 daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &newLeaf) {
//                     daf::CliqueSize subNumKeepC = 0;
//                     daf::CliqueSize subNumPovit = 0;
//                     for (const auto &node: newLeaf) {
//                         if (node.isPivot) {
//                             subNumPovit++;
//                         } else {
//                             subNumKeepC++;
//                         }
//                     }
//                     auto ncrValue = nCr[newPivotC - subNumPovit][needPivot - subNumPovit];
//                     auto cliqueIndexId = cliqueIndex.byClique(newLeaf);
//                     countingRClique[cliqueIndexId] += ncrValue;
//
//                     return true;
//                 });
//             };
//
//             for (auto leafV: leaf) {
//                 if (leafV.isPivot) {
//                     treeGraphV.removeNbr(leafV.v, {leafId, true});
//                 } else {
//                     treeGraphV.removeNbr(leafV.v, {leafId, false});
//                 }
//             }
//
//             auto removedR = removedRCliqueIdForLeaf[leafIndex];
//             auto mapped = removedR | std::views::transform(
//                               [&](const daf::Size id) {
//                                   return cliqueIndex.byId(id);
//                               }
//                           );
//
//             // DEBUG_BREAK_IF(leafId == 5);
//             bkRmClique::removeRClique(leaf, mapped, s, [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
//                 auto newLeaf = bkRmEdge::coverToVertex(c, pivots, leaf);
//                 // std::cout << " newLeaf: " << newLeaf;
//                 auto newId = tree.addNode(newLeaf);
//                 // std::cout << " newId: " << newId << std::endl;
//                 initCore(tree.adj_list[newId], newId);
//                 if (newId >= changedLeafIndex.size()) {
//                     changedLeafIndex.resize(newId * 1.5, std::numeric_limits<daf::Size>::max());
//                 }
//             });
//
//             daf::CliqueSize keepC = 0, pivotC = 0;
//             for (const auto &node: leaf) {
//                 if (node.isPivot) pivotC++;
//                 else keepC++;
//             }
//
//             daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &clique) {
//                 auto cliqueIndexId = cliqueIndex.byClique(clique);
//                 // std::cout << "cliqueIndexId: " << cliqueIndexId
//                 //           << " clique: " << clique << std::endl;
//                 // rCliqueInHeap.print("rCliqueInHeap: ");
//                 if (!rCliqueInHeap[cliqueIndexId]) return true;
//                 daf::CliqueSize subNumKeepC = 0;
//                 daf::CliqueSize subNumPovit = 0;
//                 for (const auto &node: clique) {
//                     if (node.isPivot) subNumPovit++;
//                     else subNumKeepC++;
//                 }
//                 auto ncrValue = nCr[pivotC - subNumPovit][s - keepC - subNumPovit];
//                 countingRClique[cliqueIndexId] -= ncrValue;
//                 // std::cout << "ncrValue: " << ncrValue
//                 //           << " pivotC: " << pivotC
//                 //           << " subNumPovit: " << subNumPovit
//                 //           << " keepC: " << keepC
//                 //           << " subNumKeepC: " << subNumKeepC
//                 //           << std::endl;
//                 heap.update(heapHandles[cliqueIndexId]);
//                 return true;
//             });
//
//             tree.removeNode(leafId);
//         }
//
//
// #ifndef NDEBUG
//             std::cout << "tree: ";
//             tree.printGraphPerV();
//
//             std::cout << "treeGraphV: ";
//             treeGraphV.printGraphPerV();
//
//
//             std::cout << "clique countingRClique: ";
//             for (daf::Size i = 0; i < countingRClique.size(); ++i) {
//                 std::cout << i << ": " << countingRClique[i] << " " << cliqueIndex.byId(i) << std::endl;
//             }
// #endif
//     }
//     // currentRemoveLeafIds.clear();
// #ifndef NDEBUG
//         std::cout << "tree: ";
//         tree.printGraphPerV();
//
//         std::cout << "treeGraphV: ";
//         treeGraphV.printGraphPerV();
// #endif
//
//
//     std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
//         std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
//
//     // coreE
//     // daf::printArray(coreE, edgeGraph.adj_list.size());
//
//
//     // /Users/zhangwenqian/UNSW/pivoter/a
//     // std::sort(coreE, coreE + edgeGraph.adj_list.size());
//     std::vector<std::pair<std::vector<daf::Size>, int> > sortedK;
//     sortedK.reserve(countingRClique.size());
//
//     for (daf::Size i = 0; i < cliqueIndex.size(); ++i) {
//         auto clique = cliqueIndex.byId(i);
//         std::vector<daf::Size> cliqueCopy(clique.begin(), clique.end());
//         sortedK.emplace_back(cliqueCopy, coreRClique[i]);
//     }
//
//     // if (std::accumulate(countingRClique.begin(), countingRClique.end(), 0.0) != 0) {
//     //     std::cerr << "Error: countingRClique != 0" << std::endl;
//     //     std::cerr << "countingRClique: " << countingRClique << std::endl;
//     //     std::exit(1);
//     // }
//     // std::cout << "sortedK size: " << sortedK << std::endl;
//
//     // /Users/zhangwenqian/UNSW/pivoter/a
//     std::sort(sortedK.begin(), sortedK.end(),
//               [](const auto &a, const auto &b) {
//                   return a.second < b.second; // 按照 core 值降序排序
//               });
//     auto file = fopen("/Users/zhangwenqian/UNSW/pivoter/a.out", "w");
//     for (auto i: sortedK) {
//         fprintf(file, "%d\n", (int) i.second);
//     }
//     fclose(file);
//     return sortedK;
// }
//
//
// template<class Bitset>
// void print_clique(const Bitset &bs) {
//     std::cout << '[';
//     bool first = true;
//     bkRmEdge::for_each_bit(bs, (int) bs.size(), [&](int v) {
//         if (!first) std::cout << ',';
//         first = false;
//         std::cout << v;
//         return true;
//     });
//     std::cout << "]\n";
// }
