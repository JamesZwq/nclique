// //
// // Created by _ on 25-5-26.
// //
//
// #ifndef BRONKERBOSCH_HPP
// #define BRONKERBOSCH_HPP
//
// #include <algorithm>
// #include <algorithm>
// #include <bitset>
// #include <vector>
// #include <ranges>
// #include "Global/Global.h"
// #include "graph/DynamicGraph.h"
// #include <boost/dynamic_bitset.hpp>
//
//
//
// namespace BronKerbosch {
//     using Bitset = boost::dynamic_bitset<>;
//
//     /**
//      * ， n 
//      */
//     template<class F>
//     void for_each_bit(const Bitset &bs, int n, F &&callback) {
//         //  1
//         for (size_t v = bs.find_first(); v != Bitset::npos && (int)v < n; v = bs.find_next(v)) {
//             // bs.test(v)  true，
//             if (!callback((int)v)) break;
//         }
//     }
//
//
//     [[nodiscard]] inline std::vector<TreeGraphNode>
//     coverToVertex(const Bitset &cover,
//                   const Bitset &pivots,
//                   const std::vector<TreeGraphNode> &vList) {
//         // std::cout << "coverToVertex: " << cover << std::endl;
//         // std::cout << "pivots: " << pivots << std::endl;
//         // std::cout << "vList: " << vList << std::endl;
//         std::vector<TreeGraphNode> result;
//         result.reserve(cover.count());
//
//         // cover  1 
//         auto i = cover.find_first();
//         // pivots  1 
//         auto pj = pivots.find_first();
//
//         //  cover  1 ，
//         while (i != Bitset::npos && i < vList.size()) {
//             //  pj  >= i
//             while (pj != Bitset::npos && pj < i) {
//                 pj = pivots.find_next(pj);
//             }
//             //  pj == i， pivot
//             bool isP = (pj == i);
//
//             // 
//             result.emplace_back(vList[i].v, isP);
//
//             //  cover  1 
//             i = cover.find_next(i);
//         }
//
//         return result;
//     }
//     /**
//      *  class  run() ，，
//      *  adj, n, minK ，
//      *  report 
//      */
//     template<class ReportFn>
//     void bk_run(const std::vector<Bitset> &adj,
//                 int n,
//                 int minK,
//                 Bitset R,
//                 Bitset P,
//                 Bitset pivots,
//                 ReportFn &&report) {
//         // 1)  P,X ， R
//         if (P.none()) {
//             if ((int) R.count() >= minK) {
//                 report(R, pivots);
//             }
//             return;
//         }
//
//         // 2)  pivot u ∈ P∪X， |P ∧ nbr(u)| 
//         int bestU = -1, bestCnt = -1;
//         Bitset PX = P;
//         for_each_bit(PX, n, [&](int u) {
//             Bitset nbr = adj[u] & P;
//             int cnt = (int) nbr.count();
//             if (cnt > bestCnt) {
//                 bestCnt = cnt;
//                 bestU = u;
//             }
//             return true;
//         });
//         // std::cout << "bestU: " << bestU << std::endl;
//         // 3)  P \ nbr(bestU)
//         Bitset candidates = P & ~adj[bestU];
//         for_each_bit(candidates, n, [&](int v) {
//             Bitset R2 = R;    R2.set(v);
//             Bitset P2 = P & adj[v];
//             Bitset X2 = adj[v];
//
//             // 1)  pivots
//             Bitset piv2 = pivots;
//             // 2)  v  pivot， piv2 
//             if (v == bestU) piv2.set(v);
//
//             // 3) 
//             bk_run(adj, n, minK, R2, P2, piv2, report);
//
//             P.reset(v);
//             return true;
//         });
//     }
//
//     /**
//      *  constructor ：
//      *   -  vList
//      *   -  vListMap
//      *   - ， removeEdgeList 
//      *  adj，  outN/outMinK  n/minK
//      */
//     inline std::vector<Bitset>
//         build_adj(std::vector<daf::Size> &vList,
//                   DynamicGraph<daf::Size> &edgeGraph,
//                   Bitset &staticVertex,  // ： removeEdgeList 
//                   Bitset &pivot,  // ： removeEdgeList 
//                   int &outN) {
//         std::ranges::sort(vList);
//         int n = static_cast<int>(vList.size());
//         outN = n;
//
//         // 
//         staticVertex.resize(n);
//         staticVertex.set();
//
//         pivot.resize(n);
//         pivot.reset();
//
//         // 
//         std::vector<Bitset> adj(n, Bitset(n));
//         for (int i = 0; i < n; ++i) {
//             daf::vListMap[vList[i]] = i;
//             adj[i].set();
//             adj[i].reset(i);
//             if (vList[i].isPivot) {
//                 pivot.set(i);
//             }
//         }
//         //  removeEdgeList ， staticVertex 
//         for (auto [u0, v0] : removeEdgeList) {
//             if (daf::vListMap[u0] == std::numeric_limits<daf::Size>::max() ||
//                     daf::vListMap[v0] == std::numeric_limits<daf::Size>::max()) {
//                 continue;
//             }
//             auto u = daf::vListMap[u0];
//             auto v = daf::vListMap[v0];
//             adj[u].reset(v);
//             adj[v].reset(u);
//             staticVertex.reset(u);
//             staticVertex.reset(v);
//             pivot.reset(u);
//             pivot.reset(v);
//         }
//
//         for (int i = 0; i < n; ++i) {
//             daf::vListMap[vList[i]] = std::numeric_limits<daf::Size>::max();
//         }
//         return adj;
//     }
//
//
//     /**
//      * ， main  new  BronKerbosch(...)  call run：
//      *
//      *   bronKerbosch(vList, removeEdgeList, minK, report);
//      *
//      *  report(Bitset clique) 。
//      */
//     template<class ReportFn>
//     void bronKerbosch(std::vector<TreeGraphNode> &vList,
//                       DynamicGraph<daf::Size> &edgeGraph,
//                       int minK,
//                       ReportFn &&report) {
//         int n;
//         Bitset R;
//         Bitset povit;
//         auto adj = build_adj(vList, edgeGraph, R, povit, n);
//
//         Bitset mask(n); mask.set();
//         Bitset P = (~R) & mask;
//         // X 
//         Bitset X(n);
//         // ， R
//         bk_run(adj, n, minK, R, P, povit, std::forward<ReportFn>(report));
//     }
//
// };
//
//
//
// #endif //BRONKERBOSCH_HPP
