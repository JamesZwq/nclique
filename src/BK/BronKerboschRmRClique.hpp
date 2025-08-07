// BronKerboschFns.hpp
#ifndef BKClique_HPP
#define BKClique_HPP

#include <algorithm>
#include <algorithm>
#include <bitset>
#include <vector>
#include <ranges>
#include "Global/Global.h"
#include "graph/DynamicGraph.h"
#include <boost/dynamic_bitset.hpp>

extern double nCr[1001][401];

namespace bkRmClique {
    // static constexpr int MAXN = 400;
    using Bitset = boost::dynamic_bitset<>;

    /**
     * 和原来成员版本一模一样，只是把 n 也作为参数传进来
     */
    template<class F>
    void for_each_bit(const Bitset &bs, int n, F &&callback) {
        // 首先找第一个 1
        for (size_t v = bs.find_first(); v != Bitset::npos && (int) v < n; v = bs.find_next(v)) {
            // bs.test(v) 肯定为 true，不用再测
            if (!callback((int) v)) break;
        }
    }

    inline void printBitset(const Bitset &bs, std::string name = "") {
        // BronKerboschFns.hpp
        if (!name.empty()) {
            std::cout << name << ": ";
        }
        for_each_bit(bs, (int) bs.size(), [&](int v) {
            std::cout << v << ' ';
            return true;
        });
        std::cout << std::endl;
    }

    template<typename ConflictSetsRange, typename ConflictSetsReverseRange, typename ReportFn>
    void pathSplit(int n, int r, int minK,
                   Bitset &P,
                   Bitset &pivots,
                   daf::StaticVector<daf::Size> &conflictCount,
                   const daf::StaticVector<daf::Size> &conflictMaxSize,
                   ConflictSetsRange &conflictSets,
                   ConflictSetsReverseRange &conflictSetsReverse,
                   daf::Size nextCid, ReportFn &&report) {
        // 直接调用 bk_run
        // std::cout << "===================================================================" << std::endl;
        // printBitset(P, "P");
        // printBitset(pivots, "pivots");
        // std::cout << "conflictCount: " << conflictCount << std::endl;
        // std::cout << "conflictMaxSize: " << conflictMaxSize << std::endl;
        // std::cout << "conflictSets: " ;
        // for (const auto &set : conflictSets) {
        //     std::cout << "{ ";
        //     for (auto v : set) {
        //         std::cout << v << ' ';
        //     }
        //     std::cout << " } ";
        // }
        // std::cout << std::endl;
        // std::cout << "conflictSetsReverse: " << conflictSetsReverse << std::endl;

        auto currN = P.count();
        if (currN < minK || currN - pivots.count() > minK) {
            return;
        }

        if (currN - pivots.count() == minK) {
            Bitset reportPivots(n);
            reportPivots.reset(); // 初始化为全 false
            report(P & (~pivots), reportPivots); // 直接报告当前的 P 和 pivots
            return;
        }

        // for (int i = 0; i < n; ++i) {
        //     daf::Size maxRClique = nCr[n-1][r-1];
        //     if (conflictSetsReverse[i].size() >= maxRClique) {
        //         P.reset(i); // 从 P 中删除这个点
        //         if (pivots.test(i)) {
        //             pivots.reset(i); // 从 pivots 中删除这个点
        //         }
        //         else {
        //             return;
        //         }
        //         for (auto cid: conflictSetsReverse[i]) {
        //             conflictCount[cid]--;
        //         }
        //     }
        // }
        // bkRmClique::for_each_bit(P, n, [&](int v) {
        //     auto numConflict = 0;
        //     for (auto cid: conflictSetsReverse[v]) {
        //         if (conflictCount[cid] >= conflictMaxSize[cid]) {
        //             numConflict++;
        //         }
        //     }
        //     if (numConflict >= nCr[currN - 1][r - 1]) {
        //         P.reset(v); // 从 P 中删除这个点
        //         if (pivots.test(v)) {
        //             pivots.reset(v); // 从 pivots 中删除这个点
        //             for (auto cid: conflictSetsReverse[v]) {
        //                 conflictCount[cid]--;
        //             }
        //         } else {
        //             return false; // 直接返回，后续不需要处理
        //         }
        //     }
        //     return true; // 继续处理下一个点
        // });

        bool noConflict = true;
        for (daf::Size cid = nextCid; cid < conflictCount.size(); ++cid) {
            if (conflictCount[cid] >= conflictMaxSize[cid]) {
                noConflict = false; // 有冲突
                std::vector<daf::Size> conflictInPivot;
                conflictInPivot.reserve(conflictMaxSize[cid]);
                for (auto v: conflictSets[cid]) {
                    auto idx = daf::vListMap[v];
                    if (pivots.test(idx)) {
                        conflictInPivot.push_back(idx);
                    }
                }
                for (int i = 0; i < conflictInPivot.size(); ++i) {
                    P.reset(conflictInPivot[i]);
                    pivots.reset(conflictInPivot[i]);

                    for (auto g: conflictSetsReverse[conflictInPivot[i]]) {
                        conflictCount[g]--;
                    }
                    if (i != 0) {
                        pivots.reset(conflictInPivot[i - 1]);
                        P.set(conflictInPivot[i - 1]);
                    }
                    pathSplit(n, r, minK, P, pivots,
                              conflictCount, conflictMaxSize, conflictSets, conflictSetsReverse,
                              cid + 1, std::forward<decltype(report)>(report));

                    for (auto g: conflictSetsReverse[conflictInPivot[i]]) {
                        conflictCount[g]++;
                    }

                    P.set(conflictInPivot[i]); // 恢复 P
                    pivots.set(conflictInPivot[i]); // 恢复 pivots
                }
                // reset pivots
                for (int i = 0; i + 1 < conflictInPivot.size(); ++i) {
                    pivots.set(conflictInPivot[i]);
                }
                return;
            }
        }
        if (noConflict) {
            // 没有冲突，直接运行原来的 bk_run
            if (currN >= minK) {
                report(P, pivots); // 直接报告当前的 P 和 pivots
            }
        }
    }


    /**
     * 对外主入口，等同于原来你在 main 里 new 一个 BronKerbosch(...) 然后 call run：
     *
     *   bronKerbosch(vList, removeEdgeList, minK, report);
     *
     * 其中 report(Bitset clique) 会在每次找到一个极大团时被调用。
     */
    template<
        std::ranges::input_range VListRange,
        std::ranges::input_range ConflictSetsRange,
        typename ReportFn
    >
        requires
        std::same_as<std::ranges::range_value_t<VListRange>, TreeGraphNode> &&
        std::ranges::input_range<std::ranges::range_value_t<ConflictSetsRange> > &&
        std::same_as<
            std::ranges::range_value_t<
                std::ranges::range_value_t<ConflictSetsRange>
            >,
            daf::Size
        >
    void removeRClique(VListRange &vList,
                       ConflictSetsRange &conflictSets,
                       int r,
                       int minK,
                       ReportFn &&report) {
        const int n = static_cast<int>(vList.size());
        // Bitset staticV;        // 真正的 staticVertex
        Bitset pivots(n);
        pivots.reset(); // 初始化为全 false
        std::vector<std::vector<daf::Size> > conflictSetsReverse;
        // daf::Size maxN = 0;
        conflictSetsReverse.resize(vList.size());

        for (int i = 0; i < n; ++i) {
            daf::vListMap[vList[i]] = i; // 建反向索引
            if (vList[i].isPivot) pivots.set(i);
        }


        Bitset P(n);
        P.set(); // 起始 P 为全部顶点（或 ~staticV 按需剪枝）
        // daf::StaticVector<std::pair<daf::Size, daf::Size> > removeEdgeListPairs;
        daf::StaticVector<daf::Size> conflictCount, conflictMaxSize;
        conflictCount.resize(conflictSets.size());
        conflictMaxSize.resize(conflictSets.size());
        for (int i = 0; i < conflictSets.size(); ++i) {
            for (auto j = 0; j < conflictSets[i].size(); ++j) {
                auto idx = conflictSets[i][j];
                if (conflictSetsReverse[daf::vListMap[idx]].empty()) {
                    conflictSetsReverse[daf::vListMap[idx]].reserve(conflictSets[i].size());
                }
                conflictSetsReverse[daf::vListMap[idx]].push_back(i);
            }
            conflictCount[i] = conflictSets[i].size();
            conflictMaxSize[i] = conflictSets[i].size();
        }
        // 如果一个点出现在所有 conflictSets 中, 那么则先把它从 conflictSets 中删除
        for (int i = 0; i < n; ++i) {
            daf::Size maxRClique = nCr[n-1][r-1];
            if (conflictSetsReverse[i].size() >= maxRClique) {
                P.reset(i); // 从 P 中删除这个点
                if (P.count() < minK) {
                    return; // 如果 P 的大小小于 minK，直接返回
                }
                if (pivots.test(i)) {
                    pivots.reset(i); // 从 pivots 中删除这个点
                }
                else {
                    return;
                }
                for (auto cid: conflictSetsReverse[i]) {
                    conflictCount[cid]--;
                }
            }
        }

        pathSplit(n, r, minK, P, pivots,
                  conflictCount, conflictMaxSize, conflictSets, conflictSetsReverse,
                  0, std::forward<ReportFn>(report));
    }

    inline void testBronKerbosch() {
        // 3 个点完全图：0–1, 0–2, 1–2
        // 假设 TreeGraphNode 的构造函数为 TreeGraphNode(size_t id, bool flag)
        std::vector<TreeGraphNode> vList = {
            {1, false}, // 原来是 3
            {2, true}, // 原来是 4
            {3, true}, // 原来是 6
            {4, true}, // 原来是 7
            {5, true}, // 原来是 8
            {6, true}, // 原来是 9
        };

        std::vector<std::vector<daf::Size> > conflictSets;
        conflictSets.emplace_back(std::vector<daf::Size>{1, 2, 6}); // {3,4,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{1, 3, 6});  // {3,6,9}
        conflictSets.emplace_back(std::vector<daf::Size>{4, 5, 6}); // {7,8,9}
        conflictSets.emplace_back(std::vector<daf::Size>{3, 5, 6}); // {6,8,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{1, 4, 6});  // {3,7,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{1, 5, 6});  // {3,8,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{2, 3, 6});  // {4,6,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{2, 5, 6});  // {4,8,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{2, 4, 6});  // {4,7,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{3, 4, 6});  // {6,7,9}
        int minK = 4; // 最小 clique 大小
        std::vector<double> cliqueCounts(vList.size(), 0);
        removeRClique(vList, conflictSets, 3, minK,
                      [&](const Bitset &clique, const Bitset &pivots) {
                          std::vector<int> C, P, H;
                          for (size_t i = clique.find_first(); i != Bitset::npos; i = clique.find_next(i)) {
                              C.push_back(vList[i].v);
                              (pivots.test(i) ? P : H).push_back(vList[i].v);
                          }
                          std::cout << "!!!!!!!!!!! Clique: ";
                          for (int x: C) std::cout << x << ' ';
                          std::cout << "| Pivots: ";
                          for (int x: P) std::cout << x << ' ';
                          std::cout << "| Holds: ";
                          for (int x: H) std::cout << x << ' ';
                          std::cout << '\n';

                          // use ncr do clique counting
                          for (size_t i = H.size(); i <= C.size(); ++i) {
                              auto need = C.size() - i;
                              cliqueCounts[i] += nCr[P.size()][need];
                          }
                      });

        std::cout << "Clique counts: " << cliqueCounts << std::endl;
    }
} // namespace bk


#endif // BKClique_HPP
