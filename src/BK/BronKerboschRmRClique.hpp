// BronKerboschFns.hpp
#ifndef BKClique_HPP
#define BKClique_HPP

#include <algorithm>
#include <bitset>
#include <vector>
#include <ranges>
#include "Global/Global.h"
#include "graph/DynamicGraph.h"
#include <boost/dynamic_bitset.hpp>
#include <limits>
#include <cassert>
#include <iostream>

extern double nCr[1001][401];
// ---------------- fast dynamic bitset  ( ≤ 400 bits ) ----------------
#pragma once
#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <cassert>


namespace bkRmClique {
    // static constexpr int MAXN = 400;
    // using Bitset = DynBitset;
    // using Bitset = boost::dynamic_bitset<>;
    using Bitset = DynBitset;
    using VIdx = uint16_t;
    static_assert(std::numeric_limits<VIdx>::max() >= 400, "VIdx must hold n<=400");
    static DynBitset R, P, pivots, emptyPiv, tmp2;
    // Reusable static buffers (single-threaded use)
    static std::vector<uint32_t> s_csOff, s_rsOff, s_rsCol, s_deg, s_cur;
    static std::vector<VIdx> s_csCol;
    /**
     * 和原来成员版本一模一样，只是把 n 也作为参数传进来
     */
    template<class F>
    bool for_each_bit(const Bitset &bs, int n, F &&callback) {
        // 首先找第一个 1
        for (size_t v = bs.find_first(); v != Bitset::npos && (int) v < n; v = bs.find_next(v)) {
            // bs.test(v) 肯定为 true，不用再测
            if (!callback((int) v)) {
                return false; // 如果回调返回 false，提前结束遍历
            }
        }

        return true; // 返回 true 表示遍历完成
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

    [[nodiscard]] inline std::vector<TreeGraphNode>
    coverToVertex(const Bitset &cover,
                  const Bitset &pivots,
                  const std::vector<TreeGraphNode> &vList) {
        // std::cout << "coverToVertex: " << cover << std::endl;
        // std::cout << "pivots: " << pivots << std::endl;
        // std::cout << "vList: " << vList << std::endl;
        std::vector<TreeGraphNode> result;
        result.reserve(cover.count());

        // cover 上第一个 1 位
        auto i = cover.find_first();
        // pivots 上第一个 1 位
        auto pj = pivots.find_first();

        // 只要 cover 还有 1 位，就继续
        while (i != Bitset::npos && i < vList.size()) {
            // 把 pj 移到 >= i
            while (pj != Bitset::npos && pj < i) {
                pj = pivots.find_next(pj);
            }
            // 如果 pj == i，就说明这个位置是 pivot
            bool isP = (pj == i);

            // 把这个节点加入结果
            result.emplace_back(vList[i].v, isP);

            // 移动到 cover 的下一个 1 位
            i = cover.find_next(i);
        }
        return result;
    }


    template<class ReportFn>
    void edgeSplit(const std::vector<Bitset> &adj,
                   int n,
                   int minK,
                   Bitset R,
                   Bitset P,
                   Bitset pivots,
                   ReportFn &&report) {
        // 1) 如果 P,X 都空，就报告 R
        if (P.none()) {
            if ((int) R.count() >= minK) {
                report(R, pivots);
            }
            return;
        }
        const int need = minK - (int) R.count();
        if (need > 0 && need > (int) P.count()) {
            return;
        }
        // 2) 选 pivot u ∈ P∪X，使 |P ∧ nbr(u)| 最大
        int bestU = -1, bestCnt = -1;
        for_each_bit(P, n, [&](int u) {
            int cnt = adj[u].count_and(P); // 对应实现里做 word-wise popcount(a[i] & b[i])
            if (cnt > bestCnt) {
                bestCnt = cnt;
                bestU = u;
            }
            return true;
        });
        // Bitset candidates = ;
        // std::cout << "candidates: " ;
        // std::cout << "bestU: " << bestU << std::endl;
        // printBitset(P, "P");
        // printBitset(R, "R");
        // printBitset(adj[bestU], "adj[bestU]");
        // printBitset(candidates, "candidates");

        for_each_bit(P & ~adj[bestU], n, [&](int v) {
            // std::cout << v << std::endl;
            Bitset R2 = R;
            R2.set(v);
            Bitset P2 = P & adj[v];

            // 1) 先拷一份 pivots
            Bitset piv2 = pivots;
            if (v == bestU) piv2.set(v);

            // 3) 用拷贝去递归
            edgeSplit(adj, n, minK, R2, P2, piv2, report);

            P.reset(v);
            return true;
        });
    }

    // CSR-based pathSplit with shadow counters (pSize/pivSize)
    template<typename ReportFn>
    void pathSplit(int n, int r, int minK,
                   Bitset &P,
                   Bitset &pivots,
                   int &pSize, // |P|
                   int &pivSize, // |pivots|
                   daf::StaticVector<daf::Size> &conflictCount,
                   const daf::StaticVector<daf::Size> &conflictMaxSize,
                   const std::vector<uint32_t> &csOff,
                   const std::vector<VIdx> &csCol,
                   const std::vector<uint32_t> &rsOff,
                   const std::vector<uint32_t> &rsCol,
                   daf::Size nextCid,
                   const Bitset &emptyPivotsForReport,
                   ReportFn &&report) {
        // 剪枝：规模或 pivot 约束不满足
        if (pSize < minK || (pSize - pivSize) > minK) return;
        // std::cout << "pathSplit: pSize=" << pSize << ", pivSize=" << pivSize << std::endl;
        // 正好满足 |P|-|pivots| == minK → 只保留 P 中非 pivot 的作为 clique 报告
        // if ((pSize - pivSize) == minK) {
        //     report(P & (~pivots), emptyPivotsForReport);
        //     return;
        // }
        // 寻找第一个违反的冲突集

        size_t pick = conflictCount.size();
        for (size_t cid = nextCid; cid < conflictCount.size(); ++cid) {
            if (conflictCount[cid] >= conflictMaxSize[cid]) {
                pick = cid;
                break;
            }
        }

        if (pick < conflictCount.size()) {
            // 收集该冲突集中属于 pivot 的顶点（索引）
            // {
            //     int possibleGain = 0;
            //     for (uint32_t e = csOff[pick]; e < csOff[pick + 1]; ++e) {
            //         VIdx v = csCol[e];
            //         if (P.test(v) && !pivots.test(v)) ++possibleGain;

            //     }
            //     if (pSize + possibleGain < minK) {
            //         return;
            //     }
            // }

            VIdx pivBuf[400];
            int pc = 0;
            for (uint32_t e = csOff[pick]; e < csOff[pick + 1]; ++e) {
                VIdx v = csCol[e];
                if (pivots.test(v)) pivBuf[pc++] = v;
            }

            std::sort(pivBuf, pivBuf + pc, [&](VIdx a, VIdx b) {
                uint32_t da = s_rsOff[(size_t) a + 1] - s_rsOff[(size_t) a];
                uint32_t db = s_rsOff[(size_t) b + 1] - s_rsOff[(size_t) b];
                return da > db;
            });

            for (int i = 0; i < pc; ++i) {
                VIdx v = pivBuf[i];

                bool wasP = P.test(v);
                bool wasPi = pivots.test(v);
                if (wasP) {
                    P.reset(v);
                    --pSize;
                }
                if (wasPi) {
                    pivots.reset(v);
                    --pivSize;
                }

                // v 出现于哪些冲突集：rsOff/rsCol
                for (uint32_t e = rsOff[static_cast<size_t>(v)];
                     e < rsOff[static_cast<size_t>(v) + 1]; ++e) {
                    uint32_t g = rsCol[e];
                    --conflictCount[g];
                }

                // 与原逻辑一致：除了第一个，其它分支把上一个 pivot 从 pivots 转回 P
                if (i != 0) {
                    VIdx u = pivBuf[i - 1];
                    if (pivots.test(u)) {
                        pivots.reset(u);
                        --pivSize;
                    }
                    if (!P.test(u)) {
                        P.set(u);
                        ++pSize;
                    }
                }

                pathSplit(n, r, minK, P, pivots, pSize, pivSize,
                          conflictCount, conflictMaxSize,
                          csOff, csCol, rsOff, rsCol,
                          pick + 1, emptyPivotsForReport, report);

                // 回溯
                for (uint32_t e = rsOff[static_cast<size_t>(v)];
                     e < rsOff[static_cast<size_t>(v) + 1]; ++e) {
                    uint32_t g = rsCol[e];
                    ++conflictCount[g];
                }
                if (wasPi) {
                    pivots.set(v);
                    ++pivSize;
                }
                if (wasP) {
                    P.set(v);
                    ++pSize;
                }
            }

            // 复原：除了最后一个，前面的 pivot 需要设回 pivot（匹配原实现最后的 reset pivots）
            for (int i = 0; i + 1 < pc; ++i) {
                VIdx u = pivBuf[i];
                if (!pivots.test(u)) {
                    pivots.set(u);
                    ++pivSize;
                }
            }
            return;
        }

        // 没有任何违反冲突集，直接报告（保持原逻辑）
        if (pSize >= minK) {
            // 如果P-pivots=minK, 则只报告 P 中非 pivot 的部分
            if ((pSize - pivSize) == minK) {
                report(P & (~pivots), emptyPivotsForReport);
            } else { report(P, pivots); }
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
        assert(n <= static_cast<int>(std::numeric_limits<VIdx>::max()));

        // Bitset pivots(n);  pivots.reset();
        // Bitset P(n);       P.set();meis
        pivots.setSize(n);
        pivots.reset();
        P.setSize(n);
        P.set();

        // vListMap：把原节点 id 映射到 [0..n-1]
        for (int i = 0; i < n; ++i) {
            daf::vListMap[vList[i]] = i;
            if (vList[i].isPivot) pivots.set(i);
        }

        // 初始化计数
        daf::StaticVector<daf::Size> conflictCount, conflictMaxSize;
        conflictCount.resize(conflictSets.size());
        conflictMaxSize.resize(conflictSets.size());
        size_t G = conflictSets.size();
        size_t total = 0;
        for (size_t cid = 0; cid < G; ++cid) {
            daf::Size sz = static_cast<daf::Size>(conflictSets[cid].size());
            conflictCount[cid] = sz;
            conflictMaxSize[cid] = sz;
            total += (size_t) sz;
        }

        // ---------- 构建 CSR: 冲突集 -> 顶点 (索引) ----------
        s_csOff.resize(G + 1);
        s_csOff[0] = 0;
        s_csCol.resize(total);
        s_deg.assign(static_cast<size_t>(n), 0u);

        size_t pos = 0;
        for (size_t cid = 0; cid < G; ++cid) {
            s_csOff[cid] = static_cast<uint32_t>(pos);
            // 如果发现一个全都不是 pivots 的冲突集，说明无解，直接返回
            bool hasPiv = false;
            for (auto raw: conflictSets[cid]) {
                VIdx v = static_cast<VIdx>(daf::vListMap[raw]);
#ifndef NDEBUG
                if (v >= static_cast<VIdx>(n)) { std::abort(); }
#endif
                if (pivots.test(v)) {
                    hasPiv = true;
                }
                s_csCol[pos++] = v;
                ++s_deg[static_cast<size_t>(v)];
            }
            if (!hasPiv) {
                return;
            }
        }
        s_csOff[G] = static_cast<uint32_t>(pos);

        // 前缀和得到 rsOff
        s_rsOff.resize(static_cast<size_t>(n) + 1u);
        s_rsOff[0] = 0u;
        for (int v = 0; v < n; ++v) {
            s_rsOff[static_cast<size_t>(v) + 1u] = s_rsOff[static_cast<size_t>(v)] + s_deg[static_cast<size_t>(v)];
        }

        // 反向列数组，长度 = 所有出现次数之和
        s_rsCol.resize(s_rsOff.back());

        // 游标从 rsOff 拷贝一份出来
        s_cur = s_rsOff;

        for (size_t cid = 0; cid < G; ++cid) {
            const uint32_t begin = s_csOff[cid];
            const uint32_t end = s_csOff[cid + 1];
            for (uint32_t e = begin; e < end; ++e) {
                VIdx v = s_csCol[e];
                s_rsCol[s_cur[static_cast<size_t>(v)]++] = static_cast<uint32_t>(cid);
            }
        }

        // ---------- 预剪枝：若某顶点出现在所有需要的冲突集中，则剔除（与原逻辑一致） ----------
        for (int i = 0; i < n; ++i) {
            daf::Size maxRClique = nCr[n - 1][r - 1];
            if ((daf::Size) s_deg[static_cast<size_t>(i)] >= maxRClique) {
                if (P.test(i)) P.reset(i);
                if (P.count() < (size_t) minK) return; // 保持原检查
                if (pivots.test(i)) pivots.reset(i);
                else return;
                // 同步减少相关冲突集计数
                for (uint32_t e = s_rsOff[static_cast<size_t>(i)]; e < s_rsOff[static_cast<size_t>(i) + 1]; ++e) {
                    uint32_t cid = s_rsCol[e];
                    --conflictCount[cid];
                }
            }
        }

        // 影子计数，避免重复 count()
        int pSize = (int) P.count();
        int pivSize = (int) pivots.count();

        // 复用的“空 pivots”供某个报告分支使用，避免反复分配
        // Bitset emptyPiv(n); emptyPiv.reset();
        emptyPiv.setSize(n);
        emptyPiv.reset();

        // 进入递归
        pathSplit(n, r, minK, P, pivots, pSize, pivSize,
                  conflictCount, conflictMaxSize,
                  s_csOff, s_csCol, s_rsOff, s_rsCol,
                  0, emptyPiv, std::forward<ReportFn>(report));

        conflictCount.free();
        conflictMaxSize.free();
    }

    inline void testBronKerbosch() {
        // 3 个点完全图：0–1, 0–2, 1–2
        // 假设 TreeGraphNode 的构造函数为 TreeGraphNode(size_t id, bool flag)
        std::vector<TreeGraphNode> vList = {
            {1, false}, // 原来是 3
            {2, false}, // 原来是 4
            {3, true}, // 原来是 6
            {4, true}, // 原来是 7
            {5, true}, // 原来是 8
            {6, true}, // 原来是 9
        };

        std::vector<std::vector<daf::Size> > conflictSets;
        conflictSets.emplace_back(std::vector<daf::Size>{1, 2}); // {3,4,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{1, 3, 6});  // {3,6,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{4, 5, 6}); // {7,8,9}
        // conflictSets.emplace_back(std::vector<daf::Size>{3, 5, 6}); // {6,8,9}
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