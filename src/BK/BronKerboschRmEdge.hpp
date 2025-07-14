// BronKerboschFns.hpp
#ifndef BRONKERBOSCHFNS_HPP
#define BRONKERBOSCHFNS_HPP

#include <algorithm>
#include <algorithm>
#include <bitset>
#include <vector>
#include <ranges>
#include "Global/Global.h"
#include "graph/DynamicGraph.h"
#include <boost/dynamic_bitset.hpp>


namespace bkRmEdge {
    // static constexpr int MAXN = 400;
    using Bitset = boost::dynamic_bitset<>;

    /**
     * 和原来成员版本一模一样，只是把 n 也作为参数传进来
     */
    template<class F>
    void for_each_bit(const Bitset &bs, int n, F &&callback) {
        // 首先找第一个 1
        for (size_t v = bs.find_first(); v != Bitset::npos && (int)v < n; v = bs.find_next(v)) {
            // bs.test(v) 肯定为 true，不用再测
            if (!callback((int)v)) break;
        }
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
    /**
     * 原来 class 中的 run() 方法，完全照搬逻辑，
     * 把 adj, n, minK 从成员变量变成了入参，
     * 把 report 从类成员变成了回调参数
     */
    template<class ReportFn>
    void bk_run(const std::vector<Bitset> &adj,
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

        // 2) 选 pivot u ∈ P∪X，使 |P ∧ nbr(u)| 最大
        int bestU = -1, bestCnt = -1;
        Bitset PX = P;
        for_each_bit(PX, n, [&](int u) {
            Bitset nbr = adj[u] & P;
            int cnt = (int) nbr.count();
            if (cnt > bestCnt) {
                bestCnt = cnt;
                bestU = u;
            }
            return true;
        });
        // std::cout << "bestU: " << bestU << std::endl;
        // 3) 遍历 P \ nbr(bestU)
        Bitset candidates = P & ~adj[bestU];
        for_each_bit(candidates, n, [&](int v) {
            Bitset R2 = R;    R2.set(v);
            Bitset P2 = P & adj[v];

            // 1) 先拷一份 pivots
            Bitset piv2 = pivots;
            // 2) 如果这个 v 是 pivot，就在这份 piv2 上打标
            if (v == bestU) piv2.set(v);

            // 3) 用拷贝去递归
            bk_run(adj, n, minK, R2, P2, piv2, report);

            P.reset(v);
            return true;
        });
    }

    /**
     * 等同于原来 constructor 的逻辑：
     *   - 排序 vList
     *   - 填 vListMap
     *   - 全连，再把 removeEdgeList 中的边删掉
     * 返回构造好的 adj， 并通过 outN/outMinK 传回 n/minK
     */
    inline std::vector<Bitset>
        build_adj(std::vector<TreeGraphNode> &vList,
                  daf::StaticVector<std::pair<daf::Size, daf::Size>> &removeEdgeList,
                  Bitset &staticVertex,  // 输出：那些从未在 removeEdgeList 出现过的点
                  Bitset &pivot,  // 输出：那些从未在 removeEdgeList 出现过的点
                  int &outN) {
        std::ranges::sort(vList);
        int n = (int)vList.size();
        outN = n;

        // 一开始假设所有点都是静态点
        staticVertex.resize(n);
        staticVertex.set();

        pivot.resize(n);
        pivot.reset();

        // 全连
        std::vector<Bitset> adj(n, Bitset(n));
        for (int i = 0; i < n; ++i) {
            daf::vListMap[vList[i]] = i;
            adj[i].set();
            adj[i].reset(i);
            if (vList[i].isPivot) {
                pivot.set(i);
            }
        }
        // 删除 removeEdgeList 中的边，同时把它们的两端从 staticVertex 中踢掉
        for (auto [u0, v0] : removeEdgeList) {
            if (daf::vListMap[u0] == std::numeric_limits<daf::Size>::max() ||
                    daf::vListMap[v0] == std::numeric_limits<daf::Size>::max()) {
                continue;
            }
            auto u = daf::vListMap[u0];
            auto v = daf::vListMap[v0];
            adj[u].reset(v);
            adj[v].reset(u);
            staticVertex.reset(u);
            staticVertex.reset(v);
            pivot.reset(u);
            pivot.reset(v);
        }

        for (int i = 0; i < n; ++i) {
            daf::vListMap[vList[i]] = std::numeric_limits<daf::Size>::max();
        }
        return adj;
    }


    /**
     * 对外主入口，等同于原来你在 main 里 new 一个 BronKerbosch(...) 然后 call run：
     *
     *   bronKerbosch(vList, removeEdgeList, minK, report);
     *
     * 其中 report(Bitset clique) 会在每次找到一个极大团时被调用。
     */
    template<class ReportFn>
    void bronKerbosch(std::vector<TreeGraphNode> &vList,
                      daf::StaticVector<std::pair<daf::Size, daf::Size>> &removeEdgeList,
                      int minK,
                      ReportFn &&report) {
        int n;
        Bitset R;
        Bitset povit;
        auto adj = build_adj(vList, removeEdgeList, R, povit, n);

        Bitset mask(n); mask.set();
        Bitset P = (~R) & mask;
        // X 还是空
        Bitset X(n);
        // 运行原来的递归，只不过带了预先的 R
        bk_run(adj, n, minK, R, P, povit, std::forward<ReportFn>(report));
    }

} // namespace bk

#endif // BRONKERBOSCHFNS_HPP
