// SDCT.hpp
#ifndef SDCT_HPP
#define SDCT_HPP

#include <algorithm>
#include <algorithm>
#include <bitset>
#include <vector>
#include <ranges>
#include "Global/Global.h"
#include "graph/DynamicGraph.h"
#include <boost/dynamic_bitset.hpp>

extern double nCr[1001][401];

namespace SDCT {
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
        Bitset candidates = P & ~adj[bestU];
        // std::cout << "candidates: " ;
        // // 3) 遍历 P \ nbr(bestU)

        // std::cout << "bestU: " << bestU << std::endl;
        // printBitset(P, "P");
        // printBitset(adj[bestU], "adj[bestU]");
        // printBitset(candidates, "candidates");


        for_each_bit(candidates, n, [&](int v) {
            // std::cout << v << std::endl;
            Bitset R2 = R;
            R2.set(v);
            Bitset P2 = P & adj[v];

            // 1) 先拷一份 pivots
            Bitset piv2 = pivots;
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
              daf::StaticVector<std::pair<daf::Size, daf::Size> > &removeEdgeList,
              Bitset &staticVertex, // 输出：那些从未在 removeEdgeList 出现过的点
              Bitset &pivot, // 输出：那些从未在 removeEdgeList 出现过的点
              int &outN) {
        std::ranges::sort(vList);
        int n = (int) vList.size();
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
        for (auto [u0, v0]: removeEdgeList) {
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
                      daf::StaticVector<std::pair<daf::Size, daf::Size> > &removeEdgeList,
                      int minK,
                      ReportFn &&report) {
        int n;
        Bitset R;
        Bitset povit;
        auto adj = build_adj(vList, removeEdgeList, R, povit, n);

        Bitset mask(n);
        mask.set();
        Bitset P = (~R) & mask;
        // X 还是空
        Bitset X(n);
        // 运行原来的递归，只不过带了预先的 R
        std::cout << R << " " << P << " " << povit << std::endl;
        bk_run(adj, n, minK, R, P, povit, std::forward<ReportFn>(report));
    }


    template<class ReportFn>
    void bronKerboschFromFile(const std::string &filepath,
                              int minK,
                              ReportFn &&report) {
        std::ifstream fin(filepath);
        if (!fin) {
            std::cerr << "Error: cannot open file " << filepath << std::endl;
            return;
        }
        int n, m;
        fin >> n >> m;
        // 构建邻接矩阵
        std::vector<Bitset> adj(n, Bitset(n));
        for (int i = 0; i < n; ++i) {
            adj[i].reset(); // 清空
        }
        int u, v;
        for (int i = 0; i < m; ++i) {
            fin >> u >> v;
            if (u >= 0 && u < n && v >= 0 && v < n) {
                adj[u].set(v);
                adj[v].set(u);
            }
        }

        // 初始 R 空，P 和 pivots 全 1
        Bitset R(n), P(n), pivots(n);
        R.reset();
        P.set(); // 所有顶点都在 P
        pivots = R;

        // std::cout << adj << std::endl;
        // 调用主算法

        // std::cout << R << " " << P << " " << pivots << std::endl;
        bk_run(adj, n, minK, R, P, pivots, std::forward<ReportFn>(report));
    }

    // 测试函数：从文件读取图并统计所有 k-clique 数量
    inline void testFromFile() {
        std::string filepath = "/Users/zhangwenqian/UNSW/pivoter/b";
        auto minK = 1;
        std::vector<double> cliqueCounts;
        // 先读取 n 以初始化 cliqueCounts 大小
        std::ifstream fin(filepath);
        int n, m;
        fin >> n >> m;
        cliqueCounts.assign(n + 1, 0.0);
        fin.close();

        bronKerboschFromFile(filepath, minK,
                             [&](const Bitset &clique, const Bitset &pivots) {
                                 // 计算 hold / pivot 集合
                                 std::vector<int> H, P;
                                 // std::cout << clique << " " << pivots << std::endl;
                                 for (size_t i = clique.find_first(); i != Bitset::npos; i = clique.find_next(i)) {
                                     (pivots.test(i) ? P : H).push_back((int) i);
                                 }
                                 int h = (int) H.size();
                                 int p = (int) P.size();
                                 std::cout << "!!!!!!!!!!! Clique: ";
                                 std::cout << "Findclique: H: " << H << ", P: " << P << std::endl;
                                 // 对任意 q 从 0 到 p，都会产生 h+q 大小的 clique
                                 for (int q = 0; q <= p; ++q) {
                                     int size_k = h + q;
                                     if (size_k <= n && size_k >= 0) {
                                         cliqueCounts[size_k] += nCr[p][q];
                                     }
                                 }
                             }
        );
        // std::cout << cliqueCounts << std::endl;
        for (size_t k = minK; k < cliqueCounts.size(); ++k) {
            if (cliqueCounts[k] > 0) {
                // std::cout << cliqueCounts[k] << std::endl; //output as int
                printf("%.0f\n", cliqueCounts[k]);
            }
        }
        // 输出结果
        // std::cout << "Clique counts (k: count):\n";
        // auto file = fopen("/Users/zhangwenqian/UNSW/pivoter/b.out", "w");
        // for (size_t k = minK; k < cliqueCounts.size(); ++k) {
        //     double cnt = cliqueCounts[k];
        //     if (cnt > 0) {
        //         // std::cout << cnt;
        //         fprintf(file, "%.0f\n", cnt);
        //     }
        // }
        //
        // fclose(file);
    }
} // namespace bk


#endif // SDCT_HPP
