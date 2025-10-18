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

extern double nCr[1001][401];

namespace bkRmEdge {
    // static constexpr int MAXN = 400;
    using Bitset = boost::dynamic_bitset<>;

    /**
     * ， n 
     */
    template<class F>
    void for_each_bit(const Bitset &bs, int n, F &&callback) {
        //  1
        for (size_t v = bs.find_first(); v != Bitset::npos && (int) v < n; v = bs.find_next(v)) {
            // bs.test(v)  true，
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

        // cover  1 
        auto i = cover.find_first();
        // pivots  1 
        auto pj = pivots.find_first();

        //  cover  1 ，
        while (i != Bitset::npos && i < vList.size()) {
            //  pj  >= i
            while (pj != Bitset::npos && pj < i) {
                pj = pivots.find_next(pj);
            }
            //  pj == i， pivot
            bool isP = (pj == i);

            // 
            result.emplace_back(vList[i].v, isP);

            //  cover  1 
            i = cover.find_next(i);
        }
        return result;
    }

    /**
     *  class  run() ，，
     *  adj, n, minK ，
     *  report 
     */
    template<class ReportFn>
    void bk_run(const std::vector<Bitset> &adj,
                int n,
                int minK,
                Bitset R,
                Bitset P,
                Bitset pivots,
                ReportFn &&report) {
        // 1)  P,X ， R
        // std::cout << R << " " << P << " " << pivots << std::endl;
        // 111111 000000 111110
        if (P.none()) {
            if ((int) R.count() >= minK) {
                report(R, pivots);
            }
            return;
        }

        // 2)  pivot u ∈ P∪X， |P ∧ nbr(u)| 
        int bestU = -1;
        int bestCnt = -1;
        const int Pc = (int)P.count();
        const int perfect = Pc - 1;
        //  scratch， bitset
        Bitset scratch(n);
        for_each_bit(P, n, [&](int u) {
            scratch = adj[u];   // 
            scratch &= P;       //  P 
            int cnt = (int)scratch.count();
            if (cnt > bestCnt) {
                bestCnt = cnt;
                bestU = u;
                if (bestCnt == perfect) return false; // ：
            }
            return true;
        });
        // ：candidates = P \ N(bestU)
        scratch = adj[bestU];
        scratch.flip();     // ~adj[bestU]
        scratch &= P;

        Bitset R2 = R;
        for_each_bit(scratch, n, [&](int v) {
            R2.set(v);
            Bitset P2 = P & adj[v];

            // ：/ pivot 。
            if (v == bestU) {
                pivots.set(v);
                bk_run(adj, n, minK, R2, P2, pivots, report);
                pivots.reset(v);
            } else {
                bk_run(adj, n, minK, R2, P2, pivots, report);
            }

            P.reset(v);
            R2.reset(v);
            return true;
        });
    }

    /**
     *  constructor ：
     *   -  vList
     *   -  vListMap
     *   - ， removeEdgeList 
     *  adj，  outN/outMinK  n/minK
     */
    inline std::vector<Bitset>
    build_adj(std::vector<TreeGraphNode> &vList,
              daf::StaticVector<std::pair<daf::Size, daf::Size> > &removeEdgeList,
              Bitset &staticVertex, // ： removeEdgeList 
              Bitset &pivot, // ： removeEdgeList 
              int &outN) {
        std::ranges::sort(vList);
        int n = (int) vList.size();
        outN = n;

        // 
        staticVertex.resize(n);
        staticVertex.set();

        pivot.resize(n);
        pivot.reset();

        // 
        std::vector<Bitset> adj(n, Bitset(n));
        for (int i = 0; i < n; ++i) {
            daf::vListMap[vList[i]] = i;
            adj[i].set();
            adj[i].reset(i);
            if (vList[i].isPivot) {
                pivot.set(i);
            }
        }
        //  removeEdgeList ， staticVertex 
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
     * ， main  new  BronKerbosch(...)  call run：
     *
     *   bronKerbosch(vList, removeEdgeList, minK, report);
     *
     *  report(Bitset clique) 。
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

        Bitset P = ~R;
        // X 


        // ， R
        // std::cout << "adj" << std::endl;
        // std::cout << adj << std::endl;
        // std::cout << R << " " << P << " " << povit << std::endl;
        bk_run(adj, n, minK, R, P, povit, std::forward<ReportFn>(report));
    }

    inline void testBronKerbosch() {
        // 3 ：0–1, 0–2, 1–2
        std::vector<TreeGraphNode> vList = {
            {0, false},
            {1, false},
            {2, true},
            {3, true},
            {4, true},
            {5, true}
        };
        daf::StaticVector<std::pair<daf::Size, daf::Size> > removeEdgeList;
        removeEdgeList.emplace_back(0, 1);
        // removeEdgeList.emplace_back(0, 2);
        int minK = 1;
        std::vector<double> cliqueCounts(vList.size(), 0);
        bronKerbosch(vList, removeEdgeList, minK,
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
        // 
        std::vector<Bitset> adj(n, Bitset(n));
        for (int i = 0; i < n; ++i) {
            adj[i].reset(); // 
        }
        int u, v;
        for (int i = 0; i < m; ++i) {
            fin >> u >> v;
            if (u >= 0 && u < n && v >= 0 && v < n) {
                adj[u].set(v);
                adj[v].set(u);
            }
        }

        //  R ，P  pivots  1
        Bitset R(n), P(n), pivots(n);
        R.reset();
        P.set(); //  P
        pivots = R;

        // std::cout << adj << std::endl;
        // 

        // std::cout << R << " " << P << " " << pivots << std::endl;
        bk_run(adj, n, minK, R, P, pivots, [&](const Bitset &clique, const Bitset &pivots) {
                         std::vector<int> C, P, H;
                         for (size_t i = clique.find_first(); i != Bitset::npos; i = clique.find_next(i)) {
                             C.push_back(i);
                             (pivots.test(i) ? P : H).push_back(i);
                         }
                         std::cout << "!!!!!!!!!!! Clique: ";
                         for (int x: C) std::cout << x << ' ';
                         std::cout << "| Pivots: ";
                         for (int x: P) std::cout << x << ' ';
                         std::cout << "| Holds: ";
                         for (int x: H) std::cout << x << ' ';
                         std::cout << '\n';

                     });
    }

    // ： k-clique 
    inline void testFromFile() {
        std::string filepath = "~/_/pivoter/b";
        auto minK = 1;
        std::vector<double> cliqueCounts;
        //  n  cliqueCounts 
        std::ifstream fin(filepath);
        int n, m;
        fin >> n >> m;
        cliqueCounts.assign(n + 1, 0.0);
        fin.close();

        bronKerboschFromFile(filepath, minK,
                             [&](const Bitset &clique, const Bitset &pivots) {
                                 //  hold / pivot 
                                 std::vector<int> H, P;
                                 // std::cout << clique << " " << pivots << std::endl;
                                 for (size_t i = clique.find_first(); i != Bitset::npos; i = clique.find_next(i)) {
                                     (pivots.test(i) ? P : H).push_back((int) i);
                                 }
                                 int h = (int) H.size();
                                 int p = (int) P.size();
                                 std::cout << "!!!!!!!!!!! Clique: ";
                                 std::cout << "Findclique: H: " << H << ", P: " << P << std::endl;
                                 //  q  0  p， h+q  clique
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
        // 
        // std::cout << "Clique counts (k: count):\n";
        // auto file = fopen("~/_/pivoter/b.out", "w");
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


#endif // BRONKERBOSCHFNS_HPP