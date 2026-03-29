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
    static thread_local DynBitset R, P, pivots, emptyPiv, tmp2;
    // Reusable static buffers (thread-local for parallel use)
    static thread_local std::vector<uint32_t> s_csOff, s_rsOff, s_rsCol, s_deg, s_cur;
    static thread_local std::vector<VIdx> s_csCol;
    // Recursion depth counter (shared across all pathSplit template instantiations)
    static thread_local int pathSplitDepth_ = 0;
    /**
     * ， n
     */
    template<class F>
    bool for_each_bit(const Bitset &bs, int n, F &&callback) {
        //  1
        for (size_t v = bs.find_first(); v != Bitset::npos && (int) v < n; v = bs.find_next(v)) {
            // bs.test(v)  true，
            if (!callback((int) v)) {
                return false; //  false，
            }
        }

        return true; //  true
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


    template<class ReportFn>
    void edgeSplit(const std::vector<Bitset> &adj,
                   int n,
                   int minK,
                   Bitset R,
                   Bitset P,
                   Bitset pivots,
                   ReportFn &&report) {
        // 1)  P,X ， R
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
        // 2)  pivot u ∈ P∪X， |P ∧ nbr(u)|
        int bestU = -1, bestCnt = -1;
        for_each_bit(P, n, [&](int u) {
            int cnt = adj[u].count_and(P); //  word-wise popcount(a[i] & b[i])
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

            // 1)  pivots
            Bitset piv2 = pivots;
            if (v == bestU) piv2.set(v);

            // 3)
            edgeSplit(adj, n, minK, R2, P2, piv2, report);

            P.reset(v);
            return true;
        });
    }

    // =========================================================
    // Bipartite Matching LB (König's theorem)
    // Strictly stronger than disjoint packing: μ(G) ≥ packing(F)
    // Bipartite graph: L=active conflicts, R=pivot vertices, edge iff v∈F_i
    // Returns: max matching size = LB on minimum hitting set
    // =========================================================
    inline int computeMatchingLB(
        int n,
        const Bitset &pivots,
        const daf::StaticVector<daf::Size> &conflictCount,
        const daf::StaticVector<daf::Size> &conflictMaxSize,
        const std::vector<uint32_t> &csOff,
        const std::vector<VIdx> &csCol)
    {
        // Collect active conflicts and their pivot members
        // Build compact bipartite graph for Hopcroft-Karp
        int m = 0; // number of active conflicts (left vertices)
        int activeIds[400]; // map: compact left id → original cid
        for (size_t cid = 0; cid < conflictCount.size(); ++cid) {
            if (conflictCount[cid] < conflictMaxSize[cid]) continue;
            if (m >= 400) break;
            activeIds[m++] = (int)cid;
        }
        if (m == 0) return 0;

        // Collect pivot vertices (right vertices) with compact IDs
        int a = 0;
        int pivMap[400]; // vertex → compact pivot id (-1 if not pivot)
        memset(pivMap, -1, sizeof(pivMap));
        for (int v = 0; v < n && v < 400; v++) {
            if (pivots.test(v)) pivMap[v] = a++;
        }
        if (a == 0) return 0;

        // Hungarian/Hopcroft-Karp matching on small bipartite graph
        // Use simple augmenting path algorithm (sufficient for m,a ≤ 400)
        int matchL[400], matchR[400]; // matchL[left] = right, matchR[right] = left
        memset(matchL, -1, sizeof(matchL));
        memset(matchR, -1, sizeof(matchR));

        // For each left vertex, try to find augmenting path (DFS)
        // Adjacency: left i → right pivMap[v] for each pivot v in conflict activeIds[i]
        int mu = 0;
        bool visited[400];

        // DFS augmenting path from left vertex u
        std::function<bool(int)> tryAugment = [&](int u) -> bool {
            // Iterate pivot members of conflict activeIds[u]
            int cid = activeIds[u];
            for (uint32_t e = csOff[cid]; e < csOff[cid + 1]; ++e) {
                VIdx v = csCol[e];
                if (v >= 400 || pivMap[v] < 0) continue;
                int rv = pivMap[v]; // right vertex id
                if (visited[rv]) continue;
                visited[rv] = true;
                if (matchR[rv] < 0 || tryAugment(matchR[rv])) {
                    matchL[u] = rv;
                    matchR[rv] = u;
                    return true;
                }
            }
            return false;
        };

        for (int u = 0; u < m; u++) {
            memset(visited, false, a * sizeof(bool));
            if (tryAugment(u)) mu++;
        }

        return mu;
    }

    // =========================================================
    // Greedy disjoint packing lower bound (Theorem 4 from r3bk.md)
    // Returns: max number of pivot-disjoint active conflict sets
    // Uses tmp2 as scratch bitset (safe: not used concurrently)
    // =========================================================
    inline int computePackingLB(
        int n,
        const Bitset &pivots,
        const daf::StaticVector<daf::Size> &conflictCount,
        const daf::StaticVector<daf::Size> &conflictMaxSize,
        const std::vector<uint32_t> &csOff,
        const std::vector<VIdx> &csCol,
        daf::Size startCid)
    {
        tmp2.setSize(n);
        tmp2.reset(); // tracks pivots already claimed by a packing set
        int lb = 0;
        for (size_t cid = startCid; cid < conflictCount.size(); ++cid) {
            if (conflictCount[cid] < conflictMaxSize[cid]) continue; // resolved
            // Check: all pivot members unused, and at least 1 pivot member exists
            bool disjoint = true;
            int pivCount = 0;
            for (uint32_t e = csOff[cid]; e < csOff[cid + 1]; ++e) {
                VIdx v = csCol[e];
                if (pivots.test(v)) {
                    if (tmp2.test(v)) { disjoint = false; break; }
                    pivCount++;
                }
            }
            if (!disjoint || pivCount == 0) continue;
            // Pack: mark all pivot members as used
            for (uint32_t e = csOff[cid]; e < csOff[cid + 1]; ++e) {
                VIdx v = csCol[e];
                if (pivots.test(v)) tmp2.set(v);
            }
            lb++;
        }
        return lb;
    }

    // =========================================================
    // CSR-based pathSplit with Hitting-Set pruning
    //   - Packing lower bound (Thm 4): prune if LB > budget
    //   - Smallest active conflict first (improved branching)
    //   - Most-constrained pivot ordering within conflict
    // =========================================================
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
        pathSplitDepth_++;
        struct DepthGuard_ { ~DepthGuard_() { pathSplitDepth_--; } } _dg;
        if (pathSplitDepth_ > 500) return;

        // Basic pruning: size check
        if (pSize < minK || (pSize - pivSize) > minK) return;

        // --- Find first active conflict ---
        size_t pick = conflictCount.size();
        for (size_t cid = nextCid; cid < conflictCount.size(); ++cid) {
            if (conflictCount[cid] >= conflictMaxSize[cid]) { pick = cid; break; }
        }

        if (pick < conflictCount.size()) {
            // Note: Packing LB computed in pre-processing (removeRClique).
            // Per-node LB too expensive for dense instances (O(G*|F|) per call).

            // --- Branch on pivots in the chosen conflict ---
            VIdx pivBuf[400];
            int pc = 0;
            for (uint32_t e = csOff[pick]; e < csOff[pick + 1]; ++e) {
                VIdx v = csCol[e];
                if (pivots.test(v)) pivBuf[pc++] = v;
            }

            // Most-constrained-first: sort by descending conflict degree
            std::sort(pivBuf, pivBuf + pc, [&](VIdx a, VIdx b) {
                uint32_t da = s_rsOff[(size_t) a + 1] - s_rsOff[(size_t) a];
                uint32_t db = s_rsOff[(size_t) b + 1] - s_rsOff[(size_t) b];
                return da > db;
            });

            for (int i = 0; i < pc; ++i) {
                VIdx v = pivBuf[i];

                bool wasP = P.test(v);
                bool wasPi = pivots.test(v);
                if (wasP) { P.reset(v); --pSize; }
                if (wasPi) { pivots.reset(v); --pivSize; }

                for (uint32_t e = rsOff[static_cast<size_t>(v)];
                     e < rsOff[static_cast<size_t>(v) + 1]; ++e)
                    --conflictCount[rsCol[e]];

                // Promote previous pivot to keep
                if (i != 0) {
                    VIdx u = pivBuf[i - 1];
                    if (pivots.test(u)) { pivots.reset(u); --pivSize; }
                    if (!P.test(u)) { P.set(u); ++pSize; }
                }

                // --- Recursive unit propagation ---
                // After removing v (and promoting prev), check for new singletons.
                // Queue-based: only inspect conflict sets affected by removals.
                VIdx forcedBuf[400];
                int nForced = 0;
                bool dead = false;
                {
                    // Seed queue: conflict sets affected by removing v
                    uint32_t queue[2048];
                    int qHead = 0, qTail = 0;
                    for (uint32_t e = rsOff[static_cast<size_t>(v)];
                         e < rsOff[static_cast<size_t>(v) + 1]; ++e)
                        if (qTail < 2048) queue[qTail++] = rsCol[e];

                    while (qHead < qTail && !dead) {
                        uint32_t g = queue[qHead++];
                        if (conflictCount[g] < conflictMaxSize[g]) continue;
                        int pivCnt = 0; VIdx single = 0;
                        for (uint32_t e2 = csOff[g]; e2 < csOff[g + 1]; ++e2) {
                            VIdx u = csCol[e2];
                            if (pivots.test(u)) { pivCnt++; single = u; }
                        }
                        if (pivCnt == 0) { dead = true; break; }
                        if (pivCnt == 1) {
                            P.reset(single); pivots.reset(single);
                            --pSize; --pivSize;
                            forcedBuf[nForced++] = single;
                            for (uint32_t e2 = rsOff[static_cast<size_t>(single)];
                                 e2 < rsOff[static_cast<size_t>(single) + 1]; ++e2) {
                                --conflictCount[rsCol[e2]];
                                if (qTail < 2048) queue[qTail++] = rsCol[e2];
                            }
                            if (pSize < minK) { dead = true; break; }
                        }
                    }
                }

                if (!dead) {
                    pathSplit(n, r, minK, P, pivots, pSize, pivSize,
                              conflictCount, conflictMaxSize,
                              csOff, csCol, rsOff, rsCol,
                              pick + 1, emptyPivotsForReport, report);
                }

                // Undo forced removals (reverse order)
                for (int fi = nForced - 1; fi >= 0; --fi) {
                    VIdx u = forcedBuf[fi];
                    P.set(u); pivots.set(u);
                    ++pSize; ++pivSize;
                    for (uint32_t e = rsOff[static_cast<size_t>(u)];
                         e < rsOff[static_cast<size_t>(u) + 1]; ++e)
                        ++conflictCount[rsCol[e]];
                }

                // Backtrack v
                for (uint32_t e = rsOff[static_cast<size_t>(v)];
                     e < rsOff[static_cast<size_t>(v) + 1]; ++e)
                    ++conflictCount[rsCol[e]];
                if (wasPi) { pivots.set(v); ++pivSize; }
                if (wasP) { P.set(v); ++pSize; }
            }

            // Restore promoted pivots
            for (int i = 0; i + 1 < pc; ++i) {
                VIdx u = pivBuf[i];
                if (!pivots.test(u)) { pivots.set(u); ++pivSize; }
            }
            return;
        }

        // Terminal: no active conflict — report clique
        if (pSize >= minK) {
            if ((pSize - pivSize) == minK)
                report(P & (~pivots), emptyPivotsForReport);
            else
                report(P, pivots);
        }
    }

    /**
     * ， main  new  BronKerbosch(...)  call run：
     *
     *   bronKerbosch(vList, removeEdgeList, minK, report);
     *
     *  report(Bitset clique) 。
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
                       ReportFn &&report,
                       daf::StaticVector<daf::Size>* vertexToIndex = nullptr) {
        const int n = static_cast<int>(vList.size());
        assert(n <= static_cast<int>(std::numeric_limits<VIdx>::max()));

        pivots.setSize(n);
        pivots.reset();
        P.setSize(n);
        P.set();

        // Vertex id -> index [0..n-1]; use thread-local map when provided for parallel use
        daf::StaticVector<daf::Size>* const mapPtr = vertexToIndex ? vertexToIndex : &daf::vListMap;
        for (int i = 0; i < n; ++i) {
            (*mapPtr)[vList[i].v] = static_cast<daf::Size>(i);
            if (vList[i].isPivot) pivots.set(i);
        }

        //
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

        // ----------  CSR:  ->  () ----------
        s_csOff.resize(G + 1);
        s_csOff[0] = 0;
        s_csCol.resize(total);
        s_deg.assign(static_cast<size_t>(n), 0u);

        size_t pos = 0;
        for (size_t cid = 0; cid < G; ++cid) {
            s_csOff[cid] = static_cast<uint32_t>(pos);
            //  pivots ，，
            bool hasPiv = false;
            for (auto raw: conflictSets[cid]) {
                VIdx v = static_cast<VIdx>((*mapPtr)[raw]);
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

        //  rsOff
        s_rsOff.resize(static_cast<size_t>(n) + 1u);
        s_rsOff[0] = 0u;
        for (int v = 0; v < n; ++v) {
            s_rsOff[static_cast<size_t>(v) + 1u] = s_rsOff[static_cast<size_t>(v)] + s_deg[static_cast<size_t>(v)];
        }

        // ， =
        s_rsCol.resize(s_rsOff.back());

        //  rsOff
        s_cur = s_rsOff;

        for (size_t cid = 0; cid < G; ++cid) {
            const uint32_t begin = s_csOff[cid];
            const uint32_t end = s_csOff[cid + 1];
            for (uint32_t e = begin; e < end; ++e) {
                VIdx v = s_csCol[e];
                s_rsCol[s_cur[static_cast<size_t>(v)]++] = static_cast<uint32_t>(cid);
            }
        }

        // ---------- ：，（） ----------
        for (int i = 0; i < n; ++i) {
            daf::Size maxRClique = nCr[n - 1][r - 1];
            if ((daf::Size) s_deg[static_cast<size_t>(i)] >= maxRClique) {
                if (P.test(i)) P.reset(i);
                if (P.count() < (size_t) minK) return; //
                if (pivots.test(i)) pivots.reset(i);
                else return;
                //
                for (uint32_t e = s_rsOff[static_cast<size_t>(i)]; e < s_rsOff[static_cast<size_t>(i) + 1]; ++e) {
                    uint32_t cid = s_rsCol[e];
                    --conflictCount[cid];
                }
            }
        }

        // ========== Singleton propagation (Theorem 3) ==========
        // If a conflict set has exactly 1 pivot member, that pivot MUST be removed.
        {
            bool changed = true;
            while (changed) {
                changed = false;
                for (size_t cid = 0; cid < G; ++cid) {
                    if (conflictCount[cid] < conflictMaxSize[cid]) continue; // resolved
                    int pivCount = 0;
                    VIdx singlePiv = 0;
                    for (uint32_t e = s_csOff[cid]; e < s_csOff[cid + 1]; ++e) {
                        VIdx v = s_csCol[e];
                        if (P.test(v) && pivots.test(v)) { pivCount++; singlePiv = v; }
                    }
                    if (pivCount == 0) return; // no pivot in active conflict → leaf dead
                    if (pivCount == 1) {
                        // Forced removal
                        P.reset(singlePiv);
                        pivots.reset(singlePiv);
                        for (uint32_t e = s_rsOff[singlePiv]; e < s_rsOff[singlePiv + 1]; ++e)
                            --conflictCount[s_rsCol[e]];
                        changed = true;
                        if (P.count() < (size_t)minK) return;
                    }
                }
            }
        }

        // ========== Subsumption elimination (Theorem 2) ==========
        // If F_i's pivot-members ⊆ F_j's pivot-members, F_j is redundant.
        // Only run when G is manageable to avoid O(G^2) blowup.
        if (G <= 500) {
            // Build pivot-member bitmasks for each active conflict set
            struct PivMask { size_t cid; Bitset mask; int size; };
            std::vector<PivMask> activeSets;
            activeSets.reserve(G);
            for (size_t cid = 0; cid < G; ++cid) {
                if (conflictCount[cid] < conflictMaxSize[cid]) continue;
                PivMask pm;
                pm.cid = cid;
                pm.mask.setSize(n);
                pm.mask.reset();
                pm.size = 0;
                for (uint32_t e = s_csOff[cid]; e < s_csOff[cid + 1]; ++e) {
                    VIdx v = s_csCol[e];
                    if (pivots.test(v)) { pm.mask.set(v); pm.size++; }
                }
                activeSets.push_back(std::move(pm));
            }
            // Sort by size ascending: smaller sets more likely to subsume larger ones
            std::sort(activeSets.begin(), activeSets.end(),
                      [](const PivMask &a, const PivMask &b) { return a.size < b.size; });
            // Mark subsumed sets
            for (size_t i = 0; i < activeSets.size(); ++i) {
                if (conflictCount[activeSets[i].cid] < conflictMaxSize[activeSets[i].cid]) continue;
                for (size_t j = i + 1; j < activeSets.size(); ++j) {
                    if (conflictCount[activeSets[j].cid] < conflictMaxSize[activeSets[j].cid]) continue;
                    // Check if F_i ⊆ F_j: every pivot in F_i is also in F_j
                    if (activeSets[i].size <= activeSets[j].size) {
                        bool subset = true;
                        activeSets[i].mask.for_each_bit([&](size_t bit) {
                            if (!activeSets[j].mask.test(bit)) subset = false;
                        });
                        if (subset) {
                            // F_i ⊆ F_j → F_j subsumed
                            // Set maxSize to count+1 so the set can never become active:
                            // Even if backtracking restores count, it restores to
                            // the pre-processing value which is < the inflated maxSize.
                            size_t jcid = activeSets[j].cid;
                            conflictMaxSize[jcid] = conflictCount[jcid] + 1;
                        }
                    }
                }
            }
        }

        // ========== Component decomposition LB (Theorem 5 + Corollary 5.1) ==========
        // Union-Find on pivots: pivots in same conflict set → same component.
        // Per-component packing LB sum ≥ global packing LB.
        int pSize = (int) P.count();
        int pivSize = (int) pivots.count();
        {
            int budget = pSize - minK;
            if (budget < 0) return;

            // Count active conflicts
            int numActiveConflicts = 0;
            for (size_t cid = 0; cid < G; ++cid)
                if (conflictCount[cid] >= conflictMaxSize[cid]) numActiveConflicts++;

            if (numActiveConflicts > 0) {
                // Union-Find for pivot components
                int uf[400];
                for (int i = 0; i < n; ++i) uf[i] = i;
                // Find with path compression (iterative)
                auto ufFind = [&](int x) {
                    while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
                    return x;
                };

                // Build components: union all pivots within each active conflict set
                for (size_t cid = 0; cid < G; ++cid) {
                    if (conflictCount[cid] < conflictMaxSize[cid]) continue;
                    int firstPiv = -1;
                    for (uint32_t e = s_csOff[cid]; e < s_csOff[cid + 1]; ++e) {
                        VIdx v = s_csCol[e];
                        if (pivots.test(v)) {
                            if (firstPiv >= 0) {
                                int ra = ufFind(firstPiv), rb = ufFind((int)v);
                                if (ra != rb) uf[ra] = rb;
                            } else firstPiv = (int)v;
                        }
                    }
                }

                // Per-component greedy disjoint packing LB
                int compLB[400] = {};
                tmp2.setSize(n);
                tmp2.reset(); // tracks used pivots globally

                for (size_t cid = 0; cid < G; ++cid) {
                    if (conflictCount[cid] < conflictMaxSize[cid]) continue;
                    bool disjoint = true;
                    int pivCount = 0;
                    int compRoot = -1;
                    for (uint32_t e = s_csOff[cid]; e < s_csOff[cid + 1]; ++e) {
                        VIdx v = s_csCol[e];
                        if (pivots.test(v)) {
                            if (tmp2.test(v)) { disjoint = false; break; }
                            pivCount++;
                            compRoot = ufFind((int)v);
                        }
                    }
                    if (!disjoint || pivCount == 0 || compRoot < 0) continue;
                    for (uint32_t e = s_csOff[cid]; e < s_csOff[cid + 1]; ++e) {
                        VIdx v = s_csCol[e];
                        if (pivots.test(v)) tmp2.set(v);
                    }
                    compLB[compRoot]++;
                }

                // Sum per-component LBs (Corollary 5.1)
                int totalLB = 0;
                for (int i = 0; i < n; ++i)
                    if (uf[i] == i && compLB[i] > 0) totalLB += compLB[i];

                if (totalLB > budget) return; // leaf dead

                // Note: Bipartite matching gives τ_cover (min vertex cover), NOT τ_hitting_set.
                // τ_hitting_set ≤ τ_cover = μ(matching). So matching is an UPPER bound, not lower.
                // Cannot use matching for dead-leaf detection. Disjoint packing remains the best LB.
            }
        }

        emptyPiv.setSize(n);
        emptyPiv.reset();

        //
        pathSplit(n, r, minK, P, pivots, pSize, pivSize,
                  conflictCount, conflictMaxSize,
                  s_csOff, s_csCol, s_rsOff, s_rsCol,
                  0, emptyPiv, std::forward<ReportFn>(report));

        conflictCount.free();
        conflictMaxSize.free();
    }

    inline void testBronKerbosch() {
        // 3 ：0–1, 0–2, 1–2
        //  TreeGraphNode  TreeGraphNode(size_t id, bool flag)
        std::vector<TreeGraphNode> vList = {
            {1, false}, //  3
            {2, false}, //  4
            {3, true}, //  6
            {4, true}, //  7
            {5, true}, //  8
            {6, true}, //  9
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
        int minK = 4; //  clique
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