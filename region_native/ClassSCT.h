// ClassSCT.h - orbit-aware class-weighted Succinct Clique Tree builder.
// Extracted from class_sct.cpp (validated: 20000+20000 random trials,
// sum-equality + per-pattern disjointness). buildClassSCT(G,k) returns
// CCPath leaves whose sum of support_count(leaf,0) == weighted k-clique
// count, leaves disjoint. Reuse for the region-native peel.
#ifndef CLASS_SCT_H
#define CLASS_SCT_H
#include "../src/NucleusDecomposition/CCPathCore.h"
#include <vector>
#include <algorithm>
#include <utility>
using ccpath::CCPath;
using ccpath::Vec;
struct ClassG {
    int C = 0;
    std::vector<int> w;                 // weight per class
    std::vector<std::vector<char>> A;   // A[i][j] adjacency (i!=j), A[i][i]=0
    bool adj(int i, int j) const { return A[i][j]; }
};

// An "open class": its class id and the residual weight available here.
struct OpenC { int c; int wres; };

// ======================================================================
//  buildClassSCT  — class-level Pivoter
//  ---------------------------------------------------------------------
//  EVERY leaf has empty hold (h = 0). All vertices a clique uses live in the
//  PIVOT pool, so the binomial weight C(w_c, m_c) is produced by
//  support_count's pivot term C(n_c, y_c) (support_count gives hold WEIGHT 1,
//  so weighted classes must never sit in hold).
//
//  A leaf is built from:
//    * a SPINE  : a list of (class, mult) that are FORCED. Each spine class c
//                 becomes a pivot with n_c = w_c and ell_c = u_c = mult,
//                 contributing the exact factor C(w_c, mult).
//    * a POOL   : a set of mutually-adjacent open classes (a clique), each a
//                 FREE pivot: n_c = wres, ell_c = 0, u_c = wres. They are
//                 chosen freely, summing (with the spine) to T = k.
//  T = k for every leaf (h = 0). A leaf contributes only if T is reachable.
//
//  RECURSION  gen(spine, spineSum, P):
//    spineSum = sum of spine mults. P = open classes (residual weights),
//    each adjacent to every spine class, sorted by class id.
//
//    If P is already a CLIQUE (pairwise adjacent), emit ONE leaf
//    (spine forced + P as a free pool). Otherwise pick a pivot pc in P with
//    the most P-neighbours and split, by the exact clique identity:
//
//      cliques(P) =  cliques using NO non-neighbour of pc
//                  ⊎  over each non-neighbour v of pc, in increasing class id,
//                     cliques whose LOWEST-id non-neighbour-of-pc used is v.
//
//    * "no non-neighbour of pc" : recurse on P' = {pc} ∪ N_P(pc) with the
//      SAME spine. (pc stays open; its non-neighbours are dropped.) This is
//      the compression branch — pc is now adjacent to all of P', so the next
//      pivot can only be pc itself => emit.
//    * "lowest non-neighbour used = v" : v contributes mult t in [1, wres_v].
//      Add (v, t) to the spine; new open set = classes adjacent to v that are
//      (a) still open, (b) NOT a lower-id non-neighbour of pc (canonical
//      lowest rule => those are handled by their own branch / the no-non-nbr
//      branch), and (c) v itself kept with residual wres_v - t (v can recur).
//
//  ORBIT (weight) SPLIT — the crux:
//    A class c may appear in the spine of a leaf with mult t while ALSO being
//    a free pool member (or spine again deeper) elsewhere. No vertex is ever
//    double-used: when we fix (v, t) we pass v forward with residual wres_v-t,
//    so any deeper selection from v draws from the REMAINING vertices. The
//    product of the spine factor C(w_v, t) at this level and a later factor
//    C(w_v - t, t') over the residual equals C(w_v, t) * C(w_v - t, t') =
//    C(w_v, t + t') * C(t + t', t) — i.e. it splits the (t+t')-subset of v's
//    vertices into the two roles, which is exactly the multiset the brute
//    force counts (a single class can supply many vertices to one clique;
//    they enter at different recursion depths but as one binomial overall).
//    Disjointness comes from the canonical lowest-non-neighbour branching.
// ======================================================================

// Build a leaf from a finalised spine (forced classes) + a free pool (a
// clique of universal pivots). Both are encoded as PIVOTS (empty hold) so
// the binomial weight is produced by support_count's pivot term.
//   * spine class c with mult m :  n_c = w_c, ell_c = u_c = m   (factor C(w_c,m))
//   * pool  class c             :  n_c = w_c, ell_c = 0, u_c = w_c (free)
// Each class is in AT MOST ONE of {spine, pool} for a given leaf, so there is
// no aliasing. T = k (hold is empty).
static void emitLeaf(const ClassG& G, const std::vector<std::pair<int,int>>& spine,
                     const std::vector<OpenC>& pool, int k,
                     std::vector<CCPath>& out) {
    const int C = G.C;
    Vec h = ccpath::zeros_vec(C);
    Vec n = ccpath::zeros_vec(C);
    Vec ell = ccpath::zeros_vec(C);
    Vec u = ccpath::zeros_vec(C);

    int forced = 0;
    for (auto& sp : spine) {
        int c = sp.first, mult = sp.second;
        n[c] = (int16_t)G.w[c];
        ell[c] = (int16_t)mult;
        u[c] = (int16_t)mult;
        forced += mult;
    }
    int poolCap = 0;
    for (auto& pc : pool) {
        int c = pc.c;
        n[c] = (int16_t)G.w[c];      // pc.wres == G.w[c] here (full residual)
        u[c] = (int16_t)G.w[c];      // ell stays 0
        poolCap += (int)G.w[c];
    }

    CCPath p;
    p.h = h; p.n = n; p.ell = ell; p.u = u;
    p.T = k;                          // h = 0 => T = k
    if (forced > k) return;           // spine alone overshoots
    if (forced + poolCap < k) return; // cannot reach k even filling the pool
    out.push_back(std::move(p));
}

// Recursion.
//   spine    : finalised (class, mult) decisions (forced pivots).
//   spineSum : sum of spine mults.
//   pool     : accumulated universal pivots (a clique); free pivots.
//   poolSum  : sum of pool weights.
//   P        : open classes, residual weights, each adjacent to all of pool
//              and to all spine classes; sorted by class id.
static void gen(const ClassG& G, std::vector<std::pair<int,int>>& spine,
                int spineSum, std::vector<OpenC>& pool, int poolSum,
                std::vector<OpenC> P, int k,
                std::vector<CCPath>& out) {
    if (spineSum > k) return;
    // Upper bound on reachable size: spine + pool + all open weights. If even
    // that can't reach k, prune.
    if (P.empty()) { emitLeaf(G, spine, pool, k, out); return; }

    // ---- pick pivot pc in P maximising residual weight of P adjacent to pc.
    int bestIdx = 0; long bestScore = -1;
    for (size_t i = 0; i < P.size(); ++i) {
        long sc = 0;
        for (size_t j = 0; j < P.size(); ++j)
            if (j != i && G.adj(P[i].c, P[j].c)) sc += P[j].wres;
        if (sc > bestScore) { bestScore = sc; bestIdx = (int)i; }
    }
    int pc = P[bestIdx].c;

    // Is pc universal in P (adjacent to every other open class)?
    bool universal = true;
    for (auto& q : P) if (q.c != pc && !G.adj(pc, q.c)) { universal = false; break; }

    if (universal) {
        // pc is adjacent to all remaining open classes => pull pc into the
        // free POOL and recurse on P \ {pc}. The pool stays a clique because
        // every later pool member comes from P \ {pc} ⊆ N(pc).
        std::vector<OpenC> Pr;
        Pr.reserve(P.size() - 1);
        for (auto& q : P) if (q.c != pc) Pr.push_back(q);
        pool.push_back(OpenC{pc, /*wres=*/G.w[pc]});
        gen(G, spine, spineSum, pool, poolSum + G.w[pc], std::move(Pr), k, out);
        pool.pop_back();
        return;
    }

    // pc has >=1 non-neighbour in P. Non-neighbours (EXCLUDING pc), by id.
    std::vector<int> nonNb; nonNb.reserve(P.size());   // <= |P| non-neighbours
    for (auto& q : P) if (q.c != pc && !G.adj(pc, q.c)) nonNb.push_back(q.c);
    std::sort(nonNb.begin(), nonNb.end());

    // ---------- branch A: use NO non-neighbour of pc ----------
    // Recurse on N[pc] = {pc} ∪ N_P(pc). Strictly smaller (>=1 dropped).
    {
        std::vector<OpenC> Pp;
        Pp.reserve(P.size());
        for (auto& q : P)
            if (q.c == pc || G.adj(pc, q.c)) Pp.push_back(q);
        gen(G, spine, spineSum, pool, poolSum, std::move(Pp), k, out);
    }

    // ---------- branch B: the LOWEST non-neighbour-of-pc used is v ----------
    // v's FULL multiplicity is decided here (mirrors countCanon picking class
    // v and its count j in one step). v is then removed from the open set;
    // child = open classes adjacent to v, after the canonical cut (drop pc's
    // non-neighbours with id <= v, since lower ones are other branches and v
    // itself is decided).
    for (size_t t = 0; t < nonNb.size(); ++t) {
        int v = nonNb[t];
        int wv = 0;
        for (auto& q : P) if (q.c == v) { wv = q.wres; break; }

        std::vector<char> drop(G.C, 0);
        for (size_t q = 0; q <= t; ++q) drop[nonNb[q]] = 1;  // <= : drop v too

        std::vector<OpenC> base;
        base.reserve(P.size());
        for (auto& q : P) {
            if (drop[q.c]) continue;             // v and lower pc-non-nbrs out
            if (!G.adj(v, q.c)) continue;        // must be adjacent to v
            base.push_back(q);
        }

        for (int mtt = 1; mtt <= wv; ++mtt) {
            if (spineSum + mtt > k) break;
            spine.push_back({v, mtt});
            std::vector<OpenC> child = base;     // already excludes v
            gen(G, spine, spineSum + mtt, pool, poolSum, std::move(child), k, out);
            spine.pop_back();
        }
    }
}

static std::vector<CCPath> buildClassSCT(const ClassG& G, int k) {
    std::vector<CCPath> out;
    std::vector<OpenC> P; P.reserve(G.C);   // <= G.C open classes
    for (int c = 0; c < G.C; ++c) if (G.w[c] > 0) P.push_back(OpenC{c, G.w[c]});
    std::sort(P.begin(), P.end(),
              [](const OpenC& a, const OpenC& b){ return a.c < b.c; });
    std::vector<std::pair<int,int>> spine;
    std::vector<OpenC> pool;
    gen(G, spine, 0, pool, 0, std::move(P), k, out);
    return out;
}
#endif // CLASS_SCT_H
