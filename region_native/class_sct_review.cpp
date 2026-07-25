// class_sct.cpp — build a Succinct Clique Tree (SCT) of a small WEIGHTED
// class graph as a set of DISJOINT ccpath::CCPath leaves, such that
//
//   sum_leaves support_count(leaf, b=0, nCr)  ==  total weighted-k-clique count
//
// and the leaves partition the weighted k-cliques (each represented once).
//
// A "class" c carries weight w_c >= 1 (it stands for w_c interchangeable,
// mutually-adjacent vertices). adj(i,j) is the adjacency predicate BETWEEN
// distinct classes. A weighted k-clique = a clique S of classes plus
// multiplicities m_c>=1 (sum=k); its weight = prod_c C(w_c, m_c).
//
// This is the vertex-level Pivoter (degeneracy/pivot on non-neighbours)
// LIFTED to weighted classes. See buildClassSCT for the orbit handling.
//
// Build (from region_native/):
//   g++ -O3 -std=c++17 -I../src/NucleusDecomposition -o class_sct class_sct.cpp
// Run:  ./class_sct

#include "../src/NucleusDecomposition/CCPathCore.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <random>
#include <vector>

using ccpath::CCPath;
using ccpath::Vec;

// ---------------- Pascal-table nCr (double), n up to 64 ----------------
static double PASCAL[65][65];
static void initPascal() {
    for (int n = 0; n <= 64; ++n) {
        PASCAL[n][0] = 1.0;
        for (int k = 1; k <= n; ++k)
            PASCAL[n][k] = PASCAL[n - 1][k - 1] + (k <= n - 1 ? PASCAL[n - 1][k] : 0.0);
    }
}
static double nCr_fn(int n, int k) {
    if (k < 0 || n < 0 || k > n) return 0.0;
    return PASCAL[n][k];
}

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
    std::vector<int> nonNb;
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
    std::vector<OpenC> P;
    for (int c = 0; c < G.C; ++c) if (G.w[c] > 0) P.push_back(OpenC{c, G.w[c]});
    std::sort(P.begin(), P.end(),
              [](const OpenC& a, const OpenC& b){ return a.c < b.c; });
    std::vector<std::pair<int,int>> spine;
    std::vector<OpenC> pool;
    gen(G, spine, 0, pool, 0, std::move(P), k, out);
    return out;
}

// ---------------- brute-force ground truth (countCanon) ----------------
static const ClassG* gG = nullptr;
static double countCanon(const std::vector<int>& avail, int need) {
    if (need == 0) return 1.0;
    double res = 0.0;
    for (size_t i = 0; i < avail.size(); ++i) {
        int c = avail[i];
        std::vector<int> later;
        for (size_t t = i + 1; t < avail.size(); ++t)
            if (gG->adj(c, avail[t])) later.push_back(avail[t]);
        int hi = std::min(need, gG->w[c]);
        for (int j = 1; j <= hi; ++j)
            res += nCr_fn(gG->w[c], j) * countCanon(later, need - j);
    }
    return res;
}
static double bruteTotal(const ClassG& G, int k) {
    gG = &G;
    std::vector<int> all;
    for (int c = 0; c < G.C; ++c) if (G.w[c] > 0) all.push_back(c);
    return countCanon(all, k);
}

static double sctTotal(const std::vector<CCPath>& leaves) {
    double tot = 0.0;
    for (const auto& L : leaves) {
        Vec b = ccpath::zeros_vec(L.m());
        tot += ccpath::support_count(L, b, nCr_fn);
    }
    return tot;
}

// ---------------- PER-PATTERN DISJOINTNESS verifier ----------------
// For every class-multiplicity vector m of size k whose support is a clique
// (a "weighted-clique pattern"), the SCT must represent it EXACTLY ONCE:
//   sum over leaves of [m is an admissible pivot-count y of the leaf]
//                      * prod_c C(n_c, m_c)            (the realizations)
//   ==  prod_c C(w_c, m_c)                              (the true realizations)
// Equality for EVERY pattern => no double-representation and no omission
// (this is strictly stronger than the global sum check). A leaf with empty
// forbidden admits y = m iff ell_c <= m_c <= u_c for all c and sum(m)=T.
//
// Returns true if every pattern matches.
static bool checkDisjoint(const ClassG& G, int k,
                          const std::vector<CCPath>& leaves,
                          std::vector<int>& badPattern /*out*/) {
    const int C = G.C;
    std::vector<int> m(C, 0);
    bool ok = true;

    // enumerate clique-supported multiplicity vectors summing to k.
    // recursive over classes; maintain "chosen so far form a clique".
    std::vector<int> chosen;  // classes with m_c>=1, kept as a growing clique
    std::function<void(int,int)> rec = [&](int c, int rem) {
        if (rem == 0) {
            // true weight
            double truth = 1.0;
            for (int x : chosen) truth *= nCr_fn(G.w[x], m[x]);
            // SCT represented weight for this exact pattern
            double rep = 0.0;
            for (const auto& L : leaves) {
                if (L.T != k) continue;
                // check m is admissible: ell<=m<=u for all classes, and any
                // class with m_c>0 must be in the leaf's support (n_c>0).
                bool adm = true;
                int sm = 0;
                for (int cc = 0; cc < C; ++cc) {
                    int mc = m[cc];
                    sm += mc;
                    if (mc < (int)L.ell[cc] || mc > (int)L.u[cc]) { adm = false; break; }
                }
                if (!adm || sm != k) continue;
                double prod = 1.0;
                for (int cc = 0; cc < C; ++cc)
                    if (m[cc] > 0) prod *= nCr_fn((int)L.n[cc], m[cc]);
                rep += prod;
            }
            if (std::abs(rep - truth) > 0.5) {
                ok = false;
                if (badPattern.empty()) badPattern = m;
            }
            return;
        }
        if (c >= C) return;
        // option: skip class c (m_c = 0)
        rec(c + 1, rem);
        // option: m_c in 1..min(rem, w_c), only if c is adjacent to all chosen
        bool compat = true;
        for (int x : chosen) if (!G.adj(c, x)) { compat = false; break; }
        if (compat) {
            int hi = std::min(rem, G.w[c]);
            chosen.push_back(c);
            for (int j = 1; j <= hi; ++j) {
                m[c] = j;
                rec(c + 1, rem - j);
            }
            m[c] = 0;
            chosen.pop_back();
        }
    };
    rec(0, k);
    return ok;
}

// ---- hand-built fixed cases (printed, must pass) ----
static void handCases() {
    printf("---- fixed hand cases ----\n");
    auto run = [&](const char* name, ClassG& G, int k) {
        double bf = bruteTotal(G, k);
        auto leaves = buildClassSCT(G, k);
        double sc = sctTotal(leaves);
        std::vector<int> bad;
        bool dj = checkDisjoint(G, k, leaves, bad);
        printf("  %-28s k=%d  brute=%.0f sct=%.0f  leaves=%zu  sum=%s disj=%s\n",
               name, k, bf, sc, leaves.size(),
               (std::abs(bf - sc) < 0.5 ? "OK" : "FAIL"), (dj ? "OK" : "FAIL"));
    };
    // K4, all weight 1, k=3 => C(4,3)=4 vertex-triangles
    { ClassG G; G.C = 4; G.w.assign(4, 1); G.A.assign(4, std::vector<char>(4, 1));
      for (int i = 0; i < 4; ++i) G.A[i][i] = 0; run("K4 w=1 k=3", G, 3); }
    // single class weight 4, no edges, k=2 => C(4,2)=6 (intra-class clique)
    { ClassG G; G.C = 1; G.w.assign(1, 4); G.A.assign(1, std::vector<char>(1, 0));
      run("1 class w=4 k=2", G, 2); }
    // 3 classes weight 3 each, complete, k=4 => sum over multiplicity comps
    { ClassG G; G.C = 3; G.w.assign(3, 3); G.A.assign(3, std::vector<char>(3, 1));
      for (int i = 0; i < 3; ++i) G.A[i][i] = 0; run("K3 w=3 k=4", G, 4); }
    // edgeless 3 classes w=2, k=2 => each class alone: 3*C(2,2)=3
    { ClassG G; G.C = 3; G.w.assign(3, 2); G.A.assign(3, std::vector<char>(3, 0));
      run("edgeless 3xw2 k=2", G, 2); }
    // path 0-1-2 weights {2,2,2}, k=2 => intra(3*C(2,2)=3)+edges(2*2*2=8)=11
    { ClassG G; G.C = 3; G.w = {2,2,2}; G.A.assign(3, std::vector<char>(3,0));
      G.A[0][1]=G.A[1][0]=1; G.A[1][2]=G.A[2][1]=1; run("path3 w=2 k=2", G, 2); }
    printf("\n");
}

// ============================================================================
//  INDEPENDENT GROUND TRUTH #2 — expand classes into an explicit SIMPLE graph
//  and count k-cliques by raw vertex enumeration. This shares NO logic with
//  countCanon (no "canonical multiplicity" recursion), so it guards against a
//  bug that is common to both countCanon and the SCT.
//
//  Class c (weight w_c) -> w_c vertices, mutually adjacent (intra-class clique).
//  Cross-class edge iff adj(c,c').
// ============================================================================
static long long countKCliquesExpanded(const ClassG& G, int k) {
    // build explicit vertex graph
    std::vector<int> vclass;
    for (int c = 0; c < G.C; ++c)
        for (int t = 0; t < G.w[c]; ++t) vclass.push_back(c);
    int N = (int)vclass.size();
    if (k < 0) return 0;
    if (k == 0) return 1;  // empty clique
    if (k > N) return 0;
    std::vector<std::vector<char>> Adj(N, std::vector<char>(N, 0));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            if (vclass[i] == vclass[j]) Adj[i][j] = 1;          // same class
            else Adj[i][j] = G.A[vclass[i]][vclass[j]];          // cross class
        }
    // brute force: count k-subsets that form a clique (vertices in increasing id)
    long long cnt = 0;
    std::vector<int> stk;
    std::function<void(int)> rec = [&](int start) {
        if ((int)stk.size() == k) { ++cnt; return; }
        // prune: not enough remaining vertices
        if (N - start < k - (int)stk.size()) return;
        for (int v = start; v < N; ++v) {
            bool ok = true;
            for (int u : stk) if (!Adj[u][v]) { ok = false; break; }
            if (ok) { stk.push_back(v); rec(v + 1); stk.pop_back(); }
        }
    };
    rec(0);
    return cnt;
}

// run one case through all three: brute(canon), expanded, SCT. Return true if all agree.
static bool tripleCheck(const char* name, ClassG& G, int k, bool verbose) {
    double bf = bruteTotal(G, k);
    long long ex = countKCliquesExpanded(G, k);
    auto leaves = buildClassSCT(G, k);
    double sc = sctTotal(leaves);
    std::vector<int> bad;
    bool dj = checkDisjoint(G, k, leaves, bad);
    bool sumOk = (std::abs(bf - sc) < 0.5);
    bool exOk  = (std::abs((double)ex - sc) < 0.5);
    bool canonOk = (std::abs((double)ex - bf) < 0.5);
    bool all = sumOk && exOk && dj && canonOk;
    if (verbose || !all) {
        printf("  %-26s k=%-2d canon=%.0f expanded=%lld sct=%.0f leaves=%zu | "
               "canon==exp:%s sct==exp:%s disj:%s %s\n",
               name, k, bf, ex, sc, leaves.size(),
               canonOk?"OK":"FAIL", exOk?"OK":"FAIL", dj?"OK":"FAIL",
               all?"":"  <<<<< MISMATCH");
    }
    return all;
}

static void edgeCases() {
    printf("---- targeted EDGE cases (triple check: canon / expanded-simple / SCT) ----\n");
    int fails = 0;
    auto T = [&](const char* nm, ClassG& G, int k){ if(!tripleCheck(nm,G,k,true)) ++fails; };

    // ---- k = 1 (NEVER in random test) ----
    { ClassG G; G.C=4; G.w={1,2,3,4}; G.A.assign(4,std::vector<char>(4,1));
      for(int i=0;i<4;++i)G.A[i][i]=0; T("K4 w={1,2,3,4} k=1",G,1); } // expect 1+2+3+4=10
    { ClassG G; G.C=3; G.w={2,2,2}; G.A.assign(3,std::vector<char>(3,0));
      T("edgeless3 w=2 k=1",G,1); } // expect 6

    // ---- single class (C=1), various ----
    { ClassG G; G.C=1; G.w={5}; G.A.assign(1,std::vector<char>(1,0)); T("1cls w=5 k=1",G,1);} // 5
    { ClassG G; G.C=1; G.w={5}; G.A.assign(1,std::vector<char>(1,0)); T("1cls w=5 k=3",G,3);} // C(5,3)=10
    { ClassG G; G.C=1; G.w={5}; G.A.assign(1,std::vector<char>(1,0)); T("1cls w=5 k=5",G,5);} // C(5,5)=1
    { ClassG G; G.C=1; G.w={5}; G.A.assign(1,std::vector<char>(1,0)); T("1cls w=5 k=6 (>w)",G,6);} // 0

    // ---- k > max clique size, reachable or not ----
    { ClassG G; G.C=3; G.w={1,1,1}; G.A.assign(3,std::vector<char>(3,1));
      for(int i=0;i<3;++i)G.A[i][i]=0; T("K3 w=1 k=4 (>verts)",G,4);} // 0 (only 3 verts)
    { ClassG G; G.C=4; G.w={1,1,1,1}; G.A.assign(4,std::vector<char>(4,0));
      G.A[0][1]=G.A[1][0]=1; T("1 edge +2iso w=1 k=3",G,3);} // 0 (max clique 2)

    // ---- k = total vertices (only one clique possible: whole graph clique) ----
    { ClassG G; G.C=3; G.w={2,2,2}; G.A.assign(3,std::vector<char>(3,1));
      for(int i=0;i<3;++i)G.A[i][i]=0; T("K3 w=2 k=6 (=sum, clique)",G,6);} // 1
    { ClassG G; G.C=3; G.w={2,2,2}; G.A.assign(3,std::vector<char>(3,0));
      G.A[0][1]=G.A[1][0]=1; G.A[1][2]=G.A[2][1]=1; T("path3 w=2 k=6 (=sum,notclq)",G,6);} // 0
    { ClassG G; G.C=4; G.w={3,1,2,4}; G.A.assign(4,std::vector<char>(4,1));
      for(int i=0;i<4;++i)G.A[i][i]=0; T("K4 w={3,1,2,4} k=10(=sum)",G,10);} // 1

    // ---- disconnected: two cliques sharing nothing ----
    { ClassG G; G.C=4; G.w={2,2,2,2}; G.A.assign(4,std::vector<char>(4,0));
      G.A[0][1]=G.A[1][0]=1; G.A[2][3]=G.A[3][2]=1; T("2 disjoint edges w=2 k=2",G,2);}
    { ClassG G; G.C=6; G.w={1,1,1,1,1,1}; G.A.assign(6,std::vector<char>(6,0));
      // triangle {0,1,2} + triangle {3,4,5}
      G.A[0][1]=G.A[1][0]=G.A[1][2]=G.A[2][1]=G.A[0][2]=G.A[2][0]=1;
      G.A[3][4]=G.A[4][3]=G.A[4][5]=G.A[5][4]=G.A[3][5]=G.A[5][3]=1;
      T("2 triangles w=1 k=3",G,3);} // 2

    // ---- isolated class (weight present, no edges) mixed with a clique ----
    { ClassG G; G.C=4; G.w={2,2,2,3}; G.A.assign(4,std::vector<char>(4,0));
      G.A[0][1]=G.A[1][0]=G.A[1][2]=G.A[2][1]=G.A[0][2]=G.A[2][0]=1; // tri 0-1-2, class3 isolated
      T("tri w=2 + iso w=3 k=2",G,2);}
    { ClassG G; G.C=4; G.w={2,2,2,3}; G.A.assign(4,std::vector<char>(4,0));
      G.A[0][1]=G.A[1][0]=G.A[1][2]=G.A[2][1]=G.A[0][2]=G.A[2][0]=1;
      T("tri w=2 + iso w=3 k=3",G,3);}

    // ---- LARGE weights (w up to 10, never in random test) ----
    { ClassG G; G.C=3; G.w={10,8,6}; G.A.assign(3,std::vector<char>(3,1));
      for(int i=0;i<3;++i)G.A[i][i]=0; T("K3 w={10,8,6} k=5",G,5);}
    { ClassG G; G.C=4; G.w={7,7,7,7}; G.A.assign(4,std::vector<char>(4,1));
      for(int i=0;i<4;++i)G.A[i][i]=0; T("K4 w=7 k=6",G,6);}
    { ClassG G; G.C=2; G.w={12,12}; G.A.assign(2,std::vector<char>(2,0));
      G.A[0][1]=G.A[1][0]=1; T("edge w=12 k=8",G,8);}
    { ClassG G; G.C=1; G.w={20}; G.A.assign(1,std::vector<char>(1,0)); T("1cls w=20 k=10",G,10);} // C(20,10)

    // ---- complete graph K_C with large weights, k spanning ----
    { ClassG G; G.C=5; G.w={3,3,3,3,3}; G.A.assign(5,std::vector<char>(5,1));
      for(int i=0;i<5;++i)G.A[i][i]=0; T("K5 w=3 k=7",G,7);}
    { ClassG G; G.C=6; G.w={2,2,2,2,2,2}; G.A.assign(6,std::vector<char>(6,1));
      for(int i=0;i<6;++i)G.A[i][i]=0; T("K6 w=2 k=5",G,5);}

    // ---- mixed weights incl w=1 with w large in same dense graph ----
    { ClassG G; G.C=5; G.w={1,5,1,5,1}; G.A.assign(5,std::vector<char>(5,1));
      for(int i=0;i<5;++i)G.A[i][i]=0; T("K5 w={1,5,1,5,1} k=4",G,4);}

    // ---- star graph (center adjacent to all, leaves mutually non-adj) ----
    { ClassG G; G.C=5; G.w={3,2,2,2,2}; G.A.assign(5,std::vector<char>(5,0));
      for(int j=1;j<5;++j){G.A[0][j]=G.A[j][0]=1;} T("star5 w mixed k=2",G,2);}
    { ClassG G; G.C=5; G.w={3,2,2,2,2}; G.A.assign(5,std::vector<char>(5,0));
      for(int j=1;j<5;++j){G.A[0][j]=G.A[j][0]=1;} T("star5 w mixed k=3",G,3);}

    // ---- two overlapping cliques (share a class) ----
    { ClassG G; G.C=5; G.w={2,2,2,2,2}; G.A.assign(5,std::vector<char>(5,0));
      // clique {0,1,2} and clique {2,3,4}, shared class 2
      int tri1[3]={0,1,2}, tri2[3]={2,3,4};
      for(int a=0;a<3;++a)for(int b=a+1;b<3;++b){G.A[tri1[a]][tri1[b]]=G.A[tri1[b]][tri1[a]]=1;}
      for(int a=0;a<3;++a)for(int b=a+1;b<3;++b){G.A[tri2[a]][tri2[b]]=G.A[tri2[b]][tri2[a]]=1;}
      T("2 overlap tri w=2 k=3",G,3);}
    { ClassG G; G.C=5; G.w={2,2,2,2,2}; G.A.assign(5,std::vector<char>(5,0));
      int tri1[3]={0,1,2}, tri2[3]={2,3,4};
      for(int a=0;a<3;++a)for(int b=a+1;b<3;++b){G.A[tri1[a]][tri1[b]]=G.A[tri1[b]][tri1[a]]=1;}
      for(int a=0;a<3;++a)for(int b=a+1;b<3;++b){G.A[tri2[a]][tri2[b]]=G.A[tri2[b]][tri2[a]]=1;}
      T("2 overlap tri w=2 k=4",G,4);}

    printf("  edgeCases fails=%d\n\n", fails);
}

// ============================================================================
//  EXHAUSTIVE small-graph sweep: enumerate ALL simple graphs on C classes with
//  weights drawn from a small set, for all k. Compares SCT vs the expanded
//  simple-graph k-clique count (the truly independent oracle).
// ============================================================================
static void exhaustiveSweep() {
    printf("---- EXHAUSTIVE sweep (all graphs, expanded-oracle) ----\n");
    long long cases=0, fails=0, disjFails=0;
    int firstShown=0;
    // C up to 5, weights in {1,2,3}, all edge sets, all k in 0..(sum+1)
    for (int C = 1; C <= 5; ++C) {
        int E = C*(C-1)/2;
        std::vector<std::pair<int,int>> ep;
        for (int i=0;i<C;++i) for(int j=i+1;j<C;++j) ep.push_back({i,j});
        // weight assignments: to keep it tractable, sweep a few weight vectors
        // including all-equal and mixed. We enumerate weights in {1,2,3}^C only
        // for C<=4; for C=5 restrict to {1,2}^C.
        int wmax = (C<=4)?3:2;
        int wcombos = 1; for(int i=0;i<C;++i) wcombos*=wmax;
        for (long long em = 0; em < (1LL<<E); ++em) {
            ClassG G; G.C=C; G.A.assign(C,std::vector<char>(C,0));
            for (int b=0;b<E;++b) if (em&(1LL<<b)) {
                G.A[ep[b].first][ep[b].second]=G.A[ep[b].second][ep[b].first]=1;
            }
            for (int wc = 0; wc < wcombos; ++wc) {
                G.w.assign(C,0);
                int x=wc; for(int i=0;i<C;++i){ G.w[i]=1+(x%wmax); x/=wmax; }
                int sum=0; for(int i=0;i<C;++i) sum+=G.w[i];
                for (int k=0;k<=sum+1;++k) {
                    long long ex = countKCliquesExpanded(G,k);
                    auto leaves = buildClassSCT(G,k);
                    double sc = sctTotal(leaves);
                    ++cases;
                    if (std::abs((double)ex - sc) > 0.5) {
                        ++fails;
                        if (firstShown<15){ ++firstShown;
                            printf("  EXH-FAIL C=%d k=%d expanded=%lld sct=%.1f w=[",C,k,ex,sc);
                            for(int i=0;i<C;++i)printf("%d ",G.w[i]);
                            printf("] em=%lld\n",em);
                        }
                    }
                    // disjointness only meaningful for k>=1
                    if (k>=1) {
                        std::vector<int> bad;
                        if (!checkDisjoint(G,k,leaves,bad)) ++disjFails;
                    }
                }
            }
        }
        printf("  C=%d done. cumulative cases=%lld fails=%lld disjFails=%lld\n",
               C, cases, fails, disjFails);
    }
    printf("  EXHAUSTIVE: cases=%lld  sum-fails=%lld  disj-fails=%lld\n\n",
           cases, fails, disjFails);
}

// ============================================================================
//  EFFICIENCY: leaf count and build time on dense K_C and overlapping cliques.
// ============================================================================
#include <chrono>
static long long countMaximalCliquesApprox(const ClassG& G); // fwd
static void efficiency() {
    printf("---- EFFICIENCY: leaf count vs C on dense/adversarial graphs ----\n");
    initPascal();
    // K_C, all weight w, several k. Measure leaves and time.
    for (int C = 2; C <= 18; ++C) {
        ClassG G; G.C=C; G.w.assign(C,3); G.A.assign(C,std::vector<char>(C,1));
        for(int i=0;i<C;++i)G.A[i][i]=0;
        int k = std::min(C, 6);
        auto t0=std::chrono::high_resolution_clock::now();
        auto leaves=buildClassSCT(G,k);
        auto t1=std::chrono::high_resolution_clock::now();
        double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
        printf("  K%-2d w=3 k=%d : leaves=%-8zu  time=%.2f ms\n", C, k, leaves.size(), ms);
    }
    printf("\n  -- K_C with k = C (forces all classes, single big clique) --\n");
    for (int C = 2; C <= 16; ++C) {
        ClassG G; G.C=C; G.w.assign(C,4); G.A.assign(C,std::vector<char>(C,1));
        for(int i=0;i<C;++i)G.A[i][i]=0;
        auto leaves=buildClassSCT(G,C);
        printf("  K%-2d w=4 k=%d : leaves=%zu\n", C, C, leaves.size());
    }
    printf("\n  -- adversarial: chain of overlapping triangles (sparse, many maximal cliques) --\n");
    // classes 0..C-1, triangles (i,i+1,i+2) overlapping -> path power graph
    for (int C = 4; C <= 16; C += 2) {
        ClassG G; G.C=C; G.w.assign(C,2); G.A.assign(C,std::vector<char>(C,0));
        for(int i=0;i<C;++i) for(int d=1;d<=2;++d) if(i+d<C){G.A[i][i+d]=G.A[i+d][i]=1;}
        int k=3;
        auto t0=std::chrono::high_resolution_clock::now();
        auto leaves=buildClassSCT(G,k);
        auto t1=std::chrono::high_resolution_clock::now();
        double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
        printf("  pathpow2 C=%-2d w=2 k=3 : leaves=%-6zu time=%.2f ms\n", C, leaves.size(), ms);
    }
    printf("\n  -- worst-case probe: Turán-like / complete multipartite (cocktail party) --\n");
    // cocktail party graph: C classes, complete EXCEPT a perfect matching of non-edges.
    for (int C = 4; C <= 16; C += 2) {
        ClassG G; G.C=C; G.w.assign(C,2); G.A.assign(C,std::vector<char>(C,1));
        for(int i=0;i<C;++i)G.A[i][i]=0;
        for(int i=0;i+1<C;i+=2){ G.A[i][i+1]=G.A[i+1][i]=0; } // remove matching
        int k=4;
        auto t0=std::chrono::high_resolution_clock::now();
        auto leaves=buildClassSCT(G,k);
        auto t1=std::chrono::high_resolution_clock::now();
        double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
        // also count via expanded oracle for sanity (only if not too big)
        printf("  cocktail C=%-2d w=2 k=4 : leaves=%-6zu time=%.2f ms\n", C, leaves.size(), ms);
    }
    printf("\n");
}
static long long countMaximalCliquesApprox(const ClassG&) { return 0; }

int main() {
    initPascal();
    handCases();
    edgeCases();
    exhaustiveSweep();
    efficiency();
    std::mt19937 rng(20260618u);

    const int trials = 20000;
    int passed = 0, failed = 0, firstFailShown = 0;
    int disjPass = 0, disjFail = 0, disjFailShown = 0;
    long long totalLeaves = 0;

    for (int t = 0; t < trials; ++t) {
        int C = 3 + (int)(rng() % 7);      // 3..9
        ClassG G; G.C = C; G.w.resize(C);
        for (int c = 0; c < C; ++c) G.w[c] = 1 + (int)(rng() % 4);  // 1..4
        G.A.assign(C, std::vector<char>(C, 0));
        // density spans the FULL spectrum [0,1] incl. empty & complete graphs.
        double pEdge = (double)(rng() % 101) / 100.0;
        for (int i = 0; i < C; ++i)
            for (int j = i + 1; j < C; ++j) {
                char e = ((double)(rng() % 1000) / 1000.0 < pEdge) ? 1 : 0;
                G.A[i][j] = G.A[j][i] = e;
            }
        int k = 2 + (int)(rng() % 4);      // 2..5

        double bf = bruteTotal(G, k);
        std::vector<CCPath> leaves = buildClassSCT(G, k);
        totalLeaves += (long long)leaves.size();
        double sc = sctTotal(leaves);

        if (std::abs(bf - sc) < 0.5) {
            ++passed;
        } else {
            ++failed;
            if (firstFailShown < 12) {
                ++firstFailShown;
                printf("SUM-FAIL t=%d C=%d k=%d  brute=%.1f sct=%.1f  w=[", t, C, k, bf, sc);
                for (int c = 0; c < C; ++c) printf("%d ", G.w[c]);
                printf("]  adj=\n");
                for (int i = 0; i < C; ++i) {
                    printf("    ");
                    for (int j = 0; j < C; ++j) printf("%d", (int)G.A[i][j]);
                    printf("\n");
                }
            }
        }

        // disjointness: every weighted-clique pattern represented exactly once.
        std::vector<int> bad;
        if (checkDisjoint(G, k, leaves, bad)) {
            ++disjPass;
        } else {
            ++disjFail;
            if (disjFailShown < 12) {
                ++disjFailShown;
                printf("DISJ-FAIL t=%d C=%d k=%d  pattern m=[", t, C, k);
                for (int c = 0; c < C; ++c) printf("%d ", bad[c]);
                printf("]  w=[");
                for (int c = 0; c < C; ++c) printf("%d ", G.w[c]);
                printf("]\n");
            }
        }
    }

    printf("\n==== RANDOM ADVERSARIAL TEST (class SCT vs brute force) ====\n");
    printf("[sum-equality]   trials=%d  PASS=%d  FAIL=%d\n", trials, passed, failed);
    printf("[disjointness]   trials=%d  PASS=%d  FAIL=%d  (each pattern represented exactly once)\n",
           trials, disjPass, disjFail);
    bool allOk = (failed == 0 && disjFail == 0);
    printf("[stats] avg_leaves=%.2f  overall=%s\n",
           (double)totalLeaves / trials,
           allOk ? "[ALL EXACT + DISJOINT]" : "[MISMATCH]");
    return allOk ? 0 : 1;
}
