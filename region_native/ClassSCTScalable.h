// ClassSCTScalable.h - SPARSE, degeneracy-seeded class-weighted Succinct
// Clique Tree builder. Same output contract as the dense buildClassSCT in
// ClassSCT.h:
//
//   scalableBuildClassSCT(C, w, adj, k) -> vector<ccpath::CCPath> such that
//     * leaves are DISJOINT (each weighted k-clique pattern represented once),
//     * sum over leaves of support_count(leaf, b=0, nCr) == weighted-k-clique
//       count,
//     * leaves are full-dimension-C CCPaths (h=0; spine ell=u=mult on the
//       global class id; free pool ell=0,u=w on the global class id),
// but built WITHOUT a C x C adjacency matrix and WITHOUT an O(|P|^2) pivot over
// all classes. Suitable for quotient graphs with C up to ~125000.
//
// =====================================================================
//  ALGORITHM (standard scalable Pivoter, lifted to weighted classes)
//  --------------------------------------------------------------------
//  The dense buildClassSCT runs ONE orbit-aware pivot recursion (gen()) over
//  the open set P = ALL classes. Its disjointness comes from a canonical
//  "lowest-non-neighbour" branching inside gen(); its correctness oracle is
//  the brute force countCanon (see class_sct.cpp), which enumerates cliques
//  by, for each class c in id order, taking j=1..w_c copies of c (factor
//  C(w_c,j)) and recursing on c's LATER-neighbours (neighbours after c in id
//  order). Each clique is thus generated exactly once -- by its lowest-id
//  class.
//
//  The scalable builder keeps gen() exactly (same spine / free-pool /
//  canonical-lowest-non-neighbour logic, same orbit-weight handling), but:
//
//   1. It computes a DEGENERACY ordering of the quotient graph
//      (Eppstein-Loffler-Strash: repeatedly remove a min-degree class).
//      pos[c] = position of class c in that order.
//
//   2. It enumerates cliques by SEED, where the seed of a clique is its
//      DEGENERACY-LOWEST class (replacing countCanon's id-lowest class). For
//      each seed c, in degeneracy order:
//         later(c) = { t in N(c) : pos[t] > pos[c] }            (sparse)
//      The seed contributes a SPINE entry (c, j) with j = 1..w_c (factor
//      C(w_c, j), the same role the seed plays in countCanon), and the rest
//      of the clique is built from later(c) by the SAME gen() recursion --
//      now restricted to the bounded set later(c) (|later(c)| <= degeneracy d,
//      NOT C). The seed is fully decided here (it never re-enters the open
//      set), exactly as countCanon recurses on `later` which excludes c.
//
//   3. DISJOINTNESS. Every weighted clique has a UNIQUE degeneracy-lowest
//      class, so it is produced by exactly one seed. Within a seed, gen() is
//      the validated disjoint enumerator over later(c). The seed's spine
//      multiplicity j ranges over 1..w_c and the deeper recursion draws the
//      rest of the clique from later(c) only -- so no class supplies vertices
//      to two different leaves of the same pattern. This is the same
//      guarantee Eppstein/Pivoter gives, now carrying class weights.
//
//   4. ORBIT WEIGHTS. The seed's own weight is the single spine factor
//      C(w_c, j). A NEIGHBOUR class v in later(c) is handled entirely inside
//      gen() (its multiplicity is split across spine/pool exactly as in the
//      dense version), so the per-class binomial is produced correctly. The
//      seed never appears in later(c) (later(c) is a subset of N(c), and a
//      class is not its own neighbour), so the seed's weight is never split a
//      second time -- it is one clean C(w_c, j) per leaf.
//
//   5. ADJACENCY is sparse throughout: gen() tests adj(i,j) via a hash-set
//      membership query (a per-call SparseAdj wrapper), and pivot scoring
//      intersects sorted neighbour lists. Each gen() subtree sees only a
//      bounded open set (<= d), so its internal O(|P|^2) pivot is cheap.
//
//  No C x C matrix is ever allocated; peak memory is O(C + sum of degrees +
//  output leaves).
// =====================================================================

#ifndef CLASS_SCT_SCALABLE_H
#define CLASS_SCT_SCALABLE_H

#include "../src/NucleusDecomposition/CCPathCore.h"
#include <algorithm>
#include <unordered_set>
#include <utility>
#include <vector>

namespace classsct_scalable {

using ccpath::CCPath;
using ccpath::Vec;

// An "open class": its global class id and the residual weight available here.
// Identical role to ClassSCT.h's OpenC.
struct OpenC { int c; int wres; };

// Sparse adjacency view used inside the recursion. adj[c] is c's sorted
// neighbour list (global ids). adjacency between two classes is decided by a
// hash-set membership check on the shorter neighbour list -- this keeps the
// recursion's adj(i,j) O(1) average without a C x C matrix.
struct SparseAdj {
    const std::vector<int>* w;          // weights (global)
    const std::vector<std::vector<int>>* adj;  // neighbour lists (global, sorted)
    // per-recursion neighbour membership sets (built lazily per pivot call)
    bool adjacent(int i, int j) const {
        // binary search in the shorter neighbour list
        const std::vector<int>& Ni = (*adj)[i];
        const std::vector<int>& Nj = (*adj)[j];
        if (Ni.size() <= Nj.size())
            return std::binary_search(Ni.begin(), Ni.end(), j);
        return std::binary_search(Nj.begin(), Nj.end(), i);
    }
    int weight(int c) const { return (*w)[c]; }
};

// Build a COMPACT leaf from a finalised spine (forced classes) + a free pool
// (a clique of universal pivots). The leaf's dimension is the number of
// TOUCHED classes (spine + pool), NOT C: emitting full-dimension-C leaves is
// what makes the dense builder blow up at C ~ 1e5 (4 Vecs of length C per
// leaf, millions of leaves -> hundreds of GB). The global class id of local
// slot i is recorded in p.classIds[i] (the CCPath field reserved for exactly
// this). p.classIds is SORTED by global id.
//
// Per-class encoding matches ClassSCT.h's emitLeaf EXACTLY (h=0; spine class:
// n=w_c, ell=u=mult; pool class: n=w_c, ell=0, u=w_c; T=k). support_count is
// dimension-agnostic (it only reads n/ell/u/T/forbidden over the leaf's own
// slots, with b=0 on every slot), so a compact leaf yields the IDENTICAL count
// as the full-C leaf: the C - |touched| absent classes would all carry
// n=ell=u=0 and contribute a C(0,0)=1 factor that changes nothing. Each class
// is in at most one of {spine, pool} for a given leaf, so no aliasing.
static void emitLeaf(int /*C*/, const std::vector<int>& w,
                     const std::vector<std::pair<int,int>>& spine,
                     const std::vector<OpenC>& pool, int k,
                     std::vector<CCPath>& out) {
    int forced = 0;
    for (auto& sp : spine) forced += sp.second;
    if (forced > k) return;           // spine alone overshoots
    int poolCap = 0;
    for (auto& pc : pool) poolCap += (int)w[pc.c];
    if (forced + poolCap < k) return; // cannot reach k even filling the pool

    // collect touched (globalId, role, mult) then sort by global id so the
    // compact slots are in canonical order.
    //   role 0 = spine (ell=u=mult);  role 1 = pool (ell=0, u=w)
    const int M = (int)spine.size() + (int)pool.size();
    std::vector<std::pair<int,int>> slots;   // (globalId, slotIndexIntoTemp)
    slots.reserve(M);
    std::vector<int> gid; gid.reserve(M);
    std::vector<int> nn;  nn.reserve(M);
    std::vector<int> lo;  lo.reserve(M);
    std::vector<int> up;  up.reserve(M);
    for (auto& sp : spine) {
        slots.push_back({sp.first, (int)gid.size()});
        gid.push_back(sp.first); nn.push_back(w[sp.first]);
        lo.push_back(sp.second); up.push_back(sp.second);
    }
    for (auto& pc : pool) {
        slots.push_back({pc.c, (int)gid.size()});
        gid.push_back(pc.c); nn.push_back(w[pc.c]);
        lo.push_back(0); up.push_back(w[pc.c]);
    }
    std::sort(slots.begin(), slots.end());   // by global id

    CCPath p;
    p.h = ccpath::zeros_vec(M);              // h = 0 everywhere
    p.n.resize(M); p.ell.resize(M); p.u.resize(M);
    p.classIds.resize(M);
    for (int i = 0; i < M; ++i) {
        int src = slots[i].second;
        p.n[i]   = (int16_t)nn[src];
        p.ell[i] = (int16_t)lo[src];
        p.u[i]   = (int16_t)up[src];
        p.classIds[i] = (int32_t)gid[src];
    }
    p.T = k;                                 // h = 0 => T = k
    out.push_back(std::move(p));
}

// Orbit-aware pivot recursion -- a faithful copy of ClassSCT.h's gen(), with
// the dense G.adj / G.w replaced by the SparseAdj view and C passed
// explicitly. Semantics (spine / free pool / canonical-lowest-non-neighbour
// branching / orbit-weight split) are IDENTICAL to the dense version; only the
// adjacency backing store differs. P is always a BOUNDED open set (the seed's
// later-neighbours, shrinking with depth), so the per-call O(|P|^2) pivot is
// cheap.
static void gen(const SparseAdj& G, int C,
                std::vector<std::pair<int,int>>& spine, int spineSum,
                std::vector<OpenC>& pool, int poolSum,
                std::vector<OpenC> P, int k,
                std::vector<CCPath>& out) {
    if (spineSum > k) return;
    if (P.empty()) { emitLeaf(C, *G.w, spine, pool, k, out); return; }

    // ---- pick pivot pc in P maximising residual weight of P adjacent to pc.
    int bestIdx = 0; long bestScore = -1;
    for (size_t i = 0; i < P.size(); ++i) {
        long sc = 0;
        for (size_t j = 0; j < P.size(); ++j)
            if (j != i && G.adjacent(P[i].c, P[j].c)) sc += P[j].wres;
        if (sc > bestScore) { bestScore = sc; bestIdx = (int)i; }
    }
    int pc = P[bestIdx].c;

    // Is pc universal in P (adjacent to every other open class)?
    bool universal = true;
    for (auto& q : P) if (q.c != pc && !G.adjacent(pc, q.c)) { universal = false; break; }

    if (universal) {
        // pull pc into the free POOL, recurse on P \ {pc}. Pool stays a clique
        // because every later pool member comes from P \ {pc} subset of N(pc).
        std::vector<OpenC> Pr;
        Pr.reserve(P.size() - 1);
        for (auto& q : P) if (q.c != pc) Pr.push_back(q);
        pool.push_back(OpenC{pc, /*wres=*/G.weight(pc)});
        gen(G, C, spine, spineSum, pool, poolSum + G.weight(pc), std::move(Pr), k, out);
        pool.pop_back();
        return;
    }

    // pc has >=1 non-neighbour in P. Non-neighbours (EXCLUDING pc), by id.
    std::vector<int> nonNb; nonNb.reserve(P.size());   // <= |P| non-neighbours
    for (auto& q : P) if (q.c != pc && !G.adjacent(pc, q.c)) nonNb.push_back(q.c);
    std::sort(nonNb.begin(), nonNb.end());

    // ---------- branch A: use NO non-neighbour of pc ----------
    {
        std::vector<OpenC> Pp;
        Pp.reserve(P.size());
        for (auto& q : P)
            if (q.c == pc || G.adjacent(pc, q.c)) Pp.push_back(q);
        gen(G, C, spine, spineSum, pool, poolSum, std::move(Pp), k, out);
    }

    // ---------- branch B: the LOWEST non-neighbour-of-pc used is v ----------
    for (size_t t = 0; t < nonNb.size(); ++t) {
        int v = nonNb[t];
        int wv = 0;
        for (auto& q : P) if (q.c == v) { wv = q.wres; break; }

        // canonical cut: drop pc-non-neighbours with id <= v (and v itself).
        // nonNb is small (bounded by |P| <= d), so a sorted-vector membership
        // test is cheaper than a length-C `drop` array here.
        std::vector<OpenC> base;
        base.reserve(P.size());
        for (auto& q : P) {
            // q must be adjacent to v, must not be a dropped lower non-neighbour
            if (q.c == v) continue;                       // v is decided
            bool dropped = false;
            for (size_t z = 0; z <= t; ++z) if (nonNb[z] == q.c) { dropped = true; break; }
            if (dropped) continue;
            if (!G.adjacent(v, q.c)) continue;
            base.push_back(q);
        }

        for (int mtt = 1; mtt <= wv; ++mtt) {
            if (spineSum + mtt > k) break;
            spine.push_back({v, mtt});
            std::vector<OpenC> child = base;     // already excludes v
            gen(G, C, spine, spineSum + mtt, pool, poolSum, std::move(child), k, out);
            spine.pop_back();
        }
    }
}

// Compute a degeneracy ordering of the quotient graph via the linear-time
// bucket method (Eppstein-Loffler-Strash / Matula-Beck): repeatedly remove a
// minimum-degree class. Returns pos[c] = position of class c in the order
// (0 = first removed = globally lowest in the canonical seed order). Only
// classes with w[c] > 0 participate; empty classes get pos = -1 and are never
// seeds (and never neighbours of a seed, since adj lists are over real
// classes).
static std::vector<int> degeneracyOrder(int C, const std::vector<int>& w,
                                        const std::vector<std::vector<int>>& adj) {
    std::vector<int> deg(C, 0);
    int maxd = 0;
    for (int c = 0; c < C; ++c) {
        if (w[c] <= 0) { deg[c] = -1; continue; }
        deg[c] = (int)adj[c].size();
        if (deg[c] > maxd) maxd = deg[c];
    }
    // bucket[d] = list of live classes with current degree d.
    std::vector<std::vector<int>> bucket(maxd + 1);
    for (int c = 0; c < C; ++c) if (deg[c] >= 0) bucket[deg[c]].push_back(c);
    std::vector<char> removed(C, 0);
    std::vector<int> pos(C, -1);
    int placed = 0, curd = 0;
    int nLive = 0; for (int c = 0; c < C; ++c) if (w[c] > 0) ++nLive;

    while (placed < nLive) {
        // find smallest non-empty bucket (curd only moves up between removals;
        // a removal can lower a neighbour's degree below curd, so scan down too)
        while (curd <= maxd && bucket[curd].empty()) ++curd;
        int d2 = 0; while (d2 <= maxd && bucket[d2].empty()) ++d2;
        if (d2 < curd) curd = d2;
        int v = bucket[curd].back(); bucket[curd].pop_back();
        if (removed[v]) continue;                  // stale bucket entry
        removed[v] = 1; pos[v] = placed++;
        for (int wv : adj[v]) {
            if (!removed[wv] && deg[wv] >= 0) {
                int nd = --deg[wv];
                if (nd < 0) nd = deg[wv] = 0;
                bucket[nd].push_back(wv);
                if (nd < curd) curd = nd;
            }
        }
    }
    return pos;
}

// =====================================================================
//  scalableBuildClassSCT
//  --------------------------------------------------------------------
//  C        : number of classes
//  w        : w[c] = weight of class c (>=0; 0-weight classes are ignored)
//  adj      : adj[c] = SORTED list of c's neighbour class ids (sparse; no
//             self-loops, symmetric: j in adj[i] <=> i in adj[j])
//  k        : clique size (the SCT counts weighted k-cliques)
//  returns  : DISJOINT full-dimension-C CCPath leaves (same contract as the
//             dense buildClassSCT in ClassSCT.h).
// =====================================================================
inline std::vector<CCPath>
scalableBuildClassSCT(int C, const std::vector<int>& w,
                      const std::vector<std::vector<int>>& adj, int k) {
    std::vector<CCPath> out;
    if (C <= 0 || k <= 0) return out;

    std::vector<int> pos = degeneracyOrder(C, w, adj);

    SparseAdj G;
    G.w = &w;
    G.adj = &adj;

    // Seed loop: each class c is the canonical (degeneracy-lowest) class of the
    // cliques it generates. The rest of each clique comes from c's
    // later-neighbours. The seed contributes spine (c, j) for j=1..w_c.
    std::vector<std::pair<int,int>> spine;
    std::vector<OpenC> pool;
    for (int c = 0; c < C; ++c) {
        if (w[c] <= 0) continue;
        // later(c) = neighbours after c in the degeneracy order (with weight).
        std::vector<OpenC> later;
        later.reserve(adj[c].size());
        for (int t : adj[c]) {
            if (w[t] <= 0) continue;
            if (pos[t] > pos[c]) later.push_back(OpenC{t, w[t]});
        }
        // gen() needs P sorted by global class id (its canonical-lowest
        // branching and the dense oracle both order non-neighbours by id).
        std::sort(later.begin(), later.end(),
                  [](const OpenC& a, const OpenC& b){ return a.c < b.c; });

        int hi = std::min(k, w[c]);
        for (int j = 1; j <= hi; ++j) {
            spine.clear();
            spine.push_back({c, j});
            pool.clear();
            std::vector<OpenC> P = later;     // fresh open set per seed-mult
            gen(G, C, spine, /*spineSum=*/j, pool, /*poolSum=*/0,
                std::move(P), k, out);
        }
    }
    return out;
}

}  // namespace classsct_scalable

#endif  // CLASS_SCT_SCALABLE_H
