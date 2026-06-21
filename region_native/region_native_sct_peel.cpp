// region_native_sct_peel.cpp — region-native (r,s)-nucleus peel on a
// VALIDATED class-level Succinct Clique Tree (ClassSCT). Front half
// (load/MCE/classes/pattern-enum) shared with region_native_peel.cpp; the
// region-IE support (suppOf/unionAlive/leafCount) is KEPT only as an init
// cross-check (gate G2a). The PEEL itself runs on the GLOBAL class-SCT:
//
//   * Build the quotient graph (nodes=classes, edge iff two classes co-occur
//     in some region), build its class-weighted s-clique SCT once. The SCT
//     leaves are DISJOINT CCPaths: each weighted s-clique witness lives in
//     exactly one leaf.
//   * support(pattern m) = SUM over the leaves hosting m of
//     support_count(leaf, b=m). A clean SUM, NO inclusion-exclusion across
//     regions (the whole point: leaves are disjoint so there is nothing to
//     subtract).
//   * Peel = bucket queue at PATTERN granularity. Peeling P inserts P's
//     threshold (tuple_to_threshold) into every hosting leaf's forbidden
//     antichain; a leaf whose antichain exceeds KMAX is replaced by its
//     controlled_split (a SET of CCPaths summing to the same count). The
//     leaf<->pattern maps are keyed by the ORIGINAL leaf id, which is stable
//     under splitting (children partition the parent's witnesses, so summing
//     support_count over a slot's split set is exact and the host relation
//     never changes). Affected patterns = those sharing a leaf with P.
//   * Output = the core-number distribution, weighted by pattern.mult.
#include <unordered_set>
#include <map>
#include "CCPathCore.h"
#include "ClassSCT.h"          // dense orbit-aware class-SCT (oracle)
#include "ClassSCTScalable.h"  // scalable sparse class-SCT (production path)
// counting for (r,s)-nucleus INITIAL SUPPORT. Touches no existing code.
//
// Idea (SigmodPlus section 24): compute the s-clique support of every
// region tuple directly on the Region/Class quotient, with NO CPI/SDCT
// pivot index. Pivot theory is preserved as a region union-count
// inclusion-exclusion:
//   sup(tau) = | { S : R0 subset S, |S|=s, S clique } |
//            = unionCount over host regions of  C(|cap region|-r, s-r)
// computed by the same B&B (dominance + size<s pruning) the engine uses
// for dead boxes, but on REGIONS instead of pivot boxes. Because each
// profile class is wholly in or out of a region, intersections are done
// at CLASS granularity (few classes per region), not vertex granularity.
//
// Build: g++ -O3 -std=c++17 -o region_native region_native.cpp
// Run:   ./region_native graph.edges r s [--verify N] [--mce-budget SEC]
//
// Adversarial self-check (--verify): compares sup(tau) for a sample (or
// all) tuples against direct s-clique enumeration (ground truth).
#include <algorithm>
#include <functional>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using Clock = chrono::high_resolution_clock;
static double secs(Clock::time_point a, Clock::time_point b) {
    return chrono::duration_cast<chrono::duration<double>>(b - a).count();
}

// ----- graph (CSR, undirected, 0-indexed; header "n m") -----
struct Graph {
    int n = 0;
    vector<int> off, adj;  // CSR
    bool adjacent(int u, int v) const {
        // binary search in u's neighbor list
        const int *b = &adj[off[u]], *e = &adj[off[u + 1]];
        return std::binary_search(b, e, v);
    }
    int deg(int u) const { return off[u + 1] - off[u]; }
};

static Graph load_graph(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    long n, m;
    if (fscanf(f, "%ld %ld", &n, &m) != 2) { fprintf(stderr, "bad header\n"); exit(1); }
    vector<pair<int,int>> es;
    es.reserve(m * 2);
    int u, v, maxid = (int)n - 1;
    while (fscanf(f, "%d %d", &u, &v) == 2) {
        if (u == v) continue;
        es.push_back({u, v});
        es.push_back({v, u});
        maxid = max(maxid, max(u, v));
    }
    fclose(f);
    Graph g;
    g.n = maxid + 1;
    vector<int> cnt(g.n + 1, 0);
    for (auto &e : es) cnt[e.first]++;
    g.off.assign(g.n + 1, 0);
    for (int i = 0; i < g.n; i++) g.off[i + 1] = g.off[i] + cnt[i];
    g.adj.resize(g.off[g.n]);
    vector<int> cur(g.off.begin(), g.off.end() - 1);
    for (auto &e : es) g.adj[cur[e.first]++] = e.second;
    // sort + dedup each neighbor list
    for (int i = 0; i < g.n; i++) {
        int *b = &g.adj[g.off[i]], *e = &g.adj[g.off[i + 1]];
        sort(b, e);
        int *ne = unique(b, e);
        // compact (rebuild offsets after dedup)
        (void)ne;
    }
    // rebuild with dedup compaction
    vector<int> noff(g.n + 1, 0), nadj;
    nadj.reserve(g.adj.size());
    for (int i = 0; i < g.n; i++) {
        int *b = &g.adj[g.off[i]], *e = &g.adj[g.off[i + 1]];
        int *ne = unique(b, e);
        noff[i + 1] = noff[i] + (int)(ne - b);
        for (int *p = b; p < ne; p++) nadj.push_back(*p);
    }
    g.off = move(noff);
    g.adj = move(nadj);
    return g;
}

// ----- Bron-Kerbosch with pivoting; collect maximal cliques size>=minSz -----
struct MCE {
    const Graph &g;
    int minSz;
    double budget;
    Clock::time_point t0;
    bool aborted = false;
    vector<vector<int>> cliques;
    MCE(const Graph &g_, int m, double b) : g(g_), minSz(m), budget(b) {}

    // sets as sorted vectors; intersect with a vertex's neighborhood
    static vector<int> intersect_nbr(const vector<int> &S, const Graph &g, int v) {
        vector<int> out;
        const int *b = &g.adj[g.off[v]], *e = &g.adj[g.off[v + 1]];
        out.reserve(min((size_t)(e - b), S.size()));
        size_t i = 0; const int *p = b;
        while (i < S.size() && p < e) {
            if (S[i] < *p) i++;
            else if (S[i] > *p) p++;
            else { out.push_back(S[i]); i++; p++; }
        }
        return out;
    }

    void rec(vector<int> &R, vector<int> P, vector<int> X) {
        if (aborted) return;
        if (P.empty() && X.empty()) {
            if ((int)R.size() >= minSz) cliques.push_back(R);
            return;
        }
        if (secs(t0, Clock::now()) > budget) { aborted = true; return; }
        // choose pivot u in P∪X maximizing |P ∩ N(u)|
        int bestu = -1, bestc = -1;
        auto consider = [&](const vector<int> &S) {
            for (int u : S) {
                int c = (int)intersect_nbr(P, g, u).size();
                if (c > bestc) { bestc = c; bestu = u; }
            }
        };
        consider(P); consider(X);
        // candidates = P \ N(bestu)
        vector<int> cand; cand.reserve(P.size());   // |cand| <= |P|, avoid regrowth
        {
            const int *b = &g.adj[g.off[bestu]], *e = &g.adj[g.off[bestu + 1]];
            size_t i = 0; const int *p = b;
            // P \ N(bestu)
            for (int x : P) {
                while (p < e && *p < x) p++;
                if (!(p < e && *p == x)) cand.push_back(x);
            }
        }
        for (int v : cand) {
            R.push_back(v);
            rec(R, intersect_nbr(P, g, v), intersect_nbr(X, g, v));
            R.pop_back();
            // move v from P to X
            P.erase(lower_bound(P.begin(), P.end(), v));
            X.insert(lower_bound(X.begin(), X.end(), v), v);
            if (aborted) return;
        }
    }

    bool run() {
        t0 = Clock::now();
        // Eppstein-Loffler-Strash: process vertices in degeneracy order;
        // for v, run BK-pivot on P = later-neighbours, X = earlier-
        // neighbours. This bounds |P| by the degeneracy d at every seed,
        // killing the naive root cost (P = all n vertices). Each maximal
        // clique is emitted exactly once.
        int n = g.n;
        vector<int> deg(n), order, pos(n, -1);
        for (int v = 0; v < n; v++) deg[v] = g.deg(v);
        int maxd = 0; for (int v = 0; v < n; v++) maxd = max(maxd, deg[v]);
        vector<vector<int>> bucket(maxd + 1);
        for (int v = 0; v < n; v++) bucket[deg[v]].push_back(v);
        vector<char> removed(n, 0);
        order.reserve(n);
        int curd = 0;
        for (int processed = 0; processed < n; processed++) {
            while (curd <= maxd && bucket[curd].empty()) curd++;
            // (curd may need to back down when a removal lowers a degree)
            int d2 = 0; // find smallest non-empty bucket from scratch-safe
            while (d2 <= maxd && bucket[d2].empty()) d2++;
            int v = bucket[d2].back(); bucket[d2].pop_back();
            if (removed[v]) { processed--; continue; }
            removed[v] = 1; pos[v] = (int)order.size(); order.push_back(v);
            for (int w : g_nbr(v)) {
                if (!removed[w]) {
                    deg[w]--; if (deg[w] < 0) deg[w] = 0;
                    bucket[deg[w]].push_back(w);
                }
            }
        }
        vector<int> R;
        for (int idx = 0; idx < n; idx++) {
            int v = order[idx];
            vector<int> P, X;
            P.reserve(g.deg(v)); X.reserve(g.deg(v));   // |P|+|X| == deg(v); fresh per vertex
            for (int w : g_nbr(v)) {
                if (pos[w] > idx) P.push_back(w);
                else X.push_back(w);
            }
            sort(P.begin(), P.end()); sort(X.begin(), X.end());
            R.assign(1, v);
            rec(R, move(P), move(X));
            if (aborted) return false;
        }
        return !aborted;
    }
    // neighbour range helper
    struct NbrRange { const int *b, *e; const int* begin() const{return b;} const int* end() const{return e;} };
    NbrRange g_nbr(int v) const { return {&g.adj[g.off[v]], &g.adj[g.off[v + 1]]}; }
};

// ----- nCr as double (matches engine's double-based core values) -----
static vector<vector<double>> NCR;
static void build_ncr(int N) {
    NCR.assign(N + 1, vector<double>(N + 1, 0.0));
    for (int i = 0; i <= N; i++) {
        NCR[i][0] = 1.0;
        for (int j = 1; j <= i; j++)
            NCR[i][j] = NCR[i - 1][j - 1] + NCR[i - 1][j];
    }
}
static inline double C(int n, int k) {
    if (k < 0 || n < 0 || k > n) return 0.0;
    return NCR[n][k];
}
// nCr wrapper for ccpath::support_count (forbidden-antichain IE leaf).
static double ccpath_ncr(int n, int k) { return C(n, k); }

// current resident-set, GB (Linux /proc) -- per-phase memory breakdown (MEM_DBG).
static double rssGB() {
    FILE *f = fopen("/proc/self/status", "r"); if (!f) return 0;
    char ln[256]; long kb = 0;
    while (fgets(ln, sizeof ln, f)) if (sscanf(ln, "VmRSS: %ld kB", &kb) == 1) break;
    fclose(f); return kb / 1048576.0;
}

int main(int argc, char **argv) {
    bool memDbg = getenv("MEM_DBG") != nullptr;
    auto memCk = [&](const char *tag) { if (memDbg) fprintf(stderr, "[mem] %-22s RSS=%.2fG\n", tag, rssGB()); };
    if (argc < 4) { fprintf(stderr, "usage: %s graph.edges r s [--verify N] [--mce-budget S]\n", argv[0]); return 1; }
    const char *gpath = argv[1];
    int r = atoi(argv[2]), s = atoi(argv[3]);
    int verifyN = 0; double mceBudget = 120.0;
    for (int i = 4; i < argc; i++) {
        if (!strcmp(argv[i], "--verify") && i + 1 < argc) verifyN = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--mce-budget") && i + 1 < argc) mceBudget = atof(argv[++i]);
    }

    auto T0 = Clock::now();
    Graph g = load_graph(gpath);
    auto T1 = Clock::now();
    printf("[rn] graph n=%d m=%ld  (r,s)=(%d,%d)  load=%.2fs\n",
           g.n, (long)g.adj.size() / 2, r, s, secs(T0, T1));

    // regions = maximal cliques of size >= s
    MCE mce(g, s, mceBudget);
    bool ok = mce.run();
    auto T2 = Clock::now();
    if (!ok) { printf("[rn] MCE exceeded %.0fs budget; abort.\n", mceBudget); return 0; }
    auto &regions = mce.cliques;
    for (auto &R : regions) sort(R.begin(), R.end());
    printf("[rn] regions(>=s)=%zu  MCE=%.2fs\n", regions.size(), secs(T1, T2));
    if (regions.empty()) { printf("[rn] no region >= s.\n"); return 0; }

    int maxMC = 0; for (auto &R : regions) maxMC = max(maxMC, (int)R.size());
    build_ncr(maxMC + 2);

    // ---- r-MERGEABLE DIRECT-ASSIGN (V3LM Step 1b; env SCT_NO_RMERGE to disable) ----
    // A region M is fully-mergeable iff it shares < r vertices with EVERY other
    // region. Then every r-clique in M is |host|=1 (no other region holds an
    // r-subset of it), so M behaves as an isolated clique: all its r-cliques
    // share support C(|M|-r,s-r) and peel together at that level => core =
    // C(|M|-r,s-r), count = C(|M|,r). Direct-assign these (closed form) and
    // REMOVE M from the SCT/peel. M's r-cliques never appear in any active
    // region (|host|=1) and never sit inside an active region's s-clique (an
    // active s-clique meets M in < r vertices), so this is exact + no double
    // count. The SCT peel then handles only the OVERLAPPING (non-mergeable)
    // regions -> on sparse graphs that is a tiny fraction of the patterns.
    std::map<double,double> directCoreDist;
    long long nMergeable = 0, nMergedRC = 0;
    auto Trm0 = Clock::now();
    if (!getenv("SCT_NO_RMERGE")) {
        int nRall = (int)regions.size();
        // ---- r-mergeable via r-CLIQUE sharing (replaces the O(Σ_v deg_R(v)^2) pairwise
        // overlap, which explodes on hub vertices: cit-Patents 6.0e10 ops / >128s). A
        // region M is mergeable iff every r-subset of M is unique to M (contained in no
        // other region). Enumerate (sorted r-subset, region) pairs, sort by the r-tuple,
        // and any r-subset owned by >=2 regions marks those regions NON-mergeable. Exact
        // (compares the actual r-tuples) => bit-identical. Cost O(Σ_M C(|M|,r) log). ----
        vector<char> mergeable(nRall, 1);              // assume mergeable; cleared if shared
        long long Nsub = 0; bool tooBig = false;
        for (auto &M : regions) { int W = (int)M.size();
            if (W >= r) { double c = C(W, r); if (c > 3e8 || Nsub > 300000000LL) { tooBig = true; break; } Nsub += (long long)llround(c); } }
        // pick the CHEAPER method per graph (they have opposite weaknesses):
        //   pairwise  = O(Σ_v deg_R(v)^2)   -- bad on hub-skew (cit-Patents 6e10)
        //   r-clique  = O(Σ_M C(|M|,r) log) -- bad on dense large regions (ca-HepPh r=4, 28s)
        // Estimate Σdeg^2 (cheap) and use pairwise unless r-clique has >=16x fewer ops
        // (the sort costs ~16x more per element). Either way the RESULT is identical.
        long long costPw = 0;
        { vector<int> dc(g.n, 0); for (auto &M : regions) for (int v : M) dc[v]++;
          for (int v = 0; v < g.n; v++) costPw += (long long)dc[v] * (long long)dc[v]; }
        bool rmPairwise = (getenv("SCT_RM_PAIRWISE") != nullptr) || tooBig || (costPw <= Nsub * 16);
        if (!rmPairwise) {
            vector<int> sub; sub.reserve((size_t)Nsub * (size_t)r);   // flattened sorted r-tuples
            vector<int> reg; reg.reserve((size_t)Nsub);               // owning region per tuple
            vector<int> Rs, cc(r);
            for (int M = 0; M < nRall; M++) {
                int W = (int)regions[M].size();
                if (W < r) continue;
                Rs = regions[M]; std::sort(Rs.begin(), Rs.end());     // canonical tuple order
                for (int k = 0; k < r; k++) cc[k] = k;
                while (true) {                                        // iterate r-subsets of Rs
                    for (int k = 0; k < r; k++) sub.push_back(Rs[cc[k]]);
                    reg.push_back(M);
                    int k = r - 1; while (k >= 0 && cc[k] == W - r + k) k--;
                    if (k < 0) break;
                    cc[k]++; for (int kk = k + 1; kk < r; kk++) cc[kk] = cc[kk - 1] + 1;
                }
            }
            long long N = (long long)reg.size();
            vector<int> idx((size_t)N); for (long long i = 0; i < N; i++) idx[(size_t)i] = (int)i;
            const int *sp = sub.data();
            std::sort(idx.begin(), idx.end(), [&](int a, int b) {
                const int *pa = sp + (size_t)a * r, *pb = sp + (size_t)b * r;
                for (int k = 0; k < r; k++) if (pa[k] != pb[k]) return pa[k] < pb[k];
                return false;
            });
            for (long long i = 0; i < N; ) {                         // scan groups of equal r-tuple
                long long j = i + 1; const int *pi = sp + (size_t)idx[(size_t)i] * r;
                while (j < N) { const int *pj = sp + (size_t)idx[(size_t)j] * r;
                    bool eq = true; for (int k = 0; k < r; k++) if (pi[k] != pj[k]) { eq = false; break; }
                    if (!eq) break; j++; }
                if (j - i >= 2) for (long long k = i; k < j; k++) mergeable[reg[idx[(size_t)k]]] = 0;
                i = j;
            }
            if (getenv("RM_DBG")) fprintf(stderr, "[rm-dbg] r-clique: Nsub=%lld done=%.2fs\n", N, secs(Trm0, Clock::now()));
        } else {
            // fallback (only on r-subset blowup = clique-explosion graphs, out of scope
            // for everyone incl CND): original O(Σ_v deg_R(v)^2) pairwise overlap.
            std::fill(mergeable.begin(), mergeable.end(), 0);
            vector<vector<int>> vr(g.n);
            for (int i = 0; i < nRall; i++) for (int v : regions[i]) vr[v].push_back(i);
            vector<int> cnt(nRall, 0); vector<int> dirty; dirty.reserve(256);
            for (int M = 0; M < nRall; M++) {
                for (int v : regions[M]) for (int o : vr[v]) if (o != M) { if (cnt[o] == 0) dirty.push_back(o); cnt[o]++; }
                int mx = 0; for (int o : dirty) { if (cnt[o] > mx) mx = cnt[o]; cnt[o] = 0; }
                dirty.clear();
                if (mx < r) mergeable[M] = 1;
            }
        }
        vector<vector<int>> active; active.reserve(nRall);   // <= nRall non-mergeable regions
        for (int M = 0; M < nRall; M++) {
            if (mergeable[M]) {
                int N = (int)regions[M].size();
                double cv = (N >= (int)s) ? C(N - r, s - r) : 0.0;
                directCoreDist[cv] += C(N, r);
                nMergeable++; nMergedRC += (long long)llround(C(N, r));
            } else active.push_back(std::move(regions[M]));
        }
        regions = std::move(active);
        printf("[rn] r-mergeable: %lld regions direct (%lld r-cliques); active=%zu  %.2fs\n",
               nMergeable, nMergedRC, regions.size(), secs(Trm0, Clock::now()));
        fflush(stdout);
        if (regions.empty()) {
            double mx = 0; for (auto &kv : directCoreDist) mx = max(mx, kv.first);
            printf("[sct-peel] Max core: %.0f\n", mx);
            for (auto &kv : directCoreDist) printf("core=%.0f count=%.0f\n", kv.first, kv.second);
            return 0;
        }
    }

    // vtx -> sorted region ids
    int nR = (int)regions.size();
    vector<vector<int>> vtxR(g.n);
    for (int i = 0; i < nR; i++) for (int v : regions[i]) vtxR[v].push_back(i);
    for (int v = 0; v < g.n; v++) sort(vtxR[v].begin(), vtxR[v].end());

    // classes: vertices grouped by their (sorted) region-id vector
    unordered_map<string, int> profKey; profKey.reserve(g.n);   // <= g.n distinct profiles; kill rehashing
    vector<int> classOf(g.n, -1);
    vector<vector<int>> classRegions; classRegions.reserve(g.n);   // class -> sorted region ids (its profile)
    vector<int> classSize; classSize.reserve(g.n);
    auto keyOf = [](const vector<int> &p) {
        string k; k.reserve(p.size() * 4);
        for (int x : p) { k.append((const char *)&x, 4); }
        return k;
    };
    for (int v = 0; v < g.n; v++) {
        if (vtxR[v].empty()) continue;  // not in any region
        string k = keyOf(vtxR[v]);
        auto it = profKey.find(k);
        int c;
        if (it == profKey.end()) {
            c = (int)classRegions.size();
            profKey[k] = c;
            classRegions.push_back(vtxR[v]);
            classSize.push_back(0);
        } else c = it->second;
        classOf[v] = c;
        classSize[c]++;
    }
    int nC = (int)classRegions.size();
    auto T3 = Clock::now();
    printf("[rn] classes=%d  build(vtxR+class)=%.2fs\n", nC, secs(T2, T3));
    memCk("after-classes(pats not yet)");

    // region -> sorted class ids present in it (each class wholly in/out)
    vector<vector<int>> regionClasses(nR);
    for (int c = 0; c < nC; c++)
        for (int rid : classRegions[c]) regionClasses[rid].push_back(c);
    for (int i = 0; i < nR; i++) sort(regionClasses[i].begin(), regionClasses[i].end());

    // NOTE: the peel engine enumerates ALL patterns (incl |host|=1) for
    // correctness, and leafCount/interClasses require regionClasses[i] to
    // stay SORTED by class id (binary search). So we keep the sorted order
    // here -- no non-safe/safe stable_partition (that optimization, from
    // region_native.cpp's initial-support path, would break the sort and
    // silently corrupt every support_count).

    // adversarial invariant: |region| == sum of its class sizes
    for (int i = 0; i < nR; i++) {
        long tot = 0; for (int c : regionClasses[i]) tot += classSize[c];
        if (tot != (long)regions[i].size()) {
            printf("[rn] INVARIANT FAIL: region %d |M|=%zu sum-class=%ld\n",
                   i, regions[i].size(), tot);
            return 2;
        }
    }

    // ---- region union-count: |s-cliques >=R0 in union of `regs`| ----
    // regs = list of region ids (host set). Work at class granularity:
    // intersection of regions = classes present in ALL; |.| = sum sizes.
    // B&B: pick first region M; |A1| + |rest| - |rest∩M|; prune size<s,
    // dominance (region whose class-set subsumes another). R0 occupies r
    // of the vertices inside the intersection, so each term = C(|.|-r,s-r).
    // We represent a "region" in the recursion by its class-id set
    // (sorted vector) and cache its vertex-size.
    struct Node { vector<int> classes; int vsize; };
    auto classesSize = [&](const vector<int> &cs) {
        long t = 0; for (int c : cs) t += classSize[c]; return (int)t;
    };
    auto interClasses = [&](const vector<int> &a, const vector<int> &b) {
        vector<int> out; out.reserve(std::min(a.size(), b.size()));   // |a∩b| <= min; hot in unionAlive
        size_t i = 0, j = 0;
        while (i < a.size() && j < b.size()) {
            if (a[i] < b[j]) i++;
            else if (a[i] > b[j]) j++;
            else { out.push_back(a[i]); i++; j++; }
        }
        return out;
    };
    // ===================== PEEL ENGINE (region-native) =====================
    // sup(tau) = alive s-clique witnesses >= tau in the union of tau's host
    // regions. We reuse the verified union-count B&B, but each region-
    // intersection leaf returns ccpath::support_count (alive count with the
    // peeled-tuple forbidden antichain) instead of a plain binomial. A
    // witness (composition y) is dead iff y >= some peeled tuple f (it
    // contains a peeled orbit). This is the forbidden-antichain rule,
    // validated bit-exact in test_region_forbidden.cpp.
    using ccpath::CCPath; using ccpath::Vec;
    (void)build_ncr;  // already built above

    struct Pat {
        vector<int> host;            // sorted region ids (|host|>=1)
        vector<pair<int,int>> comp;  // (classId, mult) sorted by classId, sum=r
        vector<int> classSet;        // sorted class ids of comp (subset tests) -- SCT_VERIFY-only, freed after build
        int hostSz = 0;              // |host| kept after host's region-id list is freed (peel needs only the count)
        double sup = 0; double core = -1;
        long long mult = 1;          // # actual r-cliques in this orbit
        bool alive = true; long long key = -1;
    };
    vector<Pat> pats;
    vector<int> peeled;              // pattern indices in peel order (forbidden)
    size_t maxForb = 0;              // diagnostic: largest forbidden antichain seen

    // alive witnesses >= comp inside the clique whose class-set is Cs.
    auto leafCount = [&](const vector<int> &Cs,
                         const vector<pair<int,int>> &comp) -> double {
        int M = (int)Cs.size();
        Vec nn((size_t)M), zeros((size_t)M, 0);
        for (int i = 0; i < M; i++) nn[i] = (int16_t)classSize[Cs[i]];
        CCPath p = CCPath::initial(zeros, nn, s);
        Vec b((size_t)M, 0);
        for (auto &cm : comp) {
            int pos = (int)(lower_bound(Cs.begin(), Cs.end(), cm.first) - Cs.begin());
            b[(size_t)pos] = (int16_t)cm.second;   // comp.classes ⊆ Cs guaranteed
        }
        // forbidden = peeled tuples whose class-set ⊆ Cs (their witnesses can
        // live in this clique). component-max IE handled by support_count.
        for (int pi : peeled) {
            const auto &Q = pats[pi];
            size_t i = 0, j = 0; bool sub = true;          // Q.classSet ⊆ Cs ?
            while (i < Q.classSet.size()) {
                while (j < Cs.size() && Cs[j] < Q.classSet[i]) j++;
                if (j >= Cs.size() || Cs[j] != Q.classSet[i]) { sub = false; break; }
                i++; j++;
            }
            if (!sub) continue;
            Vec fv((size_t)M, 0);
            for (auto &cm : Q.comp) {
                int pos = (int)(lower_bound(Cs.begin(), Cs.end(), cm.first) - Cs.begin());
                fv[(size_t)pos] = (int16_t)cm.second;
            }
            ccpath::insert_antichain(p.forbidden, fv);
        }
        if (p.forbidden.size() > maxForb) maxForb = p.forbidden.size();
        return ccpath::support_count(p, b, ccpath_ncr);
    };

    // forbidden-aware union over host regions. IE holds for the alive-witness
    // count (it is a set cardinality); dominance still valid because aliveness
    // is a per-witness property, so A.classes ⊆ B.classes => alive(A) ⊆ alive(B).
    std::function<double(vector<Node>, const vector<pair<int,int>> &)> unionAlive =
        [&](vector<Node> B, const vector<pair<int,int>> &comp) -> double {
        if (B.empty()) return 0.0;
        Node M = std::move(B.back()); B.pop_back();
        double here = leafCount(M.classes, comp);
        vector<Node> inter; inter.reserve(B.size());   // <= |B| surviving intersections
        for (auto &N : B) {
            vector<int> cs = interClasses(M.classes, N.classes);
            int vs = classesSize(cs);
            if (vs >= s) inter.push_back({std::move(cs), vs});
        }
        if (inter.size() > 1) {
            sort(inter.begin(), inter.end(),
                 [](const Node &a, const Node &b){ return a.classes.size() > b.classes.size(); });
            vector<Node> keep; keep.reserve(inter.size());   // <= |inter| non-dominated
            for (auto &nd : inter) {
                bool dom = false;
                for (auto &k : keep) if (nd.classes.size() <= k.classes.size()) {
                    size_t i=0,j=0; bool sub=true;
                    while (i<nd.classes.size()){
                        while(j<k.classes.size()&&k.classes[j]<nd.classes[i])j++;
                        if(j>=k.classes.size()||k.classes[j]!=nd.classes[i]){sub=false;break;}
                        i++;j++;
                    }
                    if (sub){dom=true;break;}
                }
                if (!dom) keep.push_back(std::move(nd));
            }
            inter = std::move(keep);
        }
        double interU = unionAlive(std::move(inter), comp);
        double restU  = unionAlive(std::move(B), comp);
        return here + restU - interU;
    };
    auto suppOf = [&](const Pat &P) -> double {
        vector<Node> H; H.reserve(P.host.size());
        for (int rid : P.host) H.push_back({regionClasses[rid], (int)regions[rid].size()});
        return unionAlive(std::move(H), P.comp);
    };

    // ---- patterns via GLOBAL host: emit (pattern, region) incidences, group by pattern ----
    // host(P) = {M : P's classes ⊆ regionClasses[M]} = EXACTLY the regions that enumerate P.
    // So emitting (comp-key, M) per (region, r-multiset) and grouping by comp-key yields the
    // host directly -- NO ∩classRegions[c] (kills the per-pattern K² hub blowup; the old
    // canonical-home dedup is now just the grouping). comp is padded to r (class,mult) pairs
    // with sentinel class=nC for a fixed-width sort key. Bit-identical (same pats; corehash is
    // the order-independent core distribution). SCT_PE_PERREGION forces the old per-region path.
    bool pePerRegion = getenv("SCT_PE_PERREGION") != nullptr;
    if (!pePerRegion) {
        bool peDbg = getenv("PE_DBG") != nullptr;
        const int W2 = 2 * r;
        vector<int> rec; vector<int> recReg;               // flat padded comp-keys / region per incidence
        // size the two flat buffers exactly with a cheap count-only recursion. rec
        // grows to Ninc*W2 ints (multi-GB on dense high-s graphs); reserving avoids
        // the ~2x deep-copy of geometric growth AND the transient ~1.5x peak (which
        // risks OOM on the big cells), at the cost of one extra bare recursion pass.
        {
            long long nIncEst = 0;
            auto cnt = [&](auto &&self, int idx, const vector<int> &cls, int rem) -> void {
                if (rem == 0) { nIncEst++; return; }
                for (int i = idx; i < (int)cls.size(); i++) { int c = cls[i], maxj = min(rem, classSize[c]);
                    for (int j = 1; j <= maxj; j++) self(self, i + 1, cls, rem - j); }
            };
            for (int i = 0; i < nR; i++) cnt(cnt, 0, regionClasses[i], r);
            rec.reserve((size_t)nIncEst * W2); recReg.reserve((size_t)nIncEst);
        }
        vector<pair<int,int>> cur; int curRid = 0;
        auto enumE = [&](auto &&self, int idx, const vector<int> &cls, int rem) -> void {
            if (rem == 0) {
                for (auto &cm : cur) { rec.push_back(cm.first); rec.push_back(cm.second); }
                for (int p = (int)cur.size(); p < r; p++) { rec.push_back(nC); rec.push_back(0); }   // pad
                recReg.push_back(curRid);
                return;
            }
            for (int i = idx; i < (int)cls.size(); i++) {
                int c = cls[i], maxj = min(rem, classSize[c]);
                for (int j = 1; j <= maxj; j++) { cur.push_back({c, j}); self(self, i + 1, cls, rem - j); cur.pop_back(); }
            }
        };
        for (int i = 0; i < nR; i++) { curRid = i; cur.clear(); enumE(enumE, 0, regionClasses[i], r); }
        long long Ninc = (long long)recReg.size();
        if (peDbg) fprintf(stderr, "[pe-dbg] emitted=%lld incidences t=%.2fs\n", Ninc, secs(T3, Clock::now()));
        vector<int> ord((size_t)Ninc); for (long long i = 0; i < Ninc; i++) ord[(size_t)i] = (int)i;
        const int *rp = rec.data();
        if (getenv("SCT_PE_STDSORT")) {                    // ablation: comparison sort (the old path)
            std::sort(ord.begin(), ord.end(), [&](int a, int b) {
                const int *pa = rp + (size_t)a * W2, *pb = rp + (size_t)b * W2;
                for (int k = 0; k < W2; k++) if (pa[k] != pb[k]) return pa[k] < pb[k];
                return recReg[a] < recReg[b];              // region tiebreak -> host comes out sorted
            });
        } else {
            // LSD radix sort of `ord` by the key (k0..k_{W2-1}, region). Stable
            // counting sort per column, least-significant first: region, then the
            // comp columns high..low. Produces the IDENTICAL total order to the
            // comparator above (region is the last tiebreak <=> least-significant
            // column), so the grouping and the resulting pats/host are bit-identical.
            // Replaces the O(N log N) comparison sort -- each compare did 2 random
            // W2-int gathers into the multi-GB `rec` -- with stable counting passes
            // (no comparisons, no log factor). The per-pass random gather is the cost,
            // so FEWER PASSES wins: by DEFAULT we pack each (class,mult) pair into one
            // digit class<<mb | mult (mult<=r fits in mb bits, class is the high part
            // => order-preserving, sentinel (nC,0) still sorts last) -> r+1 passes
            // instead of 2r+1. Measured: unpacked radix is 3.5x on r=3 but only ~1x at
            // r=6 (13 passes ~ comparison cost); packing halves the passes. The result
            // lands in ord2 (odd #passes) or ord (even); the pointer check swaps in O(1).
            // ABLATION: SCT_PE_RADIX_UNPACKED = one digit per int (2r+1 passes).
            int mb = 1; while ((1 << mb) <= r) mb++;        // mult in [0,r] fits in mb bits
            bool unpacked = getenv("SCT_PE_RADIX_UNPACKED") != nullptr;
            vector<int> ord2((size_t)Ninc);
            long long Mreg = (long long)nR - 1;
            long long Mcol = unpacked ? (long long)nC : ((long long)nC << mb);
            long long M = Mcol > Mreg ? Mcol : Mreg;        // widest column value
            vector<int> cnt2((size_t)M + 2);
            int *src = ord.data(), *dst = ord2.data();
            auto digit = [&](int idx, int col, bool packed) -> long long {   // col<0 => region
                if (col < 0) return recReg[idx];
                const int *e = rp + (size_t)idx * W2;
                return packed ? (((long long)e[2 * col] << mb) | e[2 * col + 1]) : e[col];
            };
            auto pass = [&](int col, bool packed) {
                std::fill(cnt2.begin(), cnt2.end(), 0);
                for (long long i = 0; i < Ninc; i++) cnt2[(size_t)digit(src[i], col, packed) + 1]++;
                for (long long b = 1; b <= M + 1; b++) cnt2[(size_t)b] += cnt2[(size_t)b - 1];   // start offsets
                for (long long i = 0; i < Ninc; i++) { int idx = src[i]; dst[cnt2[(size_t)digit(idx, col, packed)]++] = idx; }
                std::swap(src, dst);
            };
            pass(-1, false);                                // region = least significant (last tiebreak)
            if (unpacked) for (int k = W2 - 1; k >= 0; k--) pass(k, false);   // k[W2-1]..k[0]
            else          for (int c = r - 1;  c >= 0; c--) pass(c, true);    // packed pair r-1..0
            if (src == ord2.data()) ord.swap(ord2);         // land the result in ord (O(1) pointer swap)
        }
        for (long long i = 0; i < Ninc; ) {                // group consecutive equal comp-key -> one Pat
            const int *pi = rp + (size_t)ord[(size_t)i] * W2;
            long long j = i + 1;
            while (j < Ninc) { const int *pj = rp + (size_t)ord[(size_t)j] * W2;
                bool eq = true; for (int k = 0; k < W2; k++) if (pi[k] != pj[k]) { eq = false; break; }
                if (!eq) break; j++; }
            Pat P;
            for (int k = 0; k < W2; k += 2) { if (pi[k] == nC) break; P.comp.push_back({pi[k], pi[k + 1]}); }
            for (auto &cm : P.comp) P.classSet.push_back(cm.first);          // comp already class-id sorted
            P.host.reserve((size_t)(j - i));
            for (long long t = i; t < j; t++) P.host.push_back(recReg[ord[(size_t)t]]);   // sorted (region tiebreak)
            long long mu = 1; for (auto &cm : P.comp) mu *= (long long)llround(C(classSize[cm.first], cm.second));
            P.mult = mu;
            pats.push_back(std::move(P));
            i = j;
        }
        if (peDbg) fprintf(stderr, "[pe-dbg] DONE pats=%zu t=%.2fs\n", pats.size(), secs(T3, Clock::now()));
    } else {
        // old per-region path (host = ∩classRegions[c], canonical-home dedup); O(K²) on hubs.
        int curRid = 0; vector<pair<int,int>> cur;
        std::function<void(int,const vector<int>&,int)> enumAll = [&](int idx, const vector<int> &cls, int rem) {
            if (rem == 0) {
                vector<int> host = classRegions[cur[0].first];
                for (size_t i = 1; i < cur.size() && !host.empty(); i++) host = interClasses(host, classRegions[cur[i].first]);
                if (host.empty() || host[0] != curRid) return;
                Pat P; P.host = host; P.comp = cur;
                for (auto &cm : cur) P.classSet.push_back(cm.first);
                long long mu = 1; for (auto &cm : cur) mu *= (long long)llround(C(classSize[cm.first], cm.second));
                P.mult = mu; pats.push_back(std::move(P));
                return;
            }
            for (int i = idx; i < (int)cls.size(); i++) { int c = cls[i], maxj = min(rem, classSize[c]);
                for (int j = 1; j <= maxj; j++) { cur.push_back({c, j}); enumAll(i + 1, cls, rem - j); cur.pop_back(); } }
        };
        for (int i = 0; i < nR; i++) { curRid = i; enumAll(0, regionClasses[i], r); }
    }
    // EMPIRICAL TEST (env SCT_DIRECTBIN_ALL_HOST1): direct-bin ALL |host|=1
    // patterns (not only fully-mergeable regions), peel only |host|>=2. GATE
    // vs brute decides correctness. Suspected WRONG: a |host|=1 pattern's
    // witness s-clique can contain a |host|>=2 r-clique, so peeling the
    // |host|=1 pattern lowers a |host|>=2 pattern's support -> skipping it
    // corrupts the |host|>=2 cores.
    if (getenv("SCT_DIRECTBIN_ALL_HOST1")) {
        vector<Pat> keep; long long nDB=0, nDBrc=0;
        for (auto &P : pats) {
            if (P.host.size() == 1) {
                int N = (int)regions[P.host[0]].size();
                double cv = (N >= (int)s) ? C(N-r, s-r) : 0.0;
                directCoreDist[cv] += (double)P.mult; nDB++; nDBrc += P.mult;
            } else keep.push_back(std::move(P));
        }
        pats = std::move(keep);
        printf("[rn-peel] DIRECTBIN_ALL_HOST1: %lld patterns direct (%lld r-cliques); peel patterns=%zu\n",
               nDB, nDBrc, pats.size());
    }
    auto T4 = Clock::now();
    long long totalRCliques = 0; for (auto &P : pats) totalRCliques += P.mult;
    memCk("after-pattern-enum(pats)");
    printf("[rn-peel] patterns=%zu  r-cliques=%lld  enum=%.2fs\n",
           pats.size(), totalRCliques, secs(T3, T4));
    fflush(stdout);
    if (pats.empty()) { printf("[rn-peel] no patterns.\n"); return 0; }

    // ===================================================================
    //  CLASS-SCT PEEL  (replaces the region-IE peel of region_native_peel)
    // ===================================================================
    // Step 1+2: SPARSE quotient graph + SCALABLE class-SCT. Handles large nC
    // (com-dblp nC=123k, web-it 57k) where a dense C x C matrix would be ~15GB.
    // edge(i,j) iff classes i,j co-occur in some region. scalableBuildClassSCT
    // (degeneracy-seeded, sparse) returns COMPACT leaves (classIds sorted +
    // local h/n/ell/u), validated == the dense builder over 50000 trials.
    auto Tqg0 = Clock::now();
    vector<int> qw(nC); for (int c = 0; c < nC; c++) qw[c] = classSize[c];
    vector<vector<int>> qadj(nC);
    for (int M = 0; M < nR; M++) {
        const auto &rc = regionClasses[M];
        for (size_t a = 0; a < rc.size(); a++)
            for (size_t b = a + 1; b < rc.size(); b++) {
                qadj[rc[a]].push_back(rc[b]); qadj[rc[b]].push_back(rc[a]);
            }
    }
    for (int c = 0; c < nC; c++) {
        auto &v = qadj[c]; std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    auto baseLeaves = classsct_scalable::scalableBuildClassSCT(nC, qw, qadj, s);  // COMPACT
    auto Tqg1 = Clock::now();
    memCk("after-SDCT-build(slotPaths)");
    printf("[sct] quotient nC=%d  base-leaves=%zu  build=%.2fs\n",
           nC, baseLeaves.size(), secs(Tqg0, Tqg1));
    fflush(stdout);

    // m as a GLOBAL class vector (length nC): m_vec[c] = mult of class c in P.
    auto compToVec = [&](const vector<pair<int,int>> &comp) -> Vec {
        Vec v((size_t)nC, 0);
        for (auto &cm : comp) v[(size_t)cm.first] = (int16_t)cm.second;
        return v;
    };

    // region-IE init support is an INDEPENDENT reference (gate G2a); the PEEL runs on the SCT
    // support (sctSupport below). Computing region-IE suppOf per pattern was 99s on ca-AstroPh
    // (THE maps-phase bottleneck), so it is now SCT_VERIFY-only; production takes P.sup from the
    // SCT in the G2a loop below. Bit-identical (region-IE == SCT, integer-valued; verified G2a).
    bool verifyIE = getenv("SCT_VERIFY") != nullptr;
    if (verifyIE) {
        for (auto &P : pats) { P.sup = suppOf(P); }
        // sanity: SCT total s-cliques == region-IE total (each s-clique has C(s,r) r-subcliques).
        double sclSCT = 0;
        for (auto &lf : baseLeaves) { Vec zb((size_t)lf.m(), 0); sclSCT += ccpath::support_count(lf, zb, ccpath_ncr); }
        double sclIE = 0; for (auto &P : pats) sclIE += (double)P.mult * P.sup;
        sclIE /= C(s, r);
        printf("[sct] total s-cliques: class-SCT=%.0f  region-IE=%.0f  %s\n",
               sclSCT, sclIE, fabs(sclSCT - sclIE) < 0.5 ? "[OK]" : "[MISMATCH]");
        fflush(stdout);
    }
    // §96b: free build-only heavy per-pattern fields. classSet is read ONLY by suppOf (SCT_VERIFY, just above); host's
    // region-id list is needed ONLY at build (directBin) + suppOf -- the peel uses only |host| (kept as hostSz). Frees
    // ~classSet + ~host of per-pattern incidence (e.g. ca-AstroPh 4,5: 73MB + 85MB).
    for (auto &P : pats) { P.hostSz = (int)P.host.size(); vector<int>().swap(P.classSet); vector<int>().swap(P.host); }

    // Step 3: compaction + pattern<->leaf maps.
    int nLeaf = (int)baseLeaves.size();
    vector<char> hasH2(nLeaf, 0);                     // §96b: filled during enumLP (from hostSz), so leafPats can be dropped
    // NB: no full-C per-pattern vector (would be length-nC x #patterns = TB on
    // com-dblp nC=123k). Patterns stay SPARSE (comp); compToLocal maps a comp
    // straight into a leaf's local dimension via binary search on supC.

    // ---- PER-LEAF COMPACTION (FIRST, so the map-build confirm is cheap) ----
    // Each leaf touches only its support classes {c: n[c]||h[c]} (<=~10);
    // compact to that local dimension so support_count / controlled_split run
    // in the tiny local dim. CCPathCore/ClassSCT untouched.
    vector<vector<int>> supC(nLeaf);                  // leaf -> sorted support classes
    vector<vector<CCPath>> slotPaths(nLeaf);          // compact (split) path set
    for (int lid = 0; lid < nLeaf; lid++) {
        const CCPath &lf = baseLeaves[lid];
        // scalable leaves are ALREADY compact: classIds (sorted global) parallel
        // to the local h/n/ell/u dimensions. Use them directly as supC/slotPath.
        supC[lid].assign(lf.classIds.begin(), lf.classIds.end());
        slotPaths[lid].push_back(lf);
    }
    // map a global-class vector to leaf lid's local dimension (b-vector).
    auto toLocal = [&](int lid, const Vec &gv) -> Vec {
        const vector<int> &sc = supC[lid];
        Vec b((size_t)sc.size(), 0);
        for (size_t i = 0; i < sc.size(); i++) b[i] = gv[(size_t)sc[i]];
        return b;
    };
    // map a SPARSE comp straight into leaf lid's local dim (no full-C vector).
    // supC[lid] is sorted, so each comp class is binary-searched. A pattern's
    // classes are a subset of the leaf's classes when it is hosted there.
    auto compToLocal = [&](const vector<pair<int,int>> &comp, int lid) -> Vec {
        const vector<int> &sc = supC[lid];
        Vec b((size_t)sc.size(), 0);
        for (auto &cm : comp) {
            int pos = (int)(std::lower_bound(sc.begin(), sc.end(), cm.first) - sc.begin());
            if (pos < (int)sc.size() && sc[pos] == cm.first) b[(size_t)pos] = (int16_t)cm.second;
        }
        return b;
    };

    // ---- pattern<->leaf maps via PER-LEAF ENUMERATION ----
    // Every r-multiset over a leaf's classes is a valid r-clique = a registered
    // pattern (the leaf's classes form a clique). So enumerate r-multisets PER
    // LEAF and look them up in a comp->patternId index -> O(sum of pattern-leaf
    // incidences), independent of how many leaves a class sits in. (The
    // per-pattern candidate scan was 33s on ca-CondMat: 29544 leaves x a class
    // in many leaves. This is the production enumCb approach.) The host is
    // confirmed on the COMPACT leaf (local dim) so each support_count is small.
    vector<vector<int>> patLeaves(pats.size());
    vector<vector<int>> leafPats(nLeaf);
    vector<vector<Vec>> pbLocal(pats.size());         // parallel to patLeaves[pi]
    vector<vector<int16_t>> leafFlat(nLeaf);          // L3 CSR: per-leaf FLAT footprints (size leafPats[lid].size()*Mloc)
    // integer rolling-hash key (was a std::string compKey: a per-r-multiset heap
    // alloc + string hash, the bulk of the maps phase on dense graphs). Hash
    // collisions are resolved by comparing the actual comp, so lookup stays exact.
    // pattern lookup: flat open-addressing hash comp->patternId. (Was a sorted
    // (compHash, pi) array + std::lower_bound binary search: ~25 cache-cold
    // comparisons over the 23.2M-pattern array, measured as 61% of the maps phase on
    // web-Google(6,8) -- 56.5s of a 92s enumLP. An unordered_map BUILDS pathologically
    // (node-per-element, 99.7s for 848k patterns on ca-AstroPh); a single contiguous
    // table builds in one linear pass and looks up in ~1-2 probes. Collisions are
    // resolved by comparing the actual comp, so lookup stays bit-exact.)
    auto compHash = [](const vector<pair<int,int>> &comp) -> uint64_t {
        uint64_t h = 1469598103934665603ULL;
        for (auto &cm : comp) {
            h = (h ^ ((uint64_t)(uint32_t)cm.first + 1)) * 1099511628211ULL;
            h = (h ^ ((uint64_t)(uint32_t)cm.second + 1)) * 1099511628211ULL;
        }
        return h;
    };
    size_t hcap = 16; while (hcap < pats.size() * 2) hcap <<= 1;   // load factor < 0.5
    size_t hmask = hcap - 1;
    vector<int> htab(hcap, -1);
    for (int pi = 0; pi < (int)pats.size(); pi++) {
        size_t idx = compHash(pats[pi].comp) & hmask;
        while (htab[idx] != -1) idx = (idx + 1) & hmask;          // linear probe
        htab[idx] = pi;
    }
    long long mapInc = 0;                              // maps-dbg: pattern-leaf incidences
    bool mapsDbg = getenv("MAPS_DBG") != nullptr;
    // Phase 2: when set, DON'T store the ~200M pbLocal/leafPatLocB Vec payloads; each
    // blocal is recomputed on demand via localB (proven equivalent in Phase 1, gate
    // SCT_MAPS_VALIDATE). The int maps patLeaves/leafPats are always stored. Default
    // OFF keeps the stored path for A/B + corehash cross-check (both must be bit-identical).
    bool mapsRecompute = getenv("SCT_MAPS_RECOMPUTE") != nullptr;
    // ADAPTIVE memory (§79 A): recompute the COLD pbLocal (pattern->leaf footprints, read ~O(incidences)) but KEEP
    // the HOT leafPatLocB (leaf->pattern, read O(incidences*candidates) per Q-lookup) stored. Drops ~half the maps
    // memory at near-zero time cost -- vs SCT_MAPS_RECOMPUTE (full) which recomputes BOTH and pays +17% on the hot map.
    bool mapsRecomputePB = getenv("SCT_MAPS_NO_RECOMPUTE_PB") == nullptr;   // DEFAULT ON (free; §81 RSS-confirmed)
    bool recomputePB = mapsRecompute || mapsRecomputePB;   // pbLocal not stored / recomputed
    bool ondemand = getenv("SCT_ONDEMAND") != nullptr;     // §94 Stage 2: compute patLeaves on-demand, do NOT store it
    // LEVER 2 (§80): recompute the HOT leafPatLocB for the WIDEST leaves only (Mloc>=leafWmin) -- recovers the
    // OTHER ~half of the maps memory, but at a time cost on those (hot) leaves. leafWmin tunes the memory<->time
    // tradeoff (0 = off; combine with SCT_MAPS_RECOMPUTE_PB for free-half + adaptive-other-half). Mloc=supC[lid].size().
    int leafWmin = getenv("SCT_MAPS_LEAF_WMIN") ? atoi(getenv("SCT_MAPS_LEAF_WMIN")) : 0;
    vector<char> leafRecomp(nLeaf, 0);
    if (leafWmin > 0) for (int lid = 0; lid < nLeaf; lid++) leafRecomp[lid] = ((int)supC[lid].size() >= leafWmin);
    if (mapsDbg) fprintf(stderr, "[maps-dbg] patIdx-build=%.2fs pats=%zu\n", secs(Tqg1, Clock::now()), pats.size());
    auto TmapE0 = Clock::now();
    {
        vector<int> lcs, lcap; vector<pair<int,int>> cur; Vec blocal;
        cur.reserve(r);   // r-multiset depth <= r; one-time, kills startup reallocs
        // host-confirm: support_count(box,b)>0 iff the box {max(ell,b)<=y<=u, Σy=s}
        // is NONEMPTY -- every weight C(n-b,y-b) is a positive binomial, so the count
        // is >0 exactly when an integer point exists. With empty forbidden (the
        // pre-peel leaf box) this is an O(width) feasibility test, no DP. (Falls back
        // to support_count if forbidden is ever non-empty here.)
        auto hostFeasible = [](const CCPath &box, const Vec &b) -> bool {
            const int M = box.m(); int sl = 0, su = 0;
            for (int c = 0; c < M; c++) {
                int L = box.ell[c]; if ((int)b[c] > L) L = (int)b[c];
                int U = (int)box.u[c]; if (L > U) return false;
                sl += L; su += U;
            }
            return sl <= box.T && box.T <= su;
        };
        // self-recursive (Y-combinator) -> inlinable, no std::function indirection.
        // The local vector `blocal` is built INLINE during the enumeration (lcs is
        // parallel to the box's local dim, so class at enum position i sits at local
        // position i) -> no compToLocal, and the SAME blocal fills all four maps
        // (patLeaves/leafPats/pbLocal/leafPatLocB), eliminating the ~3 per-incidence
        // compToLocal remaps+allocs that dominated the maps phase on wide graphs.
        auto enumLP = [&](auto &&self, int lid, int idx, int rem) -> void {
            if (rem == 0) {
                uint64_t hk = compHash(cur);              // flat open-addressing probe (collision-safe)
                int pi = -1;
                for (size_t idx = hk & hmask; htab[idx] != -1; idx = (idx + 1) & hmask)
                    if (pats[htab[idx]].comp == cur) { pi = htab[idx]; break; }
                if (pi < 0) return;
                // confirm host on the compact leaf (filters m with no s-extension)
                const CCPath &box = slotPaths[lid][0];
                bool host = box.forbidden.empty() ? hostFeasible(box, blocal)
                                                  : (ccpath::support_count(box, blocal, ccpath_ncr) > 0.0);
                if (host) {
                    if (!ondemand) patLeaves[pi].push_back(lid);   // §94: on-demand recomputes patLeaves via class->leaves
                    if (pats[pi].hostSz >= 2) hasH2[lid] = 1;       // §96b: hasH2 computed here (host freed; leafPats droppable)
                    if (!(ondemand && (s - r) == 1)) leafPats[lid].push_back(pi);   // §96b: ondemand t=1 resolves Q via global hash
                    if (!recomputePB) pbLocal[pi].push_back(blocal);          // cold map: skip under PB or full
                    if (!mapsRecompute && !leafRecomp[lid] && !(ondemand && (s - r) == 1))   // §3b: ondemand t=1 resolves Q via global hash -> leafFlat unused
                        leafFlat[lid].insert(leafFlat[lid].end(), blocal.begin(), blocal.end());  // hot map (flat)
                    mapInc++;
                }
                return;
            }
            for (int i = idx; i < (int)lcs.size(); i++) {
                int mx = std::min(rem, lcap[i]);
                for (int j = 1; j <= mx; j++) {
                    cur.push_back({lcs[i], j}); blocal[(size_t)i] = (int16_t)j;
                    self(self, lid, i+1, rem-j);
                    cur.pop_back(); blocal[(size_t)i] = 0;
                }
            }
        };
        for (int lid = 0; lid < nLeaf; lid++) {
            const CCPath &cp = slotPaths[lid][0];          // compact leaf
            lcs = supC[lid];                               // global class ids (sorted)
            lcap.assign(cp.u.begin(), cp.u.end());         // local u, parallel to supC
            blocal.assign(lcs.size(), 0);                  // local b scratch (parallel to lcs)
            cur.clear(); enumLP(enumLP, lid, 0, r);
        }
    }
    if (mapsDbg) fprintf(stderr, "[maps-dbg] enumLP=%.2fs incidences=%lld (%.0f ns/inc)\n",
                         secs(TmapE0, Clock::now()), mapInc, mapInc ? 1e9*secs(TmapE0, Clock::now())/mapInc : 0.0);
    // ---- maps-recompute oracle: localB(pi,lid) reconstructs the stored blocal ----
    // The stored Vec payloads pbLocal[pi][k]/leafPatLocB[lid][t] are ~200M per-
    // incidence std::vector<int16_t> heap allocs (the dominant maps push_back cost).
    // But each blocal is fully determined by (pattern comp, leaf class layout):
    // blocal[i] = mult of leaf-class supC[lid][i] in pattern pi's comp -> a merge of
    // two sorted lists, recomputable on demand into a reused scratch (no alloc). This
    // lambda is the recompute; SCT_MAPS_VALIDATE proves it equals every stored payload
    // (the prerequisite for Phase 2: drop the storage and recompute at the consumers).
    auto localB = [&](int pi, int lid, Vec &out) {
        const vector<int> &sc = supC[lid];
        out.assign(sc.size(), 0);
        const auto &comp = pats[pi].comp;              // sorted (classId, mult)
        size_t i = 0, j = 0;
        while (i < sc.size() && j < comp.size()) {
            if (sc[i] < comp[j].first) i++;
            else if (sc[i] > comp[j].first) j++;
            else { out[(size_t)i] = (int16_t)comp[j].second; i++; j++; }
        }
    };
    (void)localB;
    // ---- L3 CSR: leafPatLocB is a PER-LEAF FLAT int16 array leafFlat[lid] (footprints all size Mloc=supC[lid]
    // size, concatenated in leafPats order). Footprint t = &leafFlat[lid][t*Mloc]. This removes the ~40B/footprint
    // per-Vec overhead (object + malloc header) that dominated narrow-leaf graphs, free (same data, flat layout).
    auto spanEq = [](const int16_t *p, int len, const Vec &v) -> bool {
        if ((int)v.size() != len) return false;
        for (int i = 0; i < len; i++) if (p[i] != (int16_t)v[i]) return false;
        return true;
    };
    auto hashSpan = [](const int16_t *p, int len) -> uint64_t {       // == hashVec(span) (same FNV)
        uint64_t h = 1469598103934665603ULL;
        for (int i = 0; i < len; i++) { h ^= (uint64_t)(uint16_t)p[i] + 1; h *= 1099511628211ULL; }
        return h;
    };
    // footprint of leaf-pattern t in leaf lid as a span (recompute into scr for full/wide-leaf, else CSR view).
    auto leafFP = [&](int lid, int t, int Mloc, Vec &scr) -> std::pair<const int16_t *, int> {
        if (mapsRecompute || leafRecomp[lid]) { localB(leafPats[lid][t], lid, scr); return {scr.data(), (int)scr.size()}; }
        return {leafFlat[lid].data() + (size_t)t * (size_t)Mloc, Mloc};
    };
    auto spanEqFP = [&](int lid, int t, int Mloc, Vec &scr, const Vec &v) -> bool {   // footprint(lid,t) == v ?
        auto fp = leafFP(lid, t, Mloc, scr); return spanEq(fp.first, fp.second, v);
    };
    (void)spanEq; (void)hashSpan; (void)leafFP; (void)spanEqFP;
    if (getenv("SCT_MAPS_VALIDATE") && !recomputePB) {   // prove recompute == stored (needs full storage), abort-on-mismatch
        long long chk = 0, bad = 0; Vec rb;
        for (int pi = 0; pi < (int)pats.size(); pi++)
            for (size_t k = 0; k < patLeaves[pi].size(); k++) {
                localB(pi, patLeaves[pi][k], rb);
                if (rb != pbLocal[pi][k]) { if (bad < 8) fprintf(stderr, "[maps-val] pbLocal MISMATCH pi=%d k=%d\n", pi, (int)k); bad++; }
                chk++;
            }
        for (int lid = 0; lid < nLeaf; lid++)
            for (size_t t = 0; t < leafPats[lid].size(); t++) {
                localB(leafPats[lid][t], lid, rb);
                int Mw = (int)supC[lid].size();
                if (!spanEq(leafFlat[lid].data() + (size_t)t * Mw, Mw, rb)) { if (bad < 8) fprintf(stderr, "[maps-val] leafFlat MISMATCH lid=%d t=%d\n", lid, (int)t); bad++; }
                chk++;
            }
        fprintf(stderr, "[maps-val] recompute vs stored: %lld checked, %lld bad\n", chk, bad);
        return bad ? 4 : 0;                            // validation-only run: skip the peel
    }
    // hasH2[lid]: does leaf lid host any |host|>=2 pattern? A |host|=1 pattern
    // whose leaves are ALL pure-|host|=1 removes nothing relevant when it peels
    // (its witnesses are M-exclusive, so they feed no |host|>=2 support), so its
    // source-peel can be skipped entirely. (Source-skip; SCT_NO_SKIP_H1 disables
    // both this and the target-skip.) leafPats order is enumeration order (the peel
    // accesses patterns by hash-mapped index, not order; cores bit-identical).
    // hasH2 is filled during enumLP (above) from hostSz -- correct even when leafPats is dropped (ondemand t=1).
    // The leaf->pattern maps are now built from classIds (via supC). The peel
    // works purely in local positions and NEVER reads CCPath::classIds /
    // tupleIdxs again (they are metadata; no CCPathCore algorithm touches them).
    // Strip them so the many path COPIES inside slotForbidDiff / controlled_split
    // during peeling don't drag a length-M classIds vector each time.
    for (int lid = 0; lid < nLeaf; lid++)
        for (auto &p : slotPaths[lid]) {
            p.classIds.clear(); p.classIds.shrink_to_fit();
            p.tupleIdxs.clear(); p.tupleIdxs.shrink_to_fit();
        }
    // PROBE (SCT_SLOT_REVERSE): is slot-path ORDER load-bearing for bit-identical?
    // Reverse every slot once; if corehash is unchanged, the support sums are
    // order-robust (exact-integer recovered, FP error < 0.5) -> the sub-linear slot
    // skip may reorder/mark-compact freely. If it changes, the skip must reproduce
    // the exact swap-remove order. Deterministic -> contention-insensitive.
    if (getenv("SCT_SLOT_REVERSE"))
        for (int lid = 0; lid < nLeaf; lid++) std::reverse(slotPaths[lid].begin(), slotPaths[lid].end());
    // ON-DEMAND patLeaves (§94 Stage 2): class->leaves inverted index + patLeavesOnDemand (hash-probe-smallest
    // intersection + hostFeasible). Built when SCT_ONDEMAND (or the Stage-1 verify). Every patLeaves reader goes
    // through leavesOf(), so flipping ondemand swaps stored O(pattern x leaf) for O(class x leaf) + ~1%-of-peel recompute.
    bool odBuild = ondemand || getenv("SCT_ONDEMAND_VERIFY");
    vector<vector<int>> clsLeaves; vector<int> sumEll, sumU;        // class->leaves index + per-leaf Σell / Σu (hostFeasible O(1))
    if (odBuild) {
        clsLeaves.assign((size_t)nC, {}); sumEll.assign(nLeaf, 0); sumU.assign(nLeaf, 0);
        for (int lid = 0; lid < nLeaf; lid++) {
            for (int c : supC[lid]) clsLeaves[(size_t)c].push_back(lid);
            const CCPath &box = slotPaths[lid][0]; int se = 0, su2 = 0;
            for (int c = 0; c < (int)box.u.size(); c++) { se += (int)box.ell[c]; su2 += (int)box.u[c]; }
            sumEll[lid] = se; sumU[lid] = su2;
        }
    }
    vector<int> odLeaves;
    // hosts of P = leaves with every P-class present that can extend P to an s-clique. Probe the SMALLEST class list;
    // in ONE merge pass confirm presence AND accumulate the hostFeasible delta extra=Σ max(0, P_c - ell_c). Then
    // hostFeasible == (sumEll[lid] + extra <= T <= sumU[lid]). No per-candidate alloc / no separate O(|leaf|) scan.
    auto patLeavesOnDemand = [&](int pi, vector<int> &out) {
        out.clear(); const auto &comp = pats[pi].comp; if (comp.empty()) return;
        int cstar = comp[0].first; size_t best = clsLeaves[(size_t)comp[0].first].size();   // rarest class = smallest list
        for (auto &cm : comp) { size_t z = clsLeaves[(size_t)cm.first].size(); if (z < best) { best = z; cstar = cm.first; } }
        int nc = (int)comp.size();
        for (int lid : clsLeaves[(size_t)cstar]) {
            const vector<int> &sc = supC[lid]; const CCPath &box = slotPaths[lid][0];
            size_t i = 0, j = 0; int matched = 0; long extra = 0;
            while (j < (size_t)nc && i < sc.size()) {
                if (sc[i] < comp[j].first) i++;
                else if (sc[i] > comp[j].first) j++;               // P-class comp[j] not in sc -> matched stays short
                else { if (comp[j].second > (int)box.u[i]) { matched = -1; break; }   // P_c exceeds leaf capacity -> not host
                       int d = comp[j].second - (int)box.ell[i]; if (d > 0) extra += d; matched++; i++; j++; }
            }
            if (matched != nc) continue;                           // some P-class absent / over-capacity -> not a host
            long sl = (long)sumEll[lid] + extra;                   // Σ max(ell, P) = sumEll + Σ max(0, P_c - ell_c)
            if (sl <= box.T && box.T <= sumU[lid]) out.push_back(lid);
        }
    };
    auto leavesOf = [&](int pi) -> const vector<int> & {           // on-demand recompute, else the stored list
        if (ondemand) { patLeavesOnDemand(pi, odLeaves); return odLeaves; }
        return patLeaves[pi];
    };
    (void)patLeavesOnDemand; (void)leavesOf;
    // §94 Stage 3b: affected-Q lookup WITHOUT the per-leaf footprint maps. A leaf-local composition `loc` maps to the
    // global pattern by rebuilding its global comp (leaf class supC[lid][c] with mult loc[c]) and probing the global
    // pattern hash htab (already alive, ~hcap ints). Replaces leafQ2pat + spanEqFP + leafPats -> drops leafFlat (the
    // 2.5GB footprint store) and leafPats entirely. Returns the global pattern id, or -1 if no such pattern.
    vector<pair<int,int>> glComp;
    auto globalLookup = [&](int lid, const Vec &loc, int Mloc) -> int {
        glComp.clear(); const vector<int> &sc = supC[lid];
        for (int c = 0; c < Mloc; c++) if (loc[c]) glComp.push_back({sc[c], (int)loc[c]});   // sorted by global class (sc sorted)
        uint64_t h = compHash(glComp);
        for (size_t idx = h & hmask; htab[idx] != -1; idx = (idx + 1) & hmask)
            if (pats[htab[idx]].comp == glComp) return htab[idx];
        return -1;
    };
    (void)globalLookup;
    // support(pi) = sum over hosting slots of sum over slot's paths of
    // support_count(path, b_local). Uses the pre-mapped compact b.
    Vec sctScr;                                        // reused recompute scratch for sctSupport
    auto sctSupport = [&](int pi) -> double {
        double tot = 0.0;
        const auto &ls = leavesOf(pi);
        for (size_t k = 0; k < ls.size(); k++) {
            const Vec &b = recomputePB ? (localB(pi, ls[k], sctScr), (const Vec &)sctScr)
                                       : pbLocal[pi][k];
            for (auto &p : slotPaths[ls[k]]) tot += ccpath::support_count(p, b, ccpath_ncr);
        }
        return tot;
    };
    auto T5 = Clock::now();
    printf("[sct] pattern<->leaf maps + compaction=%.2fs\n", secs(Tqg1, T5));
    if (getenv("MAPS_MEM_DBG")) {     // analytical map-payload sizes (RSS itself needs Linux); PB drops pbLocal
        long long pbB = 0, pbInc = 0;                          // pbLocal still vector<vector<Vec>> (per-footprint)
        for (auto &v : pbLocal) { pbB += 24; for (auto &fp : v) { pbB += (long long)fp.capacity() * 2 + 40; pbInc++; } }
        long long lpB = 0, lpInc = 0;                          // leafFlat is CSR: per-leaf flat int16 (no per-footprint overhead)
        for (int lid = 0; lid < nLeaf; lid++) { lpB += 24 + (long long)leafFlat[lid].capacity() * 2;
            int Mw = (int)supC[lid].size(); lpInc += Mw ? (long long)leafFlat[lid].size() / Mw : 0; }
        fprintf(stderr, "[maps-mem] pbLocal=%.1fMB(%lld inc, stored=%s) leafFlat(CSR)=%.1fMB(%lld inc) -> PB drops pbLocal, CSR drops per-Vec overhead\n",
                pbB / 1e6, pbInc, recomputePB ? "no" : "yes", lpB / 1e6, lpInc);
    }
    memCk("after-maps(patLeaves/pbLocal)");
    fflush(stdout);
    // COMP_DBG: count feasible s-COMPOSITIONS per leaf (integer pts Sum y=T, 0<=y<=u). This is
    // the peak storage of the user's DIRECT-COUNTING peel (one alive-count a_Y per composition,
    // zeroed on peel, no convolution). Viable iff #compositions stays small vs #patterns.
    if (getenv("COMP_DBG")) {
        long long totComp = 0, maxComp = 0; long long leavesCounted = 0;
        vector<long long> cdp, ndp;
        for (int lid = 0; lid < nLeaf; lid++) {
            long long lc = 0;
            for (auto &bx : slotPaths[lid]) {
                int Mb = bx.m(), Tb = bx.T;
                cdp.assign((size_t)Tb + 1, 0); cdp[0] = 1;
                for (int c = 0; c < Mb; c++) { ndp.assign((size_t)Tb + 1, 0); int uc = (int)bx.u[c];
                    for (int t = 0; t <= Tb; t++) { if (!cdp[t]) continue;
                        int my = uc; if (Tb - t < my) my = Tb - t;
                        for (int y = 0; y <= my; y++) ndp[t + y] += cdp[t]; }
                    cdp.swap(ndp); }
                lc += cdp[(size_t)Tb];
            }
            if (lc > 0) { totComp += lc; if (lc > maxComp) maxComp = lc; leavesCounted++; }
        }
        fprintf(stderr, "[comp] leaves=%lld  total s-compositions=%lld  max/leaf=%lld  avg/leaf=%.1f  vs #patterns=%zu (ratio comp/pat=%.2f)\n",
                leavesCounted, totComp, maxComp, leavesCounted ? (double)totComp / leavesCounted : 0.0, pats.size(),
                pats.size() ? (double)totComp / pats.size() : 0.0);
    }
    // CLS_LEAF_DBG (§93): size the "compute maps on-demand via class->leaves intersection" plan (user's idea).
    // (1) class-leaf incidence (the NEW small inverted-index store) vs pattern-leaf incidence (the maps we DROP).
    // (2) list-length max + per-pattern intersect cost: Σ min-list (hash-probe-smallest) / Σ sum-list (two-pointer),
    //     vs patLeaves-iter (the result size the peel walks anyway) -> the time overhead of computing vs storing.
    if (getenv("CLS_LEAF_DBG")) {
        vector<long long> clsCnt((size_t)nC + 1, 0);            // clsCnt[c] = #leaves containing class c
        long long clsLeafInc = 0;
        for (int lid = 0; lid < nLeaf; lid++) for (int c : supC[lid]) { clsCnt[c]++; clsLeafInc++; }
        long long patLeafInc = 0; for (auto &v : patLeaves) patLeafInc += (long long)v.size();
        long long maxList = 0, nClsUsed = 0;
        for (int c = 0; c <= nC; c++) if (clsCnt[c]) { nClsUsed++; if (clsCnt[c] > maxList) maxList = clsCnt[c]; }
        long long totMin = 0, totSum = 0;
        for (auto &P : pats) {
            long long mn = LLONG_MAX, sm = 0;
            for (auto &cm : P.comp) { long long l = clsCnt[cm.first]; if (l < mn) mn = l; sm += l; }
            totMin += (mn == LLONG_MAX ? 0 : mn); totSum += sm;
        }
        fprintf(stderr, "[cls-leaf] classes=%lld  class-leaf-inc=%lld (NEW store)  pattern-leaf-inc=%lld (maps DROPPED)  maps/newstore=%.1fx  maxlist=%lld avglist=%.1f\n",
                nClsUsed, clsLeafInc, patLeafInc, clsLeafInc ? (double)patLeafInc / clsLeafInc : 0.0,
                maxList, nClsUsed ? (double)clsLeafInc / nClsUsed : 0.0);
        fprintf(stderr, "[cls-leaf] intersect cost: Sum min-list(hash-probe)=%lld  Sum sum-list(two-ptr)=%lld  vs patLeaves-iter=%lld  -> overhead min=%.2fx sum=%.2fx\n",
                totMin, totSum, patLeafInc, patLeafInc ? (double)totMin / patLeafInc : 0.0, patLeafInc ? (double)totSum / patLeafInc : 0.0);
    }
    // ON-DEMAND MAPS, STAGE 1 (§94): prove patLeaves can be COMPUTED (not stored) via class->leaves intersection.
    // patLeavesOnDemand(P) = hash-probe-smallest over clsLeaves[c] (c in P), keep leaf iff every P-class is present
    // AND hostFeasible (P extends to an s-clique in the box). Asserted == the stored patLeaves on every pattern. No
    // behaviour change (gated by SCT_ONDEMAND_VERIFY). The prerequisite for Stage 2 (drop the stored patLeaves).
    if (getenv("SCT_ONDEMAND_VERIFY") && !ondemand) {              // reuses the function-scope patLeavesOnDemand/clsLeaves
        vector<int> od; long long mism = 0;
        for (int pi = 0; pi < (int)pats.size(); pi++) {
            patLeavesOnDemand(pi, od);
            if (od != patLeaves[pi]) { if (mism < 8) fprintf(stderr, "[ondemand] MISMATCH pi=%d stored=%zu od=%zu\n", pi, patLeaves[pi].size(), od.size()); mism++; }
        }
        fprintf(stderr, "[ondemand] verify: %lld/%zu patterns patLeavesOnDemand==stored %s\n",
                (long long)pats.size() - mism, pats.size(), mism ? "[FAIL]" : "[OK]");
    }

    // -------- support init: SCT (production) + optional region-IE cross-check (gate G2a) -------
    // P.sup := SCT sum-over-leaves support. Under SCT_VERIFY also compare to the region-IE init
    // (suppOf, set above) per pattern and abort on mismatch. sctSupport is computed either way
    // (it IS the production support), so production just skips the 99s region-IE entirely.
    {
        int okc = 0, badc = 0; double worst = 0;
        for (int pi = 0; pi < (int)pats.size(); pi++) {
            double sSCT = sctSupport(pi);
            if (verifyIE) {
                double sIE = pats[pi].sup;             // region-IE init
                if (fabs(sIE - sSCT) < 0.5) okc++;
                else {
                    badc++; worst = max(worst, fabs(sIE - sSCT));
                    if (badc <= 8) {
                        printf("[G2a] MISMATCH pi=%d host=%d comp=[", pi, pats[pi].hostSz);
                        for (auto &cm : pats[pi].comp) printf("(c%d:%d)", cm.first, cm.second);
                        printf("] regionIE=%.1f SCT=%.1f leaves=%zu\n", sIE, sSCT, patLeaves[pi].size());
                    }
                }
            }
            pats[pi].sup = sSCT;                        // PRODUCTION: the SCT IS the support
        }
        if (verifyIE) {
            printf("[G2a] %d/%d patterns: SCT-sum == region-IE  %s  (worst|d|=%.1f)\n",
                   okc, okc + badc, badc == 0 ? "[OK]" : "[FAIL]", worst);
            fflush(stdout);
            if (badc != 0) { printf("[G2a] FAILED — aborting before peel (correctness gate).\n"); return 3; }
        }
    }

    // ---- bucket-queue peel on the class-SCT ----
    // Peeling P: for each slot hosting P, insert tuple_to_threshold(slot,m_P)
    // into every path's forbidden, split any path whose antichain exceeds
    // KMAX. Affected patterns = those sharing a slot with P; recompute their
    // support from current slot state. Monotone: clamp new key >= curLevel.
    // KMAX bounds each path's forbidden antichain (controlled_split keeps it
    // <= KMAX), so support_count's inclusion-exclusion stays <= 2^KMAX terms.
    // With the changed-paths incremental update the split COUNT no longer
    // hurts, so a small KMAX (cheap per-call IE) wins. KMAX=1 minimises the
    // per-candidate IE/DP (<=2 terms) and measured FASTEST on every cell tested
    // (moderate ca-GrQc/ca-CondMat/dblp ~17-24% faster than KMAX=2, and even on
    // the maxSplit-stress com-dblp(3,4) 10.7s vs 12.9s) -- the affected-Q drop is
    // per-candidate-DP-bound, so the slightly larger split-set (cheap impossible-
    // skip scan) never dominates. Correctness is invariant to KMAX (verified
    // 80/80 vs brute at KMAX=1,2,12; corehash byte-identical across KMAX here).
    int KMAX = 1;
    bool kAdapt = true;                               // adaptive per-leaf KMAX (default ON)
    if (getenv("SCT_KMAX")) { KMAX = atoi(getenv("SCT_KMAX")); if (KMAX < 1) KMAX = 1; kAdapt = false; }
    // ADAPTIVE KMAX: KMAX=1 minimises per-candidate IE/DP and measured FASTEST on
    // EVERY graph tested -- moderate (~1.5-2x vs KMAX=2) AND maxSplit-blowup
    // (web-NotreDame(6,8) maxSplit~14k: 41.2s vs 53.5s for KMAX=2 on a quiet host;
    // an earlier "regression" was server contention, not real). So KMAX=1 is the
    // default everywhere observed (slot up to ~14k). The O(slot) scan in
    // slotForbidDiff is cheaply 99%-skipped (the impossible test), so it does NOT
    // overtake the DP in the measured range. The adaptive rule only HEDGES the
    // unmeasured extreme tail (slot > KTHRESH=16384, e.g. soc-pokec / web-it-2004 at
    // high (r,s)): there a leaf's KMAX rises (kml = 1 + slotSize/KTHRESH) so a runaway
    // slot self-limits. Bit-identical: support_count is invariant to the split
    // strategy (KMAX-invariance), so a per-leaf slot-dependent KMAX is exact.
    int KTHRESH = 16384; if (getenv("SCT_KMAX_THRESH")) KTHRESH = atoi(getenv("SCT_KMAX_THRESH"));
    const int KMAXCAP = 6;                            // ceiling on the adaptive KMAX
    // WITNESS-FLOOR fast path for s=r+1: when s=r+1 the dying witnesses above a
    // peeled pattern P are the single points floor_c = m_P + e_c (Σ=s), so the
    // affected Q are exactly floor_c - e_d (move one unit), grouped under ONE
    // liveness check per (c,box). The drop has a closed form (no IE/DP):
    //   drop_Q(box p) = [floor_c alive in p] * (n_p[d] - m_Q[d]).
    // This is output-sensitive (work ~ live floors, not all candidates) and was
    // validated bit-exact vs scWithTerms on 400k random boxes (incl ell>0).
    // UNIFIED witness-major affected-update (§72-75): one tail-parameterized path for every small tail t=s-r.
    // t=1 == old witness-floor, t=2 == old witness-major. The witness enumeration is O(M^t) (M=leaf width), so
    // whether it beats the general DFS is GRAPH- and LEAF-dependent (sec 75). Instead of a fixed cap we choose
    // PER LEAF: witness while its δ-enumeration is cheap relative to the general DFS's candidate space, else
    // fall back to the general DFS (both bit-identical -> pure speed). t=1,2 always win (measured across the
    // density spectrum) so they skip the gate; the gate kicks in at t>=witGateMinT. witnessTMax is a hard cap.
    int witnessTail = s - r;
    int witnessTMax = getenv("SCT_WITNESS_TMAX") ? atoi(getenv("SCT_WITNESS_TMAX")) : 8;
    int witGateMinT = getenv("SCT_WIT_GATE_MINT") ? atoi(getenv("SCT_WIT_GATE_MINT")) : 3;
    double witK = getenv("SCT_WIT_K") ? atof(getenv("SCT_WIT_K")) : 8.0;
    bool witnessActive = witnessTail >= 1 && witnessTail <= witnessTMax;
    // METHOD ROUTING by the measured crossover (sec 75-77): witness wins for the common small tails (t<=witCross),
    // the wave BATCH driver wins past it (t>witCross: batch 15.8s vs witness 79s at dblp-db t=4). t>witCross routes
    // to batch BY DEFAULT (both bit-identical). Tunable SCT_WIT_CROSS. (Per-leaf witness/general gate handles t<=cross.)
    int witCross = getenv("SCT_WIT_CROSS") ? atoi(getenv("SCT_WIT_CROSS")) : 3;
    // SKIP_H1: a |host|=1 pattern peels at EXACTLY L_M=C(|M|-r,s-r) regardless
    // of how the peel proceeds (every r-clique in its region M has support
    // >= L_M, so no witness of a |host|=1 pattern dies before curLevel=L_M --
    // verified 0 exceptions on all test graphs). Its initial key already ==
    // L_M, so recomputing its support during affected-updates is pure waste:
    // the answer never changes its bucket. We therefore SKIP |host|=1 patterns
    // as affected-update TARGETS (they still peel normally as SOURCES, removing
    // their witnesses so |host|>=2 patterns drop correctly). Since |host|=1 is
    // the majority on sparse/real graphs, this removes the bulk of the
    // per-affected scWithTerms work. Bit-identical to the full peel.
    bool skipH1 = getenv("SCT_NO_SKIP_H1") == nullptr;   // default ON
    // DFS_PRUNE: per-path AND-feasibility prune of the affected-Q DFS. The DFS's
    // cap=min(uEnv[c],rem) uses uEnv = per-coord UNION (max_z u_z) over chgOld, so
    // it admits ql that fit the union but NO single path; those score d==0 in every
    // chgOld path (applyIdx's max(ql,pl)<=p.u test) and are pure waste (~57 cands
    // generated for ~4 with a real drop). We track, per chgOld path z, whether the
    // prefix max(pl,ql) still fits u_z; once NO path is alive we prune the subtree.
    // Killing is monotone in descent (running max only grows, Σ only grows), and
    // pruned candidates are exactly applyIdx's d==0 set, so cores are BIT-IDENTICAL
    // (applyIdx stays the final oracle: an over-keep just wastes a lookup). All
    // chgOld paths share T=s, so the Σ<=T bound stays path-independent (sfp/Tcap).
    bool dfsPrune = getenv("SCT_DFS_PRUNE") != nullptr;  // default OFF (A/B flag)
    size_t maxSplit = 0;                              // diagnostic: largest split-set
    double tSFD = 0; long long slotVisits = 0;       // PROFILE: slot-scan time + path visits
    bool sfdDbg = getenv("SFD_DBG") != nullptr;       // cost-structure probe for the slot-index design
    long long sfdAff = 0, sfdCoordTests = 0, sfdFailFirst = 0;  // affected / coords-examined / failed-on-1st-coord
    // SUPP_DBG (§68 gate): histogram of |supp(a)| over forbidden-threshold insertions. |supp|==1 == axis-aligned
    // == a clean u_c reduction (the hybrid "delete a class" fast path, bit-identical). Track calls/affected-path
    // work/controlled_split triggers by |supp| -> single-class SHARE = the upper bound on the hybrid's win.
    bool suppDbg = getenv("SUPP_DBG") != nullptr;
    long long suppCalls[12] = {0}, suppAff[12] = {0}, suppSplit[12] = {0};
    long long idxMismatch = 0;                         // slot-index verify: #calls where index-find != scan-find
    // Record a pattern's LOCAL threshold into slot lid. Paths where the
    // threshold is impossible are UNCHANGED (kept in place); the rest are the
    // CHANGED paths. We collect the changed OLD paths (chgOld) and their NEW
    // replacements (chgNew) so the caller can compute affected-pattern support
    // deltas over only the changed paths (unchanged paths cancel in before-
    // after). This makes the cost independent of the total split-set size.
    // The delta-formula affected-update reads ONLY chgOld (the pre-insertion
    // snapshots); the split children cancel algebraically, so we no longer
    // materialize them into a chgNew list -- they go straight into `keep` (the
    // updated slot) by MOVE, halving the path copies in this hot routine.
    //
    // A slot can hold thousands of split paths, but P's threshold a (= bloc,
    // since every leaf path has h=0) fits in almost none of them: profiling
    // showed ~99% of path visits hit `impossible` and are skipped. The
    // impossible test only depends on P's NONZERO classes (a[c]=0 elsewhere is
    // always <= u[c]), so we pass the sparse support of bloc (plNZ = list of
    // (localpos,val)) and test just those few positions -- avoiding both the
    // per-path tuple_to_threshold allocation and the full-width O(M) scan over
    // the ~99% of paths that are unaffected.
    vector<CCPath> sfdKids;                            // reused scratch (split children)
    // ===== SLOT DOMINANCE INDEX (peel #1: find affected paths sub-linearly) =====
    // The impossible scan visits all w paths to find the ~0.2-1% AFFECTED ones (those
    // with u[c] >= bloc[c] for every plNZ coord c). Index per leaf: bkt[c][v] = slot
    // positions with u[c]==v. A path is affected iff u dominates bloc on plNZ, so the
    // candidates are bkt[pivot][v>=bloc[pivot]] for the most-selective plNZ coord; we
    // then filter survivors by the other plNZ coords. u only DECREASES on split, so
    // maxv from the initial build bounds all future values. Slot order is NOT load-
    // bearing (sec 54, SCT_SLOT_REVERSE 12/0), so find-then-act is safe. Maintained
    // via swap-remove + back-pointers bpos (cur stays dense -> consumers unchanged).
    bool slotIdx       = getenv("SCT_NO_SLOT_IDX") == nullptr;       // DEFAULT ON (escape hatch: SCT_NO_SLOT_IDX)
    bool slotIdxVerify = getenv("SCT_SLOT_IDX_VERIFY") != nullptr;   // Step 1: assert index==scan
    if (slotIdxVerify) slotIdx = false;                             // verify mode needs the full scan to drive
    bool idxOn = slotIdx || slotIdxVerify;
    bool idxDbg = getenv("IDX_DBG") != nullptr;        // localize index cost: pivot-scan vs candidate-filter vs output
    long long ixPivScan = 0, ixCand = 0, ixOut = 0;    // Σ(mv-thr) over coords / Σ pivot-candidates filtered / Σ affected
    // ===== a_Y DIRECT DEAD-SET (§87/§88): explicit per-composition alive-flag, replaces the antichain split churn =====
    // The forbidden antichain answers "is witness-composition Y alive?" via IE over the peeled-pattern generators,
    // and keeping that IE cheap forces controlled_split (measured 52% / 5.85M splits on ca-AstroPh 3,4). a_Y stores the
    // dead set EXPLICITLY: Y is dead iff Y dominates some peeled pattern, i.e. iff some addDelta enumeration already
    // marked it. Query/mark are O(1) (a hash-set insert) -- NO antichain, NO split, NO IE. Bit-identical: same dead
    // set, same witness-major drop. Witness path only (single-pattern regime); skips slotForbidDiff for the leaf.
    // DEFAULT-ON for tail t=1 (s=r+1): a_Y is faster AND leaner than the antichain on EVERY t=1 graph tested
    // (ca-GrQc/com-dblp/com-youtube/ca-AstroPh; 1.0-1.8x peel, -5..-12% RSS), bit-identical. It also turns the
    // ca-AstroPh 4,5/5,6 antichain TIMEOUTS into finishing runs. t>=2 keeps the antichain by DEFAULT (a_Y
    // over-enumerates the full δ-space on sparse t>=2; the dense-t>=2 adaptive gate is future work). Escapes:
    // SCT_AY forces a_Y for all witness tails (t<=witnessTMax); SCT_NO_AY restores the antichain everywhere.
    bool ayMode = (getenv("SCT_NO_AY") == nullptr && witnessTail == 1) || getenv("SCT_AY") != nullptr;
    // Flat open-addressing uint64 set (linear probe, power-of-2): ~5-10x faster per op than std::unordered_set
    // (no node alloc, cache-friendly) -- the per-Y dead-check/mark constant is what decides a_Y on sparse t>=2.
    struct FlatU64 {
        std::vector<uint64_t> t; size_t mask = 0, cnt = 0;
        inline bool insert(uint64_t k) {                       // true iff newly inserted (was alive)
            if (!k) k = 0x9E3779B97F4A7C15ULL;                 // remap reserved empty slot (0)
            if (t.empty()) { t.assign(16, 0); mask = 15; }
            else if ((cnt + 1) * 4 >= t.size() * 3) grow();    // load factor 0.75
            size_t i = (k * 0x9E3779B97F4A7C15ULL >> 28) & mask;
            while (t[i]) { if (t[i] == k) return false; i = (i + 1) & mask; }
            t[i] = k; ++cnt; return true;
        }
        void grow() {
            std::vector<uint64_t> old; old.swap(t);
            t.assign(old.size() * 2, 0); mask = t.size() - 1; cnt = 0;
            for (uint64_t k : old) if (k) { size_t i = (k * 0x9E3779B97F4A7C15ULL >> 28) & mask;
                while (t[i]) i = (i + 1) & mask; t[i] = k; ++cnt; }
        }
    };
    vector<FlatU64> deadY(ayMode ? nLeaf : 0);                 // per-leaf dead witness-composition hashes
    vector<char> ixBuilt(idxOn ? nLeaf : 0, 0);
    vector<int>  ixM(idxOn ? nLeaf : 0, 0), ixMaxv(idxOn ? nLeaf : 0, 0);
    vector<vector<vector<int>>> ixBkt(idxOn ? nLeaf : 0);  // [lid][c*maxv+v] -> positions
    vector<vector<int>>         ixBpos(idxOn ? nLeaf : 0); // [lid][pos*M+c] -> index in its bucket
    auto ixBuild = [&](int lid) {
        auto &cur = slotPaths[lid]; int M = (int)supC[lid].size();
        int mv = 1; for (auto &p : cur) for (int c = 0; c < M; c++) { int u = (int)p.u[c]; if (u + 1 > mv) mv = u + 1; }
        ixM[lid] = M; ixMaxv[lid] = mv;
        auto &bkt = ixBkt[lid]; bkt.assign((size_t)M * mv, {});
        auto &bp = ixBpos[lid]; bp.assign(cur.size() * (size_t)M, 0);
        for (int i = 0; i < (int)cur.size(); i++) for (int c = 0; c < M; c++) {
            int key = c * mv + (int)cur[i].u[c];
            bp[(size_t)i * M + c] = (int)bkt[key].size(); bkt[key].push_back(i);
        }
        ixBuilt[lid] = 1;
    };
    auto ixRemove = [&](int lid, int i, const Vec &uv) {            // drop position i (uv = its u)
        int M = ixM[lid], mv = ixMaxv[lid]; auto &bkt = ixBkt[lid]; auto &bp = ixBpos[lid];
        for (int c = 0; c < M; c++) { auto &B = bkt[c * mv + (int)uv[c]];
            int idx = bp[(size_t)i * M + c], last = B.back(); B[idx] = last; bp[(size_t)last * M + c] = idx; B.pop_back(); }
    };
    auto ixRelabel = [&](int lid, int from, int to, const Vec &uv) {// path uv moves from->to
        int M = ixM[lid], mv = ixMaxv[lid]; auto &bkt = ixBkt[lid]; auto &bp = ixBpos[lid];
        for (int c = 0; c < M; c++) { int idx = bp[(size_t)from * M + c]; bkt[c * mv + (int)uv[c]][idx] = to; bp[(size_t)to * M + c] = idx; }
    };
    auto ixAppend = [&](int lid, int p, const Vec &uv) {            // new path at position p
        int M = ixM[lid], mv = ixMaxv[lid]; auto &bkt = ixBkt[lid]; auto &bp = ixBpos[lid];
        if ((int)bp.size() < (p + 1) * M) bp.resize((size_t)(p + 1) * M, 0);
        for (int c = 0; c < M; c++) { int key = c * mv + (int)uv[c]; bp[(size_t)p * M + c] = (int)bkt[key].size(); bkt[key].push_back(p); }
    };
    auto ixFindAffected = [&](int lid, const vector<pair<int,int>> &plNZ, vector<int> &out) {
        out.clear(); int M = ixM[lid], mv = ixMaxv[lid]; auto &bkt = ixBkt[lid]; auto &cur = slotPaths[lid];
        int piv = -1, pivThr = 0; long best = LONG_MAX;            // pivot = min-count plNZ coord
        for (auto &pv : plNZ) { if (pv.first >= M) { out.assign(1, -1); return; }   // bloc coord beyond leaf width => no path (defensive)
            long cc = 0; for (int v = pv.second; v < mv; v++) cc += (long)bkt[pv.first * mv + v].size();
            if (idxDbg) ixPivScan += (mv - pv.second);
            if (cc < best) { best = cc; piv = pv.first; pivThr = pv.second; } }
        if (piv < 0) return;
        for (int v = pivThr; v < mv; v++) for (int pos : bkt[piv * mv + v]) {
            if (idxDbg) ixCand++;
            bool ok = true; for (auto &pv : plNZ) { if (pv.first == piv) continue;
                if ((int)cur[pos].u[pv.first] < pv.second) { ok = false; break; } }
            if (ok) out.push_back(pos);
        }
        if (idxDbg) ixOut += (long long)out.size();
    };
    vector<int> ixAff, scanAff;                        // reused: index-found / scan-found affected positions
    // IN-PLACE: unchanged paths (~99%) stay put (never moved). Only changed paths
    // are snapshotted; covers-whole / split-parent paths are swap-removed; split
    // children are appended. Same resulting slot (order-independent: leaves are a
    // disjoint set, support = sum); same chgOld set. Avoids the per-call rebuild of
    // the entire split-set (was 339M CCPath moves on com-dblp(3,4)).
    auto slotForbidDiff = [&](int lid, const Vec &bloc,
                              const vector<pair<int,int>> &plNZ, vector<CCPath> &chgOld) {
        vector<CCPath> &cur = slotPaths[lid];
        sfdKids.clear();
        chgOld.clear();
        int sdK = suppDbg ? ((int)plNZ.size() < 11 ? (int)plNZ.size() : 11) : 0;  // |supp(a)| bucket
        long long sdAff0 = sfdAff; int sdSplit = 0;
        if (suppDbg) suppCalls[sdK]++;
        if (idxOn && !ixBuilt[lid]) ixBuild(lid);
        if (slotIdxVerify) {                           // Step-1 check: index-find == scan-find (pristine state)
            ixFindAffected(lid, plNZ, ixAff);
            scanAff.clear();
            for (int i = 0; i < (int)cur.size(); i++) { bool im = false;
                for (auto &pv : plNZ) if ((int)cur[i].u[pv.first] < pv.second) { im = true; break; }
                if (!im) scanAff.push_back(i); }
            std::sort(ixAff.begin(), ixAff.end());
            if (ixAff != scanAff) {
                fprintf(stderr, "[slot-idx] MISMATCH lid=%d scan=%zu idx=%zu\n", lid, scanAff.size(), ixAff.size());
                idxMismatch++;
            }
        }
        int w = (int)cur.size();                          // live prefix [0,w)
        // effective KMAX for THIS leaf: raised as the slot grows so a blow-up slot
        // self-limits (more forbidden per path => fewer split children). Computed from
        // the entry slot size (stable within this call). Bit-identical (KMAX-invariant).
        int kml = KMAX;
        if (kAdapt) { kml = KMAX + w / KTHRESH; if (kml > KMAXCAP) kml = KMAXCAP; }
        if (slotIdx) {                                    // Step 2: INDEX-DRIVEN (skip the full scan)
            ixFindAffected(lid, plNZ, ixAff);             // affected positions (pristine slot)
            scanAff.clear();                              // reuse as deadPos
            for (int pos : ixAff) {                       // Phase A: classify (no position moves)
                sfdAff++;
                chgOld.push_back(cur[pos]);               // pre-change snapshot
                if (ccpath::covers_whole_path(cur[pos], bloc)) { scanAff.push_back(pos); continue; }
                ccpath::insert_antichain(cur[pos].forbidden, bloc);
                if ((int)cur[pos].forbidden.size() > kml) {
                    if (suppDbg) sdSplit++;
                    auto kk = ccpath::controlled_split(cur[pos], kml);
                    for (auto &k : kk) sfdKids.push_back(std::move(k));
                    scanAff.push_back(pos);               // split-parent dies (u unchanged until removed)
                }
                // else modified in place: u unchanged -> index entry stays valid.
            }
            std::sort(scanAff.begin(), scanAff.end());    // Phase B: batch swap-remove dead, descending
            for (int k = (int)scanAff.size() - 1; k >= 0; --k) {
                int p = scanAff[k]; --w;                  // (w >= p always; cur[w] is live)
                ixRemove(lid, p, cur[p].u);
                if (p != w) { cur[p] = std::move(cur[w]); ixRelabel(lid, w, p, cur[p].u); }
            }
            cur.resize(w);
            cur.reserve((size_t)w + sfdKids.size());
            for (auto &k : sfdKids) { ixAppend(lid, (int)cur.size(), k.u); cur.push_back(std::move(k)); }
            if (cur.size() > maxSplit) maxSplit = cur.size();
            if (suppDbg) { suppAff[sdK] += sfdAff - sdAff0; suppSplit[sdK] += sdSplit; }
            return;
        }
        for (int i = 0; i < w; ) {
            CCPath &p = cur[i];
            bool imposs = false;                          // impossible(p, bloc)?
            if (sfdDbg) {                                 // probe: count coord-tests + first-coord failures
                int nt = 0; for (auto &pv : plNZ) { nt++; if ((int)p.u[pv.first] < pv.second) { imposs = true; break; } }
                sfdCoordTests += nt; if (imposs && nt == 1) sfdFailFirst++;
            } else
                for (auto &pv : plNZ) if ((int)p.u[pv.first] < pv.second) { imposs = true; break; }
            if (imposs) { ++i; continue; }                // unchanged: stays in place
            sfdAff++;
            chgOld.push_back(p);                            // snapshot before change
            bool remove = false;
            if (ccpath::covers_whole_path(p, bloc)) {
                remove = true;                             // path fully dead (a==bloc)
            } else {
                ccpath::insert_antichain(p.forbidden, bloc);
                if ((int)p.forbidden.size() > kml) {
                    if (suppDbg) sdSplit++;
                    auto kk = ccpath::controlled_split(p, kml);
                    for (auto &k : kk) sfdKids.push_back(std::move(k));
                    remove = true;                         // split-parent replaced by children
                }
            }
            if (remove) {                                  // swap-remove (no self-move)
                --w;
                if (idxOn) ixRemove(lid, i, cur[i].u);     // drop i's index entries (its own u, pre-overwrite)
                if (i != w) { cur[i] = std::move(cur[w]); if (idxOn) ixRelabel(lid, w, i, cur[i].u); }
            }
            else ++i;                                      // modified in place (u unchanged -> index unchanged)
        }
        cur.resize(w);
        cur.reserve((size_t)w + sfdKids.size());          // grow once, not per child
        for (auto &k : sfdKids) { if (idxOn) ixAppend(lid, (int)cur.size(), k.u); cur.push_back(std::move(k)); }
        if (cur.size() > maxSplit) maxSplit = cur.size();
        if (suppDbg) { suppAff[sdK] += sfdAff - sdAff0; suppSplit[sdK] += sdSplit; }
    };

    // ---- fast support_count with PRECOMPUTED inclusion-exclusion terms ----
    // support_count_impl recomputes inclusion_exclusion_terms(p.forbidden) on
    // every call, even though the terms depend only on the path's forbidden
    // antichain (not on b). In the affected-update loop the SAME path is queried
    // against MANY patterns Q, so we hoist the IE-term build out of the Q loop:
    // precompute terms once per changed path, then evaluate every Q against the
    // cached terms. This is bit-identical to ccpath::support_count (same early-
    // out, same terms, same count_with_extra_lower) but removes the dominant
    // redundant IE rebuild. Reused DP scratch avoids a malloc/free per term.
    vector<double> dpA, dpB;                          // reused DP scratch
    // count_with_extra_lower (b = weight base / lower) plus an INDEPENDENT extra
    // lower `addLow` (raises the per-class lower WITHOUT changing the C(n-b,y-b)
    // weight base). addLow==nullptr reproduces ccpath::support_count exactly.
    // With addLow set, this computes the alive-witness count >= max(b, addLow)
    // still weighted by C(n-b, y-b) -- exactly the per-path peel DROP for a
    // query b=m_Q with addLow=a_p (P's threshold). The weight base staying at b
    // is what makes this correct (a plain support_count(max(b,addLow)) would
    // change the weight; see lines below).
    // CF_DBG (closed-form viability probe): per DP IE-term, count how many upper bounds
    // u_c actually TRUNCATE (Uc-Lc < total slack). The closed form replaces the DP with
    // 2^(#binding) Vandermonde binomials per term -- viable iff #binding stays small.
    bool cfDbg = getenv("CF_DBG") != nullptr;
    long long cfTerms = 0, cfBindHist[20] = {0}; double cfFormTermsVsDP = 0;
    // HYBRID closed form (sec 65): per IE-term, if the active classes give few expansion
    // terms, compute the bounded binomial-weighted composition count by Vandermonde IE
    // (pure binomials, NO convolution) instead of the DP. Bit-identical (exact integer).
    bool cfEnabled = getenv("SCT_NO_CF") == nullptr;   // DEFAULT ON (escape hatch SCT_NO_CF)
    const long long CF_CAP = 4096;                     // skip closed form if #terms exceeds this
    vector<int> cfActM, cfActLo, cfActUc;              // reused: per active class M_c, ell_c, u'_c
    auto scWithTerms = [&](const CCPath &p,
                           const vector<pair<Vec,int>> &terms,
                           const Vec &b, const Vec *addLow) -> double {
        // forbidden early-out is only valid when no addLow raises the floor; with
        // addLow we must let the IE terms account for the forbidden set instead.
        if (!addLow) {
            for (const auto &a : p.forbidden) if (ccpath::leq(a, b)) return 0.0;
        } else {
            // addLow raises the per-class floor to max(b, addLow); if ANY single
            // forbidden a <= that floor, every witness >= floor is >= a (dead) so
            // the count is 0 -- skip the IE+DP. Correct (only fires when truly 0;
            // union-covered zeros still fall through to the IE).
            const int M2 = p.m();
            for (const auto &a : p.forbidden) {
                bool below = true;
                for (int c = 0; c < M2; c++) {
                    int fl = (int)b[c] > (int)(*addLow)[c] ? (int)b[c] : (int)(*addLow)[c];
                    if ((int)a[c] > fl) { below = false; break; }
                }
                if (below) return 0.0;
            }
        }
        const int M = p.m();
        const int T = p.T;
        double total = 0.0;
        for (const auto &kv : terms) {
            const Vec &extra = kv.first;
            // ---- inlined count_with_extra_lower_impl with shared scratch ----
            bool bad = false;
            int sumL = 0, sumU = 0;
            // bounds per class (recomputed; cheap, M is tiny)
            for (int c = 0; c < M; ++c) {
                if (b[c] < 0 || (int)b[c] > (int)p.n[c]) { bad = true; break; }
                int Lc = p.ell[c];
                if ((int)b[c] > Lc) Lc = (int)b[c];
                if ((int)extra[c] > Lc) Lc = (int)extra[c];
                if (addLow && (int)(*addLow)[c] > Lc) Lc = (int)(*addLow)[c];
                int Uc = (int)p.u[c];
                if (Lc > Uc) { bad = true; break; }
                sumL += Lc; sumU += Uc;
            }
            if (bad || sumL > T || sumU < T) continue;
            if (cfDbg) {                                  // count binding upper bounds for this IE term
                int slack = T - sumL, nb = 0;
                for (int c = 0; c < M; ++c) {
                    int Lc = p.ell[c]; if ((int)b[c] > Lc) Lc = (int)b[c];
                    if ((int)extra[c] > Lc) Lc = (int)extra[c];
                    if (addLow && (int)(*addLow)[c] > Lc) Lc = (int)(*addLow)[c];
                    if ((int)p.u[c] - Lc < slack) nb++;
                }
                cfTerms++; cfBindHist[nb < 19 ? nb : 19]++;
                // closed-form work = 2^nb binomials; DP work ~ M*(T+1) cells. ratio:
                cfFormTermsVsDP += (double)(1 << (nb < 20 ? nb : 20)) / ((double)M * (T + 1) + 1);
            }
            if (cfEnabled) {                              // ---- CLOSED-FORM fast path (sec 65) ----
                int Z = T, totalM = 0; cfActM.clear(); cfActLo.clear(); cfActUc.clear();
                long long nterms = 1; bool cfBig = false;
                for (int c = 0; c < M; ++c) {
                    int bc = (int)b[c], Mc = (int)p.n[c] - bc; totalM += Mc; Z -= bc;
                    int Lc = p.ell[c]; if (bc > Lc) Lc = bc;
                    if ((int)extra[c] > Lc) Lc = (int)extra[c];
                    if (addLow && (int)(*addLow)[c] > Lc) Lc = (int)(*addLow)[c];
                    int lc = Lc - bc;                     // ell_c (low-tail length = lc terms: 0..lc-1)
                    int uc = (int)p.u[c] - bc; if (uc > Mc) uc = Mc;   // u'_c (high-tail = Mc-uc terms: uc+1..Mc)
                    int ntail = lc + (Mc - uc);
                    if (ntail > 0) { cfActM.push_back(Mc); cfActLo.push_back(lc); cfActUc.push_back(uc);
                        nterms *= (long long)(1 + ntail); if (nterms > CF_CAP) { cfBig = true; break; } }
                }
                if (!cfBig && nterms <= (long long)M * (T + 1)) {   // closed form cheaper than the DP
                    int na = (int)cfActM.size(); double cf = 0.0;
                    auto rec = [&](auto &&self, int idx, int sign, int sumJ, int sumMpick, double prodW) -> void {
                        if (idx == na) { int rem = Z - sumJ, navail = totalM - sumMpick;
                            if (rem >= 0 && rem <= navail) cf += sign * prodW * ccpath_ncr(navail, rem); return; }
                        int Mc = cfActM[idx], lc = cfActLo[idx], uc = cfActUc[idx];
                        self(self, idx + 1, sign, sumJ, sumMpick, prodW);                 // full (1+x)^Mc
                        for (int j = 0; j < lc; ++j)                                       // low tail z<ell
                            self(self, idx + 1, -sign, sumJ + j, sumMpick + Mc, prodW * ccpath_ncr(Mc, j));
                        for (int j = uc + 1; j <= Mc; ++j)                                 // high tail z>u'
                            self(self, idx + 1, -sign, sumJ + j, sumMpick + Mc, prodW * ccpath_ncr(Mc, j));
                    };
                    rec(rec, 0, 1, 0, 0, 1.0);
                    total += (double)kv.second * cf;
                    continue;                             // skip the DP for this term
                }
            }
            dpA.assign((size_t)T + 1, 0.0);
            dpA[0] = 1.0;
            for (int c = 0; c < M; ++c) {
                dpB.assign((size_t)T + 1, 0.0);
                const int bc = (int)b[c];
                const int nc = (int)p.n[c];
                int Lc = p.ell[c];
                if (bc > Lc) Lc = bc;
                if ((int)extra[c] > Lc) Lc = (int)extra[c];
                if (addLow && (int)(*addLow)[c] > Lc) Lc = (int)(*addLow)[c];
                const int Uc = (int)p.u[c];
                for (int tot = 0; tot <= T; ++tot) {
                    double w = dpA[(size_t)tot];
                    if (w == 0.0) continue;
                    int maxy = Uc;
                    if (T - tot < maxy) maxy = T - tot;
                    for (int y = Lc; y <= maxy; ++y)
                        dpB[(size_t)(tot + y)] += w * ccpath_ncr(nc - bc, y - bc);
                }
                dpA.swap(dpB);
            }
            total += (double)kv.second * dpA[(size_t)T];
        }
        return total;
    };

    // ---- lazy per-leaf lookup: local-witness HASH -> leaf-pattern indices ----
    // The affected-Q set for a peel is enumerated DIRECTLY from m_P (the DFS
    // below) instead of scanning every pattern in the leaf (~99.8% of those scans
    // failed the feasibility filter). Each enumerated candidate m_Q is looked up
    // here. Built lazily on first touch (most leaves are touched many times). Key
    // = a rolling polynomial HASH of the local m_Q (robust to any leaf width,
    // unlike a perfect mixed-radix pack which overflows for wide leaves, e.g.
    // ca-GrQc(4,5) reaches width 28). Bucket value = indices into leafPats[lid];
    // a candidate is confirmed by comparing the full local vector (only on a hash
    // hit, which is rare), so the lookup is exact.
    vector<std::unordered_map<uint64_t,vector<int>>> leafQ2pat(nLeaf);
    vector<char> leafQbuilt(nLeaf, 0);
    const uint64_t HMUL = 1099511628211ULL;          // FNV-style multiplier
    auto hashVec = [&](const Vec &v) -> uint64_t {
        uint64_t h = 1469598103934665603ULL;
        for (size_t i = 0; i < v.size(); i++) { h ^= (uint64_t)(uint16_t)v[i] + 1; h *= HMUL; }
        return h;
    };
    Vec elmScr;                                        // reused recompute scratch for ensureLeafMap
    auto ensureLeafMap = [&](int lid) {
        if (leafQbuilt[lid]) return;
        leafQbuilt[lid] = 1;
        auto &mp = leafQ2pat[lid];
        int nt = (int)leafPats[lid].size();                    // == #footprints (int map always stored)
        int Mw = (int)supC[lid].size();
        mp.reserve((size_t)nt * 2);
        for (int t = 0; t < nt; t++) {
            auto fp = leafFP(lid, t, Mw, elmScr);              // (ptr,len) span: recompute or CSR view
            mp[hashSpan(fp.first, fp.second)].push_back(t);
        }
    };

    long long npat = (long long)pats.size(), peeledN = 0, maxKey = 0;
    for (auto &P : pats) { P.key = (long long)llround(P.sup); maxKey = max(maxKey, P.key); }
    unordered_map<long long, vector<int>> bk;
    for (int pi = 0; pi < (int)pats.size(); pi++) bk[pats[pi].key].push_back(pi);
    map<double,double> coreDist;
    long long curLevel = 0;
    vector<char> seen(pats.size(), 0);
    vector<double> delta(pats.size(), 0.0);          // per-affected exact drop
    Vec uEnv, sufPl, qcand, Yscr;                     // reused per-leaf scratch (Yscr = s=r+2 witness)
    vector<long long> wdDP, cbDP2;                     // adaptive-gate scratch: #witnesses Wδ / #candidates CB
    vector<char> uLiveBuf, covBuf;                     // depth-indexed prune scratch (DFS_PRUNE)
    vector<const int16_t*> chgU, fbA;                  // u-rows / single-forbidden rows of chgOld
    vector<int> fbCrit; vector<char> fbHas;            // per-path: max critical coord / has-1-forbidden
    long long dbgGen = 0, dbgHit = 0, dbgNZ = 0;       // instrumentation: cands gen / hit / nonzero-drop
    long long witInst = 0, witMSum = 0, witMMax = 0;   // witness path: leaf-instances + leaf-width M (drives crossover)
    long long witGateW = 0, witGateG = 0;              // gated leaves (t>=minT): chose witness / fell back to general
    bool witDbg = getenv("SCT_WIT_DBG") != nullptr;
    // batch-peel kill-gate (FANIN_DBG, sec 57): fanin = total (pattern,leaf) affected-update touches /
    // distinct (leaf,level) touches = avg #patterns sharing a (leaf,curLevel). fanin>=5 => batch-peel pays.
    bool faninDbg = getenv("FANIN_DBG") != nullptr;
    long long fanA = 0, fanB = 0;
    vector<long long> leafLastLevel; if (faninDbg) leafLastLevel.assign(nLeaf, -1);
    // ===== BATCH-PEEL (§58, gate SCT_BATCH_PEEL): leaf-major affected-Q enumeration =====
    // Drain a whole curLevel wave, group its (pattern,leaf) tasks LEAF-MAJOR, and run ONE
    // affected-Q DFS per (leaf,wave) instead of once per (pattern,leaf). Each candidate Q's
    // drop = Σ over the leaf's wave-thresholds of the PROVEN single-threshold delta-formula
    // (scWithTerms with addLow=pl over that threshold's pre-image paths). This is the SAME
    // per-(P,leaf) drop math the per-pattern path uses, only reordered (enumerate-Q-once,
    // then sum-over-thresholds) so the DFS over-enumeration (gen/nz=8-14, sec 57) and the
    // candidate generation amortize across the same-level co-hosting fan-in (10-1213, sec 57).
    // EXACTNESS: cores are order-independent within a level (all level-L patterns get core L;
    // drops to higher-level Q telescope: total = pre-level-sup - post-level-sup) + slot order
    // is not load-bearing (sec 54). Intra-wave drops clamp to curLevel and never re-bucket, so
    // marking the whole wave peeled up front and skipping them is bit-identical. Cascade (a
    // higher-level Q dropping to curLevel) re-drains curLevel. Only s>r+1 (the path that owns
    // a DFS to amortize); s=r+1 keeps the proven per-pattern witness-floor path. Default OFF.
    bool batchPeel = getenv("SCT_BATCH_PEEL") != nullptr;
    if (!ayMode && !ondemand && witnessTail >= 2 && (batchPeel || witnessTail > witCross)) {
        struct BTask { int lid, pi, k; };
        vector<int> wave;
        vector<BTask> taskLL;                          // (leaf,pattern,leaf-index) tasks for the wave
        vector<CCPath> coAll;                          // accumulated changed pre-image paths (one leaf)
        vector<int> coPlIdx;                           // per pre-image -> index into coPls
        vector<Vec> coPls;                             // distinct wave-thresholds touching the leaf
        vector<vector<pair<Vec,int>>> coTerms;         // per pre-image cached IE terms
        vector<CCPath> chgTmp;                          // slotForbidDiff output (reused)
        vector<pair<int,int>> plNZ;
        Vec plScr, qlScr, uEnv, qcand;
        vector<int> aff;
        vector<int> coStart;                           // per-threshold contiguous start offset in coAll
        vector<long long> cbDP;                        // candidate-count DP scratch (saturating)
        Vec sufScr, uThrScr;                           // fallback per-threshold: pl-suffix-sum / u-envelope
        long long cbCap = getenv("SCT_BATCH_CB_CAP") ? atoll(getenv("SCT_BATCH_CB_CAP")) : 128;
        long long cbBatch = 0, cbFallback = 0;         // leaves served by single-DFS vs per-threshold fallback
        long long printMark = 0;
        while (peeledN < npat) {
            auto it0 = bk.find(curLevel);
            while (it0 == bk.end() || it0->second.empty()) {
                if (++curLevel > maxKey + 1) break;
                it0 = bk.find(curLevel);
            }
            if (curLevel > maxKey + 1) break;
            // re-drain curLevel until empty (handles within-level cascade)
            while (true) {
                auto bit = bk.find(curLevel);
                if (bit == bk.end() || bit->second.empty()) break;
                wave.clear();
                for (int pi : bit->second) if (pats[pi].alive && pats[pi].key == curLevel) wave.push_back(pi);
                bit->second.clear();
                if (wave.empty()) break;
                taskLL.clear();
                for (int pi : wave) {
                    Pat &P = pats[pi];
                    P.alive = false; P.core = (double)curLevel; peeledN++;
                    coreDist[P.core] += (double)P.mult;
                    if (skipH1 && P.hostSz == 1) {        // SOURCE-SKIP (sec: M-exclusive witnesses)
                        bool aff2 = false;
                        for (int lid : leavesOf(pi)) if (hasH2[lid]) { aff2 = true; break; }
                        if (!aff2) continue;
                    }
                    const auto &pl2 = leavesOf(pi);
                    for (int kk = 0; kk < (int)pl2.size(); kk++)
                        if (!slotPaths[pl2[kk]].empty()) taskLL.push_back({pl2[kk], pi, kk});
                }
                if ((peeledN >> 12) != printMark) { printMark = peeledN >> 12;
                    fprintf(stderr, "[peel-batch] %lld/%lld lvl=%lld maxSplit=%zu t=%.1fs\n",
                            peeledN, npat, curLevel, maxSplit, secs(T5, Clock::now())); }
                std::sort(taskLL.begin(), taskLL.end(),
                          [](const BTask &a, const BTask &b) { return a.lid < b.lid; });
                size_t ti = 0;
                while (ti < taskLL.size()) {
                    int lid = taskLL[ti].lid;
                    size_t tj = ti; while (tj < taskLL.size() && taskLL[tj].lid == lid) tj++;
                    coAll.clear(); coPlIdx.clear(); coPls.clear(); coStart.clear();
                    for (size_t t = ti; t < tj; t++) {       // apply every threshold touching this leaf
                        int pi = taskLL[t].pi, kk = taskLL[t].k;
                        const Vec &pl = recomputePB ? (localB(pi, lid, plScr), (const Vec &)plScr)
                                                    : pbLocal[pi][kk];
                        plNZ.clear();
                        for (int c = 0; c < (int)pl.size(); c++) if (pl[c]) plNZ.push_back({c, (int)pl[c]});
                        slotVisits += (long long)slotPaths[lid].size();
                        auto _sa = Clock::now();
                        slotForbidDiff(lid, pl, plNZ, chgTmp);
                        tSFD += secs(_sa, Clock::now());
                        if (chgTmp.empty()) continue;
                        coStart.push_back((int)coAll.size());     // threshold's pre-images begin here
                        int plIdx = (int)coPls.size(); coPls.push_back(pl);
                        for (auto &p : chgTmp) { coAll.push_back(std::move(p)); coPlIdx.push_back(plIdx); }
                    }
                    if (coAll.empty()) { ti = tj; continue; }
                    coStart.push_back((int)coAll.size());          // sentinel end offset
                    int Mloc = coAll.front().m(), Tcap = coAll.front().T, F = (int)coPls.size();
                    coTerms.clear(); coTerms.reserve(coAll.size());
                    for (auto &p : coAll) coTerms.push_back(ccpath::inclusion_exclusion_terms(p.forbidden, p.m()));
                    uEnv.assign((size_t)Mloc, 0);          // envelope: a Q with a drop has ql<=some changed-path u
                    for (auto &p : coAll) for (int c = 0; c < Mloc; c++) if (p.u[c] > uEnv[c]) uEnv[c] = p.u[c];
                    ensureLeafMap(lid);
                    const auto &q2p = leafQ2pat[lid];
                    const auto &qsAll = leafPats[lid];
                    qcand.assign((size_t)Mloc, 0);
                    const int16_t *uEp = uEnv.data();
                    // confirm candidate ql against the leaf-pattern map, then accumulate its drop over the
                    // pre-image entries [e0,e1). Both enumeration strategies below feed THIS, so they are
                    // bit-identical (the drop partitions exactly by threshold: Σ_j drop[e0_j,e1_j) == drop[0,all)).
                    auto confirm = [&](const Vec &ql, uint64_t h, int e0, int e1) {
                        auto itc = q2p.find(h);
                        if (itc == q2p.end()) return;
                        for (int t : itc->second)
                            if (spanEqFP(lid, t, Mloc, qlScr, ql)) {
                                int qi = qsAll[t];
                                if (!pats[qi].alive) return;       // peeled (incl. the whole wave)
                                if (skipH1 && pats[qi].hostSz == 1) return;
                                double d = 0.0;
                                for (int e = e0; e < e1; e++) {
                                    const CCPath &p = coAll[e];
                                    const Vec &plE = coPls[coPlIdx[e]];
                                    bool ok = true; int sm = 0;
                                    for (int c = 0; c < Mloc; c++) {
                                        int v = (int)ql[c] > (int)plE[c] ? (int)ql[c] : (int)plE[c];
                                        if (v > (int)p.u[c]) { ok = false; break; }
                                        sm += v;
                                    }
                                    if (ok && sm <= (int)p.T) d += scWithTerms(p, coTerms[e], ql, &plE);
                                }
                                if (d != 0.0) { if (!seen[qi]) { seen[qi] = 1; aff.push_back(qi); } delta[qi] += d; }
                                return;
                            }
                    };
                    // CB = #{ql : Σ=r, ql<=uEnv} (saturating count-DP, O(Mloc*r^2)). When CB is small the
                    // single shared DFS amortizes the same-level fan-in; when it would blow up (high r / wide
                    // leaf, no single pl to prune by) we fall back to F cheap pl_j-pruned DFS == the proven
                    // per-pattern affected-update reorganized per threshold. cbCap tunable (SCT_BATCH_CB_CAP).
                    cbDP.assign((size_t)r + 1, 0); cbDP[0] = 1;
                    for (int c = 0; c < Mloc && cbDP[r] <= cbCap; c++) {
                        int uc = (int)uEp[c]; if (uc > r) uc = r;
                        for (int t = r; t >= 1; t--) { long long s2 = 0;
                            for (int y = 1; y <= uc && y <= t; y++) s2 += cbDP[t - y];
                            cbDP[t] += s2; if (cbDP[t] > cbCap) cbDP[t] = cbCap + 1;
                        }
                    }
                    if (cbDP[r] <= cbCap) {
                        cbBatch++;
                        int allE = (int)coAll.size();
                        auto dfsB = [&](auto &&self, int c, int rem, uint64_t h) -> void {   // single shared DFS
                            if (c == Mloc) { if (rem == 0) confirm(qcand, h, 0, allE); return; }
                            int cap = (int)uEp[c]; if (cap > rem) cap = rem;
                            for (int jj = 0; jj <= cap; jj++) {
                                qcand[c] = (int16_t)jj;
                                uint64_t hc = (h ^ ((uint64_t)(uint16_t)jj + 1)) * HMUL;
                                self(self, c + 1, rem - jj, hc);
                            }
                            qcand[c] = 0;
                        };
                        dfsB(dfsB, 0, r, 1469598103934665603ULL);
                    } else {
                        cbFallback++;
                        for (int j = 0; j < F; j++) {        // per-threshold pl_j-pruned DFS (cheap O(1)/node)
                            const Vec &plj = coPls[j];
                            sufScr.assign((size_t)Mloc + 1, 0);
                            for (int c = Mloc - 1; c >= 0; c--) sufScr[c] = (int16_t)((int)sufScr[c + 1] + (int)plj[c]);
                            uThrScr.assign((size_t)Mloc, 0);     // threshold-j u-envelope over its own pre-images
                            for (int e = coStart[j]; e < coStart[j + 1]; e++) {
                                const Vec &pu = coAll[e].u;
                                for (int c = 0; c < Mloc; c++) if (pu[c] > uThrScr[c]) uThrScr[c] = pu[c];
                            }
                            const int16_t *plp = plj.data(), *ujp = uThrScr.data(), *sfp = sufScr.data();
                            int e0 = coStart[j], e1 = coStart[j + 1];
                            auto dfsT = [&](auto &&self, int c, int rem, int acc, uint64_t h) -> void {
                                if (c == Mloc) { if (rem == 0) confirm(qcand, h, e0, e1); return; }
                                if (acc + (int)sfp[c] > Tcap) return;       // Σmax(pl_j,ql)<=T (per-pattern prune)
                                int cap = (int)ujp[c]; if (cap > rem) cap = rem;
                                int plc = (int)plp[c];
                                for (int jj = 0; jj <= cap; jj++) {
                                    qcand[c] = (int16_t)jj;
                                    int mx = plc > jj ? plc : jj;
                                    uint64_t hc = (h ^ ((uint64_t)(uint16_t)jj + 1)) * HMUL;
                                    self(self, c + 1, rem - jj, acc + mx, hc);
                                }
                                qcand[c] = 0;
                            };
                            dfsT(dfsT, 0, r, 0, 1469598103934665603ULL);
                        }
                    }
                    ti = tj;
                }
                for (int qi : aff) {                        // APPLY once per wave (telescoped drop)
                    seen[qi] = 0;
                    double ns = pats[qi].sup - delta[qi];
                    delta[qi] = 0.0;
                    long long nk = (long long)llround(ns);
                    if (nk < curLevel) nk = curLevel;
                    if (nk != pats[qi].key) { pats[qi].sup = ns; pats[qi].key = nk; bk[nk].push_back(qi); }
                }
                aff.clear();
            }
        }
        fprintf(stderr, "[peel-batch] leaf-instances: single-DFS=%lld fallback=%lld (%.1f%% batched, CB-cap=%lld)\n",
                cbBatch, cbFallback, (cbBatch + cbFallback) ? 100.0 * cbBatch / (cbBatch + cbFallback) : 0.0, cbCap);
    } else {
    while (peeledN < npat) {
        auto it = bk.find(curLevel);
        while (it == bk.end() || it->second.empty()) {
            if (++curLevel > maxKey + 1) break;
            it = bk.find(curLevel);
        }
        if (curLevel > maxKey + 1) break;
        int pi = it->second.back(); it->second.pop_back();
        Pat &P = pats[pi];
        if (!P.alive || P.key != curLevel) continue;
        P.alive = false; P.core = (double)curLevel; peeledN++;
        if ((peeledN & 0x3FF) == 0)
            fprintf(stderr, "[peel] %lld/%lld lvl=%lld maxSplit=%zu t=%.1fs\n",
                    peeledN, npat, curLevel, maxSplit, secs(T5, Clock::now()));
        coreDist[P.core] += (double)P.mult;

        // SOURCE-SKIP: a |host|=1 pattern whose leaves are all pure-|host|=1
        // removes nothing that any |host|>=2 pattern depends on, so skip its
        // entire affected-update (slotForbidDiff + DFS). Correct because its
        // witnesses are M-exclusive (a |host|>=2 pattern is never hosted in
        // these leaves, so never has a witness there).
        const auto &pleaf = leavesOf(pi);   // §94: hosting leaves computed ONCE (reused by affectsH2 + the main loop)
        if (skipH1 && P.hostSz == 1) {
            bool affectsH2 = false;
            for (int lid : pleaf) if (hasH2[lid]) { affectsH2 = true; break; }
            if (!affectsH2) continue;
        }

        // INCREMENTAL, EXACT update. Only the leaves hosting P change, and
        // within a leaf only the PATHS where P's threshold applies change. For
        // an affected Q its support drop from leaf L equals
        //   sum_{changed paths p} [ SC(p_before, m_Q) - SC(p_after, m_Q) ]
        // (unchanged paths give the same count before and after, so they
        // cancel). SC = support_count carries the C(n-b,y-b) tuple reweighting
        // and the forbidden inclusion-exclusion, so this is bit-identical to a
        // full recompute over all of Q's leaves, but its cost is bounded by the
        // few paths P actually touches -- independent of the split-set size.
        // (A componentwise-max shortcut is NOT valid: the drop is not
        // SC(L, max(m_P,m_Q)) because of that C(n-b,y-b) reweighting.)
        vector<int> aff;
        Vec plScr, qlScr;                              // recompute scratch: P-side (held across k-body), Q-side (per confirm)
        vector<CCPath> chgOld;                         // pre-insertion snapshots
        vector<vector<pair<Vec,int>>> chgOldTerms;     // cached IE terms (pre-insert)
        vector<pair<int,int>> plNZ;                    // sparse nonzeros of m_P local
        for (size_t k = 0; k < pleaf.size(); k++) {
            int lid = pleaf[k];
            if (slotPaths[lid].empty()) continue;      // leaf fully peeled: no witnesses
            if (faninDbg) { fanA++; if (leafLastLevel[lid] != curLevel) { fanB++; leafLastLevel[lid] = curLevel; } }
            const Vec &pl = recomputePB ? (localB(pi, lid, plScr), (const Vec &)plScr)
                                        : pbLocal[pi][k];   // m_P local to lid (== a_p, h=0)
            int Mloc = (int)pl.size();
            // sparse support of m_P (positions where it is nonzero) -- the only
            // positions the impossible / feasibility tests depend on.
            plNZ.clear();
            for (int c = 0; c < Mloc; c++) if (pl[c]) plNZ.push_back({c, (int)pl[c]});
            // ---- a_Y DIRECT path (§88): explicit dead-set, NO slotForbidDiff / antichain / split ----
            // Witnesses dying when P peels are Y = pl + δ (Σδ=t). Y is alive iff not already marked dead by an
            // earlier-peeled sub-pattern. We enumerate those Y (same DFS as the witness path), and for each ALIVE one
            // mark it dead (one hash-set insert) and credit the witness-major drop to every Q = Y - γ. No box scan,
            // no forbidden IE, no controlled_split -> deletes the entire slotForbidDiff churn. Bit-identical: the
            // dead set is exactly {Y : Y dominates some peeled pattern}, the same set the antichain represents.
            if (ayMode && witnessActive) {
                const CCPath &box = slotPaths[lid][0];     // original leaf box (never mutated in a_Y mode)
                witInst++; witMSum += Mloc; if (Mloc > witMMax) witMMax = Mloc;
                if (!ondemand) ensureLeafMap(lid);           // §3b: ondemand resolves Q via the global hash, no per-leaf map
                const auto &q2p = leafQ2pat[lid];
                const auto &qsAll = leafPats[lid];
                const Vec &nn = box.n;                       // leaf class sizes (n_b)
                auto &dead = deadY[lid];
                Yscr = pl;                                   // scratch: pl -> +δ (Y) -> -γ (Q)
                const int16_t *uEp = box.u.data();
                const int16_t *ellp = box.ell.data();
                auto credit = [&](double w) {                // credit Q (= current Yscr); nAlive==1 by construction
                    if (w == 0.0) return;
                    int qi;
                    if (ondemand) { qi = globalLookup(lid, Yscr, Mloc); if (qi < 0) return; }
                    else {
                        auto it = q2p.find(hashVec(Yscr)); if (it == q2p.end()) return; qi = -1;
                        for (int t : it->second) if (spanEqFP(lid, t, Mloc, qlScr, Yscr)) { qi = qsAll[t]; break; }
                        if (qi < 0) return;
                    }
                    if (qi != pi && pats[qi].alive && !(skipH1 && pats[qi].hostSz == 1)) {
                        if (!seen[qi]) { seen[qi] = 1; aff.push_back(qi); }
                        delta[qi] += w;
                    }
                };
                auto remGamma = [&](auto &&self, int start, int rem, double w) -> void {
                    if (rem == 0) { credit(w); return; }
                    for (int b = start; b < Mloc; b++) {
                        int Yb = (int)Yscr[b];
                        if (Yb < 1) continue;
                        int maxm = Yb < rem ? Yb : rem;
                        for (int m = 1; m <= maxm; m++) {
                            Yscr[b] = (int16_t)(Yb - m);                      // Q_b = Y_b - m
                            int avail = (int)nn[b] - (Yb - m);                // n_b - Q_b
                            double wf = (m == 1) ? (double)avail : ccpath_ncr(avail, m);
                            self(self, b + 1, rem - m, w * wf);
                        }
                        Yscr[b] = (int16_t)Yb;
                    }
                };
                auto addDelta = [&](auto &&self, int start, int rem) -> void {
                    if (rem == 0) {
                        for (int k = 0; k < Mloc; k++)                        // feasible witness: ell <= Y <= u
                            if ((int)ellp[k] > (int)Yscr[k] || (int)Yscr[k] > (int)uEp[k]) return;
                        if (!dead.insert(hashVec(Yscr))) return;              // already dead -> no drop
                        remGamma(remGamma, 0, witnessTail, 1.0);
                        return;
                    }
                    for (int a = start; a < Mloc; a++) {
                        int room = (int)uEp[a] - (int)Yscr[a];               // Y[a] may grow (Y <= u)
                        if (room < 1) continue;
                        int maxm = room < rem ? room : rem;
                        int Ya = (int)Yscr[a];
                        for (int m = 1; m <= maxm; m++) { Yscr[a] = (int16_t)(Ya + m); self(self, a + 1, rem - m); }
                        Yscr[a] = (int16_t)Ya;
                    }
                };
                addDelta(addDelta, 0, witnessTail);
                continue;                                    // leaf done; no slotForbidDiff
            }
            // Record P (updates the stored slot via split) and capture the CHANGED
            // OLD paths (the pre-insertion snapshots where P's threshold applies).
            { auto _sa = Clock::now(); slotVisits += (long long)slotPaths[lid].size();
              slotForbidDiff(lid, pl, plNZ, chgOld); tSFD += secs(_sa, Clock::now()); }
            if (chgOld.empty()) continue;              // P touched nothing here
            // ---- ADAPTIVE PER-LEAF choice: witness vs general DFS (both bit-identical) ----
            // t=1,2 always witness (win across the density spectrum). For t>=witGateMinT estimate the two costs
            // cheaply: Wδ = #feasible witnesses (compositions of t with δ_c<=uEnv_c-pl_c) drives the witness
            // box-scan; CB = #affected-Q candidates (compositions of r with Q_c<=uEnv_c) drives the general DFS.
            // Witness per-unit (box scan) << general per-unit (scWithTerms DP), so take witness while Wδ<=CB*witK.
            bool wLeaf = witnessActive;
            if (wLeaf && witnessTail >= witGateMinT) {
                int Mw = (int)pl.size();
                uEnv.assign((size_t)Mw, 0);
                for (auto &p : chgOld) for (int c = 0; c < Mw; c++) if (p.u[c] > uEnv[c]) uEnv[c] = p.u[c];
                const long long SAT = 1LL << 40;
                wdDP.assign((size_t)witnessTail + 1, 0); wdDP[0] = 1;     // Wδ: bounded compositions of t
                for (int c = 0; c < Mw; c++) { int room = (int)uEnv[c] - (int)pl[c]; if (room < 0) room = 0; if (room > witnessTail) room = witnessTail;
                    for (int tt = witnessTail; tt >= 1; tt--) { long long s2 = 0; for (int y = 1; y <= room && y <= tt; y++) s2 += wdDP[tt - y];
                        wdDP[tt] += s2; if (wdDP[tt] > SAT) wdDP[tt] = SAT; } }
                cbDP2.assign((size_t)r + 1, 0); cbDP2[0] = 1;            // CB: bounded compositions of r
                for (int c = 0; c < Mw; c++) { int uc = (int)uEnv[c]; if (uc > r) uc = r;
                    for (int tt = r; tt >= 1; tt--) { long long s2 = 0; for (int y = 1; y <= uc && y <= tt; y++) s2 += cbDP2[tt - y];
                        cbDP2[tt] += s2; if (cbDP2[tt] > SAT) cbDP2[tt] = SAT; } }
                wLeaf = ((double)wdDP[witnessTail] <= (double)cbDP2[r] * witK);
                if (wLeaf) witGateW++; else witGateG++;     // does the gate actually split per leaf?
            }
            if (wLeaf) {
                // ---- UNIFIED WITNESS-MAJOR fast path, tail t = s-r (bit-exact vs scWithTerms, §72-74) ----
                // Dying witness Y = pl + δ (Σδ=t); each affected Q = Y - γ (Σγ=t, γ<=Y) gets drop
                //   Π_{distinct b in γ} C(n_b - Q_b, mult_b), summed over alive boxes. n is leaf-constant
                // (splits preserve n) so the weight is box-independent -> only the alive-box COUNT matters.
                // "Alive" = ell<=Y<=u AND Y not>= any forbidden (a direct dominance test, NO inclusion-
                // exclusion). δ,γ are enumerated as non-decreasing class multisets (each multiset once).
                // t=1 == the old witness-floor, t=2 == the old witness-major; one path for all tails.
                // Output-sensitive: NO DP, NO IE, NO DFS. Q==P (γ==δ) is excluded by the qi!=pi test.
                witInst++; witMSum += Mloc; if (Mloc > witMMax) witMMax = Mloc;
                ensureLeafMap(lid);
                const auto &q2p = leafQ2pat[lid];
                const auto &qsAll = leafPats[lid];
                const Vec &nn = chgOld.front().n;       // leaf class sizes (constant across the leaf's boxes)
                uEnv.assign((size_t)Mloc, 0);           // u-envelope over chgOld -> cheap δ feasibility prune
                for (auto &p : chgOld) for (int c = 0; c < Mloc; c++) if (p.u[c] > uEnv[c]) uEnv[c] = p.u[c];
                Yscr = pl;                              // scratch: pl, then +δ (build Y), then -γ (build Q)
                const int16_t *uEp = uEnv.data();
                // credit Q (= current Yscr) with the γ-weight w times the alive-box count.
                auto credit = [&](double w, int nAlive) {
                    if (w == 0.0) return;
                    auto it = q2p.find(hashVec(Yscr));
                    if (it == q2p.end()) return;
                    for (int t : it->second)
                        if (spanEqFP(lid, t, Mloc, qlScr, Yscr)) {
                            int qi = qsAll[t];
                            if (qi != pi && pats[qi].alive && !(skipH1 && pats[qi].hostSz == 1)) {
                                if (!seen[qi]) { seen[qi] = 1; aff.push_back(qi); }
                                delta[qi] += w * (double)nAlive;
                            }
                            return;
                        }
                };
                // remove t units from Y as a non-decreasing class multiset -> Q; weight = Π C(n_b-Q_b, m).
                auto remGamma = [&](auto &&self, int start, int rem, double w, int nAlive) -> void {
                    if (rem == 0) { credit(w, nAlive); return; }
                    for (int b = start; b < Mloc; b++) {
                        int Yb = (int)Yscr[b];
                        if (Yb < 1) continue;
                        int maxm = Yb < rem ? Yb : rem;
                        for (int m = 1; m <= maxm; m++) {
                            Yscr[b] = (int16_t)(Yb - m);                          // Q_b = Y_b - m
                            int avail = (int)nn[b] - (Yb - m);                    // n_b - Q_b
                            double wf = (m == 1) ? (double)avail : ccpath_ncr(avail, m);  // C(n_b-Q_b, m), m=1 fast
                            self(self, b + 1, rem - m, w * wf, nAlive);
                        }
                        Yscr[b] = (int16_t)Yb;                                    // restore
                    }
                };
                // add t units to pl as a non-decreasing class multiset -> Y; scan boxes; then enumerate γ.
                auto addDelta = [&](auto &&self, int start, int rem) -> void {
                    if (rem == 0) {
                        int nAlive = 0;                                          // chgOld boxes where Y is alive
                        for (size_t z = 0; z < chgOld.size(); z++) {
                            const CCPath &p = chgOld[z];
                            bool ok = true;
                            for (int k = 0; k < Mloc; k++)
                                if ((int)p.ell[k] > (int)Yscr[k] || (int)Yscr[k] > (int)p.u[k]) { ok = false; break; }
                            if (!ok) continue;
                            bool dead = false;
                            for (auto &a : p.forbidden) { bool le = true;
                                for (int k = 0; k < Mloc; k++) if ((int)a[k] > (int)Yscr[k]) { le = false; break; }
                                if (le) { dead = true; break; } }
                            if (!dead) nAlive++;
                        }
                        if (nAlive > 0) remGamma(remGamma, 0, witnessTail, 1.0, nAlive);
                        return;
                    }
                    for (int a = start; a < Mloc; a++) {
                        int room = (int)uEp[a] - (int)Yscr[a];                   // Y[a] may grow this much (Y<=uEnv)
                        if (room < 1) continue;
                        int maxm = room < rem ? room : rem;
                        int Ya = (int)Yscr[a];
                        for (int m = 1; m <= maxm; m++) {
                            Yscr[a] = (int16_t)(Ya + m);
                            self(self, a + 1, rem - m);
                        }
                        Yscr[a] = (int16_t)Ya;                                   // restore
                    }
                };
                addDelta(addDelta, 0, witnessTail);
                continue;                                                       // leaf done; skip DFS path
            }
            // ---- DELTA-FORMULA incremental drop (no chgNew sweep) ----
            // For a changed path p, P inserts threshold a_p (= P's tuple local to
            // p) into p.forbidden. Q's support drop from p equals the witnesses of
            // p that were alive (>= no f in A_p), >= b_Q, and are now killed by
            // a_p, i.e. also >= a_p. That count is exactly support_count of p with
            // its ORIGINAL antichain A_p evaluated at b = componentwise-max(b_Q,
            // a_p). So one scWithTerms over chgOld (with b raised to max(ql,a_p))
            // gives the drop directly -- the split children (chgNew) cancel
            // algebraically and never need to be summed per Q. This halves the
            // scWithTerms work and removes the entire chgNew x Q sweep.
            // Bit-identical to (Σ SC(chgOld) − Σ SC(chgNew)); verified vs brute.
            chgOldTerms.clear(); chgOldTerms.reserve(chgOld.size());
            for (auto &p : chgOld)
                chgOldTerms.push_back(ccpath::inclusion_exclusion_terms(p.forbidden, p.m()));
            // a_p (P's local threshold) == pl for every changed path (h=0), so we
            // use pl directly as the delta-formula addLow below -- no per-path
            // tuple_to_threshold needed.
            // envelope over the CHANGED-old paths (the only ones that can drop):
            // a Q with a nonzero drop must have a witness >= max(m_P,m_Q) in a
            // changed path. Prune Q if max(m_P,m_Q) exceeds chgOld's u-envelope.
            uEnv.assign((size_t)Mloc, 0); int Tcap = chgOld.front().T;
            for (auto &p : chgOld)
                for (int c = 0; c < Mloc; c++) if (p.u[c] > uEnv[c]) uEnv[c] = p.u[c];
            // ---- ENUMERATE affected Q directly from m_P (no full-leaf scan) ----
            // An affected Q must own a dying witness y in a changed path, i.e.
            //   y >= max(m_P, m_Q),  ell <= y <= u,  Σy = T.
            // Hence m_Q <= y <= uEnv and Σ_c max(m_P[c], m_Q[c]) <= Tcap. We DFS
            // over the leaf's local classes generating exactly the m_Q meeting
            // those bounds (each once, in canonical order), look each up in the
            // per-leaf map, and apply the delta-formula drop. This replaces the
            // O(#patterns-in-leaf) scan -- of which ~99.8% failed the bound -- by
            // an enumeration whose size is the genuinely-affected candidate set.
            ensureLeafMap(lid);
            const auto &q2p = leafQ2pat[lid];
            const auto &qsAll = leafPats[lid];
            // suffix sum of m_P: min future contribution of classes >= i to the
            // witness sum (max(pl,ql) >= pl), used to prune the DFS early.
            sufPl.assign((size_t)Mloc + 1, 0);
            for (int c = Mloc - 1; c >= 0; c--) sufPl[c] = (int16_t)(sufPl[c + 1] + (int)pl[c]);
            qcand.assign((size_t)Mloc, 0);
            // process drop for the candidate ql at leaf-pattern index t (already
            // confirmed leafPatLocB[lid][t] == ql). qi = leafPats[lid][t].
            auto applyIdx = [&](const Vec &ql, int t) {
                int qi = qsAll[t];
                if (qi == pi || !pats[qi].alive) return;
                if (skipH1 && pats[qi].hostSz == 1) return;  // peels at L_M regardless
                dbgHit++;                               // a real candidate pattern reached
                double d = 0.0;                         // drop, via delta formula
                for (size_t z = 0; z < chgOld.size(); z++) {
                    const CCPath &p = chgOld[z];
                    // a_p = pl (P's local threshold; h=0). drop_p = count(y: alive
                    // in A_p, y >= max(ql, pl)) weighted by C(n - ql, y - ql). The
                    // weight base STAYS ql (addLow=pl only raises the floor) --
                    // reweighting-correct delta. Feasible only if max(ql,pl)<=p.u.
                    bool ok = true; int sm = 0;
                    for (int c = 0; c < Mloc; c++) {
                        int v = (int)ql[c] > (int)pl[c] ? (int)ql[c] : (int)pl[c];
                        if (v > (int)p.u[c]) { ok = false; break; }
                        sm += v;
                    }
                    if (ok && sm <= (int)p.T) d += scWithTerms(p, chgOldTerms[z], ql, &pl);
                }
                if (d == 0.0) return;
                dbgNZ++;                                // candidate with a genuine nonzero drop
                if (!seen[qi]) { seen[qi] = 1; aff.push_back(qi); }
                delta[qi] += d;
            };
            // look up a complete candidate (by rolling hash) and confirm exactly.
            auto applyCand = [&](const Vec &ql, uint64_t h) {
                dbgGen++;                               // a complete ql candidate generated
                auto it = q2p.find(h);
                if (it == q2p.end()) return;
                for (int t : it->second) if (spanEqFP(lid, t, Mloc, qlScr, ql)) { applyIdx(ql, t); return; }
            };
            // DFS over leaf classes: place ql[c] in [0, min(uEnv[c],rem)], track
            // rem = r left, acc = Σ_{<c} max(pl,ql), and the running rolling hash.
            // Self-recursive lambda (Y-combinator) -> fully inlinable, no
            // std::function indirection. Each complete ql is looked up once.
            const int16_t *plp = pl.data(); const int16_t *uEp = uEnv.data();
            const int16_t *sfp = sufPl.data();
            auto dfs = [&](auto &&self, int c, int rem, int acc, uint64_t h) -> void {
                if (c == Mloc) { if (rem == 0) applyCand(qcand, h); return; }
                if (acc + (int)sfp[c] > Tcap) return;   // min future sum already too big
                int cap = (int)uEp[c]; if (cap > rem) cap = rem;
                int plc = (int)plp[c];
                for (int j = 0; j <= cap; j++) {
                    qcand[c] = (int16_t)j;
                    int mx = plc > j ? plc : j;
                    uint64_t hc = (h ^ ((uint64_t)(uint16_t)j + 1)) * HMUL;
                    self(self, c + 1, rem - j, acc + mx, hc);
                }
                qcand[c] = 0;
            };
            if (!dfsPrune) {
                dfs(dfs, 0, r, 0, 1469598103934665603ULL);
            } else {
                // ---- per-path COVERAGE + feasibility subtree prune (SCT_DFS_PRUNE) ----
                // A chgOld path z scores 0 for candidate ql (so contributes nothing to
                // the drop) iff it is u-INFEASIBLE (max(ql,pl) </=  u_z) OR COVERED
                // (some forbidden a_z <= max(ql,pl) -- exactly scWithTerms' early-out,
                // lines ~877-890). Both are MONOTONE as ql grows in descent, so once a
                // path is dead it stays dead in the whole subtree; when ALL paths are
                // dead the subtree's candidates all have d==0 and are pruned -- removing
                // exactly applyIdx's d==0 set, so cores stay BIT-IDENTICAL (applyIdx is
                // the oracle). This targets the ~90% forbidden-COVERAGE waste (vs the
                // u-only prune, which the uEnv union already made tight). KMAX=1 default
                // => <=1 forbidden/path; paths with !=1 forbidden skip coverage
                // (conservative: an over-keep only wastes a lookup, never breaks cores).
                int nChg = (int)chgOld.size();
                if ((int)chgU.size() < nChg) { chgU.resize(nChg); fbA.resize(nChg); fbCrit.resize(nChg); fbHas.resize(nChg); }
                for (int z = 0; z < nChg; z++) {
                    chgU[z] = chgOld[z].u.data();
                    if (chgOld[z].forbidden.size() == 1) {
                        const Vec &a = chgOld[z].forbidden[0];
                        fbA[z] = a.data(); fbHas[z] = 1;
                        int cm = -1;                              // max coord where a_z > pl (critical)
                        for (int c = 0; c < Mloc; c++) if ((int)a[(size_t)c] > (int)pl[c]) cm = c;
                        fbCrit[z] = cm;                           // tail a_z<=pl from coord cm+1
                    } else { fbA[z] = nullptr; fbHas[z] = 0; fbCrit[z] = Mloc; }
                }
                long long need = (long long)(Mloc + 1) * nChg;
                if ((long long)uLiveBuf.size() < need) uLiveBuf.assign((size_t)need, 1);
                if ((long long)covBuf.size()   < need) covBuf.assign((size_t)need, 1);
                for (int z = 0; z < nChg; z++) { uLiveBuf[z] = 1; covBuf[z] = fbHas[z] ? 1 : 0; }
                const int16_t **uP = chgU.data(); const int16_t **aP = fbA.data();
                const int *crit = fbCrit.data(); const char *hasF = fbHas.data();
                char *ul = uLiveBuf.data(); char *cv = covBuf.data();
                auto dfsP = [&](auto &&self, int c, int rem, int acc, uint64_t h, int base) -> void {
                    if (c == Mloc) { if (rem == 0) applyCand(qcand, h); return; }
                    if (acc + (int)sfp[c] > Tcap) return;        // Σ bound (path-independent, T=s)
                    const char *ulc = ul + base, *cvc = cv + base;
                    bool any = false;                            // any path still feasible AND uncovered?
                    for (int z = 0; z < nChg; z++) {
                        bool dead = !ulc[z] || (hasF[z] && cvc[z] && crit[z] < c);
                        if (!dead) { any = true; break; }
                    }
                    if (!any) return;                            // whole subtree scores d==0
                    int cap = (int)uEp[c]; if (cap > rem) cap = rem;
                    int plc = (int)plp[c];
                    char *uln = ul + base + nChg, *cvn = cv + base + nChg;
                    for (int j = 0; j <= cap; j++) {
                        int mx = plc > j ? plc : j;
                        for (int z = 0; z < nChg; z++) {
                            uln[z] = ulc[z] && (mx <= (int)uP[z][c]);
                            cvn[z] = hasF[z] ? (cvc[z] && ((int)aP[z][c] <= mx)) : 0;
                        }
                        qcand[c] = (int16_t)j;
                        uint64_t hc = (h ^ ((uint64_t)(uint16_t)j + 1)) * HMUL;
                        self(self, c + 1, rem - j, acc + mx, hc, base + nChg);
                    }
                    qcand[c] = 0;
                };
                dfsP(dfsP, 0, r, 0, 1469598103934665603ULL, 0);
            }
        }
        for (int qi : aff) {
            seen[qi] = 0;
            double ns = pats[qi].sup - delta[qi];     // exact incremental drop
            delta[qi] = 0.0;
            long long nk = (long long)llround(ns);
            if (nk < curLevel) nk = curLevel;          // monotone clamp
            if (nk != pats[qi].key) { pats[qi].sup = ns; pats[qi].key = nk; bk[nk].push_back(qi); }
        }
    }
    }
    auto T6 = Clock::now();
    memCk("after-peel(+index)");
    if (getenv("MEM_BREAKDOWN")) {                              // actual bytes of each major structure (post-peel, deadY full)
        double deadYB = 0, deadYent = 0;
        for (auto &d : deadY) { deadYB += (double)d.t.capacity() * 8; deadYent += (double)d.cnt; }
        double patsB = (double)pats.size() * (double)sizeof(Pat), hostInc = 0, hostB = 0, compB = 0, csB = 0;
        for (auto &P : pats) { hostInc += (double)P.host.size(); hostB += (double)P.host.capacity() * 4;
            compB += (double)P.comp.capacity() * 8; csB += (double)P.classSet.capacity() * 4; }
        double leafPatsB = 0, leafPatInc = 0, leafFlatB = 0;
        for (auto &v : leafPats) { leafPatsB += (double)v.capacity() * 4 + 24; leafPatInc += (double)v.size(); }
        for (auto &v : leafFlat) leafFlatB += (double)v.capacity() * 2 + 24;
        double slotB = 0;
        for (auto &sp : slotPaths) for (auto &b : sp) slotB += (double)(b.ell.capacity() + b.u.capacity() + b.n.capacity()) * 2 + (double)b.classIds.capacity() * 4 + 80;
        fprintf(stderr, "[mem-bd] deadY=%.0fMB(%.0f ent) | pats=%.0fMB struct + host=%.0fMB(%.0f inc) comp=%.0fMB classSet=%.0fMB | leafPats=%.0fMB(%.0f inc) leafFlat=%.0fMB | slotPaths=%.0fMB\n",
                deadYB / 1e6, deadYent, patsB / 1e6, hostB / 1e6, hostInc, compB / 1e6, csB / 1e6, leafPatsB / 1e6, leafPatInc, leafFlatB / 1e6, slotB / 1e6);
    }
    if (witDbg) fprintf(stderr, "[wit] tail=%d leaf-instances=%lld avg-M=%.1f max-M=%lld | gate: witness=%lld general=%lld (%.1f%% fell back)\n",
            witnessTail, witInst, witInst ? (double)witMSum / witInst : 0.0, witMMax,
            witGateW, witGateG, (witGateW + witGateG) ? 100.0 * witGateG / (witGateW + witGateG) : 0.0);
    fprintf(stderr, "[profile] peel=%.2fs  slotForbidDiff=%.2fs (%.0f%%)  rest(affected-update)=%.2fs  slot-path-visits=%lld\n",
            secs(T5,T6), tSFD, 100.0*tSFD/max(1e-9,secs(T5,T6)), secs(T5,T6)-tSFD, slotVisits);
    if (idxDbg) fprintf(stderr, "[idx-dbg] pivot-scan(Σ mv-thr)=%lld  candidates-filtered=%lld  affected-out=%lld | filter waste=%.1fx, pivscan/out=%.1fx\n",
            ixPivScan, ixCand, ixOut, ixOut ? (double)ixCand/ixOut : 0.0, ixOut ? (double)ixPivScan/ixOut : 0.0);
    if (sfdDbg) fprintf(stderr, "[sfd-dbg] tested=%lld affected=%lld (%.2f%%) coord-tests=%lld (%.2f/test) fail-on-1st=%lld (%.1f%% of skips)\n",
            slotVisits, sfdAff, slotVisits ? 100.0*sfdAff/slotVisits : 0.0, sfdCoordTests,
            slotVisits ? (double)sfdCoordTests/slotVisits : 0.0, sfdFailFirst,
            (slotVisits - sfdAff) ? 100.0*sfdFailFirst/(slotVisits - sfdAff) : 0.0);
    if (suppDbg) {
        long long tc = 0, ta = 0, ts = 0;
        for (int k = 0; k < 12; k++) { tc += suppCalls[k]; ta += suppAff[k]; ts += suppSplit[k]; }
        fprintf(stderr, "[supp] forbid-insert calls=%lld affected-paths=%lld splits=%lld\n", tc, ta, ts);
        for (int k = 1; k < 12; k++) if (suppCalls[k] || suppAff[k])
            fprintf(stderr, "[supp]  |supp|=%d: calls=%lld(%.1f%%) aff=%lld(%.1f%%) splits=%lld(%.1f%%)\n", k,
                    suppCalls[k], tc ? 100.0*suppCalls[k]/tc : 0.0, suppAff[k], ta ? 100.0*suppAff[k]/ta : 0.0,
                    suppSplit[k], ts ? 100.0*suppSplit[k]/ts : 0.0);
        fprintf(stderr, "[supp]  >>> SINGLE-CLASS share: calls=%.1f%%  aff-work=%.1f%%  splits=%.1f%%  (hybrid win ceiling)\n",
                tc ? 100.0*suppCalls[1]/tc : 0.0, ta ? 100.0*suppAff[1]/ta : 0.0, ts ? 100.0*suppSplit[1]/ts : 0.0);
    }
    if (slotIdxVerify) { fprintf(stderr, "[slot-idx] verify: %lld mismatched calls %s\n", idxMismatch, idxMismatch ? "FAIL" : "OK");
        if (idxMismatch) return 5; }
    if (faninDbg) fprintf(stderr, "[fanin] touches(pattern,leaf)=%lld distinct(leaf,level)=%lld  fanin=%.2f\n",
            fanA, fanB, fanB ? (double)fanA / fanB : 0.0);
    if (cfDbg) {
        fprintf(stderr, "[cf] DP IE-terms=%lld  binding-upper histogram (0..):", cfTerms);
        for (int i = 0; i <= 10; i++) if (cfBindHist[i]) fprintf(stderr, " %d:%.1f%%", i, cfTerms ? 100.0*cfBindHist[i]/cfTerms : 0.0);
        fprintf(stderr, "\n[cf] avg closed-form/DP work ratio=%.4f (closed-form wins per term if <1)\n",
                cfTerms ? cfFormTermsVsDP / cfTerms : 0.0);
    }
    printf("[sct-peel] peel=%.2fs  peeled=%lld/%lld  maxSplit(split-set)=%zu\n",
           secs(T5,T6), peeledN, npat, maxSplit);
    fprintf(stderr, "[sct-peel] dbg dfsPrune=%d cand_gen=%lld hit=%lld nz=%lld  gen/nz=%.1f\n",
            (int)dfsPrune, dbgGen, dbgHit, dbgNZ, dbgNZ ? (double)dbgGen/dbgNZ : 0.0);
    printf("[sct-peel] TIMING MCE=%.2f enum=%.2f sct-build+maps=%.2f peel=%.2f total=%.2f\n",
           secs(T1,T2), secs(T3,T4), secs(Tqg0,T5), secs(T5,T6), secs(T1,T6));
    // fold in the r-mergeable direct-assigned cores
    for (auto &kv : directCoreDist) coreDist[kv.first] += kv.second;
    double maxCore = 0; for (auto &kv : coreDist) maxCore = max(maxCore, kv.first);
    printf("[sct-peel] Max core: %.0f\n", maxCore);
    for (auto &kv : coreDist) printf("core=%.0f count=%.0f\n", kv.first, kv.second);
    return 0;
}
