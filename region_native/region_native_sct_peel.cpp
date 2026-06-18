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
        vector<int> cand;
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

int main(int argc, char **argv) {
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
        vector<vector<int>> vr(g.n);
        for (int i = 0; i < nRall; i++) for (int v : regions[i]) vr[v].push_back(i);
        vector<int> cnt(nRall, 0); vector<int> dirty; dirty.reserve(256);
        vector<char> mergeable(nRall, 0);
        for (int M = 0; M < nRall; M++) {
            for (int v : regions[M]) for (int o : vr[v]) if (o != M) {
                if (cnt[o] == 0) dirty.push_back(o); cnt[o]++;
            }
            int mx = 0; for (int o : dirty) { if (cnt[o] > mx) mx = cnt[o]; cnt[o] = 0; }
            dirty.clear();
            if (mx < r) mergeable[M] = 1;
        }
        vector<vector<int>> active;
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
    unordered_map<string, int> profKey;
    vector<int> classOf(g.n, -1);
    vector<vector<int>> classRegions;   // class -> sorted region ids (its profile)
    vector<int> classSize;
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
        vector<int> out; size_t i = 0, j = 0;
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
        vector<int> classSet;        // sorted class ids of comp (subset tests)
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
        vector<Node> inter;
        for (auto &N : B) {
            vector<int> cs = interClasses(M.classes, N.classes);
            int vs = classesSize(cs);
            if (vs >= s) inter.push_back({std::move(cs), vs});
        }
        if (inter.size() > 1) {
            sort(inter.begin(), inter.end(),
                 [](const Node &a, const Node &b){ return a.classes.size() > b.classes.size(); });
            vector<Node> keep;
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

    // ---- enumerate ALL patterns (canonical home), including |host|=1 ----
    {
        int curRid = 0;
        vector<pair<int,int>> cur;
        std::function<void(int,const vector<int>&,int)> enumAll =
            [&](int idx, const vector<int> &cls, int rem) {
            if (rem == 0) {
                vector<int> host = classRegions[cur[0].first];
                for (size_t i = 1; i < cur.size() && !host.empty(); i++)
                    host = interClasses(host, classRegions[cur[i].first]);
                if (host.empty() || host[0] != curRid) return;   // canonical home
                Pat P; P.host = host; P.comp = cur;
                for (auto &cm : cur) P.classSet.push_back(cm.first);
                sort(P.classSet.begin(), P.classSet.end());
                long long mu = 1;
                for (auto &cm : cur) mu *= (long long)llround(C(classSize[cm.first], cm.second));
                P.mult = mu;
                pats.push_back(std::move(P));
                return;
            }
            for (int i = idx; i < (int)cls.size(); i++) {
                int c = cls[i], maxj = min(rem, classSize[c]);
                for (int j = 1; j <= maxj; j++) {
                    cur.push_back({c, j}); enumAll(i+1, cls, rem-j); cur.pop_back();
                }
            }
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
    printf("[sct] quotient nC=%d  base-leaves=%zu  build=%.2fs\n",
           nC, baseLeaves.size(), secs(Tqg0, Tqg1));
    fflush(stdout);

    // m as a GLOBAL class vector (length nC): m_vec[c] = mult of class c in P.
    auto compToVec = [&](const vector<pair<int,int>> &comp) -> Vec {
        Vec v((size_t)nC, 0);
        for (auto &cm : comp) v[(size_t)cm.first] = (int16_t)cm.second;
        return v;
    };

    // region-IE init support, kept ONLY as the G2a cross-check reference.
    for (auto &P : pats) { P.sup = suppOf(P); }

    // sanity: SCT total s-cliques == region-IE total (each s-clique has
    // C(s,r) r-subcliques). Strong global gate before per-pattern work.
    {
        double sclSCT = 0;
        for (auto &lf : baseLeaves) { Vec zb((size_t)lf.m(), 0);  // compact: m()=|classIds|
            sclSCT += ccpath::support_count(lf, zb, ccpath_ncr); }
        double sclIE = 0; for (auto &P : pats) sclIE += (double)P.mult * P.sup;
        sclIE /= C(s, r);
        printf("[sct] total s-cliques: class-SCT=%.0f  region-IE=%.0f  %s\n",
               sclSCT, sclIE, fabs(sclSCT - sclIE) < 0.5 ? "[OK]" : "[MISMATCH]");
        fflush(stdout);
    }

    // Step 3: compaction + pattern<->leaf maps.
    int nLeaf = (int)baseLeaves.size();
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
    std::unordered_map<std::string,int> patIdx; patIdx.reserve(pats.size() * 2);
    auto compKey = [](const vector<pair<int,int>> &comp) {
        std::string k; k.reserve(comp.size() * 8);
        for (auto &cm : comp) { k.append((const char*)&cm.first, 4); k.append((const char*)&cm.second, 4); }
        return k;
    };
    for (int pi = 0; pi < (int)pats.size(); pi++) patIdx[compKey(pats[pi].comp)] = pi;
    {
        vector<int> lcs, lcap; vector<pair<int,int>> cur;
        std::function<void(int,int,int)> enumLP = [&](int lid, int idx, int rem) {
            if (rem == 0) {
                auto it = patIdx.find(compKey(cur));
                if (it == patIdx.end()) return;           // not a registered pattern
                int pi = it->second;
                // confirm host on the compact leaf (filters m with no s-extension)
                if (ccpath::support_count(slotPaths[lid][0], compToLocal(cur, lid), ccpath_ncr) > 0.0) {
                    patLeaves[pi].push_back(lid); leafPats[lid].push_back(pi);
                }
                return;
            }
            for (int i = idx; i < (int)lcs.size(); i++) {
                int mx = std::min(rem, lcap[i]);
                for (int j = 1; j <= mx; j++) { cur.push_back({lcs[i], j}); enumLP(lid, i+1, rem-j); cur.pop_back(); }
            }
        };
        for (int lid = 0; lid < nLeaf; lid++) {
            const CCPath &cp = slotPaths[lid][0];          // compact leaf
            lcs = supC[lid];                               // global class ids (sorted)
            lcap.assign(cp.u.begin(), cp.u.end());         // local u, parallel to supC
            cur.clear(); enumLP(lid, 0, r);
        }
    }
    for (auto &v : leafPats) std::sort(v.begin(), v.end());
    // pre-mapped compact b for every (pattern, hosting-leaf) pair.
    vector<vector<Vec>> pbLocal(pats.size());         // parallel to patLeaves[pi]
    for (int pi = 0; pi < (int)pats.size(); pi++) {
        pbLocal[pi].reserve(patLeaves[pi].size());
        for (int lid : patLeaves[pi]) pbLocal[pi].push_back(compToLocal(pats[pi].comp, lid));
    }
    // inverse: for each leaf, the local b of every pattern it hosts (parallel
    // to leafPats[lid]). Lets the incremental peel evaluate SC(slot, m_Q) for
    // every Q on a changed slot without a global remap.
    vector<vector<Vec>> leafPatLocB(nLeaf);           // parallel to leafPats[lid]
    for (int lid = 0; lid < nLeaf; lid++) {
        leafPatLocB[lid].reserve(leafPats[lid].size());
        for (int qi : leafPats[lid]) leafPatLocB[lid].push_back(compToLocal(pats[qi].comp, lid));
    }
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
    // support(pi) = sum over hosting slots of sum over slot's paths of
    // support_count(path, b_local). Uses the pre-mapped compact b.
    auto sctSupport = [&](int pi) -> double {
        double tot = 0.0;
        const auto &bl = pbLocal[pi];
        const auto &ls = patLeaves[pi];
        for (size_t k = 0; k < ls.size(); k++)
            for (auto &p : slotPaths[ls[k]]) tot += ccpath::support_count(p, bl[k], ccpath_ncr);
        return tot;
    };
    auto T5 = Clock::now();
    printf("[sct] pattern<->leaf maps + compaction=%.2fs\n", secs(Tqg1, T5));
    fflush(stdout);

    // -------------------- GATE G2a (init support equality) --------------------
    // For EVERY pattern: SCT sum-over-leaves support == region-IE suppOf.
    {
        int okc = 0, badc = 0; double worst = 0;
        for (int pi = 0; pi < (int)pats.size(); pi++) {
            double sIE = pats[pi].sup;                 // already = suppOf init
            double sSCT = sctSupport(pi);
            if (fabs(sIE - sSCT) < 0.5) okc++;
            else {
                badc++; worst = max(worst, fabs(sIE - sSCT));
                if (badc <= 8) {
                    printf("[G2a] MISMATCH pi=%d host=%zu comp=[", pi, pats[pi].host.size());
                    for (auto &cm : pats[pi].comp) printf("(c%d:%d)", cm.first, cm.second);
                    printf("] regionIE=%.1f SCT=%.1f leaves=%zu\n",
                           sIE, sSCT, patLeaves[pi].size());
                }
            }
        }
        printf("[G2a] %d/%d patterns: SCT-sum == region-IE  %s  (worst|d|=%.1f)\n",
               okc, okc + badc, badc == 0 ? "[OK]" : "[FAIL]", worst);
        fflush(stdout);
        if (badc != 0) {
            printf("[G2a] FAILED — aborting before peel (correctness gate).\n");
            return 3;
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
    // hurts, so a small KMAX (cheap per-call IE) wins: KMAX=2 finishes
    // ca-GrQc(3,4) in <1s and dblp-core30(3,4) in 0.04s vs minutes at KMAX=12.
    // Correctness is invariant to KMAX (verified 80/80 vs brute at KMAX=1,2,12).
    int KMAX = 2;
    if (getenv("SCT_KMAX")) { KMAX = atoi(getenv("SCT_KMAX")); if (KMAX < 1) KMAX = 1; }
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
    size_t maxSplit = 0;                              // diagnostic: largest split-set
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
    vector<CCPath> sfdKeep;                            // reused scratch
    auto slotForbidDiff = [&](int lid, const Vec &bloc,
                              const vector<pair<int,int>> &plNZ, vector<CCPath> &chgOld) {
        vector<CCPath> &cur = slotPaths[lid];
        vector<CCPath> &keep = sfdKeep; keep.clear(); keep.reserve(cur.size());
        chgOld.clear();
        for (auto &p : cur) {
            bool imposs = false;                          // impossible(p, bloc)?
            for (auto &pv : plNZ) if ((int)p.u[pv.first] < pv.second) { imposs = true; break; }
            if (imposs) { keep.push_back(std::move(p)); continue; }
            chgOld.push_back(p);                            // snapshot before change
            if (ccpath::covers_whole_path(p, bloc)) continue; // path fully dead (a==bloc)
            ccpath::insert_antichain(p.forbidden, bloc);
            if ((int)p.forbidden.size() > KMAX) {
                auto kids = ccpath::controlled_split(p, KMAX);
                for (auto &k : kids) keep.push_back(std::move(k));
            } else keep.push_back(std::move(p));
        }
        cur.swap(keep);
        if (cur.size() > maxSplit) maxSplit = cur.size();
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
    auto scWithTerms = [&](const CCPath &p,
                           const vector<pair<Vec,int>> &terms,
                           const Vec &b, const Vec *addLow) -> double {
        // forbidden early-out is only valid when no addLow raises the floor; with
        // addLow we must let the IE terms account for the forbidden set instead.
        if (!addLow) {
            for (const auto &a : p.forbidden) if (ccpath::leq(a, b)) return 0.0;
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
    auto ensureLeafMap = [&](int lid) {
        if (leafQbuilt[lid]) return;
        leafQbuilt[lid] = 1;
        auto &mp = leafQ2pat[lid];
        const auto &qb = leafPatLocB[lid];
        mp.reserve(qb.size() * 2);
        for (int t = 0; t < (int)qb.size(); t++) mp[hashVec(qb[t])].push_back(t);
    };

    long long npat = (long long)pats.size(), peeledN = 0, maxKey = 0;
    for (auto &P : pats) { P.key = (long long)llround(P.sup); maxKey = max(maxKey, P.key); }
    unordered_map<long long, vector<int>> bk;
    for (int pi = 0; pi < (int)pats.size(); pi++) bk[pats[pi].key].push_back(pi);
    map<double,double> coreDist;
    long long curLevel = 0;
    vector<char> seen(pats.size(), 0);
    vector<double> delta(pats.size(), 0.0);          // per-affected exact drop
    Vec uEnv, sufPl, qcand;                           // reused per-leaf scratch
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
        const auto &pleaf = patLeaves[pi];
        const auto &ploc  = pbLocal[pi];
        vector<CCPath> chgOld;                         // pre-insertion snapshots
        vector<vector<pair<Vec,int>>> chgOldTerms;     // cached IE terms (pre-insert)
        vector<pair<int,int>> plNZ;                    // sparse nonzeros of m_P local
        for (size_t k = 0; k < pleaf.size(); k++) {
            int lid = pleaf[k];
            if (slotPaths[lid].empty()) continue;      // leaf fully peeled: no witnesses
            const Vec &pl = ploc[k];                   // m_P local to lid (== a_p, h=0)
            int Mloc = (int)pl.size();
            // sparse support of m_P (positions where it is nonzero) -- the only
            // positions the impossible / feasibility tests depend on.
            plNZ.clear();
            for (int c = 0; c < Mloc; c++) if (pl[c]) plNZ.push_back({c, (int)pl[c]});
            // Record P (updates the stored slot via split) and capture the CHANGED
            // OLD paths (the pre-insertion snapshots where P's threshold applies).
            slotForbidDiff(lid, pl, plNZ, chgOld);
            if (chgOld.empty()) continue;              // P touched nothing here
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
            const auto &qbAll = leafPatLocB[lid];
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
                if (skipH1 && pats[qi].host.size() == 1) return;  // peels at L_M regardless
                double d = 0.0;                         // drop, via delta formula
                for (size_t z = 0; z < chgOld.size(); z++) {
                    const CCPath &p = chgOld[z];
                    // a_p = pl (P's local threshold; h=0). drop_p = count(y: alive
                    // in A_p, y >= max(ql, pl)) weighted by C(n - ql, y - ql). The
                    // weight base STAYS ql (addLow=pl only raises the floor) --
                    // reweighting-correct delta. Feasible only if max(ql,pl)<=p.u.
                    bool ok = true;
                    for (int c = 0; c < Mloc; c++) {
                        int v = (int)ql[c] > (int)pl[c] ? (int)ql[c] : (int)pl[c];
                        if (v > (int)p.u[c]) { ok = false; break; }
                    }
                    if (ok) d += scWithTerms(p, chgOldTerms[z], ql, &pl);
                }
                if (d == 0.0) return;
                if (!seen[qi]) { seen[qi] = 1; aff.push_back(qi); }
                delta[qi] += d;
            };
            // look up a complete candidate (by rolling hash) and confirm exactly.
            auto applyCand = [&](const Vec &ql, uint64_t h) {
                auto it = q2p.find(h);
                if (it == q2p.end()) return;
                for (int t : it->second) if (qbAll[t] == ql) { applyIdx(ql, t); return; }
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
            dfs(dfs, 0, r, 0, 1469598103934665603ULL);
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
    auto T6 = Clock::now();
    printf("[sct-peel] peel=%.2fs  peeled=%lld/%lld  maxSplit(split-set)=%zu\n",
           secs(T5,T6), peeledN, npat, maxSplit);
    printf("[sct-peel] TIMING MCE=%.2f enum=%.2f sct-build+maps=%.2f peel=%.2f total=%.2f\n",
           secs(T1,T2), secs(T3,T4), secs(Tqg0,T5), secs(T5,T6), secs(T1,T6));
    // fold in the r-mergeable direct-assigned cores
    for (auto &kv : directCoreDist) coreDist[kv.first] += kv.second;
    double maxCore = 0; for (auto &kv : coreDist) maxCore = max(maxCore, kv.first);
    printf("[sct-peel] Max core: %.0f\n", maxCore);
    for (auto &kv : coreDist) printf("core=%.0f count=%.0f\n", kv.first, kv.second);
    return 0;
}
