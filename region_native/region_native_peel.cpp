// region_native_peel.cpp — region-native END-TO-END (r,s)-nucleus peel.
// Front half (load/MCE/classes) shared with region_native.cpp. Support is
// computed via region-IE where each region-intersection leaf returns
// ccpath::support_count (alive witnesses with the peeled-tuple forbidden
// antichain). Peel = bucket queue at PATTERN granularity. Output = the
// core-number distribution, to be matched bit-for-bit against V3LM.
#include <unordered_set>
#include <map>
#include "CCPathCore.h"
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
    auto T4 = Clock::now();
    long long totalRCliques = 0; for (auto &P : pats) totalRCliques += P.mult;
    printf("[rn-peel] patterns=%zu  r-cliques=%lld  enum=%.2fs\n",
           pats.size(), totalRCliques, secs(T3, T4));
    fflush(stdout);
    if (pats.empty()) { printf("[rn-peel] no patterns.\n"); return 0; }

    // ---- initial support (empty antichains -> matches verified counting) ----
    for (auto &P : pats) { P.sup = suppOf(P); P.key = (long long)llround(P.sup); }
    auto T5 = Clock::now();
    printf("[rn-peel] initial-support=%.2fs\n", secs(T4, T5));
    fflush(stdout);

    // region -> patterns hosted there (affected-set lookup on peel)
    vector<vector<int>> regToPats(nR);
    for (int pi = 0; pi < (int)pats.size(); pi++)
        for (int rid : pats[pi].host) regToPats[rid].push_back(pi);

    // per-region peeled antichain (Pareto-min comp vectors over the region's
    // class order). Maintained incrementally; used for the single-region
    // fast delta. peeled[] (global) is kept only for the rare |host|>=2
    // fallback recompute via suppOf.
    vector<vector<Vec>> regAnti(nR);
    auto mapComp = [&](const vector<pair<int,int>> &comp, int rid) {
        const vector<int> &Cs = regionClasses[rid];
        Vec v((size_t)Cs.size(), 0);
        for (auto &cm : comp) {
            int pos = (int)(lower_bound(Cs.begin(), Cs.end(), cm.first) - Cs.begin());
            v[(size_t)pos] = (int16_t)cm.second;
        }
        return v;
    };
    // alive witnesses >= compVec in region rid, honoring regAnti[rid].
    auto regSupp = [&](int rid, const Vec &b) -> double {
        const vector<int> &Cs = regionClasses[rid];
        int M = (int)Cs.size();
        Vec nn((size_t)M, 0), zeros((size_t)M, 0);
        for (int i = 0; i < M; i++) nn[i] = (int16_t)classSize[Cs[i]];
        CCPath p = CCPath::initial(zeros, nn, s);
        p.forbidden = regAnti[rid];
        if (p.forbidden.size() > maxForb) maxForb = p.forbidden.size();
        return ccpath::support_count(p, b, ccpath_ncr);
    };

    // ---- peel: bucket queue, INCREMENTAL delta ----
    long long npat = (long long)pats.size(), peeledN = 0, maxKey = 0;
    for (auto &P : pats) maxKey = max(maxKey, P.key);
    unordered_map<long long, vector<int>> bk;
    for (int pi = 0; pi < (int)pats.size(); pi++) bk[pats[pi].key].push_back(pi);
    map<double,double> coreDist;
    long long curLevel = 0, fallbacks = 0;
    Vec maxv;                              // reused max(Q*,Q) buffer
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
        if ((peeledN & 0xFF) == 0)
            fprintf(stderr, "[peel] %lld/%lld lvl=%lld fb=%lld maxForb=%zu t=%.1fs\n",
                    peeledN, npat, curLevel, fallbacks, maxForb, secs(T5, Clock::now()));
        coreDist[P.core] += (double)P.mult;

        // collect support changes (qi -> new sup) before re-bucketing
        // SINGLE-REGION FAST PATH: every witness of P lives in its one
        // region M, so the delta to any affected Q is one regSupp(M,max).
        if (P.host.size() == 1) {
            int M = P.host[0];
            Vec pm = mapComp(P.comp, M);
            const vector<int> &Cs = regionClasses[M];
            for (int qi : regToPats[M]) {
                Pat &Q = pats[qi];
                if (qi == pi || !Q.alive) continue;
                // max(P,Q) over M's classes; feasible iff sum<=s & <=classSize
                Vec qm = mapComp(Q.comp, M);
                maxv.assign(Cs.size(), 0);
                int tot = 0; bool feas = true;
                for (size_t k = 0; k < Cs.size(); k++) {
                    int mx = max((int)pm[k], (int)qm[k]);
                    if (mx > classSize[Cs[k]]) { feas = false; break; }
                    maxv[k] = (int16_t)mx; tot += mx;
                }
                if (!feas || tot > s) continue;            // cannot co-occur
                double d = regSupp(M, maxv);               // antichain WITHOUT P yet
                if (d <= 0) continue;
                double ns = Q.sup - d;
                long long nk = (long long)llround(ns);
                if (nk < curLevel) nk = curLevel;
                if (nk != Q.key) { Q.sup = ns; Q.key = nk; bk[nk].push_back(qi); }
            }
            ccpath::insert_antichain(regAnti[M], pm);      // now record P
        } else {
            // RARE |host|>=2: record P everywhere, then full recompute of the
            // affected (shares a host region) via suppOf (global `peeled`).
            fallbacks++;
            peeled.push_back(pi);
            for (int rid : P.host) ccpath::insert_antichain(regAnti[rid], mapComp(P.comp, rid));
            // dedup affected
            static vector<char> seen2; if ((int)seen2.size()<(int)pats.size()) seen2.assign(pats.size(),0);
            vector<int> aff;
            for (int rid : P.host) for (int qi : regToPats[rid]) {
                if (qi==pi || !pats[qi].alive || seen2[qi]) continue;
                seen2[qi]=1; aff.push_back(qi);
            }
            for (int qi : aff) {
                seen2[qi]=0;
                double ns = suppOf(pats[qi]);
                long long nk = (long long)llround(ns);
                if (nk < curLevel) nk = curLevel;
                if (nk != pats[qi].key) { pats[qi].sup=ns; pats[qi].key=nk; bk[nk].push_back(qi); }
            }
        }
    }
    auto T6 = Clock::now();
    printf("[rn-peel] peel=%.2fs  peeled=%lld/%lld  multihost-fallbacks=%lld\n",
           secs(T5,T6), peeledN, npat, fallbacks);
    printf("[rn-peel] TIMING MCE=%.2f enum=%.2f initsup=%.2f peel=%.2f total=%.2f\n",
           secs(T1,T2), secs(T3,T4), secs(T4,T5), secs(T5,T6), secs(T1,T6));
    double maxCore = 0; for (auto &kv : coreDist) maxCore = max(maxCore, kv.first);
    printf("[rn-peel] Max core: %.0f\n", maxCore);
    for (auto &kv : coreDist) printf("core=%.0f count=%.0f\n", kv.first, kv.second);
    return 0;
}
