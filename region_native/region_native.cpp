// region_native.cpp — standalone prototype of region-native pivot
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
    // recursive union count over a list of class-sets (already each >= s in vsize)
    long long visited = 0;
    // NB: take B BY VALUE. An earlier by-reference version aliased B and
    // unionCount(rest) emptied it before the intersection step ran, so
    // overlaps were never subtracted (over-count on every multi-host
    // tuple). Caught by the adversarial self-check (121/500, all the
    // |Host|=1 tuples, passed; multi-host over-counted).
    std::function<double(vector<Node>)> unionCount =
        [&](vector<Node> B) -> double {
        if (B.empty()) return 0.0;
        visited++;
        Node M = std::move(B.back()); B.pop_back();
        // `rest` is exactly the post-pop B; build intersections from it
        // BEFORE recursing so nothing is consumed out from under us.
        double here = C(M.vsize - r, s - r);
        // build {M ∩ rest_i}, prune size<s and dominated
        vector<Node> inter;
        inter.reserve(B.size());
        for (auto &N : B) {
            vector<int> cs = interClasses(M.classes, N.classes);
            int vs = classesSize(cs);
            if (vs >= s) inter.push_back({move(cs), vs});
        }
        // dominance prune: drop any inter node whose class-set is a subset
        // of another (subset => its s-cliques are contained => redundant
        // in the union). Sort by size desc, keep maximal.
        if (inter.size() > 1) {
            sort(inter.begin(), inter.end(),
                 [](const Node &a, const Node &b){ return a.classes.size() > b.classes.size(); });
            vector<Node> keep;
            for (auto &nd : inter) {
                bool dom = false;
                for (auto &k : keep) {
                    // nd.classes subset of k.classes ?
                    if (nd.classes.size() <= k.classes.size()) {
                        size_t i = 0, j = 0; bool sub = true;
                        while (i < nd.classes.size()) {
                            while (j < k.classes.size() && k.classes[j] < nd.classes[i]) j++;
                            if (j >= k.classes.size() || k.classes[j] != nd.classes[i]) { sub = false; break; }
                            i++; j++;
                        }
                        if (sub) { dom = true; break; }
                    }
                }
                if (!dom) keep.push_back(move(nd));
            }
            inter = move(keep);
        }
        double interU = unionCount(std::move(inter));
        double restU = unionCount(std::move(B));   // B still intact here
        return here + restU - interU;
    };

    // enumerate region tuples and compute support.
    // tuple = class-multiplicity vector of weight r over a region's classes.
    // dedup tuples across regions by canonical key. host set computed per
    // tuple as ∩ of profiles of its classes.
    long long nTuples = 0;
    double sumSup = 0;
    // for verify
    struct Sample { vector<pair<int,int>> cj; vector<int> rep; double sup; };
    vector<Sample> samples;
    std::mt19937 rng(20260618);

    // CANONICAL-HOME tuple dedup, no hash map: enumerate a tuple's classes
    // from region M and process it ONLY if M is the lowest-id region in
    // host(tau) = ∩ profiles. Each realized tuple is thus counted exactly
    // once, at its canonical home. This is the pivot canonical-home idea
    // native to the quotient, and it removes the per-tuple string-key map
    // (millions of allocations) the first version paid.
    int curRid = 0;
    vector<pair<int,int>> cur;  // (classId, j)
    vector<Node> Bbuf;          // reused host-node buffer
    std::function<void(int,const vector<int>&,int)> enumTuple =
        [&](int idx, const vector<int> &cls, int rem) {
        if (rem == 0) {
            // host set = ∩ profiles of classes in cur
            vector<int> host = classRegions[cur[0].first];
            for (size_t i = 1; i < cur.size() && !host.empty(); i++)
                host = interClasses(host, classRegions[cur[i].first]);
            if (host.empty() || host[0] != curRid) return;  // not canonical home
            // build B&B nodes from host regions (as class sets), prune size<s
            Bbuf.clear();
            for (int rid : host) {
                int vs = (int)regions[rid].size();
                if (vs >= s) Bbuf.push_back({regionClasses[rid], vs});
            }
            double sup;
            // fast paths: |host alive|<=2 needs no recursion
            if (Bbuf.size() == 1) {
                sup = C(Bbuf[0].vsize - r, s - r);
            } else if (Bbuf.size() == 2) {
                int iv = classesSize(interClasses(Bbuf[0].classes, Bbuf[1].classes));
                sup = C(Bbuf[0].vsize - r, s - r) + C(Bbuf[1].vsize - r, s - r)
                      - (iv >= s ? C(iv - r, s - r) : 0.0);
            } else {
                sup = unionCount(Bbuf);
            }
            nTuples++; sumSup += sup;
            if (verifyN && (int)samples.size() < verifyN) {
                Sample sm; sm.cj = cur; sm.sup = sup;
                samples.push_back(move(sm));
            }
            return;
        }
        for (int i = idx; i < (int)cls.size(); i++) {
            int c = cls[i];
            int maxj = min(rem, classSize[c]);
            for (int j = 1; j <= maxj; j++) {
                cur.push_back({c, j});
                enumTuple(i + 1, cls, rem - j);
                cur.pop_back();
            }
        }
    };
    for (int i = 0; i < nR; i++) { curRid = i; enumTuple(0, regionClasses[i], r); }
    auto T4 = Clock::now();
    printf("[rn] tuples=%lld  support B&B visited=%lld  support=%.2fs\n",
           nTuples, visited, secs(T3, T4));
    printf("[rn] TIMING: MCE=%.2fs support=%.2fs (support-only, the CPI-replacement phase)\n",
           secs(T1, T2), secs(T3, T4));

    // ---- adversarial correctness: sup(tau) vs direct enumeration ----
    if (verifyN) {
        int okc = 0, bad = 0;
        for (auto &sm : samples) {
            // materialize a representative r-clique R0 of this tuple:
            // pick j vertices from each class. Need actual vertex ids:
            // build class -> member vertices lazily.
            static vector<vector<int>> classMembers;
            if (classMembers.empty()) {
                classMembers.assign(nC, {});
                for (int v = 0; v < g.n; v++) if (classOf[v] >= 0) classMembers[classOf[v]].push_back(v);
            }
            vector<int> R0;
            bool feasible = true;
            for (auto &pr : sm.cj) {
                auto &mem = classMembers[pr.first];
                if ((int)mem.size() < pr.second) { feasible = false; break; }
                for (int t = 0; t < pr.second; t++) R0.push_back(mem[t]);
            }
            if (!feasible) continue;
            // verify R0 is a clique
            bool isclique = true;
            for (size_t a = 0; a < R0.size() && isclique; a++)
                for (size_t b = a + 1; b < R0.size(); b++)
                    if (!g.adjacent(R0[a], R0[b])) { isclique = false; break; }
            if (!isclique) { printf("[rn] verify: R0 not a clique?!\n"); bad++; continue; }
            // ground truth: count (s-r)-cliques in common neighborhood of R0
            vector<int> common(g.adj.begin() + g.off[R0[0]], g.adj.begin() + g.off[R0[0] + 1]);
            for (size_t a = 1; a < R0.size(); a++) {
                vector<int> nb(g.adj.begin() + g.off[R0[a]], g.adj.begin() + g.off[R0[a] + 1]);
                vector<int> out; size_t i = 0, j = 0;
                while (i < common.size() && j < nb.size()) {
                    if (common[i] < nb[j]) i++;
                    else if (common[i] > nb[j]) j++;
                    else { out.push_back(common[i]); i++; j++; }
                }
                common = move(out);
            }
            // remove R0 vertices from common
            {
                vector<int> out;
                for (int x : common) if (find(R0.begin(), R0.end(), x) == R0.end()) out.push_back(x);
                common = move(out);
            }
            int need = s - r;
            long long gt = 0;
            vector<int> chosen;
            std::function<void(int)> rec2 = [&](int start) {
                if ((int)chosen.size() == need) { gt++; return; }
                for (int i = start; i < (int)common.size(); i++) {
                    int v = common[i];
                    bool good = true;
                    for (int c : chosen) if (!g.adjacent(v, c)) { good = false; break; }
                    if (good) { chosen.push_back(v); rec2(i + 1); chosen.pop_back(); }
                }
            };
            rec2(0);
            if (fabs((double)gt - sm.sup) < 0.5) okc++;
            else {
                bad++;
                if (bad <= 5) {
                    printf("[rn] MISMATCH tuple sup=%.0f truth=%lld  R0=", sm.sup, gt);
                    for (int x : R0) printf("%d ", x); printf("\n");
                }
            }
        }
        printf("[rn] VERIFY: %d/%d region-native == direct enumeration%s\n",
               okc, okc + bad, (bad == 0 ? "  [EXACT]" : "  [MISMATCH!]"));
    }
    return 0;
}
