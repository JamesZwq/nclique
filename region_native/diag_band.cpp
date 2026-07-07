// diag_band.cpp -- §132 t=1 DIAGONAL BAND ENGINE (prototype, separate from the production engine).
//
// Computes the (r, r+1) nucleus cell WITHOUT enumerating the full pattern space (the class-CPI
// pattern wall that makes band rows infeasible, §131). Three exact ingredients:
//
//  (1) HOST-1 CLOSED FORM (SKIP_H1 theorem, t=1): an r-clique contained in exactly one region M
//      has core = |M| - r. NEVER ENUMERATED: per region, count = C(|M|, r) minus the multi-host
//      r-cliques hosted there; the whole complement goes to the distribution in one line.
//  (2) MULTI-HOST ENUMERATION FROM PAIRWISE INTERSECTIONS: an r-clique is multi-host iff it fits
//      inside A^B for some region pair -- so candidates cost Sum_pairs(|A^B|>=r) multisets, and
//      pairs sharing >= r vertices get RARER as r grows: the lever is strongest inside the band.
//      sup0(Q) = |union of hosting regions| - r (exact at t=1; every extension vertex lies in a
//      hosting region).
//  (3) EXACT EVENT REPLAY of the interaction:
//      - NO-EARLY-DEATH LEMMA: no r-clique inside M dies before level l_M = |M| - r (every one
//        has >= l_M witnesses inside M, all alive by induction). So at l_M every witness in M is
//        alive and the REGION WAVE needs no already-dead subtraction.
//      - OWNERSHIP SPLIT: a witness W with >= 1 host-1 sub is hosted ONLY in that sub's region M
//        (any maximal clique containing W contains the sub), so it dies exactly at l_M -- owned
//        by M's wave. A witness whose EVERY r-sub is multi-host dies via multi-host deaths --
//        owned by the deadW-deduped per-death enumeration. The two classes are disjoint, so every
//        witness is killed and credited exactly once.
//      Region wave at level l_M credits each alive multi-host Q hosted in M with the closed form
//        drop = Sum_{c in classes(M)} (n_c - Q_c) * [W = Q + e_c has >= 1 host-1 r-sub].
//
// Output: the same "core=X count=Y" distribution as region_native_sct_peel (same MCE floor
// s = r+1, so the gate is a straight diff). Guards: DIAG_MAX_PATS pattern cap (clean exit 7).
//
// Build: g++ -O3 -std=c++17 -o diag_band diag_band.cpp
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>

using namespace std;
using Clock = chrono::high_resolution_clock;
static double secs(Clock::time_point a, Clock::time_point b) {
    return chrono::duration_cast<chrono::duration<double>>(b - a).count();
}

// ----- graph (CSR, undirected, 0-indexed; header "n m") -- copied from the engine -----
struct Graph {
    int n = 0;
    vector<int> off, adj;
    bool adjacent(int u, int v) const {
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
        es.push_back({u, v}); es.push_back({v, u});
        maxid = max(maxid, max(u, v));
    }
    fclose(f);
    Graph g; g.n = maxid + 1;
    vector<int> cnt(g.n + 1, 0);
    for (auto &e : es) cnt[e.first]++;
    g.off.assign(g.n + 1, 0);
    for (int i = 0; i < g.n; i++) g.off[i + 1] = g.off[i] + cnt[i];
    g.adj.resize(g.off[g.n]);
    vector<int> cur(g.off.begin(), g.off.end() - 1);
    for (auto &e : es) g.adj[cur[e.first]++] = e.second;
    vector<int> noff(g.n + 1, 0), nadj;
    nadj.reserve(g.adj.size());
    for (int i = 0; i < g.n; i++) {
        int *b = &g.adj[g.off[i]], *e = &g.adj[g.off[i + 1]];
        sort(b, e);
        int *ne = unique(b, e);
        noff[i + 1] = noff[i] + (int)(ne - b);
        for (int *p = b; p < ne; p++) nadj.push_back(*p);
    }
    g.off = std::move(noff); g.adj = std::move(nadj);
    return g;
}

// ----- Bron-Kerbosch (degeneracy-seeded) maximal cliques >= minSz -- copied from the engine -----
struct MCE {
    const Graph &g;
    int minSz;
    double budget;
    Clock::time_point t0;
    bool aborted = false;
    vector<vector<int>> cliques;
    MCE(const Graph &g_, int m, double b) : g(g_), minSz(m), budget(b) {}
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
        int bestu = -1, bestc = -1;
        auto consider = [&](const vector<int> &S) {
            for (int u : S) {
                int c = (int)intersect_nbr(P, g, u).size();
                if (c > bestc) { bestc = c; bestu = u; }
            }
        };
        consider(P); consider(X);
        vector<int> cand; cand.reserve(P.size());
        {
            const int *b = &g.adj[g.off[bestu]], *e = &g.adj[g.off[bestu + 1]];
            const int *p = b;
            for (int x : P) {
                while (p < e && *p < x) p++;
                if (!(p < e && *p == x)) cand.push_back(x);
            }
        }
        for (int v : cand) {
            R.push_back(v);
            rec(R, intersect_nbr(P, g, v), intersect_nbr(X, g, v));
            R.pop_back();
            P.erase(lower_bound(P.begin(), P.end(), v));
            X.insert(lower_bound(X.begin(), X.end(), v), v);
            if (aborted) return;
        }
    }
    struct NbrRange { const int *b, *e; const int* begin() const{return b;} const int* end() const{return e;} };
    NbrRange g_nbr(int v) const { return {&g.adj[g.off[v]], &g.adj[g.off[v + 1]]}; }
    bool run() {
        t0 = Clock::now();
        int n = g.n;
        vector<int> deg(n), order, pos(n, -1);
        for (int v = 0; v < n; v++) deg[v] = g.deg(v);
        int maxd = 0; for (int v = 0; v < n; v++) maxd = max(maxd, deg[v]);
        vector<vector<int>> bucket(maxd + 1);
        for (int v = 0; v < n; v++) bucket[deg[v]].push_back(v);
        vector<char> removed(n, 0);
        order.reserve(n);
        for (int processed = 0; processed < n; processed++) {
            int d2 = 0;
            while (d2 <= maxd && bucket[d2].empty()) d2++;
            int v = bucket[d2].back(); bucket[d2].pop_back();
            if (removed[v]) { processed--; continue; }
            removed[v] = 1; pos[v] = (int)order.size(); order.push_back(v);
            for (int w : g_nbr(v)) if (!removed[w]) { deg[w]--; if (deg[w] < 0) deg[w] = 0; bucket[deg[w]].push_back(w); }
        }
        vector<int> R;
        for (int idx = 0; idx < n; idx++) {
            int v = order[idx];
            vector<int> P, X;
            P.reserve(g.deg(v)); X.reserve(g.deg(v));
            for (int w : g_nbr(v)) { if (pos[w] > idx) P.push_back(w); else X.push_back(w); }
            sort(P.begin(), P.end()); sort(X.begin(), X.end());
            R.assign(1, v);
            rec(R, std::move(P), std::move(X));
            if (aborted) return false;
        }
        return !aborted;
    }
};

// ----- nCr as double (matches the engine's double-based counts) -----
static vector<vector<double>> NCR;
static void build_ncr(int N) {
    NCR.assign(N + 1, vector<double>(N + 1, 0.0));
    for (int i = 0; i <= N; i++) { NCR[i][0] = 1.0; for (int j = 1; j <= i; j++) NCR[i][j] = NCR[i-1][j-1] + NCR[i-1][j]; }
}
static inline double C(int n, int k) { if (k < 0 || n < 0 || k > n) return 0.0; return NCR[n][k]; }

// order-free composition fingerprint (same construction as the engine's mixCV/hashVecInc)
static inline uint64_t mixCV(int c, int v) {
    uint64_t x = ((uint64_t)(uint32_t)c << 20) ^ (uint64_t)(uint32_t)v;
    x *= 0x9E3779B97F4A7C15ULL; x ^= x >> 29; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 32;
    return x;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s graph.edges r [--mce-budget S]\n", argv[0]); return 1; }
    const char *gpath = argv[1];
    int r = atoi(argv[2]);
    int s = r + 1;
    double mceBudget = 300.0;
    for (int i = 3; i < argc; i++)
        if (!strcmp(argv[i], "--mce-budget") && i + 1 < argc) mceBudget = atof(argv[++i]);
    long long maxPats = getenv("DIAG_MAX_PATS") ? atoll(getenv("DIAG_MAX_PATS")) : 50000000LL;

    auto T0 = Clock::now();
    Graph g = load_graph(gpath);
    auto T1 = Clock::now();
    printf("[db] graph n=%d m=%ld  (r,s)=(%d,%d)  load=%.2fs\n", g.n, (long)g.adj.size() / 2, r, s, secs(T0, T1));

    MCE mce(g, s, mceBudget);
    if (!mce.run()) { printf("[db] MCE exceeded budget; abort.\n"); return 2; }
    auto &regions = mce.cliques;
    for (auto &R : regions) sort(R.begin(), R.end());
    int nR = (int)regions.size();
    auto T2 = Clock::now();
    printf("[db] regions(>=s)=%d  MCE=%.2fs\n", nR, secs(T1, T2));
    if (regions.empty()) { printf("[db] no region >= s.\n"); return 0; }
    int maxMC = 0; for (auto &R : regions) maxMC = max(maxMC, (int)R.size());
    build_ncr(maxMC + 2);

    // ----- classes (identical region-membership profiles) -- same construction as the engine -----
    vector<vector<int>> vtxR(g.n);
    for (int i = 0; i < nR; i++) for (int v : regions[i]) vtxR[v].push_back(i);
    for (int v = 0; v < g.n; v++) sort(vtxR[v].begin(), vtxR[v].end());
    unordered_map<string,int> profKey; profKey.reserve(g.n);
    vector<int> classOf(g.n, -1), classSize;
    vector<vector<int>> classRegions;                       // class -> sorted region ids
    for (int v = 0; v < g.n; v++) {
        if (vtxR[v].empty()) continue;
        string k; k.reserve(vtxR[v].size() * 4);
        for (int x : vtxR[v]) k.append((const char *)&x, 4);
        auto it = profKey.find(k);
        int c;
        if (it == profKey.end()) { c = (int)classRegions.size(); profKey[k] = c; classRegions.push_back(vtxR[v]); classSize.push_back(0); }
        else c = it->second;
        classOf[v] = c; classSize[c]++;
    }
    int nC = (int)classRegions.size();
    vector<vector<int>> regionClasses(nR);                  // region -> sorted class ids
    for (int c = 0; c < nC; c++) for (int rid : classRegions[c]) regionClasses[rid].push_back(c);
    for (int i = 0; i < nR; i++) sort(regionClasses[i].begin(), regionClasses[i].end());
    auto T3 = Clock::now();
    printf("[db] classes=%d  build=%.2fs\n", nC, secs(T2, T3));

    // ----- pair discovery: region pairs sharing >= r vertices (the only sources of multi-host) -----
    vector<pair<int,int>> pairs;
    {
        vector<int> cnt(nR, 0); vector<int> dirty; dirty.reserve(256);
        for (int M = 0; M < nR; M++) {
            for (int v : regions[M]) for (int o : vtxR[v]) if (o > M) { if (cnt[o] == 0) dirty.push_back(o); cnt[o]++; }
            for (int o : dirty) { if (cnt[o] >= r) pairs.push_back({M, o}); cnt[o] = 0; }
            dirty.clear();
        }
    }
    auto T4 = Clock::now();
    printf("[db] pairs(overlap>=r)=%zu  pair-scan=%.2fs\n", pairs.size(), secs(T3, T4));

    // ----- multi-host pattern enumeration over pair intersections (deduped) -----
    struct Pat {
        vector<pair<int,int>> comp;   // (classId, mult), sorted, sum = r
        vector<int> host;             // sorted hosting region ids (>= 2)
        vector<int> uniCls;           // sorted union of hosting regions' classes (extension space)
        double mult = 1, sup = 0;     // #r-cliques in the orbit; current support
        double core = -1; long long key = -1; bool alive = true;
    };
    vector<Pat> pats;
    unordered_map<uint64_t, vector<int>> patIdx;            // fingerprint -> candidate ids (exact compare)
    auto fpComp = [&](const vector<pair<int,int>> &comp) {
        uint64_t h = 0; for (auto &cm : comp) h ^= mixCV(cm.first, cm.second); return h;
    };
    auto lookup = [&](const vector<pair<int,int>> &comp) -> int {
        auto it = patIdx.find(fpComp(comp));
        if (it == patIdx.end()) return -1;
        for (int id : it->second) if (pats[id].comp == comp) return id;
        return -1;
    };
    {
        vector<int> cc; vector<pair<int,int>> cur;
        bool capped = false;
        for (auto &pr : pairs) {
            // common classes of the pair (class fully inside both regions)
            cc.clear();
            const auto &ca = regionClasses[pr.first];
            for (int c : ca)
                if (std::binary_search(classRegions[c].begin(), classRegions[c].end(), pr.second)) cc.push_back(c);
            // enumerate r-multisets over cc (bounded by class sizes); every one is multi-host by construction
            auto rec = [&](auto &&self, int idx, int rem) -> void {
                if (capped) return;
                if (rem == 0) {
                    if (lookup(cur) >= 0) return;                       // already found via another pair
                    if ((long long)pats.size() >= maxPats) { capped = true; return; }
                    Pat P; P.comp = cur;
                    // host = regions containing every comp class (sorted-list intersection)
                    P.host = classRegions[cur[0].first];
                    vector<int> tmp;
                    for (size_t k = 1; k < cur.size() && !P.host.empty(); k++) {
                        const auto &b = classRegions[cur[k].first];
                        tmp.clear();
                        std::set_intersection(P.host.begin(), P.host.end(), b.begin(), b.end(), back_inserter(tmp));
                        P.host.swap(tmp);
                    }
                    double mu = 1; for (auto &cm : cur) mu *= C(classSize[cm.first], cm.second);
                    P.mult = mu;
                    patIdx[fpComp(cur)].push_back((int)pats.size());
                    pats.push_back(std::move(P));
                    return;
                }
                for (int i = idx; i < (int)cc.size(); i++) {
                    int c = cc[i], mx = min(rem, classSize[c]);
                    for (int j = 1; j <= mx; j++) { cur.push_back({c, j}); self(self, i + 1, rem - j); cur.pop_back(); }
                }
            };
            rec(rec, 0, r);
            if (capped) break;
        }
        if (capped) { printf("[db] PATTERN-CAP exceeded (%lld) -- clean abort\n", maxPats); return 7; }
    }
    // union classes + sup0 = |union of hosts| - r; and per-region hosted lists for waves/aggregates
    vector<vector<int>> regPats(nR);
    {
        vector<char> mark(nC, 0);
        for (int pi = 0; pi < (int)pats.size(); pi++) {
            Pat &P = pats[pi];
            long long uni = 0;
            for (int M : P.host) {
                regPats[M].push_back(pi);
                for (int c : regionClasses[M]) if (!mark[c]) { mark[c] = 1; P.uniCls.push_back(c); uni += classSize[c]; }
            }
            sort(P.uniCls.begin(), P.uniCls.end());
            for (int c : P.uniCls) mark[c] = 0;
            P.sup = (double)uni - r;
            P.key = (long long)llround(P.sup);
        }
    }
    long long nMH = (long long)pats.size();
    double mhRC = 0; for (auto &P : pats) mhRC += P.mult;
    auto T5 = Clock::now();
    printf("[db] multi-host patterns=%lld (r-cliques=%.0f)  enum=%.2fs\n", nMH, mhRC, secs(T4, T5));

    // ----- host-1 closed-form complement per region (never enumerated) -----
    map<double,double> dist;
    {
        for (int M = 0; M < nR; M++) {
            double hosted = 0;
            for (int pi : regPats[M]) hosted += pats[pi].mult;
            double h1 = C((int)regions[M].size(), r) - hosted;
            if (h1 > 0) dist[(double)((int)regions[M].size() - r)] += h1;
        }
    }

    // ----- the multi-host peel with region-wave replay -----
    // multiHost(comp) = fits >= 2 regions (sorted-list intersection with early exit at 2).
    // MEMOIZED by the order-free composition fingerprint: sub-compositions repeat massively
    // across waves and witness enumerations (same precedent as the engine's fingerprint-only
    // deadY set; the bit-exact gates verify the aliasing risk is not realized).
    vector<int> mhA, mhB;
    unordered_map<uint64_t, char> memoMH, memoASM;
    auto multiHost = [&](const vector<pair<int,int>> &comp) -> bool {
        if (comp.empty()) return false;
        uint64_t h = 0; for (auto &cm : comp) h ^= mixCV(cm.first, cm.second);
        auto mi = memoMH.find(h);
        if (mi != memoMH.end()) return mi->second;
        mhA = classRegions[comp[0].first];
        for (size_t k = 1; k < comp.size() && mhA.size() >= 2; k++) {
            const auto &b = classRegions[comp[k].first];
            mhB.clear();
            std::set_intersection(mhA.begin(), mhA.end(), b.begin(), b.end(), back_inserter(mhB));
            mhA.swap(mhB);
        }
        bool res = mhA.size() >= 2;
        memoMH[h] = res;
        return res;
    };
    // allSubsMultiHost(W): every r-sub composition W - e_d is multi-host (else W is wave-owned).
    // Memoized by W's fingerprint.
    vector<pair<int,int>> Wc, Sc;
    auto allSubsMultiHost = [&](const vector<pair<int,int>> &W) -> bool {
        uint64_t h = 0; for (auto &cm : W) h ^= mixCV(cm.first, cm.second);
        auto mi = memoASM.find(h);
        if (mi != memoASM.end()) return mi->second;
        bool res = true;
        for (size_t d = 0; d < W.size() && res; d++) {
            Sc.clear();
            for (size_t k = 0; k < W.size(); k++) {
                int m = W[k].second - (k == d ? 1 : 0);
                if (m > 0) Sc.push_back({W[k].first, m});
            }
            if (!multiHost(Sc)) res = false;
        }
        memoASM[h] = res;
        return res;
    };
    // waves: level -> regions (only those hosting >= 1 multi-host pattern)
    unordered_map<long long, vector<int>> waves;
    long long maxKey = 0;
    for (int M = 0; M < nR; M++)
        if (!regPats[M].empty()) {
            long long lv = (long long)regions[M].size() - r;
            waves[lv].push_back(M);
            maxKey = max(maxKey, lv);
        }
    unordered_map<long long, vector<int>> bk;
    for (int pi = 0; pi < (int)pats.size(); pi++) { bk[pats[pi].key].push_back(pi); maxKey = max(maxKey, pats[pi].key); }
    unordered_set<uint64_t> deadW;                          // witnesses killed by MULTI-HOST deaths only
    vector<double> delta(pats.size(), 0.0);
    vector<char> seen(pats.size(), 0);
    vector<int> aff;
    long long peeled = 0, curLevel = 0, lastWave = -1;
    auto applyAff = [&]() {
        for (int qi : aff) {
            seen[qi] = 0;
            double ns = pats[qi].sup - delta[qi];
            delta[qi] = 0.0;
            long long nk = (long long)llround(ns);
            if (nk < curLevel) nk = curLevel;               // monotone clamp (same as the engine)
            if (nk != pats[qi].key) { pats[qi].sup = ns; pats[qi].key = nk; bk[nk].push_back(qi); }
        }
        aff.clear();
    };
    while (peeled < nMH && curLevel <= maxKey + 1) {
        if (lastWave != curLevel) {                          // region waves open the level
            lastWave = curLevel;
            auto wit = waves.find(curLevel);
            if (wit != waves.end()) {
                for (int M : wit->second) {
                    long long lv = (long long)regions[M].size() - r;
                    for (int qi : regPats[M]) {
                        Pat &Q = pats[qi];
                        if (!Q.alive) continue;
                        double drop = 0;
                        for (int c : regionClasses[M]) {
                            int qc = 0;
                            for (auto &cm : Q.comp) if (cm.first == c) { qc = cm.second; break; }
                            int avail = classSize[c] - qc;
                            if (avail <= 0) continue;
                            // W = Q + e_c dies in this wave iff it has >= 1 host-1 sub
                            Wc.clear(); bool added = false;
                            for (auto &cm : Q.comp) {
                                if (cm.first == c) { Wc.push_back({c, cm.second + 1}); added = true; }
                                else Wc.push_back(cm);
                            }
                            if (!added) { Wc.push_back({c, 1}); sort(Wc.begin(), Wc.end()); }
                            if (!allSubsMultiHost(Wc)) drop += (double)avail;
                        }
                        if (drop > 0) {
                            if (!seen[qi]) { seen[qi] = 1; aff.push_back(qi); }
                            delta[qi] += drop;
                        }
                        (void)lv;
                    }
                }
                applyAff();
            }
        }
        auto it = bk.find(curLevel);
        if (it == bk.end() || it->second.empty()) { curLevel++; continue; }
        int pi = it->second.back(); it->second.pop_back();
        Pat &P = pats[pi];
        if (!P.alive || P.key != curLevel) continue;
        P.alive = false; P.core = (double)curLevel; peeled++;
        dist[P.core] += P.mult;
        // multi-host death: enumerate witnesses W = P + e_c over the host-union classes;
        // only all-subs-multi-host W are owned here (wave-owned ones are skipped), deadW dedupes.
        for (int c : P.uniCls) {
            int pc = 0;
            for (auto &cm : P.comp) if (cm.first == c) { pc = cm.second; break; }
            int avail = classSize[c] - pc;
            if (avail <= 0) continue;
            Wc.clear(); bool added = false;
            for (auto &cm : P.comp) {
                if (cm.first == c) { Wc.push_back({c, cm.second + 1}); added = true; }
                else Wc.push_back(cm);
            }
            if (!added) { Wc.push_back({c, 1}); sort(Wc.begin(), Wc.end()); }
            if (!allSubsMultiHost(Wc)) continue;             // wave-owned witness
            uint64_t h = 0; for (auto &cm : Wc) h ^= mixCV(cm.first, cm.second);
            if (!deadW.insert(h).second) continue;           // already killed by another multi-host death
            // credit every sub Q = W - e_d (all multi-host by the test above)
            for (size_t d = 0; d < Wc.size(); d++) {
                Sc.clear();
                int freed = Wc[d].first;
                for (size_t k = 0; k < Wc.size(); k++) {
                    int m = Wc[k].second - (k == d ? 1 : 0);
                    if (m > 0) Sc.push_back({Wc[k].first, m});
                }
                int qi = lookup(Sc);
                if (qi < 0 || qi == pi || !pats[qi].alive) continue;
                int qd = 0; for (auto &cm : Sc) if (cm.first == freed) { qd = cm.second; break; }
                double w = (double)(classSize[freed] - qd);  // vertex extensions of Q covered by W
                if (!seen[qi]) { seen[qi] = 1; aff.push_back(qi); }
                delta[qi] += w;
            }
        }
        applyAff();
    }
    auto T6 = Clock::now();
    printf("[db] peel=%.2fs  peeled=%lld/%lld multi-host  deadW=%zu\n", secs(T5, T6), peeled, nMH, deadW.size());
    double maxCore = 0; for (auto &kv : dist) maxCore = max(maxCore, kv.first);
    printf("[sct-peel] Max core: %.0f\n", maxCore);
    printf("[db] TIMING MCE=%.2f classes=%.2f pairs=%.2f enum=%.2f peel=%.2f total=%.2f\n",
           secs(T1, T2), secs(T2, T3), secs(T3, T4), secs(T4, T5), secs(T5, T6), secs(T1, T6));
    for (auto &kv : dist) printf("core=%.0f count=%.0f\n", kv.first, kv.second);
    return 0;
}
