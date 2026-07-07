// nsi_query.cpp -- §136 NSI query tool: load a serialized Nucleus Spectrum Index (SCT_INDEX_OUT
// of region_native_sct_peel, format "NSI1") and answer exact spectrum queries.
//
// QUERY SEMANTICS for an r-clique R (caller guarantees R is a clique) and cell s in [s0, smax]:
//  1. Map vertices to classes. If every vertex is classed, form the class composition and look
//     it up in the pattern table. On a hit, walk the chain from the stored cold-cell core:
//     at each step, certified (kappa_prev equals the clique value at the previous cell) implies
//     the closed form C(cP - r, s - r) for ALL later cells (T3, absorbing: return immediately);
//     otherwise the cell's residue dictionary holds the exact peeled core.
//  2. On a miss (some vertex unclassed, or the composition is no pattern), R lies in no ACTIVE
//     region; if some MERGEABLE region M contains all of R, kappa = C(|M| - r, s - r) (T6:
//     mergeable regions are isolated cliques). The pattern-hit and mergeable cases are
//     structurally disjoint: a composition that matches a pattern places R inside an active
//     region (classes are wholly contained in regions).
//  3. Otherwise R is in no clique of size >= s0, so kappa_{r,s}(R) = 0 for every s >= s0.
//
// Modes:
//   ./nsi_query idx.nsi stats
//   ./nsi_query idx.nsi point  s  v1 v2 ... vr
//   ./nsi_query idx.nsi spectrum  v1 v2 ... vr
//   ./nsi_query idx.nsi count  s  k          (aggregate: # r-cliques with kappa_{r,s} >= k)
//   ./nsi_query idx.nsi bench  queries.txt   (one query per line: r vertex ids; reports latency)
//
// Build: g++ -O3 -std=c++17 -o nsi_query nsi_query.cpp
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using Clock = chrono::high_resolution_clock;

static vector<vector<double>> NCR;
static void build_ncr(int N) {
    NCR.assign(N + 1, vector<double>(N + 1, 0.0));
    for (int i = 0; i <= N; i++) { NCR[i][0] = 1.0; for (int j = 1; j <= i; j++) NCR[i][j] = NCR[i-1][j-1] + NCR[i-1][j]; }
}
static inline double C(int n, int k) { if (k < 0 || n < 0 || k > n) return 0.0; return NCR[n][k]; }

static inline uint64_t mixCV(int c, int v) {
    uint64_t x = ((uint64_t)(uint32_t)c << 20) ^ (uint64_t)(uint32_t)v;
    x *= 0x9E3779B97F4A7C15ULL; x ^= x >> 29; x *= 0xBF58476D1CE4E5B9ULL; x ^= x >> 32;
    return x;
}

struct NSI {
    int r = 0, s0 = 0, smax = 0, n = 0, nC = 0, nPats = 0, nMerg = 0;
    vector<int32_t> classOf;                       // vertex -> class (-1 = none)
    vector<int32_t> mergOff, mergAdj, mergSize;    // vertex -> mergeable-region CSR + |M|
    vector<pair<int32_t,uint8_t>> compFlat;        // pattern comps, flattened
    vector<int32_t> compOff;                       // pattern -> comp span
    vector<int32_t> cP;                            // pattern -> clique bound base c(P)
    vector<double>  k0;                            // pattern -> cold-cell core kappa_{s0}
    unordered_map<uint64_t, vector<int32_t>> patIdx;             // comp fingerprint -> ids
    vector<vector<pair<int32_t,double>>> resid;    // per cell s0+1..smax, SORTED by pid
    vector<vector<pair<double,double>>> dists;     // per cell s0..smax, (core, count)
    long long bytes = 0;
};

static bool loadNSI(const char *path, NSI &x) {
    FILE *f = fopen(path, "rb");
    if (!f) return false;
    char magic[4];
    if (fread(magic, 4, 1, f) != 1 || memcmp(magic, "NSI1", 4)) { fclose(f); return false; }
    auto r32 = [&]() { int32_t v; if (fread(&v, 4, 1, f) != 1) exit(3); return v; };
    auto r64 = [&]() { int64_t v; if (fread(&v, 8, 1, f) != 1) exit(3); return v; };
    auto rd  = [&]() { double v;  if (fread(&v, 8, 1, f) != 1) exit(3); return v; };
    auto r8  = [&]() { uint8_t v; if (fread(&v, 1, 1, f) != 1) exit(3); return v; };
    x.r = r32(); x.s0 = r32(); x.smax = r32(); x.n = r32(); x.nC = r32(); x.nPats = r32(); x.nMerg = r32();
    x.classOf.resize(x.n);
    for (int v = 0; v < x.n; v++) x.classOf[v] = r32();
    // mergeable regions -> per-vertex CSR
    {
        vector<vector<int32_t>> byV(x.n);
        x.mergSize.resize(x.nMerg);
        for (int m = 0; m < x.nMerg; m++) {
            int sz = r32(); x.mergSize[m] = sz;
            for (int i = 0; i < sz; i++) { int v = r32(); if (v >= 0 && v < x.n) byV[v].push_back(m); }
        }
        x.mergOff.assign(x.n + 1, 0);
        for (int v = 0; v < x.n; v++) x.mergOff[v + 1] = x.mergOff[v] + (int32_t)byV[v].size();
        x.mergAdj.resize(x.mergOff[x.n]);
        for (int v = 0; v < x.n; v++) copy(byV[v].begin(), byV[v].end(), x.mergAdj.begin() + x.mergOff[v]);
    }
    x.compOff.assign(x.nPats + 1, 0);
    x.cP.resize(x.nPats); x.k0.resize(x.nPats);
    x.compFlat.reserve((size_t)x.nPats * 2);
    x.patIdx.reserve((size_t)x.nPats * 2);
    for (int pi = 0; pi < x.nPats; pi++) {
        int len = r8();
        uint64_t h = 0;
        for (int i = 0; i < len; i++) { int32_t c = r32(); uint8_t m = r8(); x.compFlat.push_back({c, m}); h ^= mixCV(c, m); }
        x.compOff[pi + 1] = (int32_t)x.compFlat.size();
        x.cP[pi] = r32(); x.k0[pi] = rd();
        x.patIdx[h].push_back(pi);
    }
    x.resid.resize(x.smax - x.s0);
    for (auto &rl : x.resid) {
        int64_t cnt = r64(); rl.resize(cnt);
        for (auto &pc : rl) { pc.first = r32(); pc.second = rd(); }
        sort(rl.begin(), rl.end());                              // pid-sorted for binary search
    }
    x.dists.resize(x.smax - x.s0 + 1);
    for (auto &d : x.dists) {
        int cnt = r32(); d.resize(cnt);
        for (auto &kv : d) { kv.first = rd(); kv.second = rd(); }
    }
    x.bytes = ftell(f);
    fclose(f);
    int mx = x.s0 + 2;
    for (int m = 0; m < x.nMerg; m++) mx = max(mx, (int)x.mergSize[m]);
    for (int pi = 0; pi < x.nPats; pi++) mx = max(mx, (int)x.cP[pi]);
    build_ncr(mx + 2);
    return true;
}

// exact kappa_{r,s}(R); scratch vectors reused across calls (bench-friendly)
struct Query {
    const NSI &x;
    vector<pair<int32_t,uint8_t>> comp;
    vector<int32_t> mg, mg2;
    explicit Query(const NSI &x_) : x(x_) {}

    int lookupPattern(const int *vs, int r) {
        comp.clear();
        for (int i = 0; i < r; i++) {
            int32_t c = (vs[i] >= 0 && vs[i] < x.n) ? x.classOf[vs[i]] : -1;
            if (c < 0) return -1;
            comp.push_back({c, 1});
        }
        sort(comp.begin(), comp.end());
        int w = 0;                                              // merge duplicates -> (class, mult)
        for (int i = 0; i < (int)comp.size(); i++) {
            if (w && comp[i].first == comp[w-1].first) comp[w-1].second++;
            else comp[w++] = comp[i];
        }
        comp.resize(w);
        uint64_t h = 0; for (auto &cm : comp) h ^= mixCV(cm.first, cm.second);
        auto it = x.patIdx.find(h);
        if (it == x.patIdx.end()) return -1;
        for (int pi : it->second) {
            int a = x.compOff[pi], b = x.compOff[pi + 1];
            if (b - a != w) continue;
            bool eq = true;
            for (int i = 0; i < w; i++) if (x.compFlat[a + i] != comp[i]) { eq = false; break; }
            if (eq) return pi;
        }
        return -1;
    }
    int lookupMergeable(const int *vs, int r) {                 // region containing ALL of R, or -1
        auto span = [&](int v) { return make_pair(x.mergAdj.begin() + x.mergOff[v], x.mergAdj.begin() + x.mergOff[v + 1]); };
        if (vs[0] < 0 || vs[0] >= x.n) return -1;
        auto s0 = span(vs[0]); mg.assign(s0.first, s0.second);
        for (int i = 1; i < r && !mg.empty(); i++) {
            if (vs[i] < 0 || vs[i] >= x.n) return -1;
            auto si = span(vs[i]); mg2.clear();
            set_intersection(mg.begin(), mg.end(), si.first, si.second, back_inserter(mg2));
            mg.swap(mg2);
        }
        return mg.empty() ? -1 : mg[0];
    }
    double kappa(const int *vs, int s) {
        int pi = lookupPattern(vs, x.r);
        if (pi >= 0) {
            int c = x.cP[pi];
            double kprev = x.k0[pi];
            if (s == x.s0) return kprev;
            for (int ss = x.s0 + 1; ss <= s; ss++) {
                if (kprev == C(c - x.r, ss - 1 - x.r)) return C(c - x.r, s - x.r);   // T3, absorbing
                const auto &rl = x.resid[ss - x.s0 - 1];
                auto it = lower_bound(rl.begin(), rl.end(), make_pair((int32_t)pi, -1.0));
                if (it == rl.end() || it->first != pi) return -2.0;                  // data error
                kprev = it->second;
            }
            return kprev;
        }
        int m = lookupMergeable(vs, x.r);
        if (m >= 0) return C((int)x.mergSize[m] - x.r, s - x.r);
        return 0.0;
    }
};

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s idx.nsi stats|point|spectrum|count|bench ...\n", argv[0]); return 1; }
    NSI x;
    auto TL0 = Clock::now();
    if (!loadNSI(argv[1], x)) { fprintf(stderr, "cannot load %s\n", argv[1]); return 1; }
    double tload = chrono::duration_cast<chrono::duration<double>>(Clock::now() - TL0).count();
    const string mode = argv[2];
    Query q(x);

    if (mode == "stats") {
        long long residTot = 0; for (auto &rl : x.resid) residTot += (long long)rl.size();
        double totRC = 0;                                   // r-cliques covered (from the s0 distribution)
        if (!x.dists.empty()) for (auto &kv : x.dists[0]) totRC += kv.second;
        printf("index=%s  %.2f MB  load=%.2fs\n", argv[1], x.bytes / 1048576.0, tload);
        printf("r=%d cells s=%d..%d  n=%d classes=%d\n", x.r, x.s0, x.smax, x.n, x.nC);
        printf("patterns=%d  merg-regions=%d  residue-entries=%lld (over %d cells)\n",
               x.nPats, x.nMerg, residTot, x.smax - x.s0);
        printf("r-cliques covered at s0: %.0f  -> bytes per r-clique = %.4f (the compression claim)\n",
               totRC, totRC > 0 ? x.bytes / totRC : 0.0);
        return 0;
    }
    if (mode == "point" && argc >= 4 + x.r) {
        int s = atoi(argv[3]);
        vector<int> vs(x.r);
        for (int i = 0; i < x.r; i++) vs[i] = atoi(argv[4 + i]);
        printf("kappa_{%d,%d} = %.0f\n", x.r, s, q.kappa(vs.data(), s));
        return 0;
    }
    if (mode == "spectrum" && argc >= 3 + x.r) {
        vector<int> vs(x.r);
        for (int i = 0; i < x.r; i++) vs[i] = atoi(argv[3 + i]);
        for (int s = x.s0; s <= x.smax; s++) printf("s=%d kappa=%.0f\n", s, q.kappa(vs.data(), s));
        return 0;
    }
    if (mode == "count" && argc >= 5) {
        int s = atoi(argv[3]); double k = atof(argv[4]);
        if (s < x.s0 || s > x.smax) { fprintf(stderr, "s out of range\n"); return 1; }
        double tot = 0;
        for (auto &kv : x.dists[s - x.s0]) if (kv.first >= k) tot += kv.second;
        printf("count(kappa_{%d,%d} >= %.0f) = %.0f\n", x.r, s, k, tot);
        return 0;
    }
    if (mode == "pointfile" && argc >= 5) {                     // batch point queries at one cell (gate mode)
        int s = atoi(argv[3]);
        FILE *qf = fopen(argv[4], "r");
        if (!qf) { fprintf(stderr, "cannot open %s\n", argv[4]); return 1; }
        vector<int> flat; int v;
        while (fscanf(qf, "%d", &v) == 1) flat.push_back(v);
        fclose(qf);
        long long nq = (long long)flat.size() / x.r;
        string out; out.reserve((size_t)nq * 8);
        char buf[32];
        for (long long i = 0; i < nq; i++) {
            snprintf(buf, sizeof buf, "%.0f\n", q.kappa(&flat[i * x.r], s));
            out += buf;
        }
        fwrite(out.data(), 1, out.size(), stdout);
        return 0;
    }
    if (mode == "bench" && argc >= 4) {
        FILE *qf = fopen(argv[3], "r");
        if (!qf) { fprintf(stderr, "cannot open %s\n", argv[3]); return 1; }
        vector<int> flat; int v;
        while (fscanf(qf, "%d", &v) == 1) flat.push_back(v);
        fclose(qf);
        long long nq = (long long)flat.size() / x.r;
        if (!nq) { fprintf(stderr, "no queries\n"); return 1; }
        double sink = 0;
        auto T0 = Clock::now();
        for (long long i = 0; i < nq; i++) sink += q.kappa(&flat[i * x.r], x.smax);   // point @ deepest cell
        double tp = chrono::duration_cast<chrono::duration<double>>(Clock::now() - T0).count();
        auto T1 = Clock::now();
        for (long long i = 0; i < nq; i++)
            for (int s = x.s0; s <= x.smax; s++) sink += q.kappa(&flat[i * x.r], s);  // full spectrum
        double ts = chrono::duration_cast<chrono::duration<double>>(Clock::now() - T1).count();
        printf("bench: queries=%lld  point(s=%d)=%.0f ns/query  spectrum(%d cells)=%.0f ns/query  (sink=%.0f)\n",
               nq, x.smax, 1e9 * tp / nq, x.smax - x.s0 + 1, 1e9 * ts / nq, sink);
        return 0;
    }
    fprintf(stderr, "bad mode/args\n");
    return 1;
}
