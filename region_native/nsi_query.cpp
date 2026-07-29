// nsi_query.cpp -- serialized Nucleus Spectrum Index query tool.  It loads both
// legacy fixed-r NSI1 and the additive multi-r plane format NSI2 emitted by
// SCT_RSWEEP + SCT_INDEX_OUT, and answers exact point/row queries.
//
// LEGACY NSI1 QUERY SEMANTICS for an r-clique R (caller guarantees R is a clique) and cell s in [s0, smax]:
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
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <functional>
#include <sstream>
#include <string>
#include <sys/resource.h>
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

static int mainNSI1(int argc, char **argv) {
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

// -----------------------------------------------------------------------------
// NSI2: multi-r whole-plane index.  NSI1 above intentionally remains isolated
// so that its byte format, CLI, and query semantics stay backward compatible.

using SteadyClock = chrono::steady_clock;

struct LEReader {
    FILE *f = nullptr;
    int64_t pos = 0, size = 0, limit = 0;
    string error;

    explicit LEReader(const char *path) {
        f = fopen(path, "rb");
        if (!f) { error = string("cannot open ") + path; return; }
        if (fseek(f, 0, SEEK_END) != 0) { error = "cannot seek to end"; return; }
        long z = ftell(f);
        if (z < 0) { error = "cannot determine file size"; return; }
        size = limit = (int64_t)z;
        if (fseek(f, 0, SEEK_SET) != 0) { error = "cannot rewind file"; return; }
    }
    ~LEReader() { if (f) fclose(f); }
    bool ok() const { return f && error.empty(); }
    void fail(const string &s) { if (error.empty()) error = s; }
    bool bytes(void *dst, size_t n) {
        if (!ok()) return false;
        if ((uint64_t)n > (uint64_t)(limit - pos)) {
            fail("truncated block at byte " + to_string(pos));
            return false;
        }
        if (n && fread(dst, 1, n, f) != n) {
            fail("read failure at byte " + to_string(pos));
            return false;
        }
        pos += (int64_t)n;
        return true;
    }
    bool skip(int64_t n) {
        if (!ok()) return false;
        if (n < 0 || n > limit - pos) {
            fail("invalid skip at byte " + to_string(pos));
            return false;
        }
        if (n && fseek(f, (long)n, SEEK_CUR) != 0) {
            fail("seek failure at byte " + to_string(pos));
            return false;
        }
        pos += n;
        return true;
    }
    bool seek(int64_t p, int64_t lim) {
        if (!ok()) return false;
        if (p < 0 || lim < p || lim > size) { fail("invalid block bounds"); return false; }
        if (fseek(f, (long)p, SEEK_SET) != 0) { fail("seek failure"); return false; }
        pos = p; limit = lim;
        return true;
    }
    uint8_t u8() { uint8_t v = 0; bytes(&v, 1); return v; }
    uint16_t u16() {
        uint8_t b[2] = {0,0}; bytes(b, 2);
        return (uint16_t)b[0] | ((uint16_t)b[1] << 8);
    }
    int16_t i16() { return (int16_t)u16(); }
    uint32_t u32() {
        uint8_t b[4] = {0,0,0,0}; bytes(b, 4);
        return (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
               ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
    }
    int32_t i32() { return (int32_t)u32(); }
    uint64_t u64() {
        uint8_t b[8] = {0,0,0,0,0,0,0,0}; bytes(b, 8);
        uint64_t v = 0;
        for (int i = 0; i < 8; ++i) v |= (uint64_t)b[i] << (8 * i);
        return v;
    }
    int64_t i64() { return (int64_t)u64(); }
    double f64() {
        uint64_t bits = u64(); double v = 0;
        memcpy(&v, &bits, sizeof(v));
        return v;
    }
};

struct IntVecHash {
    size_t operator()(const vector<int32_t> &v) const noexcept {
        uint64_t h = 1469598103934665603ULL;
        for (int32_t x : v) h = (h ^ ((uint64_t)(uint32_t)x + 1)) * 1099511628211ULL;
        return (size_t)h;
    }
};

using Comp2 = vector<pair<int32_t,int16_t>>;

struct NSI2Pattern {
    Comp2 comp;
    int64_t mult = 0;
    int32_t cP = 0;
    double boundaryCore = -1;
    int32_t closedFrom = -1;
};

struct NSI2Column {
    int32_t r = 0, boundary = 0;
    vector<int32_t> mergeableRegions;
    vector<uint8_t> mergeMask;
    vector<int32_t> classRep;                 // universal class -> coarsened active representative
    vector<NSI2Pattern> patterns;
    vector<int32_t> origPid;                       // §222 NSI3: exception -> its original pattern id
    unordered_map<uint64_t, vector<int32_t>> patIdx;
    vector<vector<pair<int32_t,double>>> resid;
    vector<vector<pair<double,double>>> dists;
    int64_t fileOffset = 0, fileBytes = 0;
};

struct NSI2 {
    int32_t rmin = 0, rmax = 0, smin = 0, smax = 0;
    int32_t n = 0, nC = 0, nR = 0, nLeaf = 0, nCols = 0;
    vector<int32_t> classOf, classSize, regionSize;
    vector<vector<int32_t>> classRegions;
    vector<vector<int32_t>> classRegionsBySize;    // §222: same sets, ordered by decreasing |M|
    vector<int64_t> directory;
    bool slim = false;                             // §222: NSI3 -- only EXCEPTION patterns are stored
    vector<NSI2Column> columns;
    vector<int32_t> rToColumn;
    int64_t fileBytes = 0, sharedBytes = 0;
};

// NSI2 only needs k <= Smax-rmin.  A flat, banded Pascal table preserves the
// engine's exact recurrence/rounding while avoiding the NSI1 O(max-clique^2)
// table on a large graph with a small indexed plane.
static int NCR2N = -1, NCR2K = -1;
static vector<double> NCR2;
static bool build_ncr2(int N, int K) {
    if (N < 0 || K < 0) return false;
    const size_t stride = (size_t)K + 1;
    if ((size_t)N + 1 > numeric_limits<size_t>::max() / stride) return false;
    try { NCR2.assign(((size_t)N + 1) * stride, 0.0); }
    catch (const bad_alloc &) { return false; }
    NCR2N = N; NCR2K = K;
    for (int n = 0; n <= N; ++n) {
        NCR2[(size_t)n * stride] = 1.0;
        for (int k = 1; k <= min(n, K); ++k)
            NCR2[(size_t)n * stride + k] =
                NCR2[(size_t)(n - 1) * stride + k - 1] +
                NCR2[(size_t)(n - 1) * stride + k];
    }
    return true;
}
static inline double C2(int n, int k) {
    if (n < 0 || k < 0 || k > n || n > NCR2N || k > NCR2K) return 0.0;
    return NCR2[(size_t)n * ((size_t)NCR2K + 1) + k];
}

static bool saneCount(int64_t count, int64_t unit, int64_t remaining) {
    return count >= 0 && unit >= 0 && count <= 1000000000LL &&
           (unit == 0 || count <= remaining / unit);
}

static bool buildColumnRepresentatives(const NSI2 &x, NSI2Column &col, string &error) {
    col.mergeMask.assign((size_t)x.nR, 0);
    for (int32_t rid : col.mergeableRegions) {
        if (rid < 0 || rid >= x.nR) { error = "mergeable region id out of range"; return false; }
        if (col.mergeMask[rid]) { error = "duplicate mergeable region id"; return false; }
        if (x.regionSize[rid] < col.boundary) {
            error = "mergeable region is smaller than the column boundary"; return false;
        }
        col.mergeMask[rid] = 1;
    }

    col.classRep.assign((size_t)x.nC, -1);
    unordered_map<vector<int32_t>,int32_t,IntVecHash> byProfile;
    byProfile.reserve((size_t)x.nC * 2 + 1);
    vector<int32_t> active;
    for (int32_t c = 0; c < x.nC; ++c) {
        active.clear();
        for (int32_t rid : x.classRegions[c])
            if (x.regionSize[rid] >= col.boundary && !col.mergeMask[rid])
                active.push_back(rid);
        if (active.empty()) continue;
        auto it = byProfile.find(active);
        if (it == byProfile.end()) {
            byProfile.emplace(active, c);
            col.classRep[c] = c;
        } else {
            col.classRep[c] = it->second;
        }
    }
    return true;
}

// §222 NSI3 loader: the slim plane index. Reuses the NSI2 in-memory structures so the whole query
// layer above is shared; the differences are that region sizes are RECONSTRUCTED from the class
// profiles (the NSI2 loader already recomputed them as a cross-check, so this is the same arithmetic
// promoted to the source of truth) and that only EXCEPTION patterns exist in a column.
// §227 NSI4 loader: the PACKED slim index. It fills the SAME in-memory NSI2 structures as loadNSI3,
// so the entire query layer above is untouched -- NSI4 is an encoding of NSI3, not a second design.
// The whole file is read into memory once and decoded from a cursor: varints are not seekable, and
// a stream of one-byte reads through stdio costs far more than the file itself.
struct PackReader {
    const uint8_t *p = nullptr, *end = nullptr;
    bool bad = false;
    uint64_t var() {
        uint64_t v = 0; int sh = 0;
        while (p < end) {
            uint8_t b = *p++;
            v |= (uint64_t)(b & 127) << sh;
            if (!(b & 128)) return v;
            sh += 7;
            if (sh > 63) break;
        }
        bad = true; return 0;
    }
    int64_t zig() { uint64_t u = var(); return (int64_t)(u >> 1) ^ -(int64_t)(u & 1); }
    uint8_t u8() { if (p >= end) { bad = true; return 0; } return *p++; }
    uint64_t u64() {
        if (end - p < 8) { bad = true; return 0; }
        uint64_t v = 0; for (int i = 0; i < 8; ++i) v |= (uint64_t)p[i] << (8 * i);
        p += 8; return v;
    }
    double f64() { uint64_t u = u64(); double d; memcpy(&d, &u, 8); return d; }
};

static bool loadNSI4(const char *path, NSI2 &x, string &error) {
    FILE *f = fopen(path, "rb");
    if (!f) { error = string("cannot open ") + path; return false; }
    fseek(f, 0, SEEK_END); long z = ftell(f); fseek(f, 0, SEEK_SET);
    if (z < 8) { fclose(f); error = "file too short"; return false; }
    vector<uint8_t> raw((size_t)z);
    if (fread(raw.data(), 1, (size_t)z, f) != (size_t)z) { fclose(f); error = "short read"; return false; }
    fclose(f);
    if (memcmp(raw.data(), "NSI4", 4)) { error = "not an NSI4 file"; return false; }
    x = NSI2{};
    x.fileBytes = z; x.slim = true;
    PackReader rd{raw.data() + 4, raw.data() + z};
    x.rmin = (int32_t)rd.var(); x.rmax = (int32_t)rd.var();
    x.smin = (int32_t)rd.var(); x.smax = (int32_t)rd.var();
    x.n = (int32_t)rd.var(); x.nC = (int32_t)rd.var(); x.nR = (int32_t)rd.var();
    x.nCols = (int32_t)rd.var();
    const bool allInt = rd.u8() != 0;
    x.nLeaf = 0;
    if (rd.bad || x.rmin < 1 || x.rmax < x.rmin || x.smin != x.rmin + 1 || x.smax < x.rmax + 1 ||
        x.n < 0 || x.nC < 0 || x.nR < 0 || x.nCols != x.rmax - x.rmin + 1 || x.nC > x.n) {
        error = "invalid NSI4 header"; return false;
    }
    auto num = [&]() { return allInt ? (double)rd.zig() : rd.f64(); };

    const size_t presBytes = ((size_t)x.n + 7) / 8;
    if ((size_t)(rd.end - rd.p) < presBytes) { error = "truncated presence bitmap"; return false; }
    const uint8_t *pres = rd.p; rd.p += presBytes;
    const uint64_t live = rd.var();
    int labelBits = 1; while ((1LL << labelBits) < (long long)x.nC + 1) ++labelBits;
    const size_t labelBytes = (size_t)((live * (uint64_t)labelBits + 7) / 8);
    if (rd.bad || (size_t)(rd.end - rd.p) < labelBytes) { error = "truncated class labels"; return false; }
    const uint8_t *lab = rd.p; rd.p += labelBytes;
    x.classOf.assign((size_t)x.n, -1);
    {
        uint64_t bitPos = 0, seen = 0;
        const uint64_t mask = labelBits >= 64 ? ~0ULL : ((1ULL << labelBits) - 1);
        for (int32_t v = 0; v < x.n; ++v) {
            if (!(pres[(size_t)v >> 3] & (1u << (v & 7)))) continue;
            if (seen++ >= live) { error = "presence bitmap and label count disagree"; return false; }
            uint64_t acc = 0;
            const size_t byteAt = (size_t)(bitPos >> 3);
            const int shift = (int)(bitPos & 7);
            for (int i = 0; i < 9 && byteAt + i < labelBytes; ++i) {
                if (i * 8 - shift >= labelBits + 8) break;
                acc |= (uint64_t)lab[byteAt + i] << (8 * i);
            }
            const int32_t c = (int32_t)((acc >> shift) & mask);
            if (c < 0 || c >= x.nC) { error = "class label out of range"; return false; }
            x.classOf[v] = c;
            bitPos += (uint64_t)labelBits;
        }
        if (seen != live) { error = "presence bitmap and label count disagree"; return false; }
    }

    x.classSize.resize((size_t)x.nC);
    x.classRegions.resize((size_t)x.nC);
    x.regionSize.assign((size_t)x.nR, 0);
    for (int32_t c = 0; c < x.nC; ++c) {
        const int32_t w = (int32_t)rd.var();
        const uint64_t m = rd.var();
        if (rd.bad || w < 0 || m > (uint64_t)x.nR) { error = "invalid class record"; return false; }
        x.classSize[c] = w;
        x.classRegions[c].resize((size_t)m);
        int32_t prev = -1;
        for (uint64_t i = 0; i < m; ++i) {
            const int32_t rid = (int32_t)(prev + rd.zig());
            if (rid < 0 || rid >= x.nR || rid <= prev) { error = "class profile not sorted/in range"; return false; }
            prev = rid;
            x.classRegions[c][i] = rid;
            x.regionSize[rid] += w;                        // |M| = Sum of its class weights, as NSI3
        }
    }
    if (rd.bad) { error = "truncated class block"; return false; }

    x.directory.resize((size_t)x.nCols);
    for (int32_t i = 0; i < x.nCols; ++i) x.directory[i] = (int64_t)rd.u64();
    if (rd.bad) { error = "truncated directory"; return false; }
    x.columns.resize((size_t)x.nCols);
    x.rToColumn.assign((size_t)(x.rmax - x.rmin + 1), -1);
    for (int32_t ci = 0; ci < x.nCols; ++ci) {
        const int64_t begin = x.directory[ci];
        const int64_t fin = (ci + 1 < x.nCols) ? x.directory[ci + 1] : (int64_t)z;
        if (begin < 0 || fin < begin || fin > (int64_t)z) { error = "invalid column bounds"; return false; }
        PackReader c{raw.data() + begin, raw.data() + fin};
        NSI2Column &col = x.columns[ci];
        col.fileOffset = begin; col.fileBytes = fin - begin;
        col.r = (int32_t)c.var(); col.boundary = (int32_t)c.var();
        if (c.bad || col.r < x.rmin || col.r > x.rmax || col.boundary != col.r + 1 ||
            x.rToColumn[col.r - x.rmin] != -1) { error = "invalid NSI4 column header"; return false; }
        x.rToColumn[col.r - x.rmin] = ci;
        const uint64_t nm = c.var();
        if (c.bad || nm > (uint64_t)x.nR) { error = "invalid mergeable list"; return false; }
        col.mergeableRegions.resize((size_t)nm);
        { int32_t prev = -1;
          for (uint64_t i = 0; i < nm; ++i) { prev = (int32_t)(prev + c.zig()); col.mergeableRegions[i] = prev; } }
        if (!buildColumnRepresentatives(x, col, error)) return false;
        const uint64_t np = c.var();
        if (c.bad) { error = "invalid exception count"; return false; }
        col.patterns.resize((size_t)np);
        col.origPid.resize((size_t)np);
        int32_t prevPid = -1;
        for (uint64_t pi = 0; pi < np; ++pi) {
            prevPid = (int32_t)(prevPid + c.zig());
            col.origPid[pi] = prevPid;
            const int len = c.u8();
            if (c.bad || len <= 0 || len > col.r) { error = "invalid exception comp"; return false; }
            NSI2Pattern &P = col.patterns[pi];
            P.comp.resize((size_t)len);
            int sum = 0;
            for (int i = 0; i < len; ++i) {
                const int32_t cc = (int32_t)c.var(); const int16_t b = (int16_t)c.var();
                if (cc < 0 || cc >= x.nC || b <= 0 || b > col.r) { error = "invalid exception comp entry"; return false; }
                P.comp[i] = {cc, b}; sum += b;
            }
            if (sum != col.r) { error = "exception comp does not sum to r"; return false; }
            P.boundaryCore = allInt ? (double)c.zig() : c.f64();
            P.closedFrom = (int32_t)c.zig();
            P.cP = -1;                                     // never stored; recovered from profiles
            uint64_t h = 0;
            for (auto &cm : P.comp) h ^= mixCV(cm.first, cm.second);
            col.patIdx[h].push_back((int32_t)pi);
        }
        const uint64_t nCells = c.var();
        if (c.bad) { error = "invalid residue cell count"; return false; }
        col.resid.resize((size_t)nCells);
        for (auto &rr : col.resid) {
            const uint64_t m = c.var();
            if (c.bad) { error = "invalid residue list"; return false; }
            rr.resize((size_t)m);
            int32_t prev = -1;
            for (auto &pc : rr) {
                prev = (int32_t)(prev + c.zig()); pc.first = prev;
                pc.second = allInt ? (double)c.zig() : c.f64();
            }
        }
        const uint64_t nd = c.var();
        if (c.bad) { error = "invalid dists count"; return false; }
        col.dists.resize((size_t)nd);
        for (auto &d : col.dists) {
            const uint64_t m = c.var();
            if (c.bad) { error = "invalid dist"; return false; }
            d.resize((size_t)m);
            for (auto &kv : d) {
                kv.first = allInt ? (double)c.zig() : c.f64();
                kv.second = allInt ? (double)c.zig() : c.f64();
            }
        }
        if (c.bad) { error = "truncated column"; return false; }
    }
    {   // same sizing rule as NSI3: from the reconstructed region sizes, which upper-bound cP
        int maxCombN = x.smax;
        for (int32_t v : x.regionSize) maxCombN = max(maxCombN, v);
        if (!build_ncr2(maxCombN, x.smax - x.rmin)) { error = "binomial table build failed"; return false; }
    }
    x.classRegionsBySize = x.classRegions;
    for (auto &pr : x.classRegionsBySize)
        sort(pr.begin(), pr.end(), [&](int32_t a, int32_t b) { return x.regionSize[a] > x.regionSize[b]; });
    x.sharedBytes = x.directory.empty() ? x.fileBytes : x.directory[0];
    return true;
}

static bool loadNSI3(const char *path, NSI2 &x, string &error) {
    LEReader rd(path);
    if (!rd.ok()) { error = rd.error; return false; }
    x.fileBytes = rd.size; x.slim = true;
    char magic[4] = {};
    if (!rd.bytes(magic, 4) || memcmp(magic, "NSI3", 4) != 0) { error = "not an NSI3 file"; return false; }
    x.rmin = rd.i32(); x.rmax = rd.i32(); x.smin = rd.i32(); x.smax = rd.i32();
    x.n = rd.i32(); x.nC = rd.i32(); x.nR = rd.i32(); x.nCols = rd.i32();
    x.nLeaf = 0;
    if (!rd.ok()) { error = rd.error; return false; }
    if (x.rmin < 1 || x.rmax < x.rmin || x.smin != x.rmin + 1 || x.smax < x.rmax + 1 ||
        x.n < 0 || x.nC < 0 || x.nR < 0 || x.nCols != x.rmax - x.rmin + 1 || x.nC > x.n) {
        error = "invalid NSI3 header"; return false;
    }
    if (!saneCount(x.n, 4, rd.limit - rd.pos)) { error = "invalid vertex count"; return false; }
    x.classOf.resize((size_t)x.n);
    for (int32_t v = 0; v < x.n; ++v) {
        int32_t c = rd.i32();
        if (c < -1 || c >= x.nC) { error = "classOf id out of range"; return false; }
        x.classOf[v] = c;
    }
    x.classSize.resize((size_t)x.nC);
    x.classRegions.resize((size_t)x.nC);
    x.regionSize.assign((size_t)x.nR, 0);
    for (int32_t c = 0; c < x.nC; ++c) {
        int32_t w = rd.i32(), z = rd.i32();
        if (!rd.ok() || w < 0 || z < 0 || !saneCount(z, 4, rd.limit - rd.pos)) {
            error = "invalid class record"; return false; }
        x.classSize[c] = w;
        x.classRegions[c].resize((size_t)z);
        int32_t prev = -1;
        for (int32_t i = 0; i < z; ++i) {
            int32_t rid = rd.i32();
            if (rid < 0 || rid >= x.nR || rid <= prev) { error = "class profile not sorted/in range"; return false; }
            prev = rid;
            x.classRegions[c][i] = rid;
            x.regionSize[rid] += w;                       // reconstruct |M| = Sum of its class weights
        }
    }
    if (!rd.ok()) { error = rd.error; return false; }
    x.directory.resize((size_t)x.nCols);
    for (int32_t i = 0; i < x.nCols; ++i) x.directory[i] = (int64_t)rd.u64();
    if (!rd.ok()) { error = rd.error; return false; }
    x.columns.resize((size_t)x.nCols);
    x.rToColumn.assign((size_t)(x.rmax - x.rmin + 1), -1);
    for (int32_t ci = 0; ci < x.nCols; ++ci) {
        const int64_t begin = x.directory[ci];
        const int64_t end = (ci + 1 < x.nCols) ? x.directory[ci + 1] : x.fileBytes;
        if (!rd.seek(begin, end)) { error = rd.error; return false; }
        NSI2Column &col = x.columns[ci];
        col.fileOffset = begin; col.fileBytes = end - begin;
        col.r = rd.i32(); col.boundary = rd.i32();
        if (col.r < x.rmin || col.r > x.rmax || col.boundary != col.r + 1 ||
            x.rToColumn[col.r - x.rmin] != -1) { error = "invalid NSI3 column header"; return false; }
        x.rToColumn[col.r - x.rmin] = ci;
        int32_t nm = rd.i32();
        if (!saneCount(nm, 4, rd.limit - rd.pos)) { error = "invalid mergeable list"; return false; }
        col.mergeableRegions.resize((size_t)nm);
        for (int32_t &rid : col.mergeableRegions) rid = rd.i32();
        if (!buildColumnRepresentatives(x, col, error)) return false;
        int32_t np = rd.i32();
        if (!rd.ok() || np < 0 || !saneCount(np, 8, rd.limit - rd.pos)) {
            error = "invalid exception-pattern count"; return false; }
        col.patterns.resize((size_t)np);
        col.origPid.resize((size_t)np);
        for (int32_t pi = 0; pi < np; ++pi) {
            col.origPid[pi] = rd.i32();
            int32_t len = rd.u8();
            if (!rd.ok() || len <= 0 || len > col.r) { error = "invalid exception comp"; return false; }
            NSI2Pattern &P = col.patterns[pi];
            P.comp.resize((size_t)len);
            int sum = 0;
            for (int32_t i = 0; i < len; ++i) {
                int32_t c = rd.i32(); int16_t b = (int16_t)rd.u16();
                if (c < 0 || c >= x.nC || b <= 0 || b > col.r) { error = "invalid exception comp entry"; return false; }
                P.comp[i] = {c, b}; sum += b;
            }
            if (sum != col.r) { error = "exception comp does not sum to r"; return false; }
            P.boundaryCore = rd.f64(); P.closedFrom = rd.i32();
            P.cP = -1;                                     // never stored; recovered from profiles
            uint64_t h = 0;
            for (auto &cm : P.comp) h ^= mixCV(cm.first, cm.second);
            col.patIdx[h].push_back(pi);
        }
        int32_t nCells = rd.i32();
        if (!rd.ok() || nCells < 0) { error = "invalid residue cell count"; return false; }
        col.resid.resize((size_t)nCells);
        for (auto &rr : col.resid) {
            int64_t m = (int64_t)rd.u64();
            if (!rd.ok() || m < 0 || !saneCount((int32_t)m, 12, rd.limit - rd.pos)) {
                error = "invalid residue list"; return false; }
            rr.resize((size_t)m);
            for (auto &pc : rr) { pc.first = rd.i32(); pc.second = rd.f64(); }
        }
        int32_t nd = rd.i32();
        if (!rd.ok() || nd < 0) { error = "invalid dists count"; return false; }
        col.dists.resize((size_t)nd);
        for (auto &d : col.dists) {
            int32_t m = rd.i32();
            if (!rd.ok() || m < 0 || !saneCount(m, 16, rd.limit - rd.pos)) { error = "invalid dist"; return false; }
            d.resize((size_t)m);
            for (auto &kv : d) { kv.first = rd.f64(); kv.second = rd.f64(); }
        }
        if (!rd.ok()) { error = rd.error; return false; }
    }
    // The binomial table must cover the largest cP a query can produce. NSI2 sizes it from the
    // STORED pattern cP values; NSI3 has none, so it is sized from the reconstructed region sizes,
    // which upper-bound cP by definition. Getting this wrong silently returns 0 for EVERY certified
    // answer, because C2 is table-based and clamps out of range.
    {
        int maxCombN = x.smax;
        for (int32_t z : x.regionSize) maxCombN = max(maxCombN, z);
        if (!build_ncr2(maxCombN, x.smax - x.rmin)) { error = "binomial table build failed"; return false; }
    }
    x.classRegionsBySize = x.classRegions;
    for (auto &pr : x.classRegionsBySize)
        sort(pr.begin(), pr.end(), [&](int32_t a, int32_t b) { return x.regionSize[a] > x.regionSize[b]; });
    x.sharedBytes = x.directory.empty() ? x.fileBytes : x.directory[0];
    return true;
}

static bool loadNSI2(const char *path, NSI2 &x, string &error) {
    LEReader rd(path);
    if (!rd.ok()) { error = rd.error; return false; }
    x.fileBytes = rd.size;
    char magic[4] = {};
    if (!rd.bytes(magic, 4) || memcmp(magic, "NSI2", 4) != 0) {
        error = "not an NSI2 file"; return false;
    }
    x.rmin = rd.i32(); x.rmax = rd.i32(); x.smin = rd.i32(); x.smax = rd.i32();
    x.n = rd.i32(); x.nC = rd.i32(); x.nR = rd.i32(); x.nLeaf = rd.i32(); x.nCols = rd.i32();
    if (!rd.ok()) { error = rd.error; return false; }
    if (x.rmin < 1 || x.rmax < x.rmin || x.rmax > UINT8_MAX ||
        x.smin != x.rmin + 1 || x.smax < x.rmax + 1 ||
        x.n < 0 || x.nC < 0 || x.nR < 0 || x.nLeaf < 0 ||
        x.nCols != x.rmax - x.rmin + 1 ||
        x.nC > x.n || x.n > rd.size / 4 || x.nC > rd.size / 8 ||
        x.nR > rd.size / 4 || x.nLeaf > rd.size / 8) {
        error = "invalid NSI2 header"; return false;
    }
    if (!saneCount(x.n, 4, rd.limit - rd.pos)) { error = "invalid vertex count"; return false; }
    x.classOf.resize((size_t)x.n);
    vector<int32_t> observed((size_t)x.nC, 0);
    for (int32_t v = 0; v < x.n; ++v) {
        int32_t c = rd.i32();
        if (c < -1 || c >= x.nC) { error = "classOf id out of range"; return false; }
        x.classOf[v] = c;
        if (c >= 0) observed[c]++;
    }

    x.classSize.resize((size_t)x.nC);
    x.classRegions.resize((size_t)x.nC);
    for (int32_t c = 0; c < x.nC; ++c) {
        int32_t w = rd.i32(), z = rd.i32();
        if (w <= 0 || w != observed[c] || !saneCount(z, 4, rd.limit - rd.pos)) {
            error = "invalid class weight/profile"; return false;
        }
        x.classSize[c] = w;
        auto &p = x.classRegions[c]; p.resize((size_t)z);
        int32_t prev = -1;
        for (int32_t &rid : p) {
            rid = rd.i32();
            if (rid < 0 || rid >= x.nR || rid <= prev) {
                error = "class profile is not a sorted region-id set"; return false;
            }
            prev = rid;
        }
    }

    x.regionSize.resize((size_t)x.nR);
    for (int32_t rid = 0; rid < x.nR; ++rid) {
        int32_t z = rd.i32();
        if (z < x.smin || z > x.n || !saneCount(z, 4, rd.limit - rd.pos)) {
            error = "invalid region size"; return false;
        }
        x.regionSize[rid] = z;
        // Verify the serialized list against the equivalent shared profiles,
        // but do not duplicate it in RAM: query containment uses the profiles.
        int32_t prevVertex = -1;
        for (int32_t i = 0; i < z; ++i) {
            int32_t v = rd.i32();
            if (v < 0 || v >= x.n || v <= prevVertex || x.classOf[v] < 0 ||
                !binary_search(x.classRegions[x.classOf[v]].begin(),
                               x.classRegions[x.classOf[v]].end(), rid)) {
                error = "region vertex list disagrees with shared class profiles"; return false;
            }
            prevVertex = v;
        }
    }
    for (int32_t c = 0; c < x.nC; ++c)
        for (int32_t rid : x.classRegions[c])
            if (x.classSize[c] > x.regionSize[rid]) {
                error = "class weight exceeds a containing region size"; return false;
            }
    vector<int64_t> reconstructedRegionSize((size_t)x.nR, 0);
    for (int32_t c = 0; c < x.nC; ++c)
        for (int32_t rid : x.classRegions[c]) reconstructedRegionSize[rid] += x.classSize[c];
    for (int32_t rid = 0; rid < x.nR; ++rid)
        if (reconstructedRegionSize[rid] != x.regionSize[rid]) {
            error = "region size disagrees with shared class profiles"; return false;
        }

    for (int32_t lid = 0; lid < x.nLeaf; ++lid) {
        int32_t m = rd.i32(); (void)rd.i32();             // m, immutable leaf target T
        if (m < 0 || m > x.nC || !saneCount(m, 12, rd.limit - rd.pos)) {
            error = "invalid shared leaf"; return false;
        }
        // Each coordinate is {class:int32,h:int16,n:int16,ell:int16,u:int16}.
        if (!rd.skip((int64_t)m * 12)) { error = rd.error; return false; }
    }

    x.directory.resize((size_t)x.nCols);
    const int64_t dirStart = rd.pos;
    if (!saneCount(x.nCols, 8, rd.limit - rd.pos)) { error = "truncated r-directory"; return false; }
    for (int32_t i = 0; i < x.nCols; ++i) {
        uint64_t off = rd.u64();
        if (off > (uint64_t)numeric_limits<int64_t>::max()) { error = "column offset overflow"; return false; }
        x.directory[i] = (int64_t)off;
    }
    const int64_t dirEnd = rd.pos;
    (void)dirStart;
    x.sharedBytes = dirEnd;
    for (int32_t i = 0; i < x.nCols; ++i) {
        int64_t end = (i + 1 < x.nCols) ? x.directory[i + 1] : x.fileBytes;
        if (x.directory[i] < dirEnd || end <= x.directory[i] || end > x.fileBytes ||
            (i > 0 && x.directory[i] <= x.directory[i - 1])) {
            error = "invalid or non-monotone r-directory"; return false;
        }
    }
    if (!x.directory.empty() && x.directory.front() != dirEnd) {
        error = "unexpected bytes between the r-directory and first column"; return false;
    }

    x.columns.resize((size_t)x.nCols);
    x.rToColumn.assign((size_t)(x.rmax - x.rmin + 1), -1);
    int maxCombN = x.smax;
    for (int32_t ci = 0; ci < x.nCols; ++ci) {
        const int64_t begin = x.directory[ci];
        const int64_t end = (ci + 1 < x.nCols) ? x.directory[ci + 1] : x.fileBytes;
        if (!rd.seek(begin, end)) { error = rd.error; return false; }
        NSI2Column &col = x.columns[ci];
        col.fileOffset = begin; col.fileBytes = end - begin;
        col.r = rd.i32(); col.boundary = rd.i32();
        if (col.r < x.rmin || col.r > x.rmax || col.boundary != col.r + 1 ||
            x.rToColumn[col.r - x.rmin] != -1) {
            error = "invalid or duplicate r-column header"; return false;
        }
        x.rToColumn[col.r - x.rmin] = ci;

        int32_t nm = rd.i32();
        if (!saneCount(nm, 4, rd.limit - rd.pos)) { error = "invalid mergeable-region list"; return false; }
        col.mergeableRegions.resize((size_t)nm);
        for (int32_t &rid : col.mergeableRegions) rid = rd.i32();
        if (!buildColumnRepresentatives(x, col, error)) return false;

        int32_t np = rd.i32();
        if (np < 0 || np > 100000000 || !saneCount(np, 31, rd.limit - rd.pos)) {
            error = "invalid pattern count"; return false;
        }
        col.patterns.resize((size_t)np);
        col.patIdx.reserve((size_t)np * 2 + 1);
        for (int32_t pi = 0; pi < np; ++pi) {
            NSI2Pattern &p = col.patterns[pi];
            uint8_t len = rd.u8();
            if (len == 0 || len > col.r || !saneCount(len, 6, rd.limit - rd.pos)) {
                error = "invalid pattern composition length"; return false;
            }
            p.comp.resize(len);
            int sum = 0; int32_t prevClass = -1; uint64_t h = 0;
            for (auto &cm : p.comp) {
                cm.first = rd.i32(); cm.second = rd.i16();
                if (cm.first < 0 || cm.first >= x.nC || cm.first <= prevClass ||
                    cm.second <= 0 || sum + cm.second > col.r) {
                    error = "invalid pattern composition"; return false;
                }
                if (col.classRep[cm.first] != cm.first) {
                    error = "pattern coordinate is not the canonical active representative"; return false;
                }
                prevClass = cm.first; sum += cm.second; h ^= mixCV(cm.first, cm.second);
            }
            p.mult = rd.i64(); p.cP = rd.i32(); p.boundaryCore = rd.f64(); p.closedFrom = rd.i32();
            if (sum != col.r || p.mult <= 0 || p.cP < col.boundary || p.cP > x.n ||
                !isfinite(p.boundaryCore) || p.boundaryCore < 0 ||
                (p.closedFrom != -1 && (p.closedFrom < col.boundary || p.closedFrom > x.smax))) {
                error = "invalid pattern metadata"; return false;
            }
            col.patIdx[h].push_back(pi);
            maxCombN = max(maxCombN, p.cP);
        }

        int32_t nr = rd.i32();
        const int32_t expectedResid = x.smax - col.boundary;
        if (nr != expectedResid || !saneCount(nr, 8, rd.limit - rd.pos)) {
            error = "wrong residue-cell count"; return false;
        }
        col.resid.resize((size_t)nr);
        for (auto &cell : col.resid) {
            int64_t z = rd.i64();
            if (!saneCount(z, 12, rd.limit - rd.pos)) { error = "invalid residue dictionary"; return false; }
            cell.resize((size_t)z);
            for (auto &pc : cell) {
                pc.first = rd.i32(); pc.second = rd.f64();
                if (pc.first < 0 || pc.first >= np || !isfinite(pc.second) || pc.second < 0) {
                    error = "invalid residue entry"; return false;
                }
            }
            sort(cell.begin(), cell.end());
            for (size_t i = 1; i < cell.size(); ++i)
                if (cell[i-1].first == cell[i].first) { error = "duplicate residue pattern id"; return false; }
        }

        int32_t nd = rd.i32();
        const int32_t expectedDists = x.smax - col.boundary + 1;
        if (nd != expectedDists) { error = "wrong distribution-cell count"; return false; }
        col.dists.resize((size_t)nd);
        for (auto &dist : col.dists) {
            int32_t z = rd.i32();
            if (!saneCount(z, 16, rd.limit - rd.pos)) { error = "invalid distribution"; return false; }
            dist.resize((size_t)z);
            for (auto &kv : dist) {
                kv.first = rd.f64(); kv.second = rd.f64();
                if (!isfinite(kv.first) || !isfinite(kv.second) || kv.first < 0 || kv.second < 0) {
                    error = "invalid distribution entry"; return false;
                }
            }
        }
        if (!rd.ok() || rd.pos != end) {
            error = rd.ok() ? "column byte count disagrees with its directory extent" : rd.error;
            return false;
        }
    }
    for (int32_t ri : x.rToColumn) if (ri < 0) { error = "missing r-column"; return false; }
    for (int32_t z : x.regionSize) maxCombN = max(maxCombN, z);
    if (!build_ncr2(maxCombN, x.smax - x.rmin)) {
        error = "cannot allocate NSI2 binomial table"; return false;
    }
    return true;
}

enum class QueryCode {
    Ok,
    ROutOfRange,
    SOutOfRange,
    BadVertex,
    NotClique,
    CorruptIndex
};

static const char *queryCodeName(QueryCode c) {
    switch (c) {
        case QueryCode::Ok: return "ok";
        case QueryCode::ROutOfRange: return "r out of range";
        case QueryCode::SOutOfRange: return "s out of range";
        case QueryCode::BadVertex: return "vertex out of range";
        case QueryCode::NotClique: return "vertices do not form an r-clique";
        case QueryCode::CorruptIndex: return "missing/corrupt index entry";
    }
    return "unknown query error";
}

struct ValidationGraph {
    int n = 0;
    vector<int32_t> off, adj;
    bool adjacent(int u, int v) const {
        if (u < 0 || v < 0 || u >= n || v >= n) return false;
        return binary_search(adj.begin() + off[u], adj.begin() + off[u + 1], v);
    }
};

// macOS reports ru_maxrss in bytes while Linux reports KiB.  We use this
// process high-water mark only as an incremental retrieval measurement: the
// baseline is sampled after the index and graph loads, immediately before retrieval.
static double peakRssMiB() {
    rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
#ifdef __APPLE__
    return (double)ru.ru_maxrss / (1024.0 * 1024.0);
#else
    return (double)ru.ru_maxrss / 1024.0;
#endif
}

static bool loadValidationGraph(const char *path, ValidationGraph &g, string &error) {
    FILE *f = fopen(path, "r");
    if (!f) { error = string("cannot open graph ") + path; return false; }
    long long hn = 0, hm = 0;
    if (fscanf(f, "%lld %lld", &hn, &hm) != 2 || hn < 0 || hn > INT32_MAX || hm < 0) {
        fclose(f); error = "invalid graph header"; return false;
    }
    vector<pair<int32_t,int32_t>> es;
    if (hm <= (long long)numeric_limits<size_t>::max() / 2)
        es.reserve((size_t)hm * 2);
    int u = 0, v = 0, maxId = (int)hn - 1;
    while (fscanf(f, "%d %d", &u, &v) == 2) {
        if (u < 0 || v < 0) { fclose(f); error = "negative graph vertex id"; return false; }
        maxId = max(maxId, max(u, v));
        if (u == v) continue;
        es.push_back({u, v}); es.push_back({v, u});
    }
    if (ferror(f)) { fclose(f); error = "graph read failure"; return false; }
    fclose(f);
    if (maxId == INT32_MAX) { error = "graph vertex id overflow"; return false; }
    g.n = maxId + 1;
    sort(es.begin(), es.end());
    es.erase(unique(es.begin(), es.end()), es.end());
    g.off.assign((size_t)g.n + 1, 0);
    for (auto &e : es) g.off[(size_t)e.first + 1]++;
    for (int i = 0; i < g.n; ++i) g.off[i + 1] += g.off[i];
    g.adj.resize(es.size());
    vector<int32_t> cur(g.off.begin(), g.off.end() - 1);
    for (auto &e : es) g.adj[cur[e.first]++] = e.second;
    return true;
}

struct Query2 {
    const NSI2 &x;
    vector<int16_t> repMultiplicity;
    vector<int32_t> touchedReps;
    vector<int32_t> queryClasses;

    explicit Query2(const NSI2 &x_) : x(x_) {
        repMultiplicity.assign((size_t)x.nC, 0);
        touchedReps.reserve((size_t)x.rmax);
        queryClasses.reserve((size_t)x.rmax);
    }

    const NSI2Column *column(int r) const {
        if (r < x.rmin || r > x.rmax) return nullptr;
        int32_t ci = x.rToColumn[r - x.rmin];
        return ci < 0 ? nullptr : &x.columns[ci];
    }

    int lookupPattern(const NSI2Column &col, const int *vs) {
        touchedReps.clear();
        auto clearCounts = [&]() {
            for (int32_t c : touchedReps) repMultiplicity[c] = 0;
        };
        for (int i = 0; i < col.r; ++i) {
            int v = vs[i];
            if (v < 0 || v >= x.n) { clearCounts(); return -1; }
            int32_t c = x.classOf[v];
            if (c < 0) { clearCounts(); return -1; }
            int32_t rep = col.classRep[c];
            if (rep < 0) { clearCounts(); return -1; }
            if (repMultiplicity[rep] == 0) touchedReps.push_back(rep);
            ++repMultiplicity[rep];
        }
        uint64_t h = 0;
        for (int32_t rep : touchedReps) h ^= mixCV(rep, repMultiplicity[rep]);
        auto hit = col.patIdx.find(h);
        int found = -1;
        if (hit != col.patIdx.end()) {
            for (int32_t pi : hit->second) {
                const Comp2 &candidate = col.patterns[pi].comp;
                if (candidate.size() != touchedReps.size()) continue;
                bool same = true;
                for (const auto &cm : candidate)
                    if (repMultiplicity[cm.first] != cm.second) { same = false; break; }
                if (same) { found = pi; break; }
            }
        }
        clearCounts();
        return found;
    }

    int lookupMergeable(const NSI2Column &col, const int *vs) {
        queryClasses.clear();
        int base = -1;
        size_t best = numeric_limits<size_t>::max();
        for (int i = 0; i < col.r; ++i) {
            int v = vs[i];
            if (v < 0 || v >= x.n) return -1;
            int32_t c = x.classOf[v];
            if (c < 0) return -1;
            queryClasses.push_back(c);
            if (x.classRegions[c].size() < best) { best = x.classRegions[c].size(); base = i; }
        }
        if (base < 0) return -1;
        const auto &seed = x.classRegions[queryClasses[base]];
        for (int32_t rid : seed) {
            if (!col.mergeMask[rid]) continue;
            bool inAll = true;
            for (int i = 0; i < col.r; ++i) {
                if (i == base) continue;
                const auto &p = x.classRegions[queryClasses[i]];
                if (!binary_search(p.begin(), p.end(), rid)) { inAll = false; break; }
            }
            if (inAll) return rid;
        }
        return -1;
    }

    struct Resolved {
        const NSI2Column *col = nullptr;
        int32_t pattern = -1;
        int32_t mergeable = -1;
        int32_t cpComputed = -1;                   // §222 NSI3: cP recovered from class profiles
    };

    QueryCode resolve(const int *vs, int r, Resolved &out) {
        out = {};
        out.col = column(r);
        if (!out.col) return QueryCode::ROutOfRange;
        for (int i = 0; i < r; ++i)
            if (vs[i] < 0 || vs[i] >= x.n) return QueryCode::BadVertex;
        out.pattern = lookupPattern(*out.col, vs);
        if (x.slim) out.cpComputed = cpFromProfiles(vs, r);       // needed by both slim branches
        if (out.pattern < 0) {
            if (x.slim) { /* cP already recovered above */ }      // §222: certified => closed form
            else        out.mergeable  = lookupMergeable(*out.col, vs);
        }
        return QueryCode::Ok;
    }

    double value(const Resolved &q, int s, QueryCode &code) const {
        const NSI2Column &col = *q.col;
        if (s < col.boundary || s > x.smax) { code = QueryCode::SOutOfRange; return 0; }
        if (q.pattern >= 0) {
            const NSI2Pattern &p = col.patterns[q.pattern];
            const int32_t pid = x.slim ? col.origPid[q.pattern] : q.pattern;
            if (p.closedFrom >= 0 && s >= p.closedFrom) {
                code = QueryCode::Ok;
                // slim stores no cP: recover it exactly the same way the certified path does
                return C2((x.slim ? q.cpComputed : p.cP) - col.r, s - col.r);
            }
            if (s == col.boundary) { code = QueryCode::Ok; return p.boundaryCore; }
            const auto &cell = col.resid[s - col.boundary - 1];
            auto it = lower_bound(cell.begin(), cell.end(), make_pair(pid, -numeric_limits<double>::infinity()));
            if (it == cell.end() || it->first != pid) {
                code = QueryCode::CorruptIndex; return 0;
            }
            code = QueryCode::Ok; return it->second;
        }
        if (q.mergeable >= 0) {
            code = QueryCode::Ok;
            return C2(x.regionSize[q.mergeable] - col.r, s - col.r);
        }
        if (x.slim && q.cpComputed >= 0) {                 // §222: certified, reconstructed
            code = QueryCode::Ok;
            return C2(q.cpComputed - col.r, s - col.r);
        }
        code = QueryCode::Ok;                              // not supported by any indexed s-clique
        return 0.0;
    }

    // §222 SLIM-INDEX KERNEL: cP recovered at query time instead of read from a stored record.
    // A region M hosts an r-clique iff M contains every class of it (classes are wholly in or out of
    // a region, Thm 3.2), so the hosting set is the INTERSECTION of the vertices' class profiles and
    //     cP = max{ regionSize[M] : M in that intersection }.
    // Returns -1 when no region hosts the clique (kappa is then 0 in every cell).
    int32_t cpFromProfiles(const int *vs, int r) {
        queryClasses.clear();
        int base = -1; size_t best = numeric_limits<size_t>::max();
        for (int i = 0; i < r; ++i) {
            int v = vs[i];
            if (v < 0 || v >= x.n) return -1;
            int32_t c = x.classOf[v];
            if (c < 0) return -1;
            queryClasses.push_back(c);
            if (x.classRegions[c].size() < best) { best = x.classRegions[c].size(); base = i; }
        }
        if (base < 0) return -1;
        // Only the LARGEST hosting region matters, so walk the seed profile in decreasing size order
        // and stop at the first region present in every other profile. Membership uses the id-sorted
        // copy, so this is the same predicate with an early exit instead of a full scan.
        const auto &seed = x.classRegionsBySize[queryClasses[base]];
        for (int32_t rid : seed) {
            bool inAll = true;
            for (int i = 0; i < r && inAll; ++i) {
                if (i == base) continue;
                const auto &p2 = x.classRegions[queryClasses[i]];
                if (!binary_search(p2.begin(), p2.end(), rid)) inAll = false;
            }
            if (inAll) return x.regionSize[rid];
        }
        return -1;
    }

    double pointKernel(const int *vs, int r, int s, QueryCode &code) {
        Resolved q;
        code = resolve(vs, r, q);
        return code == QueryCode::Ok ? value(q, s, code) : 0.0;
    }

    QueryCode rowKernel(const int *vs, int r, vector<double> &out) {
        Resolved q;
        QueryCode code = resolve(vs, r, q);
        if (code != QueryCode::Ok) { out.clear(); return code; }
        out.resize((size_t)(x.smax - q.col->boundary + 1));
        for (int s = q.col->boundary; s <= x.smax; ++s) {
            out[s - q.col->boundary] = value(q, s, code);
            if (code != QueryCode::Ok) { out.clear(); return code; }
        }
        return QueryCode::Ok;
    }

    static QueryCode validateClique(const ValidationGraph &g, const int *vs, int r) {
        for (int i = 0; i < r; ++i) {
            if (vs[i] < 0 || vs[i] >= g.n) return QueryCode::BadVertex;
            for (int j = 0; j < i; ++j)
                if (vs[i] == vs[j] || !g.adjacent(vs[i], vs[j])) return QueryCode::NotClique;
        }
        return QueryCode::Ok;
    }

    double pointValidated(const ValidationGraph &g, const int *vs, int r, int s, QueryCode &code) {
        code = validateClique(g, vs, r);
        return code == QueryCode::Ok ? pointKernel(vs, r, s, code) : 0.0;
    }

    QueryCode rowValidated(const ValidationGraph &g, const int *vs, int r, vector<double> &out) {
        QueryCode code = validateClique(g, vs, r);
        if (code != QueryCode::Ok) { out.clear(); return code; }
        return rowKernel(vs, r, out);
    }
};

struct UnionFind {
    vector<int32_t> parent;
    vector<uint8_t> rank;
    explicit UnionFind(size_t n = 0) : parent(n), rank(n, 0) {
        iota(parent.begin(), parent.end(), 0);
    }
    int32_t find(int32_t x) {
        int32_t root = x;
        while (parent[root] != root) root = parent[root];
        while (parent[x] != x) { int32_t next = parent[x]; parent[x] = root; x = next; }
        return root;
    }
    void unite(int32_t a, int32_t b) {
        a = find(a); b = find(b);
        if (a == b) return;
        if (rank[a] < rank[b]) swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) ++rank[a];
    }
};

template<class Emit>
static void enumerateCliques(const ValidationGraph &g, int cliqueSize, Emit emit) {
    vector<int32_t> prefix;
    prefix.reserve((size_t)cliqueSize);
    vector<int32_t> all((size_t)g.n);
    iota(all.begin(), all.end(), 0);
    function<void(const vector<int32_t>&)> visit = [&](const vector<int32_t> &candidates) {
        if ((int)prefix.size() == cliqueSize) { emit(prefix); return; }
        const int needed = cliqueSize - (int)prefix.size();
        for (size_t i = 0; i + (size_t)needed <= candidates.size(); ++i) {
            const int32_t v = candidates[i];
            prefix.push_back(v);
            vector<int32_t> next;
            next.reserve(candidates.size() - i - 1);
            for (size_t j = i + 1; j < candidates.size(); ++j)
                if (g.adjacent(v, candidates[j])) next.push_back(candidates[j]);
            visit(next);
            prefix.pop_back();
        }
    };
    if (cliqueSize > 0) visit(all);
}

static void appendRSubcliques(const vector<int32_t> &witness, int r,
                              vector<vector<int32_t>> &out) {
    vector<int32_t> current;
    current.reserve((size_t)r);
    function<void(int)> choose = [&](int at) {
        if ((int)current.size() == r) { out.push_back(current); return; }
        const int remaining = r - (int)current.size();
        for (int i = at; i + remaining <= (int)witness.size(); ++i) {
            current.push_back(witness[i]); choose(i + 1); current.pop_back();
        }
    };
    choose(0);
}

static vector<int32_t> sizesFromUnionFind(UnionFind &uf, const vector<uint8_t> *keep = nullptr) {
    unordered_map<int32_t,int32_t> sizes;
    for (size_t i = 0; i < uf.parent.size(); ++i) {
        if (keep && !(*keep)[i]) continue;
        ++sizes[uf.find((int32_t)i)];
    }
    vector<int32_t> ans;
    ans.reserve(sizes.size());
    for (const auto &kv : sizes) ans.push_back(kv.second);
    sort(ans.begin(), ans.end());
    return ans;
}

static void printSizesRLE(const vector<int32_t> &sizes) {
    for (size_t i = 0; i < sizes.size();) {
        size_t j = i + 1;
        while (j < sizes.size() && sizes[j] == sizes[i]) ++j;
        printf("%s%d%s", i ? "," : "", sizes[i], j - i > 1 ? "x" : "");
        if (j - i > 1) printf("%zu", j - i);
        i = j;
    }
}

struct NucleiResult {
    size_t selected = 0, survivingWitnesses = 0, incidences = 0;
    double collectMs = 0, witnessMs = 0, unionMs = 0, totalMs = 0;
    double peakExtraMiB = 0;
    vector<int32_t> sizes;
};

// Retrieval path: core membership comes exclusively from NSI2; graph work is
// restricted to materializing the requested r-cliques and s-witness adjacency.
static bool retrieveConnectedNuclei(const NSI2 &x, Query2 &q, const ValidationGraph &g,
                                    int r, int s, int k, double rssAfterLoad,
                                    NucleiResult &result, string &error) {
    if (!q.column(r)) { error = "r out of range"; return false; }
    if (s < r + 1 || s > x.smax) { error = "s out of range"; return false; }
    unordered_map<vector<int32_t>,int32_t,IntVecHash> selected;
    auto total0 = SteadyClock::now();
    auto collect0 = SteadyClock::now();
    enumerateCliques(g, r, [&](const vector<int32_t> &clique) {
        QueryCode code;
        const double core = q.pointKernel(clique.data(), r, s, code);
        if (code != QueryCode::Ok) { error = queryCodeName(code); return; }
        if (core >= k) selected.emplace(clique, (int32_t)selected.size());
    });
    if (!error.empty()) return false;
    result.collectMs = chrono::duration<double,milli>(SteadyClock::now() - collect0).count();
    result.selected = selected.size();
    UnionFind uf(selected.size());
    auto witness0 = SteadyClock::now();
    enumerateCliques(g, s, [&](const vector<int32_t> &witness) {
        vector<vector<int32_t>> rs;
        appendRSubcliques(witness, r, rs);
        vector<int32_t> ids; ids.reserve(rs.size());
        for (const auto &rc : rs) {
            auto it = selected.find(rc);
            if (it == selected.end()) return;
            ids.push_back(it->second);
        }
        ++result.survivingWitnesses;
        result.incidences += ids.size();
        auto union0 = SteadyClock::now();
        for (size_t i = 1; i < ids.size(); ++i) uf.unite(ids[0], ids[i]);
        result.unionMs += chrono::duration<double,milli>(SteadyClock::now() - union0).count();
    });
    result.witnessMs = chrono::duration<double,milli>(SteadyClock::now() - witness0).count();
    result.sizes = sizesFromUnionFind(uf);
    result.totalMs = chrono::duration<double,milli>(SteadyClock::now() - total0).count();
    result.peakExtraMiB = max(0.0, peakRssMiB() - rssAfterLoad);
    return true;
}

// Independent reference: enumerate the complete r/s incidence hypergraph and
// run ordinary minimum-support peeling before constructing witness components.
static void directConnectedNuclei(const ValidationGraph &g, int r, int s, int k,
                                  size_t &selectedCount, vector<int32_t> &sizes) {
    unordered_map<vector<int32_t>,int32_t,IntVecHash> ids;
    vector<vector<int32_t>> witnesses;
    vector<vector<int32_t>> incident;
    enumerateCliques(g, s, [&](const vector<int32_t> &witness) {
        vector<vector<int32_t>> rs;
        appendRSubcliques(witness, r, rs);
        vector<int32_t> w; w.reserve(rs.size());
        for (const auto &rc : rs) {
            auto it = ids.find(rc);
            int32_t id;
            if (it == ids.end()) {
                id = (int32_t)ids.size(); ids.emplace(rc, id); incident.emplace_back();
            } else id = it->second;
            w.push_back(id);
        }
        const int32_t wi = (int32_t)witnesses.size();
        for (int32_t id : w) incident[id].push_back(wi);
        witnesses.push_back(std::move(w));
    });
    vector<int32_t> degree(ids.size()), core(ids.size());
    for (size_t i = 0; i < ids.size(); ++i) degree[i] = (int32_t)incident[i].size();
    priority_queue<pair<int32_t,int32_t>, vector<pair<int32_t,int32_t>>,
                   greater<pair<int32_t,int32_t>>> heap;
    for (size_t i = 0; i < degree.size(); ++i) heap.push({degree[i], (int32_t)i});
    vector<uint8_t> dead(ids.size(), 0), liveWitness(witnesses.size(), 1);
    while (!heap.empty()) {
        auto [d, u] = heap.top(); heap.pop();
        if (dead[u] || d != degree[u]) continue;
        dead[u] = 1; core[u] = d;
        for (int32_t wi : incident[u]) {
            if (!liveWitness[wi]) continue;
            liveWitness[wi] = 0;
            for (int32_t v : witnesses[wi])
                if (!dead[v] && degree[v] > d) {
                    --degree[v]; heap.push({degree[v], v});
                }
        }
    }
    vector<uint8_t> keep(ids.size(), 0);
    selectedCount = 0;
    for (size_t i = 0; i < core.size(); ++i)
        if (core[i] >= k) { keep[i] = 1; ++selectedCount; }
    UnionFind uf(ids.size());
    for (const auto &w : witnesses) {
        bool survives = true;
        for (int32_t id : w) if (!keep[id]) { survives = false; break; }
        if (survives) for (size_t i = 1; i < w.size(); ++i) uf.unite(w[0], w[i]);
    }
    sizes = sizesFromUnionFind(uf, &keep);
}

static bool parseIntArg(const char *s, int &v) {
    if (!s || !*s) return false;
    char *e = nullptr; errno = 0;
    long z = strtol(s, &e, 10);
    if (errno || !e || *e || z < INT32_MIN || z > INT32_MAX) return false;
    v = (int)z; return true;
}

static bool readFixedCliques(const char *path, int r, vector<int> &flat, string &error) {
    ifstream in(path);
    if (!in) { error = string("cannot open ") + path; return false; }
    string line; long long lineNo = 0;
    while (getline(in, line)) {
        ++lineNo;
        size_t hash = line.find('#'); if (hash != string::npos) line.resize(hash);
        istringstream ss(line);
        vector<long long> row; long long z;
        while (ss >> z) row.push_back(z);
        if (row.empty()) continue;
        if ((int)row.size() != r) {
            error = "query line " + to_string(lineNo) + " must contain exactly " + to_string(r) + " vertices";
            return false;
        }
        for (long long q : row) {
            if (q < INT32_MIN || q > INT32_MAX) { error = "vertex id overflow on query line " + to_string(lineNo); return false; }
            flat.push_back((int)q);
        }
    }
    return true;
}

struct BenchQuery2 { int32_t r = 0, s = 0; uint32_t off = 0; };

static bool readBenchQueries(const char *path, const NSI2 &x,
                             vector<BenchQuery2> &queries, vector<int> &flat,
                             string &error) {
    ifstream in(path);
    if (!in) { error = string("cannot open ") + path; return false; }
    string line; long long lineNo = 0;
    while (getline(in, line)) {
        ++lineNo;
        size_t hash = line.find('#'); if (hash != string::npos) line.resize(hash);
        istringstream ss(line);
        vector<long long> a; long long z;
        while (ss >> z) a.push_back(z);
        if (a.empty()) continue;
        if (a[0] < x.rmin || a[0] > x.rmax) {
            error = "r out of range on query line " + to_string(lineNo); return false;
        }
        int r = (int)a[0], s = 0, firstVertex = 0;
        if ((int)a.size() == r + 1) {
            const NSI2Column &col = x.columns[x.rToColumn[r - x.rmin]];
            s = col.boundary + (int)(queries.size() % (size_t)(x.smax - col.boundary + 1));
            firstVertex = 1;                              // r v1 ... vr; cycle point cells
        } else if ((int)a.size() == r + 2) {
            s = (int)a[1]; firstVertex = 2;               // r s v1 ... vr
        } else {
            error = "query line " + to_string(lineNo) + " must be 'r v...' or 'r s v...'";
            return false;
        }
        const NSI2Column &col = x.columns[x.rToColumn[r - x.rmin]];
        if (s < col.boundary || s > x.smax) {
            error = "s out of range on query line " + to_string(lineNo); return false;
        }
        if (flat.size() > UINT32_MAX) { error = "query file is too large"; return false; }
        BenchQuery2 q; q.r = r; q.s = s; q.off = (uint32_t)flat.size();
        for (int i = 0; i < r; ++i) {
            long long v = a[firstVertex + i];
            if (v < INT32_MIN || v > INT32_MAX) { error = "vertex id overflow on query line " + to_string(lineNo); return false; }
            flat.push_back((int)v);
        }
        queries.push_back(q);
    }
    return true;
}

struct LatencySummary { double median = 0, p95 = 0; size_t samples = 0; };

static LatencySummary summarizeLatency(vector<double> ns) {
    LatencySummary s; s.samples = ns.size();
    if (ns.empty()) return s;
    sort(ns.begin(), ns.end());
    const size_t n = ns.size();
    s.median = (n & 1) ? ns[n/2] : 0.5 * (ns[n/2 - 1] + ns[n/2]);
    size_t p = (size_t)ceil(0.95 * (double)n);
    if (p == 0) p = 1;
    s.p95 = ns[min(n - 1, p - 1)];
    return s;
}

static volatile double latencySink = 0;
static volatile uint64_t evictionSink = 0;

static void evictSoftwareCache(const vector<uint64_t> &eviction) {
    uint64_t z = evictionSink;
    // One word per 64-byte line: allocate every line in cache while avoiding
    // timing the eviction itself as part of the query.
    for (size_t i = 0; i < eviction.size(); i += 8) z += eviction[i];
    evictionSink = z;
}

template<class F>
static LatencySummary measureLatency(size_t samples, int repetitions,
                                     const vector<uint64_t> *eviction, F &&fn) {
    vector<double> times; times.reserve(samples);
    for (size_t i = 0; i < samples; ++i) {
        if (eviction) evictSoftwareCache(*eviction);
        auto a = SteadyClock::now();
        double z = 0;
        for (int rep = 0; rep < repetitions; ++rep) z += fn(i);
        auto b = SteadyClock::now();
        latencySink += z;
        times.push_back(chrono::duration<double,nano>(b - a).count() / repetitions);
    }
    return summarizeLatency(std::move(times));
}

static void printNSI2Usage(const char *prog) {
    fprintf(stderr,
        "NSI2 usage:\n"
        "  %s INDEX stats\n"
        "  %s INDEX point R S V1 ... VR                    # kernel-only\n"
        "  %s INDEX point-validated GRAPH R S V1 ... VR    # includes clique validation\n"
        "  %s INDEX row R V1 ... VR                        # all S=R+1..Smax\n"
        "  %s INDEX row-validated GRAPH R V1 ... VR\n"
        "  %s INDEX count R S K                             # K > 0\n"
        "  %s INDEX pointfile R S FILE\n"
        "  %s INDEX pointfile-validated GRAPH R S FILE\n"
        "  %s INDEX rowfile R FILE\n"
        "  %s INDEX rowfile-validated GRAPH R FILE\n"
        "  %s INDEX nuclei R S K GRAPH                     # connected k-(r,s)-nuclei + direct check\n"
        "  %s INDEX bench GRAPH QUERIES [--cold-mib N] [--warm-reps N]\n"
        "  %s INDEX archive [--vs SLIM.nsi3]               # E4 baseline: archive size accounting\n"
        "  %s INDEX archive-bench R QUERIES [--cap-gb N] [--cold-mib N]\n"
        "Bench query lines are either 'R V1 ... VR' (S cycles over the row) or\n"
        "'R S V1 ... VR'.  At least 1000 validated clique lines are required.\n"
        "archive modes need the FULL index (NSI2); the slim NSI3 drops the multiplicities.\n",
        prog, prog, prog, prog, prog, prog, prog, prog, prog, prog, prog, prog, prog, prog);
}

static int printQueryError(QueryCode code) {
    fprintf(stderr, "query failed: %s\n", queryCodeName(code));
    return 2;
}

static int mainNSI2(int argc, char **argv) {
    if (argc < 3) { printNSI2Usage(argv[0]); return 1; }
    NSI2 x; string error;
    auto load0 = SteadyClock::now();
    {
        FILE *mf = fopen(argv[1], "rb"); char mg[4] = {};
        const bool got = mf && fread(mg, 4, 1, mf) == 1;
        if (mf) fclose(mf);
        const bool isSlim = got && !memcmp(mg, "NSI3", 4);
        const bool isPack = got && !memcmp(mg, "NSI4", 4);
        if (!(isPack ? loadNSI4(argv[1], x, error)
                     : isSlim ? loadNSI3(argv[1], x, error)
                              : loadNSI2(argv[1], x, error))) {
            fprintf(stderr, "cannot load %s: %s\n", argv[1], error.c_str()); return 1;
        }
    }
    auto load1 = SteadyClock::now();
    const double loadMs = chrono::duration<double,milli>(load1 - load0).count();
    const string mode = argv[2];
    Query2 q(x);

    if (mode == "stats") {
        int64_t colTotal = 0;
        printf("index=%s format=NSI2 bytes=%lld shared-once=%lld per-column-bytes=%lld load=%.3f ms\n",
               argv[1], (long long)x.fileBytes, (long long)x.sharedBytes,
               (long long)(x.fileBytes - x.sharedBytes), loadMs);
        printf("plane r=%d..%d s<=%d n=%d classes=%d regions=%d shared-leaves=%d columns=%d\n",
               x.rmin, x.rmax, x.smax, x.n, x.nC, x.nR, x.nLeaf, x.nCols);
        printf("r  boundary  patterns  direct-regions  residues  column-bytes  offset\n");
        for (const auto &col : x.columns) {
            long long residues = 0;
            for (const auto &cell : col.resid) residues += (long long)cell.size();
            printf("%d  %d  %zu  %zu  %lld  %lld  %lld\n", col.r, col.boundary,
                   col.patterns.size(), col.mergeableRegions.size(), residues,
                   (long long)col.fileBytes, (long long)col.fileOffset);
            colTotal += col.fileBytes;
        }
        printf("byte-accounting shared-once(%lld) + per-column(%lld) = total(%lld)\n",
               (long long)x.sharedBytes, (long long)colTotal, (long long)x.fileBytes);
        return (x.sharedBytes + colTotal == x.fileBytes) ? 0 : 2;
    }

    if (mode == "nuclei") {
        int r = 0, s = 0, k = 0;
        if (argc != 7 || !parseIntArg(argv[3], r) || !parseIntArg(argv[4], s) ||
            !parseIntArg(argv[5], k) || k <= 0) {
            fprintf(stderr, "nuclei requires positive integer K: INDEX nuclei R S K GRAPH\n");
            return 1;
        }
        const NSI2Column *col = q.column(r);
        if (!col) return printQueryError(QueryCode::ROutOfRange);
        if (s < col->boundary || s > x.smax) return printQueryError(QueryCode::SOutOfRange);
        ValidationGraph g;
        if (!loadValidationGraph(argv[6], g, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        if (g.n != x.n) { fprintf(stderr, "graph/index vertex-count mismatch\n"); return 1; }
        const double rssAfterLoad = peakRssMiB();
        NucleiResult got;
        if (!retrieveConnectedNuclei(x, q, g, r, s, k, rssAfterLoad, got, error)) {
            fprintf(stderr, "nuclei retrieval failed: %s\n", error.c_str()); return 2;
        }
        auto reference0 = SteadyClock::now();
        size_t refSelected = 0; vector<int32_t> refSizes;
        directConnectedNuclei(g, r, s, k, refSelected, refSizes);
        const double referenceMs = chrono::duration<double,milli>(SteadyClock::now() - reference0).count();
        const bool exact = got.selected == refSelected && got.sizes == refSizes;
        printf("nuclei r=%d s=%d k=%d selected=%zu components=%zu sizes=", r, s, k,
               got.selected, got.sizes.size());
        printSizesRLE(got.sizes);
        printf("\n");
        printf("index_load_ms=%.3f retrieval_collect_ms=%.3f witness_scan_ms=%.3f union_find_ms=%.3f retrieval_total_ms=%.3f\n",
               loadMs, got.collectMs, got.witnessMs, got.unionMs, got.totalMs);
        printf("retrieval_peak_extra_mib=%.3f (RSS high-water after index load) surviving_witnesses=%zu incidences=%zu\n",
               got.peakExtraMiB, got.survivingWitnesses, got.incidences);
        printf("reference_selected=%zu reference_components=%zu reference_sizes=", refSelected, refSizes.size());
        printSizesRLE(refSizes);
        printf("\ncorrectness_direct_reference=%s reference_ms=%.3f\n", exact ? "EXACT-MATCH" : "MISMATCH", referenceMs);
        return exact ? 0 : 2;
    }

    if (mode == "point") {
        int r = 0, s = 0;
        if (argc < 5 || !parseIntArg(argv[3], r) || !parseIntArg(argv[4], s)) {
            printNSI2Usage(argv[0]); return 1;
        }
        if (!q.column(r)) return printQueryError(QueryCode::ROutOfRange);
        if (argc != 5 + r) { printNSI2Usage(argv[0]); return 1; }
        vector<int> vs((size_t)r);
        for (int i = 0; i < r; ++i) if (!parseIntArg(argv[5+i], vs[i])) { fprintf(stderr, "bad vertex id\n"); return 1; }
        QueryCode code; double ans = q.pointKernel(vs.data(), r, s, code);
        if (code != QueryCode::Ok) return printQueryError(code);
        printf("kappa_{%d,%d} = %.0f\n", r, s, ans); return 0;
    }

    if (mode == "point-validated") {
        int r = 0, s = 0;
        if (argc < 6 || !parseIntArg(argv[4], r) || !parseIntArg(argv[5], s)) {
            printNSI2Usage(argv[0]); return 1;
        }
        if (!q.column(r)) return printQueryError(QueryCode::ROutOfRange);
        if (argc != 6 + r) { printNSI2Usage(argv[0]); return 1; }
        ValidationGraph g;
        if (!loadValidationGraph(argv[3], g, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        if (g.n != x.n) { fprintf(stderr, "validation graph/index vertex-count mismatch\n"); return 1; }
        vector<int> vs((size_t)r);
        for (int i = 0; i < r; ++i) if (!parseIntArg(argv[6+i], vs[i])) { fprintf(stderr, "bad vertex id\n"); return 1; }
        QueryCode code; double ans = q.pointValidated(g, vs.data(), r, s, code);
        if (code != QueryCode::Ok) return printQueryError(code);
        printf("kappa_{%d,%d} = %.0f (clique validated)\n", r, s, ans); return 0;
    }

    if (mode == "row" || mode == "spectrum") {
        int r = 0;
        if (argc < 4 || !parseIntArg(argv[3], r)) { printNSI2Usage(argv[0]); return 1; }
        if (!q.column(r)) return printQueryError(QueryCode::ROutOfRange);
        if (argc != 4 + r) { printNSI2Usage(argv[0]); return 1; }
        vector<int> vs((size_t)r); vector<double> row;
        for (int i = 0; i < r; ++i) if (!parseIntArg(argv[4+i], vs[i])) { fprintf(stderr, "bad vertex id\n"); return 1; }
        QueryCode code = q.rowKernel(vs.data(), r, row);
        if (code != QueryCode::Ok) return printQueryError(code);
        int boundary = q.column(r)->boundary;
        for (size_t i = 0; i < row.size(); ++i) printf("s=%d kappa=%.0f\n", boundary + (int)i, row[i]);
        return 0;
    }

    if (mode == "row-validated") {
        int r = 0;
        if (argc < 5 || !parseIntArg(argv[4], r)) { printNSI2Usage(argv[0]); return 1; }
        if (!q.column(r)) return printQueryError(QueryCode::ROutOfRange);
        if (argc != 5 + r) { printNSI2Usage(argv[0]); return 1; }
        ValidationGraph g;
        if (!loadValidationGraph(argv[3], g, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        if (g.n != x.n) { fprintf(stderr, "validation graph/index vertex-count mismatch\n"); return 1; }
        vector<int> vs((size_t)r); vector<double> row;
        for (int i = 0; i < r; ++i) if (!parseIntArg(argv[5+i], vs[i])) { fprintf(stderr, "bad vertex id\n"); return 1; }
        QueryCode code = q.rowValidated(g, vs.data(), r, row);
        if (code != QueryCode::Ok) return printQueryError(code);
        int boundary = q.column(r)->boundary;
        for (size_t i = 0; i < row.size(); ++i) printf("s=%d kappa=%.0f\n", boundary + (int)i, row[i]);
        return 0;
    }

    if (mode == "count") {
        int r = 0, s = 0;
        if (argc != 6 || !parseIntArg(argv[3], r) || !parseIntArg(argv[4], s)) {
            printNSI2Usage(argv[0]); return 1;
        }
        char *e = nullptr; double k = strtod(argv[5], &e);
        const NSI2Column *col = q.column(r);
        if (!col) return printQueryError(QueryCode::ROutOfRange);
        if (s < col->boundary || s > x.smax) return printQueryError(QueryCode::SOutOfRange);
        if (!e || *e || !isfinite(k) || k <= 0) {
            fprintf(stderr, "K must be positive (zero-core cliques outside the indexed boundary are not materialized)\n");
            return 1;
        }
        double total = 0;
        for (auto &kv : col->dists[s - col->boundary]) if (kv.first >= k) total += kv.second;
        printf("count(kappa_{%d,%d} >= %.0f) = %.0f\n", r, s, k, total); return 0;
    }

    if (mode == "pointfile" || mode == "pointfile-validated") {
        const bool validated = mode == "pointfile-validated";
        int ai = validated ? 4 : 3, r = 0, s = 0;
        const int need = validated ? 7 : 6;
        if (argc != need || !parseIntArg(argv[ai], r) || !parseIntArg(argv[ai+1], s)) {
            printNSI2Usage(argv[0]); return 1;
        }
        const NSI2Column *pointCol = q.column(r);
        if (!pointCol) return printQueryError(QueryCode::ROutOfRange);
        if (s < pointCol->boundary || s > x.smax)
            return printQueryError(QueryCode::SOutOfRange);
        ValidationGraph g;
        if (validated && !loadValidationGraph(argv[3], g, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        if (validated && g.n != x.n) { fprintf(stderr, "validation graph/index vertex-count mismatch\n"); return 1; }
        vector<int> flat;
        if (!readFixedCliques(argv[ai+2], r, flat, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        const size_t nq = r > 0 ? flat.size() / (size_t)r : 0;
        for (size_t i = 0; i < nq; ++i) {
            QueryCode code; const int *vs = &flat[i * (size_t)r];
            double ans = validated ? q.pointValidated(g, vs, r, s, code) : q.pointKernel(vs, r, s, code);
            if (code != QueryCode::Ok) {
                fprintf(stderr, "query %zu failed: %s\n", i, queryCodeName(code)); return 2;
            }
            printf("%.0f\n", ans);
        }
        return 0;
    }

    if (mode == "rowfile" || mode == "rowfile-validated") {
        const bool validated = mode == "rowfile-validated";
        int ai = validated ? 4 : 3, r = 0;
        const int need = validated ? 6 : 5;
        if (argc != need || !parseIntArg(argv[ai], r)) { printNSI2Usage(argv[0]); return 1; }
        if (!q.column(r)) return printQueryError(QueryCode::ROutOfRange);
        ValidationGraph g;
        if (validated && !loadValidationGraph(argv[3], g, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        if (validated && g.n != x.n) { fprintf(stderr, "validation graph/index vertex-count mismatch\n"); return 1; }
        vector<int> flat; vector<double> row;
        if (!readFixedCliques(argv[ai+1], r, flat, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        const size_t nq = r > 0 ? flat.size() / (size_t)r : 0;
        for (size_t i = 0; i < nq; ++i) {
            const int *vs = &flat[i * (size_t)r];
            QueryCode code = validated ? q.rowValidated(g, vs, r, row) : q.rowKernel(vs, r, row);
            if (code != QueryCode::Ok) {
                fprintf(stderr, "query %zu failed: %s\n", i, queryCodeName(code)); return 2;
            }
            for (size_t j = 0; j < row.size(); ++j) printf("%s%.0f", j ? " " : "", row[j]);
            putchar('\n');
        }
        return 0;
    }

    if (mode == "sample") {                             // emit real r-cliques for gating/benching
        if (argc < 5) { fprintf(stderr, "usage: INDEX sample R COUNT [--by-clique]\n"); return 1; }
        int r = atoi(argv[3]); long long want = atoll(argv[4]);
        // Default: uniform over PATTERNS (what the §223 gates use, kept byte-compatible).
        // --by-clique: uniform over R-CLIQUES, i.e. patterns drawn proportionally to mult. That is
        // the honest workload for a latency benchmark, because a real query stream is over cliques,
        // not over the index's internal orbits.
        bool byClique = false;
        for (int i = 5; i < argc; ++i) if (!strcmp(argv[i], "--by-clique")) byClique = true;
        const NSI2Column *colp = q.column(r);
        if (!colp) { fprintf(stderr, "no column for r=%d\n", r); return 1; }
        const NSI2Column &col = *colp;
        if (col.patterns.empty()) { fprintf(stderr, "index has no pattern table (slim); sample from the full index\n"); return 1; }
        vector<vector<int32_t>> clsVerts((size_t)x.nC);
        for (int32_t v = 0; v < x.n; ++v) if (x.classOf[v] >= 0) clsVerts[x.classOf[v]].push_back(v);
        // a column comp is in REP space; expand a rep to any raw class that maps to it
        vector<vector<int32_t>> repRaw((size_t)x.nC);
        for (int32_t c = 0; c < x.nC; ++c) { int32_t rp = col.classRep[c]; if (rp >= 0) repRaw[rp].push_back(c); }
        unsigned seed = 987654321u;
        auto rnd = [&] { seed = seed * 1664525u + 1013904223u; return seed >> 1; };
        vector<double> cum;                                  // mult prefix sums for --by-clique
        if (byClique) {
            cum.reserve(col.patterns.size());
            double acc = 0;
            for (const auto &P : col.patterns) { acc += (double)P.mult; cum.push_back(acc); }
            if (acc <= 0) { fprintf(stderr, "column r=%d has no r-cliques\n", r); return 1; }
        }
        auto pickPattern = [&]() -> size_t {
            if (!byClique) return (size_t)(rnd() % col.patterns.size());
            const double u = (double)rnd() / 2147483648.0 * cum.back();
            return (size_t)(lower_bound(cum.begin(), cum.end(), u) - cum.begin());
        };
        long long made = 0, guard = 0;
        vector<int32_t> out;
        while (made < want && guard++ < want * 60) {
            const auto &P = col.patterns[min(pickPattern(), col.patterns.size() - 1)];
            out.clear(); bool ok = true;
            for (const auto &cm : P.comp) {
                const auto &raws = repRaw[cm.first];
                if (raws.empty()) { ok = false; break; }
                int need = cm.second, tries = 0;
                while (need > 0 && tries++ < 60) {
                    int32_t rc = raws[rnd() % raws.size()];
                    if (clsVerts[rc].empty()) continue;
                    int32_t v = clsVerts[rc][rnd() % clsVerts[rc].size()];
                    if (find(out.begin(), out.end(), v) != out.end()) continue;
                    out.push_back(v); --need;
                }
                if (need > 0) { ok = false; break; }
            }
            if (!ok || (int)out.size() != r) continue;
            sort(out.begin(), out.end());
            for (int i = 0; i < r; ++i) printf("%d%c", out[i], i + 1 == r ? '\n' : ' ');
            ++made;
        }
        fprintf(stderr, "sampled %lld r-cliques for r=%d\n", made, r);
        return 0;
    }
    if (mode == "dbg-cp") {                             // one-shot diagnostic for the slim path
        if (argc < 4 + x.rmin) { fprintf(stderr, "usage: INDEX dbg-cp R V1..VR\n"); return 1; }
        int r = atoi(argv[3]);
        vector<int> vs; for (int i = 0; i < r; ++i) vs.push_back(atoi(argv[4 + i]));
        fprintf(stderr, "slim=%d nC=%d nR=%d\n", (int)x.slim, x.nC, x.nR);
        for (int i = 0; i < r; ++i) {
            int c = vs[i] < x.n ? x.classOf[vs[i]] : -99;
            fprintf(stderr, "  v=%d class=%d |profile|=%zu\n", vs[i], c,
                    (c >= 0 && c < x.nC) ? x.classRegions[c].size() : 0);
        }
        fprintf(stderr, "  cpFromProfiles=%d\n", q.cpFromProfiles(vs.data(), r));
        long long nzRegion = 0; for (int32_t i = 0; i < x.nR; ++i) if (x.regionSize[i] > 0) nzRegion++;
        fprintf(stderr, "  regionSize: nonzero=%lld of %d, max=%d\n", nzRegion, x.nR,
                x.regionSize.empty() ? -1 : *max_element(x.regionSize.begin(), x.regionSize.end()));
        return 0;
    }
    if (mode == "verify-cp") {                          // §222: is cP recoverable from class profiles?
        if (argc < 4) { fprintf(stderr, "usage: INDEX verify-cp QUERYFILE\n"); return 1; }
        FILE *qf = fopen(argv[3], "r");
        if (!qf) { fprintf(stderr, "cannot open %s\n", argv[3]); return 1; }
        vector<int> flat; int vtmp;
        while (fscanf(qf, "%d", &vtmp) == 1) flat.push_back(vtmp);
        fclose(qf);
        // r is inferred from the workload itself: feeding an r-tuple to a column with a different r
        // tests non-cliques and produces meaningless mismatches (learned the hard way).
        int qr = 0;
        { FILE *g2 = fopen(argv[3], "r"); char line[4096];
          if (g2 && fgets(line, sizeof line, g2)) { char *tk = strtok(line, " \t\n"); while (tk) { qr++; tk = strtok(nullptr, " \t\n"); } }
          if (g2) fclose(g2); }
        if (qr <= 0) { fprintf(stderr, "cannot infer r from the workload\n"); return 1; }
        fprintf(stderr, "verify-cp: workload r=%d, testing only that column\n", qr);
        long long okPat = 0, badPat = 0, okMerg = 0, badMerg = 0, none = 0, fullCert = 0, partial = 0;
        for (int32_t ci = 0; ci < (int32_t)x.columns.size(); ++ci) {
            const NSI2Column &col = x.columns[ci];
            const int r = col.r;
            if (r != qr) continue;
            const long long nq = (long long)flat.size() / r;
            for (long long i = 0; i < nq; ++i) {
                const int *vs = &flat[i * r];
                bool bad = false;
                for (int j = 0; j < r; ++j) if (vs[j] < 0 || vs[j] >= x.n) bad = true;
                if (bad) continue;
                int pi = q.lookupPattern(col, vs);
                int32_t cpComputed = q.cpFromProfiles(vs, r);
                if (pi >= 0) {
                    const auto &P = col.patterns[pi];
                    if (P.closedFrom >= 0 && P.closedFrom <= col.boundary) fullCert++; else partial++;
                    if (cpComputed == P.cP) okPat++;
                    else { if (badPat < 5) fprintf(stderr, "  MISMATCH r=%d stored cP=%d computed=%d\n", r, P.cP, cpComputed); badPat++; }
                } else {
                    int mg = q.lookupMergeable(col, vs);
                    if (mg >= 0) { if (cpComputed == x.regionSize[mg]) okMerg++; else badMerg++; }
                    else { if (cpComputed < 0) none++; else badMerg++; }
                }
            }
        }
        // PHASE 2: the mergeable path is never reached by pattern-table samples, so generate cliques
        // that live INSIDE mergeable regions and check cP == regionSize there too. In the slim design
        // both paths collapse into the same closed form, so both must be verified.
        {
            vector<vector<int32_t>> clsVerts((size_t)x.nC);
            for (int32_t v = 0; v < x.n; ++v) if (x.classOf[v] >= 0) clsVerts[x.classOf[v]].push_back(v);
            long long mOk = 0, mBad = 0, mTried = 0;
            for (const auto &col : x.columns) {
                const int r = col.r;
                int taken = 0;
                for (int32_t rid : col.mergeableRegions) {
                    if (taken >= 200) break;
                    vector<int> verts;
                    for (int32_t c = 0; c < x.nC && (int)verts.size() < r; ++c) {
                        const auto &pr = x.classRegions[c];
                        if (!binary_search(pr.begin(), pr.end(), rid)) continue;
                        for (int32_t v : clsVerts[c]) { verts.push_back(v); if ((int)verts.size() >= r) break; }
                    }
                    if ((int)verts.size() < r) continue;
                    ++taken; ++mTried;
                    int32_t cpc = q.cpFromProfiles(verts.data(), r);
                    int32_t want = x.regionSize[rid];
                    // cP is the LARGEST hosting region, which may exceed this one; the check is that
                    // the computed value is a real hosting-region size and at least this region's.
                    if (cpc >= want) ++mOk;
                    else { if (mBad < 5) fprintf(stderr, "  MERG MISMATCH rid=%d size=%d computed=%d\n", rid, want, cpc); ++mBad; }
                }
            }
            printf("verify-cp: mergeable-region cliques tried=%lld ok=%lld MISMATCH=%lld\n", mTried, mOk, mBad);
            if (mBad) badMerg += mBad;
        }
        printf("verify-cp: patterns ok=%lld MISMATCH=%lld | mergeable ok=%lld MISMATCH=%lld | unhosted=%lld\n",
               okPat, badPat, okMerg, badMerg, none);
        printf("verify-cp: of matched patterns, fully-certified=%lld partial=%lld (%.2f%% need no record)\n",
               fullCert, partial, (fullCert + partial) ? 100.0 * fullCert / (fullCert + partial) : 0.0);
        return (badPat || badMerg) ? 2 : 0;
    }
    // ===== §226: byte anatomy of the LOADED index, and what a packed encoding would cost =====
    // §221 did this for NSI2 and it is what identified the pattern table as the only worthwhile
    // target; the same discipline decides whether slimming NSI3 further is worth building at all.
    // This mode WRITES NOTHING: it re-derives each block's byte cost from the loaded structures and
    // prices a packed alternative (bit-packed class labels, delta-varint id lists, varint integers),
    // so the payoff is known before any format change is implemented.
    if (mode == "anatomy") {
        auto vlen = [](uint64_t v) { int b = 1; while (v >= 128) { v >>= 7; ++b; } return b; };
        auto zlen = [&](int64_t v) { return vlen(((uint64_t)v << 1) ^ (uint64_t)(v >> 63)); };
        auto ilen = [&](double d, int wide) {                       // integral -> varint, else raw 8B
            double r0 = nearbyint(d);
            return (fabs(d - r0) < 1e-9 && fabs(r0) < 9e18) ? zlen((int64_t)r0) : wide;
        };
        int labelBits = 1; while ((1LL << labelBits) < (long long)x.nC + 1) ++labelBits;

        struct Blk { const char *name; double now, packed; };
        vector<Blk> blk;
        double nonIntegral = 0, integralTotal = 0;

        // A vertex whose class hosts no region can never appear in an answer, so its label is dead
        // weight. Price the alternative: a presence bitmap plus packed labels for the live ones only.
        long long liveV = 0;
        for (int32_t v = 0; v < x.n; ++v) {
            const int32_t c = x.classOf[v];
            if (c >= 0 && c < x.nC && !x.classRegions[c].empty()) ++liveV;
        }
        int liveBits = 1; while ((1LL << liveBits) < (long long)x.nC + 1) ++liveBits;
        const double classOfPacked = ceil((double)x.n * labelBits / 8.0);
        const double classOfBitmap = ceil((double)x.n / 8.0) + ceil((double)liveV * liveBits / 8.0);
        blk.push_back({"classOf", 4.0 * x.n, min(classOfPacked, classOfBitmap)});

        double profNow = 0, profPk = 0; long long profEnt = 0;
        for (int32_t c = 0; c < x.nC; ++c) {
            const auto &pr = x.classRegions[c];
            profNow += 8.0 + 4.0 * pr.size();
            profPk += vlen((uint64_t)max(0, x.classSize.empty() ? 0 : x.classSize[c])) + vlen(pr.size());
            int32_t prev = -1;
            for (int32_t rid : pr) { profPk += zlen((int64_t)rid - prev); prev = rid; }
            profEnt += (long long)pr.size();
        }
        blk.push_back({"class-profiles", profNow, profPk});

        double mgNow = 0, mgPk = 0, exNow = 0, exPk = 0, rsNow = 0, rsPk = 0, dsNow = 0, dsPk = 0;
        long long nExc = 0, nRes = 0, nDist = 0;
        for (const auto &col : x.columns) {
            mgNow += 8.0 + 4.0 + 4.0 * col.mergeableRegions.size();
            vector<int32_t> mg = col.mergeableRegions;
            sort(mg.begin(), mg.end());
            mgPk += 2.0 + vlen(mg.size());
            int32_t prev = -1;
            for (int32_t rid : mg) { mgPk += zlen((int64_t)rid - prev); prev = rid; }

            exNow += 4.0; exPk += vlen(col.patterns.size());
            int32_t prevPid = -1;
            for (size_t pi = 0; pi < col.patterns.size(); ++pi) {
                const auto &P = col.patterns[pi];
                const int32_t opid = x.slim && pi < col.origPid.size() ? col.origPid[pi] : (int32_t)pi;
                exNow += 4.0 + 1.0 + 6.0 * P.comp.size() + 8.0 + 4.0;
                exPk += zlen((int64_t)opid - prevPid) + 1.0; prevPid = opid;
                for (const auto &cm : P.comp) exPk += vlen((uint64_t)cm.first) + vlen((uint64_t)cm.second);
                exPk += ilen(P.boundaryCore, 8) + vlen((uint64_t)(P.closedFrom + 1));
                integralTotal += 1; if (ilen(P.boundaryCore, 8) == 8) nonIntegral += 1;
                ++nExc;
            }
            rsNow += 4.0; rsPk += vlen(col.resid.size());
            for (const auto &rr : col.resid) {
                rsNow += 8.0 + 12.0 * rr.size();
                rsPk += vlen(rr.size());
                int32_t pp = -1;
                for (const auto &pc : rr) {
                    rsPk += zlen((int64_t)pc.first - pp) + ilen(pc.second, 8); pp = pc.first;
                    integralTotal += 1; if (ilen(pc.second, 8) == 8) nonIntegral += 1;
                }
                nRes += (long long)rr.size();
            }
            dsNow += 4.0; dsPk += vlen(col.dists.size());
            for (const auto &d : col.dists) {
                dsNow += 4.0 + 16.0 * d.size();
                dsPk += vlen(d.size());
                for (const auto &kv : d) {
                    dsPk += ilen(kv.first, 8) + ilen(kv.second, 8);
                    integralTotal += 2;
                    if (ilen(kv.first, 8) == 8) nonIntegral += 1;
                    if (ilen(kv.second, 8) == 8) nonIntegral += 1;
                }
                nDist += (long long)d.size();
            }
        }
        blk.push_back({"mergeables", mgNow, mgPk});
        blk.push_back({"exceptions", exNow, exPk});
        blk.push_back({"residue", rsNow, rsPk});
        blk.push_back({"dists", dsNow, dsPk});
        const double hdr = 36.0 + 8.0 * x.nCols;
        blk.push_back({"header+dir", hdr, hdr});

        double now = 0, packed = 0;
        for (const auto &b : blk) { now += b.now; packed += b.packed; }
        printf("anatomy index=%s %s file=%.3f MB accounted=%.3f MB (delta %.0f B)\n",
               argv[1], x.slim ? "NSI3" : "NSI2", x.fileBytes / 1048576.0, now / 1048576.0,
               (double)x.fileBytes - now);
        printf("n=%d nC=%d nR=%d class-label-bits=%d profile-entries=%lld exceptions=%lld "
               "residue=%lld dist-entries=%lld\n",
               x.n, x.nC, x.nR, labelBits, profEnt, nExc, nRes, nDist);
        printf("block            now-B    now%%    packed-B  packed%%   shrink\n");
        for (const auto &b : blk)
            printf("%-14s %10.0f  %5.1f%%  %10.0f   %5.1f%%   %6.2fx\n", b.name, b.now,
                   now > 0 ? 100.0 * b.now / now : 0.0, b.packed,
                   packed > 0 ? 100.0 * b.packed / packed : 0.0,
                   b.packed > 0 ? b.now / b.packed : 0.0);
        printf("%-14s %10.0f  100.0%%  %10.0f   100.0%%   %6.2fx   (%.3f MB -> %.3f MB)\n",
               "TOTAL", now, packed, packed > 0 ? now / packed : 0.0,
               now / 1048576.0, packed / 1048576.0);
        printf("non-integral stored reals: %.0f of %.0f (%.4f%%) -- these are the only values a varint "
               "cannot take\n", nonIntegral, integralTotal,
               integralTotal > 0 ? 100.0 * nonIntegral / integralTotal : 0.0);
        printf("class labels: %lld of %d vertices are in a class that hosts a region (%.1f%%); "
               "flat-packed=%.0f B bitmap+live=%.0f B (chose %s)\n",
               liveV, x.n, x.n ? 100.0 * (double)liveV / x.n : 0.0, classOfPacked, classOfBitmap,
               classOfBitmap < classOfPacked ? "bitmap" : "flat");
        return 0;
    }

    // ===== §225 E4: the MATERIALIZED-ARCHIVE baseline =====
    // The archive is the only structure that answers the same queries without this index: for every
    // r-clique, kappa in every cell of its row, sorted by the clique key and probed by binary search.
    // EVERY choice below is deliberately made in the ARCHIVE's favour, so the reported gap is a
    // lower bound on the real one:
    //   * rows = SUM over patterns of mult. Distinct patterns own disjoint clique sets (a clique
    //     determines its rep-multiset), so this is a strict LOWER bound on #r-cliques. Cliques that
    //     live only in a mergeable region are reported separately and NOT charged to the archive,
    //     because they can overlap the pattern side and charging them would inflate the baseline.
    //   * 4-byte vertex ids, 4-byte kappa, and ONE key shared by the whole row's cells.
    //   * the probe workload comes from the pattern table, so every probe HITS: no early exit.
    // Where the row count exceeds the cap the archive cannot be built at all, and saying so IS the
    // result: that is the regime the index exists for.
    if (mode == "archive" || mode == "archive-bench") {
        if (x.slim) {
            fprintf(stderr,
                "archive accounting needs the FULL plane index (NSI2). The slim format drops pattern\n"
                "multiplicities on purpose, so #r-cliques cannot be read back from it. Run this mode on\n"
                "the NSI2 file and compare its output against the NSI3 file size.\n");
            return 1;
        }
        auto nCrD = [](double a, int b) {
            if (b < 0 || (double)b > a) return 0.0;
            if ((double)b > a - b) b = (int)(a - (double)b);
            double z = 1.0;
            for (int i = 0; i < b; ++i) z = z * (a - i) / (double)(i + 1);
            return z;
        };

        vector<vector<int32_t>> clsVerts((size_t)x.nC);
        for (int32_t v = 0; v < x.n; ++v) if (x.classOf[v] >= 0) clsVerts[x.classOf[v]].push_back(v);

        double archBytes = 0, patRCall = 0, mergRCall = 0;
        long long multMismatch = 0;
        printf("archive index=%s plane r=%d..%d s<=%d\n", argv[1], x.rmin, x.rmax, x.smax);
        printf("r  cells  patterns  pattern-r-cliques  mergeable-r-cliques  row-bytes  archive-MB\n");
        for (const auto &col : x.columns) {
            const int cells = x.smax - col.boundary + 1;
            vector<double> repN((size_t)x.nC, 0.0);
            for (int32_t c = 0; c < x.nC; ++c) {
                int32_t rp = col.classRep[c];
                if (rp >= 0) repN[rp] += (double)clsVerts[c].size();
            }
            double patRC = 0;
            for (const auto &P : col.patterns) {
                patRC += (double)P.mult;
                double m = 1.0;
                for (const auto &cm : P.comp) m *= nCrD(repN[cm.first], cm.second);
                if (fabs(m - (double)P.mult) > 1e-6 * max(1.0, (double)P.mult)) ++multMismatch;
            }
            double mergRC = 0;
            for (int32_t rid : col.mergeableRegions) mergRC += nCrD((double)x.regionSize[rid], col.r);
            const double rowBytes = 4.0 * col.r + 4.0 * cells;
            const double b = patRC * rowBytes;
            printf("%d  %d  %zu  %.0f  %.0f  %.0f  %.2f\n", col.r, cells, col.patterns.size(),
                   patRC, mergRC, rowBytes, b / 1048576.0);
            archBytes += b; patRCall += patRC; mergRCall += mergRC;
        }
        printf("ARCHIVE   %.0f rows  %.2f MB  (%.2f GB)\n",
               patRCall, archBytes / 1048576.0, archBytes / 1073741824.0);
        printf("INDEX     %s  %.2f MB\n", argv[1], x.fileBytes / 1048576.0);
        printf("RATIO     archive / this-index = %.1fx\n",
               x.fileBytes > 0 ? archBytes / (double)x.fileBytes : 0.0);
        printf("bytes-per-r-clique  archive=%.2f  this-index=%.4f\n",
               patRCall > 0 ? archBytes / patRCall : 0.0,
               patRCall > 0 ? (double)x.fileBytes / patRCall : 0.0);
        // an optional second index path (typically the NSI3 slim file) closes the comparison here so
        // no shell arithmetic sits between the measurement and the reported ratio
        for (int i = 3; i < argc; ++i) {
            if (strcmp(argv[i], "--vs") || i + 1 >= argc) continue;
            FILE *sf = fopen(argv[++i], "rb");
            if (!sf) { fprintf(stderr, "cannot open --vs %s\n", argv[i]); return 1; }
            fseek(sf, 0, SEEK_END); double sb = (double)ftell(sf); fclose(sf);
            printf("SLIM      %s  %.2f MB\n", argv[i], sb / 1048576.0);
            printf("RATIO     archive / slim-index = %.1fx   bytes-per-r-clique slim=%.4f\n",
                   sb > 0 ? archBytes / sb : 0.0, patRCall > 0 ? sb / patRCall : 0.0);
        }
        if (multMismatch)
            printf("WARNING   %lld patterns whose stored mult disagrees with prod C(|rep|,b);"
                   " expansion below emits the enumerated count, not mult\n", multMismatch);
        if (mode == "archive") return 0;

        // ---- materialize one row of the archive and probe it ----
        if (argc < 5) { fprintf(stderr, "usage: INDEX archive-bench R QUERYFILE [--cap-gb N] [--cold-mib N]\n"); return 1; }
        int rWant = 0;
        if (!parseIntArg(argv[3], rWant)) { fprintf(stderr, "bad R\n"); return 1; }
        double capGB = 64.0; int coldMiB = 128;
        for (int i = 5; i < argc; ++i) {
            if (!strcmp(argv[i], "--cap-gb") && i + 1 < argc) capGB = atof(argv[++i]);
            else if (!strcmp(argv[i], "--cold-mib") && i + 1 < argc) {
                if (!parseIntArg(argv[++i], coldMiB) || coldMiB < 1 || coldMiB > 4096) { fprintf(stderr, "bad --cold-mib\n"); return 1; }
            } else if (!strcmp(argv[i], "--vs") && i + 1 < argc) ++i;
            else { fprintf(stderr, "unknown archive-bench option: %s\n", argv[i]); return 1; }
        }
        const NSI2Column *colp = q.column(rWant);
        if (!colp) { fprintf(stderr, "no column for r=%d\n", rWant); return 1; }
        const NSI2Column &col = *colp;
        const int cells = x.smax - col.boundary + 1;
        double rows = 0;
        for (const auto &P : col.patterns) rows += (double)P.mult;
        // peak = keys + values + the permutation + the sorted copy of both
        const double peakBytes = rows * (4.0 * rWant + 4.0 * cells) * 2.0 + rows * 4.0;
        if (peakBytes > capGB * 1073741824.0) {
            printf("ARCHIVE-BENCH r=%d **CANNOT BE MATERIALIZED**: %.0f rows would need %.1f GB "
                   "(cap %.1f GB)\n", rWant, rows, peakBytes / 1073741824.0, capGB);
            printf("           No sorted table exists to probe, so the only alternative to the index\n"
                   "           is recomputing the cell per query. This outcome is the result.\n");
            return 0;
        }

        vector<vector<int32_t>> repV((size_t)x.nC);
        for (int32_t c = 0; c < x.nC; ++c) {
            int32_t rp = col.classRep[c];
            if (rp < 0) continue;
            auto &dst = repV[rp];
            dst.insert(dst.end(), clsVerts[c].begin(), clsVerts[c].end());
        }

        auto expand0 = SteadyClock::now();
        vector<int32_t> keys, vals;
        keys.reserve((size_t)rows * rWant); vals.reserve((size_t)rows * cells);
        // `cur` holds the tuple as the recursion built it, one block per comp entry. The emitted key
        // must be sorted, but sorting `cur` IN PLACE permutes the blocks the OUTER levels still own,
        // so every later tuple comes out corrupted. Sort into a separate buffer.
        vector<int32_t> cur((size_t)rWant), keyBuf((size_t)rWant), rowKappa((size_t)cells);
        for (size_t pi = 0; pi < col.patterns.size(); ++pi) {
            Query2::Resolved res; res.col = &col; res.pattern = (int32_t)pi;
            bool bad = false;
            for (int s = col.boundary; s <= x.smax; ++s) {
                QueryCode code = QueryCode::Ok;
                double kv = q.value(res, s, code);
                if (code != QueryCode::Ok) { bad = true; break; }
                rowKappa[s - col.boundary] = (int32_t)kv;
            }
            if (bad) { fprintf(stderr, "pattern %zu has no value at some cell\n", pi); return 2; }
            const Comp2 &comp = col.patterns[pi].comp;
            auto rec = [&](auto &&self, size_t ci, int filled) -> void {
                if (ci == comp.size()) {
                    copy(cur.begin(), cur.begin() + filled, keyBuf.begin());
                    sort(keyBuf.begin(), keyBuf.begin() + filled);
                    keys.insert(keys.end(), keyBuf.begin(), keyBuf.begin() + filled);
                    vals.insert(vals.end(), rowKappa.begin(), rowKappa.end());
                    return;
                }
                const auto &V = repV[comp[ci].first];
                const int b = comp[ci].second, m = (int)V.size();
                if (b > m) return;
                vector<int> idx((size_t)b);
                for (int i = 0; i < b; ++i) idx[i] = i;
                while (true) {
                    for (int i = 0; i < b; ++i) cur[filled + i] = V[idx[i]];
                    self(self, ci + 1, filled + b);
                    int i = b - 1;
                    while (i >= 0 && idx[i] == m - b + i) --i;
                    if (i < 0) break;
                    ++idx[i];
                    for (int j = i + 1; j < b; ++j) idx[j] = idx[j - 1] + 1;
                }
            };
            rec(rec, 0, 0);
        }
        const size_t M = vals.size() / (size_t)cells;
        const double expandS = chrono::duration<double>(SteadyClock::now() - expand0).count();

        auto sort0 = SteadyClock::now();
        vector<uint32_t> ord(M);
        for (size_t i = 0; i < M; ++i) ord[i] = (uint32_t)i;
        sort(ord.begin(), ord.end(), [&](uint32_t a, uint32_t b) {
            const int32_t *ka = &keys[(size_t)a * rWant], *kb = &keys[(size_t)b * rWant];
            for (int i = 0; i < rWant; ++i) if (ka[i] != kb[i]) return ka[i] < kb[i];
            return false;
        });
        vector<int32_t> sk((size_t)M * rWant), sv((size_t)M * cells);
        for (size_t i = 0; i < M; ++i) {
            memcpy(&sk[i * rWant], &keys[(size_t)ord[i] * rWant], sizeof(int32_t) * rWant);
            memcpy(&sv[i * cells], &vals[(size_t)ord[i] * cells], sizeof(int32_t) * cells);
        }
        keys.clear(); keys.shrink_to_fit(); vals.clear(); vals.shrink_to_fit();
        ord.clear(); ord.shrink_to_fit();
        const double sortS = chrono::duration<double>(SteadyClock::now() - sort0).count();
        const double realBytes = (double)M * (4.0 * rWant + 4.0 * cells);
        printf("ARCHIVE-BENCH r=%d materialized rows=%zu cells=%d  %.2f MB  build=%.1fs "
               "(expand %.1fs + sort %.1fs)  peak-rss=%.0f MiB\n",
               rWant, M, cells, realBytes / 1048576.0, expandS + sortS, expandS, sortS, peakRssMiB());
        if (M == 0) { fprintf(stderr, "empty archive\n"); return 2; }

        vector<BenchQuery2> allQ; vector<int> flat;
        if (!readBenchQueries(argv[4], x, allQ, flat, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        vector<BenchQuery2> queries;
        for (const auto &b : allQ) if (b.r == rWant) queries.push_back(b);
        if (queries.size() < 1000) {
            fprintf(stderr, "archive-bench needs >=1000 queries at r=%d (got %zu)\n", rWant, queries.size());
            return 1;
        }

        vector<int32_t> probe((size_t)rWant);
        long long hits = 0, misses = 0;
        auto probeOp = [&](size_t i) -> double {
            const auto &b = queries[i];
            memcpy(probe.data(), &flat[b.off], sizeof(int32_t) * rWant);
            sort(probe.begin(), probe.end());
            size_t lo = 0, hi = M;
            while (lo < hi) {
                size_t mid = lo + ((hi - lo) >> 1);
                const int32_t *km = &sk[mid * rWant];
                int c = 0;
                for (int j = 0; j < rWant; ++j) if (km[j] != probe[j]) { c = km[j] < probe[j] ? -1 : 1; break; }
                if (c < 0) lo = mid + 1; else hi = mid;
            }
            bool eq = lo < M;
            if (eq) { const int32_t *km = &sk[lo * rWant];
                      for (int j = 0; j < rWant; ++j) if (km[j] != probe[j]) { eq = false; break; } }
            if (!eq) {
                if (misses < 3) {                       // a miss means the workload and the archive
                    fprintf(stderr, "  archive MISS:");  // disagree about what an r-clique is
                    for (int j = 0; j < rWant; ++j) fprintf(stderr, " %d", probe[j]);
                    fprintf(stderr, "  (pattern=%d mergeable=%d)\n",
                            q.lookupPattern(col, probe.data()), q.lookupMergeable(col, probe.data()));
                }
                ++misses; return 0.0;
            }
            ++hits;
            return (double)sv[lo * cells + (b.s - col.boundary)];       // point query: one cell
        };
        auto rowOp = [&](size_t i) -> double {
            const auto &b = queries[i];
            memcpy(probe.data(), &flat[b.off], sizeof(int32_t) * rWant);
            sort(probe.begin(), probe.end());
            size_t lo = 0, hi = M;
            while (lo < hi) {
                size_t mid = lo + ((hi - lo) >> 1);
                const int32_t *km = &sk[mid * rWant];
                int c = 0;
                for (int j = 0; j < rWant; ++j) if (km[j] != probe[j]) { c = km[j] < probe[j] ? -1 : 1; break; }
                if (c < 0) lo = mid + 1; else hi = mid;
            }
            if (lo >= M) return 0.0;
            const int32_t *km = &sk[lo * rWant];
            for (int j = 0; j < rWant; ++j) if (km[j] != probe[j]) return 0.0;
            double z = 0;
            for (int cIdx = 0; cIdx < cells; ++cIdx) z += sv[lo * cells + cIdx];   // whole row
            return z;
        };

        const size_t warmN = queries.size(), coldN = min<size_t>(queries.size(), 1000);
        for (size_t i = 0; i < warmN; ++i) latencySink += probeOp(i) + rowOp(i);
        const long long hitsWarm = hits, missWarm = misses;
        vector<uint64_t> eviction((size_t)coldMiB * 1024 * 1024 / sizeof(uint64_t));
        for (size_t i = 0; i < eviction.size(); ++i) eviction[i] = i * 0x9E3779B97F4A7C15ULL + 1;
        LatencySummary pW = measureLatency(warmN, 1, nullptr, probeOp);
        LatencySummary rW = measureLatency(warmN, 1, nullptr, rowOp);
        LatencySummary pC = measureLatency(coldN, 1, &eviction, probeOp);
        LatencySummary rC = measureLatency(coldN, 1, &eviction, rowOp);
        printf("ARCHIVE-BENCH probes=%zu hits=%lld misses=%lld hit-rate=%.4f "
               "(a miss is an unfairly CHEAP probe for the archive)\n",
               warmN, hitsWarm, missWarm,
               (hitsWarm + missWarm) ? (double)hitsWarm / (double)(hitsWarm + missWarm) : 0.0);
        printf("operation  path       warm-med(ns)  warm-p95(ns)  cold-med(ns)  cold-p95(ns)\n");
        printf("%-9s  %-9s  %12.1f  %12.1f  %12.1f  %12.1f\n", "point", "archive", pW.median, pW.p95, pC.median, pC.p95);
        printf("%-9s  %-9s  %12.1f  %12.1f  %12.1f  %12.1f\n", "row", "archive", rW.median, rW.p95, rC.median, rC.p95);
        printf("sink=%.0f\n", (double)latencySink);
        return 0;
    }

    if (mode == "bench") {
        if (argc < 5) { printNSI2Usage(argv[0]); return 1; }
        int coldMiB = 128, warmReps = 1;
        for (int i = 5; i < argc; ++i) {
            if (!strcmp(argv[i], "--cold-mib") && i + 1 < argc) {
                if (!parseIntArg(argv[++i], coldMiB) || coldMiB < 1 || coldMiB > 4096) { fprintf(stderr, "bad --cold-mib\n"); return 1; }
            } else if (!strcmp(argv[i], "--warm-reps") && i + 1 < argc) {
                if (!parseIntArg(argv[++i], warmReps) || warmReps < 1 || warmReps > 100000) { fprintf(stderr, "bad --warm-reps\n"); return 1; }
            } else { fprintf(stderr, "unknown bench option: %s\n", argv[i]); return 1; }
        }
        ValidationGraph g;
        auto graph0 = SteadyClock::now();
        if (!loadValidationGraph(argv[3], g, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        auto graph1 = SteadyClock::now();
        if (g.n != x.n) { fprintf(stderr, "validation graph/index vertex-count mismatch\n"); return 1; }
        vector<BenchQuery2> queries; vector<int> flat;
        if (!readBenchQueries(argv[4], x, queries, flat, error)) { fprintf(stderr, "%s\n", error.c_str()); return 1; }
        if (queries.size() < 1000) { fprintf(stderr, "bench requires at least 1000 queries (got %zu)\n", queries.size()); return 1; }

        vector<double> row;
        for (size_t i = 0; i < queries.size(); ++i) {
            const auto &b = queries[i]; const int *vs = &flat[b.off]; QueryCode code;
            (void)q.pointValidated(g, vs, b.r, b.s, code);
            if (code != QueryCode::Ok) {
                fprintf(stderr, "benchmark query %zu is invalid: %s\n", i, queryCodeName(code)); return 2;
            }
            code = q.rowKernel(vs, b.r, row);
            if (code != QueryCode::Ok) {
                fprintf(stderr, "benchmark query %zu hits %s\n", i, queryCodeName(code)); return 2;
            }
        }

        const size_t warmN = queries.size();
        const size_t coldN = min<size_t>(queries.size(), 1000);
        vector<uint64_t> eviction((size_t)coldMiB * 1024 * 1024 / sizeof(uint64_t));
        for (size_t i = 0; i < eviction.size(); ++i) eviction[i] = i * 0x9E3779B97F4A7C15ULL + 1;
        bool failed = false;
        auto pointKernelOp = [&](size_t i) {
            const auto &b = queries[i]; QueryCode code;
            double z = q.pointKernel(&flat[b.off], b.r, b.s, code); failed |= code != QueryCode::Ok; return z;
        };
        auto pointValidOp = [&](size_t i) {
            const auto &b = queries[i]; QueryCode code;
            double z = q.pointValidated(g, &flat[b.off], b.r, b.s, code); failed |= code != QueryCode::Ok; return z;
        };
        auto rowKernelOp = [&](size_t i) {
            const auto &b = queries[i]; QueryCode code = q.rowKernel(&flat[b.off], b.r, row);
            failed |= code != QueryCode::Ok; return row.empty() ? 0.0 : row.back();
        };
        auto rowValidOp = [&](size_t i) {
            const auto &b = queries[i]; QueryCode code = q.rowValidated(g, &flat[b.off], b.r, row);
            failed |= code != QueryCode::Ok; return row.empty() ? 0.0 : row.back();
        };

        // One full, untimed pass establishes the warm working set.
        for (size_t i = 0; i < warmN; ++i) {
            latencySink += pointKernelOp(i) + pointValidOp(i) + rowKernelOp(i) + rowValidOp(i);
        }
        vector<double> emptyIntervals; emptyIntervals.reserve(warmN);
        for (size_t i = 0; i < warmN; ++i) {
            auto a = SteadyClock::now(); auto b = SteadyClock::now();
            emptyIntervals.push_back(chrono::duration<double,nano>(b - a).count());
        }
        LatencySummary timerCost = summarizeLatency(std::move(emptyIntervals));
        LatencySummary pkW = measureLatency(warmN, warmReps, nullptr, pointKernelOp);
        LatencySummary pvW = measureLatency(warmN, warmReps, nullptr, pointValidOp);
        LatencySummary rkW = measureLatency(warmN, warmReps, nullptr, rowKernelOp);
        LatencySummary rvW = measureLatency(warmN, warmReps, nullptr, rowValidOp);
        LatencySummary pkC = measureLatency(coldN, 1, &eviction, pointKernelOp);
        LatencySummary pvC = measureLatency(coldN, 1, &eviction, pointValidOp);
        LatencySummary rkC = measureLatency(coldN, 1, &eviction, rowKernelOp);
        LatencySummary rvC = measureLatency(coldN, 1, &eviction, rowValidOp);
        if (failed) { fprintf(stderr, "query failure during timed benchmark\n"); return 2; }

        printf("benchmark index-load=%.3f ms graph-load=%.3f ms queries=%zu warm-reps=%d cold-samples=%zu\n",
               loadMs, chrono::duration<double,milli>(graph1 - graph0).count(), warmN, warmReps, coldN);
        printf("method: warm timings amortize each steady-state call over %d repetitions; cold timings are one call\n", warmReps);
        printf("after an untimed %d MiB cache-eviction sweep (CPU-cache-cold, not disk-cold; eviction time excluded).\n", coldMiB);
        printf("row materializes all admissible s answers into a reused output buffer; validated includes O(r^2) adjacency checks.\n");
        printf("empty steady_clock interval: median=%.1f ns p95=%.1f ns (not subtracted)\n",
               timerCost.median, timerCost.p95);
        printf("operation  path       warm-med(ns)  warm-p95(ns)  cold-med(ns)  cold-p95(ns)\n");
        auto line = [&](const char *op, const char *path, const LatencySummary &w, const LatencySummary &c) {
            printf("%-9s  %-9s  %12.1f  %12.1f  %12.1f  %12.1f\n", op, path,
                   w.median, w.p95, c.median, c.p95);
        };
        line("point", "kernel", pkW, pkC); line("point", "validated", pvW, pvC);
        line("row", "kernel", rkW, rkC); line("row", "validated", rvW, rvC);
        printf("sink=%.0f\n", (double)latencySink);
        return 0;
    }

    printNSI2Usage(argv[0]); return 1;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s INDEX MODE ...\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    char magic[4] = {};
    bool got = fread(magic, 1, 4, f) == 4;
    fclose(f);
    if (!got) { fprintf(stderr, "%s is too short to be an NSI index\n", argv[1]); return 1; }
    if (!memcmp(magic, "NSI1", 4)) return mainNSI1(argc, argv);
    if (!memcmp(magic, "NSI2", 4)) return mainNSI2(argc, argv);
    if (!memcmp(magic, "NSI3", 4)) return mainNSI2(argc, argv);   // §222: slim, same query layer
    if (!memcmp(magic, "NSI4", 4)) return mainNSI2(argc, argv);   // §227: packed slim, same layer
    fprintf(stderr, "%s has unknown index magic\n", argv[1]);
    return 1;
}
