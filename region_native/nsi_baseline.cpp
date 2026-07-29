// E4 BASELINE: the materialized archive that an index must beat.
//
// The archive stores kappa for EVERY r-clique in EVERY cell of the row, sorted by the r-clique key,
// probed by binary search. That is the only alternative to our index that answers queries in
// sub-microsecond time, so it is the honest comparison target on BOTH axes:
//   SIZE  : archive = #r-cliques * (r*4 key bytes + 4 bytes per cell)   vs   our NSI1 file
//   QUERY : binary search over the archive                              vs   nsi_query's kappa()
//
// Modes:
//   size  <index.nsi>                 exact archive size arithmetic, no expansion (always works)
//   bench <index.nsi> <queries> [cap] expand, build, probe with the same workload as nsi_query bench
//
// The expansion is bounded: past <cap> r-cliques we report INFEASIBLE and print what the archive
// would have cost. That outcome is a RESULT, not a failure -- it is exactly the regime where the
// archive cannot be written and the index is the only option.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
using namespace std;
using Clock = chrono::steady_clock;

struct NSI1 {
    int r = 0, s0 = 0, smax = 0, n = 0, nC = 0, nPats = 0, nMerg = 0;
    vector<int32_t> classOf, mergSize;
    vector<vector<int32_t>> mergV;                    // mergeable region -> vertex list
    vector<pair<int32_t, uint8_t>> compFlat;          // (class, mult) pairs
    vector<int32_t> compOff, cP;
    vector<double> k0;
    vector<vector<pair<int32_t, double>>> resid;
    long long bytes = 0;
};

static bool load(const char *path, NSI1 &x) {
    FILE *f = fopen(path, "rb");
    if (!f) return false;
    char magic[4];
    if (fread(magic, 4, 1, f) != 1 || memcmp(magic, "NSI1", 4)) { fclose(f); return false; }
    auto r32 = [&] { int32_t v; if (fread(&v, 4, 1, f) != 1) exit(3); return v; };
    auto r64 = [&] { int64_t v; if (fread(&v, 8, 1, f) != 1) exit(3); return v; };
    auto rd  = [&] { double v;  if (fread(&v, 8, 1, f) != 1) exit(3); return v; };
    auto r8  = [&] { uint8_t v; if (fread(&v, 1, 1, f) != 1) exit(3); return v; };
    x.r = r32(); x.s0 = r32(); x.smax = r32(); x.n = r32(); x.nC = r32(); x.nPats = r32(); x.nMerg = r32();
    x.classOf.resize(x.n);
    for (int v = 0; v < x.n; v++) x.classOf[v] = r32();
    x.mergSize.resize(x.nMerg); x.mergV.resize(x.nMerg);
    for (int m = 0; m < x.nMerg; m++) {
        int sz = r32(); x.mergSize[m] = sz; x.mergV[m].resize(sz);
        for (int i = 0; i < sz; i++) x.mergV[m][i] = r32();
    }
    x.compOff.assign(x.nPats + 1, 0); x.cP.resize(x.nPats); x.k0.resize(x.nPats);
    for (int p = 0; p < x.nPats; p++) {
        int len = r8();                                  // writer emits the comp length as uint8
        x.compOff[p + 1] = x.compOff[p] + len;
        for (int i = 0; i < len; i++) { int32_t c = r32(); uint8_t b = r8(); x.compFlat.push_back({c, b}); }
        x.cP[p] = r32(); x.k0[p] = rd();
    }
    int nCells = x.smax - x.s0;
    x.resid.resize(nCells < 0 ? 0 : nCells);
    for (auto &rl : x.resid) {
        long long m = r64(); rl.resize(m);
        for (auto &pc : rl) { pc.first = r32(); pc.second = rd(); }
    }
    fseek(f, 0, SEEK_END); x.bytes = ftell(f);
    fclose(f);
    return true;
}

static double nCr(long long a, long long b) {           // double: only used for size arithmetic
    if (b < 0 || b > a) return 0.0;
    if (b > a - b) b = a - b;
    double r = 1.0;
    for (long long i = 0; i < b; i++) r = r * (double)(a - i) / (double)(i + 1);
    return r;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s size <index.nsi>\n       %s bench <index.nsi> <queryfile> [cap]\n", argv[0], argv[0]);
        return 1;
    }
    string mode = argv[1];
    NSI1 x;
    if (!load(argv[2], x)) { fprintf(stderr, "cannot load NSI1 %s\n", argv[2]); return 1; }
    const int r = x.r, nCells = x.smax - x.s0 + 1;

    // class -> vertices (needed both to size and to expand)
    vector<vector<int32_t>> clsV(x.nC);
    for (int v = 0; v < x.n; v++) if (x.classOf[v] >= 0 && x.classOf[v] < x.nC) clsV[x.classOf[v]].push_back(v);

    if (mode == "sample") {                            // emit random REAL r-cliques from the pattern table
        const long long want = argc >= 4 ? atoll(argv[3]) : 100000;
        unsigned seed = 12345;
        auto rnd = [&] { seed = seed * 1664525u + 1013904223u; return seed >> 1; };
        long long made = 0, guard = 0;
        vector<int32_t> out(r);
        while (made < want && guard++ < want * 50) {
            int p = (int)(rnd() % (unsigned)x.nPats);
            int filled = 0; bool ok = true;
            for (int i = x.compOff[p]; i < x.compOff[p + 1] && ok; i++) {
                const auto &V = clsV[x.compFlat[i].first];
                int b = x.compFlat[i].second;
                if ((int)V.size() < b) { ok = false; break; }
                for (int t = 0; t < b; t++) {           // sample b distinct members
                    int pick; bool dup;
                    int tries = 0;
                    do { pick = V[rnd() % V.size()]; dup = false;
                         for (int u = 0; u < filled; u++) if (out[u] == pick) { dup = true; break; }
                    } while (dup && ++tries < 50);
                    if (dup) { ok = false; break; }
                    out[filled++] = pick;
                }
            }
            if (!ok || filled != r) continue;
            sort(out.begin(), out.end());
            for (int i = 0; i < r; i++) printf("%d%c", out[i], i + 1 == r ? '\n' : ' ');
            made++;
        }
        fprintf(stderr, "sampled %lld r-cliques\n", made);
        return 0;
    }

    // exact #r-cliques: patterns carry mult = prod C(|class|, b); mergeable regions carry C(|M|, r)
    double patRC = 0.0;
    for (int p = 0; p < x.nPats; p++) {
        double m = 1.0;
        for (int i = x.compOff[p]; i < x.compOff[p + 1]; i++)
            m *= nCr((long long)clsV[x.compFlat[i].first].size(), x.compFlat[i].second);
        patRC += m;
    }
    double mergRC = 0.0;
    for (int m = 0; m < x.nMerg; m++) mergRC += nCr(x.mergSize[m], r);
    const double totRC = patRC + mergRC;

    // archive layout: r int32 key + one int32 core per cell, per r-clique
    const double archBytes = totRC * (double)(4 * r + 4 * nCells);
    printf("index    : %s  r=%d  cells=%d (s=%d..%d)  patterns=%d\n", argv[2], r, nCells, x.s0, x.smax, x.nPats);
    printf("index    : %.2f MB  (%.4f B / r-clique)\n", x.bytes / 1048576.0, x.bytes / (totRC > 0 ? totRC : 1));
    printf("r-cliques: %.0f  (pattern-side %.0f + mergeable-side %.0f)\n", totRC, patRC, mergRC);
    printf("ARCHIVE  : %.2f MB  (%.1f B / r-clique)  -> INDEX IS %.1fx SMALLER\n",
           archBytes / 1048576.0, 4.0 * r + 4.0 * nCells, x.bytes > 0 ? archBytes / x.bytes : 0.0);
    if (mode == "size") return 0;
    if (argc < 4) { fprintf(stderr, "bench needs a queryfile\n"); return 1; }

    const double cap = argc >= 5 ? atof(argv[4]) : 3e8;
    if (totRC > cap) {
        printf("ARCHIVE  : **CANNOT BE MATERIALIZED** (%.3g r-cliques > cap %.3g; would need %.1f GB)\n",
               totRC, cap, archBytes / 1073741824.0);
        printf("           This is the regime the index exists for: no sorted table can be built,\n"
               "           so the only alternative is recomputing the cell per query.\n");
        return 0;
    }

    // ---- expand patterns into the archive ----
    auto T0 = Clock::now();
    const size_t N = (size_t)totRC;
    vector<int32_t> keys; keys.reserve(N * r);
    vector<int32_t> vals; vals.reserve(N);            // one cell's core (s0) is enough to probe
    vector<int32_t> cur(r), keyBuf(r);
    // recursive product over the pattern's classes.  The emitted key must be sorted, but sorting
    // `cur` IN PLACE permutes the blocks the OUTER recursion levels still own, so every later tuple
    // comes out corrupted (it cost a 3.6% phantom miss rate before it was found).  Sort into keyBuf.
    for (int p = 0; p < x.nPats; p++) {
        const int b0 = x.compOff[p], b1 = x.compOff[p + 1];
        const int32_t core = (int32_t)x.k0[p];
        // choose comb of size b from class c, for each class in the comp
        struct Frame { int ci; vector<int32_t> pick; };
        // iterative product via recursion lambda
        auto rec = [&](auto &&self, int ci, int filled) -> void {
            if (ci == b1) {
                copy(cur.begin(), cur.begin() + filled, keyBuf.begin());
                sort(keyBuf.begin(), keyBuf.begin() + filled);
                for (int i = 0; i < filled; i++) keys.push_back(keyBuf[i]);
                vals.push_back(core);
                return;
            }
            const auto &V = clsV[x.compFlat[ci].first];
            const int b = x.compFlat[ci].second, m = (int)V.size();
            vector<int> idx(b);
            for (int i = 0; i < b; i++) idx[i] = i;
            if (b > m) return;
            while (true) {
                for (int i = 0; i < b; i++) cur[filled + i] = V[idx[i]];
                self(self, ci + 1, filled + b);
                int i = b - 1;
                while (i >= 0 && idx[i] == m - b + i) i--;
                if (i < 0) break;
                idx[i]++;
                for (int j = i + 1; j < b; j++) idx[j] = idx[j - 1] + 1;
            }
        };
        rec(rec, b0, 0);
    }
    double tExpand = chrono::duration<double>(Clock::now() - T0).count();

    // sort (key, val) together
    T0 = Clock::now();
    const size_t M = vals.size();
    vector<uint32_t> ord(M);
    for (size_t i = 0; i < M; i++) ord[i] = (uint32_t)i;
    auto less_at = [&](uint32_t a, uint32_t b) {
        const int32_t *ka = &keys[(size_t)a * r], *kb = &keys[(size_t)b * r];
        for (int i = 0; i < r; i++) if (ka[i] != kb[i]) return ka[i] < kb[i];
        return false;
    };
    sort(ord.begin(), ord.end(), less_at);
    vector<int32_t> sk((size_t)M * r), sv(M);
    for (size_t i = 0; i < M; i++) {
        memcpy(&sk[i * r], &keys[(size_t)ord[i] * r], sizeof(int32_t) * r);
        sv[i] = vals[ord[i]];
    }
    double tSort = chrono::duration<double>(Clock::now() - T0).count();
    const double realBytes = (double)M * (4.0 * r + 4.0 * nCells);
    printf("ARCHIVE  : materialized %zu entries in %.1fs (expand %.1fs + sort %.1fs), %.2f MB\n",
           M, tExpand + tSort, tExpand, tSort, realBytes / 1048576.0);

    // ---- probe with the same workload nsi_query bench uses ----
    FILE *qf = fopen(argv[3], "r");
    if (!qf) { fprintf(stderr, "cannot open %s\n", argv[3]); return 1; }
    vector<int32_t> flat; int v;
    while (fscanf(qf, "%d", &v) == 1) flat.push_back(v);
    fclose(qf);
    const long long nq = (long long)flat.size() / r;
    if (!nq) { fprintf(stderr, "no queries\n"); return 1; }
    vector<int32_t> probe(r);
    long long hits = 0; double sink = 0;
    T0 = Clock::now();
    for (long long i = 0; i < nq; i++) {
        memcpy(probe.data(), &flat[i * r], sizeof(int32_t) * r);
        sort(probe.begin(), probe.end());
        size_t lo = 0, hi = M;
        while (lo < hi) {
            size_t mid = (lo + hi) / 2;
            const int32_t *km = &sk[mid * r];
            int c = 0;
            for (int j = 0; j < r; j++) if (km[j] != probe[j]) { c = km[j] < probe[j] ? -1 : 1; break; }
            if (c < 0) lo = mid + 1; else hi = mid;
        }
        if (lo < M) {
            const int32_t *km = &sk[lo * r];
            bool eq = true;
            for (int j = 0; j < r; j++) if (km[j] != probe[j]) { eq = false; break; }
            if (eq) { hits++; sink += sv[lo]; }
        }
    }
    double tp = chrono::duration<double>(Clock::now() - T0).count();
    printf("ARCHIVE  : %lld probes, %.0f ns/probe, hits=%lld (sink=%.0f)\n", nq, 1e9 * tp / nq, hits, sink);
    return 0;
}
