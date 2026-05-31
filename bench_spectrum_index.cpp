// bench_spectrum_index.cpp
// ---------------------------------------------------------------------------
// Materialized Nucleus Spectrum Index for r=1 vertex s-clique-core values.
//
// Answers: build time, index size, query latency — the three numbers a
// cross-s index paper needs.  Compares against the no-index baseline (rebuild
// the CPI + peel from scratch for each queried s).
//
// Design
// ------
// BUILD-ONCE universal CPI: SDCT_Augmented_NoTree_Universal emits every
// maximal-clique leaf once (H=keep, P=pivot).  One SDCT traversal serves ALL
// s: leaf L contributes C(|P_L|, s-|H_L|) s-cliques, so kappa_s is obtained by
//   (1) per-leaf weights  wKeep=C(p,np), wPivot=C(p-1,np-1), np=s-|H_L|
//   (2) race-free per-vertex scatter over the vtx->leaf CSR (OpenMP)
//   (3) the validated SPIN* sparse-bucket peel (ST_V3_Peel)
// This replaces omega independent SDCT builds with ONE build + omega cheap
// recompute+peel passes.
//
// Two query backends are materialized:
//   * DENSE-by-s : per-level CSR of (vertex, kappa) for active vertices.
//                  query(s) = O(n_active_s) slice read (optimal, output-bound).
//   * ANCHOR+DELTA: per vertex (s*, kappa_{s*}) + sparse KK-shadow deltas.
//                  query(s) = reconstruct downward; smaller, slower.
//
// Correctness gate: index kappa_s is compared bit-exactly against the per-s
// brute-force path (ST_V3_Build+Peel), which bench_spectrum_sat already
// validated identical to degeneracy_cliques.
//
// Usage:
//   bench_spectrum_index <graph.edges> [--sort degen] [--repeats N]
//                        [--queries s1,s2,...] [--no-baseline] [--selftest]
// ---------------------------------------------------------------------------

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "src/graph/Graph.h"
#include "src/misc.h"
#include "src/Global/Global.h"
#include "src/SDCT_Augmented.h"
#include "src/NucleusDecomposition/NCliqueCoreDecomposition.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

using Clock = std::chrono::high_resolution_clock;
static double ms_since(Clock::time_point t0) {
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(
               Clock::now() - t0).count();
}

static constexpr double EXACT_LIMIT = 9.0e15;  // ~2^53

// ---- Kruskal-Katona lower shadow (for the anchor+delta backend) ----
static inline long double binomLD(long double n, int t) {
    if (t < 0) return 0.0L;
    if (t == 0) return 1.0L;
    if (n < (long double)t) return 0.0L;
    if (n <= 1000.0L && t <= 400) {
        long double fn = std::floor(n);
        if (fn == n) return (long double)nCr[(int)fn][t];
    }
    long double r = 1.0L;
    for (int i = 0; i < t; ++i) r = r * (n - (long double)i) / (long double)(i + 1);
    return r;
}
static inline long double largest_a(long double rem, int t) {
    long double hi = t;
    while (binomLD(hi + 1, t) <= rem) hi = hi * 2 + 2;
    long double L = t, R = hi, ans = t;
    while (L <= R) {
        long double mid = std::floor((L + R) / 2);
        if (binomLD(mid, t) <= rem) { ans = mid; L = mid + 1; } else R = mid - 1;
    }
    return ans;
}
static long double kk_shadow(long double k, int s) {
    if (k <= 0.0L) return 0.0L;
    long double shadow = 0.0L, rem = k; int j = s;
    while (rem > 0.0L && j >= 1) {
        long double a = largest_a(rem, j);
        shadow += binomLD(a, j - 1);
        rem    -= binomLD(a, j);
        j--;
    }
    return shadow;
}

// ===========================================================================
//  Universal CPI: build the dual CSR + per-leaf (keepC, pivotC) ONCE.
//  Mirrors ST_V3_Build's CSR construction but (a) drives the UNIVERSAL SDCT
//  and (b) stores keepC per leaf instead of a fixed-s countingV.
// ===========================================================================
struct UniversalCPI {
    ST_V2_Data d;                 // dual CSR; leafPivotCount filled, countingV per-s
    std::vector<int> leafKeepC;   // |H| per leaf (invariant across s)
    int omega = 0;
    size_t sigma = 0;             // total vtx-leaf incidences
};

static void buildUniversalCPI(Graph &g, UniversalCPI &U) {
    ST_V2_Data &d = U.d;
    d.numVertices = g.getGraphNodeSize();

    struct COOEntry { daf::Size vertex; daf::Size leafId; uint8_t isPivot; };
    std::vector<COOEntry> coo;
    coo.reserve(1 << 22);

    auto &leafKeepC = U.leafKeepC;

    d.numLeaves = SDCT_Augmented_NoTree_Universal(g, /*min_k=*/1,
        [&](daf::Size leafId, const daf::StaticVector<int> &keepV,
            const daf::StaticVector<int> &dropV)
        {
            int keepC = (int)keepV.size();
            int pivotC = (int)dropV.size();
            if (leafId >= d.leafPivotCount.size()) {
                size_t ns = std::max<size_t>(leafId + 1, d.leafPivotCount.size() * 2);
                d.leafPivotCount.resize(ns, 0);
                d.leafNeedPivot.resize(ns, 0);
                leafKeepC.resize(ns, 0);
            }
            d.leafPivotCount[leafId] = pivotC;
            leafKeepC[leafId] = keepC;
            for (int i = 0; i < keepC; ++i)  coo.push_back({(daf::Size)keepV[i], leafId, 0});
            for (int i = 0; i < pivotC; ++i) coo.push_back({(daf::Size)dropV[i], leafId, 1});
        });

    d.leafPivotCount.resize(d.numLeaves);
    d.leafNeedPivot.resize(d.numLeaves);
    leafKeepC.resize(d.numLeaves);

    // ---- vtx->leaf CSR (packed ids + 1 bit/incidence) ----
    d.vtxLeafOff.assign(d.numVertices + 2, 0);
    for (auto &e : coo) if (e.vertex < d.numVertices) d.vtxLeafOff[e.vertex + 1]++;
    for (daf::Size i = 1; i <= d.numVertices; ++i) d.vtxLeafOff[i] += d.vtxLeafOff[i - 1];
    const size_t total = d.vtxLeafOff[d.numVertices];
    U.sigma = total;
    d.vtxLeafIds.resize(total);
    d.vtxLeafIsPivot.assign((total + 63) >> 6, 0);
    d.leafVtxOff.assign(d.numLeaves + 1, 0);
    {
        std::vector<size_t> pos(d.numVertices, 0);
        for (auto &e : coo) {
            daf::Size v = e.vertex;
            if (v < d.numVertices) {
                size_t p = d.vtxLeafOff[v] + pos[v]++;
                d.vtxLeafIds[p] = e.leafId;
                if (e.isPivot) STV3_setBit(d.vtxLeafIsPivot, p);
                if (e.leafId < d.numLeaves) d.leafVtxOff[e.leafId + 1]++;
            }
        }
    }
    std::vector<COOEntry>().swap(coo);

    // ---- leaf->vtx CSR (transpose) ----
    for (size_t i = 1; i <= d.numLeaves; ++i) d.leafVtxOff[i] += d.leafVtxOff[i - 1];
    const size_t totalL = d.leafVtxOff[d.numLeaves];
    d.leafVtxIds.resize(totalL);
    d.leafVtxIsPivot.assign((totalL + 63) >> 6, 0);
    {
        std::vector<size_t> pos(d.numLeaves, 0);
        for (daf::Size v = 0; v < d.numVertices; ++v) {
            const size_t b = d.vtxLeafOff[v], e = d.vtxLeafOff[v + 1];
            for (size_t k = b; k < e; ++k) {
                daf::Size L = d.vtxLeafIds[k];
                if (L < d.numLeaves) {
                    size_t p = d.leafVtxOff[L] + pos[L]++;
                    d.leafVtxIds[p] = v;
                    if (STV3_getBit(d.vtxLeafIsPivot, k)) STV3_setBit(d.leafVtxIsPivot, p);
                }
            }
        }
    }

    // omega = max maximal-clique size = max_L (keepC + pivotC)
    int omega = 0;
    for (size_t L = 0; L < d.numLeaves; ++L)
        omega = std::max(omega, leafKeepC[L] + d.leafPivotCount[L]);
    U.omega = omega;
}

// ===========================================================================
//  Recompute countingV + leafNeedPivot for a target s from the universal CPI.
//  Race-free per-vertex scatter over the vtx->leaf CSR (OpenMP-parallel).
//  wKeep/wPivot are scratch arrays of size numLeaves (caller-owned, reused).
// ===========================================================================
static void recomputeCountingV(UniversalCPI &U, int s,
                               std::vector<double> &wKeep,
                               std::vector<double> &wPivot) {
    ST_V2_Data &d = U.d;
    const size_t L = d.numLeaves;
    for (size_t i = 0; i < L; ++i) {
        int p  = d.leafPivotCount[i];
        int np = s - U.leafKeepC[i];
        d.leafNeedPivot[i] = np;
        if (np >= 0 && np <= p) {
            wKeep[i]  = nCr[p][np];
            wPivot[i] = (np >= 1) ? nCr[p - 1][np - 1] : 0.0;
        } else { wKeep[i] = 0.0; wPivot[i] = 0.0; }
    }
    if (!d.countingV) d.countingV = new double[d.numVertices];
    const daf::Size n = d.numVertices;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (daf::Size v = 0; v < n; ++v) {
        const size_t b = d.vtxLeafOff[v], e = d.vtxLeafOff[v + 1];
        double sum = 0.0;
        for (size_t i = b; i < e; ++i) {
            daf::Size leaf = d.vtxLeafIds[i];
            sum += STV3_getBit(d.vtxLeafIsPivot, i) ? wPivot[leaf] : wKeep[leaf];
        }
        d.countingV[v] = sum;
    }
}

// kappa_s(v) from a peel result (coreV[v] = -1 if v never had support).
static inline double kappa_of(const double *coreV, daf::Size v) {
    return coreV[v] < 0.0 ? 0.0 : coreV[v];
}

// Resident bytes of the universal CPI ("tree" substrate): the dual CSR + the
// per-leaf metadata that must stay in memory to answer queries by recompute+
// peel.  Excludes the transient countingV (allocated/freed per peel).
template <class T> static size_t capb(const std::vector<T>& v){ return v.capacity()*sizeof(T); }
static size_t cpi_bytes(const UniversalCPI &U) {
    const ST_V2_Data &d = U.d;
    return capb(d.vtxLeafOff) + capb(d.leafVtxOff)
         + capb(d.vtxLeafIds) + capb(d.vtxLeafIsPivot)
         + capb(d.leafVtxIds) + capb(d.leafVtxIsPivot)
         + capb(d.leafPivotCount) + capb(d.leafNeedPivot)
         + capb(U.leafKeepC);
}

// ===========================================================================
//  Materialized index
// ===========================================================================
struct DenseIndex {                 // per-level CSR of active (vertex, kappa)
    std::vector<size_t> off;        // size omega+2; level s in [off[s],off[s+1])
    std::vector<uint32_t> vid;
    std::vector<double> kappa;
    size_t bytes() const {
        return off.capacity() * sizeof(size_t)
             + vid.capacity() * sizeof(uint32_t)
             + kappa.capacity() * sizeof(double);
    }
};
struct AnchorDeltaIndex {           // per-vertex (s*, kappa*) + sparse deltas
    std::vector<uint16_t> sstar;    // [n]
    std::vector<double>   kstar;    // [n]
    std::vector<size_t>   doff;     // [n+1] CSR offsets into deltaS/deltaV
    std::vector<uint16_t> deltaS;   // level of each stored delta
    std::vector<double>   deltaV;   // delta value (descending-s order per vtx)
    // query accelerator: vertices sorted by s* descending; the first prefix[s]
    // entries of `order` are exactly the vertices active at level s.
    std::vector<uint32_t> order;
    std::vector<uint32_t> prefix;   // [omega+2]; prefix[s] = #{v : s*(v) >= s}
    size_t bytes() const {          // counts the persistent index only (not order/prefix accel)
        return sstar.capacity() * sizeof(uint16_t)
             + kstar.capacity() * sizeof(double)
             + doff.capacity() * sizeof(size_t)
             + deltaS.capacity() * sizeof(uint16_t)
             + deltaV.capacity() * sizeof(double);
    }
};

// ===========================================================================
int main(int argc, char *argv[]) {
    populate_nCr();

    if (argc >= 2 && std::strcmp(argv[1], "--selftest") == 0) {
        struct { long double k; int s; long double want; } C[] = {
            {1,4,4},{4,3,6},{6,2,4},{1,3,3},{3,2,3},{10,2,5} };
        int f = 0;
        for (auto &c : C) if (std::fabsl(kk_shadow(c.k, c.s) - c.want) > 1e-9L) f++;
        std::printf("SELFTEST %s\n", f ? "FAILED" : "PASSED");
        return f ? 1 : 0;
    }
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <graph.edges> [--sort degen] [--repeats N]"
                     " [--queries s1,s2,..] [--no-baseline]\n";
        return 1;
    }
    const char *gpath = argv[1];
    std::string sortOpt = "degen";
    int repeats = 0;                     // 0 => auto
    bool do_baseline = true;
    std::vector<int> queries;
    for (int i = 2; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--sort") && i + 1 < argc) sortOpt = argv[++i];
        else if (!std::strcmp(argv[i], "--repeats") && i + 1 < argc) repeats = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--no-baseline")) do_baseline = false;
        else if (!std::strcmp(argv[i], "--queries") && i + 1 < argc) {
            char *tok = std::strtok(argv[++i], ",");
            while (tok) { queries.push_back(std::atoi(tok)); tok = std::strtok(nullptr, ","); }
        }
    }

    Graph g(gpath);
    g.printGraphInfo();
    daf::vListMap.resize(g.getGraphNodeSize() + 1);
    std::memset(daf::vListMap.data(), -1, (size_t)g.getGraphNodeSize() * sizeof(daf::Size));
    if      (sortOpt == "degen")  g.sortByDegeneracyOrder(false);
    else if (sortOpt == "degenR") g.sortByDegeneracyOrder(true);
    else if (sortOpt != "default") { std::cerr << "bad --sort\n"; return 1; }

    const int  n = (int)g.getGraphNodeSize();
    const long m = (long)g.getGraphEdgeSize();
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif

    // ===== PHASE 1: build universal CPI once =====
    UniversalCPI U;
    auto t0 = Clock::now();
    buildUniversalCPI(g, U);
    double build_universal_ms = ms_since(t0);
    const int omega = U.omega;
    std::printf("IDX_BUILD graph=%s n=%d m=%ld omega=%d leaves=%zu sigma=%zu "
                "build_universal_ms=%.1f threads=%d\n",
                gpath, n, m, omega, U.d.numLeaves, U.sigma, build_universal_ms, nthreads);
    std::fflush(stdout);

    // ===== PHASE 2: compute full spectrum from U, materialize both indexes =====
    DenseIndex dense;
    dense.off.assign(omega + 2, 0);
    AnchorDeltaIndex ad;
    ad.sstar.assign(n, 0);
    ad.kstar.assign(n, 0.0);
    std::vector<double> wKeep(U.d.numLeaves), wPivot(U.d.numLeaves);
    std::vector<double> prevK(n, 0.0), curK(n, 0.0);
    std::vector<std::vector<std::pair<uint16_t,double>>> perVtxDeltas(n);

    double recompute_ms = 0, peel_ms = 0;
    auto t_spec = Clock::now();
    for (int s = 2; s <= omega; ++s) {
        auto tr = Clock::now();
        recomputeCountingV(U, s, wKeep, wPivot);
        recompute_ms += ms_since(tr);

        auto tp = Clock::now();
        double *coreV = NCliqueVertexCoreDecomposition_ST_V3_Peel(U.d, (daf::CliqueSize)s);
        peel_ms += ms_since(tp);

        dense.off[s] = dense.vid.size();
        for (int v = 0; v < n; ++v) {
            double kv = kappa_of(coreV, v);
            curK[v] = kv;
            if (kv > 0.0) {
                dense.vid.push_back((uint32_t)v);
                dense.kappa.push_back(kv);
                ad.sstar[v] = (uint16_t)s;
                ad.kstar[v] = kv;
            }
        }
        // delta_{s-1}(v) = kappa_{s-1} - shadow_{s-1}(kappa_s); store nonzero
        if (s >= 3) {
            for (int v = 0; v < n; ++v) {
                if (curK[v] <= 0.0) continue;
                double ks = curK[v], ksm1 = prevK[v];
                if (ksm1 <= 0.0) continue;
                if (ks >= EXACT_LIMIT || ksm1 >= EXACT_LIMIT) {
                    // clamped/fuzzy: store exact residual so reconstruction stays correct
                    perVtxDeltas[v].push_back({(uint16_t)(s - 1), ksm1 - (double)kk_shadow((long double)ks, s - 1)});
                    continue;
                }
                double d = ksm1 - (double)kk_shadow((long double)ks, s - 1);
                if (std::fabs(d) > 0.5)
                    perVtxDeltas[v].push_back({(uint16_t)(s - 1), d});
            }
        }
        delete[] coreV;
        std::swap(prevK, curK);
    }
    dense.off[omega + 1] = dense.vid.size();
    double spectrum_ms = ms_since(t_spec);

    // assemble anchor+delta CSR (deltas stored descending-s for reconstruction)
    ad.doff.assign(n + 1, 0);
    for (int v = 0; v < n; ++v) ad.doff[v + 1] = ad.doff[v] + perVtxDeltas[v].size();
    ad.deltaS.resize(ad.doff[n]);
    ad.deltaV.resize(ad.doff[n]);
    for (int v = 0; v < n; ++v) {
        auto &dv = perVtxDeltas[v];
        std::sort(dv.begin(), dv.end(),
                  [](auto &a, auto &b){ return a.first > b.first; });  // descending s
        size_t base = ad.doff[v];
        for (size_t i = 0; i < dv.size(); ++i) { ad.deltaS[base+i]=dv[i].first; ad.deltaV[base+i]=dv[i].second; }
    }
    std::vector<std::vector<std::pair<uint16_t,double>>>().swap(perVtxDeltas);

    // query accelerator: bucket vertices by s* (descending) so query(s) iterates
    // only the active prefix instead of scanning all n.
    {
        std::vector<uint32_t> cnt(omega + 2, 0);
        for (int v = 0; v < n; ++v) if (ad.sstar[v] >= 2) cnt[ad.sstar[v]]++;
        ad.prefix.assign(omega + 2, 0);
        uint32_t acc = 0;
        for (int s = omega; s >= 2; --s) { acc += cnt[s]; ad.prefix[s] = acc; }
        ad.prefix[1] = ad.prefix[0] = acc;
        std::vector<uint32_t> start(omega + 2, 0), cur;
        uint32_t off = 0;
        for (int t = omega; t >= 2; --t) { start[t] = off; off += cnt[t]; }
        ad.order.resize(off);
        cur = start;
        for (int v = 0; v < n; ++v) if (ad.sstar[v] >= 2) ad.order[cur[ad.sstar[v]]++] = (uint32_t)v;
    }

    size_t n_active = 0, traj = 0; for (int v=0;v<n;++v) if (ad.sstar[v]>=2){n_active++; traj+=ad.sstar[v]-1;}
    std::printf("IDX_SIZE graph=%s n_active=%zu traj_len_sum=%zu dense_entries=%zu "
                "dense_bytes=%zu anchordelta_deltas=%zu anchordelta_bytes=%zu "
                "spectrum_ms=%.1f recompute_ms=%.1f peel_ms=%.1f\n",
                gpath, n_active, traj, dense.vid.size(), dense.bytes(),
                ad.deltaS.size(), ad.bytes(), spectrum_ms, recompute_ms, peel_ms);
    std::fflush(stdout);

    // Memory comparison: CPI/tree substrate vs the two materialized indexes.
    size_t cpib = cpi_bytes(U), dnb = dense.bytes(), adb = ad.bytes();
    std::printf("IDX_MEM graph=%s cpi_tree_bytes=%zu cpi_tree_MB=%.2f "
                "dense_bytes=%zu dense_MB=%.2f anchordelta_bytes=%zu anchordelta_MB=%.2f "
                "dense_over_cpi=%.2f anchordelta_over_cpi=%.2f bytes_per_vertex_cpi=%.1f\n",
                gpath, cpib, cpib/1e6, dnb, dnb/1e6, adb, adb/1e6,
                cpib>0 ? (double)dnb/cpib : 0.0, cpib>0 ? (double)adb/cpib : 0.0,
                (double)cpib/std::max(1,n));
    std::fflush(stdout);

    // ===== query set =====
    if (queries.empty()) {           // default: a spread across [2, omega]
        for (int s : {2,3,4,5}) if (s <= omega) queries.push_back(s);
        for (int frac : {25,50,75,90}) { int s = std::max(2, omega*frac/100); if (s<=omega) queries.push_back(s); }
        if (omega >= 2) queries.push_back(omega);
        std::sort(queries.begin(), queries.end());
        queries.erase(std::unique(queries.begin(), queries.end()), queries.end());
    }
    (void)n_active;

    // anchor+delta reconstruction of kappa_s(v); returns 0 if v inactive at s.
    auto recon = [&](int v, int s) -> double {
        if (ad.sstar[v] < s) return 0.0;
        double k = ad.kstar[v];
        int sv = ad.sstar[v];
        if (sv > s) {
            size_t e = ad.doff[v+1], di = ad.doff[v];
            for (int j = sv - 1; j >= s; --j) {
                double dval = 0.0;
                if (di < e && ad.deltaS[di] == (uint16_t)j) { dval = ad.deltaV[di]; ++di; }
                k = (double)kk_shadow((long double)k, j) + dval;
            }
        }
        return k;
    };
    // query backends return a checksum (forces the work; equal across backends)
    auto q_dense = [&](int s) -> double {
        double sum = 0; size_t b = dense.off[s], e = dense.off[s+1];
        for (size_t i = b; i < e; ++i) sum += (double)dense.vid[i] * dense.kappa[i];
        return sum;
    };
    auto q_anchordelta = [&](int s) -> double {
        double sum = 0;
        const uint32_t lim = ad.prefix[s];          // only vertices active at s
        for (uint32_t i = 0; i < lim; ++i) {
            uint32_t v = ad.order[i];
            double k = recon((int)v, s);
            if (k > 0.0) sum += (double)v * k;
        }
        return sum;
    };
    // "build tree ONCE, reuse" baseline: from the already-built universal CPI,
    // recompute support + peel for this s (NO tree rebuild).  This is the fair
    // shared-tree alternative; only the full per-s rebuild (q below) rebuilds.
    auto q_shared = [&](int s) -> double {
        recomputeCountingV(U, s, wKeep, wPivot);     // allocates U.d.countingV
        double *cv = NCliqueVertexCoreDecomposition_ST_V3_Peel(U.d, (daf::CliqueSize)s);
        double sum = 0;
        for (int v = 0; v < n; ++v) { double k = kappa_of(cv, v); if (k > 0.0) sum += (double)v * k; }
        delete[] cv;
        return sum;
    };

    // ===== correctness: index kappa_s == per-s brute force (ST_V3_Build+Peel) =====
    // Compared PER VERTEX (not just a checksum) against the path that
    // bench_spectrum_sat already validated bit-identical to degeneracy_cliques.
    int mism = 0;
    std::vector<double> idxK(n, 0.0);
    for (int s : queries) {
        ST_V2_Data bd = NCliqueVertexCoreDecomposition_ST_V3_Build(g, (daf::CliqueSize)s);
        double *bcore = (bd.numLeaves == 0) ? nullptr
                        : NCliqueVertexCoreDecomposition_ST_V3_Peel(bd, (daf::CliqueSize)s);
        if (bd.numLeaves == 0 && bd.countingV) { delete[] bd.countingV; bd.countingV = nullptr; }

        std::fill(idxK.begin(), idxK.end(), 0.0);
        for (size_t i = dense.off[s]; i < dense.off[s+1]; ++i) idxK[dense.vid[i]] = dense.kappa[i];

        int reported = 0;
        for (int v = 0; v < n; ++v) {
            double bk = bcore ? kappa_of(bcore, v) : 0.0;
            double tol = std::max(0.5, 1e-9 * std::fabs(bk));
            double dk = idxK[v], ak = recon(v, s);
            if (std::fabs(dk - bk) > tol || std::fabs(ak - bk) > tol) {
                mism++;
                if (reported < 5) {
                    std::printf("  [MISMATCH] s=%d v=%d dense=%.0f anchordelta=%.0f brute=%.0f\n",
                                s, v, dk, ak, bk);
                    reported++;
                }
            }
        }
        if (bcore) delete[] bcore;
    }
    std::printf("IDX_VERIFY graph=%s queries=%zu mismatches=%d %s\n",
                gpath, queries.size(), mism, mism ? "FAIL" : "OK");
    std::fflush(stdout);

    // ===== query latency: index vs baseline (rebuild) =====
    // Time-budgeted per backend: repeat until ~budget_ms elapse or maxrep hit,
    // then report per-op microseconds.  Auto-adapts so a slow backend (e.g.
    // anchor+delta reconstruction at omega=239 ~300ms/op) runs ~once, while a
    // ns-scale dense read runs millions of times for a stable number.
    const double budget_ms = 40.0;
    const int    maxrep    = (repeats > 0) ? repeats : 2000000;
    volatile double sink = 0;
    auto bench = [&](auto &&fn) -> double {
        // warmup
        sink += fn();
        int r = 0; auto t = Clock::now(); double el = 0;
        do { sink += fn(); ++r; el = ms_since(t); } while (el < budget_ms && r < maxrep);
        return el * 1000.0 / r;   // microseconds per op
    };
    for (int s : queries) {
        double dense_us = bench([&]{ return q_dense(s); });
        double ad_us    = bench([&]{ return q_anchordelta(s); });
        // shared-tree baseline (build universal CPI once, reuse): recompute+peel.
        // Few-shot timed (ms-scale): run until ~budget or 5 reps.
        double shared_ms = -1;
        {
            sink += q_shared(s);                      // warmup
            int r = 0; auto t = Clock::now(); double el = 0;
            do { sink += q_shared(s); ++r; el = ms_since(t); } while (el < budget_ms && r < 5);
            shared_ms = el / r;
        }
        // full per-s rebuild baseline (no shared tree): build CPI for THIS s + peel.
        double rebuild_ms = -1;
        if (do_baseline) {
            auto tb = Clock::now();
            ST_V2_Data bd = NCliqueVertexCoreDecomposition_ST_V3_Build(g, (daf::CliqueSize)s);
            if (bd.numLeaves == 0) { if (bd.countingV){delete[] bd.countingV; bd.countingV=nullptr;} }
            else { double *bc = NCliqueVertexCoreDecomposition_ST_V3_Peel(bd, (daf::CliqueSize)s); delete[] bc; }
            rebuild_ms = ms_since(tb);
        }
        size_t na_s = dense.off[s+1] - dense.off[s];
        std::printf("IDX_QUERY graph=%s s=%d n_active_s=%zu dense_us=%.4f "
                    "anchordelta_us=%.3f shared_ms=%.3f rebuild_ms=%.3f "
                    "dense_vs_shared=%.1fx dense_vs_rebuild=%.1fx\n",
                    gpath, s, na_s, dense_us, ad_us, shared_ms, rebuild_ms,
                    shared_ms > 0 ? shared_ms * 1000.0 / dense_us : 0.0,
                    rebuild_ms > 0 ? rebuild_ms * 1000.0 / dense_us : 0.0);
        std::fflush(stdout);
    }
    (void)sink;

    if (U.d.countingV) { delete[] U.d.countingV; U.d.countingV = nullptr; }
    return 0;
}
