// bench_spectrum_sat.cpp
// ---------------------------------------------------------------------------
// Measure omega (max clique size) and the cross-s "saturation" /
// compressibility of the r=1 vertex s-clique-core SPECTRUM on real graphs.
//
// For each s = 2,3,..., we build the CPI fresh (NCliqueVertexCoreDecomposition_
// ST_V3_Build) and peel it (..._ST_V3_Peel) to obtain kappa_s(v) for every
// vertex v.  We then stream each vertex's trajectory kappa_2..kappa_{s*} and
// encode it as:
//      anchor = (s*(v), kappa_{s*}(v))
//      deltas = { s : delta_s }  where  delta_s = kappa_s - KKshadow(kappa_{s+1}, s)
// A vertex is "saturated" iff all its deltas are 0 (its whole trajectory is
// reconstructible from the anchor alone via the Kruskal-Katona lower shadow).
//
// IMPORTANT (correctness): the existing SDCT_Augmented_NoTree emits a leaf only
// when the maximal-clique size cSize >= max_k (= s), and the per-leaf
// (keepV,dropV) decomposition is specific to that s.  Hence the build is NOT
// reusable across s; we rebuild per s.  This costs omega builds but is the
// trusted code path (identical to the main binary's PIVOTER_RUN_ST_V3 route),
// so the measured kappa_s values are exactly the SPIN-star results.
//
// Reported invariants (also a correctness gate, reproduces the verified
// "law_viol=0 / nest_viol=0" on real graphs):
//   * law_viol  : #(v,s) with delta_s < 0  (Shifted Containment Law violation)
//   * nest_viol : #(v,s) with kappa_{s+1}>0 but kappa_s==0 (support-nesting viol)
//
// Usage:
//   bench_spectrum_sat <graph.edges> [--smax N] [--sort degen|degenR|default]
//                      [--selftest] [--dump-kappa <file>]
//
// Output (grep-friendly):
//   SAT_S       graph=.. s=.. leaves=.. n_active=.. kappa_max=.. nz_delta=..
//               build_ms=.. peel_ms=..
//   SAT_SUMMARY graph=.. n=.. m=.. maxdeg=.. omega=.. omega_trunc=.. n_active=..
//               saturated=.. sat_frac=.. nonzero_deltas=.. traj_len_sum=..
//               naive_store=.. anchordelta_store=.. comp_ratio=.. law_viol=..
//               nest_viol=.. fuzzy=.. total_ms=..
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

#include "src/graph/Graph.h"
#include "src/misc.h"
#include "src/Global/Global.h"
#include "src/NucleusDecomposition/NCliqueCoreDecomposition.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

// Doubles represent integers exactly only below 2^53.  Above this the
// shadow/delta comparison degrades to a relative tolerance; we count such
// (v,s) pairs as "fuzzy" and report them for honesty.
static constexpr double EXACT_LIMIT = 9.0e15;  // ~2^53

// ---------------------------------------------------------------------------
//  Binomial C(n,t) in long double, with an exact nCr-table fast path.
// ---------------------------------------------------------------------------
static inline long double binomLD(long double n, int t) {
    if (t < 0) return 0.0L;
    if (t == 0) return 1.0L;
    if (n < (long double)t) return 0.0L;
    if (n <= 1000.0L && t <= 400) {
        long double fn = std::floor(n);
        if (fn == n) return (long double)nCr[(int)fn][t];
    }
    long double r = 1.0L;
    for (int i = 0; i < t; ++i)
        r = r * (n - (long double)i) / (long double)(i + 1);
    return r;
}

// Largest a >= t with C(a,t) <= rem  (assumes rem >= 1, so a=t works).
static inline long double largest_a(long double rem, int t) {
    long double hi = t;
    while (binomLD(hi + 1, t) <= rem) hi = hi * 2 + 2;   // exponential expand
    long double L = t, R = hi, ans = t;
    while (L <= R) {
        long double mid = std::floor((L + R) / 2);
        if (binomLD(mid, t) <= rem) { ans = mid; L = mid + 1; }
        else                         R = mid - 1;
    }
    return ans;
}

// Kruskal-Katona lower shadow d_s(k): write k in the s-cascade
//   k = C(a_s,s) + C(a_{s-1},s-1) + ... ,  shadow = C(a_s,s-1) + C(a_{s-1},s-2) + ...
static long double kk_shadow(long double k, int s) {
    if (k <= 0.0L) return 0.0L;
    long double shadow = 0.0L, rem = k;
    int j = s;
    while (rem > 0.0L && j >= 1) {
        long double a = largest_a(rem, j);
        shadow += binomLD(a, j - 1);
        rem    -= binomLD(a, j);
        j--;
    }
    return shadow;
}

// ---------------------------------------------------------------------------
//  Self-test for kk_shadow against hand-verified karate-club values.
// ---------------------------------------------------------------------------
static int selftest() {
    struct C { long double k; int s; long double want; const char* note; };
    C cases[] = {
        {1, 4, 4, "d_4(1)=C(4,3)=4   (karate v0: kappa4=d_4(kappa5=1))"},
        {4, 3, 6, "d_3(4)=C(4,2)=6   (karate v0: kappa3=d_3(kappa4=4))"},
        {6, 2, 4, "d_2(6)=C(4,1)=4   (karate v0: kappa2=d_2(kappa3=6))"},
        {1, 3, 3, "d_3(1)=C(3,2)=3   (karate v8: kappa3=d_3(kappa4=1))"},
        {3, 2, 3, "d_2(3)=C(3,1)=3   (karate v8: kappa2=d_2(3)+1)"},
        {0, 5, 0, "d_s(0)=0"},
        {10, 2, 5, "d_2(10)=C(5,1)=5  (10=C(5,2))"},
    };
    int fails = 0;
    for (auto& c : cases) {
        long double got = kk_shadow(c.k, c.s);
        bool ok = std::fabsl(got - c.want) < 1e-9L;
        std::printf("  d_%d(%.0Lf) = %.0Lf  (want %.0Lf)  %s   %s\n",
                    c.s, c.k, got, c.want, ok ? "OK" : "FAIL", c.note);
        if (!ok) fails++;
    }
    std::printf("SELFTEST %s (%d failures)\n", fails ? "FAILED" : "PASSED", fails);
    return fails;
}

// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // ---- argument parsing ----
    if (argc >= 2 && std::strcmp(argv[1], "--selftest") == 0) {
        populate_nCr();
        return selftest() ? 1 : 0;
    }
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <graph.edges> [--smax N] [--sort degen|degenR|default]"
                     " [--selftest] [--dump-kappa <file>]\n";
        return 1;
    }
    const char* gpath = argv[1];
    int smax_cap = 400;                 // hard cap (nCr column bound); also safety
    std::string sortOpt = "degen";
    const char* dumpKappaPath = nullptr;
    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "--smax") == 0 && i + 1 < argc)
            smax_cap = std::min(400, std::atoi(argv[++i]));
        else if (std::strcmp(argv[i], "--sort") == 0 && i + 1 < argc)
            sortOpt = argv[++i];
        else if (std::strcmp(argv[i], "--selftest") == 0) { /* handled below */ }
        else if (std::strcmp(argv[i], "--dump-kappa") == 0 && i + 1 < argc)
            dumpKappaPath = argv[++i];
    }

    populate_nCr();
    if (selftest() != 0) {
        std::cerr << "FATAL: kk_shadow self-test failed; aborting.\n";
        return 2;
    }

    // ---- load + degeneracy-sort graph (graph is never mutated by Build) ----
    Graph g(gpath);
    g.printGraphInfo();
    daf::vListMap.resize(g.getGraphNodeSize() + 1);
    std::memset(daf::vListMap.data(), -1,
                (size_t)g.getGraphNodeSize() * sizeof(daf::Size));
    if      (sortOpt == "degen")  g.sortByDegeneracyOrder(false);
    else if (sortOpt == "degenR") g.sortByDegeneracyOrder(true);
    else if (sortOpt == "default") { /* keep input order */ }
    else { std::cerr << "Unknown --sort " << sortOpt << "\n"; return 1; }

    const int  n      = (int)g.getGraphNodeSize();
    const long m      = (long)g.getGraphEdgeSize();
    const long maxdeg = (long)g.getMaxDegree();

    // ---- per-vertex streaming accumulators (O(n) memory) ----
    std::vector<double>   prevK(n, 0.0);   // kappa_{s-1}
    std::vector<double>   curK(n, 0.0);    // kappa_s
    std::vector<uint16_t> sstar(n, 0);     // max s with kappa_s>0
    std::vector<double>   kstar(n, 0.0);   // kappa_{s*}
    std::vector<uint16_t> deltaCnt(n, 0);  // # nonzero EXACT deltas on v's trajectory
    std::vector<uint8_t>  fuzzyVtx(n, 0);  // trajectory touched kappa >= 2^53 (clamped)

    uint64_t lawViol = 0, nestViol = 0, fuzzy = 0, fuzzyTrans = 0;

    // optional full-kappa dump (small graphs / correctness inspection)
    FILE* dumpf = nullptr;
    if (dumpKappaPath) {
        dumpf = std::fopen(dumpKappaPath, "w");
        if (dumpf) std::fprintf(dumpf, "# s\tvertex\tkappa_s\n");
    }

    auto t_all0 = std::chrono::high_resolution_clock::now();
    int  omega = 0;
    bool omega_trunc = false;

    for (int s = 2; s <= smax_cap; ++s) {
        auto tb0 = std::chrono::high_resolution_clock::now();
        ST_V3_Data d = NCliqueVertexCoreDecomposition_ST_V3_Build(
                           g, (daf::CliqueSize)s);
        auto tb1 = std::chrono::high_resolution_clock::now();
        long long build_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(tb1 - tb0).count();

        if (d.numLeaves == 0) {                 // no maximal clique of size >= s
            if (d.countingV) { delete[] d.countingV; d.countingV = nullptr; }
            omega = s - 1;
            break;
        }

        auto tp0 = std::chrono::high_resolution_clock::now();
        double* coreV = NCliqueVertexCoreDecomposition_ST_V3_Peel(
                            d, (daf::CliqueSize)s);   // frees d.countingV
        auto tp1 = std::chrono::high_resolution_clock::now();
        long long peel_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(tp1 - tp0).count();

        // extract kappa_s into curK (kappa = max(coreV,0))
        long   n_active_s = 0;
        double kappa_max  = 0.0;
        for (int v = 0; v < n; ++v) {
            double kv = (coreV[v] < 0.0) ? 0.0 : coreV[v];
            curK[v] = kv;
            if (kv > 0.0) { n_active_s++; if (kv > kappa_max) kappa_max = kv; }
        }
        delete[] coreV;

        if (dumpf)
            for (int v = 0; v < n; ++v)
                if (curK[v] > 0.0)
                    std::fprintf(dumpf, "%d\t%d\t%.0f\n", s, v, curK[v]);

        // update anchor (ascending => last active s is s*)
        for (int v = 0; v < n; ++v)
            if (curK[v] > 0.0) { sstar[v] = (uint16_t)s; kstar[v] = curK[v]; }

        // delta_{s-1}(v) = kappa_{s-1}(v) - d_{s-1}(kappa_s(v)) for v active at s
        uint64_t nz_this_s = 0;
        if (s >= 3) {
            for (int v = 0; v < n; ++v) {
                if (curK[v] <= 0.0) continue;        // anchor boundary, no transition
                double ks   = curK[v];               // kappa_s
                double ksm1 = prevK[v];              // kappa_{s-1}
                if (ksm1 <= 0.0) { nestViol++; fuzzyVtx[v] = 1; continue; }
                if (ks >= EXACT_LIMIT || ksm1 >= EXACT_LIMIT) {
                    // beyond double's exact-integer range the peel clamps support
                    // (kBucketKeyClamp=1e18), so the shadow/delta comparison is
                    // numerically meaningless: exclude from the trustworthy stats.
                    fuzzy++; fuzzyTrans++; fuzzyVtx[v] = 1; nz_this_s++;
                    continue;
                }
                long double sh    = kk_shadow((long double)ks, s - 1);
                double      delta = (double)((long double)ksm1 - sh);
                if (delta < -0.5) lawViol++;          // Shifted Containment Law (exact)
                if (std::fabs(delta) > 0.5) {         // nonzero exact delta
                    if (deltaCnt[v] < 65535) deltaCnt[v]++;
                    nz_this_s++;
                }
            }
        }

        std::printf("SAT_S graph=%s s=%d leaves=%zu n_active=%ld kappa_max=%.0f "
                    "nz_delta=%llu build_ms=%lld peel_ms=%lld\n",
                    gpath, s, d.numLeaves, n_active_s, kappa_max,
                    (unsigned long long)nz_this_s, build_ms, peel_ms);
        std::fflush(stdout);

        std::swap(prevK, curK);
        omega = s;
        if (s == smax_cap) omega_trunc = true;
    }

    if (dumpf) std::fclose(dumpf);

    // ---- summary / compression accounting ----
    // All-vertex (conservative): fuzzy/clamped transitions counted as
    // incompressible nonzero deltas, so this is a lower bound on compression.
    long     nActive = 0;
    uint64_t trajLenSum = 0, satCount = 0, nzDeltasAll = 0;
    // Exact-only: vertices whose ENTIRE trajectory stayed below 2^53 — the
    // trustworthy regime (equals all-vertex when fuzzy==0).
    long     nActiveExact = 0;
    uint64_t trajLenExact = 0, satCountExact = 0, nzDeltasExact = 0;
    for (int v = 0; v < n; ++v) {
        if (sstar[v] < 2) continue;
        uint64_t nzv = deltaCnt[v];
        nActive++;
        trajLenSum  += (uint64_t)(sstar[v] - 1);   // entries for s=2..s*
        nzDeltasAll += nzv;
        if (nzv == 0 && !fuzzyVtx[v]) satCount++;
        if (!fuzzyVtx[v]) {
            nActiveExact++;
            trajLenExact  += (uint64_t)(sstar[v] - 1);
            nzDeltasExact += nzv;
            if (nzv == 0) satCountExact++;
        }
    }
    nzDeltasAll += fuzzyTrans;                      // incompressible fuzzy transitions
    double satFrac    = nActive ? (double)satCount / (double)nActive : 0.0;
    double naiveStore = (double)trajLenSum;
    double adStore    = 2.0 * (double)nActive + (double)nzDeltasAll;
    double compRatio  = adStore > 0 ? naiveStore / adStore : 0.0;

    double satFracEx  = nActiveExact ? (double)satCountExact / (double)nActiveExact : 0.0;
    double naiveEx    = (double)trajLenExact;
    double adEx       = 2.0 * (double)nActiveExact + (double)nzDeltasExact;
    double compEx     = adEx > 0 ? naiveEx / adEx : 0.0;

    auto t_all1 = std::chrono::high_resolution_clock::now();
    long long total_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t_all1 - t_all0).count();

    std::printf("SAT_SUMMARY graph=%s n=%d m=%ld maxdeg=%ld omega=%d omega_trunc=%d "
                "n_active=%ld saturated=%llu sat_frac=%.4f nonzero_deltas=%llu "
                "traj_len_sum=%llu naive_store=%.0f anchordelta_store=%.0f "
                "comp_ratio=%.3f "
                "n_active_exact=%ld sat_frac_exact=%.4f comp_ratio_exact=%.3f "
                "law_viol=%llu nest_viol=%llu fuzzy=%llu total_ms=%lld\n",
                gpath, n, m, maxdeg, omega, omega_trunc ? 1 : 0,
                nActive, (unsigned long long)satCount, satFrac,
                (unsigned long long)nzDeltasAll, (unsigned long long)trajLenSum,
                naiveStore, adStore, compRatio,
                nActiveExact, satFracEx, compEx,
                (unsigned long long)lawViol, (unsigned long long)nestViol,
                (unsigned long long)fuzzy, total_ms);
    std::fflush(stdout);
    return 0;
}
