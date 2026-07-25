// test_scalable_sct.cpp — gate harness for the SPARSE degeneracy-seeded class
// SCT (ClassSCTScalable.h) against the VALIDATED dense buildClassSCT
// (ClassSCT.h, the correctness oracle).
//
//   GATE 1 (correctness): thousands of random small graphs (C 3..9, w 1..4,
//           k 2..5). Build BOTH the dense buildClassSCT (sparse->dense ClassG)
//           and scalableBuildClassSCT; assert sum-of-support_count EQUAL and
//           per-pattern disjointness (reusing class_sct.cpp's checkDisjoint)
//           holds for the scalable one. 0 fail required.
//   GATE 2 (scale): large random SPARSE quotient graphs (C 50k..125k, avg
//           degree ~10-50, planted cliques). Run scalableBuildClassSCT for
//           k=4; report build time, leaf count, peak RSS. No C x C matrix.
//   GATE 3 (degenerate consistency): a complete graph K_C returns the same
//           total as dense (binomial).
//
// Build (from region_native/):
//   g++ -O3 -std=c++17 -I../src/NucleusDecomposition -o test_scalable_sct test_scalable_sct.cpp
// Run:  ./test_scalable_sct

#include "../src/NucleusDecomposition/CCPathCore.h"
#include "ClassSCT.h"            // dense oracle: buildClassSCT(ClassG, k)
#include "ClassSCTScalable.h"    // unit under test: scalableBuildClassSCT(...)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <random>
#include <vector>

#if defined(__APPLE__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <sys/resource.h>
#endif

using ccpath::CCPath;
using ccpath::Vec;

// ---------------- Pascal-table nCr (double), n up to 256 ----------------
// Larger table than class_sct.cpp's (64) so K_C consistency with big weights
// and big planted cliques does not overflow the table.
static const int PMAX = 256;
static std::vector<std::vector<double>> PASCAL;
static void initPascal() {
    PASCAL.assign(PMAX + 1, std::vector<double>(PMAX + 1, 0.0));
    for (int n = 0; n <= PMAX; ++n) {
        PASCAL[n][0] = 1.0;
        for (int kk = 1; kk <= n; ++kk)
            PASCAL[n][kk] = PASCAL[n - 1][kk - 1] + (kk <= n - 1 ? PASCAL[n - 1][kk] : 0.0);
    }
}
static double nCr_fn(int n, int kk) {
    if (kk < 0 || n < 0 || kk > n || n > PMAX) return 0.0;
    return PASCAL[n][kk];
}

// ---------------- brute-force ground truth (countCanon), id-order ----------
// Identical to class_sct.cpp's countCanon, but parameterised on a local
// (dense) ClassG so it does not collide with the dense header's gG.
static const ClassG* gGB = nullptr;
static double countCanon(const std::vector<int>& avail, int need) {
    if (need == 0) return 1.0;
    double res = 0.0;
    for (size_t i = 0; i < avail.size(); ++i) {
        int c = avail[i];
        std::vector<int> later;
        for (size_t t = i + 1; t < avail.size(); ++t)
            if (gGB->adj(c, avail[t])) later.push_back(avail[t]);
        int hi = std::min(need, gGB->w[c]);
        for (int j = 1; j <= hi; ++j)
            res += nCr_fn(gGB->w[c], j) * countCanon(later, need - j);
    }
    return res;
}
static double bruteTotal(const ClassG& G, int k) {
    gGB = &G;
    std::vector<int> all;
    for (int c = 0; c < G.C; ++c) if (G.w[c] > 0) all.push_back(c);
    return countCanon(all, k);
}

static double sctTotal(const std::vector<CCPath>& leaves) {
    double tot = 0.0;
    for (const auto& L : leaves) {
        Vec b = ccpath::zeros_vec(L.m());
        tot += ccpath::support_count(L, b, nCr_fn);
    }
    return tot;
}

// ---------------- PER-PATTERN DISJOINTNESS verifier ----------------
// Verbatim port of class_sct.cpp's checkDisjoint: every weighted-clique
// pattern of total k must be represented EXACTLY once across the leaves.
static bool checkDisjoint(const ClassG& G, int k,
                          const std::vector<CCPath>& leaves,
                          std::vector<int>& badPattern /*out*/) {
    const int C = G.C;
    std::vector<int> m(C, 0);
    bool ok = true;
    std::vector<int> chosen;
    std::function<void(int,int)> rec = [&](int c, int rem) {
        if (rem == 0) {
            double truth = 1.0;
            for (int x : chosen) truth *= nCr_fn(G.w[x], m[x]);
            double rep = 0.0;
            for (const auto& L : leaves) {
                if (L.T != k) continue;
                bool adm = true; int sm = 0;
                for (int cc = 0; cc < C; ++cc) {
                    int mc = m[cc]; sm += mc;
                    if (mc < (int)L.ell[cc] || mc > (int)L.u[cc]) { adm = false; break; }
                }
                if (!adm || sm != k) continue;
                double prod = 1.0;
                for (int cc = 0; cc < C; ++cc)
                    if (m[cc] > 0) prod *= nCr_fn((int)L.n[cc], m[cc]);
                rep += prod;
            }
            if (std::abs(rep - truth) > 0.5) {
                ok = false;
                if (badPattern.empty()) badPattern = m;
            }
            return;
        }
        if (c >= C) return;
        rec(c + 1, rem);
        bool compat = true;
        for (int x : chosen) if (!G.adj(c, x)) { compat = false; break; }
        if (compat) {
            int hi = std::min(rem, G.w[c]);
            chosen.push_back(c);
            for (int j = 1; j <= hi; ++j) { m[c] = j; rec(c + 1, rem - j); }
            m[c] = 0; chosen.pop_back();
        }
    };
    rec(0, k);
    return ok;
}

// Expand a COMPACT scalable leaf (dimension = #touched classes, with
// classIds[i] = global id of slot i) into a FULL-dimension-C leaf, so the
// global-id-indexed dense disjointness verifier (checkDisjoint) can consume
// it. support_count is unaffected by this expansion (absent classes carry
// n=ell=u=0 => a C(0,0)=1 factor), which is exactly why the compact form is
// legal; we expand here only for the per-pattern verifier's convenience.
static CCPath expandToFull(const CCPath& cp, int C) {
    CCPath f;
    f.h = ccpath::zeros_vec(C);
    f.n = ccpath::zeros_vec(C);
    f.ell = ccpath::zeros_vec(C);
    f.u = ccpath::zeros_vec(C);
    f.T = cp.T;
    for (int i = 0; i < cp.m(); ++i) {
        int g = (int)cp.classIds[i];          // global class id of slot i
        f.h[g] = cp.h[i]; f.n[g] = cp.n[i];
        f.ell[g] = cp.ell[i]; f.u[g] = cp.u[i];
    }
    // forbidden antichains are empty for freshly-built leaves; if any existed
    // they would also need remapping, but build never sets them.
    return f;
}

// sparse adjacency -> dense ClassG (for the oracle + disjointness verifier).
static ClassG toDense(int C, const std::vector<int>& w,
                      const std::vector<std::vector<int>>& adj) {
    ClassG G; G.C = C; G.w = w;
    G.A.assign(C, std::vector<char>(C, 0));
    for (int i = 0; i < C; ++i)
        for (int j : adj[i]) if (j != i) G.A[i][j] = 1;
    return G;
}

// peak resident set size in MB (0 if unsupported).
static double peakRSS_MB() {
#if defined(__APPLE__)
    mach_task_basic_info info;
    mach_msg_type_number_t cnt = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  (task_info_t)&info, &cnt) == KERN_SUCCESS)
        return (double)info.resident_size_max / (1024.0 * 1024.0);
    return 0.0;
#elif defined(__linux__)
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_maxrss / 1024.0;   // ru_maxrss is KB on Linux
#else
    return 0.0;
#endif
}

// ===================================================================
//  GATE 1 — random small graphs, scalable vs dense oracle
// ===================================================================
static int gate1() {
    printf("==== GATE 1: scalable SCT vs dense oracle (random small graphs) ====\n");
    std::mt19937 rng(20260618u);
    const int trials = 50000;
    int sumPass = 0, sumFail = 0, disjPass = 0, disjFail = 0;
    int firstFail = 0, firstDisj = 0;
    long long totLeavesScal = 0, totLeavesDense = 0;

    for (int t = 0; t < trials; ++t) {
        int C = 3 + (int)(rng() % 7);              // 3..9
        std::vector<int> w(C);
        for (int c = 0; c < C; ++c) w[c] = 1 + (int)(rng() % 4);  // 1..4
        // full density spectrum incl. empty & complete graphs.
        double pEdge = (double)(rng() % 101) / 100.0;
        std::vector<std::vector<int>> adj(C);
        for (int i = 0; i < C; ++i)
            for (int j = i + 1; j < C; ++j)
                if ((double)(rng() % 1000) / 1000.0 < pEdge) {
                    adj[i].push_back(j); adj[j].push_back(i);
                }
        for (int i = 0; i < C; ++i) std::sort(adj[i].begin(), adj[i].end());
        int k = 2 + (int)(rng() % 4);              // 2..5

        ClassG Gd = toDense(C, w, adj);
        std::vector<CCPath> dense = buildClassSCT(Gd, k);
        std::vector<CCPath> scal  =
            classsct_scalable::scalableBuildClassSCT(C, w, adj, k);
        totLeavesDense += (long long)dense.size();
        totLeavesScal  += (long long)scal.size();

        double sDense = sctTotal(dense);
        double sScal  = sctTotal(scal);
        double bf     = bruteTotal(Gd, k);   // independent ground truth too

        // primary gate: scalable sum == dense sum (and both == brute).
        if (std::abs(sScal - sDense) < 0.5 && std::abs(sScal - bf) < 0.5) ++sumPass;
        else {
            ++sumFail;
            if (firstFail++ < 10) {
                printf("SUM-FAIL t=%d C=%d k=%d  brute=%.1f dense=%.1f scal=%.1f\n",
                       t, C, k, bf, sDense, sScal);
                printf("  w=["); for (int c = 0; c < C; ++c) printf("%d ", w[c]);
                printf("]  adj:\n");
                for (int i = 0; i < C; ++i) {
                    printf("    %d:", i);
                    for (int j : adj[i]) printf(" %d", j);
                    printf("\n");
                }
            }
        }

        // scalable leaves are COMPACT; expand to full-C for the global-id-
        // indexed per-pattern disjointness verifier.
        std::vector<CCPath> scalFull; scalFull.reserve(scal.size());
        for (const auto& L : scal) scalFull.push_back(expandToFull(L, C));
        std::vector<int> bad;
        if (checkDisjoint(Gd, k, scalFull, bad)) ++disjPass;
        else {
            ++disjFail;
            if (firstDisj++ < 10) {
                printf("DISJ-FAIL t=%d C=%d k=%d  pattern m=[", t, C, k);
                for (int c = 0; c < C; ++c) printf("%d ", bad[c]);
                printf("]  w=["); for (int c = 0; c < C; ++c) printf("%d ", w[c]);
                printf("]\n");
            }
        }
    }

    printf("[sum   scal==dense==brute] trials=%d PASS=%d FAIL=%d\n", trials, sumPass, sumFail);
    printf("[disj  per-pattern once  ] trials=%d PASS=%d FAIL=%d\n", trials, disjPass, disjFail);
    printf("[leaves] dense avg=%.2f  scalable avg=%.2f\n",
           (double)totLeavesDense / trials, (double)totLeavesScal / trials);
    bool ok = (sumFail == 0 && disjFail == 0);
    printf("GATE 1: %s\n\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}

// ===================================================================
//  GATE 3 — complete graph K_C consistency (binomial), dense vs scalable
// ===================================================================
static int gate3() {
    printf("==== GATE 3: complete-graph K_C consistency ====\n");
    int fails = 0;
    // K_C with unit weights, k-clique count must be C(C,k); also test weighted.
    struct Case { int C; int wuni; int k; };
    std::vector<Case> cases = { {6,1,3}, {8,1,4}, {10,1,5}, {5,2,3}, {7,3,4}, {12,1,6} };
    for (auto& cs : cases) {
        int C = cs.C;
        std::vector<int> w(C, cs.wuni);
        std::vector<std::vector<int>> adj(C);
        for (int i = 0; i < C; ++i)
            for (int j = 0; j < C; ++j) if (i != j) adj[i].push_back(j);
        for (int i = 0; i < C; ++i) std::sort(adj[i].begin(), adj[i].end());

        ClassG Gd = toDense(C, w, adj);
        double sDense = sctTotal(buildClassSCT(Gd, cs.k));
        double sScal  = sctTotal(classsct_scalable::scalableBuildClassSCT(C, w, adj, cs.k));
        double bf     = bruteTotal(Gd, cs.k);
        // closed form for unit weight: C(C, k).
        double closed = (cs.wuni == 1) ? nCr_fn(C, cs.k) : -1.0;
        bool ok = std::abs(sDense - sScal) < 0.5 && std::abs(sScal - bf) < 0.5
                  && (closed < 0 || std::abs(sScal - closed) < 0.5);
        if (!ok) ++fails;
        printf("  K_%d w=%d k=%d  dense=%.0f scal=%.0f brute=%.0f%s  %s\n",
               C, cs.wuni, cs.k, sDense, sScal, bf,
               closed >= 0 ? (" closed=" + std::to_string((long long)closed)).c_str() : "",
               ok ? "[OK]" : "[FAIL]");
    }
    printf("GATE 3: %s\n\n", fails == 0 ? "PASS" : "FAIL");
    return fails == 0 ? 0 : 1;
}

// ===================================================================
//  GATE 2 — SCALE: large sparse quotient graphs, scalable only
// ===================================================================
static int gate2() {
    printf("==== GATE 2: SCALE (large sparse quotient graphs, k=4) ====\n");
    printf("(scalable builder only; dense C x C is intentionally never built)\n");
    std::mt19937_64 rng(424242ull);

    struct Cfg { int C; int avgDeg; int planted; int pSize; };
    // C in 50k..125k; avg degree ~10-50; a handful of planted cliques to make
    // the SCT do real pivoting work (and exercise dense local neighbourhoods).
    std::vector<Cfg> cfgs = {
        { 50000, 10, 30, 12 },
        { 75000, 20, 40, 14 },
        {100000, 30, 50, 16 },
        {123526, 40, 60, 16 },   // ~ com-dblp nC
        {125000, 50, 50, 18 },
    };
    int k = 4;
    int fails = 0;

    for (auto& cf : cfgs) {
        int C = cf.C;
        std::vector<int> w(C);
        for (int c = 0; c < C; ++c) w[c] = 1 + (int)(rng() % 4);   // 1..4

        // Build a sparse symmetric graph as an edge set, then CSR-style adj.
        // Random Erdos-Renyi-ish backbone (avgDeg/2 random edges per node) plus
        // planted cliques. Use a set per node only transiently via final sort
        // + unique to dedup. Memory stays O(C + edges).
        std::vector<std::vector<int>> adj(C);
        long long targetEdges = (long long)C * cf.avgDeg / 2;
        std::uniform_int_distribution<int> pick(0, C - 1);
        for (long long e = 0; e < targetEdges; ++e) {
            int a = pick(rng), b = pick(rng);
            if (a == b) continue;
            adj[a].push_back(b); adj[b].push_back(a);
        }
        // planted cliques (random vertex subsets, fully connected)
        long long plantedEdges = 0;
        for (int p = 0; p < cf.planted; ++p) {
            std::vector<int> S(cf.pSize);
            for (int i = 0; i < cf.pSize; ++i) S[i] = pick(rng);
            std::sort(S.begin(), S.end()); S.erase(std::unique(S.begin(), S.end()), S.end());
            for (size_t i = 0; i < S.size(); ++i)
                for (size_t j = i + 1; j < S.size(); ++j) {
                    adj[S[i]].push_back(S[j]); adj[S[j]].push_back(S[i]);
                    ++plantedEdges;
                }
        }
        // sort+dedup neighbour lists; count edges + degeneracy proxy (maxdeg).
        long long edges = 0; int maxDeg = 0;
        for (int c = 0; c < C; ++c) {
            auto& a = adj[c];
            std::sort(a.begin(), a.end());
            a.erase(std::unique(a.begin(), a.end()), a.end());
            // drop accidental self
            a.erase(std::remove(a.begin(), a.end(), c), a.end());
            edges += (long long)a.size();
            if ((int)a.size() > maxDeg) maxDeg = (int)a.size();
        }
        edges /= 2;

        double rss0 = peakRSS_MB();
        auto t0 = std::chrono::high_resolution_clock::now();
        std::vector<CCPath> leaves =
            classsct_scalable::scalableBuildClassSCT(C, w, adj, k);
        auto t1 = std::chrono::high_resolution_clock::now();
        double build = std::chrono::duration<double>(t1 - t0).count();
        double rss1 = peakRSS_MB();

        // sanity: total weighted-4-clique count must be a finite non-negative
        // number; print it so it can be cross-checked. (Full per-pattern
        // disjointness at this scale is covered by gate 1; here we confirm it
        // runs in seconds, stays small in memory, and produces a real count.)
        double total = sctTotal(leaves);
        bool ok = (build < 60.0) && std::isfinite(total) && total >= 0.0;
        if (!ok) ++fails;
        printf("  C=%-7d edges=%-9lld maxDeg=%-4d planted=%dx%d  "
               "build=%6.2fs leaves=%-9zu total4clq=%.0f  peakRSS=%.0fMB  %s\n",
               C, edges, maxDeg, cf.planted, cf.pSize,
               build, leaves.size(), total, rss1,
               ok ? "[OK]" : "[SLOW/BAD]");
        fflush(stdout);
        (void)rss0;
    }
    printf("GATE 2: %s  (peak RSS must be far below a C x C matrix's ~%.1fGB)\n\n",
           fails == 0 ? "PASS" : "FAIL",
           (double)123526.0 * 123526.0 / (1024.0*1024.0*1024.0));
    return fails == 0 ? 0 : 1;
}

// ===================================================================
//  GATE 2b — SCALE CORRECTNESS: independent SPARSE ground-truth total.
//  Gate 2 only checks the scalable builder runs fast / small. Here we also
//  confirm its weighted-k-clique TOTAL is exactly right at a scale where the
//  dense C x C oracle would OOM, using an independent sparse counter:
//  degeneracy-seeded countCanon over sorted later-neighbour lists (no matrix).
//  This is a DIFFERENT code path (recursive count, no CCPath / pivot / pool)
//  so agreement is strong evidence of correctness, not a tautology.
// ===================================================================
// independent sparse weighted-clique counter (no dense matrix, no CCPath).
// canon over GLOBAL id order: for each class c, take j=1..w_c copies and
// recurse on c's id-later neighbours. Same model as class_sct.cpp's countCanon
// but reading the sparse adj directly. Bounded by treating later-sets via
// binary search; only feasible because the test graphs here are sparse.
static const std::vector<int>* SCW = nullptr;
static const std::vector<std::vector<int>>* SCADJ = nullptr;
static double sparseCanon(const std::vector<int>& avail, int need) {
    if (need == 0) return 1.0;
    double res = 0.0;
    for (size_t i = 0; i < avail.size(); ++i) {
        int c = avail[i];
        const std::vector<int>& Nc = (*SCADJ)[c];
        std::vector<int> later;
        for (size_t t = i + 1; t < avail.size(); ++t)
            if (std::binary_search(Nc.begin(), Nc.end(), avail[t]))
                later.push_back(avail[t]);
        int hi = std::min(need, (*SCW)[c]);
        for (int j = 1; j <= hi; ++j)
            res += nCr_fn((*SCW)[c], j) * sparseCanon(later, need - j);
    }
    return res;
}

static int gate2b() {
    printf("==== GATE 2b: SCALE CORRECTNESS (independent sparse total) ====\n");
    std::mt19937_64 rng(99991ull);
    // C large enough that a C x C oracle is out of the question, but sparse
    // enough that the independent recursive counter still finishes. Keep avg
    // degree modest and planted cliques small so sparseCanon stays tractable.
    struct Cfg { int C; int avgDeg; int planted; int pSize; int k; };
    std::vector<Cfg> cfgs = {
        { 30000,  8, 20, 8, 4 },
        { 60000, 12, 25, 9, 4 },
        { 80000, 10, 30, 7, 5 },
    };
    int fails = 0;
    for (auto& cf : cfgs) {
        int C = cf.C;
        std::vector<int> w(C);
        for (int c = 0; c < C; ++c) w[c] = 1 + (int)(rng() % 3);   // 1..3
        std::vector<std::vector<int>> adj(C);
        std::uniform_int_distribution<int> pick(0, C - 1);
        long long targetEdges = (long long)C * cf.avgDeg / 2;
        for (long long e = 0; e < targetEdges; ++e) {
            int a = pick(rng), b = pick(rng);
            if (a == b) continue;
            adj[a].push_back(b); adj[b].push_back(a);
        }
        for (int p = 0; p < cf.planted; ++p) {
            std::vector<int> S(cf.pSize);
            for (int i = 0; i < cf.pSize; ++i) S[i] = pick(rng);
            std::sort(S.begin(), S.end()); S.erase(std::unique(S.begin(), S.end()), S.end());
            for (size_t i = 0; i < S.size(); ++i)
                for (size_t j = i + 1; j < S.size(); ++j) {
                    adj[S[i]].push_back(S[j]); adj[S[j]].push_back(S[i]);
                }
        }
        for (int c = 0; c < C; ++c) {
            auto& a = adj[c];
            std::sort(a.begin(), a.end());
            a.erase(std::unique(a.begin(), a.end()), a.end());
            a.erase(std::remove(a.begin(), a.end(), c), a.end());
        }

        auto t0 = std::chrono::high_resolution_clock::now();
        double sScal = sctTotal(classsct_scalable::scalableBuildClassSCT(C, w, adj, cf.k));
        auto t1 = std::chrono::high_resolution_clock::now();

        SCW = &w; SCADJ = &adj;
        std::vector<int> all; all.reserve(C);
        for (int c = 0; c < C; ++c) if (w[c] > 0) all.push_back(c);
        auto t2 = std::chrono::high_resolution_clock::now();
        double truth = sparseCanon(all, cf.k);
        auto t3 = std::chrono::high_resolution_clock::now();

        bool ok = std::abs(sScal - truth) < 0.5;
        if (!ok) ++fails;
        printf("  C=%-6d k=%d  scalable=%.0f  sparse-truth=%.0f  "
               "(build=%.2fs truth=%.2fs)  %s\n",
               C, cf.k, sScal, truth,
               std::chrono::duration<double>(t1 - t0).count(),
               std::chrono::duration<double>(t3 - t2).count(),
               ok ? "[OK]" : "[FAIL]");
        fflush(stdout);
    }
    printf("GATE 2b: %s\n\n", fails == 0 ? "PASS" : "FAIL");
    return fails == 0 ? 0 : 1;
}

int main() {
    initPascal();
    int rc = 0;
    rc |= gate1();
    rc |= gate3();
    rc |= gate2b();
    rc |= gate2();
    printf("==== OVERALL: %s ====\n", rc == 0 ? "ALL GATES PASS" : "SOME GATE FAILED");
    return rc;
}
