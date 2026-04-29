// LowMemProbes.h — observation-only probes for RegionCPI LowMem.
// Header-only. Each probe is gated on its own env var inside the call site.
//
// The probes do not modify any algorithm state — they only print diagnostics.
// Keep this file free of production-path code so the main .cpp stays focused.
//
// Probes:
//   probe_vsafe                — PIVOTER_VSAFE_PROBE
//   probe_vsafe_gap_diag       — PIVOTER_VSAFE_GAP_DIAG (B1 vs B2 gap)
//   probe_homcc_tuple<...>     — PIVOTER_HOMCC_PROBE (tuple-level HOMCC)

#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../Global/Global.h"

// nCr table is owned by the framework; extern here so probes can read it.
extern double nCr[1001][401];

namespace pivoter::lowmem_probes {

// Sorted-list intersection size — shared utility for all probes.
inline int region_overlap(const std::vector<daf::Size> &A,
                          const std::vector<daf::Size> &B) {
    size_t i = 0, j = 0;
    int cnt = 0;
    while (i < A.size() && j < B.size()) {
        if      (A[i] < B[j]) ++i;
        else if (A[i] > B[j]) ++j;
        else                  { ++cnt; ++i; ++j; }
    }
    return cnt;
}

// V_safe probe: per-MC vertex shares (safe%, strong%, priv%, INCREMENTAL).
// For each active MC M, compute:
//   exposed(M) = vertices of M that are also in some high-overlap-with-M
//                partner M' (|M ∩ M'| ≥ r).
//   safe(M)    = M \ exposed(M).
// Reports r-clique catch counts: safeOnly (R⊆safe), strong/safeTouching
// (R∩safe≠∅), currentPrivate (frame's existing). The INCREMENTAL row is
// strong−private — what V_safe could gain over the framework's private cloud.
inline void probe_vsafe(
    daf::CliqueSize r,
    daf::Size numVertices,
    daf::Size numRegions,
    const std::vector<std::vector<daf::Size>> &regionVerts,
    const std::vector<std::vector<daf::Size>> &vtxMaxPaths) {
    auto t0 = std::chrono::high_resolution_clock::now();
    long long sum_M_size = 0, sum_exposed = 0, sum_nonPrivate = 0;
    long long sum_M_rcliques = 0, sum_safeOnly = 0, sum_safeTouching = 0;
    long long sum_currentPrivate = 0, sum_safe_size = 0;

    std::vector<int>       ovBuf(numRegions, 0);
    std::vector<bool>      hiBuf(numRegions, false);
    std::vector<daf::Size> dirty;
    dirty.reserve(256);

    auto Cnk = [&](size_t n, int k) -> long long {
        if (k < 0 || (size_t)k > n) return 0;
        if (n > 1000 || k > 400) return 0;
        return (long long)llround(nCr[n][k]);
    };

    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        for (daf::Size v : regionVerts[rid]) {
            if (v >= numVertices) continue;
            for (daf::Size other : vtxMaxPaths[v]) {
                if (other == rid) continue;
                if (ovBuf[other] == 0) dirty.push_back(other);
                ++ovBuf[other];
            }
        }
        for (auto d : dirty)
            if (ovBuf[d] >= (int)r) hiBuf[d] = true;

        size_t exposedCnt = 0, nonPrivCnt = 0;
        for (daf::Size v : regionVerts[rid]) {
            if (v >= numVertices) continue;
            bool isExp = false;
            for (daf::Size other : vtxMaxPaths[v])
                if (other != rid && hiBuf[other]) { isExp = true; break; }
            if (isExp) ++exposedCnt;
            if (vtxMaxPaths[v].size() > 1) ++nonPrivCnt;
        }
        size_t Msize    = regionVerts[rid].size();
        size_t safeSize = (Msize >= exposedCnt) ? (Msize - exposedCnt) : 0;

        long long c_M_r       = Cnk(Msize,      (int)r);
        long long c_exposed_r = Cnk(exposedCnt, (int)r);
        long long c_nonPriv_r = Cnk(nonPrivCnt, (int)r);
        long long c_safe_r    = Cnk(safeSize,   (int)r);

        sum_M_size         += (long long)Msize;
        sum_exposed        += (long long)exposedCnt;
        sum_nonPrivate     += (long long)nonPrivCnt;
        sum_safe_size      += (long long)safeSize;
        sum_M_rcliques     += c_M_r;
        sum_safeOnly       += c_safe_r;
        sum_safeTouching   += c_M_r - c_exposed_r;
        sum_currentPrivate += c_M_r - c_nonPriv_r;

        for (auto d : dirty) { ovBuf[d] = 0; hiBuf[d] = false; }
        dirty.clear();
    }

    long long incremental = sum_safeTouching - sum_currentPrivate;
    auto t1 = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    auto pct = [&](long long n, long long d) -> double {
        return d > 0 ? 100.0 * (double)n / (double)d : 0.0;
    };
    std::cout << "  [V_safe Probe] " << ms << " ms" << std::endl;
    std::cout << "    Active MCs (post Step1b):  " << numRegions << std::endl;
    std::cout << "    Σ|M|:                      " << sum_M_size << std::endl;
    std::cout << "    Σ|exposed| / Σ|M|:         " << sum_exposed
              << " (" << std::fixed << std::setprecision(1)
              << pct(sum_exposed, sum_M_size) << "%)" << std::endl;
    std::cout << "    Σ|nonPrivate| / Σ|M|:      " << sum_nonPrivate
              << " (" << pct(sum_nonPrivate, sum_M_size) << "%)" << std::endl;
    std::cout << "    Σ|safe| / Σ|M|:            " << sum_safe_size
              << " (" << pct(sum_safe_size, sum_M_size) << "%)" << std::endl;
    std::cout << "    --- r-clique counts (within active MCs) ---" << std::endl;
    std::cout << "    ΣC(|M|, r):                " << sum_M_rcliques << std::endl;
    std::cout << "    safeOnly  (R ⊆ safe):      " << sum_safeOnly
              << " (" << pct(sum_safeOnly, sum_M_rcliques) << "%)" << std::endl;
    std::cout << "    safeTouching (R∩safe≠∅):   " << sum_safeTouching
              << " (" << pct(sum_safeTouching, sum_M_rcliques) << "%)" << std::endl;
    std::cout << "    currentPrivate (frame):    " << sum_currentPrivate
              << " (" << pct(sum_currentPrivate, sum_M_rcliques) << "%)" << std::endl;
    std::cout << "    INCREMENTAL = strong - private = " << incremental
              << " (" << pct(incremental, sum_M_rcliques)
              << "% of total in M)" << std::endl;
}

// V_safe gap diagnostic: decompose Σ|M| into {private, globalsafe, permcsafe}
// vertex groups and report r-clique catch counts per group. Useful for
// quantifying B1 (globally-safe class) vs B2 (per-MC V_safe) gap.
inline void probe_vsafe_gap_diag(
    daf::CliqueSize r,
    daf::Size numVertices,
    daf::Size numRegions,
    const std::vector<std::vector<daf::Size>> &regionVerts,
    const std::vector<std::vector<daf::Size>> &vtxMaxPaths) {
    auto pairOv = [&](daf::Size r1, daf::Size r2) -> int {
        return region_overlap(regionVerts[r1], regionVerts[r2]);
    };

    // Per-vertex globally-safe flag (private OR all-pairwise-overlap < r).
    std::vector<uint8_t> isGloballySafeV(numVertices, 0);
    for (daf::Size v = 0; v < numVertices; ++v) {
        const auto &mcs = vtxMaxPaths[v];
        if (mcs.empty()) continue;
        if (mcs.size() == 1) { isGloballySafeV[v] = 1; continue; }
        bool allLow = true;
        for (size_t i = 0; i < mcs.size() && allLow; ++i)
            for (size_t j = i + 1; j < mcs.size(); ++j)
                if (pairOv(mcs[i], mcs[j]) >= (int)r) {
                    allLow = false;
                    break;
                }
        isGloballySafeV[v] = allLow ? 1 : 0;
    }

    std::vector<int>       ovBuf(numRegions, 0);
    std::vector<bool>      hiBuf(numRegions, false);
    std::vector<daf::Size> dirty;
    dirty.reserve(256);

    long long sum_private = 0, sum_globalsafe = 0, sum_permcsafe = 0;
    long long sum_M_total = 0, sum_M_rcliques = 0;
    long long catch_private = 0, catch_globalsafe = 0, catch_permcsafe = 0;

    auto Cnk = [&](size_t n, int k) -> long long {
        if (k < 0 || (size_t)k > n) return 0;
        if (n > 1000 || k > 400) return 0;
        return (long long)llround(nCr[n][k]);
    };

    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        for (daf::Size v : regionVerts[rid]) {
            if (v >= numVertices) continue;
            for (daf::Size other : vtxMaxPaths[v]) {
                if (other == rid) continue;
                if (ovBuf[other] == 0) dirty.push_back(other);
                ++ovBuf[other];
            }
        }
        for (auto d : dirty)
            if (ovBuf[d] >= (int)r) hiBuf[d] = true;

        size_t cnt_private = 0, cnt_globalsafe = 0, cnt_permcsafe = 0;
        for (daf::Size v : regionVerts[rid]) {
            if (v >= numVertices) continue;
            bool priv = (vtxMaxPaths[v].size() == 1);
            bool gs   = (isGloballySafeV[v] != 0);
            bool exposed_in_M = false;
            for (daf::Size other : vtxMaxPaths[v])
                if (other != rid && hiBuf[other]) { exposed_in_M = true; break; }
            bool permcsafe = !exposed_in_M;
            if (priv) ++cnt_private;
            if (gs) ++cnt_globalsafe;
            if (permcsafe) ++cnt_permcsafe;
        }

        size_t Msize        = regionVerts[rid].size();
        size_t exp_private  = Msize - cnt_private;
        size_t exp_global   = Msize - cnt_globalsafe;
        size_t exp_permcsafe = Msize - cnt_permcsafe;
        long long c_M       = Cnk(Msize,         (int)r);
        long long c_no_priv = Cnk(exp_private,   (int)r);
        long long c_no_gs   = Cnk(exp_global,    (int)r);
        long long c_no_pms  = Cnk(exp_permcsafe, (int)r);

        sum_M_total      += (long long)Msize;
        sum_private      += (long long)cnt_private;
        sum_globalsafe   += (long long)cnt_globalsafe;
        sum_permcsafe    += (long long)cnt_permcsafe;
        sum_M_rcliques   += c_M;
        catch_private    += c_M - c_no_priv;
        catch_globalsafe += c_M - c_no_gs;
        catch_permcsafe  += c_M - c_no_pms;

        for (auto d : dirty) { ovBuf[d] = 0; hiBuf[d] = false; }
        dirty.clear();
    }

    long long b1_extra = catch_globalsafe - catch_private;
    long long b2_extra = catch_permcsafe  - catch_globalsafe;
    auto pct = [&](long long a, long long b) -> double {
        return b > 0 ? 100.0 * a / b : 0.0;
    };
    std::cout << "  [V_safe Gap Diag] (B1 vs B2 gap source)" << std::endl;
    std::cout << "    Σ|M|=" << sum_M_total
              << " | private="    << sum_private    << " (" << pct(sum_private, sum_M_total) << "%)"
              << "  globalsafe="  << sum_globalsafe << " (" << pct(sum_globalsafe, sum_M_total) << "%)"
              << "  permcsafe="   << sum_permcsafe  << " (" << pct(sum_permcsafe, sum_M_total) << "%)"
              << std::endl;
    std::cout << "    catches: private=" << catch_private
              << " ("    << pct(catch_private, sum_M_rcliques) << "%) "
              << "globalsafe=" << catch_globalsafe << " (" << pct(catch_globalsafe, sum_M_rcliques) << "%) "
              << "permcsafe="  << catch_permcsafe  << " (" << pct(catch_permcsafe, sum_M_rcliques) << "%)"
              << std::endl;
    std::cout << "    B1 extra over private (B1 actual gain) = "
              << b1_extra << " (" << pct(b1_extra, sum_M_rcliques) << "% of M-rcliques)" << std::endl;
    std::cout << "    B2 extra over B1 (Plan B2 potential)   = "
              << b2_extra << " (" << pct(b2_extra, sum_M_rcliques) << "% of M-rcliques)" << std::endl;
}

// HOMCC tuple-level probe. Templated over Classes/RTuples to avoid leaking
// the internal struct types out of the .cpp.
//   Classes : random-access of T with `.regionIds` (sorted ascending)
//   RTuples : random-access of U with `.mult`
//   KeyOfFn : (daf::Size tidx) -> std::span<const daf::Size>
//
// Classify tuples by their MCs(t) = ∩ classes[c].regionIds:
//   singleton(|MCs|=1), vsafe(all <r), HOMCC(all ≥r), mixed.
// (By construction, vsafe and mixed are empty for active multi-MC tuples,
//  since |M_i ∩ M_j| ≥ Σ classSize ≥ r whenever both M_i,M_j ∈ MCs(t).
//  Probe still computes them so the invariant is verified empirically.)
template <typename Classes, typename RTuples, typename KeyOfFn>
void probe_homcc_tuple(
    daf::CliqueSize r,
    const std::vector<std::vector<daf::Size>> &regionVerts,
    const Classes &classes,
    const RTuples &rTuples,
    KeyOfFn keyOf) {
    auto pairOv = [&](daf::Size r1, daf::Size r2) -> int {
        return region_overlap(regionVerts[r1], regionVerts[r2]);
    };

    long long n_singleton = 0, n_vsafe = 0, n_homc = 0, n_mixed = 0;
    long long rc_singleton = 0, rc_vsafe = 0, rc_homc = 0, rc_mixed = 0;
    long long sumMCs_h = 0, maxMCs_h = 0;
    long long sumMCs_m = 0, maxMCs_m = 0;

    std::vector<daf::Size> mcsBuf, mcsTmp, distinct;

    for (daf::Size tidx = 0; tidx < rTuples.size(); ++tidx) {
        long long mult = (long long)rTuples[tidx].mult;
        auto key = keyOf(tidx);

        distinct.clear();
        distinct.reserve(key.size());
        for (auto c : key) distinct.push_back(c);
        std::sort(distinct.begin(), distinct.end());
        distinct.erase(std::unique(distinct.begin(), distinct.end()),
                       distinct.end());
        if (distinct.empty()) continue;

        // MCs(t) = ∩ classes[c].regionIds (multi-way sorted intersection).
        mcsBuf = classes[distinct[0]].regionIds;
        for (size_t i = 1; i < distinct.size() && !mcsBuf.empty(); ++i) {
            const auto &B = classes[distinct[i]].regionIds;
            mcsTmp.clear();
            size_t a = 0, b = 0;
            while (a < mcsBuf.size() && b < B.size()) {
                if      (mcsBuf[a] < B[b]) ++a;
                else if (mcsBuf[a] > B[b]) ++b;
                else { mcsTmp.push_back(mcsBuf[a]); ++a; ++b; }
            }
            mcsBuf.swap(mcsTmp);
        }

        size_t k = mcsBuf.size();
        if (k == 0) continue;
        if (k == 1) { ++n_singleton; rc_singleton += mult; continue; }

        bool anyHi = false, anyLo = false;
        for (size_t i = 0; i < k && !(anyHi && anyLo); ++i) {
            for (size_t j = i + 1; j < k; ++j) {
                int ov = pairOv(mcsBuf[i], mcsBuf[j]);
                if (ov >= (int)r) anyHi = true; else anyLo = true;
                if (anyHi && anyLo) break;
            }
        }
        if (!anyHi && anyLo) { ++n_vsafe; rc_vsafe += mult; }
        else if (anyHi && !anyLo) {
            ++n_homc; rc_homc += mult;
            sumMCs_h += k; if ((long long)k > maxMCs_h) maxMCs_h = k;
        } else {
            ++n_mixed; rc_mixed += mult;
            sumMCs_m += k; if ((long long)k > maxMCs_m) maxMCs_m = k;
        }
    }

    long long activeRC = rc_singleton + rc_vsafe + rc_homc + rc_mixed;
    auto pct = [&](long long a, long long b) -> double {
        return b > 0 ? 100.0 * a / b : 0.0;
    };
    std::cout << "  [HOMCC Probe TUPLE] |tuples|=" << rTuples.size()
              << "  singleton=" << n_singleton
              << " vsafe="     << n_vsafe
              << " HOMCC="     << n_homc
              << " mixed="     << n_mixed << std::endl;
    std::cout << "    r-cliques active=" << activeRC
              << "  singleton=" << rc_singleton << " (" << pct(rc_singleton, activeRC) << "%)"
              << "  vsafe="     << rc_vsafe     << " (" << pct(rc_vsafe, activeRC) << "%)"
              << "  HOMCC="     << rc_homc      << " (" << pct(rc_homc, activeRC) << "%)"
              << "  mixed="     << rc_mixed     << " (" << pct(rc_mixed, activeRC) << "%)"
              << std::endl;
    if (n_homc > 0) {
        std::cout << "    HOMCC tuple |MC|: avg="
                  << ((double)sumMCs_h / n_homc) << " max=" << maxMCs_h << std::endl;
    }
    if (n_mixed > 0) {
        std::cout << "    mixed tuple |MC|: avg="
                  << ((double)sumMCs_m / n_mixed) << " max=" << maxMCs_m << std::endl;
    }
}

}  // namespace pivoter::lowmem_probes
