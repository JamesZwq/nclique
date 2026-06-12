//
// Region Tuple + CPI Counting (V3Fast) — Low-Memory variant (V3LM).
//
// Identical algorithmic behavior to V3Fast but replaces the explicit
// tuple -> path incidence list (tupleToPathInfos) with a class -> path
// inverted index.  Tuple-to-path lookups at peel time are derived on the
// fly by intersecting class -> path lists for each distinct class in
// tau.key.  This cuts the bulk of the auxiliary memory (tuple-path pairs)
// by 50-100x on graphs with moderate/poor class compression, while adding
// a small amortised cost per peel for the sorted-list intersection.
//
// Only the "tuple -> paths" direction is eliminated; path -> alive-tuples
// (pathInfos[p].tupleIdxs) is kept because it shrinks monotonically as
// tuples peel and is the smallest live-state direction.
//
// === Environment-variable registry (all are off unless set) ==================
// Production-path flags:
//   PIVOTER_VSAFE_CLOUD       Phase B1: globally-safe class extension of the
//                             private-cloud (touches multi-MC classes whose
//                             pairwise overlaps are all < r). Strict superset
//                             of the existing private cloud; never a perf loss
//                             on graphs with mixed-overlap structure.
//   PIVOTER_V3_NO_PRIVATE     Disables private-cloud entirely (regression test).
// Ablation flags:
//   PIVOTER_DISABLE_RMERGE    Force-disables Step 1b r-mergeable filter — every
//                             MC stays in the active region. Used for the
//                             r-merge ablation in the paper.
// Tuning knobs:
//   PIVOTER_V3_FM_PAIR_LIMIT  Pair-count limit for fused-FM precompute.
//   PIVOTER_V3_FM_FORCE       Force fused-FM even when the precompute is
//                             estimated too expensive.
// Debug / verification:
//   PIVOTER_V3_DEBUG_VERIFY   Enables Step-5 path-bound cross-check.
// Observation-only probes (defined in LowMemProbes.h):
//   PIVOTER_VSAFE_PROBE       V_safe potential summary (per-MC).
//   PIVOTER_VSAFE_GAP_DIAG    B1 vs B2 gap diagnostic (private/globalsafe/permcsafe).
//   PIVOTER_HOMCC_PROBE       Tuple-level High-Overlap MC Cluster classification.
// ============================================================================
//

#include "NCliqueCoreDecomposition.h"
#include "LowMemProbes.h"
#include "../PhaseLogger.h"
#include "../dataStruct/robin_hood.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <map>
#include <span>
#include <unordered_map>
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#define REGNDC_HAVE_NEON 1
#elif defined(__SSE2__)
#include <emmintrin.h>
#define REGNDC_HAVE_SSE2 1
#endif
#include <unordered_set>
#include <vector>

extern double nCr[1001][401];
extern std::vector<bool> g_maxCliqueTags;
extern std::vector<std::vector<daf::Size>> g_maxCliques;

// ============================================================
// Tuple utilities (shared with V2)
// ============================================================

using TupleKey = std::vector<daf::Size>;

struct TupleHash {
    size_t operator()(const TupleKey &t) const noexcept {
        size_t h = t.size();
        for (auto x : t) h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h << 6) + (h >> 2);
        return h;
    }
};

static void enumerateMultisets(
    const std::vector<daf::Size> &classes, int size, int startIdx,
    TupleKey &current, const std::function<void()> &callback)
{
    if ((int)current.size() == size) { callback(); return; }
    for (int i = startIdx; i < (int)classes.size(); ++i) {
        current.push_back(classes[i]);
        enumerateMultisets(classes, size, i, current, callback);
        current.pop_back();
    }
}

// Enumerate s-multisets that contain tau as a sub-multiset
// classes: available class IDs (sorted), tauCounts: minimum count per class from tau
// Fills remaining s-r slots from classes with repetition allowed up to classSizes
static void enumerateSupersetsOfTau(
    const std::vector<daf::Size> &classes,
    const std::vector<daf::Size> &classSizes,
    const std::unordered_map<daf::Size, int> &tauCounts,
    int s, int startIdx, TupleKey &current,
    const std::function<void(const TupleKey &)> &callback)
{
    if ((int)current.size() == s) { callback(current); return; }
    int remaining = s - (int)current.size();
    for (int i = startIdx; i < (int)classes.size(); ++i) {
        daf::Size c = classes[i];
        // How many of c already in current?
        int alreadyUsed = 0;
        for (auto x : current) if (x == c) alreadyUsed++;
        if (alreadyUsed >= (int)classSizes[c]) continue; // can't use more
        current.push_back(c);
        enumerateSupersetsOfTau(classes, classSizes, tauCounts, s, i, current, callback);
        current.pop_back();
    }
}

static std::vector<std::pair<daf::Size, int>> getComposition(const TupleKey &t) {
    std::vector<std::pair<daf::Size, int>> comp;
    for (auto c : t) {
        if (!comp.empty() && comp.back().first == c) comp.back().second++;
        else comp.push_back({c, 1});
    }
    return comp;
}

static void enumerateSubMultisets(
    const std::vector<std::pair<daf::Size, int>> &composition, int r,
    int classIdx, TupleKey &current,
    const std::function<void(const TupleKey &)> &callback)
{
    if (classIdx == (int)composition.size()) {
        if ((int)current.size() == r) callback(current);
        return;
    }
    auto [cls, maxCnt] = composition[classIdx];
    int remaining = r - (int)current.size();
    int maxFromLater = 0;
    for (int k = classIdx + 1; k < (int)composition.size(); ++k)
        maxFromLater += composition[k].second;
    int minJ = std::max(0, remaining - maxFromLater);
    int maxJ = std::min(maxCnt, remaining);

    for (int j = minJ; j <= maxJ; ++j) {
        for (int q = 0; q < j; ++q) current.push_back(cls);
        enumerateSubMultisets(composition, r, classIdx + 1, current, callback);
        for (int q = 0; q < j; ++q) current.pop_back();
    }
}

static double computeExt(const std::vector<std::pair<daf::Size, int>> &sigmaComp,
                          const TupleKey &tau,
                          const std::vector<daf::Size> &classSizes) {
    std::unordered_map<daf::Size, int> tauCounts;
    for (auto c : tau) tauCounts[c]++;
    double ext = 1.0;
    for (auto &[cls, mi] : sigmaComp) {
        int ji = tauCounts.count(cls) ? tauCounts[cls] : 0;
        int n = (int)classSizes[cls] - ji;
        int k = mi - ji;
        if (n < k || k < 0) return 0.0;
        ext *= nCr[n][k];
    }
    return ext;
}

// ============================================================
// Main function: Region CPI (V3)
// ============================================================

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionCPI_EventPeel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    auto tStart = std::chrono::high_resolution_clock::now();
    const daf::Size numVertices = edgeGraph.n;
    const daf::Size numPaths = tree.adj_list.size();
    const daf::Size INVALID = static_cast<daf::Size>(-1);

    // ============================================================
    // Step 1: Build Regions from g_maxCliques (MaxCliqEnum, pre-mutation)
    // ============================================================

    daf::Size validPaths = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid)
        if ((int)tree.adj_list[pid].size() >= s) validPaths++;

    daf::Size numRegions = 0;
    std::vector<std::vector<daf::Size>> regionVerts;
    std::vector<std::vector<daf::Size>> vtxMaxPaths(numVertices);

    for (auto &mc : g_maxCliques) {
        if ((int)mc.size() < s) continue;
        daf::Size rid = regionVerts.size();
        regionVerts.push_back(mc); // already sorted by MaxCliqEnum
        for (daf::Size v : mc)
            if (v < numVertices) vtxMaxPaths[v].push_back(rid);
        numRegions++;
    }

    // ============================================================
    // Step 1b: r-Mergeable Region Classification
    // ============================================================
    // If ALL vertices in region M have overlap < r with every other region:
    // → all r-cliques in M are only in M → support = C(|M|-r, s-r)
    // → directly assign core value, skip all pipeline work.

    auto tStep1a = std::chrono::high_resolution_clock::now();
    daf::Size maxVtxMCs = 0;
    double totalPairs = 0;
    {
        for (daf::Size v = 0; v < numVertices; ++v) {
            maxVtxMCs = std::max(maxVtxMCs, (daf::Size)vtxMaxPaths[v].size());
            daf::Size k = vtxMaxPaths[v].size();
            totalPairs += (double)k * (k - 1) / 2;
        }
        auto step1aMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - tStep1a).count();
        std::cout << "  Step 1: " << numRegions << " MCs (≥s), max vtxMCs=" << maxVtxMCs
                  << ", total pairs=" << std::fixed << std::setprecision(0) << totalPairs
                  << ", " << step1aMs << " ms" << std::endl;
    }
    daf::phaseMark("MCEnum");

    auto tStep1b = std::chrono::high_resolution_clock::now();

    // Fused fully-mergeable check (replaces the older two-phase
    // "precompute all pair intersections, then early-exit check" scheme).
    // For each region M, maintain an overlap counter array indexed by
    // other region id; increment at each v ∈ M for each other M' ∋ v;
    // break the moment any counter reaches r. Reset only the touched
    // entries via a dirty list. Same worst-case complexity, ~10–25×
    // faster in practice (array vs hashmap, plus early-exit on the
    // accumulation phase itself — the precompute never finishes pairs
    // that are irrelevant). See docs/region_cpi_theorems.md §7.
    std::vector<bool> fullyMergeable(numRegions, true);
    std::vector<int> overlapCnt(numRegions, 0);
    std::vector<daf::Size> dirtyList;
    dirtyList.reserve(256);

    long long fmTotalIter = 0;
    daf::Size fmExits = 0;
    double fmPairLimit = 5e7;
    if (const char *env = std::getenv("PIVOTER_V3_FM_PAIR_LIMIT")) {
        fmPairLimit = std::atof(env);
    }
    const bool forceFM = std::getenv("PIVOTER_V3_FM_FORCE") != nullptr;
    // Ablation: PIVOTER_DISABLE_RMERGE=1 force-disables the r-mergeable
    // optimization (Step 1b). All MCs stay in the active region; closed-form
    // shortcut skipped. Used to measure the merging's net contribution.
    const bool disableRMerge = std::getenv("PIVOTER_DISABLE_RMERGE") != nullptr;
    const bool skipFM = !forceFM && (disableRMerge || (fmPairLimit >= 0.0 && totalPairs > fmPairLimit));

    if (skipFM) {
        std::fill(fullyMergeable.begin(), fullyMergeable.end(), false);
    } else {
        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            bool fm = true;
            for (daf::Size v : regionVerts[rid]) {
                if (v >= numVertices) continue;
                for (daf::Size other : vtxMaxPaths[v]) {
                    ++fmTotalIter;
                    if (other == rid) continue;
                    if (overlapCnt[other] == 0) dirtyList.push_back(other);
                    if (++overlapCnt[other] >= (int)r) { fm = false; ++fmExits; break; }
                }
                if (!fm) break;
            }
            fullyMergeable[rid] = fm;
            for (auto d : dirtyList) overlapCnt[d] = 0;
            dirtyList.clear();
        }
    }

    {
        auto t = std::chrono::high_resolution_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t - tStep1b).count();
        if (disableRMerge) {
            std::cout << "  Step 1b fused-FM: DISABLED (PIVOTER_DISABLE_RMERGE=1, ablation), "
                      << ms << " ms" << std::endl;
        } else if (skipFM) {
            std::cout << "  Step 1b fused-FM: skipped (estimated pairs="
                      << std::fixed << std::setprecision(0) << totalPairs
                      << " > limit=" << fmPairLimit << "), "
                      << ms << " ms" << std::endl;
        } else {
            std::cout << "  Step 1b fused-FM: " << fmTotalIter << " counter ops, "
                      << fmExits << " early exits, " << ms << " ms" << std::endl;
        }
    }

    // Directly assign fully-mergeable regions
    std::map<double, int64_t> coreDist;
    daf::Size numFullyMergeable = 0;
    uint64_t mergedRCliques = 0;

    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        if (!fullyMergeable[rid]) continue;
        numFullyMergeable++;
        int mSize = (int)regionVerts[rid].size();
        double coreVal = (mSize >= (int)s) ? nCr[mSize - r][s - r] : 0.0;
        // Cast through uint64 — C(|M|, r) can exceed UINT32_MAX on large MCs
        // (e.g. com-dblp max MC=114 → C(114, 11) ≈ 6e14).
        uint64_t numRC = (uint64_t)llround(nCr[mSize][r]);
        coreDist[coreVal] += (int64_t)numRC;
        mergedRCliques += numRC;
    }

    // Rebuild regions/vtxMaxPaths with only NON-fully-mergeable regions
    {
        std::vector<std::vector<daf::Size>> newRegionVerts;
        std::vector<daf::Size> ridMap(numRegions, INVALID); // old rid → new rid
        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            if (fullyMergeable[rid]) continue;
            ridMap[rid] = newRegionVerts.size();
            newRegionVerts.push_back(std::move(regionVerts[rid]));
        }
        regionVerts = std::move(newRegionVerts);
        numRegions = regionVerts.size();

        // Rebuild vtxMaxPaths
        for (daf::Size v = 0; v < numVertices; ++v) {
            std::vector<daf::Size> newPaths;
            for (daf::Size rid : vtxMaxPaths[v]) {
                if (ridMap[rid] != INVALID)
                    newPaths.push_back(ridMap[rid]);
            }
            vtxMaxPaths[v] = std::move(newPaths);
        }
    }

    auto tStep1bEnd = std::chrono::high_resolution_clock::now();
    auto step1bMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep1bEnd - tStep1b).count();
    std::cout << "  r-Mergeable classification: " << step1bMs << " ms" << std::endl;
    std::cout << "    Fully mergeable regions: " << numFullyMergeable
              << " (" << mergedRCliques << " r-cliques, direct)" << std::endl;
    std::cout << "    Remaining regions: " << numRegions << std::endl;

    // -------- Observation-only probes (no pipeline change). See LowMemProbes.h.
    // PIVOTER_VSAFE_GAP_DIAG: B1 (globally-safe class) vs B2 (per-MC V_safe) gap.
    // PIVOTER_VSAFE_PROBE   : V_safe potential summary (safe%/strong%/incremental%).
    if (std::getenv("PIVOTER_VSAFE_GAP_DIAG") != nullptr)
        pivoter::lowmem_probes::probe_vsafe_gap_diag(
            r, numVertices, numRegions, regionVerts, vtxMaxPaths);
    if (std::getenv("PIVOTER_VSAFE_PROBE") != nullptr)
        pivoter::lowmem_probes::probe_vsafe(
            r, numVertices, numRegions, regionVerts, vtxMaxPaths);

    // Early exit if all regions handled
    if (numRegions == 0) {
        auto tEnd = std::chrono::high_resolution_clock::now();
        auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();
        double maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
        std::cout << "\n  All regions fully r-mergeable. No peeling needed." << std::endl;
        std::cout << "  Max core: " << maxCore << std::endl;
        for (auto &[core, cnt] : coreDist)
            std::cout << "  core=" << core << " count=" << cnt << std::endl;
        std::cout << "  Total time: " << totalMs << " ms" << std::endl;
        std::vector<std::pair<std::vector<daf::Size>, double>> result;
        for (auto &[c, cnt] : coreDist)
            // Compact format: key = {lo32, hi32} encoding int64_t count
        result.push_back({{(daf::Size)(cnt & 0xFFFFFFFF), (daf::Size)((cnt >> 32) & 0xFFFFFFFF)}, (double)c});
        return result;
    }

    // ============================================================
    // Step 2: Build Overlap Classes
    // ============================================================

    struct ClassInfo { std::vector<daf::Size> regionIds; daf::Size size; };
    std::vector<ClassInfo> classes;
    std::vector<daf::Size> classOf(numVertices, INVALID);
    {
        using Profile = std::vector<daf::Size>;
        struct PH {
            size_t operator()(const Profile &p) const noexcept {
                size_t h = p.size();
                for (auto x : p) h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h<<6) + (h>>2);
                return h;
            }
        };
        std::unordered_map<Profile, daf::Size, PH> pToC;
        for (daf::Size v = 0; v < numVertices; ++v) {
            if (vtxMaxPaths[v].empty()) continue;
            auto it = pToC.find(vtxMaxPaths[v]);
            if (it == pToC.end()) {
                daf::Size cid = classes.size();
                pToC[vtxMaxPaths[v]] = cid;
                classes.push_back({vtxMaxPaths[v], 1});
                classOf[v] = cid;
            } else { classOf[v] = it->second; classes[it->second].size++; }
        }
    }
    daf::Size numClasses = classes.size();

    std::vector<daf::Size> classSizes(numClasses);
    for (daf::Size i = 0; i < numClasses; ++i) classSizes[i] = classes[i].size;

    std::vector<std::vector<daf::Size>> classesInRegion(numRegions);
    for (daf::Size cid = 0; cid < numClasses; ++cid)
        for (daf::Size rid : classes[cid].regionIds)
            classesInRegion[rid].push_back(cid);
    for (auto &v : classesInRegion) std::sort(v.begin(), v.end());

    auto tStep2 = std::chrono::high_resolution_clock::now();
    std::cout << "======= Region CPI (V3) =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Vertices: " << numVertices << ", CPI paths: " << validPaths << std::endl;
    std::cout << "  Maximal cliques (≥s): " << numRegions << std::endl;
    std::cout << "  Overlap classes: " << numClasses << std::endl;

    // Private cloud mode:
    // - do not materialize private classes as active r-tuple dimensions
    // - count all r-cliques touching private vertices directly into coreDist
    // - keep private classes only as path-side support dimensions
    const bool enablePrivateCloud = !getenv("PIVOTER_V3_NO_PRIVATE");
    const bool enableDebugVerify = getenv("PIVOTER_V3_DEBUG_VERIFY") != nullptr;

    std::vector<bool> isPrivateClass(numClasses, false);
    std::vector<daf::Size> privateClassMC(numClasses, INVALID);
    for (daf::Size cid = 0; cid < numClasses; ++cid) {
        if (classes[cid].regionIds.size() == 1) {
            isPrivateClass[cid] = true;
            privateClassMC[cid] = classes[cid].regionIds[0];
        }
    }

    // Phase B1: globally-safe class detection (PIVOTER_VSAFE_CLOUD=1).
    // A class c is globally-safe iff:
    //   k=1 (private)  OR  k≥2 AND all pairs of c's MCs have overlap < r.
    // For multi-MC class with size |c| ≥ r, automatically NOT safe (since
    // |M_i ∩ M_j| ≥ |c| ≥ r). So we only check pairwise overlap for
    // classes with |c| < r — these are candidates for "soft private".
    const bool useVsafeCloud = std::getenv("PIVOTER_VSAFE_CLOUD") != nullptr;
    std::vector<bool> isGloballySafeClass(numClasses, false);
    long long vsafe_addedClasses = 0;
    long long vsafe_addedVertices = 0;
    if (useVsafeCloud) {
        // Helper: |M_i ∩ M_j| via sorted-list intersection on regionVerts.
        // regionVerts[rid] is sorted (set-like). O(|M_i| + |M_j|) per call.
        auto pairOverlap = [&](daf::Size r1, daf::Size r2) -> int {
            const auto &A = regionVerts[r1];
            const auto &B = regionVerts[r2];
            size_t i = 0, j = 0;
            int cnt = 0;
            while (i < A.size() && j < B.size()) {
                if (A[i] < B[j]) ++i;
                else if (A[i] > B[j]) ++j;
                else { ++cnt; ++i; ++j; }
            }
            return cnt;
        };

        for (daf::Size cid = 0; cid < numClasses; ++cid) {
            const auto &mcs = classes[cid].regionIds;
            if (mcs.size() == 1) {
                isGloballySafeClass[cid] = true;  // private case
                continue;
            }
            // |c| ≥ r ⇒ pairwise overlap ≥ r ⇒ not safe
            if ((int)classSizes[cid] >= (int)r) continue;
            // Check all pairs of MCs containing c
            bool allLow = true;
            for (size_t i = 0; i < mcs.size() && allLow; ++i) {
                for (size_t j = i + 1; j < mcs.size(); ++j) {
                    if (pairOverlap(mcs[i], mcs[j]) >= (int)r) {
                        allLow = false; break;
                    }
                }
            }
            if (allLow) {
                isGloballySafeClass[cid] = true;
                if (!isPrivateClass[cid]) {
                    ++vsafe_addedClasses;
                    vsafe_addedVertices += classSizes[cid];
                }
            }
        }
        std::cout << "  [V_safe Cloud] private classes: "
                  << std::count(isPrivateClass.begin(), isPrivateClass.end(), true)
                  << ", globally-safe additions: " << vsafe_addedClasses
                  << " (extra vertices: " << vsafe_addedVertices << ")"
                  << std::endl;
    }

    // Effective "private-like" check: original private OR (V_safe enabled and globally-safe)
    auto isVsafePrivateClass = [&](daf::Size cid) -> bool {
        return isPrivateClass[cid] || (useVsafeCloud && isGloballySafeClass[cid]);
    };

    auto isActiveTupleClass = [&](daf::Size cid) -> bool {
        return !enablePrivateCloud || !isVsafePrivateClass(cid);
    };

    std::vector<std::vector<daf::Size>> activeClassesInRegion(numRegions);
    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        for (daf::Size cid : classesInRegion[rid]) {
            if (isActiveTupleClass(cid))
                activeClassesInRegion[rid].push_back(cid);
        }
    }

    std::vector<double> mcCoreVal(numRegions, 0);
    std::vector<daf::Size> privateVertexCount(numRegions, 0);
    daf::Size numPrivateClouds = 0;
    double privateRCliquesDirect = 0.0;

    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        int m = (int)regionVerts[rid].size();
        int n = m - (int)r, k = (int)s - (int)r;
        mcCoreVal[rid] = (n >= k && k >= 0) ? nCr[n][k] : 0.0;
    }
    for (daf::Size cid = 0; cid < numClasses; ++cid) {
        // Phase B1: count contribution into EVERY MC the class is in.
        // For pure private classes (k=1), this matches the original behavior.
        // For multi-MC globally-safe classes (only when useVsafeCloud=1),
        // every MC the class is in receives the |class| contribution because
        // (a) every vertex of the class is in every such MC, and
        // (b) by touching-safe theorem, an r-clique touching such a vertex
        //     in MC M_i is uniquely contained in M_i (not shared with M_j).
        // So summing over MCs counts each unique r-clique once.
        if (!isVsafePrivateClass(cid)) continue;
        for (auto rid : classes[cid].regionIds) {
            privateVertexCount[rid] += classSizes[cid];
        }
        numPrivateClouds++;
    }
    if (enablePrivateCloud) {
        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            int m = (int)regionVerts[rid].size();
            int p = (int)privateVertexCount[rid];
            if (p == 0 || m < (int)r) continue;
            int q = m - p;
            double total = nCr[m][r];
            double publicOnly = (q >= (int)r) ? nCr[q][r] : 0.0;
            // Widen to uint64: C(|M|, r) can exceed UINT32_MAX for large MCs.
            uint64_t privateTouching = (uint64_t)llround(std::max(0.0, total - publicOnly));
            if (privateTouching == 0) continue;
            coreDist[mcCoreVal[rid]] += (int64_t)privateTouching;
            privateRCliquesDirect += (double)privateTouching;
        }
    }

    // ============================================================
    // Step 3: Enumerate active r-tuples
    // ============================================================

    robin_hood::unordered_flat_map<TupleKey, daf::Size, TupleHash> rTupleIndex;
    // PERF-AUDIT MEMORY (7): tuple keys live in a single flat buffer
    // `tupleKeysFlat` of length |rTuples| × r, sliced as
    //     keyOf(tidx) = tupleKeysFlat[tidx*r .. tidx*r + r)
    // instead of giving every RTuple its own std::vector. A
    // std::vector<daf::Size> with size r occupies 24 (header) + 4r
    // bytes per tuple; the flat layout drops the 24-byte header
    // entirely. On graphs with N tuples this saves 24·N bytes — e.g.
    // ~240 MB on a 10 M-tuple instance — with no observable peel-loop
    // cost (the merge-scan we use to read the key is byte-identical
    // either way; std::span<const daf::Size> replaces the const-ref).
    //
    // mult = Π_c C(|c|, j_c) can overflow uint32 for r≥5 on graphs with
    // large class sizes (e.g. r=11 on com-dblp: C(~50, 11) ≈ 4e10).
    // Widen to uint64 so the multiplicity is represented correctly
    // throughout peeling and core aggregation.
    struct RTuple { uint64_t mult; };
    std::vector<RTuple> rTuples;
    std::vector<daf::Size> tupleKeysFlat;
    std::vector<double> tupleMinCore; // per-tuple minCore floor
    const int rPerTuple = (int)r;
    auto keyOf = [&](daf::Size tidx) -> std::span<const daf::Size> {
        return std::span<const daf::Size>(
            tupleKeysFlat.data() + (size_t)tidx * (size_t)rPerTuple,
            (size_t)rPerTuple);
    };

    {
        TupleKey cur; cur.reserve(r);
        daf::Size curRid = 0;
        auto addRTuple = [&](const TupleKey &key) {
            std::unordered_map<daf::Size, int> counts;
            for (auto c : key) counts[c]++;
            uint64_t mult = 1;
            for (auto &[c, k] : counts) {
                if ((int)classSizes[c] < k) return;
                mult *= (uint64_t)llround(nCr[classSizes[c]][k]);
            }
            if (mult == 0) return;
            auto it = rTupleIndex.find(key);
            if (it == rTupleIndex.end()) {
                rTupleIndex[key] = rTuples.size();
                rTuples.push_back({mult});
                tupleKeysFlat.insert(tupleKeysFlat.end(), key.begin(), key.end());
                tupleMinCore.push_back(enablePrivateCloud ? mcCoreVal[curRid] : 0);
            } else if (enablePrivateCloud) {
                tupleMinCore[it->second] = std::max(tupleMinCore[it->second], mcCoreVal[curRid]);
            }
        };

        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            const auto &cids = enablePrivateCloud ? activeClassesInRegion[rid] : classesInRegion[rid];
            if (cids.size() > 500) {
                std::cerr << "ABORT: MC " << rid << " has " << cids.size()
                          << " classes (>500), too large for tuple enumeration" << std::endl;
                exit(1);
            }
            curRid = rid;
            cur.clear();
            std::function<void()> cb = [&]() { addRTuple(cur); };
            enumerateMultisets(cids, r, 0, cur, cb);
        }
    }

    auto tStep3 = std::chrono::high_resolution_clock::now();
    auto step2Ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        tStep3 - tStep2).count();
    std::cout << "  Step 2+3 time (classes + tuples): " << step2Ms << " ms" << std::endl;
    daf::phaseMark("Index");
    std::cout << "  Active r-tuples: " << rTuples.size() << std::endl;
    if (enablePrivateCloud) {
        std::cout << "  Private clouds: " << numPrivateClouds << std::endl;
        std::cout << "  Private-touching r-cliques (direct): "
                  << std::fixed << std::setprecision(0) << privateRCliquesDirect << std::endl;
    } else {
        std::cout << "  Private cloud mode: disabled" << std::endl;
    }

    // ============================================================
    // Step 4: CPI Counting for initial support (THE KEY INNOVATION)
    // ============================================================
    // CORRECT formula (handles covered-but-not-encoded r-cliques):
    //
    // AggrCount(τ, P) = Σ_{b_R} [Π_R C(nh_R, j_R-b_R) C(np_R, b_R)] × C(p-Σb_R, s-h-Σb_R)
    //
    // Computed via convolution: g_c(x) per class, f = Π g_c, then sum f[t]×C(p-t, s-h-t).

    std::vector<double> support(rTuples.size(), 0.0);

    daf::Size pathsUsed = 0;
    daf::Size totalTuplesOnPaths = 0;

    // Dense per-class counters (replace the per-path unordered_map<cid, (nh, np)>).
    // Maintained across paths via a dirty list; reset at the top of each iteration.
    // See docs/region_cpi_theorems.md §7 for the same pattern used in Step 1b.
    std::vector<int> nhArr(numClasses, 0);
    std::vector<int> npArr(numClasses, 0);
    std::vector<daf::Size> classDirty;
    classDirty.reserve(64);

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        // Reset dirty counters from previous iteration (if any)
        for (auto cid : classDirty) { nhArr[cid] = 0; npArr[cid] = 0; }
        classDirty.clear();

        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        // Compute h, p and per-class (nh, np) distribution on this path
        int h = 0, p = 0;
        bool pathHasPrivateHold = false;
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            if (enablePrivateCloud && isVsafePrivateClass(cid)) {
                // Private hold → all s-cliques on this path use a private vertex
                if (!node.isPivot) pathHasPrivateHold = true;
                continue; // exclude private vertices from h/p counts
            }
            if (nhArr[cid] == 0 && npArr[cid] == 0) classDirty.push_back(cid);
            if (node.isPivot) { p++; npArr[cid]++; }
            else { h++; nhArr[cid]++; }
        }

        // Skip path if it has a private hold (all s-cliques are private-touching)
        if (pathHasPrivateHold) continue;
        if (h + p < (int)s) continue;
        pathsUsed++;

        // classDirty holds the unique classes on this path; sort for canonical
        // enumeration order (TupleKey is the sorted multiset of class ids).
        std::sort(classDirty.begin(), classDirty.end());

        // Enumerate r-multisets of this path's classes
        TupleKey cur; cur.reserve(r);
        std::vector<std::pair<daf::Size, int>> tauClasses;
        tauClasses.reserve(r);
        std::vector<double> fBuf;
        std::vector<double> nextBuf;
        std::vector<double> gcBuf;
        fBuf.reserve((size_t)r + 1);
        nextBuf.reserve((size_t)r + 1);
        gcBuf.reserve((size_t)r + 1);
        std::function<void()> cb = [&]() {
            // Build composition: (class, count) pairs
            tauClasses.clear();
            for (size_t k = 0; k < cur.size(); ) {
                daf::Size c = cur[k];
                int jc = 1;
                while (k + jc < cur.size() && cur[k + jc] == c) ++jc;
                tauClasses.push_back({c, jc});
                k += jc;
            }

            // Feasibility: each tuple class c is enumerated from classDirty, so
            // nhArr[c] + npArr[c] > 0 is already guaranteed. We still need the
            // jc <= nh+np check.
            for (auto &[c, jc] : tauClasses) {
                if (jc > nhArr[c] + npArr[c]) return;
            }

            // Look up this r-tuple
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;

            // Convolution to compute AggrCount
            // f[t] = coefficient of x^t in Π_c g_c(x)
            // where g_c(x) = Σ_{b_c} C(nh_c, j_c-b_c) × C(np_c, b_c) × x^{b_c}

            if (fBuf.empty()) fBuf.resize(1);
            fBuf[0] = 1.0;
            int fLen = 1;

            for (auto &[c, jc] : tauClasses) {
                int nhc = nhArr[c], npc = npArr[c];
                int bMin = std::max(0, jc - nhc);
                int bMax = std::min(jc, npc);
                if (bMin > bMax) return; // infeasible

                // g_c coefficients
                const int gcLen = bMax + 1;
                if ((int)gcBuf.size() < gcLen) gcBuf.resize(gcLen);
                std::fill(gcBuf.begin(), gcBuf.begin() + gcLen, 0.0);
                for (int bc = bMin; bc <= bMax; ++bc)
                    gcBuf[bc] = nCr[nhc][jc - bc] * nCr[npc][bc];

                // Convolve f with gc
                const int nextLen = fLen + gcLen - 1;
                if ((int)nextBuf.size() < nextLen) nextBuf.resize(nextLen);
                std::fill(nextBuf.begin(), nextBuf.begin() + nextLen, 0.0);
                for (int i = 0; i < fLen; ++i)
                    for (int j = 0; j < gcLen; ++j)
                        nextBuf[i + j] += fBuf[i] * gcBuf[j];
                if ((int)fBuf.size() < nextLen) fBuf.resize(nextLen);
                std::copy(nextBuf.begin(), nextBuf.begin() + nextLen, fBuf.begin());
                fLen = nextLen;
            }

            // AggrCount = Σ_t f[t] × C(p-t, s-h-t)
            double aggr = 0.0;
            for (int t = 0; t < fLen; ++t) {
                if (fBuf[t] == 0.0) continue;
                int n = p - t, k = s - h - t;
                if (n >= 0 && k >= 0 && n >= k)
                    aggr += fBuf[t] * nCr[n][k];
            }

            if (aggr > 0) {
                support[tidx] += aggr / rTuples[tidx].mult;
                totalTuplesOnPaths++;
            }
        };
        enumerateMultisets(classDirty, r, 0, cur, cb);
    }
    // Final cleanup of dirty counters so the scratch is clean for any later user
    for (auto cid : classDirty) { nhArr[cid] = 0; npArr[cid] = 0; }
    classDirty.clear();

    auto tStep4 = std::chrono::high_resolution_clock::now();
    auto step4Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep4 - tStep3).count();

    double totalSupportTuples = 0, totalRCliques = 0;
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        totalSupportTuples += rTuples[i].mult * support[i];
        totalRCliques += rTuples[i].mult;
    }
    double totalRCliquesWithPrivate = totalRCliques + privateRCliquesDirect;
    std::cout << "  CPI paths used: " << pathsUsed << std::endl;
    std::cout << "  Tuple-path pairs: " << totalTuplesOnPaths << std::endl;
    std::cout << "  Active r-cliques: " << std::fixed << std::setprecision(0) << totalRCliques << std::endl;
    if (enablePrivateCloud) {
        std::cout << "  Total r-cliques (with private cloud): "
                  << std::fixed << std::setprecision(0) << totalRCliquesWithPrivate << std::endl;
    }
    std::cout << "  Support sum (CPI): " << totalSupportTuples << std::endl;
    std::cout << "  CPI counting time: " << step4Ms << " ms" << std::endl;
    daf::phaseMark("Support");

    // ============================================================
    // Step 5+6: Constrained Path Peeling (Analytical Split)
    // ============================================================
    // Constrained Path = CPI path + per-class (min_piv, max_piv) bounds.
    // When τ is peeled on P̂: subtract old contributions, split P̂ into
    // κ disjoint sub-paths (each a ConstrainedPath), add new contributions.
    // NO s-tuple enumeration. NO inclusion-exclusion. NO BK execution.

    auto tStep6 = std::chrono::high_resolution_clock::now();

    // --- PathInfo: immutable per-path data + mutable dead boxes ---
    struct PathInfo {
        int h, T;                    // holds count, target = s - h
        std::vector<daf::Size> classIds;   // ordered class IDs on this path
        std::vector<int> nh;         // per-class hold count (parallel to classIds)
        std::vector<int> np;         // per-class pivot count (parallel to classIds)
        std::vector<daf::Size> tupleIdxs;  // alive tuple slots on this path
        // EventPeel per-slot state (parallel to tupleIdxs):
        //   liveContrib[slot] = remaining live contribution of P to the
        //                       tuple, in RAW incidence units (no /mult)
        //   tupleBase[slot]   = base factor Π_k C(nh_k, j_k); 0 exactly for
        //                       pivot-needy tuples
        std::vector<double> liveContrib;
        std::vector<double> tupleBase;
        // Lazily materialized cell (witness-orbit) set. A cell is one
        // pivot-count vector y (Σ y_k = T, 0 ≤ y_k ≤ np_k), stored as
        // packed sparse coords (k << 12 | y) in cellCoords with offsets
        // cellStart. cellsByClass[k] lists cell ids with y_k > 0.
        bool cellsBuilt = false;
        daf::Size aliveCells = 0;
        uint32_t evEpoch = 0;          // per-tuple dedupe epoch (multi-touch)
        std::vector<uint32_t> cellCoords;
        std::vector<uint32_t> cellStart;
        std::vector<uint8_t> cellAlive;
        std::vector<uint32_t> cellEpoch;   // last epoch each cell was visited
        std::vector<std::vector<uint32_t>> cellsByClass;
    };

    // LowMem: classIds is sorted (see build below), so binary search is O(log k)
    // with k = |classes on path|, typ. small.  Saves ~48-56 B/entry vs unordered_map.
    auto findClassIdx = [](const PathInfo &pi, daf::Size cid) -> int {
        auto it = std::lower_bound(pi.classIds.begin(), pi.classIds.end(), cid);
        if (it == pi.classIds.end() || *it != cid) return -1;
        return (int)(it - pi.classIds.begin());
    };

    // Compile-time bounds for stack-allocated hot-path buffers. m (classes per
    // path) is small (typ. 10-30). T (pivot target) and per-class range can
    // grow with s — for large s on graphs with big MCs, range can hit 100+.
    // Bumped bounds cover ca-HepPh (max MC 239), web-it-2004 (max MC 432)
    // and similar very dense graphs. If a path exceeds these, we fall back
    // to heap (via _heap helpers below).
    constexpr int MAX_M = 512;
    constexpr int MAX_T = 256;

    // --- Weighted feasible count: convolution DP (stack buffers, no allocation) ---
    // Caller builds a packed weight table: wtsFlat[i*wtStride + k] = weight for
    // class i at shifted allocation k (k = b_i - lower[i]), valid for 0 ≤ k < wtLen[i].
    auto countFeasibleWeighted_sb = [](const int *lower, const int *upper, int m, int T,
                                       const double *wtsFlat, int wtStride, const int *wtLen) -> double {
        int minSum = 0, maxSum = 0;
        for (int i = 0; i < m; ++i) {
            if (lower[i] > upper[i]) return 0.0;
            minSum += lower[i];
            maxSum += upper[i];
        }
        if (T < minSum || T > maxSum) return 0.0;

        int target = T - minSum;
        if (target < 0) return 0.0;

        // Fast path: m=1 — one class, k_0 = target uniquely.
        // wtLen[0] = upper[0] - lower[0] + 1, target ≤ this by the
        // T ≤ maxSum check, so the index is in range.
        if (m == 1) {
            return wtsFlat[target];
        }

        // Fast path: m=2 — single convolution loop, no DP buffers.
        if (m == 2) {
            const double *wt0 = wtsFlat;
            const double *wt1 = wtsFlat + wtStride;
            int len0 = wtLen[0];
            int len1 = wtLen[1];
            int kLo = std::max(0, target - (len1 - 1));
            int kHi = std::min(target, len0 - 1);
            double sum = 0.0;
            for (int k = kLo; k <= kHi; ++k) {
                sum += wt0[k] * wt1[target - k];
            }
            return sum;
        }

        // LowMem: dp/next buffers on heap via thread_local scratch so MC=432
        // graphs (e.g. web-it-2004) no longer abort when target > MAX_T.
        // Stack buffers were the bottleneck on MC-dominated graphs.
        thread_local std::vector<double> dpBuf;
        if ((int)dpBuf.size() < 2 * (target + 1)) dpBuf.resize(2 * (target + 1));
        double *dp = dpBuf.data();
        double *next = dp + (target + 1);
        for (int t = 0; t <= target; ++t) dp[t] = 0.0;
        dp[0] = 1.0;
        for (int i = 0; i < m; ++i) {
            int cap = upper[i] - lower[i];
            const double *wt = wtsFlat + (size_t)i * wtStride;
            int len = wtLen[i];
            int K = cap < len - 1 ? cap : len - 1;  // common kMax for the bulk

            // Initial part: t < K → kMax = t (varying length per t).
            int initEnd = K < target + 1 ? K : target + 1;
            for (int t = 0; t < initEnd; ++t) {
                double sum = 0.0;
                for (int k = 0; k <= t; ++k)
                    sum += dp[t - k] * wt[k];
                next[t] = sum;
            }

            // Bulk part: t ∈ [K, target], kMax = K constant. Process two
            // adjacent t values per iteration via explicit SIMD. The key
            // observation: a single 2-double load at &dp[t-k] yields
            // [dp[t-k], dp[t-k+1]], where lane 0 is what t needs and
            // lane 1 is what t+1 needs — natural alignment with the
            // convolution structure.
            int t = K;
#if defined(REGNDC_HAVE_NEON)
            // 2-way NEON: a single 2-double load at &dp[t-k] yields
            // [dp[t-k], dp[t-k+1]] = exactly what t and t+1 each need
            // for their k-th convolution term. Fused multiply-add
            // accumulates both lanes in one cycle.
            for (; t + 1 <= target; t += 2) {
                float64x2_t s = vdupq_n_f64(0.0);
                for (int k = 0; k <= K; ++k) {
                    float64x2_t d = vld1q_f64(&dp[t - k]);
                    float64x2_t w = vdupq_n_f64(wt[k]);
                    s = vfmaq_f64(s, d, w);
                }
                vst1q_f64(&next[t], s);
            }
#elif defined(REGNDC_HAVE_SSE2)
            for (; t + 1 <= target; t += 2) {
                __m128d s = _mm_setzero_pd();
                for (int k = 0; k <= K; ++k) {
                    __m128d d = _mm_loadu_pd(&dp[t - k]);
                    __m128d w = _mm_set1_pd(wt[k]);
                    s = _mm_add_pd(s, _mm_mul_pd(d, w));
                }
                _mm_storeu_pd(&next[t], s);
            }
#else
            // Scalar dual-accumulator fallback.
            for (; t + 1 <= target; t += 2) {
                double s0 = 0.0, s1 = 0.0;
                for (int k = 0; k <= K; ++k) {
                    double w = wt[k];
                    s0 += dp[t - k] * w;
                    s1 += dp[t + 1 - k] * w;
                }
                next[t] = s0;
                next[t + 1] = s1;
            }
#endif
            // Tail: at most one t left.
            for (; t <= target; ++t) {
                double sum = 0.0;
                for (int k = 0; k <= K; ++k)
                    sum += dp[t - k] * wt[k];
                next[t] = sum;
            }

            std::swap_ranges(dp, dp + target + 1, next);
        }
        return dp[target];
    };

    // --- Build weight tables for tuple τ' on path P ---
    // Writes into wtsFlat (caller-sized to m * wtStride). wtLen[i] = range per class.
    // Uses merge-scan on sorted key + sorted pi.classIds to avoid an unordered_map.
    // LowMem: jvec moved from stack [MAX_M] to thread_local heap so paths with
    // m > MAX_M (e.g. web-it-2004 MC=432) no longer overflow.
    auto buildWeightTables_sb = [&](daf::Size tidx, const PathInfo &pi,
                                    const int *lower, const int *upper,
                                    double *wtsFlat, int wtStride, int *wtLen) {
        int m = (int)pi.classIds.size();
        auto key = keyOf(tidx);   // sorted with repetitions
        thread_local std::vector<int> jvecBuf;
        if ((int)jvecBuf.size() < m) jvecBuf.resize(m);
        int *jvec = jvecBuf.data();
        for (int i = 0; i < m; ++i) jvec[i] = 0;
        {
            int ci = 0, ki = 0;
            while (ci < m && ki < (int)key.size()) {
                if (pi.classIds[ci] < key[ki]) { ++ci; }
                else if (pi.classIds[ci] > key[ki]) { ++ki; }
                else {
                    int jc = 1;
                    while (ki + jc < (int)key.size() && key[ki + jc] == key[ki]) ++jc;
                    jvec[ci] = jc;
                    ++ci; ki += jc;
                }
            }
        }
        for (int i = 0; i < m; ++i) {
            int lo = lower[i], hi = upper[i];
            int range = hi - lo + 1;
            double *row = wtsFlat + (size_t)i * wtStride;
            if (range <= 0) { wtLen[i] = 0; continue; }
            // LowMem: range > wtStride used to abort; callers now size
            // wtStride dynamically (max range across classes) so this
            // should never trip, but assert for safety during development.
            if (range > wtStride) { std::cerr << "buildWeightTables: range=" << range
                << " > wtStride=" << wtStride << " (bug: caller under-sized buffer)\n";
                std::abort(); }
            wtLen[i] = range;
            int jc = jvec[i];
            int nhc = pi.nh[i], npc = pi.np[i];
            for (int k = 0; k < range; ++k) {
                int tc = lo + k;
                double v = 0.0;
                if (tc >= 0 && tc <= npc) {
                    if (jc > 0) {
                        int pool = nhc + tc;
                        if (pool >= jc) v = nCr[npc][tc] * nCr[pool][jc];
                    } else {
                        v = nCr[npc][tc];
                    }
                }
                row[k] = v;
            }
        }
    };

    // --- Build PathInfo structures ---
    // LowMem-adaptive: build BOTH candidate indices during the PathInfo
    // loop, then at build-end keep only the cheaper one.
    //   * classToPaths    : Σ_p |pi.classIds|  (wins on compressible graphs
    //                                           — ca-HepPh, com-dblp)
    //   * tupleToPathInfos: Σ_p |pi.tupleIdxs| (wins on mega-clique hubs
    //                                           — web-it-2004 MC=432, where
    //                                           a few classes dominate
    //                                           thousands of paths)
    // Without this adaptive choice V3LM OOMed on web-it-2004 r=3 s=125+
    // (V3 works with 188 MB; classToPaths intersection explodes to 250 GB).
    std::vector<PathInfo> pathInfos;
    std::vector<std::vector<daf::Size>> classToPaths(numClasses);
    std::vector<std::vector<daf::Size>> tupleToPathInfos(rTuples.size());

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        for (auto cid : classDirty) { nhArr[cid] = 0; npArr[cid] = 0; }
        classDirty.clear();

        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        PathInfo pi;
        pi.h = 0;
        bool piHasPrivateHold = false;
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            if (enablePrivateCloud && isVsafePrivateClass(cid)) {
                if (!node.isPivot) piHasPrivateHold = true;
                continue;
            }
            if (nhArr[cid] == 0 && npArr[cid] == 0) classDirty.push_back(cid);
            if (node.isPivot) npArr[cid]++;
            else { nhArr[cid]++; pi.h++; }
        }
        if (piHasPrivateHold) continue;
        pi.T = s - pi.h;
        if (pi.T < 0) continue;

        // Build ordered class arrays
        pi.classIds.assign(classDirty.begin(), classDirty.end());
        std::sort(pi.classIds.begin(), pi.classIds.end());
        pi.nh.resize(pi.classIds.size());
        pi.np.resize(pi.classIds.size());
        for (int i = 0; i < (int)pi.classIds.size(); ++i) {
            pi.nh[i] = nhArr[pi.classIds[i]];
            pi.np[i] = npArr[pi.classIds[i]];
        }

        // Find tuples on this path
        TupleKey cur; cur.reserve(r);
        const std::vector<daf::Size> &pathClasses = pi.classIds;
        bool hasTuples = false;
        daf::Size piIdx = pathInfos.size();

        std::function<void()> enumCb = [&]() {
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;
            // Feasibility check
            for (size_t k = 0; k < cur.size(); ) {
                daf::Size c = cur[k];
                int jc = 1;
                while (k + jc < cur.size() && cur[k + jc] == c) ++jc;
                int idx2 = findClassIdx(pi, c);
                if (idx2 < 0) return;
                if (jc > pi.nh[idx2] + pi.np[idx2]) return;
                k += jc;
            }
            pi.tupleIdxs.push_back(tidx);
            hasTuples = true;
        };
        enumerateMultisets(pathClasses, r, 0, cur, enumCb);

        if (hasTuples) {
            pathInfos.push_back(std::move(pi));
            // Populate BOTH candidate indices. piIdx is strictly ascending,
            // so both lists are naturally sorted by path id.
            for (auto cid : pathInfos.back().classIds)
                classToPaths[cid].push_back(piIdx);
            for (auto tidx : pathInfos.back().tupleIdxs)
                tupleToPathInfos[tidx].push_back(piIdx);
        }
    }
    for (auto cid : classDirty) { nhArr[cid] = 0; npArr[cid] = 0; }
    classDirty.clear();

    // PIVOTER_HOMCC_PROBE: tuple-level High-Overlap-MC-Cluster classification.
    // See LowMemProbes.h.
    if (std::getenv("PIVOTER_HOMCC_PROBE") != nullptr)
        pivoter::lowmem_probes::probe_homcc_tuple(
            r, regionVerts, classes, rTuples, keyOf);

    // Decide which index is cheaper and free the other.
    size_t cost_class = 0, cost_tuple = 0;
    for (auto &v : classToPaths)     cost_class += v.size();
    for (auto &v : tupleToPathInfos) cost_tuple += v.size();
    bool use_class_index = (cost_class <= cost_tuple);
    std::cout << "  classToPaths entries: " << cost_class
              << "  tupleToPathInfos entries: " << cost_tuple
              << "  -> using " << (use_class_index ? "classToPaths" : "tupleToPathInfos")
              << std::endl;
    if (use_class_index) {
        std::vector<std::vector<daf::Size>>().swap(tupleToPathInfos);
    } else {
        std::vector<std::vector<daf::Size>>().swap(classToPaths);
    }

    auto tPathBuild = std::chrono::high_resolution_clock::now();
    auto pathBuildMs = std::chrono::duration_cast<std::chrono::milliseconds>(tPathBuild - tStep6).count();
    std::cout << "  PathInfo count: " << pathInfos.size() << std::endl;
    std::cout << "  PathInfo build time: " << pathBuildMs << " ms" << std::endl;

    // --- EventPeel init: per-(path, tuple) live contribution + base ---
    // liveContrib[slot] = Σ_{cells y on P} W_τ(y) in RAW incidence units,
    // computed by the same weighted convolution DP that Step 4 uses, so
    // Σ_P liveContrib/mult must reproduce the Step-4 support (checked).
    // Slots that start at zero contribution are dropped (support only
    // decreases, so they can never matter again).
    {
        auto tInit0 = std::chrono::high_resolution_clock::now();
        std::vector<double> supCheck(rTuples.size(), 0.0);
        std::vector<int> baseV, upperV;
        std::vector<double> wtsFlatV;
        std::vector<int> wtLenV;
        for (daf::Size piIdx = 0; piIdx < pathInfos.size(); ++piIdx) {
            auto &pi = pathInfos[piIdx];
            const int m = (int)pi.classIds.size();
            pi.liveContrib.resize(pi.tupleIdxs.size());
            pi.tupleBase.resize(pi.tupleIdxs.size());
            baseV.assign(m, 0); upperV.assign(m, 0);
            size_t w = 0;
            for (size_t slot = 0; slot < pi.tupleIdxs.size(); ++slot) {
                daf::Size tidx = pi.tupleIdxs[slot];
                auto key = keyOf(tidx);
                double tb = 1.0;
                for (int i = 0; i < m; ++i) { baseV[i] = 0; upperV[i] = pi.np[i]; }
                for (size_t kk = 0; kk < key.size(); ) {
                    daf::Size c = key[kk]; int jc = 1;
                    while (kk + jc < key.size() && key[kk + jc] == c) ++jc;
                    int k = findClassIdx(pi, c);
                    baseV[k] = std::max(0, jc - pi.nh[k]);
                    tb *= (jc <= pi.nh[k]) ? nCr[pi.nh[k]][jc] : 0.0;
                    kk += jc;
                }
                int wtStride = 1;
                for (int i = 0; i < m; ++i)
                    wtStride = std::max(wtStride, upperV[i] - baseV[i] + 1);
                wtsFlatV.resize((size_t)m * wtStride);
                wtLenV.resize(m);
                buildWeightTables_sb(tidx, pi, baseV.data(), upperV.data(),
                                     wtsFlatV.data(), wtStride, wtLenV.data());
                double aggr = countFeasibleWeighted_sb(
                    baseV.data(), upperV.data(), m, pi.T,
                    wtsFlatV.data(), wtStride, wtLenV.data());
                supCheck[tidx] += aggr / rTuples[tidx].mult;
                if (aggr < 0.5) continue;
                pi.tupleIdxs[w] = tidx;
                pi.liveContrib[w] = aggr;
                pi.tupleBase[w] = tb;
                ++w;
            }
            pi.tupleIdxs.resize(w);
            pi.liveContrib.resize(w);
            pi.tupleBase.resize(w);
        }
        int initMismatch = 0;
        for (daf::Size i = 0; i < rTuples.size(); ++i)
            if (std::abs(supCheck[i] - support[i]) > 0.5) initMismatch++;
        auto tInit1 = std::chrono::high_resolution_clock::now();
        std::cout << "  EventPeel init: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(tInit1 - tInit0).count()
                  << " ms, init check "
                  << (initMismatch == 0 ? "PASS" : ("MISMATCH(" + std::to_string(initMismatch) + ")"))
                  << std::endl;
    }

    // --- LowMem: free Step-1..5 scaffolding that the peel loop never reads ---
    // Everything below this comment is referenced only during PathInfo build
    // and earlier phases; the peel loop uses only rTuples, pathInfos,
    // classToPaths, support/dSup, and a handful of scalar per-tuple arrays.
    // Swap-to-empty instead of clear() to actually release the heap buffer
    // (clear keeps capacity).
    {
        robin_hood::unordered_flat_map<TupleKey, daf::Size, TupleHash>().swap(rTupleIndex);
        std::vector<std::vector<daf::Size>>().swap(regionVerts);
        std::vector<std::vector<daf::Size>>().swap(vtxMaxPaths);
        std::vector<ClassInfo>().swap(classes);
        std::vector<daf::Size>().swap(classOf);
        std::vector<daf::Size>().swap(classSizes);
        std::vector<std::vector<daf::Size>>().swap(classesInRegion);
        std::vector<bool>().swap(isPrivateClass);
        std::vector<daf::Size>().swap(privateClassMC);
        std::vector<std::vector<daf::Size>>().swap(activeClassesInRegion);
        std::vector<double>().swap(mcCoreVal);
        std::vector<daf::Size>().swap(privateVertexCount);
        // cpaths / tupleToCPaths live inside the `if (enableDebugVerify)`
        // block (lines 619-791) and are RAII-freed when that scope closes —
        // nothing to do here for production runs.
    }

    // --- LowMem: derive tuple->paths on the fly via class->path intersection ---
    //
    // For a tuple tau with distinct key classes {c_1, ..., c_k}, the paths
    // covering tau are exactly the intersection of classToPaths[c_i].  We
    // also filter out paths that lack sufficient vertices of a class to
    // realize the tuple's multiplicities (same check the original build
    // loop does).  Cost: O(sum |classToPaths[c_i]|) for the sorted-list
    // intersection plus O(|candidates|) for the multiplicity filter.
    //
    // Hot-path audit fixes (compared to the naive implementation):
    //  * `key` is already sorted with repetitions → build (class, mult)
    //    pairs in one O(k) scan, avoiding std::map allocations.
    //  * All per-call vectors (counts/cur/next/out) are captured scratch
    //    reused across calls — zero mallocs on the hot path after warmup.
    //  * Inside the filter loop the `idx < 0` branch is dead (intersection
    //    guarantees class-in-pi.classIds), removed.
    std::vector<std::pair<daf::Size, int>> pct_counts_buf;
    std::vector<daf::Size> pct_cur_buf, pct_next_buf;
    // Retire flag: set when path's alive tuple count drops to 0. Once
    // retired, the path's metadata (classIds/nh/np) is freed to reclaim
    // memory; pathsCoveringTuple's multiplicity filter skips retired
    // paths via this flag instead of dereferencing freed nh/np.
    // Sized after pathInfos build is complete (just before this lambda
    // is FIRST invoked, which happens in the peel loop below).
    std::vector<bool> pathRetired;
    // Forward-declared because pathsCoveringTuple (lambda capture by [&])
    // reads it. Sized at the same point as pathRetired below.
    std::vector<daf::Size> pathAliveCount;
    auto pathsCoveringTuple = [&](daf::Size tidx,
                                  std::vector<daf::Size> &out) {
        out.clear();
        // Fast path (adaptive): when the tuple→path direct index was
        // selected at build time, just copy that list — no intersection,
        // no multiplicity filter (tupleToPathInfos already stores only
        // paths on which tidx is feasible, by construction during build).
        // Retired paths are pre-filtered to avoid downstream readers
        // touching freed metadata.
        if (!use_class_index) {
            const auto &src = tupleToPathInfos[tidx];
            out.reserve(src.size());
            for (auto p : src)
                if (!pathRetired[p]) out.push_back(p);
            return;
        }
        const auto key = keyOf(tidx);
        if (key.empty()) return;

        // Build (class, mult) pairs from sorted key with no heap allocs.
        auto &counts = pct_counts_buf;
        counts.clear();
        for (size_t k = 0; k < key.size(); ) {
            daf::Size c = key[k];
            int jc = 1;
            while (k + jc < key.size() && key[k + jc] == c) ++jc;
            counts.emplace_back(c, jc);
            k += jc;
        }

        // Find smallest list; abort early if any class is absent.
        int smallestIdx = -1;
        size_t smallest = SIZE_MAX;
        for (int i = 0; i < (int)counts.size(); ++i) {
            daf::Size c = counts[i].first;
            if (c >= classToPaths.size() || classToPaths[c].empty()) return;
            if (classToPaths[c].size() < smallest) {
                smallest = classToPaths[c].size();
                smallestIdx = i;
            }
        }

        // Seed 'cur' with the smallest list (fewest iterations for
        // subsequent intersections). assign() reuses capacity.
        auto &cur = pct_cur_buf;
        auto &next = pct_next_buf;
        {
            const auto &first = classToPaths[counts[smallestIdx].first];
            cur.assign(first.begin(), first.end());
        }

        // Sorted-list intersection against the remaining k-1 class lists.
        for (int i = 0; i < (int)counts.size() && !cur.empty(); ++i) {
            if (i == smallestIdx) continue;
            const auto &rhs = classToPaths[counts[i].first];
            next.clear();
            size_t a = 0, b = 0;
            while (a < cur.size() && b < rhs.size()) {
                if (cur[a] < rhs[b]) ++a;
                else if (cur[a] > rhs[b]) ++b;
                else { next.push_back(cur[a]); ++a; ++b; }
            }
            std::swap(cur, next);  // preserves both capacities
        }

        // Multiplicity filter. findClassIdx is guaranteed ≥ 0 for every
        // (piIdx, c) in the intersection on alive paths. Retired paths
        // are pre-filtered: their classIds/nh/np have been freed and a
        // findClassIdx would dereference an empty vector.
        out.reserve(cur.size());
        for (auto piIdx : cur) {
            if (pathRetired[piIdx]) continue;
            const auto &pi = pathInfos[piIdx];
            bool ok = true;
            for (auto &[c, jc] : counts) {
                int idx = findClassIdx(pi, c);
                if (jc > pi.nh[idx] + pi.np[idx]) { ok = false; break; }
            }
            if (ok) out.push_back(piIdx);
        }
    };


    // --- Bucket queue setup ---
    std::vector<double> dSup = support;
    std::vector<double>().swap(support);  // LowMem: `support` is read only to seed dSup
    std::vector<bool> rPeeled(rTuples.size(), false);

    // LowMem path-retirement: track alive tuple count per path. When it drops
    // to 0, the path's deadBoxes/tupleIdxs are useless forever (any future
    // pathsCoveringTuple(t) hit implies t ∈ pi.tupleIdxs at build; if count==0,
    // all such t are already peeled, so pi is never consulted again).
    pathAliveCount.resize(pathInfos.size());
    for (daf::Size i = 0; i < pathInfos.size(); ++i)
        pathAliveCount[i] = (daf::Size)pathInfos[i].tupleIdxs.size();
    pathRetired.assign(pathInfos.size(), false);

    // Profiling counters
    long long prof_tupleUpdates = 0, prof_pathsRetired = 0;
    long long prof_bulkDeaths = 0, prof_partialDeaths = 0;
    long long prof_dyingCells = 0, prof_corrPairs = 0;
    long long prof_cellsMaterialized = 0, prof_pathsMaterialized = 0;

    // Bucket levels (support counts, core values) are stored as double
    // throughout. Rationale:
    //   - mcCoreVal = C(|M|-r, s-r) can reach 1e20+ on dense graphs (e.g.
    //     com-dblp max MC=114 gives C(111, 17) ≈ 4.5e19 — overflows uint64's
    //     full range when r grows further).
    //   - `llround` returns long long which is signed int64; any value ≥
    //     2^63 overflows it and wraps to INT64_MIN.
    //   - Step 1b and Step 3 already store core values as `double` into
    //     `coreDist`; using uint64 for peel buckets creates two different
    //     core-value spaces (clamped vs unclamped) that disagree on the
    //     labels assigned to the same r-cliques in V3 vs V3_NP.
    // Double loses precision beyond 2^53, but that is intrinsic to the
    // nCr table anyway. Using double uniformly makes V3 and V3_NP agree,
    // and agree with ST (which also uses double).
    // Bucket levels can exceed 2^32 on dense graphs (mcCoreVal ≈ 4.5e19 for
    // com-dblp max MC=114, r=3, s=20). Use uint64_t to avoid silent uint32
    // truncation. `llround` returns int64_t and wraps past 2^63; clamp via
    // INT64_MAX to avoid that trap. Core values above INT64_MAX are labelled
    // identically in V3/V3_NP (both clamp); Step 1b and Step 3 still store
    // the unclamped `double` core values into coreDist, so the output still
    // distinguishes the true large cores for directly-assigned r-cliques.
    using BigLevel = uint64_t;
    auto safeToBigLevel = [](double x) -> BigLevel {
        if (!(x > 0.0)) return 0;
        if (x >= (double)INT64_MAX) return (BigLevel)INT64_MAX;
        return (BigLevel)llround(x);
    };

    // Dense bucket array sizing. The array is indexed by support *value*, so
    // a naive cap of min(maxEffSup+2, 1000001) allocates up to 1,000,001
    // empty std::vector headers (~24 MB) whenever any support is large -- and
    // supports ARE large at high s (an active tuple's support is
    // C(|M|-r,s-r), e.g. ~1e9 at r=6,s=20, and a region floor can reach
    // ~1e19). On a cell with only a handful of tuples that dense array is
    // almost entirely empty: it is the s-independent ~24 MB memory plateau on
    // the matched-cell memory curves. Since there are at most |rTuples|
    // non-empty buckets, also cap the dense array by the tuple count; tuples
    // whose support exceeds the cap are peeled last and live in `overflow`,
    // which the peel loop and refreshAffectedPaths already consult, so the
    // output is unchanged.
    BigLevel maxRealSup = 0;
    std::vector<BigLevel> effSup(rTuples.size());
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        BigLevel sv = safeToBigLevel(dSup[i]);
        effSup[i] = std::max(sv, safeToBigLevel(tupleMinCore[i]));
        maxRealSup = std::max(maxRealSup, sv);
    }
    const BigLevel BUCKET_CAP = std::min<BigLevel>(
        std::min<BigLevel>(maxRealSup + 2, (BigLevel)1000001),
        (BigLevel)rTuples.size() + 2);

    std::vector<std::vector<daf::Size>> buckets(BUCKET_CAP);
    std::multimap<BigLevel, daf::Size> overflow;
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        if (effSup[i] < BUCKET_CAP) buckets[effSup[i]].push_back(i);
        else overflow.insert({effSup[i], i});
    }

    // Profiling counters
    long long prof_batchCount = 0;


    // Profile counters for the path-level zero-base fast path. Sparse
    // graphs hit base=0 for ~99% of tuples (the tuple's key shares no
    // class with the path, or the shared classes don't exceed nh) — so
    // caching normalize+choosePivot once per (path, batch) eliminates
    // a per-tuple ~150 ns work item across millions of tuples.
    long long prof_zeroBaseHits = 0, prof_zeroBaseFills = 0;
    // Coarser timers for the outer peel loop sections (cheap — only
    // called once per batch):
    //   pctPhase   : pathsCoveringTuple loop (build batchDead) + deadCache eviction
    //   appendPhase: appendReqBox loop (DomPrune work) + per-path append setup
    //   refreshTot : whole refreshAffectedPaths call
    long long prof_pctPhaseNs = 0, prof_appendPhaseNs = 0;
    long long prof_refreshTotalNs = 0;
    // Phase B.1 — Theorem 1 (Tuple-on-Path Saturation) detection counter.
    // Counts how many (P, τ) pairs reach alive(τ,P) ≤ 0.5, get removed
    // from the path's tupleIdxs, and stop being refreshed.
    long long prof_tupleSaturated = 0;

    // PIVOTER_TRACE_TUPLE="i,j,k": log every dSup change for those tuples.
    std::unordered_set<daf::Size> traceSet;
    if (const char *ts = std::getenv("PIVOTER_TRACE_TUPLE")) {
        std::string s2(ts); size_t pos = 0;
        while (pos < s2.size()) {
            size_t e = s2.find(',', pos);
            if (e == std::string::npos) e = s2.size();
            traceSet.insert((daf::Size)std::stoul(s2.substr(pos, e - pos)));
            pos = e + 1;
        }
        for (auto t : traceSet)
            std::cerr << "TRACE eng=EV init tidx=" << t << " dSup=" << dSup[t]
                      << " minCore=" << tupleMinCore[t] << "\n";
    }

    // PIVOTER_TRACE_PATH="p,q": log dying-tuple dispatch onto those paths.
    std::unordered_set<daf::Size> tracePathSet;
    if (const char *ps = std::getenv("PIVOTER_TRACE_PATH")) {
        std::string s2(ps); size_t pos = 0;
        while (pos < s2.size()) {
            size_t e = s2.find(',', pos);
            if (e == std::string::npos) e = s2.size();
            tracePathSet.insert((daf::Size)std::stoul(s2.substr(pos, e - pos)));
            pos = e + 1;
        }
    }

    // ====================================================================
    // EventPeel engine (SigmodPlus §10.6): factored cell-death events.
    //
    // A cell on path P is one pivot-count vector y (sum y_k = T,
    // 0 <= y_k <= np_k) — one orbit of s-cliques under within-class
    // symmetry. Peeling tau* kills the cells y >= ell(tau*|P). For a
    // surviving tuple tau', each dying cell's contribution factorizes:
    //     W(y) = g(y) * prod_{k in touched} C(nh_k + y_k, j_k),
    //     g(y) = prod_{k in S(y)} C(np_k, y_k),
    // so with the static base = prod_{k in touched} C(nh_k, j_k) the
    // batch delta is
    //     delta = G*base + sum_{dying y : S(y) cap touched != empty}
    //                        g(y)*(prod_k C(nh_k+y_k, j_k) - base),
    // where G = sum_dying g(y). base = 0 exactly for pivot-needy tuples,
    // so the per-survivor corner constraint y >= ell(tau') needs no
    // special handling. Bulk deaths (ell(tau*) == 0) never touch cells:
    // delta = liveContrib[P][tau'] and the path retires (Theorem 5).
    // Each cell dies exactly once over the whole peel.
    // ====================================================================
    auto coordK = [](uint32_t cw) -> int { return (int)(cw >> 12); };
    auto coordY = [](uint32_t cw) -> int { return (int)(cw & 0xFFFu); };
    auto CC = [](int n, int k) -> double {
        return (k < 0 || n < 0 || k > n) ? 0.0 : nCr[n][k];
    };

    auto materializeCells = [&](PathInfo &pi) {
        const int m = (int)pi.classIds.size();
        pi.cellStart.clear();
        pi.cellCoords.clear();
        pi.cellStart.push_back(0);
        std::vector<int> sufCap((size_t)m + 1, 0);
        for (int k = m - 1; k >= 0; --k) sufCap[k] = sufCap[k + 1] + pi.np[k];
        std::vector<uint32_t> cur;
        std::function<void(int, int)> rec = [&](int k, int rem) {
            if (rem == 0) {
                pi.cellCoords.insert(pi.cellCoords.end(), cur.begin(), cur.end());
                pi.cellStart.push_back((uint32_t)pi.cellCoords.size());
                return;
            }
            if (k >= m || sufCap[k] < rem) return;
            rec(k + 1, rem);
            const int cap = std::min(pi.np[k], rem);
            for (int y = 1; y <= cap; ++y) {
                cur.push_back(((uint32_t)k << 12) | (uint32_t)y);
                rec(k + 1, rem - y);
                cur.pop_back();
            }
        };
        rec(0, pi.T);
        const size_t nc = pi.cellStart.size() - 1;
        pi.cellAlive.assign(nc, (uint8_t)1);
        pi.cellEpoch.assign(nc, 0u);
        pi.evEpoch = 0;
        pi.aliveCells = (daf::Size)nc;
        pi.cellsByClass.assign((size_t)m, {});
        for (uint32_t ci = 0; ci < (uint32_t)nc; ++ci)
            for (uint32_t o = pi.cellStart[ci]; o < pi.cellStart[ci + 1]; ++o)
                pi.cellsByClass[(size_t)coordK(pi.cellCoords[o])].push_back(ci);
        pi.cellsBuilt = true;
        prof_cellsMaterialized += (long long)nc;
        prof_pathsMaterialized++;
    };

    auto retirePath = [&](daf::Size piIdx, PathInfo &pi) {
        std::vector<daf::Size>().swap(pi.tupleIdxs);
        std::vector<double>().swap(pi.liveContrib);
        std::vector<double>().swap(pi.tupleBase);
        std::vector<daf::Size>().swap(pi.classIds);
        std::vector<int>().swap(pi.nh);
        std::vector<int>().swap(pi.np);
        std::vector<uint32_t>().swap(pi.cellCoords);
        std::vector<uint32_t>().swap(pi.cellStart);
        std::vector<uint8_t>().swap(pi.cellAlive);
        std::vector<uint32_t>().swap(pi.cellEpoch);
        std::vector<std::vector<uint32_t>>().swap(pi.cellsByClass);
        pathRetired[piIdx] = true;
        prof_pathsRetired++;
    };

    // Event scratch, reused across paths and batches.
    std::vector<int> evEll;              // ell per local class idx
    std::vector<int> evEllIdx;           // local class idxs with ell > 0
    std::vector<uint32_t> evDying;       // dying cell ids this (path, batch)
    std::vector<double> evDyingG;        // g(y) per dying cell (parallel)
    std::vector<uint8_t> evMarkedG(numClasses, 0);  // global cid -> marked
    std::vector<daf::Size> evMarkedList;
    // Single-touch histogram factorization: A[k][v] = Σ_{dying y: y_k=v} g(y)
    // (exact for tuples touching exactly one marked class), and per-class
    // dying index lists for the multi-touch exact loop.
    std::vector<double> evA;                       // m × (T+1), flat
    std::vector<int> evMarkedLocal;                // local class idxs marked
    std::vector<uint8_t> evMarkedLocalG;           // local idx -> marked
    std::vector<std::vector<uint32_t>> evDyingBy;  // local k -> dying indices
    std::vector<int> evTk;               // tuple's local class idxs
    std::vector<int> evTj;               // tuple's per-class j
    // Multi-touch signature memo: tuples sharing the same
    // (marked-touched classes, j) signature share one correction core
    //   corrCore(K, j|K) = Σ_dying g·(Π_{k∈K} C(nh+y,j) − Π_{k∈K} C(nh,j)),
    // and differ only by the scalar baseOther = Π_{touched∖K} C(nh,j).
    struct SigVecHash {
        size_t operator()(const std::vector<uint32_t> &v) const {
            size_t h = 1469598103934665603ull;
            for (uint32_t x : v) { h ^= x; h *= 1099511628211ull; }
            return h;
        }
    };
    std::unordered_map<std::vector<uint32_t>, double, SigVecHash> evSigMemo;
    std::vector<uint32_t> evSigKey;
    long long prof_sigHits = 0, prof_sigMiss = 0;

    auto processPath = [&](daf::Size piIdx, const std::vector<daf::Size> &deadTuples,
                           BigLevel &scanLevel) {
        auto &pi = pathInfos[piIdx];
        if (!tracePathSet.empty() && tracePathSet.count(piIdx))
            for (auto d : deadTuples)
                std::cerr << "TRACE eng=EV PCT path=" << piIdx
                          << " dying=" << d << "\n";
        // (1) drop peeled tuples from the slot arrays (lockstep).
        {
            size_t w = 0;
            for (size_t slot = 0; slot < pi.tupleIdxs.size(); ++slot) {
                if (rPeeled[pi.tupleIdxs[slot]]) continue;
                pi.tupleIdxs[w]   = pi.tupleIdxs[slot];
                pi.liveContrib[w] = pi.liveContrib[slot];
                pi.tupleBase[w]   = pi.tupleBase[slot];
                ++w;
            }
            pi.tupleIdxs.resize(w);
            pi.liveContrib.resize(w);
            pi.tupleBase.resize(w);
        }
        if (pi.tupleIdxs.empty()) { retirePath(piIdx, pi); return; }

        // (2) classify corners; kill cells for partial corners as we go.
        bool bulk = false;
        const int mPath = (int)pi.classIds.size();
        const int Trow = pi.T + 1;
        if ((int)evMarkedLocalG.size() < mPath) evMarkedLocalG.resize(mPath, 0);
        if ((int)evDyingBy.size() < mPath) evDyingBy.resize(mPath);
        if ((size_t)evA.size() < (size_t)mPath * Trow)
            evA.assign((size_t)mPath * Trow, 0.0);
        // Event scratch is cleaned at the END of each event (with the
        // stride in effect then); see clearEventScratch below.
        auto clearEventScratch = [&]() {
            for (int k2 : evMarkedLocal) {
                std::fill(evA.begin() + (size_t)k2 * Trow,
                          evA.begin() + (size_t)(k2 + 1) * Trow, 0.0);
                evDyingBy[k2].clear();
                evMarkedLocalG[k2] = 0;
            }
            evMarkedLocal.clear();
            for (auto cid : evMarkedList) evMarkedG[cid] = 0;
            evMarkedList.clear();
            evDying.clear();
            evDyingG.clear();
            evSigMemo.clear();
        };
        double G = 0.0;
        if ((int)evEll.size() < (int)pi.classIds.size())
            evEll.resize(pi.classIds.size());

        for (auto deadIdx : deadTuples) {
            auto key = keyOf(deadIdx);
            evEllIdx.clear();
            bool feasible = true;
            for (size_t kk = 0; kk < key.size() && feasible; ) {
                daf::Size c = key[kk];
                int jc = 1;
                while (kk + jc < key.size() && key[kk + jc] == c) ++jc;
                int k = findClassIdx(pi, c);
                if (k < 0 || jc - pi.nh[k] > pi.np[k]) { feasible = false; break; }
                if (jc > pi.nh[k]) { evEll[k] = jc - pi.nh[k]; evEllIdx.push_back(k); }
                kk += jc;
            }
            if (!feasible) continue;
            if (evEllIdx.empty()) { bulk = true; break; }   // ell == 0: Theorem 5

            // partial corner: kill alive cells with y >= ell
            if (!pi.cellsBuilt) materializeCells(pi);
            int kStar = evEllIdx[0];
            for (int k2 : evEllIdx)
                if (pi.cellsByClass[k2].size() < pi.cellsByClass[kStar].size())
                    kStar = k2;
            for (uint32_t ci : pi.cellsByClass[(size_t)kStar]) {
                if (!pi.cellAlive[ci]) continue;
                bool dom = true;
                for (int kNeed : evEllIdx) {
                    int yv = 0;
                    for (uint32_t o = pi.cellStart[ci]; o < pi.cellStart[ci + 1]; ++o)
                        if (coordK(pi.cellCoords[o]) == kNeed) {
                            yv = coordY(pi.cellCoords[o]);
                            break;
                        }
                    if (yv < evEll[kNeed]) { dom = false; break; }
                }
                if (!dom) continue;
                pi.cellAlive[ci] = 0;
                pi.aliveCells--;
                double g = 1.0;
                for (uint32_t o = pi.cellStart[ci]; o < pi.cellStart[ci + 1]; ++o)
                    g *= CC(pi.np[coordK(pi.cellCoords[o])],
                            coordY(pi.cellCoords[o]));
                const uint32_t dIdx = (uint32_t)evDying.size();
                evDying.push_back(ci);
                evDyingG.push_back(g);
                for (uint32_t o = pi.cellStart[ci]; o < pi.cellStart[ci + 1]; ++o) {
                    const int k2 = coordK(pi.cellCoords[o]);
                    const int yv = coordY(pi.cellCoords[o]);
                    const daf::Size cid2 = pi.classIds[k2];
                    if (!evMarkedG[cid2]) { evMarkedG[cid2] = 1; evMarkedList.push_back(cid2); }
                    if (!evMarkedLocalG[k2]) {
                        evMarkedLocalG[k2] = 1;
                        evMarkedLocal.push_back(k2);
                    }
                    evA[(size_t)k2 * Trow + yv] += g;
                    evDyingBy[k2].push_back(dIdx);
                }
                G += g;
                prof_dyingCells++;
            }
        }

        // (3) apply the batch delta to surviving slots.
        if (bulk) {
            prof_bulkDeaths++;
            for (size_t slot = 0; slot < pi.tupleIdxs.size(); ++slot) {
                const daf::Size tidx = pi.tupleIdxs[slot];
                const double delta = pi.liveContrib[slot];
                if (delta > 0.5) {
                    dSup[tidx] -= delta / rTuples[tidx].mult;
                    if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                    if (!traceSet.empty() && traceSet.count(tidx))
                        std::cerr << "TRACE eng=EV updBULK path=" << piIdx
                                  << " tidx=" << tidx
                                  << " delta=" << (delta / rTuples[tidx].mult)
                                  << " dSup=" << dSup[tidx] << "\n";
                    prof_tupleUpdates++;
                    const BigLevel newSup = safeToBigLevel(dSup[tidx]);
                    const BigLevel newBucket =
                        std::max(newSup, safeToBigLevel(tupleMinCore[tidx]));
                    effSup[tidx] = newBucket;
                    if (newBucket < BUCKET_CAP) {
                        buckets[newBucket].push_back(tidx);
                        if (newBucket < scanLevel) scanLevel = newBucket;
                    } else {
                        overflow.insert({newBucket, tidx});
                    }
                }
            }
            clearEventScratch();
            retirePath(piIdx, pi);
            return;
        }
        prof_partialDeaths++;
        if (evDying.empty()) return;

        size_t w = 0;
        for (size_t slot = 0; slot < pi.tupleIdxs.size(); ++slot) {
            const daf::Size tidx = pi.tupleIdxs[slot];
            double delta = G * pi.tupleBase[slot];

            // corrections, only when tau' touches a marked class
            auto key = keyOf(tidx);
            bool marked = false;
            for (auto c : key)
                if (evMarkedG[c]) { marked = true; break; }
            if (marked) {
                evTk.clear(); evTj.clear();
                int nMarked = 0, kSingle = -1, jSingle = 0;
                for (size_t kk = 0; kk < key.size(); ) {
                    daf::Size c = key[kk];
                    int jc = 1;
                    while (kk + jc < key.size() && key[kk + jc] == c) ++jc;
                    const int k = findClassIdx(pi, c);
                    evTk.push_back(k);
                    evTj.push_back(jc);
                    if (evMarkedLocalG[k]) { ++nMarked; kSingle = k; jSingle = jc; }
                    kk += jc;
                }
                const double base = pi.tupleBase[slot];
                if (nMarked == 1) {
                    // Exact single-touch histogram correction:
                    //   corr = baseOther · Σ_v A[k][v]·(C(nh+v,j) − C(nh,j)),
                    // baseOther = Π_{touched≠k} C(nh,j) (no division).
                    double baseOther = 1.0;
                    for (size_t t2 = 0; t2 < evTk.size(); ++t2)
                        if (evTk[t2] != kSingle)
                            baseOther *= CC(pi.nh[evTk[t2]], evTj[t2]);
                    const double c0 = CC(pi.nh[kSingle], jSingle);
                    double acc = 0.0;
                    const double *Arow = evA.data() + (size_t)kSingle * Trow;
                    for (int v = 1; v < Trow; ++v) {
                        if (Arow[v] == 0.0) continue;
                        acc += Arow[v] * (CC(pi.nh[kSingle] + v, jSingle) - c0);
                    }
                    delta += baseOther * acc;
                    prof_corrPairs++;
                } else if (nMarked >= 2) {
                    // Multi-touch via signature memo: the correction core
                    // depends only on the marked-touched (k, j) pairs, so
                    // tuples sharing that signature share one cell loop.
                    evSigKey.clear();
                    double baseOther = 1.0;
                    for (size_t t2 = 0; t2 < evTk.size(); ++t2) {
                        const int kT = evTk[t2];
                        if (evMarkedLocalG[kT])
                            evSigKey.push_back(((uint32_t)kT << 8) |
                                               (uint32_t)evTj[t2]);
                        else
                            baseOther *= CC(pi.nh[kT], evTj[t2]);
                    }
                    auto mit = evSigMemo.find(evSigKey);
                    double corrCore;
                    if (mit != evSigMemo.end()) {
                        corrCore = mit->second;
                        prof_sigHits++;
                    } else {
                        prof_sigMiss++;
                        double cK = 1.0;
                        for (uint32_t sk : evSigKey)
                            cK *= CC(pi.nh[sk >> 8], (int)(sk & 0xFF));
                        corrCore = 0.0;
                        const uint32_t epoch = ++pi.evEpoch;
                        for (uint32_t sk : evSigKey) {
                            const int kT = (int)(sk >> 8);
                            for (uint32_t dIdx : evDyingBy[kT]) {
                                const uint32_t ci = evDying[dIdx];
                                if (pi.cellEpoch[ci] == epoch) continue;
                                pi.cellEpoch[ci] = epoch;
                                double thetaProd = 1.0;
                                for (uint32_t sk2 : evSigKey) {
                                    const int k2 = (int)(sk2 >> 8);
                                    int yv = 0;
                                    for (uint32_t o = pi.cellStart[ci];
                                         o < pi.cellStart[ci + 1]; ++o)
                                        if (coordK(pi.cellCoords[o]) == k2) {
                                            yv = coordY(pi.cellCoords[o]);
                                            break;
                                        }
                                    thetaProd *= CC(pi.nh[k2] + yv,
                                                    (int)(sk2 & 0xFF));
                                }
                                corrCore += evDyingG[dIdx] * (thetaProd - cK);
                                prof_corrPairs++;
                            }
                        }
                        evSigMemo.emplace(evSigKey, corrCore);
                    }
                    delta += baseOther * corrCore;
                }
            }

            double lc = pi.liveContrib[slot] - delta;
            if (lc < 0) lc = 0;
            if (delta > 0.5) {
                dSup[tidx] -= delta / rTuples[tidx].mult;
                if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                if (!traceSet.empty() && traceSet.count(tidx))
                    std::cerr << "TRACE eng=EV updPART path=" << piIdx
                              << " tidx=" << tidx
                              << " delta=" << (delta / rTuples[tidx].mult)
                              << " G=" << G << " lc=" << lc
                              << " dSup=" << dSup[tidx] << "\n";
                prof_tupleUpdates++;
                const BigLevel newSup = safeToBigLevel(dSup[tidx]);
                const BigLevel newBucket =
                    std::max(newSup, safeToBigLevel(tupleMinCore[tidx]));
                effSup[tidx] = newBucket;
                if (newBucket < BUCKET_CAP) {
                    buckets[newBucket].push_back(tidx);
                    if (newBucket < scanLevel) scanLevel = newBucket;
                } else {
                    overflow.insert({newBucket, tidx});
                }
            }
            if (lc < 0.5) continue;   // exhausted on P: drop the slot
            pi.tupleIdxs[w]   = tidx;
            pi.liveContrib[w] = lc;
            pi.tupleBase[w]   = pi.tupleBase[slot];
            ++w;
        }
        pi.tupleIdxs.resize(w);
        pi.liveContrib.resize(w);
        pi.tupleBase.resize(w);
        clearEventScratch();
        if (pi.tupleIdxs.empty()) retirePath(piIdx, pi);
    };

    // Per-batch path grouping (sized once; reused).
    std::vector<std::vector<daf::Size>> batchDead(pathInfos.size());
    std::vector<daf::Size> touchedPaths;

    // PIVOTER_DUMP_TUPLE_CORE=<path>: per-active-tuple final core dump
    // ("idx mult kappa" per line) for cross-engine differential debugging.
    std::vector<double> dumpTupleCore;
    if (std::getenv("PIVOTER_DUMP_TUPLE_CORE"))
        dumpTupleCore.assign(rTuples.size(), -1.0);

    // --- Batch peeling loop ---
    daf::Size numPeeled = 0;
    BigLevel currentLevel = 0, coreLevel = 0;

    std::vector<daf::Size> pctPaths;  // pathsCoveringTuple output buffer

    while (numPeeled < rTuples.size()) {
        while (currentLevel < BUCKET_CAP && buckets[currentLevel].empty())
            currentLevel++;

        if (currentLevel >= BUCKET_CAP) {
            if (overflow.empty()) break;
            currentLevel = overflow.begin()->first;
        }

        std::vector<daf::Size> batch;
        if (currentLevel < BUCKET_CAP) {
            while (!buckets[currentLevel].empty()) {
                daf::Size idx = buckets[currentLevel].back();
                buckets[currentLevel].pop_back();
                if (rPeeled[idx]) continue;
                if (effSup[idx] != currentLevel) continue;
                rPeeled[idx] = true;
                batch.push_back(idx);
            }
        } else if (!overflow.empty()) {
            auto range = overflow.equal_range(currentLevel);
            for (auto it = range.first; it != range.second; ++it) {
                daf::Size idx = it->second;
                if (!rPeeled[idx]) {
                    if (effSup[idx] == currentLevel) {
                        rPeeled[idx] = true;
                        batch.push_back(idx);
                    }
                }
            }
            overflow.erase(range.first, range.second);
        } else break;

        if (batch.empty()) { currentLevel++; continue; }
        prof_batchCount++;

        coreLevel = std::max(coreLevel, currentLevel);
        for (auto idx : batch) {
            numPeeled++;
            // minCore floor already enforced by bucket placement (effSup ≥ minCore)
            coreDist[(double)coreLevel] += (int64_t)rTuples[idx].mult;
            if (!dumpTupleCore.empty()) dumpTupleCore[idx] = (double)coreLevel;
            if (!traceSet.empty() && traceSet.count(idx))
                std::cerr << "TRACE eng=EV PEEL tidx=" << idx
                          << " level=" << coreLevel << "\n";
        }

        // Clear last batch's scratch (capacity kept).
        for (auto pi2 : touchedPaths) batchDead[pi2].clear();
        touchedPaths.clear();

        // Group this batch's peeled tuples by hosting path.
        auto _t_pct0 = std::chrono::steady_clock::now();
        for (auto idx : batch) {
            pathsCoveringTuple(idx, pctPaths);
            for (auto piIdx : pctPaths) {
                if (batchDead[piIdx].empty()) touchedPaths.push_back(piIdx);
                batchDead[piIdx].push_back(idx);
            }
        }
        prof_pctPhaseNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - _t_pct0).count();

        // Fire the death events path by path.
        auto _t_ref0 = std::chrono::steady_clock::now();
        for (auto piIdx : touchedPaths)
            if (!pathRetired[piIdx])
                processPath(piIdx, batchDead[piIdx], currentLevel);
        prof_refreshTotalNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - _t_ref0).count();

        // PERF-AUDIT MEMORY (1/3 → superseded by 7): the per-peel
        // TupleKey().swap() that used to live here has been retired in
        // favour of the flat-key layout. With std::vector<daf::Size>
        // tupleKeysFlat, individual peeled tuples cannot return their
        // r-element slice to the allocator, but the per-tuple savings
        // is already collected up front by eliminating the 24-byte
        // std::vector header. Net result is comparable peak RSS for
        // typical paper r ≤ 12 with simpler bookkeeping.


        // PERF-AUDIT MEMORY (3/3): once batchDead[piIdx] for a retired path
        // has been processed, its inner vector capacity is dead weight
        // (path's pathAliveCount==0 path was retired in the affected loop
        // and won't be touched again). Free those inner vectors.
        for (auto piIdx : touchedPaths) {
            if (pathAliveCount[piIdx] == 0) {
                std::vector<daf::Size>().swap(batchDead[piIdx]);
            }
        }
    }

    auto tStep6End = std::chrono::high_resolution_clock::now();
    auto step6Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStep6).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStart).count();

    double maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling (Batch Union) ---" << std::endl;
    std::cout << "  Peeled active tuples: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Batches: " << prof_batchCount << std::endl;
    std::cout << "  Paths retired (alive=0): " << prof_pathsRetired
              << " / " << pathInfos.size() << std::endl;
    std::cout << "  Bulk deaths (Theorem-5 retires): " << prof_bulkDeaths << std::endl;
    std::cout << "  Partial death events: " << prof_partialDeaths << std::endl;
    std::cout << "  Cells died: " << prof_dyingCells
              << "  (materialized " << prof_cellsMaterialized
              << " on " << prof_pathsMaterialized << " paths)" << std::endl;
    std::cout << "  Sparse correction pairs: " << prof_corrPairs
              << "  (sig memo " << prof_sigHits << " hits / "
              << prof_sigMiss << " miss)" << std::endl;
    {
        double peelMs = (double)step6Ms;
        auto pct = [peelMs](long long ns) -> double {
            return peelMs > 0 ? (double)ns / 1e6 / peelMs * 100.0 : 0.0;
        };
        std::cout << "  --- Peel phase breakdown (outer) ---" << std::endl;
        std::cout << "    pctPhase  : " << prof_pctPhaseNs / 1000000 << " ms ("
                  << pct(prof_pctPhaseNs) << "%)" << std::endl;
        std::cout << "    eventsTot : " << prof_refreshTotalNs / 1000000 << " ms ("
                  << pct(prof_refreshTotalNs) << "%)" << std::endl;
    }
    std::cout << "  Tuple updates: " << prof_tupleUpdates << std::endl;
    if (enablePrivateCloud) {
        std::cout << "  Direct private r-cliques: "
                  << std::fixed << std::setprecision(0) << privateRCliquesDirect << std::endl;
    }
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step6Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;
    daf::phaseMark("Peel");

    if (!dumpTupleCore.empty()) {
        std::ofstream df(std::getenv("PIVOTER_DUMP_TUPLE_CORE"));
        df << std::setprecision(17);
        for (daf::Size i = 0; i < rTuples.size(); ++i)
            df << i << " " << rTuples[i].mult << " " << dumpTupleCore[i] << "\n";
    }

    // Return compact format: one entry per core level, key[0] = count
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    for (auto &[c, cnt] : coreDist)
        // Compact format: key = {lo32, hi32} encoding int64_t count
        result.push_back({{(daf::Size)(cnt & 0xFFFFFFFF), (daf::Size)((cnt >> 32) & 0xFFFFFFFF)}, (double)c});
    return result;
}
