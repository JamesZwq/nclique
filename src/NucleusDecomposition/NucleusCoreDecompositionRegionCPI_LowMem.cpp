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
NucleusCoreDecompositionRClique_RegionCPI_LowMem(
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

    if (enableDebugVerify) {
    auto tStep5 = std::chrono::high_resolution_clock::now();

    // --- Constrained Path data structure ---
    struct ClassBounds { int nh, np; int minPiv, maxPiv; }; // per-class on a constrained path
    struct CPath {
        int h; // total holds
        std::unordered_map<daf::Size, ClassBounds> classes; // cid -> bounds
        std::vector<daf::Size> tupleIdxs; // alive tuples on this constrained path
        int totalP() const { int t=0; for(auto&[_,b]:classes) t+=b.np; return t; }
    };

    // Compute AggrCount(τ', P̂) / mult(τ') on a constrained path
    // OPTIMIZED: convolve tuple classes (O(r²)) → h[0..r].
    // Non-tuple unconstrained classes: combine into C(Pfree, T-k) lookup.
    // Total: O(r² + κ_constrained × r) instead of O(classes × T).
    auto aggrCountOnCPath = [&](daf::Size tidx, const CPath &cp) -> double {
        auto key = keyOf(tidx);
        int T = s - cp.h;
        if (T < 0) return 0.0;

        // Phase 1: Convolve tuple classes → h[0..r], and sum Pfree
        double h[32]; int hLen = 1; h[0] = 1.0;
        int Pfree = 0;

        // First pass: process tuple classes
        int ki = 0;
        while (ki < (int)key.size()) {
            daf::Size c = key[ki];
            int jc = 1;
            while (ki + jc < (int)key.size() && key[ki + jc] == c) jc++;
            ki += jc;

            auto cit = cp.classes.find(c);
            if (cit == cp.classes.end()) return 0.0;
            auto &cb = cit->second;
            int tMin = std::max(cb.minPiv, std::max(0, jc - cb.nh));
            int tMax = std::min(cb.maxPiv, cb.np);
            if (jc > cb.nh + tMax || tMin > tMax) return 0.0;

            double gc[32]; int gcLen = tMax + 1;
            for (int i = 0; i < gcLen; ++i) gc[i] = 0.0;
            for (int tc = tMin; tc <= tMax; ++tc)
                gc[tc] = nCr[cb.np][tc] * nCr[cb.nh + tc][jc];
            int nL = std::min(hLen + gcLen - 1, (int)r + 1);
            double nh2[32];
            for (int i = 0; i < nL; ++i) nh2[i] = 0.0;
            for (int i = 0; i < hLen; ++i)
                for (int j = 0; j < gcLen && i+j < nL; ++j)
                    nh2[i+j] += h[i] * gc[j];
            hLen = nL;
            for (int i = 0; i < hLen; ++i) h[i] = nh2[i];
        }

        // Phase 2: non-tuple classes
        for (auto &[c, cb] : cp.classes) {
            // Check if c is in tuple (inline, no hash set)
            bool inTuple = false;
            for (auto x : key) if (x == c) { inTuple = true; break; }
            if (inTuple) continue;

            if (cb.minPiv == 0 && cb.maxPiv >= cb.np) {
                // Unconstrained → fold into Pfree
                Pfree += cb.np;
            } else {
                // Constrained non-tuple class: convolve with h
                int tMin = cb.minPiv, tMax = std::min(cb.maxPiv, cb.np);
                if (tMin > tMax) return 0.0;
                double gc[256]; int gcLen = tMax + 1;
                for (int i = 0; i < gcLen; ++i) gc[i] = 0.0;
                for (int tc = tMin; tc <= tMax; ++tc)
                    gc[tc] = nCr[cb.np][tc];
                int nL = std::min(hLen + gcLen - 1, T + 1);
                double nh2[256];
                for (int i = 0; i < nL; ++i) nh2[i] = 0.0;
                for (int i = 0; i < hLen; ++i)
                    for (int j = 0; j < gcLen && i+j < nL; ++j)
                        nh2[i+j] += h[i] * gc[j];
                hLen = nL;
                for (int i = 0; i < hLen; ++i) h[i] = nh2[i];
            }
        }

        // Phase 3: Σ h[k] × C(Pfree, T-k)
        double aggr = 0.0;
        for (int k = 0; k < hLen; ++k) {
            if (h[k] == 0.0) continue;
            int rem = T - k;
            if (rem >= 0 && rem <= Pfree)
                aggr += h[k] * nCr[Pfree][rem];
        }
        return aggr / rTuples[tidx].mult;
    };

    // --- Initialize constrained paths from CPI paths ---
    // Use a pool of CPath objects. Each tuple tracks which cpaths it's on.
    std::vector<CPath> cpaths;
    std::vector<std::vector<daf::Size>> tupleToCPaths(rTuples.size()); // tuple -> cpath IDs

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        CPath cp;
        cp.h = 0;
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            auto &cb = cp.classes[cid];
            if (node.isPivot) { cb.np++; }
            else { cb.nh++; cp.h++; }
        }
        // Set initial bounds
        for (auto &[cid, cb] : cp.classes) { cb.minPiv = 0; cb.maxPiv = cb.np; }

        // Find tuples on this path (enumerate r-multisets of path's classes)
        std::vector<daf::Size> pathClasses;
        for (auto &[cid, _] : cp.classes) {
            if (isActiveTupleClass(cid))
                pathClasses.push_back(cid);
        }
        std::sort(pathClasses.begin(), pathClasses.end());

        daf::Size cpid = cpaths.size();
        TupleKey cur; cur.reserve(r);
        bool hasTuples = false;

        std::function<void()> enumCb = [&]() {
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;
            // Feasibility check
            std::unordered_map<daf::Size, int> counts;
            for (auto c : cur) counts[c]++;
            for (auto &[c, jc] : counts) {
                auto cit = cp.classes.find(c);
                if (cit == cp.classes.end()) return;
                if (jc > cit->second.nh + cit->second.np) return;
            }
            cp.tupleIdxs.push_back(tidx);
            hasTuples = true;
        };
        enumerateMultisets(pathClasses, r, 0, cur, enumCb);

        if (hasTuples) {
            cpaths.push_back(std::move(cp));
            for (auto tidx : cpaths.back().tupleIdxs)
                tupleToCPaths[tidx].push_back(cpid);
        }
    }

    auto tStep5End = std::chrono::high_resolution_clock::now();
    auto step5Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStep5).count();
    std::cout << "  Constrained paths: " << cpaths.size() << std::endl;
    std::cout << "  Index build time: " << step5Ms << " ms" << std::endl;

    // Verify: aggrCountOnCPath matches Step 4 support
    {
        std::vector<double> supCheck(rTuples.size(), 0.0);
        for (daf::Size cpid = 0; cpid < cpaths.size(); ++cpid)
            for (auto tidx : cpaths[cpid].tupleIdxs)
                supCheck[tidx] += aggrCountOnCPath(tidx, cpaths[cpid]);
        double checkSum = 0;
        int mismatches = 0;
        for (daf::Size i = 0; i < rTuples.size(); ++i) {
            checkSum += rTuples[i].mult * supCheck[i];
            if (std::abs(supCheck[i] - support[i]) > 0.5) mismatches++;
        }
        std::cout << "  Support sum (CPath): " << std::fixed << std::setprecision(0) << checkSum << std::endl;
        std::cout << "  CPI vs CPath match: " << (mismatches == 0 ? "PASS" : ("MISMATCH(" + std::to_string(mismatches) + ")")) << std::endl;
    }
    }

    // ============================================================
    // Step 6: Batch Union Peeling (Weighted)
    // ============================================================
    // Instead of splitting CPath objects (which explode in number),
    // we accumulate "dead boxes" per path as tuples are peeled.
    // For each alive tuple τ' on an affected path, we compute:
    //   deadCount = countUnionWeighted(base, upper, deadBoxes, T, weights)
    //   delta = deadCount - cachedOldDeadCount
    //   dSup[τ'] -= delta / mult(τ')
    //
    // WEIGHTED counting: each allocation b contributes
    //   Π_c weight_c(b_c) × C(Pfree, T - Σb_c_tuple)
    // where weight_c(b) = C(nh_c, j'_c - b) × C(np_c, b) for tuple classes
    //                    = C(np_c, b) for non-tuple classes
    //
    // countFeasibleWeighted: convolution DP, O(T × m)
    // countUnionWeighted: branch-and-bound with Pareto pruning

    auto tStep6 = std::chrono::high_resolution_clock::now();

    // --- PathInfo: immutable per-path data + mutable dead boxes ---
    struct PathInfo {
        int h, T;                    // holds count, target = s - h
        std::vector<daf::Size> classIds;   // ordered class IDs on this path
        std::vector<int> nh;         // per-class hold count (parallel to classIds)
        std::vector<int> np;         // per-class pivot count (parallel to classIds)
        std::vector<daf::Size> tupleIdxs;  // alive tuples on this path
        // LowMem: flat pool of dead boxes as int16_t. Each box has exactly
        // |classIds| elements (= m), so box i lives at
        // deadBoxesFlat[i*m .. i*m + m). Values are req[c] = max(0, jc - nh_c)
        // where jc ≤ s and nh_c ≥ 0, so req fits comfortably in int16. Peak
        // savings: kills per-box heap header + capacity rounding (typ. 3–4×
        // smaller than vector<vector<int>> for the real m=5..15 case).
        std::vector<int16_t> deadBoxesFlat;
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

    // --- normalizeBoxes (flat int16 in/out) ---
    // LowMem: boxes flow through the whole branch-and-bound as a single
    // std::vector<int16_t> of length numBoxes * m (box i at offset i*m).
    // This is the same storage format as pi.deadBoxesFlat, so the top-level
    // call uses pi.deadBoxesFlat.data() directly with zero conversion. The
    // recursion's `remaining` and normalizeBoxes's `effective`/`minimal`
    // are also int16 flat — so no std::vector<std::vector<int>> ever exists
    // on the hot path. This closes the theoretical-regression corner case
    // of the previous scratch-based design (one-dominant-path graphs).
    // PERF AUDIT REWRITE — normalizeBoxes.
    // Previous version constructed up to four std::vectors per call
    // (effective / sorted / idx+sums / minimal). On union-heavy graphs
    // we hit ~10⁸ normalizeBoxes invocations, each with ≥ 2 small heap
    // allocs. New version owns one int16_t output buffer (still per-call
    // alloc, but we keep it short-lived and the recursion uses its
    // pointer) and routes everything else through depth-aware
    // thread_local scratch.
    //
    // Recursion safety note: countUnionRec recurses with the OUTPUT of
    // normalizeBoxes (norm.boxesFlat) still live across the recursive
    // call. Therefore the OUTPUT cannot be a thread_local buffer (the
    // child's normalizeBoxes call would clobber it). It stays as a
    // per-frame std::vector. Internal scratch (cur, idx, sums, sorted
    // workspace, minimal workspace) are thread_local because they are
    // not read across the recursive boundary.
    struct NormResult { bool fullCover; int numBoxes; };
    // PERF: caller-provided outBoxes eliminates per-call std::vector
    // allocation. countUnionRec maintains a depth-indexed pool of
    // outBoxes buffers (normPool below). Recursion at depth d writes
    // to normPool[d] and reads from normPool[d-1] (caller's buffer).
    auto normalizeBoxes = [](const int *lower, const int *upper, int T,
                             const int16_t *boxesIn, int numIn, int m,
                             bool pruneDom,
                             std::vector<int16_t> &outBoxes) -> NormResult {
        outBoxes.clear();
        // Fast exit: no input boxes ⇒ empty output, fast path.
        if (numIn == 0) return {false, 0};

        // Stage 1 — clamp/feasibility. cur is one row of int16_t scratch.
        thread_local std::vector<int16_t> curBuf;
        if ((int)curBuf.size() < m) curBuf.resize(m);
        int16_t *cur = curBuf.data();

        outBoxes.reserve((size_t)numIn * m);

        for (int bi = 0; bi < numIn; ++bi) {
            const int16_t *box = boxesIn + (size_t)bi * m;
            bool impossible = false;
            int lSum = 0;
            bool equalsLower = true;
            for (int i = 0; i < m; ++i) {
                int li = lower[i];
                int v = (int)box[i] > li ? (int)box[i] : li;
                if (v > upper[i]) { impossible = true; break; }
                cur[i] = (int16_t)v;
                lSum += v;
                if (v != li) equalsLower = false;
            }
            if (impossible || lSum > T) continue;
            if (equalsLower) { outBoxes.clear(); return {true, 0}; }
            outBoxes.insert(outBoxes.end(), cur, cur + m);
        }
        int n = (int)(outBoxes.size() / (size_t)m);
        if (n == 0) return {false, 0};

        // Stage 2 — sort + dedup (thread_local idx scratch).
        thread_local std::vector<int> idxBuf;
        if ((int)idxBuf.size() < n) idxBuf.resize(n);
        int *idx = idxBuf.data();
        std::iota(idx, idx + n, 0);
        const int16_t *eff = outBoxes.data();
        auto lexLess = [eff, m](int a, int b) {
            const int16_t *pa = eff + (size_t)a * m;
            const int16_t *pb = eff + (size_t)b * m;
            for (int i = 0; i < m; ++i)
                if (pa[i] != pb[i]) return pa[i] < pb[i];
            return false;
        };
        std::sort(idx, idx + n, lexLess);

        // sortedBuf is thread_local scratch — copied back INTO outBoxes
        // (the caller-provided output) before normalizeBoxes returns.
        thread_local std::vector<int16_t> sortedBuf;
        sortedBuf.clear();
        sortedBuf.reserve(outBoxes.size());
        for (int ii = 0; ii < n; ++ii) {
            const int16_t *candidate = outBoxes.data() + (size_t)idx[ii] * m;
            if (!sortedBuf.empty()) {
                const int16_t *prev = sortedBuf.data() + (sortedBuf.size() - m);
                bool same = true;
                for (int i = 0; i < m; ++i)
                    if (prev[i] != candidate[i]) { same = false; break; }
                if (same) continue;
            }
            sortedBuf.insert(sortedBuf.end(), candidate, candidate + m);
        }
        outBoxes.assign(sortedBuf.begin(), sortedBuf.end());
        n = (int)(outBoxes.size() / (size_t)m);
        if (!pruneDom) return {false, n};

        // Stage 3 — dominance pruning (thread_local sums + minimal scratch).
        thread_local std::vector<int> sumsBuf;
        if ((int)sumsBuf.size() < n) sumsBuf.resize(n);
        int *sums = sumsBuf.data();
        for (int i = 0; i < n; ++i) {
            const int16_t *p = outBoxes.data() + (size_t)i * m;
            int s = 0;
            for (int k = 0; k < m; ++k) s += p[k];
            sums[i] = s;
        }
        if ((int)idxBuf.size() < n) idxBuf.resize(n);
        idx = idxBuf.data();
        std::iota(idx, idx + n, 0);
        const int16_t *eff2 = outBoxes.data();
        std::sort(idx, idx + n, [eff2, sums, m](int a, int b) {
            if (sums[a] != sums[b]) return sums[a] < sums[b];
            const int16_t *pa = eff2 + (size_t)a * m;
            const int16_t *pb = eff2 + (size_t)b * m;
            for (int i = 0; i < m; ++i)
                if (pa[i] != pb[i]) return pa[i] < pb[i];
            return false;
        });

        thread_local std::vector<int16_t> minimalBuf;
        minimalBuf.clear();
        minimalBuf.reserve(outBoxes.size());
        int minCount = 0;
        for (int ii = 0; ii < n; ++ii) {
            const int16_t *cand = outBoxes.data() + (size_t)idx[ii] * m;
            bool dominated = false;
            for (int k = 0; k < minCount; ++k) {
                const int16_t *kept = minimalBuf.data() + (size_t)k * m;
                bool dom = true;
                for (int i = 0; i < m; ++i)
                    if (kept[i] > cand[i]) { dom = false; break; }
                if (dom) { dominated = true; break; }
            }
            if (!dominated) {
                minimalBuf.insert(minimalBuf.end(), cand, cand + m);
                ++minCount;
            }
        }
        outBoxes.assign(minimalBuf.begin(), minimalBuf.end());
        return {false, minCount};
    };

    // Depth-indexed pool of int16_t buffers for normalizeBoxes outputs.
    // Recursion at depth d writes to normPool[d]; reads from caller's
    // (= normPool[d-1]) are safe because parent doesn't re-enter
    // normalize between recursive calls.
    constexpr int NORM_POOL_DEPTH = 64;
    thread_local std::vector<std::vector<int16_t>> normPool(NORM_POOL_DEPTH);

    // --- countUnionWeighted: branch-and-bound ---
    // choosePivotBox: pick box with fewest active dims (ties: highest sum).
    // LowMem: reads from int16 flat boxes array directly.
    // PERF: pointer-based signature avoids forcing callers to wrap their
    // working buffers (often stack arrays) in std::vector.
    auto choosePivot = [](const int *lower,
                          const int16_t *boxes, int numBoxes, int m) -> int {
        int bestIdx = 0, bestActive = m + 1, bestSum = -1;
        for (int i = 0; i < numBoxes; ++i) {
            const int16_t *b = boxes + (size_t)i * m;
            int active = 0, lsum = 0;
            for (int c = 0; c < m; ++c) {
                int bv = b[c];
                if (bv > lower[c]) ++active;
                lsum += bv;
            }
            if (active < bestActive || (active == bestActive && lsum > bestSum)) {
                bestIdx = i; bestActive = active; bestSum = lsum;
            }
        }
        return bestIdx;
    };

    struct UnionCtx {
        int m, T;
        long long recCalls;
        daf::Size tidx;          // tuple being evaluated (for weight computation)
        const PathInfo *pi;       // path being evaluated
    };

    // Helper: compute weighted feasible for given lower/upper using tuple weights
    // LowMem: wtsFlat / wtLen moved from stack [MAX_M * (MAX_T+1)] to
    // thread_local heap so graphs with np_c > MAX_T (e.g. web-it-2004 MC=432,
    // np_c can be 300+) no longer abort. Per-call cost: two size-check
    // resizes; after warmup the buffers are at peak size, zero allocation.
    // PERF: pointer-based signature (was std::vector<int>&); the caller now
    // hands in stack arrays without ever constructing a vector.
    auto feasWeighted = [&](const int *lower, const int *upper,
                            UnionCtx &ctx) -> double {
        int m = (int)ctx.pi->classIds.size();
        // wtStride = max(upper[i] - lower[i] + 1, 1) to fit any class's range.
        int wtStride = 1;
        for (int i = 0; i < m; ++i) {
            int range = upper[i] - lower[i] + 1;
            if (range > wtStride) wtStride = range;
        }
        thread_local std::vector<double> wtsFlatBuf;
        thread_local std::vector<int> wtLenBuf;
        size_t need = (size_t)m * wtStride;
        if (wtsFlatBuf.size() < need) wtsFlatBuf.resize(need);
        if ((int)wtLenBuf.size() < m) wtLenBuf.resize(m);
        double *wtsFlat = wtsFlatBuf.data();
        int *wtLen = wtLenBuf.data();
        buildWeightTables_sb(ctx.tidx, *ctx.pi, lower, upper, wtsFlat, wtStride, wtLen);
        return countFeasibleWeighted_sb(lower, upper, m, ctx.T, wtsFlat, wtStride, wtLen);
    };

    // PERF AUDIT REWRITE — countUnionRec.
    // Previous version used std::function<...> for recursion (vtable-like
    // dispatch on every call) + per-recursion std::vector allocations for
    // pivot / nextLower / nextUpper. On union-heavy graphs the function-
    // call overhead and the heap traffic dominated the peel.
    //
    // New version:
    //   - Self-recursing lambda via Y-combinator style (auto&& self),
    //     so the call site is a direct lambda invocation, no
    //     std::function indirection.
    //   - Stack-allocated arrays for pivot / nextLower / nextUpper
    //     (sized MAX_M; caller already enforces this cap on path m).
    //   - Pointer-based signature throughout, no std::vector copying
    //     across recursion levels.
    auto countUnionRec_impl = [&](auto&& self,
                                  const int *lower, const int *upper,
                                  const int16_t *boxes, int numBoxes,
                                  UnionCtx &ctx, int depth) -> double {
        ctx.recCalls++;
        // Feasibility check.
        {
            int minS = 0, maxS = 0;
            for (int i = 0; i < ctx.m; ++i) {
                if (lower[i] > upper[i]) return 0.0;
                minS += lower[i]; maxS += upper[i];
            }
            if (ctx.T < minS || ctx.T > maxS) return 0.0;
        }

        // Use depth-indexed pool slot as outBoxes for normalize.
        // Recursion at depth+1 uses normPool[depth+1], so this slot's
        // contents survive across the recursive call.
        if (depth >= NORM_POOL_DEPTH) {
            // Defensive: deeper than expected. Allocate inline (rare path).
            std::vector<int16_t> tmp;
            auto norm = normalizeBoxes(lower, upper, ctx.T, boxes, numBoxes,
                                       ctx.m, true, tmp);
            if (norm.fullCover) return feasWeighted(lower, upper, ctx);
            if (norm.numBoxes == 0) return 0.0;
            // Fall through with `tmp` (rare; performance not critical).
            int pivIdx = choosePivot(lower, tmp.data(), norm.numBoxes, ctx.m);
            int stackPivot[MAX_M], stackNL[MAX_M], stackNU[MAX_M];
            std::vector<int> heapP, heapNL, heapNU;
            int *pivot = stackPivot, *nL = stackNL, *nU = stackNU;
            if (ctx.m > MAX_M) {
                heapP.resize(ctx.m); heapNL.resize(ctx.m); heapNU.resize(ctx.m);
                pivot = heapP.data(); nL = heapNL.data(); nU = heapNU.data();
            }
            const int16_t *pv = tmp.data() + (size_t)pivIdx * ctx.m;
            for (int i = 0; i < ctx.m; ++i) pivot[i] = (int)pv[i];
            double total = feasWeighted(pivot, upper, ctx);
            int remCount = norm.numBoxes - 1;
            if (pivIdx != remCount) {
                int16_t *p = tmp.data();
                int16_t *a = p + (size_t)pivIdx * ctx.m;
                int16_t *b = p + (size_t)remCount * ctx.m;
                for (int i = 0; i < ctx.m; ++i) std::swap(a[i], b[i]);
            }
            for (int splitDim = 0; splitDim < ctx.m; ++splitDim) {
                if (pivot[splitDim] <= lower[splitDim]) continue;
                for (int i = 0; i < ctx.m; ++i) { nL[i] = lower[i]; nU[i] = upper[i]; }
                for (int earlier = 0; earlier < splitDim; ++earlier)
                    if (pivot[earlier] > nL[earlier]) nL[earlier] = pivot[earlier];
                nU[splitDim] = std::min(nU[splitDim], pivot[splitDim] - 1);
                total += self(self, nL, nU, tmp.data(), remCount, ctx, depth + 1);
            }
            return total;
        }

        std::vector<int16_t> &outBoxes = normPool[depth];
        auto norm = normalizeBoxes(lower, upper, ctx.T, boxes, numBoxes, ctx.m,
                                   true, outBoxes);
        if (norm.fullCover) return feasWeighted(lower, upper, ctx);
        if (norm.numBoxes == 0) return 0.0;

        int pivIdx = choosePivot(lower, outBoxes.data(), norm.numBoxes, ctx.m);

        int stackPivot[MAX_M], stackNextLower[MAX_M], stackNextUpper[MAX_M];
        std::vector<int> heapPivot, heapNextLower, heapNextUpper;
        int *pivot     = stackPivot;
        int *nextLower = stackNextLower;
        int *nextUpper = stackNextUpper;
        if (ctx.m > MAX_M) {
            heapPivot.resize(ctx.m);
            heapNextLower.resize(ctx.m);
            heapNextUpper.resize(ctx.m);
            pivot     = heapPivot.data();
            nextLower = heapNextLower.data();
            nextUpper = heapNextUpper.data();
        }
        {
            const int16_t *pv = outBoxes.data() + (size_t)pivIdx * ctx.m;
            for (int i = 0; i < ctx.m; ++i) pivot[i] = (int)pv[i];
        }

        double total = feasWeighted(pivot, upper, ctx);

        // Zero-copy "remove pivot": swap pivot to the last slot of outBoxes
        // and pass (numBoxes - 1) to recursion.
        int remCount = norm.numBoxes - 1;
        if (pivIdx != remCount) {
            int16_t *p = outBoxes.data();
            int16_t *a = p + (size_t)pivIdx * ctx.m;
            int16_t *b = p + (size_t)remCount * ctx.m;
            for (int i = 0; i < ctx.m; ++i) std::swap(a[i], b[i]);
        }

        for (int splitDim = 0; splitDim < ctx.m; ++splitDim) {
            if (pivot[splitDim] <= lower[splitDim]) continue;

            for (int i = 0; i < ctx.m; ++i) {
                nextLower[i] = lower[i];
                nextUpper[i] = upper[i];
            }
            for (int earlier = 0; earlier < splitDim; ++earlier)
                if (pivot[earlier] > nextLower[earlier])
                    nextLower[earlier] = pivot[earlier];
            nextUpper[splitDim] = std::min(nextUpper[splitDim], pivot[splitDim] - 1);

            total += self(self, nextLower, nextUpper,
                          outBoxes.data(), remCount, ctx, depth + 1);
        }
        return total;
    };
    auto countUnionRec = [&](const int *lower, const int *upper,
                             const int16_t *boxes, int numBoxes,
                             UnionCtx &ctx) -> double {
        return countUnionRec_impl(countUnionRec_impl, lower, upper,
                                  boxes, numBoxes, ctx, 0);
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

    // Verify: weighted feasible on full path should match aggrCountOnCPath
    if (enableDebugVerify) {
        std::vector<double> supCheck(rTuples.size(), 0.0);
        for (daf::Size piIdx = 0; piIdx < pathInfos.size(); ++piIdx) {
            auto &pi = pathInfos[piIdx];
            int m = (int)pi.classIds.size();
            for (auto tidx : pi.tupleIdxs) {
                auto key = keyOf(tidx);
                std::unordered_map<daf::Size, int> counts;
                for (auto c : key) counts[c]++;
                std::vector<int> base(m), upper(m);
                for (int i = 0; i < m; ++i) {
                    int jc = 0;
                    auto cit = counts.find(pi.classIds[i]);
                    if (cit != counts.end()) jc = cit->second;
                    base[i] = std::max(0, jc - pi.nh[i]);
                    upper[i] = pi.np[i];
                }
                int wtStride = 1;
                for (int i = 0; i < m; ++i) {
                    int range = upper[i] - base[i] + 1;
                    if (range > wtStride) wtStride = range;
                }
                std::vector<double> wtsFlatV((size_t)m * wtStride);
                std::vector<int> wtLenV(m);
                double *wtsFlat = wtsFlatV.data();
                int *wtLen = wtLenV.data();
                buildWeightTables_sb(tidx, pi, base.data(), upper.data(), wtsFlat, wtStride, wtLen);
                double aggr = countFeasibleWeighted_sb(base.data(), upper.data(), m, pi.T, wtsFlat, wtStride, wtLen);
                supCheck[tidx] += aggr / rTuples[tidx].mult;
            }
        }
        double checkSum = 0;
        int mismatches = 0;
        for (daf::Size i = 0; i < rTuples.size(); ++i) {
            checkSum += rTuples[i].mult * supCheck[i];
            if (std::abs(supCheck[i] - support[i]) > 0.5) mismatches++;
        }
        std::cout << "  Support sum (PathInfo weighted): " << std::fixed << std::setprecision(0) << checkSum << std::endl;
        std::cout << "  CPI vs PathInfo match: " << (mismatches == 0 ? "PASS" : ("MISMATCH(" + std::to_string(mismatches) + ")")) << std::endl;
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

    // --- Append requirement vector for tuple τ on path P to flat pool ---
    // req[i] = max(0, j_c - nh_c) for each class i on the path.
    // LowMem: writes m int16_t values directly to pi.deadBoxesFlat.
    // Profile counters for DomPrune (paper Alg \textsc{DomPrune}).
    long long prof_dompruneIncomingSkipped = 0;
    long long prof_dompruneOldRemoved = 0;

    // Per-batch scratch — declared BEFORE appendReqBox so the lambda's
    // [&] capture sees the name. They are read/written by both
    // appendReqBox (for oldBoxSizeByPath) and the peel loop (for the
    // others). pathInfos is fully built by this point.
    //   batchDead[piIdx]: peeled tuples on path piIdx in current batch
    //   newAppendsByPath[piIdx]: count of dead boxes still in suffix
    //                           after DomPrune accepts/removals (= length
    //                           of new-box suffix in pi.deadBoxesFlat)
    //   oldBoxSizeByPath[piIdx]: pre-batch box count, used by
    //                            appendReqBox pass (ii) as the upper
    //                            bound of "removable old boxes" so
    //                            same-batch siblings stay in the suffix
    //   touchedPaths: unique piIdx with non-empty batchDead[piIdx]
    std::vector<std::vector<daf::Size>> batchDead(pathInfos.size());
    std::vector<int> newAppendsByPath(pathInfos.size(), 0);
    std::vector<int> oldBoxSizeByPath(pathInfos.size(), 0);
    std::vector<daf::Size> touchedPaths;
    touchedPaths.reserve(64);

    // Stack scratch for the candidate dead box (avoid per-call heap alloc).
    // MAX_M is the cap enforced earlier on |classIds| for any path. If a
    // path exceeds it, we'd silently truncate, so guard with a runtime
    // check (returns false → no append, which is conservatively wrong only
    // by including a duplicate-style box; we abort to make the bug loud).
    // Returns true if the candidate box was appended (DomPrune accepted it).
    // piIdx is needed so pass (ii) can restrict its old-box-removal scan to
    // [0, oldBoxSizeByPath[piIdx]) — boxes appended earlier in the SAME
    // batch must not be displaced, otherwise the new-box suffix invariant
    // breaks (and so does newAppendsByPath bookkeeping).
    auto appendReqBox = [&](daf::Size tidx, daf::Size piIdx, PathInfo &pi) -> bool {
        int m = (int)pi.classIds.size();
        // LowMem: stack buffer regression — earlier code was using a
        // thread_local heap to support paths with m > MAX_M. We restore
        // that contract here so this code path doesn't abort on dense
        // graphs (e.g. web-it-2004 leaves with m > MAX_M=512).
        thread_local std::vector<int16_t> candBuf;
        if ((int)candBuf.size() < m) candBuf.resize(m);
        int16_t *cand = candBuf.data();
        std::memset(cand, 0, (size_t)m * sizeof(int16_t));
        auto key = keyOf(tidx);  // sorted with repetitions
        // Merge-scan key (sorted) against pi.classIds (sorted) to avoid
        // the per-call unordered_map<> construct/destruct.
        int ci = 0, ki = 0;
        int klen = (int)key.size();
        while (ci < m && ki < klen) {
            if (pi.classIds[ci] < key[ki]) { ++ci; }
            else if (pi.classIds[ci] > key[ki]) { ++ki; }
            else {
                int jc = 1;
                while (ki + jc < klen && key[ki + jc] == key[ki]) ++jc;
                int v = jc - pi.nh[ci];
                if (v > 0) cand[ci] = (int16_t)v;
                ++ci; ki += jc;
            }
        }

        // Note: when cand = 0 (peeled tuple's class usage fits within
        // hold quotas), Theorem 5 says all alive τ on P have alive=0.
        // We rely on Theorem 1 in refreshAffectedPaths to detect this
        // via newDead = totalFeas (countUnion short-circuits via
        // fullCover when D=0 is in the box list). Doing the path retire
        // here would skip required per-τ dSup updates (which T1 does
        // correctly in refresh). Net savings of doing it here are
        // small (~200 ns/τ), not worth the duplication of update logic.
        // Rely on T1 + DomPrune Pass (ii) which already cleans up old
        // boxes when cand = 0.

        // === Bidirectional DomPrune (paper Alg \textsc{DomPrune}):
        //   (i) If some existing D dominates cand (D[c] ≤ cand[c] ∀c),
        //       cand's kill region ⊆ D's, so cand contributes nothing
        //       to the union → reject cand.
        //   (ii) Conversely, every existing D' dominated by cand
        //        (cand[c] ≤ D'[c] ∀c) is now redundant — its kill
        //        region ⊆ cand's. We sweep them out before committing
        //        cand. This keeps the antichain at the smallest size
        //        the paper guarantees. The per-call cost is
        //        O(existing × m), amortized over the peel.
        // On sparse graphs DomPrune shrinks deadBoxesFlat by 1-2
        // orders of magnitude (most peeled τ★ have jc(τ★) ≤ nh on
        // every class, hence cand = 0-vector → dominated by any
        // earlier accepted box ⇒ rejected after the first append).
        // ===
        const int existing = (int)(pi.deadBoxesFlat.size() / (size_t)m);
        const int oldEnd = oldBoxSizeByPath[piIdx];  // boxes 0..oldEnd are pre-batch
        // Defensive: oldEnd should never exceed existing.
        // (existing - oldEnd) is the count of boxes already appended this
        // batch; we must NOT touch them in pass (ii).
        // Pass (i) — rejection scan: scan ALL existing (old + this-batch)
        // because a same-batch sibling might already dominate cand.
        for (int b = 0; b < existing; ++b) {
            const int16_t *D = pi.deadBoxesFlat.data() + (size_t)b * m;
            bool dominates = true;
            for (int c = 0; c < m; ++c) {
                if (D[c] > cand[c]) { dominates = false; break; }
            }
            if (dominates) {
                prof_dompruneIncomingSkipped++;
                return false;
            }
        }
        // Pass (ii) — compaction: only OLD boxes (b < oldEnd) may be
        // removed if dominated by cand. Same-batch siblings (b >= oldEnd)
        // stay untouched, preserving the new-box suffix invariant. Old
        // boxes that survive are written back compactly; same-batch
        // siblings then shift LEFT by the number of removed old boxes,
        // but they remain a contiguous suffix.
        int writeIdx = 0;
        for (int b = 0; b < oldEnd; ++b) {
            const int16_t *D = pi.deadBoxesFlat.data() + (size_t)b * m;
            bool dominated = true;
            for (int c = 0; c < m; ++c) {
                if (cand[c] > D[c]) { dominated = false; break; }
            }
            if (dominated) {
                prof_dompruneOldRemoved++;
                continue;
            }
            if (writeIdx != b) {
                std::memcpy(pi.deadBoxesFlat.data() + (size_t)writeIdx * m,
                            pi.deadBoxesFlat.data() + (size_t)b * m,
                            (size_t)m * sizeof(int16_t));
            }
            writeIdx++;
        }
        if (writeIdx != oldEnd) {
            // Shift same-batch suffix [oldEnd, existing) left by the
            // number of removed old boxes (= oldEnd - writeIdx). Use
            // memmove because source and dest overlap when removed > 0
            // and the shift is < (existing - oldEnd) bytes.
            const int removed = oldEnd - writeIdx;
            const int siblingCount = existing - oldEnd;
            if (siblingCount > 0) {
                std::memmove(pi.deadBoxesFlat.data() + (size_t)writeIdx * m,
                             pi.deadBoxesFlat.data() + (size_t)oldEnd * m,
                             (size_t)siblingCount * (size_t)m * sizeof(int16_t));
            }
            // Update the snapshot so subsequent appendReqBox calls in this
            // batch see the correct oldEnd.
            oldBoxSizeByPath[piIdx] = oldEnd - removed;
            pi.deadBoxesFlat.resize((size_t)(existing - removed) * m);
        }

        // Commit: append cand to the flat box buffer.
        size_t off = pi.deadBoxesFlat.size();
        pi.deadBoxesFlat.resize(off + (size_t)m);
        std::memcpy(pi.deadBoxesFlat.data() + off, cand, (size_t)m * sizeof(int16_t));
        return true;
    };

    std::cout << "  MinCore floor: computed inline during Step 3" << std::endl;

    // --- Per-(tuple, path) saturation cache ---
    // Stores (newDead, totalFeas) per (P, τ). totalFeas is computed
    // once at first non-trivial refresh and immutable thereafter.
    // alive(τ, P) = totalFeas - newDead. When alive ≤ 0.5, τ is
    // saturated on P (Theorem 1) and can be removed from pi_P.tupleIdxs.
    //
    // Memory: 16 bytes per entry (two doubles) vs old 8 bytes — caller
    // pays ~2× cache memory but gains the saturation pruning.
    struct CacheVal {
        double newDead;
        double totalFeas;
    };
    robin_hood::unordered_flat_map<uint64_t, CacheVal> deadCache;

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
    long long prof_unionCalls = 0, prof_totalRecCalls = 0;
    long long prof_deadBoxesAdded = 0, prof_tupleUpdates = 0;
    long long prof_pathsRetired = 0;
    long long prof_deadCacheEvicted = 0;
    daf::Size numTuplesSz = rTuples.size();

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

    BigLevel maxSup = 0;
    std::vector<BigLevel> effSup(rTuples.size());
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        BigLevel sv = safeToBigLevel(dSup[i]);
        effSup[i] = std::max(sv, safeToBigLevel(tupleMinCore[i]));
        maxSup = std::max(maxSup, effSup[i]);
    }
    const BigLevel BUCKET_CAP = std::min<BigLevel>(maxSup + 2, (BigLevel)1000001);

    std::vector<std::vector<daf::Size>> buckets(BUCKET_CAP);
    std::multimap<BigLevel, daf::Size> overflow;
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        if (effSup[i] < BUCKET_CAP) buckets[effSup[i]].push_back(i);
        else overflow.insert({effSup[i], i});
    }

    // Profiling counters
    long long prof_batchCount = 0;

    // Hoisted per-path buffers reused across tuples to avoid per-tuple heap
    // allocations. Also replace the per-tuple unordered_map<cid,int> with a
    // merge-scan over the sorted key + sorted pi.classIds.
    std::vector<int> base, upper;
    base.reserve(MAX_M); upper.reserve(MAX_M);

    // batchDead / newAppendsByPath / oldBoxSizeByPath / touchedPaths
    // are declared earlier (see the appendReqBox lambda comment) so that
    // appendReqBox's [&] capture and refreshAffectedPaths can both read
    // them. Nothing additional to declare here.

    // PERF: hoist `affected` out of the per-batch hot loop. Previous
    // version did `std::vector<daf::Size> affected; affected.reserve(...);`
    // every batch — one heap alloc per batch. With thousands of batches
    // on dense graphs this added up to noticeable pure-alloc overhead.
    // Reused with .clear() (capacity preserved).
    std::vector<daf::Size> affected;
    affected.reserve(64);

    // Profile counter for the new affected-predicate short-circuit.
    long long prof_unionSkippedByAffectedCheck = 0;
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

    auto refreshAffectedPaths = [&](const std::vector<daf::Size> &affectedPaths,
                                    BigLevel &scanLevel) {
        // Per-path cache for the all-zero base[] case. Filled lazily on
        // the first all-zero tuple of each (path, batch); reused for
        // every subsequent all-zero tuple on the same path. Stores:
        //   - affected: precomputed affected-predicate result for base=0
        //   - fullCover / numBoxes: result of normalizeBoxes(base=0,...)
        //   - boxesNoPivot: normalized boxes with pivot already removed
        //   - pivotBox: the chosen pivot, expanded to int[]
        struct PathZeroCache {
            bool valid = false;
            bool affected = false;
            bool fullCover = false;
            bool emptyAfterNorm = false;     // norm.numBoxes == 0 → newDead = 0
            int numBoxesNoPivot = 0;         // norm.numBoxes - 1 (>= 0 when not empty)
            std::vector<int16_t> boxesNoPivot;
            std::vector<int> pivotBox;
        };
        PathZeroCache zeroCache;
        // Per-path saturation collector (reused across paths via clear()).
        std::vector<daf::Size> saturatedThisPath;

        for (daf::Size piIdx : affectedPaths) {
            auto &pi = pathInfos[piIdx];
            int m = (int)pi.classIds.size();
            // Reset cache for this path.
            zeroCache.valid = false;
            // LowMem: stack-buffer cap (MAX_M/MAX_T) no longer used on the
            // weight-table hot path (moved to thread_local heap), so no
            // ceiling to enforce here.

            // In-place filter of peeled tuples. After erase the size
            // shrinks but capacity stays, so over time pi.tupleIdxs holds
            // a heap allocation sized to its historical max even when only
            // a few tuples are alive on the path. PERF-AUDIT MEMORY (2/3):
            // when capacity grows to ≥4× the live size (geometric trigger,
            // amortised O(N) over the peel), shrink via copy-and-swap to
            // give the allocator the head room back. On graphs with ~10⁵
            // paths whose tuple counts decay sharply over the peel, this
            // reclaims hundreds of MB without adding to the hot path.
            pi.tupleIdxs.erase(
                std::remove_if(pi.tupleIdxs.begin(), pi.tupleIdxs.end(),
                               [&](daf::Size t) { return rPeeled[t]; }),
                pi.tupleIdxs.end());
            if (pi.tupleIdxs.capacity() > pi.tupleIdxs.size() * 4 + 16) {
                std::vector<daf::Size>(pi.tupleIdxs).swap(pi.tupleIdxs);
            }

            if (pi.deadBoxesFlat.empty()) continue;

            // LowMem: countUnionRec reads pi.deadBoxesFlat directly.
            int nBoxes = (int)(pi.deadBoxesFlat.size() / (size_t)m);

            // Boxes appended in this batch occupy the suffix of deadBoxesFlat.
            // newAppendsByPath[piIdx] is the *actual* number of boxes added
            // (DomPrune may have rejected some of the candidates), so the
            // suffix length equals newAppendsByPath[piIdx], not
            // batchDead[piIdx].size().
            // Earlier boxes are already accounted for in deadCache[piIdx,tidx]
            // (or, if there was no prior cache hit, contribute 0 to newDead
            // because they were unaffected by τ' at the time and dead-box
            // contributions are monotone non-negative per box).
            const int numNewBoxes = newAppendsByPath[piIdx];
            const int firstNewBoxOffset = (nBoxes - numNewBoxes) * m;
            const int T = pi.T;

            // upper is tuple-independent per path: fill once.
            upper.assign(m, 0);
            for (int i = 0; i < m; ++i) upper[i] = pi.np[i];

            for (auto tidx : pi.tupleIdxs) {
                auto key = keyOf(tidx);  // sorted with repetitions
                // Merge-scan: base[i] = max(0, j_c(τ) - nh_c). Classes not
                // present in key have j_c = 0 → base[i] = 0.
                base.assign(m, 0);
                int ci = 0, ki = 0;
                bool baseAllZero = true;
                while (ci < m && ki < (int)key.size()) {
                    if (pi.classIds[ci] < key[ki]) { ++ci; }
                    else if (pi.classIds[ci] > key[ki]) { ++ki; }
                    else {
                        int jc = 1;
                        while (ki + jc < (int)key.size() && key[ki + jc] == key[ki]) ++jc;
                        int b = std::max(0, jc - pi.nh[ci]);
                        base[ci] = b;
                        if (b != 0) baseAllZero = false;
                        ++ci; ki += jc;
                    }
                }

                double newDead;
                UnionCtx ctx{m, pi.T, 0, tidx, &pi};
                bool unionSkipped = false;

                if (baseAllZero) {
                    // === Path-level fast path (all-zero base) ===
                    // base=0 makes both the affected predicate and the top
                    // normalize+choosePivot path-level (no tuple input
                    // beyond ctx). Compute once per (path, batch); reuse
                    // for every all-zero tuple on this path.
                    if (!zeroCache.valid) {
                        // Affected predicate with base=0: max(D[c], 0) = D[c].
                        if (numNewBoxes > 0) {
                            bool affected = false;
                            const int16_t *Dptr = pi.deadBoxesFlat.data() + firstNewBoxOffset;
                            for (int b = 0; b < numNewBoxes; ++b) {
                                int s = 0;
                                bool over = false;
                                for (int c = 0; c < m; ++c) {
                                    int dv = (int)Dptr[c];
                                    s += dv;
                                    if (s > T) { over = true; break; }
                                }
                                if (!over) { affected = true; break; }
                                Dptr += m;
                            }
                            zeroCache.affected = affected;
                        } else {
                            zeroCache.affected = true;
                        }
                        if (zeroCache.affected) {
                            // Use thread_local scratch as outBoxes; copy
                            // the result into zeroCache (which is the
                            // long-lived per-(path,batch) cache).
                            thread_local std::vector<int16_t> cacheFillScratch;
                            auto norm = normalizeBoxes(base.data(), upper.data(), T,
                                                       pi.deadBoxesFlat.data(), nBoxes,
                                                       m, true, cacheFillScratch);
                            zeroCache.fullCover = norm.fullCover;
                            zeroCache.emptyAfterNorm = (!norm.fullCover && norm.numBoxes == 0);
                            if (!norm.fullCover && norm.numBoxes > 0) {
                                int pivIdx = choosePivot(base.data(),
                                                         cacheFillScratch.data(),
                                                         norm.numBoxes, m);
                                zeroCache.pivotBox.resize(m);
                                const int16_t *pv = cacheFillScratch.data()
                                                  + (size_t)pivIdx * m;
                                for (int i = 0; i < m; ++i)
                                    zeroCache.pivotBox[i] = (int)pv[i];
                                zeroCache.numBoxesNoPivot = norm.numBoxes - 1;
                                zeroCache.boxesNoPivot.assign(
                                    cacheFillScratch.begin(),
                                    cacheFillScratch.begin()
                                        + (size_t)zeroCache.numBoxesNoPivot * m);
                                if (pivIdx != zeroCache.numBoxesNoPivot) {
                                    int16_t *p = zeroCache.boxesNoPivot.data();
                                    int16_t *src = cacheFillScratch.data()
                                                 + (size_t)zeroCache.numBoxesNoPivot * m;
                                    int16_t *dst = p + (size_t)pivIdx * m;
                                    for (int i = 0; i < m; ++i) dst[i] = src[i];
                                }
                            }
                        }
                        zeroCache.valid = true;
                        prof_zeroBaseFills++;
                    } else {
                        prof_zeroBaseHits++;
                    }

                    if (!zeroCache.affected) {
                        prof_unionSkippedByAffectedCheck++;
                        unionSkipped = true;
                    }
                } else {
                    // === Standard (non-zero base) path ===
                    if (numNewBoxes > 0) {
                        bool affected = false;
                        const int16_t *Dptr = pi.deadBoxesFlat.data() + firstNewBoxOffset;
                        for (int b = 0; b < numNewBoxes; ++b) {
                            int s = 0;
                            bool over = false;
                            for (int c = 0; c < m; ++c) {
                                int dv = (int)Dptr[c];
                                int v = dv > base[c] ? dv : base[c];
                                s += v;
                                if (s > T) { over = true; break; }
                            }
                            if (!over) { affected = true; break; }
                            Dptr += m;
                        }
                        if (!affected) {
                            prof_unionSkippedByAffectedCheck++;
                            unionSkipped = true;
                        }
                    }
                }

                if (unionSkipped) continue;

                if (baseAllZero) {
                    // Apply cached norm+pivot. Recursion below the top
                    // frame still calls standard countUnionRec (it
                    // re-normalizes against new bounds, which differ
                    // from the cached zeros in deeper frames).
                    ctx.recCalls++;  // top frame
                    if (zeroCache.fullCover) {
                        newDead = feasWeighted(base.data(), upper.data(), ctx);
                    } else if (zeroCache.emptyAfterNorm) {
                        newDead = 0.0;
                    } else {
                        const int *pivot = zeroCache.pivotBox.data();
                        // Pivot included in the solution: feasWeighted
                        // over [pivot, upper].
                        newDead = feasWeighted(pivot, upper.data(), ctx);

                        // Recursion only contributes if there are
                        // non-pivot boxes. When numBoxesNoPivot == 0,
                        // the for-loop below still runs but the inner
                        // call gets numIn=0 and returns 0 immediately,
                        // matching standard behaviour.
                        int stackNL[MAX_M], stackNU[MAX_M];
                        std::vector<int> heapNL, heapNU;
                        int *nL = stackNL, *nU = stackNU;
                        if (m > MAX_M) {
                            heapNL.resize(m); heapNU.resize(m);
                            nL = heapNL.data(); nU = heapNU.data();
                        }
                        for (int splitDim = 0; splitDim < m; ++splitDim) {
                            // base[splitDim] = 0 (all-zero base).
                            if (pivot[splitDim] <= 0) continue;
                            for (int i = 0; i < m; ++i) {
                                nL[i] = base[i];      // 0
                                nU[i] = upper[i];
                            }
                            for (int earlier = 0; earlier < splitDim; ++earlier)
                                if (pivot[earlier] > nL[earlier])
                                    nL[earlier] = pivot[earlier];
                            nU[splitDim] = std::min(nU[splitDim],
                                                    pivot[splitDim] - 1);
                            newDead += countUnionRec(nL, nU,
                                                     zeroCache.boxesNoPivot.data(),
                                                     zeroCache.numBoxesNoPivot, ctx);
                        }
                    }
                } else {
                    newDead = countUnionRec(base.data(), upper.data(),
                                            pi.deadBoxesFlat.data(), nBoxes, ctx);
                }

                prof_totalRecCalls += ctx.recCalls;
                prof_unionCalls++;

                uint64_t cacheKey = (uint64_t)piIdx * numTuplesSz + tidx;
                double oldDead = 0.0;
                double totalFeas = -1.0;  // sentinel: needs computation
                auto dit = deadCache.find(cacheKey);
                if (dit != deadCache.end()) {
                    oldDead = dit->second.newDead;
                    totalFeas = dit->second.totalFeas;
                }

                double delta = newDead - oldDead;
                if (delta > 0.5) {
                    // First-time cache fill: compute totalFeas (path-level
                    // feasibility for τ, no dead-box restriction). This is
                    // the constant invariant used for saturation detection
                    // (Theorem 1: alive = totalFeas - newDead).
                    if (totalFeas < 0.0) {
                        totalFeas = feasWeighted(base.data(), upper.data(), ctx);
                    }

                    double newAlive = totalFeas - newDead;
                    if (newAlive < 0.5) {
                        // Theorem 1 saturation: τ has 0 alive contribution
                        // from P. Remove τ from pi_P.tupleIdxs (deferred
                        // to batch-end), drop cache entry. dSup still
                        // updated below — the support drop is real.
                        saturatedThisPath.push_back(tidx);
                        deadCache.erase(cacheKey);
                        prof_tupleSaturated++;
                    } else {
                        deadCache[cacheKey] = {newDead, totalFeas};
                    }

                    dSup[tidx] -= delta / rTuples[tidx].mult;
                    if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                    prof_tupleUpdates++;

                    BigLevel newSup = safeToBigLevel(dSup[tidx]);
                    BigLevel newBucket = std::max(newSup, safeToBigLevel(tupleMinCore[tidx]));
                    effSup[tidx] = newBucket;
                    if (newBucket < BUCKET_CAP) {
                        buckets[newBucket].push_back(tidx);
                        if (newBucket < scanLevel) scanLevel = newBucket;
                    } else {
                        overflow.insert({newBucket, tidx});
                    }
                }
            }

            // Batch-remove saturated tuples and check for path retirement.
            // Theorem 1 says they contribute 0 to dSup going forward; safe
            // to remove from tupleIdxs, possibly retiring the path.
            if (!saturatedThisPath.empty()) {
                std::sort(saturatedThisPath.begin(), saturatedThisPath.end());
                pi.tupleIdxs.erase(
                    std::remove_if(pi.tupleIdxs.begin(), pi.tupleIdxs.end(),
                                   [&](daf::Size t) {
                                       return std::binary_search(
                                           saturatedThisPath.begin(),
                                           saturatedThisPath.end(), t);
                                   }),
                    pi.tupleIdxs.end());
                pathAliveCount[piIdx] -= (daf::Size)saturatedThisPath.size();
                saturatedThisPath.clear();

                if (pathAliveCount[piIdx] == 0) {
                    std::vector<int16_t>().swap(pi.deadBoxesFlat);
                    std::vector<daf::Size>().swap(pi.tupleIdxs);
                    std::vector<daf::Size>().swap(pi.classIds);
                    std::vector<int>().swap(pi.nh);
                    std::vector<int>().swap(pi.np);
                    pathRetired[piIdx] = true;
                    prof_pathsRetired++;
                }
            }
        }
    };

    // --- Batch peeling loop ---
    daf::Size numPeeled = 0;
    BigLevel currentLevel = 0, coreLevel = 0;

    // batchDead / newAppendsByPath / touchedPaths are declared earlier so
    // that refreshAffectedPaths can read newAppendsByPath[piIdx] (suffix
    // length of newly-appended dead boxes after DomPrune).
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
        }

        // Clear last batch's scratch (capacity kept).
        for (auto pi : touchedPaths) {
            batchDead[pi].clear();
            newAppendsByPath[pi] = 0;
        }
        touchedPaths.clear();

        auto _t_pct0 = std::chrono::steady_clock::now();
        for (auto idx : batch) {
            pathsCoveringTuple(idx, pctPaths);
            for (auto piIdx : pctPaths) {
                if (batchDead[piIdx].empty()) touchedPaths.push_back(piIdx);
                batchDead[piIdx].push_back(idx);
                // Evict stale deadCache entry: once tuple idx is peeled, no
                // future refresh queries (piIdx, idx).
                uint64_t cacheKey = (uint64_t)piIdx * numTuplesSz + idx;
                if (deadCache.erase(cacheKey)) prof_deadCacheEvicted++;
            }
        }
        prof_pctPhaseNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - _t_pct0).count();

        // Affected-path list built in-place below (retired paths dropped).
        // `affected` is declared once outside the loop and reused via
        // .clear() to preserve capacity across batches.
        auto _t_app0 = std::chrono::steady_clock::now();
        affected.clear();
        for (auto piIdx : touchedPaths) {
            auto &pi = pathInfos[piIdx];
            auto &deadTuples = batchDead[piIdx];
            pathAliveCount[piIdx] -= (daf::Size)deadTuples.size();
            if (pathAliveCount[piIdx] == 0) {
                // Path retired: no surviving tuple consults this path
                // for refresh. pathsCoveringTuple may still surface this
                // path from classToPaths / tupleToPathInfos when peeling
                // future tuples; that path now uses pathRetired[] as a
                // pre-filter so it never reads the freed nh/np.
                std::vector<int16_t>().swap(pi.deadBoxesFlat);
                std::vector<daf::Size>().swap(pi.tupleIdxs);
                std::vector<daf::Size>().swap(pi.classIds);
                std::vector<int>().swap(pi.nh);
                std::vector<int>().swap(pi.np);
                pathRetired[piIdx] = true;
                prof_pathsRetired++;
                continue;
            }
            int m = (int)pi.classIds.size();
            pi.deadBoxesFlat.reserve(pi.deadBoxesFlat.size() + deadTuples.size() * m);
            // Snapshot pre-batch box count so appendReqBox's pass (ii)
            // can restrict its old-box-removal scan to [0, oldEnd) and
            // leave any same-batch siblings in the suffix.
            const int oldBoxCountBefore = (int)(pi.deadBoxesFlat.size() / (size_t)m);
            oldBoxSizeByPath[piIdx] = oldBoxCountBefore;
            int acceptedThisBatch = 0;
            for (auto idx : deadTuples) {
                if (appendReqBox(idx, piIdx, pi)) {
                    prof_deadBoxesAdded++;
                    acceptedThisBatch++;
                }
            }
            const int curBoxCount = (int)(pi.deadBoxesFlat.size() / (size_t)m);
            // If DomPrune rejected every candidate AND removed nothing
            // from the existing buffer, the union value is unchanged →
            // skip refresh entirely on this path. (Pass-ii removals
            // would shrink curBoxCount relative to oldBoxCountBefore;
            // we re-read curBoxCount post-batch to detect them.)
            if (acceptedThisBatch == 0 && curBoxCount == oldBoxCountBefore) continue;
            // Precise new-suffix length: after the batch, the new boxes
            // form a contiguous suffix of length acceptedThisBatch
            // (pass ii only displaces them by removing some old boxes,
            // which keeps them contiguous at the end).
            newAppendsByPath[piIdx] = acceptedThisBatch;
            // PERF-AUDIT MEMORY (5): pi.deadBoxesFlat capacity grows
            // monotonically during peel via reserve()+resize(). When
            // bidirectional DomPrune (pass ii) shrinks the buffer below
            // half its capacity, the freed tail is kept reserved by
            // std::vector. Periodically shrink-to-fit when capacity ≥
            // 4× current size (geometric trigger; amortised O(N) total).
            const size_t cur = pi.deadBoxesFlat.size();
            if (pi.deadBoxesFlat.capacity() > cur * 4 + (size_t)m) {
                std::vector<int16_t>(pi.deadBoxesFlat).swap(pi.deadBoxesFlat);
            }
            affected.push_back(piIdx);
        }
        prof_appendPhaseNs += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - _t_app0).count();

        auto _t_ref0 = std::chrono::steady_clock::now();
        refreshAffectedPaths(affected, currentLevel);
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
    std::cout << "  Dead boxes added: " << prof_deadBoxesAdded << std::endl;
    std::cout << "  DomPrune incoming skipped: "
              << prof_dompruneIncomingSkipped << std::endl;
    std::cout << "  DomPrune old removed: "
              << prof_dompruneOldRemoved << std::endl;
    std::cout << "  Paths retired (alive=0): " << prof_pathsRetired
              << " / " << pathInfos.size() << std::endl;
    std::cout << "  DeadCache evictions (on peel): " << prof_deadCacheEvicted << std::endl;
    std::cout << "  DeadCache final size: " << deadCache.size() << std::endl;
    std::cout << "  Union calls: " << prof_unionCalls << std::endl;
    std::cout << "  Union skipped (unaffected by new boxes): "
              << prof_unionSkippedByAffectedCheck << std::endl;
    std::cout << "  Zero-base cache: " << prof_zeroBaseHits << " hits / "
              << prof_zeroBaseFills << " fills" << std::endl;
    std::cout << "  Total recursive calls: " << prof_totalRecCalls << std::endl;
    {
        double peelMs = (double)step6Ms;
        auto pct = [peelMs](long long ns) -> double {
            return peelMs > 0 ? (double)ns / 1e6 / peelMs * 100.0 : 0.0;
        };
        std::cout << "  --- Peel phase breakdown (outer) ---" << std::endl;
        std::cout << "    pctPhase  : " << prof_pctPhaseNs / 1000000 << " ms ("
                  << pct(prof_pctPhaseNs) << "%)" << std::endl;
        std::cout << "    appendPhz : " << prof_appendPhaseNs / 1000000 << " ms ("
                  << pct(prof_appendPhaseNs) << "%)" << std::endl;
        std::cout << "    refreshTot: " << prof_refreshTotalNs / 1000000 << " ms ("
                  << pct(prof_refreshTotalNs) << "%)" << std::endl;
        long long otherNs = (long long)peelMs * 1000000
                          - prof_pctPhaseNs - prof_appendPhaseNs - prof_refreshTotalNs;
        std::cout << "    other     : " << otherNs / 1000000 << " ms ("
                  << pct(otherNs) << "%) [bucket-pop, scratch-clear, cleanup]" << std::endl;
    }
    std::cout << "  Tuple updates: " << prof_tupleUpdates << std::endl;
    std::cout << "  Tuple saturated on path (Theorem 1): " << prof_tupleSaturated << std::endl;
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

    // Return compact format: one entry per core level, key[0] = count
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    for (auto &[c, cnt] : coreDist)
        // Compact format: key = {lo32, hi32} encoding int64_t count
        result.push_back({{(daf::Size)(cnt & 0xFFFFFFFF), (daf::Size)((cnt >> 32) & 0xFFFFFFFF)}, (double)c});
    return result;
}
