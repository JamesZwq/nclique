//
// Region Tuple + Direct Enumeration (V3Fast_NoCPI — ABLATION VARIANT)
//
// Identical to V3Fast in every phase EXCEPT Step 4 (initial support
// computation), which replaces the CPI convolution formula with direct
// enumeration of s-cliques.
//
// Purpose: ablation study for the Clique Path Index (CPI) support formula
// (Theorem 6.1). V3Fast_NoCPI keeps all other optimizations (fully-mergeable
// region isolation, tuple compression, private cloud, dead-box peeling) so
// the speed gap between V3Fast and V3Fast_NoCPI measures the CPI contribution
// in isolation.
//
// Step 4 replacement (the only algorithmic change):
//   For each non-fully-mergeable region M, enumerate every s-subset of M
//   (which is an s-clique since M is a clique), deduplicate by canonical
//   home (min-index region containing the subset), skip s-cliques using
//   any private vertex in Private Cloud mode, then enumerate r-subsets and
//   attribute +1 to the corresponding tuple's support.
//
// Complexity of Step 4:
//   V3Fast (with CPI):   O(|tuples| × |paths| × r²)
//   V3Fast_NoCPI:        O(Σ_M C(|M|, s) × C(s, r))
//   The no-CPI variant is expected to time out when max MC size is large
//   and s is moderate (e.g. com-dblp max MC = 114 + s = 20 → C(114, 20)
//   ≈ 10¹⁹ s-subsets per region). That timeout is the point of the ablation.
//
// Peeling (Step 6) is unchanged: it uses PathInfo and dead-box union
// measure, orthogonal to how initial support is computed.
//

#include "NCliqueCoreDecomposition.h"
#include "../dataStruct/robin_hood.h"
#include <algorithm>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

extern double nCr[1001][401];
extern std::vector<bool> g_maxCliqueTags;
#include "FlatCliques.h"
extern FlatCliques g_maxCliques;

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
// Main function: Region CPI (V3Fast_NoCPI ablation)
// ============================================================

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionCPI_NoCPI(
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

    for (auto mc : g_maxCliques) {
        if ((int)mc.size() < s) continue;
        daf::Size rid = regionVerts.size();
        regionVerts.emplace_back(mc.begin(), mc.end()); // already sorted by MaxCliqEnum
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
    {
        daf::Size maxVtxMCs = 0;
        double totalPairs = 0;
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

    {
        auto t = std::chrono::high_resolution_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t - tStep1b).count();
        std::cout << "  Step 1b fused-FM: " << fmTotalIter << " counter ops, "
                  << fmExits << " early exits, " << ms << " ms" << std::endl;
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

    auto isActiveTupleClass = [&](daf::Size cid) -> bool {
        return !enablePrivateCloud || !isPrivateClass[cid];
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
        if (!isPrivateClass[cid]) continue;
        daf::Size rid = privateClassMC[cid];
        if (rid == INVALID) continue;
        privateVertexCount[rid] += classSizes[cid];
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
    // mult = Π_c C(|c|, j_c) can overflow uint32 for r≥5 on graphs with
    // large class sizes (e.g. r=11 on com-dblp: C(~50, 11) ≈ 4e10). Widen
    // to uint64 so the multiplicity is represented correctly throughout
    // peeling and core aggregation.
    struct RTuple { TupleKey key; uint64_t mult; };
    std::vector<RTuple> rTuples;
    std::vector<double> tupleMinCore; // per-tuple minCore floor

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
                rTuples.push_back({key, mult});
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
    std::cout << "  Active r-tuples: " << rTuples.size() << std::endl;
    if (enablePrivateCloud) {
        std::cout << "  Private clouds: " << numPrivateClouds << std::endl;
        std::cout << "  Private-touching r-cliques (direct): "
                  << std::fixed << std::setprecision(0) << privateRCliquesDirect << std::endl;
    } else {
        std::cout << "  Private cloud mode: disabled" << std::endl;
    }

    // ============================================================
    // Step 4 (NoCPI): Direct s-clique enumeration for initial support
    // ============================================================
    // For each non-fully-mergeable region M, enumerate every s-subset of M
    // (each is an s-clique since M is a clique). Deduplicate by canonical
    // home (min-index region containing the subset) so each s-clique is
    // counted exactly once. In Private Cloud mode, skip s-cliques that use
    // any private vertex (they're already accounted for by the Step 3
    // direct assignment). Finally, enumerate r-subsets of each surviving
    // s-clique, map each to its tuple via the class-id composition, and
    // accumulate +1 into that tuple's support.
    //
    // Expected to be VERY slow / time out on graphs where max MC size is
    // large relative to s: C(|M|, s) grows exponentially in min(|M|, s).
    // That slowness is the point of this ablation — it quantifies how much
    // the CPI convolution formula (Theorem 6.1) saves.

    std::vector<double> support(rTuples.size(), 0.0);

    // Enumerator for k-subsets of a sorted vector (indices form combinations).
    auto enumerateSubsets = [&](const std::vector<daf::Size> &src, int k,
                                const std::function<void(const daf::Size *)> &cb) {
        int n = (int)src.size();
        if (k > n || k <= 0) return;
        std::vector<int> idx(k);
        for (int i = 0; i < k; ++i) idx[i] = i;
        std::vector<daf::Size> buf(k);
        while (true) {
            for (int i = 0; i < k; ++i) buf[i] = src[idx[i]];
            cb(buf.data());
            int i = k - 1;
            while (i >= 0 && idx[i] == n - k + i) --i;
            if (i < 0) break;
            ++idx[i];
            for (int j = i + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
        }
    };

    // isSubsetOfRegion: test whether sorted s_vec ⊆ sorted regionVerts[rid].
    auto isSubsetOfRegion = [&](const daf::Size *s_vec, int s_len, daf::Size rid2) -> bool {
        auto &r2 = regionVerts[rid2];
        size_t j = 0; int i = 0;
        while (i < s_len && j < r2.size()) {
            if (s_vec[i] < r2[j]) return false;
            if (s_vec[i] == r2[j]) { ++i; ++j; }
            else ++j;
        }
        return i == s_len;
    };

    // Canonical home check: rid is canonical for s_vec iff no rid2 < rid has s_vec ⊆ M_{rid2}.
    // Sufficient to look at vtxMaxPaths[s_vec[0]] (rid2 < rid, and s_vec ⊆ M_{rid2} → s_vec[0] ∈ M_{rid2}).
    auto isCanonicalHome = [&](const daf::Size *s_vec, int s_len, daf::Size rid) -> bool {
        daf::Size v0 = s_vec[0];
        for (daf::Size rid2 : vtxMaxPaths[v0]) {
            if (rid2 >= rid) break; // vtxMaxPaths is in push-order which equals region-id order
            if (isSubsetOfRegion(s_vec, s_len, rid2)) return false;
        }
        return true;
    };

    auto tStep4 = std::chrono::high_resolution_clock::now();
    long long prof_sEnumerated = 0, prof_sCanonical = 0, prof_sNonPrivate = 0;
    long long prof_rMatched = 0;

    TupleKey keyBuf; keyBuf.reserve(r);

    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        auto &M = regionVerts[rid];
        int mSize = (int)M.size();
        if (mSize < (int)s) continue;

        // Enumerate s-subsets of M (M is sorted).
        enumerateSubsets(M, (int)s, [&](const daf::Size *sVec) {
            ++prof_sEnumerated;

            // Dedup: only process if rid is the canonical (smallest-index) home.
            if (!isCanonicalHome(sVec, (int)s, rid)) return;
            ++prof_sCanonical;

            // Private-cloud filter: skip any s-clique that uses a private vertex
            // (those are accounted for by the Step 3 direct assignment block).
            if (enablePrivateCloud) {
                for (int i = 0; i < (int)s; ++i) {
                    daf::Size cid = classOf[sVec[i]];
                    if (cid != INVALID && isPrivateClass[cid]) return;
                }
            }
            ++prof_sNonPrivate;

            // Enumerate r-subsets of sVec. For each, build its tuple key
            // (sorted class-id multiset). If the tuple is registered (in
            // rTupleIndex), credit +1 to its support / mult.
            std::vector<daf::Size> sVecV(sVec, sVec + (int)s);
            enumerateSubsets(sVecV, (int)r, [&](const daf::Size *rVec) {
                keyBuf.clear();
                bool invalid = false;
                for (int i = 0; i < (int)r; ++i) {
                    daf::Size cid = classOf[rVec[i]];
                    if (cid == INVALID) { invalid = true; break; }
                    keyBuf.push_back(cid);
                }
                if (invalid) return;
                std::sort(keyBuf.begin(), keyBuf.end());
                auto it = rTupleIndex.find(keyBuf);
                if (it == rTupleIndex.end()) return;
                daf::Size tidx = it->second;
                // Each s-clique contributes 1 s-clique-containing-R event for
                // every R-realization inside it; dividing by mult gives the
                // per-realization support (same convention as CPI's formula).
                support[tidx] += 1.0 / (double)rTuples[tidx].mult;
                ++prof_rMatched;
            });
        });
    }

    auto tStep4b = std::chrono::high_resolution_clock::now();
    auto step4Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep4b - tStep4).count();

    double totalSupportTuples = 0, totalRCliques = 0;
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        totalSupportTuples += rTuples[i].mult * support[i];
        totalRCliques += rTuples[i].mult;
    }
    double totalRCliquesWithPrivate = totalRCliques + privateRCliquesDirect;
    std::cout << "  [NoCPI] s-subsets enumerated: " << prof_sEnumerated << std::endl;
    std::cout << "  [NoCPI] s-subsets canonical:  " << prof_sCanonical << std::endl;
    std::cout << "  [NoCPI] s-subsets non-private:" << prof_sNonPrivate << std::endl;
    std::cout << "  [NoCPI] r-subset matches:     " << prof_rMatched << std::endl;
    std::cout << "  Active r-cliques: " << std::fixed << std::setprecision(0) << totalRCliques << std::endl;
    if (enablePrivateCloud) {
        std::cout << "  Total r-cliques (with private cloud): "
                  << std::fixed << std::setprecision(0) << totalRCliquesWithPrivate << std::endl;
    }
    std::cout << "  Support sum (NoCPI enum): " << totalSupportTuples << std::endl;
    std::cout << "  NoCPI enumeration time: " << step4Ms << " ms" << std::endl;

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
        auto &key = rTuples[tidx].key;
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
        std::vector<std::vector<int>> deadBoxes;  // accumulated dead requirement vectors
        // classId -> index in classIds (for fast lookup)
        std::unordered_map<daf::Size, int> classToIdx;
    };

    // Compile-time bounds for stack-allocated hot-path buffers. m (classes per
    // path) is small (typ. 10-30). T (pivot target) and per-class range can
    // grow with s — for large s on graphs with big MCs, range can hit 100+.
    // Bumped bounds cover ca-HepPh (max MC 239), web-it-2004 (max MC 432).
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
        if (target >= MAX_T + 1) { std::cerr << "countFeasibleWeighted: target=" << target << " exceeds MAX_T\n"; std::abort(); }

        double dp[MAX_T + 1], next[MAX_T + 1];
        for (int t = 0; t <= target; ++t) dp[t] = 0.0;
        dp[0] = 1.0;
        for (int i = 0; i < m; ++i) {
            int cap = upper[i] - lower[i];
            const double *wt = wtsFlat + (size_t)i * wtStride;
            int len = wtLen[i];
            for (int t = 0; t <= target; ++t) next[t] = 0.0;
            for (int t = 0; t <= target; ++t) {
                double sum = 0.0;
                int kMax = std::min(cap, t);
                if (kMax >= len) kMax = len - 1;
                for (int k = 0; k <= kMax; ++k)
                    sum += dp[t - k] * wt[k];
                next[t] = sum;
            }
            std::swap_ranges(dp, dp + target + 1, next);
        }
        return dp[target];
    };

    // --- Build weight tables for tuple τ' on path P (stack buffers) ---
    // Writes into wtsFlat (expected size MAX_M × wtStride). wtLen[i] = range per class.
    // Uses merge-scan on sorted key + sorted pi.classIds to avoid an unordered_map.
    auto buildWeightTables_sb = [&](daf::Size tidx, const PathInfo &pi,
                                    const int *lower, const int *upper,
                                    double *wtsFlat, int wtStride, int *wtLen) {
        int m = (int)pi.classIds.size();
        auto &key = rTuples[tidx].key;   // sorted with repetitions
        // Merge-scan: jvec[i] = j_c(τ) for class pi.classIds[i]
        int jvec[MAX_M];
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
            if (range > wtStride) { std::cerr << "buildWeightTables: range > wtStride\n"; std::abort(); }
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

    // --- normalizeBoxes ---
    struct NormResult { bool fullCover; std::vector<std::vector<int>> boxes; };
    auto normalizeBoxes = [](const std::vector<int> &lower,
                             const std::vector<int> &upper, int T,
                             const std::vector<std::vector<int>> &boxes,
                             bool pruneDom) -> NormResult {
        int m = (int)lower.size();
        std::vector<std::vector<int>> effective;
        effective.reserve(boxes.size());
        for (auto &box : boxes) {
            std::vector<int> cur(m);
            bool impossible = false;
            int lSum = 0;
            for (int i = 0; i < m; ++i) {
                cur[i] = std::max(lower[i], box[i]);
                if (cur[i] > upper[i]) { impossible = true; break; }
                lSum += cur[i];
            }
            if (impossible || lSum > T) continue;
            if (cur == lower) return {true, {}};
            effective.push_back(std::move(cur));
        }
        std::sort(effective.begin(), effective.end());
        effective.erase(std::unique(effective.begin(), effective.end()), effective.end());
        if (!pruneDom) return {false, effective};

        // Sort by sum ascending for dominance pruning
        std::sort(effective.begin(), effective.end(), [](const std::vector<int> &a, const std::vector<int> &b) {
            int sa = 0, sb = 0;
            for (int x : a) sa += x;
            for (int x : b) sb += x;
            return sa != sb ? sa < sb : a < b;
        });
        std::vector<std::vector<int>> minimal;
        for (auto &box : effective) {
            bool dominated = false;
            for (auto &kept : minimal) {
                bool dom = true;
                for (int i = 0; i < m; ++i)
                    if (kept[i] > box[i]) { dom = false; break; }
                if (dom) { dominated = true; break; }
            }
            if (!dominated) minimal.push_back(box);
        }
        return {false, minimal};
    };

    // --- countUnionWeighted: branch-and-bound ---
    // choosePivotBox: pick box with fewest active dims (ties: highest sum)
    auto choosePivot = [](const std::vector<int> &lower,
                          const std::vector<std::vector<int>> &boxes, int m) -> int {
        int bestIdx = 0, bestActive = m + 1, bestSum = -1;
        for (int i = 0; i < (int)boxes.size(); ++i) {
            int active = 0, lsum = 0;
            for (int c = 0; c < m; ++c) {
                if (boxes[i][c] > lower[c]) ++active;
                lsum += boxes[i][c];
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

    // Forward-declare countUnionRec
    std::function<double(const std::vector<int>&, const std::vector<int>&,
                         const std::vector<std::vector<int>>&, UnionCtx&)> countUnionRec;

    // Helper: compute weighted feasible for given lower/upper using tuple weights
    // Stack-allocated wts buffer; one call per feasWeighted invocation.
    auto feasWeighted = [&](const std::vector<int> &lower, const std::vector<int> &upper,
                            UnionCtx &ctx) -> double {
        int m = (int)ctx.pi->classIds.size();
        // wtStride is the max "range" per class = max(upper[i] - lower[i] + 1).
        // Since lower[i] ≥ 0 and upper[i] ≤ np_c ≤ s, bounding by MAX_T+1 is safe.
        constexpr int wtStride = MAX_T + 1;
        double wtsFlat[MAX_M * wtStride];
        int wtLen[MAX_M];
        buildWeightTables_sb(ctx.tidx, *ctx.pi, lower.data(), upper.data(), wtsFlat, wtStride, wtLen);
        return countFeasibleWeighted_sb(lower.data(), upper.data(), m, ctx.T, wtsFlat, wtStride, wtLen);
    };

    countUnionRec = [&](const std::vector<int> &lower, const std::vector<int> &upper,
                        const std::vector<std::vector<int>> &boxes,
                        UnionCtx &ctx) -> double {
        ctx.recCalls++;
        // Feasibility check
        {
            int minS = 0, maxS = 0;
            for (int i = 0; i < ctx.m; ++i) {
                if (lower[i] > upper[i]) return 0.0;
                minS += lower[i]; maxS += upper[i];
            }
            if (ctx.T < minS || ctx.T > maxS) return 0.0;
        }

        auto norm = normalizeBoxes(lower, upper, ctx.T, boxes, true);
        if (norm.fullCover) return feasWeighted(lower, upper, ctx);
        if (norm.boxes.empty()) return 0.0;

        // Branch on pivot
        int pivIdx = choosePivot(lower, norm.boxes, ctx.m);
        std::vector<int> pivot = norm.boxes[pivIdx];

        double total = feasWeighted(pivot, upper, ctx);

        std::vector<std::vector<int>> remaining = norm.boxes;
        remaining.erase(remaining.begin() + pivIdx);

        for (int splitDim = 0; splitDim < ctx.m; ++splitDim) {
            if (pivot[splitDim] <= lower[splitDim]) continue;

            std::vector<int> nextLower = lower;
            std::vector<int> nextUpper = upper;
            for (int earlier = 0; earlier < splitDim; ++earlier)
                nextLower[earlier] = std::max(nextLower[earlier], pivot[earlier]);
            nextUpper[splitDim] = std::min(nextUpper[splitDim], pivot[splitDim] - 1);

            total += countUnionRec(nextLower, nextUpper, remaining, ctx);
        }
        return total;
    };

    // --- Build PathInfo structures ---
    std::vector<PathInfo> pathInfos;
    std::vector<std::vector<daf::Size>> tupleToPathInfos(rTuples.size());

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        PathInfo pi;
        pi.h = 0;
        bool piHasPrivateHold = false;
        std::unordered_map<daf::Size, std::pair<int,int>> cd; // cid -> (nh, np)
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            if (enablePrivateCloud && isPrivateClass[cid]) {
                if (!node.isPivot) piHasPrivateHold = true;
                continue;
            }
            if (node.isPivot) cd[cid].second++;
            else { cd[cid].first++; pi.h++; }
        }
        if (piHasPrivateHold) continue;
        pi.T = s - pi.h;
        if (pi.T < 0) continue;

        // Build ordered class arrays
        for (auto &[cid, p] : cd) pi.classIds.push_back(cid);
        std::sort(pi.classIds.begin(), pi.classIds.end());
        pi.nh.resize(pi.classIds.size());
        pi.np.resize(pi.classIds.size());
        for (int i = 0; i < (int)pi.classIds.size(); ++i) {
            pi.classToIdx[pi.classIds[i]] = i;
            pi.nh[i] = cd[pi.classIds[i]].first;
            pi.np[i] = cd[pi.classIds[i]].second;
        }

        // Find tuples on this path
        TupleKey cur; cur.reserve(r);
        std::vector<daf::Size> pathClasses;
        for (auto cid : pi.classIds) pathClasses.push_back(cid);
        bool hasTuples = false;
        daf::Size piIdx = pathInfos.size();

        std::function<void()> enumCb = [&]() {
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;
            // Feasibility check
            std::unordered_map<daf::Size, int> counts;
            for (auto c : cur) counts[c]++;
            for (auto &[c, jc] : counts) {
                auto cit = pi.classToIdx.find(c);
                if (cit == pi.classToIdx.end()) return;
                int idx2 = cit->second;
                if (jc > pi.nh[idx2] + pi.np[idx2]) return;
            }
            pi.tupleIdxs.push_back(tidx);
            hasTuples = true;
        };
        enumerateMultisets(pathClasses, r, 0, cur, enumCb);

        if (hasTuples) {
            pathInfos.push_back(std::move(pi));
            for (auto tidx : pathInfos.back().tupleIdxs)
                tupleToPathInfos[tidx].push_back(piIdx);
        }
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
                auto &key = rTuples[tidx].key;
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
                constexpr int wtStride = MAX_T + 1;
                double wtsFlat[MAX_M * wtStride];
                int wtLen[MAX_M];
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

    // --- Build requirement vector for tuple τ on path P ---
    // reqVec[i] = max(0, j_c - nh_c) for each class i on the path
    auto buildReqVec = [&](daf::Size tidx, const PathInfo &pi) -> std::vector<int> {
        int m = (int)pi.classIds.size();
        std::vector<int> req(m, 0);
        auto &key = rTuples[tidx].key;
        std::unordered_map<daf::Size, int> counts;
        for (auto c : key) counts[c]++;
        for (auto &[c, jc] : counts) {
            auto it = pi.classToIdx.find(c);
            if (it != pi.classToIdx.end()) {
                int idx2 = it->second;
                req[idx2] = std::max(0, jc - pi.nh[idx2]);
            }
        }
        return req;
    };

    std::cout << "  MinCore floor: computed inline during Step 3" << std::endl;

    // --- Per-(tuple, path) dead count cache ---
    // robin_hood flat map: ~2-3x faster than std::unordered_map for uint64 keys
    // on the peeling hot path, where every alive-tuple update queries the cache.
    robin_hood::unordered_flat_map<uint64_t, double> deadCache;

    // --- Bucket queue setup ---
    std::vector<double> dSup = support;
    std::vector<bool> rPeeled(rTuples.size(), false);

    // Profiling counters
    long long prof_unionCalls = 0, prof_totalRecCalls = 0;
    long long prof_deadBoxesAdded = 0, prof_tupleUpdates = 0;
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

    auto refreshAffectedPaths = [&](const std::unordered_set<daf::Size> &affectedPathSet,
                                    BigLevel &scanLevel) {
        for (daf::Size piIdx : affectedPathSet) {
            auto &pi = pathInfos[piIdx];
            int m = (int)pi.classIds.size();
            if (m > MAX_M) { std::cerr << "refreshAffectedPaths: m > MAX_M\n"; std::abort(); }

            std::vector<daf::Size> alive;
            alive.reserve(pi.tupleIdxs.size());
            for (auto tidx : pi.tupleIdxs) {
                if (!rPeeled[tidx])
                    alive.push_back(tidx);
            }
            pi.tupleIdxs = std::move(alive);

            if (pi.deadBoxes.empty()) continue;

            // upper is tuple-independent per path: fill once.
            upper.assign(m, 0);
            for (int i = 0; i < m; ++i) upper[i] = pi.np[i];

            for (auto tidx : pi.tupleIdxs) {
                auto &key = rTuples[tidx].key;  // sorted with repetitions
                // Merge-scan: base[i] = max(0, j_c(τ) - nh_c). Classes not
                // present in key have j_c = 0 → base[i] = 0.
                base.assign(m, 0);
                int ci = 0, ki = 0;
                while (ci < m && ki < (int)key.size()) {
                    if (pi.classIds[ci] < key[ki]) { ++ci; }
                    else if (pi.classIds[ci] > key[ki]) { ++ki; }
                    else {
                        int jc = 1;
                        while (ki + jc < (int)key.size() && key[ki + jc] == key[ki]) ++jc;
                        base[ci] = std::max(0, jc - pi.nh[ci]);
                        ++ci; ki += jc;
                    }
                }

                UnionCtx ctx{m, pi.T, 0, tidx, &pi};
                double newDead = countUnionRec(base, upper, pi.deadBoxes, ctx);
                prof_totalRecCalls += ctx.recCalls;
                prof_unionCalls++;

                uint64_t cacheKey = (uint64_t)piIdx * numTuplesSz + tidx;
                double oldDead = 0.0;
                auto dit = deadCache.find(cacheKey);
                if (dit != deadCache.end()) oldDead = dit->second;

                double delta = newDead - oldDead;
                if (delta <= 0.5) continue;

                deadCache[cacheKey] = newDead;
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
    };

    // --- Batch peeling loop ---
    daf::Size numPeeled = 0;
    BigLevel currentLevel = 0, coreLevel = 0;

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

        std::unordered_map<daf::Size, std::vector<daf::Size>> newlyDeadByPath;
        for (auto idx : batch) {
            for (auto piIdx : tupleToPathInfos[idx])
                newlyDeadByPath[piIdx].push_back(idx);
        }

        std::unordered_set<daf::Size> affectedPathSet;
        for (auto &[piIdx, deadTuples] : newlyDeadByPath) {
            auto &pi = pathInfos[piIdx];
            for (auto idx : deadTuples) {
                pi.deadBoxes.push_back(buildReqVec(idx, pi));
                prof_deadBoxesAdded++;
            }
            affectedPathSet.insert(piIdx);
        }
        refreshAffectedPaths(affectedPathSet, currentLevel);
    }

    auto tStep6End = std::chrono::high_resolution_clock::now();
    auto step6Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStep6).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStart).count();

    double maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling (Batch Union) ---" << std::endl;
    std::cout << "  Peeled active tuples: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Batches: " << prof_batchCount << std::endl;
    std::cout << "  Dead boxes added: " << prof_deadBoxesAdded << std::endl;
    std::cout << "  Union calls: " << prof_unionCalls << std::endl;
    std::cout << "  Total recursive calls: " << prof_totalRecCalls << std::endl;
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

    // Return compact format: one entry per core level, key[0] = count
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    for (auto &[c, cnt] : coreDist)
        // Compact format: key = {lo32, hi32} encoding int64_t count
        result.push_back({{(daf::Size)(cnt & 0xFFFFFFFF), (daf::Size)((cnt >> 32) & 0xFFFFFFFF)}, (double)c});
    return result;
}
