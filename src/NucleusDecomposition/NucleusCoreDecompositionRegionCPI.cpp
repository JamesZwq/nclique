//
// Region Tuple + CPI Counting (V3)
//
// Combines two orthogonal compressions:
// - Region Tuple: avoids r-clique enumeration (groups into tuples)
// - CPI Counting: avoids s-tuple enumeration (counts via path formula)
//
// Init: support(τ) = [Σ_P N(τ,P) × C(p-r+h, s-r)] / mult(τ)
// Peel: lazy s-tuple generation (only enumerate dying s-tuples during peeling)
//

#include "NCliqueCoreDecomposition.h"
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
NucleusCoreDecompositionRClique_RegionCPI(
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

    auto tStep1b = std::chrono::high_resolution_clock::now();

    // Compute pairwise intersection sizes (via vertex membership)
    std::unordered_map<uint64_t, int> interSize;
    auto pairKey = [](daf::Size a, daf::Size b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return ((uint64_t)a << 32) | b;
    };
    for (daf::Size v = 0; v < numVertices; ++v) {
        auto &rids = vtxMaxPaths[v];
        for (daf::Size i = 0; i < rids.size(); ++i)
            for (daf::Size j = i + 1; j < rids.size(); ++j) {
                auto key = pairKey(rids[i], rids[j]);
                if (!interSize.count(key)) {
                    auto &a = regionVerts[rids[i]], &b = regionVerts[rids[j]];
                    int cnt = 0;
                    auto ai = a.begin(), bi = b.begin();
                    while (ai != a.end() && bi != b.end()) {
                        if (*ai < *bi) ++ai;
                        else if (*ai > *bi) ++bi;
                        else { ++cnt; ++ai; ++bi; }
                    }
                    interSize[key] = cnt;
                }
            }
    }

    // Check each region: fully r-mergeable?
    std::vector<bool> fullyMergeable(numRegions, true);
    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        for (daf::Size v : regionVerts[rid]) {
            if (v >= numVertices) continue;
            for (daf::Size otherRid : vtxMaxPaths[v]) {
                if (otherRid == rid) continue;
                if (interSize[pairKey(rid, otherRid)] >= (int)r) {
                    fullyMergeable[rid] = false;
                    break;
                }
            }
            if (!fullyMergeable[rid]) break;
        }
    }

    // Directly assign fully-mergeable regions
    std::map<daf::Size, int64_t> coreDist;
    daf::Size numFullyMergeable = 0, mergedRCliques = 0;

    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        if (!fullyMergeable[rid]) continue;
        numFullyMergeable++;
        int mSize = (int)regionVerts[rid].size();
        daf::Size coreVal = (mSize >= (int)s) ? (daf::Size)llround(nCr[mSize - r][s - r]) : 0;
        daf::Size numRC = (daf::Size)llround(nCr[mSize][r]);
        coreDist[coreVal] += numRC;
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
        daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
        std::cout << "\n  All regions fully r-mergeable. No peeling needed." << std::endl;
        std::cout << "  Max core: " << maxCore << std::endl;
        for (auto &[core, cnt] : coreDist)
            std::cout << "  core=" << core << " count=" << cnt << std::endl;
        std::cout << "  Total time: " << totalMs << " ms" << std::endl;
        return {};
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

    // ============================================================
    // Step 3: Enumerate r-tuples (same as V2)
    // ============================================================

    std::unordered_map<TupleKey, daf::Size, TupleHash> rTupleIndex;
    struct RTuple { TupleKey key; daf::Size mult; };
    std::vector<RTuple> rTuples;

    {
        TupleKey cur; cur.reserve(r);
        auto addRTuple = [&](const TupleKey &key) {
            std::unordered_map<daf::Size, int> counts;
            for (auto c : key) counts[c]++;
            daf::Size mult = 1;
            for (auto &[c, k] : counts) {
                if ((int)classSizes[c] < k) return;
                mult *= (daf::Size)nCr[classSizes[c]][k];
            }
            if (mult == 0) return;
            if (rTupleIndex.count(key)) return;
            rTupleIndex[key] = rTuples.size();
            rTuples.push_back({key, mult});
        };

        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            auto &cids = classesInRegion[rid];
            if (cids.size() > 500) continue;
            cur.clear();
            std::function<void()> cb = [&]() { addRTuple(cur); };
            enumerateMultisets(cids, r, 0, cur, cb);
        }
    }

    auto tStep3 = std::chrono::high_resolution_clock::now();
    std::cout << "  r-tuples: " << rTuples.size() << std::endl;

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

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        // Compute h, p and class distribution on this path
        int h = 0, p = 0;
        std::unordered_map<daf::Size, std::pair<int,int>> classDistrib; // cid -> (nh, np)

        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            if (node.isPivot) { p++; classDistrib[cid].second++; }
            else { h++; classDistrib[cid].first++; }
        }

        if (h + p < (int)s) continue;
        pathsUsed++;

        // Collect unique classes on this path
        std::vector<daf::Size> pathClasses;
        for (auto &[cid, _] : classDistrib) pathClasses.push_back(cid);
        std::sort(pathClasses.begin(), pathClasses.end());

        // Enumerate r-multisets of this path's classes
        TupleKey cur; cur.reserve(r);
        std::function<void()> cb = [&]() {
            // Build composition: (class, count) pairs
            std::vector<std::pair<daf::Size, int>> tauClasses;
            {
                daf::Size prev = INVALID; int cnt = 0;
                for (auto c : cur) {
                    if (c == prev) cnt++;
                    else { if (prev != INVALID) tauClasses.push_back({prev, cnt}); prev = c; cnt = 1; }
                }
                if (prev != INVALID) tauClasses.push_back({prev, cnt});
            }

            // Feasibility: for each class c with j_c copies in tau,
            // need j_c <= n_h^c + n_p^c
            for (auto &[c, jc] : tauClasses) {
                auto it = classDistrib.find(c);
                if (it == classDistrib.end()) return;
                if (jc > it->second.first + it->second.second) return;
            }

            // Look up this r-tuple
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;

            // Convolution to compute AggrCount
            // f[t] = coefficient of x^t in Π_c g_c(x)
            // where g_c(x) = Σ_{b_c} C(nh_c, j_c-b_c) × C(np_c, b_c) × x^{b_c}

            std::vector<double> f = {1.0}; // polynomial coefficients, f[0]=1

            for (auto &[c, jc] : tauClasses) {
                auto &[nhc, npc] = classDistrib[c];
                int bMin = std::max(0, jc - nhc);
                int bMax = std::min(jc, npc);
                if (bMin > bMax) return; // infeasible

                // g_c coefficients
                std::vector<double> gc(bMax + 1, 0.0);
                for (int bc = bMin; bc <= bMax; ++bc)
                    gc[bc] = nCr[nhc][jc - bc] * nCr[npc][bc];

                // Convolve f with gc
                std::vector<double> newf(f.size() + gc.size() - 1, 0.0);
                for (int i = 0; i < (int)f.size(); ++i)
                    for (int j = 0; j < (int)gc.size(); ++j)
                        newf[i + j] += f[i] * gc[j];
                f = std::move(newf);
            }

            // AggrCount = Σ_t f[t] × C(p-t, s-h-t)
            double aggr = 0.0;
            for (int t = 0; t < (int)f.size(); ++t) {
                if (f[t] == 0.0) continue;
                int n = p - t, k = s - h - t;
                if (n >= 0 && k >= 0 && n >= k)
                    aggr += f[t] * nCr[n][k];
            }

            if (aggr > 0) {
                support[tidx] += aggr / rTuples[tidx].mult;
                totalTuplesOnPaths++;
            }
        };
        enumerateMultisets(pathClasses, r, 0, cur, cb);
    }

    auto tStep4 = std::chrono::high_resolution_clock::now();
    auto step4Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep4 - tStep3).count();

    double totalSupportTuples = 0, totalRCliques = 0;
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        totalSupportTuples += rTuples[i].mult * support[i];
        totalRCliques += rTuples[i].mult;
    }
    std::cout << "  CPI paths used: " << pathsUsed << std::endl;
    std::cout << "  Tuple-path pairs: " << totalTuplesOnPaths << std::endl;
    std::cout << "  r-cliques: " << std::fixed << std::setprecision(0) << totalRCliques << std::endl;
    std::cout << "  Support sum (CPI): " << totalSupportTuples << std::endl;
    std::cout << "  CPI counting time: " << step4Ms << " ms" << std::endl;

    // ============================================================
    // Step 5+6: Constrained Path Peeling (Analytical Split)
    // ============================================================
    // Constrained Path = CPI path + per-class (min_piv, max_piv) bounds.
    // When τ is peeled on P̂: subtract old contributions, split P̂ into
    // κ disjoint sub-paths (each a ConstrainedPath), add new contributions.
    // NO s-tuple enumeration. NO inclusion-exclusion. NO BK execution.

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
        for (auto &[cid, _] : cp.classes) pathClasses.push_back(cid);
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

    // --- Peeling with Analytical Split (double support + bucket queue) ---
    auto tStep6 = std::chrono::high_resolution_clock::now();

    // Profiling counters
    long long prof_oldContrib_ns = 0, prof_split_ns = 0, prof_newContrib_ns = 0;
    long long prof_delta_ns = 0, prof_overhead_ns = 0;
    long long prof_aggrCalls = 0, prof_splitSubpaths = 0, prof_deadPaths = 0;

    // Keep double support for exact arithmetic. Use floor() for bucket index.
    std::vector<double> dSup = support;
    std::vector<bool> rPeeled(rTuples.size(), false);

    daf::Size maxSup = 0;
    for (auto &sv : dSup) maxSup = std::max(maxSup, (daf::Size)llround(sv));
    const daf::Size BUCKET_CAP = std::min(maxSup + 2, (daf::Size)1000001);

    std::vector<std::vector<daf::Size>> buckets(BUCKET_CAP);
    std::multimap<daf::Size, daf::Size> overflow; // for large support values
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        daf::Size b = (daf::Size)llround(dSup[i]);
        if (b < BUCKET_CAP) buckets[b].push_back(i);
        else overflow.insert({b, i});
    }

    // coreDist already initialized in Step 1b (fully-mergeable regions)
    daf::Size numPeeled = 0, currentLevel = 0, coreLevel = 0;
    daf::Size totalSplits = 0;

    while (numPeeled < rTuples.size()) {
        while (currentLevel < BUCKET_CAP && buckets[currentLevel].empty())
            currentLevel++;
        daf::Size idx;
        if (currentLevel < BUCKET_CAP) {
            idx = buckets[currentLevel].back();
            buckets[currentLevel].pop_back();
        } else if (!overflow.empty()) {
            auto it = overflow.begin();
            currentLevel = it->first;
            idx = it->second;
            overflow.erase(it);
        } else break;

        if (rPeeled[idx]) continue;
        daf::Size idxBucket = (daf::Size)llround(dSup[idx]);
        if (idxBucket != currentLevel) continue;

        rPeeled[idx] = true;
        numPeeled++;
        coreLevel = std::max(coreLevel, currentLevel);
        coreDist[coreLevel] += rTuples[idx].mult;

        auto &tauKey = rTuples[idx].key;
        std::unordered_map<daf::Size, int> tauCounts;
        for (auto c : tauKey) tauCounts[c]++;

        auto cpathIds = tupleToCPaths[idx]; // copy
        for (daf::Size cpid : cpathIds) {
            auto &cp = cpaths[cpid];
            if (cp.tupleIdxs.empty()) continue;

            auto _t0 = std::chrono::high_resolution_clock::now();

            // Separate alive tuples into: affected (share class with τ) and unaffected
            std::unordered_set<daf::Size> tauClassSet(rTuples[idx].key.begin(), rTuples[idx].key.end());
            std::vector<daf::Size> affectedTuples, unaffectedTuples;
            for (auto tidx : cp.tupleIdxs) {
                if (rPeeled[tidx]) continue;
                bool shares = false;
                for (auto c : rTuples[tidx].key)
                    if (tauClassSet.count(c)) { shares = true; break; }
                if (shares) affectedTuples.push_back(tidx);
                else unaffectedTuples.push_back(tidx);
            }

            // Compute old contributions (only for AFFECTED tuples)
            std::unordered_map<daf::Size, double> oldContrib;
            for (auto tidx : affectedTuples) {
                oldContrib[tidx] = aggrCountOnCPath(tidx, cp);
                prof_aggrCalls++;
            }

            auto _t1 = std::chrono::high_resolution_clock::now();
            prof_oldContrib_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(_t1-_t0).count();

            // Find active classes
            struct ActiveClass { daf::Size cid; int mc; };
            std::vector<ActiveClass> active;
            for (auto &[c, jc] : tauCounts) {
                auto cit = cp.classes.find(c);
                if (cit == cp.classes.end()) continue;
                int mc = std::max(0, jc - cit->second.nh);
                if (mc > cit->second.minPiv)
                    active.push_back({c, mc});
            }

            // LAZY SPLIT (when s < 2r): parent survives with unaffected tuples.
            // Sub-paths only get affected tuples. Proven correct by Non-Interference Theorem.
            const bool lazy = false; // Lazy split disabled: transitive interference not handled

            // Save cp state BEFORE modifying tupleIdxs (sub-paths copy from saved state)
            CPath savedCp = cp; // copy classes/bounds (tupleIdxs will be overwritten anyway)

            if (lazy) {
                // Only remove affected + peeled from parent
                for (auto tidx : affectedTuples) {
                    auto &vec = tupleToCPaths[tidx];
                    vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end());
                }
                { auto &vec = tupleToCPaths[idx];
                  vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end()); }
                cp.tupleIdxs = unaffectedTuples; // parent keeps unaffected
            } else {
                for (auto tidx : cp.tupleIdxs) {
                    auto &vec = tupleToCPaths[tidx];
                    vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end());
                }
                cp.tupleIdxs.clear();
            }

            auto _t2 = std::chrono::high_resolution_clock::now();
            prof_overhead_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(_t2-_t1).count();

            // Compute new contributions from sub-paths
            std::unordered_map<daf::Size, double> newContrib;
            if (active.empty()) {
                prof_deadPaths++;
                // Path fully dead: unaffected on parent also lose support
                if (lazy) {
                    for (auto tidx : cp.tupleIdxs) {
                        double oc = aggrCountOnCPath(tidx, savedCp);
                        prof_aggrCalls++;
                        dSup[tidx] -= oc;
                        if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                        daf::Size nb = (daf::Size)llround(dSup[tidx]);
                        if (nb < BUCKET_CAP) { buckets[nb].push_back(tidx); if (nb < currentLevel) currentLevel = nb; }
                        else overflow.insert({nb, tidx});
                        auto &vec = tupleToCPaths[tidx];
                        vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end());
                    }
                    cp.tupleIdxs.clear();
                }
            } else {
                for (int i = 0; i < (int)active.size(); ++i) {
                    CPath subCp = savedCp; // use SAVED state (original bounds)
                    bool feasible = true;

                    for (int j = 0; j < i; ++j) {
                        auto &cb = subCp.classes[active[j].cid];
                        cb.minPiv = std::max(cb.minPiv, active[j].mc);
                        if (cb.minPiv > cb.maxPiv) { feasible = false; break; }
                    }
                    if (!feasible) continue;

                    {
                        auto &cb = subCp.classes[active[i].cid];
                        cb.maxPiv = std::min(cb.maxPiv, active[i].mc - 1);
                        if (cb.minPiv > cb.maxPiv) continue;
                    }

                    int minTotal = 0, maxTotal = 0;
                    for (auto &[_, cb] : subCp.classes) {
                        minTotal += cb.minPiv;
                        maxTotal += cb.maxPiv;
                    }
                    if ((s - subCp.h) < minTotal || (s - subCp.h) > maxTotal) continue;

                    subCp.tupleIdxs.clear();
                    for (auto tidx : affectedTuples) {
                        double c = aggrCountOnCPath(tidx, subCp);
                        prof_aggrCalls++;
                        if (c > 1e-9) {
                            subCp.tupleIdxs.push_back(tidx);
                            newContrib[tidx] += c;
                        }
                    }
                    // Lazy: unaffected stay on parent. Full: copy to sub-path.
                    if (!lazy) {
                        for (auto tidx : unaffectedTuples)
                            subCp.tupleIdxs.push_back(tidx);
                    }
                    if (subCp.tupleIdxs.empty()) continue;

                    daf::Size newCpid = cpaths.size();
                    cpaths.push_back(std::move(subCp));
                    totalSplits++;
                    prof_splitSubpaths++;

                    for (auto tidx : cpaths.back().tupleIdxs)
                        tupleToCPaths[tidx].push_back(newCpid);
                }
            }

            auto _t3 = std::chrono::high_resolution_clock::now();
            prof_split_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(_t3-_t2).count();

            // Apply net delta in double, re-bucket
            for (auto &[tidx, oc] : oldContrib) {
                double nc = newContrib.count(tidx) ? newContrib[tidx] : 0.0;
                dSup[tidx] += (nc - oc);
                if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                daf::Size newBucket = (daf::Size)llround(dSup[tidx]);
                if (newBucket < BUCKET_CAP) {
                    buckets[newBucket].push_back(tidx);
                    if (newBucket < currentLevel) currentLevel = newBucket;
                } else {
                    overflow.insert({newBucket, tidx});
                }
            }
        }
    }

    auto tStep6End = std::chrono::high_resolution_clock::now();
    auto step6Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStep6).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStart).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling (Analytical Split) ---" << std::endl;
    std::cout << "  Peeled: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Total splits: " << totalSplits << " (sub-paths: " << prof_splitSubpaths
              << ", dead: " << prof_deadPaths << ")" << std::endl;
    std::cout << "  Final cpaths: " << cpaths.size() << std::endl;
    std::cout << "  AggrCount calls: " << prof_aggrCalls << std::endl;
    std::cout << "  Time breakdown:" << std::endl;
    std::cout << "    oldContrib: " << prof_oldContrib_ns/1000000 << " ms" << std::endl;
    std::cout << "    split+new:  " << prof_split_ns/1000000 << " ms" << std::endl;
    std::cout << "    overhead:   " << prof_overhead_ns/1000000 << " ms" << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step6Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;

    return {};
}
