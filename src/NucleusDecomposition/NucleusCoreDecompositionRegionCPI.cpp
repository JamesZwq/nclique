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

    // --- Weighted feasible count: convolution DP ---
    // For each class i, weight_i(b) depends on lower[i]..upper[i] and the
    // class's nh/np and tuple's j_c. We pass per-class weight tables.
    // weightTables[i][b - lower[i]] = weight for allocation b on class i.
    // Returns Σ_{b: lower<=b<=upper, Σb=T} Π_i weightTables[i][b_i - lower[i]]
    auto countFeasibleWeighted = [](const std::vector<int> &lower,
                                    const std::vector<int> &upper, int T,
                                    const std::vector<std::vector<double>> &weightTables) -> double {
        int m = (int)lower.size();
        int minSum = 0, maxSum = 0;
        for (int i = 0; i < m; ++i) {
            if (lower[i] > upper[i]) return 0.0;
            minSum += lower[i];
            maxSum += upper[i];
        }
        if (T < minSum || T > maxSum) return 0.0;

        int target = T - minSum;
        if (target < 0) return 0.0;

        // dp[t] = weighted count of allocations with shifted-sum = t
        // Process each class: convolve dp with class's weight polynomial
        std::vector<double> dp(target + 1, 0.0);
        dp[0] = 1.0;
        for (int i = 0; i < m; ++i) {
            int cap = upper[i] - lower[i]; // range: 0..cap (shifted)
            auto &wt = weightTables[i]; // wt[k] = weight for b_i = lower[i] + k
            std::vector<double> next(target + 1, 0.0);
            // For each target sum t, next[t] = Σ_{k=0..min(cap,t)} dp[t-k] * wt[k]
            for (int t = 0; t <= target; ++t) {
                double sum = 0.0;
                int kMax = std::min(cap, t);
                for (int k = 0; k <= kMax; ++k) {
                    if (k < (int)wt.size())
                        sum += dp[t - k] * wt[k];
                }
                next[t] = sum;
            }
            dp.swap(next);
        }
        return dp[target];
    };

    // --- Build weight tables for tuple τ' on path P with given lower/upper ---
    // tc = total pivot allocation from class c in the s-clique.
    // For tuple classes (j'_c > 0): w(tc) = C(np_c, tc) × C(nh_c + tc, j'_c)
    //   (choose tc pivots, then choose j'_c of the first nh_c + tc vertices for the tuple)
    // For non-tuple classes (j'_c = 0): w(tc) = C(np_c, tc)
    // Returns weightTables[i][tc - lower[i]] for each class i
    auto buildWeightTables = [&](daf::Size tidx, const PathInfo &pi,
                                 const std::vector<int> &lower,
                                 const std::vector<int> &upper) -> std::vector<std::vector<double>> {
        int m = (int)pi.classIds.size();
        auto &key = rTuples[tidx].key;
        std::unordered_map<daf::Size, int> counts;
        for (auto c : key) counts[c]++;

        std::vector<std::vector<double>> wts(m);
        for (int i = 0; i < m; ++i) {
            int lo = lower[i], hi = upper[i];
            int range = hi - lo + 1;
            if (range <= 0) { wts[i] = {}; continue; }
            wts[i].resize(range, 0.0);

            auto cit = counts.find(pi.classIds[i]);
            int jc = (cit != counts.end()) ? cit->second : 0;
            int nhc = pi.nh[i], npc = pi.np[i];

            for (int k = 0; k < range; ++k) {
                int tc = lo + k; // total pivots allocated from class i
                if (tc < 0 || tc > npc) continue;
                if (jc > 0) {
                    // tuple class: C(np_c, tc) × C(nh_c + tc, j_c)
                    int pool = nhc + tc;
                    if (pool >= jc)
                        wts[i][k] = nCr[npc][tc] * nCr[pool][jc];
                } else {
                    // non-tuple class: C(np_c, tc)
                    wts[i][k] = nCr[npc][tc];
                }
            }
        }
        return wts;
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
    auto feasWeighted = [&](const std::vector<int> &lower, const std::vector<int> &upper,
                            UnionCtx &ctx) -> double {
        auto wts = buildWeightTables(ctx.tidx, *ctx.pi, lower, upper);
        return countFeasibleWeighted(lower, upper, ctx.T, wts);
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
        std::unordered_map<daf::Size, std::pair<int,int>> cd; // cid -> (nh, np)
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            if (node.isPivot) cd[cid].second++;
            else { cd[cid].first++; pi.h++; }
        }
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
    {
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
                auto wts = buildWeightTables(tidx, pi, base, upper);
                double aggr = countFeasibleWeighted(base, upper, pi.T, wts);
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

    // --- Per-(tuple, path) dead count cache ---
    daf::Size numTuplesSz = rTuples.size();
    std::unordered_map<uint64_t, double> deadCache;

    // --- Bucket queue setup ---
    std::vector<double> dSup = support;
    std::vector<bool> rPeeled(rTuples.size(), false);

    daf::Size maxSup = 0;
    for (auto &sv : dSup) maxSup = std::max(maxSup, (daf::Size)llround(sv));
    const daf::Size BUCKET_CAP = std::min(maxSup + 2, (daf::Size)1000001);

    std::vector<std::vector<daf::Size>> buckets(BUCKET_CAP);
    std::multimap<daf::Size, daf::Size> overflow;
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        daf::Size b = (daf::Size)llround(dSup[i]);
        if (b < BUCKET_CAP) buckets[b].push_back(i);
        else overflow.insert({b, i});
    }

    // Profiling counters
    long long prof_unionCalls = 0, prof_totalRecCalls = 0;
    long long prof_deadBoxesAdded = 0, prof_tupleUpdates = 0;
    long long prof_batchCount = 0;

    // --- Batch peeling loop ---
    daf::Size numPeeled = 0, currentLevel = 0, coreLevel = 0;

    while (numPeeled < rTuples.size()) {
        // Find next non-empty level
        while (currentLevel < BUCKET_CAP && buckets[currentLevel].empty())
            currentLevel++;

        // Drain ALL tuples at currentLevel (batch)
        // Use rPeeled to deduplicate (mark during drain, not after)
        std::vector<daf::Size> batch;
        if (currentLevel < BUCKET_CAP) {
            while (!buckets[currentLevel].empty()) {
                daf::Size idx = buckets[currentLevel].back();
                buckets[currentLevel].pop_back();
                if (rPeeled[idx]) continue;
                daf::Size idxBucket = (daf::Size)llround(dSup[idx]);
                if (idxBucket != currentLevel) continue;
                rPeeled[idx] = true; // mark immediately to avoid duplicates
                batch.push_back(idx);
            }
        } else if (!overflow.empty()) {
            currentLevel = overflow.begin()->first;
            auto range = overflow.equal_range(currentLevel);
            for (auto it = range.first; it != range.second; ++it) {
                daf::Size idx = it->second;
                if (!rPeeled[idx]) {
                    daf::Size idxBucket = (daf::Size)llround(dSup[idx]);
                    if (idxBucket == currentLevel) {
                        rPeeled[idx] = true;
                        batch.push_back(idx);
                    }
                }
            }
            overflow.erase(range.first, range.second);
        } else break;

        if (batch.empty()) { currentLevel++; continue; }
        prof_batchCount++;

        // Assign core level
        for (auto idx : batch) {
            numPeeled++;
            coreLevel = std::max(coreLevel, currentLevel);
            coreDist[coreLevel] += rTuples[idx].mult;
        }

        // Collect affected paths
        std::unordered_set<daf::Size> affectedPathSet;
        for (auto idx : batch)
            for (auto piIdx : tupleToPathInfos[idx])
                affectedPathSet.insert(piIdx);

        // Process each affected path ONCE
        for (daf::Size piIdx : affectedPathSet) {
            auto &pi = pathInfos[piIdx];
            int m = (int)pi.classIds.size();

            // Add new dead boxes from batch tuples on this path
            for (auto idx : batch) {
                // Check if this tuple is actually on this path
                bool onPath = false;
                for (auto p : tupleToPathInfos[idx])
                    if (p == piIdx) { onPath = true; break; }
                if (!onPath) continue;

                std::vector<int> req = buildReqVec(idx, pi);
                pi.deadBoxes.push_back(std::move(req));
                prof_deadBoxesAdded++;
            }

            // Collect alive tuples (remove peeled ones)
            std::vector<daf::Size> alive;
            alive.reserve(pi.tupleIdxs.size());
            for (auto tidx : pi.tupleIdxs)
                if (!rPeeled[tidx])
                    alive.push_back(tidx);
            pi.tupleIdxs = std::move(alive);

            if (pi.deadBoxes.empty()) continue;

            for (auto tidx : pi.tupleIdxs) {
                // Build base (lower bound) and upper for this tuple on this path
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

                // Compute new dead count via countUnionWeighted
                UnionCtx ctx{m, pi.T, 0, tidx, &pi};
                double newDead = countUnionRec(base, upper, pi.deadBoxes, ctx);
                prof_totalRecCalls += ctx.recCalls;
                prof_unionCalls++;

                // Look up cached old dead count
                uint64_t cacheKey = (uint64_t)piIdx * numTuplesSz + tidx;
                double oldDead = 0.0;
                auto dit = deadCache.find(cacheKey);
                if (dit != deadCache.end()) oldDead = dit->second;

                double delta = newDead - oldDead;
                if (delta > 0.5) {
                    deadCache[cacheKey] = newDead;
                    dSup[tidx] -= delta / rTuples[tidx].mult;
                    if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                    prof_tupleUpdates++;

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
    }

    auto tStep6End = std::chrono::high_resolution_clock::now();
    auto step6Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStep6).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStart).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling (Batch Union) ---" << std::endl;
    std::cout << "  Peeled: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Batches: " << prof_batchCount << std::endl;
    std::cout << "  Dead boxes added: " << prof_deadBoxesAdded << std::endl;
    std::cout << "  Union calls: " << prof_unionCalls << std::endl;
    std::cout << "  Total recursive calls: " << prof_totalRecCalls << std::endl;
    std::cout << "  Tuple updates: " << prof_tupleUpdates << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step6Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;

    return {};
}
