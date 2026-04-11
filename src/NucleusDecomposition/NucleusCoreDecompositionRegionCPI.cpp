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
    // Step 2: Build Overlap Classes (same as V2)
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
    // Step 5: Build tuple-path index for Strategy A peeling
    // ============================================================
    // For each path P: list of (tuple_idx, AggrCount_contribution) pairs.
    // For each tuple: list of paths containing it.
    // NO s-tuple enumeration.

    auto tStep5 = std::chrono::high_resolution_clock::now();

    // Path info: per path, store class distribution and tuple list
    struct PathTupleInfo {
        daf::Size tidx;
        double aggrContrib; // AggrCount(τ,P) / mult(τ) = per-r-clique support from this path
    };

    struct PathData {
        int h, p;
        std::unordered_map<daf::Size, std::pair<int,int>> classDistrib;
        std::vector<PathTupleInfo> tuples; // tuples on this path
    };

    std::vector<PathData> pathDataVec(numPaths);
    std::vector<std::vector<daf::Size>> tupleToPathIds(rTuples.size()); // tuple -> path IDs

    // Rebuild path class distributions (same as Step 4, but store them)
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        auto &pd = pathDataVec[pid];
        pd.h = 0; pd.p = 0;
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            if (node.isPivot) { pd.p++; pd.classDistrib[cid].second++; }
            else { pd.h++; pd.classDistrib[cid].first++; }
        }
    }

    // For each tuple-path pair found during Step 4: record it
    // (Re-enumerate, this time storing the mapping)
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &pd = pathDataVec[pid];
        if (pd.h + pd.p < (int)s) continue;

        std::vector<daf::Size> pathClasses;
        for (auto &[cid, _] : pd.classDistrib) pathClasses.push_back(cid);
        std::sort(pathClasses.begin(), pathClasses.end());

        TupleKey cur; cur.reserve(r);
        std::function<void()> cb = [&]() {
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;

            // Check feasibility and compute AggrCount (same convolution as Step 4)
            std::vector<std::pair<daf::Size, int>> tauClasses;
            {
                daf::Size prev = INVALID; int cnt = 0;
                for (auto c : cur) {
                    if (c == prev) cnt++;
                    else { if (prev != INVALID) tauClasses.push_back({prev, cnt}); prev = c; cnt = 1; }
                }
                if (prev != INVALID) tauClasses.push_back({prev, cnt});
            }
            for (auto &[c, jc] : tauClasses) {
                auto dit = pd.classDistrib.find(c);
                if (dit == pd.classDistrib.end()) return;
                if (jc > dit->second.first + dit->second.second) return;
            }

            std::vector<double> f = {1.0};
            for (auto &[c, jc] : tauClasses) {
                auto &[nhc, npc] = pd.classDistrib[c];
                int bMin = std::max(0, jc - nhc);
                int bMax = std::min(jc, npc);
                if (bMin > bMax) return;
                std::vector<double> gc(bMax + 1, 0.0);
                for (int bc = bMin; bc <= bMax; ++bc)
                    gc[bc] = nCr[nhc][jc - bc] * nCr[npc][bc];
                std::vector<double> newf(f.size() + gc.size() - 1, 0.0);
                for (int i = 0; i < (int)f.size(); ++i)
                    for (int j = 0; j < (int)gc.size(); ++j)
                        newf[i + j] += f[i] * gc[j];
                f = std::move(newf);
            }
            double aggr = 0.0;
            for (int t = 0; t < (int)f.size(); ++t) {
                if (f[t] == 0.0) continue;
                int n = pd.p - t, k = s - pd.h - t;
                if (n >= 0 && k >= 0 && n >= k) aggr += f[t] * nCr[n][k];
            }
            if (aggr <= 0) return;

            double contrib = aggr / rTuples[tidx].mult;
            pd.tuples.push_back({tidx, contrib});
            tupleToPathIds[tidx].push_back(pid);
        };
        enumerateMultisets(pathClasses, r, 0, cur, cb);
    }

    auto tStep5End = std::chrono::high_resolution_clock::now();
    auto step5Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStep5).count();
    daf::Size totalPathTuplePairs = 0;
    for (auto &pd : pathDataVec) totalPathTuplePairs += pd.tuples.size();
    std::cout << "  Tuple-path pairs (index): " << totalPathTuplePairs << std::endl;
    std::cout << "  Index build time: " << step5Ms << " ms" << std::endl;

    // ============================================================
    // Step 6: Strategy A Peeling (CPI + class removal, NO s-tuples)
    // ============================================================
    // Key insight: when all tuples using class c on path P are peeled,
    // remove c's vertices from P (update h, p). Future SharedCount
    // computations use updated parameters → dead s-cliques are
    // automatically excluded. NO inclusion-exclusion needed.

    auto tStep6 = std::chrono::high_resolution_clock::now();

    daf::Size maxSup = 0;
    for (auto &sv : support) maxSup = std::max(maxSup, (daf::Size)sv);

    std::vector<std::vector<daf::Size>> buckets(maxSup + 2);
    std::vector<daf::Size> curSupport(rTuples.size());
    std::vector<bool> rPeeled(rTuples.size(), false);

    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        curSupport[i] = (daf::Size)support[i];
        buckets[curSupport[i]].push_back(i);
    }

    // Per-path: alive tuple count per class (for class removal tracking)
    // aliveClassCount[pid][cid] = number of alive tuples on pid using class cid
    std::vector<std::unordered_map<daf::Size, int>> aliveClassCount(numPaths);
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        for (auto &pti : pathDataVec[pid].tuples) {
            auto &key = rTuples[pti.tidx].key;
            std::unordered_set<daf::Size> seen;
            for (auto c : key) {
                if (seen.insert(c).second) // count each class once per tuple
                    aliveClassCount[pid][c]++;
            }
        }
    }

    // Helper: compute SharedCount(τ, τ', P) using CURRENT path parameters
    // = number of s-cliques on P containing both τ and τ' r-cliques
    // Combined requirement: for each class c, need max(j_c^τ, j_c^{τ'}) vertices
    // → pivot need = max(0, max(j_c^τ, j_c^{τ'}) - current_nh_c)
    auto computeSharedCount = [&](daf::Size pid, daf::Size tidx1, daf::Size tidx2) -> double {
        auto &pd = pathDataVec[pid];
        auto &key1 = rTuples[tidx1].key;
        auto &key2 = rTuples[tidx2].key;

        // Build class counts for both tuples
        std::unordered_map<daf::Size, int> counts1, counts2;
        for (auto c : key1) counts1[c]++;
        for (auto c : key2) counts2[c]++;

        // Merged classes
        std::unordered_set<daf::Size> allClasses;
        for (auto &[c, _] : counts1) allClasses.insert(c);
        for (auto &[c, _] : counts2) allClasses.insert(c);

        // Combined requirement: max(j1_c, j2_c) per class
        // Constrained allocation: distribute s-h pivots with minimum per class
        int totalMinPiv = 0;
        struct ClassAlloc { daf::Size cid; int minPiv; int maxExtra; int npc; };
        std::vector<ClassAlloc> allocs;

        for (daf::Size c : allClasses) {
            auto dit = pd.classDistrib.find(c);
            if (dit == pd.classDistrib.end()) return 0.0;
            int nhc = dit->second.first, npc = dit->second.second;
            int j1 = counts1.count(c) ? counts1[c] : 0;
            int j2 = counts2.count(c) ? counts2[c] : 0;
            int need = std::max(j1, j2);
            int pivNeed = std::max(0, need - nhc);
            if (pivNeed > npc) return 0.0;
            totalMinPiv += pivNeed;
            allocs.push_back({c, pivNeed, npc - pivNeed, npc});
        }

        // Add classes on P not in either tuple (their pivots are freely available)
        for (auto &[c, dp] : pd.classDistrib) {
            if (allClasses.count(c)) continue;
            if (dp.second > 0) allocs.push_back({c, 0, dp.second, dp.second});
        }

        int extraNeeded = (s - pd.h) - totalMinPiv;
        if (extraNeeded < 0) return 0.0;

        // Convolution: distribute extraNeeded among classes
        std::vector<double> f = {1.0};
        for (auto &ac : allocs) {
            int maxE = std::min(ac.maxExtra, extraNeeded);
            if (maxE < 0) continue;
            std::vector<double> gc(maxE + 1, 0.0);
            for (int e = 0; e <= maxE; ++e)
                gc[e] = nCr[ac.npc][ac.minPiv + e];
            std::vector<double> newf(f.size() + gc.size() - 1, 0.0);
            for (int i = 0; i < (int)f.size(); ++i)
                for (int j = 0; j < (int)gc.size(); ++j)
                    newf[i + j] += f[i] * gc[j];
            f = std::move(newf);
        }
        return (extraNeeded < (int)f.size()) ? f[extraNeeded] : 0.0;
    };

    std::map<daf::Size, int64_t> coreDist;
    daf::Size numPeeled = 0, currentLevel = 0, coreLevel = 0;

    while (numPeeled < rTuples.size()) {
        while (currentLevel <= maxSup && buckets[currentLevel].empty())
            currentLevel++;
        if (currentLevel > maxSup) break;

        daf::Size idx = buckets[currentLevel].back();
        buckets[currentLevel].pop_back();
        if (rPeeled[idx] || curSupport[idx] != currentLevel) continue;

        rPeeled[idx] = true;
        numPeeled++;
        coreLevel = std::max(coreLevel, currentLevel);
        coreDist[coreLevel] += rTuples[idx].mult;

        auto &tauKey = rTuples[idx].key;

        for (daf::Size pid : tupleToPathIds[idx]) {
            auto &pd = pathDataVec[pid];

            // Step 1: Compute SharedCount(τ, τ', P) for each alive τ' on P
            // using CURRENT h/p (before class removal)
            for (auto &pti : pd.tuples) {
                if (pti.tidx == idx || rPeeled[pti.tidx]) continue;

                double shared = computeSharedCount(pid, idx, pti.tidx);
                double delta = shared / rTuples[pti.tidx].mult;

                if (delta > 0.5) {
                    daf::Size ridx2 = pti.tidx;
                    daf::Size oldSup = curSupport[ridx2];
                    daf::Size red = (daf::Size)(delta + 0.5);
                    daf::Size newSup = oldSup > red ? oldSup - red : 0;
                    if (newSup < oldSup) {
                        curSupport[ridx2] = newSup;
                        buckets[newSup].push_back(ridx2);
                        if (newSup < currentLevel) currentLevel = newSup;
                    }
                }
            }

            // Step 2: Update alive class counts and remove dead classes
            std::unordered_set<daf::Size> tauClasses;
            for (auto c : tauKey) tauClasses.insert(c);

            for (daf::Size c : tauClasses) {
                auto &cnt = aliveClassCount[pid][c];
                cnt--;
                if (cnt == 0) {
                    // All tuples using class c on this path are peeled
                    // Remove c's vertices from path: update h and p
                    auto dit = pd.classDistrib.find(c);
                    if (dit != pd.classDistrib.end()) {
                        pd.h -= dit->second.first;
                        pd.p -= dit->second.second;
                        pd.classDistrib.erase(dit);
                    }
                }
            }
        }
    }

    auto tStep6End = std::chrono::high_resolution_clock::now();
    auto step6Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStep6).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStart).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling ---" << std::endl;
    std::cout << "  Peeled: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step6Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;

    return {};
}
