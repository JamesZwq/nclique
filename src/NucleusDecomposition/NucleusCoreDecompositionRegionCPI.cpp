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
    // Step 5+6: ST-style Tree Mutation with Region Compression
    // ============================================================
    // Key idea: maintain the CPI tree (mutable, like ST). When a class
    // is fully peeled on a path → remove its vertices → path shrinks.
    // If path becomes infeasible → LeafDeath (subtract full contribution).
    // Support recomputed via AggrCount on the modified tree.

    auto tStep5 = std::chrono::high_resolution_clock::now();

    // --- Per-path mutable state ---
    struct PathState {
        std::unordered_map<daf::Size, std::pair<int,int>> classDistrib; // cid -> (nh, np)
        int h = 0, p = 0;
        std::vector<daf::Size> tupleIdxs; // alive tuples on this path
        bool alive = true;
    };

    std::vector<PathState> pathStates(numPaths);
    std::vector<std::vector<daf::Size>> tupleToPathIds(rTuples.size());

    // Per (path, class): alive tuple count
    std::vector<std::unordered_map<daf::Size, int>> aliveClassCount(numPaths);

    // Initialize path states from CPI tree
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        auto &ps = pathStates[pid];
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            if (node.isPivot) { ps.p++; ps.classDistrib[cid].second++; }
            else { ps.h++; ps.classDistrib[cid].first++; }
        }
        if (ps.h + ps.p < (int)s) continue;

        // Enumerate tuples on this path
        std::vector<daf::Size> pathClasses;
        for (auto &[cid, _] : ps.classDistrib) pathClasses.push_back(cid);
        std::sort(pathClasses.begin(), pathClasses.end());

        TupleKey cur; cur.reserve(r);
        std::function<void()> cb = [&]() {
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;
            std::unordered_map<daf::Size, int> counts;
            for (auto c : cur) counts[c]++;
            for (auto &[c, jc] : counts) {
                auto cit = ps.classDistrib.find(c);
                if (cit == ps.classDistrib.end()) return;
                if (jc > cit->second.first + cit->second.second) return;
            }
            ps.tupleIdxs.push_back(tidx);
            tupleToPathIds[tidx].push_back(pid);
            // Update alive class count
            std::unordered_set<daf::Size> seen;
            for (auto c : cur)
                if (seen.insert(c).second) aliveClassCount[pid][c]++;
        };
        enumerateMultisets(pathClasses, r, 0, cur, cb);
        ps.alive = !ps.tupleIdxs.empty();
    }

    // AggrCount on mutable path (same Vandermonde formula, no min/max bounds)
    auto aggrCount = [&](daf::Size tidx, const PathState &ps) -> double {
        auto &key = rTuples[tidx].key;
        int targetTotal = s - ps.h;
        if (targetTotal < 0) return 0.0;

        double f[512]; int fLen = 1; f[0] = 1.0;

        int ki = 0;
        while (ki < (int)key.size()) {
            daf::Size c = key[ki];
            int jc = 1;
            while (ki + jc < (int)key.size() && key[ki + jc] == c) jc++;
            ki += jc;

            auto cit = ps.classDistrib.find(c);
            if (cit == ps.classDistrib.end()) return 0.0;
            int nhc = cit->second.first, npc = cit->second.second;
            int tMin = std::max(0, jc - nhc);
            int tMax = std::min({jc, npc, (int)500});
            if (jc > nhc + npc || tMin > tMax) return 0.0;

            double gc[512]; int gcLen = tMax + 1;
            for (int i = 0; i < gcLen; ++i) gc[i] = 0.0;
            for (int tc = tMin; tc <= tMax; ++tc)
                gc[tc] = nCr[npc][tc] * nCr[nhc + tc][jc];
            int newLen = fLen + gcLen - 1;
            double newf[512];
            for (int i = 0; i < newLen; ++i) newf[i] = 0.0;
            for (int i = 0; i < fLen; ++i)
                for (int j = 0; j < gcLen; ++j)
                    newf[i + j] += f[i] * gc[j];
            fLen = std::min(newLen, targetTotal + 1);
            for (int i = 0; i < fLen; ++i) f[i] = newf[i];
        }

        for (auto &[c, dp] : ps.classDistrib) {
            bool inTuple = false;
            for (auto x : key) if (x == c) { inTuple = true; break; }
            if (inTuple) continue;
            int npc = dp.second;
            if (npc == 0) continue; // skip empty classes
            double gc[512]; int gcLen = std::min(npc + 1, 512);
            for (int i = 0; i < gcLen; ++i) gc[i] = nCr[npc][i];
            int newLen = fLen + gcLen - 1;
            double newf[512];
            for (int i = 0; i < newLen; ++i) newf[i] = 0.0;
            for (int i = 0; i < fLen; ++i)
                for (int j = 0; j < gcLen; ++j)
                    newf[i + j] += f[i] * gc[j];
            fLen = std::min(newLen, targetTotal + 1);
            for (int i = 0; i < fLen; ++i) f[i] = newf[i];
        }

        return (targetTotal < fLen) ? f[targetTotal] / rTuples[tidx].mult : 0.0;
    };

    auto tStep5End = std::chrono::high_resolution_clock::now();
    auto step5Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStep5).count();
    daf::Size numAlivePaths = 0;
    for (auto &ps : pathStates) if (ps.alive) numAlivePaths++;
    std::cout << "  Alive paths: " << numAlivePaths << std::endl;
    std::cout << "  Index build time: " << step5Ms << " ms" << std::endl;

    // --- Peeling ---
    auto tStep6 = std::chrono::high_resolution_clock::now();

    std::vector<double> dSup = support;
    std::vector<bool> rPeeled(rTuples.size(), false);

    // Per-(tuple, path) contribution cache for subtract/add updates
    std::unordered_map<daf::Size, double> cachedContrib;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &ps = pathStates[pid];
        if (!ps.alive) continue;
        for (auto tidx : ps.tupleIdxs) {
            daf::Size cacheKey = (daf::Size)pid * rTuples.size() + tidx;
            cachedContrib[cacheKey] = aggrCount(tidx, ps);
        }
    }

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

    std::map<daf::Size, int64_t> coreDist;
    daf::Size numPeeled = 0, currentLevel = 0, coreLevel = 0;
    long long totalClassRemovals = 0, totalLeafDeaths = 0, totalRecomputes = 0;

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
        if ((daf::Size)llround(dSup[idx]) != currentLevel) continue;

        rPeeled[idx] = true;
        numPeeled++;
        coreLevel = std::max(coreLevel, currentLevel);
        coreDist[coreLevel] += rTuples[idx].mult;

        auto &tauKey = rTuples[idx].key;
        std::unordered_set<daf::Size> tauClasses(tauKey.begin(), tauKey.end());

        // Process each path containing this tuple
        for (daf::Size pid : tupleToPathIds[idx]) {
            auto &ps = pathStates[pid];
            if (!ps.alive) continue;

            // Step 1: Decrement alive class counts, check for class removal
            bool treeModified = false;
            for (daf::Size c : tauClasses) {
                auto &cnt = aliveClassCount[pid][c];
                cnt--;
                if (cnt == 0) {
                    // Class c fully dead on path pid → remove vertices
                    auto dit = ps.classDistrib.find(c);
                    if (dit != ps.classDistrib.end()) {
                        ps.h -= dit->second.first;
                        ps.p -= dit->second.second;
                        ps.classDistrib.erase(dit);
                        treeModified = true;
                        totalClassRemovals++;
                    }
                }
            }

            // Step 2: Check LeafDeath (path too small for s-cliques)
            if (treeModified && ps.h + ps.p < (int)s) {
                // LeafDeath: subtract ALL remaining contributions, kill path
                totalLeafDeaths++;
                for (auto tidx : ps.tupleIdxs) {
                    if (rPeeled[tidx]) continue;
                    daf::Size ck = (daf::Size)pid * rTuples.size() + tidx;
                    double oldC = cachedContrib.count(ck) ? cachedContrib[ck] : 0.0;
                    dSup[tidx] -= oldC;
                    if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                    cachedContrib.erase(ck);
                    daf::Size newB = (daf::Size)llround(dSup[tidx]);
                    if (newB < BUCKET_CAP) {
                        buckets[newB].push_back(tidx);
                        if (newB < currentLevel) currentLevel = newB;
                    } else overflow.insert({newB, tidx});
                }
                for (auto tidx : ps.tupleIdxs) {
                    auto &vec = tupleToPathIds[tidx];
                    vec.erase(std::remove(vec.begin(), vec.end(), pid), vec.end());
                }
                ps.tupleIdxs.clear();
                ps.alive = false;
                continue;
            }

            // Step 3: If tree was modified (class removed but path survived):
            // Recompute AggrCount for alive tuples and update dSup (subtract old, add new)
            if (treeModified) {
                totalRecomputes++;
                // Remove peeled tuples
                std::vector<daf::Size> alive;
                for (auto tidx : ps.tupleIdxs)
                    if (!rPeeled[tidx]) alive.push_back(tidx);
                ps.tupleIdxs.clear();

                for (auto tidx : alive) {
                    // Subtract old cached contribution
                    daf::Size cacheKey = (daf::Size)pid * rTuples.size() + tidx;
                    double oldC = cachedContrib.count(cacheKey) ? cachedContrib[cacheKey] : 0.0;
                    // Compute new contribution on modified path
                    double newC = aggrCount(tidx, ps);
                    if (newC > 1e-9) {
                        ps.tupleIdxs.push_back(tidx);
                        cachedContrib[cacheKey] = newC;
                    } else {
                        cachedContrib.erase(cacheKey);
                        // Remove path from tuple's path list
                        auto &vec = tupleToPathIds[tidx];
                        vec.erase(std::remove(vec.begin(), vec.end(), pid), vec.end());
                    }
                    // Update support
                    dSup[tidx] += (newC - oldC);
                    if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                    daf::Size newB = (daf::Size)llround(dSup[tidx]);
                    if (newB < BUCKET_CAP) {
                        buckets[newB].push_back(tidx);
                        if (newB < currentLevel) currentLevel = newB;
                    } else overflow.insert({newB, tidx});
                }
            }

            // Step 4: If NO class was removed on this path: use Theorem 2 delta
            // to update neighbors' support (approximation, corrected on next class removal)
            if (!treeModified) {
                auto &tauKeyRef = rTuples[idx].key;
                std::unordered_map<daf::Size, int> tCounts;
                for (auto c : tauKeyRef) tCounts[c]++;

                for (auto tidx : ps.tupleIdxs) {
                    if (rPeeled[tidx] || tidx == idx) continue;
                    // Theorem 2: Δ = constrained allocation for specific C_r'
                    auto &tpKey = rTuples[tidx].key;
                    std::unordered_map<daf::Size, int> tpCounts;
                    for (auto c : tpKey) tpCounts[c]++;

                    // For each class: compute l_c, u_c
                    int totalMand = 0;
                    bool feasible = true;
                    struct AI { int lc, uc; };
                    std::vector<AI> allocs;
                    for (auto &[c, dp] : ps.classDistrib) {
                        int nhc = dp.first, npc = dp.second;
                        int jTau = tCounts.count(c) ? tCounts[c] : 0;
                        int jTauP = tpCounts.count(c) ? tpCounts[c] : 0;
                        int qpc = std::max(0, jTauP - nhc);
                        int mc = std::max(0, jTau - nhc);
                        int lc = std::max(0, mc - qpc);
                        int uc = npc - qpc;
                        if (uc < 0 || uc < lc) { feasible = false; break; }
                        totalMand += qpc + lc;
                        allocs.push_back({lc, uc});
                    }
                    if (!feasible) continue;
                    int extra = (s - ps.h) - totalMand;
                    if (extra < 0) continue;

                    double f[64]; int fLen = 1; f[0] = 1.0;
                    for (auto &ai : allocs) {
                        int maxE = std::min(ai.uc - ai.lc, extra);
                        if (maxE < 0) { feasible = false; break; }
                        double gc[64]; int gcLen = maxE + 1;
                        for (int e = 0; e <= maxE; ++e)
                            gc[e] = nCr[ai.uc][ai.lc + e];
                        int newLen = fLen + gcLen - 1;
                        double nf[128];
                        for (int i = 0; i < newLen; ++i) nf[i] = 0.0;
                        for (int i = 0; i < fLen; ++i)
                            for (int j = 0; j < gcLen; ++j)
                                nf[i+j] += f[i] * gc[j];
                        fLen = std::min(newLen, extra + 1);
                        for (int i = 0; i < fLen; ++i) f[i] = nf[i];
                    }
                    if (!feasible) continue;
                    double delta = (extra < fLen) ? f[extra] : 0.0;

                    if (delta > 0.5) {
                        dSup[tidx] -= delta;
                        if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                        daf::Size newB = (daf::Size)llround(dSup[tidx]);
                        if (newB < BUCKET_CAP) {
                            buckets[newB].push_back(tidx);
                            if (newB < currentLevel) currentLevel = newB;
                        } else overflow.insert({newB, tidx});
                    }
                }
            }

            // Remove peeled tuple from path
            {
                auto &vec = ps.tupleIdxs;
                vec.erase(std::remove(vec.begin(), vec.end(), idx), vec.end());
            }
        }
    }

    auto tStep6End = std::chrono::high_resolution_clock::now();
    auto step6Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStep6).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStart).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling (Tree Mutation) ---" << std::endl;
    std::cout << "  Peeled: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Class removals: " << totalClassRemovals
              << ", LeafDeaths: " << totalLeafDeaths
              << ", Recomputes: " << totalRecomputes << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step6Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;

    return {};
}
