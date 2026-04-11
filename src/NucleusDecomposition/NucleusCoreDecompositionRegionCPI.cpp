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
    // Step 1: Build Maximal-Clique Regions (same as V2)
    // ============================================================

    std::vector<bool> pathValid(numPaths, false);
    daf::Size validPaths = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if ((int)tree.adj_list[pid].size() >= s) { pathValid[pid] = true; validPaths++; }
    }

    std::vector<daf::Size> regionOf(numPaths, INVALID);
    std::vector<daf::Size> maximalPathIds;
    daf::Size numMaximal = 0, numSub = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (!pathValid[pid]) continue;
        bool isMax = (pid < g_maxCliqueTags.size()) ? g_maxCliqueTags[pid] : true;
        if (isMax) { regionOf[pid] = maximalPathIds.size(); maximalPathIds.push_back(pid); numMaximal++; }
    }
    daf::Size numRegions = maximalPathIds.size();

    std::vector<std::vector<daf::Size>> regionVerts(numRegions);
    std::vector<std::vector<daf::Size>> vtxMaxPaths(numVertices);
    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        auto &leaf = tree.adj_list[maximalPathIds[rid]];
        regionVerts[rid].reserve(leaf.size());
        for (const auto &node : leaf) regionVerts[rid].push_back(node.v);
        std::sort(regionVerts[rid].begin(), regionVerts[rid].end());
        for (daf::Size v : regionVerts[rid])
            if (v < numVertices) vtxMaxPaths[v].push_back(rid);
    }

    daf::Size orphans = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (!pathValid[pid] || regionOf[pid] != INVALID) continue;
        numSub++;
        auto &leaf = tree.adj_list[pid];
        std::vector<daf::Size> pVerts;
        pVerts.reserve(leaf.size());
        for (const auto &node : leaf) pVerts.push_back(node.v);
        std::sort(pVerts.begin(), pVerts.end());
        daf::Size rareV = pVerts[0], rareCount = vtxMaxPaths[pVerts[0]].size();
        for (daf::Size v : pVerts)
            if (vtxMaxPaths[v].size() < rareCount) { rareV = v; rareCount = vtxMaxPaths[v].size(); }
        for (daf::Size rid : vtxMaxPaths[rareV]) {
            auto &qv = regionVerts[rid];
            if (qv.size() <= pVerts.size()) continue;
            if (std::includes(qv.begin(), qv.end(), pVerts.begin(), pVerts.end())) {
                regionOf[pid] = rid; break;
            }
        }
        if (regionOf[pid] == INVALID) {
            regionOf[pid] = maximalPathIds.size();
            maximalPathIds.push_back(pid);
            regionVerts.push_back(pVerts);
            for (daf::Size v : pVerts)
                if (v < numVertices) vtxMaxPaths[v].push_back(regionVerts.size() - 1);
            numRegions++; orphans++;
        }
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
    std::cout << "  Vertices: " << numVertices << ", paths: " << validPaths << std::endl;
    std::cout << "  Maximal cliques: " << numMaximal << ", sub merged: " << numSub
              << " (orphans: " << orphans << ")" << std::endl;
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
    // support(τ) = [Σ_P N(τ,P) × C(p-r+h, s-r)] / mult(τ)
    // where N(τ,P) = Π_c C(n_p^c, j_c - n_h^c)

    std::vector<double> support(rTuples.size(), 0.0);

    // Per-path info for CPI counting
    struct PathClassInfo { daf::Size cid; int nh, np; };

    // Also build: for each r-tuple, list of regions containing it (for lazy s-tuple enum)
    // (reuse classesInRegion for this)

    daf::Size pathsUsed = 0;
    daf::Size totalTuplesOnPaths = 0;

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue; // path too small for s-cliques

        // Compute h, p and class distribution on this path
        int h = 0, p = 0;
        std::unordered_map<daf::Size, std::pair<int,int>> classDistrib; // cid -> (nh, np)
        std::vector<daf::Size> pathClasses; // unique sorted class IDs on this path

        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            if (node.isPivot) { p++; classDistrib[cid].second++; }
            else { h++; classDistrib[cid].first++; }
        }

        if (h + p < (int)s) continue; // not enough vertices
        pathsUsed++;

        // Per-path support contribution (same for ALL r-cliques on this path)
        double perPathSup = (p - r + h >= 0 && s - r >= 0) ? nCr[p - r + h][s - r] : 0.0;
        if (perPathSup <= 0) continue;

        // Collect unique classes on this path
        for (auto &[cid, _] : classDistrib) pathClasses.push_back(cid);
        std::sort(pathClasses.begin(), pathClasses.end());

        // Enumerate r-multisets of this path's classes
        TupleKey cur; cur.reserve(r);
        std::function<void()> cb = [&]() {
            // Check feasibility: for each class c in cur with count j_c:
            // need n_h^c <= j_c and j_c - n_h^c <= n_p^c
            std::unordered_map<daf::Size, int> counts;
            for (auto c : cur) counts[c]++;
            double N = 1.0;
            for (auto &[c, jc] : counts) {
                auto it = classDistrib.find(c);
                if (it == classDistrib.end()) return;
                int nhc = it->second.first, npc = it->second.second;
                if (nhc > jc) return; // too many holds for this class
                int pivNeeded = jc - nhc;
                if (pivNeeded > npc) return; // not enough pivots
                N *= nCr[npc][pivNeeded];
            }
            // Also check classes in classDistrib but NOT in counts: nh must be 0
            // (holds from classes not used by tau are mandatory in every r-clique on P,
            //  but tau doesn't use them => conflict if nh > 0)
            // Actually: ALL holds must be in the r-clique (encoding condition).
            // So for any class c with nh_c > 0: tau must use at least nh_c from c.
            for (auto &[c, dp] : classDistrib) {
                if (counts.count(c)) continue;
                if (dp.first > 0) return; // class c has holds but tau doesn't use it
            }

            if (N <= 0) return;

            // Look up this r-tuple
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;

            daf::Size tidx = it->second;
            support[tidx] += N * perPathSup / rTuples[tidx].mult;
            totalTuplesOnPaths++;
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
    // Step 5: Peeling with lazy s-tuple generation
    // ============================================================

    auto tStep5 = std::chrono::high_resolution_clock::now();

    // Build incidence lazily: for each r-tuple, list of regions containing all its classes
    std::vector<std::vector<daf::Size>> tupleRegions(rTuples.size());
    for (daf::Size tidx = 0; tidx < rTuples.size(); ++tidx) {
        auto &key = rTuples[tidx].key;
        // Find regions containing ALL classes in the tuple
        // A region R contains all classes iff classesInRegion[R] ⊇ unique classes of key
        std::unordered_set<daf::Size> neededClasses(key.begin(), key.end());

        // Start with regions of first class, intersect
        if (neededClasses.empty()) continue;
        daf::Size firstClass = *neededClasses.begin();
        for (daf::Size rid : classes[firstClass].regionIds) {
            bool allPresent = true;
            for (daf::Size c : neededClasses) {
                auto &cr = classesInRegion[rid];
                if (!std::binary_search(cr.begin(), cr.end(), c)) {
                    allPresent = false; break;
                }
            }
            if (allPresent) tupleRegions[tidx].push_back(rid);
        }
    }

    // Bucket queue peeling
    daf::Size maxSup = 0;
    for (auto &s : support) maxSup = std::max(maxSup, (daf::Size)s);

    std::vector<std::vector<daf::Size>> buckets(maxSup + 2);
    std::vector<daf::Size> curSupport(rTuples.size());
    std::vector<bool> rPeeled(rTuples.size(), false);

    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        curSupport[i] = (daf::Size)support[i];
        buckets[curSupport[i]].push_back(i);
    }

    // Dead s-tuple tracking
    std::unordered_set<TupleKey, TupleHash> deadSTuples;

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

        // Lazy s-tuple cascade: for each region containing tau, enumerate dying s-tuples
        auto &tauKey = rTuples[idx].key;
        std::unordered_map<daf::Size, int> tauCounts;
        for (auto c : tauKey) tauCounts[c]++;

        for (daf::Size rid : tupleRegions[idx]) {
            auto &cids = classesInRegion[rid];
            if (cids.size() > 500) continue;

            // Enumerate s-multisets from cids that contain tauKey as sub-multiset
            TupleKey cur; cur.reserve(s);

            // Pre-fill with tau's classes
            cur = tauKey;

            // Fill remaining s-r slots from cids (with repetition up to class size)
            std::function<void(int)> fillRemaining = [&](int startIdx) {
                if ((int)cur.size() == s) {
                    // Validate: for each class, count ≤ classSize
                    std::unordered_map<daf::Size, int> counts;
                    for (auto c : cur) counts[c]++;
                    for (auto &[c, k] : counts)
                        if ((int)classSizes[c] < k) return;

                    TupleKey key = cur;
                    std::sort(key.begin(), key.end());

                    // Check if already dead
                    if (deadSTuples.count(key)) return;

                    // Mark dead
                    deadSTuples.insert(key);

                    // Find r-sub-tuples and update support
                    auto comp = getComposition(key);
                    TupleKey subKey; subKey.reserve(r);
                    std::function<void(const TupleKey &)> subCb = [&](const TupleKey &sub) {
                        auto it = rTupleIndex.find(sub);
                        if (it == rTupleIndex.end()) return;
                        daf::Size ridx2 = it->second;
                        if (ridx2 == idx || rPeeled[ridx2]) return;
                        double ext = computeExt(comp, sub, classSizes);
                        if (ext <= 0) return;
                        daf::Size oldSup = curSupport[ridx2];
                        daf::Size red = (daf::Size)ext;
                        daf::Size newSup = oldSup > red ? oldSup - red : 0;
                        if (newSup < oldSup) {
                            curSupport[ridx2] = newSup;
                            buckets[newSup].push_back(ridx2);
                            if (newSup < currentLevel) currentLevel = newSup;
                        }
                    };
                    enumerateSubMultisets(comp, r, 0, subKey, subCb);
                    return;
                }
                for (int i = startIdx; i < (int)cids.size(); ++i) {
                    cur.push_back(cids[i]);
                    fillRemaining(i);
                    cur.pop_back();
                }
            };

            // Sort cur (tau's classes) to maintain sorted order for enumeration
            std::sort(cur.begin(), cur.end());

            // We need to enumerate s-multisets that contain tauKey
            // Approach: enumerate ALL s-multisets from cids, filter those ⊇ tauKey
            // More efficient: start from tauKey, add s-r more classes
            cur.clear();
            // Enumerate s-multisets from cids that are ⊇ tauKey (as multiset)
            // Build by picking exactly jc copies of class c (for c in tau), then s-r more
            // Use simpler approach: enumerate all s-multisets, check containment
            std::function<void()> cb = [&]() {
                // Check if cur ⊇ tauKey (as multiset)
                std::unordered_map<daf::Size, int> curCounts;
                for (auto c : cur) curCounts[c]++;
                for (auto &[c, k] : tauCounts)
                    if (curCounts[c] < k) return;
                // Valid: check class sizes
                for (auto &[c, k] : curCounts)
                    if ((int)classSizes[c] < k) return;

                TupleKey key = cur;
                if (deadSTuples.count(key)) return;
                deadSTuples.insert(key);

                auto comp = getComposition(key);
                TupleKey subKey; subKey.reserve(r);
                std::function<void(const TupleKey &)> subCb = [&](const TupleKey &sub) {
                    auto it = rTupleIndex.find(sub);
                    if (it == rTupleIndex.end()) return;
                    daf::Size ridx2 = it->second;
                    if (ridx2 == idx || rPeeled[ridx2]) return;
                    double ext = computeExt(comp, sub, classSizes);
                    if (ext <= 0) return;
                    daf::Size oldSup = curSupport[ridx2];
                    daf::Size red = (daf::Size)ext;
                    daf::Size newSup = oldSup > red ? oldSup - red : 0;
                    if (newSup < oldSup) {
                        curSupport[ridx2] = newSup;
                        buckets[newSup].push_back(ridx2);
                        if (newSup < currentLevel) currentLevel = newSup;
                    }
                };
                enumerateSubMultisets(comp, r, 0, subKey, subCb);
            };
            enumerateMultisets(cids, s, 0, cur, cb);
        }
    }

    auto tStep5End = std::chrono::high_resolution_clock::now();
    auto step5Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStep5).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStart).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling ---" << std::endl;
    std::cout << "  Peeled: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Dead s-tuples: " << deadSTuples.size() << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step5Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;

    return {};
}
