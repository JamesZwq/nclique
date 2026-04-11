//
// Region Tuple V2: General (r,s)-Nucleus Decomposition
//
// V1 (NucleusCoreDecompositionRegion.cpp): correct for s=r+1 only.
// V2: correct for ANY s via bipartite incidence between r-tuples and s-tuples.
//
// Key idea: enumerate BOTH r-tuples and s-tuples (compressed by regions).
// Build incidence. Peel r-tuples; when one is peeled, mark incident s-tuples
// dead and reduce support of other incident r-tuples.
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
// Tuple utilities
// ============================================================

using TupleKey = std::vector<daf::Size>;

struct TupleHash {
    size_t operator()(const TupleKey &t) const noexcept {
        size_t h = t.size();
        for (auto x : t) h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h << 6) + (h >> 2);
        return h;
    }
};

// Enumerate all multisets of given size from sorted class list
// callback(current_multiset) called for each valid multiset
static void enumerateMultisets(
    const std::vector<daf::Size> &classes, int size, int startIdx,
    TupleKey &current, const std::function<void()> &callback)
{
    if ((int)current.size() == size) {
        callback();
        return;
    }
    for (int i = startIdx; i < (int)classes.size(); ++i) {
        current.push_back(classes[i]);
        enumerateMultisets(classes, size, i, current, callback);
        current.pop_back();
    }
}

// Enumerate all r-sub-multisets of an s-multiset sigma
// sigma given as (class, count) pairs. Pick j_i from each class (Σj_i = r).
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
    int remainingClasses = (int)composition.size() - classIdx - 1;
    // j_i can be 0 to min(maxCnt, remaining), but must leave enough for later
    int maxJ = std::min(maxCnt, remaining);
    int minJ = std::max(0, remaining - /* max from remaining classes */ 0);
    // compute max possible from remaining classes
    int maxFromLater = 0;
    for (int k = classIdx + 1; k < (int)composition.size(); ++k)
        maxFromLater += composition[k].second;
    minJ = std::max(0, remaining - maxFromLater);

    for (int j = minJ; j <= maxJ; ++j) {
        for (int q = 0; q < j; ++q) current.push_back(cls);
        enumerateSubMultisets(composition, r, classIdx + 1, current, callback);
        for (int q = 0; q < j; ++q) current.pop_back();
    }
}

// Compute composition of a sorted tuple: (class, count) pairs
static std::vector<std::pair<daf::Size, int>> getComposition(const TupleKey &t) {
    std::vector<std::pair<daf::Size, int>> comp;
    for (auto c : t) {
        if (!comp.empty() && comp.back().first == c) comp.back().second++;
        else comp.push_back({c, 1});
    }
    return comp;
}

// Compute ext(sigma, tau) = Product_i C(|c_i| - j_i, m_i - j_i)
static double computeExt(const std::vector<std::pair<daf::Size, int>> &sigmaComp,
                          const TupleKey &tau,
                          const std::vector<daf::Size> &classSizes) {
    // Build tau composition
    std::unordered_map<daf::Size, int> tauCounts;
    for (auto c : tau) tauCounts[c]++;

    double ext = 1.0;
    for (auto &[cls, mi] : sigmaComp) {
        int ji = tauCounts.count(cls) ? tauCounts[cls] : 0;
        int classSize = classSizes[cls];
        int n = classSize - ji;
        int k = mi - ji;
        if (n < k || k < 0) return 0.0;
        ext *= nCr[n][k];
    }
    return ext;
}

// ============================================================
// Main function
// ============================================================

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionV2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    auto tStart = std::chrono::high_resolution_clock::now();
    const daf::Size numVertices = edgeGraph.n;
    const daf::Size numPaths = tree.adj_list.size();
    const daf::Size INVALID = static_cast<daf::Size>(-1);

    // ============================================================
    // Step 1: Build Maximal-Clique Regions (same as V1)
    // ============================================================

    std::vector<bool> pathValid(numPaths, false);
    daf::Size validPaths = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (tree.adj_list[pid].size() >= r) { pathValid[pid] = true; validPaths++; }
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

    auto tStep1 = std::chrono::high_resolution_clock::now();
    std::cout << "======= Region Tuple V2 (general s) =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Vertices: " << numVertices << ", paths: " << validPaths << std::endl;
    std::cout << "  Maximal cliques: " << numMaximal << ", sub merged: " << numSub
              << " (orphans: " << orphans << ")" << std::endl;

    // ============================================================
    // Step 2: Build Overlap Classes (same as V1)
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
    std::cout << "  Overlap classes: " << numClasses << std::endl;

    // ============================================================
    // Step 3: Enumerate r-tuples and s-tuples
    // ============================================================

    // 3a. Enumerate r-tuples
    std::unordered_map<TupleKey, daf::Size, TupleHash> rTupleIndex;
    struct RTuple { TupleKey key; daf::Size mult; };
    std::vector<RTuple> rTuples;

    {
        TupleKey cur; cur.reserve(r);
        auto addRTuple = [&](daf::Size rid, const TupleKey &key) {
            std::unordered_map<daf::Size, int> counts;
            for (auto c : key) counts[c]++;
            daf::Size mult = 1;
            for (auto &[c, k] : counts) {
                if ((int)classSizes[c] < k) return;
                mult *= (daf::Size)nCr[classSizes[c]][k];
            }
            if (mult == 0) return;
            auto it = rTupleIndex.find(key);
            if (it == rTupleIndex.end()) {
                rTupleIndex[key] = rTuples.size();
                rTuples.push_back({key, mult});
            }
        };

        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            auto &cids = classesInRegion[rid];
            if (cids.size() > 500) continue;
            cur.clear();
            std::function<void()> cb = [&]() { addRTuple(rid, cur); };
            enumerateMultisets(cids, r, 0, cur, cb);
        }
    }

    // 3b. Enumerate s-tuples
    std::unordered_map<TupleKey, daf::Size, TupleHash> sTupleIndex;
    struct STuple {
        TupleKey key;
        std::vector<std::pair<daf::Size, double>> incidentRTuples; // (rTupleIdx, ext)
        bool alive;
    };
    std::vector<STuple> sTuples;

    // Per r-tuple: incident s-tuples
    std::vector<std::vector<daf::Size>> rToS(rTuples.size());

    {
        TupleKey cur; cur.reserve(s);
        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            auto &cids = classesInRegion[rid];
            if (cids.size() > 500) continue;
            cur.clear();
            std::function<void()> cb = [&]() {
                // Check validity: for each class c appearing k times, need |c| >= k
                std::unordered_map<daf::Size, int> counts;
                for (auto c : cur) counts[c]++;
                for (auto &[c, k] : counts)
                    if ((int)classSizes[c] < k) return;

                TupleKey key = cur;
                if (sTupleIndex.count(key)) return; // already enumerated from another region

                // Find all r-sub-tuples
                auto comp = getComposition(key);
                TupleKey subKey; subKey.reserve(r);
                std::vector<std::pair<daf::Size, double>> incidents;

                std::function<void(const TupleKey &)> subCb = [&](const TupleKey &sub) {
                    auto it = rTupleIndex.find(sub);
                    if (it == rTupleIndex.end()) return;
                    double ext = computeExt(comp, sub, classSizes);
                    if (ext > 0)
                        incidents.push_back({it->second, ext});
                };
                enumerateSubMultisets(comp, r, 0, subKey, subCb);

                if (!incidents.empty()) {
                    daf::Size sid = sTuples.size();
                    sTupleIndex[key] = sid;
                    sTuples.push_back({key, std::move(incidents), true});
                    for (auto &[ridx, ext] : sTuples.back().incidentRTuples)
                        rToS[ridx].push_back(sid);
                }
            };
            enumerateMultisets(cids, s, 0, cur, cb);
        }
    }

    auto tStep3 = std::chrono::high_resolution_clock::now();
    auto step3Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep3 - tStep2).count();
    std::cout << "  r-tuples: " << rTuples.size() << ", s-tuples: " << sTuples.size() << std::endl;
    std::cout << "  Enumeration time: " << step3Ms << " ms" << std::endl;

    // ============================================================
    // Step 4: Compute initial support
    // ============================================================

    std::vector<double> support(rTuples.size(), 0.0);
    for (daf::Size sid = 0; sid < sTuples.size(); ++sid) {
        for (auto &[ridx, ext] : sTuples[sid].incidentRTuples)
            support[ridx] += ext;
    }

    // Verification: total support sum
    double totalSupportCPI = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (!pathValid[pid]) continue;
        auto &leaf = tree.adj_list[pid];
        daf::Size h = 0, p = 0;
        for (const auto &node : leaf) { if (node.isPivot) p++; else h++; }
        if ((int)(h + p) >= s && (int)p >= (int)(s - h))
            totalSupportCPI += nCr[p][s - h];
    }
    totalSupportCPI *= nCr[s][r];

    double totalSupportTuples = 0;
    double totalRCliques = 0;
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        totalSupportTuples += rTuples[i].mult * support[i];
        totalRCliques += rTuples[i].mult;
    }

    double relErr = std::abs(totalSupportTuples - totalSupportCPI) / std::max(1.0, totalSupportCPI);
    std::cout << "  r-cliques: " << std::fixed << std::setprecision(0) << totalRCliques << std::endl;
    std::cout << "  Support sum (tuples): " << totalSupportTuples
              << "  (CPI): " << totalSupportCPI << std::endl;
    std::cout << "  VERIFICATION: " << (relErr < 1e-6 ? "PASS" : "MISMATCH") << std::endl;

    // ============================================================
    // Step 5: Cascade Peeling
    // ============================================================

    auto tStep5 = std::chrono::high_resolution_clock::now();

    daf::Size maxSup = 0;
    for (auto &s : support) maxSup = std::max(maxSup, (daf::Size)s);

    std::vector<std::vector<daf::Size>> buckets(maxSup + 2);
    std::vector<daf::Size> curSupport(rTuples.size());
    std::vector<bool> rPeeled(rTuples.size(), false);

    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        curSupport[i] = (daf::Size)support[i];
        buckets[curSupport[i]].push_back(i);
    }

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

        // Cascade: mark incident s-tuples dead, update neighbors
        for (daf::Size sid : rToS[idx]) {
            if (!sTuples[sid].alive) continue;
            sTuples[sid].alive = false;

            for (auto &[ridx2, ext] : sTuples[sid].incidentRTuples) {
                if (ridx2 == idx) continue; // skip self
                if (rPeeled[ridx2]) continue;

                daf::Size oldSup = curSupport[ridx2];
                daf::Size red = (daf::Size)ext;
                daf::Size newSup = oldSup > red ? oldSup - red : 0;
                if (newSup < oldSup) {
                    curSupport[ridx2] = newSup;
                    buckets[newSup].push_back(ridx2);
                    if (newSup < currentLevel) currentLevel = newSup;
                }
            }
        }
    }

    auto tStep5End = std::chrono::high_resolution_clock::now();
    auto step5Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStep5).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStart).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling ---" << std::endl;
    std::cout << "  Peeled: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step5Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;

    return {};
}
