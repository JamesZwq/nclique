//
// Region V2 Fast: V2 with MC-accelerated peeling
//
// Same as V2 but during peeling: when a popped r-tuple belongs to
// only ONE alive MC, delete that entire MC (batch peel all its
// unique tuples + cascade). This avoids per-tuple s-tuple cascade
// for large MCs.
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
extern std::vector<std::vector<daf::Size>> g_maxCliques;

// ============================================================
// Tuple utilities (same as V2)
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

static std::vector<std::pair<daf::Size, int>> getComposition(const TupleKey &t) {
    std::vector<std::pair<daf::Size, int>> comp;
    for (auto c : t) {
        if (!comp.empty() && comp.back().first == c) comp.back().second++;
        else comp.push_back({c, 1});
    }
    return comp;
}

static double computeExt(const std::vector<std::pair<daf::Size, int>> &sigmaComp,
                          const TupleKey &tau,
                          const std::vector<daf::Size> &classSizes) {
    std::unordered_map<daf::Size, int> tauCounts;
    for (auto c : tau) tauCounts[c]++;
    double ext = 1.0;
    for (auto &[cls, mi] : sigmaComp) {
        int ji = tauCounts.count(cls) ? tauCounts[cls] : 0;
        int n = classSizes[cls] - ji, k = mi - ji;
        if (n < k || k < 0) return 0.0;
        ext *= nCr[n][k];
    }
    return ext;
}

// ============================================================
// Main function
// ============================================================

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionV2_Fast(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    auto tStart = std::chrono::high_resolution_clock::now();
    const daf::Size numVertices = edgeGraph.n;
    const daf::Size INVALID = static_cast<daf::Size>(-1);

    // ============================================================
    // Step 1: Build Regions from MaxCliqEnum
    // ============================================================
    std::vector<std::vector<daf::Size>> regionVerts;
    std::vector<std::vector<daf::Size>> vtxMaxPaths(numVertices);
    for (auto &mc : g_maxCliques) {
        if ((int)mc.size() < s) continue;
        daf::Size rid = regionVerts.size();
        regionVerts.push_back(mc);
        for (daf::Size v : mc)
            if (v < numVertices) vtxMaxPaths[v].push_back(rid);
    }
    daf::Size numRegions = regionVerts.size();

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

    // ============================================================
    // Step 3: Enumerate r-tuples and s-tuples (same as V2)
    // ============================================================
    std::unordered_map<TupleKey, daf::Size, TupleHash> rTupleIndex;
    struct RTuple { TupleKey key; daf::Size mult; };
    std::vector<RTuple> rTuples;

    // NEW: track which MCs each r-tuple belongs to
    std::vector<std::vector<daf::Size>> rTupleMCs; // rTupleMCs[tid] = list of MC indices

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
                rTupleMCs.push_back({rid});
            } else {
                // Already exists: add this MC to its list (if not already there)
                auto &mcs = rTupleMCs[it->second];
                if (mcs.empty() || mcs.back() != rid) mcs.push_back(rid);
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

    // NEW: per-MC → list of r-tuples
    std::vector<std::vector<daf::Size>> mcRTuples(numRegions);
    for (daf::Size tid = 0; tid < rTuples.size(); ++tid)
        for (auto mi : rTupleMCs[tid])
            mcRTuples[mi].push_back(tid);

    // Enumerate s-tuples (same as V2)
    std::unordered_map<TupleKey, daf::Size, TupleHash> sTupleIndex;
    struct STuple {
        TupleKey key;
        std::vector<std::pair<daf::Size, double>> incidentRTuples;
        bool alive;
    };
    std::vector<STuple> sTuples;
    std::vector<std::vector<daf::Size>> rToS(rTuples.size());

    // NEW: per s-tuple → which MC it came from
    std::vector<daf::Size> sTupleMC; // which MC this s-tuple belongs to

    {
        TupleKey cur; cur.reserve(s);
        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            auto &cids = classesInRegion[rid];
            if (cids.size() > 500) continue;
            cur.clear();
            std::function<void()> cb = [&]() {
                std::unordered_map<daf::Size, int> counts;
                for (auto c : cur) counts[c]++;
                for (auto &[c, k] : counts)
                    if ((int)classSizes[c] < k) return;

                TupleKey key = cur;
                if (sTupleIndex.count(key)) return;

                auto comp = getComposition(key);
                TupleKey subKey; subKey.reserve(r);
                std::vector<std::pair<daf::Size, double>> incidents;

                std::function<void(const TupleKey &)> subCb = [&](const TupleKey &sub) {
                    auto it = rTupleIndex.find(sub);
                    if (it == rTupleIndex.end()) return;
                    double ext = computeExt(comp, sub, classSizes);
                    if (ext > 0) incidents.push_back({it->second, ext});
                };
                enumerateSubMultisets(comp, r, 0, subKey, subCb);

                if (!incidents.empty()) {
                    daf::Size sid = sTuples.size();
                    sTupleIndex[key] = sid;
                    sTuples.push_back({key, std::move(incidents), true});
                    sTupleMC.push_back(rid);
                    for (auto &[ridx, ext] : sTuples.back().incidentRTuples)
                        rToS[ridx].push_back(sid);
                }
            };
            enumerateMultisets(cids, s, 0, cur, cb);
        }
    }

    auto tStep3 = std::chrono::high_resolution_clock::now();

    // ============================================================
    // Step 4: Compute initial support (same as V2)
    // ============================================================
    std::vector<double> support(rTuples.size(), 0.0);
    for (daf::Size sid = 0; sid < sTuples.size(); ++sid)
        for (auto &[ridx, ext] : sTuples[sid].incidentRTuples)
            support[ridx] += ext;

    double totalRCliques = 0;
    for (daf::Size i = 0; i < rTuples.size(); ++i) totalRCliques += rTuples[i].mult;

    std::cout << "======= Region V2 Fast =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  MCs: " << numRegions << ", Classes: " << numClasses << std::endl;
    std::cout << "  r-tuples: " << rTuples.size() << ", s-tuples: " << sTuples.size() << std::endl;
    std::cout << "  r-cliques: " << std::fixed << std::setprecision(0) << totalRCliques << std::endl;
    std::cout << "  Enum time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(tStep3 - tStart).count()
              << " ms" << std::endl;

    // ============================================================
    // Step 5: MC-Accelerated Cascade Peeling
    // ============================================================
    auto tStep5 = std::chrono::high_resolution_clock::now();

    daf::Size maxSup = 0;
    for (auto &sv : support) maxSup = std::max(maxSup, (daf::Size)sv);

    std::vector<std::vector<daf::Size>> buckets(maxSup + 2);
    std::vector<daf::Size> curSupport(rTuples.size());
    std::vector<bool> rPeeled(rTuples.size(), false);

    // NEW: MC alive tracking
    std::vector<bool> mcAlive(numRegions, true);
    // NEW: per r-tuple, count of alive MCs
    std::vector<int> rTupleAliveMCCount(rTuples.size());
    for (daf::Size tid = 0; tid < rTuples.size(); ++tid)
        rTupleAliveMCCount[tid] = (int)rTupleMCs[tid].size();

    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        curSupport[i] = (daf::Size)support[i];
        buckets[curSupport[i]].push_back(i);
    }

    std::map<daf::Size, int64_t> coreDist;
    daf::Size numPeeled = 0, currentLevel = 0, coreLevel = 0;
    long long batchKills = 0;

    // Helper: peel a single r-tuple (standard V2 cascade)
    auto peelOne = [&](daf::Size idx) {
        if (rPeeled[idx]) return;
        rPeeled[idx] = true;
        numPeeled++;
        coreDist[coreLevel] += rTuples[idx].mult;

        for (daf::Size sid : rToS[idx]) {
            if (!sTuples[sid].alive) continue;
            sTuples[sid].alive = false;
            for (auto &[ridx2, ext] : sTuples[sid].incidentRTuples) {
                if (ridx2 == idx || rPeeled[ridx2]) continue;
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
    };

    // Helper: kill entire MC (batch peel all its unique tuples)
    auto killMC = [&](daf::Size mi) {
        if (!mcAlive[mi]) return;
        mcAlive[mi] = false;
        batchKills++;

        // Mark all of M's s-tuples as dead + cascade
        for (daf::Size tid : mcRTuples[mi]) {
            // Decrement alive MC count for this r-tuple
            rTupleAliveMCCount[tid]--;
        }

        // Peel all r-tuples whose only alive MC was this one
        for (daf::Size tid : mcRTuples[mi]) {
            if (rPeeled[tid]) continue;
            if (rTupleAliveMCCount[tid] == 0) {
                // This r-tuple has no alive MCs left → peel it
                peelOne(tid);
            }
        }
    };

    while (numPeeled < rTuples.size()) {
        while (currentLevel <= maxSup && buckets[currentLevel].empty())
            currentLevel++;
        if (currentLevel > maxSup) break;

        daf::Size idx = buckets[currentLevel].back();
        buckets[currentLevel].pop_back();
        if (rPeeled[idx] || curSupport[idx] != currentLevel) continue;

        coreLevel = std::max(coreLevel, currentLevel);

        // Check: does this r-tuple belong to only 1 alive MC?
        if (rTupleAliveMCCount[idx] == 1) {
            // Find which MC
            daf::Size targetMC = INVALID;
            for (auto mi : rTupleMCs[idx])
                if (mcAlive[mi]) { targetMC = mi; break; }

            if (targetMC != INVALID) {
                // Kill the entire MC
                killMC(targetMC);
                // The popped tuple was peeled inside killMC (if it had 0 alive MCs)
                // If not already peeled (shouldn't happen), peel it now
                if (!rPeeled[idx]) peelOne(idx);
                continue;
            }
        }

        // Normal path: peel single tuple via s-tuple cascade
        peelOne(idx);
    }

    auto tStep5End = std::chrono::high_resolution_clock::now();
    auto step5Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStep5).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStart).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- MC-Accelerated Peeling ---" << std::endl;
    std::cout << "  Peeled: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  MC batch kills: " << batchKills << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    std::cout << "  Peeling time: " << step5Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;

    // Build result
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    for (auto &[c, cnt] : coreDist)
        for (int64_t i = 0; i < cnt; ++i) result.push_back({{}, (double)c});
    return result;
}
