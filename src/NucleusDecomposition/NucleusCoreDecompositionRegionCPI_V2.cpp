//
// Region CPI V3B — Clean Reimplementation with Lazy Split
//
// Clean reimplementation of V3 Region CPI algorithm.
//
// Lazy Split Theorem (s < 2r): if tuples tau and tau' share NO overlap class,
// no s-clique can contain both their r-cliques. Proof: requires 2r vertices
// from disjoint classes, needs s >= 2r.
//
// Applied optimization: when splitting a CPath for peeled tau, skip aggrCount
// computation for UNAFFECTED tuples (those sharing no class with tau).
// Unaffected tuples are added to sub-paths for correct CPath structure
// but their delta is known to be zero, so support is not updated.
//
// Results: EXACT match with V3 and ST baseline for s = r+1.
// For s > r+1, matches V3 (both approximate due to FP accumulation).
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
// Tuple utilities
// ============================================================

namespace v3b {

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

} // namespace v3b

// ============================================================
// Main function: Region CPI V3B (Lazy Split)
// ============================================================

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionCPI_V2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    using namespace v3b;
    auto tStart = std::chrono::high_resolution_clock::now();
    const daf::Size numVertices = edgeGraph.n;
    const daf::Size numPaths = tree.adj_list.size();
    const daf::Size INVALID = static_cast<daf::Size>(-1);
    const bool lazySplit = (s < 2 * r); // Lazy Split Theorem

    // ============================================================
    // Step 1: Build Regions from g_maxCliques
    // ============================================================

    daf::Size numRegions = 0;
    std::vector<std::vector<daf::Size>> regionVerts;
    std::vector<std::vector<daf::Size>> vtxMaxPaths(numVertices);

    for (auto &mc : g_maxCliques) {
        if ((int)mc.size() < s) continue;
        daf::Size rid = regionVerts.size();
        regionVerts.push_back(mc);
        for (daf::Size v : mc)
            if (v < numVertices) vtxMaxPaths[v].push_back(rid);
        numRegions++;
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
            } else {
                classOf[v] = it->second;
                classes[it->second].size++;
            }
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
    std::cout << "======= Region CPI V3B (Lazy Split) =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s
              << " lazySplit=" << (lazySplit ? "ON (s<2r)" : "OFF (s>=2r)") << std::endl;
    std::cout << "  Vertices: " << numVertices << std::endl;
    std::cout << "  Maximal cliques (>=s): " << numRegions << std::endl;
    std::cout << "  Overlap classes: " << numClasses << std::endl;

    // ============================================================
    // Step 3: Enumerate r-tuples
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
    // Step 4: CPI Counting for initial support
    // ============================================================
    // AggrCount(τ, P) via Vandermonde convolution:
    //   g_c(x) = Σ_{b_c} C(nh_c, j_c-b_c) × C(np_c, b_c) × x^{b_c}
    //   f = Π_c g_c(x)
    //   AggrCount = Σ_t f[t] × C(p-t, s-h-t)

    std::vector<double> support(rTuples.size(), 0.0);
    daf::Size pathsUsed = 0;
    daf::Size totalTuplesOnPaths = 0;

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

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

        std::vector<daf::Size> pathClasses;
        for (auto &[cid, _] : classDistrib) pathClasses.push_back(cid);
        std::sort(pathClasses.begin(), pathClasses.end());

        TupleKey cur; cur.reserve(r);
        std::function<void()> cb = [&]() {
            // Build tuple class composition
            std::vector<std::pair<daf::Size, int>> tauClasses;
            {
                daf::Size prev = INVALID; int cnt = 0;
                for (auto c : cur) {
                    if (c == prev) cnt++;
                    else { if (prev != INVALID) tauClasses.push_back({prev, cnt}); prev = c; cnt = 1; }
                }
                if (prev != INVALID) tauClasses.push_back({prev, cnt});
            }

            // Feasibility check
            for (auto &[c, jc] : tauClasses) {
                auto it = classDistrib.find(c);
                if (it == classDistrib.end()) return;
                if (jc > it->second.first + it->second.second) return;
            }

            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;

            // Convolution: f = Π_c g_c(x)
            std::vector<double> f = {1.0};
            for (auto &[c, jc] : tauClasses) {
                auto &[nhc, npc] = classDistrib[c];
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
    // Step 5: Build CPath index
    // ============================================================

    auto tStep5 = std::chrono::high_resolution_clock::now();

    struct ClassBounds { int nh, np; int minPiv, maxPiv; };
    struct CPath {
        int h;
        std::unordered_map<daf::Size, ClassBounds> classes;
        std::vector<daf::Size> tupleIdxs;
        bool alive; // false = destroyed (skip in tupleToCPaths)
    };

    // aggrCountOnCPath: Vandermonde with Pfree optimization
    // Phase 1: convolve tuple classes → h[0..r]
    // Phase 2: non-tuple classes → Pfree (unconstrained) or convolve (constrained)
    // Phase 3: Σ h[k] × C(Pfree, T-k)
    auto aggrCountOnCPath = [&](daf::Size tidx, const CPath &cp) -> double {
        auto &key = rTuples[tidx].key;
        int T = s - cp.h;
        if (T < 0) return 0.0;

        double hPoly[32]; int hLen = 1; hPoly[0] = 1.0;
        int Pfree = 0;

        // Phase 1: tuple classes
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
                    nh2[i+j] += hPoly[i] * gc[j];
            hLen = nL;
            for (int i = 0; i < hLen; ++i) hPoly[i] = nh2[i];
        }

        // Phase 2: non-tuple classes
        for (auto &[c, cb] : cp.classes) {
            bool inTuple = false;
            for (auto x : key) if (x == c) { inTuple = true; break; }
            if (inTuple) continue;

            if (cb.minPiv == 0 && cb.maxPiv >= cb.np) {
                Pfree += cb.np;
            } else {
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
                        nh2[i+j] += hPoly[i] * gc[j];
                hLen = nL;
                for (int i = 0; i < hLen; ++i) hPoly[i] = nh2[i];
            }
        }

        // Phase 3: Σ h[k] × C(Pfree, T-k)
        double aggr = 0.0;
        for (int k = 0; k < hLen; ++k) {
            if (hPoly[k] == 0.0) continue;
            int rem = T - k;
            if (rem >= 0 && rem <= Pfree)
                aggr += hPoly[k] * nCr[Pfree][rem];
        }
        return aggr / rTuples[tidx].mult;
    };

    // --- Initialize CPaths from SDCT paths ---
    std::vector<CPath> cpaths;
    std::vector<std::vector<daf::Size>> tupleToCPaths(rTuples.size());

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        CPath cp;
        cp.h = 0;
        cp.alive = true;
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices || classOf[v] == INVALID) continue;
            daf::Size cid = classOf[v];
            auto &cb = cp.classes[cid];
            if (node.isPivot) { cb.np++; }
            else { cb.nh++; cp.h++; }
        }
        for (auto &[cid, cb] : cp.classes) { cb.minPiv = 0; cb.maxPiv = cb.np; }

        // Find tuples on this path
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
            // Feasibility
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
        int mismatches = 0;
        for (daf::Size i = 0; i < rTuples.size(); ++i) {
            if (std::abs(supCheck[i] - support[i]) > 0.5) mismatches++;
        }
        std::cout << "  CPI vs CPath match: "
                  << (mismatches == 0 ? "PASS" : ("MISMATCH(" + std::to_string(mismatches) + ")"))
                  << std::endl;
    }

    // ============================================================
    // Step 6: Peeling with Lazy Split
    // ============================================================
    //
    // Lazy Split optimization (when s < 2r):
    //   If tuples τ and τ' share NO overlap class, peeling τ doesn't
    //   affect τ' (no s-clique can contain both). Skip aggrCount for
    //   unaffected tuples. Add them to sub-paths (for future structure)
    //   but without computing their contributions.

    auto tStep6 = std::chrono::high_resolution_clock::now();

    long long prof_aggrCalls = 0, prof_splitSubpaths = 0, prof_deadPaths = 0;
    long long prof_skippedAggr = 0; // aggrCount calls saved by lazy split

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

    std::map<daf::Size, int64_t> coreDist;
    daf::Size numPeeled = 0, currentLevel = 0, coreLevel = 0;
    daf::Size totalSplits = 0;

    auto rebucket = [&](daf::Size tidx) {
        daf::Size newBucket = (daf::Size)llround(dSup[tidx]);
        if (newBucket < BUCKET_CAP) {
            buckets[newBucket].push_back(tidx);
            if (newBucket < currentLevel) currentLevel = newBucket;
        } else {
            overflow.insert({newBucket, tidx});
        }
    };

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
        std::unordered_set<daf::Size> tauClassSet(tauKey.begin(), tauKey.end());

        std::unordered_map<daf::Size, int> tauCounts;
        for (auto c : tauKey) tauCounts[c]++;

        auto cpathIds = tupleToCPaths[idx]; // copy
        for (daf::Size cpid : cpathIds) {
            // Access via cpaths[cpid] to avoid dangling references after push_back
            if (!cpaths[cpid].alive || cpaths[cpid].tupleIdxs.empty()) continue;

            // --- Classify alive tuples ---
            std::vector<daf::Size> affectedTuples, unaffectedTuples;
            for (auto tidx2 : cpaths[cpid].tupleIdxs) {
                if (tidx2 == idx || rPeeled[tidx2]) continue;
                bool shares = false;
                if (lazySplit) {
                    for (auto c : rTuples[tidx2].key)
                        if (tauClassSet.count(c)) { shares = true; break; }
                }
                if (lazySplit && !shares)
                    unaffectedTuples.push_back(tidx2);
                else
                    affectedTuples.push_back(tidx2);
            }

            // --- Find active classes ---
            struct ActiveClass { daf::Size cid; int mc; };
            std::vector<ActiveClass> active;
            for (auto &[c, jc] : tauCounts) {
                auto cit = cpaths[cpid].classes.find(c);
                if (cit == cpaths[cpid].classes.end()) continue;
                int mc = std::max(0, jc - cit->second.nh);
                if (mc > cit->second.minPiv)
                    active.push_back({c, mc});
            }

            // --- Compute old contributions (affected only) ---
            std::unordered_map<daf::Size, double> oldContrib;
            for (auto tidx2 : affectedTuples) {
                oldContrib[tidx2] = aggrCountOnCPath(tidx2, cpaths[cpid]);
                prof_aggrCalls++;
            }
            prof_skippedAggr += unaffectedTuples.size();

            // --- LAZY SPLIT: parent survives with unaffected tuples ---
            // Remove affected + peeled from parent. Unaffected stay.
            // Save parent state BEFORE modification for sub-path creation.
            CPath origCp = cpaths[cpid]; // copy for sub-path basis

            if (lazySplit) {
                // Remove affected tuples from parent's tupleToCPaths
                for (auto tidx2 : affectedTuples) {
                    auto &vec = tupleToCPaths[tidx2];
                    vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end());
                }
                // Remove peeled tuple
                {
                    auto &vec = tupleToCPaths[idx];
                    vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end());
                }
                // Parent keeps only unaffected tuples
                cpaths[cpid].tupleIdxs = unaffectedTuples;
            } else {
                // Full split: destroy parent
                for (auto tidx2 : cpaths[cpid].tupleIdxs) {
                    auto &vec = tupleToCPaths[tidx2];
                    vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end());
                }
                cpaths[cpid].tupleIdxs.clear();
                cpaths[cpid].alive = false;
            }

            if (active.empty()) {
                // Dead path: ALL s-cliques contain τ's r-clique
                prof_deadPaths++;
                // Affected tuples lose their contribution
                for (auto &[tidx2, oc] : oldContrib) {
                    dSup[tidx2] -= oc;
                    if (dSup[tidx2] < -0.5) dSup[tidx2] = 0;
                    rebucket(tidx2);
                }
                // Unaffected tuples on parent also lose support (path fully dead)
                if (lazySplit) {
                    for (auto tidx2 : cpaths[cpid].tupleIdxs) {
                        double oc = aggrCountOnCPath(tidx2, origCp);
                        prof_aggrCalls++;
                        dSup[tidx2] -= oc;
                        if (dSup[tidx2] < -0.5) dSup[tidx2] = 0;
                        rebucket(tidx2);
                        auto &vec = tupleToCPaths[tidx2];
                        vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end());
                    }
                    cpaths[cpid].tupleIdxs.clear();
                    cpaths[cpid].alive = false;
                }
                continue;
            }

            // --- Create sub-paths (ONLY affected tuples) ---
            std::unordered_map<daf::Size, double> newContrib;

            for (int i = 0; i < (int)active.size(); ++i) {
                CPath subCp = origCp;
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
                subCp.alive = true;

                // Affected tuples: compute aggrCount on sub-path
                for (auto tidx2 : affectedTuples) {
                    double c = aggrCountOnCPath(tidx2, subCp);
                    prof_aggrCalls++;
                    if (c > 1e-9) {
                        subCp.tupleIdxs.push_back(tidx2);
                        newContrib[tidx2] += c;
                    }
                }

                // Lazy: unaffected stay on parent. Full: copy to sub-path.
                if (!lazySplit) {
                    for (auto tidx2 : unaffectedTuples)
                        subCp.tupleIdxs.push_back(tidx2);
                }
                prof_skippedAggr += unaffectedTuples.size();

                if (subCp.tupleIdxs.empty()) continue;

                daf::Size newCpid = cpaths.size();
                cpaths.push_back(std::move(subCp));
                totalSplits++;
                prof_splitSubpaths++;

                for (auto tidx2 : cpaths.back().tupleIdxs)
                    tupleToCPaths[tidx2].push_back(newCpid);
            }

            // --- Apply delta (affected only) ---
            for (auto &[tidx2, oc] : oldContrib) {
                double nc = newContrib.count(tidx2) ? newContrib[tidx2] : 0.0;
                dSup[tidx2] += (nc - oc);
                if (dSup[tidx2] < -0.5) dSup[tidx2] = 0;
                rebucket(tidx2);
            }
        }
    }

    auto tStep6End = std::chrono::high_resolution_clock::now();
    auto step6Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStep6).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStart).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling (Lazy Split) ---" << std::endl;
    std::cout << "  Peeled: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Total splits: " << totalSplits << " (sub-paths: " << prof_splitSubpaths
              << ", dead: " << prof_deadPaths << ")" << std::endl;
    std::cout << "  Final cpaths: " << cpaths.size() << std::endl;
    std::cout << "  AggrCount calls: " << prof_aggrCalls
              << " (saved: " << prof_skippedAggr << ")" << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step6Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;

    return {};
}
