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
    // Step 6: Strategy A Peeling (CPI-based, NO s-tuples)
    // ============================================================
    // When τ is peeled: for each path P containing τ, for each other τ' on P:
    // Δ = "new dying s-cliques on P containing τ' due to τ"
    //
    // For correctness: we track which tuples have been peeled per path.
    // A s-clique on P is dead iff it contains an r-clique from any peeled tuple on P.
    //
    // Simple approach: maintain per-tuple alive_support (initially from CPI).
    // When τ is peeled: iterate over paths, compute delta for each neighbor τ'.
    //
    // Delta computation: on path P, the dying s-cliques for τ' are those
    // containing both τ and τ' r-cliques. We compute this via the
    // constrained allocation formula (AggrCount for the "union requirement").

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

    // Per-path: track peeled tuples' requirement vectors for dead s-clique tracking
    // For now: simple approach — recompute alive support per (tuple, path) on demand
    // using the set of peeled tuples on that path.
    //
    // Track per path: list of peeled requirement vectors
    struct ReqVec { std::unordered_map<daf::Size, int> mR; }; // m_c = max(0, j_c - nh_c)
    std::vector<std::vector<ReqVec>> pathDeadReqs(numPaths); // per path: dead reqs

    std::map<daf::Size, int64_t> coreDist;
    daf::Size numPeeled = 0, currentLevel = 0, coreLevel = 0;

    // Helper: compute number of s-cliques on path P containing C_r' from τ'
    // that ALSO satisfy a set of requirement vectors (i.e., contain r-cliques from dead tuples)
    // Uses inclusion-exclusion over the requirement vectors.
    auto countSatisfying = [&](daf::Size pid, daf::Size tidx,
                               const std::vector<ReqVec> &reqs) -> double {
        auto &pd = pathDataVec[pid];
        auto &tauKey = rTuples[tidx].key;

        // Build tau's class distribution on this path
        std::unordered_map<daf::Size, int> tauCounts;
        for (auto c : tauKey) tauCounts[c]++;

        // For each subset S of reqs: compute the number of s-cliques containing τ'
        // AND satisfying all requirements in S.
        // Combined requirement: for each class c, need max(q'_c, max_{i∈S} m_c^(i)) pivots
        // where q'_c = max(0, j'_c - nh_c) = min pivots τ' needs from c

        double total = 0;
        int nReqs = reqs.size();
        if (nReqs > 20) nReqs = 20; // cap for safety

        for (int mask = 0; mask < (1 << nReqs); ++mask) {
            // Combined requirement: for each class c, the minimum pivot count
            std::unordered_map<daf::Size, int> combined;

            // Start with τ' requirements
            for (auto &[c, jc] : tauCounts) {
                int nhc = pd.classDistrib.count(c) ? pd.classDistrib[c].first : 0;
                combined[c] = std::max(0, jc - nhc);
            }

            // Add requirements from selected dead tuples
            int bits = __builtin_popcount(mask);
            for (int b = 0; b < nReqs; ++b) {
                if (!(mask & (1 << b))) continue;
                for (auto &[c, mc] : reqs[b].mR) {
                    combined[c] = std::max(combined[c], mc);
                }
            }

            // Count s-cliques with n_Q^c ≥ combined[c] for all c,
            // and Σ n_Q^c = s - h, and n_Q^c ≤ np_c
            // This is the constrained allocation formula
            // Using convolution: for each class c, choices = np_c - combined[c] extra pivots
            // (after fixing combined[c] mandatory pivots)

            int mandatoryTotal = 0;
            bool feasible = true;
            std::vector<std::pair<daf::Size, int>> allocClasses; // (class, max_extra)
            for (auto &[c, minPiv] : combined) {
                int npc = pd.classDistrib.count(c) ? pd.classDistrib[c].second : 0;
                if (minPiv > npc) { feasible = false; break; }
                mandatoryTotal += minPiv;
                allocClasses.push_back({c, npc - minPiv});
            }
            // Also include classes on P not in combined (they contribute 0 mandatory but have available pivots)
            if (feasible) {
                for (auto &[c, dp] : pd.classDistrib) {
                    if (combined.count(c)) continue;
                    allocClasses.push_back({c, dp.second}); // all pivots available as extra
                }
            }

            int extraNeeded = (s - pd.h) - mandatoryTotal;
            if (!feasible || extraNeeded < 0) continue;

            // Count allocations: distribute extraNeeded among allocClasses
            // Each class can contribute 0..max_extra
            // Count = coefficient of x^extraNeeded in Π (1 + x + ... + x^{max_extra_c})
            // = coefficient of x^extraNeeded in Π Σ_{e=0}^{max_extra_c} C(1, 1)^e ... nah

            // Actually: we're choosing specific pivot VERTICES, not just counts.
            // For class c: we need combined[c] + e_c pivots total, choosing from np_c.
            // Ways = C(np_c, combined[c] + e_c). With e_c extra.
            // Total = Σ_{e_c ≥ 0, Σe_c = extraNeeded} Π C(np_c, combined_c + e_c)

            // Compute via convolution
            std::vector<double> f = {1.0};
            for (auto &[c, maxE] : allocClasses) {
                int npc = pd.classDistrib.count(c) ? pd.classDistrib[c].second : 0;
                int minPiv = combined.count(c) ? combined[c] : 0;
                std::vector<double> gc(std::min(maxE, extraNeeded) + 1, 0.0);
                for (int e = 0; e < (int)gc.size(); ++e) {
                    int choose = minPiv + e;
                    if (choose <= npc)
                        gc[e] = nCr[npc][choose];
                }
                std::vector<double> newf(f.size() + gc.size() - 1, 0.0);
                for (int i = 0; i < (int)f.size(); ++i)
                    for (int j = 0; j < (int)gc.size(); ++j)
                        newf[i + j] += f[i] * gc[j];
                f = std::move(newf);
            }

            double count = (extraNeeded < (int)f.size()) ? f[extraNeeded] : 0.0;

            // Inclusion-exclusion sign
            total += (bits % 2 == 0 ? 1 : -1) * count;
        }

        return total;
    };

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

        // Build requirement vector for τ
        auto &tauKey = rTuples[idx].key;

        for (daf::Size pid : tupleToPathIds[idx]) {
            auto &pd = pathDataVec[pid];

            // Compute τ's requirement vector on this path
            ReqVec newReq;
            std::unordered_map<daf::Size, int> tauCounts;
            for (auto c : tauKey) tauCounts[c]++;
            for (auto &[c, jc] : tauCounts) {
                int nhc = pd.classDistrib.count(c) ? pd.classDistrib[c].first : 0;
                int mc = std::max(0, jc - nhc);
                if (mc > 0) newReq.mR[c] = mc;
            }

            // For each other tuple τ' on this path: compute Δ
            // Δ = (s-cliques containing τ' AND satisfying newReq)
            //   - (s-cliques containing τ' AND satisfying newReq AND any previous dead req)
            // = countSatisfying(pid, tidx', {newReq})
            //   - countSatisfying(pid, tidx', {newReq, prev_reqs...})
            // But this is just: countSatisfying with {all_reqs_including_new}
            //                  - countSatisfying with {prev_reqs_only}
            // = Δ(alive before) - Δ(alive after) = net new deaths for τ'

            // Simpler: alive_P(τ') was maintained. After peeling τ:
            // new alive_P(τ') = countSatisfying(pid, tidx', {} means "not satisfying ANY dead req")
            // Hmm, countSatisfying with empty reqs = total s-cliques containing τ' = initial support from P

            // Actually the simplest correct approach:
            // alive_after = count of s-cliques on P containing τ' and NOT containing any dead tuple's r-clique
            // Δ = alive_before - alive_after

            auto &deadReqs = pathDeadReqs[pid];

            for (auto &pti : pd.tuples) {
                if (pti.tidx == idx || rPeeled[pti.tidx]) continue;

                // alive_before = support from this path (before peeling τ)
                // = countSatisfying with complement of deadReqs (alive = NOT satisfying any dead req)
                // = total - countSatisfying(deadReqs)

                // After adding newReq to deadReqs:
                // alive_after = total - countSatisfying(deadReqs ∪ {newReq})

                // Δ = alive_before - alive_after
                //   = countSatisfying(deadReqs ∪ {newReq}) - countSatisfying(deadReqs)

                // This is exactly: (s-cliques satisfying newReq) that DON'T satisfy any previous dead req
                // = countSatisfying({newReq}) - countSatisfying({newReq} ∩ some dead req)

                // Using inclusion-exclusion over deadReqs ∪ {newReq}:

                std::vector<ReqVec> allReqs = deadReqs;
                allReqs.push_back(newReq);

                double afterCount = countSatisfying(pid, pti.tidx, allReqs);
                double beforeCount = countSatisfying(pid, pti.tidx, deadReqs);
                double delta = afterCount - beforeCount;

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

            // Add newReq to this path's dead list
            deadReqs.push_back(newReq);
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
