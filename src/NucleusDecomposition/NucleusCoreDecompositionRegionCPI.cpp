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
    auto aggrCountOnCPath = [&](daf::Size tidx, const CPath &cp) -> double {
        auto &key = rTuples[tidx].key;
        std::unordered_map<daf::Size, int> counts;
        for (auto c : key) counts[c]++;

        int p = cp.totalP();

        // Convolution: for each class c, generate g_c(x)
        // b_c = pivots chosen from c for the r-clique
        // b_c ranges: max(0, j_c - nh_c) ≤ b_c ≤ min(j_c, maxPiv_c)
        // AND b_c ≤ np_c (available pivots)
        // AND j_c - b_c ≤ nh_c
        // Ways per b_c: C(nh_c, j_c - b_c) × C(min(np_c, maxPiv_c), b_c) ... hmm

        // Actually: constrained path restricts HOW MANY pivots from c are in Q_s.
        // For the s-clique: minPiv_c ≤ n_Q^c ≤ maxPiv_c.
        // For the r-clique C_r' from τ': C_r' uses j_c vertices from c,
        //   of which (j_c - b_c) from holds and b_c from pivots.
        // C_r' ⊆ S means the b_c pivots are among Q_s's pivots from c.

        // The s-clique has n_Q^c pivots from c (with minPiv ≤ n_Q^c ≤ maxPiv).
        // C_r' uses b_c of those pivots.
        // Ways: C(nh_c, j_c-b_c) × C(n_Q^c, b_c) × (remaining pivots chosen freely)

        // We need to sum over b_c AND n_Q^c. But AggrCount sums over all s-cliques
        // containing all possible C_r' from τ'. Let me reformulate.

        // Per specific C_r' (with b_c pivots from each c):
        //   q'_c = b_c (pivots from c in C_r')
        //   extra pivots from c: e_c = n_Q^c - b_c, with max(minPiv_c - b_c, 0) ≤ e_c ≤ maxPiv_c - b_c
        //   ways to choose extra: C(np_c - b_c, e_c)
        //   total extra: Σ e_c = (s - h) - Σ b_c

        // To get AggrCount (total over all C_r' in τ'):
        //   Σ_{b_c} [Π C(nh_c, j_c-b_c) C(np_c, b_c)] × [s-cliques containing this C_r' on cp]
        //   where s-cliques containing C_r' = constrained allocation with bounds

        // This is a double summation: over b_c (for the r-clique) and e_c (for extra pivots).
        // Total pivots from c in s-clique = b_c + e_c, constrained by [minPiv, maxPiv].

        // Combined: for each class c, total pivots t_c = b_c + e_c.
        // minPiv_c ≤ t_c ≤ maxPiv_c
        // Ways to split t_c into (b_c from r-clique, e_c extra): C(nh_c, j_c-b_c) × C(np_c, t_c)
        //   ... wait, this isn't right. The r-clique picks b_c specific pivots, then the
        //   remaining t_c - b_c are chosen from the remaining np_c - b_c.
        //   So ways = C(nh_c, j_c-b_c) × C(np_c, b_c) × C(np_c-b_c, t_c-b_c)
        //   Hmm, but C(np_c, b_c) × C(np_c-b_c, t_c-b_c) = C(np_c, t_c) × C(t_c, b_c)

        // Actually for AggrCount: we're counting (C_r', S) pairs where C_r' ∈ τ' and S is
        // an s-clique on cp containing C_r'. Each S counted once per C_r' it contains.
        // Then divide by mult(τ') to get per-C_r' count.

        // Per C_r' from τ' (with b_c pivots from c):
        //   s-cliques containing C_r' on cp = Σ_{e_c} Π C(np_c - b_c, e_c)
        //   where for each c: max(0, minPiv_c - b_c) ≤ e_c ≤ maxPiv_c - b_c
        //   and Σ(b_c + e_c) = s - h, i.e., Σe_c = (s-h) - Σb_c

        // Total AggrCount = Σ_{b_c} [Π C(nh_c, j_c-b_c) × C(np_c, b_c)] × count(e_c's)
        // Summed over all valid b_c and e_c.

        // Merge b_c + e_c into t_c: t_c = total pivots from c in s-clique.
        // For each c: t_c ranges from max(minPiv_c, max(0, j_c - nh_c)) to min(maxPiv_c, np_c)
        // And within t_c: b_c ranges from max(0, j_c-nh_c) to min(j_c, t_c)
        // Ways = Σ_{b_c} C(nh_c, j_c-b_c) × C(np_c, b_c) × C(np_c-b_c, t_c-b_c)

        // Using identity: C(np_c, b_c) × C(np_c-b_c, t_c-b_c) = C(np_c, t_c) × C(t_c, b_c)
        // So: ways = C(np_c, t_c) × Σ_{b_c} C(nh_c, j_c-b_c) × C(t_c, b_c)
        //         = C(np_c, t_c) × C(nh_c + t_c, j_c)  [Vandermonde convolution]

        // So for each class c:
        //   g_c[t_c] = C(np_c, t_c) × C(nh_c + t_c, j_c)
        //   for t_c in [max(minPiv_c, max(0, j_c-nh_c)), min(maxPiv_c, np_c)]

        // AggrCount = Σ_{Σt_c = s-h} Π g_c[t_c]
        // = [x^{s-h}] Π_c (Σ_{t_c} g_c[t_c] x^{t_c})

        // Then divide by mult(τ') for per-C_r' count.

        // This is clean! The generating function per class is g_c(x).

        std::vector<double> f = {1.0};
        int targetTotal = s - cp.h;

        for (auto &[c, cb] : cp.classes) {
            int jc = counts.count(c) ? counts[c] : 0;
            int tMin = std::max(cb.minPiv, std::max(0, jc - cb.nh));
            int tMax = std::min(cb.maxPiv, cb.np);
            if (jc > cb.nh + tMax) return 0.0; // infeasible: not enough vertices for j_c
            if (tMin > tMax) {
                // This class can't contribute enough pivots
                // If jc > 0 and tMax < tMin: tuple can't be on this path
                if (jc > 0) return 0.0;
                // If jc == 0: class just contributes to total pivots
                tMin = std::max(tMin, 0); // but bounded by minPiv
                if (tMin > tMax) return 0.0;
            }

            std::vector<double> gc(tMax + 1, 0.0);
            for (int tc = tMin; tc <= tMax; ++tc) {
                if (cb.nh + tc < jc) continue; // C(nh+tc, jc) invalid
                gc[tc] = nCr[cb.np][tc] * nCr[cb.nh + tc][jc];
            }

            std::vector<double> newf(f.size() + gc.size() - 1, 0.0);
            for (int i = 0; i < (int)f.size(); ++i)
                for (int j = 0; j < (int)gc.size(); ++j)
                    newf[i + j] += f[i] * gc[j];
            f = std::move(newf);
        }

        double aggr = (targetTotal >= 0 && targetTotal < (int)f.size()) ? f[targetTotal] : 0.0;
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

    // --- Peeling with Analytical Split ---
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

    std::map<daf::Size, int64_t> coreDist;
    daf::Size numPeeled = 0, currentLevel = 0, coreLevel = 0;
    daf::Size totalSplits = 0;

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
        std::unordered_map<daf::Size, int> tauCounts;
        for (auto c : tauKey) tauCounts[c]++;

        // For each constrained path containing τ:
        auto cpathIds = tupleToCPaths[idx]; // copy (will be modified)
        for (daf::Size cpid : cpathIds) {
            auto &cp = cpaths[cpid];
            if (cp.tupleIdxs.empty()) continue;

            // Step A: Subtract old contributions from this cpath
            for (auto tidx : cp.tupleIdxs) {
                if (rPeeled[tidx] && tidx != idx) continue;
                if (tidx == idx) continue;
                double oldContrib = aggrCountOnCPath(tidx, cp);
                if (oldContrib > 0.5) {
                    daf::Size oldSup = curSupport[tidx];
                    daf::Size red = (daf::Size)(oldContrib + 0.5);
                    curSupport[tidx] = oldSup > red ? oldSup - red : 0;
                }
            }

            // Step B: Compute τ's requirement and find active classes
            struct ActiveClass { daf::Size cid; int mc; };
            std::vector<ActiveClass> active;
            for (auto &[c, jc] : tauCounts) {
                auto cit = cp.classes.find(c);
                if (cit == cp.classes.end()) continue;
                int mc = std::max(0, jc - cit->second.nh);
                if (mc > cit->second.minPiv) // m_c exceeds current lower bound
                    active.push_back({c, mc});
            }

            // Build list of alive tuples (excluding peeled ones)
            std::vector<daf::Size> aliveTuples;
            for (auto tidx : cp.tupleIdxs)
                if (!rPeeled[tidx]) aliveTuples.push_back(tidx);

            // Remove old cpath from all tuples' cpath lists
            for (auto tidx : cp.tupleIdxs) {
                auto &vec = tupleToCPaths[tidx];
                vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end());
            }
            cp.tupleIdxs.clear();

            if (active.empty()) {
                // All s-cliques on cp contain τ's r-clique → entire cp is dead
                // (old contributions already subtracted; don't add anything back)
            } else {
                // Step C: Split into |active| disjoint sub-paths
                for (int i = 0; i < (int)active.size(); ++i) {
                    CPath subCp = cp; // copy parent's bounds
                    bool feasible = true;

                    // For j < i: tighten lower bound to m_{c_j}
                    for (int j = 0; j < i; ++j) {
                        auto &cb = subCp.classes[active[j].cid];
                        cb.minPiv = std::max(cb.minPiv, active[j].mc);
                        if (cb.minPiv > cb.maxPiv) { feasible = false; break; }
                    }
                    if (!feasible) continue;

                    // For c_i: tighten upper bound to m_{c_i} - 1
                    {
                        auto &cb = subCp.classes[active[i].cid];
                        cb.maxPiv = std::min(cb.maxPiv, active[i].mc - 1);
                        if (cb.minPiv > cb.maxPiv) continue; // infeasible
                    }

                    // Check total pivots feasibility
                    int minTotal = 0, maxTotal = 0;
                    for (auto &[_, cb] : subCp.classes) {
                        minTotal += cb.minPiv;
                        maxTotal += cb.maxPiv;
                    }
                    int need = s - subCp.h;
                    if (need < minTotal || need > maxTotal) continue;

                    // Find alive tuples on this sub-path
                    subCp.tupleIdxs.clear();
                    for (auto tidx : aliveTuples) {
                        double contrib = aggrCountOnCPath(tidx, subCp);
                        if (contrib > 1e-9) subCp.tupleIdxs.push_back(tidx);
                    }
                    if (subCp.tupleIdxs.empty()) continue;

                    // Add sub-path and register with tuples
                    daf::Size newCpid = cpaths.size();
                    cpaths.push_back(std::move(subCp));
                    totalSplits++;

                    // Step D: Add new contributions from sub-path
                    for (auto tidx : cpaths.back().tupleIdxs) {
                        tupleToCPaths[tidx].push_back(newCpid);
                        double newContrib = aggrCountOnCPath(tidx, cpaths.back());
                        if (newContrib > 0.5) {
                            curSupport[tidx] += (daf::Size)(newContrib + 0.5);
                            buckets[curSupport[tidx]].push_back(tidx);
                        }
                    }
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
    std::cout << "  Total splits: " << totalSplits << std::endl;
    std::cout << "  Final cpaths: " << cpaths.size() << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step6Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;

    return {};
}
