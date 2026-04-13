//
// V4: Region-based (r,s)-Nucleus Decomposition
//
// Phase 1: Iteratively delete MCs with private vertices.
//          core(unique tuples) = C(|M|-r, s-r). Update MC lists.
// Phase 2: Remaining tuples — compute IE support, assign core.
//
// No CPI. No cliqueIndex. No BK. Pure MC/region logic.
//

#include "NCliqueCoreDecomposition.h"
#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <unordered_map>
#include <vector>
#include <set>

extern double nCr[1001][401];
extern std::vector<std::vector<daf::Size>> g_maxCliques;

using TupleKey = std::vector<daf::Size>;
struct TupleHash {
    size_t operator()(const TupleKey &t) const noexcept {
        size_t h = t.size();
        for (auto x : t) h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h << 6) + (h >> 2);
        return h;
    }
};

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex)
{
    auto tStart = std::chrono::high_resolution_clock::now();
    const daf::Size numVertices = edgeGraph.getGraphNodeSize();
    const daf::Size INVALID = std::numeric_limits<daf::Size>::max();

    // ================================================================
    // Build MCs + vertex MC membership
    // ================================================================
    std::vector<std::vector<daf::Size>> mcs;
    for (auto &mc : g_maxCliques)
        if ((int)mc.size() >= s) mcs.push_back(mc);
    daf::Size numMC = mcs.size();

    std::vector<bool> mcAlive(numMC, true);
    // Per-vertex: set of alive MC indices
    std::vector<std::vector<daf::Size>> vtxMCs(numVertices);
    for (daf::Size mi = 0; mi < numMC; ++mi)
        for (auto v : mcs[mi])
            if (v < numVertices) vtxMCs[v].push_back(mi);

    auto mcCoreVal = [&](daf::Size mi) -> double {
        int n = (int)mcs[mi].size() - (int)r, k = (int)s - (int)r;
        return (n >= k && k >= 0) ? nCr[n][k] : 0.0;
    };

    // ================================================================
    // Phase 1: Iteratively delete MCs with private vertices
    // ================================================================
    // Private vertex of MC M = vertex v where vtxMCs[v] == {M}
    //
    // When M is deleted:
    //   - Tuples containing M's private class: core = C(|M|-r, s-r)
    //   - Count = C(|M|, r) - C(|M| - #priv, r)
    //   - Remove M from all vertices' MC lists
    //   - Some vertices may now become private to another MC → cascade

    std::map<double, double> coreDist; // core → r-clique count

    auto getPrivCount = [&](daf::Size mi) -> int {
        int p = 0;
        for (auto v : mcs[mi])
            if (v < numVertices && vtxMCs[v].size() == 1 && vtxMCs[v][0] == mi) p++;
        return p;
    };

    // Priority queue: smallest C(|M|-r, s-r) first
    using PQE = std::pair<double, daf::Size>;
    std::priority_queue<PQE, std::vector<PQE>, std::greater<>> pq;
    for (daf::Size mi = 0; mi < numMC; ++mi)
        if (getPrivCount(mi) > 0) pq.push({mcCoreVal(mi), mi});

    daf::Size phase1MCs = 0;
    double phase1RCliques = 0;

    while (!pq.empty()) {
        auto [val, mi] = pq.top(); pq.pop();
        if (!mcAlive[mi]) continue;
        int priv = getPrivCount(mi);
        if (priv == 0) continue;

        // Assign core
        double core = val;
        int M = (int)mcs[mi].size();
        double withPriv = nCr[M][r] - nCr[M - priv][r];
        if (withPriv > 0) {
            coreDist[core] += withPriv;
            phase1RCliques += withPriv;
        }

        // Delete MC
        mcAlive[mi] = false;
        phase1MCs++;

        // Update vertex MC lists + cascade
        std::set<daf::Size> checkMCs;
        for (auto v : mcs[mi]) {
            if (v >= numVertices) continue;
            auto &ml = vtxMCs[v];
            ml.erase(std::remove(ml.begin(), ml.end(), mi), ml.end());
            // If vertex now has exactly 1 MC: that MC gains a private vertex
            if (ml.size() == 1 && mcAlive[ml[0]])
                checkMCs.insert(ml[0]);
        }
        for (auto m2 : checkMCs)
            if (getPrivCount(m2) > 0)
                pq.push({mcCoreVal(m2), m2});
    }

    auto tPhase1 = std::chrono::high_resolution_clock::now();

    // ================================================================
    // Phase 2: Remaining MCs (no private vertices)
    // ================================================================
    // Build classes + tuples from remaining alive MCs.
    // Compute support via IE. Assign core = support.
    // (Simple: no cascading. Correct for symmetric cases like K_{2,2,2,2}.)

    std::vector<daf::Size> remainMCs;
    for (daf::Size mi = 0; mi < numMC; ++mi)
        if (mcAlive[mi]) remainMCs.push_back(mi);

    double phase2RCliques = 0;

    if (!remainMCs.empty()) {
        // Build classes from remaining alive profiles
        using Profile = std::vector<daf::Size>;
        struct PH { size_t operator()(const Profile &p) const noexcept {
            size_t h=p.size(); for(auto x:p) h^=std::hash<daf::Size>()(x)+0x9e3779b9ULL+(h<<6)+(h>>2); return h; }};
        std::unordered_map<Profile, daf::Size, PH> pToC;
        std::vector<daf::Size> classOf(numVertices, INVALID);
        std::vector<daf::Size> classSizes;

        for (daf::Size v = 0; v < numVertices; ++v) {
            if (vtxMCs[v].empty()) continue;
            auto it = pToC.find(vtxMCs[v]);
            if (it == pToC.end()) { daf::Size c = classSizes.size(); pToC[vtxMCs[v]] = c; classSizes.push_back(0); }
            classOf[v] = pToC[vtxMCs[v]]; classSizes[classOf[v]]++;
        }

        // MC class sizes
        std::vector<std::unordered_map<daf::Size, int>> mcCS(numMC);
        std::vector<std::vector<daf::Size>> mcCL(numMC);
        for (auto mi : remainMCs) {
            for (auto v : mcs[mi]) if (v < numVertices && classOf[v] != INVALID) mcCS[mi][classOf[v]]++;
            for (auto &[c, _] : mcCS[mi]) mcCL[mi].push_back(c);
            std::sort(mcCL[mi].begin(), mcCL[mi].end());
        }

        // Enumerate tuples
        struct TI { TupleKey key; daf::Size mult; };
        std::vector<TI> tuples;
        std::unordered_map<TupleKey, daf::Size, TupleHash> tidx;
        { TupleKey cur; cur.reserve(r);
          std::function<void(const std::vector<daf::Size>&, int)> en;
          en = [&](const std::vector<daf::Size> &cl, int st) {
            if ((int)cur.size() == r) { if (tidx.count(cur)) return;
              std::unordered_map<daf::Size, int> cnt; for (auto c : cur) cnt[c]++;
              daf::Size mult = 1; for (auto &[c, k] : cnt) { if ((int)classSizes[c] < k) return;
                mult *= (daf::Size)(nCr[classSizes[c]][k] + 0.5); }
              if (!mult) return; tidx[cur] = tuples.size(); tuples.push_back({cur, mult}); return; }
            for (int i = st; i < (int)cl.size(); ++i) { cur.push_back(cl[i]); en(cl, i); cur.pop_back(); }
          };
          for (auto mi : remainMCs) { cur.clear(); en(mcCL[mi], 0); }
        }

        // Tuple → MCs + IE support
        auto mcInterSz = [&](const std::vector<daf::Size> &ms) -> int {
            if (ms.empty()) return 0; if (ms.size() == 1) return (int)mcs[ms[0]].size();
            auto cur = mcs[ms[0]]; for (size_t i = 1; i < ms.size(); ++i) {
                std::vector<daf::Size> nxt; std::set_intersection(cur.begin(), cur.end(),
                    mcs[ms[i]].begin(), mcs[ms[i]].end(), std::back_inserter(nxt)); cur = std::move(nxt); }
            return (int)cur.size();
        };

        for (auto &t : tuples) {
            // Find MCs containing this tuple
            std::unordered_map<daf::Size, int> cnt; for (auto c : t.key) cnt[c]++;
            std::vector<daf::Size> tMCs;
            for (auto mi : remainMCs) {
                bool ok = true;
                for (auto &[c, k] : cnt) { auto it = mcCS[mi].find(c); if (it == mcCS[mi].end() || it->second < k) { ok = false; break; } }
                if (ok) tMCs.push_back(mi);
            }

            // IE support
            int p = (int)tMCs.size(); double sup = 0;
            for (int mask = 1; mask < (1 << p); ++mask) {
                std::vector<daf::Size> sub; for (int i = 0; i < p; ++i) if (mask & (1 << i)) sub.push_back(tMCs[i]);
                int isz = mcInterSz(sub); int n = isz - (int)r, k = (int)s - (int)r;
                double v = (n >= k && k >= 0) ? nCr[n][k] : 0.0;
                sup += (__builtin_popcount(mask) % 2 == 1 ? 1 : -1) * v;
            }
            sup = std::max(0.0, sup);

            coreDist[sup] += t.mult;
            phase2RCliques += t.mult;
        }
    }

    auto tEnd = std::chrono::high_resolution_clock::now();

    // Output
    printf("======= V4 Region-Based =======\n");
    printf("  r=%d s=%d, MCs: %zu\n", (int)r, (int)s, (size_t)numMC);
    printf("  Phase 1: %zu MCs peeled, %.0f r-cliques, %lld ms\n",
        (size_t)phase1MCs, phase1RCliques,
        std::chrono::duration_cast<std::chrono::milliseconds>(tPhase1 - tStart).count());
    printf("  Phase 2: %zu MCs remain, %.0f r-cliques, %lld ms\n",
        remainMCs.size(), phase2RCliques,
        std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tPhase1).count());
    printf("  Total: %lld ms\n",
        std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count());

    double maxCore = 0;
    for (auto &[c, _] : coreDist) maxCore = std::max(maxCore, c);
    printf("  Max core: %.0f\n", maxCore);
    for (auto &[c, cnt] : coreDist) printf("  core=%.0f count=%.0f\n", c, cnt);

    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    for (auto &[c, cnt] : coreDist) result.push_back({{}, c});
    return result;
}
