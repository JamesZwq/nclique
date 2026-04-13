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
    // Phase 2: Remaining MCs — brute force peeling
    // ================================================================
    // Enumerate ALL s-cliques in remaining MCs.
    // Track alive status. Standard peeling with cascading.

    std::vector<daf::Size> remainMCs;
    for (daf::Size mi = 0; mi < numMC; ++mi)
        if (mcAlive[mi]) remainMCs.push_back(mi);

    double phase2RCliques = 0;

    if (!remainMCs.empty()) {
        // Enumerate all s-cliques from remaining MCs (as sorted vertex vectors)
        std::set<std::vector<daf::Size>> allSCliques;
        for (auto mi : remainMCs) {
            auto &mc = mcs[mi];
            int n = (int)mc.size();
            // Enumerate all s-subsets of mc
            std::vector<int> idx(s);
            for (int i = 0; i < s; ++i) idx[i] = i;
            while (true) {
                std::vector<daf::Size> sc(s);
                for (int i = 0; i < s; ++i) sc[i] = mc[idx[i]];
                std::sort(sc.begin(), sc.end());
                allSCliques.insert(sc);
                // Next combination
                int i = s - 1;
                while (i >= 0 && idx[i] == n - s + i) i--;
                if (i < 0) break;
                idx[i]++;
                for (int j = i + 1; j < s; ++j) idx[j] = idx[j-1] + 1;
            }
        }

        // Enumerate all r-cliques from s-cliques
        std::map<std::vector<daf::Size>, int> rCliqueSupport; // r-clique → support count
        std::map<std::vector<daf::Size>, std::vector<int>> rCliqueToSCliques; // r-clique → s-clique indices

        std::vector<std::vector<daf::Size>> sCliqueList(allSCliques.begin(), allSCliques.end());
        std::vector<bool> sAlive(sCliqueList.size(), true);

        for (int si = 0; si < (int)sCliqueList.size(); ++si) {
            auto &sc = sCliqueList[si];
            // Enumerate r-subsets
            std::vector<int> idx(r);
            for (int i = 0; i < r; ++i) idx[i] = i;
            while (true) {
                std::vector<daf::Size> rc(r);
                for (int i = 0; i < r; ++i) rc[i] = sc[idx[i]];
                rCliqueSupport[rc]++;
                rCliqueToSCliques[rc].push_back(si);
                int i = r - 1;
                while (i >= 0 && idx[i] == s - r + i) i--;
                if (i < 0) break;
                idx[i]++;
                for (int j = i + 1; j < r; ++j) idx[j] = idx[j-1] + 1;
            }
        }

        // Peeling
        std::map<std::vector<daf::Size>, double> rCliqueCore;
        double minCore2 = 0;

        while (!rCliqueSupport.empty()) {
            // Find min support
            double minSup = 1e18;
            for (auto &[rc, sup] : rCliqueSupport)
                if (sup < minSup) minSup = sup;
            minCore2 = std::max(minCore2, minSup);

            // Collect all r-cliques at min support
            std::vector<std::vector<daf::Size>> batch;
            for (auto &[rc, sup] : rCliqueSupport)
                if (sup <= minCore2) batch.push_back(rc);

            // Peel batch
            for (auto &rc : batch) {
                rCliqueCore[rc] = minCore2;
                // Kill s-cliques containing this r-clique
                for (auto si : rCliqueToSCliques[rc]) {
                    if (!sAlive[si]) continue;
                    sAlive[si] = false;
                    // Decrement support of other r-cliques in this s-clique
                    auto &sc = sCliqueList[si];
                    std::vector<int> idx(r);
                    for (int i = 0; i < r; ++i) idx[i] = i;
                    while (true) {
                        std::vector<daf::Size> rc2(r);
                        for (int i = 0; i < r; ++i) rc2[i] = sc[idx[i]];
                        if (rc2 != rc && rCliqueSupport.count(rc2))
                            rCliqueSupport[rc2]--;
                        int i = r - 1;
                        while (i >= 0 && idx[i] == s - r + i) i--;
                        if (i < 0) break;
                        idx[i]++;
                        for (int j = i + 1; j < r; ++j) idx[j] = idx[j-1] + 1;
                    }
                }
                rCliqueSupport.erase(rc);
            }
        }

        // Aggregate to core distribution
        for (auto &[rc, core] : rCliqueCore) {
            coreDist[core]++;
            phase2RCliques++;
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

    // Build result: one entry per r-clique for correct dispatcher output
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    for (auto &[c, cnt] : coreDist)
        for (int i = 0; i < (int)cnt; ++i)
            result.push_back({{}, c});
    return result;
}
