//
// V4: Two-Phase (r,s)-Nucleus Decomposition
//
// Phase 1: Iterative MC deletion — MCs with private vertices.
//   Count: C(|M|,r) - C(|M|-|priv|,r) r-cliques at level C(|M|-r,s-r).
//   No tuples, no classes needed.
//
// Phase 2: Remaining MCs (no private vertices).
//   Build classes+tuples, compute IE support, peel.
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

    // ================================================================
    // Build MCs
    // ================================================================
    std::vector<std::vector<daf::Size>> mcs;
    for (auto &mc : g_maxCliques)
        if ((int)mc.size() >= s) mcs.push_back(mc);
    daf::Size numMC = mcs.size();

    // Vertex → alive MC list
    std::vector<std::vector<daf::Size>> vtxAliveMCs(numVertices);
    for (daf::Size mi = 0; mi < numMC; ++mi)
        for (auto v : mcs[mi])
            if (v < numVertices) vtxAliveMCs[v].push_back(mi);

    std::vector<bool> mcAlive(numMC, true);

    // ================================================================
    // Phase 1: Iterative MC deletion
    // ================================================================
    // MC has "private vertex" iff some vertex v has vtxAliveMCs[v] == {mi}
    // When M is deleted: r-cliques with private vertices get core = C(|M|-r, s-r)
    // Count = C(|M|, r) - C(|M| - #private, r)

    auto mcCoreValue = [&](daf::Size mi) -> double {
        int n = (int)mcs[mi].size() - (int)r, k = (int)s - (int)r;
        return (n >= k && k >= 0) ? nCr[n][k] : 0.0;
    };

    auto getPrivateCount = [&](daf::Size mi) -> int {
        int priv = 0;
        for (auto v : mcs[mi])
            if (v < numVertices && vtxAliveMCs[v].size() == 1 && vtxAliveMCs[v][0] == mi)
                priv++;
        return priv;
    };

    // Priority queue: (C-value, MC index), smallest first
    using PQE = std::pair<double, daf::Size>;
    std::priority_queue<PQE, std::vector<PQE>, std::greater<>> pq;
    for (daf::Size mi = 0; mi < numMC; ++mi)
        if (getPrivateCount(mi) > 0) pq.push({mcCoreValue(mi), mi});

    std::map<double, double> phase1Dist; // core -> r-clique count (double for large nCr)
    double minCore = 0;
    daf::Size phase1MCs = 0;

    while (!pq.empty()) {
        auto [val, mi] = pq.top(); pq.pop();
        if (!mcAlive[mi]) continue;
        int priv = getPrivateCount(mi);
        if (priv == 0) continue;

        minCore = std::max(minCore, val);

        // Count r-cliques with at least one private vertex
        int M = (int)mcs[mi].size();
        double total = nCr[M][r];
        double noPriv = nCr[M - priv][r]; // r-cliques using NO private vertex
        double withPriv = total - noPriv;
        if (withPriv > 0) phase1Dist[minCore] += withPriv;

        // Delete MC
        mcAlive[mi] = false;
        phase1MCs++;

        // Update profiles, find new private vertices
        std::set<daf::Size> check;
        for (auto v : mcs[mi]) {
            if (v >= numVertices) continue;
            auto &am = vtxAliveMCs[v];
            am.erase(std::remove(am.begin(), am.end(), mi), am.end());
            if (am.size() == 1 && mcAlive[am[0]]) check.insert(am[0]);
        }
        for (auto m2 : check)
            if (getPrivateCount(m2) > 0) pq.push({mcCoreValue(m2), m2});
    }

    auto tP1 = std::chrono::high_resolution_clock::now();

    // ================================================================
    // Phase 2: Remaining MCs — build classes, tuples, IE support
    // ================================================================
    std::vector<daf::Size> remainMCs;
    for (daf::Size mi = 0; mi < numMC; ++mi)
        if (mcAlive[mi]) remainMCs.push_back(mi);

    std::map<double, double> phase2Dist;
    daf::Size phase2Tuples = 0;

    if (!remainMCs.empty()) {
        // Build classes from remaining alive profiles
        const daf::Size INVALID = std::numeric_limits<daf::Size>::max();
        using Profile = std::vector<daf::Size>;
        struct PH {
            size_t operator()(const Profile &p) const noexcept {
                size_t h = p.size();
                for (auto x : p) h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h<<6) + (h>>2);
                return h;
            }
        };
        std::unordered_map<Profile, daf::Size, PH> pToC;
        std::vector<daf::Size> classOf(numVertices, INVALID);
        std::vector<daf::Size> classSizes;
        std::vector<std::vector<daf::Size>> classVerts;

        for (daf::Size v = 0; v < numVertices; ++v) {
            if (vtxAliveMCs[v].empty()) continue;
            auto it = pToC.find(vtxAliveMCs[v]);
            if (it == pToC.end()) {
                daf::Size cid = classSizes.size();
                pToC[vtxAliveMCs[v]] = cid;
                classSizes.push_back(0);
                classVerts.emplace_back();
            }
            daf::Size cid = pToC[vtxAliveMCs[v]];
            classOf[v] = cid;
            classSizes[cid]++;
            classVerts[cid].push_back(v);
        }

        // MC class info
        std::vector<std::unordered_map<daf::Size, int>> mcCS(numMC);
        std::vector<std::vector<daf::Size>> mcCL(numMC);
        for (auto mi : remainMCs) {
            for (auto v : mcs[mi])
                if (v < numVertices && classOf[v] != INVALID) mcCS[mi][classOf[v]]++;
            for (auto &[c, _] : mcCS[mi]) mcCL[mi].push_back(c);
            std::sort(mcCL[mi].begin(), mcCL[mi].end());
        }

        // Enumerate tuples
        struct TI { TupleKey key; daf::Size mult; };
        std::vector<TI> tuples;
        std::unordered_map<TupleKey, daf::Size, TupleHash> tidx;
        {
            TupleKey cur; cur.reserve(r);
            std::function<void(const std::vector<daf::Size>&, int)> en;
            en = [&](const std::vector<daf::Size> &cl, int st) {
                if ((int)cur.size() == r) {
                    if (tidx.count(cur)) return;
                    std::unordered_map<daf::Size, int> cnt;
                    for (auto c : cur) cnt[c]++;
                    daf::Size mult = 1;
                    for (auto &[c, k] : cnt) {
                        if ((int)classSizes[c] < k) return;
                        mult *= (daf::Size)(nCr[classSizes[c]][k] + 0.5);
                    }
                    if (!mult) return;
                    tidx[cur] = tuples.size();
                    tuples.push_back({cur, mult});
                    return;
                }
                for (int i = st; i < (int)cl.size(); ++i) {
                    cur.push_back(cl[i]);
                    en(cl, i);
                    cur.pop_back();
                }
            };
            for (auto mi : remainMCs) {
                cur.clear();
                en(mcCL[mi], 0);
            }
        }

        // Tuple → MCs
        std::vector<std::vector<daf::Size>> tMCs(tuples.size());
        for (daf::Size ti = 0; ti < tuples.size(); ++ti) {
            std::unordered_map<daf::Size, int> cnt;
            for (auto c : tuples[ti].key) cnt[c]++;
            for (auto mi : remainMCs) {
                bool ok = true;
                for (auto &[c, k] : cnt) {
                    auto it = mcCS[mi].find(c);
                    if (it == mcCS[mi].end() || it->second < k) { ok = false; break; }
                }
                if (ok) tMCs[ti].push_back(mi);
            }
        }

        // MC intersection
        auto interSz = [&](const std::vector<daf::Size> &ms) -> int {
            if (ms.empty()) return 0;
            if (ms.size() == 1) return (int)mcs[ms[0]].size();
            auto cur = mcs[ms[0]];
            for (size_t i = 1; i < ms.size(); ++i) {
                std::vector<daf::Size> nxt;
                std::set_intersection(cur.begin(), cur.end(), mcs[ms[i]].begin(), mcs[ms[i]].end(), std::back_inserter(nxt));
                cur = std::move(nxt);
            }
            return (int)cur.size();
        };

        // IE support
        for (daf::Size ti = 0; ti < tuples.size(); ++ti) {
            auto &tm = tMCs[ti];
            int p = (int)tm.size();
            double sup = 0;
            for (int mask = 1; mask < (1 << p); ++mask) {
                std::vector<daf::Size> sub;
                for (int i = 0; i < p; ++i) if (mask & (1 << i)) sub.push_back(tm[i]);
                int isz = interSz(sub);
                int n = isz - (int)r, k = (int)s - (int)r;
                double v = (n >= k && k >= 0) ? nCr[n][k] : 0.0;
                sup += (__builtin_popcount(mask) % 2 == 1 ? 1 : -1) * v;
            }
            sup = std::max(0.0, sup);
            double core = sup;  // Phase 2 tuples get their own support as core
            phase2Dist[core] += tuples[ti].mult;
            phase2Tuples++;
        }
    }

    auto tEnd = std::chrono::high_resolution_clock::now();

    // Merge
    std::map<double, double> coreDist;
    for (auto &[c, n] : phase1Dist) coreDist[c] += n;
    for (auto &[c, n] : phase2Dist) coreDist[c] += n;

    printf("======= V4 Two-Phase =======\n");
    printf("  r=%d s=%d, MCs: %zu\n", (int)r, (int)s, (size_t)numMC);
    printf("  Phase 1: %zu MCs peeled, %lld ms\n", (size_t)phase1MCs,
        std::chrono::duration_cast<std::chrono::milliseconds>(tP1 - tStart).count());
    printf("  Phase 2: %zu MCs, %zu tuples, %lld ms\n", remainMCs.size(), (size_t)phase2Tuples,
        std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tP1).count());
    printf("  Total time: %lld ms\n",
        std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count());

    double maxCore = 0;
    for (auto &[c, _] : coreDist) maxCore = std::max(maxCore, c);
    printf("  Max core: %.0f\n", maxCore);
    for (auto &[c, cnt] : coreDist)
        printf("  core=%.0f count=%.0f\n", c, cnt);

    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    for (auto &[c, cnt] : coreDist) result.push_back({{}, c});
    return result;
}
