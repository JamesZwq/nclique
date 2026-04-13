//
// V4: Direct Core Assignment — no peeling needed.
//
// Theorem: core(T) = max_{M ⊇ T} C(|M|-r, s-r)
// where M ranges over maximal cliques of size ≥ s.
//
// Algorithm: O(MaxCliqEnum + Σ|M|)
//

#include "NCliqueCoreDecomposition.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <functional>

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
    // Step 1: MCs + Classes
    // ================================================================
    std::vector<std::pair<int, daf::Size>> mcSizeIdx; // (size, index)
    std::vector<std::vector<daf::Size>> mcs;
    for (auto &mc : g_maxCliques) {
        if ((int)mc.size() >= s) {
            mcSizeIdx.push_back({(int)mc.size(), mcs.size()});
            mcs.push_back(mc);
        }
    }
    daf::Size numMC = mcs.size();

    // Vertex → MC membership
    std::vector<std::vector<daf::Size>> vtxMCs(numVertices);
    for (daf::Size mi = 0; mi < numMC; ++mi)
        for (auto v : mcs[mi])
            if (v < numVertices) vtxMCs[v].push_back(mi);

    // Classes
    std::vector<daf::Size> classOf(numVertices, INVALID);
    std::vector<daf::Size> classSizes;
    std::vector<std::vector<daf::Size>> classVertices;
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
            if (vtxMCs[v].empty()) continue;
            auto &prof = vtxMCs[v];
            auto it = pToC.find(prof);
            if (it == pToC.end()) {
                daf::Size cid = classSizes.size();
                pToC[prof] = cid;
                classSizes.push_back(0);
                classVertices.emplace_back();
            }
            daf::Size cid = pToC[vtxMCs[v]];
            classOf[v] = cid;
            classSizes[cid]++;
            classVertices[cid].push_back(v);
        }
    }

    // Classes per MC + class sizes within MC
    std::vector<std::vector<daf::Size>> mcClasses(numMC);
    std::vector<std::unordered_map<daf::Size, int>> mcClassSize(numMC);
    for (daf::Size mi = 0; mi < numMC; ++mi) {
        for (auto v : mcs[mi]) {
            if (classOf[v] == INVALID) continue;
            mcClassSize[mi][classOf[v]]++;
        }
        for (auto &[c, _] : mcClassSize[mi])
            mcClasses[mi].push_back(c);
        std::sort(mcClasses[mi].begin(), mcClasses[mi].end());
    }

    // ================================================================
    // Step 2: Enumerate tuples + assign core directly
    // ================================================================
    // core(τ) = max over MCs containing τ of C(|M|-r, s-r)
    //         = C(max_MC_size - r, s - r)

    struct TupleInfo {
        TupleKey key;
        daf::Size mult;
        std::vector<daf::Size> representative;
        double core;
    };
    std::vector<TupleInfo> tuples;
    std::unordered_map<TupleKey, daf::Size, TupleHash> tupleIndex;
    {
        TupleKey cur; cur.reserve(r);
        std::function<void(const std::vector<daf::Size>&, int, daf::Size)> enumerate;
        // Pass mc index to track which MC we're enumerating from
        enumerate = [&](const std::vector<daf::Size> &classes, int start, daf::Size mcIdx) {
            if ((int)cur.size() == r) {
                auto it = tupleIndex.find(cur);
                if (it != tupleIndex.end()) {
                    // Already exists: update core if this MC is larger
                    double newCore = nCr[(int)mcs[mcIdx].size() - (int)r][(int)s - (int)r];
                    if (newCore > tuples[it->second].core)
                        tuples[it->second].core = newCore;
                    return;
                }
                std::unordered_map<daf::Size, int> counts;
                for (auto c : cur) counts[c]++;
                daf::Size mult = 1;
                std::vector<daf::Size> rep;
                for (auto &[c, k] : counts) {
                    if ((int)classSizes[c] < k) return;
                    mult *= (daf::Size)(nCr[classSizes[c]][k] + 0.5);
                    for (int i = 0; i < k; ++i) rep.push_back(classVertices[c][i]);
                }
                if (mult == 0) return;
                std::sort(rep.begin(), rep.end());
                double core = nCr[(int)mcs[mcIdx].size() - (int)r][(int)s - (int)r];
                tupleIndex[cur] = tuples.size();
                tuples.push_back({cur, mult, std::move(rep), core});
                return;
            }
            for (int i = start; i < (int)classes.size(); ++i) {
                cur.push_back(classes[i]);
                enumerate(classes, i, mcIdx);
                cur.pop_back();
            }
        };
        // Enumerate from LARGEST MCs first (so first insertion has highest core)
        // But we update anyway, so order doesn't matter
        for (daf::Size mi = 0; mi < numMC; ++mi) {
            if (mcClasses[mi].size() > 500) continue;
            cur.clear();
            enumerate(mcClasses[mi], 0, mi);
        }
    }

    auto tEnd = std::chrono::high_resolution_clock::now();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();

    daf::Size totalRCliques = 0;
    for (auto &t : tuples) totalRCliques += t.mult;

    printf("======= Region + ST (V4 Direct) =======\n");
    printf("  r=%d s=%d, MCs: %zu, Tuples: %zu, r-cliques: %zu (%.1fx)\n",
           (int)r, (int)s, (size_t)numMC, tuples.size(), (size_t)totalRCliques,
           tuples.empty() ? 0.0 : (double)totalRCliques / tuples.size());
    printf("  Total time: %lld ms\n", totalMs);

    // ================================================================
    // Result
    // ================================================================
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    std::map<double, daf::Size> coreDist;
    for (auto &t : tuples) {
        result.push_back({t.representative, t.core});
        coreDist[t.core] += t.mult;
    }
    double maxCore = 0;
    for (auto &[c, cnt] : coreDist) maxCore = std::max(maxCore, c);
    printf("  Max core: %.0f\n", maxCore);
    for (auto &[c, cnt] : coreDist)
        printf("  core=%.0f count=%zu\n", c, (size_t)cnt);

    return result;
}
