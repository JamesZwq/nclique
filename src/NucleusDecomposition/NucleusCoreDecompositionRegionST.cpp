//
// Region + ST (V4): Maximal-Clique-Level Peeling
//
// No CPI, no cliqueIndex, no BK, no convolution.
// Support computed via inclusion-exclusion over maximal cliques.
// Each IE term = ONE C(n,k). Typically 2-8 terms per tuple.
//

#include "NCliqueCoreDecomposition.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>
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
    // Step 1: MCs, Classes, Tuples
    // ================================================================
    // Filter maximal cliques ≥ s
    std::vector<std::vector<daf::Size>> mcs; // mc[i] = sorted vertex list
    for (auto &mc : g_maxCliques)
        if ((int)mc.size() >= s) mcs.push_back(mc);
    daf::Size numMC = mcs.size();

    // Vertex → MC membership
    std::vector<std::vector<daf::Size>> vtxMCs(numVertices);
    for (daf::Size mi = 0; mi < numMC; ++mi)
        for (auto v : mcs[mi])
            if (v < numVertices) vtxMCs[v].push_back(mi);

    // Overlap classes
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
            auto it = pToC.find(vtxMCs[v]);
            if (it == pToC.end()) {
                daf::Size cid = classSizes.size();
                pToC[vtxMCs[v]] = cid;
                classSizes.push_back(0);
                classVertices.emplace_back();
            }
            daf::Size cid = pToC[vtxMCs[v]];
            classOf[v] = cid;
            classSizes[cid]++;
            classVertices[cid].push_back(v);
        }
    }

    // Classes per MC
    std::vector<std::vector<daf::Size>> mcClasses(numMC);
    for (daf::Size mi = 0; mi < numMC; ++mi) {
        robin_hood::unordered_flat_set<daf::Size> seen;
        for (auto v : mcs[mi])
            if (classOf[v] != INVALID && seen.insert(classOf[v]).second)
                mcClasses[mi].push_back(classOf[v]);
        std::sort(mcClasses[mi].begin(), mcClasses[mi].end());
    }

    // Tuples + representative
    struct TupleInfo {
        TupleKey key;
        daf::Size mult;
        std::vector<daf::Size> representative;
        std::vector<std::pair<daf::Size, int>> composition; // (classId, count)
    };
    std::vector<TupleInfo> tuples;
    std::unordered_map<TupleKey, daf::Size, TupleHash> tupleIndex;
    {
        TupleKey cur; cur.reserve(r);
        std::function<void(const std::vector<daf::Size>&, int)> enumerate;
        enumerate = [&](const std::vector<daf::Size> &classes, int start) {
            if ((int)cur.size() == r) {
                if (tupleIndex.count(cur)) return;
                std::unordered_map<daf::Size, int> counts;
                for (auto c : cur) counts[c]++;
                daf::Size mult = 1;
                std::vector<daf::Size> rep;
                std::vector<std::pair<daf::Size, int>> comp;
                for (auto &[c, k] : counts) {
                    if ((int)classSizes[c] < k) return;
                    mult *= (daf::Size)(nCr[classSizes[c]][k] + 0.5);
                    for (int i = 0; i < k; ++i) rep.push_back(classVertices[c][i]);
                    comp.push_back({c, k});
                }
                if (mult == 0) return;
                std::sort(rep.begin(), rep.end());
                tupleIndex[cur] = tuples.size();
                tuples.push_back({cur, mult, std::move(rep), std::move(comp)});
                return;
            }
            for (int i = start; i < (int)classes.size(); ++i) {
                cur.push_back(classes[i]);
                enumerate(classes, i);
                cur.pop_back();
            }
        };
        for (daf::Size mi = 0; mi < numMC; ++mi) {
            if (mcClasses[mi].size() > 500) continue;
            cur.clear();
            enumerate(mcClasses[mi], 0);
        }
    }

    daf::Size totalRCliques = 0;
    for (auto &t : tuples) totalRCliques += t.mult;

    // ================================================================
    // Step 2: Tuple → MCs mapping
    // ================================================================
    // For each tuple: which MCs contain it?
    // τ feasible in MC m iff for each class c in τ: |c ∩ m| ≥ j_c
    std::vector<std::vector<daf::Size>> tupleMCs(tuples.size());
    std::vector<std::vector<daf::Size>> mcTuples(numMC);

    // Per MC: class sizes within MC
    std::vector<std::unordered_map<daf::Size, int>> mcClassSize(numMC);
    for (daf::Size mi = 0; mi < numMC; ++mi)
        for (auto v : mcs[mi])
            if (classOf[v] != INVALID) mcClassSize[mi][classOf[v]]++;

    for (daf::Size tid = 0; tid < tuples.size(); ++tid) {
        for (daf::Size mi = 0; mi < numMC; ++mi) {
            bool feasible = true;
            for (auto &[c, k] : tuples[tid].composition) {
                auto it = mcClassSize[mi].find(c);
                if (it == mcClassSize[mi].end() || it->second < k) { feasible = false; break; }
            }
            if (feasible) {
                tupleMCs[tid].push_back(mi);
                mcTuples[mi].push_back(tid);
            }
        }
    }

    // ================================================================
    // Step 3: Precompute MC pairwise intersection sizes
    // ================================================================
    // mcInter[i][j] = |M_i ∩ M_j| (only for overlapping pairs)
    std::vector<std::unordered_map<daf::Size, int>> mcInter(numMC);
    for (daf::Size v = 0; v < numVertices; ++v) {
        auto &ms = vtxMCs[v];
        for (int i = 0; i < (int)ms.size(); ++i)
            for (int j = i+1; j < (int)ms.size(); ++j) {
                mcInter[ms[i]][ms[j]]++;
                mcInter[ms[j]][ms[i]]++;
            }
    }

    auto tInit = std::chrono::high_resolution_clock::now();
    printf("======= Region + ST (V4 MC-level) =======\n");
    printf("  r=%d s=%d, MCs: %zu, Tuples: %zu, r-cliques: %zu (%.1fx)\n",
           (int)r, (int)s, (size_t)numMC, tuples.size(), (size_t)totalRCliques,
           tuples.empty() ? 0.0 : (double)totalRCliques / tuples.size());
    printf("  Init time: %lld ms\n",
        std::chrono::duration_cast<std::chrono::milliseconds>(tInit - tStart).count());

    // ================================================================
    // Step 4: Compute initial support via IE over MCs
    // ================================================================
    // support(τ) = Σ_{∅≠A ⊆ MCs(τ)} (-1)^{|A|+1} C(|∩A|-r, s-r)
    // |∩A| = intersection size of all MCs in A

    std::vector<bool> mcAlive(numMC, true);

    // Compute intersection size of a set of MCs
    auto mcIntersectionSize = [&](const std::vector<daf::Size> &mcSet) -> int {
        if (mcSet.empty()) return 0;
        if (mcSet.size() == 1) return (int)mcs[mcSet[0]].size();
        // Intersect sorted vertex lists
        std::vector<daf::Size> cur = mcs[mcSet[0]];
        for (int i = 1; i < (int)mcSet.size(); ++i) {
            std::vector<daf::Size> next;
            std::set_intersection(cur.begin(), cur.end(),
                                  mcs[mcSet[i]].begin(), mcs[mcSet[i]].end(),
                                  std::back_inserter(next));
            cur = std::move(next);
        }
        return (int)cur.size();
    };

    // Compute support of tuple using IE over its alive MCs
    auto computeSupport = [&](daf::Size tid) -> double {
        std::vector<daf::Size> aliveMCs;
        for (auto mi : tupleMCs[tid])
            if (mcAlive[mi]) aliveMCs.push_back(mi);
        if (aliveMCs.empty()) return 0.0;

        int p = (int)aliveMCs.size();
        double result = 0;
        for (int mask = 1; mask < (1 << p); ++mask) {
            std::vector<daf::Size> subset;
            for (int i = 0; i < p; ++i)
                if (mask & (1 << i)) subset.push_back(aliveMCs[i]);
            int interSize = mcIntersectionSize(subset);
            int n = interSize - (int)r, k = (int)s - (int)r;
            double val = (n >= k && k >= 0) ? nCr[n][k] : 0.0;
            int sign = (__builtin_popcount(mask) % 2 == 1) ? 1 : -1;
            result += sign * val;
        }
        return std::max(0.0, result);
    };

    std::vector<double> support(tuples.size());
    for (daf::Size tid = 0; tid < tuples.size(); ++tid)
        support[tid] = computeSupport(tid);

    auto tSupport = std::chrono::high_resolution_clock::now();
    double totalSupport = 0;
    for (daf::Size i = 0; i < tuples.size(); ++i) totalSupport += tuples[i].mult * support[i];
    printf("  Support sum: %.0f, compute: %lld ms\n", totalSupport,
        std::chrono::duration_cast<std::chrono::milliseconds>(tSupport - tInit).count());

    // ================================================================
    // Step 5: Peeling
    // ================================================================
    std::vector<double> coreTuple(tuples.size(), 0.0);
    std::vector<bool> tuplePeeled(tuples.size(), false);

    // Bucket queue on tuples
    constexpr int BUCKET_MAX = 5000000;
    double rawMax = 0;
    for (auto &sv : support) rawMax = std::max(rawMax, sv);
    int maxBucket = (int)std::min(rawMax, (double)BUCKET_MAX);

    int curBucket = 0;
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::set<std::pair<double, daf::Size>> overflowSet;
    std::vector<int> bucket_of(tuples.size(), -1);
    std::vector<daf::Size> pos_in_bucket(tuples.size());
    std::vector<double> overflowStoredVal(tuples.size(), -1);

    auto bucketInsert = [&](daf::Size tid) {
        double val = std::max(0.0, support[tid]);
        int b = std::min((int)val, maxBucket + 1);
        if (b <= maxBucket) { pos_in_bucket[tid] = buckets[b].size(); buckets[b].push_back(tid); }
        else { overflowSet.insert({val, tid}); overflowStoredVal[tid] = val; }
        bucket_of[tid] = b;
    };
    auto bucketMove = [&](daf::Size tid) {
        if (tuplePeeled[tid]) return;
        double val = std::max(0.0, support[tid]);
        int oldB = bucket_of[tid];
        int newB = std::min((int)val, maxBucket + 1);
        if (oldB == newB && (newB <= maxBucket || val == overflowStoredVal[tid])) return;
        if (oldB >= 0 && oldB <= maxBucket) {
            auto &bk = buckets[oldB]; daf::Size p = pos_in_bucket[tid];
            if (p < bk.size()-1) { daf::Size l = bk.back(); bk[p] = l; pos_in_bucket[l] = p; }
            bk.pop_back();
        } else if (oldB == maxBucket+1) { overflowSet.erase({overflowStoredVal[tid], tid}); }
        if (newB <= maxBucket) {
            pos_in_bucket[tid] = buckets[newB].size(); buckets[newB].push_back(tid);
            if (newB < curBucket) curBucket = newB;
        } else { overflowSet.insert({val, tid}); overflowStoredVal[tid] = val; }
        bucket_of[tid] = newB;
    };
    for (daf::Size tid = 0; tid < tuples.size(); ++tid) bucketInsert(tid);

    daf::Size remainingTuples = tuples.size();
    double minCore = 0;
    long long totalIters = 0;

    std::vector<daf::Size> batchTids;
    robin_hood::unordered_flat_set<daf::Size> affectedTuples;
    robin_hood::unordered_flat_set<daf::Size> deadMCs; // MCs killed in this iteration

    auto tPeel = std::chrono::high_resolution_clock::now();

    while (remainingTuples > 0) {
        batchTids.clear();
        affectedTuples.clear();
        deadMCs.clear();

        // Pop min-support tuples
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) {
            if (overflowSet.empty()) break;
            auto it = overflowSet.begin();
            daf::Size tid = it->second; overflowSet.erase(it);
            if (tuplePeeled[tid]) continue;
            minCore = std::max(support[tid], minCore);
            tuplePeeled[tid] = true; coreTuple[tid] = minCore; remainingTuples--;
            batchTids.push_back(tid);
            goto process;
        }
        minCore = std::max((double)curBucket, minCore);
        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty() && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                daf::Size tid = buckets[curBucket].back(); buckets[curBucket].pop_back();
                if (tuplePeeled[tid]) continue;
                tuplePeeled[tid] = true; coreTuple[tid] = minCore; remainingTuples--;
                batchTids.push_back(tid);
            }
            if (curBucket+1 < (int)buckets.size() && !buckets[curBucket+1].empty() && (curBucket+1) <= (int)minCore) curBucket++;
            else break;
        }

        process:
        if (batchTids.empty() || remainingTuples == 0) { if (remainingTuples == 0) break; continue; }
        totalIters++;

        // For each peeled tuple: check if any MC should die
        for (auto tid : batchTids) {
            for (auto mi : tupleMCs[tid]) {
                if (!mcAlive[mi]) continue;
                // MC dies if ALL its tuples are peeled or have support ≤ minCore
                // Simplified: MC dies when C(|M|-r, s-r) ≤ minCore
                // (its unique tuples have this static support)
                if (nCr[(int)mcs[mi].size() - (int)r][(int)s - (int)r] <= minCore) {
                    mcAlive[mi] = false;
                    deadMCs.insert(mi);
                }
            }
        }

        // Update support for tuples affected by dead MCs
        if (!deadMCs.empty()) {
            for (auto mi : deadMCs) {
                for (auto tid : mcTuples[mi]) {
                    if (!tuplePeeled[tid]) affectedTuples.insert(tid);
                }
            }
            for (auto tid : affectedTuples) {
                double newSup = computeSupport(tid);
                support[tid] = newSup;
                bucketMove(tid);
            }
        }
    }

    auto tEnd = std::chrono::high_resolution_clock::now();
    auto peelMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tPeel).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();
    printf("  Peeling: %lld iter\n", totalIters);
    printf("  Peeling time: %lld ms\n", peelMs);
    printf("  Total time: %lld ms\n", totalMs);

    // ================================================================
    // Result
    // ================================================================
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    std::map<double, daf::Size> coreDist;
    for (daf::Size tid = 0; tid < tuples.size(); ++tid) {
        result.push_back({tuples[tid].representative, coreTuple[tid]});
        coreDist[coreTuple[tid]] += tuples[tid].mult;
    }
    double maxCore = 0;
    for (auto &[c, cnt] : coreDist) maxCore = std::max(maxCore, c);
    printf("  Max core: %.0f\n", maxCore);
    for (auto &[c, cnt] : coreDist)
        printf("  core=%.0f count=%zu\n", c, (size_t)cnt);

    return result;
}
