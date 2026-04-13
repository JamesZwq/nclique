//
// Region + ST (V4): Tuple-level peeling with O(r) per-update.
//
// No cliqueIndex, no BK, no C(n,r) enumeration.
// Each support update = ONE C(n,k) per (tuple, leaf).
// Dead tuples tracked per leaf; inclusion-exclusion for corrections.
//

#include "NCliqueCoreDecomposition.h"
#include "graph/DynamicGraphSet.h"
#include <algorithm>
#include <chrono>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
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
    const daf::Size numLeaves = tree.adj_list.size();
    const daf::Size INVALID = std::numeric_limits<daf::Size>::max();

    // ================================================================
    // Step 1: Classes + Tuples + Representatives
    // ================================================================
    std::vector<std::vector<daf::Size>> vtxRegions(numVertices);
    daf::Size numRegions = 0;
    for (auto &mc : g_maxCliques) {
        if ((int)mc.size() < s) continue;
        daf::Size rid = numRegions++;
        for (daf::Size v : mc)
            if (v < numVertices) vtxRegions[v].push_back(rid);
    }

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
            if (vtxRegions[v].empty()) continue;
            auto it = pToC.find(vtxRegions[v]);
            if (it == pToC.end()) {
                daf::Size cid = classSizes.size();
                pToC[vtxRegions[v]] = cid;
                classSizes.push_back(0);
                classVertices.emplace_back();
            }
            daf::Size cid = pToC[vtxRegions[v]];
            classOf[v] = cid;
            classSizes[cid]++;
            classVertices[cid].push_back(v);
        }
    }
    daf::Size numClasses = classSizes.size();

    std::vector<std::vector<daf::Size>> classesInRegion(numRegions);
    for (daf::Size v = 0; v < numVertices; ++v) {
        if (classOf[v] == INVALID) continue;
        for (auto rid : vtxRegions[v]) {
            auto &cr = classesInRegion[rid];
            if (std::find(cr.begin(), cr.end(), classOf[v]) == cr.end())
                cr.push_back(classOf[v]);
        }
    }
    for (auto &v : classesInRegion) std::sort(v.begin(), v.end());

    struct TupleInfo {
        TupleKey key;
        daf::Size mult;
        std::vector<daf::Size> representative;
        // Class composition: classId → count
        std::vector<std::pair<daf::Size, int>> composition;
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
        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            if (classesInRegion[rid].size() > 500) continue;
            cur.clear();
            enumerate(classesInRegion[rid], 0);
        }
    }

    daf::Size totalRCliques = 0;
    for (auto &t : tuples) totalRCliques += t.mult;

    auto tTuples = std::chrono::high_resolution_clock::now();
    printf("======= Region + ST (V4) =======\n");
    printf("  r=%d s=%d, Tuples: %zu, r-cliques: %zu (%.1fx)\n",
           (int)r, (int)s, tuples.size(), (size_t)totalRCliques,
           tuples.empty() ? 0.0 : (double)totalRCliques / tuples.size());

    // ================================================================
    // Step 2: Per-leaf data + initial support
    // ================================================================
    // Per-leaf: hold/pivot per class
    struct LeafClassInfo {
        int nh = 0, np = 0; // hold count, pivot count for this class on this leaf
    };
    // Sparse: only store classes present on each leaf
    struct LeafInfo {
        int h = 0, p = 0;          // total holds, pivots
        std::unordered_map<daf::Size, LeafClassInfo> classes;
        bool alive = true;
    };
    std::vector<LeafInfo> leafInfo(numLeaves);

    for (daf::Size L = 0; L < numLeaves; ++L) {
        auto &li = leafInfo[L];
        for (const auto &node : tree.adj_list[L]) {
            if (node.v >= numVertices || classOf[node.v] == INVALID) continue;
            daf::Size c = classOf[node.v];
            if (node.isPivot) { li.p++; li.classes[c].np++; }
            else { li.h++; li.classes[c].nh++; }
        }
    }

    // Per-leaf: dead tuples (for inclusion-exclusion)
    std::vector<std::vector<daf::Size>> leafDeadTuples(numLeaves);

    // Compute representative's subP on a leaf
    // O(r) per call via treeGraphV set lookup (not O(|leaf|) scan)
    auto repSubP = [&](daf::Size tid, daf::Size L) -> int {
        auto &rep = tuples[tid].representative;
        int subP = 0;
        for (auto v : rep)
            if (treeGraphV.adj_list[v].count(TreeGraphNode{(daf::Size)L, true})) subP++;
        return subP;
    };

    auto repClassPivots = [&](daf::Size tid, daf::Size L) -> std::unordered_map<daf::Size, int> {
        auto &rep = tuples[tid].representative;
        std::unordered_map<daf::Size, int> result;
        for (auto v : rep)
            if (treeGraphV.adj_list[v].count(TreeGraphNode{(daf::Size)L, true}))
                result[classOf[v]]++;
        return result;
    };

    auto repOnLeaf = [&](daf::Size tid, daf::Size L) -> bool {
        auto &rep = tuples[tid].representative;
        for (auto v : rep)
            if (!treeGraphV.adj_list[v].count(TreeGraphNode{(daf::Size)L, true}) &&
                !treeGraphV.adj_list[v].count(TreeGraphNode{(daf::Size)L, false}))
                return false;
        return true;
    };

    // Compute |{S on L : T_rep ⊆ S, S ⊇ ALL tuples in requiredTids}|
    // bRep and subPRep are pre-computed for the representative
    auto computeJoint = [&](const std::unordered_map<daf::Size, int> &bRep, int subPRep,
                            const std::vector<daf::Size> &requiredTids, daf::Size L) -> double {
        auto &li = leafInfo[L];

        int totalMandatory = subPRep;

        // For each class: compute max pivot requirement across all required tuples
        std::unordered_map<daf::Size, int> maxReq;
        for (auto &[c, b] : bRep) maxReq[c] = b;

        for (auto rtid : requiredTids) {
            for (auto &[c, jc] : tuples[rtid].composition) {
                auto cit = li.classes.find(c);
                if (cit == li.classes.end()) return 0.0;
                int dc = std::max(0, jc - cit->second.nh);
                if (maxReq.find(c) == maxReq.end()) maxReq[c] = dc;
                else maxReq[c] = std::max(maxReq[c], dc);
            }
        }

        for (auto &[c, req] : maxReq) {
            int brc = bRep.count(c) ? bRep.at(c) : 0;
            int extra = std::max(0, req - brc);
            auto cit = li.classes.find(c);
            int availExtra = (cit != li.classes.end()) ? cit->second.np - brc : 0;
            if (extra > availExtra) return 0.0;
            totalMandatory += extra;
        }

        int F = li.p - totalMandatory;
        int N = (s - li.h) - totalMandatory;
        if (F < 0 || N < 0 || F < N) return 0.0;
        return nCr[F][N];
    };

    // Alive delta with inclusion-exclusion over dead tuples on L
    // bRep and subPRep are pre-computed for the representative
    auto computeAliveDelta = [&](const std::unordered_map<daf::Size, int> &bRep, int subPRep,
                                 daf::Size peeledTid, daf::Size L) -> double {
        auto &deadList = leafDeadTuples[L];

        // Filter: only dead tuples feasible on this leaf
        std::vector<daf::Size> relevant;
        for (auto dtid : deadList) {
            bool feasible = true;
            for (auto &[c, k] : tuples[dtid].composition) {
                auto it = leafInfo[L].classes.find(c);
                if (it == leafInfo[L].classes.end() || it->second.nh + it->second.np < k)
                    { feasible = false; break; }
            }
            if (feasible) relevant.push_back(dtid);
        }

        if (relevant.empty()) {
            return computeJoint(bRep, subPRep, {peeledTid}, L);
        }

        int p = std::min((int)relevant.size(), 20);
        double result = 0;
        for (int mask = 0; mask < (1 << p); ++mask) {
            std::vector<daf::Size> required = {peeledTid};
            for (int i = 0; i < p; ++i)
                if (mask & (1 << i)) required.push_back(relevant[i]);
            int sign = (__builtin_popcount(mask) % 2 == 0) ? 1 : -1;
            result += sign * computeJoint(bRep, subPRep, required, L);
        }
        return std::max(0.0, result);
    };

    // Initial support via representative
    std::vector<double> support(tuples.size(), 0.0);
    std::vector<std::vector<std::pair<daf::Size, double>>> repLeafContrib(tuples.size());
    std::vector<std::vector<daf::Size>> leafRepTuples(numLeaves);

    for (daf::Size tid = 0; tid < tuples.size(); ++tid) {
        auto &rep = tuples[tid].representative;
        int bestV = 0; size_t bestSz = treeGraphV.adj_list[rep[0]].size();
        for (int vi = 1; vi < (int)rep.size(); ++vi) {
            size_t sz = treeGraphV.adj_list[rep[vi]].size();
            if (sz < bestSz) { bestSz = sz; bestV = vi; }
        }
        for (const auto &node0 : treeGraphV.adj_list[rep[bestV]]) {
            daf::Size L = node0.v;
            if (!leafInfo[L].alive || (int)tree.adj_list[L].size() < s) continue;
            if (!repOnLeaf(tid, L)) continue;
            int subP = repSubP(tid, L);
            int rp = leafInfo[L].p - subP, rn = (s - leafInfo[L].h) - subP;
            if (rp >= 0 && rn >= 0 && rp >= rn) {
                double c = nCr[rp][rn];
                support[tid] += c;
                repLeafContrib[tid].push_back({L, c});
                leafRepTuples[L].push_back(tid);
            }
        }
    }

    auto tSupport = std::chrono::high_resolution_clock::now();
    double totalSupport = 0;
    for (daf::Size i = 0; i < tuples.size(); ++i) totalSupport += tuples[i].mult * support[i];
    printf("  Support sum: %.0f, init: %lld ms\n", totalSupport,
        std::chrono::duration_cast<std::chrono::milliseconds>(tSupport - tTuples).count());

    // ================================================================
    // Step 3: Peeling — no BK, no cliqueIndex
    // ================================================================
    std::vector<double> coreTuple(tuples.size(), 0.0);
    std::vector<bool> tuplePeeled(tuples.size(), false);

    // Bucket queue (per-tuple)
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
    robin_hood::unordered_flat_set<daf::Size> affectedTupleSet;

    auto tPeel = std::chrono::high_resolution_clock::now();

    while (remainingTuples > 0) {
        batchTids.clear();
        affectedTupleSet.clear();

        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) {
            if (overflowSet.empty()) break;
            auto it = overflowSet.begin();
            daf::Size tid = it->second; overflowSet.erase(it);
            if (tuplePeeled[tid]) continue;
            minCore = std::max(support[tid], minCore);
            tuplePeeled[tid] = true; coreTuple[tid] = minCore; remainingTuples--;
            batchTids.push_back(tid);
            goto phase1;
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

        phase1:
        if (batchTids.empty() || remainingTuples == 0) { if (remainingTuples == 0) break; continue; }
        totalIters++;

        // ============================================================
        // Phase 1+2: For each peeled tuple, find affected leaves,
        // compute delta = C(F, N) per (surviving_tuple, leaf)
        // ============================================================

        // Collect all (peeledTid, leaf) pairs first, then process
        // This ensures dead list is consistent within a batch
        std::vector<std::pair<daf::Size, daf::Size>> peelLeafPairs; // (peeledTid, leafId)

        for (auto peeledTid : batchTids) {
            auto &pComp = tuples[peeledTid].composition;

            daf::Size rareC = pComp[0].first;
            for (auto &[c, k] : pComp)
                if (classSizes[c] < classSizes[rareC]) rareC = c;

            robin_hood::unordered_flat_set<daf::Size> candidates;
            for (auto v : classVertices[rareC])
                for (const auto &nd : treeGraphV.adj_list[v])
                    candidates.insert(nd.v);

            for (auto L : candidates) {
                if (!leafInfo[L].alive) continue;
                bool feasible = true;
                for (auto &[c, k] : pComp) {
                    auto it = leafInfo[L].classes.find(c);
                    if (it == leafInfo[L].classes.end() || it->second.nh + it->second.np < k)
                        { feasible = false; break; }
                }
                if (feasible) peelLeafPairs.push_back({peeledTid, L});
            }
        }

        // Compute deltas (dead list does NOT yet include batch tuples)
        constexpr int IE_THRESHOLD = 8; // max dead tuples for inclusion-exclusion (2^8 = 256)

        for (auto &[peeledTid, L] : peelLeafPairs) {
            if (!leafInfo[L].alive) continue;
            int existingDead = (int)leafDeadTuples[L].size();

            if (existingDead < IE_THRESHOLD) {
                // Fast path: inclusion-exclusion (exact, O(2^p × r))
                for (auto tid : leafRepTuples[L]) {
                    if (tuplePeeled[tid]) continue;
                    auto bRep = repClassPivots(tid, L);
                    int subPRep = 0;
                    for (auto &[c, b] : bRep) subPRep += b;
                    double delta = computeAliveDelta(bRep, subPRep, peeledTid, L);
                    if (delta > 0) {
                        support[tid] -= delta;
                        if (support[tid] < 0) support[tid] = 0;
                        auto &contribs = repLeafContrib[tid];
                        for (auto &pr : contribs)
                            if (pr.first == L) { pr.second -= delta; break; }
                        affectedTupleSet.insert(tid);
                    }
                }
            } else {
                // Fallback: kill leaf entirely (subtract remaining contribution)
                // This is exact: remaining contrib was tracked in repLeafContrib
                for (auto tid : leafRepTuples[L]) {
                    if (tuplePeeled[tid]) continue;
                    auto &contribs = repLeafContrib[tid];
                    for (int i = 0; i < (int)contribs.size(); ++i) {
                        if (contribs[i].first == L) {
                            support[tid] -= contribs[i].second;
                            if (support[tid] < 0) support[tid] = 0;
                            contribs[i] = contribs.back();
                            contribs.pop_back();
                            affectedTupleSet.insert(tid);
                            break;
                        }
                    }
                }
                leafRepTuples[L].clear();
                leafInfo[L].alive = false;
                for (const auto &i : tree.adj_list[L])
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(L));
                tree.adj_list[L].clear();
                tree.recycleNode(L);
            }
        }

        // Add batch tuples to dead lists (only for alive leaves)
        for (auto &[peeledTid, L] : peelLeafPairs) {
            if (leafInfo[L].alive)
                leafDeadTuples[L].push_back(peeledTid);
        }

        // Deferred bucketMove
        for (auto tid : affectedTupleSet) {
            if (!tuplePeeled[tid]) bucketMove(tid);
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
