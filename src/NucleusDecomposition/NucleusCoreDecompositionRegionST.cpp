//
// Region + ST (V4): Tuple-level peeling WITHOUT per-r-clique enumeration.
//
// Key innovation: Phase 2 uses class-level LeafDeath check + per-tuple
// support update, avoiding the O(C(n,r)) r-clique enumeration.
//
// No cliqueIndex. No per-r-clique support tracking.
// Support tracked per-tuple via representative r-clique.
//

#include "NCliqueCoreDecomposition.h"
#include "../BK/BronKerboschRmRClique.hpp"
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

    // Classes per region
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

    // Tuples + representative
    struct TupleInfo {
        TupleKey key;
        daf::Size mult;
        std::vector<daf::Size> representative; // specific vertex IDs
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
                for (auto &[c, k] : counts) {
                    if ((int)classSizes[c] < k) return;
                    mult *= (daf::Size)(nCr[classSizes[c]][k] + 0.5);
                    for (int i = 0; i < k; ++i) rep.push_back(classVertices[c][i]);
                }
                if (mult == 0) return;
                std::sort(rep.begin(), rep.end());
                tupleIndex[cur] = tuples.size();
                tuples.push_back({cur, mult, std::move(rep)});
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
    // Step 2: Initial support via representative + build mappings
    // ================================================================
    // Per-leaf: pivot/keep counts
    std::vector<int> leafPivotC(numLeaves, 0);
    std::vector<int> leafNeedPivot(numLeaves, 0);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        int pC = 0, kC = 0;
        for (const auto &node : tree.adj_list[L]) {
            if (node.isPivot) pC++; else kC++;
        }
        leafPivotC[L] = pC;
        leafNeedPivot[L] = s - kC;
    }

    std::vector<double> support(tuples.size(), 0.0);
    // repLeafContrib[tid] = {(leafId, contrib), ...}
    std::vector<std::vector<std::pair<daf::Size, double>>> repLeafContrib(tuples.size());
    // leafRepTuples[L] = {tid, ...} — tuples whose representative is on leaf L
    std::vector<std::vector<daf::Size>> leafRepTuples(numLeaves);

    auto repSubP = [&](const std::vector<daf::Size> &rep, daf::Size L) -> int {
        const auto &leaf = tree.adj_list[L];
        int subP = 0;
        for (auto v : rep)
            for (const auto &ln : leaf)
                if (ln.v == v) { if (ln.isPivot) subP++; break; }
        return subP;
    };

    auto repOnLeaf = [&](const std::vector<daf::Size> &rep, daf::Size L) -> bool {
        const auto &leaf = tree.adj_list[L];
        for (auto v : rep) {
            bool found = false;
            for (const auto &ln : leaf) if (ln.v == v) { found = true; break; }
            if (!found) return false;
        }
        return true;
    };

    for (daf::Size tid = 0; tid < tuples.size(); ++tid) {
        auto &rep = tuples[tid].representative;
        // Find leaves containing ALL rep vertices
        int bestV = 0; size_t bestSz = treeGraphV.adj_list[rep[0]].size();
        for (int vi = 1; vi < (int)rep.size(); ++vi) {
            size_t sz = treeGraphV.adj_list[rep[vi]].size();
            if (sz < bestSz) { bestSz = sz; bestV = vi; }
        }
        for (const auto &node0 : treeGraphV.adj_list[rep[bestV]]) {
            daf::Size L = node0.v;
            if (tree.adj_list[L].empty() || (int)tree.adj_list[L].size() < s) continue;
            if (!repOnLeaf(rep, L)) continue;
            int subP = repSubP(rep, L);
            int rp = leafPivotC[L] - subP, rn = leafNeedPivot[L] - subP;
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
    // Step 3: Peeling
    // ================================================================
    std::vector<double> coreTuple(tuples.size(), 0.0);
    std::vector<bool> tuplePeeled(tuples.size(), false);

    // Per-tuple bucket queue
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
        if (b <= maxBucket) {
            pos_in_bucket[tid] = buckets[b].size();
            buckets[b].push_back(tid);
        } else {
            overflowSet.insert({val, tid});
            overflowStoredVal[tid] = val;
        }
        bucket_of[tid] = b;
    };

    auto bucketMove = [&](daf::Size tid) {
        if (tuplePeeled[tid]) return;
        double val = std::max(0.0, support[tid]);
        int oldB = bucket_of[tid];
        int newB = std::min((int)val, maxBucket + 1);
        if (oldB == newB && (newB <= maxBucket || val == overflowStoredVal[tid])) return;
        // Remove from old
        if (oldB >= 0 && oldB <= maxBucket) {
            auto &bk = buckets[oldB];
            daf::Size p = pos_in_bucket[tid];
            if (p < bk.size()-1) { daf::Size l = bk.back(); bk[p] = l; pos_in_bucket[l] = p; }
            bk.pop_back();
        } else if (oldB == maxBucket+1) {
            overflowSet.erase({overflowStoredVal[tid], tid});
        }
        // Insert to new
        if (newB <= maxBucket) {
            pos_in_bucket[tid] = buckets[newB].size();
            buckets[newB].push_back(tid);
            if (newB < curBucket) curBucket = newB;
        } else {
            overflowSet.insert({val, tid});
            overflowStoredVal[tid] = val;
        }
        bucket_of[tid] = newB;
    };

    for (daf::Size tid = 0; tid < tuples.size(); ++tid) bucketInsert(tid);

    daf::Size remainingTuples = tuples.size();
    double minCore = 0;
    long long totalIters = 0, cntLeafDeath = 0, cntBK = 0;

    // Batch of peeled tuples per iteration
    std::vector<daf::Size> batchTids;
    // Affected tuples for deferred bucketMove
    robin_hood::unordered_flat_set<daf::Size> affectedTupleSet;

    auto tPeel = std::chrono::high_resolution_clock::now();

    while (remainingTuples > 0) {
        batchTids.clear();
        affectedTupleSet.clear();

        // Pop tuples at current level
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) {
            if (overflowSet.empty()) break;
            auto it = overflowSet.begin();
            daf::Size tid = it->second;
            overflowSet.erase(it);
            if (tuplePeeled[tid]) continue;
            minCore = std::max(support[tid], minCore);
            tuplePeeled[tid] = true;
            coreTuple[tid] = minCore;
            remainingTuples--;
            batchTids.push_back(tid);
            goto phase1;
        }

        minCore = std::max((double)curBucket, minCore);
        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                daf::Size tid = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                if (tuplePeeled[tid]) continue;
                tuplePeeled[tid] = true;
                coreTuple[tid] = minCore;
                remainingTuples--;
                batchTids.push_back(tid);
            }
            if (curBucket+1 < (int)buckets.size() && !buckets[curBucket+1].empty()
                && (curBucket+1) <= (int)minCore) curBucket++;
            else break;
        }

        phase1:
        if (batchTids.empty() || remainingTuples == 0) {
            if (remainingTuples == 0) break;
            continue;
        }
        totalIters++;

        // ============================================================
        // Phase 1+2: Find affected leaves + class-level update
        // ============================================================
        // For each peeled tuple, find leaves containing its r-cliques.
        // Use class → vertex → treeGraphV.

        robin_hood::unordered_flat_set<daf::Size> processedLeaves;

        for (auto peeledTid : batchTids) {
            auto &key = tuples[peeledTid].key;
            std::unordered_map<daf::Size, int> classNeeds;
            for (auto c : key) classNeeds[c]++;

            // Find rarest class for iteration
            daf::Size rareC = key[0]; int rareN = classNeeds[key[0]];
            for (auto &[c, k] : classNeeds)
                if ((int)classSizes[c] < (int)classSizes[rareC]) { rareC = c; rareN = k; }

            // Collect candidate leaves
            robin_hood::unordered_flat_set<daf::Size> candidates;
            for (auto v : classVertices[rareC])
                for (const auto &nd : treeGraphV.adj_list[v])
                    candidates.insert(nd.v);

            for (auto L : candidates) {
                if (processedLeaves.count(L)) continue;
                const auto &leaf = tree.adj_list[L];
                if (leaf.empty() || (int)leaf.size() < s) continue;

                // Compute class counts on leaf
                std::unordered_map<daf::Size, int> leafCC; // total per class
                std::unordered_map<daf::Size, int> leafHold, leafPiv;
                for (const auto &ln : leaf) {
                    if (ln.v >= numVertices || classOf[ln.v] == INVALID) continue;
                    daf::Size c = classOf[ln.v];
                    leafCC[c]++;
                    if (ln.isPivot) leafPiv[c]++; else leafHold[c]++;
                }

                // Check feasibility: does this leaf contain any r-clique of peeled tuple?
                bool feasible = true;
                for (auto &[c, k] : classNeeds)
                    if (leafCC[c] < k) { feasible = false; break; }
                if (!feasible) continue;

                processedLeaves.insert(L);
                int n = (int)leaf.size();
                int pivotC = leafPivotC[L], keepC = n - pivotC;
                int needPivot = s - keepC;

                // Class-level LeafDeath check
                daf::Size maxRCliquePerVertex = (daf::Size)(nCr[n-1][r-1] + 0.5);

                // Compute per-class conflict degree from ALL peeled tuples in batch
                bool leafDead = false;
                int forcedPivotRemove = 0, forcedHoldRemove = 0;

                // For each class on this leaf: compute total conflict from all peeled tuples
                for (const auto &ln : leaf) {
                    if (ln.v >= numVertices || classOf[ln.v] == INVALID) continue;
                    daf::Size ci = classOf[ln.v];
                    daf::Size totalConflict = 0;
                    for (auto pt : batchTids) {
                        auto &pkey = tuples[pt].key;
                        std::unordered_map<daf::Size, int> pCounts;
                        for (auto c : pkey) pCounts[c]++;
                        if (pCounts.find(ci) == pCounts.end() || pCounts[ci] == 0) continue;
                        // Check feasibility of this peeled tuple on this leaf
                        bool pFeas = true;
                        for (auto &[c, k] : pCounts)
                            if (leafCC[c] < k) { pFeas = false; break; }
                        if (!pFeas) continue;
                        // Conflict from this peeled tuple for vertex in class ci
                        daf::Size conf = (daf::Size)(nCr[leafCC[ci]-1][pCounts[ci]-1] + 0.5);
                        for (auto &[c, k] : pCounts) {
                            if (c == ci) continue;
                            conf *= (daf::Size)(nCr[leafCC[c]][k] + 0.5);
                        }
                        totalConflict += conf;
                    }
                    if (totalConflict >= maxRCliquePerVertex) {
                        if (!ln.isPivot) { leafDead = true; break; }
                        forcedPivotRemove++;
                    }
                }
                if (!leafDead) {
                    int remPivots = pivotC - forcedPivotRemove;
                    int remTotal = n - forcedPivotRemove - forcedHoldRemove;
                    if (remTotal < (int)s || remPivots < needPivot) leafDead = true;
                }

                if (leafDead) {
                    cntLeafDeath++;
                    // Subtract representative contributions for all tuples on this leaf
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
                    // Remove leaf from treeGraphV
                    for (const auto &i : leaf)
                        treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(L));
                    tree.adj_list[L].clear();
                    tree.recycleNode(L);
                } else {
                    cntBK++;
                    // BK case: enumerate peeled r-cliques on this leaf
                    std::vector<std::vector<daf::Size>> removedRC;
                    for (auto pt : batchTids) {
                        auto &pkey = tuples[pt].key;
                        std::unordered_map<daf::Size, int> pCounts;
                        for (auto c : pkey) pCounts[c]++;
                        // Check feasibility
                        bool pFeas = true;
                        for (auto &[c, k] : pCounts)
                            if (leafCC[c] < k) { pFeas = false; break; }
                        if (!pFeas) continue;
                        // Enumerate r-cliques of this tuple on this leaf
                        std::unordered_map<daf::Size, std::vector<daf::Size>> cverts;
                        for (const auto &ln : leaf)
                            if (ln.v < numVertices && classOf[ln.v] != INVALID)
                                cverts[classOf[ln.v]].push_back(ln.v);
                        std::vector<std::pair<daf::Size, int>> classList(pCounts.begin(), pCounts.end());
                        std::vector<daf::Size> current;
                        std::function<void(int)> gen = [&](int ci) {
                            if (ci == (int)classList.size()) {
                                auto sorted = current;
                                std::sort(sorted.begin(), sorted.end());
                                removedRC.push_back(sorted);
                                return;
                            }
                            auto [cls, need] = classList[ci];
                            auto &verts = cverts[cls];
                            std::function<void(int,int)> choose = [&](int start, int rem) {
                                if (rem == 0) { gen(ci+1); return; }
                                for (int i = start; i <= (int)verts.size()-rem; ++i) {
                                    current.push_back(verts[i]);
                                    choose(i+1, rem-1);
                                    current.pop_back();
                                }
                            };
                            choose(0, need);
                        };
                        gen(0);
                    }

                    // Subtract old rep contributions
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

                    // Remove old leaf from treeGraphV
                    for (const auto &lv : leaf) {
                        if (lv.isPivot) treeGraphV.removeNbr(lv.v, {L, true});
                        else treeGraphV.removeNbr(lv.v, {L, false});
                    }

                    // Run BK
                    auto oldRepTuples = leafRepTuples[L];
                    leafRepTuples[L].clear();

                    bkRmClique::removeRClique(leaf, removedRC, r, s,
                        [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                            auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                            auto newId = tree.addNode(newLeaf);
                            const auto &stored = tree.adj_list[newId];
                            int newPC = 0, newKC = 0;
                            for (const auto &i : stored) {
                                if (i.isPivot) { treeGraphV.addNbr(i.v, {newId, true}); newPC++; }
                                else { treeGraphV.addNbr(i.v, {newId, false}); newKC++; }
                            }
                            // Extend arrays
                            if (newId >= leafPivotC.size()) {
                                leafPivotC.resize(newId+1, 0);
                                leafNeedPivot.resize(newId+1, 0);
                                leafRepTuples.resize(newId+1);
                            }
                            leafPivotC[newId] = newPC;
                            leafNeedPivot[newId] = s - newKC;

                            // Add rep contributions for tuples whose rep is on new sub-leaf
                            for (auto tid : oldRepTuples) {
                                if (tuplePeeled[tid]) continue;
                                if (repOnLeaf(tuples[tid].representative, newId)) {
                                    int subP = repSubP(tuples[tid].representative, newId);
                                    int rp = newPC - subP, rn = (s - newKC) - subP;
                                    if (rp >= 0 && rn >= 0 && rp >= rn) {
                                        double contrib = nCr[rp][rn];
                                        support[tid] += contrib;
                                        repLeafContrib[tid].push_back({newId, contrib});
                                        leafRepTuples[newId].push_back(tid);
                                        affectedTupleSet.insert(tid);
                                    }
                                }
                            }
                        });

                    tree.removeNode(L);
                }
            }
        }

        // Deferred bucketMove for all affected tuples
        for (auto tid : affectedTupleSet) {
            if (!tuplePeeled[tid]) bucketMove(tid);
        }
    }

    auto tEnd = std::chrono::high_resolution_clock::now();
    auto peelMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tPeel).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();
    printf("  Peeling: %lld iter, %lld LeafDeath, %lld BK\n", totalIters, cntLeafDeath, cntBK);
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
