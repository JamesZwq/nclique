//
// Region + ST (V4): ST peeling with tuple-level compression.
//
// Same correctness as ST. Key optimization: tuple batching.
// When popping an r-clique, also pop all r-cliques of its tuple
// (they share the same core value by Support Equality Theorem).
//
// Peeling: per-r-clique support tracking (same as ST).
// Bucket queue: per-r-clique (same as ST), but batch pop by tuple.
//

#include "NCliqueCoreDecomposition.h"
#include "../BK/BronKerboschRmRClique.hpp"
#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"
#include <algorithm>
#include <chrono>
#include <cstring>
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
    const daf::Size INVALID = std::numeric_limits<daf::Size>::max();

    // ================================================================
    // Step 1: Build overlap classes + tuples from maximal cliques
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
            }
            classOf[v] = pToC[vtxRegions[v]];
            classSizes[classOf[v]]++;
        }
    }

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

    // Enumerate tuples
    struct TupleInfo { TupleKey key; daf::Size mult; };
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
                for (auto &[c, k] : counts) {
                    if ((int)classSizes[c] < k) return;
                    mult *= (daf::Size)(nCr[classSizes[c]][k] + 0.5);
                }
                if (mult == 0) return;
                tupleIndex[cur] = tuples.size();
                tuples.push_back({cur, mult});
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

    printf("======= Region + ST (V4) =======\n");
    printf("  r=%d s=%d, Tuples: %zu, r-cliques: %zu (%.1fx compression)\n",
           (int)r, (int)s, tuples.size(), (size_t)totalRCliques,
           tuples.empty() ? 0.0 : (double)totalRCliques / tuples.size());

    // ================================================================
    // Step 2: Build cliqueIndex + support + r-clique → tuple mapping
    // ================================================================
    const daf::Size numLeaves = tree.adj_list.size();
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

    StaticCliqueIndex cliqueIndex(r);
    std::vector<double> countingRClique;

    daf::timeCount("fused build+counting (V4)", [&]() {
        cliqueIndex.buildWithFullEnum(tree, edgeGraph.adj_list.size(),
            [&](daf::Size leafIdx, StaticCliqueIndex::Id cliqueId,
                daf::CliqueSize subNumPivot, const uint8_t*) {
                if (cliqueId >= countingRClique.size())
                    countingRClique.resize(cliqueId + 1, 0);
                int rp = leafPivotC[leafIdx] - (int)subNumPivot;
                int rn = leafNeedPivot[leafIdx] - (int)subNumPivot;
                if (rp >= 0 && rn >= 0 && rp >= rn)
                    countingRClique[cliqueId] += nCr[rp][rn];
            });
    });

    // Map r-clique → tuple + build tuple → r-cliques index
    std::vector<daf::Size> rCliqueToTuple(cliqueIndex.size(), INVALID);
    std::vector<std::vector<daf::Size>> tupleRCliques(tuples.size());
    {
        TupleKey keyBuf; keyBuf.reserve(r);
        for (daf::Size id = 0; id < cliqueIndex.size(); ++id) {
            auto rc = cliqueIndex.byId(id);
            keyBuf.clear();
            bool valid = true;
            for (auto v : rc) {
                if (v >= numVertices || classOf[v] == INVALID) { valid = false; break; }
                keyBuf.push_back(classOf[v]);
            }
            if (!valid) continue;
            std::sort(keyBuf.begin(), keyBuf.end());
            auto it = tupleIndex.find(keyBuf);
            if (it != tupleIndex.end()) {
                rCliqueToTuple[id] = it->second;
                tupleRCliques[it->second].push_back(id);
            }
        }
    }

    auto tIndex = std::chrono::high_resolution_clock::now();
    printf("  cliqueIndex: %zu r-cliques, build+map: %lld ms\n",
        (size_t)cliqueIndex.size(),
        std::chrono::duration_cast<std::chrono::milliseconds>(tIndex - tStart).count());

    // ================================================================
    // Step 3: Peeling (ST mechanism + tuple batch pop)
    // ================================================================
    std::vector<double> coreTuple(tuples.size(), 0.0);
    std::vector<bool> tuplePeeled(tuples.size(), false);
    std::vector<bool> rCliqueAlive(cliqueIndex.size(), true);

    // Bucket queue (per-r-clique, same structure as ST)
    constexpr int BUCKET_MAX = 5000000;
    double rawMax = 0;
    for (daf::Size id = 0; id < cliqueIndex.size(); ++id)
        rawMax = std::max(rawMax, countingRClique[id]);
    int maxBucket = (int)std::min(rawMax, (double)BUCKET_MAX);

    int curBucket = 0;  // declared before lambdas that reference it
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::set<std::pair<double, daf::Size>> overflowSet;
    std::vector<int> bucket_of(cliqueIndex.size(), -1);
    std::vector<daf::Size> pos_in_bucket(cliqueIndex.size());
    std::vector<double> overflowStoredVal(cliqueIndex.size(), -1);

    auto bucketInsert = [&](daf::Size id) {
        double val = std::max(0.0, countingRClique[id]);
        int b = std::min((int)val, maxBucket + 1);
        if (b <= maxBucket) {
            pos_in_bucket[id] = buckets[b].size();
            buckets[b].push_back(id);
        } else {
            overflowSet.insert({val, id});
            overflowStoredVal[id] = val;
        }
        bucket_of[id] = b;
    };

    auto bucketMove = [&](daf::Size id) {
        if (!rCliqueAlive[id]) return;
        double val = std::max(0.0, countingRClique[id]);
        int oldB = bucket_of[id];
        int newB = std::min((int)val, maxBucket + 1);

        if (oldB == newB && (newB <= maxBucket || val == overflowStoredVal[id])) return;

        // Remove from old bucket
        if (oldB >= 0 && oldB <= maxBucket) {
            auto &bk = buckets[oldB];
            daf::Size myPos = pos_in_bucket[id];
            if (myPos < bk.size() - 1) {
                daf::Size last = bk.back();
                bk[myPos] = last;
                pos_in_bucket[last] = myPos;
            }
            bk.pop_back();
        } else if (oldB == maxBucket + 1) {
            overflowSet.erase({overflowStoredVal[id], id});
        }

        // Insert to new bucket
        if (newB <= maxBucket) {
            pos_in_bucket[id] = buckets[newB].size();
            buckets[newB].push_back(id);
            // CRITICAL: pull curBucket back if item moved to lower bucket
            if (newB < curBucket) curBucket = newB;
        } else {
            overflowSet.insert({val, id});
            overflowStoredVal[id] = val;
        }
        bucket_of[id] = newB;
    };

    for (daf::Size id = 0; id < cliqueIndex.size(); ++id) bucketInsert(id);

    daf::Size remainingRCliques = cliqueIndex.size();
    double minCore = 0;
    long long totalIters = 0, cntLeafDeath = 0, cntBK = 0;

    std::vector<daf::Size> changedLeaf;
    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), INVALID);
    std::vector<std::vector<daf::Size>> removedRCliqueIdForLeaf;
    std::vector<daf::Size> currentRemoveRcliqueIds;
    std::vector<daf::Size> vertexConflictDeg;

    auto tPeel = std::chrono::high_resolution_clock::now();

    while (remainingRCliques > 0) {
        for (auto &lid : changedLeaf) changedLeafIndex[lid] = INVALID;
        changedLeaf.clear();
        removedRCliqueIdForLeaf.clear();
        currentRemoveRcliqueIds.clear();

        // --- Pop r-cliques at current bucket level ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) {
            if (overflowSet.empty()) break;
            // Handle overflow
            auto it = overflowSet.begin();
            daf::Size id = it->second;
            overflowSet.erase(it);
            if (!rCliqueAlive[id]) continue;
            minCore = std::max(countingRClique[id], minCore);
            rCliqueAlive[id] = false;
            currentRemoveRcliqueIds.push_back(id);
            remainingRCliques--;
            // Batch pop: also pop all r-cliques of the same tuple
            daf::Size tid = rCliqueToTuple[id];
            if (tid != INVALID && !tuplePeeled[tid]) {
                tuplePeeled[tid] = true;
                coreTuple[tid] = minCore;
                for (auto rid : tupleRCliques[tid]) {
                    if (rCliqueAlive[rid]) {
                        rCliqueAlive[rid] = false;
                        // Remove from bucket/overflow
                        int b = bucket_of[rid];
                        if (b >= 0 && b <= maxBucket) {
                            auto &bk = buckets[b];
                            daf::Size p = pos_in_bucket[rid];
                            if (p < bk.size()-1) { daf::Size l = bk.back(); bk[p] = l; pos_in_bucket[l] = p; }
                            bk.pop_back();
                        } else if (b == maxBucket+1) {
                            overflowSet.erase({overflowStoredVal[rid], rid});
                        }
                        bucket_of[rid] = -1;
                        currentRemoveRcliqueIds.push_back(rid);
                        remainingRCliques--;
                    }
                }
            }
            goto phase1;
        }

        minCore = std::max((double)curBucket, minCore);
        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                daf::Size id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                if (!rCliqueAlive[id]) continue;
                rCliqueAlive[id] = false;
                currentRemoveRcliqueIds.push_back(id);
                remainingRCliques--;
                // Batch pop: also pop entire tuple
                daf::Size tid = rCliqueToTuple[id];
                if (tid != INVALID && !tuplePeeled[tid]) {
                    tuplePeeled[tid] = true;
                    coreTuple[tid] = minCore;
                    for (auto rid : tupleRCliques[tid]) {
                        if (rCliqueAlive[rid]) {
                            rCliqueAlive[rid] = false;
                            int b = bucket_of[rid];
                            if (b >= 0 && b <= maxBucket) {
                                auto &bk = buckets[b];
                                daf::Size p = pos_in_bucket[rid];
                                if (p < bk.size()-1) { daf::Size l = bk.back(); bk[p] = l; pos_in_bucket[l] = p; }
                                bk.pop_back();
                            } else if (b == maxBucket+1) {
                                overflowSet.erase({overflowStoredVal[rid], rid});
                            }
                            bucket_of[rid] = -1;
                            currentRemoveRcliqueIds.push_back(rid);
                            remainingRCliques--;
                        }
                    }
                }
            }
            if (curBucket+1 < (int)buckets.size() && !buckets[curBucket+1].empty()
                && (curBucket+1) <= (int)minCore) curBucket++;
            else break;
        }

        phase1:
        if (remainingRCliques == 0) break;
        totalIters++;

        // === Phase 1: Find affected leaves (same as ST) ===
        for (auto rmId : currentRemoveRcliqueIds) {
            auto rClique = cliqueIndex.byId(rmId);
            daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                [&](const TreeGraphNode &uClique) {
                    auto &idx = changedLeafIndex[uClique.v];
                    if (idx == INVALID) {
                        idx = removedRCliqueIdForLeaf.size();
                        removedRCliqueIdForLeaf.emplace_back();
                        changedLeaf.push_back(uClique.v);
                    }
                    removedRCliqueIdForLeaf[idx].push_back(rmId);
                });
        }

        // === Phase 2: LeafDeath / BK (same as ST) ===
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            auto leafId = changedLeaf[idx];
            const auto &leaf = tree.adj_list[leafId];
            if (leaf.empty()) continue;
            auto &removedR = removedRCliqueIdForLeaf[changedLeafIndex[leafId]];

            int keepC = 0, pivotC = 0;
            for (const auto &node : leaf) {
                if (node.isPivot) pivotC++; else keepC++;
            }
            int needPivot = s - keepC;
            int n = (int)leaf.size();

            daf::Size maxRCliquePerVertex = (daf::Size)(nCr[n-1][r-1] + 0.5);
            vertexConflictDeg.assign(n, 0);
            for (int i = 0; i < n; ++i) daf::vListMap[leaf[i].v] = (daf::Size)i;
            for (auto rmId : removedR) {
                auto rc = cliqueIndex.byId(rmId);
                for (auto v : rc) {
                    daf::Size pos = daf::vListMap[v];
                    if (pos < (daf::Size)n) vertexConflictDeg[pos]++;
                }
            }

            bool leafDead = false;
            int forcedPivotRemove = 0;
            for (int i = 0; i < n; ++i) {
                if (vertexConflictDeg[i] >= maxRCliquePerVertex) {
                    if (!leaf[i].isPivot) { leafDead = true; break; }
                    forcedPivotRemove++;
                }
            }
            if (!leafDead) {
                int rp = pivotC - forcedPivotRemove;
                int rt = n - forcedPivotRemove;
                if (rt < (int)s || rp < needPivot) leafDead = true;
            }

            if (leafDead) {
                cntLeafDeath++;
                daf::enumerateCombinations(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &clique) {
                        auto cid = cliqueIndex.byClique(clique);
                        if (!rCliqueAlive[cid]) return true;
                        daf::CliqueSize subP = 0;
                        for (const auto &nd : clique) if (nd.isPivot) subP++;
                        countingRClique[cid] -= nCr[pivotC - subP][needPivot - subP];
                        if (countingRClique[cid] < 0) countingRClique[cid] = 0;
                        bucketMove(cid);
                        return true;
                    });
                for (const auto &i : leaf)
                    treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                tree.adj_list[leafId].clear();
                tree.recycleNode(leafId);
            } else {
                cntBK++;
                for (const auto &lv : leaf) {
                    if (lv.isPivot) treeGraphV.removeNbr(lv.v, {leafId, true});
                    else treeGraphV.removeNbr(lv.v, {leafId, false});
                }

                auto mapped = removedR | std::views::transform(
                    [&](daf::Size id) { return cliqueIndex.byId(id); });

                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                        auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                        auto newId = tree.addNode(newLeaf);
                        const auto &stored = tree.adj_list[newId];
                        daf::CliqueSize newPC = 0, newKC = 0;
                        for (const auto &i : stored) {
                            if (i.isPivot) { treeGraphV.addNbr(i.v, {newId, true}); newPC++; }
                            else { treeGraphV.addNbr(i.v, {newId, false}); newKC++; }
                        }
                        int np = s - newKC;
                        daf::enumerateCombinations(stored, r,
                            [&](const daf::StaticVector<TreeGraphNode> &rclique) {
                                daf::CliqueSize subP = 0;
                                for (const auto &nd : rclique) if (nd.isPivot) subP++;
                                if (subP <= np && newPC - subP < 1001 && np - subP < 401)
                                    countingRClique[cliqueIndex.byClique(rclique)] += nCr[newPC - subP][np - subP];
                                return true;
                            });
                        if (newId >= changedLeafIndex.size())
                            changedLeafIndex.resize(newId * 2, INVALID);
                        // Extend leaf arrays
                        if (newId >= leafPivotC.size()) {
                            leafPivotC.resize(newId+1, 0);
                            leafNeedPivot.resize(newId+1, 0);
                        }
                        leafPivotC[newId] = newPC;
                        leafNeedPivot[newId] = s - newKC;
                    });

                // Subtract old leaf
                daf::enumerateCombinations(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &clique) {
                        auto cid = cliqueIndex.byClique(clique);
                        if (!rCliqueAlive[cid]) return true;
                        daf::CliqueSize subP = 0;
                        for (const auto &nd : clique) if (nd.isPivot) subP++;
                        countingRClique[cid] -= nCr[pivotC - subP][needPivot - subP];
                        if (countingRClique[cid] < 0) countingRClique[cid] = 0;
                        bucketMove(cid);
                        return true;
                    });

                tree.removeNode(leafId);
            }
        }
    }

    auto tEnd = std::chrono::high_resolution_clock::now();
    auto peelMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tPeel).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();
    printf("  Peeling: %lld iter, %lld LeafDeath, %lld BK\n", totalIters, cntLeafDeath, cntBK);
    printf("  Peeling time: %lld ms\n", peelMs);
    printf("  Total time: %lld ms\n", totalMs);

    // ================================================================
    // Result: per-r-clique core values via tuple
    // ================================================================
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    std::map<double, daf::Size> coreDist;
    for (daf::Size id = 0; id < cliqueIndex.size(); ++id) {
        daf::Size tid = rCliqueToTuple[id];
        double core = (tid != INVALID) ? coreTuple[tid] : 0.0;
        auto rc = cliqueIndex.byId(id);
        result.push_back({std::vector<daf::Size>(rc.begin(), rc.end()), core});
        coreDist[core]++;
    }
    double maxCore = 0;
    for (auto &[c, cnt] : coreDist) maxCore = std::max(maxCore, c);
    printf("  Max core: %.0f\n", maxCore);
    for (auto &[c, cnt] : coreDist)
        printf("  core=%.0f count=%zu\n", c, (size_t)cnt);

    return result;
}
