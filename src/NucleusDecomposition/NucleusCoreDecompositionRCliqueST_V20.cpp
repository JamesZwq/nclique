//
// ST V20: Sequential d=1 Split — NO BK recursion
//
// Key insight: Process removals ONE AT A TIME. Each removal is d=1:
//   - Delta: ONE nCr per r-clique (analytical closed-form)
//   - Split: Theorem 1 deterministic construction (t ≤ r subpaths, no recursion)
// After splitting, future removals hit SMALLER subpaths → cheaper.
//
// Compare pathSplit: O(Π|F_i|) recursive subpaths with backtracking.
// V20: O(k × C(n_avg, r)) where n_avg shrinks after each split.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <cstring>
#include <set>

#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V20(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    long long duration_init = 0, duration_pop = 0, duration_intersect = 0;
    long long duration_split = 0, duration_support = 0;
    long long cntLeafDeath = 0, cntSplit = 0;

    auto time_start = std::chrono::high_resolution_clock::now();

    // ========== INIT ==========
    StaticCliqueIndex localIndex(r);
    StaticCliqueIndex &cliqueIndex = prebuiltIndex ? *prebuiltIndex : localIndex;
    if (!prebuiltIndex) {
        daf::timeCount("clique Index build", [&]() {
            localIndex.build(tree, edgeGraph.adj_list.size());
        });
    }

    std::vector<double> countingRClique;
    daf::timeCount("countingPerRClique (V20)", [&]() {
        countingRClique.assign(cliqueIndex.size(), 0.0);
        for (const auto &leaf : tree.adj_list) {
            if (leaf.size() < r) continue;
            daf::CliqueSize pivotC = 0, keepC = 0;
            for (const auto &i : leaf) { if (i.isPivot) pivotC++; else keepC++; }
            int needPivot = s - (int)keepC;
            daf::enumerateCombinations(leaf, r,
                [&](const daf::StaticVector<TreeGraphNode> &rc) {
                daf::CliqueSize sp = 0;
                for (const auto &nd : rc) if (nd.isPivot) sp++;
                if (sp <= needPivot) {
                    int row = (int)pivotC-(int)sp, col = needPivot-(int)sp;
                    if (row>=0 && row<1001 && col>=0 && col<401)
                        countingRClique[cliqueIndex.byClique(rc)] += nCr[row][col];
                }
                return true;
            });
        }
    });

    const daf::Size nClique = cliqueIndex.size();
    daf::log_memory("r-clique index + counting");

    std::vector<double> coreRClique(nClique, 0);
    daf::StaticVector<bool> rCliqueInHeap(nClique);
    rCliqueInHeap.resize(nClique);
    memset(rCliqueInHeap.getData(), true, nClique * sizeof(bool));

    // ========== Hybrid bucket+set PQ ==========
    constexpr double BUCKET_THRESHOLD = 5000000.0;
    double rawMax = 0;
    for (daf::Size i = 0; i < nClique; ++i) rawMax = std::max(rawMax, countingRClique[i]);
    int maxBucket = (int)std::min(rawMax, BUCKET_THRESHOLD);
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::set<std::pair<double, daf::Size>> overflowSet;
    std::vector<int> bucket_of(nClique, -1);
    std::vector<daf::Size> pos_in_bucket(nClique);
    std::vector<double> overflowStoredVal(nClique, -1);
    for (daf::Size i = 0; i < nClique; ++i) {
        if (countingRClique[i] <= BUCKET_THRESHOLD) {
            int b=(int)countingRClique[i]; bucket_of[i]=b;
            pos_in_bucket[i]=buckets[b].size(); buckets[b].push_back(i);
        } else {
            overflowSet.insert({countingRClique[i],i}); overflowStoredVal[i]=countingRClique[i];
        }
    }
    int curBucket = 0;
    daf::Size remainingInHeap = nClique;

    auto bucketMove = [&](daf::Size id) {
        if (!rCliqueInHeap[id]) return;
        double val = std::max(0.0, countingRClique[id]);
        int oldB = bucket_of[id];
        if (oldB==-1) overflowSet.erase({overflowStoredVal[id],id});
        if (val<=BUCKET_THRESHOLD) {
            int newB=(int)val;
            if (oldB>=0&&newB==oldB) return;
            if (oldB>=0) { auto&v=buckets[oldB]; auto p=pos_in_bucket[id];
                if(p<v.size()-1){auto l=v.back();v[p]=l;pos_in_bucket[l]=p;} v.pop_back(); }
            bucket_of[id]=newB; pos_in_bucket[id]=buckets[newB].size();
            buckets[newB].push_back(id); if(newB<curBucket) curBucket=newB;
        } else {
            if (oldB>=0) { auto&v=buckets[oldB]; auto p=pos_in_bucket[id];
                if(p<v.size()){if(p<v.size()-1){auto l=v.back();v[p]=l;pos_in_bucket[l]=p;} v.pop_back();} }
            overflowSet.insert({val,id}); overflowStoredVal[id]=val; bucket_of[id]=-1;
        }
    };
    auto drainOverflow = [&]() {
        while (!overflowSet.empty()) {
            auto id=overflowSet.begin()->second;
            if (!rCliqueInHeap[id]) { overflowSet.erase(overflowSet.begin()); continue; }
            if (countingRClique[id]<=BUCKET_THRESHOLD) {
                overflowSet.erase(overflowSet.begin());
                int b=(int)countingRClique[id]; bucket_of[id]=b;
                pos_in_bucket[id]=buckets[b].size(); buckets[b].push_back(id);
            } else break;
        }
    };

    daf::log_memory("Other index");
    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::cout << "=========================begin (r>=3 ST_V20)===========================" << std::endl;
    double minCore = 0;
    long long totalIters = 0;

    // Reusable buffers
    std::vector<daf::Size> currentRemoveIds;
    currentRemoveIds.reserve(nClique);

    while (remainingInHeap > 0) {
        auto t0 = std::chrono::high_resolution_clock::now();
        currentRemoveIds.clear();

        // --- Pop ---
        drainOverflow();
        while (curBucket<(int)buckets.size()&&buckets[curBucket].empty()) curBucket++;
        if (curBucket>=(int)buckets.size()) {
            if (!overflowSet.empty()) {
                while (!overflowSet.empty()) {
                    auto id=overflowSet.begin()->second; overflowSet.erase(overflowSet.begin());
                    if (!rCliqueInHeap[id]) continue;
                    minCore=std::max(countingRClique[id],minCore);
                    rCliqueInHeap[id]=false; currentRemoveIds.push_back(id);
                    coreRClique[id]=minCore; remainingInHeap--;
                    while (!overflowSet.empty()) {
                        auto nxt=overflowSet.begin()->second;
                        if (!rCliqueInHeap[nxt]){overflowSet.erase(overflowSet.begin());continue;}
                        if (countingRClique[nxt]<=minCore) {
                            overflowSet.erase(overflowSet.begin());
                            rCliqueInHeap[nxt]=false; currentRemoveIds.push_back(nxt);
                            coreRClique[nxt]=minCore; remainingInHeap--;
                        } else break;
                    }
                    break;
                }
                if (currentRemoveIds.empty()) break;
                goto pop_done20;
            }
            break;
        }
        minCore=std::max((double)curBucket,minCore);
        while (curBucket<(int)buckets.size()&&!buckets[curBucket].empty()&&curBucket<=(int)minCore) {
            while (!buckets[curBucket].empty()) {
                auto id=buckets[curBucket].back(); buckets[curBucket].pop_back();
                rCliqueInHeap[id]=false; currentRemoveIds.push_back(id);
                coreRClique[id]=minCore; remainingInHeap--;
            }
            if (curBucket+1<(int)buckets.size()&&!buckets[curBucket+1].empty()&&(curBucket+1)<=(int)minCore)
                curBucket++; else break;
        }
        pop_done20:
        duration_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now()-t0).count();
        if (remainingInHeap==0) break;
        totalIters++;

        // --- Process each removed r-clique SEQUENTIALLY (d=1 each) ---
        for (auto rmId : currentRemoveIds) {
            auto rmClique = cliqueIndex.byId(rmId); // the removed r-clique's vertices

            // Find all leaves containing this r-clique via treeGraphV
            // (intersect the adjacency lists of all vertices in rmClique)
            auto t1 = std::chrono::high_resolution_clock::now();
            std::vector<daf::Size> affectedLeaves;
            daf::intersect_dense_sets_multi(rmClique, treeGraphV.adj_list,
                [&](const TreeGraphNode &u) {
                    affectedLeaves.push_back(u.v);
                });
            duration_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now()-t1).count();

            // For each affected leaf: d=1 analytical delta + Theorem 1 split
            for (auto leafId : affectedLeaves) {
                const auto &leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;

                int n = (int)leaf.size();
                if (n < (int)s) continue;

                daf::CliqueSize pivotC = 0, keepC = 0;
                for (const auto &nd : leaf) { if (nd.isPivot) pivotC++; else keepC++; }
                int needPivot = s - (int)keepC;
                int a = (int)pivotC;

                daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
                for (int i = 0; i < n; ++i) mapRef[leaf[i].v] = (daf::Size)i;

                // Build pivot mask for the removed r-clique on this leaf
                uint64_t rmMask = 0;
                int rmPivotCount = 0;
                std::vector<int> rmPivotPositions; // leaf-local positions of rm's pivots
                for (auto v : rmClique) {
                    daf::Size pos = mapRef[v];
                    if (pos < (daf::Size)n && leaf[pos].isPivot) {
                        rmMask |= (1ULL << pos);
                        rmPivotPositions.push_back((int)pos);
                        rmPivotCount++;
                    }
                }

                if (rmPivotCount == 0) {
                    // All vertices of rm are holds → leaf must die (can't avoid rm)
                    cntLeafDeath++;
                    auto t_s = std::chrono::high_resolution_clock::now();
                    for (const auto &v : leaf) treeGraphV.removeNbr(v.v, {leafId, v.isPivot});
                    // Subtract ALL support from this leaf
                    daf::enumerateCombinations(leaf, r,
                        [&](const daf::StaticVector<TreeGraphNode> &rc) {
                        auto cid = cliqueIndex.byClique(rc);
                        if (!rCliqueInHeap[cid]) return true;
                        daf::CliqueSize sp=0;
                        for (const auto &nd : rc) if(nd.isPivot) sp++;
                        countingRClique[cid] -= nCr[a-sp][needPivot-sp];
                        if (countingRClique[cid]<0) countingRClique[cid]=0;
                        bucketMove(cid);
                        return true;
                    });
                    tree.removeNode(leafId);
                    duration_support += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now()-t_s).count();
                    continue;
                }

                auto t_sp = std::chrono::high_resolution_clock::now();
                cntSplit++;

                // ===== STEP 1: Analytical d=1 delta =====
                // For each r-clique q' on this leaf:
                //   delta(q') = -nCr(a - |piv(q' ∪ rm)|, needPivot - |piv(q' ∪ rm)|)
                daf::enumerateCombinations(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &rc) {
                    auto cid = cliqueIndex.byClique(rc);
                    if (!rCliqueInHeap[cid]) return true;
                    uint64_t qMask = 0;
                    for (const auto &nd : rc)
                        if (nd.isPivot) qMask |= (1ULL << mapRef[nd.v]);
                    uint64_t combined = qMask | rmMask;
                    int u = __builtin_popcountll(combined);
                    int row = a - u, col = needPivot - u;
                    if (row >= 0 && col >= 0 && row < 1001 && col < 401) {
                        countingRClique[cid] -= nCr[row][col];
                        if (countingRClique[cid] < 0) countingRClique[cid] = 0;
                        bucketMove(cid);
                    }
                    return true;
                });

                // ===== STEP 2: Theorem 1 deterministic split =====
                // Remove old leaf from treeGraphV
                for (const auto &v : leaf) treeGraphV.removeNbr(v.v, {leafId, v.isPivot});

                // Generate t = rmPivotCount subpaths
                std::sort(rmPivotPositions.begin(), rmPivotPositions.end());
                int t = rmPivotCount;

                for (int i = 0; i < t; ++i) {
                    // Subpath i: promote positions [0..i-1] to hold, remove position [i]
                    std::vector<TreeGraphNode> subpath;
                    daf::CliqueSize newKeepC = 0, newPivotC = 0;
                    for (int p = 0; p < n; ++p) {
                        if (p == rmPivotPositions[i]) continue; // removed
                        bool isPiv = leaf[p].isPivot;
                        // Check if this position should be promoted to hold
                        for (int j = 0; j < i; ++j) {
                            if (p == rmPivotPositions[j]) { isPiv = false; break; }
                        }
                        subpath.push_back({leaf[p].v, isPiv});
                        if (isPiv) newPivotC++; else newKeepC++;
                    }

                    // Prune: subpath too small?
                    if ((int)(newKeepC + newPivotC) < (int)s) continue;
                    if ((int)newKeepC > (int)s) continue;
                    int newNeedPivot = s - newKeepC;
                    if ((int)newPivotC < newNeedPivot) continue;

                    auto newId = tree.addNode(subpath);
                    for (const auto &v : tree.adj_list[newId])
                        treeGraphV.addNbr(v.v, {newId, v.isPivot});
                }

                tree.removeNode(leafId);

                duration_split += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now()-t_sp).count();
            }
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "time: " << elapsed << " ms" << std::endl;
    std::cout << "Time Breakdown (ms):" << std::endl;
    std::cout << "  Init:      " << duration_init / 1e6 << std::endl;
    std::cout << "  Pop:       " << duration_pop / 1e6 << std::endl;
    std::cout << "  Intersect: " << duration_intersect / 1e6 << std::endl;
    std::cout << "  Split:     " << duration_split / 1e6 << std::endl;
    std::cout << "  Support:   " << duration_support / 1e6 << std::endl;
    std::cout << "  Cases: LeafDeath=" << cntLeafDeath << " Split=" << cntSplit
              << " iters=" << totalIters << std::endl;

    std::vector<std::pair<std::vector<daf::Size>, double>> sortedK;
    sortedK.reserve(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        auto c = cliqueIndex.byId(i);
        sortedK.emplace_back(std::vector<daf::Size>(c.begin(), c.end()), coreRClique[i]);
    }
    return sortedK;
}
