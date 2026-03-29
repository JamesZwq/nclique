//
// Single-thread r=1 nucleus decomposition — immutable tree + per-leaf counters.
//
// Tree and treeGraphV are NEVER modified. Each leaf maintains only two counters
// (remainPivots, alive). Batch peel + per-leaf delta via nCr lookup.
//
// vs old ST: eliminates leafRmInfo, removedPivots sort/dedup, tree.removeNbrs(),
//            tree.recycleNode(). Just integer counters + nCr table lookups.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <boost/heap/d_ary_heap.hpp>

#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace VCD_ST {
    double * countingPerVertex(const DynamicGraph<TreeGraphNode> &treeGraph,
                               const Graph &edgeGraph,
                               const daf::CliqueSize k) {
        const daf::Size n = edgeGraph.adj_list_offsets.size();
        auto *countingV = new double[n];
        memset(countingV, 0, n * sizeof(double));
        daf::StaticVector<daf::Size> povit;
        daf::StaticVector<daf::Size> keepC;
        for (const auto &clique : treeGraph.adj_list) {
            povit.clear();
            keepC.clear();
            if (clique.size() < k) continue;
            for (auto &i : clique) {
                if (i.isPivot) povit.push_back(i.v);
                else keepC.push_back(i.v);
            }
            int needPivot = int(k) - int(keepC.size());
            double keepVal = nCr[povit.size()][needPivot];
            for (const auto &v : keepC)
                countingV[v] += keepVal;
            double pivotVal = 0;
            const int npv = needPivot - 1;
            if (npv >= 0 && npv <= static_cast<int>(povit.size()) - 1)
                pivotVal = nCr[povit.size() - 1][npv];
            for (const auto &v : povit)
                countingV[v] += pivotVal;
        }
        povit.free();
        keepC.free();
        return countingV;
    }
}

double * NCliqueVertexCoreDecomposition_ST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    // --- Precompute per-leaf metadata ---
    const daf::Size numLeaves = tree.adj_list.size();
    std::vector<int> leafRemainPivots(numLeaves);
    std::vector<int> leafNeedPivot(numLeaves);
    std::vector<uint8_t> leafAlive(numLeaves);

    for (daf::Size L = 0; L < numLeaves; ++L) {
        int keeps = 0, pivots = 0;
        for (auto &node : tree.adj_list[L]) {
            if (node.isPivot) pivots++;
            else keeps++;
        }
        leafRemainPivots[L] = pivots;
        leafNeedPivot[L] = (int)k - keeps;
        leafAlive[L] = (leafNeedPivot[L] >= 0 && leafNeedPivot[L] <= pivots) ? 1 : 0;
    }

    // --- H3: Build flat CSR vertex→leaf index (replaces treeGraphV hash set iteration) ---
    // TreeGraphNode is 8B (63-bit v + 1-bit isPivot). Store as compact pair<Size,uint8_t>.
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;
    struct VLeafEntry { daf::Size leafId; uint8_t isPivot; };
    std::vector<daf::Size> vtxLeafOff(numVertices + 2, 0); // CSR offsets
    // Count phase
    for (daf::Size L = 0; L < numLeaves; ++L) {
        for (const auto &node : tree.adj_list[L]) {
            if (node.v < numVertices) vtxLeafOff[node.v + 1]++;
        }
    }
    // Prefix sum
    for (daf::Size i = 1; i <= numVertices + 1; ++i)
        vtxLeafOff[i] += vtxLeafOff[i - 1];
    // Fill phase
    std::vector<VLeafEntry> vtxLeafData(vtxLeafOff[numVertices]);
    std::vector<daf::Size> vtxLeafPos(numVertices + 1, 0); // write cursor
    for (daf::Size L = 0; L < numLeaves; ++L) {
        for (const auto &node : tree.adj_list[L]) {
            daf::Size v = node.v;
            if (v < numVertices) {
                daf::Size pos = vtxLeafOff[v] + vtxLeafPos[v]++;
                vtxLeafData[pos] = {L, (uint8_t)node.isPivot};
            }
        }
    }
    vtxLeafPos.clear();
    vtxLeafPos.shrink_to_fit();

    // --- Initial per-vertex support ---
    auto countingV = VCD_ST::countingPerVertex(tree, edgeGraph, k);
    auto coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;

    // --- Hybrid Bucket + d-ary heap overflow ---
    constexpr double BUCKET_THRESHOLD = 5000000.0;
    int maxBucket = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] > 0 && countingV[i] <= BUCKET_THRESHOLD)
            maxBucket = std::max(maxBucket, (int)countingV[i]);
    }
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(numVertices, -1);
    std::vector<daf::Size> pos_in_bucket(numVertices);

    struct CmpVtx {
        const double *cnt;
        bool operator()(daf::Size a, daf::Size b) const { return cnt[a] > cnt[b]; }
    };
    using OverflowHeap = boost::heap::d_ary_heap<daf::Size, boost::heap::arity<8>,
        boost::heap::mutable_<true>, boost::heap::compare<CmpVtx>>;
    OverflowHeap overflowHeap{CmpVtx{countingV}};
    std::vector<OverflowHeap::handle_type> overflowHandles(numVertices);
    std::vector<uint8_t> inOverflow(numVertices, 0);

    std::vector<uint8_t> vertexInHeap(numVertices, 0);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] <= 0) continue;
        if (countingV[i] <= BUCKET_THRESHOLD) {
            int b = (int)countingV[i];
            bucket_of[i] = b;
            pos_in_bucket[i] = buckets[b].size();
            buckets[b].push_back(i);
        } else {
            overflowHandles[i] = overflowHeap.push(i);
            inOverflow[i] = 1;
        }
        vertexInHeap[i] = 1;
        remainingInHeap++;
    }
    int curBucket = 0;

    auto bucketMove = [&](daf::Size id) {
        if (!vertexInHeap[id]) return;
        double val = std::max(0.0, countingV[id]);
        int oldB = bucket_of[id];
        if (oldB == -1 && inOverflow[id]) {
            overflowHeap.update(overflowHandles[id]);
            if (val <= BUCKET_THRESHOLD) {
                overflowHeap.erase(overflowHandles[id]);
                inOverflow[id] = 0;
                int newB = (int)val;
                if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
                bucket_of[id] = newB; pos_in_bucket[id] = buckets[newB].size();
                buckets[newB].push_back(id); if (newB < curBucket) curBucket = newB;
            }
            return;
        }
        if (val <= BUCKET_THRESHOLD) {
            int newB = (int)val;
            if (oldB == newB) return;
            if (oldB >= 0) {
                auto &oldVec = buckets[oldB];
                daf::Size myPos = pos_in_bucket[id];
                if (myPos < oldVec.size() - 1) {
                    daf::Size last = oldVec.back();
                    oldVec[myPos] = last;
                    pos_in_bucket[last] = myPos;
                }
                oldVec.pop_back();
            }
            if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
            bucket_of[id] = newB;
            pos_in_bucket[id] = buckets[newB].size();
            buckets[newB].push_back(id);
            if (newB < curBucket) curBucket = newB;
        } else {
            if (oldB >= 0) {
                auto &oldVec = buckets[oldB];
                daf::Size myPos = pos_in_bucket[id];
                if (myPos < oldVec.size() - 1) {
                    daf::Size last = oldVec.back();
                    oldVec[myPos] = last;
                    pos_in_bucket[last] = myPos;
                }
                oldVec.pop_back();
            }
            if (!inOverflow[id]) {
                overflowHandles[id] = overflowHeap.push(id);
                inOverflow[id] = 1;
            } else {
                overflowHeap.update(overflowHandles[id]);
            }
            bucket_of[id] = -1;
        }
    };

    // --- Per-leaf batch tracking (reused each iteration) ---
    // leafRemovedPivots[L] = how many pivots removed this round
    // leafDies[L] = 1 if any keepC removed this round
    std::vector<int> leafRemovedPivots(numLeaves, 0);
    std::vector<uint8_t> leafDies(numLeaves, 0);
    std::vector<uint8_t> leafAffected(numLeaves, 0);
    std::vector<daf::Size> affectedLeaves;
    affectedLeaves.reserve(4096);

    // Dirty vertex tracking for deferred bucket moves
    std::vector<uint8_t> dirtyMark(numVertices, 0);
    std::vector<daf::Size> dirtyVertices;
    dirtyVertices.reserve(8192);

    daf::StaticVector<daf::Size> currentRemoveVertexIds(numVertices);

    std::cout << "=========================begin=========================" << std::endl;
    std::cout << "vertices in heap: " << remainingInHeap << std::endl;

    double minCore = 0;

    while (remainingInHeap > 0) {
        // --- Bucket pop: batch all vertices at current level ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;

        bool fromOverflow = false;
        if (curBucket >= (int)buckets.size()) {
            // Drain overflow heap → buckets
            while (!overflowHeap.empty()) {
                daf::Size oid = overflowHeap.top();
                if (!vertexInHeap[oid]) { overflowHeap.pop(); inOverflow[oid]=0; continue; }
                if (countingV[oid] <= BUCKET_THRESHOLD) {
                    overflowHeap.pop(); inOverflow[oid]=0;
                    int b = std::max(0, (int)countingV[oid]);
                    if (b >= (int)buckets.size()) buckets.resize(b + 1);
                    bucket_of[oid] = b; pos_in_bucket[oid] = buckets[b].size();
                    buckets[b].push_back(oid);
                    continue;
                }
                break;
            }
            // Try overflow min
            if (!overflowHeap.empty()) {
                daf::Size oid = overflowHeap.top();
                overflowHeap.pop(); inOverflow[oid]=0;
                minCore = std::max(countingV[oid], minCore);
                vertexInHeap[oid] = 0;
                currentRemoveVertexIds.push_back(oid);
                coreV[oid] = minCore;
                remainingInHeap--;
                // Drain more at same level
                while (!overflowHeap.empty()) {
                    daf::Size nid = overflowHeap.top();
                    if (!vertexInHeap[nid]) { overflowHeap.pop(); inOverflow[nid]=0; continue; }
                    if (countingV[nid] <= minCore) {
                        overflowHeap.pop(); inOverflow[nid]=0;
                        vertexInHeap[nid] = 0;
                        currentRemoveVertexIds.push_back(nid);
                        coreV[nid] = minCore;
                        remainingInHeap--;
                    } else break;
                }
                fromOverflow = true;
            }
            if (!fromOverflow) break;
            // Re-check buckets
            curBucket = 0;
            while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        }

        if (!fromOverflow) {
            minCore = std::max((double)curBucket, minCore);

            while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
                   && (double)curBucket <= minCore) {
                while (!buckets[curBucket].empty()) {
                    auto id = buckets[curBucket].back();
                    buckets[curBucket].pop_back();
                    vertexInHeap[id] = 0;
                    currentRemoveVertexIds.push_back(id);
                    coreV[id] = minCore;
                    remainingInHeap--;
                }
                if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                    && (double)(curBucket + 1) <= minCore) {
                    curBucket++;
                } else break;
            }
        }

        if (remainingInHeap == 0) break;

        // --- Phase 1: find affected leaves via flat CSR index (H3) ---
        for (int vi = 0; vi < (int)currentRemoveVertexIds.size(); ++vi) {
            auto v = currentRemoveVertexIds[vi];
            const daf::Size begin = vtxLeafOff[v];
            const daf::Size end = vtxLeafOff[v + 1];
            // H1: prefetch next vertex's CSR region
            if (vi + 1 < (int)currentRemoveVertexIds.size()) {
                auto nextV = currentRemoveVertexIds[vi + 1];
                __builtin_prefetch(&vtxLeafData[vtxLeafOff[nextV]], 0, 1);
            }
            for (daf::Size ei = begin; ei < end; ++ei) {
                const auto &entry = vtxLeafData[ei];
                daf::Size leafId = entry.leafId;
                if (!leafAlive[leafId]) continue;
                if (!leafAffected[leafId]) {
                    leafAffected[leafId] = 1;
                    affectedLeaves.push_back(leafId);
                }
                if (!entry.isPivot) {
                    leafDies[leafId] = 1;
                } else {
                    leafRemovedPivots[leafId]++;
                }
            }
        }

        // --- Phase 2: process each affected leaf ONCE ---
        for (auto leafId : affectedLeaves) {
            int old_rp = leafRemainPivots[leafId];
            int np = leafNeedPivot[leafId];
            int new_rp = old_rp - leafRemovedPivots[leafId];
            bool dies = leafDies[leafId] || np > new_rp || new_rp < 0;

            double deltaKeep, deltaPivot;
            if (dies) {
                deltaKeep = (np >= 0 && np <= old_rp) ? nCr[old_rp][np] : 0.0;
                deltaPivot = (np >= 1 && old_rp >= 1) ? nCr[old_rp - 1][np - 1] : 0.0;
            } else {
                deltaKeep = nCr[old_rp][np] - nCr[new_rp][np];
                deltaPivot = (np >= 1) ? nCr[old_rp - 1][np - 1] - nCr[new_rp - 1][np - 1] : 0.0;
            }

            for (auto &node : tree.adj_list[leafId]) {
                if (!vertexInHeap[node.v]) continue;
                double delta = node.isPivot ? deltaPivot : deltaKeep;
                if (delta > 0) {
                    countingV[node.v] -= delta;
                    if (countingV[node.v] < 0) countingV[node.v] = 0;
                    if (!dirtyMark[node.v]) {
                        dirtyMark[node.v] = 1;
                        dirtyVertices.push_back(node.v);
                    }
                }
            }

            leafRemainPivots[leafId] = new_rp;
            if (dies) leafAlive[leafId] = 0;
        }

        // --- Batch bucket moves ---
        for (auto v : dirtyVertices) {
            if (vertexInHeap[v]) bucketMove(v);
            dirtyMark[v] = 0;
        }
        dirtyVertices.clear();

        // --- Cleanup ---
        for (auto leafId : affectedLeaves) {
            leafRemovedPivots[leafId] = 0;
            leafDies[leafId] = 0;
            leafAffected[leafId] = 0;
        }
        affectedLeaves.clear();
        currentRemoveVertexIds.clear();
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    delete[] countingV;
    currentRemoveVertexIds.free();
    return coreV;
}
