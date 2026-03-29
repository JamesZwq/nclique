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
        for (daf::Size i = 0; i < n; ++i)
            countingV[i] = std::round(countingV[i]);
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

    // --- Pure d-ary heap PQ (same as REF) ---
    struct CmpVtx {
        const double *cnt;
        bool operator()(daf::Size a, daf::Size b) const { return cnt[a] > cnt[b]; }
    };
    using VtxHeap = boost::heap::d_ary_heap<daf::Size, boost::heap::arity<8>,
        boost::heap::mutable_<true>, boost::heap::compare<CmpVtx>>;
    VtxHeap pq{CmpVtx{countingV}};
    pq.reserve(numVertices);
    std::vector<VtxHeap::handle_type> pqHandles(numVertices);

    std::vector<uint8_t> vertexInHeap(numVertices, 0);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] <= 0) continue;
        pqHandles[i] = pq.push(i);
        vertexInHeap[i] = 1;
        remainingInHeap++;
    }

    auto bucketMove = [&](daf::Size id) {
        if (!vertexInHeap[id]) return;
        pq.update(pqHandles[id]);
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
        // --- Heap pop: batch all vertices at current min level ---
        while (!pq.empty() && !vertexInHeap[pq.top()]) pq.pop();
        if (pq.empty()) break;

        minCore = std::max(countingV[pq.top()], minCore);
        while (!pq.empty() && countingV[pq.top()] <= minCore) {
            auto id = pq.top(); pq.pop();
            if (!vertexInHeap[id]) continue;
            vertexInHeap[id] = 0;
            currentRemoveVertexIds.push_back(id);
            coreV[id] = minCore;
            remainingInHeap--;
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
                    countingV[node.v] = std::max(std::round(countingV[node.v]), 0.0);
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
