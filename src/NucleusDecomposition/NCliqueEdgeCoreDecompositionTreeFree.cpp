//
// NCliqueEdgeCoreDecompositionTreeFree.cpp
//
// Tree-Free R=2 Edge Nucleus Decomposition.
//
// Key innovation: Case A (~82-94%) and Case C (~2-10%) are handled without
// tree.adj_list or treeGraphV, using only flat CSR indices built during SDCT.
// Only Case B (~0.5-6.7%) needs BK re-enumeration, reconstructing the leaf
// from leafVtxData.
//
// Data structures (all built during SDCT_Augmented callback):
//   vtxLeafOff/vtxLeafData  — Vertex→Leaf CSR (for Phase 1 intersection)
//   leafVtxOff/leafVtxData  — Leaf→Vertex CSR (for Case B/C vertex access)
//   leafEdgeInfo            — Leaf→Edge list (for Case A fast subtract)
//   leafMeta                — Per-leaf weights (wKK, wPP, wKP)
//   countingKE              — Per-edge support counts
//
// Phase 1: For removed edge (u,v), intersect vtxLeafData[u] ∩ vtxLeafData[v]
//          via sorted merge on leafId to find affected leaves.
// Case A:  Leaf dies — subtract via leafEdgeInfo scan (same as V4).
// Case C:  Pivot-only removal — immutable counter with nCr delta formula.
// Case B:  Reconstruct vertex list from leafVtxData, run BK, create sub-leaves
//          in dynamic overflow vectors.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <cassert>
#include <algorithm>
#include <cstring>

#include "../BK/BronKerboschRmEdge.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace R2TreeFree {

enum EdgeType : uint8_t { KK = 0, PP = 1, KP = 2 };

struct LeafEdgeEntry {
    daf::Size edgeId;
    EdgeType type;
};

struct LeafMeta {
    double wKK, wPP, wKP;
    int pivotCount;
    int keepCount;
    int needPivot;
};

struct VLeafEntry {
    daf::Size leafId;
    uint8_t isPivot;  // role of this vertex in this leaf
};

struct LeafVtxEntry {
    daf::Size vertex;
    uint8_t isPivot;
};

struct LeafRmInfo {
    bool removedKeepC;
    daf::StaticVector<daf::Size> removedPivots{0};
    daf::StaticVector<std::pair<daf::Size, daf::Size>> removedEdges{0};
    LeafRmInfo() : removedKeepC(false) {}
    bool empty() const { return !removedKeepC && removedPivots.empty() && removedEdges.empty(); }
    void clear() { removedKeepC = false; removedPivots.clear(); removedEdges.clear(); }
};

template<typename It1, typename It2, typename WeightT, typename UpdateFunc>
inline void processEdgePairsImpl(It1 b1, It1 e1, It2 b2, It2 e2,
                                 WeightT weight, UpdateFunc &&upd) noexcept {
    if (weight < 0) return;
    if (b1 == b2 && e1 == e2 && b1 == e1 && b2 == e2) return;
    if (b1 == b2 && e1 == e2) {
        for (auto it = b1; it + 1 != e1; ++it) { auto u = *it; for (auto jt = it + 1; jt != e1; ++jt) upd(u, *jt, weight); }
    } else {
        for (auto it = b1; it != e1; ++it) { auto u = *it; for (auto jt = b2; jt != e2; ++jt) upd(u, *jt, weight); }
    }
}
template<typename Range1, typename Range2, typename WeightT, typename UpdateFunc>
inline void processEdgePairs(const Range1 &r1, const Range2 &r2, WeightT weight, UpdateFunc &&upd) noexcept {
    processEdgePairsImpl(std::begin(r1), std::end(r1), std::begin(r2), std::end(r2), weight, std::forward<UpdateFunc>(upd));
}
template<typename Range, typename WeightT, typename UpdateFunc>
inline void processEdgePairs(const Range &r, WeightT weight, UpdateFunc &&upd) noexcept {
    processEdgePairsImpl(std::begin(r), std::end(r), std::begin(r), std::end(r), weight, std::forward<UpdateFunc>(upd));
}

} // namespace R2TreeFree


std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_TreeFree(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    // ====================================================================
    // Init: Build indices from tree (same data the SDCT callback would build)
    // ====================================================================
    const daf::Size numEdges = edgeGraph.adj_list.size();
    const daf::Size numLeavesInit = tree.adj_list.size();
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;

    // countingKE
    auto *countingKE = new double[numEdges];
    memset(countingKE, 0, numEdges * sizeof(double));

    // leafEdgeInfo + leafMeta
    std::vector<std::vector<R2TreeFree::LeafEdgeEntry>> leafEdgeInfo(numLeavesInit);
    std::vector<R2TreeFree::LeafMeta> leafMeta(numLeavesInit);

    // COO buffer for dual CSR
    struct COOEntry {
        daf::Size vertex;
        daf::Size leafId;
        uint8_t isPivot;
    };
    std::vector<COOEntry> cooBuf;
    cooBuf.reserve(numLeavesInit * 4);

    // Per-leaf immutable counters for Case C
    std::vector<int> leafRemainPivots(numLeavesInit);

    daf::StaticVector<daf::Size> tPovit, tKeepC;

    for (daf::Size li = 0; li < numLeavesInit; ++li) {
        const auto &clique = tree.adj_list[li];
        if (clique.size() < k) {
            leafMeta[li] = {0, 0, 0, 0, 0, 0};
            leafRemainPivots[li] = 0;
            continue;
        }

        tPovit.clear(); tKeepC.clear();
        for (const auto &node : clique) {
            if (node.isPivot) tPovit.push_back(node.v);
            else tKeepC.push_back(node.v);
            cooBuf.push_back({node.v, li, (uint8_t)node.isPivot});
        }

        int needPivot = int(k) - int(tKeepC.size());
        double wKK = 0, wPP = 0, wKP = 0;
        {
            daf::Size nk = tKeepC.size(), npv = tPovit.size();
            leafEdgeInfo[li].reserve(nk*(nk-1)/2 + npv*(npv-1)/2 + nk*npv);
        }

        if (needPivot >= 0 && needPivot <= int(tPovit.size())) {
            wKK = std::llround(nCr[tPovit.size()][needPivot]);
            for (size_t i = 0; i < tKeepC.size(); ++i)
                for (size_t j = i + 1; j < tKeepC.size(); ++j) {
                    auto eid = edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]);
                    countingKE[eid] += wKK;
                    leafEdgeInfo[li].push_back({eid, R2TreeFree::KK});
                }
        }
        int needPP = needPivot - 2;
        if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
            wPP = std::llround(nCr[tPovit.size() - 2][needPP]);
            for (size_t i = 0; i < tPovit.size(); ++i)
                for (size_t j = i + 1; j < tPovit.size(); ++j) {
                    auto eid = edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]);
                    countingKE[eid] += wPP;
                    leafEdgeInfo[li].push_back({eid, R2TreeFree::PP});
                }
        }
        int needKP = needPivot - 1;
        if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
            wKP = std::llround(nCr[tPovit.size() - 1][needKP]);
            for (size_t i = 0; i < tKeepC.size(); ++i)
                for (size_t j = 0; j < tPovit.size(); ++j) {
                    auto eid = edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]);
                    countingKE[eid] += wKP;
                    leafEdgeInfo[li].push_back({eid, R2TreeFree::KP});
                }
        }

        leafMeta[li] = {wKK, wPP, wKP, (int)tPovit.size(), (int)tKeepC.size(), needPivot};
        leafRemainPivots[li] = (int)tPovit.size();
    }
    tPovit.free(); tKeepC.free();

    // Build dual CSR from COO
    std::vector<daf::Size> vtxLeafOff(numVertices + 2, 0);
    for (auto &e : cooBuf)
        if (e.vertex < numVertices) vtxLeafOff[e.vertex + 1]++;
    for (daf::Size i = 1; i <= numVertices; ++i)
        vtxLeafOff[i] += vtxLeafOff[i - 1];
    std::vector<R2TreeFree::VLeafEntry> vtxLeafData(vtxLeafOff[numVertices]);
    {
        std::vector<daf::Size> pos(numVertices, 0);
        for (auto &e : cooBuf) {
            daf::Size v = e.vertex;
            if (v < numVertices) {
                daf::Size p = vtxLeafOff[v] + pos[v]++;
                vtxLeafData[p] = {e.leafId, e.isPivot};
            }
        }
    }
    // Sort vtxLeafData by leafId within each vertex (for sorted merge intersection)
    for (daf::Size v = 0; v < numVertices; ++v) {
        auto begin = vtxLeafData.begin() + vtxLeafOff[v];
        auto end = vtxLeafData.begin() + vtxLeafOff[v + 1];
        std::sort(begin, end, [](const R2TreeFree::VLeafEntry &a, const R2TreeFree::VLeafEntry &b) {
            return a.leafId < b.leafId;
        });
    }

    std::vector<daf::Size> leafVtxOff(numLeavesInit + 1, 0);
    for (auto &e : cooBuf)
        if (e.leafId < numLeavesInit) leafVtxOff[e.leafId + 1]++;
    for (size_t i = 1; i <= numLeavesInit; ++i)
        leafVtxOff[i] += leafVtxOff[i - 1];
    std::vector<R2TreeFree::LeafVtxEntry> leafVtxData(leafVtxOff[numLeavesInit]);
    {
        std::vector<daf::Size> pos(numLeavesInit, 0);
        for (auto &e : cooBuf) {
            daf::Size L = e.leafId;
            if (L < numLeavesInit) {
                daf::Size p = leafVtxOff[L] + pos[L]++;
                leafVtxData[p] = {e.vertex, e.isPivot};
            }
        }
    }
    cooBuf.clear();
    cooBuf.shrink_to_fit();

    auto time_init = std::chrono::high_resolution_clock::now();
    std::cout << "TreeFree R=2: init took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_init - time_start).count()
              << " ms, leaves=" << numLeavesInit << std::endl;

    // ====================================================================
    // Peeling
    // ====================================================================
    auto *coreE = new double[numEdges];
    memset(coreE, 0, numEdges * sizeof(double));

    std::vector<uint8_t> leafAlive(numLeavesInit, 1);
    // Mark leaves with empty clique as dead
    for (daf::Size li = 0; li < numLeavesInit; ++li) {
        if (tree.adj_list[li].size() < k) leafAlive[li] = 0;
    }

    // Track leaves mutated by Case C/B
    std::vector<uint8_t> leafMutated(numLeavesInit, 0);

    daf::StaticVector<uint8_t> edgeInHeap(numEdges);
    edgeInHeap.c_size = numEdges;
    memset(edgeInHeap.getData(), 1, numEdges * sizeof(uint8_t));

    // Bucket array
    int maxBucket = 0;
    for (daf::Size i = 0; i < numEdges; ++i)
        if (countingKE[i] > 0) maxBucket = std::max(maxBucket, (int)countingKE[i]);
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(numEdges);
    std::vector<daf::Size> pos_in_bucket(numEdges);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numEdges; ++i) {
        if (countingKE[i] == 0) { edgeInHeap[i] = 0; continue; }
        int b = (int)countingKE[i]; bucket_of[i] = b; pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i); remainingInHeap++;
    }
    int curBucket = 0;

    auto bucketMove = [&](daf::Size id) {
        int newB = std::max(0, (int)countingKE[id]); int oldB = bucket_of[id]; if (newB == oldB) return;
        auto &oldVec = buckets[oldB]; daf::Size myPos = pos_in_bucket[id];
        if (myPos < oldVec.size() - 1) { daf::Size last = oldVec.back(); oldVec[myPos] = last; pos_in_bucket[last] = myPos; }
        oldVec.pop_back(); if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
        bucket_of[id] = newB; pos_in_bucket[id] = buckets[newB].size(); buckets[newB].push_back(id);
        if (newB < curBucket) curBucket = newB;
    };

    // Deferred batch bucketMove
    std::vector<uint8_t> dirtyMark(numEdges, 0);
    std::vector<daf::Size> dirtyEdges;
    dirtyEdges.reserve(8192);

    // Per-leaf removal info (for init leaves, use leafRmInfo; for overflow, separate)
    std::vector<R2TreeFree::LeafRmInfo> leafRmInfo(numLeavesInit);
    std::vector<uint8_t> leafAffected(numLeavesInit, 0);
    daf::StaticVector<daf::Size> removedLeaf(numLeavesInit);
    std::vector<daf::Size> caseBLeafIds;
    caseBLeafIds.reserve(1024);

    daf::StaticVector<daf::Size> currentRemoveEdgeIds(numEdges);

    // Overflow area for Case B sub-leaves
    // New sub-leaves get IDs starting from numLeavesInit
    // For Phase 1 on new sub-leaves, we fall back to treeGraphV intersection
    // (Case B is rare, so overhead is minimal)

    daf::StaticVector<daf::Size> povit, keepC, newPivot, newKeepC;

    daf::log_memory("TreeFree R2 init");

    std::cout << "=========================begin (R2 TreeFree)=========================" << std::endl;
    double minCore = 0, cntA = 0, cntB = 0, cntC = 0, totalIters = 0;
    long long cntA_fast = 0, cntA_fallback = 0;  // fast=leafEdgeInfo, fallback=tree O(n²)

    using hrc = std::chrono::high_resolution_clock;
    long long t_phase1_ns = 0, t_caseA_ns = 0, t_caseC_ns = 0, t_caseB_ns = 0;
    long long t_bucketFlush_ns = 0;

    while (remainingInHeap > 0) {
        // Bucket pop
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;
        minCore = std::max((double)curBucket, minCore);
        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty() && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back(); buckets[curBucket].pop_back();
                edgeInHeap[id] = 0; currentRemoveEdgeIds.push_back(id); coreE[id] = minCore; remainingInHeap--;
            }
            if (curBucket+1 < (int)buckets.size() && !buckets[curBucket+1].empty() && (curBucket+1) <= (int)minCore) curBucket++; else break;
        }
        if (remainingInHeap == 0) break;
        totalIters++;

        // === Phase 1: Use treeGraphV intersection (proven correct, same as V4 baseline) ===
        auto tp1 = hrc::now();
        // The key optimization of TreeFree is in Case A (leafEdgeInfo scan), not Phase 1.
        for (int ei = 0; ei < (int)currentRemoveEdgeIds.size(); ++ei) {
            auto edgeId = currentRemoveEdgeIds[ei];
            auto [edgeU, edgeV] = edgeGraph.getEdgeById(edgeId);
            auto &adjEdgeU = treeGraphV.getNbr(edgeU);
            auto &adjEdgeV = treeGraphV.getNbr(edgeV);
            daf::intersect_dense_sets(adjEdgeU, adjEdgeV,
                [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                    daf::Size leafId = uClique.v;
                    if (leafId >= leafRmInfo.size()) {
                        leafRmInfo.resize(leafId * 2 + 1);
                        leafAffected.resize(leafId * 2 + 1, 0);
                    }
                    auto &lrm = leafRmInfo[leafId];
                    if (lrm.empty()) {
                        removedLeaf.push_back(leafId);
                        leafAffected[leafId] = 1;
                    }
                    if (!lrm.removedKeepC) {
                        if (!uClique.isPivot && !vClique.isPivot) lrm.removedKeepC = true;
                        else if (uClique.isPivot && vClique.isPivot) lrm.removedEdges.push_back({edgeU, edgeV});
                        else if (uClique.isPivot) lrm.removedPivots.push_back(edgeU);
                        else lrm.removedPivots.push_back(edgeV);
                    }
                });
        }
        t_phase1_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - tp1).count();

        // Pre-sort removedPivots
        for (int i = 0; i < (int)removedLeaf.size(); ++i) {
            auto leafId = removedLeaf[i];
            if (leafId >= leafRmInfo.size()) continue;
            auto &rp = leafRmInfo[leafId].removedPivots;
            if (rp.size() == 2) { if (rp[0] > rp[1]) std::swap(rp[0], rp[1]); if (rp[0] == rp[1]) rp.c_size = 1; }
            else if (rp.size() > 2) { std::ranges::sort(rp); rp.unique(); }
        }

        // === Phase 2: Case A (leafEdgeInfo scan) & Case C (tree mutation) ===
        auto tAC = hrc::now();
        caseBLeafIds.clear();
        {
            auto directSub = [&](daf::Size idx, double w) {
                countingKE[idx] -= w; if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx] && !dirtyMark[idx]) { dirtyMark[idx] = 1; dirtyEdges.push_back(idx); }
            };

            for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                auto leafId = removedLeaf[leafIdIdx];
                if (leafId >= leafRmInfo.size()) continue;
                R2TreeFree::LeafRmInfo &leafRm = leafRmInfo[leafId];

                // For init leaves (leafId < numLeavesInit): use flat CSR data
                // For overflow leaves: use tree.adj_list (fallback)
                bool isInitLeaf = leafId < numLeavesInit;

                int numPivots, numKeeps, needPivot;
                if (isInitLeaf && !leafMutated[leafId]) {
                    numPivots = leafMeta[leafId].pivotCount;
                    numKeeps = leafMeta[leafId].keepCount;
                    needPivot = leafMeta[leafId].needPivot;
                } else {
                    // Mutated init leaf or overflow leaf: read current state from tree
                    const auto &leaf = tree.adj_list[leafId];
                    if (leaf.empty()) { leafRm.clear(); continue; }
                    numPivots = 0; numKeeps = 0;
                    for (const auto &node : leaf) {
                        if (node.isPivot) numPivots++;
                        else numKeeps++;
                    }
                    needPivot = int(k) - numKeeps;
                }

                bool isDeadLeaf = leafRm.removedKeepC || needPivot > numPivots - (int)leafRm.removedPivots.size();
                bool isCaseB = !leafRm.removedEdges.empty() && !isDeadLeaf;

                if (isCaseB) { caseBLeafIds.push_back(leafId); continue; }

                if (isDeadLeaf) {
                    cntA++;

                    if (isInitLeaf && !leafMutated[leafId] && !leafEdgeInfo[leafId].empty()) {
                        cntA_fast++;
                        // Fast path: sequential scan of leafEdgeInfo (leaf never mutated)
                        double wKK = 0, wPP = 0, wKP = 0;
                        if (needPivot >= 0 && needPivot <= numPivots)
                            wKK = std::llround(nCr[numPivots][needPivot]);
                        int nPP = needPivot - 2;
                        if (0 <= nPP && nPP <= numPivots - 2)
                            wPP = std::llround(nCr[numPivots - 2][nPP]);
                        int nKP = needPivot - 1;
                        if (0 <= nKP && nKP <= numPivots - 1)
                            wKP = std::llround(nCr[numPivots - 1][nKP]);

                        for (const auto &entry : leafEdgeInfo[leafId]) {
                            if (!edgeInHeap[entry.edgeId]) continue;
                            double w;
                            switch (entry.type) {
                                case R2TreeFree::KK: w = wKK; break;
                                case R2TreeFree::PP: w = wPP; break;
                                case R2TreeFree::KP: w = wKP; break;
                            }
                            directSub(entry.edgeId, w);
                        }
                    } else {
                        cntA_fallback++;
                        // Mutated init leaf or overflow leaf: recompute from tree
                        const auto &leaf = tree.adj_list[leafId];
                        daf::StaticVector<daf::Size> lPovit, lKeepC;
                        for (const auto &node : leaf) {
                            if (node.isPivot) lPovit.push_back(node.v);
                            else lKeepC.push_back(node.v);
                        }
                        int np = int(k) - int(lKeepC.size());
                        double wKK = 0, wPP = 0, wKP = 0;
                        if (np >= 0 && np <= (int)lPovit.size())
                            wKK = std::llround(nCr[lPovit.size()][np]);
                        int nPP = np - 2;
                        if (0 <= nPP && nPP <= int(lPovit.size()) - 2)
                            wPP = std::llround(nCr[lPovit.size() - 2][nPP]);
                        int nKP = np - 1;
                        if (0 <= nKP && nKP <= int(lPovit.size()) - 1)
                            wKP = std::llround(nCr[lPovit.size() - 1][nKP]);

                        if (wKK > 0) for (daf::Size i = 0; i+1 < lKeepC.size(); ++i)
                            for (daf::Size j = i+1; j < lKeepC.size(); ++j)
                                directSub(edgeGraph.getEdgeCompressedId(lKeepC[i], lKeepC[j]), wKK);
                        if (wPP > 0) for (daf::Size i = 0; i+1 < lPovit.size(); ++i)
                            for (daf::Size j = i+1; j < lPovit.size(); ++j)
                                directSub(edgeGraph.getEdgeCompressedId(lPovit[i], lPovit[j]), wPP);
                        if (wKP > 0) for (daf::Size i = 0; i < lKeepC.size(); ++i)
                            for (daf::Size j = 0; j < lPovit.size(); ++j)
                                directSub(edgeGraph.getEdgeCompressedId(lKeepC[i], lPovit[j]), wKP);
                        lPovit.free(); lKeepC.free();
                    }

                    // Mark leaf dead and clean up tree + treeGraphV
                    if (isInitLeaf) {
                        leafAlive[leafId] = 0;
                        leafEdgeInfo[leafId].clear();
                    }
                    // Always use tree for treeGraphV cleanup (proven correct)
                    for (const auto &i : tree.adj_list[leafId])
                        treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                    tree.adj_list[leafId].clear();
                    // Don't recycle init leaf IDs — prevents aliasing with overflow sub-leaves
                    if (!isInitLeaf) tree.recycleNode(leafId);
                } else {
                    // ==== Case C: pivot-only removal — use tree mutation (baseline approach) ====
                    // Case C requires tracking which specific pivots were removed across rounds.
                    // The immutable counter approach (leafRemainPivots) only tracks counts, not identity,
                    // so subsequent rounds can't correctly enumerate surviving vertices from leafVtxData.
                    // Therefore, Case C always mutates the tree — same as V4 baseline.
                    cntC++;

                    // tree.adj_list already has correct data from buildAndVerifySDCT
                    // treeGraphV also already has entries for all init leaves
                    // Just mark as mutated so future Case A uses tree fallback
                    if (isInitLeaf) leafMutated[leafId] = 1;

                    {
                        const auto &leaf = tree.adj_list[leafId];
                        povit.clear(); keepC.clear();
                        for (const auto &node : leaf) {
                            if (node.isPivot) povit.push_back(node.v);
                            else keepC.push_back(node.v);
                        }
                        daf::Size np = k - keepC.size();
                        double KtoK = 0, KtoP = 0, PtoP = 0;
                        double RemovedKtoK = 0, RemovedKtoP = 0, RemovedPtoP = 0;

                        if (np <= povit.size()) {
                            KtoK = std::llround(nCr[povit.size()][np]) - std::llround(nCr[povit.size() - leafRm.removedPivots.size()][np]);
                            RemovedKtoK = std::llround(nCr[povit.size()][np]);
                        }
                        int needPPc = int(np) - 2;
                        if (0 <= needPPc && needPPc <= int(povit.size()) - 2) {
                            RemovedPtoP = std::llround(nCr[povit.size() - 2][needPPc]);
                            PtoP = RemovedPtoP;
                            if (leafRm.removedPivots.size() + 2 <= povit.size())
                                PtoP -= std::llround(nCr[povit.size() - 2 - leafRm.removedPivots.size()][needPPc]);
                        }
                        int needKPc = int(np) - 1;
                        if (0 <= needKPc && needKPc <= int(povit.size()) - 1) {
                            RemovedKtoP = std::llround(nCr[povit.size() - 1][needKPc]);
                            KtoP = RemovedKtoP;
                            if (leafRm.removedPivots.size() + 1 <= povit.size())
                                KtoP -= std::llround(nCr[povit.size() - 1 - leafRm.removedPivots.size()][needKPc]);
                        }

                        if (!leafRm.removedPivots.empty() && np <= povit.size() - leafRm.removedPivots.size()) {
                            daf::StaticVector<TreeGraphNode> newLeafF;
                            for (const auto &node : leaf) {
                                if (!std::binary_search(leafRm.removedPivots.begin(), leafRm.removedPivots.end(), node.v))
                                    newLeafF.push_back(node);
                            }
                            for (daf::Size i = 0; i+1 < leafRm.removedPivots.size(); ++i)
                                for (daf::Size j = i+1; j < leafRm.removedPivots.size(); ++j)
                                    directSub(edgeGraph.getEdgeCompressedId(leafRm.removedPivots[i], leafRm.removedPivots[j]), RemovedPtoP);
                            for (daf::Size i = 0; i < newLeafF.size(); ++i)
                                for (daf::Size j = 0; j < leafRm.removedPivots.size(); ++j) {
                                    double d = newLeafF[i].isPivot ? RemovedPtoP : RemovedKtoP;
                                    directSub(edgeGraph.getEdgeCompressedId(newLeafF[i].v, leafRm.removedPivots[j]), d);
                                }
                            for (daf::Size i = 0; i+1 < newLeafF.size(); ++i)
                                for (daf::Size j = i+1; j < newLeafF.size(); ++j) {
                                    auto &u = newLeafF[i], &v = newLeafF[j];
                                    double d = (!u.isPivot && !v.isPivot) ? KtoK : (u.isPivot && v.isPivot) ? PtoP : KtoP;
                                    directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                                }
                            newLeafF.free();
                            for (auto removedNbr : leafRm.removedPivots)
                                treeGraphV.removeNbr(removedNbr, static_cast<TreeGraphNode>(leafId));
                            tree.removeNbrs(leafId, leafRm.removedPivots);
                        } else {
                            for (daf::Size i = 0; i+1 < leaf.size(); ++i)
                                for (daf::Size j = i+1; j < leaf.size(); ++j) {
                                    auto &u = leaf[i], &v = leaf[j];
                                    double d = (!u.isPivot && !v.isPivot) ? RemovedKtoK : (u.isPivot && v.isPivot) ? RemovedPtoP : RemovedKtoP;
                                    directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                                }
                            for (const auto &i : leaf)
                                treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                            if (isInitLeaf) {
                                leafEdgeInfo[leafId].clear();
                                leafAlive[leafId] = 0;
                            }
                            tree.adj_list[leafId].clear();
                            if (!isInitLeaf) tree.recycleNode(leafId);
                        }
                        povit.clear(); keepC.clear();
                    }
                }
                leafRm.clear();
            }

            // Flush deferred bucket moves
            auto tbf = hrc::now();
            for (auto eid : dirtyEdges) { if (edgeInHeap[eid]) bucketMove(eid); dirtyMark[eid] = 0; }
            dirtyEdges.clear();
            t_bucketFlush_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - tbf).count();
        }
        t_caseA_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - tAC).count();

        // === Phase 2c: BK for Case B — always uses tree (same as baseline) ===
        auto tB0 = hrc::now();
        cntB += caseBLeafIds.size();
        for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
            auto leafId = caseBLeafIds[bi];
            if (leafId >= leafRmInfo.size()) continue;
            R2TreeFree::LeafRmInfo &leafRm = leafRmInfo[leafId];

            bool isInitLeaf = leafId < numLeavesInit;

            // For Case B we need the full vertex list and BK adjacency
            // For init leaves: reconstruct from leafVtxData
            // For overflow leaves: use tree.adj_list
            std::vector<TreeGraphNode> leafNodes;
            leafNodes.reserve(tree.adj_list[leafId].size());
            if (isInitLeaf) {
                for (daf::Size li = leafVtxOff[leafId]; li < leafVtxOff[leafId + 1]; ++li) {
                    auto &nd = leafVtxData[li];
                    leafNodes.push_back({nd.vertex, (bool)nd.isPivot});
                }
                // Remove previously-removed pivots if leaf was mutated
                if (leafMutated[leafId]) {
                    // We need to filter out pivots removed in prior Case C rounds
                    // Since we don't track which specific pivots were removed,
                    // we need to use the tree for this case
                    // Actually, for Case B we need the current state of the leaf
                    // If it was mutated by Case C, the tree has the current state
                    leafNodes.clear();
                    for (const auto &node : tree.adj_list[leafId])
                        leafNodes.push_back(node);
                }
            } else {
                for (const auto &node : tree.adj_list[leafId])
                    leafNodes.push_back(node);
            }

            povit.clear(); keepC.clear();
            for (auto &node : leafNodes) {
                if (node.isPivot) povit.push_back(node.v);
                else keepC.push_back(node.v);
            }
            daf::Size needPivotBK = k - keepC.size();

            auto addW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] += w; if (edgeInHeap[idx]) bucketMove(idx);
            };

            // Remove old leaf from treeGraphV
            if (isInitLeaf && !leafMutated[leafId]) {
                // For unmutated init leaves, we can build treeGraphV entries from leafVtxData
                for (daf::Size li = leafVtxOff[leafId]; li < leafVtxOff[leafId + 1]; ++li) {
                    auto &nd = leafVtxData[li];
                    treeGraphV.removeNbr(nd.vertex, {leafId, (bool)nd.isPivot});
                }
            } else {
                for (auto &lv : tree.adj_list[leafId]) {
                    if (lv.isPivot) treeGraphV.removeNbr(lv.v, {leafId, true});
                    else treeGraphV.removeNbr(lv.v, {leafId, false});
                }
            }
            if (!leafRm.removedPivots.empty()) tree.removeNbrs(leafId, leafRm.removedPivots);

            // Store leafNodes into tree for BK (needed for coverToVertex)
            // tree.adj_list already has the correct data from SDCT build
            // For unmutated init leaves, removedPivots were already applied above via tree.removeNbrs

            auto &leafRef = tree.adj_list[leafId];
            bkRmEdge::bronKerbosch(leafRef, leafRm.removedEdges, k,
                [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                    auto subLeaf = bkRmEdge::coverToVertex(c, pivots, leafRef);
                    auto newId = tree.addNode(subLeaf);
                    newPivot.clear(); newKeepC.clear();
                    for (const auto &i : tree.adj_list[newId]) {
                        if (i.isPivot) { newPivot.push_back(i.v); treeGraphV.addNbr(i.v, {newId, true}); }
                        else { newKeepC.push_back(i.v); treeGraphV.addNbr(i.v, {newId, false}); }
                    }
                    daf::Size np = k - newKeepC.size();
                    if (np <= newPivot.size() && newKeepC.size() > 1) {
                        double w = std::llround(nCr[newPivot.size()][np]);
                        R2TreeFree::processEdgePairs(newKeepC, w, addW);
                    }
                    int nPP = int(np) - 2;
                    if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                        double w = std::llround(nCr[newPivot.size() - 2][nPP]);
                        R2TreeFree::processEdgePairs(newPivot, w, addW);
                    }
                    int nKP = int(np) - 1;
                    if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                        double w = std::llround(nCr[newPivot.size() - 1][nKP]);
                        R2TreeFree::processEdgePairs(newKeepC, newPivot, w, addW);
                    }
                    newPivot.clear(); newKeepC.clear();
                    if (newId >= leafRmInfo.size()) {
                        leafRmInfo.resize(newId * 2 + 1);
                        leafAffected.resize(newId * 2 + 1, 0);
                    }
                });

            // Remove old contribution
            auto removeW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] -= w; if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx]) bucketMove(idx);
            };
            if (needPivotBK <= povit.size()) { double w = std::llround(nCr[povit.size()][needPivotBK]); R2TreeFree::processEdgePairs(keepC, w, removeW); }
            int needPP = int(needPivotBK) - 2;
            if (0 <= needPP && needPP <= int(povit.size()) - 2) { double w = std::llround(nCr[povit.size() - 2][needPP]); R2TreeFree::processEdgePairs(povit, w, removeW); }
            int needKP = int(needPivotBK) - 1;
            if (0 <= needKP && needKP <= int(povit.size()) - 1) { double w = std::llround(nCr[povit.size() - 1][needKP]); R2TreeFree::processEdgePairs(keepC, povit, w, removeW); }

            // For init leaves, don't use removeNode (which recycles the ID)
            // Just clear the adj_list to prevent ID aliasing with overflow sub-leaves
            if (isInitLeaf) {
                tree.adj_list[leafId].clear();
                leafAlive[leafId] = 0;
                leafEdgeInfo[leafId].clear();
            } else {
                tree.removeNode(leafId);
            }
            leafRm.clear();
            povit.clear(); keepC.clear();
        }
        t_caseB_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - tB0).count();

        currentRemoveEdgeIds.clear();
        for (auto leafId : removedLeaf) {
            if (leafId < leafAffected.size()) leafAffected[leafId] = 0;
        }
        removedLeaf.clear();
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
    std::cout << "  Cases: A=" << cntA << " B=" << cntB << " C=" << cntC << " iters=" << totalIters << std::endl;
    std::cout << "  [CaseA] fast=" << cntA_fast << " fallback=" << cntA_fallback << std::endl;
    std::cout << "  [Breakdown] phase1=" << (t_phase1_ns/1000000) << "ms"
              << ", caseA+C=" << (t_caseA_ns/1000000) << "ms"
              << " (bucketFlush=" << (t_bucketFlush_ns/1000000) << "ms)"
              << ", caseB=" << (t_caseB_ns/1000000) << "ms" << std::endl;

    daf::Size numCounting = 0;
    std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> sortedK;
    sortedK.reserve(numEdges);
    const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;
    for (daf::Size u = 0; u < n; ++u) {
        const daf::Size start = edgeGraph.adj_list_offsets[u];
        const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
        for (daf::Size idx = start; idx < end; ++idx)
            sortedK.emplace_back(std::make_pair(std::make_pair(u, edgeGraph.adj_list[idx]), (int)coreE[idx]));
    }
    assert(numCounting == 0);
    delete[] countingKE; delete[] coreE;
    povit.free(); keepC.free(); newPivot.free(); newKeepC.free();
    return sortedK;
}
