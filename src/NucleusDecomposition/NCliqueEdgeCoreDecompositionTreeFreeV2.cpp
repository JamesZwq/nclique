//
// NCliqueEdgeCoreDecompositionTreeFreeV2.cpp
//
// Optimized Tree-Free R=2 Edge Nucleus Decomposition.
//
// Improvements over TreeFree.cpp:
//   Opt 4: Remove unused dual CSR (vtxLeafOff/vtxLeafData/leafVtxOff/leafVtxData/COO)
//          — these were never used during peeling, saving ~30-40% init time
//   Opt 5: Edge→Leaf reverse index for Phase 1 — replaces treeGraphV hash-set
//          intersection with direct CSR scan for init leaves
//   Opt 6: Rebuild leafEdgeInfo after Case C — re-enables fast path for future Case A
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <cassert>
#include <algorithm>
#include <cstring>

#include "../BK/BronKerboschRmEdge.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace R2TreeFreeV2 {

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

// Opt 5: Edge→Leaf CSR entry
struct EdgeLeafEntry {
    daf::Size leafId;
    EdgeType edgeType;
    daf::Size pivotVertex;  // for KP: which vertex is the pivot; unused for KK/PP
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

} // namespace R2TreeFreeV2


std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_TreeFreeV2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    // ====================================================================
    // Init: Build indices from tree
    // ====================================================================
    const daf::Size numEdges = edgeGraph.adj_list.size();
    const daf::Size numLeavesInit = tree.adj_list.size();
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;

    // countingKE
    auto *countingKE = new double[numEdges];
    memset(countingKE, 0, numEdges * sizeof(double));

    // leafEdgeInfo + leafMeta
    std::vector<std::vector<R2TreeFreeV2::LeafEdgeEntry>> leafEdgeInfo(numLeavesInit);
    std::vector<R2TreeFreeV2::LeafMeta> leafMeta(numLeavesInit);

    // Opt 4: No COO buffer, no dual CSR — not needed during peeling

    // Opt 5: Edge→Leaf COO buffer (built during init, converted to CSR after)
    struct EdgeLeafCOO {
        daf::Size edgeId;
        daf::Size leafId;
        R2TreeFreeV2::EdgeType edgeType;
        daf::Size pivotVertex;
    };
    std::vector<EdgeLeafCOO> edgeLeafCOO;
    edgeLeafCOO.reserve(numLeavesInit * 4);

    daf::StaticVector<daf::Size> tPovit, tKeepC;

    for (daf::Size li = 0; li < numLeavesInit; ++li) {
        const auto &clique = tree.adj_list[li];
        if (clique.size() < k) {
            leafMeta[li] = {0, 0, 0, 0, 0, 0};
            continue;
        }

        tPovit.clear(); tKeepC.clear();
        for (const auto &node : clique) {
            if (node.isPivot) tPovit.push_back(node.v);
            else tKeepC.push_back(node.v);
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
                    leafEdgeInfo[li].push_back({eid, R2TreeFreeV2::KK});
                    edgeLeafCOO.push_back({eid, li, R2TreeFreeV2::KK, 0});
                }
        }
        int needPP = needPivot - 2;
        if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
            wPP = std::llround(nCr[tPovit.size() - 2][needPP]);
            for (size_t i = 0; i < tPovit.size(); ++i)
                for (size_t j = i + 1; j < tPovit.size(); ++j) {
                    auto eid = edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]);
                    countingKE[eid] += wPP;
                    leafEdgeInfo[li].push_back({eid, R2TreeFreeV2::PP});
                    edgeLeafCOO.push_back({eid, li, R2TreeFreeV2::PP, 0});
                }
        }
        int needKP = needPivot - 1;
        if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
            wKP = std::llround(nCr[tPovit.size() - 1][needKP]);
            for (size_t i = 0; i < tKeepC.size(); ++i)
                for (size_t j = 0; j < tPovit.size(); ++j) {
                    auto eid = edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]);
                    countingKE[eid] += wKP;
                    leafEdgeInfo[li].push_back({eid, R2TreeFreeV2::KP});
                    // For KP: tPovit[j] is the pivot vertex
                    edgeLeafCOO.push_back({eid, li, R2TreeFreeV2::KP, tPovit[j]});
                }
        }

        leafMeta[li] = {wKK, wPP, wKP, (int)tPovit.size(), (int)tKeepC.size(), needPivot};
    }
    tPovit.free(); tKeepC.free();

    // Opt 5: Build Edge→Leaf CSR from COO
    std::vector<daf::Size> edgeLeafOff(numEdges + 1, 0);
    for (auto &e : edgeLeafCOO)
        edgeLeafOff[e.edgeId + 1]++;
    for (daf::Size i = 1; i <= numEdges; ++i)
        edgeLeafOff[i] += edgeLeafOff[i - 1];
    std::vector<R2TreeFreeV2::EdgeLeafEntry> edgeLeafData(edgeLeafOff[numEdges]);
    {
        std::vector<daf::Size> pos(numEdges, 0);
        for (auto &e : edgeLeafCOO) {
            daf::Size p = edgeLeafOff[e.edgeId] + pos[e.edgeId]++;
            edgeLeafData[p] = {e.leafId, e.edgeType, e.pivotVertex};
        }
    }
    edgeLeafCOO.clear();
    edgeLeafCOO.shrink_to_fit();

    auto time_init = std::chrono::high_resolution_clock::now();
    std::cout << "TreeFreeV2 R=2: init took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(time_init - time_start).count()
              << " ms, leaves=" << numLeavesInit
              << ", edgeLeafEntries=" << edgeLeafData.size() << std::endl;

    // ====================================================================
    // Peeling
    // ====================================================================
    auto *coreE = new double[numEdges];
    memset(coreE, 0, numEdges * sizeof(double));

    std::vector<uint8_t> leafAlive(numLeavesInit, 1);
    for (daf::Size li = 0; li < numLeavesInit; ++li) {
        if (tree.adj_list[li].size() < k) leafAlive[li] = 0;
    }

    std::vector<uint8_t> leafMutated(numLeavesInit, 0);
    // Opt 6: tracks whether leafEdgeInfo was rebuilt after Case C
    // leafMutated stays 1 forever (Edge→Leaf CSR is stale), but this flag
    // re-enables the Case A fast path via the rebuilt leafEdgeInfo.
    std::vector<uint8_t> leafEdgeInfoRebuilt(numLeavesInit, 0);

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

    std::vector<R2TreeFreeV2::LeafRmInfo> leafRmInfo(numLeavesInit);
    std::vector<uint8_t> leafAffected(numLeavesInit, 0);
    daf::StaticVector<daf::Size> removedLeaf(numLeavesInit);
    std::vector<daf::Size> caseBLeafIds;
    caseBLeafIds.reserve(1024);

    daf::StaticVector<daf::Size> currentRemoveEdgeIds(numEdges);

    daf::StaticVector<daf::Size> povit, keepC, newPivot, newKeepC;

    daf::log_memory("TreeFreeV2 R2 init");

    std::cout << "=========================begin (R2 TreeFreeV2)=========================" << std::endl;
    double minCore = 0, cntA = 0, cntB = 0, cntC = 0, totalIters = 0;
    long long cntA_fast = 0, cntA_fallback = 0;
    long long cntC_rebuild = 0;  // Opt 6: how many Case C rebuilds

    using hrc = std::chrono::high_resolution_clock;
    long long t_phase1_ns = 0, t_caseA_ns = 0, t_caseC_ns = 0, t_caseB_ns = 0;
    long long t_bucketFlush_ns = 0;
    long long t_phase1_edgeLeaf_ns = 0, t_phase1_overflow_ns = 0;

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

        // === Phase 1: Opt 5 — Edge→Leaf CSR scan for init leaves, treeGraphV fallback for overflow ===
        auto tp1 = hrc::now();
        for (int ei = 0; ei < (int)currentRemoveEdgeIds.size(); ++ei) {
            auto edgeId = currentRemoveEdgeIds[ei];
            auto [edgeU, edgeV] = edgeGraph.getEdgeById(edgeId);

            // Opt 5: Scan Edge→Leaf CSR for unmutated init leaves only
            // Mutated leaves have stale CSR entries (Case C changed edge types)
            // so they must fall through to the treeGraphV intersection below.
            auto tp1a = hrc::now();
            for (daf::Size idx = edgeLeafOff[edgeId]; idx < edgeLeafOff[edgeId + 1]; ++idx) {
                auto &entry = edgeLeafData[idx];
                daf::Size leafId = entry.leafId;
                if (leafId == (daf::Size)-1) continue;
                if (leafId >= numLeavesInit) continue;  // overflow → treeGraphV below
                if (!leafAlive[leafId]) continue;
                if (leafMutated[leafId]) continue;  // stale CSR → treeGraphV below

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
                    switch (entry.edgeType) {
                        case R2TreeFreeV2::KK:
                            lrm.removedKeepC = true;
                            break;
                        case R2TreeFreeV2::PP:
                            lrm.removedEdges.push_back({edgeU, edgeV});
                            break;
                        case R2TreeFreeV2::KP:
                            lrm.removedPivots.push_back(entry.pivotVertex);
                            break;
                    }
                }
            }
            t_phase1_edgeLeaf_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - tp1a).count();

            // Overflow leaves (from Case B): treeGraphV intersection
            auto tp1b = hrc::now();
            auto &adjEdgeU = treeGraphV.getNbr(edgeU);
            auto &adjEdgeV = treeGraphV.getNbr(edgeV);
            daf::intersect_dense_sets(adjEdgeU, adjEdgeV,
                [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                    daf::Size leafId = uClique.v;
                    // Skip unmutated init leaves (already handled by Edge→Leaf CSR)
                    // Mutated init leaves must be handled here (stale CSR entries)
                    if (leafId < numLeavesInit && !leafMutated[leafId]) return;
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
            t_phase1_overflow_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(hrc::now() - tp1b).count();
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

        // === Phase 2: Case A & Case C ===
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
                R2TreeFreeV2::LeafRmInfo &leafRm = leafRmInfo[leafId];

                bool isInitLeaf = leafId < numLeavesInit;

                int numPivots, numKeeps, needPivot;
                if (isInitLeaf && (!leafMutated[leafId] || leafEdgeInfoRebuilt[leafId])) {
                    numPivots = leafMeta[leafId].pivotCount;
                    numKeeps = leafMeta[leafId].keepCount;
                    needPivot = leafMeta[leafId].needPivot;
                } else {
                    const auto& leaf = tree.adj_list[leafId];
                    if (leaf.empty()) { leafRm.clear(); continue; }
                    numPivots = 0; numKeeps = 0;
                    for (const auto& node : leaf) {
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

                    if (isInitLeaf && (!leafMutated[leafId] || leafEdgeInfoRebuilt[leafId]) && !leafEdgeInfo[leafId].empty()) {
                        cntA_fast++;
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
                                case R2TreeFreeV2::KK: w = wKK; break;
                                case R2TreeFreeV2::PP: w = wPP; break;
                                case R2TreeFreeV2::KP: w = wKP; break;
                            }
                            directSub(entry.edgeId, w);
                        }
                    } else {
                        cntA_fallback++;
                        const auto& leaf = tree.adj_list[leafId];
                        daf::StaticVector<daf::Size> lPovit, lKeepC;
                        for (const auto& node : leaf) {
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

                    // Mark leaf dead + invalidate Edge→Leaf CSR entries
                    if (isInitLeaf) {
                        leafAlive[leafId] = 0;
                        // Opt 5: Invalidate Edge→Leaf entries for this leaf
                        // (leafEdgeInfo gives us the edge IDs, so we can find them in the CSR)
                        // No need — the leafAlive check in Phase 1 will skip them
                        leafEdgeInfo[leafId].clear();
                    }
                    for (const auto& i : tree.adj_list[leafId])
                        treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                    tree.adj_list[leafId].clear();
                    if (!isInitLeaf) tree.recycleNode(leafId);
                } else {
                    // ==== Case C: pivot-only removal ====
                    cntC++;

                    if (isInitLeaf) leafMutated[leafId] = 1;

                    {
                        const auto& leaf = tree.adj_list[leafId];
                        povit.clear(); keepC.clear();
                        for (const auto& node : leaf) {
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
                            for (const auto& node : leaf) {
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

                            // Opt 6: Rebuild leafEdgeInfo after Case C so future Case A uses fast path
                            if (isInitLeaf) {
                                cntC_rebuild++;
                                leafEdgeInfo[leafId].clear();
                                auto &newLeaf = tree.adj_list[leafId];
                                daf::StaticVector<daf::Size> rPovit, rKeepC;
                                for (const auto& node : newLeaf) {
                                    if (node.isPivot) rPovit.push_back(node.v);
                                    else rKeepC.push_back(node.v);
                                }
                                int rNp = int(k) - int(rKeepC.size());
                                double rWKK = 0, rWPP = 0, rWKP = 0;
                                {
                                    daf::Size nk = rKeepC.size(), npv = rPovit.size();
                                    leafEdgeInfo[leafId].reserve(nk*(nk-1)/2 + npv*(npv-1)/2 + nk*npv);
                                }

                                if (rNp >= 0 && rNp <= (int)rPovit.size()) {
                                    rWKK = std::llround(nCr[rPovit.size()][rNp]);
                                    for (size_t ri = 0; ri < rKeepC.size(); ++ri)
                                        for (size_t rj = ri + 1; rj < rKeepC.size(); ++rj)
                                            leafEdgeInfo[leafId].push_back({edgeGraph.getEdgeCompressedId(rKeepC[ri], rKeepC[rj]), R2TreeFreeV2::KK});
                                }
                                int rNPP = rNp - 2;
                                if (0 <= rNPP && rNPP <= int(rPovit.size()) - 2) {
                                    rWPP = std::llround(nCr[rPovit.size() - 2][rNPP]);
                                    for (size_t ri = 0; ri < rPovit.size(); ++ri)
                                        for (size_t rj = ri + 1; rj < rPovit.size(); ++rj)
                                            leafEdgeInfo[leafId].push_back({edgeGraph.getEdgeCompressedId(rPovit[ri], rPovit[rj]), R2TreeFreeV2::PP});
                                }
                                int rNKP = rNp - 1;
                                if (0 <= rNKP && rNKP <= int(rPovit.size()) - 1) {
                                    rWKP = std::llround(nCr[rPovit.size() - 1][rNKP]);
                                    for (size_t ri = 0; ri < rKeepC.size(); ++ri)
                                        for (size_t rj = 0; rj < rPovit.size(); ++rj)
                                            leafEdgeInfo[leafId].push_back({edgeGraph.getEdgeCompressedId(rKeepC[ri], rPovit[rj]), R2TreeFreeV2::KP});
                                }
                                // Update leafMeta
                                leafMeta[leafId] = {rWKK, rWPP, rWKP, (int)rPovit.size(), (int)rKeepC.size(), rNp};
                                leafMutated[leafId] = 1;  // keep mutated (Edge→Leaf CSR is stale)
                                leafEdgeInfoRebuilt[leafId] = 1;  // re-enable Case A fast path
                                rPovit.free(); rKeepC.free();
                            }
                        } else {
                            for (daf::Size i = 0; i+1 < leaf.size(); ++i)
                                for (daf::Size j = i+1; j < leaf.size(); ++j) {
                                    auto &u = leaf[i], &v = leaf[j];
                                    double d = (!u.isPivot && !v.isPivot) ? RemovedKtoK : (u.isPivot && v.isPivot) ? RemovedPtoP : RemovedKtoP;
                                    directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                                }
                            for (const auto& i : leaf)
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

        // === Phase 2c: BK for Case B ===
        auto tB0 = hrc::now();
        cntB += caseBLeafIds.size();
        for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
            auto leafId = caseBLeafIds[bi];
            if (leafId >= leafRmInfo.size()) continue;
            R2TreeFreeV2::LeafRmInfo &leafRm = leafRmInfo[leafId];

            bool isInitLeaf = leafId < numLeavesInit;

            // Opt 4: Always use tree.adj_list for Case B (no leafVtxData CSR)
            std::vector<TreeGraphNode> leafNodes;
            leafNodes.reserve(tree.adj_list[leafId].size());
            for (const auto& node : tree.adj_list[leafId])
                leafNodes.push_back(node);

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
            for (auto &lv : tree.adj_list[leafId]) {
                if (lv.isPivot) treeGraphV.removeNbr(lv.v, {leafId, true});
                else treeGraphV.removeNbr(lv.v, {leafId, false});
            }
            if (!leafRm.removedPivots.empty()) tree.removeNbrs(leafId, leafRm.removedPivots);

            auto &leafRef = tree.adj_list[leafId];
            bkRmEdge::bronKerbosch(leafRef, leafRm.removedEdges, k,
                [&](const bkRmEdge::Bitset &c, const bkRmEdge::Bitset &pivots) {
                    auto subLeaf = bkRmEdge::coverToVertex(c, pivots, leafRef);
                    auto newId = tree.addNode(subLeaf);
                    newPivot.clear(); newKeepC.clear();
                    for (const auto& i : tree.adj_list[newId]) {
                        if (i.isPivot) { newPivot.push_back(i.v); treeGraphV.addNbr(i.v, {newId, true}); }
                        else { newKeepC.push_back(i.v); treeGraphV.addNbr(i.v, {newId, false}); }
                    }
                    daf::Size np = k - newKeepC.size();
                    if (np <= newPivot.size() && newKeepC.size() > 1) {
                        double w = std::llround(nCr[newPivot.size()][np]);
                        R2TreeFreeV2::processEdgePairs(newKeepC, w, addW);
                    }
                    int nPP = int(np) - 2;
                    if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                        double w = std::llround(nCr[newPivot.size() - 2][nPP]);
                        R2TreeFreeV2::processEdgePairs(newPivot, w, addW);
                    }
                    int nKP = int(np) - 1;
                    if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                        double w = std::llround(nCr[newPivot.size() - 1][nKP]);
                        R2TreeFreeV2::processEdgePairs(newKeepC, newPivot, w, addW);
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
            if (needPivotBK <= povit.size()) { double w = std::llround(nCr[povit.size()][needPivotBK]); R2TreeFreeV2::processEdgePairs(keepC, w, removeW); }
            int needPP = int(needPivotBK) - 2;
            if (0 <= needPP && needPP <= int(povit.size()) - 2) { double w = std::llround(nCr[povit.size() - 2][needPP]); R2TreeFreeV2::processEdgePairs(povit, w, removeW); }
            int needKP = int(needPivotBK) - 1;
            if (0 <= needKP && needKP <= int(povit.size()) - 1) { double w = std::llround(nCr[povit.size() - 1][needKP]); R2TreeFreeV2::processEdgePairs(keepC, povit, w, removeW); }

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
    std::cout << "  [CaseC] rebuilds=" << cntC_rebuild << std::endl;
    std::cout << "  [Breakdown] phase1=" << (t_phase1_ns/1000000) << "ms"
              << " (edgeLeaf=" << (t_phase1_edgeLeaf_ns/1000000) << "ms"
              << ", overflow=" << (t_phase1_overflow_ns/1000000) << "ms)"
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
