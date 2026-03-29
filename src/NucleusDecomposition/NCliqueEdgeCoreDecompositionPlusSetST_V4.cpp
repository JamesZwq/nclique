//
// R2 ST V4: Edge-Leaf Dual Index for Immutable Case A
//
// Key optimization: Case A (leaf death, 94.6% of events) uses a pre-built
// leafEdgeInfo[leafId] scan instead of O(n²) vertex-pair enumeration.
// Each entry stores {edgeId, edgeType} so the weight lookup is O(1).
//
// Phase 1, Case B, Case C: unchanged from baseline (proven correct).
// Deferred batch bucketMove from V1.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <cassert>

#include "../BK/BronKerboschRmEdge.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace R2STV4 {

enum EdgeType : uint8_t { KK = 0, PP = 1, KP = 2 };

struct LeafEdgeEntry {
    daf::Size edgeId;
    EdgeType type;
};

// Per-leaf metadata for Case A fast path
struct LeafMeta {
    double wKK, wPP, wKP;
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

// Fused countingPerEdge + leafEdgeInfo construction
struct InitResult {
    double *countingKE;
    std::vector<std::vector<LeafEdgeEntry>> leafEdgeInfo;
    std::vector<LeafMeta> leafMeta;
};

InitResult countingPerEdgeWithIndex(
    const DynamicGraph<TreeGraphNode> &treeGraph,
    const Graph &edgeGraph,
    const daf::CliqueSize k) {

    const daf::Size numEdges = edgeGraph.adj_list.size();
    const daf::Size numLeaves = treeGraph.adj_list.size();

    auto *countingE = new double[numEdges];
    memset(countingE, 0, numEdges * sizeof(double));

    std::vector<std::vector<LeafEdgeEntry>> leafEdgeInfo(numLeaves);
    std::vector<LeafMeta> leafMeta(numLeaves, {0, 0, 0});

    daf::StaticVector<daf::Size> tPovit, tKeepC;

    for (daf::Size li = 0; li < numLeaves; ++li) {
        const auto &clique = treeGraph.adj_list[li];
        if (clique.size() < k) continue;

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
                    countingE[eid] += wKK;
                    leafEdgeInfo[li].push_back({eid, KK});
                }
        }

        int needPP = needPivot - 2;
        if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
            wPP = std::llround(nCr[tPovit.size() - 2][needPP]);
            for (size_t i = 0; i < tPovit.size(); ++i)
                for (size_t j = i + 1; j < tPovit.size(); ++j) {
                    auto eid = edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]);
                    countingE[eid] += wPP;
                    leafEdgeInfo[li].push_back({eid, PP});
                }
        }

        int needKP = needPivot - 1;
        if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
            wKP = std::llround(nCr[tPovit.size() - 1][needKP]);
            for (size_t i = 0; i < tKeepC.size(); ++i)
                for (size_t j = 0; j < tPovit.size(); ++j) {
                    auto eid = edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]);
                    countingE[eid] += wKP;
                    leafEdgeInfo[li].push_back({eid, KP});
                }
        }

        leafMeta[li] = {wKK, wPP, wKP};
    }

    tPovit.free(); tKeepC.free();
    return {countingE, std::move(leafEdgeInfo), std::move(leafMeta)};
}

} // namespace R2STV4


std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V4(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    // === Init: fused countingPerEdge + leafEdgeInfo ===
    auto initResult = R2STV4::countingPerEdgeWithIndex(tree, edgeGraph, k);
    auto *countingKE = initResult.countingKE;
    auto &leafEdgeInfo = initResult.leafEdgeInfo;
    auto &leafMeta = initResult.leafMeta;

    const daf::Size numEdgesInit = edgeGraph.adj_list.size();

    auto *coreE = new double[numEdgesInit];
    memset(coreE, 0, numEdgesInit * sizeof(double));

    // Track leaves mutated by Case C/B — fast path only safe for unmutated leaves
    const daf::Size numLeavesInit = tree.adj_list.size();
    std::vector<uint8_t> leafMutated(numLeavesInit, 0);

    daf::StaticVector<daf::Size> povit, keepC, newPivot, newKeepC;
    daf::StaticVector<daf::Size> currentRemoveEdgeIds(numEdgesInit);

    daf::StaticVector<uint8_t> edgeInHeap(numEdgesInit);
    edgeInHeap.c_size = numEdgesInit;
    memset(edgeInHeap.getData(), 1, numEdgesInit * sizeof(uint8_t));

    daf::StaticVector<daf::Size> removedLeaf(tree.adj_list.size());
    daf::StaticVector<R2STV4::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();

    // --- Bucket array ---
    const daf::Size numEdges = numEdgesInit;
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
        auto& oldVec = buckets[oldB]; daf::Size myPos = pos_in_bucket[id];
        if (myPos < oldVec.size() - 1) { daf::Size last = oldVec.back(); oldVec[myPos] = last; pos_in_bucket[last] = myPos; }
        oldVec.pop_back(); if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
        bucket_of[id] = newB; pos_in_bucket[id] = buckets[newB].size(); buckets[newB].push_back(id);
        if (newB < curBucket) curBucket = newB;
    };

    // Deferred batch bucketMove
    std::vector<uint8_t> dirtyMark(numEdges, 0);
    std::vector<daf::Size> dirtyEdges;
    dirtyEdges.reserve(8192);

    std::vector<uint8_t> leafAffected(leafRmInfo.size(), 0);
    std::vector<daf::Size> caseBLeafIds;
    caseBLeafIds.reserve(1024);

    daf::log_memory("R2 ST V4 init (with leaf-edge index)");

    std::cout << "=========================begin (R2 ST V4)=========================" << std::endl;
    double minCore = 0, cntA = 0, cntB = 0, cntC = 0, totalIters = 0;

    while (remainingInHeap > 0) {
        // --- Bucket pop ---
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

        // === Phase 1: baseline intersect_dense_sets (proven correct) ===
        for (int ei = 0; ei < (int)currentRemoveEdgeIds.size(); ++ei) {
            auto edgeId = currentRemoveEdgeIds[ei];
            auto [edgeU, edgeV] = edgeGraph.getEdgeById(edgeId);
            auto &adjEdgeU = treeGraphV.getNbr(edgeU);
            auto &adjEdgeV = treeGraphV.getNbr(edgeV);
            daf::intersect_dense_sets(adjEdgeU, adjEdgeV,
                [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
                    daf::Size leafId = uClique.v;
                    auto &lrm = leafRmInfo[leafId];
                    if (lrm.empty()) { removedLeaf.push_back(leafId); leafAffected[leafId] = 1; }
                    if (!lrm.removedKeepC) {
                        if (!uClique.isPivot && !vClique.isPivot) lrm.removedKeepC = true;
                        else if (uClique.isPivot && vClique.isPivot) lrm.removedEdges.push_back({edgeU, edgeV});
                        else if (uClique.isPivot) lrm.removedPivots.push_back(edgeU);
                        else lrm.removedPivots.push_back(edgeV);
                    }
                });
        }
        // Pre-sort removedPivots
        for (int i = 0; i < (int)removedLeaf.size(); ++i) {
            auto &rp = leafRmInfo[removedLeaf[i]].removedPivots;
            if (rp.size() == 2) { if (rp[0] > rp[1]) std::swap(rp[0], rp[1]); if (rp[0] == rp[1]) rp.c_size = 1; }
            else if (rp.size() > 2) { std::ranges::sort(rp); rp.unique(); }
        }

        // === Phase 2: Case A (leafEdgeInfo scan) & Case C (baseline O(n²)) ===
        caseBLeafIds.clear();
        {
            daf::StaticVector<daf::Size> tPovit, tKeepC;

            auto directSub = [&](daf::Size idx, double w) {
                countingKE[idx] -= w; if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx] && !dirtyMark[idx]) { dirtyMark[idx] = 1; dirtyEdges.push_back(idx); }
            };

            for (int leafIdIdx = 0; leafIdIdx < (int)removedLeaf.size(); ++leafIdIdx) {
                auto leafId = removedLeaf[leafIdIdx];
                R2STV4::LeafRmInfo &leafRm = leafRmInfo[leafId];

                const auto& leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;

                tPovit.clear(); tKeepC.clear();
                daf::Size numKeeps = 0;
                for (const auto& node : leaf) {
                    if (node.isPivot) tPovit.push_back(node.v);
                    else { tKeepC.push_back(node.v); numKeeps++; }
                }
                daf::Size needPivot = k - numKeeps;
                daf::Size numPivots = tPovit.size();

                bool isDeadLeaf = leafRm.removedKeepC || needPivot > numPivots - leafRm.removedPivots.size();
                bool isCaseB = !leafRm.removedEdges.empty() && !isDeadLeaf;

                if (isCaseB) { caseBLeafIds.push_back(leafId); continue; }

                if (isDeadLeaf) {
                    // ==== Case A: leaf dies — subtract full contribution ====
                    cntA++;

                    // Compute current weights from live tree state
                    double wKK = 0, wPP = 0, wKP = 0;
                    if (needPivot <= tPovit.size())
                        wKK = std::llround(nCr[tPovit.size()][needPivot]);
                    int nPP = int(needPivot) - 2;
                    if (0 <= nPP && nPP <= int(tPovit.size()) - 2)
                        wPP = std::llround(nCr[tPovit.size() - 2][nPP]);
                    int nKP = int(needPivot) - 1;
                    if (0 <= nKP && nKP <= int(tPovit.size()) - 1)
                        wKP = std::llround(nCr[tPovit.size() - 1][nKP]);

                    // Use leafEdgeInfo scan only for unmutated init leaves
                    if (leafId < numLeavesInit && !leafMutated[leafId]
                        && leafId < leafEdgeInfo.size() && !leafEdgeInfo[leafId].empty()) {
                        // Fast path: sequential scan — leaf was never mutated
                        for (const auto &entry : leafEdgeInfo[leafId]) {
                            if (!edgeInHeap[entry.edgeId]) continue;
                            double w;
                            switch (entry.type) {
                                case R2STV4::KK: w = wKK; break;
                                case R2STV4::PP: w = wPP; break;
                                case R2STV4::KP: w = wKP; break;
                            }
                            directSub(entry.edgeId, w);
                        }
                    } else {
                        // Fallback: O(n²) pair enumeration
                        if (wKK > 0) {
                            for (daf::Size i = 0; i+1 < tKeepC.size(); ++i)
                                for (daf::Size j = i+1; j < tKeepC.size(); ++j)
                                    directSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j]), wKK);
                        }
                        if (wPP > 0) {
                            for (daf::Size i = 0; i+1 < tPovit.size(); ++i)
                                for (daf::Size j = i+1; j < tPovit.size(); ++j)
                                    directSub(edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j]), wPP);
                        }
                        if (wKP > 0) {
                            for (daf::Size i = 0; i < tKeepC.size(); ++i)
                                for (daf::Size j = 0; j < tPovit.size(); ++j)
                                    directSub(edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j]), wKP);
                        }
                    }
                    // Tree mutation (same as baseline)
                    // Clear leafEdgeInfo to prevent stale data after recycleNode reuses this ID
                    if (leafId < leafEdgeInfo.size())
                        leafEdgeInfo[leafId].clear();
                    for (const auto& i : leaf)
                        treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                    tree.adj_list[leafId].clear();
                    tree.recycleNode(leafId);
                } else {
                    // ==== Case C: baseline O(n²) (unchanged, 2.2% of events) ====
                    cntC++;
                    if (leafId < numLeavesInit) leafMutated[leafId] = 1;
                    double KtoK = 0, KtoP = 0, PtoP = 0;
                    double RemovedKtoK = 0, RemovedKtoP = 0, RemovedPtoP = 0;

                    if (needPivot <= tPovit.size()) {
                        KtoK = std::llround(nCr[tPovit.size()][needPivot]) - std::llround(nCr[tPovit.size() - leafRm.removedPivots.size()][needPivot]);
                        RemovedKtoK = std::llround(nCr[tPovit.size()][needPivot]);
                    }
                    int needPP = int(needPivot) - 2;
                    if (0 <= needPP && needPP <= int(tPovit.size()) - 2) {
                        RemovedPtoP = std::llround(nCr[tPovit.size() - 2][needPP]);
                        PtoP = RemovedPtoP;
                        if (leafRm.removedPivots.size() + 2 <= tPovit.size())
                            PtoP -= std::llround(nCr[tPovit.size() - 2 - leafRm.removedPivots.size()][needPP]);
                    }
                    int needKP = int(needPivot) - 1;
                    if (0 <= needKP && needKP <= int(tPovit.size()) - 1) {
                        RemovedKtoP = std::llround(nCr[tPovit.size() - 1][needKP]);
                        KtoP = RemovedKtoP;
                        if (leafRm.removedPivots.size() + 1 <= tPovit.size())
                            KtoP -= std::llround(nCr[tPovit.size() - 1 - leafRm.removedPivots.size()][needKP]);
                    }

                    if (!leafRm.removedPivots.empty() && needPivot <= tPovit.size() - leafRm.removedPivots.size()) {
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
                    } else {
                        for (daf::Size i = 0; i+1 < leaf.size(); ++i)
                            for (daf::Size j = i+1; j < leaf.size(); ++j) {
                                auto &u = leaf[i], &v = leaf[j];
                                double d = (!u.isPivot && !v.isPivot) ? RemovedKtoK : (u.isPivot && v.isPivot) ? RemovedPtoP : RemovedKtoP;
                                directSub(edgeGraph.getEdgeCompressedId(u.v, v.v), d);
                            }
                        for (const auto& i : leaf)
                            treeGraphV.removeNbr(i.v, static_cast<TreeGraphNode>(leafId));
                        if (leafId < leafEdgeInfo.size()) leafEdgeInfo[leafId].clear();
                        tree.adj_list[leafId].clear();
                        tree.recycleNode(leafId);
                    }
                }
                leafRmInfo[leafId].clear();
            }
            tPovit.free(); tKeepC.free();

            // Flush deferred bucket moves
            for (auto eid : dirtyEdges) { if (edgeInHeap[eid]) bucketMove(eid); dirtyMark[eid] = 0; }
            dirtyEdges.clear();
        }

        // === Phase 2c: BK for Case B (unchanged from baseline) ===
        cntB += caseBLeafIds.size();
        for (int bi = 0; bi < (int)caseBLeafIds.size(); ++bi) {
            auto leafId = caseBLeafIds[bi];
            R2STV4::LeafRmInfo &leafRm = leafRmInfo[leafId];
            const auto& leaf = tree.adj_list[leafId];

            povit.clear(); keepC.clear();
            for (const auto& node : leaf) { if (node.isPivot) povit.push_back(node.v); else keepC.push_back(node.v); }
            daf::Size needPivot = k - keepC.size();

            auto addW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] += w; if (edgeInHeap[idx]) bucketMove(idx);
            };

            for (const auto& leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
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
                        R2STV4::processEdgePairs(newKeepC, w, addW);
                    }
                    int nPP = int(np) - 2;
                    if (0 <= nPP && nPP <= int(newPivot.size()) - 2) {
                        double w = std::llround(nCr[newPivot.size() - 2][nPP]);
                        R2STV4::processEdgePairs(newPivot, w, addW);
                    }
                    int nKP = int(np) - 1;
                    if (0 <= nKP && nKP <= int(newPivot.size()) - 1) {
                        double w = std::llround(nCr[newPivot.size() - 1][nKP]);
                        R2STV4::processEdgePairs(newKeepC, newPivot, w, addW);
                    }
                    newPivot.clear(); newKeepC.clear();
                    if (newId >= leafRmInfo.size()) {
                        removedLeaf.reserve(newId * 1.5);
                        leafRmInfo.resize(newId * 1.5);
                        leafAffected.resize(newId * 1.5, 0);
                    }
                });

            auto removeW = [&](daf::Size u, daf::Size v, double w) {
                auto idx = edgeGraph.getEdgeCompressedId(u, v);
                countingKE[idx] -= w; if (countingKE[idx] < 0) countingKE[idx] = 0;
                if (edgeInHeap[idx]) bucketMove(idx);
            };
            if (needPivot <= povit.size()) { double w = std::llround(nCr[povit.size()][needPivot]); R2STV4::processEdgePairs(keepC, w, removeW); }
            int needPP = int(needPivot) - 2;
            if (0 <= needPP && needPP <= int(povit.size()) - 2) { double w = std::llround(nCr[povit.size() - 2][needPP]); R2STV4::processEdgePairs(povit, w, removeW); }
            int needKP = int(needPivot) - 1;
            if (0 <= needKP && needKP <= int(povit.size()) - 1) { double w = std::llround(nCr[povit.size() - 1][needKP]); R2STV4::processEdgePairs(keepC, povit, w, removeW); }

            tree.removeNode(leafId);
            if (leafId < leafEdgeInfo.size()) leafEdgeInfo[leafId].clear();
            leafRmInfo[leafId].clear();
            povit.clear(); keepC.clear();
        }

        currentRemoveEdgeIds.clear();
        for (auto leafId : removedLeaf) leafAffected[leafId] = 0;
        removedLeaf.clear();
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;
    std::cout << "  Cases: A=" << cntA << " B=" << cntB << " C=" << cntC << " iters=" << totalIters << std::endl;

    daf::Size numCounting = 0;
    std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> sortedK;
    sortedK.reserve(numEdgesInit);
    const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;
    for (daf::Size u = 0; u < n; ++u) {
        const daf::Size start = edgeGraph.adj_list_offsets[u];
        const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
        for (daf::Size idx = start; idx < end; ++idx)
            sortedK.emplace_back(std::make_pair(std::make_pair(u, edgeGraph.adj_list[idx]), (int) coreE[idx]));
    }
    assert(numCounting == 0);
    delete[] countingKE; delete[] coreE;
    povit.free(); keepC.free(); newPivot.free(); newKeepC.free();
    return sortedK;
}
