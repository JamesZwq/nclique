//
// Local H-index iterative convergence for r=2 (edge) nucleus decomposition.
//
// Tree is immutable. Each edge independently computes its core value
// via H-index over its neighboring leaves, iterating until convergence.
// Uses exact s-clique enumeration: for each s-clique containing edge e,
// level = min(coreE[e'] for all other edges e' in the s-clique).
//
// Optimizations:
//   1. Integer bucket H-index (no sort)
//   2. Dirty queue with propagation
//   3. Early termination when H-index >= currentCore
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <vector>
#include <cstring>
#include <queue>
#include <cmath>

#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace ECD_Local {

// Initial support counting per edge (same as PlusECDSet_ST)
double * countingPerEdge(const DynamicGraph<TreeGraphNode> &treeGraph,
                            const Graph &edgeGraph,
                            const daf::CliqueSize k) {
    auto *countingE = new double[edgeGraph.adj_list.size()];
    memset(countingE, 0, edgeGraph.adj_list.size() * sizeof(double));

    daf::StaticVector<daf::Size> tPovit, tKeepC;
    for (const auto &clique : treeGraph.adj_list) {
        if (clique.size() < k) continue;
        tPovit.clear(); tKeepC.clear();
        for (const auto &node : clique) {
            if (node.isPivot) tPovit.push_back(node.v);
            else tKeepC.push_back(node.v);
        }
        int needPivot = int(k) - int(tKeepC.size());

        // K-K edges
        if (needPivot >= 0 && needPivot <= int(tPovit.size())) {
            double val = (nCr[tPovit.size()][needPivot]);
            for (size_t i = 0; i + 1 < tKeepC.size(); ++i)
                for (size_t j = i + 1; j < tKeepC.size(); ++j)
                    countingE[edgeGraph.getEdgeCompressedId(tKeepC[i], tKeepC[j])] += val;
        }
        // P-P edges
        int needPP = needPivot - 2;
        if (needPP >= 0 && needPP <= int(tPovit.size()) - 2) {
            double val = (nCr[tPovit.size() - 2][needPP]);
            for (size_t i = 0; i + 1 < tPovit.size(); ++i)
                for (size_t j = i + 1; j < tPovit.size(); ++j)
                    countingE[edgeGraph.getEdgeCompressedId(tPovit[i], tPovit[j])] += val;
        }
        // K-P edges
        int needKP = needPivot - 1;
        if (needKP >= 0 && needKP <= int(tPovit.size()) - 1) {
            double val = (nCr[tPovit.size() - 1][needKP]);
            for (size_t i = 0; i < tKeepC.size(); ++i)
                for (size_t j = 0; j < tPovit.size(); ++j)
                    countingE[edgeGraph.getEdgeCompressedId(tKeepC[i], tPovit[j])] += val;
        }
    }
    tPovit.free(); tKeepC.free();
    return countingE;
}

} // namespace ECD_Local


// ============================================================
// Compute H-index for edge e=(u,v) using exact s-clique enumeration.
// For each leaf containing both u and v:
//   - others = leaf \ {u,v}, classified as keeps/pivots
//   - enumerate all (s-2)-subsets from others (all keeps + choose pivots)
//   - for each s-clique, level = min(coreE[e'] for all other edges e')
//   - bucket[level] += 1
// ============================================================
static double computeEdgeHIndex(
    daf::Size edgeU, daf::Size edgeV, daf::Size selfEdgeId,
    const DynamicGraph<TreeGraphNode> &tree,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    const double *coreE,
    const Graph &edgeGraph,
    const daf::CliqueSize s,
    const std::vector<int> &leafNeedPivot,
    const std::vector<uint8_t> &leafValid,
    double currentCore,
    double &totalSupportOut,
    std::vector<TreeGraphNode> &otherKeeps,
    std::vector<TreeGraphNode> &otherPivots,
    std::vector<double> &buckets)
{
    totalSupportOut = 0;
    if (currentCore <= 0) return 0;

    int bucketSize = (int)currentCore;
    if (bucketSize < 1) bucketSize = 1;
    if ((int)buckets.size() < bucketSize + 1)
        buckets.resize(bucketSize + 1, 0);
    for (int i = 0; i <= bucketSize; ++i) buckets[i] = 0;

    double rawTotalSupport = 0;

    auto &adjU = treeGraphV.getNbr(edgeU);
    auto &adjV = treeGraphV.getNbr(edgeV);

    daf::intersect_dense_sets(adjU, adjV,
        [&](const TreeGraphNode &uClique, const TreeGraphNode &vClique) {
            daf::Size leafId = uClique.v;
            if (!leafValid[leafId]) return;

            const auto &leaf = tree.adj_list[leafId];
            int needPivot = leafNeedPivot[leafId];

            bool edgeUIsPivot = uClique.isPivot;
            bool edgeVIsPivot = vClique.isPivot;

            // Classify other vertices (leaf \ {u,v})
            otherKeeps.clear();
            otherPivots.clear();
            for (const auto &node : leaf) {
                if (node.v == edgeU || node.v == edgeV) continue;
                if (node.isPivot) otherPivots.push_back(node);
                else otherKeeps.push_back(node);
            }

            // effectiveNeedPivot for the (s-2) additional vertices
            int effectiveNeedPivot = needPivot;
            if (edgeUIsPivot) effectiveNeedPivot--;
            if (edgeVIsPivot) effectiveNeedPivot--;

            int numOtherKeeps = (int)otherKeeps.size();
            int numOtherPivots = (int)otherPivots.size();
            int need = (int)s - 2; // how many more vertices

            if (effectiveNeedPivot < 0 || effectiveNeedPivot > numOtherPivots) return;
            if (numOtherKeeps + numOtherPivots < need) return;
            if (numOtherKeeps > need) return;

            // Helper lambda: given chosen pivots, compute level and add to bucket
            auto processSClique = [&](const TreeGraphNode *chosenPivotData, int chosenPivotCount) {
                daf::StaticVector<daf::Size> sVerts;
                sVerts.push_back(edgeU);
                sVerts.push_back(edgeV);
                for (const auto &k : otherKeeps) sVerts.push_back(k.v);
                for (int i = 0; i < chosenPivotCount; ++i) sVerts.push_back(chosenPivotData[i].v);

                double minOtherCore = 1e18;
                int nv = (int)sVerts.size();
                for (int i = 0; i < nv; ++i) {
                    for (int j = i + 1; j < nv; ++j) {
                        auto eid = edgeGraph.getEdgeCompressedId(sVerts[i], sVerts[j]);
                        if (eid == selfEdgeId) continue;
                        if (eid == std::numeric_limits<daf::Size>::max()) continue;
                        double c = coreE[eid];
                        if (c < minOtherCore) minOtherCore = c;
                    }
                }
                sVerts.free();

                if (minOtherCore > 9e17) return;
                if (minOtherCore < 1) return;

                int level = (minOtherCore > bucketSize) ? bucketSize : (int)minOtherCore;
                buckets[level] += 1;
                rawTotalSupport += 1;
            };

            if (effectiveNeedPivot == 0) {
                // s-clique = {u, v} + all otherKeeps, no pivots needed
                processSClique(nullptr, 0);
            } else {
                // Enumerate all combinations of effectiveNeedPivot pivots from otherPivots
                daf::enumerateCombinations(otherPivots, effectiveNeedPivot,
                    [&](const daf::StaticVector<TreeGraphNode> &chosenPivots) {
                        processSClique(chosenPivots.getData(), (int)chosenPivots.size());
                        return true;
                    });
            }
        });

    totalSupportOut = rawTotalSupport;
    if (rawTotalSupport < 1) return 0;

    // Bucket H-index scan (high to low)
    double accumulated = 0;
    for (int c = bucketSize; c >= 1; --c) {
        accumulated += buckets[c];
        if (accumulated >= c) return c;
    }
    return 0;
}


// ============================================================
// r=2 Local H-index with dirty-queue propagation
// ============================================================
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>>
NCliqueEdgeCoreDecomposition_Local(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k) {

    auto time_start = std::chrono::high_resolution_clock::now();

    const daf::Size numLeaves = tree.adj_list.size();
    const daf::Size numEdges = edgeGraph.adj_list.size();

    // Precompute per-leaf metadata (immutable)
    std::vector<int> leafNeedPivot(numLeaves);
    std::vector<uint8_t> leafValid(numLeaves);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        int keeps = 0, pivots = 0;
        for (const auto &node : tree.adj_list[L]) {
            if (node.isPivot) pivots++;
            else keeps++;
        }
        leafNeedPivot[L] = (int)k - keeps;
        leafValid[L] = (leafNeedPivot[L] >= 0 && leafNeedPivot[L] <= pivots) ? 1 : 0;
    }

    // Initial support
    auto *countingE = ECD_Local::countingPerEdge(tree, edgeGraph, k);
    auto *coreE = new double[numEdges];
    for (daf::Size i = 0; i < numEdges; ++i)
        coreE[i] = countingE[i];
    delete[] countingE;

    // Scratch buffers
    std::vector<TreeGraphNode> otherKeeps, otherPivots;
    std::vector<double> buckets;
    otherKeeps.reserve(64);
    otherPivots.reserve(64);
    buckets.reserve(4096);

    // Dirty queue
    std::vector<uint8_t> inQueue(numEdges, 0);
    std::queue<daf::Size> dirtyQueue;

    const daf::Size n = edgeGraph.adj_list_offsets.size() - 1;

    std::cout << "=========================begin (Local H-index r=2)=========================" << std::endl;
    std::cout << "edges: " << numEdges << ", leaves: " << numLeaves << std::endl;

    // Phase 1: initial full scan
    for (daf::Size u = 0; u < n; ++u) {
        const daf::Size start = edgeGraph.adj_list_offsets[u];
        const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
        for (daf::Size idx = start; idx < end; ++idx) {
            if (coreE[idx] <= 0) continue;
            daf::Size v = edgeGraph.adj_list[idx];

            double totalSup = 0;
            double newCore = computeEdgeHIndex(
                u, v, idx, tree, treeGraphV, coreE, edgeGraph, k,
                leafNeedPivot, leafValid, coreE[idx], totalSup,
                otherKeeps, otherPivots, buckets);
            newCore = std::min(newCore, coreE[idx]);

            if (newCore != coreE[idx]) {
                double oldCore = coreE[idx];
                coreE[idx] = newCore;
                // Enqueue co-leaf edges
                auto &adjU = treeGraphV.getNbr(u);
                auto &adjV = treeGraphV.getNbr(v);
                daf::intersect_dense_sets(adjU, adjV,
                    [&](const TreeGraphNode &uC, const TreeGraphNode &vC) {
                        daf::Size leafId = uC.v;
                        if (!leafValid[leafId]) return;
                        const auto &leaf = tree.adj_list[leafId];
                        for (size_t i = 0; i < leaf.size(); ++i) {
                            for (size_t j = i + 1; j < leaf.size(); ++j) {
                                auto eid = edgeGraph.getEdgeCompressedId(leaf[i].v, leaf[j].v);
                                if (eid != std::numeric_limits<daf::Size>::max() &&
                                    coreE[eid] > 0 && oldCore >= coreE[eid] && !inQueue[eid]) {
                                    inQueue[eid] = 1;
                                    dirtyQueue.push(eid);
                                }
                            }
                        }
                    });
            }
        }
    }

    // Phase 2: dirty queue propagation
    long long recomputeCount = 0;
    int iteration = 1;

    while (!dirtyQueue.empty()) {
        iteration++;
        size_t roundSize = dirtyQueue.size();
        for (size_t qi = 0; qi < roundSize; ++qi) {
            daf::Size edgeId = dirtyQueue.front();
            dirtyQueue.pop();
            inQueue[edgeId] = 0;
            recomputeCount++;

            if (coreE[edgeId] <= 0) continue;

            auto [u, v] = edgeGraph.getEdgeById(edgeId);

            double totalSup = 0;
            double newCore = computeEdgeHIndex(
                u, v, edgeId, tree, treeGraphV, coreE, edgeGraph, k,
                leafNeedPivot, leafValid, coreE[edgeId], totalSup,
                otherKeeps, otherPivots, buckets);
            newCore = std::min(newCore, coreE[edgeId]);

            if (newCore != coreE[edgeId]) {
                double oldCore = coreE[edgeId];
                coreE[edgeId] = newCore;
                // Propagate
                auto &adjU = treeGraphV.getNbr(u);
                auto &adjV = treeGraphV.getNbr(v);
                daf::intersect_dense_sets(adjU, adjV,
                    [&](const TreeGraphNode &uC, const TreeGraphNode &vC) {
                        daf::Size leafId = uC.v;
                        if (!leafValid[leafId]) return;
                        const auto &leaf = tree.adj_list[leafId];
                        for (size_t i = 0; i < leaf.size(); ++i) {
                            for (size_t j = i + 1; j < leaf.size(); ++j) {
                                auto eid = edgeGraph.getEdgeCompressedId(leaf[i].v, leaf[j].v);
                                if (eid != std::numeric_limits<daf::Size>::max() &&
                                    coreE[eid] > 0 && oldCore >= coreE[eid] && !inQueue[eid]) {
                                    inQueue[eid] = 1;
                                    dirtyQueue.push(eid);
                                }
                            }
                        }
                    });
            }
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "Local H-index r=2 converged in " << iteration << " iterations, "
              << recomputeCount << " edge re-evaluations" << std::endl;
    std::cout << "time: " << elapsed << " ms" << std::endl;

    // Build output
    std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> result;
    result.reserve(numEdges);
    for (daf::Size u = 0; u < n; ++u) {
        const daf::Size start = edgeGraph.adj_list_offsets[u];
        const daf::Size end = edgeGraph.adj_list_offsets[u + 1];
        for (daf::Size idx = start; idx < end; ++idx) {
            result.emplace_back(std::make_pair(u, edgeGraph.adj_list[idx]), (int)coreE[idx]);
        }
    }

    delete[] coreE;
    return result;
}
