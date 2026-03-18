//
// Single-thread optimized r=1 nucleus decomposition.
// Stripped of all OMP infrastructure: no locks, no per-thread vectors,
// no deferred dirty tracking, no atomic directives.
// Direct countingV updates + immediate bucketMove.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>

#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace VCD_ST {
    struct LeafRmInfo {
        bool removedKeepC;
        daf::StaticVector<daf::Size> removedPivots{0};

        LeafRmInfo() : removedKeepC(false) {}
        bool empty() const {
            return !removedKeepC && removedPivots.empty();
        }

        void init(auto capacity = 400) {
            removedKeepC = false;
            removedPivots.reserve(capacity);
        }

        void clear() {
            removedKeepC = false;
            removedPivots.clear();
        }
    };

    double * countingPerVertex(const DynamicGraph<TreeGraphNode> &treeGraph,
                                                     const Graph &edgeGraph,
                                                     const daf::CliqueSize k) {
        double *countingV = new double[edgeGraph.adj_list_offsets.size()];
        memset(countingV, 0, edgeGraph.adj_list_offsets.size() * sizeof(double));
        daf::StaticVector<daf::Size> povit;
        daf::StaticVector<daf::Size> keepC;
        for (const auto &clique: treeGraph.adj_list) {
            povit.clear();
            keepC.clear();
            if (clique.size() < k) {
                continue;
            }
            for (auto &i: clique) {
                if (i.isPivot) {
                    povit.push_back(i.v);
                } else {
                    keepC.push_back(i.v);
                }
            }

            int needPivot = int(k) - int(keepC.size());
            for (const auto &v: keepC) {
                countingV[v] += nCr[povit.size()][needPivot];
            }
            double eachPivotKcliques = 0;
            const int needPivotWithV = needPivot - 1;
            if (needPivotWithV >= 0 && needPivotWithV <= static_cast<int>(povit.size()) - 1) {
                eachPivotKcliques = nCr[povit.size() - 1][needPivotWithV];
            }
            for (const auto &v: povit) {
                countingV[v] += eachPivotKcliques;
            }
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

    auto countingV = VCD_ST::countingPerVertex(tree, edgeGraph, k);
    const daf::Size numVertices = edgeGraph.adj_list_offsets.size() - 1;
    auto coreV = new double[numVertices + 1];
    for (daf::Size i = 0; i <= numVertices; ++i) coreV[i] = -1.0;

    // --- Bucket array ---
    int maxBucket = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] > 0)
            maxBucket = std::max(maxBucket, (int)countingV[i]);
    }
    std::vector<std::vector<daf::Size>> buckets(maxBucket + 2);
    std::vector<int> bucket_of(numVertices);
    std::vector<daf::Size> pos_in_bucket(numVertices);
    std::vector<bool> vertexInHeap(numVertices, false);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] <= 0) continue;
        int b = (int)countingV[i];
        bucket_of[i] = b;
        pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i);
        vertexInHeap[i] = true;
        remainingInHeap++;
    }
    int curBucket = 0;

    // Bucket move helper — called immediately inline (no deferred tracking)
    auto bucketMove = [&](daf::Size id) {
        int newB = std::max(0, (int)countingV[id]);
        int oldB = bucket_of[id];
        if (newB == oldB) return;
        auto& oldVec = buckets[oldB];
        daf::Size myPos = pos_in_bucket[id];
        if (myPos < oldVec.size() - 1) {
            daf::Size last = oldVec.back();
            oldVec[myPos] = last;
            pos_in_bucket[last] = myPos;
        }
        oldVec.pop_back();
        if (newB >= (int)buckets.size()) buckets.resize(newB + 1);
        bucket_of[id] = newB;
        pos_in_bucket[id] = buckets[newB].size();
        buckets[newB].push_back(id);
        if (newB < curBucket) curBucket = newB;
    };

    // Leaf tracking
    daf::StaticVector<daf::Size> removedLeaf(tree.adj_list.size());
    daf::StaticVector<VCD_ST::LeafRmInfo> leafRmInfo(tree.adj_list.size());
    leafRmInfo.c_size = tree.adj_list.size();
    std::vector<bool> leafAffected(tree.adj_list.size(), false);

    daf::StaticVector<daf::Size> currentRemoveVertexIds(numVertices);

    std::cout << "=========================begin=========================" << std::endl;
    std::cout << "vertices in heap: " << remainingInHeap << std::endl;

    double minCore = 0;

    while (remainingInHeap > 0) {
        // --- Bucket pop ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((double)curBucket, minCore);

        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && (double)curBucket <= minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                vertexInHeap[id] = false;
                currentRemoveVertexIds.push_back(id);
                coreV[id] = minCore;
                remainingInHeap--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (double)(curBucket + 1) <= minCore) {
                curBucket++;
            } else break;
        }

        if (remainingInHeap == 0) break;

        // --- Phase 1: identify affected leaves (serial, no locks) ---
        for (int vi = 0; vi < (int)currentRemoveVertexIds.size(); ++vi) {
            auto v = currentRemoveVertexIds[vi];
            auto &adjClique = treeGraphV.getNbr(v);
            for (const auto &clique : adjClique) {
                daf::Size leafId = clique.v;
                if (tree.adj_list[leafId].empty()) continue;
                auto &lrm = leafRmInfo[leafId];
                if (lrm.empty()) {
                    removedLeaf.push_back(leafId);
                    leafAffected[leafId] = true;
                }
                if (!lrm.removedKeepC) {
                    if (!clique.isPivot) {
                        lrm.removedKeepC = true;
                    } else {
                        lrm.removedPivots.push_back(v);
                    }
                }
            }
        }

        // --- Phase 2: countingV deltas + tree mutations + immediate bucket moves ---
        for (int li = 0; li < (int)removedLeaf.size(); ++li) {
            auto leafId = removedLeaf[li];
            auto &leaf = tree.adj_list[leafId];
            VCD_ST::LeafRmInfo &leafRm = leafRmInfo[leafId];
            if (leaf.empty()) continue;

            // Sort removedPivots inline
            std::ranges::sort(leafRm.removedPivots);
            leafRm.removedPivots.unique();

            int numPivots = 0, numKeeps = 0;
            for (auto &node : leaf) {
                if (node.isPivot) numPivots++;
                else numKeeps++;
            }

            int needPivot = (int)k - numKeeps;
            int remainPivots = numPivots - (int)leafRm.removedPivots.size();
            bool leafDies = leafRm.removedKeepC || needPivot < 0 || needPivot > remainPivots;

            if (leafDies) {
                double keepValue = (needPivot >= 0 && needPivot <= numPivots)
                                   ? nCr[numPivots][needPivot] : 0.0;
                double pivotValue = (needPivot >= 1 && needPivot - 1 <= numPivots - 1)
                                    ? nCr[numPivots - 1][needPivot - 1] : 0.0;

                for (auto &node : leaf) {
                    double delta = node.isPivot ? pivotValue : keepValue;
                    if (delta > 0) {
                        countingV[node.v] -= delta;
                        countingV[node.v] = std::max(countingV[node.v], 0.0);
                        if (vertexInHeap[node.v])
                            bucketMove(node.v);
                    }
                }
                leaf.clear();
                tree.recycleNode(leafId);
            } else if (!leafRm.removedPivots.empty() && needPivot <= remainPivots) {
                double KeepDelta = nCr[numPivots][needPivot] - nCr[remainPivots][needPivot];
                double RemovedPivotFull = (needPivot >= 1) ? nCr[numPivots - 1][needPivot - 1] : 0.0;
                double PivotDelta = (needPivot >= 1)
                    ? nCr[numPivots - 1][needPivot - 1] - nCr[remainPivots - 1][needPivot - 1]
                    : 0.0;

                for (auto rp : leafRm.removedPivots) {
                    if (RemovedPivotFull > 0) {
                        countingV[rp] -= RemovedPivotFull;
                        countingV[rp] = std::max(countingV[rp], 0.0);
                        if (vertexInHeap[rp])
                            bucketMove(rp);
                    }
                }

                for (auto &node : leaf) {
                    if (node.isPivot && std::binary_search(
                            leafRm.removedPivots.begin(),
                            leafRm.removedPivots.end(), node.v))
                        continue;
                    double delta = node.isPivot ? PivotDelta : KeepDelta;
                    if (delta > 0) {
                        countingV[node.v] -= delta;
                        countingV[node.v] = std::max(countingV[node.v], 0.0);
                        if (vertexInHeap[node.v])
                            bucketMove(node.v);
                    }
                }
                tree.removeNbrs(leafId, leafRm.removedPivots);
            }
        }

        // --- Cleanup ---
        for (int li = 0; li < (int)removedLeaf.size(); ++li) {
            leafRmInfo[removedLeaf[li]].clear();
            leafAffected[removedLeaf[li]] = false;
        }
        currentRemoveVertexIds.clear();
        removedLeaf.clear();
    }

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count() << " ms" << std::endl;

    delete[] countingV;
    currentRemoveVertexIds.free();
    removedLeaf.free();
    leafRmInfo.free();
    return coreV;
}
