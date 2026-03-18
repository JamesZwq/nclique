//
// Single-thread optimized r=1 nucleus decomposition.
// Optimizations over parallel version:
//   - No OMP: no locks, no per-thread vectors, no atomic directives
//   - Immediate bucketMove (no deferred dirty tracking)
//   - Integer arithmetic (long long countingV) eliminates float overhead
//   - vector<uint8_t> instead of vector<bool> for flags
//   - Small-array fast path: skip sort for ≤2 removedPivots
//   - Merged Phase 2 (delta + tree mutation in single pass)
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

    long long * countingPerVertex(const DynamicGraph<TreeGraphNode> &treeGraph,
                                  const Graph &edgeGraph,
                                  const daf::CliqueSize k) {
        const daf::Size n = edgeGraph.adj_list_offsets.size();
        auto *countingV = new long long[n];
        memset(countingV, 0, n * sizeof(long long));
        daf::StaticVector<daf::Size> povit;
        daf::StaticVector<daf::Size> keepC;
        for (const auto &clique: treeGraph.adj_list) {
            povit.clear();
            keepC.clear();
            if (clique.size() < k) continue;
            for (auto &i: clique) {
                if (i.isPivot) povit.push_back(i.v);
                else keepC.push_back(i.v);
            }
            int needPivot = int(k) - int(keepC.size());
            long long keepVal = std::llround(nCr[povit.size()][needPivot]);
            for (const auto &v: keepC)
                countingV[v] += keepVal;
            long long pivotVal = 0;
            const int needPivotWithV = needPivot - 1;
            if (needPivotWithV >= 0 && needPivotWithV <= static_cast<int>(povit.size()) - 1)
                pivotVal = std::llround(nCr[povit.size() - 1][needPivotWithV]);
            for (const auto &v: povit)
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
    std::vector<uint8_t> vertexInHeap(numVertices, 0);
    daf::Size remainingInHeap = 0;
    for (daf::Size i = 0; i < numVertices; ++i) {
        if (countingV[i] <= 0) continue;
        int b = (int)countingV[i];
        bucket_of[i] = b;
        pos_in_bucket[i] = buckets[b].size();
        buckets[b].push_back(i);
        vertexInHeap[i] = 1;
        remainingInHeap++;
    }
    int curBucket = 0;

    // Bucket move helper — called immediately inline
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
    std::vector<uint8_t> leafAffected(tree.adj_list.size(), 0);

    daf::StaticVector<daf::Size> currentRemoveVertexIds(numVertices);

    std::cout << "=========================begin=========================" << std::endl;
    std::cout << "vertices in heap: " << remainingInHeap << std::endl;

    long long minCore = 0;

    while (remainingInHeap > 0) {
        // --- Bucket pop ---
        while (curBucket < (int)buckets.size() && buckets[curBucket].empty()) curBucket++;
        if (curBucket >= (int)buckets.size()) break;

        minCore = std::max((long long)curBucket, minCore);

        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                vertexInHeap[id] = 0;
                currentRemoveVertexIds.push_back(id);
                coreV[id] = (double)minCore;
                remainingInHeap--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (curBucket + 1) <= (int)minCore) {
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
                    leafAffected[leafId] = 1;
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

        // --- Phase 2: merged delta + tree mutation + immediate bucket moves ---
        for (int li = 0; li < (int)removedLeaf.size(); ++li) {
            auto leafId = removedLeaf[li];
            auto &leaf = tree.adj_list[leafId];
            VCD_ST::LeafRmInfo &leafRm = leafRmInfo[leafId];
            if (leaf.empty()) continue;

            // Sort removedPivots — fast path for small arrays
            const auto rpSize = leafRm.removedPivots.size();
            if (rpSize == 2) {
                if (leafRm.removedPivots[0] > leafRm.removedPivots[1])
                    std::swap(leafRm.removedPivots[0], leafRm.removedPivots[1]);
                if (leafRm.removedPivots[0] == leafRm.removedPivots[1])
                    leafRm.removedPivots.c_size = 1;
            } else if (rpSize > 2) {
                std::ranges::sort(leafRm.removedPivots);
                leafRm.removedPivots.unique();
            }

            int numPivots = 0, numKeeps = 0;
            for (auto &node : leaf) {
                if (node.isPivot) numPivots++;
                else numKeeps++;
            }

            int needPivot = (int)k - numKeeps;
            int remainPivots = numPivots - (int)leafRm.removedPivots.size();
            bool leafDies = leafRm.removedKeepC || needPivot < 0 || needPivot > remainPivots;

            if (leafDies) {
                long long keepValue = (needPivot >= 0 && needPivot <= numPivots)
                                   ? std::llround(nCr[numPivots][needPivot]) : 0;
                long long pivotValue = (needPivot >= 1 && needPivot - 1 <= numPivots - 1)
                                    ? std::llround(nCr[numPivots - 1][needPivot - 1]) : 0;

                for (auto &node : leaf) {
                    long long delta = node.isPivot ? pivotValue : keepValue;
                    if (delta > 0) {
                        countingV[node.v] -= delta;
                        if (countingV[node.v] < 0) countingV[node.v] = 0;
                        if (vertexInHeap[node.v])
                            bucketMove(node.v);
                    }
                }
                leaf.clear();
                tree.recycleNode(leafId);
            } else if (!leafRm.removedPivots.empty() && needPivot <= remainPivots) {
                long long KeepDelta = std::llround(nCr[numPivots][needPivot]) - std::llround(nCr[remainPivots][needPivot]);
                long long RemovedPivotFull = (needPivot >= 1) ? std::llround(nCr[numPivots - 1][needPivot - 1]) : 0;
                long long PivotDelta = (needPivot >= 1)
                    ? std::llround(nCr[numPivots - 1][needPivot - 1]) - std::llround(nCr[remainPivots - 1][needPivot - 1])
                    : 0;

                for (auto rp : leafRm.removedPivots) {
                    if (RemovedPivotFull > 0) {
                        countingV[rp] -= RemovedPivotFull;
                        if (countingV[rp] < 0) countingV[rp] = 0;
                        if (vertexInHeap[rp])
                            bucketMove(rp);
                    }
                }

                for (auto &node : leaf) {
                    if (node.isPivot) {
                        // Linear scan for small removedPivots (usually 1-3 elements)
                        bool isRemoved = false;
                        for (daf::Size ri = 0; ri < leafRm.removedPivots.size(); ++ri) {
                            if (leafRm.removedPivots[ri] == node.v) { isRemoved = true; break; }
                            if (leafRm.removedPivots[ri] > node.v) break; // sorted, can early exit
                        }
                        if (isRemoved) continue;
                    }
                    long long delta = node.isPivot ? PivotDelta : KeepDelta;
                    if (delta > 0) {
                        countingV[node.v] -= delta;
                        if (countingV[node.v] < 0) countingV[node.v] = 0;
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
            leafAffected[removedLeaf[li]] = 0;
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
