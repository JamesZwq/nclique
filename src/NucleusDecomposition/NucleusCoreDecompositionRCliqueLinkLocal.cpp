//
// Framework 2: Exact Local H-Index via LeafCliqueInfo Containment
//
// Correct exact H-index: for each r-clique C, for each s-clique S containing C,
// compute q(S,C) = min{core(D) : D ⊂ S, D ≠ C, |D|=r}. Then
// H-index(C) = largest k such that |{S : q(S,C) ≥ k}| ≥ k.
//
// Key optimization over CPI_Exact:
//   Instead of enumerating C(s,r) r-subsets per s-clique and doing hash lookups,
//   iterate leafCliqueInfo entries and check containment via positional bitset.
//   This replaces ~C(s,r) hash lookups (~30-60 cycles each) with
//   ~|leafCliqueInfo| × r position tests (~2 cycles each).
//
// Tree is IMMUTABLE. No BK, no tree mutation.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <vector>
#include <cstring>
#include <queue>
#include <cmath>

#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace LinkLocal {

struct LeafCliqueEntry {
    daf::Size cliqueId;
    double ncrValue;
    uint8_t positions[8]; // positions of this r-clique's vertices in the leaf (max r=8)
    uint8_t numPos;       // = r
};

struct DualIndex {
    std::vector<std::vector<LeafCliqueEntry>> leafCliqueInfo;
    std::vector<std::vector<daf::Size>> cliqueLeafIds;
    std::vector<double> counting;
};

DualIndex buildDualIndex(
    const DynamicGraph<TreeGraphNode> &treeGraph,
    const StaticCliqueIndex &cliqueIndex,
    daf::CliqueSize r, daf::CliqueSize s) {

    const daf::Size nClique = cliqueIndex.size();
    const daf::Size numLeaves = treeGraph.adj_list.size();

    DualIndex result;
    result.leafCliqueInfo.resize(numLeaves);
    result.cliqueLeafIds.resize(nClique);
    result.counting.assign(nClique, 0);

    daf::StaticVector<daf::Size> &mapRef = daf::vListMap;

    for (daf::Size leafId = 0; leafId < numLeaves; ++leafId) {
        const auto &leaf = treeGraph.adj_list[leafId];
        if (leaf.size() < r) continue;

        daf::CliqueSize pivotC = 0, keepC = 0;
        for (const auto &i : leaf) {
            if (i.isPivot) pivotC++;
            else keepC++;
        }
        int needPivot = s - static_cast<int>(keepC);

        // Set up position map for this leaf
        int n = (int)leaf.size();
        for (int i = 0; i < n; ++i)
            mapRef[leaf[i].v] = (daf::Size)i;

        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            daf::CliqueSize subP = 0;
            for (const auto &node : rClique) if (node.isPivot) subP++;
            double ncrValue = nCr[pivotC - subP][needPivot - subP];
            auto id = cliqueIndex.byClique(rClique);
            if (id < nClique) {
                result.counting[id] += ncrValue;

                LeafCliqueEntry entry;
                entry.cliqueId = id;
                entry.ncrValue = ncrValue;
                entry.numPos = (uint8_t)r;
                for (int j = 0; j < (int)r; ++j)
                    entry.positions[j] = (uint8_t)mapRef[rClique[j].v];

                result.leafCliqueInfo[leafId].push_back(entry);
                result.cliqueLeafIds[id].push_back(leafId);
            }
            return true;
        });
    }
    return result;
}

} // namespace LinkLocal


// ============================================================
// Compute EXACT H-index for r-clique C using leafCliqueInfo containment.
//
// For each leaf containing C:
//   For each s-clique S from that leaf containing C:
//     For each OTHER r-clique D in leafCliqueInfo[leaf]:
//       Check if D ⊂ S using positional membership array
//       Track min core among contained D's
//     bucket[minCore] += 1
//
// H-index = largest k where accumulated(bucket[k..]) >= k
//
// Uses uint8_t membership array instead of uint64_t bitset to support
// leaves with up to 256 vertices (com-dblp has leaves up to 114).
// ============================================================
static long long computeHIndexExactViaContainment(
    daf::Size cliqueId,
    const StaticCliqueIndex &cliqueIndex,
    const DynamicGraph<TreeGraphNode> &tree,
    const long long *coreC,
    const std::vector<std::vector<LinkLocal::LeafCliqueEntry>> &leafCliqueInfo,
    const std::vector<std::vector<daf::Size>> &cliqueLeafIds,
    daf::CliqueSize r, daf::CliqueSize s,
    long long currentCore,
    std::vector<long long> &buckets,
    // per-leaf metadata
    const std::vector<int> &leafPivotC,
    const std::vector<int> &leafNeedPivot,
    const std::vector<uint8_t> &leafValid,
    // scratch buffers (caller-owned, avoid repeated allocation)
    std::vector<uint8_t> &inSClique,
    std::vector<uint16_t> &otherPivotPos)
{
    if (currentCore <= 0) return 0;

    int bucketSize = (int)currentCore;
    if (bucketSize < 1) bucketSize = 1;
    if ((int)buckets.size() < bucketSize + 1)
        buckets.resize(bucketSize + 1, 0);
    for (int i = 0; i <= bucketSize; ++i) buckets[i] = 0;

    if (cliqueId >= cliqueLeafIds.size()) return 0;

    auto cliqueVerts = cliqueIndex.byId(cliqueId);
    daf::StaticVector<daf::Size> &mapRef = daf::vListMap;

    long long rawTotalSupport = 0;

    for (auto leafId : cliqueLeafIds[cliqueId]) {
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.empty()) continue;
        if (!leafValid[leafId]) continue;

        int n = (int)leaf.size();
        int needPivot = leafNeedPivot[leafId];

        // Set up position map for this leaf
        for (int i = 0; i < n; ++i)
            mapRef[leaf[i].v] = (daf::Size)i;

        // Find C's pivots in this leaf
        int cPivotsInLeaf = 0;
        for (auto v : cliqueVerts) {
            if (leaf[mapRef[v]].isPivot) cPivotsInLeaf++;
        }

        // Mark C's positions in inSClique (base membership)
        if ((int)inSClique.size() < n) inSClique.resize(n, 0);
        // Clear relevant positions
        for (int i = 0; i < n; ++i) inSClique[i] = 0;

        // Base membership: C's vertices + all other keeps
        for (auto v : cliqueVerts)
            inSClique[mapRef[v]] = 1;

        otherPivotPos.clear();
        int numOtherKeeps = 0;
        for (int i = 0; i < n; ++i) {
            if (inSClique[i]) continue; // C's vertex
            if (leaf[i].isPivot) {
                otherPivotPos.push_back((uint16_t)i);
            } else {
                inSClique[i] = 1; // other keep — always in s-clique
                numOtherKeeps++;
            }
        }
        int numOtherPivots = (int)otherPivotPos.size();

        int effectiveNeedPivot = needPivot - cPivotsInLeaf;
        if (effectiveNeedPivot < 0 || effectiveNeedPivot > numOtherPivots) continue;
        int need = (int)s - (int)r;
        if (numOtherKeeps + numOtherPivots < need) continue;
        if (numOtherKeeps > need) continue;

        const auto &leafEntries = leafCliqueInfo[leafId];

        // Lambda to process one s-clique. inSClique[pos]==1 for all positions in S.
        auto processSClique = [&]() {
            long long minOtherCore = std::numeric_limits<long long>::max();
            bool hasOther = false;

            for (const auto &entry : leafEntries) {
                if (entry.cliqueId == cliqueId) continue;

                // Check if ALL of entry's vertices are in inSClique
                bool allIn = true;
                for (int j = 0; j < entry.numPos; ++j) {
                    if (!inSClique[entry.positions[j]]) {
                        allIn = false;
                        break;
                    }
                }
                if (!allIn) continue;

                long long c = coreC[entry.cliqueId];
                if (c < minOtherCore) minOtherCore = c;
                hasOther = true;
                if (minOtherCore <= 0) break;
            }

            if (!hasOther || minOtherCore <= 0) return;

            int level = (minOtherCore > bucketSize) ? bucketSize : (int)minOtherCore;
            buckets[level] += 1;
            rawTotalSupport += 1;
        };

        if (effectiveNeedPivot == 0) {
            processSClique();
        } else if (effectiveNeedPivot == 1) {
            for (int p = 0; p < numOtherPivots; ++p) {
                inSClique[otherPivotPos[p]] = 1;
                processSClique();
                inSClique[otherPivotPos[p]] = 0;
            }
        } else if (effectiveNeedPivot == 2) {
            for (int p1 = 0; p1 < numOtherPivots; ++p1) {
                inSClique[otherPivotPos[p1]] = 1;
                for (int p2 = p1 + 1; p2 < numOtherPivots; ++p2) {
                    inSClique[otherPivotPos[p2]] = 1;
                    processSClique();
                    inSClique[otherPivotPos[p2]] = 0;
                }
                inSClique[otherPivotPos[p1]] = 0;
            }
        } else {
            // General case: recursive enumeration with set/unset on inSClique
            auto enumPivots = [&](auto &&self, int start, int remaining) -> void {
                if (remaining == 0) {
                    processSClique();
                    return;
                }
                for (int p = start; p <= numOtherPivots - remaining; ++p) {
                    inSClique[otherPivotPos[p]] = 1;
                    self(self, p + 1, remaining - 1);
                    inSClique[otherPivotPos[p]] = 0;
                }
            };
            enumPivots(enumPivots, 0, effectiveNeedPivot);
        }

        // Clean up: unmark C's vertices and other keeps from inSClique
        // (already all 0 for pivot positions since we unset them)
        for (auto v : cliqueVerts) inSClique[mapRef[v]] = 0;
        for (int i = 0; i < n; ++i) inSClique[i] = 0;
    }

    if (rawTotalSupport < 1) return 0;

    // Bucket H-index scan
    long long accumulated = 0;
    for (int c = bucketSize; c >= 1; --c) {
        accumulated += buckets[c];
        if (accumulated >= c) return c;
    }
    return 0;
}


// ============================================================
// Main entry: Link-Graph Local H-Index for r≥3 (corrected)
// ============================================================
std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLinkLocal(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    auto time_start = std::chrono::high_resolution_clock::now();

    // Build clique index
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build (LinkLocal)", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });
    const daf::Size numCliques = cliqueIndex.size();
    const daf::Size numLeaves = tree.adj_list.size();

    // Build dual index (with positional encoding)
    auto dualIndex = daf::timeCount("buildDualIndex (LinkLocal)", [&]() {
        return LinkLocal::buildDualIndex(tree, cliqueIndex, r, s);
    });
    auto &leafCliqueInfo = dualIndex.leafCliqueInfo;
    auto &cliqueLeafIds = dualIndex.cliqueLeafIds;

    // Precompute per-leaf metadata
    std::vector<int> leafPivotC(numLeaves), leafNeedPivot(numLeaves);
    std::vector<uint8_t> leafValid(numLeaves);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        int keeps = 0, pivots = 0;
        for (const auto &node : tree.adj_list[L]) {
            if (node.isPivot) pivots++;
            else keeps++;
        }
        leafPivotC[L] = pivots;
        leafNeedPivot[L] = (int)s - keeps;
        leafValid[L] = (leafNeedPivot[L] >= 0 && leafNeedPivot[L] <= pivots) ? 1 : 0;
    }

    // Initial support → core estimates
    auto *coreC = new long long[numCliques];
    for (daf::Size i = 0; i < numCliques; ++i)
        coreC[i] = dualIndex.counting[i];

    // Scratch
    std::vector<long long> buckets;
    buckets.reserve(4096);
    std::vector<uint8_t> inSClique;
    inSClique.reserve(256);
    std::vector<uint16_t> otherPivotPosScratch;
    otherPivotPosScratch.reserve(256);

    // Dirty queue
    std::vector<uint8_t> inQueue(numCliques, 0);
    std::queue<daf::Size> dirtyQueue;

    auto initEnd = std::chrono::high_resolution_clock::now();
    long long initMs = std::chrono::duration_cast<std::chrono::milliseconds>(initEnd - time_start).count();

    std::cout << "=========================begin (Link-Graph Local H-index r>=" << (int)r << ")=========================" << std::endl;
    std::cout << "r-cliques: " << numCliques << ", leaves: " << numLeaves << std::endl;
    std::cout << "Init time: " << initMs << " ms" << std::endl;

    // Phase 1: initial full scan
    for (daf::Size cid = 0; cid < numCliques; ++cid) {
        if (coreC[cid] <= 0) continue;

        long long newCore = computeHIndexExactViaContainment(
            cid, cliqueIndex, tree, coreC,
            leafCliqueInfo, cliqueLeafIds, r, s,
            coreC[cid], buckets, leafPivotC, leafNeedPivot, leafValid,
            inSClique, otherPivotPosScratch);
        newCore = std::min(newCore, coreC[cid]);

        if (newCore != coreC[cid]) {
            long long oldCore = coreC[cid];
            coreC[cid] = newCore;
            // Enqueue co-leaf r-cliques via leafCliqueInfo
            for (auto leafId : cliqueLeafIds[cid]) {
                if (!leafValid[leafId]) continue;
                for (const auto &entry : leafCliqueInfo[leafId]) {
                    if (entry.cliqueId == cid) continue;
                    if (coreC[entry.cliqueId] > 0 && oldCore >= coreC[entry.cliqueId] && !inQueue[entry.cliqueId]) {
                        inQueue[entry.cliqueId] = 1;
                        dirtyQueue.push(entry.cliqueId);
                    }
                }
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
            daf::Size cid = dirtyQueue.front();
            dirtyQueue.pop();
            inQueue[cid] = 0;
            recomputeCount++;

            if (coreC[cid] <= 0) continue;

            long long newCore = computeHIndexExactViaContainment(
                cid, cliqueIndex, tree, coreC,
                leafCliqueInfo, cliqueLeafIds, r, s,
                coreC[cid], buckets, leafPivotC, leafNeedPivot, leafValid,
                inSClique, otherPivotPosScratch);
            newCore = std::min(newCore, coreC[cid]);

            if (newCore != coreC[cid]) {
                long long oldCore = coreC[cid];
                coreC[cid] = newCore;
                for (auto leafId : cliqueLeafIds[cid]) {
                    if (!leafValid[leafId]) continue;
                    for (const auto &entry : leafCliqueInfo[leafId]) {
                        if (entry.cliqueId == cid) continue;
                        if (coreC[entry.cliqueId] > 0 && oldCore >= coreC[entry.cliqueId] && !inQueue[entry.cliqueId]) {
                            inQueue[entry.cliqueId] = 1;
                            dirtyQueue.push(entry.cliqueId);
                        }
                    }
                }
            }
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();
    std::cout << "Link-Graph Local r>=" << (int)r << " converged in " << iteration
              << " iterations, " << recomputeCount << " re-evaluations" << std::endl;
    std::cout << "time: " << elapsed << " ms" << std::endl;

    // Build output
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    result.reserve(numCliques);
    for (daf::Size i = 0; i < numCliques; ++i) {
        auto verts = cliqueIndex.byId(i);
        std::vector<daf::Size> copy(verts.begin(), verts.end());
        result.emplace_back(std::move(copy), (int)coreC[i]);
    }

    delete[] coreC;
    return result;
}
