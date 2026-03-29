//
// ST V13: Path-Fragility Peeling (PFP) — No r-Clique IDs
//
// Core idea: Peel by path fragility instead of r-clique support.
// For isolated paths (~85%+), compute core numbers in O(r) per path
// without enumerating any r-clique. For non-isolated paths, use
// on-demand vertex-path intersection to compute r-clique support.
//
// Eliminates: StaticCliqueIndex, cliqueLeafIds, countingRClique, leafCliqueInfo
// Space: O(V + Σ|P|) — same as CPI itself
//

#include "NCliqueCoreDecomposition.h"
#include "dataStruct/CliqueHashMap.h"
#include <chrono>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <numeric>

#include "../BK/BronKerboschRmRClique.hpp"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace RCliqueSTv13 {

struct VtxPathEntry {
    daf::Size pathId;
    uint8_t isPivot;
};

struct PathInfo {
    daf::CliqueSize pivotC;
    daf::CliqueSize keepC;
    int needPivot;
    int fragility;
    bool alive;
};

static int computeFragility(int pivotC, int needPivot, int r) {
    int minRP = std::min(r, pivotC);
    int remP = pivotC - minRP;
    int remN = needPivot - minRP;
    if (remP < 0 || remN < 0 || remP < remN) return 0;
    return (int)(nCr[remP][remN] + 0.5);
}

static long long localSupport(int pivotC, int needPivot, int b) {
    int remP = pivotC - b;
    int remN = needPivot - b;
    if (remP < 0 || remN < 0 || remP < remN) return 0;
    return (double)(nCr[remP][remN] + 0.5);
}

} // namespace RCliqueSTv13


std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V13(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    using namespace RCliqueSTv13;

    long long duration_init = 0, duration_pop = 0;
    long long duration_isolated = 0, duration_nonisolated = 0, duration_output = 0;
    long long cntIsolatedDeath = 0, cntBK = 0, cntRebucket = 0;

    auto time_start = std::chrono::high_resolution_clock::now();
    const daf::Size numLeaves = tree.adj_list.size();

    // ===== Build vtxPaths + pathInfo =====
    std::vector<std::vector<VtxPathEntry>> vtxPaths;
    std::vector<PathInfo> pathInfo(numLeaves);
    {
        daf::Size maxV = 0;
        for (daf::Size L = 0; L < numLeaves; ++L)
            for (const auto &node : tree.adj_list[L])
                maxV = std::max(maxV, node.v + 1);
        vtxPaths.resize(maxV);

        for (daf::Size L = 0; L < numLeaves; ++L) {
            const auto &leaf = tree.adj_list[L];
            if (leaf.size() < (size_t)r) {
                pathInfo[L] = {0, 0, 0, 0, false};
                continue;
            }
            daf::CliqueSize pC = 0, kC = 0;
            for (const auto &node : leaf) {
                if (node.isPivot) pC++; else kC++;
                vtxPaths[node.v].push_back({L, (uint8_t)node.isPivot});
            }
            int np = s - (int)kC;
            int frag = computeFragility(pC, np, r);
            pathInfo[L] = {pC, kC, np, frag, true};
        }
    }

    // ===== DIAGNOSTIC: Isolation + vertex sharing =====
    {
        daf::Size totalAlive = 0, isolated = 0, nonIsolated = 0;
        // vtxPaths degree distribution: how many paths does each vertex appear in?
        std::map<int, int> vtxDegDist; // degree -> count
        size_t totalVtxPathEntries = 0;
        for (daf::Size v = 0; v < vtxPaths.size(); ++v) {
            if (!vtxPaths[v].empty()) {
                vtxDegDist[(int)vtxPaths[v].size()]++;
                totalVtxPathEntries += vtxPaths[v].size();
            }
        }

        // Per-path: how many shared vertices (appear in >1 path)?
        long long totalSharedVerts = 0;
        int maxShared = 0;
        for (daf::Size L = 0; L < numLeaves; ++L) {
            if (!pathInfo[L].alive) continue;
            totalAlive++;
            bool iso = true;
            int sharedCount = 0;
            for (const auto &node : tree.adj_list[L]) {
                if (vtxPaths[node.v].size() > 1) {
                    iso = false;
                    sharedCount++;
                }
            }
            if (iso) isolated++; else nonIsolated++;
            totalSharedVerts += sharedCount;
            maxShared = std::max(maxShared, sharedCount);
        }
        double pct = totalAlive > 0 ? 100.0 * isolated / totalAlive : 0.0;
        double avgShared = totalAlive > 0 ? (double)totalSharedVerts / totalAlive : 0.0;

        // Keep vertices (appear in ALL paths via hold set) distribution
        size_t keepOnlyVtx = 0, pivotOnlyVtx = 0, bothVtx = 0;
        for (daf::Size v = 0; v < vtxPaths.size(); ++v) {
            if (vtxPaths[v].empty()) continue;
            bool hasKeep = false, hasPivot = false;
            for (const auto &e : vtxPaths[v]) {
                if (e.isPivot) hasPivot = true; else hasKeep = true;
            }
            if (hasKeep && hasPivot) bothVtx++;
            else if (hasKeep) keepOnlyVtx++;
            else pivotOnlyVtx++;
        }

        std::cout << "===== PATH ISOLATION DIAGNOSTIC (r=" << r << ", s=" << s << ") =====" << std::endl;
        std::cout << "  Total alive paths:   " << totalAlive << std::endl;
        std::cout << "  Isolated paths:      " << isolated << " (" << pct << "%)" << std::endl;
        std::cout << "  Non-isolated paths:  " << nonIsolated << std::endl;
        std::cout << "  Avg shared verts/path:" << avgShared << ", max:" << maxShared << std::endl;
        std::cout << "  Vertex degree distribution (degree -> #vertices):" << std::endl;
        for (const auto &[deg, cnt] : vtxDegDist) {
            if (deg <= 10 || cnt > 100)
                std::cout << "    deg=" << deg << " -> " << cnt << " vertices" << std::endl;
        }
        std::cout << "  Total vtxPath entries: " << totalVtxPathEntries << std::endl;
        std::cout << "  Vertex roles: keepOnly=" << keepOnlyVtx << " pivotOnly=" << pivotOnlyVtx << " both=" << bothVtx << std::endl;

        // Leaf-death analysis: how many paths would die at curLevel=0?
        // A path dies if ANY hold vertex has ALL its r-cliques removed.
        // At level 0: fragility=0 means the weakest r-clique has support 0.
        // Count paths with fragility=0 (immediate death candidates)
        daf::Size frag0 = 0, frag1 = 0, fragGt1 = 0;
        for (daf::Size L = 0; L < numLeaves; ++L) {
            if (!pathInfo[L].alive) continue;
            if (pathInfo[L].fragility == 0) frag0++;
            else if (pathInfo[L].fragility == 1) frag1++;
            else fragGt1++;
        }
        std::cout << "  Fragility: frag=0:" << frag0 << " frag=1:" << frag1 << " frag>1:" << fragGt1 << std::endl;
        std::cout << "===== END DIAGNOSTIC =====" << std::endl;
        return {};  // skip peeling, return empty result
    }

    // ===== pathBucket =====
    int maxBucket = 0;
    daf::Size aliveCount = 0;
    for (daf::Size L = 0; L < numLeaves; ++L) {
        if (pathInfo[L].alive) {
            maxBucket = std::max(maxBucket, pathInfo[L].fragility);
            aliveCount++;
        }
    }
    std::vector<std::vector<daf::Size>> pathBucket(maxBucket + 2);
    std::vector<int> path_bucket_of(numLeaves, -1);
    for (daf::Size L = 0; L < numLeaves; ++L) {
        if (pathInfo[L].alive) {
            int b = pathInfo[L].fragility;
            path_bucket_of[L] = b;
            pathBucket[b].push_back(L);
        }
    }

    struct DeadPathRecord {
        std::vector<TreeGraphNode> vertices;
        int coreLevel;
    };
    std::vector<DeadPathRecord> deadPaths;
    deadPaths.reserve(numLeaves);

    std::vector<std::pair<std::vector<daf::Size>, double>> bkOutput;

    daf::log_memory("V13 index (vtxPaths + pathBucket)");
    duration_init = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::cout << "=========================begin (r>=3 ST_V13 PFP)===========================" << std::endl;
    std::cout << "Paths: " << aliveCount << ", maxBucket: " << maxBucket << std::endl;

    // ===== Helpers =====
    auto ensureSize = [&](daf::Size id) {
        if (id >= pathInfo.size()) {
            pathInfo.resize(id + 1);
            path_bucket_of.resize(id + 1, -1);
        }
    };

    auto deactivatePath = [&](daf::Size pathId) {
        pathInfo[pathId].alive = false;
        for (const auto &node : tree.adj_list[pathId]) {
            auto &vp = vtxPaths[node.v];
            for (size_t i = 0; i < vp.size(); ) {
                if (vp[i].pathId == pathId) {
                    vp[i] = vp.back();
                    vp.pop_back();
                } else ++i;
            }
        }
        aliveCount--;
    };

    auto isIsolated = [&](daf::Size pathId) -> bool {
        for (const auto &node : tree.adj_list[pathId])
            for (const auto &entry : vtxPaths[node.v])
                if (entry.pathId != pathId) return false;
        return true;
    };

    auto addNewPath = [&](std::vector<TreeGraphNode> &newLeaf, int floorLevel) -> daf::Size {
        daf::Size newId = tree.addNodePresorted(newLeaf);
        ensureSize(newId);
        const auto &leaf = tree.adj_list[newId];
        daf::CliqueSize pC = 0, kC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pC++; else kC++;
            if (node.v >= vtxPaths.size()) vtxPaths.resize(node.v + 1);
            vtxPaths[node.v].push_back({newId, (uint8_t)node.isPivot});
        }
        int np = s - (int)kC;
        int frag = std::max(computeFragility(pC, np, r), floorLevel);
        pathInfo[newId] = {pC, kC, (int)(s - (int)kC), frag, true};
        if (frag >= (int)pathBucket.size()) pathBucket.resize(frag + 1);
        path_bucket_of[newId] = frag;
        pathBucket[frag].push_back(newId);
        aliveCount++;
        return newId;
    };

    auto rebucketPath = [&](daf::Size pathId, int newLevel) {
        int oldB = path_bucket_of[pathId];
        if (newLevel == oldB) return;
        auto &oldVec = pathBucket[oldB];
        for (size_t i = 0; i < oldVec.size(); ++i) {
            if (oldVec[i] == pathId) {
                oldVec[i] = oldVec.back();
                oldVec.pop_back();
                break;
            }
        }
        if (newLevel >= (int)pathBucket.size()) pathBucket.resize(newLevel + 1);
        path_bucket_of[pathId] = newLevel;
        pathBucket[newLevel].push_back(pathId);
    };

    // On-demand globalSup for one r-clique (returns the actual value)
    auto computeGlobalSup = [&](const daf::StaticVector<TreeGraphNode> &rClique) -> long long {
        // Find vertex with fewest covering paths
        daf::Size minVtx = rClique[0].v;
        size_t minSize = vtxPaths[rClique[0].v].size();
        for (daf::Size j = 1; j < r; ++j) {
            if (vtxPaths[rClique[j].v].size() < minSize) {
                minSize = vtxPaths[rClique[j].v].size();
                minVtx = rClique[j].v;
            }
        }

        long long globalSup = 0;
        for (const auto &entry : vtxPaths[minVtx]) {
            daf::Size qId = entry.pathId;
            if (!pathInfo[qId].alive) continue;

            bool allIn = true;
            int bInQ = 0;
            for (daf::Size j = 0; j < r; ++j) {
                daf::Size v = rClique[j].v;
                if (v == minVtx) {
                    if (entry.isPivot) bInQ++;
                    continue;
                }
                bool found = false;
                for (const auto &e2 : vtxPaths[v]) {
                    if (e2.pathId == qId) {
                        found = true;
                        if (e2.isPivot) bInQ++;
                        break;
                    }
                }
                if (!found) { allIn = false; break; }
            }
            if (!allIn) continue;

            globalSup += localSupport(pathInfo[qId].pivotC,
                                       pathInfo[qId].needPivot, bInQ);
        }
        return globalSup;
    };

    // ===== Main peeling loop =====
    int curLevel = 0;
    long long totalIters = 0;

    while (aliveCount > 0) {
        auto t_pop = std::chrono::high_resolution_clock::now();
        while (curLevel < (int)pathBucket.size() && pathBucket[curLevel].empty())
            curLevel++;
        if (curLevel >= (int)pathBucket.size()) break;

        daf::Size pathId = pathBucket[curLevel].back();
        pathBucket[curLevel].pop_back();

        if (!pathInfo[pathId].alive) continue;

        duration_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::high_resolution_clock::now() - t_pop).count();

        totalIters++;

        if (isIsolated(pathId)) {
            // ===== ISOLATED PATH =====
            auto t_iso = std::chrono::high_resolution_clock::now();
            cntIsolatedDeath++;

            const auto &leaf = tree.adj_list[pathId];
            DeadPathRecord rec;
            rec.vertices.assign(leaf.begin(), leaf.end());
            rec.coreLevel = curLevel;
            deadPaths.push_back(std::move(rec));

            deactivatePath(pathId);
            tree.adj_list[pathId].clear();

            duration_isolated += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - t_iso).count();
        } else {
            // ===== NON-ISOLATED PATH =====
            auto t_noniso = std::chrono::high_resolution_clock::now();
            cntBK++;

            const auto &leaf = tree.adj_list[pathId];
            if (leaf.empty()) continue;
            int n = (int)leaf.size();
            daf::CliqueSize pivotC = pathInfo[pathId].pivotC;
            int needPivot = pathInfo[pathId].needPivot;

            daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
            for (int i = 0; i < n; ++i)
                mapRef[leaf[i].v] = (daf::Size)i;

            // Enumerate ALL r-cliques in this path to find:
            //   1) Candidates with globalSup ≤ curLevel → remove
            //   2) Minimum globalSup across all r-cliques → re-bucket target
            struct CandidateClique {
                std::vector<daf::Size> vertices;
            };
            std::vector<CandidateClique> toRemove;
            long long minGlobalSup = std::numeric_limits<long long>::max();

            daf::enumerateCombinationsWithIdx(leaf, r,
                [&](const daf::StaticVector<TreeGraphNode> &rClique, const size_t*) {
                // First quick filter: localSup on this path
                daf::CliqueSize subP = 0;
                for (daf::Size j = 0; j < r; ++j)
                    if (rClique[j].isPivot) subP++;

                long long ls = localSupport(pivotC, needPivot, subP);

                // For r-cliques with localSup > curLevel, we know
                // globalSup >= localSup > curLevel, so they won't be removed.
                // But we still need their globalSup for re-bucketing.
                // Optimization: use localSup as lower bound for minGlobalSup tracking.
                if (ls > curLevel) {
                    // Don't compute full globalSup, but track localSup as lower bound
                    minGlobalSup = std::min(minGlobalSup, ls);
                    return true;
                }

                // For r-cliques with localSup ≤ curLevel, compute full globalSup
                long long gs = computeGlobalSup(rClique);

                if (gs <= curLevel) {
                    std::vector<daf::Size> verts(r);
                    for (daf::Size j = 0; j < r; ++j) verts[j] = rClique[j].v;
                    std::sort(verts.begin(), verts.end());
                    toRemove.push_back({std::move(verts)});
                } else {
                    minGlobalSup = std::min(minGlobalSup, gs);
                }
                return true;
            });

            if (toRemove.empty()) {
                // No r-cliques qualify — skip-level re-bucket
                cntRebucket++;
                int newLevel = (minGlobalSup == std::numeric_limits<long long>::max())
                    ? curLevel + 1
                    : (int)minGlobalSup;
                newLevel = std::max(newLevel, curLevel + 1);
                if (newLevel >= (int)pathBucket.size())
                    pathBucket.resize(newLevel + 1);
                path_bucket_of[pathId] = newLevel;
                pathBucket[newLevel].push_back(pathId);

                duration_nonisolated += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - t_noniso).count();
                continue;
            }

            // Record removed r-cliques
            for (auto &cand : toRemove)
                bkOutput.emplace_back(std::move(cand.vertices), curLevel);

            // Check leaf death
            daf::Size maxRCliquePerVertex = (daf::Size)((double)(nCr[n - 1][r - 1] + 0.5));
            std::vector<daf::Size> vertexConflictDeg(n, 0);
            size_t bkStart = bkOutput.size() - toRemove.size();
            for (size_t ci = bkStart; ci < bkOutput.size(); ++ci)
                for (auto v : bkOutput[ci].first) {
                    daf::Size pos = mapRef[v];
                    if (pos < (daf::Size)n)
                        vertexConflictDeg[pos]++;
                }

            bool leafDead = false;
            int forcedPivotRemove = 0;
            for (int i = 0; i < n; ++i) {
                if (vertexConflictDeg[i] >= maxRCliquePerVertex) {
                    if (!leaf[i].isPivot) { leafDead = true; break; }
                    forcedPivotRemove++;
                }
            }
            if (!leafDead) {
                int remainingPivots = (int)pivotC - forcedPivotRemove;
                int remainingTotal = n - forcedPivotRemove;
                if (remainingTotal < (int)s || remainingPivots < needPivot)
                    leafDead = true;
            }

            if (leafDead) {
                // Path dies — collect neighbors, record, deactivate
                std::vector<daf::Size> neighborPaths;
                neighborPaths.reserve(leaf.size() * 4);
                for (const auto &node : leaf)
                    for (const auto &entry : vtxPaths[node.v])
                        if (entry.pathId != pathId && pathInfo[entry.pathId].alive)
                            neighborPaths.push_back(entry.pathId);
                std::sort(neighborPaths.begin(), neighborPaths.end());
                neighborPaths.erase(std::unique(neighborPaths.begin(), neighborPaths.end()),
                                    neighborPaths.end());

                DeadPathRecord rec;
                rec.vertices.assign(leaf.begin(), leaf.end());
                rec.coreLevel = curLevel;
                deadPaths.push_back(std::move(rec));

                deactivatePath(pathId);
                tree.adj_list[pathId].clear();

                for (auto npId : neighborPaths) {
                    if (!pathInfo[npId].alive) continue;
                    int oldB = path_bucket_of[npId];
                    if (oldB > curLevel)
                        rebucketPath(npId, curLevel);
                }
            } else {
                // BK split
                std::vector<std::vector<daf::Size>> conflictSets;
                conflictSets.reserve(toRemove.size());
                for (size_t ci = bkStart; ci < bkOutput.size(); ++ci)
                    conflictSets.push_back(bkOutput[ci].first);

                auto leafCopy = tree.adj_list[pathId];
                deactivatePath(pathId);
                tree.adj_list[pathId].clear();

                std::vector<TreeGraphNode> reusableLeaf;
                reusableLeaf.reserve(400);

                bkRmClique::removeRClique(leafCopy, conflictSets, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &bkPivots) {
                        reusableLeaf.clear();
                        c.for_each_bit([&](size_t parentPos) {
                            bool isP = bkPivots.test(parentPos);
                            reusableLeaf.push_back({leafCopy[parentPos].v, isP});
                        });
                        if ((int)reusableLeaf.size() >= (int)s)
                            addNewPath(reusableLeaf, curLevel);
                    });
            }

            duration_nonisolated += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - t_noniso).count();
        }
    }

    // ===== OUTPUT PHASE =====
    auto t_out = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<std::vector<daf::Size>, double>> sortedK;
    size_t estTotal = bkOutput.size();
    for (const auto &dp : deadPaths) {
        int nn = (int)dp.vertices.size();
        if (nn >= r) estTotal += (size_t)(nCr[nn][r] + 0.5);
    }
    sortedK.reserve(estTotal);

    for (auto &entry : bkOutput)
        sortedK.push_back(std::move(entry));

    for (const auto &dp : deadPaths) {
        int nn = (int)dp.vertices.size();
        if (nn < r) continue;
        int pC = 0, kC = 0;
        for (const auto &node : dp.vertices) {
            if (node.isPivot) pC++; else kC++;
        }
        int np = s - kC;

        daf::enumerateCombinationsWithIdx(dp.vertices, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique, const size_t*) {
            daf::CliqueSize subP = 0;
            std::vector<daf::Size> verts(r);
            for (daf::Size j = 0; j < r; ++j) {
                verts[j] = rClique[j].v;
                if (rClique[j].isPivot) subP++;
            }
            std::sort(verts.begin(), verts.end());
            long long ls = localSupport(pC, np, subP);
            int core = std::max((int)ls, dp.coreLevel);
            sortedK.emplace_back(std::move(verts), core);
            return true;
        });
    }

    // Deduplicate: keep highest core per clique
    std::sort(sortedK.begin(), sortedK.end(),
              [](const auto &a, const auto &b) {
                  if (a.first != b.first) return a.first < b.first;
                  return a.second < b.second;
              });
    {
        std::vector<std::pair<std::vector<daf::Size>, double>> deduped;
        deduped.reserve(sortedK.size());
        for (size_t i = 0; i < sortedK.size(); ++i) {
            if (i + 1 < sortedK.size() && sortedK[i].first == sortedK[i + 1].first)
                continue;
            deduped.push_back(std::move(sortedK[i]));
        }
        sortedK = std::move(deduped);
    }

    std::sort(sortedK.begin(), sortedK.end(),
              [](const auto &a, const auto &b) { return a.second < b.second; });

    duration_output = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::high_resolution_clock::now() - t_out).count();

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - time_start).count();

    std::cout << "time: " << elapsed << " ms" << std::endl;
    std::cout << "Time Breakdown (ms):" << std::endl;
    std::cout << "  Init:         " << duration_init / 1000000.0 << std::endl;
    std::cout << "  Pop:          " << duration_pop / 1000000.0 << std::endl;
    std::cout << "  Isolated:     " << duration_isolated / 1000000.0 << std::endl;
    std::cout << "  NonIsolated:  " << duration_nonisolated / 1000000.0 << std::endl;
    std::cout << "  Output:       " << duration_output / 1000000.0 << std::endl;
    std::cout << "  Cases: IsolatedDeath=" << cntIsolatedDeath
              << " BK=" << cntBK
              << " Rebucket=" << cntRebucket
              << " iters=" << totalIters << std::endl;
    std::cout << "  Total r-cliques output: " << sortedK.size() << std::endl;

    return sortedK;
}
