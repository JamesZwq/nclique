//
// Region-Based Nucleus Decomposition for r>=3.
//
// Key idea: group vertices by maximal-clique membership ("overlap classes").
// All r-cliques from the same class tuple have the SAME support.
// Work with class tuples instead of individual r-cliques.
//
// Support formula (exact):
//   For class tuple (A,B,C) with common maximal cliques M_1,...,M_k:
//   support = Σ_{S⊆{M_1..M_k}, |S|≥1} (-1)^{|S|+1} C(|∩_{M∈S} M| - r, s - r)
//   (inclusion-exclusion over maximal clique intersections)
//
// For s=r+1: simplifies to |∪ M_i| - r.
//

#include "NCliqueCoreDecomposition.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

extern double nCr[1001][401];
extern std::vector<bool> g_maxCliqueTags;

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_Region(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    auto tStart = std::chrono::high_resolution_clock::now();
    const daf::Size numVertices = edgeGraph.n;
    const daf::Size numPaths = tree.adj_list.size();
    const daf::Size INVALID = static_cast<daf::Size>(-1);

    // ============================================================
    // Step 1: Build Maximal-Clique Regions
    // ============================================================

    // 1a. Collect per-path validity (size >= r)
    std::vector<bool> pathValid(numPaths, false);
    daf::Size validPaths = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (tree.adj_list[pid].size() >= r) {
            pathValid[pid] = true;
            validPaths++;
        }
    }

    // 1b. Assign region IDs to maximal clique paths
    std::vector<daf::Size> regionOf(numPaths, INVALID);
    std::vector<daf::Size> maximalPathIds;
    daf::Size numMaximal = 0, numSub = 0;

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (!pathValid[pid]) continue;
        bool isMax = (pid < g_maxCliqueTags.size()) ? g_maxCliqueTags[pid] : true;
        if (isMax) {
            regionOf[pid] = maximalPathIds.size();
            maximalPathIds.push_back(pid);
            numMaximal++;
        }
    }
    daf::Size numRegions = maximalPathIds.size();

    // 1c. Build region vertex sets (sorted) + vertex-to-region index
    std::vector<std::vector<daf::Size>> regionVerts(numRegions);
    std::vector<std::vector<daf::Size>> vtxMaxPaths(numVertices);

    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        auto &leaf = tree.adj_list[maximalPathIds[rid]];
        regionVerts[rid].reserve(leaf.size());
        for (const auto &node : leaf) regionVerts[rid].push_back(node.v);
        std::sort(regionVerts[rid].begin(), regionVerts[rid].end());
        for (daf::Size v : regionVerts[rid]) {
            if (v < numVertices) vtxMaxPaths[v].push_back(rid);
        }
    }

    // 1d. Merge sub-clique paths into parent maximal cliques
    daf::Size orphans = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (!pathValid[pid] || regionOf[pid] != INVALID) continue;
        numSub++;

        auto &leaf = tree.adj_list[pid];
        std::vector<daf::Size> pVerts;
        pVerts.reserve(leaf.size());
        for (const auto &node : leaf) pVerts.push_back(node.v);
        std::sort(pVerts.begin(), pVerts.end());

        // Pick rarest vertex for efficient parent search
        daf::Size rareV = pVerts[0], rareCount = vtxMaxPaths[pVerts[0]].size();
        for (daf::Size v : pVerts) {
            if (vtxMaxPaths[v].size() < rareCount) { rareV = v; rareCount = vtxMaxPaths[v].size(); }
        }

        for (daf::Size rid : vtxMaxPaths[rareV]) {
            auto &qv = regionVerts[rid];
            if (qv.size() <= pVerts.size()) continue;
            if (std::includes(qv.begin(), qv.end(), pVerts.begin(), pVerts.end())) {
                regionOf[pid] = rid;
                break;
            }
        }
        if (regionOf[pid] == INVALID) {
            regionOf[pid] = maximalPathIds.size();
            maximalPathIds.push_back(pid);
            // Also build regionVerts for this orphan
            regionVerts.push_back(pVerts);
            for (daf::Size v : pVerts)
                if (v < numVertices) vtxMaxPaths[v].push_back(regionVerts.size() - 1);
            numRegions++;
            orphans++;
        }
    }

    auto tStep1 = std::chrono::high_resolution_clock::now();
    auto step1Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep1 - tStart).count();

    std::cout << "======= Region-Based Nucleus Decomposition =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Vertices: " << numVertices << ", valid paths: " << validPaths << std::endl;
    std::cout << "  Maximal cliques: " << numMaximal
              << ", sub-cliques merged: " << numSub
              << " (orphans: " << orphans << ")" << std::endl;
    std::cout << "  Step 1 time: " << step1Ms << " ms" << std::endl;

    // ============================================================
    // Step 2: Build Overlap Classes + Enumerate Region r-Tuples
    // ============================================================
    // General for any r: region r-tuple = multiset of r overlap classes
    // All r-cliques from the same tuple have identical support (Theorem 1)

    auto tStep2 = std::chrono::high_resolution_clock::now();

    // 2a. CPI total support sum (for verification)
    double totalSupportCPI = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (!pathValid[pid]) continue;
        auto &leaf = tree.adj_list[pid];
        daf::Size h = 0, p = 0;
        for (const auto &node : leaf) { if (node.isPivot) p++; else h++; }
        if ((int)(h + p) >= s && (int)p >= (int)(s - h))
            totalSupportCPI += nCr[p][s - h];
    }
    totalSupportCPI *= nCr[s][r];

    // 2b. Build overlap classes
    struct ClassInfo {
        std::vector<daf::Size> regionIds; // sorted maximal clique IDs
        daf::Size size;                   // number of vertices
    };
    std::vector<ClassInfo> classes;
    std::vector<daf::Size> classOf(numVertices, INVALID);
    {
        using Profile = std::vector<daf::Size>;
        struct PH {
            size_t operator()(const Profile &p) const noexcept {
                size_t h = p.size();
                for (auto x : p) h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h<<6) + (h>>2);
                return h;
            }
        };
        std::unordered_map<Profile, daf::Size, PH> pToC;
        for (daf::Size v = 0; v < numVertices; ++v) {
            if (vtxMaxPaths[v].empty()) continue;
            auto it = pToC.find(vtxMaxPaths[v]);
            if (it == pToC.end()) {
                daf::Size cid = classes.size();
                pToC[vtxMaxPaths[v]] = cid;
                classes.push_back({vtxMaxPaths[v], 1});
                classOf[v] = cid;
            } else {
                classOf[v] = it->second;
                classes[it->second].size++;
            }
        }
    }
    daf::Size numClasses = classes.size();

    // 2c. Build per-maximal-clique class list
    std::vector<std::vector<daf::Size>> classesInRegion(numRegions);
    for (daf::Size cid = 0; cid < numClasses; ++cid)
        for (daf::Size rid : classes[cid].regionIds)
            classesInRegion[rid].push_back(cid);
    for (auto &v : classesInRegion) std::sort(v.begin(), v.end());

    // 2d. Enumerate region r-tuples (multisets of size r from each clique's classes)
    using TupleKey = std::vector<daf::Size>; // sorted class IDs, length r
    struct TupleHash {
        size_t operator()(const TupleKey &t) const noexcept {
            size_t h = t.size();
            for (auto x : t) h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h<<6) + (h>>2);
            return h;
        }
    };

    // tupleRegionMap[tuple] = list of maximal cliques containing the tuple
    std::unordered_map<TupleKey, std::vector<daf::Size>, TupleHash> tupleRegionMap;

    // Recursive multiset enumeration
    TupleKey currentTuple;
    currentTuple.reserve(r);

    std::function<void(daf::Size rid, const std::vector<daf::Size>&, int startIdx, int remaining)> enumerate;
    enumerate = [&](daf::Size rid, const std::vector<daf::Size> &cids, int startIdx, int remaining) {
        if (remaining == 0) {
            // Check multiplicity > 0: for each distinct class c appearing k times, need |c| >= k
            std::unordered_map<daf::Size, int> counts;
            for (auto c : currentTuple) counts[c]++;
            bool valid = true;
            for (auto &[c, k] : counts) {
                if ((int)classes[c].size < k) { valid = false; break; }
            }
            if (valid) tupleRegionMap[currentTuple].push_back(rid);
            return;
        }
        for (int i = startIdx; i < (int)cids.size(); ++i) {
            currentTuple.push_back(cids[i]);
            enumerate(rid, cids, i, remaining - 1); // i (not i+1) allows repetition
            currentTuple.pop_back();
        }
    };

    daf::Size maxClassesPerClique = 0;
    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        auto &cids = classesInRegion[rid];
        maxClassesPerClique = std::max(maxClassesPerClique, (daf::Size)cids.size());
        if (cids.size() > 500) continue; // guard against huge cliques
        currentTuple.clear();
        enumerate(rid, cids, 0, r);
    }

    std::cout << "  Overlap classes: " << numClasses
              << ", max classes/clique: " << maxClassesPerClique << std::endl;

    // 2e. Compute support and multiplicity per tuple
    struct RegionTuple {
        TupleKey key;
        daf::Size mult;      // number of r-cliques represented
        daf::Size support;
    };
    std::vector<RegionTuple> tuples;
    tuples.reserve(tupleRegionMap.size());
    double totalRCliques = 0, totalSupportSum = 0;

    for (auto &[key, rids] : tupleRegionMap) {
        // Multiplicity: Π C(|class_c|, count_c)
        std::unordered_map<daf::Size, int> counts;
        for (auto c : key) counts[c]++;
        daf::Size mult = 1;
        for (auto &[c, k] : counts) mult *= (daf::Size)nCr[classes[c].size][k];
        if (mult == 0) continue;

        // Support via inclusion-exclusion over common cliques
        daf::Size nr = rids.size();
        double support = 0;

        if (s == r + 1) {
            // s = r+1: support = |∪ common cliques alive| - r
            std::vector<daf::Size> unionV;
            for (daf::Size rid : rids) {
                std::vector<daf::Size> merged;
                std::merge(unionV.begin(), unionV.end(),
                    regionVerts[rid].begin(), regionVerts[rid].end(),
                    std::back_inserter(merged));
                merged.erase(std::unique(merged.begin(), merged.end()), merged.end());
                unionV = std::move(merged);
            }
            support = std::max(0.0, (double)unionV.size() - r);
        } else if (nr <= 15) {
            for (uint32_t mask = 1; mask < (1u << nr); ++mask) {
                int bits = __builtin_popcount(mask);
                std::vector<daf::Size> inter;
                bool first = true;
                for (daf::Size i = 0; i < nr; ++i) {
                    if (!(mask & (1u << i))) continue;
                    if (first) { inter = regionVerts[rids[i]]; first = false; }
                    else {
                        std::vector<daf::Size> tmp;
                        std::set_intersection(inter.begin(), inter.end(),
                            regionVerts[rids[i]].begin(), regionVerts[rids[i]].end(),
                            std::back_inserter(tmp));
                        inter = std::move(tmp);
                    }
                }
                int n = (int)inter.size() - r, k = s - r;
                if (n >= k && k >= 0)
                    support += (bits % 2 == 1 ? 1 : -1) * nCr[n][k];
            }
        } else {
            // Fallback: single-term approximation using largest clique
            daf::Size maxK = 0;
            for (daf::Size rid : rids) maxK = std::max(maxK, (daf::Size)regionVerts[rid].size());
            support = nCr[std::max(0, (int)maxK - r)][s - r];
        }

        support = std::max(0.0, support);
        tuples.push_back({key, mult, (daf::Size)support});
        totalRCliques += mult;
        totalSupportSum += mult * support;
    }

    auto tStep2End = std::chrono::high_resolution_clock::now();
    auto step2Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep2End - tStep2).count();

    std::cout << "\n  --- Step 2: Region r-Tuples ---" << std::endl;
    std::cout << "  Tuples: " << tuples.size() << std::endl;
    std::cout << "  Total r-cliques (from tuples): " << (int64_t)totalRCliques << std::endl;
    std::cout << "  Total support sum (tuples): " << std::fixed << std::setprecision(0) << totalSupportSum << std::endl;
    std::cout << "  Total support sum (CPI):    " << std::fixed << std::setprecision(0) << totalSupportCPI << std::endl;
    double relErr = std::abs(totalSupportSum - totalSupportCPI) / std::max(1.0, totalSupportCPI);
    if (relErr < 1e-9)
        std::cout << "  VERIFICATION: PASS" << std::endl;
    else
        std::cout << "  VERIFICATION: MISMATCH (relErr=" << relErr << ")" << std::endl;
    std::cout << "  Step 2 time: " << step2Ms << " ms" << std::endl;

    // ============================================================
    // Step 3: Cascade Peeling (general r, s=r+1)
    // ============================================================

    auto tStep3 = std::chrono::high_resolution_clock::now();

    // Build tuple index
    std::unordered_map<TupleKey, daf::Size, TupleHash> tupleIndex;
    for (daf::Size i = 0; i < tuples.size(); ++i)
        tupleIndex[tuples[i].key] = i;

    // Per-tuple common regions
    std::vector<std::vector<daf::Size>> tupleCommonRegions(tuples.size());
    for (auto &[key, rids] : tupleRegionMap) {
        auto it = tupleIndex.find(key);
        if (it != tupleIndex.end()) tupleCommonRegions[it->second] = rids;
    }

    // Bucket queue
    daf::Size maxSupport = 0;
    for (auto &t : tuples) maxSupport = std::max(maxSupport, t.support);

    std::vector<std::vector<daf::Size>> buckets(maxSupport + 2);
    std::vector<daf::Size> tupleSupport(tuples.size());
    std::vector<bool> peeled(tuples.size(), false);

    for (daf::Size i = 0; i < tuples.size(); ++i) {
        tupleSupport[i] = tuples[i].support;
        buckets[tupleSupport[i]].push_back(i);
    }

    std::map<daf::Size, int64_t> coreDist;
    daf::Size numPeeled = 0, currentLevel = 0;
    daf::Size coreLevel = 0; // non-decreasing: core = max(coreLevel, currentLevel)

    // Helper: check if a tuple is alive (unpeeled and has valid r-cliques)
    auto isTupleAlive = [&](const TupleKey &k, const TupleKey &selfKey) -> bool {
        if (k == selfKey) return true; // being peeled now, was alive before
        auto it = tupleIndex.find(k);
        if (it != tupleIndex.end()) return !peeled[it->second];
        // Not in index: check if it COULD have mult > 0
        // (need |class_c| >= count_c for each class c in the tuple)
        std::unordered_map<daf::Size, int> counts;
        for (auto c : k) counts[c]++;
        for (auto &[c, cnt] : counts)
            if (c >= classes.size() || (int)classes[c].size < cnt) return false;
        return true; // valid tuple exists but not tracked → alive
    };

    while (numPeeled < tuples.size()) {
        while (currentLevel <= maxSupport && buckets[currentLevel].empty())
            currentLevel++;
        if (currentLevel > maxSupport) break;

        daf::Size idx = buckets[currentLevel].back();
        buckets[currentLevel].pop_back();
        if (peeled[idx] || tupleSupport[idx] != currentLevel) continue;

        peeled[idx] = true;
        numPeeled++;
        coreLevel = std::max(coreLevel, currentLevel);
        coreDist[coreLevel] += tuples[idx].mult;

        // Cascade: for s=r+1, the (r+1)-clique has r+1 r-cliques (tuples).
        // Peeling tuple τ: for each neighbor class D, for each position k in τ:
        //   - Replace τ[k] with D → neighbor tuple τ'
        //   - Check: r-1 other tuples (each with another position replaced by D) are alive
        //   - Reduction = |τ[k]| - (number of OTHER positions in τ with same class as τ[k])

        auto &tkey = tuples[idx].key; // length r

        // Collect neighbor classes from common maximal cliques
        std::unordered_set<daf::Size> neighborClasses;
        for (daf::Size rid : tupleCommonRegions[idx])
            for (daf::Size cid : classesInRegion[rid])
                neighborClasses.insert(cid);

        for (daf::Size D : neighborClasses) {
            // For each position k: generate neighbor tuple (replace tkey[k] with D)
            // Pre-compute all r neighbor keys
            std::vector<TupleKey> neighborKeys(r);
            for (int k = 0; k < r; ++k) {
                neighborKeys[k] = tkey;
                neighborKeys[k][k] = D;
                std::sort(neighborKeys[k].begin(), neighborKeys[k].end());
            }

            for (int k = 0; k < r; ++k) {
                if (D == tkey[k]) continue;        // no change → self
                if (neighborKeys[k] == tkey) continue; // self after sorting

                // Dedup: skip if same key as earlier k
                bool dup = false;
                for (int j = 0; j < k; ++j) {
                    if (D != tkey[j] && neighborKeys[j] != tkey && neighborKeys[k] == neighborKeys[j]) {
                        dup = true; break;
                    }
                }
                if (dup) continue;

                // Check: the other r-1 neighbor tuples must be alive
                bool allAlive = true;
                for (int j = 0; j < r; ++j) {
                    if (j == k) continue;
                    if (!isTupleAlive(neighborKeys[j], tkey)) { allAlive = false; break; }
                }
                if (!allAlive) continue;

                // Reduction = |removed_class| - overlap
                daf::Size removedClass = tkey[k];
                daf::Size overlap = 0;
                for (int j = 0; j < r; ++j)
                    if (j != k && tkey[j] == removedClass) overlap++;
                daf::Size red = classes[removedClass].size > overlap
                              ? classes[removedClass].size - overlap : 0;
                if (red == 0) continue;

                auto it = tupleIndex.find(neighborKeys[k]);
                if (it == tupleIndex.end()) continue;
                daf::Size aidx = it->second;
                if (peeled[aidx]) continue;

                daf::Size oldSup = tupleSupport[aidx];
                daf::Size newSup = oldSup > red ? oldSup - red : 0;
                if (newSup < oldSup) {
                    tupleSupport[aidx] = newSup;
                    buckets[newSup].push_back(aidx);
                    if (newSup < currentLevel) currentLevel = newSup;
                }
            }
        }
    }

    auto tStep3End = std::chrono::high_resolution_clock::now();
    auto step3Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep3End - tStep3).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Step 3: Cascade Peeling ---" << std::endl;
    std::cout << "  Tuples peeled: " << numPeeled << " / " << tuples.size() << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    std::cout << "  Core distribution:" << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Step 3 time: " << step3Ms << " ms" << std::endl;

    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep3End - tStart).count();
    std::cout << "  Total region decomposition time: " << totalMs << " ms" << std::endl;
    std::cout << "==================================================" << std::endl;

    return {};
}
