//
// Region Tuple + CPI Counting (V3)
//
// Combines two orthogonal compressions:
// - Region Tuple: avoids r-clique enumeration (groups into tuples)
// - CPI Counting: avoids s-tuple enumeration (counts via path formula)
//
// Init: support(τ) = [Σ_P N(τ,P) × C(p-r+h, s-r)] / mult(τ)
// Peel: lazy s-tuple generation (only enumerate dying s-tuples during peeling)
//

#include "NCliqueCoreDecomposition.h"
#include <algorithm>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

extern double nCr[1001][401];
extern std::vector<bool> g_maxCliqueTags;
extern std::vector<std::vector<daf::Size>> g_maxCliques;

// ============================================================
// Tuple utilities (shared with V2)
// ============================================================

using TupleKey = std::vector<daf::Size>;

struct TupleHash {
    size_t operator()(const TupleKey &t) const noexcept {
        size_t h = t.size();
        for (auto x : t) h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h << 6) + (h >> 2);
        return h;
    }
};

static void enumerateMultisets(
    const std::vector<daf::Size> &classes, int size, int startIdx,
    TupleKey &current, const std::function<void()> &callback)
{
    if ((int)current.size() == size) { callback(); return; }
    for (int i = startIdx; i < (int)classes.size(); ++i) {
        current.push_back(classes[i]);
        enumerateMultisets(classes, size, i, current, callback);
        current.pop_back();
    }
}

// Enumerate s-multisets that contain tau as a sub-multiset
// classes: available class IDs (sorted), tauCounts: minimum count per class from tau
// Fills remaining s-r slots from classes with repetition allowed up to classSizes
static void enumerateSupersetsOfTau(
    const std::vector<daf::Size> &classes,
    const std::vector<daf::Size> &classSizes,
    const std::unordered_map<daf::Size, int> &tauCounts,
    int s, int startIdx, TupleKey &current,
    const std::function<void(const TupleKey &)> &callback)
{
    if ((int)current.size() == s) { callback(current); return; }
    int remaining = s - (int)current.size();
    for (int i = startIdx; i < (int)classes.size(); ++i) {
        daf::Size c = classes[i];
        // How many of c already in current?
        int alreadyUsed = 0;
        for (auto x : current) if (x == c) alreadyUsed++;
        if (alreadyUsed >= (int)classSizes[c]) continue; // can't use more
        current.push_back(c);
        enumerateSupersetsOfTau(classes, classSizes, tauCounts, s, i, current, callback);
        current.pop_back();
    }
}

static std::vector<std::pair<daf::Size, int>> getComposition(const TupleKey &t) {
    std::vector<std::pair<daf::Size, int>> comp;
    for (auto c : t) {
        if (!comp.empty() && comp.back().first == c) comp.back().second++;
        else comp.push_back({c, 1});
    }
    return comp;
}

static void enumerateSubMultisets(
    const std::vector<std::pair<daf::Size, int>> &composition, int r,
    int classIdx, TupleKey &current,
    const std::function<void(const TupleKey &)> &callback)
{
    if (classIdx == (int)composition.size()) {
        if ((int)current.size() == r) callback(current);
        return;
    }
    auto [cls, maxCnt] = composition[classIdx];
    int remaining = r - (int)current.size();
    int maxFromLater = 0;
    for (int k = classIdx + 1; k < (int)composition.size(); ++k)
        maxFromLater += composition[k].second;
    int minJ = std::max(0, remaining - maxFromLater);
    int maxJ = std::min(maxCnt, remaining);

    for (int j = minJ; j <= maxJ; ++j) {
        for (int q = 0; q < j; ++q) current.push_back(cls);
        enumerateSubMultisets(composition, r, classIdx + 1, current, callback);
        for (int q = 0; q < j; ++q) current.pop_back();
    }
}

static double computeExt(const std::vector<std::pair<daf::Size, int>> &sigmaComp,
                          const TupleKey &tau,
                          const std::vector<daf::Size> &classSizes) {
    std::unordered_map<daf::Size, int> tauCounts;
    for (auto c : tau) tauCounts[c]++;
    double ext = 1.0;
    for (auto &[cls, mi] : sigmaComp) {
        int ji = tauCounts.count(cls) ? tauCounts[cls] : 0;
        int n = (int)classSizes[cls] - ji;
        int k = mi - ji;
        if (n < k || k < 0) return 0.0;
        ext *= nCr[n][k];
    }
    return ext;
}

// ============================================================
// Main function: Region CPI (V3)
// ============================================================

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionCPI(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    auto tStart = std::chrono::high_resolution_clock::now();
    const daf::Size numVertices = edgeGraph.n;
    const daf::Size numPaths = tree.adj_list.size();
    const daf::Size INVALID = static_cast<daf::Size>(-1);

    // ============================================================
    // Step 0: Connected Component Decomposition of Clique Overlap Graph
    // ============================================================
    // Key insight: maximal cliques that don't overlap can be processed
    // independently. We decompose into connected components and handle:
    //   - Isolated cliques (no overlap): closed-form core = C(|M|-r, s-r)
    //   - Multi-clique components: full V3 pipeline (Steps 1-6)

    // First, collect valid maximal cliques (size >= s)
    std::vector<daf::Size> mcIndices; // indices into g_maxCliques
    for (daf::Size i = 0; i < g_maxCliques.size(); ++i)
        if ((int)g_maxCliques[i].size() >= s) mcIndices.push_back(i);
    daf::Size numMC = mcIndices.size();

    // Build vertex -> clique membership (for valid cliques only)
    // mcLocalId: position in mcIndices (0..numMC-1)
    std::vector<std::vector<daf::Size>> vtxToMC(numVertices);
    for (daf::Size lid = 0; lid < numMC; ++lid) {
        for (daf::Size v : g_maxCliques[mcIndices[lid]])
            if (v < numVertices) vtxToMC[v].push_back(lid);
    }

    // Union-Find for clique overlap graph
    std::vector<daf::Size> uf_parent(numMC), uf_rank(numMC, 0);
    for (daf::Size i = 0; i < numMC; ++i) uf_parent[i] = i;
    std::function<daf::Size(daf::Size)> uf_find = [&](daf::Size x) -> daf::Size {
        while (uf_parent[x] != x) { uf_parent[x] = uf_parent[uf_parent[x]]; x = uf_parent[x]; }
        return x;
    };
    auto uf_union = [&](daf::Size a, daf::Size b) {
        a = uf_find(a); b = uf_find(b);
        if (a == b) return;
        if (uf_rank[a] < uf_rank[b]) std::swap(a, b);
        uf_parent[b] = a;
        if (uf_rank[a] == uf_rank[b]) uf_rank[a]++;
    };

    // Two cliques overlap if they share a vertex
    for (daf::Size v = 0; v < numVertices; ++v) {
        auto &mcs = vtxToMC[v];
        for (daf::Size i = 1; i < mcs.size(); ++i)
            uf_union(mcs[0], mcs[i]);
    }

    // Group cliques by component
    std::unordered_map<daf::Size, std::vector<daf::Size>> components; // root -> list of local IDs
    for (daf::Size lid = 0; lid < numMC; ++lid)
        components[uf_find(lid)].push_back(lid);

    // Classify components
    daf::Size numIsolated = 0, numLargeComp = 0;
    int64_t isolatedRCliques = 0;
    std::map<daf::Size, int64_t> coreDist; // accumulated core distribution
    std::vector<daf::Size> largePoolCliques; // local IDs for full V3

    for (auto &[root, members] : components) {
        if (members.size() == 1) {
            // --- Isolated clique: core = C(|M|-r, s-r) for all r-cliques ---
            daf::Size lid = members[0];
            int mSize = (int)g_maxCliques[mcIndices[lid]].size();
            daf::Size coreVal = (daf::Size)llround(nCr[mSize - r][s - r]);
            daf::Size numRCliquesInMC = (daf::Size)llround(nCr[mSize][r]);
            coreDist[coreVal] += numRCliquesInMC;
            isolatedRCliques += numRCliquesInMC;
            numIsolated++;
            continue;
        }

        // Multi-clique component: send to full V3 pipeline
        numLargeComp++;
        for (daf::Size lid : members) largePoolCliques.push_back(lid);
    }

    auto tStep0 = std::chrono::high_resolution_clock::now();
    auto step0Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep0 - tStart).count();

    std::cout << "======= Region CPI (V3) =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Step 0: Connected Component Decomposition" << std::endl;
    std::cout << "    Valid maximal cliques (>= s): " << numMC << std::endl;
    std::cout << "    Components: " << components.size() << std::endl;
    std::cout << "      Isolated (1 clique):  " << numIsolated
              << " (" << isolatedRCliques << " r-cliques, closed-form)" << std::endl;
    std::cout << "      Multi-clique (full V3): " << numLargeComp
              << " (" << largePoolCliques.size() << " cliques)" << std::endl;
    std::cout << "    Step 0 time: " << step0Ms << " ms" << std::endl;

    // If all cliques handled by Step 0, skip Steps 1-6
    if (largePoolCliques.empty()) {
        auto tEnd = std::chrono::high_resolution_clock::now();
        auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();
        daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
        std::cout << "\n  All cliques handled by Step 0 decomposition." << std::endl;
        std::cout << "  Max core: " << maxCore << std::endl;
        for (auto &[core, cnt] : coreDist)
            std::cout << "  core=" << core << " count=" << cnt << std::endl;
        std::cout << "  Total time: " << totalMs << " ms" << std::endl;
        std::cout << "==============================================" << std::endl;
        return {};
    }

    // ============================================================
    // Step 1: Build Regions from large-component cliques only
    // ============================================================

    daf::Size validPaths = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid)
        if ((int)tree.adj_list[pid].size() >= s) validPaths++;

    daf::Size numRegions = 0;
    std::vector<std::vector<daf::Size>> regionVerts;
    std::vector<std::vector<daf::Size>> vtxMaxPaths(numVertices);

    // Only include cliques from large components
    for (daf::Size lid : largePoolCliques) {
        auto &mc = g_maxCliques[mcIndices[lid]];
        daf::Size rid = regionVerts.size();
        regionVerts.push_back(mc);
        for (daf::Size v : mc)
            if (v < numVertices) vtxMaxPaths[v].push_back(rid);
        numRegions++;
    }

    // ============================================================
    // Step 2: Build Overlap Classes (standard)
    // ============================================================

    struct ClassInfo { std::vector<daf::Size> regionIds; daf::Size size; };
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
            } else { classOf[v] = it->second; classes[it->second].size++; }
        }
    }
    daf::Size numClasses = classes.size();

    std::vector<daf::Size> classSizes(numClasses);
    for (daf::Size i = 0; i < numClasses; ++i) classSizes[i] = classes[i].size;

    std::vector<std::vector<daf::Size>> classesInRegion(numRegions);
    for (daf::Size cid = 0; cid < numClasses; ++cid)
        for (daf::Size rid : classes[cid].regionIds)
            classesInRegion[rid].push_back(cid);
    for (auto &v : classesInRegion) std::sort(v.begin(), v.end());

    // ============================================================
    // Step 2b: r-Mergeable class analysis (for tuple merging in Step 3b)
    // ============================================================
    // A class C is r-mergeable in region R if for ALL other regions R' in C's
    // profile: |R ∩ R'| < r. This means any r-clique using vertices of C within
    // R cannot extend to another region.
    //
    // Consequence: In region R, if ALL classes of a tuple T are r-mergeable
    // (or private to R), then ALL r-cliques of type T within R have the same
    // support = C(|R|-r, s-r). Such tuples can be MERGED into a single
    // "super-tuple" with combined multiplicity.

    // For class C with profile P_C = {R1,..,Rk}:
    //   C is r-mergeable in Ri iff for ALL j != i: |Ri ∩ Rj| < r
    //
    // Key insight: we don't need ALL pairwise intersection sizes. We only need
    // to check pairs from the same class profile. And we can use a threshold
    // approach: a class is NOT r-mergeable in Ri if ANY Rj in its profile has
    // |Ri ∩ Rj| >= r.
    //
    // Fast approach: for each multi-region class C, check if its profile regions
    // have pairwise intersection >= r. Use per-region neighbor lists to quickly
    // compute intersection sizes for just the needed pairs.

    // classRMergeable[cid][rid] = true if class cid is r-mergeable in region rid
    std::vector<std::unordered_map<daf::Size, bool>> classRMergeable(numClasses);
    daf::Size mergeableClassRegionPairs = 0;

    // Build intersection data lazily: for each multi-region class, compute
    // pairwise intersection sizes between its profile regions.
    // Use the fact: |Ri ∩ Rj| = number of vertices in both Ri and Rj
    //             = Σ classes C' with both Ri and Rj in C'.profile : |C'|

    // Approach: Build per-region -> set of (class, size) for multi-region classes.
    // Then for a pair (Ri, Rj): |Ri ∩ Rj| = Σ of sizes of classes appearing in both.

    // regionMultiClasses[rid] = list of (cid, size) for multi-region classes in rid
    std::vector<std::vector<std::pair<daf::Size, daf::Size>>> regionMultiClasses(numRegions);
    for (daf::Size cid = 0; cid < numClasses; ++cid) {
        if (classes[cid].regionIds.size() <= 1) continue;
        for (daf::Size rid : classes[cid].regionIds)
            regionMultiClasses[rid].push_back({cid, classes[cid].size});
    }

    // Cache for pairwise intersections (computed on demand)
    std::unordered_map<uint64_t, int> pairIntersectionCache;
    auto pairKey = [](daf::Size a, daf::Size b) -> uint64_t {
        return ((uint64_t)std::min(a,b) << 32) | (uint64_t)std::max(a,b);
    };

    auto getIntersection = [&](daf::Size r1, daf::Size r2) -> int {
        uint64_t key = pairKey(r1, r2);
        auto it = pairIntersectionCache.find(key);
        if (it != pairIntersectionCache.end()) return it->second;

        // Compute: iterate multi-region classes of the smaller region list
        int isz = 0;
        auto &listA = regionMultiClasses[r1];
        auto &listB = regionMultiClasses[r2];
        auto &smaller = (listA.size() <= listB.size()) ? listA : listB;
        daf::Size otherRid = (listA.size() <= listB.size()) ? r2 : r1;

        for (auto &[cid, sz] : smaller) {
            // Check if this class is also in otherRid
            bool inOther = false;
            for (daf::Size rid3 : classes[cid].regionIds)
                if (rid3 == otherRid) { inOther = true; break; }
            if (inOther) isz += (int)sz;
        }

        pairIntersectionCache[key] = isz;
        return isz;
    };

    daf::Size skippedLargeClass = 0, skippedLargeProfile = 0;
    for (daf::Size cid = 0; cid < numClasses; ++cid) {
        auto &regs = classes[cid].regionIds;
        if (regs.size() <= 1) {
            for (daf::Size rid : regs)
                classRMergeable[cid][rid] = true;
            continue;
        }

        // Quick filter 1: |C| >= r → |Ri ∩ Rj| >= |C| >= r → not mergeable
        if ((int)classes[cid].size >= (int)r) {
            skippedLargeClass++;
            continue; // classRMergeable entries default to absent (= non-mergeable)
        }

        // Quick filter 2: profile size > 500 → skip (too expensive to check,
        // and almost certainly non-mergeable due to large intersections)
        if (regs.size() > 500) {
            skippedLargeProfile++;
            continue;
        }

        for (daf::Size rid : regs) {
            bool mergeable = true;
            for (daf::Size rid2 : regs) {
                if (rid == rid2) continue;
                int isz = getIntersection(rid, rid2);
                if (isz >= (int)r) { mergeable = false; break; }
            }
            classRMergeable[cid][rid] = mergeable;
            if (mergeable) mergeableClassRegionPairs++;
        }
    }

    daf::Size pairsCached = pairIntersectionCache.size();
    // Free data structures no longer needed
    regionMultiClasses.clear();
    regionMultiClasses.shrink_to_fit();
    pairIntersectionCache.clear();

    auto tStep2 = std::chrono::high_resolution_clock::now();
    auto step2Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep2 - tStep0).count();
    std::cout << "  Steps 1-2 (large components):" << std::endl;
    std::cout << "    Regions (large pool): " << numRegions << std::endl;
    std::cout << "    CPI paths: " << validPaths << std::endl;
    std::cout << "    Overlap classes (before merge): " << numClasses << std::endl;
    std::cout << "    r-mergeable multi-region (class,region) pairs: " << mergeableClassRegionPairs << std::endl;
    std::cout << "    Pair intersection cache: " << pairsCached << " pairs" << std::endl;
    std::cout << "    Step 1-2 time: " << step2Ms << " ms" << std::endl;

    // ============================================================
    // Step 2c: r-Mergeable Class Merging
    // ============================================================
    // For each region R, merge all r-mergeable classes into a SINGLE merged
    // class M_R. This dramatically reduces class count and tuple count on
    // sparse graphs.
    //
    // Key invariant: each merged class M_R belongs to exactly ONE region R,
    // so classSizes[M_R] is globally unique (= count of mergeable vertices in R).
    //
    // Challenge: vertex v in multiple regions may be mergeable in all of them,
    // so it belongs to M_R1 in R1 and M_R2 in R2. We can't store both in
    // classOf[v]. Solution: for CPI paths, determine the path's region and use
    // a per-region class lookup for mergeable vertices.
    //
    // Data structures:
    //   mergedClassForRegion[rid] = merged class ID for region rid (INVALID if none)
    //   vertexIsMergeable[v]      = true if v is mergeable in ANY region
    //   pathRegion[pid]           = region ID that contains CPI path pid
    //   For path processing: if v is mergeable, use mergedClassForRegion[pathRegion]
    //                        otherwise use classOf[v] (unchanged)

    // A class C is "fully mergeable" if C is r-mergeable in ALL regions of its
    // profile. Only fully-mergeable classes can be safely merged: their vertices
    // are removed from C (size→0) and placed into per-region merged classes.
    //
    // Partial mergeability (r-mergeable in SOME regions) would require per-region
    // class sizes, which conflicts with the global tuple mult framework.
    daf::Size numOrigClasses = numClasses;
    // r-merging: disabled by default (increases split count for large s).
    // Enable with PIVOTER_MERGE=1 for experiments.
    const bool mergeDisabled = (getenv("PIVOTER_MERGE") == nullptr);

    std::vector<bool> classFullyMergeable(numOrigClasses, false);
    if (!mergeDisabled) {
      for (daf::Size cid = 0; cid < numOrigClasses; ++cid) {
        auto &regs = classes[cid].regionIds;
        if (regs.empty()) continue;
        bool allMerge = true;
        for (daf::Size rid : regs) {
            auto it = classRMergeable[cid].find(rid);
            if (it == classRMergeable[cid].end() || !it->second) {
                allMerge = false; break;
            }
        }
        classFullyMergeable[cid] = allMerge;
      }
    }

    // Identify fully-mergeable vertices and count per region
    std::vector<bool> vertexIsMergeable(numVertices, false);
    std::vector<daf::Size> mergeableCount(numRegions, 0);
    daf::Size fullyMergeableVertices = 0;
    for (daf::Size v = 0; v < numVertices; ++v) {
        if (classOf[v] == INVALID) continue;
        if (!classFullyMergeable[classOf[v]]) continue;
        vertexIsMergeable[v] = true;
        fullyMergeableVertices++;
        for (daf::Size rid : classes[classOf[v]].regionIds)
            mergeableCount[rid]++;
    }

    // Create merged classes only for regions with fully-mergeable vertices
    std::vector<daf::Size> mergedClassForRegion(numRegions, INVALID);
    daf::Size numMergedClasses = 0;
    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        if (mergeableCount[rid] == 0) continue;
        daf::Size mcid = classes.size();
        classes.push_back({{rid}, mergeableCount[rid]});
        mergedClassForRegion[rid] = mcid;
        numMergedClasses++;
    }
    numClasses = classes.size();
    classSizes.resize(numClasses);
    for (daf::Size i = numOrigClasses; i < numClasses; ++i)
        classSizes[i] = classes[i].size;

    // Zero out fully-mergeable original classes (their vertices are now in merged classes)
    for (daf::Size cid = 0; cid < numOrigClasses; ++cid) {
        if (!classFullyMergeable[cid]) continue;
        classSizes[cid] = 0;
        classes[cid].size = 0;
    }

    // Rebuild classesInRegion with merged classes
    classesInRegion.assign(numRegions, {});
    for (daf::Size cid = 0; cid < numClasses; ++cid) {
        if (classSizes[cid] == 0) continue; // skip empty classes
        for (daf::Size rid : classes[cid].regionIds)
            classesInRegion[rid].push_back(cid);
    }
    for (auto &v : classesInRegion) std::sort(v.begin(), v.end());

    // Build path-to-region mapping: for each CPI path, find which region
    // contains all its vertices. Use the first vertex's region set and intersect.
    std::vector<daf::Size> pathRegion(numPaths, INVALID);
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        // Find a vertex in the path that's in the large pool
        daf::Size bestRid = INVALID;
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            if (v >= numVertices) continue;
            auto &vregs = vtxMaxPaths[v];
            if (vregs.empty()) continue;
            if (bestRid == INVALID) {
                // First vertex with regions: pick its first region as candidate
                // (any region containing this vertex)
                bestRid = vregs[0];
            }
            // Verify: is bestRid still valid for this vertex?
            bool found = false;
            for (daf::Size rid : vregs)
                if (rid == bestRid) { found = true; break; }
            if (!found) {
                // Need to intersect: find a region common to all vertices so far
                // For correctness, try all regions of this vertex
                bestRid = INVALID;
                for (daf::Size candRid : vregs) {
                    // Check if candRid contains all previous vertices
                    // (expensive but paths are small)
                    bool ok = true;
                    for (const auto &prevNode : leaf) {
                        if (&prevNode == &node) break;
                        daf::Size pv = prevNode.v;
                        if (pv >= numVertices) continue;
                        auto &pvregs = vtxMaxPaths[pv];
                        if (pvregs.empty()) continue;
                        bool inCand = false;
                        for (daf::Size rid2 : pvregs)
                            if (rid2 == candRid) { inCand = true; break; }
                        if (!inCand) { ok = false; break; }
                    }
                    if (ok) { bestRid = candRid; break; }
                }
                if (bestRid == INVALID) break;
            }
        }
        pathRegion[pid] = bestRid;
    }

    // Helper: get the effective class of vertex v on a path in region rid.
    // If v is fully mergeable and rid has a merged class, return that merged class.
    // Otherwise return classOf[v] (which may have size 0 if fully merged; caller skips).
    auto effectiveClass = [&](daf::Size v, daf::Size rid) -> daf::Size {
        if (v >= numVertices || classOf[v] == INVALID) return INVALID;
        if (!vertexIsMergeable[v] || rid == INVALID) return classOf[v];
        daf::Size mcid = mergedClassForRegion[rid];
        return (mcid != INVALID) ? mcid : classOf[v];
    };

    // Count original classes-in-region entries (for comparison, before merge)
    // All original classes with non-zero original size are counted
    daf::Size origClassesInRegionEntries = 0;
    for (daf::Size cid = 0; cid < numOrigClasses; ++cid) {
        // classes[cid].size is now 0 for merged classes; use classFullyMergeable to detect
        if (classes[cid].size > 0 || classFullyMergeable[cid])
            origClassesInRegionEntries += classes[cid].regionIds.size();
    }
    // Count single-region classes that were merged (they reduce entries)
    daf::Size singleRegMerged = 0, multiRegMerged = 0;
    for (daf::Size cid = 0; cid < numOrigClasses; ++cid) {
        if (!classFullyMergeable[cid]) continue;
        if (classes[cid].regionIds.size() == 1) singleRegMerged++;
        else multiRegMerged++;
    }

    auto tStep2c = std::chrono::high_resolution_clock::now();
    auto step2cMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep2c - tStep2).count();
    std::cout << "  Step 2c: r-Mergeable Class Merging" << std::endl;
    std::cout << "    Fully-mergeable: " << singleRegMerged << " single-region + "
              << multiRegMerged << " multi-region = "
              << (singleRegMerged + multiRegMerged) << " / " << numOrigClasses << std::endl;
    std::cout << "    Fully-mergeable vertices: " << fullyMergeableVertices << std::endl;
    std::cout << "    Merged classes created: " << numMergedClasses << std::endl;
    std::cout << "    Active classes (after merge): "
              << std::count_if(classSizes.begin(), classSizes.end(), [](daf::Size s){return s>0;})
              << std::endl;
    {
        daf::Size totalClassesInRegions = 0;
        for (auto &v : classesInRegion) totalClassesInRegions += v.size();
        std::cout << "    Classes-in-region entries: " << origClassesInRegionEntries
                  << " -> " << totalClassesInRegions << std::endl;
    }
    std::cout << "    Merge time: " << step2cMs << " ms" << std::endl;

    // ============================================================
    // Step 3: Enumerate r-tuples (using merged classes)
    // ============================================================

    std::unordered_map<TupleKey, daf::Size, TupleHash> rTupleIndex;
    struct RTuple { TupleKey key; daf::Size mult; };
    std::vector<RTuple> rTuples;

    {
        TupleKey cur; cur.reserve(r);
        auto addRTuple = [&](const TupleKey &key) {
            std::unordered_map<daf::Size, int> counts;
            for (auto c : key) counts[c]++;
            daf::Size mult = 1;
            for (auto &[c, k] : counts) {
                if ((int)classSizes[c] < k) return;
                mult *= (daf::Size)nCr[classSizes[c]][k];
            }
            if (mult == 0) return;
            if (rTupleIndex.count(key)) return;
            rTupleIndex[key] = rTuples.size();
            rTuples.push_back({key, mult});
        };

        for (daf::Size rid = 0; rid < numRegions; ++rid) {
            auto &cids = classesInRegion[rid];
            if (cids.size() > 500) continue;
            cur.clear();
            std::function<void()> cb = [&]() { addRTuple(cur); };
            enumerateMultisets(cids, r, 0, cur, cb);
        }
    }

    auto tStep3 = std::chrono::high_resolution_clock::now();
    daf::Size numRawTuples = rTuples.size();
    std::cout << "  Raw r-tuples (after merge): " << numRawTuples << std::endl;

    // ============================================================
    // Step 3b: r-Mergeable Tuple Classification
    // ============================================================
    // A tuple T is "closed-form" in region R if ALL its classes are r-mergeable
    // in R (including single-region classes which are trivially private).
    // Such tuples have support = C(|R|-r, s-r) for all their r-cliques.
    //
    // We merge all closed-form tuples for the same region into a SUPER-TUPLE:
    //   super_mult = sum of mults of constituent tuples = C(mergeablePoolSize, r)
    //   support = C(|R|-r, s-r) (closed form, no CPI needed)
    //
    // Non-closed-form tuples keep individual support computed via CPI.

    // For each tuple, determine if it's closed-form and in which region
    std::vector<bool> tupleIsClosedForm(rTuples.size(), false);
    std::vector<daf::Size> tupleClosedFormRegion(rTuples.size(), INVALID);

    // For each region, compute the "mergeable pool" size = sum of sizes of
    // all r-mergeable classes in that region (using post-merge class set)
    // Note: merged classes M_R are single-region, size = mergeableCount[R],
    // and are trivially "mergeable" in their region.
    std::vector<daf::Size> mergeablePoolSizeV(numRegions, 0);
    for (daf::Size rid = 0; rid < numRegions; ++rid) {
        for (daf::Size cid : classesInRegion[rid]) {
            if (classSizes[cid] == 0) continue;
            if (cid >= numOrigClasses) {
                // Merged class: always mergeable in its region
                mergeablePoolSizeV[rid] += classSizes[cid];
            } else if (classes[cid].regionIds.size() == 1) {
                // Single-region original class (not fully-mergeable, else size=0): mergeable
                mergeablePoolSizeV[rid] += classSizes[cid];
            } else {
                auto it = classRMergeable[cid].find(rid);
                if (it != classRMergeable[cid].end() && it->second)
                    mergeablePoolSizeV[rid] += classSizes[cid];
            }
        }
    }

    // Pre-compute: for each class, the set of regions where it's r-mergeable
    // (includes single-region classes and merged classes which are trivially mergeable)
    std::vector<std::vector<daf::Size>> classMergeableRegions(numClasses);
    for (daf::Size cid = 0; cid < numClasses; ++cid) {
        if (classSizes[cid] == 0) continue;
        if (cid >= numOrigClasses || classes[cid].regionIds.size() == 1) {
            // Merged class or single-region: trivially mergeable in its region(s)
            classMergeableRegions[cid] = classes[cid].regionIds;
        } else {
            for (auto &[rid, mg] : classRMergeable[cid])
                if (mg) classMergeableRegions[cid].push_back(rid);
            std::sort(classMergeableRegions[cid].begin(), classMergeableRegions[cid].end());
        }
    }

    for (daf::Size tidx = 0; tidx < rTuples.size(); ++tidx) {
        auto &key = rTuples[tidx].key;
        // Collect unique classes in this tuple
        std::vector<daf::Size> tupleCids;
        {
            daf::Size prev = INVALID;
            for (auto c : key) if (c != prev) { tupleCids.push_back(c); prev = c; }
        }

        // Find a region where ALL classes are r-mergeable:
        // Intersect the mergeable-region sets of all classes
        // Start with the smallest set for efficiency
        int smallestIdx = 0;
        size_t smallestSize = classMergeableRegions[tupleCids[0]].size();
        for (int i = 1; i < (int)tupleCids.size(); ++i) {
            if (classMergeableRegions[tupleCids[i]].size() < smallestSize) {
                smallestSize = classMergeableRegions[tupleCids[i]].size();
                smallestIdx = i;
            }
        }

        // Check each region in the smallest set
        for (daf::Size rid : classMergeableRegions[tupleCids[smallestIdx]]) {
            bool allMergeable = true;
            for (int i = 0; i < (int)tupleCids.size(); ++i) {
                if (i == smallestIdx) continue;
                auto &regs = classMergeableRegions[tupleCids[i]];
                if (!std::binary_search(regs.begin(), regs.end(), rid)) {
                    allMergeable = false; break;
                }
            }
            if (allMergeable) {
                tupleIsClosedForm[tidx] = true;
                tupleClosedFormRegion[tidx] = rid;
                break;
            }
        }
    }

    // Build super-tuples: one per region that has closed-form tuples
    // superTuples[rid] = {combined mult, support}
    struct SuperTuple { daf::Size mult; double support; daf::Size region; };
    std::vector<SuperTuple> superTuples;
    std::unordered_map<daf::Size, daf::Size> regionToSuperIdx; // rid -> super-tuple index

    daf::Size closedFormTuples = 0, closedFormRCliques = 0;
    daf::Size cpiTuples = 0;

    for (daf::Size tidx = 0; tidx < rTuples.size(); ++tidx) {
        if (!tupleIsClosedForm[tidx]) { cpiTuples++; continue; }
        closedFormTuples++;
        closedFormRCliques += rTuples[tidx].mult;
        daf::Size rid = tupleClosedFormRegion[tidx];
        auto it = regionToSuperIdx.find(rid);
        if (it == regionToSuperIdx.end()) {
            regionToSuperIdx[rid] = superTuples.size();
            int regSize = (int)regionVerts[rid].size();
            double sup = (regSize >= (int)s) ? nCr[regSize - r][s - r] : 0.0;
            superTuples.push_back({rTuples[tidx].mult, sup, rid});
        } else {
            superTuples[it->second].mult += rTuples[tidx].mult;
        }
    }

    // Verify: for each region with a super-tuple, the combined mult should
    // equal C(mergeablePoolSizeV[rid], r)
    for (auto &[rid, sidx] : regionToSuperIdx) {
        daf::Size expected = (daf::Size)llround(nCr[mergeablePoolSizeV[rid]][r]);
        if (superTuples[sidx].mult != expected) {
            std::cout << "  WARNING: super-tuple region " << rid
                      << " mult=" << superTuples[sidx].mult
                      << " expected C(" << mergeablePoolSizeV[rid] << "," << r
                      << ")=" << expected << std::endl;
        }
    }

    auto tStep3b = std::chrono::high_resolution_clock::now();
    auto step3bMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep3b - tStep3).count();
    std::cout << "  Closed-form tuples: " << closedFormTuples << " / " << numRawTuples
              << " (" << closedFormRCliques << " r-cliques)" << std::endl;
    std::cout << "  Super-tuples (regions): " << superTuples.size() << std::endl;
    std::cout << "  CPI-needed tuples: " << cpiTuples << std::endl;
    std::cout << "  Merge classification time: " << step3bMs << " ms" << std::endl;

    // ============================================================
    // Step 4: CPI Counting for initial support
    // ============================================================
    // Only compute CPI support for non-closed-form tuples.
    // Closed-form tuples already have support = C(|R|-r, s-r).

    std::vector<double> support(rTuples.size(), 0.0);

    // Pre-assign closed-form support
    for (daf::Size tidx = 0; tidx < rTuples.size(); ++tidx) {
        if (!tupleIsClosedForm[tidx]) continue;
        daf::Size rid = tupleClosedFormRegion[tidx];
        int regSize = (int)regionVerts[rid].size();
        support[tidx] = (regSize >= (int)s) ? nCr[regSize - r][s - r] : 0.0;
    }

    daf::Size pathsUsed = 0;
    daf::Size totalTuplesOnPaths = 0;
    daf::Size pathsWithRegion = 0, pathsWithoutRegion = 0;

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        daf::Size pRid = pathRegion[pid]; // region for this path

        // Compute h, p and class distribution on this path
        // Use effectiveClass to map mergeable vertices to merged classes
        int h = 0, p = 0;
        std::unordered_map<daf::Size, std::pair<int,int>> classDistrib; // cid -> (nh, np)

        for (const auto &node : leaf) {
            daf::Size v = node.v;
            daf::Size cid = effectiveClass(v, pRid);
            if (cid == INVALID) continue;
            if (classSizes[cid] == 0) continue; // class was fully merged away
            if (node.isPivot) { p++; classDistrib[cid].second++; }
            else { h++; classDistrib[cid].first++; }
        }

        if (h + p < (int)s) continue;
        pathsUsed++;
        if (pRid != INVALID) pathsWithRegion++;
        else pathsWithoutRegion++;

        // Collect unique classes on this path
        std::vector<daf::Size> pathClasses;
        for (auto &[cid, _] : classDistrib) pathClasses.push_back(cid);
        std::sort(pathClasses.begin(), pathClasses.end());

        // Enumerate r-multisets of this path's classes
        TupleKey cur; cur.reserve(r);
        std::function<void()> cb = [&]() {
            // Look up this r-tuple
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;

            // Skip closed-form tuples (support already assigned)
            if (tupleIsClosedForm[tidx]) return;

            // Build composition: (class, count) pairs
            std::vector<std::pair<daf::Size, int>> tauClasses;
            {
                daf::Size prev = INVALID; int cnt = 0;
                for (auto c : cur) {
                    if (c == prev) cnt++;
                    else { if (prev != INVALID) tauClasses.push_back({prev, cnt}); prev = c; cnt = 1; }
                }
                if (prev != INVALID) tauClasses.push_back({prev, cnt});
            }

            // Feasibility: for each class c with j_c copies in tau,
            // need j_c <= n_h^c + n_p^c
            for (auto &[c, jc] : tauClasses) {
                auto dit = classDistrib.find(c);
                if (dit == classDistrib.end()) return;
                if (jc > dit->second.first + dit->second.second) return;
            }

            // Convolution to compute AggrCount
            // f[t] = coefficient of x^t in Π_c g_c(x)
            // where g_c(x) = Σ_{b_c} C(nh_c, j_c-b_c) × C(np_c, b_c) × x^{b_c}

            std::vector<double> f = {1.0}; // polynomial coefficients, f[0]=1

            for (auto &[c, jc] : tauClasses) {
                auto &[nhc, npc] = classDistrib[c];
                int bMin = std::max(0, jc - nhc);
                int bMax = std::min(jc, npc);
                if (bMin > bMax) return; // infeasible

                // g_c coefficients
                std::vector<double> gc(bMax + 1, 0.0);
                for (int bc = bMin; bc <= bMax; ++bc)
                    gc[bc] = nCr[nhc][jc - bc] * nCr[npc][bc];

                // Convolve f with gc
                std::vector<double> newf(f.size() + gc.size() - 1, 0.0);
                for (int i = 0; i < (int)f.size(); ++i)
                    for (int j = 0; j < (int)gc.size(); ++j)
                        newf[i + j] += f[i] * gc[j];
                f = std::move(newf);
            }

            // AggrCount = Σ_t f[t] × C(p-t, s-h-t)
            double aggr = 0.0;
            for (int t = 0; t < (int)f.size(); ++t) {
                if (f[t] == 0.0) continue;
                int n = p - t, k = s - h - t;
                if (n >= 0 && k >= 0 && n >= k)
                    aggr += f[t] * nCr[n][k];
            }

            if (aggr > 0) {
                support[tidx] += aggr / rTuples[tidx].mult;
                totalTuplesOnPaths++;
            }
        };
        enumerateMultisets(pathClasses, r, 0, cur, cb);
    }

    auto tStep4 = std::chrono::high_resolution_clock::now();
    auto step4Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep4 - tStep3b).count();

    double totalSupportTuples = 0, totalRCliques = 0;
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        totalSupportTuples += rTuples[i].mult * support[i];
        totalRCliques += rTuples[i].mult;
    }
    std::cout << "  CPI paths used: " << pathsUsed
              << " (region: " << pathsWithRegion << ", no-region: " << pathsWithoutRegion << ")" << std::endl;
    std::cout << "  Tuple-path pairs (CPI only): " << totalTuplesOnPaths << std::endl;
    std::cout << "  r-cliques: " << std::fixed << std::setprecision(0) << totalRCliques << std::endl;
    std::cout << "  Support sum: " << totalSupportTuples << std::endl;
    std::cout << "  CPI counting time: " << step4Ms << " ms" << std::endl;

    // ============================================================
    // Step 5+6: Constrained Path Peeling (Analytical Split)
    // ============================================================
    // Constrained Path = CPI path + per-class (min_piv, max_piv) bounds.
    // When τ is peeled on P̂: subtract old contributions, split P̂ into
    // κ disjoint sub-paths (each a ConstrainedPath), add new contributions.
    // NO s-tuple enumeration. NO inclusion-exclusion. NO BK execution.

    auto tStep5 = std::chrono::high_resolution_clock::now();

    // --- Constrained Path data structure ---
    struct ClassBounds { int nh, np; int minPiv, maxPiv; }; // per-class on a constrained path
    struct CPath {
        int h; // total holds
        std::unordered_map<daf::Size, ClassBounds> classes; // cid -> bounds
        std::vector<daf::Size> tupleIdxs; // alive tuples on this constrained path
        int totalP() const { int t=0; for(auto&[_,b]:classes) t+=b.np; return t; }
    };

    // Compute AggrCount(τ', P̂) / mult(τ') on a constrained path
    // OPTIMIZED: convolve tuple classes (O(r²)) → h[0..r].
    // Non-tuple unconstrained classes: combine into C(Pfree, T-k) lookup.
    // Total: O(r² + κ_constrained × r) instead of O(classes × T).
    auto aggrCountOnCPath = [&](daf::Size tidx, const CPath &cp) -> double {
        auto &key = rTuples[tidx].key;
        int T = s - cp.h;
        if (T < 0) return 0.0;

        // Phase 1: Convolve tuple classes → h[0..r], and sum Pfree
        double h[32]; int hLen = 1; h[0] = 1.0;
        int Pfree = 0;

        // First pass: process tuple classes
        int ki = 0;
        while (ki < (int)key.size()) {
            daf::Size c = key[ki];
            int jc = 1;
            while (ki + jc < (int)key.size() && key[ki + jc] == c) jc++;
            ki += jc;

            auto cit = cp.classes.find(c);
            if (cit == cp.classes.end()) return 0.0;
            auto &cb = cit->second;
            int tMin = std::max(cb.minPiv, std::max(0, jc - cb.nh));
            int tMax = std::min(cb.maxPiv, cb.np);
            if (jc > cb.nh + tMax || tMin > tMax) return 0.0;

            double gc[32]; int gcLen = tMax + 1;
            for (int i = 0; i < gcLen; ++i) gc[i] = 0.0;
            for (int tc = tMin; tc <= tMax; ++tc)
                gc[tc] = nCr[cb.np][tc] * nCr[cb.nh + tc][jc];
            int nL = std::min(hLen + gcLen - 1, (int)r + 1);
            double nh2[32];
            for (int i = 0; i < nL; ++i) nh2[i] = 0.0;
            for (int i = 0; i < hLen; ++i)
                for (int j = 0; j < gcLen && i+j < nL; ++j)
                    nh2[i+j] += h[i] * gc[j];
            hLen = nL;
            for (int i = 0; i < hLen; ++i) h[i] = nh2[i];
        }

        // Phase 2: non-tuple classes
        for (auto &[c, cb] : cp.classes) {
            // Check if c is in tuple (inline, no hash set)
            bool inTuple = false;
            for (auto x : key) if (x == c) { inTuple = true; break; }
            if (inTuple) continue;

            if (cb.minPiv == 0 && cb.maxPiv >= cb.np) {
                // Unconstrained → fold into Pfree
                Pfree += cb.np;
            } else {
                // Constrained non-tuple class: convolve with h
                int tMin = cb.minPiv, tMax = std::min(cb.maxPiv, cb.np);
                if (tMin > tMax) return 0.0;
                double gc[256]; int gcLen = tMax + 1;
                for (int i = 0; i < gcLen; ++i) gc[i] = 0.0;
                for (int tc = tMin; tc <= tMax; ++tc)
                    gc[tc] = nCr[cb.np][tc];
                int nL = std::min(hLen + gcLen - 1, T + 1);
                double nh2[256];
                for (int i = 0; i < nL; ++i) nh2[i] = 0.0;
                for (int i = 0; i < hLen; ++i)
                    for (int j = 0; j < gcLen && i+j < nL; ++j)
                        nh2[i+j] += h[i] * gc[j];
                hLen = nL;
                for (int i = 0; i < hLen; ++i) h[i] = nh2[i];
            }
        }

        // Phase 3: Σ h[k] × C(Pfree, T-k)
        double aggr = 0.0;
        for (int k = 0; k < hLen; ++k) {
            if (h[k] == 0.0) continue;
            int rem = T - k;
            if (rem >= 0 && rem <= Pfree)
                aggr += h[k] * nCr[Pfree][rem];
        }
        return aggr / rTuples[tidx].mult;
    };

    // --- Initialize constrained paths from CPI paths ---
    // Use a pool of CPath objects. Each tuple tracks which cpaths it's on.
    std::vector<CPath> cpaths;
    std::vector<std::vector<daf::Size>> tupleToCPaths(rTuples.size()); // tuple -> cpath IDs

    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        auto &leaf = tree.adj_list[pid];
        if ((int)leaf.size() < (int)s) continue;

        daf::Size pRid = pathRegion[pid]; // region for this path

        CPath cp;
        cp.h = 0;
        for (const auto &node : leaf) {
            daf::Size v = node.v;
            daf::Size cid = effectiveClass(v, pRid);
            if (cid == INVALID) continue;
            if (classSizes[cid] == 0) continue; // merged away
            auto &cb = cp.classes[cid];
            if (node.isPivot) { cb.np++; }
            else { cb.nh++; cp.h++; }
        }
        // Set initial bounds
        for (auto &[cid, cb] : cp.classes) { cb.minPiv = 0; cb.maxPiv = cb.np; }

        // Find tuples on this path (enumerate r-multisets of path's classes)
        std::vector<daf::Size> pathClasses;
        for (auto &[cid, _] : cp.classes) pathClasses.push_back(cid);
        std::sort(pathClasses.begin(), pathClasses.end());

        daf::Size cpid = cpaths.size();
        TupleKey cur; cur.reserve(r);
        bool hasTuples = false;

        std::function<void()> enumCb = [&]() {
            auto it = rTupleIndex.find(cur);
            if (it == rTupleIndex.end()) return;
            daf::Size tidx = it->second;
            // Feasibility check
            std::unordered_map<daf::Size, int> counts;
            for (auto c : cur) counts[c]++;
            for (auto &[c, jc] : counts) {
                auto cit = cp.classes.find(c);
                if (cit == cp.classes.end()) return;
                if (jc > cit->second.nh + cit->second.np) return;
            }
            cp.tupleIdxs.push_back(tidx);
            hasTuples = true;
        };
        enumerateMultisets(pathClasses, r, 0, cur, enumCb);

        if (hasTuples) {
            cpaths.push_back(std::move(cp));
            for (auto tidx : cpaths.back().tupleIdxs)
                tupleToCPaths[tidx].push_back(cpid);
        }
    }

    auto tStep5End = std::chrono::high_resolution_clock::now();
    auto step5Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep5End - tStep5).count();
    std::cout << "  Constrained paths: " << cpaths.size() << std::endl;
    std::cout << "  Index build time: " << step5Ms << " ms" << std::endl;

    // Verify: aggrCountOnCPath matches Step 4 support
    {
        std::vector<double> supCheck(rTuples.size(), 0.0);
        for (daf::Size cpid = 0; cpid < cpaths.size(); ++cpid)
            for (auto tidx : cpaths[cpid].tupleIdxs)
                supCheck[tidx] += aggrCountOnCPath(tidx, cpaths[cpid]);
        double checkSum = 0;
        int mismatches = 0;
        for (daf::Size i = 0; i < rTuples.size(); ++i) {
            checkSum += rTuples[i].mult * supCheck[i];
            if (std::abs(supCheck[i] - support[i]) > 0.5) mismatches++;
        }
        std::cout << "  Support sum (CPath): " << std::fixed << std::setprecision(0) << checkSum << std::endl;
        std::cout << "  CPI vs CPath match: " << (mismatches == 0 ? "PASS" : ("MISMATCH(" + std::to_string(mismatches) + ")")) << std::endl;
    }

    // --- Peeling with Analytical Split (double support + bucket queue) ---
    auto tStep6 = std::chrono::high_resolution_clock::now();

    // Profiling counters
    long long prof_oldContrib_ns = 0, prof_split_ns = 0, prof_newContrib_ns = 0;
    long long prof_delta_ns = 0, prof_overhead_ns = 0;
    long long prof_aggrCalls = 0, prof_splitSubpaths = 0, prof_deadPaths = 0;

    // Keep double support for exact arithmetic. Use floor() for bucket index.
    std::vector<double> dSup = support;
    std::vector<bool> rPeeled(rTuples.size(), false);

    daf::Size maxSup = 0;
    for (auto &sv : dSup) maxSup = std::max(maxSup, (daf::Size)llround(sv));
    const daf::Size BUCKET_CAP = std::min(maxSup + 2, (daf::Size)1000001);

    std::vector<std::vector<daf::Size>> buckets(BUCKET_CAP);
    std::multimap<daf::Size, daf::Size> overflow; // for large support values
    for (daf::Size i = 0; i < rTuples.size(); ++i) {
        daf::Size b = (daf::Size)llround(dSup[i]);
        if (b < BUCKET_CAP) buckets[b].push_back(i);
        else overflow.insert({b, i});
    }

    // coreDist already initialized in Step 0 (may have isolated/small component entries)
    daf::Size numPeeled = 0, currentLevel = 0, coreLevel = 0;
    daf::Size totalSplits = 0;

    while (numPeeled < rTuples.size()) {
        while (currentLevel < BUCKET_CAP && buckets[currentLevel].empty())
            currentLevel++;
        daf::Size idx;
        if (currentLevel < BUCKET_CAP) {
            idx = buckets[currentLevel].back();
            buckets[currentLevel].pop_back();
        } else if (!overflow.empty()) {
            auto it = overflow.begin();
            currentLevel = it->first;
            idx = it->second;
            overflow.erase(it);
        } else break;

        if (rPeeled[idx]) continue;
        daf::Size idxBucket = (daf::Size)llround(dSup[idx]);
        if (idxBucket != currentLevel) continue;

        rPeeled[idx] = true;
        numPeeled++;
        coreLevel = std::max(coreLevel, currentLevel);
        coreDist[coreLevel] += rTuples[idx].mult;

        auto &tauKey = rTuples[idx].key;
        std::unordered_map<daf::Size, int> tauCounts;
        for (auto c : tauKey) tauCounts[c]++;

        auto cpathIds = tupleToCPaths[idx]; // copy
        for (daf::Size cpid : cpathIds) {
            auto &cp = cpaths[cpid];
            if (cp.tupleIdxs.empty()) continue;

            auto _t0 = std::chrono::high_resolution_clock::now();

            // Separate alive tuples into: affected (share class with τ) and unaffected
            std::unordered_set<daf::Size> tauClassSet(rTuples[idx].key.begin(), rTuples[idx].key.end());
            std::vector<daf::Size> affectedTuples, unaffectedTuples;
            for (auto tidx : cp.tupleIdxs) {
                if (rPeeled[tidx]) continue;
                bool shares = false;
                for (auto c : rTuples[tidx].key)
                    if (tauClassSet.count(c)) { shares = true; break; }
                if (shares) affectedTuples.push_back(tidx);
                else unaffectedTuples.push_back(tidx);
            }

            // Compute old contributions (only for AFFECTED tuples)
            std::unordered_map<daf::Size, double> oldContrib;
            for (auto tidx : affectedTuples) {
                oldContrib[tidx] = aggrCountOnCPath(tidx, cp);
                prof_aggrCalls++;
            }

            auto _t1 = std::chrono::high_resolution_clock::now();
            prof_oldContrib_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(_t1-_t0).count();

            // Find active classes
            struct ActiveClass { daf::Size cid; int mc; };
            std::vector<ActiveClass> active;
            for (auto &[c, jc] : tauCounts) {
                auto cit = cp.classes.find(c);
                if (cit == cp.classes.end()) continue;
                int mc = std::max(0, jc - cit->second.nh);
                if (mc > cit->second.minPiv)
                    active.push_back({c, mc});
            }

            for (auto tidx : cp.tupleIdxs) {
                auto &vec = tupleToCPaths[tidx];
                vec.erase(std::remove(vec.begin(), vec.end(), cpid), vec.end());
            }
            cp.tupleIdxs.clear();

            auto _t2 = std::chrono::high_resolution_clock::now();
            prof_overhead_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(_t2-_t1).count();

            // Compute new contributions from sub-paths
            std::unordered_map<daf::Size, double> newContrib;
            if (active.empty()) {
                prof_deadPaths++;
            } else {
                for (int i = 0; i < (int)active.size(); ++i) {
                    CPath subCp = cp;
                    bool feasible = true;

                    for (int j = 0; j < i; ++j) {
                        auto &cb = subCp.classes[active[j].cid];
                        cb.minPiv = std::max(cb.minPiv, active[j].mc);
                        if (cb.minPiv > cb.maxPiv) { feasible = false; break; }
                    }
                    if (!feasible) continue;

                    {
                        auto &cb = subCp.classes[active[i].cid];
                        cb.maxPiv = std::min(cb.maxPiv, active[i].mc - 1);
                        if (cb.minPiv > cb.maxPiv) continue;
                    }

                    int minTotal = 0, maxTotal = 0;
                    for (auto &[_, cb] : subCp.classes) {
                        minTotal += cb.minPiv;
                        maxTotal += cb.maxPiv;
                    }
                    if ((s - subCp.h) < minTotal || (s - subCp.h) > maxTotal) continue;

                    subCp.tupleIdxs.clear();
                    // Only compute aggrCount for AFFECTED tuples (share class with τ)
                    for (auto tidx : affectedTuples) {
                        double c = aggrCountOnCPath(tidx, subCp);
                        prof_aggrCalls++;
                        if (c > 1e-9) {
                            subCp.tupleIdxs.push_back(tidx);
                            newContrib[tidx] += c;
                        }
                    }
                    // Unaffected tuples: same contribution as parent → add all
                    for (auto tidx : unaffectedTuples)
                        subCp.tupleIdxs.push_back(tidx);
                    if (subCp.tupleIdxs.empty()) continue;

                    daf::Size newCpid = cpaths.size();
                    cpaths.push_back(std::move(subCp));
                    totalSplits++;
                    prof_splitSubpaths++;

                    for (auto tidx : cpaths.back().tupleIdxs)
                        tupleToCPaths[tidx].push_back(newCpid);
                }
            }

            auto _t3 = std::chrono::high_resolution_clock::now();
            prof_split_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(_t3-_t2).count();

            // Apply net delta in double, re-bucket
            for (auto &[tidx, oc] : oldContrib) {
                double nc = newContrib.count(tidx) ? newContrib[tidx] : 0.0;
                dSup[tidx] += (nc - oc);
                if (dSup[tidx] < -0.5) dSup[tidx] = 0;
                daf::Size newBucket = (daf::Size)llround(dSup[tidx]);
                if (newBucket < BUCKET_CAP) {
                    buckets[newBucket].push_back(tidx);
                    if (newBucket < currentLevel) currentLevel = newBucket;
                } else {
                    overflow.insert({newBucket, tidx});
                }
            }
        }
    }

    auto tStep6End = std::chrono::high_resolution_clock::now();
    auto step6Ms = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStep6).count();
    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(tStep6End - tStart).count();

    daf::Size maxCore = coreDist.empty() ? 0 : coreDist.rbegin()->first;
    std::cout << "\n  --- Cascade Peeling (Analytical Split) ---" << std::endl;
    std::cout << "  Peeled: " << numPeeled << " / " << rTuples.size() << std::endl;
    std::cout << "  Total splits: " << totalSplits << " (sub-paths: " << prof_splitSubpaths
              << ", dead: " << prof_deadPaths << ")" << std::endl;
    std::cout << "  Final cpaths: " << cpaths.size() << std::endl;
    std::cout << "  AggrCount calls: " << prof_aggrCalls << std::endl;
    std::cout << "  Time breakdown:" << std::endl;
    std::cout << "    oldContrib: " << prof_oldContrib_ns/1000000 << " ms" << std::endl;
    std::cout << "    split+new:  " << prof_split_ns/1000000 << " ms" << std::endl;
    std::cout << "    overhead:   " << prof_overhead_ns/1000000 << " ms" << std::endl;
    std::cout << "  Max core: " << maxCore << std::endl;
    for (auto &[core, cnt] : coreDist)
        std::cout << "  core=" << core << " count=" << cnt << std::endl;
    std::cout << "  Peeling time: " << step6Ms << " ms" << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;

    return {};
}
