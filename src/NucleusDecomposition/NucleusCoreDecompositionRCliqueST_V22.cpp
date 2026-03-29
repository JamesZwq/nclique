//
// ST V22: Lazy Bipartite Peeling — ZERO pathSplit
//
// No pre-enumeration of s-cliques. During peeling:
//   When r-clique q is removed from leaf P:
//     1. Enumerate s-cliques containing q on P: C(a-t, need-t) cliques
//     2. For each: check alive (via hash set), if alive → destroy
//     3. On destruction: find other r-cliques (via leaf cache), decrement support
//
// Each s-clique is destroyed exactly once. Amortized cost: O(|destroyed| × C(s,r)).
// No pathSplit, no tree mutation, no BK recursion.
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <cstring>
#include <set>
#include <unordered_set>

#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace V22 {
// Compact s-clique identifier: leafId + pivot bitmask
struct SCliqueKey {
    daf::Size leafId;
    uint64_t pivMask; // which pivots are chosen (among leaf's pivot positions)
    bool operator==(const SCliqueKey &o) const { return leafId == o.leafId && pivMask == o.pivMask; }
};
struct SCliqueHash {
    size_t operator()(const SCliqueKey &k) const {
        return std::hash<uint64_t>()(k.pivMask) ^ (std::hash<daf::Size>()(k.leafId) << 32);
    }
};

struct CacheEntry {
    daf::Size cliqueId;
    uint64_t pivMask; // pivot positions among leaf's pivots
};
} // namespace V22

std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V22(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    long long dur_init = 0, dur_pop = 0, dur_intersect = 0, dur_peel = 0;
    long long cntDestroyed = 0;
    auto T0 = std::chrono::high_resolution_clock::now();

    // ========== INIT ==========
    StaticCliqueIndex localCI(r);
    StaticCliqueIndex &CI = prebuiltIndex ? *prebuiltIndex : localCI;
    if (!prebuiltIndex)
        daf::timeCount("clique Index build (V22)", [&]() { localCI.build(tree, edgeGraph.adj_list.size()); });
    const daf::Size nCl = CI.size();

    // Count support via standard counting
    std::vector<double> support;
    daf::timeCount("countingPerRClique (V22)", [&]() {
        support.assign(nCl, 0.0);
        for (const auto &leaf : tree.adj_list) {
            if (leaf.size() < r) continue;
            daf::CliqueSize pC = 0, kC = 0;
            for (const auto &v : leaf) { if (v.isPivot) pC++; else kC++; }
            int need = s - (int)kC;
            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rc) {
                daf::CliqueSize sp = 0;
                for (const auto &v : rc) if (v.isPivot) sp++;
                if (sp <= need) {
                    int R = (int)pC - (int)sp, C = need - (int)sp;
                    if (R >= 0 && R < 1001 && C >= 0 && C < 401)
                        support[CI.byClique(rc)] += nCr[R][C];
                }
                return true;
            });
        }
    });

    // Build per-leaf r-clique cache (lazy — build on first access)
    std::vector<std::vector<V22::CacheEntry>> leafCache(tree.adj_list.size());
    // Per-leaf pivot position list (for mapping)
    std::vector<std::vector<int>> leafPivPos(tree.adj_list.size());

    auto ensureCache = [&](daf::Size leafId) {
        if (!leafCache[leafId].empty()) return;
        const auto &leaf = tree.adj_list[leafId];
        int n = (int)leaf.size();
        if (n < (int)r) return;
        // Build pivot position list
        auto &pivPos = leafPivPos[leafId];
        pivPos.clear();
        for (int i = 0; i < n; i++)
            if (leaf[i].isPivot) pivPos.push_back(i);
        // Map vertices
        daf::StaticVector<daf::Size> &M = daf::vListMap;
        for (int i = 0; i < n; i++) M[leaf[i].v] = (daf::Size)i;
        // Build cache using lookupRaw (sorted vertex array)
        auto &cache = leafCache[leafId];
        daf::Size vertBuf[8];
        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rc) {
            uint64_t mask = 0;
            for (int j = 0; j < (int)r; j++) {
                vertBuf[j] = rc[j].v;
                if (rc[j].isPivot) {
                    daf::Size pos = M[rc[j].v];
                    for (int pi = 0; pi < (int)pivPos.size(); pi++) {
                        if ((daf::Size)pivPos[pi] == pos) { mask |= (1ULL << pi); break; }
                    }
                }
            }
            cache.push_back({CI.lookupRaw(vertBuf), mask});
            return true;
        });
    };

    // Dead s-clique tracker
    std::unordered_set<uint64_t> deadSCliques; // packed key: (leafId << 32) | pivMask
    // For leaves with >32 pivots, use a separate set with full SCliqueKey
    std::unordered_set<V22::SCliqueKey, V22::SCliqueHash> deadSCliquesLarge;

    auto packKey = [](daf::Size leafId, uint64_t pivMask) -> uint64_t {
        return ((uint64_t)leafId << 32) | (pivMask & 0xFFFFFFFF);
    };

    auto isDead = [&](daf::Size leafId, uint64_t pivMask, int pivotC) -> bool {
        if (pivotC <= 32)
            return deadSCliques.count(packKey(leafId, pivMask)) > 0;
        else
            return deadSCliquesLarge.count({leafId, pivMask}) > 0;
    };

    auto markDead = [&](daf::Size leafId, uint64_t pivMask, int pivotC) {
        if (pivotC <= 32)
            deadSCliques.insert(packKey(leafId, pivMask));
        else
            deadSCliquesLarge.insert({leafId, pivMask});
    };

    std::vector<double> core(nCl, 0);
    std::vector<bool> inHeap(nCl, true);

    // ========== Bucket PQ ==========
    constexpr double BT = 5e6;
    double rawMax = 0;
    for (daf::Size i = 0; i < nCl; i++) rawMax = std::max(rawMax, support[i]);
    int maxB = (int)std::min(rawMax, BT);
    std::vector<std::vector<daf::Size>> bkts(maxB + 2);
    std::set<std::pair<double, daf::Size>> ovf;
    std::vector<int> bof(nCl, -1);
    std::vector<daf::Size> pib(nCl);
    std::vector<double> osv(nCl, -1);
    for (daf::Size i = 0; i < nCl; i++) {
        if (support[i] <= BT) { int b = (int)support[i]; bof[i] = b; pib[i] = bkts[b].size(); bkts[b].push_back(i); }
        else { ovf.insert({support[i], i}); osv[i] = support[i]; }
    }
    int curB = 0; daf::Size rem = nCl;
    auto bmove = [&](daf::Size id) {
        if (!inHeap[id]) return;
        double v = std::max(0.0, support[id]); int ob = bof[id];
        if (ob == -1) ovf.erase({osv[id], id});
        if (v <= BT) { int nb = (int)v; if (ob >= 0 && nb == ob) return;
            if (ob >= 0) { auto &bk = bkts[ob]; auto p = pib[id]; if (p < bk.size() - 1) { auto l = bk.back(); bk[p] = l; pib[l] = p; } bk.pop_back(); }
            bof[id] = nb; pib[id] = bkts[nb].size(); bkts[nb].push_back(id); if (nb < curB) curB = nb;
        } else { if (ob >= 0) { auto &bk = bkts[ob]; auto p = pib[id]; if (p < bk.size()) { if (p < bk.size() - 1) { auto l = bk.back(); bk[p] = l; pib[l] = p; } bk.pop_back(); } }
            ovf.insert({v, id}); osv[id] = v; bof[id] = -1; }
    };
    auto drain = [&]() { while (!ovf.empty()) { auto id = ovf.begin()->second;
        if (!inHeap[id]) { ovf.erase(ovf.begin()); continue; }
        if (support[id] <= BT) { ovf.erase(ovf.begin()); int b = (int)support[id]; bof[id] = b; pib[id] = bkts[b].size(); bkts[b].push_back(id); } else break; } };

    dur_init = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - T0).count();
    daf::log_memory("V22 init");

    std::cout << "=========================begin (r>=3 ST_V22)===========================" << std::endl;
    double minCore = 0; long long iters = 0;
    std::vector<daf::Size> popIds;

    while (rem > 0) {
        auto t0 = std::chrono::high_resolution_clock::now();
        popIds.clear();
        drain();
        while (curB < (int)bkts.size() && bkts[curB].empty()) curB++;
        if (curB >= (int)bkts.size()) {
            if (!ovf.empty()) { while (!ovf.empty()) {
                auto id = ovf.begin()->second; ovf.erase(ovf.begin());
                if (!inHeap[id]) continue;
                minCore = std::max(support[id], minCore); inHeap[id] = false;
                popIds.push_back(id); core[id] = minCore; rem--;
                while (!ovf.empty()) { auto nx = ovf.begin()->second;
                    if (!inHeap[nx]) { ovf.erase(ovf.begin()); continue; }
                    if (support[nx] <= minCore) { ovf.erase(ovf.begin()); inHeap[nx] = false; popIds.push_back(nx); core[nx] = minCore; rem--; } else break; }
                break; } if (popIds.empty()) break; goto pd; }
            break; }
        minCore = std::max((double)curB, minCore);
        while (curB < (int)bkts.size() && !bkts[curB].empty() && curB <= (int)minCore) {
            while (!bkts[curB].empty()) { auto id = bkts[curB].back(); bkts[curB].pop_back();
                inHeap[id] = false; popIds.push_back(id); core[id] = minCore; rem--; }
            if (curB + 1 < (int)bkts.size() && !bkts[curB + 1].empty() && (curB + 1) <= (int)minCore) curB++; else break; }
        pd:
        dur_pop += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t0).count();
        if (rem == 0) break;
        iters++;

        // ===== Lazy bipartite peeling =====
        for (auto rmId : popIds) {
            auto rmVerts = CI.byId(rmId);

            // Find leaves containing this r-clique
            auto t1 = std::chrono::high_resolution_clock::now();
            std::vector<daf::Size> leaves;
            daf::intersect_dense_sets_multi(rmVerts, treeGraphV.adj_list,
                [&](const TreeGraphNode &u) { leaves.push_back(u.v); });
            dur_intersect += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t1).count();

            auto tp = std::chrono::high_resolution_clock::now();
            for (auto leafId : leaves) {
                const auto &leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;
                int n = (int)leaf.size();
                if (n < (int)s) continue;

                daf::CliqueSize pivotC = 0, keepC = 0;
                for (const auto &v : leaf) { if (v.isPivot) pivotC++; else keepC++; }
                int need = s - (int)keepC;
                if (need < 0 || need > (int)pivotC) continue;

                // Ensure cache + mapping + mask buckets
                if (leafId >= leafCache.size()) { leafCache.resize(leafId + 1); leafPivPos.resize(leafId + 1); }
                ensureCache(leafId);
                const auto &cache = leafCache[leafId];
                const auto &pivPos = leafPivPos[leafId];

                // Build mask→entries index for fast submask lookup
                // (rebuilt per leaf access — small cost vs scanning full cache)
                std::unordered_map<uint64_t, std::vector<daf::Size>> maskIdx;
                for (size_t ci = 0; ci < cache.size(); ci++)
                    maskIdx[cache[ci].pivMask].push_back(cache[ci].cliqueId);

                // Build removed r-clique's pivot mask on this leaf
                daf::StaticVector<daf::Size> &M = daf::vListMap;
                for (int i = 0; i < n; i++) M[leaf[i].v] = (daf::Size)i;

                uint64_t rmPivMask = 0;
                int rmPivCnt = 0;
                for (auto v : rmVerts) {
                    daf::Size pos = M[v];
                    if (pos < (daf::Size)n && leaf[pos].isPivot) {
                        for (int pi = 0; pi < (int)pivPos.size(); pi++)
                            if ((daf::Size)pivPos[pi] == pos) { rmPivMask |= (1ULL << pi); rmPivCnt++; break; }
                    }
                }

                // Enumerate s-cliques containing rm on this leaf:
                // Choose (need - rmPivCnt) more pivots from the (pivotC - rmPivCnt) remaining pivots
                int extraNeed = need - rmPivCnt;
                int extraAvail = (int)pivotC - rmPivCnt;
                if (extraNeed < 0 || extraNeed > extraAvail) continue;

                // Collect non-rm pivot indices
                std::vector<int> freePivIdx;
                freePivIdx.reserve(extraAvail);
                for (int pi = 0; pi < (int)pivotC; pi++)
                    if (!(rmPivMask & (1ULL << pi))) freePivIdx.push_back(pi);

                // Enumerate C(extraAvail, extraNeed) combinations
                std::vector<int> combo(extraNeed);
                for (int i = 0; i < extraNeed; i++) combo[i] = i;

                auto processCombo = [&]() {
                    // Build s-clique pivot mask = rmPivMask | chosen extras
                    uint64_t sMask = rmPivMask;
                    for (int i = 0; i < extraNeed; i++) sMask |= (1ULL << freePivIdx[combo[i]]);

                    // Check alive
                    if (isDead(leafId, sMask, (int)pivotC)) return;

                    // Mark dead
                    markDead(leafId, sMask, (int)pivotC);
                    cntDestroyed++;

                    // Find r-cliques via submask enumeration (O(2^need) lookups vs O(|cache|) scan)
                    uint64_t sub = sMask;
                    do {
                        auto it = maskIdx.find(sub);
                        if (it != maskIdx.end()) {
                            for (auto cid : it->second) {
                                if (cid == rmId || !inHeap[cid]) continue;
                                support[cid] -= 1.0;
                                if (support[cid] < 0) support[cid] = 0;
                                bmove(cid);
                            }
                        }
                        if (sub == 0) break;
                        sub = (sub - 1) & sMask;
                    } while (true);
                };

                if (extraNeed == 0) {
                    processCombo();
                } else {
                    while (true) {
                        processCombo();
                        int i = extraNeed - 1;
                        while (i >= 0 && combo[i] == extraAvail - extraNeed + i) i--;
                        if (i < 0) break;
                        combo[i]++;
                        for (int j = i + 1; j < extraNeed; j++) combo[j] = combo[j - 1] + 1;
                    }
                }
            }
            dur_peel += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - tp).count();
        }
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - T0).count();
    std::cout << "time: " << elapsed << " ms" << std::endl;
    std::cout << "Time Breakdown (ms):" << std::endl;
    std::cout << "  Init:      " << dur_init / 1e6 << std::endl;
    std::cout << "  Pop:       " << dur_pop / 1e6 << std::endl;
    std::cout << "  Intersect: " << dur_intersect / 1e6 << std::endl;
    std::cout << "  Peel:      " << dur_peel / 1e6 << std::endl;
    std::cout << "  Destroyed: " << cntDestroyed << ", iters: " << iters << std::endl;

    std::vector<std::pair<std::vector<daf::Size>, double>> out;
    out.reserve(nCl);
    for (daf::Size i = 0; i < nCl; i++) { auto c = CI.byId(i); out.emplace_back(std::vector<daf::Size>(c.begin(), c.end()), core[i]); }
    return out;
}
