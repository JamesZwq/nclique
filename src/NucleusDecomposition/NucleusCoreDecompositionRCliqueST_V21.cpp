//
// ST V21: Sequential d=1 with cache propagation — ZERO BK recursion
//
// For each removal q from leaf P:
//   1. Analytical delta: Δ(q') = -nCr(a-|piv(q'∪q)|, need-|piv(q'∪q)|)
//   2. Theorem-1 deterministic split: t ≤ r subpaths
//   3. Cache propagation: parent entries → children via containment
//
// Cost: O(m × r² × C(n,r))  vs  pathSplit: O(Π|F_i| × C(n',r) × r)
//

#include "NCliqueCoreDecomposition.h"
#include <chrono>
#include <algorithm>
#include <cstring>
#include <set>

#include "dataStruct/CliqueHashMap.h"
#include "graph/DynamicGraphSet.h"

extern double nCr[1001][401];

namespace V21 {
struct CacheEntry {
    daf::Size cliqueId;
    uint64_t pivotMask; // leaf-local pivot positions as bitmask
};
} // namespace V21

std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V21(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    long long dur_init=0, dur_pop=0, dur_intersect=0, dur_delta=0, dur_split=0;
    long long cntDeath=0, cntSplit=0;
    auto T0 = std::chrono::high_resolution_clock::now();

    // ========== INIT ==========
    StaticCliqueIndex localIndex(r);
    StaticCliqueIndex &CI = prebuiltIndex ? *prebuiltIndex : localIndex;
    if (!prebuiltIndex)
        daf::timeCount("clique Index build", [&]() { localIndex.build(tree, edgeGraph.adj_list.size()); });

    std::vector<double> cnt;
    daf::timeCount("countingPerRClique (V21)", [&]() {
        cnt.assign(CI.size(), 0.0);
        for (const auto &leaf : tree.adj_list) {
            if (leaf.size() < r) continue;
            daf::CliqueSize pC=0, kC=0;
            for (const auto &v : leaf) { if(v.isPivot)pC++;else kC++; }
            int need = s-(int)kC;
            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rc) {
                daf::CliqueSize sp=0;
                for (const auto &v : rc) if(v.isPivot)sp++;
                if (sp<=need) { int R=(int)pC-(int)sp,C=need-(int)sp;
                    if(R>=0&&R<1001&&C>=0&&C<401) cnt[CI.byClique(rc)]+=nCr[R][C]; }
                return true;
            });
        }
    });

    const daf::Size nCl = CI.size();
    daf::log_memory("index+counting");
    std::vector<double> core(nCl, 0);
    daf::StaticVector<bool> inHeap(nCl); inHeap.resize(nCl);
    memset(inHeap.getData(), true, nCl*sizeof(bool));

    // Per-leaf cache (lazy)
    std::vector<std::vector<V21::CacheEntry>> leafCache(tree.adj_list.size());

    // Build cache for a leaf
    auto buildCache = [&](daf::Size leafId) {
        const auto &leaf = tree.adj_list[leafId];
        auto &cache = leafCache[leafId];
        if (!cache.empty() || leaf.size() < r) return;
        daf::StaticVector<daf::Size> &M = daf::vListMap;
        int n=(int)leaf.size();
        for (int i=0;i<n;i++) M[leaf[i].v]=(daf::Size)i;
        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rc) {
            uint64_t mask=0;
            for (const auto &v : rc) if(v.isPivot) mask|=(1ULL<<M[v.v]);
            cache.push_back({CI.byClique(rc), mask});
            return true;
        });
    };

    // ========== Bucket+Set PQ ==========
    constexpr double BT=5e6;
    double rawMax=0;
    for (daf::Size i=0;i<nCl;i++) rawMax=std::max(rawMax,cnt[i]);
    int maxB=(int)std::min(rawMax,BT);
    std::vector<std::vector<daf::Size>> bkts(maxB+2);
    std::set<std::pair<double,daf::Size>> ovf;
    std::vector<int> bof(nCl,-1);
    std::vector<daf::Size> pib(nCl);
    std::vector<double> osv(nCl,-1);
    for (daf::Size i=0;i<nCl;i++) {
        if(cnt[i]<=BT){int b=(int)cnt[i];bof[i]=b;pib[i]=bkts[b].size();bkts[b].push_back(i);}
        else{ovf.insert({cnt[i],i});osv[i]=cnt[i];}
    }
    int curB=0; daf::Size rem=nCl;

    auto bmove=[&](daf::Size id){
        if(!inHeap[id])return;
        double v=std::max(0.0,cnt[id]); int ob=bof[id];
        if(ob==-1)ovf.erase({osv[id],id});
        if(v<=BT){int nb=(int)v;if(ob>=0&&nb==ob)return;
            if(ob>=0){auto&bk=bkts[ob];auto p=pib[id];if(p<bk.size()-1){auto l=bk.back();bk[p]=l;pib[l]=p;}bk.pop_back();}
            bof[id]=nb;pib[id]=bkts[nb].size();bkts[nb].push_back(id);if(nb<curB)curB=nb;
        }else{if(ob>=0){auto&bk=bkts[ob];auto p=pib[id];if(p<bk.size()){if(p<bk.size()-1){auto l=bk.back();bk[p]=l;pib[l]=p;}bk.pop_back();}}
            ovf.insert({v,id});osv[id]=v;bof[id]=-1;}
    };
    auto drain=[&](){while(!ovf.empty()){auto id=ovf.begin()->second;
        if(!inHeap[id]){ovf.erase(ovf.begin());continue;}
        if(cnt[id]<=BT){ovf.erase(ovf.begin());int b=(int)cnt[id];bof[id]=b;pib[id]=bkts[b].size();bkts[b].push_back(id);}else break;}};

    dur_init=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-T0).count();

    std::cout<<"=========================begin (r>=3 ST_V21)==========================="<<std::endl;
    double minCore=0; long long iters=0;
    std::vector<daf::Size> popIds;
    popIds.reserve(nCl);

    while (rem>0) {
        auto t0=std::chrono::high_resolution_clock::now();
        popIds.clear();
        drain();
        while(curB<(int)bkts.size()&&bkts[curB].empty())curB++;
        if(curB>=(int)bkts.size()){
            if(!ovf.empty()){while(!ovf.empty()){
                auto id=ovf.begin()->second;ovf.erase(ovf.begin());
                if(!inHeap[id])continue;
                minCore=std::max(cnt[id],minCore);inHeap[id]=false;popIds.push_back(id);core[id]=minCore;rem--;
                while(!ovf.empty()){auto nx=ovf.begin()->second;
                    if(!inHeap[nx]){ovf.erase(ovf.begin());continue;}
                    if(cnt[nx]<=minCore){ovf.erase(ovf.begin());inHeap[nx]=false;popIds.push_back(nx);core[nx]=minCore;rem--;}else break;}
                break;}if(popIds.empty())break;goto pdone;}
            break;}
        minCore=std::max((double)curB,minCore);
        while(curB<(int)bkts.size()&&!bkts[curB].empty()&&curB<=(int)minCore){
            while(!bkts[curB].empty()){auto id=bkts[curB].back();bkts[curB].pop_back();
                inHeap[id]=false;popIds.push_back(id);core[id]=minCore;rem--;}
            if(curB+1<(int)bkts.size()&&!bkts[curB+1].empty()&&(curB+1)<=(int)minCore)curB++;else break;}
        pdone:
        dur_pop+=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-t0).count();
        if(rem==0)break;
        iters++;

        // ===== Process each removed r-clique SEQUENTIALLY (d=1) =====
        for (auto rmId : popIds) {
            auto rmVerts = CI.byId(rmId);

            // Find leaves containing this r-clique
            auto t1=std::chrono::high_resolution_clock::now();
            std::vector<daf::Size> leaves;
            daf::intersect_dense_sets_multi(rmVerts, treeGraphV.adj_list,
                [&](const TreeGraphNode &u) { leaves.push_back(u.v); });
            dur_intersect+=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-t1).count();

            for (auto leafId : leaves) {
                const auto &leaf = tree.adj_list[leafId];
                if (leaf.empty()) continue;
                int n=(int)leaf.size();
                if (n<(int)s) continue;

                daf::CliqueSize pivotC=0, keepC=0;
                for (const auto &v : leaf) { if(v.isPivot)pivotC++;else keepC++; }
                int need = s-(int)keepC;
                int a=(int)pivotC;

                // Map vertices to leaf-local positions
                daf::StaticVector<daf::Size> &M = daf::vListMap;
                for (int i=0;i<n;i++) M[leaf[i].v]=(daf::Size)i;

                // Build removed r-clique's pivot mask on this leaf
                uint64_t rmMask=0;
                int rmPivCnt=0;
                int rmPivPos[8];
                for (auto v : rmVerts) {
                    daf::Size pos=M[v];
                    if (pos<(daf::Size)n && leaf[pos].isPivot) {
                        rmMask|=(1ULL<<pos);
                        rmPivPos[rmPivCnt++]=(int)pos;
                    }
                }

                if (rmPivCnt==0) {
                    // All rm vertices are holds → leaf dies
                    cntDeath++;
                    auto td=std::chrono::high_resolution_clock::now();
                    // Ensure cache exists for subtraction
                    if (leafId>=leafCache.size()) leafCache.resize(leafId+1);
                    buildCache(leafId);
                    for (const auto &e : leafCache[leafId]) {
                        if(!inHeap[e.cliqueId])continue;
                        int p=__builtin_popcountll(e.pivotMask);
                        cnt[e.cliqueId]-=nCr[a-p][need-p];
                        if(cnt[e.cliqueId]<0)cnt[e.cliqueId]=0;
                        bmove(e.cliqueId);
                    }
                    for (const auto &v : leaf) treeGraphV.removeNbr(v.v,{leafId,v.isPivot});
                    tree.removeNode(leafId);
                    leafCache[leafId].clear(); leafCache[leafId].shrink_to_fit();
                    dur_delta+=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-td).count();
                    continue;
                }

                cntSplit++;
                auto td=std::chrono::high_resolution_clock::now();

                // Ensure cache exists
                if (leafId>=leafCache.size()) leafCache.resize(leafId+1);
                buildCache(leafId);
                auto &parentCache = leafCache[leafId];

                // ===== STEP 1: Analytical d=1 delta =====
                for (const auto &e : parentCache) {
                    if (!inHeap[e.cliqueId]) continue;
                    uint64_t combined = e.pivotMask | rmMask;
                    int u = __builtin_popcountll(combined);
                    int R = a-u, C = need-u;
                    if (R>=0 && C>=0 && R<1001 && C<401) {
                        cnt[e.cliqueId] -= nCr[R][C];
                        if (cnt[e.cliqueId]<0) cnt[e.cliqueId]=0;
                        bmove(e.cliqueId);
                    }
                }
                dur_delta+=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-td).count();

                // ===== STEP 2: Theorem-1 split + cache propagation =====
                auto ts=std::chrono::high_resolution_clock::now();

                // Remove old leaf from treeGraphV
                for (const auto &v : leaf) treeGraphV.removeNbr(v.v,{leafId,v.isPivot});

                // Sort pivot positions for deterministic ordering
                std::sort(rmPivPos, rmPivPos+rmPivCnt);

                // Generate t = rmPivCnt subpaths
                for (int i=0; i<rmPivCnt; i++) {
                    // Subpath i: promote rmPivPos[0..i-1] to hold, remove rmPivPos[i]
                    int excludePos = rmPivPos[i];

                    // Build subpath vertex list + pivot bitmask for containment
                    std::vector<TreeGraphNode> subLeaf;
                    subLeaf.reserve(n-1);
                    uint64_t childPivMask=0; // which PARENT positions are pivots in child
                    daf::CliqueSize newPivC=0, newKeepC=0;

                    for (int p=0; p<n; p++) {
                        if (p==excludePos) continue;
                        bool isPiv = leaf[p].isPivot;
                        // Promote positions [0..i-1] to hold
                        for (int j=0; j<i; j++) {
                            if (p==rmPivPos[j]) { isPiv=false; break; }
                        }
                        subLeaf.push_back({leaf[p].v, isPiv});
                        if (isPiv) { childPivMask|=(1ULL<<p); newPivC++; }
                        else newKeepC++;
                    }

                    // Prune infeasible subpaths
                    int newNeed = s-newKeepC;
                    if ((int)(newKeepC+newPivC)<(int)s || (int)newKeepC>(int)s || (int)newPivC<newNeed)
                        continue;

                    auto newId = tree.addNode(subLeaf);
                    for (const auto &v : tree.adj_list[newId])
                        treeGraphV.addNbr(v.v, {newId, v.isPivot});

                    // Cache propagation: parent → child via containment
                    if (newId>=leafCache.size()) leafCache.resize(newId+1);
                    auto &childCache = leafCache[newId];
                    childCache.clear();

                    // Build position remap: parent pos → child pos
                    uint8_t remap[400];
                    memset(remap, 255, sizeof(remap));
                    uint8_t ci=0;
                    for (int p=0; p<n; p++) {
                        if (p==excludePos) continue;
                        remap[p]=ci++;
                    }

                    for (const auto &e : parentCache) {
                        // Check: all pivot positions of this r-clique must be in childPivMask or promoted holds
                        // Equivalently: all positions in e.pivotMask must NOT be excludePos,
                        // and the r-clique's vertices must all be in the child
                        // Since child = parent minus excludePos, check if excludePos is NOT one of e's positions
                        // Actually: e.pivotMask has the parent-local positions of pivots of this r-clique.
                        // The r-clique is in the child iff none of its vertex positions == excludePos.
                        // But e.pivotMask only tracks PIVOT positions. What about hold positions?
                        // An r-clique consists of both pivot and hold vertices. The hold positions are always in the child
                        // (holds are never removed). The pivot positions: check if they're all still present.
                        // A pivot position is present in child iff it's not excludePos.
                        // Note: positions promoted to hold (rmPivPos[0..i-1]) are still in child, just as holds.

                        // So: r-clique is in child iff excludePos is NOT one of its VERTEX positions.
                        // e.pivotMask has the r-clique's PIVOT positions. But the r-clique might also have
                        // holds at other positions. We need to check ALL positions.
                        // Problem: we only store pivotMask, not all positions!

                        // Fix: the r-clique is absent from child ONLY if one of its vertices is at excludePos.
                        // But we don't know which positions the holds are at (only pivotMask stored).
                        // We need to store the full position set, or check differently.

                        // Alternative: check if the r-clique shares a vertex with the excluded position.
                        // The vertex at excludePos is leaf[excludePos].v. The r-clique's vertices can be
                        // recovered from CI.byId(e.cliqueId). Check if any vertex == leaf[excludePos].v.

                        // This is O(r) per entry. Total: O(|cache| × r) per subpath.

                        auto rc = CI.byId(e.cliqueId);
                        bool inChild = true;
                        for (auto v : rc) {
                            if (M[v]==(daf::Size)excludePos) { inChild=false; break; }
                        }
                        if (!inChild) continue;

                        // Build new pivot mask in CHILD-local positions
                        uint64_t newMask=0;
                        for (auto v : rc) {
                            daf::Size parentPos = M[v];
                            if (parentPos<(daf::Size)n && remap[parentPos]!=255) {
                                // Check if this position is a pivot in the child
                                if (childPivMask & (1ULL<<parentPos))
                                    newMask |= (1ULL << remap[parentPos]);
                            }
                        }
                        childCache.push_back({e.cliqueId, newMask});
                    }
                }

                tree.removeNode(leafId);
                parentCache.clear(); parentCache.shrink_to_fit();

                dur_split+=std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-ts).count();
            }
        }
    }

    auto elapsed=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-T0).count();
    std::cout<<"time: "<<elapsed<<" ms"<<std::endl;
    std::cout<<"Time Breakdown (ms):"<<std::endl;
    std::cout<<"  Init:      "<<dur_init/1e6<<std::endl;
    std::cout<<"  Pop:       "<<dur_pop/1e6<<std::endl;
    std::cout<<"  Intersect: "<<dur_intersect/1e6<<std::endl;
    std::cout<<"  Delta:     "<<dur_delta/1e6<<std::endl;
    std::cout<<"  Split:     "<<dur_split/1e6<<std::endl;
    std::cout<<"  Cases: Death="<<cntDeath<<" Split="<<cntSplit<<" iters="<<iters<<std::endl;

    std::vector<std::pair<std::vector<daf::Size>,double>> out;
    out.reserve(nCl);
    for (daf::Size i=0;i<nCl;i++){auto c=CI.byId(i);out.emplace_back(std::vector<daf::Size>(c.begin(),c.end()),core[i]);}
    return out;
}
