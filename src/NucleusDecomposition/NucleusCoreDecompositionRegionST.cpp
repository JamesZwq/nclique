//
// V4: Unified MC + Tuple Peeling
//
// Single loop interleaving MC deletion (private vertex theorem)
// with shared tuple peeling (IE support).
//

#include "NCliqueCoreDecomposition.h"
#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <unordered_map>
#include <vector>
#include <set>

extern double nCr[1001][401];
extern std::vector<std::vector<daf::Size>> g_maxCliques;

using TupleKey = std::vector<daf::Size>;
struct TupleHash {
    size_t operator()(const TupleKey &t) const noexcept {
        size_t h = t.size();
        for (auto x : t) h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h << 6) + (h >> 2);
        return h;
    }
};

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex)
{
    auto tStart = std::chrono::high_resolution_clock::now();
    const daf::Size numVertices = edgeGraph.getGraphNodeSize();
    const daf::Size INVALID = std::numeric_limits<daf::Size>::max();

    // ================================================================
    // Build MCs + Classes + Tuples
    // ================================================================
    std::vector<std::vector<daf::Size>> mcs;
    for (auto &mc : g_maxCliques)
        if ((int)mc.size() >= s) mcs.push_back(mc);
    daf::Size numMC = mcs.size();

    std::vector<bool> mcAlive(numMC, true);
    std::vector<std::vector<daf::Size>> vtxAliveMCs(numVertices);
    for (daf::Size mi = 0; mi < numMC; ++mi)
        for (auto v : mcs[mi])
            if (v < numVertices) vtxAliveMCs[v].push_back(mi);

    // Classes (stable: built from initial graph, valid throughout)
    using Profile = std::vector<daf::Size>;
    struct PH { size_t operator()(const Profile &p) const noexcept {
        size_t h=p.size(); for(auto x:p) h^=std::hash<daf::Size>()(x)+0x9e3779b9ULL+(h<<6)+(h>>2); return h; }};
    std::unordered_map<Profile, daf::Size, PH> pToC;
    std::vector<daf::Size> classOf(numVertices, INVALID);
    std::vector<daf::Size> classSizes;
    std::vector<std::vector<daf::Size>> classVerts;
    for (daf::Size v = 0; v < numVertices; ++v) {
        if (vtxAliveMCs[v].empty()) continue;
        auto it = pToC.find(vtxAliveMCs[v]);
        if (it == pToC.end()) { daf::Size c=classSizes.size(); pToC[vtxAliveMCs[v]]=c;
            classSizes.push_back(0); classVerts.emplace_back(); }
        daf::Size c = pToC[vtxAliveMCs[v]];
        classOf[v]=c; classSizes[c]++; classVerts[c].push_back(v);
    }

    // MC class info
    std::vector<std::unordered_map<daf::Size,int>> mcCS(numMC);
    std::vector<std::vector<daf::Size>> mcCL(numMC);
    for (daf::Size mi=0; mi<numMC; ++mi) {
        for (auto v:mcs[mi]) if(v<numVertices && classOf[v]!=INVALID) mcCS[mi][classOf[v]]++;
        for (auto&[c,_]:mcCS[mi]) mcCL[mi].push_back(c);
        std::sort(mcCL[mi].begin(), mcCL[mi].end());
    }

    // Tuples
    struct TupleInfo { TupleKey key; daf::Size mult; };
    std::vector<TupleInfo> tuples;
    std::unordered_map<TupleKey,daf::Size,TupleHash> tidx;
    { TupleKey cur; cur.reserve(r);
      std::function<void(const std::vector<daf::Size>&,int)> en;
      en=[&](const std::vector<daf::Size>&cl,int st){
        if((int)cur.size()==r){ if(tidx.count(cur))return;
          std::unordered_map<daf::Size,int> cnt; for(auto c:cur) cnt[c]++;
          daf::Size mult=1; for(auto&[c,k]:cnt){ if((int)classSizes[c]<k)return;
            mult*=(daf::Size)(nCr[classSizes[c]][k]+0.5); }
          if(!mult)return; tidx[cur]=tuples.size(); tuples.push_back({cur,mult}); return; }
        for(int i=st;i<(int)cl.size();++i){ cur.push_back(cl[i]); en(cl,i); cur.pop_back(); }
      };
      for(daf::Size mi=0;mi<numMC;++mi){ if(mcCL[mi].size()>500)continue; cur.clear(); en(mcCL[mi],0); }
    }

    // Tuple → MCs mapping
    std::vector<std::vector<daf::Size>> tupleMCs(tuples.size());
    std::vector<std::vector<daf::Size>> mcTuples(numMC);
    for (daf::Size ti=0; ti<tuples.size(); ++ti) {
        std::unordered_map<daf::Size,int> cnt; for(auto c:tuples[ti].key) cnt[c]++;
        for (daf::Size mi=0; mi<numMC; ++mi) {
            bool ok=true; for(auto&[c,k]:cnt){ auto it=mcCS[mi].find(c);
                if(it==mcCS[mi].end()||it->second<k){ok=false;break;} }
            if(ok){ tupleMCs[ti].push_back(mi); mcTuples[mi].push_back(ti); }
        }
    }

    // ================================================================
    // IE Support computation
    // ================================================================
    auto mcInterSize = [&](const std::vector<daf::Size>&ms) -> int {
        if(ms.empty())return 0; if(ms.size()==1)return(int)mcs[ms[0]].size();
        auto cur=mcs[ms[0]]; for(size_t i=1;i<ms.size();++i){
            std::vector<daf::Size> nxt; std::set_intersection(cur.begin(),cur.end(),
                mcs[ms[i]].begin(),mcs[ms[i]].end(),std::back_inserter(nxt)); cur=std::move(nxt); }
        return(int)cur.size(); };

    auto computeSupport = [&](daf::Size ti) -> double {
        std::vector<daf::Size> alive; for(auto mi:tupleMCs[ti]) if(mcAlive[mi]) alive.push_back(mi);
        int p=(int)alive.size(); if(!p)return 0;
        double res=0; for(int mask=1;mask<(1<<p);++mask){
            std::vector<daf::Size> sub; for(int i=0;i<p;++i) if(mask&(1<<i)) sub.push_back(alive[i]);
            int isz=mcInterSize(sub); int n=isz-(int)r,k=(int)s-(int)r;
            double v=(n>=k&&k>=0)?nCr[n][k]:0.0;
            res+=(__builtin_popcount(mask)%2==1?1:-1)*v; }
        return std::max(0.0,res); };

    std::vector<double> support(tuples.size());
    for (daf::Size ti=0; ti<tuples.size(); ++ti) support[ti]=computeSupport(ti);

    // ================================================================
    // Unified peeling loop
    // ================================================================
    std::vector<bool> tuplePeeled(tuples.size(), false);
    std::vector<double> coreTuple(tuples.size(), 0);

    // Bucket queue on tuples
    constexpr int BMAX=5000000;
    double rawMax=0; for(auto&sv:support) rawMax=std::max(rawMax,sv);
    int maxB=(int)std::min(rawMax,(double)BMAX);
    int curBucket=0;
    std::vector<std::vector<daf::Size>> buckets(maxB+2);
    std::set<std::pair<double,daf::Size>> overflow;
    std::vector<int> bucket_of(tuples.size(),-1);
    std::vector<daf::Size> pos_in(tuples.size());
    std::vector<double> ovVal(tuples.size(),-1);

    auto bIns=[&](daf::Size ti){ double v=std::max(0.0,support[ti]);
        int b=std::min((int)v,maxB+1);
        if(b<=maxB){pos_in[ti]=buckets[b].size();buckets[b].push_back(ti);}
        else{overflow.insert({v,ti});ovVal[ti]=v;} bucket_of[ti]=b; };
    auto bMov=[&](daf::Size ti){ if(tuplePeeled[ti])return;
        double v=std::max(0.0,support[ti]); int ob=bucket_of[ti]; int nb=std::min((int)v,maxB+1);
        if(ob==nb&&(nb<=maxB||v==ovVal[ti]))return;
        if(ob>=0&&ob<=maxB){auto&bk=buckets[ob];daf::Size p=pos_in[ti];
            if(p<bk.size()-1){auto l=bk.back();bk[p]=l;pos_in[l]=p;}bk.pop_back();}
        else if(ob==maxB+1)overflow.erase({ovVal[ti],ti});
        if(nb<=maxB){pos_in[ti]=buckets[nb].size();buckets[nb].push_back(ti);if(nb<curBucket)curBucket=nb;}
        else{overflow.insert({v,ti});ovVal[ti]=v;} bucket_of[ti]=nb; };

    for(daf::Size ti=0;ti<tuples.size();++ti) bIns(ti);
    daf::Size remaining=tuples.size(); double minCore=0;

    // MC priority queue: (C value, mc index)
    using PQE=std::pair<double,daf::Size>;
    std::priority_queue<PQE,std::vector<PQE>,std::greater<>> mcPQ;
    auto mcVal=[&](daf::Size mi)->double{ int n=(int)mcs[mi].size()-(int)r,k=(int)s-(int)r;
        return(n>=k&&k>=0)?nCr[n][k]:0.0; };
    auto hasPriv=[&](daf::Size mi)->bool{ if(!mcAlive[mi])return false;
        for(auto v:mcs[mi]) if(v<numVertices&&vtxAliveMCs[v].size()==1&&vtxAliveMCs[v][0]==mi) return true;
        return false; };
    for(daf::Size mi=0;mi<numMC;++mi) if(hasPriv(mi)) mcPQ.push({mcVal(mi),mi});

    // Kill MC + cascade + recompute
    auto killMC=[&](daf::Size mi){
        if(!mcAlive[mi])return;
        mcAlive[mi]=false;
        // Recompute support for all tuples in this MC
        for(auto ti:mcTuples[mi]){
            if(tuplePeeled[ti])continue;
            support[ti]=computeSupport(ti);
            bMov(ti);
        }
        // Update vertex profiles, find new private vertices
        std::set<daf::Size> check;
        for(auto v:mcs[mi]){ if(v>=numVertices)continue;
            auto&am=vtxAliveMCs[v]; am.erase(std::remove(am.begin(),am.end(),mi),am.end());
            if(am.size()==1&&mcAlive[am[0]]) check.insert(am[0]); }
        for(auto m2:check) if(hasPriv(m2)) mcPQ.push({mcVal(m2),m2});
    };

    auto tPeel = std::chrono::high_resolution_clock::now();

    while(remaining>0){
        // Kill MCs whose C value ≤ current level
        while(!mcPQ.empty()){
            auto[v,mi]=mcPQ.top();
            if(!mcAlive[mi]||!hasPriv(mi)){mcPQ.pop();continue;}
            if(v>minCore+0.5) break; // MC too strong for current level
            mcPQ.pop(); killMC(mi);
        }

        // Pop tuple
        while(curBucket<(int)buckets.size()&&buckets[curBucket].empty()) curBucket++;
        if(curBucket>=(int)buckets.size()){
            if(overflow.empty())break;
            auto it=overflow.begin(); daf::Size ti=it->second; overflow.erase(it);
            if(tuplePeeled[ti])continue;
            minCore=std::max(support[ti],minCore);
            tuplePeeled[ti]=true; coreTuple[ti]=minCore; remaining--;
            // Check if any MC should die now
            while(!mcPQ.empty()){ auto[v,mi]=mcPQ.top();
                if(!mcAlive[mi]||!hasPriv(mi)){mcPQ.pop();continue;}
                if(v>minCore+0.5)break; mcPQ.pop(); killMC(mi); }
            continue;
        }

        minCore=std::max((double)curBucket,minCore);

        // Kill MCs at current level
        while(!mcPQ.empty()){ auto[v,mi]=mcPQ.top();
            if(!mcAlive[mi]||!hasPriv(mi)){mcPQ.pop();continue;}
            if(v>minCore+0.5)break; mcPQ.pop(); killMC(mi); }

        // Pop all tuples at current bucket
        while(curBucket<(int)buckets.size()&&!buckets[curBucket].empty()&&curBucket<=(int)minCore){
            while(!buckets[curBucket].empty()){
                daf::Size ti=buckets[curBucket].back(); buckets[curBucket].pop_back();
                if(tuplePeeled[ti])continue;
                tuplePeeled[ti]=true; coreTuple[ti]=minCore; remaining--;
            }
            if(curBucket+1<(int)buckets.size()&&!buckets[curBucket+1].empty()&&(curBucket+1)<=(int)minCore) curBucket++;
            else break;
        }

        // After popping: check MC deaths again
        while(!mcPQ.empty()){ auto[v,mi]=mcPQ.top();
            if(!mcAlive[mi]||!hasPriv(mi)){mcPQ.pop();continue;}
            if(v>minCore+0.5)break; mcPQ.pop(); killMC(mi); }
    }

    auto tEnd = std::chrono::high_resolution_clock::now();

    // ================================================================
    // Output
    // ================================================================
    std::map<double,daf::Size> coreDist;
    for(daf::Size ti=0;ti<tuples.size();++ti) coreDist[coreTuple[ti]]+=tuples[ti].mult;

    printf("======= V4 Unified =======\n");
    printf("  r=%d s=%d, MCs: %zu, Tuples: %zu\n",(int)r,(int)s,(size_t)numMC,tuples.size());
    printf("  Total time: %lld ms\n",
        std::chrono::duration_cast<std::chrono::milliseconds>(tEnd-tStart).count());
    printf("  Peel time: %lld ms\n",
        std::chrono::duration_cast<std::chrono::milliseconds>(tEnd-tPeel).count());
    double maxCore=0; for(auto&[c,_]:coreDist) maxCore=std::max(maxCore,c);
    printf("  Max core: %.0f\n",maxCore);
    for(auto&[c,cnt]:coreDist) printf("  core=%.0f count=%zu\n",c,(size_t)cnt);

    std::vector<std::pair<std::vector<daf::Size>,double>> result;
    for(auto&[c,cnt]:coreDist) result.push_back({{},c});
    return result;
}
