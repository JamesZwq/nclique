// NewFastLocal V9: Pre-allocated memory pool for leafCls
// Eliminates malloc contention in parallel precompute
// Use nClique as upper bound for total storage
#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include <algorithm>
#include <numeric>
#include <cstring>
#include <chrono>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
extern double nCr[1001][401];
namespace NewFastLocal {
std::vector<std::pair<std::vector<daf::Size>,int>>
NucleusCoreDecompositionNewFastLocal(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s)
{
    const bool verbose=std::getenv("PIVOTER_VERBOSE")!=nullptr;
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build",[&](){cliqueIndex.build(tree,edgeGraph.adj_list.size());});
    const daf::Size nClique=cliqueIndex.size(), nLeaves=tree.adj_list.size();
#ifdef _OPENMP
    const int nt=omp_get_max_threads();
#else
    const int nt=1;
#endif
    auto t0=std::chrono::high_resolution_clock::now();

    // Pass 1 (cheap): count r-cliques per leaf to build CSR offsets
    // Only enumeration, no storage
    std::vector<double> sup(nClique,0.0);
    std::vector<std::vector<double>> tl(nt,std::vector<double>(nClique,0.0));
    std::vector<uint32_t> lsz(nLeaves,0);

#ifdef _OPENMP
    #pragma omp parallel num_threads(nt)
    { int tid=omp_get_thread_num(); auto& loc=tl[tid];
      #pragma omp for schedule(dynamic,32) nowait
      for(daf::Size li=0;li<nLeaves;++li){
        const auto& leaf=tree.adj_list[li]; if((int)leaf.size()<r) continue;
        daf::CliqueSize pC=0,kC=0;
        for(const auto& n:leaf){pC+=n.isPivot;kC+=!n.isPivot;}
        int nP=s-(int)kC; if(nP<0||nP>(int)pC) continue;
        uint32_t cnt=0;
        daf::enumerateCombinations(leaf,r,[&](const daf::StaticVector<TreeGraphNode>& rc){
          daf::CliqueSize sp=0; for(const auto& n:rc) sp+=n.isPivot;
          int p1=pC-sp,p2=nP-sp; daf::Size cid=cliqueIndex.byClique(rc);
          ++cnt;
          if(p1>=0&&p1<1001&&p2>=0&&p2<401) loc[cid]+=nCr[p1][p2];
          return true;
        }); lsz[li]=cnt;
      }
    }
    #pragma omp parallel for schedule(static,1024) num_threads(nt)
    for(daf::Size i=0;i<nClique;++i) for(int t=0;t<nt;++t) sup[i]+=tl[t][i];
#else
    for(daf::Size li=0;li<nLeaves;++li){
      const auto& leaf=tree.adj_list[li]; if((int)leaf.size()<r) continue;
      daf::CliqueSize pC=0,kC=0;
      for(const auto& n:leaf){pC+=n.isPivot;kC+=!n.isPivot;}
      int nP=s-(int)kC; if(nP<0||nP>(int)pC) continue;
      uint32_t cnt=0;
      daf::enumerateCombinations(leaf,r,[&](const daf::StaticVector<TreeGraphNode>& rc){
        daf::CliqueSize sp=0; for(const auto& n:rc) sp+=n.isPivot;
        int p1=pC-sp,p2=nP-sp; daf::Size cid=cliqueIndex.byClique(rc);
        ++cnt;
        if(p1>=0&&p1<1001&&p2>=0&&p2<401) sup[cid]+=nCr[p1][p2];
        return true;
      }); lsz[li]=cnt;
    }
#endif

    // Build CSR offsets
    std::vector<uint32_t> loff(nLeaves+1,0);
    for(daf::Size li=0;li<nLeaves;++li) loff[li+1]=loff[li]+lsz[li];
    uint32_t tot=loff[nLeaves];
    // Pre-allocate flat array (single malloc)
    std::vector<uint32_t> flat(tot);

    // Pass 2: fill flat array in parallel (no malloc, direct writes to pre-allocated)
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,32) num_threads(nt)
#endif
    for(daf::Size li=0;li<nLeaves;++li){
      if(!lsz[li]) continue;
      const auto& leaf=tree.adj_list[li];
      uint32_t pos=loff[li];
      daf::enumerateCombinations(leaf,r,[&](const daf::StaticVector<TreeGraphNode>& rc){
        flat[pos++]=(uint32_t)cliqueIndex.byClique(rc); return true;
      });
    }

    // Active leaves
    std::vector<daf::Size> act; act.reserve(nLeaves);
    for(daf::Size li=0;li<nLeaves;++li) if(lsz[li]>=2) act.push_back(li);

    if(verbose){
      auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now()-t0).count();
      std::cout<<"[Local] Precompute: "<<ms<<"ms act="<<act.size()<<" tot="<<tot<<std::endl;
    }

    // Peeling - pure array ops, zero enumeration
    auto tp=std::chrono::high_resolution_clock::now();
    std::vector<double> core(sup),newCore(nClique);
    const int pnt=std::min(nt,16);
    std::vector<std::vector<double>> tlC(pnt,std::vector<double>(nClique));
    std::size_t iters=0;
    while(iters<100){
      iters++;
#ifdef _OPENMP
      #pragma omp parallel num_threads(pnt)
      { int tid=omp_get_thread_num(); auto& my=tlC[tid];
        std::fill(my.begin(),my.end(),1e18);
        #pragma omp for schedule(dynamic,128) nowait
        for(daf::Size ai=0;ai<act.size();++ai){
          daf::Size li=act[ai];
          const uint32_t*cls=flat.data()+loff[li]; int n=(int)lsz[li];
          // Sort once O(n log n), then h-index O(n), without-self O(1) per element
          std::vector<double> sc(n);
          for(int i=0;i<n;++i) sc[i]=core[cls[i]];
          std::sort(sc.rbegin(),sc.rend());
          // Compute hFull
          int hF=0; for(int k=0;k<n;++k){if(sc[k]>=(double)(k+1))hF=k+1;else break;}
          if(!hF) continue;
          // Count elements with core >= hF (prefix of sorted array)
          int cntHF=0; for(int k=0;k<n;++k){if(sc[k]>=(double)hF)cntHF++;else break;}
          for(int i=0;i<n;++i){
            uint32_t ci=cls[i]; int h;
            if(core[ci]<(double)hF){
              h=hF; // not in top group, removing ci doesn't reduce hFull
            } else {
              // ci is in top group; after removal: (cntHF-1) >= hF ?
              h=((cntHF-1)>=hF)?hF:hF-1;
            }
            if((double)h<my[ci]) my[ci]=(double)h;
          }
        }
      }
#else
      std::fill(tlC[0].begin(),tlC[0].end(),1e18);
      for(daf::Size ai=0;ai<act.size();++ai){
        daf::Size li=act[ai];
        const uint32_t*cls=flat.data()+loff[li]; int n=(int)lsz[li];
        std::vector<double> sc(n);
        for(int i=0;i<n;++i) sc[i]=core[cls[i]];
        std::sort(sc.rbegin(),sc.rend());
        int hF=0; for(int k=0;k<n;++k){if(sc[k]>=(double)(k+1))hF=k+1;else break;}
        if(!hF) continue;
        int cntHF=0; for(int k=0;k<n;++k){if(sc[k]>=(double)hF)cntHF++;else break;}
        for(int i=0;i<n;++i){
          uint32_t ci=cls[i]; int h;
          if(core[ci]<(double)hF) h=hF;
          else h=((cntHF-1)>=hF)?hF:hF-1;
          if((double)h<tlC[0][ci]) tlC[0][ci]=(double)h;
        }
      }
#endif
      // Min-reduction
#ifdef _OPENMP
      #pragma omp parallel for schedule(static,1024) num_threads(pnt)
#endif
      for(daf::Size i=0;i<nClique;++i){
        double v=sup[i]; for(int t=0;t<pnt;++t) if(tlC[t][i]<v) v=tlC[t][i];
        newCore[i]=v;
      }
      bool conv=true;
#ifdef _OPENMP
      #pragma omp parallel for reduction(&&:conv) num_threads(pnt)
#endif
      for(daf::Size i=0;i<nClique;++i) if(newCore[i]!=core[i]) conv=false;
      std::swap(core,newCore);
      if(conv) break;
    }
    if(verbose){
      auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now()-tp).count();
      std::cout<<"[Local] Peeling: "<<ms<<"ms ("<<iters<<" iters)"<<std::endl;
    }
    std::vector<std::pair<std::vector<daf::Size>,int>> result;
    for(daf::Size i=0;i<nClique;++i)
      result.push_back({std::vector<daf::Size>(cliqueIndex.byId(i).begin(),cliqueIndex.byId(i).end()),(int)core[i]});
    return result;
}
} // namespace NewFastLocal
