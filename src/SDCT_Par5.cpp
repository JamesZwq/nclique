/*
 * SDCT_Par5: Prefetch + flat leaf arena
 * New over Par4:
 * 1. __builtin_prefetch on adj_list traversal
 * 2. Flat leaf storage: contiguous LeafArena avoids per-leaf heap alloc
 */

#include <sys/mman.h>
#include <omp.h>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <iostream>

#include "degeneracy_algorithm_cliques_V.h"
#include "degeneracy_helper.h"
#include "LinkedList.h"
#include "MemoryManager.h"
#include "misc.h"
#include "tree/MultiBranchTree.h"
#include "graph/DynamicGraph.h"
#include "graph/Graph.h"
#include "dataStruct/CliqueCSR.hpp"

extern double nCr[1001][401];

struct Arena5 {
    int* base=nullptr; int top=0; int cap=0; bool owned=true;
    void init(int n){
        cap=n+16;
        base=(int*)std::malloc(cap*sizeof(int));
        top=0; owned=true;
    }
    ~Arena5(){if(owned&&base)std::free(base);}
    int* alloc_raw(int n){if(n<=0)n=1;int*p=base+top;top+=n;return p;}
    int save(){return top;}
    void restore(int t){top=t;}
};
static thread_local Arena5 g_arena5;

struct MarkArray5 {
    std::vector<int> mark; int gen=0;
    void init(int n){mark.assign(n,0);}
    void next(){++gen;}
    void set(int v){mark[v]=gen;}
    bool get(int v)const{return mark[v]==gen;}
};
static thread_local MarkArray5 g_mark5;

// Flat leaf arena: [total, v0, drop0, v1, drop1, ...] per leaf
// Stores data as flat int array (no per-leaf malloc during computation)
// flush() reconstructs vectors only once at the end
struct LeafArena5 {
    std::vector<int> buf;
    std::vector<int> offsets;
    void reserve(int n){buf.reserve(n*10);offsets.reserve(n);}
    void add(int* keepV,int keepSz,int* dropV,int dropSz){
        offsets.push_back((int)buf.size());
        int total=keepSz+dropSz;
        buf.push_back(total);
        TreeGraphNode tmp[512];
        for(int i=0;i<keepSz;i++)tmp[i]={(daf::Size)keepV[i],false};
        for(int i=0;i<dropSz;i++)tmp[keepSz+i]={(daf::Size)dropV[i],true};
        if(total<=32){
            for(int i=1;i<total;i++){TreeGraphNode k=tmp[i];int j=i-1;
                while(j>=0&&tmp[j]>k){tmp[j+1]=tmp[j];j--;}tmp[j+1]=k;}
        }else{std::sort(tmp,tmp+total);}
        for(int i=0;i<total;i++){buf.push_back((int)tmp[i].v);buf.push_back((int)tmp[i].isPivot);}
    }
    void flush(DynamicGraph<TreeGraphNode>& out){
        for(int oi=0;oi<(int)offsets.size();oi++){
            int pos=offsets[oi],sz=buf[pos++];
            std::vector<TreeGraphNode> leaf; leaf.reserve(sz);
            for(int i=0;i<sz;i++){leaf.emplace_back((daf::Size)buf[pos],(bool)buf[pos+1]);pos+=2;}
            out.adj_list.push_back(std::move(leaf));
        }
    }
    void flush_range(DynamicGraph<TreeGraphNode>& out, size_t base_idx){
        for(int oi=0;oi<(int)offsets.size();oi++){
            int pos=offsets[oi],sz=buf[pos++];
            std::vector<TreeGraphNode> leaf; leaf.reserve(sz);
            for(int i=0;i<sz;i++){leaf.emplace_back((daf::Size)buf[pos],(bool)buf[pos+1]);pos+=2;}
            out.adj_list[base_idx+oi]=std::move(leaf);
        }
    }
    void flush_csr(CliqueCSR<int>& out){
        for(int oi=0;oi<(int)offsets.size();oi++){
            int pos=offsets[oi],sz=buf[pos++];
            // Extract vertex IDs directly (skip isPivot flags)
            std::vector<int> vertices;
            vertices.reserve(sz);
            for(int i=0;i<sz;i++){
                vertices.push_back((int)buf[pos]);  // only store vertex id, skip isPivot
                pos+=2;
            }
            out.add_clique(vertices);
        }
    }
    size_t size() const { return offsets.size(); }
};
static thread_local LeafArena5 g_leafarena5;

static int findPivotMark5(
    int** pivotNonNeighbors,int* numNonNeighbors,
    int* vertexSets,int* vertexLookup,
    int** neighborsInP,int* numNeighbors,
    int beginP,int beginR)
{
    int pivot=-1,maxInP=-1,sz=beginR-beginP;
    for(int j=beginP;j<beginR;j++){
        int v=vertexSets[j],numP=MY_MIN(sz,numNeighbors[v]),inP=0;
        for(int k=0;k<numP;k++){int l=vertexLookup[neighborsInP[v][k]];if(l>=beginP&&l<beginR)inP++;else break;}
        if(inP>maxInP){maxInP=inP;pivot=v;}
    }
    g_mark5.next();
    int nPN=MY_MIN(sz,numNeighbors[pivot]);
    for(int j=0;j<nPN;j++){int nb=neighborsInP[pivot][j],l=vertexLookup[nb];if(l>=beginP&&l<beginR)g_mark5.set(nb);else break;}
    int*buf=g_arena5.alloc_raw(sz),nc=0;
    for(int j=beginP;j<beginR;j++){int v=vertexSets[j];if(!g_mark5.get(v))buf[nc++]=v;}
    *pivotNonNeighbors=buf;*numNonNeighbors=nc;
    return pivot;
}

static void moveToRFast5(
    int vertex,int* vertexSets,int* vertexLookup,
    int** neighborsInP,int* numNeighbors,
    int* pBeginP,int* pBeginR,int* pNewBeginP,int* pNewBeginR)
{
    int vl=vertexLookup[vertex];
    (*pBeginR)--;vertexSets[vl]=vertexSets[*pBeginR];vertexLookup[vertexSets[*pBeginR]]=vl;
    vertexSets[*pBeginR]=vertex;vertexLookup[vertex]=*pBeginR;
    *pNewBeginP=*pBeginP;*pNewBeginR=*pBeginP;
    int sizeOfP=*pBeginR-*pBeginP;
    g_mark5.next();
    int numNbr=MY_MIN(sizeOfP,numNeighbors[vertex]);
    for(int k=0;k<numNbr;k++)g_mark5.set(neighborsInP[vertex][k]);
    for(int j=*pBeginP;j<*pBeginR;j++){
        int nb=vertexSets[j];
        if(g_mark5.get(nb)){
            vertexSets[j]=vertexSets[*pNewBeginR];vertexLookup[vertexSets[*pNewBeginR]]=j;
            vertexSets[*pNewBeginR]=nb;vertexLookup[nb]=*pNewBeginR;(*pNewBeginR)++;
        }
    }
    for(int j2=*pNewBeginP;j2<*pNewBeginR;j2++){
        int tv=vertexSets[j2],numP=MY_MIN(sizeOfP,numNeighbors[tv]),cnt=0;
        for(int k=0;k<numP;k++){int nbr=neighborsInP[tv][k],nl=vertexLookup[nbr];
            if(nl>=*pNewBeginP&&nl<*pNewBeginR){neighborsInP[tv][k]=neighborsInP[tv][cnt];neighborsInP[tv][cnt++]=nbr;}}
    }
}

static void fillInPandXArena5(
    int vertex,int* vertexSets,int* vertexLookup,
    Graph& edgeGraph,int** neighborsInP,int* numNeighbors,
    int* pBeginX,int* pBeginP,int* pBeginR,
    int* pNewBeginX,int* pNewBeginP,int* pNewBeginR)
{
    int vl=vertexLookup[vertex];
    (*pBeginR)--;vertexSets[vl]=vertexSets[*pBeginR];vertexLookup[vertexSets[*pBeginR]]=vl;
    vertexSets[*pBeginR]=vertex;vertexLookup[vertex]=*pBeginR;
    *pNewBeginR=*pBeginR;*pNewBeginP=*pBeginR;
    auto [nb,ne]=edgeGraph.getNbr(vertex);
    if(nb<ne)__builtin_prefetch(&edgeGraph.adj_list[nb],0,1);
    for(int idx=nb;idx<ne;++idx){
        if(idx+4<ne)__builtin_prefetch(&edgeGraph.adj_list[idx+4],0,0);
        int neighbor=edgeGraph.adj_list[idx];
        if(neighbor<=vertex)continue;
        int nl=vertexLookup[neighbor];(*pNewBeginP)--;
        vertexSets[nl]=vertexSets[*pNewBeginP];vertexLookup[vertexSets[*pNewBeginP]]=nl;
        vertexSets[*pNewBeginP]=neighbor;vertexLookup[neighbor]=*pNewBeginP;
    }
    *pNewBeginX=*pNewBeginP;
    int pSize=*pNewBeginR-*pNewBeginP;
    for(int j=*pNewBeginP;j<*pNewBeginR;++j){
        int v=vertexSets[j];numNeighbors[v]=0;
        int as=MY_MIN(pSize,edgeGraph.getNbrCount(v));
        neighborsInP[v]=g_arena5.alloc_raw(as>0?as:1);
    }
    for(int j=*pNewBeginP;j<*pNewBeginR;++j){
        int v=vertexSets[j];
        auto [nb2,ne2]=edgeGraph.getNbr(v);
        if(nb2<ne2)__builtin_prefetch(&edgeGraph.adj_list[nb2],0,1);
        for(auto idx=nb2;idx<ne2;++idx){
            if(idx+4<ne2)__builtin_prefetch(&edgeGraph.adj_list[idx+4],0,0);
            int w=edgeGraph.adj_list[idx];
            if(w<=v)continue;
            int wl=vertexLookup[w];
            if(wl>=*pNewBeginP&&wl<*pNewBeginR){
                neighborsInP[v][numNeighbors[v]++]=w;
                neighborsInP[w][numNeighbors[w]++]=v;
            }
        }
    }
}

static void recurse5(
    int* vertexSets,int* vertexLookup,
    int** neighborsInP,int* numNeighbors,
    int beginP,int beginR,
    int* keepV,int keepSz,int* dropV,int dropSz,
    int max_k,int min_k)
{
    if(keepSz>max_k)return;
    if(beginP>=beginR){
        int cSize=keepSz+dropSz;
        if(cSize<min_k)return;
        if(keepSz==max_k) g_leafarena5.add(keepV,keepSz,dropV,0);
        else              g_leafarena5.add(keepV,keepSz,dropV,dropSz);
        return;
    }
    int*cands;int nCands;
    int arenaTop=g_arena5.save();
    int pivot=findPivotMark5(&cands,&nCands,vertexSets,vertexLookup,
                             neighborsInP,numNeighbors,beginP,beginR);
    for(int ci=0;ci<nCands;ci++){
        int vertex=cands[ci],newBeginP,newBeginR;
        moveToRFast5(vertex,vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     &beginP,&beginR,&newBeginP,&newBeginR);
        if(vertex==pivot){dropV[dropSz]=vertex;
            recurse5(vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     newBeginP,newBeginR,keepV,keepSz,dropV,dropSz+1,max_k,min_k);
        }else{keepV[keepSz]=vertex;
            recurse5(vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     newBeginP,newBeginR,keepV,keepSz+1,dropV,dropSz,max_k,min_k);
        }
    }
    g_arena5.restore(arenaTop);
}

DynamicGraph<TreeGraphNode> SDCT_Par5(Graph& edgeGraph,int max_k,int min_k){
    auto size=(int)edgeGraph.getGraphNodeSize();
    int nthreads=omp_get_max_threads();
    std::cout<<"SDCT_Par5 "<<nthreads<<" threads"<<std::endl;
    std::vector<LeafArena5> thread_leaves(nthreads);

    // Pre-allocate all per-thread buffers SERIALLY before parallel region
    // This eliminates malloc contention (32 threads competing for glibc malloc lock)
    int arena_cap = size * 8;  // Balanced: enough for large neighborhoods
    struct ThreadBufs {
        int* vertexSets=nullptr; int* vertexLookup=nullptr;
        int* numNeighbors=nullptr; int** neighborsInP=nullptr;
        int* arena=nullptr; int* mark=nullptr;
    };
    std::vector<ThreadBufs> bufs(nthreads);
    double t_alloc0 = omp_get_wtime();
    for(int t=0;t<nthreads;t++){
        bufs[t].vertexSets   = (int*)std::malloc(size*sizeof(int));
        bufs[t].vertexLookup = (int*)std::malloc(size*sizeof(int));
        bufs[t].numNeighbors = (int*)std::malloc(size*sizeof(int));
        bufs[t].neighborsInP = (int**)std::malloc(size*sizeof(int*));
        bufs[t].arena        = (int*)std::malloc(arena_cap*sizeof(int));
        bufs[t].mark         = (int*)std::calloc(size, sizeof(int));
    }
    double t_alloc1 = omp_get_wtime();
    printf("Serial pre-alloc took: %.1f ms\n", (t_alloc1-t_alloc0)*1000);

    #pragma omp parallel
    {
        int tid=omp_get_thread_num();

        // Point thread-local structures at pre-allocated buffers
        g_arena5.base = bufs[tid].arena;
        g_arena5.cap  = arena_cap;
        g_arena5.top  = size;
        g_arena5.owned = false;
        g_mark5.mark.assign(bufs[tid].mark, bufs[tid].mark + size);
        g_mark5.gen   = 0;
        g_leafarena5.buf.clear(); g_leafarena5.offsets.clear();
        g_leafarena5.reserve(std::max(1,size/nthreads)*20);

        int* vertexSets    = bufs[tid].vertexSets;
        int* vertexLookup  = bufs[tid].vertexLookup;
        int* numNeighbors  = bufs[tid].numNeighbors;
        int** neighborsInP = bufs[tid].neighborsInP;

        // First-touch initialization: each thread initializes ONLY its own buffers
        // This achieves NUMA locality without redundant work
        for(int i=0;i<size;i++){
            vertexSets[i]=i; 
            vertexLookup[i]=i; 
            numNeighbors[i]=1; 
            neighborsInP[i]=bufs[tid].arena+i;
        }

        // Barrier to synchronize all threads after initialization
        #pragma omp barrier

        int beginX=0,beginP=0,beginR=size;
        int keepV[MAX_CSIZE],dropV[MAX_CSIZE];

        double t_work=0;
        double t_par_start = omp_get_wtime();

        #pragma omp for schedule(guided) nowait
        for(int vertex=0;vertex<size;vertex++){
            double tw0=omp_get_wtime();
            // beginR is deterministic: after processing vertices 0..vertex-1,
            // each one increments beginR by 1, so beginR = size - vertex
            int localBeginR = size - vertex;
            int localBeginP = 0;
            int localBeginX = 0;
            (void)beginX; (void)beginP; (void)beginR;
            int arenaBase=g_arena5.save();
            int newBeginX,newBeginP,newBeginR;
            fillInPandXArena5(vertex,vertexSets,vertexLookup,edgeGraph,
                              neighborsInP,numNeighbors,
                              &localBeginX,&localBeginP,&localBeginR,&newBeginX,&newBeginP,&newBeginR);
            keepV[0]=vertex;
            recurse5(vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     newBeginP,newBeginR,keepV,1,dropV,0,max_k,min_k);
            g_arena5.restore(arenaBase);
            t_work += omp_get_wtime()-tw0;
        }
        double t_par_end = omp_get_wtime();
        #pragma omp critical
        printf("Thread %2d: work=%.1fms total=%.1fms efficiency=%.0f%%\n",
               tid, t_work*1000, (t_par_end-t_par_start)*1000,
               100.0*t_work/(t_par_end-t_par_start));
        #pragma omp critical
        thread_leaves[tid]=std::move(g_leafarena5);
        // No free here - bufs are freed after parallel region
    }

    // Free all pre-allocated buffers
    for(int t=0;t<nthreads;t++){
        std::free(bufs[t].vertexSets); std::free(bufs[t].vertexLookup);
        std::free(bufs[t].numNeighbors); std::free(bufs[t].neighborsInP);
        std::free(bufs[t].arena); std::free(bufs[t].mark);
    }

    size_t total=0;
    double t_merge0 = omp_get_wtime();
    for(auto&tl:thread_leaves)total+=tl.size();
    
    DynamicGraph<TreeGraphNode> treeGraph(size);
    treeGraph.adj_list.reserve(total);
    
    // Serial flush with pre-reserved space (no reallocation)
    for(int t=0;t<nthreads;t++) thread_leaves[t].flush(treeGraph);
    
    double t_merge1 = omp_get_wtime();
    printf("Result merge took: %.1f ms (total cliques: %zu)\n", (t_merge1-t_merge0)*1000, total);
    return treeGraph;
}

// CSR version: much faster merge using flat data structure
CliqueCSR<int> SDCT_Par5_CSR(Graph& edgeGraph,int max_k,int min_k){
    auto size=(int)edgeGraph.getGraphNodeSize();
    int nthreads=omp_get_max_threads();
    std::cout<<"SDCT_Par5_CSR "<<nthreads<<" threads"<<std::endl;
    std::vector<LeafArena5> thread_leaves(nthreads);

    // Pre-allocate all per-thread buffers SERIALLY before parallel region
    int arena_cap = size * 8;
    struct ThreadBufs {
        int* vertexSets=nullptr; int* vertexLookup=nullptr;
        int* numNeighbors=nullptr; int** neighborsInP=nullptr;
        int* arena=nullptr; int* mark=nullptr;
    };
    std::vector<ThreadBufs> bufs(nthreads);
    double t_alloc0 = omp_get_wtime();
    for(int t=0;t<nthreads;t++){
        bufs[t].vertexSets   = (int*)std::malloc(size*sizeof(int));
        bufs[t].vertexLookup = (int*)std::malloc(size*sizeof(int));
        bufs[t].numNeighbors = (int*)std::malloc(size*sizeof(int));
        bufs[t].neighborsInP = (int**)std::malloc(size*sizeof(int*));
        bufs[t].arena        = (int*)std::malloc(arena_cap*sizeof(int));
        bufs[t].mark         = (int*)std::calloc(size, sizeof(int));
    }
    double t_alloc1 = omp_get_wtime();
    printf("Serial pre-alloc took: %.1f ms\n", (t_alloc1-t_alloc0)*1000);

    #pragma omp parallel
    {
        int tid=omp_get_thread_num();

        g_arena5.base = bufs[tid].arena;
        g_arena5.cap  = arena_cap;
        g_arena5.top  = size;
        g_arena5.owned = false;
        g_mark5.mark.assign(bufs[tid].mark, bufs[tid].mark + size);
        g_mark5.gen   = 0;
        g_leafarena5.buf.clear(); g_leafarena5.offsets.clear();
        g_leafarena5.reserve(std::max(1,size/nthreads)*20);

        int* vertexSets    = bufs[tid].vertexSets;
        int* vertexLookup  = bufs[tid].vertexLookup;
        int* numNeighbors  = bufs[tid].numNeighbors;
        int** neighborsInP = bufs[tid].neighborsInP;

        // First-touch initialization: each thread initializes ONLY its own buffers
        for(int i=0;i<size;i++){
            vertexSets[i]=i; 
            vertexLookup[i]=i; 
            numNeighbors[i]=1; 
            neighborsInP[i]=bufs[tid].arena+i;
        }

        #pragma omp barrier

        int beginX=0,beginP=0,beginR=size;
        int keepV[MAX_CSIZE],dropV[MAX_CSIZE];

        double t_work=0;
        double t_par_start = omp_get_wtime();

        #pragma omp for schedule(guided) nowait
        for(int vertex=0;vertex<size;vertex++){
            double tw0=omp_get_wtime();
            int localBeginR = size - vertex;
            int localBeginP = 0;
            int localBeginX = 0;
            (void)beginX; (void)beginP; (void)beginR;
            int arenaBase=g_arena5.save();
            int newBeginX,newBeginP,newBeginR;
            fillInPandXArena5(vertex,vertexSets,vertexLookup,edgeGraph,
                              neighborsInP,numNeighbors,
                              &localBeginX,&localBeginP,&localBeginR,&newBeginX,&newBeginP,&newBeginR);
            keepV[0]=vertex;
            recurse5(vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     newBeginP,newBeginR,keepV,1,dropV,0,max_k,min_k);
            g_arena5.restore(arenaBase);
            t_work += omp_get_wtime()-tw0;
        }
        double t_par_end = omp_get_wtime();
        #pragma omp critical
        printf("Thread %2d: work=%.1fms total=%.1fms efficiency=%.0f%%\n",
               tid, t_work*1000, (t_par_end-t_par_start)*1000,
               100.0*t_work/(t_par_end-t_par_start));
        #pragma omp critical
        thread_leaves[tid]=std::move(g_leafarena5);
    }

    // Free all pre-allocated buffers
    for(int t=0;t<nthreads;t++){
        std::free(bufs[t].vertexSets); std::free(bufs[t].vertexLookup);
        std::free(bufs[t].numNeighbors); std::free(bufs[t].neighborsInP);
        std::free(bufs[t].arena); std::free(bufs[t].mark);
    }

    size_t total=0;
    double t_merge0 = omp_get_wtime();
    for(auto&tl:thread_leaves)total+=tl.size();
    
    CliqueCSR<int> result;
    result.reserve_cliques(total);
    
    // Serial flush with CSR format (much faster!)
    for(int t=0;t<nthreads;t++) thread_leaves[t].flush_csr(result);
    
    double t_merge1 = omp_get_wtime();
    printf("Result merge took: %.1f ms (total cliques: %zu)\n", (t_merge1-t_merge0)*1000, total);
    return result;
}
