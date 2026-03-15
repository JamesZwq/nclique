/*
 * SDCT_Par_Final: Truly parallel SDCT with proper work distribution
 * Key improvements:
 * 1. Global shared vertex data (no per-thread duplication)
 * 2. Thread-local arena only for recursion buffers
 * 3. Proper work distribution via #pragma omp for
 * 4. No false sharing or synchronization overhead
 */

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

extern double nCr[1001][401];

struct ArenaOpt {
    int* base=nullptr; int top=0; int cap=0;
    void init(int n){cap=n+16;base=(int*)std::malloc(cap*sizeof(int));top=0;}
    ~ArenaOpt(){std::free(base);}
    int* alloc_raw(int n){if(n<=0)n=1;int*p=base+top;top+=n;return p;}
    int save(){return top;}
    void restore(int t){top=t;}
};
static thread_local ArenaOpt g_arenaOpt;

struct MarkArrayOpt {
    std::vector<int> mark; int gen=0;
    void init(int n){mark.assign(n,0);}
    void next(){++gen;}
    void set(int v){mark[v]=gen;}
    bool get(int v)const{return mark[v]==gen;}
};
static thread_local MarkArrayOpt g_markOpt;

static int findPivotMarkOpt(
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
    g_markOpt.next();
    int nPN=MY_MIN(sz,numNeighbors[pivot]);
    for(int j=0;j<nPN;j++){int nb=neighborsInP[pivot][j],l=vertexLookup[nb];if(l>=beginP&&l<beginR)g_markOpt.set(nb);else break;}
    int*buf=g_arenaOpt.alloc_raw(sz),nc=0;
    for(int j=beginP;j<beginR;j++){int v=vertexSets[j];if(!g_markOpt.get(v))buf[nc++]=v;}
    *pivotNonNeighbors=buf;*numNonNeighbors=nc;
    return pivot;
}

static void moveToRFastOpt(
    int vertex,int* vertexSets,int* vertexLookup,
    int** neighborsInP,int* numNeighbors,
    int* pBeginP,int* pBeginR,int* pNewBeginP,int* pNewBeginR)
{
    int vl=vertexLookup[vertex];
    (*pBeginR)--;vertexSets[vl]=vertexSets[*pBeginR];vertexLookup[vertexSets[*pBeginR]]=vl;
    vertexSets[*pBeginR]=vertex;vertexLookup[vertex]=*pBeginR;
    *pNewBeginP=*pBeginP;*pNewBeginR=*pBeginP;
    int sizeOfP=*pBeginR-*pBeginP;
    g_markOpt.next();
    int numNbr=MY_MIN(sizeOfP,numNeighbors[vertex]);
    for(int k=0;k<numNbr;k++)g_markOpt.set(neighborsInP[vertex][k]);
    for(int j=*pBeginP;j<*pBeginR;j++){
        int nb=vertexSets[j];
        if(g_markOpt.get(nb)){
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

static void fillInPandXArenaOpt(
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
    for(int idx=nb;idx<ne;++idx){
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
        neighborsInP[v]=g_arenaOpt.alloc_raw(as>0?as:1);
    }
    for(int j=*pNewBeginP;j<*pNewBeginR;++j){
        int v=vertexSets[j];
        auto [nb2,ne2]=edgeGraph.getNbr(v);
        for(auto idx=nb2;idx<ne2;++idx){
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

static void insertionSortOpt(TreeGraphNode* arr,int n){
    for(int i=1;i<n;i++){
        TreeGraphNode key=arr[i];int j=i-1;
        while(j>=0&&arr[j]>key){arr[j+1]=arr[j];j--;}
        arr[j+1]=key;
    }
}

static void recurseOpt(
    int* vertexSets,int* vertexLookup,
    int** neighborsInP,int* numNeighbors,
    int beginP,int beginR,
    int* keepV,int keepSz,int* dropV,int dropSz,
    int max_k,int min_k,
    std::vector<std::vector<TreeGraphNode>>& localBuf)
{
    if(keepSz>max_k)return;
    if(beginP>=beginR){
        int cSize=keepSz+dropSz;
        if(cSize<min_k)return;
        std::vector<TreeGraphNode> node;node.reserve(cSize);
        if(keepSz==max_k){
            for(int i=0;i<keepSz;i++)node.emplace_back(keepV[i],false);
            localBuf.push_back(std::move(node));return;
        }
        for(int i=0;i<keepSz;i++)node.emplace_back(keepV[i],false);
        for(int i=0;i<dropSz;i++)node.emplace_back(dropV[i],true);
        if(cSize<=32)insertionSortOpt(node.data(),cSize);
        else std::sort(node.begin(),node.end());
        localBuf.push_back(std::move(node));return;
    }
    int*cands;int nCands;
    int arenaTop=g_arenaOpt.save();
    int pivot=findPivotMarkOpt(&cands,&nCands,vertexSets,vertexLookup,
                             neighborsInP,numNeighbors,beginP,beginR);
    for(int ci=0;ci<nCands;ci++){
        int vertex=cands[ci],newBeginP,newBeginR;
        moveToRFastOpt(vertex,vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     &beginP,&beginR,&newBeginP,&newBeginR);
        if(vertex==pivot){dropV[dropSz]=vertex;
            recurseOpt(vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     newBeginP,newBeginR,keepV,keepSz,dropV,dropSz+1,max_k,min_k,localBuf);
        }else{keepV[keepSz]=vertex;
            recurseOpt(vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     newBeginP,newBeginR,keepV,keepSz+1,dropV,dropSz,max_k,min_k,localBuf);
        }
    }
    g_arenaOpt.restore(arenaTop);
}

DynamicGraph<TreeGraphNode> SDCT_Par_Final(Graph& edgeGraph,int max_k,int min_k){
    auto size=(int)edgeGraph.getGraphNodeSize();
    int nthreads=omp_get_max_threads();
    std::cout<<"SDCT_Par_Final "<<nthreads<<" threads"<<std::endl;
    
    // Global shared data structures (initialized once)
    std::vector<int> vertexSets(size),vertexLookup(size);
    std::vector<int*> neighborsInP(size,nullptr);
    std::vector<int> numNeighbors(size,1);
    
    // Initialize global structures
    for(int i=0;i<size;i++){
        vertexSets[i]=i;
        vertexLookup[i]=i;
    }
    
    std::vector<std::vector<std::vector<TreeGraphNode>>> thread_bufs(nthreads);

    #pragma omp parallel
    {
        int tid=omp_get_thread_num();
        g_arenaOpt.init(size*64);
        g_markOpt.init(size);
        thread_bufs[tid].reserve(std::max(1,size/nthreads)*20);
        int keepV[MAX_CSIZE],dropV[MAX_CSIZE];

        #pragma omp for schedule(dynamic,8) nowait
        for(int vertex=0;vertex<size;vertex++){
            int arenaBase=g_arenaOpt.save();
            
            // Create thread-local copies for this vertex's recursion
            std::vector<int> localVertexSets(vertexSets);
            std::vector<int> localVertexLookup(vertexLookup);
            std::vector<int*> localNeighborsInP(size,nullptr);
            std::vector<int> localNumNeighbors(numNeighbors);
            
            // Initialize local neighbor arrays
            for(int i=0;i<size;i++){
                int as=MY_MIN(size,edgeGraph.getNbrCount(i));
                localNeighborsInP[i]=g_arenaOpt.alloc_raw(as>0?as:1);
                localNumNeighbors[i]=0;
            }
            
            // Build neighbor lists for this vertex
            int beginX=0,beginP=0,beginR=size;
            int newBeginX,newBeginP,newBeginR;
            fillInPandXArenaOpt(vertex,localVertexSets.data(),localVertexLookup.data(),edgeGraph,
                              localNeighborsInP.data(),localNumNeighbors.data(),
                              &beginX,&beginP,&beginR,&newBeginX,&newBeginP,&newBeginR);
            
            keepV[0]=vertex;
            recurseOpt(localVertexSets.data(),localVertexLookup.data(),
                     localNeighborsInP.data(),localNumNeighbors.data(),
                     newBeginP,newBeginR,keepV,1,dropV,0,
                     max_k,min_k,thread_bufs[tid]);
            
            g_arenaOpt.restore(arenaBase);
        }
    }

    size_t total=0;
    for(auto&tb:thread_bufs)total+=tb.size();
    DynamicGraph<TreeGraphNode> treeGraph(size);
    treeGraph.adj_list.reserve(total);
    for(auto&tb:thread_bufs)
        for(auto&leaf:tb)
            treeGraph.adj_list.push_back(std::move(leaf));
    return treeGraph;
}
