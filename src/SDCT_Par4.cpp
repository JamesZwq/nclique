/*
 * SDCT_Par4: Fast moveToR with mark array
 * moveToR O(P*deg) -> O(P + deg) using generation-stamped mark
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

struct Arena4 {
    int* base=nullptr; int top=0; int cap=0;
    void init(int n){cap=n+16;base=(int*)std::malloc(cap*sizeof(int));top=0;}
    ~Arena4(){std::free(base);}
    int* alloc_raw(int n){if(n<=0)n=1;int*p=base+top;top+=n;return p;}
    int save(){return top;}
    void restore(int t){top=t;}
};
static thread_local Arena4 g_arena4;

struct MarkArray4 {
    std::vector<int> mark;
    int gen=0;
    void init(int n){mark.assign(n,0);}
    void next(){++gen;}
    void set(int v){mark[v]=gen;}
    bool get(int v)const{return mark[v]==gen;}
};
static thread_local MarkArray4 g_mark4;

// Pivot finder with bitmap
static int findPivotMark4(
    int** pivotNonNeighbors, int* numNonNeighbors,
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int beginP, int beginR)
{
    int pivot=-1, maxInP=-1, sz=beginR-beginP;
    for(int j=beginP;j<beginR;j++){
        int v=vertexSets[j];
        int numPotential=MY_MIN(sz,numNeighbors[v]);
        int numInP=0;
        for(int k=0;k<numPotential;k++){
            int loc=vertexLookup[neighborsInP[v][k]];
            if(loc>=beginP&&loc<beginR)numInP++; else break;
        }
        if(numInP>maxInP){maxInP=numInP;pivot=v;}
    }
    g_mark4.next();
    int nPivNbr=MY_MIN(sz,numNeighbors[pivot]);
    for(int j=0;j<nPivNbr;j++){
        int nb=neighborsInP[pivot][j];
        int loc=vertexLookup[nb];
        if(loc>=beginP&&loc<beginR)g_mark4.set(nb); else break;
    }
    int*buf=g_arena4.alloc_raw(sz);
    int nc=0;
    for(int j=beginP;j<beginR;j++){
        int v=vertexSets[j];
        if(!g_mark4.get(v))buf[nc++]=v;
    }
    *pivotNonNeighbors=buf; *numNonNeighbors=nc;
    return pivot;
}

// Fast moveToR: mark vertex's neighbors O(deg), then partition P O(|P|)
// Replaces original O(|P|*deg) double loop
static void moveToRFast4(
    int vertex,
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int* pBeginP, int* pBeginR,
    int* pNewBeginP, int* pNewBeginR)
{
    int vertexLocation=vertexLookup[vertex];
    (*pBeginR)--;
    vertexSets[vertexLocation]=vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]]=vertexLocation;
    vertexSets[*pBeginR]=vertex;
    vertexLookup[vertex]=*pBeginR;

    *pNewBeginP=*pBeginP;
    *pNewBeginR=*pBeginP;

    int sizeOfP=*pBeginR-*pBeginP;

    // Mark all neighbors of vertex in P using mark array (new round)
    g_mark4.next();
    int numNbr=MY_MIN(sizeOfP,numNeighbors[vertex]);
    for(int k=0;k<numNbr;k++)
        g_mark4.set(neighborsInP[vertex][k]);

    // Partition P: neighbors of vertex go into newP
    // Must scan from pBeginP and track carefully
    int j = *pBeginP;
    while (j < *pBeginR) {
        int neighbor = vertexSets[j];
        int neighborLocation = j;
        int numPotentialNeighbors = MY_MIN(sizeOfP, numNeighbors[neighbor]);
        // Check if neighbor is adjacent to vertex using mark (O(1))
        bool isAdj = false;
        // We need to check: is 'neighbor' in neighborsInP[vertex]?
        // We already marked all neighbors of vertex, check via mark
        if (g_mark4.get(neighbor)) {
            isAdj = true;
        }
        if (isAdj) {
            vertexSets[neighborLocation] = vertexSets[*pNewBeginR];
            vertexLookup[vertexSets[*pNewBeginR]] = neighborLocation;
            vertexSets[*pNewBeginR] = neighbor;
            vertexLookup[neighbor] = *pNewBeginR;
            (*pNewBeginR)++;
        }
        j++;
    }

    // Compact neighborsInP (identical to original: reorder but don't update numNeighbors)
    for (int j2 = *pNewBeginP; j2 < *pNewBeginR; j2++) {
        int tv = vertexSets[j2];
        int numPotential = MY_MIN(sizeOfP, numNeighbors[tv]);
        int cnt = 0;
        for (int k = 0; k < numPotential; k++) {
            int nbr = neighborsInP[tv][k];
            int nloc = vertexLookup[nbr];
            if (nloc >= *pNewBeginP && nloc < *pNewBeginR) {
                // swap to front, same as original
                neighborsInP[tv][k] = neighborsInP[tv][cnt];
                neighborsInP[tv][cnt] = nbr;
                cnt++;
            }
        }
        // NOTE: do NOT update numNeighbors[tv] — original doesn't either
    }
}

// fillInPandX with raw alloc
static void fillInPandXArena4(
    int vertex, int* vertexSets, int* vertexLookup,
    Graph& edgeGraph, int** neighborsInP, int* numNeighbors,
    int* pBeginX, int* pBeginP, int* pBeginR,
    int* pNewBeginX, int* pNewBeginP, int* pNewBeginR)
{
    int vl=vertexLookup[vertex];
    (*pBeginR)--;
    vertexSets[vl]=vertexSets[*pBeginR]; vertexLookup[vertexSets[*pBeginR]]=vl;
    vertexSets[*pBeginR]=vertex; vertexLookup[vertex]=*pBeginR;
    *pNewBeginR=*pBeginR; *pNewBeginP=*pBeginR;

    auto [nb,ne]=edgeGraph.getNbr(vertex);
    for(int idx=nb;idx<ne;++idx){
        int neighbor=edgeGraph.adj_list[idx];
        if(neighbor<=vertex)continue;
        int nl=vertexLookup[neighbor]; (*pNewBeginP)--;
        vertexSets[nl]=vertexSets[*pNewBeginP]; vertexLookup[vertexSets[*pNewBeginP]]=nl;
        vertexSets[*pNewBeginP]=neighbor; vertexLookup[neighbor]=*pNewBeginP;
    }
    *pNewBeginX=*pNewBeginP;

    int pSize=*pNewBeginR-*pNewBeginP;
    for(int j=*pNewBeginP;j<*pNewBeginR;++j){
        int v=vertexSets[j]; numNeighbors[v]=0;
        int as=MY_MIN(pSize,edgeGraph.getNbrCount(v));
        neighborsInP[v]=g_arena4.alloc_raw(as>0?as:1);
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

static void insertionSort4(TreeGraphNode* arr,int n){
    for(int i=1;i<n;i++){
        TreeGraphNode key=arr[i]; int j=i-1;
        while(j>=0&&arr[j]>key){arr[j+1]=arr[j];j--;}
        arr[j+1]=key;
    }
}

static void recurse4(
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int beginP, int beginR,
    int* keepV, int keepSz, int* dropV, int dropSz,
    int max_k, int min_k,
    std::vector<std::vector<TreeGraphNode>>& localBuf)
{
    if(keepSz>max_k)return;
    if(beginP>=beginR){
        int cSize=keepSz+dropSz;
        if(cSize<min_k)return;
        std::vector<TreeGraphNode> node; node.reserve(cSize);
        if(keepSz==max_k){
            for(int i=0;i<keepSz;i++)node.emplace_back(keepV[i],false);
            localBuf.push_back(std::move(node));return;
        }
        for(int i=0;i<keepSz;i++)node.emplace_back(keepV[i],false);
        for(int i=0;i<dropSz;i++)node.emplace_back(dropV[i],true);
        if(cSize<=32)insertionSort4(node.data(),cSize);
        else std::sort(node.begin(),node.end());
        localBuf.push_back(std::move(node));return;
    }

    int*cands; int nCands;
    int arenaTop=g_arena4.save();
    int pivot=findPivotMark4(&cands,&nCands,
                             vertexSets,vertexLookup,
                             neighborsInP,numNeighbors,beginP,beginR);
    for(int ci=0;ci<nCands;ci++){
        int vertex=cands[ci];
        int newBeginP,newBeginR;
        moveToRFast4(vertex,vertexSets,vertexLookup,
                     neighborsInP,numNeighbors,
                     &beginP,&beginR,&newBeginP,&newBeginR);
        if(vertex==pivot){
            dropV[dropSz]=vertex;
            recurse4(vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     newBeginP,newBeginR,keepV,keepSz,dropV,dropSz+1,
                     max_k,min_k,localBuf);
        }else{
            keepV[keepSz]=vertex;
            recurse4(vertexSets,vertexLookup,neighborsInP,numNeighbors,
                     newBeginP,newBeginR,keepV,keepSz+1,dropV,dropSz,
                     max_k,min_k,localBuf);
        }
    }
    g_arena4.restore(arenaTop);
}

DynamicGraph<TreeGraphNode> SDCT_Par4(Graph& edgeGraph,int max_k,int min_k){
    auto size=(int)edgeGraph.getGraphNodeSize();
    int nthreads=omp_get_max_threads();
    std::cout<<"SDCT_Par4 "<<nthreads<<" threads"<<std::endl;

    std::vector<std::vector<std::vector<TreeGraphNode>>> thread_bufs(nthreads);

    #pragma omp parallel
    {
        int tid=omp_get_thread_num();
        g_arena4.init(size*64);
        g_mark4.init(size);

        std::vector<int> vertexSets(size),vertexLookup(size);
        std::vector<int*> neighborsInP(size,nullptr);
        std::vector<int> numNeighbors(size,1);
        for(int i=0;i<size;i++){
            vertexSets[i]=i; vertexLookup[i]=i;
            neighborsInP[i]=g_arena4.alloc_raw(1);
            numNeighbors[i]=1;
        }
        int beginX=0,beginP=0,beginR=size;
        thread_bufs[tid].reserve(std::max(1,size/nthreads)*20);
        int keepV[MAX_CSIZE],dropV[MAX_CSIZE];

        #pragma omp for schedule(dynamic,1) nowait
        for(int vertex=0;vertex<size;vertex++){
            int arenaBase=g_arena4.save();
            int newBeginX,newBeginP,newBeginR;
            fillInPandXArena4(vertex,vertexSets.data(),vertexLookup.data(),edgeGraph,
                              neighborsInP.data(),numNeighbors.data(),
                              &beginX,&beginP,&beginR,&newBeginX,&newBeginP,&newBeginR);
            keepV[0]=vertex;
            recurse4(vertexSets.data(),vertexLookup.data(),
                     neighborsInP.data(),numNeighbors.data(),
                     newBeginP,newBeginR,keepV,1,dropV,0,
                     max_k,min_k,thread_bufs[tid]);
            g_arena4.restore(arenaBase);
            beginR++;
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
