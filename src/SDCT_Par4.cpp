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

// Arena holding both permanent per-thread state (vertexSets/vertexLookup/
// neighborsInP/numNeighbors at the bottom) and BK scratch above. A single
// malloc per thread eliminates allocator-lock contention at high T.
struct Arena4 {
    int* base = nullptr;
    size_t top = 0;
    size_t cap = 0;
    size_t perm_end = 0;          // restore() never drops top below this
    void init(size_t n_ints) {
        if (base) std::free(base);                     // re-init safe
        cap = n_ints + 16;
        base = (int*)std::malloc(cap * sizeof(int));
        top = 0; perm_end = 0;
    }
    ~Arena4(){ std::free(base); }
    int* alloc_raw(int n) { if (n <= 0) n = 1; int* p = base + top; top += n; return p; }
    void mark_perm() { perm_end = top; }
    size_t save() { return top; }
    void restore(size_t t) { top = (t < perm_end ? perm_end : t); }
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
    size_t arenaTop=g_arena4.save();
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

    // Per-thread arena layout (single malloc, no std::vectors → no allocator-lock
    // contention at high T, no zero-fill of unused state):
    //   [neighborsInP ptrs (2*size ints) | vertexSets (size) | vertexLookup (size)
    //    | numNeighbors (size) | BK scratch ...]
    // numNeighbors / neighborsInP are written by fillInPandXArena4 before being
    // read inside recurse4, so no init is needed.
    const size_t perm_ints    = 5 * (size_t)size;
    const size_t scratch_ints = std::min<size_t>(
        (size_t)32 * 1024 * 1024,                                        // 128 MB hard cap
        std::max<size_t>(2u * 1024u * 1024u, (size_t)size * 8u));        // ≥ 8 MB, otherwise 32·size

    #pragma omp parallel
    {
        int tid=omp_get_thread_num();
        g_arena4.init(perm_ints + scratch_ints);
        g_mark4.init(size);

        // neighborsInP first → 16-byte aligned (malloc returns 16-aligned, top=0).
        int** neighborsInP = reinterpret_cast<int**>(g_arena4.alloc_raw(2 * size));
        int*  vertexSets   = g_arena4.alloc_raw(size);
        int*  vertexLookup = g_arena4.alloc_raw(size);
        int*  numNeighbors = g_arena4.alloc_raw(size);
        // Identity permutation: only state read before fillInPandXArena4 touches a vertex.
        for (int i = 0; i < size; ++i) { vertexSets[i] = i; vertexLookup[i] = i; }
        g_arena4.mark_perm();   // restore() will not pop below this point

        int beginX=0,beginP=0,beginR=size;
        thread_bufs[tid].reserve(std::max(1,size/nthreads)*20);
        int keepV[MAX_CSIZE],dropV[MAX_CSIZE];

        #pragma omp for schedule(dynamic,8) nowait
        for(int vertex=0;vertex<size;vertex++){
            size_t arenaBase=g_arena4.save();
            int newBeginX,newBeginP,newBeginR;
            fillInPandXArena4(vertex,vertexSets,vertexLookup,edgeGraph,
                              neighborsInP,numNeighbors,
                              &beginX,&beginP,&beginR,&newBeginX,&newBeginP,&newBeginR);
            keepV[0]=vertex;
            recurse4(vertexSets,vertexLookup,
                     neighborsInP,numNeighbors,
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
