/*
 * SDCT_Par3: Further optimizations over SDCT_Par2
 * 1. No memset in arena: raw alloc for neighborsInP
 * 2. Bitmap-based pivot neighbor marking: O(degree) instead of O(P*degree)
 * 3. Insertion sort for small leaf nodes
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

// Arena without zero-init
struct Arena3 {
    int* base = nullptr;
    int  top  = 0;
    int  cap  = 0;
    void init(int n) {
        cap  = n + 16;
        base = static_cast<int*>(std::malloc(cap * sizeof(int)));
        top  = 0;
    }
    ~Arena3() { std::free(base); }
    int* alloc_raw(int n)  { if (n<=0) n=1; int* p=base+top; top+=n; return p; }
    int* alloc_zero(int n) { if (n<=0) n=1; int* p=base+top; top+=n; std::memset(p,0,n*sizeof(int)); return p; }
    int  save()         { return top; }
    void restore(int t) { top = t; }
};
static thread_local Arena3 g_arena3;

// Generation-stamped bitmap — O(1) mark/clear, no reset needed
struct BitMark {
    std::vector<int> mark;
    int gen = 0;
    void init(int n) { mark.assign(n, 0); }
    void next()      { ++gen; }
    void set(int v)  { mark[v] = gen; }
    bool get(int v)  { return mark[v] == gen; }
};
static thread_local BitMark g_bitmark;

// Pivot finder with bitmap marking
static int findPivotBitmap(
    int** pivotNonNeighbors, int* numNonNeighbors,
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int beginP, int beginR)
{
    int pivot = -1, maxInP = -1;
    int sz = beginR - beginP;

    for (int j = beginP; j < beginR; j++) {
        int v = vertexSets[j];
        int numPotential = MY_MIN(sz, numNeighbors[v]);
        int numInP = 0;
        for (int k = 0; k < numPotential; k++) {
            int loc = vertexLookup[neighborsInP[v][k]];
            if (loc >= beginP && loc < beginR) numInP++;
            else break;
        }
        if (numInP > maxInP) { maxInP = numInP; pivot = v; }
    }

    // Mark pivot neighbors with bitmap
    g_bitmark.next();
    int numPivotNbr = MY_MIN(sz, numNeighbors[pivot]);
    for (int j = 0; j < numPivotNbr; j++) {
        int nb  = neighborsInP[pivot][j];
        int loc = vertexLookup[nb];
        if (loc >= beginP && loc < beginR) g_bitmark.set(nb);
        else break;
    }

    // Build P \ N(pivot) without memset
    int* buf = g_arena3.alloc_raw(sz);
    int numCands = 0;
    for (int j = beginP; j < beginR; j++) {
        int v = vertexSets[j];
        if (!g_bitmark.get(v)) buf[numCands++] = v;
    }
    *pivotNonNeighbors = buf;
    *numNonNeighbors   = numCands;
    return pivot;
}

// fillInPandX with raw alloc
static void fillInPandXArena3(
    int vertex,
    int* vertexSets, int* vertexLookup,
    Graph& edgeGraph,
    int** neighborsInP, int* numNeighbors,
    int* pBeginX, int* pBeginP, int* pBeginR,
    int* pNewBeginX, int* pNewBeginP, int* pNewBeginR)
{
    int vertexLocation = vertexLookup[vertex];
    (*pBeginR)--;
    vertexSets[vertexLocation] = vertexSets[*pBeginR];
    vertexLookup[vertexSets[*pBeginR]] = vertexLocation;
    vertexSets[*pBeginR] = vertex;
    vertexLookup[vertex] = *pBeginR;
    *pNewBeginR = *pBeginR;
    *pNewBeginP = *pBeginR;

    auto [nbr_begin, nbr_end] = edgeGraph.getNbr(vertex);
    for (int idx = nbr_begin; idx < nbr_end; ++idx) {
        int neighbor = edgeGraph.adj_list[idx];
        if (neighbor <= vertex) continue;
        int neighborLocation = vertexLookup[neighbor];
        (*pNewBeginP)--;
        vertexSets[neighborLocation] = vertexSets[*pNewBeginP];
        vertexLookup[vertexSets[*pNewBeginP]] = neighborLocation;
        vertexSets[*pNewBeginP] = neighbor;
        vertexLookup[neighbor] = *pNewBeginP;
    }
    *pNewBeginX = *pNewBeginP;

    int pSize = *pNewBeginR - *pNewBeginP;
    for (int j = *pNewBeginP; j < *pNewBeginR; ++j) {
        int v = vertexSets[j];
        numNeighbors[v] = 0;
        int allocSz = MY_MIN(pSize, edgeGraph.getNbrCount(v));
        neighborsInP[v] = g_arena3.alloc_raw(allocSz > 0 ? allocSz : 1);
    }
    for (int j = *pNewBeginP; j < *pNewBeginR; ++j) {
        int v = vertexSets[j];
        auto [nb, ne] = edgeGraph.getNbr(v);
        for (auto idx = nb; idx < ne; ++idx) {
            int w = edgeGraph.adj_list[idx];
            if (w <= v) continue;
            int wloc = vertexLookup[w];
            if (wloc >= *pNewBeginP && wloc < *pNewBeginR) {
                neighborsInP[v][numNeighbors[v]++] = w;
                neighborsInP[w][numNeighbors[w]++] = v;
            }
        }
    }
}

// Insertion sort for tiny arrays
static void insertionSort(TreeGraphNode* arr, int n) {
    for (int i = 1; i < n; i++) {
        TreeGraphNode key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > key) { arr[j+1] = arr[j]; j--; }
        arr[j+1] = key;
    }
}

// Recursive BK kernel
static void recurse3(
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int beginP, int beginR,
    int* keepV, int keepSz,
    int* dropV, int dropSz,
    int max_k, int min_k,
    std::vector<std::vector<TreeGraphNode>>& localBuf)
{
    if (keepSz > max_k) return;
    if (beginP >= beginR) {
        int cSize = keepSz + dropSz;
        if (cSize < min_k) return;
        std::vector<TreeGraphNode> node;
        node.reserve(cSize);
        if (keepSz == max_k) {
            for (int i = 0; i < keepSz; i++) node.emplace_back(keepV[i], false);
            localBuf.push_back(std::move(node));
            return;
        }
        for (int i = 0; i < keepSz; i++) node.emplace_back(keepV[i], false);
        for (int i = 0; i < dropSz; i++) node.emplace_back(dropV[i], true);
        if (cSize <= 32) insertionSort(node.data(), cSize);
        else std::sort(node.begin(), node.end());
        localBuf.push_back(std::move(node));
        return;
    }

    int* cands; int nCands;
    int arenaTop = g_arena3.save();
    int pivot = findPivotBitmap(&cands, &nCands,
                                vertexSets, vertexLookup,
                                neighborsInP, numNeighbors,
                                beginP, beginR);
    for (int ci = 0; ci < nCands; ci++) {
        int vertex = cands[ci];
        int newBeginP, newBeginR;
        moveToRDegeneracyCliques(vertex, vertexSets, vertexLookup,
                                 neighborsInP, numNeighbors,
                                 &beginP, &beginR, &newBeginP, &newBeginR);
        if (vertex == pivot) {
            dropV[dropSz] = vertex;
            recurse3(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR, keepV, keepSz, dropV, dropSz+1,
                     max_k, min_k, localBuf);
        } else {
            keepV[keepSz] = vertex;
            recurse3(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR, keepV, keepSz+1, dropV, dropSz,
                     max_k, min_k, localBuf);
        }
    }
    g_arena3.restore(arenaTop);
}

// Main entry: SDCT_Par3
DynamicGraph<TreeGraphNode> SDCT_Par3(Graph& edgeGraph, int max_k, int min_k) {
    auto size = (int)edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    std::cout << "SDCT_Par3 " << nthreads << " threads" << std::endl;

    std::vector<std::vector<std::vector<TreeGraphNode>>> thread_bufs(nthreads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        g_arena3.init(size * 64);
        g_bitmark.init(size);

        std::vector<int>  vertexSets(size);
        std::vector<int>  vertexLookup(size);
        std::vector<int*> neighborsInP(size, nullptr);
        std::vector<int>  numNeighbors(size, 1);
        for (int i = 0; i < size; i++) {
            vertexSets[i]   = i;
            vertexLookup[i] = i;
            neighborsInP[i] = g_arena3.alloc_raw(1);
            numNeighbors[i] = 1;
        }

        int beginX = 0, beginP = 0, beginR = size;
        thread_bufs[tid].reserve(std::max(1, size / nthreads) * 20);
        int keepV[MAX_CSIZE], dropV[MAX_CSIZE];

        #pragma omp for schedule(dynamic, 1) nowait
        for (int vertex = 0; vertex < size; vertex++) {
            int arenaBase = g_arena3.save();
            int newBeginX, newBeginP, newBeginR;
            fillInPandXArena3(vertex,
                              vertexSets.data(), vertexLookup.data(), edgeGraph,
                              neighborsInP.data(), numNeighbors.data(),
                              &beginX, &beginP, &beginR,
                              &newBeginX, &newBeginP, &newBeginR);
            keepV[0] = vertex;
            recurse3(vertexSets.data(), vertexLookup.data(),
                     neighborsInP.data(), numNeighbors.data(),
                     newBeginP, newBeginR, keepV, 1, dropV, 0,
                     max_k, min_k, thread_bufs[tid]);
            g_arena3.restore(arenaBase);
            beginR++;
        }
    }

    size_t total = 0;
    for (auto& tb : thread_bufs) total += tb.size();
    DynamicGraph<TreeGraphNode> treeGraph(size);
    treeGraph.adj_list.reserve(total);
    for (auto& tb : thread_bufs)
        for (auto& leaf : tb)
            treeGraph.adj_list.push_back(std::move(leaf));
    return treeGraph;
}
