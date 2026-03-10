/*
 * SDCT_Par2: High-performance parallel SDCT
 * Key optimizations:
 * 1. Thread-local arena allocator: eliminates calloc/free lock contention
 * 2. Dynamic scheduling (dynamic,1): fixes load imbalance
 * 3. Stack arrays for keepV/dropV: no heap in hot path
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

// ============================================================
// Thread-local arena allocator
// Stack discipline: save top, alloc, restore top = free all
// ============================================================
struct Arena {
    int* base = nullptr;
    int  top  = 0;
    int  cap  = 0;

    void init(int n) {
        cap  = n + 16;
        base = static_cast<int*>(std::malloc(cap * sizeof(int)));
        top  = 0;
    }
    ~Arena() { std::free(base); }

    int* alloc(int n) {
        if (n <= 0) n = 1;
        int* p = base + top;
        top += n;
        std::memset(p, 0, n * sizeof(int));
        return p;
    }
    int  save()        { return top; }
    void restore(int t){ top = t; }
};

static thread_local Arena g_arena;

// ============================================================
// Arena-based pivot: same logic as findBestPivotNonNeighborsDegeneracyCliques
// but uses arena instead of calloc
// ============================================================
static int findPivotArena(
    int** pivotNonNeighbors, int* numNonNeighbors,
    int* vertexSets, int* vertexLookup,
    int** neighborsInP, int* numNeighbors,
    int beginP, int beginR)
{
    int pivot = -1;
    int maxIntersectionSize = -1;

    for (int j = beginP; j < beginR; j++) {
        int vertex = vertexSets[j];
        int numPotential = MY_MIN(beginR - beginP, numNeighbors[vertex]);
        int numInP = 0;
        for (int k = 0; k < numPotential; k++) {
            int loc = vertexLookup[neighborsInP[vertex][k]];
            if (loc >= beginP && loc < beginR) numInP++;
            else break;
        }
        if (numInP > maxIntersectionSize) {
            maxIntersectionSize = numInP;
            pivot = vertex;
        }
    }

    int sz = beginR - beginP;
    int* buf = g_arena.alloc(sz);
    std::memcpy(buf, &vertexSets[beginP], sz * sizeof(int));
    *numNonNeighbors = sz;

    int numPivotNbr = MY_MIN(sz, numNeighbors[pivot]);
    for (int j = 0; j < numPivotNbr; j++) {
        int loc = vertexLookup[neighborsInP[pivot][j]];
        if (loc >= beginP && loc < beginR)
            buf[loc - beginP] = -1;
        else break;
    }

    for (int j = 0; j < *numNonNeighbors; ) {
        if (buf[j] == -1) {
            buf[j] = buf[--(*numNonNeighbors)];
        } else {
            j++;
        }
    }

    *pivotNonNeighbors = buf;
    return pivot;
}

// ============================================================
// Arena-based fillInPandX
// ============================================================
static void fillInPandXArena(
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
        neighborsInP[v] = g_arena.alloc(allocSz > 0 ? allocSz : 1);
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

// ============================================================
// Recursive BK kernel
// ============================================================
static void recurse2(
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
        std::sort(node.begin(), node.end());
        localBuf.push_back(std::move(node));
        return;
    }

    int* cands;
    int  nCands;
    int arenaTop = g_arena.save();
    int pivot = findPivotArena(&cands, &nCands,
                               vertexSets, vertexLookup,
                               neighborsInP, numNeighbors,
                               beginP, beginR);

    for (int ci = 0; ci < nCands; ci++) {
        int vertex = cands[ci];
        int newBeginP, newBeginR;
        moveToRDegeneracyCliques(vertex,
                                 vertexSets, vertexLookup,
                                 neighborsInP, numNeighbors,
                                 &beginP, &beginR,
                                 &newBeginP, &newBeginR);
        if (vertex == pivot) {
            dropV[dropSz] = vertex;
            recurse2(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR, keepV, keepSz, dropV, dropSz+1,
                     max_k, min_k, localBuf);
        } else {
            keepV[keepSz] = vertex;
            recurse2(vertexSets, vertexLookup, neighborsInP, numNeighbors,
                     newBeginP, newBeginR, keepV, keepSz+1, dropV, dropSz,
                     max_k, min_k, localBuf);
        }
    }
    g_arena.restore(arenaTop);
}

// ============================================================
// SDCT_Par2: main entry
// ============================================================
DynamicGraph<TreeGraphNode> SDCT_Par2(Graph& edgeGraph, int max_k, int min_k) {
    auto size = (int)edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    std::cout << "SDCT_Par2 " << nthreads << " threads" << std::endl;

    std::vector<std::vector<std::vector<TreeGraphNode>>> thread_bufs(nthreads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        // Arena size: each vertex can have up to size neighbours;
        // each recursion level allocs up to 'size' ints for pivot buf
        // + up to degree entries for neighborsInP.
        // Use 64x size as a safe upper bound.
        g_arena.init(size * 64);

        std::vector<int>  vertexSets(size);
        std::vector<int>  vertexLookup(size);
        std::vector<int*> neighborsInP(size, nullptr);
        std::vector<int>  numNeighbors(size, 1);
        for (int i = 0; i < size; i++) {
            vertexSets[i]   = i;
            vertexLookup[i] = i;
            neighborsInP[i] = g_arena.alloc(1);
            numNeighbors[i] = 1;
        }

        int beginX = 0, beginP = 0, beginR = size;
        thread_bufs[tid].reserve(std::max(1, size / nthreads) * 20);

        int keepV[MAX_CSIZE], dropV[MAX_CSIZE];

        #pragma omp for schedule(dynamic, 1) nowait
        for (int vertex = 0; vertex < size; vertex++) {
            int arenaBase = g_arena.save();

            int newBeginX, newBeginP, newBeginR;
            fillInPandXArena(vertex,
                             vertexSets.data(), vertexLookup.data(),
                             edgeGraph,
                             neighborsInP.data(), numNeighbors.data(),
                             &beginX, &beginP, &beginR,
                             &newBeginX, &newBeginP, &newBeginR);

            keepV[0] = vertex;
            recurse2(vertexSets.data(), vertexLookup.data(),
                     neighborsInP.data(), numNeighbors.data(),
                     newBeginP, newBeginR,
                     keepV, 1, dropV, 0,
                     max_k, min_k,
                     thread_bufs[tid]);

            g_arena.restore(arenaBase);
            beginR++;
        }
    } // end parallel

    // Merge thread buffers
    size_t total = 0;
    for (auto& tb : thread_bufs) total += tb.size();
    DynamicGraph<TreeGraphNode> treeGraph(size);
    treeGraph.adj_list.reserve(total);
    for (auto& tb : thread_bufs)
        for (auto& leaf : tb)
            treeGraph.adj_list.push_back(std::move(leaf));
    return treeGraph;
}
