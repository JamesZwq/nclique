/*
 * SDCT_Parallel_Lightweight: Lightweight parallel SDCT
 * Key insight: Minimize per-thread memory allocation
 * Use stack-based allocation instead of heap
 */

#include<cassert>
#include<climits>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include <iostream>
#include <omp.h>
#include <cstring>

#include"degeneracy_algorithm_cliques_V.h"
#include"degeneracy_helper.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"misc.h"
#include"tree/MultiBranchTree.h"

extern double nCr[1001][401];

void listAllCliquesDegeneracyRecursive_VedgeGraphSDCT(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR, daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int min_k, std::vector<std::vector<TreeGraphNode>> &tree_buffer);

DynamicGraph<TreeGraphNode> SDCT_Parallel_Lightweight(Graph &edgeGraph, int max_k, int min_k) {
    auto size = edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    std::cout << "SDCT_Parallel_Lightweight with " << nthreads << " threads" << std::endl;
    
    // Global shared data (read-only after initialization)
    std::vector<int> globalVertexSets(size);
    std::vector<int> globalVertexLookup(size);
    std::vector<int *> globalNeighborsInP(size);
    std::vector<int> globalNumNeighbors(size, 1);

    for (int i = 0; i < size; ++i) {
        globalVertexLookup[i] = i;
        globalVertexSets[i] = i;
        globalNeighborsInP[i] = static_cast<int *>(Calloc(1, sizeof(int)));
    }

    // Thread-local results
    std::vector<std::vector<std::vector<TreeGraphNode>>> thread_results(nthreads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        thread_results[tid].reserve(size / nthreads * 10);

        // Allocate thread-local data on stack (minimal)
        int* vertexSets = (int*)alloca(size * sizeof(int));
        int* vertexLookup = (int*)alloca(size * sizeof(int));
        int** neighborsInP = (int**)alloca(size * sizeof(int*));
        int* numNeighbors = (int*)alloca(size * sizeof(int));

        // Copy global data
        memcpy(vertexSets, globalVertexSets.data(), size * sizeof(int));
        memcpy(vertexLookup, globalVertexLookup.data(), size * sizeof(int));
        memcpy(numNeighbors, globalNumNeighbors.data(), size * sizeof(int));

        for (int i = 0; i < size; ++i) {
            neighborsInP[i] = static_cast<int *>(Calloc(1, sizeof(int)));
        }

        daf::StaticVector<int> dropV(MAX_CSIZE);
        daf::StaticVector<int> keepV(MAX_CSIZE);

        int beginX = 0;
        int beginP = 0;
        int beginR = size;

        // Process vertices with large chunk size
        #pragma omp for schedule(dynamic, 64) nowait
        for (int vertex = 0; vertex < size; ++vertex) {
            int newBeginX, newBeginP, newBeginR;

            fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(
                vertex,
                vertexSets, vertexLookup,
                edgeGraph,
                neighborsInP, numNeighbors,
                &beginX, &beginP, &beginR,
                &newBeginX, &newBeginP, &newBeginR);

            dropV.clear();
            keepV.clear();
            keepV.push_back(vertex);

            listAllCliquesDegeneracyRecursive_VedgeGraphSDCT(
                vertexSets,
                vertexLookup, neighborsInP,
                numNeighbors, newBeginP,
                newBeginR, keepV, dropV, max_k, min_k, thread_results[tid]);

            beginR = beginR + 1;
        }

        // Cleanup
        for (int i = 0; i < size; ++i) {
            Free(neighborsInP[i]);
        }
    }

    // Merge results from all threads
    DynamicGraph<TreeGraphNode> treeGraph(size);
    size_t total_cliques = 0;
    for (int tid = 0; tid < nthreads; ++tid) {
        total_cliques += thread_results[tid].size();
    }
    treeGraph.adj_list.reserve(total_cliques);

    for (int tid = 0; tid < nthreads; ++tid) {
        for (auto &clique : thread_results[tid]) {
            treeGraph.adj_list.push_back(std::move(clique));
        }
    }

    // Cleanup global data
    for (int i = 0; i < size; ++i) {
        Free(globalNeighborsInP[i]);
    }

    return treeGraph;
}
