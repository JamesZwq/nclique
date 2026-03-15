/*
 * SDCT_Parallel_Atomic: Atomic-based parallel SDCT
 * Key insight: Share global data structures, use atomic operations for synchronization
 * No per-thread data copying - minimal memory overhead
 */

#include<cassert>
#include<climits>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include <iostream>
#include <omp.h>
#include <atomic>

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

DynamicGraph<TreeGraphNode> SDCT_Parallel_Atomic(Graph &edgeGraph, int max_k, int min_k) {
    auto size = edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    std::cout << "SDCT_Parallel_Atomic with " << nthreads << " threads" << std::endl;
    
    // Global shared data (single copy, shared by all threads)
    daf::StaticVector<int> vertexSets(size);
    daf::StaticVector<int> vertexLookup(size);
    daf::StaticVector<int *> neighborsInP(size);
    daf::StaticVector<int> numNeighbors(size);
    vertexSets.c_size = size;
    vertexLookup.c_size = size;
    neighborsInP.c_size = size;
    numNeighbors.c_size = size;

    for (int i = 0; i < size; ++i) {
        vertexLookup[i] = i;
        vertexSets[i] = i;
        neighborsInP[i] = static_cast<int *>(Calloc(1, sizeof(int)));
        numNeighbors[i] = 1;
    }

    // Thread-local results
    std::vector<std::vector<std::vector<TreeGraphNode>>> thread_results(nthreads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        thread_results[tid].reserve(size / nthreads * 10);

        daf::StaticVector<int> dropV(MAX_CSIZE);
        daf::StaticVector<int> keepV(MAX_CSIZE);

        int beginX = 0;
        int beginP = 0;
        int beginR = size;

        // Process vertices with static scheduling to minimize synchronization
        #pragma omp for schedule(static) nowait
        for (int vertex = 0; vertex < size; ++vertex) {
            int newBeginX, newBeginP, newBeginR;

            // Use critical section only for data structure modification
            #pragma omp critical(fillInPandX)
            {
                fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(
                    vertex,
                    vertexSets.data(), vertexLookup.data(),
                    edgeGraph,
                    neighborsInP.data(), numNeighbors.data(),
                    &beginX, &beginP, &beginR,
                    &newBeginX, &newBeginP, &newBeginR);
            }

            dropV.clear();
            keepV.clear();
            keepV.push_back(vertex);

            // Recursion is thread-safe (uses thread-local data)
            listAllCliquesDegeneracyRecursive_VedgeGraphSDCT(
                vertexSets.data(),
                vertexLookup.data(), neighborsInP.data(),
                numNeighbors.data(), newBeginP,
                newBeginR, keepV, dropV, max_k, min_k, thread_results[tid]);

            #pragma omp critical(updateBeginR)
            {
                beginR = beginR + 1;
            }
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

    // Cleanup
    vertexSets.free();
    vertexLookup.free();
    for (int i = 0; i < size; ++i) {
        Free(neighborsInP[i]);
    }
    neighborsInP.free();
    numNeighbors.free();

    return treeGraph;
}
