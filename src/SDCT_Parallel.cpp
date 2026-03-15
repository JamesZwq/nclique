/*
 * SDCT_Parallel: True parallel version of SDCT
 * Key insight: The outer loop over vertices can be parallelized
 * Each thread processes different vertices independently
 */

#include<cassert>
#include<climits>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <mutex>

#include"degeneracy_algorithm_cliques_V.h"
#include"degeneracy_helper.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"misc.h"
#include"tree/MultiBranchTree.h"
#include <cstring>

extern double nCr[1001][401];

void listAllCliquesDegeneracyRecursive_VedgeGraphSDCT(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR, daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int min_k, std::vector<std::vector<TreeGraphNode>> &tree_buffer);

DynamicGraph<TreeGraphNode> SDCT_Parallel(Graph &edgeGraph, int max_k, int min_k) {
    auto size = edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    std::cout << "SDCT_Parallel with " << nthreads << " threads" << std::endl;
    
    // Global shared data structures (read-only after initialization)
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

    // Thread-local buffers for results
    std::vector<std::vector<std::vector<TreeGraphNode>>> thread_results(nthreads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        thread_results[tid].reserve(size / nthreads * 10);

        // Thread-local copies of mutable data
        daf::StaticVector<int> local_vertexSets(size);
        daf::StaticVector<int> local_vertexLookup(size);
        daf::StaticVector<int *> local_neighborsInP(size);
        daf::StaticVector<int> local_numNeighbors(size);
        local_vertexSets.c_size = size;
        local_vertexLookup.c_size = size;
        local_neighborsInP.c_size = size;
        local_numNeighbors.c_size = size;

        for (int i = 0; i < size; ++i) {
            local_vertexLookup[i] = i;
            local_vertexSets[i] = i;
            local_neighborsInP[i] = static_cast<int *>(Calloc(1, sizeof(int)));
            local_numNeighbors[i] = 1;
        }

        daf::StaticVector<int> dropV(MAX_CSIZE);
        daf::StaticVector<int> keepV(MAX_CSIZE);

        int beginX = 0;
        int beginP = 0;
        int beginR = size;

        #pragma omp for schedule(dynamic, 4) nowait
        for (int vertex = 0; vertex < size; ++vertex) {
            int newBeginX, newBeginP, newBeginR;

            fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(
                vertex,
                local_vertexSets.data(), local_vertexLookup.data(),
                edgeGraph,
                local_neighborsInP.data(), local_numNeighbors.data(),
                &beginX, &beginP, &beginR,
                &newBeginX, &newBeginP, &newBeginR);

            dropV.clear();
            keepV.clear();
            keepV.push_back(vertex);

            listAllCliquesDegeneracyRecursive_VedgeGraphSDCT(
                local_vertexSets.data(),
                local_vertexLookup.data(), local_neighborsInP.data(),
                local_numNeighbors.data(), newBeginP,
                newBeginR, keepV, dropV, max_k, min_k, thread_results[tid]);

            beginR = beginR + 1;
        }

        // Cleanup thread-local data
        for (int i = 0; i < size; ++i) {
            Free(local_neighborsInP[i]);
        }
        local_vertexSets.free();
        local_vertexLookup.free();
        local_neighborsInP.free();
        local_numNeighbors.free();
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
    vertexSets.free();
    vertexLookup.free();
    for (int i = 0; i < size; ++i) {
        Free(neighborsInP[i]);
    }
    neighborsInP.free();
    numNeighbors.free();

    return treeGraph;
}
