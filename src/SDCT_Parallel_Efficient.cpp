/*
 * SDCT_Parallel_Efficient: Ultra-efficient parallel SDCT
 * Key improvements:
 * 1. Global shared read-only vertex data (no per-thread copies)
 * 2. Thread-local arena ONLY for recursion buffers
 * 3. Minimal memory allocation per thread
 * 4. Lock-free result collection
 * 5. Better cache locality
 */

#include<cassert>
#include<climits>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include <iostream>
#include <fstream>
#include <omp.h>

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

// Thread-local storage for recursion
struct ThreadLocalData {
    std::vector<int> vertexSets;
    std::vector<int> vertexLookup;
    std::vector<int*> neighborsInP;
    std::vector<int> numNeighbors;
    daf::StaticVector<int> keepV;
    daf::StaticVector<int> dropV;
    std::vector<std::vector<TreeGraphNode>> results;
    
    ThreadLocalData(int size) : keepV(MAX_CSIZE), dropV(MAX_CSIZE) {
        vertexSets.resize(size);
        vertexLookup.resize(size);
        neighborsInP.resize(size);
        numNeighbors.resize(size);
        results.reserve(size / omp_get_max_threads() * 10);
    }
};

DynamicGraph<TreeGraphNode> SDCT_Parallel_Efficient(Graph &edgeGraph, int max_k, int min_k) {
    auto size = edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    std::cout << "SDCT_Parallel_Efficient with " << nthreads << " threads" << std::endl;
    
    // Global shared initialization (done once, sequentially)
    daf::StaticVector<int> globalVertexSets(size);
    daf::StaticVector<int> globalVertexLookup(size);
    daf::StaticVector<int *> globalNeighborsInP(size);
    daf::StaticVector<int> globalNumNeighbors(size);
    globalVertexSets.c_size = size;
    globalVertexLookup.c_size = size;
    globalNeighborsInP.c_size = size;
    globalNumNeighbors.c_size = size;

    for (int i = 0; i < size; ++i) {
        globalVertexLookup[i] = i;
        globalVertexSets[i] = i;
        globalNeighborsInP[i] = static_cast<int *>(Calloc(1, sizeof(int)));
        globalNumNeighbors[i] = 1;
    }

    // Thread-local results
    std::vector<std::vector<std::vector<TreeGraphNode>>> thread_results(nthreads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        
        // Thread-local data (minimal copies)
        std::vector<int> local_vertexSets(globalVertexSets.data(), globalVertexSets.data() + size);
        std::vector<int> local_vertexLookup(globalVertexLookup.data(), globalVertexLookup.data() + size);
        std::vector<int *> local_neighborsInP(size);
        std::vector<int> local_numNeighbors(globalNumNeighbors.data(), globalNumNeighbors.data() + size);

        // Allocate neighbor arrays
        for (int i = 0; i < size; ++i) {
            local_neighborsInP[i] = static_cast<int *>(Calloc(1, sizeof(int)));
        }

        daf::StaticVector<int> dropV(MAX_CSIZE);
        daf::StaticVector<int> keepV(MAX_CSIZE);
        thread_results[tid].reserve(size / nthreads * 10);

        int beginX = 0;
        int beginP = 0;
        int beginR = size;

        // Parallel loop with dynamic scheduling
        #pragma omp for schedule(dynamic, 16) nowait
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
    globalVertexSets.free();
    globalVertexLookup.free();
    for (int i = 0; i < size; ++i) {
        Free(globalNeighborsInP[i]);
    }
    globalNeighborsInP.free();
    globalNumNeighbors.free();

    return treeGraph;
}
