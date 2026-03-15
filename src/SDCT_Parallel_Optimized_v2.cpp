/*
 * SDCT_Parallel_Optimized_v2: Work-stealing parallel SDCT
 * Key insight: Use global shared data with minimal synchronization
 * Process vertices in batches, each batch is independent
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

DynamicGraph<TreeGraphNode> SDCT_Parallel_Optimized_v2(Graph &edgeGraph, int max_k, int min_k) {
    auto size = edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    std::cout << "SDCT_Parallel_Optimized_v2 with " << nthreads << " threads" << std::endl;
    
    // Global shared data (initialized once)
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
        thread_results[tid].reserve(size / nthreads * 10);

        // Thread-local copies (minimal)
        std::vector<int> local_vertexSets(globalVertexSets.data(), globalVertexSets.data() + size);
        std::vector<int> local_vertexLookup(globalVertexLookup.data(), globalVertexLookup.data() + size);
        std::vector<int *> local_neighborsInP(size);
        std::vector<int> local_numNeighbors(globalNumNeighbors.data(), globalNumNeighbors.data() + size);

        for (int i = 0; i < size; ++i) {
            local_neighborsInP[i] = static_cast<int *>(Calloc(1, sizeof(int)));
        }

        daf::StaticVector<int> dropV(MAX_CSIZE);
        daf::StaticVector<int> keepV(MAX_CSIZE);

        int beginX = 0;
        int beginP = 0;
        int beginR = size;

        // Process vertices with guided scheduling for better load balancing
        #pragma omp for schedule(guided) nowait
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

        // Cleanup
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
