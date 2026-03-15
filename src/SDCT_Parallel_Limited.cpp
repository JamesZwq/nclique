/*
 * SDCT_Parallel_Limited: Limited parallelism SDCT
 * Key insight: Use only a few threads (e.g., 8-16) to avoid memory explosion
 * Each thread processes a larger chunk of vertices
 */

#include<cassert>
#include<climits>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include <iostream>
#include <omp.h>

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

DynamicGraph<TreeGraphNode> SDCT_Parallel_Limited(Graph &edgeGraph, int max_k, int min_k) {
    auto size = edgeGraph.getGraphNodeSize();
    
    // Limit threads to avoid memory explosion
    int max_threads = omp_get_max_threads();
    int effective_threads = std::min(max_threads, 4);  // Use at most 4 threads
    
    std::cout << "SDCT_Parallel_Limited with " << effective_threads << " threads (max available: " << max_threads << ")" << std::endl;
    
    // Thread-local results
    std::vector<std::vector<std::vector<TreeGraphNode>>> thread_results(effective_threads);

    #pragma omp parallel num_threads(effective_threads)
    {
        int tid = omp_get_thread_num();
        thread_results[tid].reserve(size / effective_threads * 10);

        // Thread-local copies
        std::vector<int> local_vertexSets(size);
        std::vector<int> local_vertexLookup(size);
        std::vector<int *> local_neighborsInP(size);
        std::vector<int> local_numNeighbors(size, 1);

        for (int i = 0; i < size; ++i) {
            local_vertexLookup[i] = i;
            local_vertexSets[i] = i;
            local_neighborsInP[i] = static_cast<int *>(Calloc(1, sizeof(int)));
        }

        daf::StaticVector<int> dropV(MAX_CSIZE);
        daf::StaticVector<int> keepV(MAX_CSIZE);

        int beginX = 0;
        int beginP = 0;
        int beginR = size;

        // Use guided scheduling for better load balancing
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
    for (int tid = 0; tid < effective_threads; ++tid) {
        total_cliques += thread_results[tid].size();
    }
    treeGraph.adj_list.reserve(total_cliques);

    for (int tid = 0; tid < effective_threads; ++tid) {
        for (auto &clique : thread_results[tid]) {
            treeGraph.adj_list.push_back(std::move(clique));
        }
    }

    return treeGraph;
}
