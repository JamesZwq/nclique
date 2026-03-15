/*
 * SDCT_Parallel_Batch: Batch-based parallel SDCT
 * Key insight: Process vertices in large batches to reduce thread overhead
 * Each batch has its own data, but batches are processed in parallel
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

DynamicGraph<TreeGraphNode> SDCT_Parallel_Batch(Graph &edgeGraph, int max_k, int min_k) {
    auto size = edgeGraph.getGraphNodeSize();
    int nthreads = omp_get_max_threads();
    std::cout << "SDCT_Parallel_Batch with " << nthreads << " threads" << std::endl;
    
    // Determine batch size: process vertices in batches to reduce overhead
    int batch_size = std::max(1, (int)size / (nthreads * 4));  // 4 batches per thread
    int num_batches = ((int)size + batch_size - 1) / batch_size;
    
    std::cout << "Batch size: " << batch_size << ", Num batches: " << num_batches << std::endl;
    
    // Results for each batch
    std::vector<std::vector<std::vector<TreeGraphNode>>> batch_results(num_batches);

    #pragma omp parallel for schedule(dynamic, 1) collapse(1)
    for (int batch_id = 0; batch_id < num_batches; ++batch_id) {
        int start_vertex = batch_id * batch_size;
        int end_vertex = std::min(start_vertex + batch_size, (int)size);
        
        // Each batch has its own data structures
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

        daf::StaticVector<int> dropV(MAX_CSIZE);
        daf::StaticVector<int> keepV(MAX_CSIZE);
        batch_results[batch_id].reserve(batch_size * 10);

        int beginX = 0;
        int beginP = 0;
        int beginR = size;

        // Process vertices in this batch
        for (int vertex = start_vertex; vertex < end_vertex; ++vertex) {
            int newBeginX, newBeginP, newBeginR;

            fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(
                vertex,
                vertexSets.data(), vertexLookup.data(),
                edgeGraph,
                neighborsInP.data(), numNeighbors.data(),
                &beginX, &beginP, &beginR,
                &newBeginX, &newBeginP, &newBeginR);

            dropV.clear();
            keepV.clear();
            keepV.push_back(vertex);

            listAllCliquesDegeneracyRecursive_VedgeGraphSDCT(
                vertexSets.data(),
                vertexLookup.data(), neighborsInP.data(),
                numNeighbors.data(), newBeginP,
                newBeginR, keepV, dropV, max_k, min_k, batch_results[batch_id]);

            beginR = beginR + 1;
        }

        // Cleanup
        for (int i = 0; i < size; ++i) {
            Free(neighborsInP[i]);
        }
        vertexSets.free();
        vertexLookup.free();
        neighborsInP.free();
        numNeighbors.free();
    }

    // Merge results from all batches
    DynamicGraph<TreeGraphNode> treeGraph(size);
    size_t total_cliques = 0;
    for (int bid = 0; bid < num_batches; ++bid) {
        total_cliques += batch_results[bid].size();
    }
    treeGraph.adj_list.reserve(total_cliques);

    for (int bid = 0; bid < num_batches; ++bid) {
        for (auto &clique : batch_results[bid]) {
            treeGraph.adj_list.push_back(std::move(clique));
        }
    }

    return treeGraph;
}
