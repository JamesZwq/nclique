/*
 * SDCT_Parallel_Recursive: Task-based parallel SDCT
 * Key insight: Parallelize the recursion tree, not the vertex loop
 * This allows sharing of data structures and automatic load balancing
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

// Forward declarations
void listAllCliquesDegeneracyRecursive_VedgeGraphSDCT(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR, daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int min_k, std::vector<std::vector<TreeGraphNode>> &tree_buffer);

// Parallel recursive function
void listAllCliquesDegeneracyRecursive_Parallel(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR, daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int min_k, std::vector<std::vector<TreeGraphNode>> &tree_buffer,
    int depth) {
    
    if (keepV.size() > max_k) return;
    
    if (beginP >= beginR) {
        int cSize = keepV.size() + dropV.size();
        if (cSize < min_k) return;
        
        std::vector<TreeGraphNode> node;
        node.reserve(cSize);
        
        if (keepV.size() == max_k) {
            for (int i = 0; i < keepV.size(); i++)
                node.emplace_back(keepV[i], false);
            tree_buffer.push_back(std::move(node));
            return;
        }
        
        for (int i = 0; i < keepV.size(); i++)
            node.emplace_back(keepV[i], false);
        for (int i = 0; i < dropV.size(); i++)
            node.emplace_back(dropV[i], true);
        
        if (cSize <= 32)
            std::sort(node.begin(), node.end());
        else
            std::sort(node.begin(), node.end());
        
        tree_buffer.push_back(std::move(node));
        return;
    }
    
    // Find pivot
    int pivot = -1, maxInP = -1, sz = beginR - beginP;
    for (int j = beginP; j < beginR; j++) {
        int v = vertexSets[j];
        int numP = MY_MIN(sz, numNeighbors[v]);
        int inP = 0;
        for (int k = 0; k < numP; k++) {
            int l = vertexLookup[neighborsInP[v][k]];
            if (l >= beginP && l < beginR) inP++;
            else break;
        }
        if (inP > maxInP) {
            maxInP = inP;
            pivot = v;
        }
    }
    
    // Mark pivot non-neighbors
    int nPN = MY_MIN(sz, numNeighbors[pivot]);
    std::vector<bool> marked(beginR, false);
    for (int j = 0; j < nPN; j++) {
        int nb = neighborsInP[pivot][j];
        int l = vertexLookup[nb];
        if (l >= beginP && l < beginR) marked[l] = true;
        else break;
    }
    
    // Get candidates
    std::vector<int> cands;
    for (int j = beginP; j < beginR; j++) {
        if (!marked[j]) cands.push_back(vertexSets[j]);
    }
    
    // Process candidates in parallel (only at shallow depths)
    if (depth < 3 && cands.size() > 1) {
        #pragma omp taskgroup
        {
            for (int ci = 0; ci < cands.size(); ci++) {
                #pragma omp task
                {
                    int vertex = cands[ci];
                    int vl = vertexLookup[vertex];
                    
                    // Save state
                    int saved_beginR = beginR;
                    int saved_beginP = beginP;
                    
                    // Move to R
                    beginR--;
                    vertexSets[vl] = vertexSets[beginR];
                    vertexLookup[vertexSets[beginR]] = vl;
                    vertexSets[beginR] = vertex;
                    vertexLookup[vertex] = beginR;
                    
                    int newBeginP = beginP;
                    int newBeginR = beginP;
                    int sizeOfP = beginR - beginP;
                    
                    // Mark neighbors
                    std::vector<bool> neighbor_marked(beginR, false);
                    int numNbr = MY_MIN(sizeOfP, numNeighbors[vertex]);
                    for (int k = 0; k < numNbr; k++)
                        neighbor_marked[vertexLookup[neighborsInP[vertex][k]]] = true;
                    
                    // Partition P
                    for (int j = beginP; j < beginR; j++) {
                        int nb = vertexSets[j];
                        if (neighbor_marked[j]) {
                            vertexSets[j] = vertexSets[newBeginR];
                            vertexLookup[vertexSets[newBeginR]] = j;
                            vertexSets[newBeginR] = nb;
                            vertexLookup[nb] = newBeginR;
                            newBeginR++;
                        }
                    }
                    
                    // Recurse
                    if (vertex == pivot) {
                        dropV.push_back(vertex);
                        listAllCliquesDegeneracyRecursive_Parallel(
                            vertexSets, vertexLookup, neighborsInP, numNeighbors,
                            newBeginP, newBeginR, keepV, dropV, max_k, min_k, tree_buffer, depth + 1);
                        dropV.pop_back();
                    } else {
                        keepV.push_back(vertex);
                        listAllCliquesDegeneracyRecursive_Parallel(
                            vertexSets, vertexLookup, neighborsInP, numNeighbors,
                            newBeginP, newBeginR, keepV, dropV, max_k, min_k, tree_buffer, depth + 1);
                        keepV.pop_back();
                    }
                    
                    // Restore state
                    beginR = saved_beginR;
                    beginP = saved_beginP;
                }
            }
        }
    } else {
        // Sequential processing for deep recursion
        for (int ci = 0; ci < cands.size(); ci++) {
            int vertex = cands[ci];
            int vl = vertexLookup[vertex];
            
            beginR--;
            vertexSets[vl] = vertexSets[beginR];
            vertexLookup[vertexSets[beginR]] = vl;
            vertexSets[beginR] = vertex;
            vertexLookup[vertex] = beginR;
            
            int newBeginP = beginP;
            int newBeginR = beginP;
            int sizeOfP = beginR - beginP;
            
            std::vector<bool> neighbor_marked(beginR, false);
            int numNbr = MY_MIN(sizeOfP, numNeighbors[vertex]);
            for (int k = 0; k < numNbr; k++)
                neighbor_marked[vertexLookup[neighborsInP[vertex][k]]] = true;
            
            for (int j = beginP; j < beginR; j++) {
                int nb = vertexSets[j];
                if (neighbor_marked[j]) {
                    vertexSets[j] = vertexSets[newBeginR];
                    vertexLookup[vertexSets[newBeginR]] = j;
                    vertexSets[newBeginR] = nb;
                    vertexLookup[nb] = newBeginR;
                    newBeginR++;
                }
            }
            
            if (vertex == pivot) {
                dropV.push_back(vertex);
                listAllCliquesDegeneracyRecursive_Parallel(
                    vertexSets, vertexLookup, neighborsInP, numNeighbors,
                    newBeginP, newBeginR, keepV, dropV, max_k, min_k, tree_buffer, depth + 1);
                dropV.pop_back();
            } else {
                keepV.push_back(vertex);
                listAllCliquesDegeneracyRecursive_Parallel(
                    vertexSets, vertexLookup, neighborsInP, numNeighbors,
                    newBeginP, newBeginR, keepV, dropV, max_k, min_k, tree_buffer, depth + 1);
                keepV.pop_back();
            }
        }
    }
}

DynamicGraph<TreeGraphNode> SDCT_Parallel_Recursive(Graph &edgeGraph, int max_k, int min_k) {
    auto size = edgeGraph.getGraphNodeSize();
    std::cout << "SDCT_Parallel_Recursive" << std::endl;
    
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

    std::vector<std::vector<TreeGraphNode>> treeBuffer;
    daf::StaticVector<int> dropV(MAX_CSIZE);
    daf::StaticVector<int> keepV(MAX_CSIZE);

    int beginX = 0;
    int beginP = 0;
    int beginR = size;

    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int vertex = 0; vertex < size; ++vertex) {
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

                #pragma omp task
                {
                    listAllCliquesDegeneracyRecursive_Parallel(
                        vertexSets.data(),
                        vertexLookup.data(), neighborsInP.data(),
                        numNeighbors.data(), newBeginP,
                        newBeginR, keepV, dropV, max_k, min_k, treeBuffer, 0);
                }

                beginR = beginR + 1;
            }
        }
    }

    DynamicGraph<TreeGraphNode> treeGraph(size);
    treeGraph.adj_list.reserve(treeBuffer.size());
    for (auto &clique : treeBuffer) {
        treeGraph.adj_list.push_back(std::move(clique));
    }

    vertexSets.free();
    vertexLookup.free();
    for (int i = 0; i < size; ++i) {
        Free(neighborsInP[i]);
    }
    neighborsInP.free();
    numNeighbors.free();

    return treeGraph;
}
