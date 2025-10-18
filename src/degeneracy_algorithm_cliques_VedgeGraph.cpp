/*
    This file contains the algorithm for listing all cliques
    according to the algorithm of Jain et al. specified in
    "The power of pivoting for exact clique counting." (WSDM 2020).

    This code is a modified version of the code of quick-cliques-1.0 library for counting
    maximal cliques by Darren Strash (first name DOT last name AT gmail DOT com).

    Original author: Darren Strash (first name DOT last name AT gmail DOT com)

    Copyright (c) 2011 Darren Strash. This code is released under the GNU Public License (GPL) 3.0.

    Modifications Copyright (c) 2020 Shweta Jain

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#include<cassert>
#include<climits>
#include<cstdio>
#include<cstdlib>
#include<vector>
#include <iostream>
#include <fstream>

#include"degeneracy_algorithm_cliques_V.h"
#include"degeneracy_helper.h"
#include"LinkedList.h"
#include"MemoryManager.h"
#include"misc.h"
#include"tree/MultiBranchTree.h"
#include <cstring>
#include <NucleusDecomposition/NCliqueCoreDecomposition.h>
// #include"nCr.h"

extern double nCr[1001][401];

void listAllCliquesDegeneracyRecursive_VedgeGraph(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR, daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int min_k, DynamicGraph<TreeGraphNode> &tree) ;
/*! \brief Computes the vertex v in P union X that has the most neighbors in P,
           and places P \ {neighborhood of v} in an array. These are the 
           vertices to consider adding to the partial clique during the current
           recursive call of the algorithm.

    \param pivotNonNeighbors  An intially unallocated pointer, which will contain the set 
                              P \ {neighborhood of v} when this function completes.

    \param numNonNeighbors A pointer to a single integer, which has been preallocated,
                           which will contain the number of elements in pivotNonNeighbors.
 
    \param vertexSets An array containing sets of vertices divided into sets X, P, R, and other.
 
    \param vertexLookup A lookup table indexed by vertex number, storing the index of that 
                        vertex in vertexSets.

    \param neighborsInP Maps vertices to arrays of neighbors such that 
                        neighbors in P fill the first cells

    \param numNeighbors An the neighbor of neighbors a vertex had in P,
                        the first time this function is called, this bound is 
                        used to keep us from allocating more than linear space.
 
    \param beginX The index where set X begins in vertexSets.
 
    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.
    \param min_k

*/


/*!

    \param adjList An array of linked lists, representing the input graph in the
                   "typical" adjacency list format.
 
    \param adjacencyList an array of arrays, representing the input graph in a more
                         compact and cache-friendly adjacency list format. (not currently used)

    \param cliques A linked list of cliques to return. <b>(only available when compiled 
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \param degree An array, indexed by vertex, containing the degree of that vertex. (not currently used)

    \param size The number of vertices in the graph.

    \return the number of maximal cliques of the input graph.
*/

DynamicGraph<TreeGraphNode> listAllCliquesDegeneracy_VedgeGraph(Graph &edgeGraph, int max_k, int min_k) {
    // vertex sets are stored in an array like this:
    // |--X--|--P--|
    auto size = edgeGraph.getGraphNodeSize();
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

    int beginX = 0;
    int beginP = 0;
    int beginR = size;

    daf::StaticVector<int> dropV(MAX_CSIZE);
    daf::StaticVector<int> keepV(MAX_CSIZE);

    DynamicGraph<TreeGraphNode> treeGraph(edgeGraph.getGraphNodeSize());
    for (int vertex = 0; vertex < edgeGraph.getGraphNodeSize(); ++vertex) {
        // if (edgeGraph.coreV[vertex] < min_k) {
        //     continue; // skip vertices with core less than min_k
        // }
        int newBeginX, newBeginP, newBeginR;

        // of vertex
        fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(vertex,
                                                     vertexSets.data(), vertexLookup.data(),
                                                     edgeGraph,
                                                     neighborsInP.data(), numNeighbors.data(),
                                                     &beginX, &beginP, &beginR,
                                                     &newBeginX, &newBeginP, &newBeginR);


        dropV.clear();
        keepV.clear();
        keepV.push_back(vertex);
        // keepV.push_back(vertex);
        listAllCliquesDegeneracyRecursive_VedgeGraph(vertexSets.data(),
                                                     vertexLookup.data(), neighborsInP.data(),
                                                     numNeighbors.data(), newBeginP,
                                                     newBeginR, keepV, dropV, max_k, min_k, treeGraph);

        beginR = beginR + 1;
    }
    // tree.printTree();
    // tree.initLeafsParent();
    // tree.cliqueCount().print();

    // auto file1 = fopen("~/_/pivoter/outB.txt", "w");
    // for (auto &i: treeGraph.cliqueCount()) {
    //     printf("%lf\n", i);
    //     fprintf(file1, "%lf\n", i);
    // }
    // fclose(file1);
    vertexSets.free();
    vertexLookup.free();
    for (int i = 0; i < size; ++i) {
        Free(neighborsInP[i]);
    }
    neighborsInP.free();
    numNeighbors.free();


    return treeGraph;
}

/*! \brief Recursively list all maximal cliques containing all of
           all vertices in R, some vertices in P and no vertices in X.

    \param cliqueCount A pointer to the number of maximal cliques computed
                       thus far.

    \param cliques A linked list of cliques to return. <b>(only available when compiled
                   with RETURN_CLIQUES_ONE_BY_ONE defined)</b>

    \param partialClique A linked list storing R, the partial clique for this
                         recursive call.

    \param vertexSets An array containing sets of vertices divided into sets X, P, R and other.

    \param vertexLookup A lookup table indexed by vertex number, storing the index of that
                        vertex in vertexSets.

    \param neighborsInP Maps vertices to arrays of neighbors such that
                        neighbors in P fill the first cells

    \param numNeighbors An the neighbor of neighbors a vertex had in P,
                        the first time this function is called, this bound is
                        used to keep us from allocating more than linear space.


    \param beginP The index where set P begins in vertexSets.

    \param beginR The index where set R begins in vertexSets.

*/

void listAllCliquesDegeneracyRecursive_VedgeGraph(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR, daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int min_k, DynamicGraph<TreeGraphNode> &tree) {
    // std::cout << "max_k: " << max_k << std::endl;
    // if (keep > 3) {
    //     std::cerr << "keep should be less than 4" << std::endl;
    // }
    if (keepV.size() > max_k) {
        return;
    }
    if ((beginP >= beginR)) {
        auto cSize = keepV.size() + dropV.size();
        std::vector<TreeGraphNode> newNode;
        if (keepV.size() == max_k) {
            newNode.reserve(cSize);
            for (uint64_t i: keepV) {
                newNode.emplace_back(i, false);
            }
            tree.addNode(newNode);
            return;
        }
        newNode.reserve(cSize);
        for (uint64_t i: keepV) {
            newNode.emplace_back(i, false);
        }
        for (uint64_t i: dropV) {
            newNode.emplace_back(i, true);
        }
        std::ranges::sort(newNode);
        tree.addNode(newNode);
        return;
    }

    int *myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough;

    // get the candidates to add to R to make a maximal clique
    int pivot = findBestPivotNonNeighborsDegeneracyCliques(&myCandidatesToIterateThrough,
                                                           &numCandidatesToIterateThrough,
                                                           vertexSets, vertexLookup,
                                                           neighborsInP, numNeighbors,
                                                           beginP, beginR);


    // add candiate vertices to the partial clique one at a time and
    // search for maximal cliquesRazer Viper V3 Pro
    if (numCandidatesToIterateThrough != 0) {
        int iterator = 0;
        while (iterator < numCandidatesToIterateThrough) {
            // vertex to be added to the partial clique
            int vertex = myCandidatesToIterateThrough[iterator];

            int newBeginP, newBeginR;

            // swap vertex into R and update all data structures
            moveToRDegeneracyCliques(vertex,
                                     vertexSets, vertexLookup,
                                     neighborsInP, numNeighbors,
                                     &beginP, &beginR, &newBeginP,
                                     &newBeginR);


            // recursively compute maximal cliques with new sets R, P and X
            if (vertex == pivot) {
                // dropV[drop] = vertex;
                dropV.push_back(vertex);
                listAllCliquesDegeneracyRecursive_VedgeGraph(vertexSets,
                                                             vertexLookup, neighborsInP,
                                                             numNeighbors, newBeginP,
                                                             newBeginR, keepV, dropV, max_k, min_k, tree);
                dropV.pop_back();
            } else {
                keepV.push_back(vertex);
                listAllCliquesDegeneracyRecursive_VedgeGraph(vertexSets,
                                                             vertexLookup, neighborsInP,
                                                             numNeighbors, newBeginP,
                                                             newBeginR, keepV, dropV, max_k, min_k, tree);
                keepV.pop_back();
            }

            moveFromRToXDegeneracyCliques(vertex,
                                          vertexSets, vertexLookup,
                                          &beginP, &beginR);

            iterator++;
        }

        // swap vertices that were moved to X back into P, for higher recursive calls.
        iterator = 0;
        while (iterator < numCandidatesToIterateThrough) {
            int vertex = myCandidatesToIterateThrough[iterator];
            int vertexLocation = vertexLookup[vertex];

            beginP--;
            vertexSets[vertexLocation] = vertexSets[beginP];
            vertexSets[beginP] = vertex;
            vertexLookup[vertex] = beginP;
            vertexLookup[vertexSets[vertexLocation]] = vertexLocation;

            iterator++;
        }
    }


    // don't need to check for emptiness before freeing, since
    // something will always be there (we allocated enough memory
    // for all of P, which is nonempty)
    Free(myCandidatesToIterateThrough);

    return;
}
