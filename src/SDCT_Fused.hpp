// SDCT_Fused: builds tree + invokes onLeaf callback for each leaf during BK enumeration.
//
// Key feature: separate emit_min_k and store_min_k:
//   - emit_min_k (=r): callback fires for ALL leaves with size ≥ r
//   - store_min_k (=s): only leaves with size ≥ s are stored in the tree
//
// Callback signature:
//   void onLeaf(daf::Size leafId, const std::vector<TreeGraphNode> &leaf, bool stored)
//   - leafId: valid tree node ID if stored=true; undefined if stored=false
//   - leaf: the vertex list (keep + pivot)
//   - stored: whether this leaf was added to the tree
//
// This allows R≥3 to build r-clique indices from ALL leaves while only storing
// s-clique-relevant leaves in the tree.

#pragma once
// Included after all other headers in degeneracy_cliques.cpp

template<typename OnLeafFunc>
static void listAllCliquesDegeneracyRecursive_Fused(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR, daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int emit_min_k, int store_min_k,
    DynamicGraph<TreeGraphNode> &tree, OnLeafFunc &onLeaf) {

    if ((int)keepV.size() > max_k) return;

    if (beginP >= beginR) {
        int cSize = (int)keepV.size() + (int)dropV.size();
        if (cSize < emit_min_k) return;

        std::vector<TreeGraphNode> newNode;
        newNode.reserve(cSize);
        for (uint64_t i: keepV) newNode.emplace_back(i, false);
        if ((int)keepV.size() < max_k) {
            for (uint64_t i: dropV) newNode.emplace_back(i, true);
            std::ranges::sort(newNode);
        }

        if (cSize >= store_min_k) {
            // Store in tree + callback with stored=true
            auto leafId = tree.addNode(newNode);
            onLeaf(leafId, tree.adj_list[leafId], true);
        } else {
            // Don't store in tree — callback with stored=false
            onLeaf(daf::Size(-1), newNode, false);
        }
        return;
    }

    int *myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough;

    int pivot = findBestPivotNonNeighborsDegeneracyCliques(&myCandidatesToIterateThrough,
                                                           &numCandidatesToIterateThrough,
                                                           vertexSets, vertexLookup,
                                                           neighborsInP, numNeighbors,
                                                           beginP, beginR);

    if (numCandidatesToIterateThrough != 0) {
        int iterator = 0;
        while (iterator < numCandidatesToIterateThrough) {
            int vertex = myCandidatesToIterateThrough[iterator];
            int newBeginP, newBeginR;

            moveToRDegeneracyCliques(vertex, vertexSets, vertexLookup,
                                     neighborsInP, numNeighbors,
                                     &beginP, &beginR, &newBeginP, &newBeginR);

            if (vertex == pivot) {
                dropV.push_back(vertex);
                listAllCliquesDegeneracyRecursive_Fused(vertexSets, vertexLookup,
                    neighborsInP, numNeighbors, newBeginP, newBeginR,
                    keepV, dropV, max_k, emit_min_k, store_min_k, tree, onLeaf);
                dropV.pop_back();
            } else {
                keepV.push_back(vertex);
                listAllCliquesDegeneracyRecursive_Fused(vertexSets, vertexLookup,
                    neighborsInP, numNeighbors, newBeginP, newBeginR,
                    keepV, dropV, max_k, emit_min_k, store_min_k, tree, onLeaf);
                keepV.pop_back();
            }
            iterator++;
        }
    }
    Free(myCandidatesToIterateThrough);
}

// SDCT_Fused: main entry point
// emit_min_k: minimum leaf size to invoke callback (typically r)
// store_min_k: minimum leaf size to store in tree (typically s)
template<typename OnLeafFunc>
DynamicGraph<TreeGraphNode> SDCT_Fused(Graph &edgeGraph, int max_k,
                                        int emit_min_k, int store_min_k,
                                        OnLeafFunc &&onLeaf) {
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

    int beginX = 0, beginP = 0, beginR = size;
    daf::StaticVector<int> dropV(MAX_CSIZE);
    daf::StaticVector<int> keepV(MAX_CSIZE);

    DynamicGraph<TreeGraphNode> treeGraph(edgeGraph.getGraphNodeSize());
    for (int vertex = 0; vertex < edgeGraph.getGraphNodeSize(); ++vertex) {
        int newBeginX, newBeginP, newBeginR;
        fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(vertex,
            vertexSets.data(), vertexLookup.data(), edgeGraph,
            neighborsInP.data(), numNeighbors.data(),
            &beginX, &beginP, &beginR, &newBeginX, &newBeginP, &newBeginR);

        dropV.clear();
        keepV.clear();
        keepV.push_back(vertex);
        listAllCliquesDegeneracyRecursive_Fused(vertexSets.data(),
            vertexLookup.data(), neighborsInP.data(), numNeighbors.data(),
            newBeginP, newBeginR, keepV, dropV, max_k,
            emit_min_k, store_min_k, treeGraph, onLeaf);

        beginR = beginR + 1;
    }

    vertexSets.free();
    vertexLookup.free();
    for (int i = 0; i < size; ++i) Free(neighborsInP[i]);
    neighborsInP.free();
    numNeighbors.free();

    return treeGraph;
}
