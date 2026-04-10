// SDCT_MaxClique: identical to SDCT_Fused but tracks X to tag maximal cliques.
// X does NOT affect path construction. It only determines isMaximalClique per path.

#pragma once

template<typename OnLeafFunc>
static void listAllCliquesDegeneracyRecursive_MaxClique(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR, daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int emit_min_k, int store_min_k,
    DynamicGraph<TreeGraphNode> &tree, OnLeafFunc &onLeaf,
    Graph &edgeGraph, std::vector<int> &xVerts) {

    if ((int)keepV.size() > max_k) return;

    if (beginP >= beginR) {
        // P is empty. X empty → maximal clique.
        bool isMaximal = xVerts.empty();
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
            auto leafId = tree.addNode(newNode);
            onLeaf(leafId, tree.adj_list[leafId], true, isMaximal);
        } else {
            onLeaf(daf::Size(-1), newNode, false, isMaximal);
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

            // Compute newX = {u ∈ xVerts : u is neighbor of vertex}
            std::vector<int> newX;
            for (int u : xVerts) {
                bool adj = false;
                // Check u → vertex
                auto [uB, uE] = edgeGraph.getNbr(u);
                for (int ei = uB; ei < uE; ++ei) {
                    if ((int)edgeGraph.adj_list[ei] == vertex) { adj = true; break; }
                }
                // Check vertex → u (bidirectional)
                if (!adj) {
                    auto [vB, vE] = edgeGraph.getNbr(vertex);
                    for (int ei = vB; ei < vE; ++ei) {
                        if ((int)edgeGraph.adj_list[ei] == u) { adj = true; break; }
                    }
                }
                if (adj) newX.push_back(u);
            }

            if (vertex == pivot) {
                dropV.push_back(vertex);
                listAllCliquesDegeneracyRecursive_MaxClique(vertexSets, vertexLookup,
                    neighborsInP, numNeighbors, newBeginP, newBeginR,
                    keepV, dropV, max_k, emit_min_k, store_min_k, tree, onLeaf,
                    edgeGraph, newX);
                dropV.pop_back();
            } else {
                keepV.push_back(vertex);
                listAllCliquesDegeneracyRecursive_MaxClique(vertexSets, vertexLookup,
                    neighborsInP, numNeighbors, newBeginP, newBeginR,
                    keepV, dropV, max_k, emit_min_k, store_min_k, tree, onLeaf,
                    edgeGraph, newX);
                keepV.pop_back();
            }

            // After recursion: vertex has been processed → add to X
            xVerts.push_back(vertex);

            iterator++;
        }

        // Restore X: remove the vertices we added in this call
        for (int i = 0; i < numCandidatesToIterateThrough; ++i)
            xVerts.pop_back();
    }
    Free(myCandidatesToIterateThrough);
}

template<typename OnLeafFunc>
DynamicGraph<TreeGraphNode> SDCT_MaxClique(Graph &edgeGraph, int max_k,
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

        // Build initial X = earlier neighbors of vertex (neighbor < vertex in degeneracy order)
        // fillInPandX doesn't fill X, so we build it manually from the graph.
        std::vector<int> xVerts;
        {
            auto [nbBegin, nbEnd] = edgeGraph.getNbr(vertex);
            for (int idx = nbBegin; idx < nbEnd; ++idx) {
                int nb = edgeGraph.adj_list[idx];
                if (nb < vertex) xVerts.push_back(nb); // earlier neighbor
            }
        }

        dropV.clear();
        keepV.clear();
        keepV.push_back(vertex);
        listAllCliquesDegeneracyRecursive_MaxClique(vertexSets.data(),
            vertexLookup.data(), neighborsInP.data(), numNeighbors.data(),
            newBeginP, newBeginR, keepV, dropV, max_k,
            emit_min_k, store_min_k, treeGraph, onLeaf,
            edgeGraph, xVerts);

        beginR = beginR + 1;
    }

    vertexSets.free();
    vertexLookup.free();
    for (int i = 0; i < size; ++i) Free(neighborsInP[i]);
    neighborsInP.free();
    numNeighbors.free();

    return treeGraph;
}
