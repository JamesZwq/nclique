//
// SDCT_Augmented.inl — template implementation for callback-based SDCT.
//
// Included at the bottom of SDCT_Augmented.h. Do NOT compile separately.
//
#pragma once

#include "degeneracy_algorithm_cliques_V.h"
#include "misc.h"
#include "LinkedList.h"
#include "MemoryManager.h"
#include <vector>
#include <algorithm>

extern double nCr[1001][401];

// ---------------------------------------------------------------------------
//  Shared context struct to avoid passing many pointers through recursion.
//  The callback and emptyDrop live here, accessed via a single pointer.
// ---------------------------------------------------------------------------
struct AugCtx {
    using CallbackFn = void(*)(void *userData, daf::Size leafId,
        const daf::StaticVector<int> &keepV, const daf::StaticVector<int> &dropV);

    CallbackFn callback;
    void *userData;
    daf::StaticVector<int> emptyDrop{1}; // capacity 1 to avoid null data_

    void invoke(daf::Size leafId, const daf::StaticVector<int> &keepV,
                const daf::StaticVector<int> &dropV) const {
        callback(userData, leafId, keepV, dropV);
    }
};

// ---------------------------------------------------------------------------
//  Recursive BK — NO tree storage
// ---------------------------------------------------------------------------
static void bkRecurse_NoTree(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR,
    daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int min_k,
    size_t &leafCounter,
    AugCtx &ctx)
{
    if ((int)keepV.size() > max_k) return;

    if (beginP >= beginR) {
        int cSize = (int)keepV.size() + (int)dropV.size();
        if (cSize < min_k || cSize < max_k) return;

        daf::Size leafId = static_cast<daf::Size>(leafCounter++);
        if ((int)keepV.size() == max_k) {
            ctx.invoke(leafId, keepV, ctx.emptyDrop);
        } else {
            ctx.invoke(leafId, keepV, dropV);
        }
        return;
    }

    int *myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough;

    int pivot = findBestPivotNonNeighborsDegeneracyCliques(
        &myCandidatesToIterateThrough, &numCandidatesToIterateThrough,
        vertexSets, vertexLookup, neighborsInP, numNeighbors,
        beginP, beginR);

    if (numCandidatesToIterateThrough != 0) {
        int iterator = 0;
        while (iterator < numCandidatesToIterateThrough) {
            int vertex = myCandidatesToIterateThrough[iterator];
            int newBeginP, newBeginR;

            moveToRDegeneracyCliques(vertex,
                vertexSets, vertexLookup, neighborsInP, numNeighbors,
                &beginP, &beginR, &newBeginP, &newBeginR);

            if (vertex == pivot) {
                dropV.push_back(vertex);
                bkRecurse_NoTree(
                    vertexSets, vertexLookup, neighborsInP, numNeighbors,
                    newBeginP, newBeginR, keepV, dropV, max_k, min_k,
                    leafCounter, ctx);
                dropV.pop_back();
            } else {
                keepV.push_back(vertex);
                bkRecurse_NoTree(
                    vertexSets, vertexLookup, neighborsInP, numNeighbors,
                    newBeginP, newBeginR, keepV, dropV, max_k, min_k,
                    leafCounter, ctx);
                keepV.pop_back();
            }
            iterator++;
        }
    }

    Free(myCandidatesToIterateThrough);
}

// ---------------------------------------------------------------------------
//  Recursive BK — with tree storage
// ---------------------------------------------------------------------------
static void bkRecurse_WithTree(
    int *vertexSets, int *vertexLookup,
    int **neighborsInP, int *numNeighbors,
    int beginP, int beginR,
    daf::StaticVector<int> &keepV, daf::StaticVector<int> &dropV,
    int max_k, int min_k,
    DynamicGraph<TreeGraphNode> &tree,
    AugCtx &ctx)
{
    if ((int)keepV.size() > max_k) return;

    if (beginP >= beginR) {
        int cSize = (int)keepV.size() + (int)dropV.size();
        if (cSize < min_k) return;

        std::vector<TreeGraphNode> newNode;
        if ((int)keepV.size() == max_k) {
            newNode.reserve(keepV.size());
            for (uint64_t i : keepV) newNode.emplace_back(i, false);
            daf::Size leafId = tree.addNode(newNode);
            ctx.invoke(leafId, keepV, ctx.emptyDrop);
            return;
        }
        newNode.reserve(cSize);
        for (uint64_t i : keepV) newNode.emplace_back(i, false);
        for (uint64_t i : dropV) newNode.emplace_back(i, true);
        std::ranges::sort(newNode);
        daf::Size leafId = tree.addNode(newNode);
        ctx.invoke(leafId, keepV, dropV);
        return;
    }

    int *myCandidatesToIterateThrough;
    int numCandidatesToIterateThrough;

    int pivot = findBestPivotNonNeighborsDegeneracyCliques(
        &myCandidatesToIterateThrough, &numCandidatesToIterateThrough,
        vertexSets, vertexLookup, neighborsInP, numNeighbors,
        beginP, beginR);

    if (numCandidatesToIterateThrough != 0) {
        int iterator = 0;
        while (iterator < numCandidatesToIterateThrough) {
            int vertex = myCandidatesToIterateThrough[iterator];
            int newBeginP, newBeginR;

            moveToRDegeneracyCliques(vertex,
                vertexSets, vertexLookup, neighborsInP, numNeighbors,
                &beginP, &beginR, &newBeginP, &newBeginR);

            if (vertex == pivot) {
                dropV.push_back(vertex);
                bkRecurse_WithTree(
                    vertexSets, vertexLookup, neighborsInP, numNeighbors,
                    newBeginP, newBeginR, keepV, dropV, max_k, min_k,
                    tree, ctx);
                dropV.pop_back();
            } else {
                keepV.push_back(vertex);
                bkRecurse_WithTree(
                    vertexSets, vertexLookup, neighborsInP, numNeighbors,
                    newBeginP, newBeginR, keepV, dropV, max_k, min_k,
                    tree, ctx);
                keepV.pop_back();
            }
            iterator++;
        }
    }

    Free(myCandidatesToIterateThrough);
}

// ---------------------------------------------------------------------------
//  Top-level: SDCT_Augmented (with tree)
// ---------------------------------------------------------------------------
template<typename OnLeafFn>
DynamicGraph<TreeGraphNode> SDCT_Augmented(
    Graph &edgeGraph, int max_k, int min_k, OnLeafFn &&onLeaf)
{
    // Wrap the lambda into AugCtx via a type-erased trampoline
    auto trampoline = [](void *ud, daf::Size leafId,
        const daf::StaticVector<int> &keepV, const daf::StaticVector<int> &dropV) {
        (*static_cast<OnLeafFn*>(ud))(leafId, keepV, dropV);
    };
    AugCtx ctx;
    ctx.callback = trampoline;
    ctx.userData = &onLeaf;

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
        int newBeginX, newBeginP, newBeginR;

        fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(vertex,
            vertexSets.data(), vertexLookup.data(), edgeGraph,
            neighborsInP.data(), numNeighbors.data(),
            &beginX, &beginP, &beginR,
            &newBeginX, &newBeginP, &newBeginR);

        dropV.clear();
        keepV.clear();
        keepV.push_back(vertex);
        bkRecurse_WithTree(
            vertexSets.data(), vertexLookup.data(),
            neighborsInP.data(), numNeighbors.data(),
            newBeginP, newBeginR, keepV, dropV, max_k, min_k,
            treeGraph, ctx);

        beginR = beginR + 1;
    }

    vertexSets.free();
    vertexLookup.free();
    for (int i = 0; i < size; ++i) Free(neighborsInP[i]);
    neighborsInP.free();
    numNeighbors.free();

    return treeGraph;
}

// ---------------------------------------------------------------------------
//  Top-level: SDCT_Augmented_NoTree (tree-free)
// ---------------------------------------------------------------------------
template<typename OnLeafFn>
size_t SDCT_Augmented_NoTree(
    Graph &edgeGraph, int max_k, int min_k, OnLeafFn &&onLeaf)
{
    // Wrap the lambda into AugCtx via a type-erased trampoline
    auto trampoline = [](void *ud, daf::Size leafId,
        const daf::StaticVector<int> &keepV, const daf::StaticVector<int> &dropV) {
        (*static_cast<OnLeafFn*>(ud))(leafId, keepV, dropV);
    };
    AugCtx ctx;
    ctx.callback = trampoline;
    ctx.userData = &onLeaf;

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

    size_t leafCounter = 0;
    for (int vertex = 0; vertex < edgeGraph.getGraphNodeSize(); ++vertex) {
        int newBeginX, newBeginP, newBeginR;

        fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(vertex,
            vertexSets.data(), vertexLookup.data(), edgeGraph,
            neighborsInP.data(), numNeighbors.data(),
            &beginX, &beginP, &beginR,
            &newBeginX, &newBeginP, &newBeginR);

        dropV.clear();
        keepV.clear();
        keepV.push_back(vertex);
        bkRecurse_NoTree(
            vertexSets.data(), vertexLookup.data(),
            neighborsInP.data(), numNeighbors.data(),
            newBeginP, newBeginR, keepV, dropV, max_k, min_k,
            leafCounter, ctx);

        beginR = beginR + 1;
    }

    vertexSets.free();
    vertexLookup.free();
    for (int i = 0; i < size; ++i) Free(neighborsInP[i]);
    neighborsInP.free();
    numNeighbors.free();

    return leafCounter;
}

// ---------------------------------------------------------------------------
//  Top-level: SDCT_Augmented_NoTree_Interleaved (tree-free + dual callback)
//  onLeaf(leafId, keepV, dropV) — called at each leaf emission
//  onVertexDone(vertex)         — called after each top-level vertex's subtree
// ---------------------------------------------------------------------------
template<typename OnLeafFn, typename OnVertexDoneFn>
size_t SDCT_Augmented_NoTree_Interleaved(
    Graph &edgeGraph, int max_k, int min_k,
    OnLeafFn &&onLeaf, OnVertexDoneFn &&onVertexDone)
{
    // Wrap the onLeaf lambda into AugCtx via a type-erased trampoline
    auto trampoline = [](void *ud, daf::Size leafId,
        const daf::StaticVector<int> &keepV, const daf::StaticVector<int> &dropV) {
        (*static_cast<OnLeafFn*>(ud))(leafId, keepV, dropV);
    };
    AugCtx ctx;
    ctx.callback = trampoline;
    ctx.userData = &onLeaf;

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

    size_t leafCounter = 0;
    for (int vertex = 0; vertex < edgeGraph.getGraphNodeSize(); ++vertex) {
        int newBeginX, newBeginP, newBeginR;

        fillInPandXForRecursiveCallDegeneracyCliquesEdgeGraph(vertex,
            vertexSets.data(), vertexLookup.data(), edgeGraph,
            neighborsInP.data(), numNeighbors.data(),
            &beginX, &beginP, &beginR,
            &newBeginX, &newBeginP, &newBeginR);

        dropV.clear();
        keepV.clear();
        keepV.push_back(vertex);
        bkRecurse_NoTree(
            vertexSets.data(), vertexLookup.data(),
            neighborsInP.data(), numNeighbors.data(),
            newBeginP, newBeginR, keepV, dropV, max_k, min_k,
            leafCounter, ctx);

        beginR = beginR + 1;

        // Signal that vertex's subtree is complete — countingV[vertex] is finalized
        onVertexDone(vertex);
    }

    vertexSets.free();
    vertexLookup.free();
    for (int i = 0; i < size; ++i) Free(neighborsInP[i]);
    neighborsInP.free();
    numNeighbors.free();

    return leafCounter;
}
