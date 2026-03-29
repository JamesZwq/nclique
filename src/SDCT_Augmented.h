//
// SDCT_Augmented: Callback-based SDCT construction.
//
// At each leaf emission, invokes onLeaf(leafId, keepV, dropV) so that
// downstream index-building (support counting, CSR construction, etc.)
// can be fused into the single SDCT traversal — eliminating redundant
// post-hoc passes over the tree.
//
// Two entry points:
//   SDCT_Augmented()        — builds DynamicGraph<TreeGraphNode> AND calls onLeaf
//   SDCT_Augmented_NoTree() — calls onLeaf only, never stores tree (for R=1 tree-free peeling)
//
#pragma once
#ifndef SDCT_AUGMENTED_H
#define SDCT_AUGMENTED_H

#include "graph/Graph.h"
#include "graph/DynamicGraph.h"
#include "Global/Global.h"
#include <span>
#include <cstddef>

// ---------------------------------------------------------------------------
//  SDCT with tree storage + callback
// ---------------------------------------------------------------------------
//  OnLeafFn signature:
//    void(daf::Size leafId,
//         const daf::StaticVector<int>& keepV,
//         const daf::StaticVector<int>& dropV)
//
template<typename OnLeafFn>
DynamicGraph<TreeGraphNode> SDCT_Augmented(
    Graph &edgeGraph, int max_k, int min_k, OnLeafFn &&onLeaf);

// ---------------------------------------------------------------------------
//  SDCT without tree storage (tree-free) + callback
// ---------------------------------------------------------------------------
//  Returns the total number of leaves emitted.
//
template<typename OnLeafFn>
size_t SDCT_Augmented_NoTree(
    Graph &edgeGraph, int max_k, int min_k, OnLeafFn &&onLeaf);

// ---------------------------------------------------------------------------
//  SDCT without tree storage (tree-free) + dual callback
//  onLeaf called at each leaf; onVertexDone called after each top-level
//  vertex's subtree completes (vertex is then finalized).
//  Returns the total number of leaves emitted.
// ---------------------------------------------------------------------------
template<typename OnLeafFn, typename OnVertexDoneFn>
size_t SDCT_Augmented_NoTree_Interleaved(
    Graph &edgeGraph, int max_k, int min_k,
    OnLeafFn &&onLeaf, OnVertexDoneFn &&onVertexDone);

// Implementation is in SDCT_Augmented.cpp (template — explicit instantiation
// is not needed because callers include this header and the .cpp is compiled
// as part of the translation unit via the GLOB_RECURSE in CMakeLists.txt).
// We keep the implementation in a separate .inl file included at the bottom.
#include "SDCT_Augmented.inl"

#endif // SDCT_AUGMENTED_H
