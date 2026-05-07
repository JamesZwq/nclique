// (1,s)-Nucleus hierarchy via CPI-leaf-as-DSU-connector.  Strictly
// matches paper Def 4 (s-connectivity through admissible s-cliques)
// and paper §5 "hierarchy as peel by-product".
//
// Algorithm (post-peel, single pass):
//   1.  Reconstruct leaf-first vertex lists from V3's inverse CSR
//       (vtxLeafOff / vtxLeafIds).  Compute each leaf L's birth level
//       kP_L = min_{v ∈ L} core[v].  O(Σ) time, O(Σ) temp memory.
//   2.  Sort leaves by descending kP.
//   3.  DSU over (vertices ∪ leaves).  For each leaf L in sorted
//       order, union L with every v ∈ L.  This makes L a "connector
//       node": two vertices are in the same DSU component iff they
//       are s-connected via admissible CPI cliques (paper Def 4),
//       provably equivalent because every admissible s-clique of L
//       at level kP_L witnesses pairwise membership through V_h.
//   4.  Each merge of two non-trivial components fires an elder-rule
//       event.  Output the resulting forest as CSV.
//
// Memory: V3 peel itself stays at 10Σ — this routine runs *after*
// peel, takes O(Σ) temporary memory for the leaf-first scratch and
// frees it on return.

#pragma once

#include <cstddef>
#include <vector>

#include "graph/DynamicGraph.h"
#include "Global/Global.h"

struct ST_V2_Data;  // forward

namespace nucleus_hier {

// Build the hierarchy directly from V3's dual CSR (no leaf-first tree
// required).  `coreV[v]` is κ_s(v) (-1 if v not in any nucleus).
// `min_size` filters trivial branches; pass 1 to keep all.
void buildAndDumpHierarchyFromCSR(
    const ST_V2_Data& v3data,
    const double*     coreV,
    daf::Size         numVertices,
    daf::Size         s,
    int               min_size,
    const char*       out_path);

}  // namespace nucleus_hier
