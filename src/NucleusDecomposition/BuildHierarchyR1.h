// (1,s)-Nucleus hierarchy from the preserved static CPI.
//
// The default implementation follows the paper BuildHier algorithm:
// each CPI leaf emits the clique-entry events induced by its hold set,
// pivot set, and pivot requirement, then an elder-rule DSU scan builds
// the join tree.
//
// The old leaf-as-connector builder is kept as an explicit comparison
// mode (`PIVOTER_HIER_VERSION=leaf`) for local A/B benchmarking.
//
// Algorithm (post-peel):
//   1.  Read each CPI leaf's hold vertices, pivot vertices, and eta
//       from V3's packed leaf CSR.
//   2.  Emit the BuildHier clique-entry events for that leaf.
//   3.  Sort all events by descending level.
//   4.  Run an elder-rule DSU scan over graph vertices and output CSV.
//
// Memory: V3 peel itself stays at the paper's peel-state budget. This
// routine runs after peel and takes O(total emitted event size) scratch.

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

void buildAndDumpHierarchyFromCSRExact(
    const ST_V2_Data& v3data,
    const double*     coreV,
    daf::Size         numVertices,
    daf::Size         s,
    int               min_size,
    const char*       out_path);

void buildAndDumpHierarchyFromCSRLeafConnector(
    const ST_V2_Data& v3data,
    const double*     coreV,
    daf::Size         numVertices,
    daf::Size         s,
    int               min_size,
    const char*       out_path);

}  // namespace nucleus_hier
