//
// Created by _ on 25-3-4.
//

#ifndef NCLIQUECOREDECOMPOSITION_H
#define NCLIQUECOREDECOMPOSITION_H
#include "../tree/MultiBranchTree.h"
#include <tbb/spin_mutex.h>

#include <ranges>
#include <map>

#include "graph/DynamicGraph.h"
#include "graph/DynamicGraphSet.h"
#include "graph/Graph.h"
#include "dataStruct/CliqueHashMap.h"

template<typename T>
class DynamicGraphSet;

extern double nCr[1001][401];

template<typename T_Key, typename T_Value>
struct ThreadSafeMap {
    std::map<T_Key, T_Value, std::greater<> > map;
    tbb::spin_mutex mutex;

    void insert(T_Key key, T_Value value) {
        // tbb::spin_mutex::scoped_lock lock(mutex);
        map[key] = value;
    }

    T_Value add(T_Key key, T_Value value) {
        // tbb::spin_mutex::scoped_lock lock(mutex);
        map[key] += value;
        return map[key];
    }

    T_Value get(T_Key key) {
        // tbb::spin_mutex::scoped_lock lock(mutex);
        return map[key];
    }

    friend std::ostream &operator<<(std::ostream &os, const ThreadSafeMap &map) {
        os << "{";
        bool first = true;
        for (const auto &pair: map.map) {
            if (!first) {
                os << ", ";
            }
            os << pair.first << ": " << pair.second;
            first = false;
        }
        os << "}";
        return os;
    }
};

inline double getCore(daf::CliqueSize povit, daf::CliqueSize keep, daf::CliqueSize r, daf::CliqueSize s) {
    return nCr[povit + keep - r][s - r];
}

void baseNucleusCoreDecomposition(const MultiBranchTree &tree, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> baseNucleusEdgeCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraph<daf::Size> &treeGraphV, daf::CliqueSize k);


double * baseNucleusCoreDecompositionLeaf(const MultiBranchTree &tree, daf::CliqueSize k);

void baseNucleusCoreDecompositionPar(const MultiBranchTree &tree, daf::CliqueSize k);

void baseNucleusCoreDecompositionParHash(const MultiBranchTree &tree, daf::CliqueSize k);
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraph<TreeGraphNode> &treeGraphV, daf::CliqueSize k);


std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);


std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > baseNucleusEdgeCoreDecompositionSet(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSetKCore(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);


double * NCliqueVertexCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > NucleusCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionHierarchy(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

void NucleusCoreDecompositionRemoveSclique(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);


std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionRClique(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// §105 M3: tuple-batched r-clique nucleus decomposition (PIVOTER_RUN_TUPLE_BATCH).
// Requires g_m1ClassOf (PIVOTER_M1 class computation) to be populated.
std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionRCliqueTupleBatch(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionRCliqueRef(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Single-thread optimized r≥3 (no OMP, integer arithmetic, Case A/C fast paths)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V2: Opt 1 — Leaf-Clique Reverse Index (eliminates Support hash lookups)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V3: Opt 2 — Clique-Leaf Mapping (eliminates Intersect hash set intersection)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V3(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V4: Opt 1+2 Combined (eliminates both Intersect and Support bottlenecks)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V4(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V5: Δ-Support with Positional Containment (bitset containment, no hash lookups in BK)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V5(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V6: Case C Extraction — r=2-style delta formula for pivot-only removal
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V6(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V7: Relaxed Case C — pivot-only conflict avoids BK entirely
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V7(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V8: Merged Init + Direct Bitset Containment (eliminates dual-index pass + robin_hood in BK)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V8(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V9: Combined — no treeGraphV, positional LeafCliqueEntry, no coverToVertex alloc
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V9(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V10: Single-pass init via buildWithCallback, reuse newEntries vector, addNodePresorted
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V10(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V11: Aggressive — single-pass init, batch bucketMove, compact cliqueLeafIds, skip dead in BK
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V11(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V12: On-demand r-clique peeling — eliminates leafCliqueInfo for memory savings
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V12(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V13: Path-Fragility Peeling — No r-Clique IDs, O(V + Σ|P|) space
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V13(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V16: Adaptive Hybrid — analytical d=1 + BK with MAX_SUBLEAVES cap + LeafDeath fallback
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V16(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V23: Incremental IE — analytical delta + pending list, no pathSplit
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V23(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V22: Explicit bipartite peeling — enumerate s-cliques, ZERO pathSplit
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V22(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V21: Sequential d=1 — analytical delta + Theorem-1 split + cache propagation, ZERO BK
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V21(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V18: Adaptive — lazy leafCliqueInfo (dead=V17 direct enum, BK=V11 cached), bucket+set PQ
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V18(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V20: V18 + sparse bucket PQ (map<double, vector<Size>>)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V20(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Region-based nucleus decomposition: peel regions (vertex profile groups), not individual r-cliques
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_Region(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionV2: general (r,s) via bipartite incidence between r-tuples and s-tuples
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionV2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionV2 Fast: MC deletion + V2 cascade (Phase 1 + Phase 2)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionV2_Fast(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Region + ST (V4): Region tuple compression + ST peeling mechanism
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionCPI (V3): Region Tuple + CPI Counting (no s-tuple enumeration)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionCPI NoCPI (ablation): replaces Step 4 with direct s-clique enumeration.
// All other optimizations (region FM, tuple, private cloud, dead-box peeling)
// unchanged. Used to measure the CPI formula's contribution in isolation.
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_NoCPI(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionCPI Cloud: Vandermonde-Chu + IE support, analytical cloud events
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_Cloud(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionCPI V3 + Hierarchy: outputs per-level connected component counts
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_Hierarchy(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionCPI V3 + Class-based Hierarchy (V3HC): DSU nodes = classes + FM MCs
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_HierarchyClass(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionCPI V3 low-memory variant (V3LM): replaces tuple->path index with
// class->path inverted index, reducing aux memory by 50-100x on graphs
// with poor class compression.
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_LowMem(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionCPI V3 low-memory + NoCPI Step 4 ablation.
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_LowMem_NoCPI(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// EventPeel: factored cell-death event engine (quotient peel). Replaces the
// dead-box B&B/union machinery with per-path alive-cell sets + rank-1
// factored death events (SigmodPlus §10.6). Same Steps 1-4 as V3LM.
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_EventPeel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// 4-tier RegND ablation dispatcher.  PIVOTER_TIER={1,2,3,4} selects:
//   T1 \regnd, T2 \regndplus, T3 \regndplusplus, T4 \regndstar (headline).
// Routes to LowMem / LowMem_NoCPI with PIVOTER_RECOMPUTE_PEEL,
// PIVOTER_V3_NO_PRIVATE, PIVOTER_VSAFE_CLOUD toggled as appropriate.
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_Tier(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionCPI V3 low-memory + class-based hierarchy post-processing.
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_LowMem_Hier(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionCPI V3B: Lazy Split optimization (s < 2r: unaffected tuples stay on parent)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionCPI_V2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// CCPath: lazy + eager hybrid via threshold antichain and ell/u bound
// restriction. Faithful port of solved.py (verified by 2000 random tests
// against brute-force enumeration).
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_CCPath(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// RegionExact: exact compressed peeling on r-class / s-class tuples
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_RegionExact(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Quotient Lab: report clean-leaf quotient compression, then run exact V20
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_QuotientLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V19: V18 + pure d-ary heap PQ (no overflow)
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V19(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V17: Hybrid — Reference's lightweight init + V11's bucket PQ
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V17(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V15: Incremental Exception Tracking (IET) — BK-free path splitting
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V15(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// ST V14: Analytical Sub-Leaf Construction (d=1 all-pivot) + V11 bitmask containment
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_V14(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Local CPI Vertex-Proxy H-index for r>=3 (no peeling, no BK, immutable tree)
std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLocalCPI_VP(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Local CPI Exact H-index for r>=3 (no peeling, no BK, immutable tree, exact s-clique enum)
std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLocalCPI_Exact(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Framework 1: Link-Graph Peeling (no BK, no tree mutation, pairwise weight subtraction)
std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLinkPeel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Framework 2: Link-Graph Local H-Index (no BK, no tree mutation, weighted H-index on leafCliqueInfo)
std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLinkLocal(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Batch Parallel version (new aggressive optimization)
namespace BatchParallel {
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRCliqueBatchParallel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);
}


std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionRCliquePar(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);


#endif //NCLIQUECOREDECOMPOSITION_H
// Level-based Parallel Algorithm
namespace LevelParallel {
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRCliqueLevelPar(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);
}

// Advanced Parallel Algorithm (Target: 3x+ speedup on 8 threads)
namespace AdvancedParallel {
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionAdvancedParallel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);
}

// Ultra Parallel Algorithm (Target: 3x+ speedup on 8 threads, optimized version)
namespace UltraParallel {
std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionUltraParallel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);
}


// Local H-index iterative convergence (r=1, no peeling)
double * NCliqueVertexCoreDecomposition_Local(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index naive version (full scan every iteration, no dirty queue)
double * NCliqueVertexCoreDecomposition_LocalNaive(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index V2: core-level enqueue filter + timestamp skip
double * NCliqueVertexCoreDecomposition_LocalV2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index V3: OpenMP parallel round-based convergence
double * NCliqueVertexCoreDecomposition_LocalV3(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index V4 declaration moved below ST_V2_Data definition.

// Single-thread optimized versions (no OMP overhead)
double * NCliqueVertexCoreDecomposition_ST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// ST V2: Tree-free R1 peeling via Augmented SDCT + dual CSR
// Does not need tree or treeGraphV — builds indices inline during SDCT.
double * NCliqueVertexCoreDecomposition_ST_V2(
    Graph &edgeGraph, daf::CliqueSize k);

// ST V2 split phases: Build must run BEFORE edgeGraph.beSingleEdge()
struct ST_V2_Data {
    struct VLeafEntry { daf::Size leafId; uint8_t isPivot; };
    struct LeafVtxEntry { daf::Size vertex; uint8_t isPivot; };

    daf::Size numVertices = 0;
    size_t numLeaves = 0;

    // Offsets are size_t (uint64) because cumulative offset = total Σ can
    // exceed uint32 max (4.29B) on billion-edge graphs at moderate s.
    // Friendster: Σ reaches 18.96B at s=5; uint32 wraps to ~700M, causing
    // OOB writes during dual-CSR fill. Element vectors stay 32-bit since
    // individual leaf/vertex IDs fit in uint32.
    std::vector<size_t> vtxLeafOff;
    std::vector<size_t> leafVtxOff;

    // Legacy unpacked layout (8 bytes/incidence including padding).
    // V2 uses these; V3 leaves them empty and uses the packed layout below.
    std::vector<VLeafEntry> vtxLeafData;
    std::vector<LeafVtxEntry> leafVtxData;

    // Packed layout used by V3 (4 bytes/incidence + 1 bit/incidence).
    // Splits the (id, isPivot) pair into two parallel arrays:
    //   - {vtx,leaf}LeafIds: one uint32 per incidence (the id)
    //   - {vtx,leaf}LeafIsPivot: 1 bit per incidence, packed into uint64 words
    // Saves ~50% memory on the dominant CSR arrays at billion-edge scale.
    std::vector<daf::Size> vtxLeafIds;
    std::vector<uint64_t>  vtxLeafIsPivot;
    std::vector<daf::Size> leafVtxIds;
    std::vector<uint64_t>  leafVtxIsPivot;

    std::vector<int> leafPivotCount;
    std::vector<int> leafNeedPivot;

    double *countingV = nullptr;
};

// Bit-array helpers: 1 bit per index, packed into uint64 words.
inline void STV3_setBit(std::vector<uint64_t>& bits, size_t i) {
    bits[i >> 6] |= (uint64_t(1) << (i & 63));
}
inline bool STV3_getBit(const std::vector<uint64_t>& bits, size_t i) {
    return (bits[i >> 6] >> (i & 63)) & 1;
}

ST_V2_Data NCliqueVertexCoreDecomposition_ST_V2_Build(
    Graph &edgeGraph, daf::CliqueSize k);
double * NCliqueVertexCoreDecomposition_ST_V2_Peel(
    ST_V2_Data &d, daf::CliqueSize k);

// ST V2 interleaved probe: measures feasibility of peeling during construction
void NCliqueVertexCoreDecomposition_ST_V2_InterleavedProbe(
    Graph &edgeGraph, daf::CliqueSize k);

// ST V3: V2 + segfault fix (sparse bucket queue keyed by int64_t).  Same
// data layout, same algorithm; fixes the dense-bucket overflow that
// crashed V2 on web-Stanford at s>=8.  Build phase is identical to V2,
// so we typedef ST_V3_Data = ST_V2_Data for binary compatibility.
using ST_V3_Data = ST_V2_Data;
ST_V3_Data NCliqueVertexCoreDecomposition_ST_V3_Build(
    Graph &edgeGraph, daf::CliqueSize k);
double * NCliqueVertexCoreDecomposition_ST_V3_Peel(
    ST_V3_Data &d, daf::CliqueSize k);
double * NCliqueVertexCoreDecomposition_ST_V3(
    Graph &edgeGraph, daf::CliqueSize k);

// Local H-index V4: async parallel with in-place updates.
// Reads paths from ST_V2_Data dual CSR (same layout as SPIN★) so the
// SPIN vs SPIN★ comparison reflects algorithm differences, not data
// structure overhead.
double * NCliqueVertexCoreDecomposition_LocalV4_Peel(
    ST_V2_Data &data, daf::CliqueSize k);
double * NCliqueVertexCoreDecomposition_LocalV4(
    Graph &edgeGraph, daf::CliqueSize k);

// ST V3 Lean: memory-tighter SPIN★ — drops per-leaf persistent state arrays
// (leafPivotCount, leafNeedPivot, leafAlive, leafRemainPivots) and recomputes
// np / old_rp from per-event leaf scans. Saves ~13 bytes/leaf permanent at
// the cost of ~1.5-2x peel CPU.
ST_V2_Data NCliqueVertexCoreDecomposition_ST_V3_Lean_Build(
    Graph &edgeGraph, daf::CliqueSize k);
double * NCliqueVertexCoreDecomposition_ST_V3_Lean_Peel(
    ST_V2_Data &d, daf::CliqueSize k);
double * NCliqueVertexCoreDecomposition_ST_V3_Lean(
    Graph &edgeGraph, daf::CliqueSize k);

// Interleaved construction-decomposition: peels vertices during SDCT construction
// Must be called BEFORE edgeGraph.beSingleEdge() (needs original graph for SDCT)
double * NCliqueVertexCoreDecomposition_Interleaved(
    Graph &edgeGraph, daf::CliqueSize k);

// InterleavedV2: Opt1 (fix drainHeap) + Opt2 (guard drain) + Opt3 (pre-reserve)
double * NCliqueVertexCoreDecomposition_InterleavedV2(
    Graph &edgeGraph, daf::CliqueSize k);

// OnDemand R=1: minimal-storage, CSR-based peeling with on-demand leaf state
// Must be called BEFORE edgeGraph.beSingleEdge() (needs original graph)
double * NCliqueVertexCoreDecomposition_OnDemand(
    Graph &edgeGraph, daf::CliqueSize k);

// OnDemand R=2: CSR-based init + Case A/C closed-form delta + Case B BK fallback
// Same interface as other R=2 variants (needs tree for Case B)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>>
PlusNucleusEdgeCoreDecompositionSet_OnDemand(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Tree-free R=2 edge nucleus decomposition via Augmented SDCT + dual CSR
// Must be called with tree for Case B BK fallback
struct ST_V2_EdgeData {
    struct VLeafEntry { daf::Size leafId; uint8_t isPivot; };
    struct LeafVtxEntry { daf::Size vertex; uint8_t isPivot; };
    struct LeafEdgeEntry { daf::Size edgeId; uint8_t edgeType; }; // 0=KK, 1=PP, 2=KP

    daf::Size numVertices = 0;
    daf::Size numEdges = 0;
    size_t numLeaves = 0;

    std::vector<daf::Size> vtxLeafOff;
    std::vector<VLeafEntry> vtxLeafData;
    std::vector<daf::Size> leafVtxOff;
    std::vector<LeafVtxEntry> leafVtxData;
    std::vector<std::vector<LeafEdgeEntry>> leafEdgeInfo;

    std::vector<int> leafPivotCount;
    std::vector<int> leafKeepCount;
    std::vector<int> leafNeedPivot;
    std::vector<long long> leafWKK, leafWPP, leafWKP;

    long long *countingKE = nullptr;
};

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_TreeFree(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// TreeFreeV2: Opt4 (no dual CSR) + Opt5 (Edge→Leaf index) + Opt6 (rebuild leafEdgeInfo after Case C)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_TreeFreeV2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_ST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V1: Edge-Leaf Dual Index (eliminates getEdgeCompressedId from peeling)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V1(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V2: V1 + Immutable Tree for Case A/C (eliminates removeNbr calls)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V3: V2 + Deferred Batch BucketMove (eliminates redundant bucket moves)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V3(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V4: Edge-Leaf Dual Index + Immutable Case A/C + Flat Phase 1
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V4(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V5: Sorted Vector Phase 1 (merge-scan intersection, ~8x less memory)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V5(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V6: Immutable Case A — leafAlive bitset skips removeNbr for dead leaves
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V6(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V7: Lazy Purge — staleCount-guided cleanup of dead leaf entries
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_ST_V7(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V7b: Periodic Compaction — purge dirty vertices every K iterations
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_ST_V7b(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V7c: Always-Purge — brute force cleanup before every intersection
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_ST_V7c(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V8: Case A LeafEdgeInfo fast-path (skip O(|leaf|^2) enumeration)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_ST_V8(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V8b: CSR init + deltaAccum batch bucketMove
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_ST_V8b(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 DCLP: Dual-CSR Lazy Peeling (edge-path transpose index + CSR-based Case C)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DCLP(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 DCLP2: DCLP + Complement BK for Case B d>=2 (BK on affected vertices only)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DCLP2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 Hybrid: DCLP + static edge->leaf transpose for original leaves + dynamic leaf graph
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_Hybrid(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 Hybrid Lab: isolated experimental copy; future hybrid experiments should go here
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_HybridLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 HardLeaf Lab: isolated DCLP-profile variant for hard-leaf sparsity/dominance experiments
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_HardLeafLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_HardLeafHybridLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DefectLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DefectRoutedLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DefectD2Lab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_DefectD2OrbitLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V9: Vertex→Leaf CSR + leafAlive + Vertex-Event Grouping
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V9(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 ST V10: Case C CSR fast-path + overflow CSR for Case B sub-leaves
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_ST_V10(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// R2 Lazy Verification: O(|P|) vertex updates + on-demand edge support recomputation
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>> PlusNucleusEdgeCoreDecompositionSet_Lazy(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionCorrect(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Local H-index r=2 (edge, no peeling, vertex-proxy)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>>
NCliqueEdgeCoreDecomposition_Local(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index r=2 V2: same exact H-index, but batch propagation by dirty leaf
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>>
NCliqueEdgeCoreDecomposition_LocalV2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index r=2 V3: batched propagation + edge->leaf CSR + s=3 fast path
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>>
NCliqueEdgeCoreDecomposition_LocalV3(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index r=2 V4: V3 + full edge-id hash lookup in the H-index hot path
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>>
NCliqueEdgeCoreDecomposition_LocalV4(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index r=2 V5: V3 + queue-time exact support-at-core fast check for s=3
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>>
NCliqueEdgeCoreDecomposition_LocalV5(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index r≥3 vertex-proxy (no peeling)
std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLocal(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);

// Local H-index r≥3 with exact s-clique enumeration (Online BK, no peeling)
std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRCliqueLocal_BK(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex = nullptr);
