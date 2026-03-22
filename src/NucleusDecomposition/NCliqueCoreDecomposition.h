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
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

std::vector<std::pair<std::vector<daf::Size>, double> > NucleusCoreDecompositionHierarchy(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

void NucleusCoreDecompositionRemoveSclique(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);


std::vector<std::pair<std::vector<daf::Size>, int> > NucleusCoreDecompositionRClique(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

std::vector<std::pair<std::vector<daf::Size>, int> > NucleusCoreDecompositionRCliqueRef(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// Single-thread optimized r≥3 (no OMP, integer arithmetic, Case A/C fast paths)
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V2: Opt 1 — Leaf-Clique Reverse Index (eliminates Support hash lookups)
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V2(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V3: Opt 2 — Clique-Leaf Mapping (eliminates Intersect hash set intersection)
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V3(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V4: Opt 1+2 Combined (eliminates both Intersect and Support bottlenecks)
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V4(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V5: Δ-Support with Positional Containment (bitset containment, no hash lookups in BK)
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V5(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V6: Case C Extraction — r=2-style delta formula for pivot-only removal
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V6(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V7: Relaxed Case C — pivot-only conflict avoids BK entirely
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V7(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V8: Merged Init + Direct Bitset Containment (eliminates dual-index pass + robin_hood in BK)
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V8(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V9: Combined — no treeGraphV, positional LeafCliqueEntry, no coverToVertex alloc
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V9(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V10: Single-pass init via buildWithCallback, reuse newEntries vector, addNodePresorted
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V10(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// ST V11: Aggressive — single-pass init, batch bucketMove, compact cliqueLeafIds, skip dead in BK
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRClique_ST_V11(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// Local CPI Vertex-Proxy H-index for r>=3 (no peeling, no BK, immutable tree)
std::vector<std::pair<std::vector<daf::Size>, int>>
NucleusCoreDecompositionRCliqueLocalCPI_VP(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// Local CPI Exact H-index for r>=3 (no peeling, no BK, immutable tree, exact s-clique enum)
std::vector<std::pair<std::vector<daf::Size>, int>>
NucleusCoreDecompositionRCliqueLocalCPI_Exact(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// Framework 1: Link-Graph Peeling (no BK, no tree mutation, pairwise weight subtraction)
std::vector<std::pair<std::vector<daf::Size>, int>>
NucleusCoreDecompositionRCliqueLinkPeel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// Framework 2: Link-Graph Local H-Index (no BK, no tree mutation, weighted H-index on leafCliqueInfo)
std::vector<std::pair<std::vector<daf::Size>, int>>
NucleusCoreDecompositionRCliqueLinkLocal(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// Batch Parallel version (new aggressive optimization)
namespace BatchParallel {
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRCliqueBatchParallel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);
}


std::vector<std::pair<std::vector<daf::Size>, int> > NucleusCoreDecompositionRCliquePar(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);


#endif //NCLIQUECOREDECOMPOSITION_H
// Level-based Parallel Algorithm
namespace LevelParallel {
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRCliqueLevelPar(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);
}

// Advanced Parallel Algorithm (Target: 3x+ speedup on 8 threads)
namespace AdvancedParallel {
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionAdvancedParallel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);
}

// Ultra Parallel Algorithm (Target: 3x+ speedup on 8 threads, optimized version)
namespace UltraParallel {
std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionUltraParallel(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);
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

// Local H-index V4: async parallel with in-place updates
double * NCliqueVertexCoreDecomposition_LocalV4(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Single-thread optimized versions (no OMP overhead)
double * NCliqueVertexCoreDecomposition_ST(
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

std::vector<std::pair<std::vector<daf::Size>, int> > NucleusCoreDecompositionCorrect(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// Local H-index r=2 (edge, no peeling, vertex-proxy)
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int>>
NCliqueEdgeCoreDecomposition_Local(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

// Local H-index r≥3 vertex-proxy (no peeling)
std::vector<std::pair<std::vector<daf::Size>, int>>
NucleusCoreDecompositionRCliqueLocal(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);

// Local H-index r≥3 with exact s-clique enumeration (Online BK, no peeling)
std::vector<std::pair<std::vector<daf::Size>, int>>
NucleusCoreDecompositionRCliqueLocal_BK(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);