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


// Single-thread optimized versions (no OMP overhead)
double * NCliqueVertexCoreDecomposition_ST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecompositionSet_ST(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize k);

std::vector<std::pair<std::vector<daf::Size>, int> > NucleusCoreDecompositionCorrect(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s);