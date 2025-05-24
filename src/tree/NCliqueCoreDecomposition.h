//
// Created by 张文谦 on 25-3-4.
//

#ifndef NCLIQUECOREDECOMPOSITION_H
#define NCLIQUECOREDECOMPOSITION_H
#include "MultiBranchTree.h"
#include <tbb/spin_mutex.h>

#include <ranges>
#include <map>

#include "graph/DynamicGraph.h"

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

void baseNucleusCoreDecompositionLeaf(const MultiBranchTree &tree, daf::CliqueSize k);

void baseNucleusCoreDecompositionPar(const MultiBranchTree &tree, daf::CliqueSize k);

void baseNucleusCoreDecompositionParHash(const MultiBranchTree &tree, daf::CliqueSize k);
std::vector<std::pair<std::pair<daf::Size, daf::Size>, int> > PlusNucleusEdgeCoreDecomposition(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraph<TreeGraphNode> &treeGraphV, daf::CliqueSize k);
#endif //NCLIQUECOREDECOMPOSITION_H
