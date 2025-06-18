// Created by 张文谦 on 25-4-24.
// Updated 2025-06-15: switched underlying adjacency container
// from google::dense_hash_set to robin_hood::unordered_flat_set.

#include "DynamicGraphSet.h"

#include "DynamicBipartiteGraph.hpp"
#include "DynamicGraph.h"
#include "tree/MultiBranchTree.h"

// #include <robin_hood.h>
#include <limits>
#include <cstring>

// -----------------------------------------------------------------------------
//  Helper: reserve a robin_hood set with ~30% slack
// -----------------------------------------------------------------------------

template <typename SetT>
inline void initFlatSet(SetT& set, daf::Size expected) {
    set.reserve(static_cast<size_t>(expected * 1.3));
}

// -----------------------------------------------------------------------------
//  Template specialisation: leaf → node (TreeGraphNode ↔ leafId)
// -----------------------------------------------------------------------------

template<>
DynamicGraphSet<TreeGraphNode>::DynamicGraphSet(
        const daf::StaticVector<TreeNode*>& leafList, daf::Size minK) {
    adj_list.resize(leafList.size());

    for (auto* leaf : leafList) {
        if (leaf->MaxDeep < minK) continue;
        initFlatSet(adj_list[leaf->leafId], leaf->MaxDeep);
        auto* node = leaf;
        while (node->v != ROOTID) {
            adj_list[leaf->leafId].insert({node->v, node->isPivot});
            node = node->parent;
        }
    }
    removedNodes.reserve(leafList.size() / 2);
}

// -----------------------------------------------------------------------------
//  Template specialisation: node → leaf  (daf::Size ↔ leafId)
// -----------------------------------------------------------------------------

template<>
DynamicGraphSet<daf::Size>::DynamicGraphSet(
        const daf::StaticVector<TreeNode*>& leafList,
        daf::Size n,
        daf::Size minK) {
    adj_list.resize(n);

    // 1) count degrees first
    std::vector<unsigned int> degree(n, 0);
    for (auto* leaf : leafList) {
        if (leaf->MaxDeep < minK) continue;
        auto* node = leaf;
        while (node->v != ROOTID) {
            degree[node->v]++;
            node = node->parent;
        }
    }

    // 2) pre-reserve each adjacency set
    for (daf::Size i = 0; i < n; ++i) {
        initFlatSet(adj_list[i], degree[i]);
    }

    // 3) populate
    for (auto* leaf : leafList) {
        if (leaf->MaxDeep < minK) continue;
        auto* node = leaf;
        while (node->v != ROOTID) {
            adj_list[node->v].insert(leaf->leafId);
            node = node->parent;
        }
    }

    removedNodes.reserve(n / 2);
}

// -----------------------------------------------------------------------------
//  Template specialisation: node → leaf from DynamicGraph<TreeGraphNode>
// -----------------------------------------------------------------------------

template<>
DynamicGraphSet<daf::Size>::DynamicGraphSet(
        const DynamicGraph<TreeGraphNode>& treeGraph,
        daf::Size n,
        daf::Size minK) {
    adj_list.resize(n);

    // 1) degree counting
    std::vector<unsigned int> degree(n, 0);
    for (const auto& leaf : treeGraph.adj_list) {
        if (leaf.size() < minK) continue;
        for (const auto& i : leaf) {
            degree[i.v]++;
        }
    }

    // 2) reserve
    for (daf::Size i = 0; i < n; ++i) {
        initFlatSet(adj_list[i], degree[i]);
    }

    // 3) populate
    for (daf::Size leafId = 0; leafId < treeGraph.adj_list.size(); ++leafId) {
        const auto& leaf = treeGraph.adj_list[leafId];
        if (leaf.size() < minK) continue;
        for (const auto& i : leaf) {
            adj_list[i.v].insert(leafId);
        }
    }
    removedNodes.reserve(n / 2);
}

template<>
DynamicGraphSet<TreeGraphNode>::DynamicGraphSet(
        const DynamicGraph<TreeGraphNode>& treeGraph,
        daf::Size n,
        daf::Size minK) {
    adj_list.resize(n);
    std::vector<unsigned int> degree(n, 0);

    // 1) degree counting
    for (const auto& leaf : treeGraph.adj_list) {
        if (leaf.size() < minK) continue;
        for (const auto& nd : leaf) degree[nd.v]++;
    }

    // 2) reserve storage
    for (daf::Size i = 0; i < n; ++i) initFlatSet(adj_list[i], degree[i]);

    // 3) populate adjacency
    for (daf::Size leafId = 0; leafId < treeGraph.adj_list.size(); ++leafId) {
        const auto& leaf = treeGraph.adj_list[leafId];
        if (leaf.size() < minK) continue;
        for (const auto& nd : leaf) adj_list[nd.v].insert({leafId, nd.isPivot});
    }
    removedNodes.reserve(n / 2);
}
// -----------------------------------------------------------------------------
//  Template specialisation: TreeGraphNode variant using edgeCore filtering
// -----------------------------------------------------------------------------

template<>
DynamicGraphSet<TreeGraphNode>::DynamicGraphSet(
        const DynamicGraph<TreeGraphNode>& treeGraph,
        const Graph& edgeGraph,
        daf::Size n,
        daf::Size minK) {
    adj_list.resize(n);

    // 1) degree counting with core filtering
    std::vector<unsigned int> degree(n, 0);
    for (const auto& leaf : treeGraph.adj_list) {
        if (leaf.size() < minK) continue;
        daf::Size leafCore = std::numeric_limits<daf::Size>::max();
        for (const auto& i : leaf) {
            if (!i.isPivot) leafCore = std::min(leafCore, edgeGraph.coreV[i.v]);
        }
        for (const auto& i : leaf) {
            if (edgeGraph.coreV[i.v] > leafCore) break;
            degree[i.v]++;
        }
    }

    // 2) reserve
    for (daf::Size i = 0; i < n; ++i) {
        initFlatSet(adj_list[i], degree[i]);
    }

    // 3) populate
    for (daf::Size leafId = 0; leafId < treeGraph.adj_list.size(); ++leafId) {
        const auto& leaf = treeGraph.adj_list[leafId];
        if (leaf.size() < minK) continue;
        daf::Size leafCore = std::numeric_limits<daf::Size>::max();
        for (const auto& i : leaf) {
            if (!i.isPivot) leafCore = std::min(leafCore, edgeGraph.coreV[i.v]);
        }
        for (const auto& i : leaf) {
            if (edgeGraph.coreV[i.v] > leafCore) break;
            adj_list[i.v].insert({leafId, i.isPivot});
        }
    }
    removedNodes.reserve(n / 2);
}

// -----------------------------------------------------------------------------
//  Counting helpers remain unchanged
// -----------------------------------------------------------------------------

template<>
daf::StaticVector<double> DynamicGraphSet<TreeGraphNode>::cliqueCount() {
    auto maxdegree = this->maxDegree();
    daf::StaticVector<double> counts(maxdegree + 1);
    counts.c_size = maxdegree + 1;
    std::memset(counts.data, 0, (maxdegree + 1) * sizeof(double));

    for (const auto& leaf : this->adj_list) {
        int pivotCount = 0, nonPivotCount = 0;
        for (const auto& node : leaf) {
            (node.isPivot ? pivotCount : nonPivotCount)++;
        }
        const daf::CliqueSize rsize = pivotCount + nonPivotCount;
        for (daf::CliqueSize i = 0; i <= pivotCount; ++i) {
            const daf::Size k = rsize - i;
            counts[k] += nCr[pivotCount][i];
        }
    }
    return counts;
}


template<>
double DynamicGraphSet<TreeGraphNode>::cliqueCount(daf::Size k) {
    double counts = 0;
    for (const auto& leaf : this->adj_list) {
        int pivotCount = 0, nonPivotCount = 0;
        for (const auto& node : leaf) {
            (node.isPivot ? pivotCount : nonPivotCount)++;
        }
        if (k > pivotCount + nonPivotCount || nonPivotCount > k) continue;
        counts += nCr[pivotCount][k - nonPivotCount];
    }
    return counts;
}


template<>
daf::Size DynamicGraphSet<TreeGraphNode>::numBipartEdge(const daf::Size k) {
    daf::Size counts = 0;
    for (const auto& leaf : this->adj_list) {
        if (leaf.size() < k) continue;
        counts += nCr[leaf.size()][k];
    }
    return counts;
}
