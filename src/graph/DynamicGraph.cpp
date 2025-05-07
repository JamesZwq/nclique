//
// Created by 张文谦 on 25-4-24.
//

#include "DynamicGraph.h"

#include "tree/MultiBranchTree.h"

// template<typename T>
// leaf to node
template<>
DynamicGraph<TreeGraphNode>::DynamicGraph(const daf::StaticVector<TreeNode *> &leafList,
                                daf::Size minK) {
    adj_list.resize(leafList.size());
    for (auto leaf: leafList) {
        if (leaf->MaxDeep < minK) continue;
        adj_list[leaf->leafId].reserve(leaf->MaxDeep);
        auto node = leaf;
        while (node->v != ROOTID) {
            adj_list[leaf->leafId].emplace_back(node->v, node->isPivot);
            node = node->parent;
        }
    }
    removedNodes.reserve(leafList.size()/2);
    for (auto &i: adj_list) {
        std::ranges::sort(i);
    }
}

// node to leaf
template<>
DynamicGraph<daf::Size>::DynamicGraph(const daf::StaticVector<TreeNode *> &leafList,
                                daf::Size n,
                                daf::Size minK) {
    adj_list.resize(n);
    unsigned int *degree = new unsigned int[n];
    memset(degree, 0, n * sizeof(unsigned int));
    for (auto leaf: leafList) {
        if (leaf->MaxDeep < minK) continue;
        auto node = leaf;
        while (node->v != ROOTID) {
            degree[node->v]++;
            node = node->parent;
        }
    }

    for (daf::Size i = 0; i < n; ++i) {
        adj_list[i].reserve(degree[i]);
    }

    for (auto leaf: leafList) {
        if (leaf->MaxDeep < minK) continue;
        auto node = leaf;
        while (node->v != ROOTID) {
            adj_list[node->v].push_back(leaf->leafId);
            node = node->parent;
        }
    }
    for (auto &i: adj_list) {
        std::ranges::sort(i);
    }
    removedNodes.reserve(n/2);
    delete[] degree;
}

