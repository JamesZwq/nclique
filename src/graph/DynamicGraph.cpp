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


template<>
DynamicGraph<daf::Size>::DynamicGraph(const DynamicGraph<TreeGraphNode> &treeGraph,
                                daf::Size n,
                                daf::Size minK) {
    adj_list.resize(n);
    unsigned int *degree = new unsigned int[n];
    memset(degree, 0, n * sizeof(unsigned int));
    for (auto leaf: treeGraph.adj_list) {
        if (leaf.size() < minK) continue;
        auto node = leaf;
        for (const auto &i: leaf) {
            degree[i.v]++;
        }
    }

    for (daf::Size i = 0; i < n; ++i) {
        adj_list[i].reserve(degree[i]);
    }

    // for (auto leaf: treeGraph.adj_list) {
    //     if (leaf.size() < minK) continue;
    //     // while (node->v != ROOTID) {
    //     //     adj_list[node->v].push_back(leaf->leafId);
    //     //     node = node->parent;
    //     // }
    //     for (const auto &i: leaf) {
    //         adj_list[i.v].push_back(leaf[0].leafId);
    //     }
    // }
    for (daf::Size leafId = 0; leafId < treeGraph.adj_list.size(); ++leafId) {
        const auto &leaf = treeGraph.adj_list[leafId];
        if (leaf.size() < minK) continue;
        for (const auto &i: leaf) {
            adj_list[i.v].push_back(leafId);
        }
    }
    for (auto &i: adj_list) {
        std::ranges::sort(i);
    }
    removedNodes.reserve(n/2);
    delete[] degree;
}


template<>
DynamicGraph<TreeGraphNode>::DynamicGraph(const DynamicGraph<TreeGraphNode> &treeGraph,
                                daf::Size n,
                                daf::Size minK) {
    adj_list.resize(n);
    unsigned int *degree = new unsigned int[n];
    memset(degree, 0, n * sizeof(unsigned int));
    for (auto leaf: treeGraph.adj_list) {
        if (leaf.size() < minK) continue;
        auto node = leaf;
        for (const auto &i: leaf) {
            degree[i.v]++;
        }
    }

    for (daf::Size i = 0; i < n; ++i) {
        adj_list[i].reserve(degree[i]);
    }

    // for (auto leaf: treeGraph.adj_list) {
    //     if (leaf.size() < minK) continue;
    //     // while (node->v != ROOTID) {
    //     //     adj_list[node->v].push_back(leaf->leafId);
    //     //     node = node->parent;
    //     // }
    //     for (const auto &i: leaf) {
    //         adj_list[i.v].push_back(leaf[0].leafId);
    //     }
    // }
    for (daf::Size leafId = 0; leafId < treeGraph.adj_list.size(); ++leafId) {
        const auto &leaf = treeGraph.adj_list[leafId];
        if (leaf.size() < minK) continue;
        for (const auto &i: leaf) {
            adj_list[i.v].push_back({leafId, i.isPivot});
        }
    }
    for (auto &i: adj_list) {
        std::ranges::sort(i);
    }
    removedNodes.reserve(n/2);
    delete[] degree;
}

template<>
daf::StaticVector<double> DynamicGraph<TreeGraphNode>::cliqueCount() {
    auto maxdegree = this->maxDegree();
    daf::StaticVector<double> counts(maxdegree + 1);
    counts.c_size = maxdegree + 1;
    memset(counts.data, 0, (maxdegree + 1) * sizeof(double));
    for (auto leaf: this->adj_list) {
        int pivotCount = 0, nonPivotCount = 0;
        for (const auto &node: leaf) {
            if (node.isPivot) {
                pivotCount++;
            } else {
                nonPivotCount++;
            }
        }
        const daf::CliqueSize rsize = pivotCount + nonPivotCount;
        for (daf::CliqueSize i = 0; i <= pivotCount; i++) {
            const daf::Size k = rsize - i;
            counts[k] += nCr[pivotCount][i];
        }
    }
    return counts;
}


template<>
double DynamicGraph<TreeGraphNode>::cliqueCount(daf::Size k) {
    double counts = 0;
    for (const auto& leaf: this->adj_list) {
        int pivotCount = 0, nonPivotCount = 0;
        for (const auto &node: leaf) {
            if (node.isPivot) {
                pivotCount++;
            } else {
                nonPivotCount++;
            }
        }
        if (k > pivotCount + nonPivotCount || nonPivotCount > k) {
            continue;
        }
        counts += nCr[pivotCount][k - nonPivotCount];
    }
    return counts;
}


template<>
daf::Size DynamicGraph<TreeGraphNode>::numBipartEdge(const daf::Size k) {
    daf::Size counts = 0;
    for (const auto& leaf: this->adj_list) {
        if (leaf.size() < k) continue;
        counts += nCr[leaf.size()][k];
    }
    return counts;
}

