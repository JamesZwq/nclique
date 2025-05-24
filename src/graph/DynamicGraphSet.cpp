//
// Created by 张文谦 on 25-4-24.
//

#include "DynamicGraphSet.h"

#include "DynamicBipartiteGraph.hpp"
#include "DynamicGraph.h"

#include "tree/MultiBranchTree.h"

// template<typename T>
// leaf to node

void initDenseHashSet(google::dense_hash_set<TreeGraphNode> &set, daf::Size n) {
    set.set_empty_key(TreeGraphNode::EMPTYKEY);
    set.set_deleted_key(TreeGraphNode::DELETEDKEY);
    set.resize(n * 1.3);
}
void initDenseHashSet(google::dense_hash_set<daf::Size> &set, daf::Size n) {
    set.set_empty_key(std::numeric_limits<daf::Size>::max());
    set.set_deleted_key(std::numeric_limits<daf::Size>::max() - 1);
    set.resize(n * 1.3);
}

template<>
DynamicGraphSet<TreeGraphNode>::DynamicGraphSet(const daf::StaticVector<TreeNode *> &leafList,
                                                daf::Size minK) {
    adj_list.resize(leafList.size());
    for (auto leaf: leafList) {
        if (leaf->MaxDeep < minK) continue;
        initDenseHashSet(adj_list[leaf->leafId], leaf->MaxDeep);
        auto node = leaf;
        while (node->v != ROOTID) {
            adj_list[leaf->leafId].insert({node->v, node->isPivot});
            node = node->parent;
        }
    }
    removedNodes.reserve(leafList.size()/2);

}

// node to leaf
template<>
DynamicGraphSet<daf::Size>::DynamicGraphSet(const daf::StaticVector<TreeNode *> &leafList,
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
        initDenseHashSet(adj_list[i], degree[i]);
    }

    for (auto leaf: leafList) {
        if (leaf->MaxDeep < minK) continue;
        auto node = leaf;
        while (node->v != ROOTID) {
            adj_list[node->v].insert(leaf->leafId);
            node = node->parent;
        }
    }

    removedNodes.reserve(n/2);
    delete[] degree;
}


template<>
DynamicGraphSet<daf::Size>::DynamicGraphSet(const DynamicGraphSet<TreeGraphNode> &treeGraph,
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
        // adj_list[i].reserve(degree[i]);
        initDenseHashSet(adj_list[i], degree[i]);
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
            adj_list[i.v].insert(leafId);
        }
    }
    removedNodes.reserve(n/2);
    delete[] degree;
}


template<>
DynamicGraphSet<TreeGraphNode>::DynamicGraphSet(const DynamicGraphSet<TreeGraphNode> &treeGraph,
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
        // adj_list[i].reserve(degree[i]);
        initDenseHashSet(adj_list[i], degree[i]);
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
            adj_list[i.v].insert({leafId, i.isPivot});
        }
    }
    removedNodes.reserve(n/2);
    delete[] degree;
}

template<>
daf::StaticVector<double> DynamicGraphSet<TreeGraphNode>::cliqueCount() {
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
double DynamicGraphSet<TreeGraphNode>::cliqueCount(daf::Size k) {
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
daf::Size DynamicGraphSet<TreeGraphNode>::numBipartEdge(const daf::Size k) {
    daf::Size counts = 0;
    for (const auto& leaf: this->adj_list) {
        if (leaf.size() < k) continue;
        counts += nCr[leaf.size()][k];
    }
    return counts;
}

