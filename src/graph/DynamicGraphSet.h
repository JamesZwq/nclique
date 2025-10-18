//
// Created by _ on 25-4-24.
//
#pragma once
#ifndef DynamicGraphSet_H
#define DynamicGraphSet_H
#include "tree/MultiBranchTree.h"
#include <functional>
#include <cstddef>
#include <utility>

// #include <sparsehash/dense_hash_set>

#include "DynamicGraph.h"
#include "dataStruct/CliqueHashMap.h"
#include "Global/ostreamOverload.hpp"

template<typename T>
class DynamicGraphSet {
public:
    explicit DynamicGraphSet(const daf::StaticVector<TreeNode *> &leafList, daf::Size minK);

    explicit DynamicGraphSet(const daf::StaticVector<TreeNode *> &leafList, daf::Size n, daf::Size minK);

    explicit DynamicGraphSet(const DynamicGraph<TreeGraphNode> &treeGraph,
                             daf::Size n,
                             daf::Size minK);

    explicit DynamicGraphSet() = default;

    explicit DynamicGraphSet(daf::Size n) {
        adj_list.reserve(n);
    }

    explicit DynamicGraphSet(const DynamicGraph<TreeGraphNode> &treeGraph, const Graph &edgeGraph, daf::Size n, daf::Size minK);

    void addNbr(daf::Size node_id, T nbr) {
        // if (adj_list[node_id].empty()) {
        //     adj_list[node_id].insert(nbr);
        //     return;
        // }
        // if (nbr < adj_list[node_id].back()) {
        //     adj_list[node_id].insert(nbr);
        //     std::sort(adj_list[node_id].begin(), adj_list[node_id].end());
        // } else {
        //     adj_list[node_id].insert(nbr);
        // }
        adj_list[node_id].insert(nbr);
    }

    robin_hood::unordered_flat_set<T> &getNbr(daf::Size node_id) {
        return adj_list[node_id];
    }

    void removeNbr(daf::Size node_id, T nbr) {
        auto &nbrs = adj_list[node_id];
        nbrs.erase(nbr);
    }

    robin_hood::unordered_flat_set<T> &removeNbrs(daf::Size node_id, robin_hood::unordered_flat_set<daf::Size> &nbrs) {
        auto &s = adj_list[node_id];
        // s.erase(nbrs);
        for (const auto &nbr: nbrs) {
            s.erase(nbr);
        }
        return s;
    }

    [[nodiscard]] daf::Size addNode(std::vector<T> &nbrs) {
        if (removedNodes.empty()) {
            adj_list.emplace_back(std::move(nbrs));
            return adj_list.size() - 1;
        }
        daf::Size node_id = removedNodes.back();
        removedNodes.pop_back();
        adj_list[node_id] = std::move(nbrs);
        return node_id;
    }

    void removeNode(daf::Size node_id) {
#ifndef NDEBUG
            if (node_id >= adj_list.size()) {
                throw std::out_of_range("Node ID out of range.");
            }
#endif
        removedNodes.push_back(node_id);
        adj_list[node_id].clear();
    }

    void printGraphPerV() {
        // sort

        for (daf::Size i = 0; i < adj_list.size(); ++i) {
            std::cout << i << ": ";
            for (const auto &nbr: adj_list[i]) {
                std::cout << nbr << " ";
            }
            std::cout << std::endl;
        }
    }


    daf::StaticVector<double> cliqueCount();

    [[nodiscard]] daf::Size maxDegree() const {
        daf::Size max = 0;
        for (const auto &i: adj_list) {
            if (i.size() > max) {
                max = i.size();
            }
        }
        return max;
    }

    DynamicGraphSet clone() const {
        DynamicGraphSet clone;
        clone.adj_list = adj_list;
        clone.removedNodes = removedNodes;
        return clone;
    }

    double cliqueCount(daf::Size k);

    daf::Size numBipartEdge(daf::Size k);

    std::vector<robin_hood::unordered_flat_set<T> > adj_list;
    // std::vector<robin_hood::unordered_flat_set<T> > adj_list;
    // robin_hood::unordered_flat_set<><T>

    //  ID
private:
    std::vector<daf::Size> removedNodes;
};


// return r clique list, and the graph, r clique is the node, connect to treeGraph, like bipartite graph.
// std::pair<DynamicGraphSet<daf::Size>, std::vector<Clique> > getCliqueGraph(
//     const DynamicGraphSet<TreeGraphNode> &treeGraph,
//     daf::CliqueSize r, daf::CliqueSize s);
#endif //DynamicGraphSet_H
