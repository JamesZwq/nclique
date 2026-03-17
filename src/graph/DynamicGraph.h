//
// Created by _ on 25-4-24.
//
#pragma once
#ifndef DynamicGraph_H
#define DynamicGraph_H
#include "Graph.h"
#include "tree/MultiBranchTree.h"
#include <functional>
#include <cstddef>
#include <utility>
#include <type_traits>

// #include "dataStruct/CliqueHashMap.h"


template<typename T>
class DynamicGraph {
    public:
        explicit DynamicGraph(const daf::StaticVector<TreeNode *> &leafList, daf::Size minK);
        explicit DynamicGraph(const daf::StaticVector<TreeNode *> &leafList, daf::Size n, daf::Size minK);
        explicit DynamicGraph(const DynamicGraph<TreeGraphNode> &treeGraph,
                                daf::Size n,
                                daf::Size minK);

        explicit DynamicGraph() = default;
        explicit DynamicGraph(daf::Size n) {
            adj_list.reserve(n);
        }
    void addNbr(daf::Size node_id, T nbr) {
            if (adj_list[node_id].empty()) {
                adj_list[node_id].push_back(nbr);
                return;
            }
            if (nbr < adj_list[node_id].back()) {
                adj_list[node_id].push_back(nbr);
                std::sort(adj_list[node_id].begin(), adj_list[node_id].end());
            } else {
                adj_list[node_id].push_back(nbr);
            }
        }

        std::vector<T> &getNbr(daf::Size node_id) {
            return adj_list[node_id];
        }

        void removeNbr(daf::Size node_id, T nbr) {
            auto &nbrs = adj_list[node_id];
            auto it = std::lower_bound(nbrs.begin(), nbrs.end(), nbr);
            if (it != nbrs.end() && *it == nbr) {
                nbrs.erase(it);
            }
        }

    std::vector<T> &removeNbrs(daf::Size node_id, daf::StaticVector<daf::Size> &nbrs) {
            auto &lst = adj_list[node_id];
            std::size_t N = lst.size(), M = nbrs.size();
            std::size_t i = 0, j = 0, k = 0;

            // ：i  lst，j  nbrs；k
            while (i < N && j < M) {
                if (lst[i] == nbrs[j]) {
                    // ：（i++、j++）， k
                    ++i; ++j;
                }
                else {
                    // lst[i] ，
                    lst[k++] = std::move(lst[i++]);
                }
            }
            //  lst （nbrs ），
            while (i < N) {
                lst[k++] = std::move(lst[i++]);
            }
            //
            lst.resize(k);
            return lst;
        }

        daf::Size addNode(std::vector<T> &nbrs) {
            std::sort(nbrs.begin(), nbrs.end());
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

        // Push node_id to free list without clearing adj_list (caller already cleared it)
        void recycleNode(daf::Size node_id) {
            removedNodes.push_back(node_id);
        }

        void printGraphPerV() {
            // sort
            // for (auto &i: adj_list) {
            //     std::sort(i.begin(), i.end());
            // }
            for (daf::Size i = 0; i < adj_list.size(); ++i) {
                std::cout << i << ": ";
                for (const auto &nbr: adj_list[i]) {
                    std::cout << nbr << " ";
                }
                std::cout << std::endl;
            }
        }

        void printAdjStats() {
            if constexpr (std::is_same_v<T, TreeGraphNode>) {
                size_t maxLen = 0;
                for (const auto& adj : adj_list) {
                    if (adj.size() > maxLen) {
                        maxLen = adj.size();
                    }
                }

                std::cout << "Adj stats list (size " << maxLen + 1 << "): ";
                for (size_t n = 0; n <= maxLen; ++n) {
                    size_t count = 0;
                    for (const auto& adj : adj_list) {
                        if (adj.size() >= n) {
                            size_t nonPivotCount = 0;
                            for (const auto& treeNode : adj) {
                                if (!treeNode.isPivot) {
                                    nonPivotCount++;
                                }
                            }
                            if (nonPivotCount < n) {
                                count++;
                            }
                        }
                    }
                    std::cout << count << ";";
                }
                std::cout << std::endl;
            } else {
                std::cout << "printAdjStats only supported for DynamicGraph<TreeNode*>" << std::endl;
            }
        }


    [[nodiscard]] daf::StaticVector<double> cliqueCount() const;

    [[nodiscard]] daf::Size maxDegree() const {
        daf::Size max = 0;
        for (const auto &i: adj_list) {
            if (i.size() > max) {
                max = i.size();
            }
        }
        return max;
    }

    bool operator==(const DynamicGraph &other) const {
        auto sorted_this = adj_list;
        auto sorted_other = other.adj_list;
        for (auto &neighbors : sorted_this) {
            std::sort(neighbors.begin(), neighbors.end());
        }
        for (auto &neighbors : sorted_other) {
            std::sort(neighbors.begin(), neighbors.end());
        }
        std::sort(sorted_this.begin(), sorted_this.end());
        std::sort(sorted_other.begin(), sorted_other.end());
        return sorted_this == sorted_other;
    }

    bool operator!=(const DynamicGraph &other) const {
        return !(*this == other);
    }

    DynamicGraph clone() const {
        DynamicGraph clone;
        clone.adj_list = adj_list;
        clone.removedNodes = removedNodes;
        return clone;
    }

    [[nodiscard]] daf::Size size() const {
        return adj_list.size() - removedNodes.size();
    }

    [[nodiscard]]  double cliqueCount(daf::Size k) const;

        daf::Size numBipartEdge(daf::Size k);

        daf::StaticVector<double> cliqueCountPerV(daf::Size maxV, daf::Size k) const;
        daf::StaticVector<double> cliqueCountPerVAcc(daf::Size maxV, daf::Size k) const;

        std::vector<std::vector<T> > adj_list;
private:
        std::vector<daf::Size> removedNodes;

};


// // return r clique list, and the graph, r clique is the node, connect to treeGraph, like bipartite graph.
// std::pair<DynamicGraph<daf::Size>,std::vector<Clique>> getCliqueGraph(const DynamicGraph<TreeGraphNode> &treeGraph,
//         daf::CliqueSize r, daf::CliqueSize s);
#endif //DynamicGraph_H