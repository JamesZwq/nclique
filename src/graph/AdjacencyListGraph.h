#pragma once
#ifndef SUBGRAPHMATCHING_ADJACENCYLIST_GRAPH_H
#define SUBGRAPHMATCHING_ADJACENCYLIST_GRAPH_H

#include <vector>
#include <string>
#include "Global/Global.h"

class AdjacencyListGraph {
public:
    // --- Constructors and Destructor ---
    explicit AdjacencyListGraph(daf::Size n = 0);
    explicit AdjacencyListGraph(const std::string &file_path);
    AdjacencyListGraph(const AdjacencyListGraph &) = delete;
    AdjacencyListGraph &operator=(const AdjacencyListGraph &) = delete;
    virtual ~AdjacencyListGraph() = default;

    // --- Dynamic Operations ---
    void addEdge(daf::Size u, daf::Size v);
    void removeEdge(daf::Size u, daf::Size v);

    // --- Core Graph API ---
    [[nodiscard]] const std::vector<daf::Size>& getNbrs(daf::Size node_id) const;
    [[nodiscard]] daf::Size getDegree(daf::Size node_id) const;
    [[nodiscard]] daf::Size getGraphEdgeSize() const;
    [[nodiscard]] daf::Size getGraphNodeSize() const;
    [[nodiscard]] daf::Size getMaxDegree() const;
    [[nodiscard]] bool hasEdge(daf::Size u, daf::Size v) const;
    [[nodiscard]] double getAvgDegree() const;

    // --- Utility Functions ---
    void printGraph() const;
    void resize(daf::Size new_n);

private:
    void _remove_edge_directed(daf::Size u, daf::Size v);
    void _update_max_degree();

    daf::Size n_; // Number of vertices
    daf::Size m_; // Number of edges
    std::vector<std::vector<daf::Size>> adj_;
    daf::Size max_degree_;
};

#endif //SUBGRAPHMATCHING_ADJACENCYLIST_GRAPH_H