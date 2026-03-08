#include "AdjacencyListGraph.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

AdjacencyListGraph::AdjacencyListGraph(daf::Size n) : n_(n), m_(0), adj_(n), max_degree_(0) {}

AdjacencyListGraph::AdjacencyListGraph(const std::string &file_path) : n_(0), m_(0), max_degree_(0) {
    std::ifstream input_file(file_path);
    if (!input_file.is_open()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        exit(1);
    }

    daf::Size max_node_id = 0;
    std::string line;
    std::vector<std::pair<daf::Size, daf::Size>> edges;

    while (std::getline(input_file, line)) {
        if (line[0] == '#' || line.empty()) {
            continue;
        }
        std::istringstream iss(line);
        daf::Size u, v;
        iss >> u >> v;
        if (u == v) continue;
        edges.emplace_back(u, v);
        max_node_id = std::max({max_node_id, u, v});
    }

    n_ = max_node_id + 1;
    adj_.resize(n_);

    for (const auto& edge : edges) {
        addEdge(edge.first, edge.second);
    }
}

void AdjacencyListGraph::addEdge(daf::Size u, daf::Size v) {
    if (u >= n_ || v >= n_) {
        daf::Size new_size = std::max(u, v) + 1;
        resize(new_size);
    }

    // Avoid parallel edges
    if (hasEdge(u, v)) return;

    adj_[u].push_back(v);
    adj_[v].push_back(u);
    std::sort(adj_[u].begin(), adj_[u].end());
    std::sort(adj_[v].begin(), adj_[v].end());

    m_++;
    max_degree_ = std::max({max_degree_, (daf::Size)adj_[u].size(), (daf::Size)adj_[v].size()});
}

void AdjacencyListGraph::removeEdge(daf::Size u, daf::Size v) {
    if (u >= n_ || v >= n_ || !hasEdge(u,v)) {
        return;
    }

    _remove_edge_directed(u, v);
    _remove_edge_directed(v, u);
    m_--;

    if (getDegree(u) == max_degree_ - 1 || getDegree(v) == max_degree_ - 1) {
        _update_max_degree();
    }
}

void AdjacencyListGraph::_remove_edge_directed(daf::Size u, daf::Size v) {
    auto& nbrs = adj_[u];
    auto it = std::lower_bound(nbrs.begin(), nbrs.end(), v);
    if (it != nbrs.end() && *it == v) {
        nbrs.erase(it);
    }
}

const std::vector<daf::Size>& AdjacencyListGraph::getNbrs(daf::Size node_id) const {
    return adj_[node_id];
}

daf::Size AdjacencyListGraph::getDegree(daf::Size node_id) const {
    return adj_[node_id].size();
}

daf::Size AdjacencyListGraph::getGraphEdgeSize() const {
    return m_;
}

daf::Size AdjacencyListGraph::getGraphNodeSize() const {
    return n_;
}

daf::Size AdjacencyListGraph::getMaxDegree() const {
    return max_degree_;
}

bool AdjacencyListGraph::hasEdge(daf::Size u, daf::Size v) const {
    if (u >= n_ || v >= n_) return false;
    const auto& nbrs = (getDegree(u) < getDegree(v)) ? adj_[u] : adj_[v];
    auto it = std::lower_bound(nbrs.begin(), nbrs.end(), (getDegree(u) < getDegree(v)) ? v : u);
    return (it != nbrs.end() && *it == ((getDegree(u) < getDegree(v)) ? v : u));
}

double AdjacencyListGraph::getAvgDegree() const {
    return n_ > 0 ? (2.0 * m_) / n_ : 0.0;
}

void AdjacencyListGraph::printGraph() const {
    for (daf::Size i = 0; i < n_; ++i) {
        std::cout << i << ": ";
        for (const auto& nbr : adj_[i]) {
            std::cout << nbr << " ";
        }
        std::cout << std::endl;
    }
}

void AdjacencyListGraph::resize(daf::Size new_n) {
    if (new_n > n_) {
        adj_.resize(new_n);
        n_ = new_n;
    }
}

void AdjacencyListGraph::_update_max_degree() {
    max_degree_ = 0;
    for (daf::Size i = 0; i < n_; ++i) {
        if (adj_[i].size() > max_degree_) {
            max_degree_ = adj_[i].size();
        }
    }
}