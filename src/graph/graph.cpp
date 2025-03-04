//
// Created by 张文谦 on 24-7-24.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include "Graph.h"
#include "Global.h"

Graph::Graph(const std::string &file_path) {
    max_degree = 0;
    std::ifstream input_file(file_path);
    if (!input_file.is_open()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        exit(1);
    }

    daf::Size max_node_id = 0;
    std::string line;

    //    read the file to get the max_node_id
    while (std::getline(input_file, line)) {
        if (line[0] == '#' || line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        daf::Size node_id, nbr_id;
        iss >> node_id >> nbr_id;
        if (node_id == nbr_id) continue;
        max_node_id = std::max(max_node_id, std::max(node_id, nbr_id));
    }
    n = max_node_id + 1;

    std::vector<daf::Size> degree(max_node_id + 1, 0);

    input_file.clear();
    input_file.seekg(0, std::ios::beg);

    //    read the file to get the degree of each node
    while (std::getline(input_file, line)) {
        if (line[0] == '#' || line.empty()) {
            continue;
        }
        std::istringstream iss(line);
        daf::Size node_id, nbr_id;
        iss >> node_id >> nbr_id;

        if (node_id == nbr_id) continue;
        if (std::max(node_id, nbr_id) > max_node_id) {
            std::cout << "skip the first line: " << line << std::endl;
            continue;
        }

        ++degree[node_id];
        ++degree[nbr_id];
    }
    adj_list_offsets.resize(max_node_id + 2);
    adj_list_offsets[0] = 0;
    for (daf::Size i = 1; i < adj_list_offsets.size(); ++i) {
        adj_list_offsets[i] = adj_list_offsets[i - 1] + degree[i - 1];
        max_degree = std::max(max_degree, degree[i - 1]);
    }

    adj_list.resize(adj_list_offsets.back());
    input_file.clear();
    std::fill(degree.begin(), degree.end(), 0);
    input_file.seekg(0, std::ios::beg);

    //    read the file to get the adj_list
    while (std::getline(input_file, line)) {
        if (line[0] == '#' || line.empty()) {
            continue;
        }

        std::istringstream iss(line);
        daf::Size node_id, nbr_id;
        iss >> node_id >> nbr_id;
        if (node_id == nbr_id) continue;
        if (std::max(node_id, nbr_id) > max_node_id) {
            std::cout << "skip the first line: " << line << std::endl;
            continue;
        }
        adj_list[adj_list_offsets[node_id] + degree[node_id]++] = nbr_id;
        adj_list[adj_list_offsets[nbr_id] + degree[nbr_id]++] = node_id;
    }
    for (daf::Size i = 0; i < adj_list_offsets.size() - 1; ++i) {
        std::sort(adj_list.begin() + adj_list_offsets[i], adj_list.begin() + adj_list_offsets[i + 1]);
    }
    input_file.close();
    //    std::vector<std::pair<daf::Size, daf::Size>> edges;
    //
    //    edges.reserve(getGraphEdgeSize()-1);
    //    for (daf::Size i : vertex) {
    //        auto v = vertex[i];
    //        auto [begin, end] = getNbrVec(v);
    //        for (auto j = begin; j != end; ++j) {
    ////            if (v < *j) {
    //            edges.emplace_back(v, *j);
    ////            }
    //        }
    //    }
    //    std::sort(edges.begin(), edges.end());
    //    std::cout << edges << std::endl;
    //    std::cout << adj_list << std::endl;
}

std::ostream &operator<<(std::ostream &os, const Graph &g) {
    // std::vector<std::pair<daf::Size, daf::Size>> edges;

    // edges.reserve(g.getGraphEdgeSize()-1);
    // for (daf::Size v = 0; v < g.getGraphNodeSize(); ++v) {
    //     for (daf::Size j = g.adj_list_offsets[v]; j < g.adj_list_offsets[v + 1]; ++j) {
    //         if (v < g.adj_list[j]) {
    //             edges.emplace_back(v, g.adj_list[j]);
    //         }
    //     }
    // }
    //
    // for (const auto &edge : edges) {
    //     os << edge.first << " " << edge.second << std::endl;
    // }

    // for (daf::Size i = 0; i < g.getGraphNodeSize(); ++i) {
    //     os << i << ": ";
    //     for (daf::Size j = g.getNbr(i).first; j < g.getNbr(i).second; ++j) {
    //         os << g.adj_list[j] << " ";
    //     }
    //     os << std::endl;
    // }

    std::cout << g.adj_list_offsets << std::endl;
    std::cout << g.adj_list << std::endl;

    return os;
}


std::pair<daf::Size, daf::Size> Graph::getNbr(const daf::Size node_id) const {
    return {adj_list_offsets[node_id], adj_list_offsets[node_id + 1]};
}


daf::Size Graph::getNbrCount(daf::Size node_id) const {
    return adj_list_offsets[node_id + 1] - adj_list_offsets[node_id];
}

daf::Size Graph::getGraphEdgeSize() const {
    return adj_list.size() / 2;
}

daf::Size Graph::getGraphNodeSize() const {
    return adj_list_offsets.size() - 1;
}

daf::Size Graph::getEdgeSize(daf::Size edge_id) const {
    return adj_list_offsets[edge_id + 1] - adj_list_offsets[edge_id];
}

daf::Size Graph::getAvgDegree() const {
    return adj_list.size() / (adj_list_offsets.size() - 1);
}


daf::Size Graph::getMaxDegree() const {
    return max_degree;
}

void Graph::sortByDegree(bool reverse) {
    std::vector<daf::Size> order(getGraphNodeSize());
    for (daf::Size i = 0; i < getGraphNodeSize(); ++i) {
        order[i] = i;
    }
    std::ranges::sort(order, [this, &reverse](daf::Size a, daf::Size b) {
        return reverse ? getNbrCount(a) > getNbrCount(b) : getNbrCount(a) < getNbrCount(b);
    });
    sortVertexByGivenOrder(order);
}


void Graph::sortByBFSTraversal() {
    std::vector<bool> visited(getGraphNodeSize(), false);
    daf::Size order = 0;
    std::vector<daf::Size> ordered_vertex(getGraphNodeSize());
    for (daf::Size i = 0; i < getGraphNodeSize(); ++i) {
        if (visited[i]) {
            continue;
        }
        std::queue<daf::Size> q;
        q.push(i);
        visited[i] = true;
        while (!q.empty()) {
            daf::Size node_id = q.front();
            q.pop();
            ordered_vertex[order++] = node_id;
            for (daf::Size j = getNbr(node_id).first; j < getNbr(node_id).second; ++j) {
                if (!visited[adj_list[j]]) {
                    visited[adj_list[j]] = true;
                    q.push(adj_list[j]);
                }
            }
        }
    }
    std::ranges::reverse(ordered_vertex);
    sortVertexByGivenOrder(ordered_vertex);
}

void Graph::sortByDegeneracyOrdering() {
    // std::vector<daf::Size> vertex_order;
    // kCore::FlatArrayCoreDecomposition(*this, vertex_order);
    // std::ranges::reverse(vertex_order);
    // std::cout << "vertex_order: " << vertex_order << std::endl;
    // sortVertexByGivenOrder(vertex_order);
}

void Graph::printGraphInfo() const {
    std::cout << "Graph Info: " << std::endl;
    std::cout << "Node Size: " << getGraphNodeSize() << std::endl;
    std::cout << "Edge Size: " << getGraphEdgeSize() << std::endl;
    std::cout << "Max Degree: " << getMaxDegree() << std::endl;
    std::cout << "Avg Degree: " << getAvgDegree() << std::endl;
    std::cout << std::endl;
}

void Graph::sortVertexByGivenOrder(const std::vector<daf::Size> &vertex_order) {
    std::vector<daf::Size> vertex_order_index(getGraphNodeSize());
    for (daf::Size i = 0; i < getGraphNodeSize(); ++i) {
        vertex_order_index[vertex_order[i]] = i;
    }

    std::vector<std::vector<daf::Size> > adj_list_copy(getGraphNodeSize());
    for (daf::Size i = 0; i < getGraphNodeSize(); ++i) {
        adj_list_copy[i] = std::vector<daf::Size>(getNbrCount(i));
        for (daf::Size j = getNbr(i).first; j < getNbr(i).second; ++j) {
            adj_list_copy[i][j - getNbr(i).first] = vertex_order_index[adj_list[j]];
        }
    }

    for (daf::Size i = 0; i < getGraphNodeSize(); ++i) {
        std::ranges::sort(adj_list_copy[i]);
    }

    const auto n = getGraphNodeSize();


    for (daf::Size v = 0; v < n; ++v) {
        adj_list_offsets[v + 1] = adj_list_copy[vertex_order[v]].size();
    }

    for (daf::Size i = 1; i <= n; ++i) {
        adj_list_offsets[i] += adj_list_offsets[i - 1];
    }
    //    std::cout << adj_list_offsets << std::endl;

    for (daf::Size i = 0; i < n; ++i) {
        auto v = vertex_order[i];
        for (daf::Size j = 0; j < adj_list_copy[v].size(); ++j) {
            adj_list[adj_list_offsets[i] + j] = adj_list_copy[v][j];
        }
    }
}
