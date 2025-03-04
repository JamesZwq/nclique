//
// Created by 张文谦 on 24-7-24.
//

#ifndef SUBGRAPHMATCHING_GRAPH_H
#define SUBGRAPHMATCHING_GRAPH_H


#include <vector>
#include "Global.h"
#include <span>

class Graph {
public:
    explicit Graph(const std::string &file_path);

    //    print the Graph with <<
    [[nodiscard]] std::pair<daf::Size, daf::Size> getNbr(daf::Size node_id) const;

    [[nodiscard]] auto getNbrVec(const daf::Size node_id) const {
        const auto [fst, snd] = getNbr(node_id);
        return std::ranges::subrange(adj_list.begin() + fst, adj_list.begin() + snd);
    }

    [[nodiscard]] daf::Size getNbrCount(daf::Size node_id) const;

    [[nodiscard]] daf::Size getGraphEdgeSize() const;

    [[nodiscard]] daf::Size getGraphNodeSize() const;

    [[nodiscard]] daf::Size getEdgeSize(daf::Size edge_id) const;

    [[nodiscard]] daf::Size getAvgDegree() const;

    [[nodiscard]] daf::Size getMaxDegree() const;

    void sortVertexByGivenOrder(const std::vector<daf::Size> &order);

    void sortByDegeneracyOrdering();

    void printGraphInfo() const;

    void sortByDegree(bool reverse = true);

    void sortByBFSTraversal();

    friend std::ostream &operator<<(std::ostream &os, const Graph &g);

    void enumerateKCliques() const;

    void CliqueCounting();

    std::vector<daf::Size> adj_list_offsets;
    std::vector<daf::Size> adj_list;
    daf::Size n;
    daf::Size max_degree;
};


#endif //SUBGRAPHMATCHING_GRAPH_H
