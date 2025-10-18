//
// Created by _ on 24-7-24.
//
#pragma once
#ifndef SUBGRAPHMATCHING_GRAPH_H
#define SUBGRAPHMATCHING_GRAPH_H


#include <vector>
#include "Global/Global.h"
#include <span>
// #include <google/dense_hash_map>

class MultiBranchTree;
class TreeNode;

class Graph {
public:


    explicit Graph() = default;
    // @brief constructor
    // @param file_path the path of the graph file
    // @param singleEdge if true, 
    explicit Graph(const std::string &file_path, bool singleEdge = false);

    // std::vector<daf::Size> adj_list_offsets;
    // std::vector<daf::Size> adj_list;
    // daf::Size n;
    // daf::Size max_degree;
    explicit Graph(const MultiBranchTree &tree, const daf::StaticVector<TreeNode *> &leafList, daf::Size minK = 0);

    explicit Graph(const MultiBranchTree &tree,
        const daf::StaticVector<TreeNode *> &leafList,
        const Graph &treeGraphV,
        daf::Size minK = 0);


    explicit Graph(const MultiBranchTree &tree, const Graph &leafGraph, const daf::StaticVector<TreeNode *> &leafList, daf::Size r, daf::Size s);



    // explicit Graph(const MultiBranchTree &tree, const Graph &leafGraph, const daf::StaticVector<TreeNode *> &leafList, daf::Size r, daf::Size s);

    // remove copy and move
    Graph(const Graph &) = delete;
    Graph(Graph &&) {
        std::cout << "graph move constructor" << std::endl;
        // adj_list_offsets = std::move(adj_list_offsets);
        // adj_list = std::move(adj_list);
        // n = std::move(n);
        // max_degree = std::move(max_degree);
    }
    Graph &operator=(const Graph &) = delete;
    Graph &operator=(Graph &&) = delete;
    virtual ~Graph() {
        adj_list.free();
        adj_list_offsets.free();
    }
    //    print the Graph with <<
    [[nodiscard]] std::pair<daf::Size, daf::Size> getNbr(daf::Size node_id) const;

    std::pair<std::vector<daf::Size>::const_iterator, std::vector<daf::Size>::const_iterator> getNbrInter(
        daf::Size node_id) const;

    [[nodiscard]] auto getNbrVec(const daf::Size node_id) const {
        const auto [fst, snd] = getNbr(node_id);
        using diff_t = std::iterator_traits<decltype(adj_list.begin())>::difference_type;
        return std::ranges::subrange(
            adj_list.begin() + static_cast<diff_t>(fst),
            adj_list.begin() + static_cast<diff_t>(snd)
        );
    }

    [[nodiscard]] daf::Size getNbrCount(daf::Size node_id) const;

    [[nodiscard]] daf::Size getDegree(daf::Size node_id) const {
        return getNbrCount(node_id);
    }

    [[nodiscard]] daf::Size getGraphEdgeSize() const;

    [[nodiscard]] daf::Size getGraphNodeSize() const;

    [[nodiscard]] daf::Size getEdgeSize(daf::Size edge_id) const;

    [[nodiscard]] double getAvgDegree() const;

    [[nodiscard]] daf::Size getMaxDegree() const;

    void sortVertexByGivenOrder(const std::vector<daf::Size> &order);

    daf::Size getDegreeCount(daf::Size d) const;

    std::vector<daf::Size> sortByDegeneracyOrder();

    void printGraphInfo() const;

    void sortByDegree(bool reverse = true);

    void sortByBFSTraversal();

    friend std::ostream &operator<<(std::ostream &os, const Graph &g);

    void enumerateKCliques() const;

    void CliqueCounting();

    std::vector<std::vector<daf::Size>> BronKerboschPivot();

    void printGraph() const {
        std::cout << *this;
    }


    void printGraphPerV() const {
        for (daf::Size i = 0; i < n; ++i) {
            std::cout << i << ": ";
            for (daf::Size j = adj_list_offsets[i]; j < adj_list_offsets[i + 1]; ++j) {
                std::cout << adj_list[j] << " ";
            }
            std::cout << std::endl;
        }
    }

    // daf::Size getEdgeIndex(daf::Size u, daf::Size v) const {
    //     // u < v
    //     if (u > v) std::swap(u, v);
    //     // std::cout << "u: " << u << " v: " << v << std::endl;
    //     auto [begin, end] = getNbr(u);
    //     const auto it = std::lower_bound(adj_list.begin() + begin, adj_list.begin() + end, v);
    //     if (it != adj_list.end() && *it == v) {
    //         return it - adj_list.begin();
    //     }
    //     return -1;
    // }
    void beSingleEdge();

    void buildEdgeIdMap() {
        //  beSingleEdge()，adj_list_offsets/new_adj  u<v 
        // 1) 
        // edgeIdMap_.set_empty_key( std::numeric_limits<uint64_t>::max() );
        // edgeIdMap_.set_deleted_key( std::numeric_limits<uint64_t>::max() - 1 );
        // 
        // edgeIdMap_.resize( static_cast<size_t>(adj_list.c_size * 1.3) );
        edgeIdMap_.reserve(static_cast<size_t>(adj_list.c_size * 1.3));
        eidToNode.resize(adj_list.c_size);
        // 2)  CSR  map
        for (daf::Size u = 0; u < n; ++u) {
            bool useMap = (getNbrCount(u) >= 1000);
            for (daf::Size eid = adj_list_offsets[u]; eid < adj_list_offsets[u + 1]; ++eid) {
                daf::Size v = adj_list[eid];
                // beSingleEdge  v>u 
                if (useMap) {
                    uint64_t key = (static_cast<uint64_t>(u) << 32) | static_cast<uint64_t>(v);
                    edgeIdMap_[key] = eid;
                }
                eidToNode[eid] = u;
            }
        }
    }


    inline daf::Size getEdgeIndex(daf::Size u, daf::Size v) const noexcept {
        // 1. 
        if (u > v) {
            auto t = u; u = v; v = t;
        }

        // 2.  [b,e)
        auto [b, e] = getNbr(u);
        const daf::Size* data = adj_list.data();

        // 3.  v
        std::size_t lo = b, hi = e;
        while (lo < hi) {
            std::size_t mid = lo + ((hi - lo) >> 1);
            if (data[mid] < v) lo = mid + 1;
            else               hi = mid;
        }

        // 4. 
        if (lo < e && data[lo] == v) {
            return lo;
        }
        // error here
        std::cerr << "Error: [" << u << ", " << v << "] not found" << std::endl;
        return static_cast<daf::Size>(-1);
    }

    [[nodiscard]] daf::Size getEdgeCompressedId(daf::Size u, daf::Size v) const noexcept {
        if (u > v) std::swap(u,v);
        if (getNbrCount(u) < 1000) {
            return getEdgeIndex(u, v);
        }
        const uint64_t key = (u << 32) | v;
        const auto it = edgeIdMap_.find(key);
        if (it != edgeIdMap_.end()) return it->second;
        return static_cast<daf::Size>(-1);
    }

    std::pair<daf::Size, daf::Size> getEdgeById(daf::Size eid) const {
        if (eid >= adj_list.c_size) {
            std::cerr << "Error: eid out of range" << std::endl;
            return {static_cast<daf::Size>(-1), static_cast<daf::Size>(-1)};
        }
        return {eidToNode[eid], adj_list[eid]};
    }

    void initCore();

    daf::StaticVector<daf::Size> adj_list_offsets;
    daf::StaticVector<daf::Size> adj_list;
    daf::StaticVector<daf::Size> eidToNode;
    daf::StaticVector<daf::Size> coreV;
    daf::Size n;
    daf::Size max_degree;

    // key = (u<<32)|v, value =  ID
    // google::dense_hash_map<uint64_t, daf::Size> edgeIdMap_;
    robin_hood::unordered_flat_map<uint64_t, daf::Size> edgeIdMap_;
private:
    void BronKerboschPivotHelp(std::vector<daf::Size>& R,
                           std::vector<daf::Size>& P,
                           std::vector<daf::Size>& X,
                           std::vector<std::vector<daf::Size>>& cliques);
};


#endif //SUBGRAPHMATCHING_GRAPH_H
