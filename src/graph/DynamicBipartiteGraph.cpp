//
// Created by _ on 25-5-18.
//

#include "DynamicBipartiteGraph.hpp"
extern double nCr[1001][401];
// DynamicBipartiteGraph::DynamicBipartiteGraph(DynamicGraph<TreeGraphNode> &treeGraph, Graph &edgeGraph) {
//     // 
//     nLeft_ = edgeGraph.adj_list.size();
//     nRight_ = treeGraph.adj_list.size();
//
//     // （）
//     daf::Size numEdge = 0;
//     for (const auto &clique: treeGraph.adj_list) {
//         numEdge += nCr[clique.size()][2];
//     }
//     bucketSize_ = numEdge / (nLeft_ + nRight_) + 1;
//
//     // 、
//     leftAdj_.reserve(nLeft_);
//     rightAdj_.reserve(nRight_);
//     for (int i = 0; i < nLeft_; ++i) {
//         leftAdj_.emplace_back(createEmptySet());
//     }
//     for (int i = 0; i < nRight_; ++i) {
//         rightAdj_.emplace_back(createEmptySet());
//     }
//
//     daf::Size cliqueId = 0;
//     for (const auto &clique: treeGraph.adj_list) {
//         for (daf::Size i = 0; i < clique.size(); ++i) {
//             for (daf::Size j = i + 1; j < clique.size(); ++j) {
//                 auto u = clique[i].v;
//                 auto v = clique[j].v;
//                 auto index = edgeGraph.getEdgeCompressedId(u, v);
//                 leftAdj_[index].insert(cliqueId);
//                 rightAdj_[cliqueId].insert(index);
//             }
//         }
//         ++cliqueId;
//     }
//
// }