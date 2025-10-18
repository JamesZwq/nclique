// //
// // Created by _ on 25-5-18.
// //
//
// #ifndef DYNAMICBIPARTITEGRAPH_HPP
// #define DYNAMICBIPARTITEGRAPH_HPP
//
//
//
// #pragma once
// #include <vector>
// #include <stdexcept>
// #include <google/dense_hash_set>
//
// #include "DynamicGraph.h"
// #include "Global/Global.h"
// using google::dense_hash_set;
// extern double nCr[1001][401];
// class DynamicBipartiteGraph {
// public:
//     //  dense_hash_set    
//     static constexpr int EMPTY_KEY   = -1;
//     static constexpr int DELETED_KEY = -2;
//
//     // ：nLeft, nRight ；maxEdges 
//     DynamicBipartiteGraph(daf::Size nLeft, daf::Size nRight, size_t maxEdges)
//         : nLeft_(nLeft), nRight_(nRight)
//     {
//         if (nLeft_ <= 0 || nRight_ <= 0)
//             throw std::invalid_argument("Node counts must be positive");
//
//         // （）
//         bucketSize_ = maxEdges / (nLeft_ + nRight_) + 1;
//
//         // 、
//         leftAdj_.reserve(nLeft_);
//         rightAdj_.reserve(nRight_);
//         for (int i = 0; i < nLeft_; ++i) {
//             leftAdj_.emplace_back(createEmptySet());
//         }
//         for (int i = 0; i < nRight_; ++i) {
//             rightAdj_.emplace_back(createEmptySet());
//         }
//     }
//
//     // DynamicBipartiteGraph(DynamicGraph<TreeGraphNode> &treeGraph, Graph &edgeGraph);
//
//     // ， ID（）
//     int addLeftNode() {
//         leftAdj_.emplace_back(createEmptySet());
//         return nLeft_++;
//     }
//
//     // ， ID（）
//     int addRightNode() {
//         rightAdj_.emplace_back(createEmptySet());
//         return nRight_++;
//     }
//
//     //  (u->v)， u/v ， addLeftNode/addRightNode
//     bool addEdge(int u, int v) {
//         if (u > nLeft_) {
//             if (u == nLeft_) {
//                 addLeftNode();
//             } else {
//                 throw std::out_of_range("Left node ID out of range");
//             }
//         }
//
//         if (v > nRight_) {
//             if (v == nRight_) {
//                 addRightNode();
//             } else {
//                 throw std::out_of_range("Right node ID out of range");
//             }
//         }
//         auto &Ls = leftAdj_[u];
//         if (Ls.find(v) != Ls.end()) return false;
//         Ls.insert(v);
//         rightAdj_[v].insert(u);
//         return true;
//     }
//
//     // 
//     bool removeEdge(int u, int v) {
//         auto &Ls = leftAdj_[u];
//         if (Ls.find(v) == Ls.end()) return false;
//         Ls.erase(v);
//         rightAdj_[v].erase(u);
//         return true;
//     }
//
//     // 
//     bool hasEdge(int u, int v) const {
//         return leftAdj_[u].find(v) != leftAdj_[u].end();
//     }
//
//     // /
//     int numLeft() const  { return nLeft_; }
//     int numRight() const { return nRight_; }
//
//     dense_hash_set<int> &getLeftAdj(int u) {
//         if (u < 0 || u >= nLeft_)
//             throw std::out_of_range("Left node ID out of range");
//         return leftAdj_[u];
//     }
//
//     dense_hash_set<int> &getRightAdj(int v) {
//         if (v < 0 || v >= nRight_)
//             throw std::out_of_range("Right node ID out of range");
//         return rightAdj_[v];
//     }
//
// private:
//     int nLeft_, nRight_;
//     size_t bucketSize_;
//     std::vector<dense_hash_set<int>> leftAdj_, rightAdj_;
//
//     //  dense_hash_set
//     dense_hash_set<int> createEmptySet() const {
//         dense_hash_set<int> s;
//         s.set_empty_key(EMPTY_KEY);
//         s.set_deleted_key(DELETED_KEY);
//         s.resize(bucketSize_ * 2);
//         return s;
//     }
// };
//
// #endif //DYNAMICBIPARTITEGRAPH_HPP
