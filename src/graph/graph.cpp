//
// Created by 张文谦 on 24-7-24.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include "Graph.h"
#include "../Global/Global.h"
#include "tree/MultiBranchTree.h"

Graph::Graph(const std::string &file_path, bool singleEdge) {
    max_degree = 0;
    std::ifstream input_file(file_path);
    if (!input_file.is_open()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        exit(1);
    }

    daf::Size max_node_id = 0;
    std::string line;

    if (std::getline(input_file, line)) {
        std::istringstream iss(line);
        iss >> n;
        max_node_id = n - 1;
        if (n == 0 || n == 1) {
            input_file.clear();
            input_file.seekg(0, std::ios::beg);
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
        }
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
        if (!singleEdge) {
            ++degree[node_id];
            ++degree[nbr_id];
        } else {
            ++degree[std::min(node_id, nbr_id)];
        }
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
        if (!singleEdge) {
            adj_list[adj_list_offsets[node_id] + degree[node_id]++] = nbr_id;
            adj_list[adj_list_offsets[nbr_id] + degree[nbr_id]++] = node_id;
        } else {
            auto min_id = std::min(node_id, nbr_id);
            adj_list[adj_list_offsets[min_id] + degree[min_id]++] = std::max(node_id, nbr_id);
        }
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


std::pair<std::vector<daf::Size>::const_iterator, std::vector<daf::Size>::const_iterator>
Graph::getNbrInter(const daf::Size node_id) const {
    auto [nbr_start, nbr_end] = getNbr(node_id);
    return {adj_list.cbegin() + nbr_start, adj_list.cbegin() + nbr_end};
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

double Graph::getAvgDegree() const {
    return adj_list.size() / (adj_list_offsets.size() - 1.0);
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

// void Graph::sortByDegeneracyOrdering() {
// std::vector<daf::Size> vertex_order;
// kCore::FlatArrayCoreDecomposition(*this, vertex_order);
// std::ranges::reverse(vertex_order);
// std::cout << "vertex_order: " << vertex_order << std::endl;
// sortVertexByGivenOrder(vertex_order);
// }

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

daf::Size Graph::getDegreeCount(daf::Size d) const {
    if (d > max_degree) {
        return 0;
    }
    daf::Size count = 0;
    for (daf::Size i = 0; i < getGraphNodeSize(); ++i) {
        if (getNbrCount(i) > d) {
            count++;
        }
    }
    return count;
}

Graph::Graph(const MultiBranchTree &tree, const daf::StaticVector<TreeNode *> &leafList, daf::Size minK) {
    max_degree = 0;
    n = tree.root->children.size();
    adj_list_offsets.resize(n + 2);
    adj_list_offsets[0] = 0;
    adj_list_offsets[1] = 0;
    // std::vector<daf::Size> degree(n, 0);
    for (const auto leaf: leafList) {
        if (leaf->MaxDeep < minK) continue;
        auto node = leaf;
        while (node->v != ROOTID) {
            adj_list_offsets[node->v + 2]++;
            node = node->parent;
        }
    }
    std::ofstream out("/Users/zhangwenqian/UNSW/pivoter/adj_list.txt");
    for (const auto &i: adj_list_offsets) {
        out << i << "\n";
    }
    out.close();
    // degree 分布已经输出到文件中
    std::cout << "The degree distribution has been output to the file: /Users/zhangwenqian/UNSW/pivoter/adj_list.txt" <<
            std::endl;
    for (daf::Size i = 1; i < n + 2; ++i) {
        max_degree = std::max(max_degree, adj_list_offsets[i]);
        adj_list_offsets[i] = adj_list_offsets[i - 1] + adj_list_offsets[i];
    }
    // std::cout << adj_list_offsets << std::endl;
    // adj_list_offsets[n] = adj_list_offsets[n - 1];

    adj_list.resize(adj_list_offsets.back());

    for (const auto leaf: leafList) {
        if (leaf->MaxDeep < minK) continue;
        auto node = leaf;
        while (node->v != ROOTID) {
            adj_list[adj_list_offsets[node->v + 1]++] = leaf->leafId;
            node = node->parent;
        }
    }

    // std::cout << "adj_list_offsets: " << adj_list_offsets << std::endl;
    // std::cout << "adj_list: " << adj_list << std::endl;
    // print adj List to file
    // /Users/zhangwenqian/UNSW/pivoter/adj_list.txt
}

Graph::Graph(const MultiBranchTree &tree, const Graph &leafGraph, const daf::StaticVector<TreeNode *> &leafList,
             daf::Size r, daf::Size s) {
    if (r < 2) {
        std::cerr << "minK must be greater than 2" << std::endl;
        exit(1);
    }
    daf::Size numNbr = 0;
    daf::StaticVector<daf::Size> candidateNbr(leafList.size());
    daf::StaticVector<daf::Size> selectNbr(leafList.size());
    candidateNbr.c_size = leafList.size();
    n = 0;
    max_degree = 0;
    adj_list_offsets.resize(leafList.size() + 1);
    adj_list.reserve(leafList.size() * 6);
    std::function<void(TreeNode *)> rec = [&](const TreeNode *node) {
        if (node->children.size() == 0) {
            if (n++ != node->leafId) {
                std::cout << "Error node: " << node->leafId << " count: " << n << std::endl;
                exit(1);
            }
            adj_list_offsets[node->leafId + 1] += numNbr + adj_list_offsets[node->leafId] - 1;
            max_degree = std::max(max_degree, numNbr);
            for (auto nbr: selectNbr) {
                if (nbr == node->leafId) continue;
                adj_list.push_back(nbr);
            }
            // candidateNbr.print();
            std::sort(adj_list.begin() + adj_list_offsets[node->leafId], adj_list.end());
            return;
        }
        for (const auto &child: node->children) {
            if (child->MaxDeep < s) continue;
            auto [nbr_begin, nbr_end] = leafGraph.getNbr(child->v);
            for (daf::Size i = nbr_begin; i < nbr_end; ++i) {
                const auto leafId = leafGraph.adj_list[i];
                candidateNbr[leafId]++;
                if (candidateNbr[leafId] == r) {
                    numNbr++;
                    selectNbr.push_back(leafId);
                }
            }
            rec(child);
            for (daf::Size i = nbr_begin; i < nbr_end; ++i) {
                const auto leafId = leafGraph.adj_list[i];
                if (candidateNbr[leafId] == r) {
                    numNbr--;
                    selectNbr.pop_back();
                }
                candidateNbr[leafId]--;
            }
        }
    };
    for (const auto child: tree.root->children) {
        if (child->MaxDeep < s) continue;
        auto [nbr_begin, nbr_end] = leafGraph.getNbr(child->v);
        daf::CliqueSize numKeep = 0;
        std::vector<daf::Size> keep;
        for (daf::Size i = nbr_begin; i < nbr_end; ++i) {
            const auto leafId = leafGraph.adj_list[i];
            candidateNbr[leafId] = 1;
        }
        rec(child);
        // candidateNbr.clear();
        for (daf::Size i = nbr_begin; i < nbr_end; ++i) {
            const auto leafId = leafGraph.adj_list[i];
            candidateNbr[leafId] = 0;
        }
    }
    // tmpAdj_list.free();
    candidateNbr.free();
    selectNbr.free();
}

void Graph::BronKerboschPivotHelp(std::vector<daf::Size> &R,
                                  std::vector<daf::Size> &P,
                                  std::vector<daf::Size> &X,
                                  std::vector<std::vector<daf::Size> > &cliques) {
    // 若 P 和 X 均为空，则 R 是一个极大团
    if (P.empty() && X.empty()) {
        cliques.push_back(R);
        return;
    }

    // 选择枢轴 u：从 P ∪ X 中选一个邻居最多的顶点
    int u = -1;
    int maxDegree = -1;
    std::vector<daf::Size> unionPX = P;
    unionPX.insert(unionPX.end(), X.begin(), X.end());
    for (int vertex: unionPX) {
        auto [nbegin, nend] = getNbr(vertex);
        int degree = nend - nbegin;
        if (degree > maxDegree) {
            maxDegree = degree;
            u = vertex;
        }
    }

    // 利用枢轴 u 计算 P\N(u)
    std::vector<daf::Size> diff;
    auto [nbegin, nend] = getNbr(u);
    std::set_difference(P.begin(), P.end(),
                        // adj_list[u].begin(), adj_list[u].end(),
                        adj_list.begin() + nbegin, adj_list.begin() + nend,
                        std::back_inserter(diff));

    // 对 diff 中的每个顶点 v 进行递归搜索
    for (daf::Size v: diff) {
        std::vector<daf::Size> R_new = R;
        R_new.push_back(v);

        auto [nbegin, nend] = getNbr(v);
        // P_new = P ∩ N(v)
        std::vector<daf::Size> P_new;
        std::set_intersection(P.begin(), P.end(),
                              // adj_list[v].begin(), adj_list[v].end(),
                              adj_list.begin() + nbegin, adj_list.begin() + nend,
                              std::back_inserter(P_new));

        // X_new = X ∩ N(v)
        std::vector<daf::Size> X_new;
        std::set_intersection(X.begin(), X.end(),
                              // adj_list[v].begin(), adj_list[v].end(),
                              adj_list.begin() + nbegin, adj_list.begin() + nend,
                              std::back_inserter(X_new));

        BronKerboschPivotHelp(R_new, P_new, X_new, cliques);

        // 从 P 中移除 v，并加入 X（保持集合的不重复性）
        P.erase(std::remove(P.begin(), P.end(), v), P.end());
        X.push_back(v);
        // 保持 X 有序以确保下一次 set_intersection 正确
        std::sort(X.begin(), X.end());
    }
}

std::vector<std::vector<daf::Size> > Graph::BronKerboschPivot() {
    std::vector<daf::Size> R; // 当前团（起始为空）
    std::vector<daf::Size> P; // 可拓展的顶点，初始为所有顶点
    std::vector<daf::Size> X; // 已处理的顶点集合，初始为空
    std::vector<std::vector<daf::Size> > cliques; // 保存所有极大团

    // 初始化 P 为所有顶点（假定 0,1,2,...,n-1 已经有序）
    for (daf::Size i = 0; i < n; i++) {
        P.push_back(i);
    }

    // 执行递归算法
    BronKerboschPivotHelp(R, P, X, cliques);

    // 输出所有极大团
    std::cout << "所有极大团：" << std::endl;
    // std::cout << cliques << std::endl;
    // /Users/zhangwenqian/UNSW/pivoter/a
    auto out = fopen("/Users/zhangwenqian/UNSW/pivoter/b", "w");
    for (const auto &clique: cliques) {
        for (int vertex: clique) {
            fprintf(out, "%d ", vertex);
        }
        fprintf(out, "\n");
    }

    return cliques;
}


// cover tree to graph
Graph::Graph(const MultiBranchTree &tree,
             const daf::StaticVector<TreeNode *> &leafList,
             const Graph &treeGraphV,
             daf::Size minK) {
    // treeGraphV.getAvgDegree();
    adj_list_offsets.resize(treeGraphV.adj_list_offsets.size());
    adj_list.reserve(treeGraphV.adj_list.size());
}
