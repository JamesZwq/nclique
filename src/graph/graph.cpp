//
// Created by _ on 24-7-24.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include "graph/Graph.h"
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


// std::pair<std::vector<daf::Size>::const_iterator, std::vector<daf::Size>::const_iterator>
// Graph::getNbrInter(const daf::Size node_id) const {
//     auto [nbr_start, nbr_end] = getNbr(node_id);
//     return {adj_list.cbegin() + nbr_start, adj_list.cbegin() + nbr_end};
// }


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

std::vector<daf::Size> Graph::sortByDegeneracyOrder() {
    const daf::Size n = getGraphNodeSize();
    if (n == 0) return {};

    /*------------------------------------------------------------
     * 1. 预处理：统计度数并记录最大度
     *-----------------------------------------------------------*/
    std::vector<int> deg(n);          // 当前度数
    int max_deg = 0;
    for (daf::Size v = 0; v < n; ++v) {
        auto [l, r] = getNbr(v);      // 邻接区间 [l, r)
        deg[v] = static_cast<int>(r - l);
        max_deg = std::max(max_deg, deg[v]);
    }

    /*------------------------------------------------------------
     * 2. 建立 bin 数组：bin[d] 保存“度 = d 的顶点段”的起始下标
     *-----------------------------------------------------------*/
    std::vector<int> bin(max_deg + 1, 0);
    for (int d : deg) ++bin[d];               // 先统计每个度出现次数

    // 前缀和 → 起始位
    int start = 0;
    for (int d = 0; d <= max_deg; ++d) {
        int cnt = bin[d];
        bin[d] = start;
        start += cnt;
    }

    /*------------------------------------------------------------
     * 3. vert[pos] = 顶点编号，pos[v] = v 在 vert 中的位置
     *-----------------------------------------------------------*/
    std::vector<int>          vert(n);
    std::vector<int>          pos(n);

    for (daf::Size v = 0; v < n; ++v) {
        int d = deg[v];
        int p = bin[d]++;          // bin[d] 现在是插入点，插完再 +1
        vert[p] = static_cast<int>(v);
        pos[v]  = p;
    }

    // 把 bin 恢复成段首下标（上一步多 +1 了）
    for (int d = max_deg; d > 0; --d) bin[d] = bin[d - 1];
    bin[0] = 0;

    /*------------------------------------------------------------
     * 4. 剥皮主循环：O(|V| + |E|)
     *-----------------------------------------------------------*/
    std::vector<daf::Size> ordered_vertex;
    ordered_vertex.reserve(n);

    int degeneracy = 0;       // 可选：记录图的 k-core 退化度

    for (int i = 0; i < static_cast<int>(n); ++i) {
        int v = vert[i];
        ordered_vertex.push_back(static_cast<daf::Size>(v));

        degeneracy = std::max(degeneracy, deg[v]);

        // 所有度数严格大于 deg[v] 的邻居度数减一，并在 vert 中前移
        for (auto idx = getNbr(v).first; idx < getNbr(v).second; ++idx) {
            int u = static_cast<int>(adj_list[idx]);
            if (deg[u] > deg[v]) {
                int du  = deg[u];
                int pu  = pos[u];        // u 当前在 vert 的位置
                int pw  = bin[du];       // du 段的首元素
                int w   = vert[pw];      // 该段当前首顶点

                if (u != w) {
                    // 交换 u 和 w 的位置，并更新 pos[]
                    std::swap(vert[pu], vert[pw]);
                    pos[u] = pw;
                    pos[w] = pu;
                }
                ++bin[du];               // du 段首右移一格
                --deg[u];                // 度数减一
            }
        }
    }

    /*------------------------------------------------------------
     * 5. 可选：将高退化度顶点排在前面，保持与 BFS 实现一致
     *-----------------------------------------------------------*/
    // std::ranges::reverse(ordered_vertex);

    /*------------------------------------------------------------
     * 6. 调用框架内已有的重排函数
     *-----------------------------------------------------------*/
    sortVertexByGivenOrder(ordered_vertex);

    // 如需使用 degeneracy，可在此返回或存成员变量
    return ordered_vertex;
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
    std::ofstream out("~/_/pivoter/adj_list.txt");
    for (const auto &i: adj_list_offsets) {
        out << i << "\n";
    }
    out.close();
    // degree 分布已经输出到文件中
    std::cout << "The degree distribution has been output to the file: ~/_/pivoter/adj_list.txt" <<
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
    // ~/_/pivoter/adj_list.txt
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
            auto [nbr_begin, nbr_end] = this->getNbr(child->v);
            for (daf::Size i = nbr_begin; i < nbr_end; ++i) {
                const auto leafId = this->adj_list[i];
                candidateNbr[leafId]++;
                if (candidateNbr[leafId] == r) {
                    numNbr++;
                    selectNbr.push_back(leafId);
                }
            }
            rec(child);
            for (daf::Size i = nbr_begin; i < nbr_end; ++i) {
                const auto leafId = this->adj_list[i];
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
        auto [nbr_begin, nbr_end] = this->getNbr(child->v);
        daf::CliqueSize numKeep = 0;
        std::vector<daf::Size> keep;
        for (daf::Size i = nbr_begin; i < nbr_end; ++i) {
            const auto leafId = this->adj_list[i];
            candidateNbr[leafId] = 1;
        }
        rec(child);
        // candidateNbr.clear();
        for (daf::Size i = nbr_begin; i < nbr_end; ++i) {
            const auto leafId = this->adj_list[i];
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
    // ~/_/pivoter/a
    auto out = fopen("~/_/pivoter/b", "w");
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


void Graph::beSingleEdge() {
    // 1) 第一遍：统计每个 u 的“新度数” deg[u] = |{v in N(u) | v>u}|
    std::vector<daf::Size> deg(n, 0);
    daf::Size new_max_deg = 0;
    for (daf::Size u = 0; u < n; ++u) {
        for (daf::Size ei = adj_list_offsets[u]; ei < adj_list_offsets[u+1]; ++ei) {
            daf::Size v = adj_list[ei];
            if (v > u) {
                ++deg[u];
            }
        }
        if (deg[u] > new_max_deg) new_max_deg = deg[u];
    }

    // 2) 构建新的 offsets 前缀和
    daf::StaticVector<daf::Size> new_offsets(n+1);
    new_offsets.c_size = n + 1;
    new_offsets[0] = 0;
    for (daf::Size u = 0; u < n; ++u) {
        new_offsets[u+1] = new_offsets[u] + deg[u];
    }
    daf::Size m_new = new_offsets[n];

    // 3) 分配并填充新的邻接数组；重用 deg 作为写入指针
    daf::StaticVector<daf::Size> new_adj(m_new);
    new_adj.c_size = m_new;
    std::fill(deg.begin(), deg.end(), 0);
    for (daf::Size u = 0; u < n; ++u) {
        daf::Size base = new_offsets[u];
        for (daf::Size ei = adj_list_offsets[u]; ei < adj_list_offsets[u+1]; ++ei) {
            daf::Size v = adj_list[ei];
            if (v > u) {
                new_adj[ base + deg[u] ] = v;
                ++deg[u];
            }
        }
    }

    // 4) 用新数组替换旧数组，并更新 max_degree 与边数 m
    adj_list_offsets.swap(new_offsets);
    adj_list.swap(new_adj);
    max_degree = new_max_deg;
    new_offsets.free();
    new_adj.free();
    // （如果你有成员 m，记得同时更新 m = m_new;）
}

void Graph::initCore() {
    auto *bin = new daf::Size[this->getMaxDegree() + 1];
    auto *pos = new daf::Size[this->getGraphNodeSize() + 1];
    auto *vert = new daf::Size[this->getGraphNodeSize() + 1];
    // auto *coreV = new daf::Size[this->getGraphNodeSize() + 1];
    coreV.resize(this->getGraphNodeSize() + 1);
    for (daf::Size i = 0; i < this->getGraphNodeSize(); ++i) {
        coreV[i] = this->getNbrCount(i);
    }

    std::fill_n(bin, (this->getMaxDegree() + 1), 0);

    for (auto v = 0; v < this->getGraphNodeSize(); ++v) {
        bin[coreV[v]] += 1;
    }

    daf::Size start = 0;

    for (daf::Size d = 0; d <= this->getMaxDegree(); ++d) {
        const daf::Size num = bin[d];
        bin[d] = start;
        start += num;
    }

    for (auto v = 0; v < this->getGraphNodeSize(); ++v) {
        pos[v] = bin[coreV[v]];
        vert[pos[v]] = v;
        bin[coreV[v]] += 1;
    }

    for (daf::Size d = this->getMaxDegree(); d--;) {
        bin[d + 1] = bin[d];
    }
    bin[0] = 0;

    for (daf::Size i = 0; i < this->getGraphNodeSize(); ++i) {
        auto v = vert[i];
        //        std::cout << v << " ";
        std::pair<daf::Size, daf::Size> nbr = this->getNbr(v);
        for (daf::Size j = nbr.first; j < nbr.second; j++) {
            const daf::Size u = this->adj_list[j];
            if (coreV[u] > coreV[v]) {
                const daf::Size du = coreV[u];
                const daf::Size pu = pos[u];
                const daf::Size pw = bin[du];
                const daf::Size w = vert[pw];

                if (u != w) {
                    pos[u] = pw;
                    pos[w] = pu;
                    vert[pu] = w;
                    vert[pw] = u;
                }

                bin[du]++;
                coreV[u]--;
            }
        }
    }

    delete[] bin;
    delete[] pos;
    delete[] vert;
}