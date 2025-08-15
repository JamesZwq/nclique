// filename: compare_csr_adj.cpp
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <vector>
#include <unordered_set>
using namespace std;
using namespace std::chrono;

// 随机生成无向图的边列表
vector<pair<int,int>> make_edges(int n, int m) {
    mt19937_64 rng(12345);
    uniform_int_distribution<int> dist(0, n-1);
    unordered_set<unsigned long long> seen;
    vector<pair<int,int>> edges;
    edges.reserve(m);
    while ((int)edges.size() < m) {
        int u = dist(rng), v = dist(rng);
        if (u == v) continue;
        if (u > v) swap(u,v);
        unsigned long long code = (unsigned long long)u<<32 | (unsigned long long)v;
        if (seen.insert(code).second) {
            edges.emplace_back(u,v);
        }
    }
    return edges;
}

// 构造 vector<vector<int>> 邻接表
vector<vector<int>> build_adj(int n, const vector<pair<int,int>>& edges) {
    vector<vector<int>> adj(n);
    for (auto &e: edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    return adj;
}

// 构造 CSR (offsets + neighbors)
void build_csr(int n,
               const vector<pair<int,int>>& edges,
               vector<int>& offsets,
               vector<int>& nbrs)
{
    offsets.assign(n+1, 0);
    for (auto &e: edges) {
        offsets[e.first+1]++;
        offsets[e.second+1]++;
    }
    for (int i = 1; i <= n; i++)
        offsets[i] += offsets[i-1];

    nbrs.resize(offsets[n]);
    vector<int> pos = offsets;
    for (auto &e: edges) {
        nbrs[pos[e.first]++] = e.second;
        nbrs[pos[e.second]++] = e.first;
    }
}

int main(int argc, char** argv){
    int n = 10000, m = 50000, q = 10000000;
    if (argc >= 4) {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        q = atoi(argv[3]);
    }
    cout << "Nodes="<<n<<" Edges="<<m<<" Queries="<<q<<"\n";

    // 1) 生成随机图
    auto edges = make_edges(n, m);

    // 2) 构造两种结构
    auto adj = build_adj(n, edges);
    vector<int> csr_off, csr_nbr;
    build_csr(n, edges, csr_off, csr_nbr);

    // 3) 两个相同种子、同步生成随机序列的 RNG
    mt19937 rng_list(123), rng_csr(123);
    uniform_int_distribution<int> dist_node(0, n-1);

    // 4) 直接在循环里生成随机 (u,i) 访问，避免大内存
    //    List-of-lists
    auto t1 = high_resolution_clock::now();
    long long sum1 = 0;
    for (int k = 0; k < q; ++k) {
        int u = dist_node(rng_list);
        int deg = adj[u].size();
        int i = deg ? (int)(rng_list() % deg) : 0;
        sum1 += adj[u][i];
    }
    auto t2 = high_resolution_clock::now();

    //    CSR
    auto t3 = high_resolution_clock::now();
    long long sum2 = 0;
    for (int k = 0; k < q; ++k) {
        int u = dist_node(rng_csr);
        int deg = adj[u].size();
        int i = deg ? (int)(rng_csr() % deg) : 0;
        sum2 += csr_nbr[ csr_off[u] + i ];
    }
    auto t4 = high_resolution_clock::now();

    double tl = duration<double>(t2 - t1).count();
    double tc = duration<double>(t4 - t3).count();

    cout << fixed << setprecision(6);
    cout << "List-of-lists: " << tl << " s, checksum=" << sum1 << "\n";
    cout << "CSR         : " << tc << " s, checksum=" << sum2 << "\n";
    cout << "Speedup CSR / list = " << (tl / tc) << "x\n";

    return 0;
}