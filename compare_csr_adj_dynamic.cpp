// filename: compare_csr_adj_dynamic.cpp
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

// 随机生成无向图的边列表（允许重复，节省内存）
vector<pair<int,int>> make_edges(int n, int m) {
    mt19937_64 rng(12345);
    uniform_int_distribution<int> dist(0, n-1);
    vector<pair<int,int>> edges;
    edges.reserve(m);
    while ((int)edges.size() < m) {
        int u = dist(rng), v = dist(rng);
        if (u != v) edges.emplace_back(u, v);
    }
    return edges;
}

// 构造 vector<vector<int>> 邻接表并测试随机访问
double test_adj(int n,
                const vector<pair<int,int>>& edges,
                long long q,
                long long &checksum)
{
    // build adj
    vector<vector<int>> adj(n);
    adj.reserve(n);
    for (auto &e : edges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
    }
    // random access
    mt19937 rng(123);
    uniform_int_distribution<int> dist_node(0, n-1);
    auto t0 = high_resolution_clock::now();
    checksum = 0;
    for (long long i = 0; i < q; ++i) {
        int u = dist_node(rng);
        auto &neigh = adj[u];
        int d = neigh.size();
        if (d) checksum += neigh[rng() % d];
    }
    auto t1 = high_resolution_clock::now();
    return duration<double>(t1 - t0).count();
}

// 构造 CSR 并测试随机访问
double test_csr(int n,
                const vector<pair<int,int>>& edges,
                long long q,
                long long &checksum)
{
    // build csr
    vector<int> offsets(n+1);
    for (auto &e: edges) {
        offsets[e.first + 1]++;
        offsets[e.second + 1]++;
    }
    for (int i = 1; i <= n; ++i) offsets[i] += offsets[i-1];
    vector<int> nbrs(offsets[n]);
    vector<int> pos = offsets;
    for (auto &e: edges) {
        nbrs[pos[e.first]++] = e.second;
        nbrs[pos[e.second]++] = e.first;
    }
    // free edges? they remain outside

    // random access
    mt19937 rng(123);
    uniform_int_distribution<int> dist_node(0, n-1);
    auto t0 = high_resolution_clock::now();
    checksum = 0;
    for (long long i = 0; i < q; ++i) {
        int u = dist_node(rng);
        int d = offsets[u+1] - offsets[u];
        if (d) checksum += nbrs[ offsets[u] + (rng() % d) ];
    }
    auto t1 = high_resolution_clock::now();
    return duration<double>(t1 - t0).count();
}

int main(int argc, char** argv) {
    int n = 100000, m = 500000;
    long long q = 10000000;
    if (argc >= 4) {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        q = atoll(argv[3]);
    }
    cout << "Nodes="<<n<<" Edges="<<m<<" Queries="<<q<<"\n";

    // generate
    auto edges = make_edges(n, m);

    long long sum1=0, sum2=0;
    double t_list = test_adj(n, edges, q, sum1);
    // free adj by scope exit inside function

    double t_csr  = test_csr(n, edges, q, sum2);

    cout << fixed << setprecision(6);
    cout << "List-of-lists: " << t_list << " s, checksum="<<sum1<<"\n";
    cout << "CSR         : " << t_csr  << " s, checksum="<<sum2<<"\n";
    cout << "Speedup CSR / list = " << (t_list / t_csr) << "x\n";
    return 0;
}