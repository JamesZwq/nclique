// filename: compare_value_vs_pointer.cpp
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

struct EdgeData {
    size_t nbr;
    double weight;
};

static vector<pair<int,int>> make_edges(int n, int m) {
    mt19937_64 rng(12345);
    uniform_int_distribution<int> dist(0, n-1);
    unordered_set<unsigned long long> seen;
    vector<pair<int,int>> edges;
    edges.reserve(m);
    while ((int)edges.size() < m) {
        int u = dist(rng), v = dist(rng);
        if (u == v) continue;
        if (u > v) swap(u,v);
        unsigned long long code = ((unsigned long long)u<<32) | (unsigned long long)v;
        if (seen.insert(code).second)
            edges.emplace_back(u,v);
    }
    return edges;
}

int main(int argc, char** argv){
    int n = 100000;         // 节点数
    int m = 200000;         // 边数
    long long q = 20000000; // 随机查询次数
    if (argc >= 4) {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        q = atoll(argv[3]);
    }
    cout << "Nodes="<<n<<" Edges="<<m<<" Queries="<<q<<"\n";

    // 1) 随机图
    auto edges = make_edges(n,m);

    // 2) 构造 value 版邻接表
    vector<vector<EdgeData>> adj_val(n);
    adj_val.reserve(n);
    for (auto &e: edges) {
        double w = (double)(e.first ^ e.second) / m;  // 任意 weight
        adj_val[e.first].push_back({(size_t)e.second, w});
        adj_val[e.second].push_back({(size_t)e.first, w});
    }

    // 3) 构造 pointer 版邻接表（用一个大池子分配所有 EdgeData）
    vector<EdgeData> pool;
    pool.reserve(2*m);
    vector<vector<EdgeData*>> adj_ptr(n);
    adj_ptr.reserve(n);
    for (auto &e: edges) {
        double w = (double)(e.first ^ e.second) / m;
        pool.push_back({(size_t)e.second, w});
        adj_ptr[e.first].push_back(&pool.back());
        pool.push_back({(size_t)e.first, w});
        adj_ptr[e.second].push_back(&pool.back());
    }

    // 4) 随机访问 benchmark
    mt19937_64 rng(2024);
    uniform_int_distribution<size_t> dist_node(0, n-1);

    // value
    auto t0 = high_resolution_clock::now();
    long long sum_val = 0;
    for (long long i = 0; i < q; ++i) {
        size_t u = dist_node(rng);
        auto &L = adj_val[u];
        size_t d = L.size();
        if (d) {
            size_t idx = (size_t)(rng() % d);
            sum_val += L[idx].nbr;  // 读 nbr 字段
        }
    }
    auto t1 = high_resolution_clock::now();

    // pointer
    // 重置 RNG 以保证访问序列完全一样
    rng.seed(2024);
    auto t2 = high_resolution_clock::now();
    long long sum_ptr = 0;
    for (long long i = 0; i < q; ++i) {
        size_t u = dist_node(rng);
        auto &L = adj_ptr[u];
        size_t d = L.size();
        if (d) {
            size_t idx = (size_t)(rng() % d);
            sum_ptr += L[idx]->nbr;  // 额外一次指针解引用
        }
    }
    auto t3 = high_resolution_clock::now();

    double tv = duration<double>(t1-t0).count();
    double tp = duration<double>(t3-t2).count();
    cout << fixed << setprecision(6);
    cout << "  Value  adj: " << tv << " s, checksum="<<sum_val<<"\n";
    cout << "Pointer adj: " << tp << " s, checksum="<<sum_ptr<<"\n";
    cout << "Speedup Value/Pointer = " << (tv/tp) << "x\n";

    return 0;
}