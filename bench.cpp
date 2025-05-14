#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
#include <tuple>

// 对单次测试取多次平均，减少抖动
double single_run(size_t N, size_t M, int repeat) {
    std::vector<int> vec(N);
    std::iota(vec.begin(), vec.end(), 0);
    std::unordered_map<int,int> mp;
    mp.reserve(N);
    for (int i = 0; i < (int)N; ++i) mp[i] = i;

    std::mt19937_64 rng(12345);
    std::uniform_int_distribution<int> dist(0, (int)N - 1);
    std::vector<int> queries(M);
    for (auto &q : queries) q = dist(rng);

    volatile int sink = 0;
    double t_bin = 0, t_map = 0;
    for (int r = 0; r < repeat; ++r) {
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int x : queries) {
            auto it = std::lower_bound(vec.begin(), vec.end(), x);
            sink ^= *it;
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        t_bin += std::chrono::duration<double, std::micro>(t1 - t0).count();

        t0 = std::chrono::high_resolution_clock::now();
        for (int x : queries) {
            sink ^= mp[x];
        }
        t1 = std::chrono::high_resolution_clock::now();
        t_map += std::chrono::duration<double, std::micro>(t1 - t0).count();
    }
    return (t_bin/repeat) - (t_map/repeat); // 正：二分慢；负：二分快
}

int main() {
    std::vector<size_t> Ns = { 10, 20, 50, 100, 200, 500, 1000, 5000 };
    std::vector<size_t> Ms = { 1, 10, 100, 500, 1000, 5000 };

    const int REPEAT = 20;
    std::cout << "N\tM\tΔμs (binary - hashmap)\n";
    std::cout << "----------------------------------\n";
    for (auto N : Ns) {
        for (auto M : Ms) {
            if (M > N * 100) continue; // 避免超时
            double diff = single_run(N, M, REPEAT);
            std::cout
              << N << "\t"
              << M << "\t"
              << (diff) << "\n";
        }
    }
    return 0;
}