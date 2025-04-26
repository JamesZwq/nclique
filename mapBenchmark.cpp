#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <random>
#include <chrono>
#include <numeric>    // iota
#include <iomanip>    // fixed, setprecision

int main() {
    // 要测试的 N 值
    const std::vector<size_t> sizes = {
        100,   // 1e6
    };
    // 每个 N 下的查询次数
    const size_t num_queries = 100'00000;

    std::mt19937 rng(2025);
    std::uniform_int_distribution<int> dist;

    std::cout << "      N"
              << "\t" << "lower_bound (ms)"
              << "\t" << "unordered_map (ms)"
              << "\t" << "found_lb"
              << "\t" << "found_um"
              << "\n";

    for (size_t N : sizes) {
        // ——— 1. 构造数据结构（不计时） ———
        std::vector<int> data(N);
        std::iota(data.begin(), data.end(), 0);


        std::unordered_map<int,int> mapdata;
        mapdata.reserve(N);

        for (size_t i = 0; i < N; ++i) {
            mapdata[(int)i] = (int)i;
        }

        // 生成查询序列（半命中半未命中）
        std::vector<int> queries(num_queries);
        for (size_t i = 0; i < num_queries; ++i) {
            if ((i & 1) == 0)       queries[i] = dist(rng) % (int)N;
            else                    queries[i] = (int)N + (dist(rng) % (int)N);
        }

        // ——— 2. 测试 lower_bound ———
        volatile size_t found_lb = 0;  // volatile 防止优化掉循环
        auto t1 = std::chrono::steady_clock::now();
        for (int key : queries) {
            auto it = std::lower_bound(data.begin(), data.end(), key);
            if (it != data.end() && *it == key) ++found_lb;
        }
        auto t2 = std::chrono::steady_clock::now();
        double ms_lb = std::chrono::duration<double, std::milli>(t2 - t1).count();

        // ——— 3. 测试 unordered_map::find ———
        volatile size_t found_um = 0;
        auto t3 = std::chrono::steady_clock::now();
        for (int key : queries) {
            if (mapdata.find(key) != mapdata.end()) ++found_um;
        }
        auto t4 = std::chrono::steady_clock::now();
        double ms_um = std::chrono::duration<double, std::milli>(t4 - t3).count();

        // ——— 4. 打印结果 ———
        std::cout << std::setw(8) << N << "\t"
                  << std::fixed << std::setprecision(3)
                  << std::setw(12) << ms_lb << "\t"
                  << std::setw(15) << ms_um << "\t"
                  << std::setw(8) << found_lb << "\t"
                  << std::setw(8) << found_um
                  << "\n";
    }
    return 0;
}