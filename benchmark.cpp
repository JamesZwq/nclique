/**********************************************************************
 * benchmark_sets_vector_google.cpp   –   C++20 micro-benchmark
 *
 *  比较 6 类容器在 5 种操作上的性能：
 *      1) std::unordered_set
 *      2) absl::flat_hash_set         (Abseil)
 *      3) tsl::robin_set              (tsl-robin-map)
 *      4) folly::F14FastSet           (Facebook Folly)
 *      5) google::dense_hash_set      (Google SparseHash)
 *      6) 已排序 std::vector          (二分插入 / 查找 / 删除；双指针交集)
 *
 *  • 每类容器重复 R 次，输出 3 种操作的总耗时（毫秒）。
 *    参数：N=1000，R=1000
 *  • 操作：Insert+Lookup+Erase、Intersection、Iteration。
 *  • 修正：Intersection 中的循环结果累加到 volatile 变量，防止编译器优化掉。
 *  • 可用宏 -DNO_ABSL / -DNO_TSL / -DNO_FOLLY / -DNO_GOOGLE 屏蔽。
 *********************************************************************/

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#ifndef NO_ABSL
  #include <absl/container/flat_hash_set.h>
#endif
#ifndef NO_TSL
  #include <tsl/robin_set.h>
#endif
#ifndef NO_FOLLY
  #include <folly/container/F14Set.h>
#endif
#ifndef NO_GOOGLE
  #include <google/dense_hash_set>
#endif

using Clock = std::chrono::high_resolution_clock;
struct Timer {
    Clock::time_point t0;
    void   tic() { t0 = Clock::now(); }
    double toc() const {
        return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    }
};

std::vector<std::uint64_t> makeRandom(std::size_t N, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::uint64_t> dist(0, std::numeric_limits<std::uint64_t>::max() - 2);
    std::vector<std::uint64_t> v(N);
    for (auto& x : v) x = dist(rng);
    return v;
}

struct SortedVector {
    std::vector<std::uint64_t> v;
    void reserve(std::size_t n) { v.reserve(n); }
    void insert(std::uint64_t x) { auto it = std::lower_bound(v.begin(), v.end(), x);
                                   if (it == v.end() || *it != x) v.insert(it, x);
    }
    bool contains(std::uint64_t x) const { return std::binary_search(v.begin(), v.end(), x); }
    void erase(std::uint64_t x) { auto it = std::lower_bound(v.begin(), v.end(), x);
                                  if (it != v.end() && *it == x) v.erase(it);
    }
    std::size_t size() const { return v.size(); }
    auto begin()  const { return v.begin(); }
    auto end()    const { return v.end();   }
};

std::size_t intersectionSize(const SortedVector& A, const SortedVector& B) {
    const auto& a = A.v, &b = B.v;
    std::size_t i = 0, j = 0, cnt = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j])      ++i;
        else if (b[j] < a[i]) ++j;
        else { ++cnt; ++i; ++j; }
    }
    return cnt;
}

#ifndef NO_GOOGLE
struct DenseHashSet {
    using SetT = google::dense_hash_set<std::uint64_t>;
    SetT s;
    DenseHashSet() {
        auto empty_key   = std::numeric_limits<std::uint64_t>::max();
        auto deleted_key = empty_key - 1;
        s.set_empty_key(empty_key);
        s.set_deleted_key(deleted_key);
    }
    template<class It>
    DenseHashSet(It first, It last) : DenseHashSet() {
        s.resize(std::distance(first, last) * 2);
        for (auto it = first; it != last; ++it) s.insert(*it);
    }
    void reserve(std::size_t n)          { s.resize(n * 2); }
    void insert(std::uint64_t x)         { s.insert(x); }
    bool contains(std::uint64_t x) const { return s.find(x) != s.end(); }
    void erase(std::uint64_t x)          { s.erase(x); }
    std::size_t size() const { return s.size(); }
    auto begin() const { return s.begin(); }
    auto end()   const { return s.end();   }
};
#endif

// Insert + Lookup + Erase 基准
template<class Set>
double benchInsertLookupErase(const std::vector<std::uint64_t>& data, std::size_t lookups) {
    Set s; s.reserve(data.size());
    Timer t; t.tic();
    for (auto x: data) s.insert(x);
    for (std::size_t i = 0; i < lookups; ++i) s.contains(data[i % data.size()]);
    for (auto x: data) s.erase(x);
    return t.toc();
}
double benchInsertLookupErase(const std::vector<std::uint64_t>& data, std::size_t lookups, SortedVector*) {
    SortedVector s; s.reserve(data.size());
    Timer t; t.tic();
    for (auto x: data) s.insert(x);
    for (std::size_t i = 0; i < lookups; ++i) s.contains(data[i % data.size()]);
    for (auto x: data) s.erase(x);
    return t.toc();
}

// Intersection 基准，累加结果到 volatile 阻止优化
template<class Set>
double benchIntersection(const std::vector<std::uint64_t>& A, const std::vector<std::uint64_t>& B) {
    Set s1(A.begin(), A.end()), s2(B.begin(), B.end());
    const Set& small = (s1.size() < s2.size() ? s1 : s2);
    const Set& large = (s1.size() < s2.size() ? s2 : s1);
    volatile std::size_t cnt = 0;
    Timer t; t.tic();
    for (auto x: small) if (large.contains(x)) ++cnt;
    (void)cnt;
    return t.toc();
}
double benchIntersection(const std::vector<std::uint64_t>& A, const std::vector<std::uint64_t>& B, SortedVector*) {
    SortedVector s1, s2; s1.v = A; s2.v = B;
    std::sort(s1.v.begin(), s1.v.end());
    std::sort(s2.v.begin(), s2.v.end());
    volatile std::size_t cnt = 0;
    Timer t; t.tic();
    cnt = intersectionSize(s1, s2);
    (void)cnt;
    return t.toc();
}

// Iteration 基准
template<class Set>
double benchIteration(const std::vector<std::uint64_t>& data) {
    Set s(data.begin(), data.end()); volatile std::uint64_t sum = 0;
    Timer t; t.tic(); for (auto x: s) sum += x;
    (void)sum;
    return t.toc();
}
double benchIteration(const std::vector<std::uint64_t>& data, SortedVector*) {
    SortedVector s; s.v = data;
    std::sort(s.v.begin(), s.v.end()); volatile std::uint64_t sum = 0;
    Timer t; t.tic(); for (auto x: s) sum += x;
    return t.toc();
}

// 统一调用：总耗时输出
template<class Set>
void runCase(const std::string& name, std::size_t N, int R = 1000, std::size_t w = 14) {
    std::vector<double> tIns(R), tInt(R), tIt(R);
    for (int r = 0; r < R; ++r) {
        auto d1 = makeRandom(N, 100 + r);
        auto d2 = makeRandom(N, 200 + r);
        if constexpr(std::is_same_v<Set, SortedVector>) {
            tIns[r] = benchInsertLookupErase(d1, N, static_cast<Set*>(nullptr));
            tInt [r] = benchIntersection       (d1, d2, static_cast<Set*>(nullptr));
            tIt  [r] = benchIteration          (d1, static_cast<Set*>(nullptr));
        } else {
            tIns[r] = benchInsertLookupErase<Set>(d1, N);
            tInt [r] = benchIntersection    <Set>(d1, d2);
            tIt  [r] = benchIteration       <Set>(d1);
        }
    }
    auto sum = [&](auto& v){ return std::accumulate(v.begin(), v.end(), 0.0); };
    std::cout << std::left << std::setw(w) << name
              << " | Sum(I+L+E): " << std::setw(10) << sum(tIns) << "ms"
              << " | Sum(Int): "    << std::setw(10) << sum(tInt) << "ms"
              << " | Sum(It): "     << sum(tIt)  << "ms\n";
}

int main() {
    const std::size_t N = 1000;
    const int R = 10000;
    std::cout << "Bench N = " << N << "  (total of " << R << " runs)\n";
    runCase<std::unordered_set<std::uint64_t>>("unordered", N);
#ifndef NO_ABSL
    runCase<absl::flat_hash_set<std::uint64_t>>("absl_flat", N);
#endif
#ifndef NO_TSL
    runCase<tsl::robin_set<std::uint64_t>>("tsl_robin", N);
#endif
#ifndef NO_FOLLY
    runCase<folly::F14FastSet<std::uint64_t>>("folly_F14", N);
#endif
#ifndef NO_GOOGLE
    runCase<DenseHashSet>("google_dense", N);
#endif
    runCase<SortedVector>("vec_sorted", N);
    return 0;
}
