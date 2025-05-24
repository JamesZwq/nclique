/**********************************************************************
 * benchmark_sets_vector_google.cpp   –   C++20 micro-benchmark
 *
 *  比较 6 类容器在 4 种操作上的性能：
 *      1) std::unordered_set
 *      2) absl::flat_hash_set         (Abseil)
 *      3) tsl::robin_set              (tsl-robin-map)
 *      4) folly::F14FastSet           (Facebook Folly)
 *      5) google::dense_hash_set      (Google SparseHash)
 *      6) 已排序 std::vector          (二分插入 / 查找 / 删除；双指针交集)
 *
 *  • 每类容器重复 R 次，输出平均耗时（毫秒）。
 *  • 可用宏 -DNO_ABSL / -DNO_TSL / -DNO_FOLLY / -DNO_GOOGLE 屏蔽某库。
 *  • dense_hash_set 需先 set_empty_key 及 set_deleted_key；随机数据避开这两个键。
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

/* ---------- 计时器 ---------- */
using Clock = std::chrono::high_resolution_clock;
struct Timer {
    Clock::time_point t0;
    void   tic() { t0 = Clock::now(); }
    double toc() const {
        return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    }
};

/* ---------- 生成随机 uint64，排除 UINT64_MAX 和 UINT64_MAX-1 ---------- */
std::vector<std::uint64_t> makeRandom(std::size_t N, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    // 上界设为 max()-2，以避开 empty_key = max() 和 deleted_key = max()-1
    std::uniform_int_distribution<std::uint64_t> dist(
        0, std::numeric_limits<std::uint64_t>::max() - 2);
    std::vector<std::uint64_t> v(N);
    for (auto& x : v) x = dist(rng);
    return v;
}

/* ======================================================================
 *  Sorted Vector “Set” 适配器
 * ====================================================================*/
struct SortedVector {
    std::vector<std::uint64_t> v;

    void reserve(std::size_t n) { v.reserve(n); }

    void insert(std::uint64_t x) {
        auto it = std::lower_bound(v.begin(), v.end(), x);
        if (it == v.end() || *it != x) v.insert(it, x);
    }
    bool contains(std::uint64_t x) const {
        return std::binary_search(v.begin(), v.end(), x);
    }
    void erase(std::uint64_t x) {
        auto it = std::lower_bound(v.begin(), v.end(), x);
        if (it != v.end() && *it == x) v.erase(it);
    }

    std::size_t size() const { return v.size(); }
    auto begin()  const { return v.begin(); }
    auto end()    const { return v.end();   }
};

std::size_t intersectionSize(const SortedVector& A, const SortedVector& B) {
    const auto& a = A.v; const auto& b = B.v;
    std::size_t i = 0, j = 0, cnt = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j]) ++i;
        else if (b[j] < a[i]) ++j;
        else { ++cnt; ++i; ++j; }
    }
    return cnt;
}

/* ======================================================================
 *  Google dense_hash_set 适配器（已添加 deleted_key）
 * ====================================================================*/
#ifndef NO_GOOGLE
struct DenseHashSet {
    using SetT = google::dense_hash_set<std::uint64_t>;
    SetT s;

    DenseHashSet() {
        // 选择两个不在数据域内的 sentinel 值
        auto empty_key   = std::numeric_limits<std::uint64_t>::max();     // 2^64-1
        auto deleted_key = empty_key - 1;                                 // 2^64-2

        s.set_empty_key(empty_key);
        s.set_deleted_key(deleted_key);
    }

    /* 区间构造 */
    template<class It>
    DenseHashSet(It first, It last) : DenseHashSet() {
        std::size_t n = std::distance(first, last);
        s.resize(n * 2);  // 负载 < 0.5
        for (auto it = first; it != last; ++it) s.insert(*it);
    }

    /* 常规接口 */
    void reserve(std::size_t n)          { s.resize(n * 2); }
    void insert(std::uint64_t x)         { s.insert(x); }
    bool contains(std::uint64_t x) const { return s.find(x) != s.end(); }
    void erase(std::uint64_t x)          { s.erase(x); }

    std::size_t size() const { return s.size(); }
    auto begin() const { return s.begin(); }
    auto end()   const { return s.end();   }
};
#endif

/* ======================================================================
 *  通用基准：Insert + Lookup + Erase
 * ====================================================================*/
template<class Set>
double benchInsertLookupErase(const std::vector<std::uint64_t>& data,
                              std::size_t lookups)
{
    Set s;
    s.reserve(data.size());

    Timer tm; tm.tic();
    for (auto x : data) s.insert(x);

    std::size_t hit = 0;
    for (std::size_t i = 0; i < lookups; ++i)
        hit += s.contains(data[i % data.size()]);

    for (auto x : data) s.erase(x);
    (void)hit;
    return tm.toc();
}

/* vector 专用重载 */
double benchInsertLookupErase(const std::vector<std::uint64_t>& data,
                              std::size_t lookups,
                              SortedVector*)
{
    SortedVector s; s.reserve(data.size());
    Timer tm; tm.tic();
    for (auto x : data) s.insert(x);
    std::size_t hit = 0;
    for (std::size_t i = 0; i < lookups; ++i)
        hit += s.contains(data[i % data.size()]);
    for (auto x : data) s.erase(x);
    (void)hit;
    return tm.toc();
}

/* ======================================================================
 *  交集基准
 * ====================================================================*/
template<class Set>
double benchIntersection(const std::vector<std::uint64_t>& A,
                         const std::vector<std::uint64_t>& B)
{
    Set s1(A.begin(), A.end()), s2(B.begin(), B.end());
    const Set& small = (s1.size() < s2.size()) ? s1 : s2;
    const Set& large = (s1.size() < s2.size()) ? s2 : s1;

    Timer tm; tm.tic();
    std::size_t cnt = 0;
    for (auto x : small) if (large.contains(x)) ++cnt;
    (void)cnt;
    return tm.toc();
}

/* vector 专用 */
double benchIntersection(const std::vector<std::uint64_t>& A,
                         const std::vector<std::uint64_t>& B,
                         SortedVector*)
{
    SortedVector s1, s2;
    s1.v.assign(A.begin(), A.end());
    s2.v.assign(B.begin(), B.end());
    std::sort(s1.v.begin(), s1.v.end());
    std::sort(s2.v.begin(), s2.v.end());

    Timer tm; tm.tic();
    auto sz = intersectionSize(s1, s2);
    (void)sz;
    return tm.toc();
}

/* ======================================================================
 *  统一调用（重复 R 次取平均）
 * ====================================================================*/
template<class Set>
void runCase(const std::string& name,
             std::size_t N,
             int R = 5,
             std::size_t w = 13)
{
    std::vector<double> tIns(R), tInt(R);
    for (int r = 0; r < R; ++r) {
        auto d1 = makeRandom(N, 100 + r);
        auto d2 = makeRandom(N, 200 + r);

        if constexpr (std::is_same_v<Set, SortedVector>) {
            tIns[r] = benchInsertLookupErase(d1, N, static_cast<Set*>(nullptr));
            tInt[r] = benchIntersection    (d1, d2, static_cast<Set*>(nullptr));
        } else {
            tIns[r] = benchInsertLookupErase<Set>(d1, N);
            tInt[r] = benchIntersection    <Set>(d1, d2);
        }
    }
    auto avg = [](const std::vector<double>& v){
        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    };

    std::cout << std::left << std::setw(w) << name
              << " | Insert+Lookup+Erase: "
              << std::setw(8) << std::fixed << std::setprecision(3) << avg(tIns) << " ms"
              << " | Intersection: "
              << avg(tInt) << " ms\n";
}

/* ======================================================================*/
int main()
{
    const std::size_t N = 1'000'000;  // 元素数
    const int R = 5;                  // 重复次数

    std::cout << "Bench N = " << N << "  (avg of " << R << " runs)\n";

    runCase<std::unordered_set<std::uint64_t>>("unordered", N, R);
#ifndef NO_ABSL
    runCase<absl::flat_hash_set<std::uint64_t>>("absl_flat", N, R);
#endif
#ifndef NO_TSL
    runCase<tsl::robin_set<std::uint64_t>>("tsl_robin", N, R);
#endif
#ifndef NO_FOLLY
    runCase<folly::F14FastSet<std::uint64_t>>("folly_F14", N, R);
#endif
#ifndef NO_GOOGLE
    runCase<DenseHashSet>("google_dense", N, R);
#endif
    runCase<SortedVector>("vec_sorted", N, R);

    return 0;
}