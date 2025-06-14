#include <vector>
#include <cassert>
#include <algorithm>
#include <chrono>
#include <functional>
#include <iterator>

// 泛用：从随机访问容器 C 中选 size 个元素
// cb(const std::vector<T>& comb)  返回 false 可提前终止
template<typename Container, typename Fn>
bool enumComb(const Container &C, size_t size, Fn cb) {
    size_t n = C.size();
    if (size > n) return true;
    std::vector<typename Container::value_type> comb(size);
    // dfs(pos, left): 从 C[pos..] 还要选 left 个
    std::function<bool(size_t,size_t)> dfs = [&](size_t pos, size_t left) -> bool {
        if (left == 0) {
            return cb(comb);
        }
        // 剪枝：剩余元素不足
        if (n - pos < left) return true;
        // 枚举：在 [pos .. n-left] 中取一个
        for (size_t i = pos; i + left <= n; ++i) {
            comb[size - left] = C[i];
            if (!dfs(i + 1, left - 1)) return false;
        }
        return true;
    };
    return dfs(0, size);
}

// 主函数：从 S ∪ R 中选 k 个，且至少含 1 个来自 S
// S, R 都必须是随机访问容器，value_type 相同
template<typename S, typename R, typename Fn>
bool enumAtLeastOneFromTwo(const S &s,
                            const R &r,
                            size_t k,
                            Fn cb)
{
    size_t sz = s.size(), zr = r.size();
    assert(k > 0 && k <= sz + zr);

    // 一次性分配好 merged buffer
    std::vector<typename S::value_type> merged(k);

    // t = 1..min(k, sz) 从 s 中选 t，从 r 中选 k-t
    for (size_t t = 1; t <= std::min(k, sz); ++t) {
        bool contS = enumComb(s, t, [&](const std::vector<typename S::value_type>& c1){
            // 枚举 R 中补齐 k-t
            return enumComb(r, k - t, [&](const std::vector<typename R::value_type>& c2){
                // 一次两段 memcpy
                std::copy_n(c1.data(), c1.size(), merged.data());
                std::copy_n(c2.data(), c2.size(), merged.data() + c1.size());
                // merged 已升序（两段各自升序且分区不重叠）
                return cb(merged);
            });
        });
        if (!contS) return false;  // S 层提前终止
    }
    return true;
}

#include <iostream>
int main(){
    std::vector<int> S = {1,2,3,4};
    std::vector<int> R = {5,6,7,8,9};
    int k = 3;

    const int TRIALS = 1;
    auto start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < TRIALS; ++i){
        enumAtLeastOneFromTwo(S, R, k, [&](const std::vector<int>& comb){
            std::cout << "Combination: ";
            for (const auto& x : comb) {
                std::cout << x << ' ';
            }
            std::cout << '\n';
            return true;
        });
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Ran " << TRIALS << " trials in " << ms << " ms\n";
    return 0;
}