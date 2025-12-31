#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <chrono>
#include <random>
#include <cstdint>
#include <iomanip>
#include <cstring>

// ================== 基础数学工具 ==================
template<typename T>
static inline unsigned __int128 binom_u128(T n, T r) noexcept {
    if (r > n) return 0;
    if (r == 0 || r == n) return 1;
    if (r > n - r) r = n - r;
    unsigned __int128 res = 1;
    for (T i = 1; i <= r; ++i) res = res * (n - r + i) / i;
    return res;
}

// 模拟 Rank 计算
uint64_t mock_rank(const std::vector<uint32_t>& c) {
    unsigned __int128 acc = 0;
    for (size_t i = 0; i < c.size(); ++i) {
        // 强制转换 i+1 为 uint32_t，与 c[i] 匹配
        acc += binom_u128(c[i], (uint32_t)(i + 1));
    }
    return (uint64_t)acc;
}

// Unrank 辅助二分
uint32_t find_max_v(unsigned __int128 rank, uint32_t r, uint32_t max_n) {
    uint32_t low = r, high = max_n, ans = r;
    while (low <= high) {
        uint32_t mid = low + (high - low) / 2;
        if (binom_u128(mid, r) <= rank) { ans = mid; low = mid + 1; }
        else { high = mid - 1; }
    }
    return ans;
}

// 模拟 Unrank 计算
void mock_unrank(uint64_t key, int k, uint32_t maxV, std::vector<uint32_t>& out) {
    out.resize(k);
    unsigned __int128 curr = key;
    uint32_t search_limit = maxV;
    
    // 【修复点】：这里 i 是 int，v 是 uint32_t，需要统一类型
    for (int i = k; i >= 1; --i) {
        // 将 i 转为 uint32_t
        uint32_t r = (uint32_t)i; 
        uint32_t v = find_max_v(curr, r, search_limit);
        out[i-1] = v;
        
        // 【修复点】：显式转换类型
        curr -= binom_u128(v, r); 
        
        if (v > 0) search_limit = v - 1;
    }
}

// ================== 方案 A: Hash Index ==================
class HashIndex {
public:
    int k_;
    std::vector<uint32_t> pool_; 
    std::unordered_map<uint64_t, uint32_t> map_;

    HashIndex(int k) : k_(k) {}

    void build(const std::vector<std::vector<uint32_t>>& cliques) {
        pool_.reserve(cliques.size() * k_);
        map_.reserve(cliques.size());
        uint32_t id = 0;
        for (const auto& c : cliques) {
            uint64_t key = mock_rank(c);
            pool_.insert(pool_.end(), c.begin(), c.end());
            map_[key] = id++;
        }
    }

    uint32_t byClique(const std::vector<uint32_t>& c) {
        uint64_t key = mock_rank(c);
        return map_[key];
    }

    std::vector<uint32_t> byId(uint32_t id) {
        std::vector<uint32_t> res;
        res.reserve(k_);
        size_t start = id * k_;
        for(int i=0; i<k_; ++i) res.push_back(pool_[start + i]);
        return res;
    }
    
    size_t memory_usage() {
        // 粗略估算：Pool + HashMap (Entry overhead + Buckets)
        return pool_.capacity() * 4 + map_.size() * (8 + 4 + 16) + map_.bucket_count() * 8;
    }
};

// ================== 方案 B: Sorted Vector Index ==================
class SortedIndex {
public:
    int k_;
    uint32_t maxV_;
    std::vector<uint64_t> keys_;

    SortedIndex(int k, uint32_t maxV) : k_(k), maxV_(maxV) {}

    void build(const std::vector<std::vector<uint32_t>>& cliques) {
        keys_.reserve(cliques.size());
        for (const auto& c : cliques) {
            keys_.push_back(mock_rank(c));
        }
        std::sort(keys_.begin(), keys_.end());
    }

    uint32_t byClique(const std::vector<uint32_t>& c) {
        uint64_t key = mock_rank(c);
        auto it = std::lower_bound(keys_.begin(), keys_.end(), key);
        return (uint32_t)std::distance(keys_.begin(), it);
    }

    std::vector<uint32_t> byId(uint32_t id) {
        std::vector<uint32_t> res;
        mock_unrank(keys_[id], k_, maxV_, res);
        return res;
    }

    size_t memory_usage() {
        return keys_.capacity() * 8;
    }
};

// ================== 测试工具 ==================
class Timer {
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point start_;
public:
    Timer() : start_(Clock::now()) {}
    double elapsed_ms() {
        auto end = Clock::now();
        return std::chrono::duration<double, std::milli>(end - start_).count();
    }
};

std::vector<std::vector<uint32_t>> generate_data(size_t N, int k, uint32_t maxV) {
    std::vector<std::vector<uint32_t>> data;
    data.reserve(N);
    std::mt19937 gen(42);
    std::uniform_int_distribution<uint32_t> dist(0, maxV);
    
    for(size_t i=0; i<N; ++i) {
        std::vector<uint32_t> c;
        while(c.size() < k) {
            uint32_t v = dist(gen);
            bool exist = false;
            for(auto x : c) if(x==v) exist=true;
            if(!exist) c.push_back(v);
        }
        std::sort(c.begin(), c.end());
        data.push_back(c);
    }
    return data;
}

// ================== Main ==================
int main() {
    int k = 4;
    uint32_t maxV = 1000000;
    
    std::vector<size_t> N_list = {100000, 1000000, 5000000, 10000000}; 

    std::cout << std::left << std::setw(10) << "N" 
              << std::setw(10) << "Method"
              << std::setw(15) << "Build(ms)"
              << std::setw(15) << "Query(ID->C)"
              << std::setw(15) << "Find(C->ID)"
              << std::setw(15) << "Mem(MB)" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    for (size_t N : N_list) {
        auto data = generate_data(N, k, maxV);

        // --- Test Hash ---
        {
            HashIndex idx(k);
            Timer t1;
            idx.build(data);
            double t_build = t1.elapsed_ms();

            Timer t2;
            // 只要跑一部分 id->c，否则 1000w 次太慢
            size_t sample = std::min(N, (size_t)100000);
            for(size_t i=0; i<sample; ++i) { volatile auto v = idx.byId(i); }
            double t_byId = t2.elapsed_ms() * (double)N / sample; 

            Timer t3;
            for(size_t i=0; i<sample; ++i) { volatile auto id = idx.byClique(data[i]); }
            double t_byClique = t3.elapsed_ms() * (double)N / sample;

            std::cout << std::left << std::setw(10) << N 
                      << std::setw(10) << "HASH"
                      << std::setw(15) << t_build
                      << std::setw(15) << t_byId
                      << std::setw(15) << t_byClique
                      << std::setw(15) << (idx.memory_usage() / 1024.0 / 1024.0) << std::endl;
        }

        // --- Test Sorted Vector ---
        {
            SortedIndex idx(k, maxV);
            Timer t1;
            idx.build(data);
            double t_build = t1.elapsed_ms();

            Timer t2;
            size_t sample = std::min(N, (size_t)100000);
            for(size_t i=0; i<sample; ++i) { volatile auto v = idx.byId(i); }
            double t_byId = t2.elapsed_ms() * (double)N / sample;

            Timer t3;
            for(size_t i=0; i<sample; ++i) { volatile auto id = idx.byClique(data[i]); }
            double t_byClique = t3.elapsed_ms() * (double)N / sample;

            std::cout << std::left << std::setw(10) << "" 
                      << std::setw(10) << "SORTED"
                      << std::setw(15) << t_build
                      << std::setw(15) << t_byId
                      << std::setw(15) << t_byClique
                      << std::setw(15) << (idx.memory_usage() / 1024.0 / 1024.0) << std::endl;
        }
        std::cout << std::string(80, '-') << std::endl;
    }

    return 0;
}
