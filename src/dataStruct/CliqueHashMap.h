#pragma once

#include <vector>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <type_traits>

// ... (保留之前的辅助 struct has_v_member 不变) ...
template<typename T, typename = void>
struct has_v_member : std::false_type {};
template<typename T>
struct has_v_member<T, std::void_t<decltype(std::declval<T>().v)>> : std::true_type {};

// 辅助 Traits: 检测是否有 .size() 方法
template<typename T, typename = void>
struct has_size_method : std::false_type {};

template<typename T>
struct has_size_method<T, std::void_t<decltype(std::declval<T>().size())>> : std::true_type {};
class StaticCliqueIndex {
public:
    using Id = uint32_t;
    using Size = uint32_t;

private:
    Size k_;
    Size maxV_ = 0;
    std::vector<uint64_t> keys_;
    bool is_finalized_ = false;

    // 【新优化】组合数查找表 C[r][n]
    // 扁平化存储：index = (r-1) * (maxV_ + 1) + n
    // 使用 uint128 存储以防溢出
    std::vector<unsigned __int128> C_table_;
    bool use_lookup_table_ = false;

    // ------------------------------------------------
    // 基础计算 (Fallback)
    // ------------------------------------------------
    static inline unsigned __int128 binom_u128_calc(Size n, Size r) noexcept {
        if (r > n) return 0;
        if (r == 0 || r == n) return 1;
        if (r > n - r) r = n - r;
        unsigned __int128 res = 1;
        for (Size i = 1; i <= r; ++i) res = res * (n - r + i) / i;
        return res;
    }

    // ------------------------------------------------
    // 查表计算 (Fast Path)
    // ------------------------------------------------
    inline unsigned __int128 get_binom_fast(Size n, Size r) const {
        // r 从 1 开始，存放在 r-1 行
        // 表的大小是 maxV_ + 1
        return C_table_[(size_t)(r - 1) * (maxV_ + 1) + n];
    }

    // ------------------------------------------------
    // Rank (统一接口)
    // ------------------------------------------------
    uint64_t rank64(const std::vector<Size> &c) const noexcept {
        unsigned __int128 acc = 0;
        if (use_lookup_table_) {
            for (size_t i = 0; i < c.size(); ++i) {
                // c[i] 是 n, (i+1) 是 r
                acc += get_binom_fast(c[i], (Size)(i + 1));
            }
        } else {
            for (size_t i = 0; i < c.size(); ++i) {
                acc += binom_u128_calc(c[i], (Size)(i + 1));
            }
        }
        return static_cast<uint64_t>(acc);
    }

    // ------------------------------------------------
    // Unrank (优化重点)
    // ------------------------------------------------
    void unrank64(uint64_t key, std::vector<Size>& out_nodes) const {
        out_nodes.resize(k_);
        unsigned __int128 current_rank = key;
        Size search_limit = maxV_;

        for (int i = k_; i >= 1; --i) {
            Size r = (Size)i;
            Size v;

            if (use_lookup_table_) {
                // 【极速版】在预计算的数组中二分查找
                // 这一行的起始位置
                const auto* row_start = &C_table_[(size_t)(r - 1) * (maxV_ + 1)];
                // 我们要在 [r, search_limit] 范围内找最大的 v 使得 C(v, r) <= rank
                // 等价于找第一个 > rank 的位置，然后减 1

                // std::upper_bound 在有序数组上极其快
                const auto* it = std::upper_bound(row_start + r, row_start + search_limit + 1, current_rank);
                v = (Size)(std::distance(row_start, it) - 1);

                current_rank -= *(row_start + v); // 直接减去查表值
            } else {
                // 【慢速版】
                v = find_max_v_fallback(current_rank, r, search_limit);
                current_rank -= binom_u128_calc(v, r);
            }

            out_nodes[i - 1] = v;
            if (v > 0) search_limit = v - 1;
        }
    }

    Size find_max_v_fallback(unsigned __int128 rank, Size r, Size max_n) const {
        Size low = r, high = max_n, ans = r;
        while (low <= high) {
            Size mid = low + (high - low) / 2;
            if (binom_u128_calc(mid, r) <= rank) { ans = mid; low = mid + 1; }
            else { high = mid - 1; }
        }
        return ans;
    }

public:
    StaticCliqueIndex(Size K) : k_(K) {}

    // 预分配
    void reserve(size_t estimated_count) { keys_.reserve(estimated_count); }

    // 插入
    template<typename Iterable>
    void insert(const Iterable& c) {
        if (is_finalized_) throw std::runtime_error("Cannot insert after finish()");
        // 临时计算 rank，这里还没启用查表，因为表需要 maxV 才能建
        // 建议 insert 阶段只用 calc
        std::vector<Size> temp; temp.reserve(k_);
        for(auto v : c) temp.push_back(static_cast<Size>(v));
        std::sort(temp.begin(), temp.end());

        // 临时使用 calc 计算 rank，因为此时还没 finish，不知道 maxV
        unsigned __int128 acc = 0;
        for (size_t i = 0; i < temp.size(); ++i)
            acc += binom_u128_calc(temp[i], (Size)(i + 1));

        keys_.push_back(static_cast<uint64_t>(acc));
    }

    // 兼容接口
    // =========================================================
    // 通用插入接口：支持任意可迭代容器 (Vector, Set, Span, List...)
    // =========================================================
    template<typename PivotIterable, typename KeepIterable>
    void insert(const PivotIterable& pivots, const KeepIterable& keeps) {
        if (is_finalized_) throw std::runtime_error("Cannot insert after finish()");

        std::vector<Size> c;
        // 预估大小，避免多次分配 (假设 pivots 和 keeps 都有 size() 方法)
        // 如果容器没有 size() (比如 forward_list)，这行可能需要 SFINAE 保护，
        // 但对于常见图算法容器 (vector/set) 都是安全的。
        if constexpr (has_size_method<PivotIterable>::value && has_size_method<KeepIterable>::value) {
            c.reserve(pivots.size() + keeps.size());
        } else {
            c.reserve(k_);
        }

        // 1. 处理 Pivots (假设都是整数类型)
        for (const auto& v : pivots) {
            c.push_back(static_cast<Size>(v));
        }

        // 2. 处理 Keeps (自动适配 .v 或 整数)
        for (const auto& v : keeps) {
            // 获取 v 的实际类型 (去除引用和 const)
            using ElementType = std::decay_t<decltype(v)>;

            if constexpr (has_v_member<ElementType>::value) {
                // 如果是结构体 (如 TreeGraphNode)，取 .v
                c.push_back(static_cast<Size>(v.v));
            } else {
                // 如果是普通整数
                c.push_back(static_cast<Size>(v));
            }
        }

        // 3. 排序 (Rank 计算的前提)
        std::sort(c.begin(), c.end());

        // 4. 计算 Rank Key
        // 注意：这里必须用 calc 版本，因为此时 Lookup Table 可能还未构建
        unsigned __int128 acc = 0;
        for (size_t i = 0; i < c.size(); ++i) {
            acc += binom_u128_calc(c[i], (Size)(i + 1));
        }
        keys_.push_back(static_cast<uint64_t>(acc));
    }

    // 【核心】完成构建并初始化查找表
    void finish(Size maxV) {
        if (is_finalized_) return;
        maxV_ = maxV;

        // 1. 排序 Keys
        std::sort(keys_.begin(), keys_.end());
        keys_.shrink_to_fit();

        // 2. 【新】构建 Lookup Table
        // 内存消耗：k * maxV * 16 bytes.
        // 例如 k=5, maxV=300k => 24MB. 非常小！
        try {
            size_t table_size = (size_t)k_ * (maxV_ + 1);
            C_table_.resize(table_size);

            // 填充表格
            for (Size r = 1; r <= k_; ++r) {
                for (Size n = 0; n <= maxV_; ++n) {
                    C_table_[(size_t)(r - 1) * (maxV_ + 1) + n] = binom_u128_calc(n, r);
                }
            }
            use_lookup_table_ = true;

            double mb = (double)table_size * 16.0 / 1024.0 / 1024.0;
            std::cout << "Lookup Table Built. Mem: " << mb << " MB" << std::endl;

        } catch (const std::bad_alloc& e) {
            std::cerr << "Warning: Failed to alloc lookup table. Fallback to calculation." << std::endl;
            use_lookup_table_ = false;
        }

        is_finalized_ = true;

        double indexMB = (double)keys_.size() * 8.0 / 1024.0 / 1024.0;
        std::cout << "Clique Index Finalized: " << keys_.size() << " cliques. Index Mem: " << indexMB << " MB" << std::endl;
    }

    // ------- 查询阶段 -------

    std::vector<Size> byId(Id id) const {
        if (id >= keys_.size()) throw std::out_of_range("byId id out of range");
        uint64_t key = keys_[id];
        std::vector<Size> nodes;
        unrank64(key, nodes); // 内部会自动使用 Fast Path
        return nodes;
    }

    template<typename Iterable>
    Id byClique(const Iterable &c) const {
        // ... (同前，略，内部调用 rank64 会自动用 Fast Path) ...
        std::vector<Size> sorted_c; sorted_c.reserve(k_);
        for(auto v : c) sorted_c.push_back(static_cast<Size>(v));
        std::sort(sorted_c.begin(), sorted_c.end());

        uint64_t key = rank64(sorted_c);
        auto it = std::lower_bound(keys_.begin(), keys_.end(), key);
        if (it == keys_.end() || *it != key) throw std::runtime_error("clique not found");
        return static_cast<Id>(std::distance(keys_.begin(), it));
    }

    // ... (保留其他 byClique 重载，逻辑不变) ...
    // 指针接口
    template<typename T_Pivot, typename T_Keep>
    Id byClique(const T_Pivot *pivots, Size pivotCount,
                const T_Keep *keeps, Size keepCount) const {
        std::vector<Size> c; c.reserve(k_);
        for(size_t i=0; i<pivotCount; ++i) c.push_back(static_cast<Size>(pivots[i]));
        for(size_t i=0; i<keepCount; ++i) {
            if constexpr (has_v_member<T_Keep>::value) c.push_back(static_cast<Size>(keeps[i].v));
            else c.push_back(static_cast<Size>(keeps[i]));
        }
        return byClique(c);
    }

    template<typename PivotCont, typename KeepCont>
    Id byClique(const PivotCont &pivots, const KeepCont &keeps) const {
        std::vector<Size> c; c.reserve(k_);
        for(auto v : pivots) c.push_back(static_cast<Size>(v));
        for(auto v : keeps) {
             if constexpr (has_v_member<typename KeepCont::value_type>::value) c.push_back(static_cast<Size>(v.v));
             else c.push_back(static_cast<Size>(v));
        }
        return byClique(c);
    }

    size_t size() const noexcept { return keys_.size(); }
    void print(std::ostream &os = std::cout) const { os << "Index Size: " << keys_.size() << "\n"; }
};