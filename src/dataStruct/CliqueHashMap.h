#pragma once
#include <iterator>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <limits>
#include <iostream>
#include "robin_hood.h"
#include "graph/DynamicGraph.h"

// ===== 小工具：逐项计算 C(n,r)，用 128bit 防溢出 =====
template<typename T>
static inline unsigned __int128 binom_u128(T n, T r) noexcept {
    if (r > n) return 0;
    if (r == 0 || r == n) return 1;
    if (r > n - r) r = n - r;
    unsigned __int128 res = 1;
    for (T i = 1; i <= r; ++i) {
        res = res * (n - r + i) / i;
    }
    return res;
}

// ========== 索引类 ==========
class StaticCliqueIndex {
public:
    using Id = uint32_t; // 可换 size_t

private:
    size_t k_;
    // id -> k_ 连排节点
    std::vector<daf::Size> pool_;
    // rank/hash -> id
    std::vector<robin_hood::unordered_flat_map<uint64_t, Id> > mapList_;
    daf::Size numClique = 0; // 用于统计总 clique 数量
    // ---- rank 计算（组合数排名，保证无碰撞）----
    uint64_t rank64(const daf::Size *data, size_t len) const noexcept {
        unsigned __int128 acc = 0;
        for (size_t i = 0; i < len; ++i) {
            acc += binom_u128(data[i], static_cast<daf::Size>(i + 1));
        }
        return static_cast<uint64_t>(acc);
    }

    uint64_t rank64(const std::vector<daf::Size> &c) const noexcept {
        unsigned __int128 acc = 0;
        for (size_t i = 0; i < c.size(); ++i) {
            acc += binom_u128(c[i], static_cast<daf::Size>(i + 1));
        }
        return static_cast<uint64_t>(acc); // 足够 k<=8 & n<=2^60；否则可返回 128bit 并哈希
    }

    static uint64_t rank64_from_two_sorted(
        const daf::Size *pivots, size_t pCount,
        const TreeGraphNode *keeps, size_t keepCount) noexcept {
        size_t k = pCount + keepCount;
        unsigned __int128 acc = 0;
        size_t i = 0, j = 0;
        // idx 从 1 开始，对应 binom(c[idx-1], idx)
        for (size_t idx = 1; idx <= k; ++idx) {
            daf::Size v;
            if (i < pCount && (j >= keepCount || pivots[i] < keeps[j].v)) {
                v = pivots[i++];
            } else {
                v = keeps[j++].v;
            }
            acc += binom_u128(v, static_cast<daf::Size>(idx));
        }
        return static_cast<uint64_t>(acc);
    }

public:
    StaticCliqueIndex(size_t K) : k_(K) {
    }

    StaticCliqueIndex(const DynamicGraph<TreeGraphNode> &treeGraph, const daf::Size maxV, size_t K) : k_(K) {
        build(treeGraph, maxV);
    }

    // ------- 批量构建 -------
    void build(const DynamicGraph<TreeGraphNode> &treeGraph, const daf::Size maxV) {
        daf::StaticVector<daf::Size> countV = treeGraph.cliqueCountPerVAcc(maxV, k_);
        daf::Size N = treeGraph.cliqueCount(k_);
        mapList_.resize(countV.size());
        for (daf::Size i = 0; i < countV.size(); ++i) {
            if (countV[i] > 0) {
                mapList_[i].reserve(static_cast<size_t>(countV[i] * 1.1));
            }
        }
        std::cout << "Clique Index: " << N << " cliques, " << countV.size() << " vertices." << std::endl;
        pool_.reserve(N * k_);
        numClique = 0;
        daf::StaticVector<daf::Size> keep;
        daf::StaticVector<daf::Size> drop;

        for (const auto &clique: treeGraph.adj_list) {
            // Clique c = raw;
            // std::sort(c.begin(), c.end());          // 保升序，确保 rank 唯一
            // auto key = rank64(c);
            for (const auto &node: clique) {
                if (node.isPivot) drop.push_back(node.v);
                else keep.push_back(node.v);
            }
            daf::enumerateCombinations(keep, drop, k_,
                                       [&](const daf::StaticVector<daf::Size> &keep,
                                           const daf::StaticVector<daf::Size> &combination) {
                                           // 1) 扩容并准备输出区间
                                           size_t start = pool_.size();
                                           pool_.resize(start + k_);
                                           daf::Size *out = pool_.data() + start;

                                           // 2) 线性归并 combination (长度 m) 和 keep (长度 n),  m+n == k_
                                           const daf::Size *a = combination.data;
                                           const daf::Size *a_end = a + combination.size();
                                           const daf::Size *b = keep.data;
                                           const daf::Size *b_end = b + keep.size();

                                           while (a < a_end && b < b_end) {
                                               *out++ = (*a < *b ? *a++ : *b++);
                                           }
                                           if (a < a_end) {
                                               // 剩余 a 段
                                               size_t rem = a_end - a;
                                               std::memcpy(out, a, rem * sizeof(daf::Size));
                                           } else if (b < b_end) {
                                               // 剩余 b 段
                                               size_t rem = b_end - b;
                                               std::memcpy(out, b, rem * sizeof(daf::Size));
                                           }

                                           // 3) 计算组合排名：对刚刚写入的 k_ 个元素做 rank
                                           unsigned __int128 acc = 0;
                                           for (size_t i = 0; i < k_; ++i) {
                                               daf::Size v = pool_[start + i];
                                               acc += binom_u128(v, static_cast<daf::Size>(i + 1));
                                           }
                                           uint64_t key = static_cast<uint64_t>(acc);

                                           // 4) 插入到 map（假设永不失败）
                                           daf::Size vmin = pool_[start]; // merge 后第一位即最小值
                                           mapList_[vmin].emplace(key, numClique++);

                                           return true;
                                       }
            );
            keep.clear();
            drop.clear();
        }

        keep.free();
        drop.free();
        countV.free();
        pool_.shrink_to_fit();
    }

    std::span<const daf::Size> byId(Id id) const {
        size_t start = static_cast<size_t>(id) * k_;
        if (start + k_ > pool_.size()) {
            throw std::out_of_range("bad clique id");
        }
        return {pool_.data() + start, k_};
    }

    template<typename Iterable>
    Id byClique(const Iterable &c) const {
        // Assume c is sorted and has length k_
        const auto r = c.size();
        if (static_cast<size_t>(r) != k_) {
            throw std::invalid_argument("byClique: wrong clique size");
        }
        // First element is minimum
        auto it_begin = std::begin(c);
        daf::Size vmin = *it_begin;
        // Compute rank directly from the iterable
        unsigned __int128 acc = 0;
        size_t idx = 0;
        for (auto it = it_begin; it != std::end(c); ++it, ++idx) {
            acc += binom_u128(static_cast<daf::Size>(*it),
                              static_cast<daf::Size>(idx + 1));
        }
        uint64_t key = static_cast<uint64_t>(acc);
        // Lookup in mapList_
        const auto &bucket = mapList_[vmin];
        auto it = bucket.find(key);
        if (it == bucket.end()) {
            std::cerr << "Error: clique not found for key=" << key
                    << " vmin=" << vmin << "\n";
            throw std::runtime_error("byClique: clique not found");
        }
        return it->second;
    }

    template<typename Iterable>
    std::pair<Id, bool> byNewClique(const Iterable &c) {
        // Assume c is sorted and has length k_
        if (c.size() != k_) {
            throw std::invalid_argument("byNewClique: wrong clique size");
        }

        // 2) 计算最小顶点，用作桶索引
        daf::Size vmin = *std::begin(c);

        // 3) 计算 rank
        unsigned __int128 acc = 0;
        size_t idx = 0;
        for (auto it = std::begin(c); it != std::end(c); ++it, ++idx) {
            acc += binom_u128(static_cast<daf::Size>(*it),
                              static_cast<daf::Size>(idx + 1));
        }
        uint64_t key = static_cast<uint64_t>(acc);

        // 4) 在 mapList_ 中查找
        auto &bucket = mapList_[vmin];
        auto itMap = bucket.find(key);
        if (itMap != bucket.end()) {
            // 已存在，直接返回 ID
            return {itMap->second, false};
        }

        // 5) 不存在：插入新 clique
        Id newId = static_cast<Id>(numClique);
        // 5.1 扩容 pool_: 依次追加 c 中的每一个元素
        pool_.insert(pool_.end(), std::begin(c), std::end(c));

        // 5.2 向桶里插入 (key -> newId)
        bucket.emplace(key, newId);
        ++numClique;

        return {newId, true};
    }

    Id byClique(const daf::Size *pivots, size_t pivotCount,
                const TreeGraphNode *keeps, size_t keepCount) const {
        if (pivotCount == 0 || pivotCount + keepCount != k_) {
            throw std::invalid_argument(
                "byClique: pivotCount must be >=1 and pivotCount+keepCount==k");
        }
        // 取两段首元素的最小值作为桶索引
        daf::Size vmin = keeps && keepCount > 0
                             ? std::min(pivots[0], keeps[0].v)
                             : pivots[0];

        // 一次归并计算 rank
        uint64_t key = rank64_from_two_sorted(
            pivots, pivotCount,
            keeps, keepCount
        );

        // 在桶里查找
        const auto &bucket = mapList_[vmin];
        auto it = bucket.find(key);
        if (it == bucket.end()) {
            throw std::runtime_error("byClique: clique not found");
        }
        return it->second;
    }

    template<
        typename PivotContainer,
        typename KeepContainer,
        // 回退到 C++17/14 时可把下面两行换成 static_assert
        typename = std::enable_if_t<
            std::is_convertible_v<typename PivotContainer::value_type, daf::Size> &&
            std::is_convertible_v<typename KeepContainer::value_type, daf::Size>
        >
    >
    Id byClique(const PivotContainer &pivots,
                const KeepContainer &keeps) const {
        const size_t pCount = pivots.size();
        const size_t kCount = keeps.size();
        if (pCount == 0 || pCount + kCount != k_) {
            throw std::invalid_argument(
                "byClique: pivotCount must be >=1 and pivotCount+keepCount==k");
        }

        // 1) 选出桶索引 vmin = min(pivots.front(), keeps.front())
        daf::Size firstP = static_cast<daf::Size>(*std::begin(pivots));
        daf::Size vmin = firstP;
        if (kCount > 0) {
            daf::Size firstK = static_cast<daf::Size>(*std::begin(keeps));
            vmin = std::min(firstP, firstK);
        }

        // 2) 一次归并同时计算 rank
        unsigned __int128 acc = 0;
        size_t idx = 1; // binom(_, idx)
        auto itP = std::begin(pivots), endP = std::end(pivots);
        auto itK = std::begin(keeps), endK = std::end(keeps);

        while (itP != endP && itK != endK) {
            daf::Size vP = static_cast<daf::Size>(*itP);
            daf::Size vK = static_cast<daf::Size>(*itK);
            if (vP < vK) {
                acc += binom_u128(vP, static_cast<daf::Size>(idx++));
                ++itP;
            } else {
                acc += binom_u128(vK, static_cast<daf::Size>(idx++));
                ++itK;
            }
        }
        // 把剩下的都吃掉
        while (itP != endP) {
            daf::Size vP = static_cast<daf::Size>(*itP++);
            acc += binom_u128(vP, static_cast<daf::Size>(idx++));
        }
        while (itK != endK) {
            daf::Size vK = static_cast<daf::Size>(*itK++);
            acc += binom_u128(vK, static_cast<daf::Size>(idx++));
        }

        uint64_t key = static_cast<uint64_t>(acc);
        const auto &bucket = mapList_[vmin];
        auto found = bucket.find(key);
        if (found == bucket.end()) {
            throw std::runtime_error("byClique: clique not found");
        }
        return found->second;
    }

    auto size() const noexcept { return numClique; }

    void verify() const {
        // 1) 验证 byId / byClique 正向、反向一致
        for (Id id = 0; id < numClique; ++id) {
            // 拿回 CSR 存储的 clique
            auto span_c = byId(id);
            // 转成可排序的 vector
            std::vector<daf::Size> c(span_c.begin(), span_c.end());
            // 最小值检查
            daf::Size vmin = *std::min_element(c.begin(), c.end());
            if (vmin != c.front()) {
                std::cerr << "Error: id=" << id
                        << " pool_ clique not sorted, min=" << vmin
                        << " but front=" << c.front() << "\n";
                throw std::runtime_error("verify failed: unsorted byId");
            }
            // 反查
            Id id2 = byClique(c);
            if (id2 != id) {
                std::cerr << "Error: byClique returned " << id2
                        << " for id=" << id << "\n";
                throw std::runtime_error("verify failed: byClique mismatch");
            }
        }

        // 2) 验证所有 mapList_ 桶中的 key/id
        for (size_t i = 0; i < mapList_.size(); ++i) {
            const auto &bucket = mapList_[i];
            for (auto const &kv: bucket) {
                uint64_t key = kv.first;
                Id id = kv.second;
                if (id >= numClique) {
                    std::cerr << "Error: mapList_[" << i << "] has id="
                            << id << " >= size()\n";
                    throw std::runtime_error("verify failed: id out of range");
                }
                // 拿出 clique，并检查最小值
                auto span_c = byId(id);
                if (span_c.front() != static_cast<daf::Size>(i)) {
                    std::cerr << "Error: mapList_[" << i
                            << "] entry id=" << id
                            << " has min=" << span_c.front() << "\n";
                    throw std::runtime_error("verify failed: bucket index mismatch");
                }
                // 重新计算 rank
                uint64_t key2 = rank64(span_c.data(), span_c.size());
                if (key2 != key) {
                    std::cerr << "Error: mapList_[" << i
                            << "] entry id=" << id
                            << " key mismatch, stored=" << key
                            << " computed=" << key2 << "\n";
                    throw std::runtime_error("verify failed: key mismatch");
                }
            }
        }

        std::cout << "StaticCliqueIndex verification passed: "
                << numClique << " cliques, "
                << mapList_.size() << " buckets.\n";
    }

    // 输出索引内容：key: clique nodes
    friend std::ostream &operator<<(std::ostream &os, const StaticCliqueIndex &idx) {
        for (Id id = 0; id < idx.numClique; ++id) {
            auto span_c = idx.byId(id);
            uint64_t key = idx.rank64(span_c.data(), span_c.size());
            os << key << ":";
            for (auto v: span_c) {
                os << " " << v;
            }
            os << "\n";
        }
        return os;
    }

    // 打印到指定输出流（默认 std::cout）
    void print(std::ostream &os = std::cout) const {
        os << *this;
    }
};
