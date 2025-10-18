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

// ===== ： C(n,r)， 128bit  =====
template<typename T>
static inline unsigned __int128 binom_u128(T n, T r) noexcept {
    if (r > n) return 0;
    if (r == 0 || r == n) return 1;
    if (r > n - r) r = n - r;

    switch (r) {        //  r = 4 
        case 1:  return n;
        case 2:  return (unsigned __int128)n * (n - 1) / 2;
        case 3:  return (unsigned __int128)n * (n - 1) * (n - 2) / 6;
        case 4:  return (unsigned __int128)n * (n - 1) * (n - 2) * (n - 3) / 24;
        case 5:  return (unsigned __int128)n * (n - 1) * (n - 2) * (n - 3) * (n - 4) / 120;
        case 6:  return (unsigned __int128)n * (n - 1) * (n - 2) * (n - 3) * (n - 4) * (n - 5) / 720;
        default: {
            unsigned __int128 res = 1;
            for (T i = 1; i <= r; ++i)
                res = res * (n - r + i) / i;
            return res;
        }
    }
}

// ==========  ==========
class StaticCliqueIndex {
public:
    using Id = uint32_t; //  daf::Size

private:
    daf::Size k_;
    // id -> k_ 
    std::vector<daf::Size> pool_;
    // rank/hash -> id
    std::vector<robin_hood::unordered_flat_map<uint64_t, Id> > mapList_;
    daf::Size numClique = 0; //  clique 
    // ---- rank （，）----
    uint64_t rank64(const daf::Size *data, daf::Size len) const noexcept {
        unsigned __int128 acc = 0;
        for (daf::Size i = 0; i < len; ++i) {
            acc += binom_u128(data[i], static_cast<daf::Size>(i + 1));
        }
        return static_cast<uint64_t>(acc);
    }

    uint64_t rank64(const std::vector<daf::Size> &c) const noexcept {
        unsigned __int128 acc = 0;
        for (daf::Size i = 0; i < c.size(); ++i) {
            acc += binom_u128(c[i], static_cast<daf::Size>(i + 1));
        }
        return static_cast<uint64_t>(acc); //  k<=8 & n<=2^60； 128bit 
    }

    static uint64_t rank64_from_two_sorted(
        const daf::Size *pivots, daf::Size pCount,
        const TreeGraphNode *keeps, daf::Size keepCount) noexcept {
        daf::Size k = pCount + keepCount;
        unsigned __int128 acc = 0;
        daf::Size i = 0, j = 0;
        // idx  1 ， binom(c[idx-1], idx)
        for (daf::Size idx = 1; idx <= k; ++idx) {
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
    StaticCliqueIndex(daf::Size K) : k_(K) {
    }

    StaticCliqueIndex(const DynamicGraph<TreeGraphNode> &treeGraph, const daf::Size maxV, daf::Size K) : k_(K) {
        build(treeGraph, maxV);
    }

    // -------  -------
    void build(const DynamicGraph<TreeGraphNode> &treeGraph, const daf::Size maxV) {
        daf::StaticVector<double> countV = treeGraph.cliqueCountPerVAcc(maxV, k_);
        daf::Size N = treeGraph.cliqueCount(k_);
        mapList_.resize(countV.size());
        for (daf::Size i = 0; i < countV.size(); ++i) {
            if (countV[i] > 0) {
                mapList_[i].reserve(static_cast<daf::Size>(countV[i] * 1.1));
            }
        }
        std::cout << "Clique Index: " << N << " cliques, " << countV.size() << " vertices." << std::endl;
        pool_.reserve(N * k_);
        numClique = 0;
        daf::StaticVector<daf::Size> keep;
        daf::StaticVector<daf::Size> drop;

        for (const auto &clique: treeGraph.adj_list) {
            // Clique c = raw;
            // std::sort(c.begin(), c.end());          // ， rank 
            // auto key = rank64(c);
            for (const auto &node: clique) {
                if (node.isPivot) drop.push_back(node.v);
                else keep.push_back(node.v);
            }
            daf::enumerateCombinations(keep, drop, k_,
                                       [&](const daf::StaticVector<daf::Size> &keep,
                                           const daf::StaticVector<daf::Size> &combination) {
                                           // 1) 
                                           daf::Size start = pool_.size();
                                           pool_.resize(start + k_);
                                           daf::Size *out = pool_.data() + start;

                                           // 2)  combination ( m)  keep ( n),  m+n == k_
                                           const daf::Size *a = combination.data();
                                           const daf::Size *a_end = a + combination.size();
                                           const daf::Size *b = keep.data();
                                           const daf::Size *b_end = b + keep.size();

                                           while (a < a_end && b < b_end) {
                                               *out++ = (*a < *b ? *a++ : *b++);
                                           }
                                           if (a < a_end) {
                                               //  a 
                                               daf::Size rem = a_end - a;
                                               std::memcpy(out, a, rem * sizeof(daf::Size));
                                           } else if (b < b_end) {
                                               //  b 
                                               daf::Size rem = b_end - b;
                                               std::memcpy(out, b, rem * sizeof(daf::Size));
                                           }

                                           // 3) ： k_  rank
                                           unsigned __int128 acc = 0;
                                           for (daf::Size i = 0; i < k_; ++i) {
                                               daf::Size v = pool_[start + i];
                                               acc += binom_u128(v, static_cast<daf::Size>(i + 1));
                                           }
                                           uint64_t key = static_cast<uint64_t>(acc);

                                           // 4)  map（）
                                           daf::Size vmin = pool_[start]; // merge 
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
        daf::Size start = static_cast<daf::Size>(id) * k_;
        if (start + k_ > pool_.size()) {
            throw std::out_of_range("bad clique id");
        }
        return {pool_.data() + start, k_};
    }

    template<typename Iterable>
    Id byClique(const Iterable &c) const {
        // Assume c is sorted and has length k_
        const auto r = c.size();
        if (static_cast<daf::Size>(r) != k_) {
            throw std::invalid_argument("byClique: wrong clique size");
        }
        // First element is minimum
        auto it_begin = std::begin(c);
        daf::Size vmin = *it_begin;
        // Compute rank directly from the iterable
        unsigned __int128 acc = 0;
        daf::Size idx = 0;
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

        // 2) ，
        daf::Size vmin = *std::begin(c);

        // 3)  rank
        unsigned __int128 acc = 0;
        daf::Size idx = 0;
        for (auto it = std::begin(c); it != std::end(c); ++it, ++idx) {
            acc += binom_u128(static_cast<daf::Size>(*it),
                              static_cast<daf::Size>(idx + 1));
        }
        uint64_t key = static_cast<uint64_t>(acc);

        // 4)  mapList_ 
        auto &bucket = mapList_[vmin];
        auto itMap = bucket.find(key);
        if (itMap != bucket.end()) {
            // ， ID
            return {itMap->second, false};
        }

        // 5) ： clique
        Id newId = static_cast<Id>(numClique);
        // 5.1  pool_:  c 
        pool_.insert(pool_.end(), std::begin(c), std::end(c));

        // 5.2  (key -> newId)
        bucket.emplace(key, newId);
        ++numClique;

        return {newId, true};
    }

    Id byClique(const daf::Size *pivots, daf::Size pivotCount,
                const TreeGraphNode *keeps, daf::Size keepCount) const {
        if (pivotCount == 0 || pivotCount + keepCount != k_) {
            throw std::invalid_argument(
                "byClique: pivotCount must be >=1 and pivotCount+keepCount==k");
        }
        // 
        daf::Size vmin = keeps && keepCount > 0
                             ? std::min(pivots[0], keeps[0].v)
                             : pivots[0];

        //  rank
        uint64_t key = rank64_from_two_sorted(
            pivots, pivotCount,
            keeps, keepCount
        );

        // 
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
        //  C++17/14  static_assert
        typename = std::enable_if_t<
            std::is_convertible_v<typename PivotContainer::value_type, daf::Size> &&
            std::is_convertible_v<typename KeepContainer::value_type, daf::Size>
        > >
    Id byClique(const PivotContainer &pivots,
                const KeepContainer &keeps) const {
        const daf::Size pCount = pivots.size();
        const daf::Size kCount = keeps.size();
        if (pCount == 0 || pCount + kCount != k_) {
            throw std::invalid_argument(
                "byClique: pivotCount must be >=1 and pivotCount+keepCount==k");
        }

        // 1)  vmin = min(pivots.front(), keeps.front())
        daf::Size firstP = static_cast<daf::Size>(*std::begin(pivots));
        daf::Size vmin = firstP;
        if (kCount > 0) {
            daf::Size firstK = static_cast<daf::Size>(*std::begin(keeps));
            vmin = std::min(firstP, firstK);
        }

        // 2)  rank
        unsigned __int128 acc = 0;
        daf::Size idx = 1; // binom(_, idx)
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
        // 
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
        // 1)  byId / byClique 、
        for (Id id = 0; id < numClique; ++id) {
            //  CSR  clique
            auto span_c = byId(id);
            //  vector
            std::vector<daf::Size> c(span_c.begin(), span_c.end());
            // 
            daf::Size vmin = *std::min_element(c.begin(), c.end());
            if (vmin != c.front()) {
                std::cerr << "Error: id=" << id
                        << " pool_ clique not sorted, min=" << vmin
                        << " but front=" << c.front() << "\n";
                throw std::runtime_error("verify failed: unsorted byId");
            }
            // 
            Id id2 = byClique(c);
            if (id2 != id) {
                std::cerr << "Error: byClique returned " << id2
                        << " for id=" << id << "\n";
                throw std::runtime_error("verify failed: byClique mismatch");
            }
        }

        // 2)  mapList_  key/id
        for (daf::Size i = 0; i < mapList_.size(); ++i) {
            const auto &bucket = mapList_[i];
            for (auto const &kv: bucket) {
                uint64_t key = kv.first;
                Id id = kv.second;
                if (id >= numClique) {
                    std::cerr << "Error: mapList_[" << i << "] has id="
                            << id << " >= size()\n";
                    throw std::runtime_error("verify failed: id out of range");
                }
                //  clique，
                auto span_c = byId(id);
                if (span_c.front() != static_cast<daf::Size>(i)) {
                    std::cerr << "Error: mapList_[" << i
                            << "] entry id=" << id
                            << " has min=" << span_c.front() << "\n";
                    throw std::runtime_error("verify failed: bucket index mismatch");
                }
                //  rank
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

    // ：key: clique nodes
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

    // （ std::cout）
    void print(std::ostream &os = std::cout) const {
        os << *this;
    }
};
