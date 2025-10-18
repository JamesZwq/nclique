//
// Created by _ on 2025/8/12.
//

#ifndef PIVOTER_COREDISJOIN_HPP
#define PIVOTER_COREDISJOIN_HPP
#include <vector>
#include <numeric>
#include <cstdint>
#include <utility>
#include <iostream>
#include <unordered_map>
#include <fstream>

#include "Global/Global.h"


class CoreDisJoin {
public:
    using idx_t = daf::Size;

    explicit CoreDisJoin(daf::Size n, daf::Size numK) noexcept {
        codeDisjointSets.resize(numK);
        for (auto &ds: codeDisjointSets) {
            ds.reset(n);
        }
        this->numK = numK;
        currKIndex = 0;
        this->n = n;
        // coreList.resize(numK);
    }

    inline void unite(idx_t a, idx_t b) noexcept {
        // std::cout << "Unite " << a << " and " << b << " at level " << currKIndex << "\n";
        for (int k = currKIndex; k >= 0; --k) {
            if (!codeDisjointSets[k].unite(a, b)) {
                return;
            }
        }
    }

    // Component 26: [26: ([4, 5, 7]) 27: ([4, 6, 7])]
    inline void addK() noexcept {
        currKIndex++;
        if (currKIndex >= codeDisjointSets.size()) {
            codeDisjointSets.emplace_back(n);
        }
    }

    friend std::ostream &operator<<(std::ostream &os, CoreDisJoin &cdj) {
        os << "coreDisJoin(n=" << cdj.n << ", numK=" << cdj.numK << ", currKIndex=" << cdj.currKIndex << ")\n";
        // print by component
        for (size_t i = 0; i <= cdj.currKIndex; ++i) {
            auto &ds = cdj.codeDisjointSets[i];
            // Only print components with size > 1
            std::unordered_map<idx_t, std::vector<idx_t> > components;
            for (idx_t v = 0; v < ds.n(); ++v) {
                idx_t p = ds.find(v);
                components[p].push_back(v);
            }
            os << "DisjointSet(n=" << ds.n() << ", comp_cnt=" << ds.count() << ")\n";
            for (const auto &[rep, members]: components) {
                if (members.size() == 1) continue;
                os << "  Component " << rep << ": [";
                for (size_t i = 0; i < members.size(); ++i) {
                    if (i > 0) os << " ";
                    os << members[i];
                }
                os << "]\n";
            }
        }
        return os;
    }

    void print() {
        std::cout << *this;
    }

    // template <class Indexable>
    // void print(const Indexable &data) {
    //     std::cout << "coreDisJoin(n=" << n << ", numK=" << numK << ", currKIndex=" << currKIndex << ")\n";
    //     // print by component with user-provided labels for component representative
    //     for (size_t i = 0; i < static_cast<size_t>(currKIndex); ++i) {
    //         codeDisjointSets[i].print_with_data(std::cout, data);
    //     }
    // }

    // ：mapper(id) -> 
    template<class Mapper>
    void print(const Mapper &mapper) {
        std::cout << "coreDisJoin(n=" << n << ", numK=" << numK << ", currKIndex=" << currKIndex << ")\n";
        for (size_t i = 0; i <= static_cast<size_t>(currKIndex); ++i) {
            std::cout << "K = " << i << ":\n";
            auto &ds = codeDisjointSets[i];
            std::unordered_map<idx_t, std::vector<idx_t> > components;
            for (idx_t v = 0; v < ds.n(); ++v) {
                idx_t p = ds.find(v);
                components[p].push_back(v);
            }
            std::cout << "DisjointSet(n=" << ds.n() << ", comp_cnt=" << ds.count() << ")\n";
            for (const auto &kv: components) {
                idx_t rep = kv.first;
                const auto &members = kv.second;
                if (members.size() == 1) continue;
                std::cout << "  Component " << rep << ": [";
                for (size_t i = 0; i < members.size(); ++i) {
                    if (i) std::cout << ' ';
                    std::cout << members[i] << ": (" << mapper(members[i]) << ")";
                }
                std::cout << "]\n";
            }
        }
    }

    // 
    template<class Mapper, class PathLike>
    bool print_to_file(const PathLike &filePath, const Mapper &mapper) {
        std::ofstream ofs(filePath);
        if (!ofs) {
            std::cerr << "[CoreDisJoin] Failed to open file: " << filePath << "\n";
            return false;
        }
        ofs << "coreDisJoin(n=" << n << ", numK=" << numK << ", currKIndex=" << currKIndex << ")\n";
        for (size_t i = 0; i <= static_cast<size_t>(currKIndex); ++i) {
            ofs << "K = " << i << ":\n";
            auto &ds = codeDisjointSets[i];
            std::unordered_map<idx_t, std::vector<idx_t> > components;
            for (idx_t v = 0; v < ds.n(); ++v) {
                idx_t p = ds.find(v);
                components[p].push_back(v);
            }
            ofs << "DisjointSet(n=" << ds.n() << ", comp_cnt=" << ds.count() << ")\n";
            for (const auto &kv: components) {
                idx_t rep = kv.first;
                const auto &members = kv.second;
                if (members.size() == 1) continue;
                ofs << "  Component " << rep << ": [";
                for (size_t i = 0; i < members.size(); ++i) {
                    if (i) ofs << ' ';
                    ofs << members[i] << ": (" << mapper(members[i]) << ")";
                }
                ofs << "]\n";
            }
        }
        ofs.flush();
        return true;
    }

    bool validate_top_layer(daf::Size minSize) noexcept {
        if (codeDisjointSets.empty()) return true;
        size_t top = static_cast<size_t>(currKIndex);
        if (top >= codeDisjointSets.size()) top = codeDisjointSets.size() - 1; // 
        return codeDisjointSets[top].validate_min_component_size(static_cast<idx_t>(minSize));
    }


    // private:
    class DisjointSet {
    public:
        explicit DisjointSet(idx_t n = 0) noexcept { reset(n); }

        //  n 
        void reset(const idx_t n) noexcept {
            parent_.resize(n);
            size_.assign(n, 1);
            std::iota(parent_.begin(), parent_.end(), 0);
            comp_cnt_ = n;
            // beCall_.assign(n, false); //  false
        }

        // ，
        void reserve(idx_t n) {
            parent_.reserve(n);
            size_.reserve(n);
            // beCall_.reserve(n);
        }

        // （ + ）
        inline idx_t find(idx_t x) noexcept {
            while (parent_[x] != x) {
                parent_[x] = parent_[parent_[x]]; // path halving
                x = parent_[x];
            }
            return x;
        }

        //  a  b ； false
        inline bool unite(idx_t a, idx_t b) noexcept {
            a = find(a);
            b = find(b);
            // beCall_[a] = true;
            // beCall_[b] = true;
            if (a == b) return false;
            // （）
            if (size_[a] < size_[b]) std::swap(a, b);
            parent_[b] = a;
            size_[a] += size_[b];
            --comp_cnt_;
            return true;
        }

        // 
        inline bool same(idx_t a, idx_t b) noexcept { return find(a) == find(b); }

        //  x 
        inline idx_t size(idx_t x) noexcept { return size_[find(x)]; }

        // 
        inline idx_t count() const noexcept { return comp_cnt_; }

        // 
        inline idx_t n() const noexcept { return static_cast<idx_t>(parent_.size()); }

        friend std::ostream &operator<<(std::ostream &os, DisjointSet &cdj) {
            // print by component
            std::unordered_map<idx_t, std::vector<idx_t> > components;
            for (idx_t v = 0; v < cdj.n(); ++v) {
                idx_t p = cdj.find(v);
                components[p].push_back(v);
            }
            os << "DisjointSet(n=" << cdj.n() << ", comp_cnt=" << cdj.comp_cnt_ << ")\n";
            for (const auto &[rep, members]: components) {
                if (members.size() == 1) continue;
                os << "  Component " << rep << ": [";
                for (size_t i = 0; i < members.size(); ++i) {
                    if (i > 0) os << " ";
                    os << members[i];
                }
                os << "]\n";
            }
            return os;
        }

        template<class Indexable>
        void print_with_data(std::ostream &os, const Indexable &data) {
            // Build components map: representative -> members
            std::unordered_map<idx_t, std::vector<idx_t> > components;
            for (idx_t v = 0; v < n(); ++v) {
                idx_t p = find(v); // path halving allowed here
                components[p].push_back(v);
            }
            os << "DisjointSet(n=" << n() << ", comp_cnt=" << comp_cnt_ << ")\n";
            for (const auto &kv: components) {
                idx_t rep = kv.first;
                const auto &members = kv.second;
                if (members.size() == 1) continue;
                os << "  Component " << rep << " (" << data[rep] << "): [";
                for (size_t i = 0; i < members.size(); ++i) {
                    if (i) os << ' ';
                    os << members[i];
                }
                os << "]\n";
            }
        }

        // ： rep  mapper(rep) 
        template<class Mapper>
        void print_with_mapper(std::ostream &os, const Mapper &mapper) {
            std::unordered_map<idx_t, std::vector<idx_t> > components;
            for (idx_t v = 0; v < n(); ++v) {
                idx_t p = find(v); // 
                components[p].push_back(v);
            }
            os << "DisjointSet(n=" << n() << ", comp_cnt=" << comp_cnt_ << ")\n";
            for (const auto &kv: components) {
                idx_t rep = kv.first;
                const auto &members = kv.second;
                if (members.size() == 1) continue;
                os << "  Component " << rep << ": [";
                for (size_t i = 0; i < members.size(); ++i) {
                    if (i) os << ' ';
                    // os << mapper(members[i]);
                    os << members[i] << ": (" << mapper(members[i]) << ")";
                }
                os << "]\n";
            }
        }


        bool validate_min_component_size(idx_t minSize) noexcept {
            if (n() == 0) return true;
            std::unordered_map<idx_t, idx_t> comp_sz;
            comp_sz.reserve(static_cast<size_t>(n()));
            for (idx_t v = 0; v < n(); ++v) {
                idx_t r = find(v); // ，
                ++comp_sz[r];
            }
            for (const auto &kv: comp_sz) {
                idx_t sz = kv.second;
                if (sz == 1) continue; // 
                if (sz < minSize) {
                    std::cout << "DisjointSet validation failed: "
                            << "component size " << sz << " < " << minSize
                            << " for representative " << kv.first << "\n";
                    return false; // 2..minSize-1 
                }
            }
            return true;
        }

        // private:
        std::vector<idx_t> parent_;
        std::vector<idx_t> size_;
        // std::vector<bool> beCall_;
        idx_t comp_cnt_ = 0;
    };

    std::vector<DisjointSet> codeDisjointSets;
    // std::vector<daf::Size> coreList;
    daf::Size n;
    daf::Size numK; //  k 
    daf::Size currKIndex; //  k 
};

#endif //PIVOTER_COREDISJOIN_HPP
