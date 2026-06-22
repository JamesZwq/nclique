#pragma once
//
// FlatCliques: CSR layout for a set of maximal cliques.
// Replaces vector<vector<daf::Size>> to avoid 100M+ tiny per-clique heap
// allocations (glibc arena fragmentation ~2x bloat + slowness).
// clique i = data[off[i] .. off[i+1]).  operator[] yields a span (no copy).
//
#include "Global/Global.h"
#include <vector>
#include <span>
#include <algorithm>
#include <cstddef>

struct FlatCliques {
    std::vector<daf::Size> data;
    std::vector<size_t> off{0};            // clique i = data[off[i] .. off[i+1])
    size_t size() const { return off.size() - 1; }
    bool empty() const { return off.size() <= 1; }
    std::span<const daf::Size> operator[](size_t i) const {
        return std::span<const daf::Size>(data.data() + off[i], off[i+1] - off[i]);
    }
    void appendSorted(const std::vector<int>& clique) {   // matches the old per-clique sort
        size_t base = data.size();
        for (int v : clique) data.push_back((daf::Size)v);
        std::sort(data.begin() + base, data.end());
        off.push_back(data.size());
    }
    struct Iter { const FlatCliques* fc; size_t i;
        bool operator!=(const Iter& o) const { return i != o.i; }
        void operator++() { ++i; }
        std::span<const daf::Size> operator*() const { return (*fc)[i]; } };
    Iter begin() const { return {this, 0}; }
    Iter end()   const { return {this, size()}; }
    void reserveData(size_t nCliques, size_t nVerts) { off.reserve(nCliques+1); data.reserve(nVerts); }
};
