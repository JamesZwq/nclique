#pragma once
#ifndef CLIQUE_CSR_HPP
#define CLIQUE_CSR_HPP

#include <vector>
#include <cstdint>
#include <span>
#include <cstring>
#include "CSR.hpp"

extern double nCr[1001][401];

/**
 * CliqueCSR: CSR-based storage for clique results
 *
 * Stores cliques (sets of vertices) in compressed sparse row format:
 * - offset_: clique start positions in data_
 * - data_: concatenated vertex lists for all cliques
 * - pivot_: per-vertex isPivot flags (optional, for cliqueCount)
 *
 * This is much more efficient than std::vector<std::vector<T>> for:
 * - Memory layout (single contiguous allocation)
 * - Cache locality
 * - Fast merge/flush operations (just memcpy)
 */
template<typename NodeType = std::uint32_t>
class CliqueCSR {
    using index_t = NodeType;

    std::vector<index_t> offset_;  // size = num_cliques + 1
    std::vector<index_t> data_;    // concatenated clique vertex lists
    std::vector<uint8_t> pivot_;   // isPivot flag per vertex (parallel to data_)

public:
    CliqueCSR() {
        offset_.push_back(0);  // Initialize with starting offset
    }
    
    explicit CliqueCSR(size_t num_cliques) {
        offset_.reserve(num_cliques + 1);
        offset_.push_back(0);
    }
    
    /**
     * Add a clique (vertex list)
     */
    void add_clique(const std::vector<index_t>& vertices) {
        for (index_t v : vertices) {
            data_.push_back(v);
        }
        offset_.push_back(static_cast<index_t>(data_.size()));
    }
    
    /**
     * Add a clique from raw pointer
     */
    void add_clique(const index_t* vertices, size_t count) {
        for (size_t i = 0; i < count; ++i) {
            data_.push_back(vertices[i]);
        }
        offset_.push_back(static_cast<index_t>(data_.size()));
    }
    
    /**
     * Get number of cliques
     */
    [[nodiscard]] size_t num_cliques() const noexcept {
        return offset_.empty() ? 0 : offset_.size() - 1;
    }
    
    /**
     * Get clique i as a span
     */
    [[nodiscard]] std::span<const index_t> clique(size_t i) const noexcept {
        return std::span<const index_t>(&data_[offset_[i]], offset_[i+1] - offset_[i]);
    }
    
    /**
     * Get total number of vertices across all cliques
     */
    [[nodiscard]] size_t total_vertices() const noexcept {
        return data_.size();
    }
    
    /**
     * Reserve space for cliques
     */
    void reserve_cliques(size_t n) {
        offset_.reserve(n + 1);
    }
    
    /**
     * Reserve space for vertex data
     */
    void reserve_vertices(size_t n) {
        data_.reserve(n);
    }

    /**
     * Bulk init from pre-built flat arrays (zero-copy move).
     * offsets: size = total_cliques+1, offsets[i] = start of clique i in data.
     * data: concatenated vertex IDs.
     */
    void init_from_flat(std::vector<index_t>&& offsets_in, std::vector<index_t>&& data_in) {
        offset_ = std::move(offsets_in);
        data_   = std::move(data_in);
    }

    /**
     * Bulk init with pivot flags (zero-copy move).
     */
    void init_from_flat(std::vector<index_t>&& offsets_in, std::vector<index_t>&& data_in,
                        std::vector<uint8_t>&& pivot_in) {
        offset_ = std::move(offsets_in);
        data_   = std::move(data_in);
        pivot_  = std::move(pivot_in);
    }

    /**
     * Whether pivot flags are available.
     */
    [[nodiscard]] bool has_pivot() const noexcept { return !pivot_.empty(); }

    /**
     * Compute clique counts per size k.
     * Returns vector where result[k] = number of k-cliques encoded in this SDCT.
     * Requires pivot flags to be present.
     */
    [[nodiscard]] std::vector<double> cliqueCount() const {
        size_t maxSz = 0;
        size_t nc = num_cliques();
        for (size_t i = 0; i < nc; i++) {
            size_t sz = (size_t)(offset_[i+1] - offset_[i]);
            if (sz > maxSz) maxSz = sz;
        }
        std::vector<double> counts(maxSz + 1, 0.0);
        for (size_t i = 0; i < nc; i++) {
            size_t beg = (size_t)offset_[i], end = (size_t)offset_[i+1];
            int pivotCount = 0, nonPivotCount = 0;
            for (size_t j = beg; j < end; j++) {
                if (pivot_[j]) pivotCount++;
                else nonPivotCount++;
            }
            int rsize = pivotCount + nonPivotCount;
            for (int q = 0; q <= pivotCount; q++) {
                int k = rsize - q;
                counts[(size_t)k] += nCr[pivotCount][q];
            }
        }
        return counts;
    }

    /**
     * Clear all cliques
     */
    void clear() {
        offset_.clear();
        data_.clear();
        offset_.push_back(0);
    }
    
    /**
     * Get raw offset array (for advanced usage)
     */
    [[nodiscard]] const std::vector<index_t>& offsets() const noexcept {
        return offset_;
    }
    
    /**
     * Get raw data array (for advanced usage)
     */
    [[nodiscard]] const std::vector<index_t>& vertices() const noexcept {
        return data_;
    }
};

#endif // CLIQUE_CSR_HPP
