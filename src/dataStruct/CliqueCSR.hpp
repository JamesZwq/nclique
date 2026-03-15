#pragma once
#ifndef CLIQUE_CSR_HPP
#define CLIQUE_CSR_HPP

#include <vector>
#include <cstdint>
#include <span>
#include "CSR.hpp"

/**
 * CliqueCSR: CSR-based storage for clique results
 * 
 * Stores cliques (sets of vertices) in compressed sparse row format:
 * - offset_: clique start positions in data_
 * - data_: concatenated vertex lists for all cliques
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
