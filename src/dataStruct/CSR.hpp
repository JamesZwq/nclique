//
// Compact‑Sparse‑Row graph container
// Author: 张文谦   Date: 2025‑08‑05
//
#ifndef CSR_HPP
#define CSR_HPP

#include <vector>
#include <cstdint>
#include <type_traits>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <span>            // C++20

/**
 * A minimal, header‑only CSR adjacency structure.
 *
 *  * IndexType     – unsigned integral type used for vertex ids and
 *                    edge array offsets.  uint8_t / uint16_t / uint32_t /
 *                    uint64_t are all valid.
 *  * The graph is stored as **undirected by default**; if you set
 *    directed==true in the constructor, every (u,v) in the input edge list
 *    will be inserted exactly once.
 *
 *  Interface (all O(1) unless otherwise stated)
 *  -------------------------------------------------------
 *   size()                    – number of vertices
 *   edges()                   – number of directed edges stored
 *   degree(v)                 – out‑degree of v
 *   neighbours(v) / neighbors – std::span with the adjacency of v
 *   neighbour_begin/end       – raw pointer iterators
 *   has_edge(u,v)             – O(deg(u)) check (binary search because rows
 *                                are sorted)
 */
template<typename IndexType = std::uint32_t>
class CSR {
    static_assert(std::is_unsigned_v<IndexType>,
                  "IndexType must be an unsigned integral type.");

    using index_t = IndexType;

    std::vector<index_t> offset_;   // size = n + 1
    std::vector<index_t> adj_;      // concatenated neighbour lists
    bool                 directed_{false};

public:
    /*------------------------------------------------------------
     *  Construction helpers
     *----------------------------------------------------------*/
    CSR() = default;

    /**
     * Build from an edge list.
     *
     *  @param n          number of vertices (ids must be in [0,n) )
     *  @param edges      list of (u,v) pairs
     *  @param directed   whether to keep the edge list as‑is (true) or store
     *                    an undirected graph (false, default)
     *
     *  ***Throws*** std::invalid_argument if an id >= n is encountered or
     *              ids are not strictly less than std::numeric_limits<index_t>::max().
     */
    CSR(index_t n,
        const std::vector<std::pair<index_t,index_t>>& edges,
        bool directed = false)
        : offset_(n + 1, 0),
          directed_(directed)
    {
        // 1) degree counting
        for (auto [u,v] : edges) {
            if (u >= n || v >= n)
                throw std::invalid_argument("vertex id out of range");
            if constexpr (!std::is_same_v<index_t, std::uint64_t>) {
                constexpr auto lim = std::numeric_limits<index_t>::max();
                if (u >= lim || v >= lim)
                    throw std::invalid_argument("vertex id exceeds IndexType range");
            }
            ++offset_[u + 1];
            if (!directed_) ++offset_[v + 1];
        }

        // 2) prefix sum → offsets
        for (index_t i = 1; i <= n; ++i) offset_[i] += offset_[i - 1];

        adj_.resize(offset_.back());

        // 3) scatter
        std::vector<index_t> cursor(offset_.begin(), offset_.end() - 1);
        for (auto [u,v] : edges) {
            adj_[cursor[u]++] = v;
            if (!directed_) adj_[cursor[v]++] = u;
        }

        // 4) sort rows for fast queries
        for (index_t v = 0; v < n; ++v) {
            auto beg = adj_.begin() + offset_[v];
            auto end = adj_.begin() + offset_[v + 1];
            std::sort(beg, end);
        }
    }

    /**
     * Build a *directed clique* from a vertex-presence Bit‑set.
     *
     *  Every bit position represents a possible vertex id; bits that are `0`
     *  are treated as *absent* (holes).  For every ordered pair (u,v) with
     *  u < v and both vertices present, we store a single directed edge
     *  u → v.  (If you need the symmetric direction too, call the undirected
     *  edge‑list constructor with both (u,v) and (v,u).)
     *
     *  Works with:
     *    • boost::dynamic_bitset
     *    • std::bitset
     *    • any custom bitset that provides .size() and .test(pos)
     *
     *  Complexity:  O(N²) in the number of present vertices (clique size),
     *  but still very fast for the small cliques this path is intended for.
     */
    template<class Bitset>
    explicit CSR(const Bitset& bs)
        : offset_(bs.size() + 1, 0),
          directed_(true)
    {
        const index_t N = static_cast<index_t>(bs.size());

        /* 1) prefix degrees (only vertices whose bit == 1 participate) */
        for (index_t u = 0; u < N; ++u) {
            offset_[u + 1] = offset_[u];                // carry‑over
            if (!bs.test(u)) continue;

            index_t deg = 0;
            for (index_t v = u + 1; v < N; ++v)
                if (bs.test(v)) ++deg;

            offset_[u + 1] += deg;                      // u has ‘deg’ out‑edges
        }

        /* 2) fill adjacency array */
        adj_.resize(offset_.back());
        index_t cursor = 0;
        for (index_t u = 0; u < N; ++u) {
            if (!bs.test(u)) continue;
            for (index_t v = u + 1; v < N; ++v)
                if (bs.test(v)) adj_[cursor++] = v;
        }
        // rows are naturally sorted because v increases monotonically
    }

    /*------------------------------------------------------------
     *  Basic accessors
     *----------------------------------------------------------*/
    [[nodiscard]] index_t size()  const noexcept { return offset_.empty() ? 0 : static_cast<index_t>(offset_.size() - 1); }
    [[nodiscard]] index_t edges() const noexcept { return static_cast<index_t>(adj_.size()); }
    [[nodiscard]] bool    directed() const noexcept { return directed_; }

    [[nodiscard]] index_t degree(index_t v) const noexcept {
        return offset_[v + 1] - offset_[v];
    }

    [[nodiscard]] std::span<const index_t> neighbours(index_t v) const noexcept {
        return std::span<const index_t>(&adj_[offset_[v]], degree(v));
    }
    // US spelling alias
    [[nodiscard]] std::span<const index_t> neighbors(index_t v) const noexcept {
        return neighbours(v);
    }

    [[nodiscard]] const index_t* neighbour_begin(index_t v) const noexcept {
        return &adj_[offset_[v]];
    }
    [[nodiscard]] const index_t* neighbour_end(index_t v) const noexcept {
        return &adj_[offset_[v + 1]];
    }

    /*------------------------------------------------------------
     *  Convenience algorithms
     *----------------------------------------------------------*/
    [[nodiscard]] bool has_edge(index_t u, index_t v) const noexcept {
        auto row = neighbours(u);
        return std::binary_search(row.begin(), row.end(), v);
    }

    /**
     * Iterate every neighbour of v, calling F(index_t) -> void/bool.
     * If F returns bool and returns false, iteration stops early.
     */
    template<typename F>
    void for_each_neighbour(index_t v, F&& f) const {
        for (index_t w : neighbours(v))
            if constexpr (std::is_same_v<std::invoke_result_t<F,index_t>, bool>) {
                if(!f(w)) break;
            } else {
                f(w);
            }
    }
};

#endif // CSR_HPP
