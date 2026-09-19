#pragma once

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace partial {
using Vertex = uint32_t;
using Count = uint64_t;
using Vertices = std::vector<Vertex>;
constexpr Vertex absent = std::numeric_limits<Vertex>::max();
constexpr Count infinity = std::numeric_limits<Count>::max();
using Clock = std::chrono::steady_clock;
inline double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
inline void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
inline void checked_add(Count& a, Count b) {
    if (b >= infinity - a) throw std::overflow_error("64-bit count overflow");
    a += b;
}

struct Graph {
    Vertex n = 0;
    size_t m = 0;
    std::vector<size_t> off;
    Vertices neighbors;
    std::span<const Vertex> row(Vertex v) const {
        return std::span<const Vertex>(neighbors).subspan(off[v], off[v + 1] - off[v]);
    }
    bool adjacent(Vertex u, Vertex v) const {
        auto r = row(u);
        return std::binary_search(r.begin(), r.end(), v);
    }
    static Graph from_edges(Vertex n, std::vector<std::pair<Vertex, Vertex>> edges) {
        for (auto& [u, v] : edges) {
            require(u < n && v < n && u != v, "invalid simple-graph edge");
            if (v < u) std::swap(u, v);
        }
        std::sort(edges.begin(), edges.end());
        require(std::adjacent_find(edges.begin(), edges.end()) == edges.end(), "duplicate edge");
        Graph g;
        g.n = n;
        g.m = edges.size();
        g.off.assign(static_cast<size_t>(n) + 1, 0);
        for (auto [u, v] : edges) { ++g.off[u + 1]; ++g.off[v + 1]; }
        std::partial_sum(g.off.begin(), g.off.end(), g.off.begin());
        g.neighbors.resize(2 * g.m);
        auto cursor = g.off;
        for (auto [u, v] : edges) {
            g.neighbors[cursor[u]++] = v;
            g.neighbors[cursor[v]++] = u;
        }
        for (Vertex v = 0; v < n; ++v)
            std::sort(g.neighbors.begin() + static_cast<ptrdiff_t>(g.off[v]),
                      g.neighbors.begin() + static_cast<ptrdiff_t>(g.off[v + 1]));
        return g;
    }
    static Graph read(const std::string& path) {
        std::ifstream in(path);
        require(bool(in), "cannot open graph");
        std::string line;
        size_t nn = 0, mm = 0;
        bool header = false;
        std::vector<std::pair<Vertex, Vertex>> edges;
        while (std::getline(in, line)) {
            const auto pos = line.find_first_not_of(" \t\r");
            if (pos == std::string::npos || line[pos] == '#') continue;
            std::istringstream stream(line);
            uint64_t a, b;
            require(bool(stream >> a >> b), "invalid graph line");
            if (!header) {
                require(a < absent, "too many vertices");
                nn = a; mm = b; header = true; edges.reserve(mm);
            } else {
                require(a < nn && b < nn, "edge endpoint out of range");
                edges.emplace_back(static_cast<Vertex>(a), static_cast<Vertex>(b));
            }
        }
        require(header && edges.size() == mm, "graph header count differs from edges");
        return from_edges(static_cast<Vertex>(nn), std::move(edges));
    }
};

inline Vertices intersect(const Vertices& candidates, std::span<const Vertex> neighbors) {
    Vertices result;
    result.reserve(std::min(candidates.size(), neighbors.size()));
    if (neighbors.size() > 8 * candidates.size()) {
        for (Vertex v : candidates)
            if (std::binary_search(neighbors.begin(), neighbors.end(), v)) result.push_back(v);
    } else {
        std::set_intersection(candidates.begin(), candidates.end(), neighbors.begin(), neighbors.end(),
                              std::back_inserter(result));
    }
    return result;
}
inline size_t intersection_size(const Vertices& candidates, std::span<const Vertex> neighbors) {
    size_t count = 0;
    if (neighbors.size() > 8 * candidates.size()) {
        for (Vertex v : candidates)
            count += std::binary_search(neighbors.begin(), neighbors.end(), v);
    } else {
        size_t i = 0, j = 0;
        while (i < candidates.size() && j < neighbors.size()) {
            if (candidates[i] < neighbors[j]) ++i;
            else if (neighbors[j] < candidates[i]) ++j;
            else { ++count; ++i; ++j; }
        }
    }
    return count;
}

class Heap {
    Vertices heap_, position_;
    std::vector<Count>& keys_;
    bool less(Vertex a, Vertex b) const {
        return keys_[a] < keys_[b] || (keys_[a] == keys_[b] && a < b);
    }
    void swap_at(size_t a, size_t b) {
        std::swap(heap_[a], heap_[b]);
        position_[heap_[a]] = static_cast<Vertex>(a);
        position_[heap_[b]] = static_cast<Vertex>(b);
    }
    void down(size_t at) {
        while (2 * at + 1 < heap_.size()) {
            size_t child = 2 * at + 1;
            if (child + 1 < heap_.size() && less(heap_[child + 1], heap_[child])) ++child;
            if (!less(heap_[child], heap_[at])) break;
            swap_at(at, child); at = child;
        }
    }
public:
    explicit Heap(std::vector<Count>& keys): heap_(keys.size()), position_(keys.size()), keys_(keys) {
        std::iota(heap_.begin(), heap_.end(), 0);
        std::iota(position_.begin(), position_.end(), 0);
        for (size_t i = heap_.size() / 2; i; --i) down(i - 1);
    }
    bool empty() const { return heap_.empty(); }
    Vertex first() const { return heap_.front(); }
    bool contains(Vertex v) const { return position_[v] != absent; }
    Vertex pop() {
        Vertex v = heap_.front();
        swap_at(0, heap_.size() - 1);
        heap_.pop_back(); position_[v] = absent;
        if (!heap_.empty()) down(0);
        return v;
    }
    void decrease(Vertex v, Count value) {
        require(contains(v) && value <= keys_[v], "invalid heap decrease");
        keys_[v] = value;
        size_t at = position_[v];
        while (at && less(heap_[at], heap_[(at - 1) / 2])) {
            size_t parent = (at - 1) / 2;
            swap_at(at, parent); at = parent;
        }
    }
};

struct Seeds {
    Vertices order, rank, ordinary;
    Vertex maximum = 0;
    explicit Seeds(const Graph& g): rank(g.n), ordinary(g.n) {
        std::vector<Count> degree(g.n);
        for (Vertex v = 0; v < g.n; ++v) degree[v] = g.row(v).size();
        Heap heap(degree);
        order.reserve(g.n);
        while (!heap.empty()) {
            Vertex v = heap.pop();
            maximum = std::max(maximum, static_cast<Vertex>(degree[v]));
            ordinary[v] = maximum;
            rank[v] = static_cast<Vertex>(order.size()); order.push_back(v);
            for (Vertex u : g.row(v)) if (heap.contains(u)) heap.decrease(u, degree[u] - 1);
        }
    }
};

class Combinations {
    int columns_;
    std::vector<Count> values_;
public:
    Combinations(Vertex maximum, int s): columns_(s + 1),
        values_((static_cast<size_t>(maximum) + 1) * static_cast<size_t>(s + 1), 0) {
        const size_t columns = static_cast<size_t>(columns_);
        for (Vertex n = 0; n <= maximum; ++n) {
            values_[static_cast<size_t>(n) * columns] = 1;
            for (int k = 1; k <= s && static_cast<Vertex>(k) <= n; ++k) {
                const size_t col = static_cast<size_t>(k);
                Count a = values_[static_cast<size_t>(n - 1) * columns + col - 1];
                Count b = values_[static_cast<size_t>(n - 1) * columns + col];
                values_[static_cast<size_t>(n) * columns + col] =
                    (a >= infinity - b) ? infinity : a + b;
            }
        }
    }
    Count operator()(int n, int k) const {
        if (k < 0 || k > n || n < 0) return 0;
        require(k < columns_, "combination column out of range");
        size_t i = static_cast<size_t>(n) * static_cast<size_t>(columns_) + static_cast<size_t>(k);
        require(i < values_.size(), "combination row out of range");
        if (values_[i] == infinity) throw std::overflow_error("binomial coefficient exceeds uint64");
        return values_[i];
    }
};

struct Index {
    std::vector<size_t> off{0};
    Vertices holds, vertices;
    std::vector<size_t> incident_off;
    Vertices incident;
    size_t size() const { return holds.size(); }
    std::span<const Vertex> row(size_t p) const {
        return std::span<const Vertex>(vertices).subspan(off[p], off[p + 1] - off[p]);
    }
    void append(const Vertices& h, const Vertices& p) {
        require(size() < (static_cast<size_t>(absent) / 2), "too many paths for packed incidence");
        holds.push_back(static_cast<Vertex>(h.size()));
        vertices.insert(vertices.end(), h.begin(), h.end());
        vertices.insert(vertices.end(), p.begin(), p.end());
        off.push_back(vertices.size());
    }
    void transpose(Vertex n) {
        incident_off.assign(static_cast<size_t>(n) + 1, 0);
        for (Vertex v : vertices) ++incident_off[v + 1];
        std::partial_sum(incident_off.begin(), incident_off.end(), incident_off.begin());
        incident.resize(vertices.size());
        auto cursor = incident_off;
        for (size_t p = 0; p < size(); ++p) {
            auto members = row(p);
            for (size_t i = 0; i < members.size(); ++i)
                incident[cursor[members[i]]++] = static_cast<Vertex>(2 * p + (i >= holds[p]));
        }
    }
    std::span<const Vertex> touching(Vertex v) const {
        return std::span<const Vertex>(incident).subspan(incident_off[v], incident_off[v + 1] - incident_off[v]);
    }
    size_t capacity_bytes() const {
        return (off.capacity() + incident_off.capacity()) * sizeof(size_t) +
               (holds.capacity() + vertices.capacity() + incident.capacity()) * sizeof(Vertex);
    }
};

struct Edge { Count level; Vertex u, v; };
using Forest = std::vector<Edge>;
struct DSU {
    Vertices parent, size;
    explicit DSU(Vertex n): parent(n), size(n, 1) { std::iota(parent.begin(), parent.end(), 0); }
    Vertex find(Vertex v) {
        while (parent[v] != v) { parent[v] = parent[parent[v]]; v = parent[v]; }
        return v;
    }
    bool join(Vertex u, Vertex v) {
        u = find(u); v = find(v);
        if (u == v) return false;
        if (size[u] < size[v]) std::swap(u, v);
        parent[v] = u; size[u] += size[v]; return true;
    }
};
inline bool descending_edge(const Edge& a, const Edge& b) { return a.level > b.level; }

inline Forest hierarchy(Vertex n, const Index& index, int s, const std::vector<Count>& labels,
                        const Vertices* supplied_order = nullptr) {
    Vertices sorted;
    if (!supplied_order) {
        sorted.reserve(n);
        for (Vertex v = 0; v < n; ++v) if (labels[v]) sorted.push_back(v);
        std::sort(sorted.begin(), sorted.end(), [&](Vertex u, Vertex v) {
            return labels[u] > labels[v] || (labels[u] == labels[v] && u < v);
        });
    }
    const Vertices& order = supplied_order ? *supplied_order : sorted;
    std::vector<uint8_t> active(n, 0);
    Vertices active_h(index.size(), 0), active_p(index.size(), 0), representative(index.size(), absent);
    DSU dsu(n);
    Forest forest; forest.reserve(n);
    auto connect = [&](Count level, Vertex u, Vertex v) {
        if (dsu.join(u, v)) forest.push_back({level, u, v});
    };
    for (Vertex v : order) {
        Count level = labels[v];
        if (!level) break;
        active[v] = 1;
        for (Vertex code : index.touching(v)) {
            size_t p = code / 2;
            if (representative[p] != absent) { connect(level, representative[p], v); continue; }
            if (code & 1) ++active_p[p]; else ++active_h[p];
            if (active_h[p] == index.holds[p] && active_p[p] >= static_cast<Vertex>(s) - index.holds[p]) {
                representative[p] = v;
                for (Vertex u : index.row(p)) if (active[u]) connect(level, v, u);
            }
        }
    }
    return forest;
}

class ThresholdForest {
    Vertices parent_;
    std::vector<Count> level_;
public:
    ThresholdForest() = default;
    ThresholdForest(Vertex n, Forest forest): parent_(n), level_(n, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
        Vertices sizes(n, 1);
        std::sort(forest.begin(), forest.end(), descending_edge);
        for (const auto& edge : forest) {
            Vertex u = edge.u, v = edge.v;
            while (parent_[u] != u) u = parent_[u];
            while (parent_[v] != v) v = parent_[v];
            if (u == v) continue;
            if (sizes[u] < sizes[v]) std::swap(u, v);
            parent_[v] = u; level_[v] = edge.level; sizes[u] += sizes[v];
        }
    }
    bool covers(Count level, Vertex u, Vertex v) const {
        if (!level || u == v) return true;
        while (parent_[u] != u && level_[u] >= level) u = parent_[u];
        while (parent_[v] != v && level_[v] >= level) v = parent_[v];
        return u == v;
    }
};

struct Stats {
    double seed_s = 0, witness_s = 0, witness_hierarchy_s = 0, search_s = 0;
    double transpose_s = 0, peel_s = 0, hierarchy_s = 0, total_s = 0;
    uint64_t expanded = 0, tested = 0, minimum_pass = 0, discarded = 0, coverage_queries = 0;
    uint64_t witness_cliques = 0, witness_members = 0, settled = 0;
    uint64_t paths = 0, members = 0, index_bytes = 0, broadcasts = 0, queue_updates = 0;
};

inline Index witnesses(const Graph& g, int s, const Combinations& choose,
                       const std::vector<Count>& upper, std::vector<Count>& lower, Stats& stats,
                       size_t seed_budget = 0) {
    Vertices order(g.n); std::iota(order.begin(), order.end(), 0);
    auto better = [&](Vertex u, Vertex v) {
        if (upper[u] != upper[v]) return upper[u] > upper[v];
        if (g.row(u).size() != g.row(v).size()) return g.row(u).size() > g.row(v).size();
        return u < v;
    };
    if (seed_budget && seed_budget < order.size()) {
        std::partial_sort(order.begin(), order.begin() + static_cast<ptrdiff_t>(seed_budget), order.end(), better);
        order.resize(seed_budget);
    } else std::sort(order.begin(), order.end(), better);
    Index index;
    std::set<Vertices> seen;
    for (Vertex v : order) {
        if (lower[v] == upper[v]) continue;
        Vertices clique{v};
        auto r = g.row(v);
        Vertices candidates(r.begin(), r.end());
        while (!candidates.empty()) {
            Vertex u = *std::min_element(candidates.begin(), candidates.end(), better);
            clique.push_back(u);
            candidates = intersect(candidates, g.row(u));
        }
        if (clique.size() < static_cast<size_t>(s)) continue;
        std::sort(clique.begin(), clique.end());
        if (!seen.insert(clique).second) continue;
        Count value = choose(static_cast<int>(clique.size()) - 1, s - 1);
        for (Vertex u : clique) {
            require(value <= upper[u], "witness violates upper bound");
            lower[u] = std::max(lower[u], value);
        }
        index.append({}, clique);
    }
    stats.witness_cliques = index.size(); stats.witness_members = index.vertices.size();
    for (Vertex v = 0; v < g.n; ++v) stats.settled += lower[v] == upper[v];
    index.transpose(g.n);
    return index;
}

inline bool minima_certified(const Vertices& h, const Vertices& p, int s,
                             const std::vector<Count>& lower, const std::vector<Count>& upper) {
    int q = s - static_cast<int>(h.size());
    if (q < 0 || static_cast<size_t>(q) > p.size()) return true;
    Count low = infinity, up = infinity;
    for (Vertex v : h) { low = std::min(low, lower[v]); up = std::min(up, upper[v]); }
    if (!q) return low == up;
    if (!h.empty() && low == up && std::all_of(p.begin(), p.end(), [&](Vertex v) { return lower[v] >= low; })) return true;
    if (!h.empty() && up > low) {
        size_t eligible = 0;
        for (Vertex v : p) eligible += lower[v] >= low && upper[v] > low;
        if (eligible >= static_cast<size_t>(q)) return false;
    }
    Vertices ordered = p;
    std::sort(ordered.begin(), ordered.end(), [&](Vertex u, Vertex v) { return lower[u] > lower[v]; });
    for (size_t start = 0; start < ordered.size();) {
        Count level = lower[ordered[start]];
        size_t end = start, uncertain = 0;
        while (end < ordered.size() && lower[ordered[end]] == level) {
            uncertain += upper[ordered[end]] > level; ++end;
        }
        bool allowed = h.empty() || (low >= level && up > level);
        if (allowed && uncertain && start + uncertain >= static_cast<size_t>(q)) return false;
        start = end;
    }
    return true;
}

template<class Check>
bool path_events(const Vertices& h, const Vertices& p, int s,
                 const std::vector<Count>& labels, Check check) {
    int q = s - static_cast<int>(h.size());
    if (q < 0 || static_cast<size_t>(q) > p.size()) return true;
    Count hold_level = infinity;
    for (Vertex v : h) hold_level = std::min(hold_level, labels[v]);
    if (!q) {
        for (Vertex v : h) if (!check(hold_level, h[0], v)) return false;
        return true;
    }
    if (!h.empty() && std::all_of(p.begin(), p.end(), [&](Vertex v) { return labels[v] >= hold_level; })) {
        for (Vertex v : h) if (!check(hold_level, h[0], v)) return false;
        for (Vertex v : p) if (!check(hold_level, h[0], v)) return false;
        return true;
    }
    Vertices ordered = p;
    std::sort(ordered.begin(), ordered.end(), [&](Vertex u, Vertex v) {
        return labels[u] > labels[v] || (labels[u] == labels[v] && u < v);
    });
    Count level = std::min(hold_level, labels[ordered[static_cast<size_t>(q - 1)]]);
    Vertex anchor = h.empty() ? ordered[0] : h[0];
    for (Vertex v : h) if (!check(level, anchor, v)) return false;
    for (int i = 0; i < q; ++i) if (!check(level, anchor, ordered[static_cast<size_t>(i)])) return false;
    for (size_t i = static_cast<size_t>(q); i < ordered.size(); ++i)
        if (!check(std::min(hold_level, labels[ordered[i]]), ordered[0], ordered[i])) return false;
    return true;
}

class Search {
    const Graph& g_;
    int s_;
    const std::vector<Count>& lower_;
    const std::vector<Count>& upper_;
    const ThresholdForest* cover_;
    Stats& stats_;
    Index& index_;
    Vertices holds_, pivots_;
    void expand(Vertices candidates) {
        int q = s_ - static_cast<int>(holds_.size());
        if (q < 0 || pivots_.size() + candidates.size() < static_cast<size_t>(q)) return;
        if (cover_) {
            ++stats_.tested;
            Vertices optional;
            if (q) { optional = pivots_; optional.insert(optional.end(), candidates.begin(), candidates.end()); }
            if (minima_certified(holds_, optional, s_, lower_, upper_)) {
                ++stats_.minimum_pass;
                bool covered = path_events(holds_, optional, s_, lower_, [&](Count level, Vertex u, Vertex v) {
                    ++stats_.coverage_queries;
                    return cover_->covers(level, u, v);
                });
                if (covered) { ++stats_.discarded; return; }
            }
        }
        ++stats_.expanded;
        if (!q) { index_.append(holds_, {}); return; }
        if (candidates.empty()) { index_.append(holds_, pivots_); return; }
        auto close = [&]() {
            Vertices p = pivots_; p.insert(p.end(), candidates.begin(), candidates.end());
            index_.append(holds_, p);
        };
        if (q == 1) { close(); return; }
        Vertex pivot = candidates[0];
        size_t best = 0, minimum = candidates.size();
        for (Vertex u : candidates) {
            size_t degree = intersection_size(candidates, g_.row(u));
            minimum = std::min(minimum, degree);
            if (degree > best) { pivot = u; best = degree; }
        }
        if (minimum + 1 == candidates.size()) { close(); return; }
        Vertices branches{pivot};
        for (Vertex u : candidates) if (u != pivot && !g_.adjacent(pivot, u)) branches.push_back(u);
        for (Vertex u : branches) {
            candidates.erase(std::lower_bound(candidates.begin(), candidates.end(), u));
            auto next = intersect(candidates, g_.row(u));
            auto& added = (u == pivot) ? pivots_ : holds_;
            added.push_back(u); expand(std::move(next)); added.pop_back();
        }
    }
public:
    Search(const Graph& g, int s, const std::vector<Count>& lower, const std::vector<Count>& upper,
           const ThresholdForest* cover, Stats& stats, Index& index):
        g_(g), s_(s), lower_(lower), upper_(upper), cover_(cover), stats_(stats), index_(index) {}
    void run(const Seeds& seeds) {
        for (Vertex v : seeds.order) {
            Vertices candidates;
            for (Vertex u : g_.row(v)) if (seeds.rank[u] > seeds.rank[v]) candidates.push_back(u);
            if (candidates.size() + 1 < static_cast<size_t>(s_)) continue;
            holds_.assign(1, v); pivots_.clear(); expand(std::move(candidates));
        }
    }
};

struct PeelResult { std::vector<Count> core; Vertices order; };
inline PeelResult peel(Vertex n, const Index& index, int s, const Combinations& choose,
                       const std::vector<Count>* lower, const std::vector<Count>* upper, Stats& stats) {
    std::vector<Count> support(n, 0), keys(n), wh(index.size()), wp(index.size());
    Vertices count(index.size());
    std::vector<uint8_t> fixed(n, 0), touched(index.size(), 0), dead(index.size(), 0), dirty(n, 0);
    for (Vertex v = 0; v < n; ++v) fixed[v] = upper && (*lower)[v] == (*upper)[v];
    for (size_t p = 0; p < index.size(); ++p) {
        auto members = index.row(p);
        count[p] = static_cast<Vertex>(members.size()) - index.holds[p];
        int q = s - static_cast<int>(index.holds[p]);
        wh[p] = choose(static_cast<int>(count[p]), q);
        wp[p] = choose(static_cast<int>(count[p]) - 1, q - 1);
        for (size_t i = 0; i < members.size(); ++i)
            if (!fixed[members[i]]) checked_add(support[members[i]], i < index.holds[p] ? wh[p] : wp[p]);
    }
    auto priority = [&](Vertex v) {
        return upper ? std::min((*upper)[v], std::max((*lower)[v], support[v])) : support[v];
    };
    for (Vertex v = 0; v < n; ++v) keys[v] = priority(v);
    Heap heap(keys);
    PeelResult result; result.core.resize(n); result.order.reserve(n);
    Vertices batch, affected, changed;
    Count level = 0;
    while (!heap.empty()) {
        level = std::max(level, keys[heap.first()]); batch.clear();
        while (!heap.empty() && keys[heap.first()] <= level) {
            Vertex v = heap.pop(); result.core[v] = level;
            result.order.push_back(v); batch.push_back(v);
        }
        if (heap.empty()) break;
        affected.clear(); changed.clear();
        for (Vertex v : batch) for (Vertex code : index.touching(v)) {
            Vertex p = code / 2;
            if (dead[p]) continue;
            if (!touched[p]) { touched[p] = 1; affected.push_back(p); }
            if (code & 1) --count[p]; else dead[p] = 1;
        }
        for (Vertex p : affected) {
            int q = s - static_cast<int>(index.holds[p]);
            if (count[p] < static_cast<Vertex>(q)) dead[p] = 1;
            Count new_h = dead[p] ? 0 : choose(static_cast<int>(count[p]), q);
            Count new_p = dead[p] ? 0 : choose(static_cast<int>(count[p]) - 1, q - 1);
            require(new_h <= wh[p] && new_p <= wp[p], "negative path loss");
            Count loss_h = wh[p] - new_h, loss_p = wp[p] - new_p;
            wh[p] = new_h; wp[p] = new_p; touched[p] = 0;
            auto members = index.row(p);
            for (size_t i = 0; i < members.size(); ++i) {
                Count delta = i < index.holds[p] ? loss_h : loss_p;
                if (!delta) continue;
                ++stats.broadcasts;
                Vertex v = members[i];
                if (!heap.contains(v) || fixed[v]) continue;
                require(support[v] >= delta, "support underflow"); support[v] -= delta;
                if (!dirty[v]) { dirty[v] = 1; changed.push_back(v); }
            }
        }
        for (Vertex v : changed) {
            dirty[v] = 0; Count key = priority(v);
            if (key != keys[v]) { heap.decrease(v, key); ++stats.queue_updates; }
        }
    }
    return result;
}

struct Result { PeelResult peel; Forest forest; Index index; Stats stats; };
inline Result solve(const Graph& g, int s, const std::string& mode) {
    require(s >= 2, "s must be at least two");
    require(mode == "full" || mode == "bounds" || mode == "partial" || mode == "partial-budget", "unknown mode");
    Result result;
    Stats& stats = result.stats;
    auto started = Clock::now(), phase = started;
    Seeds seeds(g);
    Combinations choose(seeds.maximum + 1, std::min(s, static_cast<int>(seeds.maximum) + 2));
    stats.seed_s = elapsed(phase);
    std::vector<Count> lower, upper;
    Forest witness_forest;
    ThresholdForest cover;
    if (mode != "full") {
        phase = Clock::now(); lower.assign(g.n, 0); upper.resize(g.n);
        for (Vertex v = 0; v < g.n; ++v) upper[v] = choose(static_cast<int>(seeds.ordinary[v]), s - 1);
        Index witness = witnesses(g, s, choose, upper, lower, stats, mode == "partial-budget" ? 256 : 0);
        stats.witness_s = elapsed(phase);
        phase = Clock::now(); witness_forest = hierarchy(g.n, witness, s, lower);
        cover = ThresholdForest(g.n, witness_forest);
        stats.witness_hierarchy_s = elapsed(phase);
    }
    phase = Clock::now();
    Search search(g, s, lower, upper, mode == "full" || mode == "bounds" ? nullptr : &cover, stats, result.index);
    search.run(seeds);
    stats.search_s = elapsed(phase);
    cover = ThresholdForest();
    phase = Clock::now(); result.index.transpose(g.n); stats.transpose_s = elapsed(phase);
    phase = Clock::now();
    result.peel = peel(g.n, result.index, s, choose, mode == "full" ? nullptr : &lower,
                       mode == "full" ? nullptr : &upper, stats);
    stats.peel_s = elapsed(phase);
    phase = Clock::now();
    Vertices order(result.peel.order.rbegin(), result.peel.order.rend());
    result.forest = hierarchy(g.n, result.index, s, result.peel.core, &order);
    if (!witness_forest.empty()) {
        result.forest.insert(result.forest.end(), witness_forest.begin(), witness_forest.end());
        std::sort(result.forest.begin(), result.forest.end(), descending_edge);
        DSU dsu(g.n); size_t length = 0;
        for (auto edge : result.forest) if (dsu.join(edge.u, edge.v)) result.forest[length++] = edge;
        result.forest.resize(length);
    }
    stats.hierarchy_s = elapsed(phase);
    stats.paths = result.index.size(); stats.members = result.index.vertices.size();
    stats.index_bytes = result.index.capacity_bytes();
    stats.total_s = elapsed(started);
    return result;
}

inline bool forest_covers(Vertex n, Forest lower, Forest upper) {
    std::sort(lower.begin(), lower.end(), descending_edge);
    std::sort(upper.begin(), upper.end(), descending_edge);
    DSU dsu(n); size_t i = 0;
    for (auto edge : upper) {
        while (i < lower.size() && lower[i].level >= edge.level) {
            dsu.join(lower[i].u, lower[i].v); ++i;
        }
        if (edge.level && dsu.find(edge.u) != dsu.find(edge.v)) return false;
    }
    return true;
}
} // namespace partial
