#pragma once

#include "../r1_theory_20260905/cpp/partial_index.hpp"
#include <tuple>

namespace allsize {
using namespace partial;
inline double ms(Clock::time_point start) { return 1000 * elapsed(start); }

struct Work {
    uint64_t events = 0, source_reads = 0, target_reads = 0;
    uint64_t count_reads = 0, updates = 0, rank_calls = 0, list_moves = 0;
};

// This is the existing hold/pivot partition, with stopping rules valid for
// every size in 2..maximum, rather than for one fixed size.
class Search {
    const Graph& graph_;
    int maximum_;
    Index& index_;
    Vertices holds_, pivots_;
    void visit(Vertices candidates) {
        if (holds_.size() > static_cast<size_t>(maximum_)) return;
        if (holds_.size() + pivots_.size() + candidates.size() < 2) return;
        if (holds_.size() == static_cast<size_t>(maximum_)) {
            index_.append(holds_, {});
            return;
        }
        auto close = [&] {
            Vertices optional = pivots_;
            optional.insert(optional.end(), candidates.begin(), candidates.end());
            index_.append(holds_, optional);
        };
        if (candidates.empty() || holds_.size() + 1 == static_cast<size_t>(maximum_)) {
            close();
            return;
        }
        Vertices universal;
        Vertex pivot = absent;
        size_t best = 0;
        for (Vertex u : candidates) {
            const size_t degree = intersection_size(candidates, graph_.row(u));
            if (degree + 1 == candidates.size()) universal.push_back(u);
            else if (pivot == absent || degree > best) { pivot = u; best = degree; }
        }
        if (universal.size() == candidates.size()) { close(); return; }
        const size_t saved = pivots_.size();
        Vertices remaining;
        std::set_difference(candidates.begin(), candidates.end(), universal.begin(), universal.end(),
                            std::back_inserter(remaining));
        candidates.swap(remaining);
        pivots_.insert(pivots_.end(), universal.begin(), universal.end());
        Vertices branches{pivot};
        for (Vertex u : candidates) if (u != pivot && !graph_.adjacent(pivot, u)) branches.push_back(u);
        for (Vertex u : branches) {
            auto at = std::lower_bound(candidates.begin(), candidates.end(), u);
            require(at != candidates.end() && *at == u, "missing search branch");
            candidates.erase(at);
            auto& selected = u == pivot ? pivots_ : holds_;
            selected.push_back(u);
            visit(intersect(candidates, graph_.row(u)));
            selected.pop_back();
        }
        pivots_.resize(saved);
    }
public:
    Search(const Graph& graph, int maximum, Index& index): graph_(graph), maximum_(maximum), index_(index) {}
    void run() {
        for (Vertex v = 0; v < graph_.n; ++v) {
            Vertices later;
            for (Vertex u : graph_.row(v)) if (u > v) later.push_back(u);
            holds_.assign(1, v);
            pivots_.clear();
            visit(std::move(later));
        }
    }
};

struct Layout {
    Index paths;
    std::vector<size_t> reverse_off, channel{0};
    Vertices reverse, owner, initial_top, lo, hi;
    int maximum;
    Layout(const Graph& graph, int bound): maximum(bound) { Search(graph, bound, paths).run(); }
    void prepare(Vertex n) {
        require(paths.vertices.size() < absent, "prototype occurrence IDs exceed uint32");
        initial_top.assign(n, 1);
        reverse_off.assign(static_cast<size_t>(n) + 1, 0);
        owner.resize(paths.vertices.size());
        lo.resize(paths.size()); hi.resize(paths.size());
        for (size_t p = 0; p < paths.size(); ++p) {
            lo[p] = std::max<Vertex>(2, paths.holds[p]);
            hi[p] = std::min<size_t>(maximum, paths.row(p).size());
            require(lo[p] <= hi[p], "empty size range");
            channel.push_back(channel.back() + hi[p] - lo[p] + 1);
            for (size_t i = paths.off[p]; i < paths.off[p + 1]; ++i) {
                const Vertex v = paths.vertices[i];
                owner[i] = static_cast<Vertex>(p);
                ++reverse_off[v + 1];
                initial_top[v] = std::max(initial_top[v], hi[p]);
            }
        }
        std::partial_sum(reverse_off.begin(), reverse_off.end(), reverse_off.begin());
        reverse.resize(paths.vertices.size());
        auto cursor = reverse_off;
        for (size_t i = 0; i < paths.vertices.size(); ++i)
            reverse[cursor[paths.vertices[i]]++] = static_cast<Vertex>(i);
    }
    std::span<const Vertex> touching(Vertex v) const {
        return std::span<const Vertex>(reverse).subspan(reverse_off[v], reverse_off[v + 1] - reverse_off[v]);
    }
    bool pivot(size_t i) const { const size_t p = owner[i]; return i >= paths.off[p] + paths.holds[p]; }
    bool valid(size_t p, int s) const { return static_cast<int>(lo[p]) <= s && s <= static_cast<int>(hi[p]); }
    size_t slot(size_t p, int s) const { return channel[p] + s - lo[p]; }
    size_t bytes() const {
        return paths.capacity_bytes() + (reverse_off.capacity() + channel.capacity()) * sizeof(size_t)
            + (reverse.capacity() + owner.capacity() + initial_top.capacity() + lo.capacity() + hi.capacity()) * sizeof(Vertex);
    }
};

inline Count contribution(const Layout& index, size_t p, bool pivot, int s, int live,
                          const Combinations& choose) {
    return choose(live - static_cast<int>(pivot), s - static_cast<int>(index.paths.holds[p]) - static_cast<int>(pivot));
}

struct Result {
    std::vector<Count> core;
    Work work;
    size_t state_bytes = 0;
    double peel_ms = 0;
};

inline Result fixed(const Layout& index, Vertex n, const Combinations& choose) {
    const auto start = Clock::now();
    Result result;
    result.core.assign((index.maximum + 1) * static_cast<size_t>(n), 0);
    for (int s = 2; s <= index.maximum; ++s) {
        std::vector<Count> degree(n, 0);
        Vertices pivots(index.paths.size());
        std::vector<uint8_t> alive(index.paths.size(), 0);
        for (size_t p = 0; p < index.paths.size(); ++p) if (index.valid(p, s)) {
            alive[p] = 1;
            pivots[p] = static_cast<Vertex>(index.paths.row(p).size() - index.paths.holds[p]);
            const Count h = contribution(index, p, false, s, pivots[p], choose);
            const Count z = contribution(index, p, true, s, pivots[p], choose);
            for (size_t i = index.paths.off[p]; i < index.paths.off[p + 1]; ++i) {
                ++result.work.count_reads;
                checked_add(degree[index.paths.vertices[i]], index.pivot(i) ? z : h);
            }
        }
        Heap queue(degree);
        result.state_bytes = std::max(result.state_bytes,
            degree.capacity() * sizeof(Count) + pivots.capacity() * sizeof(Vertex) + alive.capacity()
            + 2 * static_cast<size_t>(n) * sizeof(Vertex));
        while (!queue.empty()) {
            const Vertex v = queue.pop();
            const Count k = degree[v];
            result.core[static_cast<size_t>(s) * n + v] = k;
            if (!k) continue;
            ++result.work.events;
            for (Vertex occurrence : index.touching(v)) {
                ++result.work.source_reads;
                const size_t p = index.owner[occurrence];
                if (!alive[p]) continue;
                const bool removed_pivot = index.pivot(occurrence);
                const int q = pivots[p];
                Count h, z;
                if (removed_pivot) {
                    h = contribution(index, p, false, s, q, choose) - contribution(index, p, false, s, q - 1, choose);
                    z = contribution(index, p, true, s, q, choose) - contribution(index, p, true, s, q - 1, choose);
                    --pivots[p];
                } else {
                    h = contribution(index, p, false, s, q, choose);
                    z = contribution(index, p, true, s, q, choose);
                    alive[p] = 0;
                }
                if (!h && !z) continue;
                for (size_t i = index.paths.off[p]; i < index.paths.off[p + 1]; ++i) {
                    ++result.work.target_reads;
                    const Vertex u = index.paths.vertices[i];
                    const Count loss = index.pivot(i) ? z : h;
                    if (queue.contains(u) && degree[u] > k && loss) {
                        const Count next = loss >= degree[u] - k ? k : degree[u] - loss;
                        queue.decrease(u, next);
                        ++result.work.updates;
                    }
                }
            }
        }
    }
    result.peel_ms = ms(start);
    return result;
}

inline bool choose_at_most(uint64_t n, int r, Count limit) {
    if (n < static_cast<uint64_t>(r)) return true;
    r = static_cast<int>(std::min<uint64_t>(r, n - r));
    unsigned __int128 value = 1;
    for (int j = 1; j <= r; ++j) {
        value = value * (n - j + 1) / j;
        if (value > limit) return false;
    }
    return value <= limit;
}

inline uint64_t band(Count degree, int s, Vertex n, Work& work) {
    ++work.rank_calls;
    if (!degree) return 0;
    if (s == 2) return degree + 1;
    uint64_t low = s, high = static_cast<uint64_t>(n) + 1;
    while (high - low > 1) {
        const uint64_t middle = (low + high) / 2;
        if (choose_at_most(middle - 1, s - 1, degree)) low = middle;
        else high = middle;
    }
    return low;
}

// An indexed heap with one handle per vertex. Unlike the fixed-s heap, a
// removed handle can be reinserted when its next size becomes active.
class Queue {
    Vertices heap_, position_;
    const Vertices& top_;
    const std::vector<Count>& degree_;
    const std::vector<uint64_t>& band_;
    auto key(Vertex v) const { return std::tuple(band_[v], -static_cast<int>(top_[v]), degree_[v], v); }
    void swap_at(size_t a, size_t b) {
        std::swap(heap_[a], heap_[b]);
        position_[heap_[a]] = static_cast<Vertex>(a);
        position_[heap_[b]] = static_cast<Vertex>(b);
    }
    void up(size_t i) {
        while (i && key(heap_[i]) < key(heap_[(i - 1) / 2])) {
            size_t parent = (i - 1) / 2; swap_at(i, parent); i = parent;
        }
    }
    void down(size_t i) {
        while (2 * i + 1 < heap_.size()) {
            size_t child = 2 * i + 1;
            if (child + 1 < heap_.size() && key(heap_[child + 1]) < key(heap_[child])) ++child;
            if (key(heap_[i]) <= key(heap_[child])) break;
            swap_at(i, child); i = child;
        }
    }
public:
    Queue(const Vertices& top, const std::vector<Count>& degree, const std::vector<uint64_t>& bands)
        : position_(top.size(), absent), top_(top), degree_(degree), band_(bands) { heap_.reserve(top.size()); }
    bool empty() const { return heap_.empty(); }
    void insert(Vertex v) { require(position_[v] == absent, "duplicate queue handle"); position_[v] = heap_.size(); heap_.push_back(v); up(heap_.size() - 1); }
    Vertex pop() {
        Vertex v = heap_[0]; swap_at(0, heap_.size() - 1); heap_.pop_back(); position_[v] = absent;
        if (!heap_.empty()) down(0);
        return v;
    }
    void decrease(Vertex v) { require(position_[v] != absent, "missing queue handle"); up(position_[v]); }
    size_t bytes() const { return (heap_.capacity() + position_.capacity()) * sizeof(Vertex); }
};

inline Result joint(const Layout& index, Vertex n, const Combinations& choose, bool bucketed) {
    const auto start = Clock::now();
    Result result;
    result.core.assign((index.maximum + 1) * static_cast<size_t>(n), 0);
    Vertices top = index.initial_top;
    Vertices live_pivots(index.channel.back()), hold_min(index.paths.size(), index.maximum);
    Vertices heads, next, previous;
    if (bucketed) {
        heads.assign(index.channel.back(), absent);
        next.assign(index.paths.vertices.size(), absent);
        previous.assign(index.paths.vertices.size(), absent);
    }
    std::vector<Count> degree(n, 0), level(index.maximum + 1, 0);
    std::vector<uint64_t> bands(n, 0);
    Queue queue(top, degree, bands);
    auto insert = [&](Vertex occurrence, int s) {
        if (!bucketed) return;
        const size_t p = index.owner[occurrence];
        if (!index.valid(p, s)) return;
        const size_t slot = index.slot(p, s);
        const Vertex first = heads[slot];
        next[occurrence] = first;
        previous[occurrence] = absent;
        if (first != absent) previous[first] = occurrence;
        heads[slot] = occurrence;
        ++result.work.list_moves;
    };
    auto erase = [&](Vertex occurrence, int s) {
        if (!bucketed) return;
        const size_t p = index.owner[occurrence];
        if (!index.valid(p, s)) return;
        if (previous[occurrence] == absent) heads[index.slot(p, s)] = next[occurrence];
        else next[previous[occurrence]] = next[occurrence];
        if (next[occurrence] != absent) previous[next[occurrence]] = previous[occurrence];
        next[occurrence] = previous[occurrence] = absent;
        ++result.work.list_moves;
    };
    for (size_t p = 0; p < index.paths.size(); ++p) {
        const Vertex q = static_cast<Vertex>(index.paths.row(p).size() - index.paths.holds[p]);
        std::fill(live_pivots.begin() + index.channel[p], live_pivots.begin() + index.channel[p + 1], q);
        for (size_t i = index.paths.off[p]; i < index.paths.off[p + 1]; ++i) insert(static_cast<Vertex>(i), top[index.paths.vertices[i]]);
    }
    auto enter = [&](Vertex v) {
        const int s = top[v];
        if (s < 2) return;
        degree[v] = 0;
        for (Vertex occurrence : index.touching(v)) {
            ++result.work.count_reads;
            const size_t p = index.owner[occurrence];
            if (!index.valid(p, s) || static_cast<int>(hold_min[p]) < s) continue;
            checked_add(degree[v], contribution(index, p, index.pivot(occurrence), s, live_pivots[index.slot(p, s)], choose));
        }
        degree[v] = std::max(degree[v], level[s]);
        bands[v] = band(degree[v], s, n, result.work);
        queue.insert(v);
    };
    for (Vertex v = 0; v < n; ++v) enter(v);
    std::tuple<uint64_t, int, Count> last{0, -index.maximum, 0};
    while (!queue.empty()) {
        const Vertex v = queue.pop();
        const int s = top[v];
        const Count k = degree[v];
        const auto selected = std::tuple(bands[v], -s, k);
        require(selected >= last, "joint priority moved backward");
        last = selected;
        level[s] = k;
        result.core[static_cast<size_t>(s) * n + v] = k;
        ++result.work.events;
        for (Vertex occurrence : index.touching(v)) {
            ++result.work.source_reads;
            const size_t p = index.owner[occurrence];
            erase(occurrence, s);
            if (!index.valid(p, s)) { insert(occurrence, s - 1); continue; }
            const bool removed_pivot = index.pivot(occurrence);
            const size_t slot = index.slot(p, s);
            const int q = live_pivots[slot];
            Count h = 0, z = 0;
            if (static_cast<int>(hold_min[p]) >= s) {
                if (removed_pivot) {
                    h = contribution(index, p, false, s, q, choose) - contribution(index, p, false, s, q - 1, choose);
                    z = contribution(index, p, true, s, q, choose) - contribution(index, p, true, s, q - 1, choose);
                } else {
                    h = contribution(index, p, false, s, q, choose);
                    z = contribution(index, p, true, s, q, choose);
                }
            }
            if (removed_pivot) { require(q > 0, "negative live pivot count"); --live_pivots[slot]; }
            else hold_min[p] = std::min<Vertex>(hold_min[p], s - 1);
            auto update = [&](Vertex i) {
                ++result.work.target_reads;
                const Vertex u = index.paths.vertices[i];
                const Count loss = index.pivot(i) ? z : h;
                if (u == v || static_cast<int>(top[u]) != s || degree[u] <= k || !loss) return;
                degree[u] = loss >= degree[u] - k ? k : degree[u] - loss;
                if (!degree[u] || !choose_at_most(bands[u] - 1, s - 1, degree[u]))
                    bands[u] = band(degree[u], s, n, result.work);
                queue.decrease(u);
                ++result.work.updates;
            };
            if (h || z) {
                if (bucketed) for (Vertex i = heads[slot]; i != absent; i = next[i]) update(i);
                else for (size_t i = index.paths.off[p]; i < index.paths.off[p + 1]; ++i) update(static_cast<Vertex>(i));
            }
            insert(occurrence, s - 1);
        }
        --top[v];
        enter(v);
    }
    result.state_bytes = (top.capacity() + live_pivots.capacity() + hold_min.capacity() + heads.capacity()
        + next.capacity() + previous.capacity()) * sizeof(Vertex)
        + (degree.capacity() + level.capacity()) * sizeof(Count) + bands.capacity() * sizeof(uint64_t) + queue.bytes();
    result.peel_ms = ms(start);
    return result;
}
}
