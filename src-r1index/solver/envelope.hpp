#pragma once

#include "ceiling.hpp"
#include <boost/multiprecision/cpp_int.hpp>

namespace envelope {
using namespace allsize;
using Big = boost::multiprecision::cpp_int;

inline Vertex coloring_bound(const Graph& graph) {
    Vertices color(graph.n, absent), seen(static_cast<size_t>(graph.n) + 1, absent);
    Vertex colors = 0;
    for (Vertex i = graph.n; i > 0; --i) {
        const Vertex v = i - 1;
        for (Vertex u : graph.row(v)) if (color[u] != absent) seen[color[u]] = v;
        Vertex c = 0;
        while (seen[c] == v) ++c;
        color[v] = c;
        colors = std::max(colors, c + 1);
    }
    for (Vertex v = 0; v < graph.n; ++v) for (Vertex u : graph.row(v))
        require(color[u] != color[v], "invalid coloring certificate");
    return std::max<Vertex>(1, colors ? colors - 1 : 0);
}

inline Big binomial(Vertex n, int k) {
    if (k < 0 || static_cast<Vertex>(k) > n) return 0;
    k = std::min<int>(k, n - k);
    Big value = 1;
    for (int j = 1; j <= k; ++j) { value *= n - j + 1; value /= j; }
    return value;
}

struct Curve {
    size_t stride;
    std::vector<Count> rounded;
    std::vector<uint8_t> flags;
    Curve(Vertex maximum_core, int maximum_size, Vertex b, bool hybrid, bool turan = false): stride(static_cast<size_t>(maximum_core) + 1),
        rounded(stride * (maximum_size + 1)), flags(rounded.size()) {
        if (turan) {
            std::vector<Big> polynomial(maximum_size), without(maximum_size);
            polynomial[0] = 1;
            for (Vertex q = 0; q <= maximum_core; ++q) {
                for (int s = 2; s <= maximum_size; ++s) {
                    const size_t i = slot(s, q);
                    if (static_cast<Vertex>(s - 1) > b || polynomial[s - 1] >= infinity) flags[i] = 1;
                    else rounded[i] = polynomial[s - 1].convert_to<Count>();
                }
                if (q == maximum_core) break;
                // Divide by one smallest-part factor, then increase that part.
                const Vertex size = q / b;
                without[0] = 1;
                for (int j = 1; j < maximum_size; ++j) without[j] = polynomial[j] - size * without[j - 1];
                for (int j = 1; j < maximum_size; ++j) polynomial[j] += without[j - 1];
            }
            return;
        }
        for (int s = 2; s <= maximum_size; ++s) {
            const int a = s - 1;
            Big denominator = 1;
            for (int j = 0; j < a; ++j) denominator *= b;
            const Big coefficient = binomial(b, a);
            for (Vertex q = 0; q <= maximum_core; ++q) {
                const size_t i = slot(s, q);
                if (hybrid && static_cast<Vertex>(a) > b) { flags[i] = 1; continue; }
                Big numerator, divisor = 1;
                if (hybrid && q > b) {
                    numerator = coefficient;
                    for (int j = 0; j < a; ++j) numerator *= q;
                    divisor = denominator;
                } else numerator = binomial(q, a);
                const Big value = (numerator + divisor - 1) / divisor;
                if (value >= infinity) flags[i] |= 1;
                else rounded[i] = value.convert_to<Count>();
                if (numerator % divisor != 0) flags[i] |= 2;
            }
        }
    }
    size_t slot(int s, Vertex q) const { return static_cast<size_t>(s) * stride + q; }
    bool fits(int s, Vertex q, Count d) const {
        const size_t i = slot(s, q);
        return !(flags[i] & 1) && rounded[i] <= d;
    }
    Count lower(int s, Vertex q) const {
        const size_t i = slot(s, q);
        require(!(flags[i] & 1), "active bound exceeds the count domain");
        return rounded[i];
    }
    Count exact(int s, Vertex q) const {
        require(flags[slot(s, q)] == 0, "certificate is nonintegral or exceeds count domain");
        return rounded[slot(s, q)];
    }
    Vertex rank(int s, Count d, Vertex cap, Work& work) const {
        ++work.rank_calls;
        if (!d) return 0;
        Vertex low = s - 1;
        require(low <= cap && fits(s, low, d), "invalid positive frontier");
        uint64_t high = static_cast<uint64_t>(cap) + 1;
        while (high - low > 1) {
            const Vertex middle = (static_cast<uint64_t>(low) + high) / 2;
            if (fits(s, middle, d)) low = middle;
            else high = middle;
        }
        return low;
    }
    size_t bytes() const { return rounded.capacity() * sizeof(Count) + flags.capacity(); }
};

struct EnvelopeResult {
    CeilingResult result;
    uint64_t nonbinomial_events = 0, nonbinomial_outputs = 0;
};

inline EnvelopeResult peel(const Layout& index, Vertex n, const Combinations& choose,
                           const Vertices& ordinary, const Curve& curve, Vertex b, bool hybrid) {
    const auto start = Clock::now();
    EnvelopeResult output;
    auto& result = output.result.data;
    result.core.assign((index.maximum + 1) * static_cast<size_t>(n), 0);
    std::copy(ordinary.begin(), ordinary.end(), result.core.begin() + 2 * static_cast<size_t>(n));
    Vertices top = index.initial_top;
    Vertices live_pivots(index.channel.back()), hold_min(index.paths.size(), index.maximum);
    std::vector<Count> degree(n, 0), level(index.maximum + 1, 0);
    std::vector<uint64_t> bands(n, 0);
    Queue queue(top, degree, bands);
    Vertex current_band = 0;
    for (size_t p = 0; p < index.paths.size(); ++p) {
        const Vertex q = index.paths.row(p).size() - index.paths.holds[p];
        std::fill(live_pivots.begin() + index.channel[p], live_pivots.begin() + index.channel[p + 1], q);
    }
    auto floor = [&](int s) { return std::max(level[s], curve.lower(s, current_band)); };
    auto update_band = [&](Vertex v) { bands[v] = curve.rank(top[v], degree[v], ordinary[v], result.work); };
    auto enter = [&](Vertex v) {
        const int s = top[v];
        if (s < 3) return;
        degree[v] = 0;
        for (Vertex occurrence : index.touching(v)) {
            ++result.work.count_reads;
            const size_t p = index.owner[occurrence];
            if (!index.valid(p, s) || static_cast<int>(hold_min[p]) < s) continue;
            checked_add(degree[v], contribution(index, p, index.pivot(occurrence), s, live_pivots[index.slot(p, s)], choose));
        }
        degree[v] = std::max(degree[v], floor(s));
        update_band(v);
        queue.insert(v);
    };
    for (Vertex v = 0; v < n; ++v) enter(v);
    while (!queue.empty()) {
        const Vertex v = queue.pop();
        require(bands[v] >= current_band, "envelope band moved backward");
        current_band = bands[v];
        const int old_top = top[v];
        const bool certified = current_band == ordinary[v];
        const int new_top = certified ? 2 : old_top - 1;
        if (certified) {
            ++output.result.ceiling_events;
            output.result.ceiling_outputs += old_top - new_top;
            if (hybrid && current_band > b) {
                ++output.nonbinomial_events;
                output.nonbinomial_outputs += old_top - new_top;
            }
        }
        for (int s = new_top + 1; s <= old_top; ++s) {
            const Count value = certified ? curve.exact(s, ordinary[v]) : degree[v];
            require(value >= level[s], "envelope output below earlier level");
            level[s] = value;
            result.core[static_cast<size_t>(s) * n + v] = value;
        }
        ++result.work.events;
        for (Vertex occurrence : index.touching(v)) {
            ++result.work.source_reads;
            const size_t p = index.owner[occurrence];
            const int low = std::max<int>(index.lo[p], new_top + 1);
            const int high = std::min<int>(index.hi[p], old_top);
            if (low > high) continue;
            const bool removed_pivot = index.pivot(occurrence);
            if (static_cast<int>(hold_min[p]) >= low) {
                for (size_t i = index.paths.off[p]; i < index.paths.off[p + 1]; ++i) {
                    ++result.work.target_reads;
                    const Vertex u = index.paths.vertices[i];
                    const int s = top[u];
                    if (u == v || s < low || s > high || static_cast<int>(hold_min[p]) < s) continue;
                    const Count minimum = floor(s);
                    if (degree[u] <= minimum) continue;
                    const int q = live_pivots[index.slot(p, s)];
                    const bool target_pivot = index.pivot(i);
                    Count loss = contribution(index, p, target_pivot, s, q, choose);
                    if (removed_pivot) loss -= contribution(index, p, target_pivot, s, q - 1, choose);
                    if (!loss) continue;
                    degree[u] = loss >= degree[u] - minimum ? minimum : degree[u] - loss;
                    if (!degree[u] || !curve.fits(s, bands[u], degree[u])) update_band(u);
                    queue.decrease(u);
                    ++result.work.updates;
                }
            }
            if (removed_pivot) {
                for (int s = low; s <= high; ++s) {
                    Vertex& q = live_pivots[index.slot(p, s)];
                    require(q > 0, "negative envelope pivot count");
                    --q;
                }
            } else hold_min[p] = std::min<Vertex>(hold_min[p], new_top);
        }
        top[v] = new_top;
        enter(v);
    }
    result.state_bytes = (top.capacity() + live_pivots.capacity() + hold_min.capacity()) * sizeof(Vertex)
        + (degree.capacity() + level.capacity()) * sizeof(Count) + bands.capacity() * sizeof(uint64_t) + queue.bytes();
    result.peel_ms = ms(start);
    return output;
}
}
