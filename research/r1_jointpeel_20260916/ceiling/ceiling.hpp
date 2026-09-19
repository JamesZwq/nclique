#pragma once

#include "../joint.hpp"

namespace allsize {
struct CeilingResult {
    Result data;
    uint64_t ceiling_events = 0, ceiling_outputs = 0;
};

inline Result fixed_sparse(const Layout& index, Vertex n, const Combinations& choose, const Vertices& ordinary) {
    const auto start = Clock::now();
    Result result;
    result.core.assign((index.maximum + 1) * static_cast<size_t>(n), 0);
    std::copy(ordinary.begin(), ordinary.end(), result.core.begin() + 2 * static_cast<size_t>(n));
    for (int s = 3; s <= index.maximum; ++s) {
        std::vector<Count> initial(n, 0);
        Vertices pivots(index.paths.size());
        std::vector<uint8_t> alive(index.paths.size(), 0);
        for (size_t p = 0; p < index.paths.size(); ++p) if (index.valid(p, s)) {
            alive[p] = 1;
            pivots[p] = static_cast<Vertex>(index.paths.row(p).size() - index.paths.holds[p]);
            const Count h = contribution(index, p, false, s, pivots[p], choose);
            const Count z = contribution(index, p, true, s, pivots[p], choose);
            for (size_t i = index.paths.off[p]; i < index.paths.off[p + 1]; ++i) {
                ++result.work.count_reads;
                checked_add(initial[index.paths.vertices[i]], index.pivot(i) ? z : h);
            }
        }
        Vertices local(n, absent), original;
        std::vector<Count> degree;
        for (Vertex v = 0; v < n; ++v) if (initial[v]) {
            local[v] = static_cast<Vertex>(degree.size());
            original.push_back(v);
            degree.push_back(initial[v]);
        }
        const size_t preparation_bytes = initial.capacity() * sizeof(Count);
        std::vector<Count>().swap(initial);
        Heap queue(degree);
        result.state_bytes = std::max(result.state_bytes,
            degree.capacity() * sizeof(Count) + pivots.capacity() * sizeof(Vertex) + alive.capacity()
            + (local.capacity() + original.capacity()) * sizeof(Vertex)
            + std::max(preparation_bytes, 2 * degree.size() * sizeof(Vertex)));
        while (!queue.empty()) {
            const Vertex local_v = queue.pop(), v = original[local_v];
            const Count k = degree[local_v];
            result.core[static_cast<size_t>(s) * n + v] = k;
            ++result.work.events;
            for (Vertex occurrence : index.touching(v)) {
                ++result.work.source_reads;
                const size_t p = index.owner[occurrence];
                if (!alive[p]) continue;
                const int q = pivots[p];
                Count h, z;
                if (index.pivot(occurrence)) {
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
                    const Vertex u = local[index.paths.vertices[i]];
                    const Count loss = index.pivot(i) ? z : h;
                    if (u != absent && queue.contains(u) && degree[u] > k && loss) {
                        queue.decrease(u, loss >= degree[u] - k ? k : degree[u] - loss);
                        ++result.work.updates;
                    }
                }
            }
        }
    }
    result.peel_ms = ms(start);
    return result;
}

inline CeilingResult ceiling(const Layout& index, Vertex n, const Combinations& choose,
                             const Vertices& ordinary, bool batch) {
    const auto start = Clock::now();
    CeilingResult output;
    auto& result = output.data;
    result.core.assign((index.maximum + 1) * static_cast<size_t>(n), 0);
    std::copy(ordinary.begin(), ordinary.end(), result.core.begin() + 2 * static_cast<size_t>(n));
    Vertices top = index.initial_top;
    Vertices live_pivots(index.channel.back()), hold_min(index.paths.size(), index.maximum);
    std::vector<Count> degree(n, 0), level(index.maximum + 1, 0);
    std::vector<uint64_t> bands(n, 0);
    Queue queue(top, degree, bands);
    uint64_t current_band = 0;
    for (size_t p = 0; p < index.paths.size(); ++p) {
        const Vertex q = static_cast<Vertex>(index.paths.row(p).size() - index.paths.holds[p]);
        std::fill(live_pivots.begin() + index.channel[p], live_pivots.begin() + index.channel[p + 1], q);
    }
    auto floor = [&](int s) { return std::max(level[s], choose(static_cast<int>(current_band) - 1, s - 1)); };
    auto update_band = [&](Vertex v) {
        bands[v] = std::min<uint64_t>(ordinary[v] + 1, band(degree[v], top[v], ordinary[v] + 1, result.work));
    };
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
        require(bands[v] >= current_band, "ceiling band moved backward");
        current_band = bands[v];
        const int old_top = top[v];
        const bool certified = current_band == static_cast<uint64_t>(ordinary[v]) + 1;
        const int new_top = certified && batch ? 2 : old_top - 1;
        if (certified) {
            ++output.ceiling_events;
            output.ceiling_outputs += old_top - new_top;
        }
        for (int s = new_top + 1; s <= old_top; ++s) {
            const Count value = certified ? choose(ordinary[v], s - 1) : degree[v];
            require(value >= level[s], "ceiling output below earlier level");
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
                    if (!degree[u] || !choose_at_most(bands[u] - 1, s - 1, degree[u])) update_band(u);
                    queue.decrease(u);
                    ++result.work.updates;
                }
            }
            if (removed_pivot) {
                for (int s = low; s <= high; ++s) {
                    Vertex& q = live_pivots[index.slot(p, s)];
                    require(q > 0, "negative ceiling pivot count");
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
