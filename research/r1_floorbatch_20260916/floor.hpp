#pragma once

#include "../r1_envelopecore_20260916/envelope.hpp"

namespace floorbatch {
using namespace envelope;

struct Report {
    allsize::Result data;
    uint64_t ceiling_events = 0, ceiling_outputs = 0;
    uint64_t floor_events = 0, floor_outputs = 0, seed_reads = 0;
    uint64_t audited_states = 0;
    double seed_ms = 0;
};

template<bool Audit = false>
Report peel(const Layout& index, Vertex n, const Combinations& choose,
            const Vertices& ordinary, const Curve& curve, bool seeded, bool batch,
            const std::vector<Count>* expected = nullptr) {
    const auto start = Clock::now();
    Report output;
    auto& r = output.data;
    r.core.assign((index.maximum + 1) * static_cast<size_t>(n), 0);
    std::copy(ordinary.begin(), ordinary.end(), r.core.begin() + 2 * static_cast<size_t>(n));
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
    size_t seed_scratch_bytes = 0;
    if (seeded) {
        const auto seed_start = Clock::now();
        std::vector<Count> counts(n);
        seed_scratch_bytes = counts.capacity() * sizeof(Count);
        for (int s = 3; s <= index.maximum; ++s) {
            std::fill(counts.begin(), counts.end(), 0);
            for (size_t p = 0; p < index.paths.size(); ++p) if (index.valid(p, s)) {
                const int q = index.paths.row(p).size() - index.paths.holds[p];
                const Count h = contribution(index, p, false, s, q, choose);
                const Count z = contribution(index, p, true, s, q, choose);
                for (size_t i = index.paths.off[p]; i < index.paths.off[p + 1]; ++i) {
                    ++output.seed_reads;
                    checked_add(counts[index.paths.vertices[i]], index.pivot(i) ? z : h);
                }
            }
            Count minimum = infinity;
            for (Vertex v = 0; v < n; ++v) if (static_cast<int>(top[v]) >= s) {
                require(counts[v] > 0, "positive initial layer has zero degree");
                minimum = std::min(minimum, counts[v]);
                if (static_cast<int>(top[v]) == s) degree[v] = counts[v];
            }
            level[s] = minimum == infinity ? 0 : minimum;
        }
        output.seed_ms = ms(seed_start);
    }
    auto floor = [&](int s) { return std::max(level[s], curve.lower(s, current_band)); };
    auto count = [&](Vertex v, int s) {
        Count value = 0;
        for (Vertex occurrence : index.touching(v)) {
            ++r.work.count_reads;
            const size_t p = index.owner[occurrence];
            if (index.valid(p, s) && static_cast<int>(hold_min[p]) >= s)
                checked_add(value, contribution(index, p, index.pivot(occurrence), s, live_pivots[index.slot(p, s)], choose));
        }
        return std::max(value, floor(s));
    };
    auto insert = [&](Vertex v) {
        bands[v] = curve.rank(top[v], degree[v], ordinary[v], r.work);
        queue.insert(v);
    };
    for (Vertex v = 0; v < n; ++v) if (top[v] >= 3) {
        if (!seeded) degree[v] = count(v, top[v]);
        insert(v);
    }
    auto audit = [&] {
        if constexpr (Audit) {
            require(expected != nullptr, "audit requires oracle outputs");
            ++output.audited_states;
            for (size_t p = 0; p < index.paths.size(); ++p) {
                const auto row = index.paths.row(p);
                for (int s = std::max<int>(3, index.lo[p]); s <= static_cast<int>(index.hi[p]); ++s) {
                    bool holds_live = true;
                    Vertex pivots = 0;
                    for (size_t j = 0; j < row.size(); ++j) {
                        if (j < index.paths.holds[p]) holds_live &= top[row[j]] >= static_cast<Vertex>(s);
                        else pivots += top[row[j]] >= static_cast<Vertex>(s);
                    }
                    require(holds_live == (hold_min[p] >= static_cast<Vertex>(s)), "hold state differs from live set");
                    require(pivots == live_pivots[index.slot(p, s)], "pivot channel differs from live set");
                }
            }
            for (Vertex v = 0; v < n; ++v) {
                for (int s = 3; s <= index.maximum; ++s) {
                    const size_t at = static_cast<size_t>(s) * n + v;
                    if (static_cast<int>(top[v]) >= s) require((*expected)[at] >= floor(s), "floor exceeds an unfinished answer");
                    else require(r.core[at] == (*expected)[at], "completed answer differs");
                }
                if (top[v] >= 3) {
                    const auto saved = r.work.count_reads;
                    require(degree[v] == count(v, top[v]), "maintained clipped degree differs");
                    r.work.count_reads = saved;
                }
            }
        }
    };
    audit();
    while (!queue.empty()) {
        const Vertex v = queue.pop();
        require(bands[v] >= current_band, "floor frontier moved backward");
        current_band = bands[v];
        const int old_top = top[v];
        int new_top = old_top - 1;
        Count cached_entry = 0;
        auto assign = [&](int s, Count value) {
            require(value >= level[s], "assigned value below layer floor");
            if constexpr (Audit) require(value == (*expected)[static_cast<size_t>(s) * n + v], "assignment differs from oracle");
            level[s] = value;
            r.core[static_cast<size_t>(s) * n + v] = value;
        };
        if (current_band == ordinary[v]) {
            new_top = 2;
            ++output.ceiling_events;
            output.ceiling_outputs += old_top - new_top;
            for (int s = 3; s <= old_top; ++s) assign(s, curve.exact(s, ordinary[v]));
        } else {
            assign(old_top, degree[v]);
            const auto before = output.floor_outputs;
            while (new_top >= 3) {
                cached_entry = count(v, new_top);
                if (!batch || cached_entry > floor(new_top)) break;
                assign(new_top, cached_entry);
                ++output.floor_outputs;
                --new_top;
            }
            output.floor_events += output.floor_outputs != before;
        }
        ++r.work.events;
        // One target counter uses only its own size, even for a size interval.
        for (Vertex occurrence : index.touching(v)) {
            ++r.work.source_reads;
            const size_t p = index.owner[occurrence];
            const int low = std::max<int>(index.lo[p], new_top + 1), high = std::min<int>(index.hi[p], old_top);
            if (low > high) continue;
            const bool removed_pivot = index.pivot(occurrence);
            if (static_cast<int>(hold_min[p]) >= low) {
                for (size_t i = index.paths.off[p]; i < index.paths.off[p + 1]; ++i) {
                    ++r.work.target_reads;
                    const Vertex u = index.paths.vertices[i];
                    const int s = top[u];
                    if (u == v || s < low || s > high || static_cast<int>(hold_min[p]) < s) continue;
                    const Count minimum = floor(s);
                    if (degree[u] <= minimum) continue;
                    const int q = live_pivots[index.slot(p, s)];
                    const bool pivot = index.pivot(i);
                    Count loss = contribution(index, p, pivot, s, q, choose);
                    if (removed_pivot) loss -= contribution(index, p, pivot, s, q - 1, choose);
                    if (!loss) continue;
                    degree[u] = loss >= degree[u] - minimum ? minimum : degree[u] - loss;
                    if (!degree[u] || !curve.fits(s, bands[u], degree[u]))
                        bands[u] = curve.rank(s, degree[u], ordinary[u], r.work);
                    queue.decrease(u);
                    ++r.work.updates;
                }
            }
            if (removed_pivot) {
                for (int s = low; s <= high; ++s) {
                    auto& q = live_pivots[index.slot(p, s)];
                    require(q > 0, "negative live pivot count");
                    --q;
                }
            } else hold_min[p] = std::min<Vertex>(hold_min[p], new_top);
        }
        top[v] = new_top;
        if (new_top >= 3) { degree[v] = cached_entry; insert(v); }
        audit();
    }
    r.state_bytes = (top.capacity() + live_pivots.capacity() + hold_min.capacity()) * sizeof(Vertex)
        + (degree.capacity() + level.capacity()) * sizeof(Count) + bands.capacity() * sizeof(uint64_t)
        + queue.bytes() + seed_scratch_bytes;
    r.peel_ms = ms(start);
    return output;
}
}
