#pragma once
// terminal::Solver<Count>::solve (sink form, no audit, no replay), statement for statement, with per-phase timers:
// a profiling aid for tail_solver.hpp (2026-09-23).  Only the timers are added.
#include "../r1_terminal_20260918/terminal.hpp"
#include <span>
namespace terminalprobe {
using namespace allsize;
struct Phases { double alloc_ms = 0, init_ms = 0, bounds_ms = 0, heap_ms = 0, peel_ms = 0, total_ms = 0; };
template<class Count> struct Solver {
    using Base = fullrange::Kernel<Count>; using Choose = typename Base::Combinations; using Queue = typename Base::Queue;
    using TS = terminal::Solver<Count>; using Value = typename TS::Value;
    template<class Sink> static Phases solve(const Graph& graph, const terminal::Index& index, const Choose& choose, const Vertices& ordinary, Sink sink) {
        using terminal::RowId; using terminal::Offset;
        Phases ph; const auto start = Clock::now(); const Vertex n = graph.n;
        tworoads::Statistics stats; orderdp::Extra extra; terminal::Metrics metrics; uint64_t events = 0;
        std::vector<Count> core(2 * static_cast<size_t>(n), 0);
        auto slot = [&](int s) { return static_cast<size_t>(s & 1) * n; };
        std::copy(ordinary.begin(), ordinary.end(), core.begin() + slot(2));
        sink(2, std::span<const Count>(core).subspan(slot(2), n));
        Vertices order(n), next_order; std::iota(order.begin(), order.end(), 0);
        if (!std::is_sorted(ordinary.begin(), ordinary.end()))
            std::sort(order.begin(), order.end(), [&](Vertex a, Vertex b) { return ordinary[a] < ordinary[b] || (ordinary[a] == ordinary[b] && a < b); });
        next_order.reserve(n);
        const Vertex maximum = ordinary.empty() ? 0 : *std::max_element(ordinary.begin(), ordinary.end());
        struct Cache { Count key = Base::infinity, value = 0; };
        std::vector<Cache> cache(4096);
        for (int s = 3; s <= index.maximum; ++s) {
            auto upper = std::span<Count>(core).subspan(slot(s), n);
            const auto previous = std::span<const Count>(core).subspan(slot(s - 1), n);
            std::fill(upper.begin(), upper.end(), 0);
            auto t0 = Clock::now();
            std::vector<Count> support(n, 0), wh(index.rows.size(), 0), wp(index.rows.size(), 0), wx(index.group_row.size(), 0);
            Vertices count(index.rows.size(), 0), choices(index.group_row.size(), 0);
            std::vector<Offset> choice_off(index.group_row.size(), 0); Vertices choice_size(index.group_row.size(), 0), scratch;
            std::vector<uint8_t> live(n, 0), touched(index.rows.size(), 0), dead(index.rows.size(), 1), dirty(n, 0);
            ph.alloc_ms += ms(t0); t0 = Clock::now();
            for (RowId p = 0; p < index.rows.size(); ++p) {
                const auto& row = index.rows[p]; if (!row.valid(s)) continue;
                dead[p] = 0; count[p] = row.pivots();
                const Vertex z = row.end - row.pivot_end;
                const auto value = TS::coefficients(index, p, s, count[p], z, choose, metrics);
                wh[p] = value.h; wp[p] = value.q;
                auto initialize = [&](Offset begin, Offset end, Count weight) {
                    if (!weight) return;
                    for (Offset i = begin; i < end; ++i) Base::checked_add(support[index.members[i]], weight);
                };
                initialize(row.begin, row.hold_end, value.h); initialize(row.hold_end, row.pivot_end, value.q);
                if (row.group != absent) {
                    const Vertex g = row.group; wx[g] = value.x; choices[g] = z;
                    initialize(row.pivot_end, row.end, value.x);
                    choice_off[g] = scratch.size();
                    if (value.x) { choice_size[g] = z; scratch.insert(scratch.end(), index.members.begin() + row.pivot_end, index.members.begin() + row.end); }
                }
            }
            ph.init_ms += ms(t0); t0 = Clock::now();
            Vertex remaining = 0; next_order.clear(); std::fill(cache.begin(), cache.end(), Cache{});
            for (Vertex v = 0; v < n; ++v) {
                if (!support[v]) { next_order.push_back(v); continue; }
                live[v] = 1; ++remaining; const Count a = previous[v];
                auto& entry = cache[static_cast<size_t>((a ^ (a >> 17) ^ (a >> 37)) & (cache.size() - 1))];
                if (entry.key != a) entry = {a, Base::integer_upper(a, s - 2, maximum, stats)};
                upper[v] = entry.value;
                require(upper[v] > 0, "zero bound for a clique member");
            }
            ph.bounds_ms += ms(t0); t0 = Clock::now();
            Queue heap(support, upper, false, extra);
            ph.heap_ms += ms(t0); t0 = Clock::now();
            Vertices batch, changed; std::vector<RowId> affected; size_t cursor = 0;
            auto stream_key = [&]() -> Count {
                while (cursor < order.size() && !live[order[cursor]]) ++cursor;
                return cursor < order.size() ? upper[order[cursor]] : Base::infinity;
            };
            if (!remaining) { sink(s, upper); ph.peel_ms += ms(t0); break; }
            Count level = 0;
            while (remaining) {
                level = std::max(level, std::min(stream_key(), heap.first_key()));
                require(level != Base::infinity, "missing factored minimum"); batch.clear();
                while (std::min(stream_key(), heap.first_key()) <= level) {
                    Vertex v;
                    if (heap.first_key() <= level) v = heap.pop();
                    else { v = order[cursor++]; require(!heap.contains(v), "implicit heap duplicate"); }
                    require(live[v], "duplicate factored removal"); live[v] = 0; --remaining;
                    upper[v] = level; next_order.push_back(v); batch.push_back(v); ++events;
                }
                if (!remaining) break;
                affected.clear(); changed.clear();
                for (Vertex v : batch) for (uint64_t code : index.touching(v)) {
                    const RowId p = code >> 2; const unsigned role = code & 3;
                    if (dead[p]) continue;
                    if (!touched[p]) { touched[p] = 1; affected.push_back(p); }
                    if (role == 0) dead[p] = 1;
                    else if (role == 1) { require(count[p] > 0, "pivot counter underflow"); --count[p]; }
                    else { auto& z = choices[index.rows[p].group]; require(z > 0, "choice counter underflow"); --z; }
                }
                auto subtract = [&](Vertex v, Count loss) {
                    if (!live[v]) return;
                    require(support[v] >= loss, "factored degree underflow"); support[v] -= loss;
                    if (!dirty[v]) { dirty[v] = 1; changed.push_back(v); }
                };
                auto scan = [&](Offset begin, Offset end, Count loss) { if (!loss) return; for (Offset i = begin; i < end; ++i) subtract(index.members[i], loss); };
                for (RowId p : affected) {
                    const auto& row = index.rows[p]; const Vertex g = row.group;
                    const Value value = dead[p] ? Value{} : TS::coefficients(index, p, s, count[p], g == absent ? 0 : choices[g], choose, metrics);
                    require(value.h <= wh[p] && value.q <= wp[p], "negative factored common loss");
                    const Count lh = wh[p] - value.h, lq = wp[p] - value.q;
                    wh[p] = value.h; wp[p] = value.q; touched[p] = 0; if (!value.h) dead[p] = 1;
                    scan(row.begin, row.hold_end, lh); scan(row.hold_end, row.pivot_end, lq);
                    if (g == absent) continue;
                    require(value.x <= wx[g], "negative choice loss");
                    const Count lx = wx[g] - value.x; wx[g] = value.x;
                    if (!lx) continue;
                    const Offset begin = choice_off[g]; Vertex length = choice_size[g], at = 0;
                    while (at < length) {
                        Vertex v = scratch[begin + at]; bool replacement = false;
                        while (!live[v]) { --length; if (at == length) break; v = scratch[begin + length]; replacement = true; }
                        if (at == length) break;
                        if (replacement) scratch[begin + at] = v;
                        subtract(v, lx); ++at;
                    }
                    choice_size[g] = length;
                }
                for (Vertex v : changed) {
                    dirty[v] = 0; const Count next = std::min(upper[v], support[v]);
                    if (heap.contains(v)) { if (next < heap.key(v)) heap.decrease(v, next); }
                    else if (support[v] < upper[v]) heap.insert(v, next);
                }
            }
            ph.peel_ms += ms(t0);
            require(next_order.size() == n, "incomplete factored order"); order.swap(next_order);
            sink(s, upper);
        }
        ph.total_ms = ms(start); return ph;
    }
};
}
