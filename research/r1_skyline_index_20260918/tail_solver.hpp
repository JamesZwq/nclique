#pragma once
// Tail-certified all-size solver (2026-09-23).
//
// Contract: the same as terminal::Solver<Count>::solve with a row sink.  Rows s = 2, 3, ... of the vertex s-clique core
// matrix are delivered one at a time in increasing s, the first all-zero row included, and the loop stops there.
//
// Idea.  At size s every active vertex v (clique number omega(v) >= s) has the clique floor C(omega-1, s-1) as a lower
// bound and the integer Kruskal-Katona bound from kappa_{s-1}(v) as an upper bound.  When the two coincide, kappa_s(v)
// is known before any peeling; the tail theorem then keeps it on the floor at every larger size ("settled").  Only the
// residue R_s (active, not settled) needs a peel, and a residue count only moves through rows that hold a residue
// member.  A size with a residue is peeled in one of two ways.
//  * Residue peel: over the rows valid at s that touch R_s and over their members (the relevant vertices).  A settled
//    relevant vertex is never keyed by a count; it leaves the stream exactly when the level reaches its known value
//    (asserted), and a row with no residue member is never initialised or updated.
//  * Full peel: the per-size peel of terminal::Solver, statement for statement, over every valid row and every active
//    vertex; settled vertices are peeled like the others and their result is asserted equal to the floor.  It is chosen
//    when the residue touches most of the valid rows: there the residue peel saves little and pays for the marking.
// Correctness of the residue peel: for every level k, {kappa_s >= k} = D_k is the largest set whose members each lie in
// >= k s-cliques inside it; removing residue vertices with fewer than k cliques inside (current residue) + {settled :
// kappa >= k} ends exactly at D_k minus the settled part, and a clique through a residue vertex lies on a relevant row,
// so every other vertex of that clique is relevant (RESULTS_FINAL 17.10).
//
// Tail = false puts every active vertex into the residue; with the default policy every size is then a full peel.
#include "../r1_terminal_20260918/terminal.hpp"
#include <span>
#include <type_traits>

namespace tailpeel {
using namespace allsize;

struct Policy {
    double dense = 0.30;     // residue touching volume / valid incidences above which the size takes the full peel
    double relevant = 0.60;  // relevant incidences / valid incidences above which the marking stops for the full peel
    int force = 0;           // 0 adaptive, 1 residue peel whenever there is a residue, 2 full peel whenever there is one
};
struct SizeStats {
    int s = 0, mode = 0;     // mode 0: closed form (no residue) or past the clique number, 1: residue peel, 2: full peel
    uint64_t active = 0, residue = 0, settled_new = 0, residue_degree = 0, valid_incidences = 0;
    uint64_t settled_early = 0;   // settled by the integer cascade although kappa_{s-1} is above the floor of size s-1
    uint64_t relevant_vertices = 0, relevant_rows = 0, relevant_incidences = 0;
    double ms = 0;
};
struct Stats {
    std::vector<SizeStats> sizes;
    double omega_ms = 0, upper_ms = 0, mark_ms = 0, init_ms = 0, peel_ms = 0, order_ms = 0, total_ms = 0;
    uint64_t events = 0, subtracts = 0, heap_pops = 0, last_residue_size = 0;
    uint64_t residue_sizes = 0, full_sizes = 0, aborted_marks = 0;
};

// What the tail solver reads besides the index: the clique number of every vertex (capped at the index maximum: the
// largest row through it, 0 for a vertex in no edge), the largest of them, and the total length of the rows valid at
// each size.
struct Prepared {
    Vertices omega;
    Vertex max_omega = 0;
    std::vector<uint64_t> valid_incidences;
};
inline void scan_rows(const terminal::Index& index, Prepared& out) {       // the row part of Prepared
    out.valid_incidences.assign(static_cast<size_t>(index.maximum) + 2, 0);
    for (const auto& row : index.rows) { const uint64_t len = row.end - row.begin; out.valid_incidences[row.lo] += len; out.valid_incidences[static_cast<size_t>(row.hi) + 1] -= len; }
    for (size_t i = 1; i < out.valid_incidences.size(); ++i) out.valid_incidences[i] += out.valid_incidences[i - 1];
    for (Vertex w : out.omega) out.max_omega = std::max(out.max_omega, w);
}
// Index::prepare (the reverse lists, identical) fused with Prepared: one pass over the members serves both.
inline Prepared prepare(terminal::Index& index, Vertex n) {
    using terminal::RowId; using terminal::Offset;
    index.reverse_off.assign(static_cast<size_t>(n) + 1, 0);
    for (Vertex v : index.members) ++index.reverse_off[v + 1];
    std::partial_sum(index.reverse_off.begin(), index.reverse_off.end(), index.reverse_off.begin());
    index.reverse.resize(index.members.size());
    auto cursor = index.reverse_off;
    Prepared out; out.omega.assign(n, 0);
    out.valid_incidences.assign(static_cast<size_t>(index.maximum) + 2, 0);
    for (RowId p = 0; p < index.rows.size(); ++p) {
        const auto& row = index.rows[p];
        for (Offset i = row.begin; i < row.end; ++i) {
            const Vertex u = index.members[i];
            const uint64_t role = i < row.hold_end ? 0 : (i < row.pivot_end ? 1 : 2);
            index.reverse[cursor[u]++] = static_cast<terminal::Code>((static_cast<uint64_t>(p) << 2) | role);
            out.omega[u] = std::max<Vertex>(out.omega[u], row.hi);
        }
        const uint64_t len = row.end - row.begin;
        out.valid_incidences[row.lo] += len; out.valid_incidences[static_cast<size_t>(row.hi) + 1] -= len;
        out.max_omega = std::max(out.max_omega, row.hi);
    }
    for (size_t i = 1; i < out.valid_incidences.size(); ++i) out.valid_incidences[i] += out.valid_incidences[i - 1];
    return out;
}

template<class Count, bool Tail = true> struct Solver {
    using Base = fullrange::Kernel<Count>;
    using Choose = typename Base::Combinations;
    using Queue = typename Base::Queue;
    using TS = terminal::Solver<Count>;
    using Value = typename TS::Value;

    template<class Sink> static Stats solve(const Graph& graph, const terminal::Index& index, const Choose& choose,
                                            const Vertices& ordinary, Sink sink, const Prepared* prepared = nullptr,
                                            const Policy policy = {}) {
        using terminal::RowId; using terminal::Offset;
        const auto t_all = Clock::now();
        Stats st;
        const Vertex n = graph.n; const int S = index.maximum;
        // clique numbers and valid-row volumes: from tailpeel::prepare when given, else one pass over the rows
        auto t0 = Clock::now();
        Prepared own;
        if (!prepared) {
            own.omega.assign(n, 0);
            for (const auto& row : index.rows)
                for (Offset i = row.begin; i < row.end; ++i) { Vertex& w = own.omega[index.members[i]]; w = std::max<Vertex>(w, row.hi); }
            scan_rows(index, own);
        }
        const Prepared& pre = prepared ? *prepared : own;
        const Vertices& omega = pre.omega; const Vertex max_omega = pre.max_omega; const auto& valid_incid = pre.valid_incidences;
        require(omega.size() == n && valid_incid.size() == static_cast<size_t>(S) + 2, "prepared data of another index");
        st.omega_ms = ms(t0);
        // one value row, updated in place: row s-1 at the start of size s, keys during the peel, row s at its end
        std::vector<Count> val(n);
        for (Vertex v = 0; v < n; ++v) val[v] = ordinary[v];
        sink(2, std::span<const Count>(val));
        // the stream order: the active vertices of the previous size, by nondecreasing value (inactive ones never return)
        Vertices order, next_order;
        for (Vertex v = 0; v < n; ++v) {
            require((ordinary[v] >= 1) == (omega[v] >= 2), "clique number disagrees with the ordinary core");
            if (ordinary[v]) order.push_back(v);
        }
        if (!std::is_sorted(order.begin(), order.end(), [&](Vertex a, Vertex b) { return ordinary[a] < ordinary[b]; }))
            std::stable_sort(order.begin(), order.end(), [&](Vertex a, Vertex b) { return ordinary[a] < ordinary[b]; });
        const Vertex maximum = ordinary.empty() ? 0 : *std::max_element(ordinary.begin(), ordinary.end());

        // state kept across sizes; only the touched entries are reset after a size
        const size_t R = index.rows.size(), G = index.group_row.size();
        std::vector<uint8_t> settled(n, 0), in_r(n, 0), live(n, 0), relevant(n, 0), dirty(n, 0);
        std::vector<Count> support(n, 0);
        std::vector<uint8_t> dead(R, 1), touched(R, 0), marked(R, 0);
        Vertices count(R, 0);
        std::vector<Count> wh(R, 0), wp(R, 0), wx(G, 0);
        Vertices choices(G, 0), choice_size(G, 0), scratch;
        std::vector<Offset> choice_off(G, 0);
        std::vector<RowId> rel_rows, affected;
        Vertices rel_verts, residue, batch, changed, removal;
        struct Cache { Count key = Base::infinity, value = 0; };
        std::vector<Cache> cache(4096);
        tworoads::Statistics bstats; orderdp::Extra extra; terminal::Metrics metrics;

        for (int s = 3; s <= S; ++s) {
            const auto t_size = Clock::now();
            SizeStats ss; ss.s = s; ss.valid_incidences = valid_incid[static_cast<size_t>(s)];
            const Vertex sv = static_cast<Vertex>(s);
            if (max_omega < sv) {                                      // first all-zero row: deliver it and stop
                std::fill(val.begin(), val.end(), Count{0});
                sink(s, std::span<const Count>(val));
                st.sizes.push_back(ss); break;
            }
            // values known before peeling: settled vertices, and residue vertices with their upper bound.  The order
            // loses the vertices that stop being active (their value becomes 0) and keeps its sort: the key is a
            // nondecreasing function of the previous value.
            t0 = Clock::now();
            std::fill(cache.begin(), cache.end(), Cache{});
            residue.clear();
            uint64_t residue_degree = 0;
            {
                Count* const value = val.data(); const Vertex* const om = omega.data();
                uint8_t* const setl = settled.data(); uint8_t* const res = in_r.data();
                Vertex* const ord = order.data(); const size_t size = order.size(); size_t kept = 0;
                for (size_t i = 0; i < size; ++i) {
                    const Vertex v = ord[i];
                    if (om[v] < sv) { value[v] = 0; continue; }
                    ord[kept++] = v;
                    const Count floor_value = choose(static_cast<int>(om[v]) - 1, s - 1);
                    if (setl[v]) { value[v] = floor_value; continue; }
                    const Count a = value[v];
                    require(a > 0, "active vertex without a value one size below");
                    auto& entry = cache[static_cast<size_t>((a ^ (a >> 17) ^ (a >> 37)) & (cache.size() - 1))];
                    if (entry.key != a) entry = {a, Base::integer_upper(a, s - 2, maximum, bstats)};
                    const Count upper = entry.value;
                    require(upper >= floor_value, "upper bound below the clique floor");
                    if (Tail && upper == floor_value) {
                        setl[v] = 1; value[v] = floor_value; ++ss.settled_new;
                        if (s > 3 && a != choose(static_cast<int>(om[v]) - 1, s - 2)) ++ss.settled_early;
                    }
                    else { value[v] = upper; res[v] = 1; residue.push_back(v); residue_degree += index.reverse_off[v + 1] - index.reverse_off[v]; }
                }
                order.resize(kept); ss.active = kept;
            }
            ss.residue = residue.size(); ss.residue_degree = residue_degree;
            st.upper_ms += ms(t0);
            if (residue.empty()) {                                     // closed form, and so is every larger size
                sink(s, std::span<const Count>(val));
                ss.ms = ms(t_size); st.sizes.push_back(ss); continue;
            }
            st.last_residue_size = static_cast<uint64_t>(s);
            // the relevant rows (valid at s, holding a residue member), unless the residue is too wide to pay for them
            t0 = Clock::now();
            bool full = policy.force == 2 || (policy.force == 0 && static_cast<double>(residue_degree) > policy.dense * static_cast<double>(ss.valid_incidences));
            if (!full) {
                const double cap = policy.force == 1 ? std::numeric_limits<double>::infinity() : policy.relevant * static_cast<double>(ss.valid_incidences);
                uint8_t* const mark = marked.data(); const terminal::Row* const rows = index.rows.data();
                uint64_t incid = 0;
                for (Vertex v : residue) {
                    for (uint64_t code : index.touching(v)) {
                        const RowId p = code >> 2;
                        if (mark[p]) continue;
                        const auto& row = rows[p];
                        if (!row.valid(s)) continue;
                        mark[p] = 1; rel_rows.push_back(p); incid += row.end - row.begin;
                    }
                    if (static_cast<double>(incid) > cap) break;
                }
                if (static_cast<double>(incid) > cap) {                // too wide: undo the marks and peel in full
                    for (RowId p : rel_rows) mark[p] = 0;
                    rel_rows.clear(); full = true; ++st.aborted_marks;
                } else { std::sort(rel_rows.begin(), rel_rows.end()); ss.relevant_rows = rel_rows.size(); ss.relevant_incidences = incid; }
            }
            ss.mode = full ? 2 : 1; ++(full ? st.full_sizes : st.residue_sizes);
            st.mark_ms += ms(t0);
            // initial counts: every valid row and every active vertex (full), or the relevant rows with counts kept for
            // the residue members only and every member collected as relevant (residue)
            t0 = Clock::now();
            scratch.clear();
            Vertex remaining = 0;
            {
                Count* const sup = support.data(); uint8_t* const alive = live.data(); uint8_t* const res = in_r.data();
                uint8_t* const rel = relevant.data(); uint8_t* const gone = dead.data();
                Vertex* const cnt = count.data(); Count* const hw = wh.data(); Count* const pw = wp.data();
                const Vertex* const members = index.members.data(); const terminal::Row* const rows = index.rows.data();
                auto open_row = [&](RowId p) -> Value {
                    const auto& row = rows[p];
                    gone[p] = 0; cnt[p] = row.pivots();
                    const Vertex z = static_cast<Vertex>(row.end - row.pivot_end);
                    const Value value = TS::coefficients(index, p, s, cnt[p], z, choose, metrics);
                    hw[p] = value.h; pw[p] = value.q;
                    if (row.group != absent) {
                        const Vertex g = row.group;
                        wx[g] = value.x; choices[g] = z; choice_off[g] = scratch.size(); choice_size[g] = 0;
                        if (value.x) { choice_size[g] = z; scratch.insert(scratch.end(), index.members.begin() + row.pivot_end, index.members.begin() + row.end); }
                    }
                    return value;
                };
                if (full) {
                    auto add = [&](Offset b, Offset e, Count w) { if (!w) return; for (Offset i = b; i < e; ++i) Base::checked_add(sup[members[i]], w); };
                    for (RowId p = 0; p < R; ++p) {
                        const auto& row = rows[p];
                        if (!row.valid(s)) continue;
                        const Value value = open_row(p);
                        add(row.begin, row.hold_end, value.h); add(row.hold_end, row.pivot_end, value.q);
                        if (row.group != absent) add(row.pivot_end, row.end, value.x);
                    }
                    for (Vertex v : order) alive[v] = 1;
                    remaining = static_cast<Vertex>(order.size());
                } else {
                    auto add = [&](Offset b, Offset e, Count w) {
                        for (Offset i = b; i < e; ++i) {
                            const Vertex u = members[i];
                            if (!rel[u]) { rel[u] = 1; alive[u] = 1; rel_verts.push_back(u); }
                            if (res[u]) Base::checked_add(sup[u], w);
                        }
                    };
                    for (RowId p : rel_rows) {
                        const auto& row = rows[p];
                        const Value value = open_row(p);
                        add(row.begin, row.hold_end, value.h); add(row.hold_end, row.pivot_end, value.q);
                        if (row.group != absent) add(row.pivot_end, row.end, value.x);
                    }
                    remaining = static_cast<Vertex>(rel_verts.size());
                    ss.relevant_vertices = rel_verts.size();
                }
                for (Vertex v : residue) require(sup[v] > 0, "residue vertex without a clique");
            }
            st.init_ms += ms(t0);
            // the level-batch peel; keys are min(count, val) in the heap, val in the stream (the order).
            // Residue peel: settled vertices carry no count and leave the stream at their known value.
            t0 = Clock::now();
            removal.clear();
            auto peel = [&](auto residue_tag) {
                constexpr bool Residue = decltype(residue_tag)::value;
                // local pointers and counters: byte stores in the loop cannot then force reloads of the arrays' bases
                Count* const sup = support.data(); Count* const key = val.data();
                uint8_t* const alive = live.data(); uint8_t* const res = in_r.data(); uint8_t* const dirt = dirty.data();
                uint8_t* const gone = dead.data(); uint8_t* const seen = touched.data();
                Vertex* const cnt = count.data(); Vertex* const ch = choices.data(); Vertex* const chsize = choice_size.data();
                Count* const hw = wh.data(); Count* const pw = wp.data(); Count* const xw = wx.data();
                Vertex* const pool = scratch.data(); const Offset* const choff = choice_off.data();
                const Vertex* const members = index.members.data(); const terminal::Row* const rows = index.rows.data();
                const Vertex* const ord = order.data(); const size_t order_size = order.size();
                uint64_t events = 0, subtracts = 0, pops = 0;
                Queue heap(std::span<const Count>(support), std::span<const Count>(val), false, extra);
                size_t cursor = 0;
                auto stream_key = [&]() -> Count {
                    while (cursor < order_size && !alive[ord[cursor]]) ++cursor;
                    return cursor < order_size ? key[ord[cursor]] : Base::infinity;
                };
                Count level = 0;
                while (remaining) {
                    level = std::max(level, std::min(stream_key(), heap.first_key()));
                    require(level != Base::infinity, "missing minimum");
                    batch.clear();
                    while (std::min(stream_key(), heap.first_key()) <= level) {
                        Vertex v;
                        if (heap.first_key() <= level) { v = heap.pop(); ++pops; }
                        else { v = ord[cursor++]; require(!heap.contains(v), "implicit heap duplicate"); }
                        require(alive[v], "duplicate removal"); alive[v] = 0; --remaining;
                        if (!res[v]) require(key[v] == level, "a settled vertex left at a level other than its value");
                        key[v] = level; removal.push_back(v); batch.push_back(v); ++events;
                    }
                    if (!remaining) break;
                    affected.clear(); changed.clear();
                    for (Vertex v : batch)
                        for (uint64_t code : index.touching(v)) {
                            const RowId p = code >> 2; const unsigned role = code & 3;
                            if (gone[p]) continue;
                            if (!seen[p]) { seen[p] = 1; affected.push_back(p); }
                            if (role == 0) gone[p] = 1;
                            else if (role == 1) { require(cnt[p] > 0, "pivot counter underflow"); --cnt[p]; }
                            else { auto& z = ch[rows[p].group]; require(z > 0, "choice counter underflow"); --z; }
                        }
                    auto subtract = [&](Vertex v, Count loss) {
                        if (!alive[v]) return;
                        if constexpr (Residue) { if (!res[v]) return; }
                        require(sup[v] >= loss, "count underflow"); sup[v] -= loss; ++subtracts;
                        if (!dirt[v]) { dirt[v] = 1; changed.push_back(v); }
                    };
                    auto scan = [&](Offset b, Offset e, Count loss) { if (!loss) return; for (Offset i = b; i < e; ++i) subtract(members[i], loss); };
                    for (RowId p : affected) {
                        const auto& row = rows[p]; const Vertex g = row.group;
                        const Value value = gone[p] ? Value{} : TS::coefficients(index, p, s, cnt[p], g == absent ? 0 : ch[g], choose, metrics);
                        require(value.h <= hw[p] && value.q <= pw[p], "negative loss");
                        const Count lh = hw[p] - value.h, lq = pw[p] - value.q;
                        hw[p] = value.h; pw[p] = value.q; seen[p] = 0; if (!value.h) gone[p] = 1;
                        scan(row.begin, row.hold_end, lh); scan(row.hold_end, row.pivot_end, lq);
                        if (g == absent) continue;
                        require(value.x <= xw[g], "negative choice loss");
                        const Count lx = xw[g] - value.x; xw[g] = value.x;
                        if (!lx) continue;
                        const Offset begin = choff[g]; Vertex length = chsize[g], at = 0;
                        while (at < length) {
                            Vertex v = pool[begin + at]; bool replacement = false;
                            while (!alive[v]) { --length; if (at == length) break; v = pool[begin + length]; replacement = true; }
                            if (at == length) break;
                            if (replacement) pool[begin + at] = v;
                            subtract(v, lx); ++at;
                        }
                        chsize[g] = length;
                    }
                    for (Vertex v : changed) {
                        dirt[v] = 0; const Count next = std::min(key[v], sup[v]);
                        if (heap.contains(v)) { if (next < heap.key(v)) heap.decrease(v, next); }
                        else if (sup[v] < key[v]) heap.insert(v, next);
                    }
                }
                st.events += events; st.subtracts += subtracts; st.heap_pops += pops;
            };
            if (full) peel(std::false_type{}); else peel(std::true_type{});
            st.peel_ms += ms(t0);
            // the order of the next size: the removal sequence (nondecreasing levels), merged in the residue peel with
            // the order restricted to the non-relevant vertices (settled: their value is a nondecreasing function of the
            // previous one, so the old order sorts them)
            t0 = Clock::now();
            if (full) order.swap(removal);
            else {
                next_order.clear();
                size_t j = 0;
                for (Vertex v : order) {
                    if (relevant[v]) continue;
                    while (j < removal.size() && val[removal[j]] < val[v]) next_order.push_back(removal[j++]);
                    next_order.push_back(v);
                }
                while (j < removal.size()) next_order.push_back(removal[j++]);
                require(next_order.size() == order.size(), "incomplete order");
                order.swap(next_order);
            }
            st.order_ms += ms(t0);
            // reset the touched state
            if (full) {
                std::fill(dead.begin(), dead.end(), uint8_t{1});
                for (Vertex v : order) support[v] = 0;
            } else {
                for (RowId p : rel_rows) { marked[p] = 0; dead[p] = 1; touched[p] = 0; }
                rel_rows.clear();
                for (Vertex v : rel_verts) { relevant[v] = 0; support[v] = 0; live[v] = 0; dirty[v] = 0; }
                rel_verts.clear();
            }
            for (Vertex v : residue) in_r[v] = 0;
            sink(s, std::span<const Count>(val));
            ss.ms = ms(t_size); st.sizes.push_back(ss);
        }
        st.total_ms = ms(t_all);
        return st;
    }
};
}
