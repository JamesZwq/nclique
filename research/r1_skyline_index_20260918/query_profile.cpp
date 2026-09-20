// Query-latency profile of a chain index file, against the S trees baseline built in the same process from the same
// decomposition (one tree and one DFS array per size over vertices, parent-pointer climb, a community is one memory copy;
// the layout of stages/index.cpp `vertices` mode and of the paper's baseline).  Workloads: the fixed one (seed 20260918:
// 20,000 own-level, 20,000 half-level and 1,000 root-level community queries over active vertices, s uniform in
// [2, omega(v)]), then (a) the own-level queries split by answer size into deciles, (b) a stratified workload (seed
// 20260921): for every clique size s in a fixed list, 5,000 own-level queries with v uniform among the vertices active at
// s.  Every latency is the per-query average of the median of five passes after one warm-up pass, one thread; the chain
// index is the loaded form (the form the tool measures).  Usage: query_profile <index.cx>   -> one JSON line.
#include "../../src-r1index/chain_index.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <random>
using namespace chainindex;
using Index = ChainIndex<double>; using Clock = std::chrono::steady_clock;
struct Q { uint32_t v, u; int s; double k; };

// ---- the baseline: S trees over vertices, materialized from the index (same nodes, same DFS order, vertex labels)
struct STrees {
    struct Size { std::vector<double> top; std::vector<uint32_t> parent, lo, hi, dfs; };
    std::vector<Size> sizes;                 // indexed by s
    std::vector<uint32_t> own_off, own;      // per vertex: own node at sizes 2 .. omega(v), ragged
    uint64_t bytes = 0;
    static STrees build(const Index& ix) {
        STrees st; const uint32_t n = ix.n; st.sizes.resize(ix.max_size + 1);
        for (int s = 2; s <= ix.max_size; ++s) {
            const Index::Layer& L = ix.layers[s]; const uint32_t N = static_cast<uint32_t>(L.size.size()); const size_t R = L.runs.size() / 2; Size& S = st.sizes[s];
            std::vector<uint64_t> pre(R + 1, 0); for (size_t j = 0; j < R; ++j) pre[j + 1] = pre[j] + (L.runs[2 * j + 1] - L.runs[2 * j]);
            auto pos = [&](uint32_t r, uint32_t v) { return r < R ? pre[r] + (v - L.runs[2 * static_cast<size_t>(r)]) : pre[R]; };
            S.top.resize(N); S.parent.resize(N); S.lo.resize(N); S.hi.resize(N);
            for (uint32_t x = 0; x < N; ++x) { const size_t e = static_cast<size_t>(x) + L.size[x];
                S.top[x] = ix.top_at(L, x); S.parent[x] = L.parent[x];
                S.lo[x] = static_cast<uint32_t>(pos(L.entry[2 * static_cast<size_t>(x)], L.entry[2 * static_cast<size_t>(x) + 1]));
                S.hi[x] = static_cast<uint32_t>(pos(L.entry[2 * e], L.entry[2 * e + 1])); }
            for (size_t j = 0; j < R; ++j) for (uint32_t v = L.runs[2 * j]; v < L.runs[2 * j + 1]; ++v) S.dfs.push_back(v);
            st.bytes += 20ull * N + 4ull * S.dfs.size();
        }
        st.own_off.assign(n + 1, 0);
        for (uint32_t v = 0; v < n; ++v) { const uint32_t c = ix.chain_of(v); st.own_off[v + 1] = st.own_off[v] + (ix.omega[c] >= 2 ? ix.omega[c] - 1 : 0); }
        st.own.resize(st.own_off[n]);
        for (uint32_t v = 0; v < n; ++v) { const uint32_t c = ix.chain_of(v); for (int s = 2; s <= ix.omega[c]; ++s) st.own[st.own_off[v] + s - 2] = ix.own_node(c, s); }
        st.bytes += 4ull * (n + 1) + 4ull * st.own.size();
        return st;
    }
    inline uint32_t locate(uint32_t v, int s, double k, uint32_t& lo, uint32_t& hi) const {
        const Size& S = sizes[s]; uint32_t x = own[own_off[v] + s - 2];
        while (S.parent[x] != kNone && S.top[S.parent[x]] >= k) x = S.parent[x];
        lo = S.lo[x]; hi = S.hi[x]; return x;
    }
};

static double median5(std::array<double, 5> t) { std::sort(t.begin(), t.end()); return t[2]; }
struct Timed { double locate_ns, list_ns, st_locate_ns, st_list_ns; double output, ranges; };
static Timed measure(const Index& ix, const STrees& st, const std::vector<Q>& qs, std::vector<uint32_t>& ids) {
    std::array<double, 5> tl{}, te{}, sl{}, se{}; uint64_t z = 0, outputs = 0, nr = 0, st_outputs = 0;
    for (int pass = 0; pass < 6; ++pass) { const auto b = Clock::now();
        for (const auto& q : qs) { Index::Runs r; uint32_t nd = 0; ix.community_runs(q.v, q.s, q.k, r, nd); z += r.nmid + r.lo0 + nd; }
        const double el = std::chrono::duration<double, std::nano>(Clock::now() - b).count() / qs.size(); if (pass) tl[pass - 1] = el; }
    for (int pass = 0; pass < 6; ++pass) { uint64_t cnt = 0; const auto b = Clock::now();
        for (const auto& q : qs) { Index::Runs r; uint32_t nd = 0; ix.community_runs(q.v, q.s, q.k, r, nd); cnt += static_cast<uint64_t>(Index::expand(r, ids.data()) - ids.data()); if (pass == 0) nr += r.count(); }
        const double el = std::chrono::duration<double, std::nano>(Clock::now() - b).count() / qs.size(); if (pass == 0) outputs = cnt; else te[pass - 1] = el; }
    for (int pass = 0; pass < 6; ++pass) { const auto b = Clock::now();
        for (const auto& q : qs) { uint32_t lo, hi; z += st.locate(q.v, q.s, q.k, lo, hi) + lo + hi; }
        const double el = std::chrono::duration<double, std::nano>(Clock::now() - b).count() / qs.size(); if (pass) sl[pass - 1] = el; }
    for (int pass = 0; pass < 6; ++pass) { uint64_t cnt = 0; const auto b = Clock::now();
        for (const auto& q : qs) { uint32_t lo, hi; st.locate(q.v, q.s, q.k, lo, hi); std::memcpy(ids.data(), st.sizes[q.s].dfs.data() + lo, (hi - lo) * sizeof(uint32_t)); cnt += hi - lo; }
        const double el = std::chrono::duration<double, std::nano>(Clock::now() - b).count() / qs.size(); if (pass == 0) st_outputs = cnt; else se[pass - 1] = el; }
    if (z == 42) std::cerr << "";   // keep the locate loops observable
    if (st_outputs != outputs) { std::cerr << "S trees and index disagree on the answer size\n"; std::exit(1); }
    return {median5(tl), median5(te), median5(sl), median5(se), double(outputs) / qs.size(), double(nr) / qs.size()};
}
static uint64_t output_of(const Index& ix, const Q& q) { Index::Runs r; uint32_t nd = 0; ix.community_runs(q.v, q.s, q.k, r, nd); return Index::total(r); }
static void put(const Timed& t) { std::cout << "\"locate_ns\":" << t.locate_ns << ",\"list_ns\":" << t.list_ns << ",\"st_locate_ns\":" << t.st_locate_ns << ",\"st_list_ns\":" << t.st_list_ns << ",\"output\":" << t.output << ",\"ranges\":" << t.ranges; }

int main(int argc, char** argv) {
    if (argc != 2) { std::cerr << "usage: query_profile <index.cx>\n"; return 2; }
    const Index ix = Index::load(argv[1]); const uint32_t n = ix.n; const STrees st = STrees::build(ix);
    std::vector<uint32_t> active; for (uint32_t v = 0; v < n; ++v) if (ix.omega[ix.chain_of(v)] >= 2) active.push_back(v);
    std::vector<uint32_t> ids(static_cast<size_t>(n) + Index::kSlack);
    std::mt19937_64 rng(20260918);
    auto draw = [&](int regime, int count) { std::vector<Q> qs; for (int i = 0; i < count; ++i) { const uint32_t v = active[rng() % active.size()]; const uint32_t c = ix.chain_of(v);
        const int s = 2 + static_cast<int>(rng() % static_cast<uint64_t>(ix.omega[c] - 1)); const double x = ix.value(v, s);
        const double k = regime == 0 ? x : (regime == 1 ? std::max<double>(1.0, std::floor(x / 2)) : 1.0); qs.push_back({v, static_cast<uint32_t>(rng() % n), s, k}); } return qs; };
    const std::vector<Q> own = draw(0, 20000), half = draw(1, 20000), root = draw(2, 1000);
    std::cout << std::fixed << std::setprecision(3) << "{\"n\":" << n << ",\"s_max\":" << ix.max_size << ",\"active\":" << active.size() << ",\"index_bytes\":" << ix.bytes_total() << ",\"strees_struct_bytes\":" << st.bytes << ",\"regimes\":{";
    const char* names[3] = {"own", "half", "root"}; const std::vector<Q>* sets[3] = {&own, &half, &root};
    for (int r = 0; r < 3; ++r) { const Timed t = measure(ix, st, *sets[r], ids); std::cout << (r ? "," : "") << "\"" << names[r] << "\":{\"queries\":" << sets[r]->size() << ","; put(t); std::cout << "}"; }
    // (a) own-level queries by answer size: ten deciles of the 20,000 queries sorted by community size
    std::vector<Q> sorted = own; std::stable_sort(sorted.begin(), sorted.end(), [&](const Q& a, const Q& b) { return output_of(ix, a) < output_of(ix, b); });
    std::cout << "},\"own_by_output\":[";
    for (int d = 0; d < 10; ++d) { std::vector<Q> part(sorted.begin() + sorted.size() * d / 10, sorted.begin() + sorted.size() * (d + 1) / 10);
        const Timed t = measure(ix, st, part, ids); const uint64_t lo = output_of(ix, part.front()), hi = output_of(ix, part.back());
        std::cout << (d ? "," : "") << "{\"decile\":" << d << ",\"queries\":" << part.size() << ",\"output_min\":" << lo << ",\"output_max\":" << hi << ","; put(t); std::cout << "}"; }
    // (b) stratified by clique size: own level, v uniform among the vertices active at s
    std::cout << "],\"by_size\":[";
    std::mt19937_64 rng2(20260921); bool first = true;
    const int sizes[] = {2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128, 160, 192, 256, 320, 384, 448, 500};
    for (int s : sizes) { if (s > ix.max_size) break;
        std::vector<uint32_t> act; for (uint32_t v : active) if (ix.omega[ix.chain_of(v)] >= s) act.push_back(v);
        if (act.empty()) continue;
        std::vector<Q> qs; for (int i = 0; i < 5000; ++i) { const uint32_t v = act[rng2() % act.size()]; qs.push_back({v, 0, s, ix.value(v, s)}); }
        const Timed t = measure(ix, st, qs, ids); const Index::Layer& L = ix.layers[s];
        std::cout << (first ? "" : ",") << "{\"s\":" << s << ",\"active\":" << act.size() << ",\"nodes\":" << L.size.size() << ",\"runs\":" << L.runs.size() / 2 << ",\"queries\":" << qs.size() << ","; put(t); std::cout << "}"; first = false; }
    std::cout << "]}\n";
    return 0;
}
