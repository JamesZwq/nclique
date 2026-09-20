// Query-latency profile of a chain index file: the fixed workload (seed 20260918: 20,000 own-level, 20,000 half-level and
// 1,000 root-level community queries over active vertices, s uniform in [2, omega(v)]) reproduced from the loaded index,
// then (a) the own-level queries split by output size into deciles, (b) a stratified workload (seed 20260921): for every
// clique size s in a fixed bucket list, 5,000 own-level queries with v uniform among the vertices active at s.  Every
// latency is the per-query average of the median of five passes after one warm-up pass, one thread, on the loaded index
// (the same form the tool measures).  Usage: query_profile <index.cx>   -> one JSON line.
#include "../../src-r1index/chain_index.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
using namespace chainindex;
using Index = ChainIndex<double>; using Clock = std::chrono::steady_clock;
struct Q { uint32_t v, u; int s; double k; };

static double median5(std::array<double, 5> t) { std::sort(t.begin(), t.end()); return t[2]; }
struct Timed { double locate_ns, list_ns; double output, ranges; };
static Timed measure(const Index& ix, const std::vector<Q>& qs, std::vector<uint32_t>& ids) {
    std::array<double, 5> tl{}, te{}; uint64_t z = 0, outputs = 0, nr = 0;
    for (int pass = 0; pass < 6; ++pass) { const auto st = Clock::now();
        for (const auto& q : qs) { Index::Runs r; uint32_t nd = 0; ix.community_runs(q.v, q.s, q.k, r, nd); z += r.nmid + r.lo0 + nd; }
        const double el = std::chrono::duration<double, std::nano>(Clock::now() - st).count() / qs.size(); if (pass) tl[pass - 1] = el; }
    for (int pass = 0; pass < 6; ++pass) { uint64_t cnt = 0; const auto st = Clock::now();
        for (const auto& q : qs) { Index::Runs r; uint32_t nd = 0; ix.community_runs(q.v, q.s, q.k, r, nd); cnt += static_cast<uint64_t>(Index::expand(r, ids.data()) - ids.data()); if (pass == 0) nr += r.count(); }
        const double el = std::chrono::duration<double, std::nano>(Clock::now() - st).count() / qs.size(); if (pass == 0) outputs = cnt; else te[pass - 1] = el; }
    if (z == 42) std::cerr << "";   // keep the locate loop observable
    return {median5(tl), median5(te), double(outputs) / qs.size(), double(nr) / qs.size()};
}
static uint64_t output_of(const Index& ix, const Q& q) { Index::Runs r; uint32_t nd = 0; ix.community_runs(q.v, q.s, q.k, r, nd); return Index::total(r); }

int main(int argc, char** argv) {
    if (argc != 2) { std::cerr << "usage: query_profile <index.cx>\n"; return 2; }
    const Index ix = Index::load(argv[1]); const uint32_t n = ix.n;
    std::vector<uint32_t> active; for (uint32_t v = 0; v < n; ++v) if (ix.omega[ix.chain_of(v)] >= 2) active.push_back(v);
    std::vector<uint32_t> ids(static_cast<size_t>(n) + Index::kSlack);
    std::mt19937_64 rng(20260918);
    auto draw = [&](int regime, int count) { std::vector<Q> qs; for (int i = 0; i < count; ++i) { const uint32_t v = active[rng() % active.size()]; const uint32_t c = ix.chain_of(v);
        const int s = 2 + static_cast<int>(rng() % static_cast<uint64_t>(ix.omega[c] - 1)); const double x = ix.value(v, s);
        const double k = regime == 0 ? x : (regime == 1 ? std::max<double>(1.0, std::floor(x / 2)) : 1.0); qs.push_back({v, static_cast<uint32_t>(rng() % n), s, k}); } return qs; };
    const std::vector<Q> own = draw(0, 20000), half = draw(1, 20000), root = draw(2, 1000);
    std::cout << std::fixed << std::setprecision(3) << "{\"n\":" << n << ",\"s_max\":" << ix.max_size << ",\"active\":" << active.size() << ",\"regimes\":{";
    const char* names[3] = {"own", "half", "root"}; const std::vector<Q>* sets[3] = {&own, &half, &root};
    for (int r = 0; r < 3; ++r) { const Timed t = measure(ix, *sets[r], ids);
        std::cout << (r ? "," : "") << "\"" << names[r] << "\":{\"queries\":" << sets[r]->size() << ",\"locate_ns\":" << t.locate_ns << ",\"list_ns\":" << t.list_ns << ",\"output\":" << t.output << ",\"ranges\":" << t.ranges << "}"; }
    // (a) own-level queries by output size: ten deciles of the 20,000 queries sorted by community size
    std::vector<Q> sorted = own; std::stable_sort(sorted.begin(), sorted.end(), [&](const Q& a, const Q& b) { return output_of(ix, a) < output_of(ix, b); });
    std::cout << "},\"own_by_output\":[";
    for (int d = 0; d < 10; ++d) { std::vector<Q> part(sorted.begin() + sorted.size() * d / 10, sorted.begin() + sorted.size() * (d + 1) / 10);
        const Timed t = measure(ix, part, ids); uint64_t lo = output_of(ix, part.front()), hi = output_of(ix, part.back());
        std::cout << (d ? "," : "") << "{\"decile\":" << d << ",\"queries\":" << part.size() << ",\"output_min\":" << lo << ",\"output_max\":" << hi << ",\"output\":" << t.output << ",\"ranges\":" << t.ranges << ",\"locate_ns\":" << t.locate_ns << ",\"list_ns\":" << t.list_ns << "}"; }
    // (b) stratified by clique size: own level, v uniform among the vertices active at s
    std::cout << "],\"by_size\":[";
    std::mt19937_64 rng2(20260921); bool first = true;
    const int sizes[] = {2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128, 160, 192, 256, 320, 384, 448, 500};
    for (int s : sizes) { if (s > ix.max_size) break;
        std::vector<uint32_t> act; for (uint32_t v : active) if (ix.omega[ix.chain_of(v)] >= s) act.push_back(v);
        if (act.empty()) continue;
        std::vector<Q> qs; for (int i = 0; i < 5000; ++i) { const uint32_t v = act[rng2() % act.size()]; qs.push_back({v, 0, s, ix.value(v, s)}); }
        const Timed t = measure(ix, qs, ids); const Index::Layer& L = ix.layers[s];
        std::cout << (first ? "" : ",") << "{\"s\":" << s << ",\"active\":" << act.size() << ",\"nodes\":" << L.size.size() << ",\"runs\":" << L.runs.size() / 2 << ",\"queries\":" << qs.size() << ",\"output\":" << t.output << ",\"ranges\":" << t.ranges << ",\"locate_ns\":" << t.locate_ns << ",\"list_ns\":" << t.list_ns << "}"; first = false; }
    std::cout << "]}\n";
    return 0;
}
