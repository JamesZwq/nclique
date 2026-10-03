// Clique-tree census (2026-10-04): runs the terminal builder with a counting sink instead of the row index, so the size
// of the clique tree of a graph too large to materialise can be measured.  Per builder mode: rows, member incidences
// (hold / pivot / choice), groups, builder states, the rows and incidences valid at each size, and the bytes the
// narrow and wide row indexes would take.  Usage: tree_census <graph.edges> [mode ...]
#include "common.hpp"
#include <iostream>
using namespace orderreplay;
#include "../r1_orderreplay_20260917/shared_harness.inc"

struct Census {
    int maximum;
    terminal::BuildWork work;
    uint64_t rows = 0, groups = 0, hold = 0, pivot = 0, choice = 0, max_row = 0;
    std::vector<uint64_t> valid_rows, valid_inc;      // difference arrays over sizes
    explicit Census(int s): maximum(s), valid_rows(s + 3, 0), valid_inc(s + 3, 0) {}
    void append(const Vertices& h, const Vertices& q, const Vertices& x = {}, bool zero = false) {
        require(!h.empty(), "terminal requires a hold");
        const Vertex low = std::max<size_t>(2, h.size() + (!x.empty() && !zero));
        const Vertex high = std::min<size_t>(maximum, h.size() + q.size() + !x.empty());
        require(low <= high, "invalid terminal size interval");
        if (!x.empty()) { ++groups; work.choices += x.size(); }
        ++rows; hold += h.size(); pivot += q.size(); choice += x.size();
        const uint64_t len = h.size() + q.size() + x.size(); max_row = std::max(max_row, len);
        valid_rows[low] += 1; valid_rows[high + 1] -= 1; valid_inc[low] += len; valid_inc[high + 1] -= len;
    }
};

template<int Mode> static void run(const Input& in, int mode) {
    const int S = std::max(2, static_cast<int>(in.d) + 1);
    Census c(S); const auto t = Clock::now();
    terminal::Builder<Mode, Census>(in.graph, c).run();
    const double build_ms = ms(t);
    const uint64_t inc = c.hold + c.pivot + c.choice;
    // narrow row: 4 offsets + group, lo, hi (4 bytes each); wide: 8-byte offsets.  Members 4 bytes, reverse one code each.
    const uint64_t narrow = c.rows * 28 + inc * 8 + c.groups * 5, wide = c.rows * 48 + inc * 12 + c.groups * 9;
    std::cout << "{\"mode\":" << mode << ",\"n\":" << in.graph.n << ",\"m\":" << in.graph.m << ",\"degeneracy\":" << in.d
              << ",\"rows\":" << c.rows << ",\"groups\":" << c.groups << ",\"hold\":" << c.hold << ",\"pivot\":" << c.pivot
              << ",\"choice\":" << c.choice << ",\"incidences\":" << inc << ",\"max_row\":" << c.max_row
              << ",\"states\":" << c.work.states << ",\"narrow_bytes\":" << narrow << ",\"wide_bytes\":" << wide
              << ",\"build_ms\":" << build_ms << ",\"valid\":[";
    uint64_t r = 0, i = 0; bool first = true;
    for (int s = 2; s <= S; ++s) { r += c.valid_rows[s]; i += c.valid_inc[s];
        if (s == 2) { r = 0; i = 0; for (int t2 = 0; t2 <= 2; ++t2) { r += c.valid_rows[t2]; i += c.valid_inc[t2]; } }
        if (!r) break;
        if (s > 12 && s % 10 && s != S) continue;              // sizes 2..12, then every tenth
        std::cout << (first ? "" : ",") << "[" << s << "," << r << "," << i << "]"; first = false; }
    std::cout << "]}" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: tree_census <graph.edges> [mode ...]\n"; return 2; }
    const auto t = Clock::now(); const Input in = prepare(argv[1]);
    std::cerr << "loaded n=" << in.graph.n << " m=" << in.graph.m << " d=" << in.d << " in " << ms(t) << " ms\n";
    std::vector<int> modes; for (int a = 2; a < argc; ++a) modes.push_back(std::stoi(argv[a]));
    if (modes.empty()) modes = {0, 2};
    for (int mode : modes) { if (mode == 0) run<0>(in, 0); else if (mode == 1) run<1>(in, 1); else run<2>(in, 2); }
}
