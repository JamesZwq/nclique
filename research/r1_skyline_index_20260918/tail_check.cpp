// tail_check: the tail-certified solver (tail_solver.hpp) against the frozen terminal solver.
//   --random <count>   random, split, complete and planted-clique graphs; every row compared value by value for
//                      terminal::Solver and four tail-solver configurations: adaptive, forced residue peel, forced full
//                      peel, and Tail=false
//   --graph <path>     one real graph: every row hashed and compared across the five, their solve times, and the
//                      per-size statistics of the adaptive and the two forced runs (one JSON line)
//   --time <solver> <graph> <repeats>   one solver, repeated (profiling aid)
// Timing runs are sequential in one process, each solver once, in the order listed.  TAIL_DENSE / TAIL_RELEVANT
// override the adaptive thresholds.
#define main skyline_stage1_count_main
#include "count.cpp"
#undef main
#include "tail_solver.hpp"
#include "terminal_probe.hpp"
#include <iostream>
#include <array>
#include <cstdlib>

static uint64_t cov_sizes = 0, cov_residue_sizes = 0, cov_mixed_sizes = 0, cov_closed_sizes = 0, cov_settled_new = 0;   // coverage of the random suite
static cpp_int bound_rows(const terminal::Index& ti, uint32_t n, int d) {   // as count_bound_terminal in chain_index_tool
    int qmax = d + 1; for (const auto& row : ti.rows) qmax = std::max(qmax, static_cast<int>(row.end - row.hold_end));   // rows truncated at the size bound
    std::vector<cpp_int> peaks(qmax + 1); for (int q = 0; q <= qmax; ++q) peaks[q] = choose_int(q, q / 2);
    std::vector<cpp_int> bound(n);
    for (const auto& row : ti.rows) { const int q = static_cast<int>(row.end - row.hold_end);
        for (terminal::Offset i = row.begin; i < row.end; ++i) bound[ti.members[i]] += peaks[q]; }
    cpp_int maximum = d; for (const auto& b : bound) if (b > maximum) maximum = b; return maximum * 2;
}
template<class F> static void by_width(const cpp_int& bound, F&& f) {
    for (unsigned bits : {64u, 128u, 256u, 512u}) {
        if (bound >= (cpp_int(1) << bits) - 1) continue;
        try { if (bits == 64) f(uint64_t{}, bits); else if (bits == 128) f((unsigned __int128){}, bits); else if (bits == 256) f(boost::multiprecision::uint256_t{}, bits); else f(boost::multiprecision::uint512_t{}, bits); return; }
        catch (const std::overflow_error& e) { std::cerr << "[width] " << bits << " bits overflowed (" << e.what() << "); retrying wider\n"; if (bits == 512) throw; }
    }
    throw std::overflow_error("count bound exceeds 512 bits");
}
template<class T> static uint64_t hash_value(uint64_t h, T value) {
    if constexpr (std::is_same_v<T, uint64_t>) { h ^= value; h *= 1099511628211ULL; }
    else if constexpr (std::is_same_v<T, unsigned __int128>) { h ^= static_cast<uint64_t>(value); h *= 1099511628211ULL; h ^= static_cast<uint64_t>(value >> 64); h *= 1099511628211ULL; }
    else { for (int i = 0; i < std::numeric_limits<T>::digits; i += 64) { h ^= static_cast<uint64_t>(value & T(std::numeric_limits<uint64_t>::max())); h *= 1099511628211ULL; value >>= 64; } }
    return h;
}
template<class T> struct Recorder {
    bool keep = false; std::vector<uint64_t> hashes; std::vector<std::vector<T>> rows;
    void add(int s, std::span<const T> row) {
        require(s == static_cast<int>(hashes.size()) + 2, "rows out of order");
        uint64_t h = 1469598103934665603ULL; for (const T& x : row) h = hash_value(h, x); hashes.push_back(h);
        if (keep) rows.emplace_back(row.begin(), row.end());
    }
};
using Clock2 = std::chrono::steady_clock;
static double since(Clock2::time_point a) { return std::chrono::duration<double, std::milli>(Clock2::now() - a).count(); }

static tailpeel::Policy env_policy(int force) {                     // the adaptive thresholds can be overridden for tuning
    tailpeel::Policy p; p.force = force;
    if (const char* x = std::getenv("TAIL_DENSE")) p.dense = std::stod(x);
    if (const char* x = std::getenv("TAIL_RELEVANT")) p.relevant = std::stod(x);
    return p;
}
static void print_sizes(const char* key, const tailpeel::Stats& st) {   // [s, mode, active, residue, residue_degree, valid_incid, rel_incid, rel_vertices, ms]
    std::cout << ",\"" << key << "\":[";
    for (size_t i = 0; i < st.sizes.size(); ++i) { const auto& x = st.sizes[i];
        std::cout << (i ? "," : "") << "[" << x.s << "," << x.mode << "," << x.active << "," << x.residue << "," << x.residue_degree << "," << x.valid_incidences
                  << "," << x.relevant_incidences << "," << x.relevant_vertices << "," << x.ms << "]"; }
    std::cout << "]";
}
template<class T> static void compare_one(const Input& in, const terminal::Index& ti, const tailpeel::Prepared& omega, double prepare_ms, double fused_ms,
                                          bool keep, bool print, unsigned bits, const std::string& name) {
    typename Kernel<T>::Combinations choose(in.d + 1, ti.maximum);
    constexpr int K = 5;   // terminal, tail (adaptive), tail forced to the residue peel, tail forced to the full peel, Tail=false
    std::array<Recorder<T>, K> rec; for (auto& r : rec) r.keep = keep;
    std::array<double, K> t{}; std::array<tailpeel::Stats, K> st;
    auto t0 = Clock2::now();
    terminal::Solver<T>::solve(in.graph, ti, choose, in.ordinary, nullptr, [&](int s, std::span<const T> r) { rec[0].add(s, r); });
    t[0] = since(t0);
    for (int k = 1; k < K; ++k) {
        auto sink = [&](int s, std::span<const T> r) { rec[k].add(s, r); };
        t0 = Clock2::now();
        if (k == 4) st[k] = tailpeel::Solver<T, false>::solve(in.graph, ti, choose, in.ordinary, sink, nullptr, env_policy(0));   // also covers the unprepared path
        else st[k] = tailpeel::Solver<T, true>::solve(in.graph, ti, choose, in.ordinary, sink, &omega, env_policy(k == 1 ? 0 : k == 2 ? 1 : 2));
        t[k] = since(t0);
    }
    for (int k = 1; k < K; ++k) {
        require(rec[k].hashes.size() == rec[0].hashes.size(), "different number of rows");
        for (size_t i = 0; i < rec[0].hashes.size(); ++i) {
            if (keep) require(rec[k].rows[i] == rec[0].rows[i], "row differs from the terminal solver");
            require(rec[k].hashes[i] == rec[0].hashes[i], "row hash differs from the terminal solver");
        }
    }
    for (int k : {1, 2}) for (const auto& x : st[k].sizes) { if (!x.active) continue; ++cov_sizes; cov_settled_new += x.settled_new;
        if (x.residue) { ++cov_residue_sizes; if (x.relevant_vertices > x.residue) ++cov_mixed_sizes; } else ++cov_closed_sizes; }
    if (!print) return;
    const auto& tail = st[1];
    uint64_t act = 0, res = 0, closed = 0, valid = 0, rel = 0;
    for (const auto& x : st[2].sizes) { act += x.active; res += x.residue; closed += x.active && !x.residue; valid += x.valid_incidences * (x.active > 0); rel += x.relevant_incidences; }
    std::cout << std::fixed << std::setprecision(3)
              << "{\"graph\":\"" << name << "\",\"n\":" << in.graph.n << ",\"m\":" << in.graph.m << ",\"d\":" << in.d << ",\"rows\":" << ti.rows.size() << ",\"incidences\":" << ti.members.size()
              << ",\"bits\":" << bits << ",\"sizes\":" << rec[0].hashes.size() << ",\"identical\":true"
              << ",\"prepare_ms\":" << prepare_ms << ",\"fused_prepare_ms\":" << fused_ms << ",\"terminal_ms\":" << t[0] << ",\"tail_ms\":" << t[1] << ",\"residue_ms\":" << t[2] << ",\"fullpeel_ms\":" << t[3] << ",\"notail_ms\":" << t[4]
              << ",\"tail_omega_ms\":" << tail.omega_ms << ",\"tail_total_ms\":" << tail.total_ms << ",\"tail_upper_ms\":" << tail.upper_ms << ",\"tail_mark_ms\":" << tail.mark_ms << ",\"tail_init_ms\":" << tail.init_ms << ",\"tail_peel_ms\":" << tail.peel_ms << ",\"tail_order_ms\":" << tail.order_ms
              << ",\"tail_residue_sizes\":" << tail.residue_sizes << ",\"tail_full_sizes\":" << tail.full_sizes << ",\"tail_aborted_marks\":" << tail.aborted_marks
              << ",\"active_pairs\":" << act << ",\"residue_pairs\":" << res << ",\"closed_sizes\":" << closed << ",\"last_residue_size\":" << tail.last_residue_size
              << ",\"valid_incidences_active_sizes\":" << valid << ",\"relevant_incidences\":" << rel
              << ",\"tail_events\":" << tail.events << ",\"fullpeel_events\":" << st[3].events;
    print_sizes("per_size_tail", st[1]); print_sizes("per_size_residue", st[2]); print_sizes("per_size_full", st[3]);
    std::cout << "}\n";
}
static void run_graph(const Input& in, bool keep, bool print, const std::string& name) {
    const int S = std::max(2, static_cast<int>(in.d) + 1);
    terminal::Index ti(S); terminal::build(in.graph, ti, 0);
    auto t0 = Clock2::now(); ti.prepare(in.graph.n); const double prepare_ms = since(t0);
    const auto plain = ti.reverse; const auto plain_off = ti.reverse_off;
    t0 = Clock2::now(); const tailpeel::Prepared omega = tailpeel::prepare(ti, in.graph.n); const double fused_ms = since(t0);
    require(ti.reverse == plain && ti.reverse_off == plain_off, "fused prepare differs from Index::prepare");
    by_width(bound_rows(ti, in.graph.n, static_cast<int>(in.d)), [&](auto tag, unsigned bits) { using T = decltype(tag); compare_one<T>(in, ti, omega, prepare_ms, fused_ms, keep, print, bits, name); });
}
static void check_small(const Graph& g, uint64_t& graphs) {
    try { Seeds z(g); Input in{g, z.ordinary, z.maximum}; run_graph(in, true, false, "small"); ++graphs; }
    catch (const std::exception& e) { std::cerr << "failure on n=" << g.n << " edges:"; for (Vertex u = 0; u < g.n; ++u) for (Vertex w : g.row(u)) if (u < w) std::cerr << ' ' << u << '-' << w; std::cerr << '\n'; throw; }
}
static void random_suite(uint64_t count) {
    std::mt19937_64 rng(20260923); uint64_t graphs = 0;
    for (Vertex n = 0; n <= 6; ++n) {                                  // every labelled graph up to six vertices
        std::vector<std::pair<Vertex, Vertex>> p; for (Vertex a = 0; a < n; ++a) for (Vertex b = a + 1; b < n; ++b) p.emplace_back(a, b);
        for (uint64_t mask = 0; mask < (uint64_t{1} << p.size()); ++mask) { std::vector<std::pair<Vertex, Vertex>> e; for (size_t i = 0; i < p.size(); ++i) if (mask >> i & 1) e.push_back(p[i]); check_small(Graph::from_edges(n, std::move(e)), graphs); }
    }
    for (Vertex k = 2; k <= 14; ++k) check_small(complete(k), graphs);
    for (Vertex h = 1; h <= 8; ++h) for (Vertex c = 1; c <= 8; ++c) { check_small(split_graph(h, c, false), graphs); check_small(split_graph(h, c, true), graphs); }
    for (uint64_t t = 0; t < count; ++t) {                             // random: G(n,p), and planted cliques with noise
        const Vertex n = 7 + static_cast<Vertex>(rng() % 34);
        std::vector<std::pair<Vertex, Vertex>> e; std::set<std::pair<Vertex, Vertex>> seen;
        auto add = [&](Vertex a, Vertex b) { if (a == b) return; if (a > b) std::swap(a, b); if (seen.insert({a, b}).second) e.emplace_back(a, b); };
        const int kind = static_cast<int>(t % 3);
        const double p = kind == 0 ? 0.05 + 0.9 * std::uniform_real_distribution<double>(0, 1)(rng) : 0.02 + 0.2 * std::uniform_real_distribution<double>(0, 1)(rng);
        for (Vertex a = 0; a < n; ++a) for (Vertex b = a + 1; b < n; ++b) if (std::uniform_real_distribution<double>(0, 1)(rng) < p) add(a, b);
        if (kind >= 1) {                                               // planted cliques of sizes 3..12, overlapping at random
            const int cliques = 1 + static_cast<int>(rng() % 5);
            for (int q = 0; q < cliques; ++q) { const Vertex size = std::min<Vertex>(n, 3 + static_cast<Vertex>(rng() % 10)); std::vector<Vertex> pick(n); std::iota(pick.begin(), pick.end(), 0); std::shuffle(pick.begin(), pick.end(), rng);
                for (Vertex i = 0; i < size; ++i) for (Vertex j = i + 1; j < size; ++j) add(pick[i], pick[j]); }
        }
        check_small(Graph::from_edges(n, std::move(e)), graphs);
    }
    std::cout << "{\"passed\":true,\"graphs\":" << graphs << ",\"active_sizes\":" << cov_sizes << ",\"sizes_with_residue\":" << cov_residue_sizes << ",\"sizes_with_settled_in_the_peel\":" << cov_mixed_sizes << ",\"closed_form_sizes\":" << cov_closed_sizes << ",\"settlements\":" << cov_settled_new << "}\n";
}

template<class T> static void time_one(const Input& in, terminal::Index& ti, const std::string& which, int repeats) {
    typename Kernel<T>::Combinations choose(in.d + 1, ti.maximum); uint64_t h = 0;
    for (int r = 0; r < repeats; ++r) {
        const auto t0 = Clock2::now();
        auto sink = [&](int s, std::span<const T> row) { h ^= static_cast<uint64_t>(s) * row.size(); };
        if (which == "terminal") { ti.prepare(in.graph.n); terminal::Solver<T>::solve(in.graph, ti, choose, in.ordinary, nullptr, sink); }
        else if (which == "probe") { ti.prepare(in.graph.n); const auto ph = terminalprobe::Solver<T>::solve(in.graph, ti, choose, in.ordinary, sink);
            std::cout << "phases alloc " << ph.alloc_ms << " init " << ph.init_ms << " bounds " << ph.bounds_ms << " heap " << ph.heap_ms << " peel " << ph.peel_ms << " total " << ph.total_ms << "\n"; }
        else {
            const tailpeel::Prepared omega = tailpeel::prepare(ti, in.graph.n);
            if (which == "notail") tailpeel::Solver<T, false>::solve(in.graph, ti, choose, in.ordinary, sink, &omega, env_policy(0));
            else if (which == "residue") tailpeel::Solver<T, true>::solve(in.graph, ti, choose, in.ordinary, sink, &omega, env_policy(1));
            else if (which == "fullpeel") tailpeel::Solver<T, true>::solve(in.graph, ti, choose, in.ordinary, sink, &omega, env_policy(2));
            else { require(which == "tail", "unknown solver"); const auto st = tailpeel::Solver<T, true>::solve(in.graph, ti, choose, in.ordinary, sink, &omega, env_policy(0));
                std::cout << "phases omega " << st.omega_ms << " upper " << st.upper_ms << " mark " << st.mark_ms << " init " << st.init_ms << " peel " << st.peel_ms << " order " << st.order_ms << " total " << st.total_ms << "\n"; }
        }
        std::cout << which << " run " << r << ": " << since(t0) << " ms\n" << std::flush;
    }
    std::cout << "checksum " << h << "\n";
}

int main(int argc, char** argv) {
    try {
        if (argc == 5 && std::string(argv[1]) == "--time") {           // --time <terminal|tail|residue|fullpeel|notail> <graph> <repeats>: profiling aid
            const Input in = prepare(argv[3]); const int S = std::max(2, static_cast<int>(in.d) + 1);
            terminal::Index ti(S); terminal::build(in.graph, ti, 0); ti.prepare(in.graph.n);   // bound_rows reads the reverse lists; each timed run prepares again
            by_width(bound_rows(ti, in.graph.n, static_cast<int>(in.d)), [&](auto tag, unsigned) { using T = decltype(tag); time_one<T>(in, ti, argv[2], std::stoi(argv[4])); });
            return 0;
        }
        if (argc == 3 && std::string(argv[1]) == "--random") { random_suite(std::stoull(argv[2])); return 0; }
        if (argc == 3 && std::string(argv[1]) == "--scan") {           // row visits of the tree pass: all rows, rows not yet too small, valid rows
            const Input in = prepare(argv[2]); const int S = std::max(2, static_cast<int>(in.d) + 1);
            terminal::Index ti(S); terminal::build(in.graph, ti, 0); const auto pre = tailpeel::prepare(ti, in.graph.n);
            uint64_t all = 0, alive = 0, valid = 0, incid = ti.members.size();
            for (Vertex v = 0; v < in.graph.n; ++v) { const uint64_t w = pre.omega[v]; if (w < 2) continue;
                for (uint64_t code : ti.touching(v)) { const auto& row = ti.rows[code >> 2];
                    all += w - 1; alive += row.hi - 1; valid += row.hi - std::max<Vertex>(row.lo, 2) + 1; } }
            std::cout << "{\"incidences\":" << incid << ",\"rows\":" << ti.rows.size() << ",\"visits_all\":" << all << ",\"visits_hi_ge_s\":" << alive << ",\"visits_valid\":" << valid << "}\n";
            return 0;
        }
        require(argc == 3 && std::string(argv[1]) == "--graph", "usage: tail_check --random <count> | --graph <path>");
        const Input in = prepare(argv[2]); std::string name = argv[2]; name = name.substr(name.find_last_of('/') + 1);
        run_graph(in, false, true, name);
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
