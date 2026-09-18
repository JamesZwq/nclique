// ChainIndex tool: build from a graph (all-size engine + per-size trees + chains + aligned labels),
// save/load, brute-force selftest with a disk round trip, and the benchmark of RESULTS_FINAL.md.
#define main skyline_stage1_count_main
#include "count.cpp"
#undef main
#include "chain_index.hpp"
#include <chrono>
#include <filesystem>
#include <functional>

using chainindex::ChainIndex; using chainindex::kNone;

struct BuildTimes { double solve_ms = 0, trees_ms = 0, chains_ms = 0, layout_ms = 0; };

// Build the index from an Input; perm[v_input] = internal id.  Returns the index (finished) and perm.
template<class T> static ChainIndex<T> build_chain_index(const Input& in, std::vector<uint32_t>& perm, BuildTimes& bt) {
    using Clock = std::chrono::steady_clock; auto ms = [](Clock::time_point a) { return std::chrono::duration<double, std::milli>(Clock::now() - a).count(); };
    const Graph& g = in.graph; const uint32_t n = g.n; const int S = std::max(2, static_cast<int>(in.d) + 1);
    auto t0 = Clock::now();
    Layout layout(g, S); layout.prepare(n); typename Kernel<T>::Combinations choose(in.d + 1, S);
    terminal::Index ti(S); terminal::build(g, ti, 0); ti.prepare(n);
    auto out = terminal::Solver<T>::solve(g, ti, choose, in.ordinary); const auto& core = out.common.data.core;
    bt.solve_ms = ms(t0); t0 = Clock::now();
    // per-size trees: preorder node ids, per-vertex own node
    ChainIndex<T> ix; ix.n = n; ix.max_size = S; ix.layers.resize(S + 1);
    std::vector<std::vector<uint32_t>> own(S + 1);
    for (int s = 2; s <= S; ++s) {
        auto tr = make_tree(g, ti, core, s); auto& L = ix.layers[s]; std::vector<uint32_t> renum(tr.nodes.size(), kNone);
        std::function<void(int, uint32_t)> dfs = [&](int x, uint32_t p) { const uint32_t id = static_cast<uint32_t>(L.top.size()); renum[x] = id;
            L.top.push_back(static_cast<T>(tr.nodes[x].hi)); L.parent.push_back(p); L.size.push_back(0);
            for (int y : tr.nodes[x].children) dfs(y, id); L.size[id] = static_cast<uint32_t>(L.top.size()) - id; };
        for (size_t i = 0; i < tr.nodes.size(); ++i) if (tr.nodes[i].parent < 0) dfs(static_cast<int>(i), kNone);
        own[s].assign(n, kNone); for (uint32_t v = 0; v < n; ++v) if (tr.leaf[v] >= 0) own[s][v] = renum[tr.leaf[v]];
    }
    bt.trees_ms = ms(t0); t0 = Clock::now();
    // chains: group vertices by their own-node tuple; order chains by the size-2 node (preorder => own-first contiguity at s = 2), inactive last
    std::vector<int> om(n, 0); for (uint32_t v = 0; v < n; ++v) for (int s = 2; s <= S; ++s) if (core[static_cast<size_t>(s) * n + v] > T{0}) om[v] = s;
    std::map<std::vector<uint32_t>, uint32_t> ids; std::vector<uint32_t> tmpchain(n); std::vector<std::vector<uint32_t>> members;
    for (uint32_t v = 0; v < n; ++v) { std::vector<uint32_t> key; for (int s = 2; s <= om[v]; ++s) key.push_back(own[s][v]);
        auto [it, fresh] = ids.emplace(std::move(key), static_cast<uint32_t>(ids.size())); if (fresh) members.emplace_back(); tmpchain[v] = it->second; members[it->second].push_back(v); }
    const uint32_t C = static_cast<uint32_t>(members.size());
    std::vector<uint32_t> order(C); std::iota(order.begin(), order.end(), 0);
    auto key2 = [&](uint32_t c) { const uint32_t v = members[c][0]; return om[v] >= 2 ? own[2][v] : kNone; };
    std::stable_sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) { return key2(a) < key2(b); });
    std::vector<uint32_t> rank_of(C); for (uint32_t r = 0; r < C; ++r) rank_of[order[r]] = r;
    perm.assign(n, kNone); ix.chains = C; ix.start_pos.assign(C + 1, 0); uint32_t next = 0;
    for (uint32_t r = 0; r < C; ++r) { ix.start_pos[r] = next; for (uint32_t v : members[order[r]]) perm[v] = next++; }
    ix.start_pos[C] = next; require(next == n, "relabel incomplete");
    ix.start_bits.assign((n + 63) / 64, 0); for (uint32_t r = 0; r < C; ++r) { const uint32_t p = ix.start_pos[r]; ix.start_bits[p >> 6] |= 1ull << (p & 63); }
    ix.start_cum.assign(ix.start_bits.size() + 1, 0); for (size_t w = 0; w < ix.start_bits.size(); ++w) ix.start_cum[w + 1] = ix.start_cum[w] + static_cast<uint32_t>(__builtin_popcountll(ix.start_bits[w]));
    bt.chains_ms = ms(t0); t0 = Clock::now();
    // per chain: omega, sigma, trajectory, residue (chain r's representative vertex)
    ix.omega.assign(C, 0); ix.sigma.assign(C, 0); ix.traj_off.assign(C + 1, 0); ix.residue_off.assign(C + 1, 0);
    for (uint32_t r = 0; r < C; ++r) { const uint32_t v = members[order[r]][0]; const int o = om[v]; require(o < 256, "omega exceeds a byte");
        ix.omega[r] = static_cast<uint8_t>(o); int sg = o + 1;
        for (int s = 2; s <= o; ++s) if (sg == o + 1 && cpp_int(core[static_cast<size_t>(s) * n + v]) == choose_int(o - 1, s - 1)) sg = s;
        ix.sigma[r] = static_cast<uint8_t>(sg);
        ix.traj_off[r + 1] = ix.traj_off[r] + (o >= 2 ? o - 1 : 0); ix.residue_off[r + 1] = ix.residue_off[r] + (o >= 2 ? sg - 2 : 0); }
    ix.traj_node.assign(ix.traj_off[C], kNone); ix.residue.assign(ix.residue_off[C], T{0});
    for (uint32_t r = 0; r < C; ++r) { const uint32_t v = members[order[r]][0]; const int o = om[v];
        for (int s = 2; s <= o; ++s) ix.traj_node[ix.traj_off[r] + (s - 2)] = own[s][v];
        for (int s = 2; s < ix.sigma[r]; ++s) ix.residue[ix.residue_off[r] + (s - 2)] = core[static_cast<size_t>(s) * n + v]; }
    // per size: slice of chain ids in preorder of their own node (own-first), bucket offsets
    for (int s = 2; s <= S; ++s) { auto& L = ix.layers[s]; const size_t N = L.top.size(); std::vector<uint32_t> cnt(N + 1, 0);
        for (uint32_t r = 0; r < C; ++r) if (s <= ix.omega[r]) ++cnt[ix.traj_node[ix.traj_off[r] + (s - 2)] + 1];
        for (size_t i = 0; i < N; ++i) cnt[i + 1] += cnt[i];
        L.bucket.assign(cnt.begin(), cnt.end() - 1); L.slice.assign(cnt[N], kNone); std::vector<uint32_t> fill(cnt.begin(), cnt.end() - 1);
        for (uint32_t r = 0; r < C; ++r) if (s <= ix.omega[r]) { const uint32_t x = ix.traj_node[ix.traj_off[r] + (s - 2)]; L.slice[fill[x]++] = r; } }
    ix.finish(); bt.layout_ms = ms(t0);
    return ix;
}

// ------------------------------------------------------------ selftest: brute force + disk round trip
template<class T> static void selftest_graph(const Graph& g, const std::string& tmp, uint64_t& queries, uint64_t& members, uint64_t& values, uint64_t& ladders) {
    Seeds z(g); Input in{g, z.ordinary, z.maximum}; std::vector<uint32_t> perm; BuildTimes bt; auto built = build_chain_index<T>(in, perm, bt);
    built.save(tmp); auto ix = ChainIndex<T>::load(tmp);
    require(ix.n == built.n && ix.chains == built.chains && ix.max_size == built.max_size && ix.bytes_total() == built.bytes_total(), "round trip header");
    std::vector<uint32_t> inv(g.n); for (uint32_t v = 0; v < g.n; ++v) inv[perm[v]] = v;
    const int S = ix.max_size; const uint32_t n = g.n;
    // reference core matrix from the frozen control
    Layout layout(g, S); layout.prepare(n); typename Kernel<T>::Combinations choose(in.d + 1, S);
    const auto core = Kernel<T>{}.template fixed_sparse<true>(layout, n, choose, in.ordinary).core;
    std::vector<uint32_t> ranges, ids; std::vector<std::pair<T, uint64_t>> lad;
    for (int s = 2; s <= S + 1; ++s) for (uint32_t v = 0; v < n; ++v) { const T truth = s <= S ? core[static_cast<size_t>(s) * n + v] : T{0}; require(ix.value(perm[v], s) == truth, "value"); ++values; }
    for (int s = 2; s <= S; ++s) { auto cl = bottomup::clique_masks(g, s);
        for (uint32_t v = 0; v < n; ++v) { const T kv = core[static_cast<size_t>(s) * n + v]; if (kv == T{0}) continue;
            for (T k = 1; k <= kv; ++k) {
                std::vector<uint8_t> inside(n); for (uint32_t u = 0; u < n; ++u) inside[u] = core[static_cast<size_t>(s) * n + u] >= k;
                CountDSU brute(n);
                for (auto mask : cl) { bool ok = true; for (uint32_t u = 0; u < n; ++u) if ((mask >> u) & 1) ok &= inside[u];
                    if (ok) { uint32_t first = kNone; for (uint32_t u = 0; u < n; ++u) if ((mask >> u) & 1) { if (first == kNone) first = u; else brute.join(first, u); } } }
                std::vector<uint32_t> truth; for (uint32_t u = 0; u < n; ++u) if (inside[u] && brute.find(u) == brute.find(v)) truth.push_back(u);
                ranges.clear(); require(ix.community_ranges(perm[v], s, k, ranges) != kNone, "community missing"); ChainIndex<T>::expand(ranges, ids);
                for (auto& x : ids) x = inv[x]; std::sort(ids.begin(), ids.end()); require(ids == truth, "community differs from brute force"); ++queries;
                for (uint32_t u = 0; u < n; ++u) { const bool t = inside[u] && brute.find(u) == brute.find(v); require(ix.member(perm[u], perm[v], s, k) == t, "membership"); ++members; }
                if (k == kv) { ix.ladder(perm[v], s, lad); require(!lad.empty() && lad.front().first == kv, "ladder start");
                    for (size_t i = 0; i < lad.size(); ++i) { ranges.clear(); ix.community_ranges(perm[v], s, lad[i].first, ranges); uint64_t cnt = 0; for (size_t j = 0; j < ranges.size(); j += 2) cnt += ranges[j + 1] - ranges[j];
                        require(cnt == lad[i].second, "ladder count"); if (i) require(lad[i].first < lad[i - 1].first, "ladder order"); } ++ladders; }
            } } }
}
static void tool_selftest() {
    const std::string tmp = (std::filesystem::temp_directory_path() / "chainindex_selftest.cx").string();
    uint64_t graphs = 0, queries = 0, members = 0, values = 0, ladders = 0; std::mt19937_64 rng(20260919);
    auto one = [&](const Graph& g) { try { selftest_graph<uint64_t>(g, tmp, queries, members, values, ladders); } catch (const std::exception& e) {
        std::cerr << "selftest failure on n=" << g.n << " edges:"; for (Vertex u = 0; u < g.n; ++u) for (Vertex w : g.row(u)) if (u < w) std::cerr << ' ' << u << '-' << w; std::cerr << '\n'; throw; } ++graphs; };
    for (Vertex n = 0; n <= 6; ++n) { std::vector<std::pair<Vertex, Vertex>> p; for (Vertex a = 0; a < n; ++a) for (Vertex c = a + 1; c < n; ++c) p.emplace_back(a, c);
        for (uint64_t mask = 0; mask < (uint64_t{1} << p.size()); ++mask) { std::vector<std::pair<Vertex, Vertex>> e; for (size_t i = 0; i < p.size(); ++i) if (mask >> i & 1) e.push_back(p[i]); one(Graph::from_edges(n, std::move(e))); } }
    for (int t = 0; t < 200; ++t) { Vertex n = 7 + rng() % 4; std::vector<std::pair<Vertex, Vertex>> e; for (Vertex a = 0; a < n; ++a) for (Vertex c = a + 1; c < n; ++c) if (rng() % 2) e.emplace_back(a, c); one(Graph::from_edges(n, std::move(e))); }
    for (bool x : {false, true}) for (Vertex h : {1, 2, 4}) one(split_graph(h, 4, x));
    one(complete(8)); std::filesystem::remove(tmp);
    std::cout << "{\"passed\":true,\"graphs\":" << graphs << ",\"community_queries\":" << queries << ",\"membership_checks\":" << members << ",\"value_checks\":" << values << ",\"ladder_checks\":" << ladders << "}\n";
}

// ------------------------------------------------------------ benchmark
template<class T> static void bench(const Input& in, unsigned bits, const std::string& outpath) {
    using Clock = std::chrono::steady_clock; auto ms = [](Clock::time_point a) { return std::chrono::duration<double, std::milli>(Clock::now() - a).count(); };
    std::vector<uint32_t> perm; BuildTimes bt; auto t0 = Clock::now(); auto built = build_chain_index<T>(in, perm, bt); const double build_ms = ms(t0);
    t0 = Clock::now(); built.save(outpath); const double save_ms = ms(t0);
    { std::ofstream pf(outpath + ".perm", std::ios::binary); pf.write(reinterpret_cast<const char*>(perm.data()), perm.size() * 4); }
    t0 = Clock::now(); auto ix = ChainIndex<T>::load(outpath); const double load_ms = ms(t0);
    require(ix.bytes_total() == built.bytes_total() && ix.chains == built.chains, "load mismatch");
    const uint64_t file_bytes = std::filesystem::file_size(outpath);
    const uint32_t n = ix.n; std::vector<uint32_t> active; for (uint32_t v = 0; v < n; ++v) if (ix.omega[ix.chain_of(v)] >= 2) active.push_back(v);
    struct Q { uint32_t v, u; int s; T k; }; std::mt19937_64 rng(20260918);
    auto draw = [&](int regime, int count) { std::vector<Q> qs; for (int i = 0; i < count; ++i) { const uint32_t v = active[rng() % active.size()]; const uint32_t c = ix.chain_of(v);
        const int s = 2 + static_cast<int>(rng() % static_cast<uint64_t>(ix.omega[c] - 1)); const T x = ix.value(v, s);
        const T k = regime == 0 ? x : (regime == 1 ? std::max<T>(T{1}, x / 2) : T{1}); qs.push_back({v, static_cast<uint32_t>(rng() % n), s, k}); } return qs; };
    const std::vector<Q> own = draw(0, 20000), half = draw(1, 20000), root = draw(2, 1000);
    const std::vector<Q> mq = [&] { std::vector<Q> m; for (int r = 0; r < 3; ++r) { auto q = draw(r, 6667); m.insert(m.end(), q.begin(), q.end()); } return m; }();
    std::vector<uint32_t> ranges, ids; ranges.reserve(1 << 20); ids.reserve(1 << 22); std::vector<std::pair<T, uint64_t>> lad;
    auto median5 = [](std::array<double, 5> t) { std::sort(t.begin(), t.end()); return t[2]; };
    auto time_ranges = [&](const std::vector<Q>& qs, bool explicit_ids) { std::array<double, 5> ts{}; uint64_t outputs = 0;
        for (int pass = 0; pass < 6; ++pass) { uint64_t cnt = 0; const auto st = Clock::now();
            for (const auto& q : qs) { ranges.clear(); ix.community_ranges(q.v, q.s, q.k, ranges);
                if (explicit_ids) { ChainIndex<T>::expand(ranges, ids); cnt += ids.size(); } else { uint64_t m = 0; for (size_t j = 0; j < ranges.size(); j += 2) m += ranges[j + 1] - ranges[j]; cnt += m; } }
            const double el = std::chrono::duration<double, std::nano>(Clock::now() - st).count() / qs.size(); if (pass == 0) outputs = cnt; else ts[pass - 1] = el; }
        return std::pair<double, double>{median5(ts), double(outputs) / qs.size()}; };
    auto time_member = [&]() { std::array<double, 5> ts{}; uint64_t sum = 0;
        for (int pass = 0; pass < 6; ++pass) { uint64_t z = 0; const auto st = Clock::now(); for (const auto& q : mq) z += ix.member(q.u, q.v, q.s, q.k);
            const double el = std::chrono::duration<double, std::nano>(Clock::now() - st).count() / mq.size(); sum += z; if (pass) ts[pass - 1] = el; }
        return std::pair<double, uint64_t>{median5(ts), sum}; };
    auto time_value = [&]() { std::array<double, 5> ts{}; uint64_t h = 0; std::vector<std::pair<uint32_t, int>> vq;
        for (int i = 0; i < 200000; ++i) { const uint32_t v = static_cast<uint32_t>(rng() % n); vq.emplace_back(v, 2 + static_cast<int>(rng() % static_cast<uint64_t>(ix.omega[ix.chain_of(v)] + 1))); }
        for (int pass = 0; pass < 6; ++pass) { const auto st = Clock::now(); for (const auto& [v, s] : vq) h ^= static_cast<uint64_t>(ix.value(v, s)) + s;
            const double el = std::chrono::duration<double, std::nano>(Clock::now() - st).count() / vq.size(); if (pass) ts[pass - 1] = el; }
        return std::pair<double, uint64_t>{median5(ts), h}; };
    auto time_ladder = [&]() { std::array<double, 5> ts{}; uint64_t steps = 0;
        for (int pass = 0; pass < 6; ++pass) { uint64_t z = 0; const auto st = Clock::now(); for (const auto& q : own) { ix.ladder(q.v, q.s, lad); z += lad.size(); }
            const double el = std::chrono::duration<double, std::nano>(Clock::now() - st).count() / own.size(); if (pass == 0) steps = z; else ts[pass - 1] = el; }
        return std::pair<double, double>{median5(ts), double(steps) / own.size()}; };
    const auto ro = time_ranges(own, false), rh = time_ranges(half, false), rr = time_ranges(root, false), eo = time_ranges(own, true), eh = time_ranges(half, true), er = time_ranges(root, true);
    const auto mm = time_member(); const auto vv = time_value(); const auto ll = time_ladder();
    uint64_t depth_max = 0; for (const auto& L : ix.layers) { std::vector<uint32_t> d(L.top.size(), 0); for (uint32_t x = 0; x < L.top.size(); ++x) { if (L.parent[x] != kNone) d[x] = d[L.parent[x]] + 1; depth_max = std::max<uint64_t>(depth_max, d[x]); } }
    std::cout << std::fixed << std::setprecision(3) << "{\"passed\":true,\"n\":" << n << ",\"m\":" << in.graph.m << ",\"s_max\":" << ix.max_size << ",\"count_bits\":" << bits << ",\"chains\":" << ix.chains
        << ",\"canonical_nodes\":" << ix.node_count() << ",\"max_depth\":" << depth_max
        << ",\"bytes_map\":" << ix.bytes_map() << ",\"bytes_chains\":" << ix.bytes_chains() << ",\"bytes_layers\":" << ix.bytes_layers() << ",\"bytes_total\":" << ix.bytes_total() << ",\"file_bytes\":" << file_bytes << ",\"perm_bytes\":" << 4ull * n
        << ",\"solve_ms\":" << bt.solve_ms << ",\"trees_ms\":" << bt.trees_ms << ",\"chains_ms\":" << bt.chains_ms << ",\"layout_ms\":" << bt.layout_ms << ",\"build_ms\":" << build_ms << ",\"save_ms\":" << save_ms << ",\"load_ms\":" << load_ms
        << ",\"range_own_ns\":" << ro.first << ",\"range_half_ns\":" << rh.first << ",\"range_root_ns\":" << rr.first
        << ",\"explicit_own_ns\":" << eo.first << ",\"explicit_half_ns\":" << eh.first << ",\"explicit_root_ns\":" << er.first
        << ",\"own_output\":" << eo.second << ",\"half_output\":" << eh.second << ",\"root_output\":" << er.second
        << ",\"member_ns\":" << mm.first << ",\"member_checksum\":" << mm.second << ",\"value_ns\":" << vv.first << ",\"value_checksum\":" << vv.second
        << ",\"ladder_ns\":" << ll.first << ",\"ladder_steps\":" << ll.second << "}\n";
}

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--selftest") { tool_selftest(); return 0; }
        require(argc == 4 && std::string(argv[1]) == "--bench", "usage: chain_index_tool --selftest | --bench <graph> <out.cx>");
        Input in = prepare(argv[2]); Layout l(in.graph, std::max(2, static_cast<int>(in.d) + 1)); l.prepare(in.graph.n);
        const unsigned w = width(count_bound(in.graph, l, in.d)); const std::string out = argv[3];
        if (w == 64) bench<uint64_t>(in, w, out); else if (w == 128) bench<unsigned __int128>(in, w, out);
        else if (w == 256) bench<boost::multiprecision::uint256_t>(in, w, out); else bench<boost::multiprecision::uint512_t>(in, w, out);
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
