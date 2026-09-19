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

#if defined(__APPLE__)
#include <mach/mach.h>
static uint64_t rss_now() {   // current resident set size in bytes (macOS)
    mach_task_basic_info info; mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &count) != KERN_SUCCESS) return 0;
    return info.resident_size;
}
#else
#include <unistd.h>
static uint64_t rss_now() {   // current resident set size in bytes (Linux: /proc/self/statm, second field in pages)
    std::ifstream f("/proc/self/statm"); uint64_t size = 0, resident = 0; if (!(f >> size >> resident)) return 0;
    return resident * static_cast<uint64_t>(sysconf(_SC_PAGESIZE));
}
#endif
static double g_ti_ms = 0;   // time to build the shared row index (outside build_chain_index)
struct BuildTimes { double solve_ms = 0, trees_ms = 0, chains_ms = 0, layout_ms = 0;
                    uint64_t rss_start = 0, rss_solve = 0, rss_trees = 0, rss_chains = 0, rss_layout = 0, ti_bytes = 0, core_bytes = 0, own_bytes = 0; };

// The clique-tree row index of the terminal solver, built once per graph and shared by every count width.
static terminal::Index build_terminal_index(const Input& in) {
    const int S = std::max(2, static_cast<int>(in.d) + 1);
    terminal::Index ti(S); terminal::build(in.graph, ti, 0); ti.prepare(in.graph.n); return ti;
}
static uint64_t terminal_index_bytes(const terminal::Index& ti) {
    return ti.rows.size() * sizeof(terminal::Row) + 4ull * (ti.members.size() + ti.reverse.size() + ti.group_row.size()) + 8ull * ti.reverse_off.size() + ti.zero_choice.size();
}
// Upper bound on every count the solver accumulates, from the rows themselves: a row whose free part (pivots and
// choices) has q members holds at most C(q, floor(q/2)) cliques of any one size through any member, so a member's
// support is at most the sum of those peaks over its rows.  One spare bit; the solver still checks every addition.
static cpp_int count_bound_terminal(const terminal::Index& ti, uint32_t n, int d) {
    std::vector<cpp_int> peaks(d + 2); for (int q = 0; q <= d + 1; ++q) peaks[q] = choose_int(q, q / 2);
    std::vector<cpp_int> bound(n);
    for (const auto& row : ti.rows) { const int q = static_cast<int>(row.end - row.hold_end); require(q >= 0 && q < static_cast<int>(peaks.size()), "row free part exceeds the clique bound");
        for (Vertex i = row.begin; i < row.end; ++i) bound[ti.members[i]] += peaks[q]; }
    cpp_int maximum = d; for (const auto& b : bound) if (b > maximum) maximum = b; return maximum * 2;
}
// Run f(T{}, bits) at the first width that holds `bound`; if the solver reports an overflow, retry one width up.
template<class F> static void dispatch_width(const cpp_int& bound, F&& f) {
    for (unsigned bits : {64u, 128u, 256u, 512u}) {
        if (bound >= (cpp_int(1) << bits) - 1) continue;
        try { if (bits == 64) f(uint64_t{}, bits); else if (bits == 128) f((unsigned __int128){}, bits); else if (bits == 256) f(boost::multiprecision::uint256_t{}, bits); else f(boost::multiprecision::uint512_t{}, bits); return; }
        catch (const std::overflow_error& e) { std::cerr << "[width] " << bits << " bits overflowed (" << e.what() << "); retrying wider\n"; if (bits == 512) throw; }
    }
    throw std::overflow_error("count bound exceeds 512 bits");
}

// Build the index from an Input and the shared row index; perm[v_input] = internal id.  Memory O(n + chains x sizes):
// the solver streams one core row at a time (two-row window), each row is turned into that size's canonical tree at
// once, and the chains are refined size by size in a prefix trie (one node per distinct own-node prefix), so neither
// the s_max x n core matrix nor per-size own-node arrays exist.  Chain ranks follow the lexicographic order of the
// preorder-id tuple (X_2, X_3, ...): a trie DFS that emits a node's terminating chain before its children, which are
// created in key order.  Each DFS array visits a node's own chains and child subtrees by ascending smallest rank.
template<class T> static ChainIndex<T> build_chain_index(const Input& in, const terminal::Index& ti, std::vector<uint32_t>& perm, BuildTimes& bt) {
    using Clock = std::chrono::steady_clock; auto ms = [](Clock::time_point a) { return std::chrono::duration<double, std::milli>(Clock::now() - a).count(); };
    const Graph& g = in.graph; const uint32_t n = g.n; const int S = ti.maximum; require(S >= 2, "size bound");
    bt.rss_start = rss_now(); bt.ti_bytes = terminal_index_bytes(ti); bt.core_bytes = 2ull * n * sizeof(T); bt.own_bytes = 0;
    typename Kernel<T>::Combinations choose(in.d + 1, S);
    ChainIndex<T> ix; ix.n = n; ix.max_size = S; ix.layers.resize(S + 1);
    // per-size trees in creation ids: parent, top, children, preorder id (creation child order)
    struct Tree0 { std::vector<int> parent; std::vector<T> hi; std::vector<std::vector<int>> children; std::vector<uint32_t> pre; };
    std::vector<Tree0> trees(S + 1);
    // prefix trie: node = (chain prefix up to some size); level, parent, own tree node (creation id) at that size, kappa there,
    // the chain that terminates here (-1 if none), children in key order
    struct Trie { std::vector<int32_t> parent; std::vector<uint8_t> level; std::vector<uint32_t> own; std::vector<T> kappa; std::vector<int32_t> chain; std::vector<std::vector<int32_t>> children; };
    Trie trie; std::vector<int32_t> cls(n, -1); std::vector<uint8_t> active(n, 0); std::vector<int32_t> chain_of(n, -1);
    std::vector<int32_t> chain_node; uint32_t next_chain = 0;             // per chain in creation order: its trie node
    auto terminate = [&](Vertex v) { const int32_t t = cls[v]; int32_t& c = trie.chain[t];
        if (c < 0) { c = static_cast<int32_t>(next_chain++); chain_node.push_back(t); } chain_of[v] = c; active[v] = 0; };
    double trees_ms = 0; std::vector<std::pair<uint64_t, Vertex>> keys;
    auto on_row = [&](int s, std::span<const T> row) {
        if (s < 2 || s > S) return; const auto t0 = Clock::now();
        auto tr = make_tree_row<T>(g, ti, row, s); auto& t = trees[s]; const size_t N = tr.nodes.size();
        t.parent.resize(N); t.hi.resize(N); t.children.resize(N); t.pre.assign(N, kNone);
        for (size_t i = 0; i < N; ++i) { t.parent[i] = tr.nodes[i].parent; t.hi[i] = static_cast<T>(tr.nodes[i].hi); t.children[i] = tr.nodes[i].children; }
        uint32_t next = 0; std::function<void(int)> dfs0 = [&](int x) { t.pre[x] = next++; for (int y : t.children[x]) dfs0(y); };
        for (size_t i = 0; i < N; ++i) if (t.parent[i] < 0) dfs0(static_cast<int>(i));
        // refine: a vertex active at s extends its prefix by pre_s(own node); one active at s-1 only terminates (omega = s-1)
        keys.clear();
        for (Vertex v = 0; v < n; ++v) {
            if (tr.leaf[v] >= 0) { require(s == 2 || active[v], "support nesting: active at s but not at s-1");
                const uint64_t parent = s == 2 ? 0xFFFFFFFFull : static_cast<uint32_t>(cls[v]); keys.emplace_back((parent << 32) | t.pre[tr.leaf[v]], v); }
            else if (active[v]) terminate(v);
        }
        std::sort(keys.begin(), keys.end());
        for (size_t i = 0; i < keys.size();) { size_t j = i; while (j < keys.size() && keys[j].first == keys[i].first) ++j;
            const Vertex v0 = keys[i].second; const int32_t id = static_cast<int32_t>(trie.parent.size()); const int32_t parent = s == 2 ? -1 : cls[v0];
            trie.parent.push_back(parent); trie.level.push_back(static_cast<uint8_t>(s)); trie.own.push_back(static_cast<uint32_t>(tr.leaf[v0])); trie.kappa.push_back(row[v0]); trie.chain.push_back(-1); trie.children.emplace_back();
            if (parent >= 0) trie.children[parent].push_back(id);
            for (size_t k = i; k < j; ++k) { cls[keys[k].second] = id; active[keys[k].second] = 1; }
            i = j; }
        trees_ms += ms(t0);
    };
    const auto tsolve = Clock::now();
    terminal::Solver<T>::solve(g, ti, choose, in.ordinary, nullptr, on_row);
    for (Vertex v = 0; v < n; ++v) if (active[v]) terminate(v);           // rows past the last delivered one are zero
    bt.solve_ms = ms(tsolve) - trees_ms; bt.trees_ms = trees_ms; bt.rss_solve = rss_now(); bt.rss_trees = bt.rss_solve; auto t0 = Clock::now();
    keys.clear(); keys.shrink_to_fit();
    // chain ranks: trie DFS, a node's terminating chain before its children (children were created in key order); inactive vertices last
    std::vector<uint32_t> rank_of(next_chain, kNone); uint32_t r = 0;
    std::function<void(int32_t)> dfs = [&](int32_t t) { if (trie.chain[t] >= 0) rank_of[static_cast<size_t>(trie.chain[t])] = r++; for (int32_t c : trie.children[t]) dfs(c); };
    for (size_t t = 0; t < trie.parent.size() && trie.level[t] == 2; ++t) dfs(static_cast<int32_t>(t));
    require(r == next_chain, "chain ranking incomplete");
    uint32_t inactive = 0; for (Vertex v = 0; v < n; ++v) if (chain_of[v] < 0) ++inactive;
    const uint32_t C = next_chain + (inactive ? 1 : 0); ix.chains = C;
    std::vector<uint32_t> crank(n); for (Vertex v = 0; v < n; ++v) crank[v] = chain_of[v] < 0 ? next_chain : rank_of[static_cast<size_t>(chain_of[v])];
    ix.start_pos.assign(C + 1, 0); for (Vertex v = 0; v < n; ++v) ++ix.start_pos[crank[v] + 1];
    for (uint32_t c = 0; c < C; ++c) ix.start_pos[c + 1] += ix.start_pos[c];
    { std::vector<uint32_t> fill(ix.start_pos.begin(), ix.start_pos.end() - 1); perm.assign(n, kNone); for (Vertex v = 0; v < n; ++v) perm[v] = fill[crank[v]]++; }
    ix.start_bits.assign((n + 63) / 64, 0); for (uint32_t c = 0; c < C; ++c) { const uint32_t p = ix.start_pos[c]; ix.start_bits[p >> 6] |= 1ull << (p & 63); }
    ix.start_cum.assign(ix.start_bits.size() + 1, 0); for (size_t w = 0; w < ix.start_bits.size(); ++w) ix.start_cum[w + 1] = ix.start_cum[w] + static_cast<uint32_t>(__builtin_popcountll(ix.start_bits[w]));
    std::vector<int32_t> node_of_rank(C, -1); for (uint32_t c = 0; c < next_chain; ++c) node_of_rank[rank_of[c]] = chain_node[c];
    bt.chains_ms = ms(t0); bt.rss_chains = rss_now(); t0 = Clock::now();
    // per chain (by rank): omega = trie level, kappa along the trie path, sigma, residues, trajectory in creation ids
    ix.omega.assign(C, 0); ix.sigma.assign(C, 0); ix.traj_off.assign(C + 1, 0); ix.residue_off.assign(C + 1, 0);
    std::vector<T> path; std::vector<uint32_t> ownpath;
    auto walk = [&](int32_t t) { path.clear(); ownpath.clear(); for (int32_t x = t; x >= 0; x = trie.parent[x]) { path.push_back(trie.kappa[x]); ownpath.push_back(trie.own[x]); }
        std::reverse(path.begin(), path.end()); std::reverse(ownpath.begin(), ownpath.end()); };   // index s-2
    for (uint32_t c = 0; c < C; ++c) { const int32_t t = node_of_rank[c]; const int o = t < 0 ? 0 : trie.level[t]; require(o < 256, "omega exceeds a byte");
        ix.omega[c] = static_cast<uint8_t>(o); int sg = o + 1;
        if (t >= 0) { walk(t); require(static_cast<int>(path.size()) == o - 1, "trie path length");
            for (int s = 2; s <= o; ++s) if (sg == o + 1 && cpp_int(path[s - 2]) == choose_int(o - 1, s - 1)) sg = s; }
        ix.sigma[c] = static_cast<uint8_t>(sg);
        ix.traj_off[c + 1] = ix.traj_off[c] + (o >= 2 ? o - 1 : 0); ix.residue_off[c + 1] = ix.residue_off[c] + (o >= 2 ? sg - 2 : 0); }
    ix.traj_node.assign(ix.traj_off[C], kNone); ix.residue.assign(ix.residue_off[C], T{0});
    for (uint32_t c = 0; c < C; ++c) { const int32_t t = node_of_rank[c]; if (t < 0) continue; walk(t); const int o = ix.omega[c];
        for (int s = 2; s <= o; ++s) ix.traj_node[ix.traj_off[c] + (s - 2)] = ownpath[s - 2];
        for (int s = 2; s < ix.sigma[c]; ++s) ix.residue[ix.residue_off[c] + (s - 2)] = path[s - 2]; }
    // per size: own chains per node (ranks ascending), smallest rank per subtree, DFS by ascending smallest rank => final ids, DFS array, buckets
    for (int s = 2; s <= S; ++s) { const auto& t = trees[s]; const size_t N = t.parent.size(); auto& L = ix.layers[s];
        auto own0 = [&](uint32_t c) { return ix.traj_node[ix.traj_off[c] + (s - 2)]; };
        std::vector<uint32_t> cnt(N + 1, 0); for (uint32_t c = 0; c < C; ++c) if (s <= ix.omega[c]) ++cnt[own0(c) + 1];
        for (size_t i = 0; i < N; ++i) cnt[i + 1] += cnt[i];
        std::vector<uint32_t> ownlist(cnt[N]), fill(cnt.begin(), cnt.end() - 1);
        for (uint32_t c = 0; c < C; ++c) if (s <= ix.omega[c]) ownlist[fill[own0(c)]++] = c;
        std::vector<uint32_t> minrank(N, kNone), bypre(N); for (size_t x = 0; x < N; ++x) bypre[t.pre[x]] = static_cast<uint32_t>(x);
        for (size_t i = N; i-- > 0;) { const uint32_t x = bypre[i]; if (cnt[x + 1] > cnt[x]) minrank[x] = std::min(minrank[x], ownlist[cnt[x]]);
            if (t.parent[x] >= 0) minrank[t.parent[x]] = std::min(minrank[t.parent[x]], minrank[x]); }
        std::vector<uint32_t> renum(N, kNone);
        std::function<void(uint32_t, uint32_t)> dfs1 = [&](uint32_t x, uint32_t p) {
            const uint32_t id = static_cast<uint32_t>(L.top.size()); renum[x] = id;
            L.top.push_back(t.hi[x]); L.parent.push_back(p); L.size.push_back(0); L.bucket.push_back(static_cast<uint32_t>(L.slice.size()));
            std::vector<std::pair<uint32_t, int64_t>> items;   // (key, item): item < 0 encodes own chain rank -(r + 1), item >= 0 a child node
            for (uint32_t j = cnt[x]; j < cnt[x + 1]; ++j) items.emplace_back(ownlist[j], -static_cast<int64_t>(ownlist[j]) - 1);
            for (int y : t.children[x]) items.emplace_back(minrank[y], static_cast<int64_t>(y));
            std::sort(items.begin(), items.end());
            for (const auto& [key, item] : items) { if (item < 0) L.slice.push_back(static_cast<uint32_t>(-item - 1)); else dfs1(static_cast<uint32_t>(item), id); }
            L.size[id] = static_cast<uint32_t>(L.top.size()) - id; };
        std::vector<std::pair<uint32_t, uint32_t>> roots; for (size_t x = 0; x < N; ++x) if (t.parent[x] < 0) roots.emplace_back(minrank[x], static_cast<uint32_t>(x));
        std::sort(roots.begin(), roots.end()); for (const auto& [key, x] : roots) dfs1(x, kNone);
        require(L.slice.size() == cnt[N] && L.top.size() == N, "dfs array incomplete");
        for (uint32_t c = 0; c < C; ++c) if (s <= ix.omega[c]) ix.traj_node[ix.traj_off[c] + (s - 2)] = renum[own0(c)]; }
    ix.finish(); bt.layout_ms = ms(t0); bt.rss_layout = rss_now();
    return ix;
}

// ------------------------------------------------------------ selftest: brute force + disk round trip
template<class T> static void selftest_graph(const Graph& g, const std::string& tmp, uint64_t& queries, uint64_t& members, uint64_t& values, uint64_t& ladders) {
    Seeds z(g); Input in{g, z.ordinary, z.maximum}; std::vector<uint32_t> perm; BuildTimes bt; const terminal::Index ti = build_terminal_index(in); auto built = build_chain_index<T>(in, ti, perm, bt);
    std::vector<uint32_t> inv(g.n); for (uint32_t v = 0; v < g.n; ++v) inv[perm[v]] = v;
    const int S = built.max_size; const uint32_t n = g.n;
    // reference core matrix from the frozen control
    Layout layout(g, S); layout.prepare(n); typename Kernel<T>::Combinations choose(in.d + 1, S);
    const auto core = Kernel<T>{}.template fixed_sparse<true>(layout, n, choose, in.ordinary).core;
    std::vector<uint32_t> ranges, ids; std::vector<std::pair<T, uint64_t>> lad;
    auto check = [&](const ChainIndex<T>& ix) {
        for (int s = 2; s <= S + 1; ++s) for (uint32_t v = 0; v < n; ++v) { const T truth = s <= S ? core[static_cast<size_t>(s) * n + v] : T{0}; require(ix.value(perm[v], s) == truth, "value"); ++values; }
        for (int s = 2; s <= S; ++s) { auto cl = bottomup::clique_masks(g, s);
            for (uint32_t v = 0; v < n; ++v) { const T kv = core[static_cast<size_t>(s) * n + v]; if (kv == T{0}) continue;
                for (T k = 1; k <= kv; ++k) {
                    std::vector<uint8_t> inside(n); for (uint32_t u = 0; u < n; ++u) inside[u] = core[static_cast<size_t>(s) * n + u] >= k;
                    CountDSU brute(n);
                    for (auto mask : cl) { bool ok = true; for (uint32_t u = 0; u < n; ++u) if ((mask >> u) & 1) ok &= inside[u];
                        if (ok) { uint32_t first = kNone; for (uint32_t u = 0; u < n; ++u) if ((mask >> u) & 1) { if (first == kNone) first = u; else brute.join(first, u); } } }
                    std::vector<uint32_t> truth; for (uint32_t u = 0; u < n; ++u) if (inside[u] && brute.find(u) == brute.find(v)) truth.push_back(u);
                    const uint32_t node = ix.community_ranges(perm[v], s, k, ranges); require(node != kNone, "community missing");
                    for (size_t j = 0; j < ranges.size(); j += 2) require(ranges[j] < ranges[j + 1] && (j == 0 || ranges[j] != ranges[j - 1]), "ranges well formed and fully merged");
                    ChainIndex<T>::expand(ranges, ids);
                    for (auto& x : ids) x = inv[x]; std::sort(ids.begin(), ids.end()); require(ids == truth, "community differs from brute force"); ++queries;
                    if (ix.compact) { typename ChainIndex<T>::Runs r; uint32_t nd = kNone; require(ix.community_runs(perm[v], s, k, r, nd) && nd == node && r.count() * 2 == ranges.size(), "pointer form");
                        std::vector<uint32_t> ex(ChainIndex<T>::total(r) + ChainIndex<T>::kSlack); require(ChainIndex<T>::expand(r, ex.data()) == ex.data() + ex.size() - ChainIndex<T>::kSlack, "pointer expand"); ex.resize(ex.size() - ChainIndex<T>::kSlack);
                        for (auto& x : ex) x = inv[x]; std::sort(ex.begin(), ex.end()); require(ex == truth, "pointer form differs from brute force"); }
                    for (uint32_t u = 0; u < n; ++u) { const bool t = inside[u] && brute.find(u) == brute.find(v); require(ix.member(perm[u], perm[v], s, k) == t, "membership"); ++members; }
                    if (k == kv) { ix.ladder(perm[v], s, lad); require(!lad.empty() && lad.front().first == kv, "ladder start");
                        for (size_t i = 0; i < lad.size(); ++i) { ix.community_ranges(perm[v], s, lad[i].first, ranges); const uint64_t cnt = ChainIndex<T>::total(ranges);
                            require(cnt == lad[i].second, "ladder count"); if (i) require(lad[i].first < lad[i - 1].first, "ladder order"); } ++ladders; }
                } } }
        // s = 2: every community is one range (chains ranked by the size-2 preorder)
        for (uint32_t v = 0; v < n; ++v) { const T kv = core[static_cast<size_t>(2) * n + v]; for (T k = 1; k <= kv; ++k) { ix.community_ranges(perm[v], 2, k, ranges); require(ranges.size() == 2, "size-2 community is one range"); } }
    };
    check(built);                                            // build form (chain ids)
    { auto full = built; full.compact_runs(false); require(full.compact && !full.packed_tops, "compaction (full tops)"); check(full); }   // compact form, tops as T
    const uint64_t pairs = built.pairs_total(); built.compact_runs(true); require(built.compact && built.packed_tops && built.runs_total() <= pairs, "compaction");
    check(built);                                            // compact form (runs, packed tops)
    built.save(tmp); auto ix = ChainIndex<T>::load(tmp);
    require(ix.n == built.n && ix.chains == built.chains && ix.max_size == built.max_size && ix.bytes_total() == built.bytes_total() && ix.runs_total() == built.runs_total(), "round trip header");
    check(ix);                                               // loaded
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
template<class T> static void bench(const Input& in, const terminal::Index& ti, unsigned bits, const std::string& outpath) {
    using Clock = std::chrono::steady_clock; auto ms = [](Clock::time_point a) { return std::chrono::duration<double, std::milli>(Clock::now() - a).count(); };
    std::vector<uint32_t> perm; BuildTimes bt; auto t0 = Clock::now(); auto built = build_chain_index<T>(in, ti, perm, bt); const double build_ms = ms(t0);
    const uint64_t slice_bytes_layers = built.bytes_layers(), slice_bytes_total = built.bytes_total(), pairs = built.pairs_total();
    const uint32_t n = built.n; std::vector<uint32_t> active; for (uint32_t v = 0; v < n; ++v) if (built.omega[built.chain_of(v)] >= 2) active.push_back(v);
    const ChainIndex<T>* px = &built;   // the index under measurement: build form, then the compact form with T tops, then the loaded packed form
#define ix (*px)
    struct Q { uint32_t v, u; int s; T k; }; std::mt19937_64 rng(20260918);
    auto draw = [&](int regime, int count) { std::vector<Q> qs; for (int i = 0; i < count; ++i) { const uint32_t v = active[rng() % active.size()]; const uint32_t c = ix.chain_of(v);
        const int s = 2 + static_cast<int>(rng() % static_cast<uint64_t>(ix.omega[c] - 1)); const T x = ix.value(v, s);
        const T k = regime == 0 ? x : (regime == 1 ? std::max<T>(T{1}, x / 2) : T{1}); qs.push_back({v, static_cast<uint32_t>(rng() % n), s, k}); } return qs; };
    const std::vector<Q> own = draw(0, 20000), half = draw(1, 20000), root = draw(2, 1000);
    const std::vector<Q> mq = [&] { std::vector<Q> m; for (int r = 0; r < 3; ++r) { auto q = draw(r, 6667); m.insert(m.end(), q.begin(), q.end()); } return m; }();
    std::vector<uint32_t> ranges; ranges.reserve(1 << 20); std::vector<uint32_t> ids(static_cast<size_t>(n) + ChainIndex<T>::kSlack); std::vector<std::pair<T, uint64_t>> lad;   // ids: caller-owned output buffer with slack
    auto median5 = [](std::array<double, 5> t) { std::sort(t.begin(), t.end()); return t[2]; };
    auto time_climb = [&](const std::vector<Q>& qs) { std::array<double, 5> ts{}; uint64_t z = 0;   // own node lookup + climb only (the part the top encoding touches)
        for (int pass = 0; pass < 6; ++pass) { const auto st = Clock::now(); for (const auto& q : qs) { const uint32_t x = ix.own_node(ix.chain_of(q.v), q.s); z += ix.climb(q.s, x, q.k); }
            const double el = std::chrono::duration<double, std::nano>(Clock::now() - st).count() / qs.size(); if (pass) ts[pass - 1] = el; }
        return std::pair<double, uint64_t>{median5(ts), z}; };
    auto time_ptr = [&](const std::vector<Q>& qs) { std::array<double, 5> ts{}; uint64_t z = 0;   // compact form: climb + pointer, no copy
        for (int pass = 0; pass < 6; ++pass) { const auto st = Clock::now(); for (const auto& q : qs) { typename ChainIndex<T>::Runs r; uint32_t nd = 0; ix.community_runs(q.v, q.s, q.k, r, nd); z += r.nmid + r.lo0 + nd; }
            const double el = std::chrono::duration<double, std::nano>(Clock::now() - st).count() / qs.size(); if (pass) ts[pass - 1] = el; }
        return std::pair<double, uint64_t>{median5(ts), z}; };
    struct RangeStat { double ns, vertices, ranges, chains; };   // per query: latency, output vertices, merged ranges, chains in the segment
    auto time_ranges = [&](const std::vector<Q>& qs, bool explicit_ids) { std::array<double, 5> ts{}; uint64_t outputs = 0, nr = 0, nc = 0;
        for (int pass = 0; pass < 6; ++pass) { uint64_t cnt = 0; const auto st = Clock::now();
            for (const auto& q : qs) { const uint32_t node = ix.community_ranges(q.v, q.s, q.k, ranges);
                if (explicit_ids) cnt += static_cast<uint64_t>(ChainIndex<T>::expand(ranges, ids.data()) - ids.data()); else cnt += ChainIndex<T>::total(ranges);
                if (pass == 0) { nr += ranges.size() / 2; if (ix.compact) nc += ranges.size() / 2; else { uint32_t b, e; ChainIndex<T>::slice_bounds(ix.layers[q.s], node, b, e); nc += e - b; } } }
            const double el = std::chrono::duration<double, std::nano>(Clock::now() - st).count() / qs.size(); if (pass == 0) outputs = cnt; else ts[pass - 1] = el; }
        return RangeStat{median5(ts), double(outputs) / qs.size(), double(nr) / qs.size(), double(nc) / qs.size()}; };
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
    // build form (chain-id DFS arrays): the ablation
    const auto sro = time_ranges(own, false), srh = time_ranges(half, false), srr = time_ranges(root, false), seo = time_ranges(own, true), seh = time_ranges(half, true), ser = time_ranges(root, true);
    const auto sll = time_ladder(); const auto sco = time_climb(own), sch = time_climb(half), scr = time_climb(root);
    // compact form with tops kept as T (isolates the top encoding), measured on a copy
    ChainIndex<T> full = built; full.compact_runs(false); px = &full;
    const auto fro = time_ranges(own, false), frh = time_ranges(half, false), frr = time_ranges(root, false), feo = time_ranges(own, true);
    const auto fco = time_climb(own), fch = time_climb(half), fcr = time_climb(root); const auto fpo = time_ptr(own), fph = time_ptr(half), fpr = time_ptr(root);
    const auto fmm = time_member(); const auto fvv = time_value(); const auto fll = time_ladder(); const uint64_t full_bytes_total = full.bytes_total(), full_bytes_layers = full.bytes_layers();
    full = ChainIndex<T>{}; px = &built;
    // compact form with packed tops (the file format): convert, save, load, measure on the loaded index
    t0 = Clock::now(); built.compact_runs(true); const double compact_ms = ms(t0);
    t0 = Clock::now(); built.save(outpath); const double save_ms = ms(t0);
    { std::ofstream pf(outpath + ".perm", std::ios::binary); pf.write(reinterpret_cast<const char*>(perm.data()), perm.size() * 4); }
    t0 = Clock::now(); auto loaded = ChainIndex<T>::load(outpath); const double load_ms = ms(t0);
    require(loaded.bytes_total() == built.bytes_total() && loaded.chains == built.chains && loaded.runs_total() == built.runs_total(), "load mismatch");
    const uint64_t file_bytes = std::filesystem::file_size(outpath); px = &loaded;
    const auto ro = time_ranges(own, false), rh = time_ranges(half, false), rr = time_ranges(root, false), eo = time_ranges(own, true), eh = time_ranges(half, true), er = time_ranges(root, true);
    require(ro.vertices == sro.vertices && rh.vertices == srh.vertices && rr.vertices == srr.vertices && eo.vertices == seo.vertices, "forms disagree on output size");
    const auto po = time_ptr(own), ph = time_ptr(half), pr = time_ptr(root); const auto co = time_climb(own), ch = time_climb(half), cr = time_climb(root);
    require(co.second == sco.second && ch.second == sch.second && cr.second == scr.second, "climbs differ between forms");
    const auto mm = time_member(); const auto vv = time_value(); const auto ll = time_ladder();
#undef ix
    const ChainIndex<T>& ix = loaded;
    // per-vertex S trees with values, stage-2 `vertices` accounting (index.cpp): nodes (W + 12) each, 8 bytes per (vertex, size) pair
    // (DFS array entry + own-node pointer), 4 (n + 1) offsets, 2 n omega/sigma, 8 (n + 1) residue offsets, W per residue cell
    uint64_t pairs_v = 0, residue_v = 0; for (uint32_t c = 0; c < ix.chains; ++c) { const uint64_t sz = ix.start_pos[c + 1] - ix.start_pos[c]; if (ix.omega[c] >= 2) { pairs_v += sz * (ix.omega[c] - 1); residue_v += sz * (ix.sigma[c] - 2); } }
    const uint64_t baseline_vertex_bytes = ix.node_count() * (chainindex::Traits<T>::W + 12) + 8ull * pairs_v + 4ull * (n + 1) + 2ull * n + 8ull * (n + 1) + chainindex::Traits<T>::W * residue_v;
    uint64_t depth_max = 0; for (const auto& L : ix.layers) { std::vector<uint32_t> d(L.size.size(), 0); for (uint32_t x = 0; x < L.size.size(); ++x) { if (L.parent[x] != kNone) d[x] = d[L.parent[x]] + 1; depth_max = std::max<uint64_t>(depth_max, d[x]); } }
    std::cout << std::fixed << std::setprecision(3) << "{\"passed\":true,\"n\":" << n << ",\"m\":" << in.graph.m << ",\"s_max\":" << ix.max_size << ",\"count_bits\":" << bits << ",\"chains\":" << ix.chains
        << ",\"canonical_nodes\":" << ix.node_count() << ",\"max_depth\":" << depth_max
        << ",\"pairs_total\":" << pairs << ",\"runs_total\":" << ix.runs_total() << ",\"vertex_pairs\":" << pairs_v << ",\"vertex_residue_cells\":" << residue_v << ",\"baseline_vertex_bytes\":" << baseline_vertex_bytes
        << ",\"bytes_map\":" << ix.bytes_map() << ",\"bytes_chains\":" << ix.bytes_chains() << ",\"bytes_layers\":" << ix.bytes_layers() << ",\"bytes_total\":" << ix.bytes_total() << ",\"file_bytes\":" << file_bytes << ",\"perm_bytes\":" << 4ull * n
        << ",\"slice_bytes_layers\":" << slice_bytes_layers << ",\"slice_bytes_total\":" << slice_bytes_total
        << ",\"solve_ms\":" << bt.solve_ms << ",\"trees_ms\":" << bt.trees_ms << ",\"chains_ms\":" << bt.chains_ms << ",\"layout_ms\":" << bt.layout_ms << ",\"build_ms\":" << build_ms << ",\"compact_ms\":" << compact_ms << ",\"save_ms\":" << save_ms << ",\"load_ms\":" << load_ms
        << ",\"ti_ms\":" << g_ti_ms << ",\"ti_bytes\":" << bt.ti_bytes << ",\"rss_start\":" << bt.rss_start << ",\"rss_after_solve\":" << bt.rss_solve << ",\"rss_after_chains\":" << bt.rss_chains << ",\"rss_after_layout\":" << bt.rss_layout
        << ",\"slice_range_own_ns\":" << sro.ns << ",\"slice_range_half_ns\":" << srh.ns << ",\"slice_range_root_ns\":" << srr.ns
        << ",\"slice_explicit_own_ns\":" << seo.ns << ",\"slice_explicit_half_ns\":" << seh.ns << ",\"slice_explicit_root_ns\":" << ser.ns << ",\"slice_ladder_ns\":" << sll.first
        << ",\"slice_own_ranges\":" << sro.ranges << ",\"slice_half_ranges\":" << srh.ranges << ",\"slice_root_ranges\":" << srr.ranges
        << ",\"full_bytes_total\":" << full_bytes_total << ",\"full_bytes_layers\":" << full_bytes_layers
        << ",\"full_range_own_ns\":" << fro.ns << ",\"full_range_half_ns\":" << frh.ns << ",\"full_range_root_ns\":" << frr.ns << ",\"full_explicit_own_ns\":" << feo.ns
        << ",\"full_climb_own_ns\":" << fco.first << ",\"full_climb_half_ns\":" << fch.first << ",\"full_climb_root_ns\":" << fcr.first
        << ",\"full_ptr_own_ns\":" << fpo.first << ",\"full_ptr_half_ns\":" << fph.first << ",\"full_ptr_root_ns\":" << fpr.first
        << ",\"full_member_ns\":" << fmm.first << ",\"full_value_ns\":" << fvv.first << ",\"full_ladder_ns\":" << fll.first
        << ",\"climb_own_ns\":" << co.first << ",\"climb_half_ns\":" << ch.first << ",\"climb_root_ns\":" << cr.first
        << ",\"slice_climb_own_ns\":" << sco.first << ",\"slice_climb_half_ns\":" << sch.first << ",\"slice_climb_root_ns\":" << scr.first
        << ",\"ptr_own_ns\":" << po.first << ",\"ptr_half_ns\":" << ph.first << ",\"ptr_root_ns\":" << pr.first << ",\"ptr_checksum\":" << (po.second ^ ph.second ^ pr.second)
        << ",\"range_own_ns\":" << ro.ns << ",\"range_half_ns\":" << rh.ns << ",\"range_root_ns\":" << rr.ns
        << ",\"explicit_own_ns\":" << eo.ns << ",\"explicit_half_ns\":" << eh.ns << ",\"explicit_root_ns\":" << er.ns
        << ",\"own_output\":" << eo.vertices << ",\"half_output\":" << eh.vertices << ",\"root_output\":" << er.vertices
        << ",\"own_ranges\":" << ro.ranges << ",\"half_ranges\":" << rh.ranges << ",\"root_ranges\":" << rr.ranges
        << ",\"own_chains\":" << ro.chains << ",\"half_chains\":" << rh.chains << ",\"root_chains\":" << rr.chains
        << ",\"member_ns\":" << mm.first << ",\"member_checksum\":" << mm.second << ",\"value_ns\":" << vv.first << ",\"value_checksum\":" << vv.second
        << ",\"ladder_ns\":" << ll.first << ",\"ladder_steps\":" << ll.second << "}\n";
}

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--selftest") { tool_selftest(); return 0; }
        require((argc == 3 && std::string(argv[1]) == "--build-only") || (argc == 4 && std::string(argv[1]) == "--bench"), "usage: chain_index_tool --selftest | --build-only <graph> | --bench <graph> <out.cx>");
        Input in = prepare(argv[2]); const uint64_t rss_loaded = rss_now(); const auto tti = std::chrono::steady_clock::now();
        const terminal::Index ti = build_terminal_index(in); const cpp_int bound = count_bound_terminal(ti, in.graph.n, in.d); const uint64_t rss_ti = rss_now();
        g_ti_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tti).count();
        if (std::string(argv[1]) == "--build-only") {
            dispatch_width(bound, [&](auto tag, unsigned bits) { using T = decltype(tag); std::vector<uint32_t> perm; BuildTimes bt; auto ix = build_chain_index<T>(in, ti, perm, bt);
                const uint64_t build_form_bytes = ix.bytes_total(); const auto tc = std::chrono::steady_clock::now(); ix.compact_runs(true); const double compact_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tc).count();
                std::cout << "{\"n\":" << in.graph.n << ",\"count_bits\":" << bits << ",\"chains\":" << ix.chains << ",\"canonical_nodes\":" << ix.node_count() << ",\"index_bytes\":" << ix.bytes_total() << ",\"build_form_bytes\":" << build_form_bytes << ",\"compact_ms\":" << compact_ms
                    << ",\"bytes_map\":" << ix.bytes_map() << ",\"bytes_chains\":" << ix.bytes_chains() << ",\"bytes_layers\":" << ix.bytes_layers()
                    << ",\"ti_ms\":" << g_ti_ms << ",\"ti_bytes\":" << bt.ti_bytes << ",\"rss_loaded\":" << rss_loaded << ",\"rss_with_ti\":" << rss_ti
                    << ",\"rss_start\":" << bt.rss_start << ",\"rss_after_solve\":" << bt.rss_solve << ",\"rss_after_chains\":" << bt.rss_chains << ",\"rss_after_layout\":" << bt.rss_layout
                    << ",\"solve_ms\":" << bt.solve_ms << ",\"trees_ms\":" << bt.trees_ms << ",\"chains_ms\":" << bt.chains_ms << ",\"layout_ms\":" << bt.layout_ms << "}\n"; });
            return 0;
        }
        const std::string out = argv[3];
        dispatch_width(bound, [&](auto tag, unsigned bits) { using T = decltype(tag); bench<T>(in, ti, bits, out); });
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
