// Stage 2: build the S-tree baseline and the skyline index (THEORY.md
// Sections 4-5, IMPLEMENTATION2.md) from the stage-1 canonical trees,
// check both against brute force and against each other, measure exact
// bytes and query latency.  The stage-1 source is included verbatim (its
// main renamed) so that make_tree, shadow, twins and the selftest graphs
// are shared byte-for-byte with `count`.
#define main skyline_stage1_count_main
#include "count.cpp"
#undef main

#include <chrono>
#include <functional>
#include <optional>

static constexpr uint32_t none = UINT32_MAX;

template<class T> struct IdxNode { T top{}; uint32_t parent = none, size = 0, bucket = 0, ahi = none; };

template<class T> struct Built {
    const Graph* g = nullptr; int maximum = 2; unsigned bits = 64;
    std::vector<int> cls; std::vector<std::vector<Vertex>> groups;      // shared class maps
    std::vector<int> omega, sigma; std::vector<std::vector<uint8_t>> sky; // per class
    std::vector<std::vector<T>> core;                                     // [s][v], reference only
    std::vector<std::vector<IdxNode<T>>> nodes;                           // [s] DFS-numbered canonical nodes
    // baseline
    std::vector<std::vector<uint32_t>> base_bucket, base_classes;         // [s] own-first DFS class array
    std::vector<uint64_t> leaf_off; std::vector<uint32_t> leaf_ids;      // CSR: class -> X_s(c) for s = 2..omega(c)
    // skyline
    std::vector<std::vector<uint32_t>> sky_classes; std::vector<std::vector<uint8_t>> sky_gamma; // [s] entries
    std::vector<uint64_t> loc_off; std::vector<uint32_t> loc_ids; std::vector<uint8_t> loc_size; // CSR: class -> skyline nodes
    std::vector<std::vector<uint32_t>> rev_off, rev_ids;                   // [s] node of T_s -> nodes of T_{s+1} with ahi there
    // values (Block D, shared)
    std::vector<uint64_t> residue_off; std::vector<T> residue;
    std::optional<typename Kernel<T>::Combinations> choose;
    // byte accounting
    uint64_t bytes_shared = 0, bytes_base_nodes = 0, bytes_base_pairs = 0, bytes_sky_nodes = 0, bytes_sky_entries = 0,
             bytes_sky_location = 0, bytes_sky_cross = 0, bytes_block_d = 0, bytes_block_d_variant_b = 0;
    uint64_t active_pairs = 0, skyline_entries = 0, canonical_nodes = 0;

    T value(uint32_t c, int s) const {
        if (s > omega[c]) return T{0};
        if (s >= sigma[c]) return (*choose)(omega[c] - 1, s - 1);
        return residue[residue_off[c] + static_cast<size_t>(s - 2)];
    }
    uint32_t base_leaf(uint32_t c, int s) const { return s > omega[c] ? none : leaf_ids[leaf_off[c] + static_cast<size_t>(s - 2)]; }
    uint32_t climb(uint32_t x, int s, const T& k) const {
        const auto& ns = nodes[s];
        while (ns[x].parent != none && ns[ns[x].parent].top >= k) x = ns[x].parent;
        return x;
    }
    uint32_t locate(uint32_t c, int s) const {                     // X_s(c) through Location and the ahi chain
        if (s > omega[c]) return none;
        size_t i = loc_off[c], e = loc_off[c + 1];
        while (i < e && loc_size[i] < s) ++i;
        require(i < e, "location has no size >= s for an active class");
        uint32_t x = loc_ids[i];
        for (int t = loc_size[i]; t > s; --t) { x = nodes[t][x].ahi; require(x != none, "missing ahi"); }
        return x;
    }
    // community queries append class ids to `out` (caller-preallocated), return the node
    uint32_t base_community(Vertex v, int s, const T& k, std::vector<uint32_t>& out) const {
        const uint32_t c = cls[v]; const uint32_t l = base_leaf(c, s); if (l == none) return none;
        const uint32_t n = climb(l, s, k); const auto& ns = nodes[s]; const auto& a = base_classes[s];
        const uint32_t b = base_bucket[s][n], e = n + ns[n].size < ns.size() ? base_bucket[s][n + ns[n].size] : static_cast<uint32_t>(a.size());
        out.insert(out.end(), a.begin() + b, a.begin() + e);
        return n;
    }
    uint32_t sky_community(Vertex v, int s, const T& k, std::vector<uint32_t>& out, std::vector<uint32_t>& adm, std::vector<uint32_t>& next) const {
        const uint32_t c = cls[v]; uint32_t n = locate(c, s); if (n == none) return none;
        n = climb(n, s, k);
        const auto& ns = nodes[s]; const auto& a = sky_classes[s];
        const uint32_t hi = n + ns[n].size;
        const uint32_t b = ns[n].bucket, e = hi < ns.size() ? ns[hi].bucket : static_cast<uint32_t>(a.size());
        out.insert(out.end(), a.begin() + b, a.begin() + e);   // size s: every entry qualifies (gamma < s)
        adm.clear(); for (uint32_t x = n; x < hi; ++x) adm.push_back(x);
        for (int t = s; t < maximum && !adm.empty(); ++t) {
            next.clear();
            const auto& ro = rev_off[t]; const auto& ri = rev_ids[t]; const auto& nn = nodes[t + 1];
            const auto& na = sky_classes[t + 1]; const auto& ng = sky_gamma[t + 1];
            for (uint32_t x : adm) for (uint32_t j = ro[x]; j < ro[x + 1]; ++j) {
                const uint32_t m = ri[j]; next.push_back(m);
                const uint32_t mb = nn[m].bucket, me = m + 1 < nn.size() ? nn[m + 1].bucket : static_cast<uint32_t>(na.size());
                for (uint32_t q = mb; q < me && ng[q] < s; ++q) out.push_back(na[q]);
            }
            adm.swap(next);
        }
        return n;
    }
    // aligned labels: internal ids where every chain is one id range; rank over a bitmap of chain starts
    std::vector<Vertex> perm, inv;                 // external -> internal, internal -> external (perm is input translation, inv is test-only)
    std::vector<uint64_t> start_bits; std::vector<uint32_t> start_cum;   // bitmap + cumulative popcount per word
    std::vector<uint32_t> start_pos, chain_at_rank, rank_of_chain;       // range start per rank, rank -> class id, class id -> rank
    uint32_t rank_of(Vertex vint) const { const uint64_t w = vint >> 6; const uint64_t mask = (vint & 63) == 63 ? ~0ull : ((2ull << (vint & 63)) - 1);
        return start_cum[w] + static_cast<uint32_t>(__builtin_popcountll(start_bits[w] & mask)) - 1; }
    uint32_t chain_of_internal(Vertex vint) const { return chain_at_rank[rank_of(vint)]; }
    // community as internal-id ranges (start, end) per chain; returns the node
    uint32_t aligned_community_ranges(Vertex vint, int s, const T& k, std::vector<uint32_t>& ranges) const {
        const uint32_t c = chain_of_internal(vint); const uint32_t l = base_leaf(c, s); if (l == none) return none;
        const uint32_t n = climb(l, s, k); const auto& ns = nodes[s]; const auto& a = base_classes[s];
        const uint32_t b = base_bucket[s][n], e = n + ns[n].size < ns.size() ? base_bucket[s][n + ns[n].size] : static_cast<uint32_t>(a.size());
        for (uint32_t j = b; j < e; ++j) { const uint32_t r = rank_of_chain[a[j]]; ranges.push_back(start_pos[r]); ranges.push_back(start_pos[r + 1]); }
        return n;
    }
    static void ranges_to_ids(const std::vector<uint32_t>& ranges, std::vector<Vertex>& out) {
        size_t total = 0; for (size_t j = 0; j < ranges.size(); j += 2) total += ranges[j + 1] - ranges[j];
        out.resize(total); Vertex* w = out.data();
        for (size_t j = 0; j < ranges.size(); j += 2) { const uint32_t lo = ranges[j], hi = ranges[j + 1]; for (uint32_t x = lo; x < hi; ++x) *w++ = x; } }
    bool aligned_member(Vertex uint_, Vertex vint, int s, const T& k) const {
        const uint32_t lv = base_leaf(chain_of_internal(vint), s), lu = base_leaf(chain_of_internal(uint_), s); if (lv == none || lu == none) return false;
        const uint32_t n = climb(lv, s, k); return lu >= n && lu < n + nodes[s][n].size;
    }
    bool base_member(Vertex u, Vertex v, int s, const T& k) const {
        const uint32_t lv = base_leaf(cls[v], s), lu = base_leaf(cls[u], s); if (lv == none || lu == none) return false;
        const uint32_t n = climb(lv, s, k); return lu >= n && lu < n + nodes[s][n].size;
    }
    bool sky_member(Vertex u, Vertex v, int s, const T& k) const {
        const uint32_t lv = locate(cls[v], s), lu = locate(cls[u], s); if (lv == none || lu == none) return false;
        const uint32_t n = climb(lv, s, k); return lu >= n && lu < n + nodes[s][n].size;
    }
};

enum class ClassMode { twins, chains, vertices, aligned };
template<class T> static Built<T> build_index(const Input& in, unsigned bits, ClassMode mode) {
    Built<T> b; const Graph& g = in.graph; const Vertex n = g.n;
    b.g = &g; b.maximum = std::max(2, static_cast<int>(in.d) + 1); b.bits = bits; b.choose.emplace(in.d + 1, b.maximum);
    const int S = b.maximum;
    Layout layout(g, S); layout.prepare(n);
    terminal::Index ti(S); terminal::build(g, ti, 0); ti.prepare(n);
    auto out = terminal::Solver<T>::solve(g, ti, *b.choose, in.ordinary);
    auto fixed = Kernel<T>{}.template fixed_sparse<true>(layout, n, *b.choose, in.ordinary);
    require(out.common.data.core == fixed.core, "terminal core differs from frozen fixed_sparse control");
    const auto& flat = out.common.data.core;
    b.core.assign(S + 1, std::vector<T>(n));
    for (int s = 2; s <= S; ++s) for (Vertex v = 0; v < n; ++v) b.core[s][v] = flat[static_cast<size_t>(s) * n + v];
    // Phase A: per-size trees, DFS numbering, per-vertex own nodes
    b.nodes.resize(S + 1); b.base_bucket.resize(S + 1); b.base_classes.resize(S + 1); b.sky_classes.resize(S + 1); b.sky_gamma.resize(S + 1);
    std::vector<std::vector<uint32_t>> vertex_leaf(S + 1);        // [s][v] -> node id, none if inactive
    std::vector<std::vector<Vertex>> creator(S + 1);              // [s][node] -> a vertex attached at that node (from make_tree)
    for (int s = 2; s <= S; ++s) {
        auto tr = make_tree(g, ti, flat, s);
        std::vector<uint32_t> renum(tr.nodes.size(), none);
        std::function<void(int, uint32_t)> dfs = [&](int x, uint32_t p) {
            const uint32_t id = static_cast<uint32_t>(b.nodes[s].size()); renum[x] = id;
            b.nodes[s].push_back({static_cast<T>(tr.nodes[x].hi), p, 0, 0, none}); creator[s].push_back(static_cast<Vertex>(tr.nodes[x].creator));
            for (int y : tr.nodes[x].children) dfs(y, id);
            b.nodes[s][id].size = static_cast<uint32_t>(b.nodes[s].size()) - id; };
        for (size_t i = 0; i < tr.nodes.size(); ++i) if (tr.nodes[i].parent < 0) dfs(static_cast<int>(i), none);
        require(b.nodes[s].size() == tr.nodes.size(), "DFS lost nodes");
        vertex_leaf[s].assign(n, none); for (Vertex v = 0; v < n; ++v) if (tr.leaf[v] >= 0) vertex_leaf[s][v] = renum[tr.leaf[v]];
    }
    // Phase B: classes = twins (closed neighbourhoods) or chains (tuple of own nodes over all sizes)
    if (mode == ClassMode::twins) b.cls = twins(g, b.groups);
    else if (mode == ClassMode::vertices) { b.cls.resize(n); b.groups.assign(n, {}); for (Vertex v = 0; v < n; ++v) { b.cls[v] = static_cast<int>(v); b.groups[v] = {v}; } }
    else {
        std::vector<int> om(n, 0); for (Vertex v = 0; v < n; ++v) for (int s = 2; s <= S; ++s) if (b.core[s][v] > T{0}) om[v] = s;
        std::map<std::vector<uint32_t>, uint32_t> ids; b.cls.assign(n, 0); b.groups.clear();
        for (Vertex v = 0; v < n; ++v) { std::vector<uint32_t> key; for (int s = 2; s <= om[v]; ++s) { require(vertex_leaf[s][v] != none, "active vertex without own node"); key.push_back(vertex_leaf[s][v]); }
            auto [it, fresh] = ids.emplace(std::move(key), static_cast<uint32_t>(ids.size())); if (fresh) b.groups.emplace_back(); b.cls[v] = static_cast<int>(it->second); b.groups[it->second].push_back(v); }
    }
    const size_t nc = b.groups.size();
    b.omega.assign(nc, 0); b.sigma.assign(nc, 0); b.sky.assign(nc, std::vector<uint8_t>(S + 1, 0));
    std::map<std::pair<int, cpp_int>, cpp_int> cache;
    auto sig = [&](int s, const T& x) { auto q = std::make_pair(s, cpp_int(x)); auto it = cache.find(q);
        return it == cache.end() ? cache.emplace(q, shadow(s, q.second)).first->second : it->second; };
    for (size_t c = 0; c < nc; ++c) {
        const Vertex v = b.groups[c][0];
        for (int s = 2; s <= S; ++s) if (b.core[s][v] > T{0}) b.omega[c] = s;
        for (Vertex u : b.groups[c]) for (int s = 2; s <= S; ++s) require(b.core[s][u] == b.core[s][v], "F5 twin core equality");
        b.sigma[c] = b.omega[c] + 1;
        for (int s = 2; s <= b.omega[c]; ++s) if (b.sigma[c] == b.omega[c] + 1 && cpp_int(b.core[s][v]) == choose_int(b.omega[c] - 1, s - 1)) b.sigma[c] = s;
        for (int s = 2; s < b.omega[c]; ++s) { const cpp_int a = cpp_int(b.core[s][v]), d = sig(s, b.core[s + 1][v]);
            require(a >= d, "F2 shadow bound"); b.sky[c][s] = a != d; }
        if (b.omega[c] >= 2) b.sky[c][b.omega[c]] = 1;
    }
    // Phase C: per-size buckets over classes
    std::vector<std::vector<uint32_t>> leaf(S + 1);                // [s][class] -> node id, none if inactive (temporary)
    for (int s = 2; s <= S; ++s) {
        leaf[s].assign(nc, none);
        for (size_t c = 0; c < nc; ++c) if (s <= b.omega[c]) { const uint32_t l = vertex_leaf[s][b.groups[c][0]]; require(l != none, "active class without leaf"); leaf[s][c] = l;
            for (Vertex u : b.groups[c]) require(vertex_leaf[s][u] == l, "class members with different own nodes"); }
        const size_t N = b.nodes[s].size();
        // own classes per node (baseline) and skyline entries per node, then packed in node id order (= DFS preorder, own first)
        std::vector<uint32_t> own_cnt(N + 1, 0), sky_cnt(N + 1, 0);
        for (size_t c = 0; c < nc; ++c) if (leaf[s][c] != none) { ++own_cnt[leaf[s][c] + 1]; if (b.sky[c][s]) ++sky_cnt[leaf[s][c] + 1]; }
        for (size_t i = 0; i < N; ++i) { own_cnt[i + 1] += own_cnt[i]; sky_cnt[i + 1] += sky_cnt[i]; }
        // Preorder ids with own-first placement: a node's own classes are placed at own_cnt[id] and its subtree's classes are
        // exactly those of ids in [id, id+size), which is contiguous because ids are preorder.
        b.base_bucket[s].assign(N, 0); b.base_classes[s].assign(own_cnt[N], 0); b.sky_classes[s].assign(sky_cnt[N], 0); b.sky_gamma[s].assign(sky_cnt[N], 0);
        std::vector<uint32_t> own_fill(own_cnt.begin(), own_cnt.end() - 1);
        std::vector<std::vector<std::pair<uint8_t, uint32_t>>> sky_tmp(N);
        for (size_t c = 0; c < nc; ++c) if (leaf[s][c] != none) {
            const uint32_t x = leaf[s][c]; b.base_classes[s][own_fill[x]++] = static_cast<uint32_t>(c);
            if (b.sky[c][s]) { int p = s - 1; while (p >= 2 && !b.sky[c][p]) --p; require(p < 256, "gamma exceeds a byte");
                sky_tmp[x].emplace_back(static_cast<uint8_t>(p >= 2 ? p : 0), static_cast<uint32_t>(c)); } }
        for (uint32_t x = 0; x < N; ++x) {
            b.base_bucket[s][x] = own_cnt[x]; b.nodes[s][x].bucket = sky_cnt[x];
            std::sort(sky_tmp[x].begin(), sky_tmp[x].end());
            for (size_t j = 0; j < sky_tmp[x].size(); ++j) { b.sky_gamma[s][sky_cnt[x] + j] = sky_tmp[x][j].first; b.sky_classes[s][sky_cnt[x] + j] = sky_tmp[x][j].second; } }
        b.canonical_nodes += N; b.active_pairs += own_cnt[N]; b.skyline_entries += sky_cnt[N];
    }
    // baseline leaf CSR
    b.leaf_off.assign(nc + 1, 0);
    for (size_t c = 0; c < nc; ++c) b.leaf_off[c + 1] = b.leaf_off[c] + (b.omega[c] >= 2 ? b.omega[c] - 1 : 0);
    b.leaf_ids.assign(b.leaf_off[nc], none);
    for (size_t c = 0; c < nc; ++c) for (int s = 2; s <= b.omega[c]; ++s) b.leaf_ids[b.leaf_off[c] + (s - 2)] = leaf[s][c];
    // skyline location CSR
    b.loc_off.assign(nc + 1, 0);
    for (size_t c = 0; c < nc; ++c) { uint64_t k = 0; for (int s = 2; s <= b.omega[c]; ++s) k += b.sky[c][s]; b.loc_off[c + 1] = b.loc_off[c] + k; }
    b.loc_ids.assign(b.loc_off[nc], none); b.loc_size.assign(b.loc_off[nc], 0);
    for (size_t c = 0; c < nc; ++c) { uint64_t j = b.loc_off[c]; for (int s = 2; s <= b.omega[c]; ++s) if (b.sky[c][s]) { b.loc_ids[j] = leaf[s][c]; b.loc_size[j] = static_cast<uint8_t>(s); ++j; } }
    // cross-size pointers (THEORY.md Section 6 Step 3) and reverse lists
    b.rev_off.resize(S + 1); b.rev_ids.resize(S + 1);
    for (int s = 2; s < S; ++s) {
        auto& lower = b.nodes[s]; auto& upper = b.nodes[s + 1];
        std::vector<uint32_t> cnt(lower.size() + 1, 0);
        for (uint32_t m = 0; m < upper.size(); ++m) {
            const Vertex u = creator[s + 1][m]; const cpp_int lv = sig(s, upper[m].top);
            uint32_t x = leaf[s][b.cls[u]]; require(x != none, "creator inactive at the smaller size");
            require(cpp_int(lower[x].top) >= lv, "F2 at the container search");
            while (lower[x].parent != none && cpp_int(lower[lower[x].parent].top) >= lv) x = lower[x].parent;
            upper[m].ahi = x; ++cnt[x + 1]; }
        for (size_t i = 0; i < lower.size(); ++i) cnt[i + 1] += cnt[i];
        b.rev_off[s] = cnt; b.rev_ids[s].assign(cnt.back(), none); std::vector<uint32_t> fill(cnt.begin(), cnt.end() - 1);
        for (uint32_t m = 0; m < upper.size(); ++m) b.rev_ids[s][fill[upper[m].ahi]++] = m;
    }
    // Block D
    b.residue_off.assign(nc + 1, 0);
    for (size_t c = 0; c < nc; ++c) { b.residue_off[c] = b.residue.size(); for (int s = 2; s < b.sigma[c]; ++s) b.residue.push_back(b.core[s][b.groups[c][0]]); }
    b.residue_off[nc] = b.residue.size();
    // bytes: W = count width in bytes; node records W+12 (baseline) and W+16 (skyline); pairs/entries 4 B ids; gamma and size bytes 1 B
    const uint64_t W = bits / 8;
    b.bytes_shared = 4ull * n + 4ull * (nc + 1) + 4ull * n;                          // class label per vertex, class->vertex CSR
    b.bytes_base_nodes = b.canonical_nodes * (W + 12); b.bytes_base_pairs = 4ull * b.active_pairs + 4ull * (nc + 1) + 4ull * b.active_pairs;
    b.bytes_sky_nodes = b.canonical_nodes * (W + 16); b.bytes_sky_entries = 5ull * b.skyline_entries;
    b.bytes_sky_location = 5ull * b.skyline_entries + 8ull * (nc + 1);
    for (int s = 2; s < S; ++s) b.bytes_sky_cross += 4ull * (b.rev_off[s].size() + b.rev_ids[s].size());
    uint64_t bentries = 0; for (size_t c = 0; c < nc; ++c) for (int s = 2; s < b.sigma[c]; ++s) if (b.sky[c][s]) ++bentries;
    b.bytes_block_d = 2ull * nc + 8ull * (nc + 1) + W * b.residue.size(); b.bytes_block_d_variant_b = 2ull * nc + 8ull * (nc + 1) + W * bentries;
    if (mode == ClassMode::vertices) b.bytes_shared = 0;   // classes are the vertices themselves: no map needed
    if (mode == ClassMode::aligned) {
        // chain order: DFS order of the size-2 array first (so size-2 communities are single ranges), then the rest
        std::vector<uint32_t> order; std::vector<uint8_t> seen(nc, 0);
        for (uint32_t c : b.base_classes[2]) if (!seen[c]) { seen[c] = 1; order.push_back(c); }
        for (uint32_t c = 0; c < nc; ++c) if (!seen[c]) { seen[c] = 1; order.push_back(c); }
        b.perm.assign(n, absent); b.inv.assign(n, absent); b.start_pos.assign(nc + 1, 0); b.chain_at_rank.assign(nc, 0); b.rank_of_chain.assign(nc, 0);
        Vertex next = 0;
        for (uint32_t r = 0; r < nc; ++r) { const uint32_t c = order[r]; b.chain_at_rank[r] = c; b.rank_of_chain[c] = r; b.start_pos[r] = next;
            for (Vertex v : b.groups[c]) { b.perm[v] = next; b.inv[next] = v; ++next; } }
        b.start_pos[nc] = next; require(next == n, "aligned relabel incomplete");
        b.start_bits.assign((n + 63) / 64, 0); for (uint32_t r = 0; r < nc; ++r) { const Vertex p = b.start_pos[r]; b.start_bits[p >> 6] |= 1ull << (p & 63); }
        b.start_cum.assign(b.start_bits.size() + 1, 0); for (size_t w = 0; w < b.start_bits.size(); ++w) b.start_cum[w + 1] = b.start_cum[w] + static_cast<uint32_t>(__builtin_popcountll(b.start_bits[w]));
        for (Vertex v = 0; v < n; ++v) require(b.chain_of_internal(b.perm[v]) == static_cast<uint32_t>(b.cls[v]), "rank disagrees with the class map");
        b.bytes_shared = 8ull * b.start_bits.size() + 4ull * b.start_cum.size() + 4ull * b.start_pos.size() + 4ull * b.chain_at_rank.size() + 4ull * b.rank_of_chain.size();
    }
    return b;
}

template<class T> static void expand(const Built<T>& b, const std::vector<uint32_t>& classes, std::vector<Vertex>& out) {
    out.clear(); for (uint32_t c : classes) out.insert(out.end(), b.groups[c].begin(), b.groups[c].end());
}

// ---------------------------------------------------------------- selftest against brute force
template<class T> static void selftest_graph(const Graph& g, uint64_t& queries, uint64_t& members, uint64_t& values, ClassMode mode) {
    Seeds z(g); Input in{g, z.ordinary, z.maximum}; auto b = build_index<T>(in, 64, mode);
    const int S = b.maximum; std::vector<uint32_t> cls_out, adm, next; std::vector<Vertex> va, vb;
    for (int s = 2; s <= S; ++s) {
        auto cl = bottomup::clique_masks(g, s);
        for (Vertex v = 0; v < g.n; ++v) {
            for (int t = 2; t <= S + 1; ++t) { const T truth = t <= S ? b.core[t][v] : T{0}; require(b.value(b.cls[v], t) == truth, "selftest value"); ++values; }
            const T kv = b.core[s][v]; if (kv == T{0}) continue;
            for (T k = 1; k <= kv; ++k) {
                std::vector<uint8_t> inside(g.n); for (Vertex u = 0; u < g.n; ++u) inside[u] = b.core[s][u] >= k;
                CountDSU brute(g.n);
                for (auto mask : cl) { bool ok = true; for (Vertex u = 0; u < g.n; ++u) if ((mask >> u) & 1) ok &= inside[u];
                    if (ok) { Vertex first = absent; for (Vertex u = 0; u < g.n; ++u) if ((mask >> u) & 1) { if (first == absent) first = u; else brute.join(first, u); } } }
                std::vector<Vertex> truth; for (Vertex u = 0; u < g.n; ++u) if (inside[u] && brute.find(u) == brute.find(v)) truth.push_back(u);
                cls_out.clear(); require(b.base_community(v, s, k, cls_out) != none, "baseline community missing"); expand(b, cls_out, va); std::sort(va.begin(), va.end());
                cls_out.clear(); require(b.sky_community(v, s, k, cls_out, adm, next) != none, "skyline community missing"); expand(b, cls_out, vb); std::sort(vb.begin(), vb.end());
                require(va == truth, "baseline community differs from brute force"); require(vb == truth, "skyline community differs from brute force"); ++queries;
                if (!b.perm.empty()) { std::vector<uint32_t> rg; require(b.aligned_community_ranges(b.perm[v], s, k, rg) != none, "aligned community missing");
                    std::vector<Vertex> ids; Built<T>::ranges_to_ids(rg, ids); for (auto& x : ids) x = b.inv[x]; std::sort(ids.begin(), ids.end()); require(ids == truth, "aligned community differs from brute force");
                    require(b.value(b.chain_of_internal(b.perm[v]), s) == b.core[s][v], "aligned value"); }
                for (Vertex u = 0; u < g.n; ++u) { const bool t = inside[u] && brute.find(u) == brute.find(v);
                    require(b.base_member(u, v, s, k) == t, "baseline membership"); require(b.sky_member(u, v, s, k) == t, "skyline membership"); ++members;
                    if (!b.perm.empty()) require(b.aligned_member(b.perm[u], b.perm[v], s, k) == t, "aligned membership"); }
            }
        }
    }
}
static void index_selftest() {
    uint64_t graphs = 0, queries = 0, members = 0, values = 0; std::mt19937_64 rng(20260918);
    auto one = [&](const Graph& g) { try { selftest_graph<uint64_t>(g, queries, members, values, ClassMode::twins); selftest_graph<uint64_t>(g, queries, members, values, ClassMode::chains); selftest_graph<uint64_t>(g, queries, members, values, ClassMode::vertices); selftest_graph<uint64_t>(g, queries, members, values, ClassMode::aligned); } catch (const std::exception& e) {
        std::cerr << "selftest failure on n=" << g.n << " edges:"; for (Vertex u = 0; u < g.n; ++u) for (Vertex w : g.row(u)) if (u < w) std::cerr << ' ' << u << '-' << w; std::cerr << '\n'; throw; } ++graphs; };
    for (Vertex n = 0; n <= 6; ++n) { std::vector<std::pair<Vertex, Vertex>> p; for (Vertex a = 0; a < n; ++a) for (Vertex c = a + 1; c < n; ++c) p.emplace_back(a, c);
        for (uint64_t mask = 0; mask < (uint64_t{1} << p.size()); ++mask) { std::vector<std::pair<Vertex, Vertex>> e; for (size_t i = 0; i < p.size(); ++i) if (mask >> i & 1) e.push_back(p[i]); one(Graph::from_edges(n, std::move(e))); } }
    for (int t = 0; t < 200; ++t) { Vertex n = 7 + rng() % 4; std::vector<std::pair<Vertex, Vertex>> e; for (Vertex a = 0; a < n; ++a) for (Vertex c = a + 1; c < n; ++c) if (rng() % 2) e.emplace_back(a, c); one(Graph::from_edges(n, std::move(e))); }
    for (bool x : {false, true}) for (Vertex h : {1, 2, 4}) one(split_graph(h, 4, x));
    one(complete(8));
    std::cout << "{\"passed\":true,\"graphs\":" << graphs << ",\"community_queries\":" << queries << ",\"membership_checks\":" << members << ",\"value_checks\":" << values << "}\n";
}

// ---------------------------------------------------------------- measurement
template<class T> static void run_graph(const Input& in, unsigned bits, ClassMode mode) {
    auto b = build_index<T>(in, bits, mode); const Graph& g = in.graph;
    std::vector<Vertex> active; for (Vertex v = 0; v < g.n; ++v) if (b.omega[b.cls[v]] >= 2) active.push_back(v);
    struct Q { Vertex v, u; int s; T k; };
    std::mt19937_64 rng(20260918);
    auto draw = [&](int regime, int count) { std::vector<Q> qs; for (int i = 0; i < count; ++i) { const Vertex v = active[rng() % active.size()];
        const int s = 2 + static_cast<int>(rng() % static_cast<uint64_t>(b.omega[b.cls[v]] - 1)); const T x = b.value(b.cls[v], s);
        const T k = regime == 0 ? x : (regime == 1 ? std::max<T>(T{1}, x / 2) : T{1}); qs.push_back({v, static_cast<Vertex>(rng() % g.n), s, k}); } return qs; };
    const std::vector<Q> own = draw(0, 20000), half = draw(1, 20000), root = draw(2, 1000), mq = [&] { std::vector<Q> m; for (int r = 0; r < 3; ++r) { auto q = draw(r, 6667); m.insert(m.end(), q.begin(), q.end()); } return m; }();
    std::vector<uint32_t> co, adm, next; std::vector<Vertex> va, vb; co.reserve(1 << 20); va.reserve(1 << 22); vb.reserve(1 << 22);
    std::vector<uint32_t> rg; rg.reserve(1 << 20);
    // correctness pass: both designs agree on every timed query
    for (const auto* qs : {&own, &half, &root}) for (const auto& q : *qs) {
        co.clear(); b.base_community(q.v, q.s, q.k, co); expand(b, co, va); std::sort(va.begin(), va.end());
        co.clear(); b.sky_community(q.v, q.s, q.k, co, adm, next); expand(b, co, vb); std::sort(vb.begin(), vb.end());
        require(va == vb, "timed community cross-check");
        if (!b.perm.empty()) { rg.clear(); b.aligned_community_ranges(b.perm[q.v], q.s, q.k, rg); std::vector<Vertex> ids; Built<T>::ranges_to_ids(rg, ids); for (auto& x : ids) x = b.inv[x]; std::sort(ids.begin(), ids.end()); require(ids == va, "timed aligned cross-check"); } }
    for (const auto& q : mq) { require(b.base_member(q.u, q.v, q.s, q.k) == b.sky_member(q.u, q.v, q.s, q.k), "timed membership cross-check");
        if (!b.perm.empty()) require(b.aligned_member(b.perm[q.u], b.perm[q.v], q.s, q.k) == b.base_member(q.u, q.v, q.s, q.k), "timed aligned membership cross-check"); }
    auto median5 = [](std::array<double, 5> t) { std::sort(t.begin(), t.end()); return t[2]; };
    const bool direct = mode == ClassMode::vertices;           // class ids are vertex ids: the slice copy is the answer
    auto time_aligned = [&](const std::vector<Q>& qs, bool explicit_ids) { std::array<double, 5> ts{}; uint64_t outputs = 0;
        for (int pass = 0; pass < 6; ++pass) { uint64_t n = 0; const auto st = std::chrono::steady_clock::now();
            for (const auto& q : qs) { rg.clear(); b.aligned_community_ranges(b.perm[q.v], q.s, q.k, rg);
                if (explicit_ids) { Built<T>::ranges_to_ids(rg, va); n += va.size(); } else { uint64_t m = 0; for (size_t j = 0; j < rg.size(); j += 2) m += rg[j + 1] - rg[j]; n += m; } }
            const auto en = std::chrono::steady_clock::now(); if (pass == 0) outputs = n; else ts[pass - 1] = std::chrono::duration<double, std::nano>(en - st).count() / qs.size(); }
        return std::pair<double, double>{median5(ts), double(outputs) / qs.size()}; };
    auto time_community = [&](const std::vector<Q>& qs, bool sky) { std::array<double, 5> ts{}; uint64_t outputs = 0;
        for (int pass = 0; pass < 6; ++pass) { uint64_t n = 0; const auto st = std::chrono::steady_clock::now();
            for (const auto& q : qs) { co.clear(); if (sky) b.sky_community(q.v, q.s, q.k, co, adm, next); else b.base_community(q.v, q.s, q.k, co);
                if (direct) { va.assign(co.begin(), co.end()); n += va.size(); } else { expand(b, co, va); n += va.size(); } }
            const auto en = std::chrono::steady_clock::now(); if (pass == 0) outputs = n; else ts[pass - 1] = std::chrono::duration<double, std::nano>(en - st).count() / qs.size(); }
        return std::pair<double, double>{median5(ts), double(outputs) / qs.size()}; };
    auto time_member_aligned = [&]() { std::array<double, 5> ts{}; uint64_t sum = 0; std::vector<std::pair<Vertex, Vertex>> pv; for (const auto& q : mq) pv.emplace_back(b.perm[q.u], b.perm[q.v]);
        for (int pass = 0; pass < 6; ++pass) { uint64_t z = 0; const auto st = std::chrono::steady_clock::now();
            for (size_t i = 0; i < mq.size(); ++i) z += b.aligned_member(pv[i].first, pv[i].second, mq[i].s, mq[i].k);
            const auto en = std::chrono::steady_clock::now(); sum += z; if (pass) ts[pass - 1] = std::chrono::duration<double, std::nano>(en - st).count() / mq.size(); }
        return std::pair<double, uint64_t>{median5(ts), sum}; };
    auto time_member = [&](bool sky) { std::array<double, 5> ts{}; uint64_t sum = 0;
        for (int pass = 0; pass < 6; ++pass) { uint64_t z = 0; const auto st = std::chrono::steady_clock::now();
            for (const auto& q : mq) z += sky ? b.sky_member(q.u, q.v, q.s, q.k) : b.base_member(q.u, q.v, q.s, q.k);
            const auto en = std::chrono::steady_clock::now(); sum += z; if (pass) ts[pass - 1] = std::chrono::duration<double, std::nano>(en - st).count() / mq.size(); }
        return std::pair<double, uint64_t>{median5(ts), sum}; };
    auto time_value = [&]() { std::array<double, 5> ts{}; uint64_t h = 0; std::vector<std::pair<Vertex, int>> vq;
        for (int i = 0; i < 200000; ++i) { const Vertex v = static_cast<Vertex>(rng() % g.n); vq.emplace_back(v, 2 + static_cast<int>(rng() % static_cast<uint64_t>(b.omega[b.cls[v]] + 1))); }
        for (int pass = 0; pass < 6; ++pass) { const auto st = std::chrono::steady_clock::now();
            for (const auto& [v, s] : vq) h ^= static_cast<uint64_t>(b.value(b.cls[v], s)) + s;
            const auto en = std::chrono::steady_clock::now(); if (pass) ts[pass - 1] = std::chrono::duration<double, std::nano>(en - st).count() / vq.size(); }
        return std::pair<double, uint64_t>{median5(ts), h}; };
    const auto bo = time_community(own, false), so = time_community(own, true), bh = time_community(half, false), sh = time_community(half, true),
               br = time_community(root, false), sr = time_community(root, true); const auto bm = time_member(false), sm = time_member(true); const auto val = time_value();
    std::pair<double, double> ao{0, 0}, ah{0, 0}, ar{0, 0}, ro{0, 0}, rh{0, 0}, rr{0, 0}; std::pair<double, uint64_t> am{0, 0};
    if (!b.perm.empty()) { ao = time_aligned(own, true); ah = time_aligned(half, true); ar = time_aligned(root, true); ro = time_aligned(own, false); rh = time_aligned(half, false); rr = time_aligned(root, false); am = time_member_aligned(); }
    const uint64_t base_no_d = b.bytes_shared + b.bytes_base_nodes + b.bytes_base_pairs, sky_no_d = b.bytes_shared + b.bytes_sky_nodes + b.bytes_sky_entries + b.bytes_sky_location + b.bytes_sky_cross;
    const uint64_t aligned_map = (g.n + 7) / 8 + (g.n + 63) / 64 * 4 + 8ull * b.groups.size();   // bitmap + rank directory + one range per class (projection)
    std::cout << std::fixed << std::setprecision(3) << "{\"passed\":true,\"mode\":\"" << (mode == ClassMode::twins ? "twins" : mode == ClassMode::chains ? "chains" : mode == ClassMode::vertices ? "vertices" : "aligned") << "\",\"bytes_shared_aligned\":" << aligned_map << ",\"n\":" << g.n << ",\"m\":" << g.m << ",\"s_max\":" << b.maximum << ",\"count_bits\":" << bits
        << ",\"classes\":" << b.groups.size() << ",\"active_pairs\":" << b.active_pairs << ",\"skyline_entries\":" << b.skyline_entries << ",\"canonical_nodes\":" << b.canonical_nodes
        << ",\"bytes_shared\":" << b.bytes_shared << ",\"bytes_base_nodes\":" << b.bytes_base_nodes << ",\"bytes_base_pairs\":" << b.bytes_base_pairs
        << ",\"bytes_sky_nodes\":" << b.bytes_sky_nodes << ",\"bytes_sky_entries\":" << b.bytes_sky_entries << ",\"bytes_sky_location\":" << b.bytes_sky_location << ",\"bytes_sky_cross\":" << b.bytes_sky_cross
        << ",\"bytes_block_d\":" << b.bytes_block_d << ",\"bytes_block_d_variant_b\":" << b.bytes_block_d_variant_b
        << ",\"base_without_d\":" << base_no_d << ",\"sky_without_d\":" << sky_no_d << ",\"base_with_d\":" << base_no_d + b.bytes_block_d << ",\"sky_with_d\":" << sky_no_d + b.bytes_block_d
        << ",\"own_base_ns\":" << bo.first << ",\"own_sky_ns\":" << so.first << ",\"own_output\":" << bo.second
        << ",\"half_base_ns\":" << bh.first << ",\"half_sky_ns\":" << sh.first << ",\"half_output\":" << bh.second
        << ",\"root_base_ns\":" << br.first << ",\"root_sky_ns\":" << sr.first << ",\"root_output\":" << br.second
        << ",\"aligned_own_ns\":" << ao.first << ",\"aligned_half_ns\":" << ah.first << ",\"aligned_root_ns\":" << ar.first
        << ",\"aligned_range_own_ns\":" << ro.first << ",\"aligned_range_half_ns\":" << rh.first << ",\"aligned_range_root_ns\":" << rr.first << ",\"aligned_member_ns\":" << am.first << ",\"aligned_member_checksum\":" << am.second
        << ",\"member_base_ns\":" << bm.first << ",\"member_sky_ns\":" << sm.first << ",\"member_checksum\":" << bm.second << ",\"member_checksum_sky\":" << sm.second
        << ",\"value_ns\":" << val.first << ",\"value_checksum\":" << val.second << "}\n";
}

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--selftest") { index_selftest(); return 0; }
        require((argc == 3 || argc == 4) && std::string(argv[1]) == "--graph", "usage: index --selftest | --graph path [twins|chains]");
        const std::string m = argc == 4 ? argv[3] : "twins";
        const ClassMode mode = m == "chains" ? ClassMode::chains : m == "vertices" ? ClassMode::vertices : m == "aligned" ? ClassMode::aligned : ClassMode::twins;
        require(m == "twins" || m == "chains" || m == "vertices" || m == "aligned", "unknown class mode");
        Input in = prepare(argv[2]); Layout l(in.graph, std::max(2, static_cast<int>(in.d) + 1)); l.prepare(in.graph.n);
        const unsigned w = width(count_bound(in.graph, l, in.d));
        if (w == 64) run_graph<uint64_t>(in, w, mode); else if (w == 128) run_graph<unsigned __int128>(in, w, mode);
        else if (w == 256) run_graph<boost::multiprecision::uint256_t>(in, w, mode); else run_graph<boost::multiprecision::uint512_t>(in, w, mode);
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
