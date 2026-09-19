// ChainIndex: the r = 1 all-size community/value index (THEORY.md, CHAINS.md, RESULTS_FINAL.md).
//
// Vertices are partitioned into chains (equal own canonical node at every size, Lemma C1) and relabelled so
// that every chain is one id range; the vertex-to-chain map is a bitmap with a rank directory.  Per size s the
// canonical merge tree is stored over chains with node tops, parents, subtree sizes, and its DFS array as
// maximal runs of consecutive labels plus a per-node entry point (run, label): a community is located in O(1)
// after the climb and reported as head range + whole runs + tail range (Lemma C6); every (2,k)-community is
// one range (Lemma C4).  Values: per chain omega, sigma and the residues for s < sigma; s >= sigma is
// C(omega-1, s-1).  Node tops and residues are stored at per-size byte widths (1, 2, 4, 8, 16 or 32 bytes;
// aligned constant-width loads, one width dispatch per query).  Flat arrays only; file format CHAINX03.
// Header-only, templated on the count type T (64/128-bit builtins or boost fixed-width integers).
#pragma once
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

namespace chainindex {

static constexpr uint32_t kNone = UINT32_MAX;

template<class T> struct Traits {                       // little-endian bytes <-> T, at the full width W or at a given width w <= W
    static constexpr size_t W = sizeof(T);
    static void put(const T& x, unsigned char* p, size_t w = W) { std::memcpy(p, &x, w); }
    static T get(const unsigned char* p, size_t w = W) { T x = 0; std::memcpy(&x, p, w); return x; }
    template<size_t Wb> static T load(const unsigned char* p) { T x = 0; std::memcpy(&x, p, Wb); return x; }   // constant-size load, Wb <= W
    static unsigned width(const T& x) { unsigned w = 1; T y = x >> 8; while (y != 0) { ++w; y >>= 8; } return w; }   // bytes needed
};
template<unsigned Bits> struct Traits<boost::multiprecision::number<boost::multiprecision::cpp_int_backend<Bits, Bits, boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>> {
    using T = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<Bits, Bits, boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>;
    static constexpr size_t W = Bits / 8;
    static void put(const T& x, unsigned char* p, size_t w = W) { std::memset(p, 0, w); std::vector<unsigned char> tmp; boost::multiprecision::export_bits(x, std::back_inserter(tmp), 8, false); std::memcpy(p, tmp.data(), std::min(tmp.size(), w)); }
    static T get(const unsigned char* p, size_t w = W) { if (w <= 8) { uint64_t v = 0; std::memcpy(&v, p, w); return T(v); } T x; boost::multiprecision::import_bits(x, p, p + w, 8, false); return x; }
    template<size_t Wb> static T load(const unsigned char* p) { if constexpr (Wb <= 8) { uint64_t v = 0; std::memcpy(&v, p, Wb); return T(v); } else { T x; boost::multiprecision::import_bits(x, p, p + Wb, 8, false); return x; } }
    static unsigned width(const T& x) { return x == 0 ? 1u : static_cast<unsigned>(boost::multiprecision::msb(x) / 8 + 1); }
};
// stored widths are rounded up to 1, 2, 4, 8, 16, 32 or 64 bytes: a field of w bytes at index y sits at offset y * w, aligned
// to w, so a constant-size load never crosses a cache line; the width is dispatched once per query (with_width)
static inline unsigned round_width(unsigned w) { unsigned r = 1; while (r < w) r <<= 1; return r; }
template<class F> static inline auto with_width(unsigned w, F&& f) {
    switch (w) { case 1: return f(std::integral_constant<size_t, 1>{}); case 2: return f(std::integral_constant<size_t, 2>{}); case 4: return f(std::integral_constant<size_t, 4>{});
        case 8: return f(std::integral_constant<size_t, 8>{}); case 16: return f(std::integral_constant<size_t, 16>{}); case 32: return f(std::integral_constant<size_t, 32>{});
        default: return f(std::integral_constant<size_t, 64>{}); }
}

template<class T> struct ChainIndex {
    // header
    uint32_t n = 0, chains = 0; int max_size = 2;
    // block 1: vertex -> chain (bitmap of chain starts + rank directory), chain -> vertex range
    std::vector<uint64_t> start_bits; std::vector<uint32_t> start_cum; std::vector<uint32_t> start_pos;   // start_pos[chains] = n
    // block 2: per chain
    std::vector<uint8_t> omega, sigma; std::vector<uint32_t> traj_off, traj_node;   // traj_node[traj_off[c] + (s-2)] = own node of chain c at size s
    std::vector<uint32_t> residue_off; std::vector<T> residue;                     // build form: residue values for 2 <= s < sigma(c), element offsets
    // compact form: per size s a byte width res_w[s] (the widest residue stored at s), res_prefix[s] = sum of res_w[2..s-1];
    // chain c's residues are the bytes residue_bytes[residue_boff[c] + res_prefix[s], +res_w[s]) for 2 <= s < sigma(c)
    std::vector<uint8_t> res_w; std::vector<uint32_t> res_prefix, residue_boff; std::vector<uint8_t> residue_bytes;
    // block 3: per size s (index s, 0 and 1 unused).  Nodes in DFS preorder; node x's subtree is [x, x + size[x]).
    // build form:   slice = DFS array of chain ids; node x holds slice[bucket[x], bucket[x + size[x]]).
    // compact form: runs = the DFS array as maximal runs of consecutive vertex ids (lo, hi pairs);
    //               entry[2x], entry[2x+1] = (run, vertex) where node x's segment begins; entry[2N] = (#runs, 0) is the sentinel.
    //               compact_runs() converts; files hold the compact form.
    // top (build form): node top values as T; compact form: top_bytes at top_w bytes per node (the widest top of the layer)
    struct Layer { std::vector<T> top; std::vector<uint8_t> top_bytes; uint8_t top_w = 1; std::vector<uint32_t> parent, size, jump, bucket, slice, entry, runs; };
    std::vector<Layer> layers;
    bool compact = false;            // run form (entry points + runs) instead of chain-id arrays
    bool packed_tops = true;         // compact form: node tops at per-layer byte widths (else kept as T)
    // binomials for the certified tail, C(a, b) for a <= max_size
    std::vector<std::vector<T>> binom;

    // ---- derived accessors
    uint32_t chain_of(uint32_t v) const { const uint32_t w = v >> 6; const uint64_t mask = (v & 63) == 63 ? ~0ull : ((2ull << (v & 63)) - 1);
        return start_cum[w] + static_cast<uint32_t>(__builtin_popcountll(start_bits[w] & mask)) - 1; }
    uint32_t own_node(uint32_t c, int s) const { return s > omega[c] || s < 2 ? kNone : traj_node[traj_off[c] + static_cast<size_t>(s - 2)]; }
    T top_at(const Layer& L, uint32_t x) const { return compact && packed_tops ? with_width(L.top_w, [&](auto wb) { return Traits<T>::template load<decltype(wb)::value>(L.top_bytes.data() + static_cast<size_t>(x) * decltype(wb)::value); }) : L.top[x]; }
    template<size_t Wb> uint32_t climb_bytes(const Layer& L, uint32_t x, const T& k) const {   // compact form: tops at Wb bytes, aligned constant-size loads
        const unsigned char* tb = L.top_bytes.data();
        auto top = [&](uint32_t y) -> T { return Traits<T>::template load<Wb>(tb + static_cast<size_t>(y) * Wb); };
        while (L.parent[x] != kNone && top(L.parent[x]) >= k) { const uint32_t j = L.jump[x]; if (j != kNone && top(j) >= k) x = j; else x = L.parent[x]; }
        return x;
    }
    uint32_t climb_values(const Layer& L, uint32_t x, const T& k) const {                   // build form: tops as T
        while (L.parent[x] != kNone && L.top[L.parent[x]] >= k) { const uint32_t j = L.jump[x]; if (j != kNone && L.top[j] >= k) x = j; else x = L.parent[x]; }
        return x;
    }
    uint32_t climb(int s, uint32_t x, const T& k) const {   // highest ancestor with top >= k; tops decrease upward
        const Layer& L = layers[s]; if (!compact || !packed_tops) return climb_values(L, x, k);
        return with_width(L.top_w, [&](auto wb) { return climb_bytes<decltype(wb)::value>(L, x, k); });
    }
    // ---- queries (internal labels)
    T value(uint32_t v, int s) const {
        const uint32_t c = chain_of(v);
        if (s < 2 || s > omega[c]) return T{0};
        if (s >= sigma[c]) return binom[omega[c] - 1][s - 1];
        if (compact) { const unsigned char* q = residue_bytes.data() + residue_boff[c] + res_prefix[s]; return with_width(res_w[s], [&](auto wb) { return Traits<T>::template load<decltype(wb)::value>(q); }); }
        return residue[residue_off[c] + static_cast<size_t>(s - 2)];
    }
    // build form: segment of the DFS array at size s holding the chains of `node` and its descendants: slice[b, e)
    static void slice_bounds(const Layer& L, uint32_t node, uint32_t& b, uint32_t& e) {
        b = L.bucket[node]; const uint32_t hi = node + L.size[node]; e = hi < L.bucket.size() ? L.bucket[hi] : static_cast<uint32_t>(L.slice.size()); }
    // compact form: a community as head range [lo0, hi0), nmid whole runs at mid (lo, hi pairs), tail range [lo1, hi1); empty parts have lo == hi
    struct Runs { uint32_t lo0, hi0; const uint32_t* mid; uint32_t nmid; uint32_t lo1, hi1; uint32_t count() const { return (lo0 < hi0) + nmid + (lo1 < hi1); } };
    static void node_runs(const Layer& L, uint32_t x, Runs& out) {
        const size_t e = static_cast<size_t>(x) + L.size[x];                    // exit = entry of the next preorder node (sentinel at N)
        const uint32_t rx = L.entry[2 * static_cast<size_t>(x)], vx = L.entry[2 * static_cast<size_t>(x) + 1], re = L.entry[2 * e], ve = L.entry[2 * e + 1];
        if (rx == re) { out = {vx, ve, nullptr, 0, 0, 0}; return; }
        out.lo0 = vx; out.hi0 = L.runs[2 * static_cast<size_t>(rx) + 1]; out.mid = L.runs.data() + 2 * (static_cast<size_t>(rx) + 1); out.nmid = re - rx - 1;
        if (re < L.runs.size() / 2) { out.lo1 = L.runs[2 * static_cast<size_t>(re)]; out.hi1 = ve; } else { out.lo1 = out.hi1 = 0; }
    }
    // compact form only: locate the community of v at (s, k) in O(1); false if v is inactive at s
    bool community_runs(uint32_t v, int s, const T& k, Runs& out, uint32_t& node) const {
        const uint32_t x = own_node(chain_of(v), s); if (x == kNone) { node = kNone; return false; }
        node = climb(s, x, k); node_runs(layers[s], node, out); return true;
    }
    // community of v at (s, k) as vertex-id ranges (pairs lo, hi; fully merged); `ranges` is replaced; returns the node or kNone
    uint32_t community_ranges(uint32_t v, int s, const T& k, std::vector<uint32_t>& ranges) const {
        ranges.clear(); const uint32_t c = chain_of(v); const uint32_t x = own_node(c, s); if (x == kNone) return kNone;
        const Layer& L = layers[s]; const uint32_t node = climb(s, x, k);
        if (compact) { Runs r; node_runs(L, node, r);
            if (r.lo0 < r.hi0) { ranges.push_back(r.lo0); ranges.push_back(r.hi0); }
            ranges.insert(ranges.end(), r.mid, r.mid + 2 * static_cast<size_t>(r.nmid));
            if (r.lo1 < r.hi1) { ranges.push_back(r.lo1); ranges.push_back(r.hi1); }
            return node; }
        uint32_t b, e; slice_bounds(L, node, b, e);
        for (uint32_t j = b; j < e; ++j) { const uint32_t cc = L.slice[j]; const uint32_t lo = start_pos[cc], hi = start_pos[cc + 1];
            if (!ranges.empty() && ranges.back() == lo) ranges.back() = hi; else { ranges.push_back(lo); ranges.push_back(hi); } }
        return node;
    }
    static uint64_t total(const std::vector<uint32_t>& ranges) { uint64_t t = 0; for (size_t j = 0; j < ranges.size(); j += 2) t += ranges[j + 1] - ranges[j]; return t; }
    static uint64_t total(const Runs& r) { uint64_t t = (r.hi0 - r.lo0) + (r.hi1 - r.lo1); for (uint32_t j = 0; j < r.nmid; ++j) t += r.mid[2 * static_cast<size_t>(j) + 1] - r.mid[2 * static_cast<size_t>(j)]; return t; }
    // explicit vertex ids of `ranges` written to out[0, total); returns one past the last id written.  Branchless per range:
    // every block of eight ids is stored unconditionally and the pointer advances by the true length, so the caller's
    // buffer must have kSlack spare slots after `total`.
    static constexpr size_t kSlack = 8;
    typedef uint32_t V4 __attribute__((vector_size(16)));
    static inline void store8(uint32_t x, uint32_t* w) { const V4 i0 = {0, 1, 2, 3}, i1 = {4, 5, 6, 7}; const V4 xs = {x, x, x, x}; const V4 a = i0 + xs, b = i1 + xs; std::memcpy(w, &a, 16); std::memcpy(w + 4, &b, 16); }
    static uint32_t* fill(uint32_t x, const uint32_t hi, uint32_t* w) {
        const uint32_t len = hi - x; store8(x, w);
        if (len <= 8) return w + len;
        uint32_t* const end = w + len; x += 8; w += 8;
        for (; w + 8 <= end; x += 8, w += 8) store8(x, w);
        store8(x, w);                      // tail, over-writes into the slack
        return end;
    }
    static uint32_t* expand(const std::vector<uint32_t>& ranges, uint32_t* out) { uint32_t* w = out; for (size_t j = 0; j < ranges.size(); j += 2) w = fill(ranges[j], ranges[j + 1], w); return w; }
    static uint32_t* expand(const Runs& r, uint32_t* out) { uint32_t* w = fill(r.lo0, r.hi0, out); for (uint32_t j = 0; j < r.nmid; ++j) w = fill(r.mid[2 * static_cast<size_t>(j)], r.mid[2 * static_cast<size_t>(j) + 1], w); return fill(r.lo1, r.hi1, w); }
    static void expand(const std::vector<uint32_t>& ranges, std::vector<uint32_t>& out) { const size_t t = total(ranges); out.resize(t + kSlack); expand(ranges, out.data()); out.resize(t); }
    bool member(uint32_t u, uint32_t v, int s, const T& k) const {
        const uint32_t xv = own_node(chain_of(v), s), xu = own_node(chain_of(u), s); if (xv == kNone || xu == kNone) return false;
        const uint32_t node = climb(s, xv, k); return xu >= node && xu < node + layers[s].size[node];
    }
    // ladder of v at size s: (level, community vertex count) for every nucleus containing v, from its own level upward
    void ladder(uint32_t v, int s, std::vector<std::pair<T, uint64_t>>& out) const {
        out.clear(); const uint32_t x = own_node(chain_of(v), s); if (x == kNone) return;
        const Layer& L = layers[s];
        for (uint32_t node = x; node != kNone; node = L.parent[node]) {
            uint64_t count = 0;
            if (compact) { Runs r; node_runs(L, node, r); count = total(r); }
            else { uint32_t b, e; slice_bounds(L, node, b, e); for (uint32_t j = b; j < e; ++j) { const uint32_t cc = L.slice[j]; count += start_pos[cc + 1] - start_pos[cc]; } }
            out.emplace_back(top_at(L, node), count);
        }
    }
    // build form -> compact form: (1) merge consecutive vertex ranges along each DFS array into maximal runs, record every node's
    // entry point, drop the chain ids; (2) pack node tops and residues at per-size byte widths
    void compact_runs(bool pack_tops = true) {
        if (compact) return; packed_tops = pack_tops;
        for (int s = 2; s <= max_size; ++s) { Layer& L = layers[s]; const size_t N = L.top.size(); L.entry.assign(2 * (N + 1), 0); L.runs.clear();
            size_t x = 0;   // preorder => bucket non-decreasing; a node's entry is the run and vertex at its segment start
            for (size_t j = 0; j < L.slice.size(); ++j) {
                const uint32_t cc = L.slice[j], lo = start_pos[cc], hi = start_pos[cc + 1];
                if (!L.runs.empty() && L.runs.back() == lo) L.runs.back() = hi; else { L.runs.push_back(lo); L.runs.push_back(hi); }
                const uint32_t r = static_cast<uint32_t>(L.runs.size() / 2 - 1);
                while (x < N && L.bucket[x] == j) { L.entry[2 * x] = r; L.entry[2 * x + 1] = lo; ++x; } }
            const uint32_t M = static_cast<uint32_t>(L.runs.size() / 2);
            for (; x <= N; ++x) { L.entry[2 * x] = M; L.entry[2 * x + 1] = 0; }   // sentinel (and any empty trailing segment)
            L.slice.clear(); L.slice.shrink_to_fit(); L.bucket.clear(); L.bucket.shrink_to_fit();
            if (pack_tops) { unsigned w = 1; for (const T& t : L.top) w = std::max(w, Traits<T>::width(t)); w = round_width(w); L.top_w = static_cast<uint8_t>(w);
                L.top_bytes.assign(N * w, 0); for (size_t i = 0; i < N; ++i) Traits<T>::put(L.top[i], L.top_bytes.data() + i * w, w);
                L.top.clear(); L.top.shrink_to_fit(); } }
        // residues: the widest stored value per size, then one packed block per chain
        res_w.assign(max_size + 1, 1); res_prefix.assign(max_size + 2, 0);
        for (uint32_t c = 0; c < chains; ++c) for (int s = 2; s < sigma[c]; ++s) res_w[s] = static_cast<uint8_t>(std::max<unsigned>(res_w[s], Traits<T>::width(residue[residue_off[c] + (s - 2)])));
        for (int s = 2; s <= max_size; ++s) res_w[s] = static_cast<uint8_t>(round_width(res_w[s]));
        for (int s = 2; s <= max_size; ++s) res_prefix[s + 1] = res_prefix[s] + res_w[s];
        residue_boff.assign(chains + 1, 0);
        for (uint32_t c = 0; c < chains; ++c) { const uint64_t bytes = residue_boff[c] + (sigma[c] > 2 ? res_prefix[sigma[c]] : 0); if (bytes > UINT32_MAX) throw std::overflow_error("residue block exceeds 4 GB"); residue_boff[c + 1] = static_cast<uint32_t>(bytes); }
        residue_bytes.assign(residue_boff[chains], 0);
        for (uint32_t c = 0; c < chains; ++c) for (int s = 2; s < sigma[c]; ++s) Traits<T>::put(residue[residue_off[c] + (s - 2)], residue_bytes.data() + residue_boff[c] + res_prefix[s], res_w[s]);
        residue.clear(); residue.shrink_to_fit(); residue_off.clear(); residue_off.shrink_to_fit();
        compact = true;
    }
    uint64_t runs_total() const { uint64_t r = 0; for (const auto& L : layers) r += L.runs.size() / 2; return r; }
    uint64_t pairs_total() const { uint64_t r = 0; for (const auto& L : layers) r += L.slice.size(); return r; }
    // ---- bytes (in memory; jump pointers are derived and not stored on disk)
    uint64_t bytes_map() const { return 8ull * start_bits.size() + 4ull * start_cum.size() + 4ull * start_pos.size(); }
    uint64_t bytes_chains() const { return omega.size() + sigma.size() + 4ull * traj_off.size() + 4ull * traj_node.size() + 4ull * residue_off.size() + Traits<T>::W * residue.size()
        + res_w.size() + 4ull * res_prefix.size() + 4ull * residue_boff.size() + residue_bytes.size(); }
    uint64_t bytes_layers() const { uint64_t b = 0; for (const auto& L : layers) b += Traits<T>::W * L.top.size() + L.top_bytes.size() + 4ull * (L.parent.size() + L.size.size() + L.jump.size() + L.bucket.size() + L.slice.size() + L.entry.size() + L.runs.size()); return b; }
    uint64_t bytes_total() const { return bytes_map() + bytes_chains() + bytes_layers(); }
    uint64_t node_count() const { uint64_t c = 0; for (const auto& L : layers) c += L.size.size(); return c; }

    // ---- jump pointers (Myers' single skip pointer per node): built from parent links, nodes in DFS preorder (parents before children)
    static void build_jumps(Layer& L) {
        const size_t N = L.size.size(); L.jump.assign(N, kNone); std::vector<uint32_t> depth(N, 0);
        for (uint32_t x = 0; x < N; ++x) {
            const uint32_t p = L.parent[x]; if (p == kNone) { depth[x] = 0; L.jump[x] = kNone; continue; }
            depth[x] = depth[p] + 1;
            const uint32_t jp = L.jump[p];
            if (jp != kNone && L.jump[jp] != kNone && depth[p] - depth[jp] == depth[jp] - depth[L.jump[jp]]) L.jump[x] = L.jump[jp];
            else L.jump[x] = p;
        }
    }
    void finish() {
        for (auto& L : layers) if (!L.size.empty()) build_jumps(L);
        binom.assign(max_size + 1, std::vector<T>(max_size + 1, T{0}));
        for (int a = 0; a <= max_size; ++a) { binom[a][0] = T{1}; for (int b = 1; b <= a; ++b) binom[a][b] = binom[a - 1][b - 1] + (b <= a - 1 ? binom[a - 1][b] : T{0}); }
    }

    // ---- disk format (compact form): magic, header, then arrays as (u64 count, raw bytes)
    template<class V> static void wv(std::ofstream& f, const std::vector<V>& v) { const uint64_t c = v.size(); f.write(reinterpret_cast<const char*>(&c), 8); if (c) f.write(reinterpret_cast<const char*>(v.data()), c * sizeof(V)); }
    template<class V> static void rv(std::ifstream& f, std::vector<V>& v) { uint64_t c = 0; f.read(reinterpret_cast<char*>(&c), 8); v.resize(c); if (c) f.read(reinterpret_cast<char*>(v.data()), c * sizeof(V)); }
    static void wt(std::ofstream& f, const std::vector<T>& v) { const uint64_t c = v.size(); f.write(reinterpret_cast<const char*>(&c), 8); std::vector<unsigned char> buf(c * Traits<T>::W); for (size_t i = 0; i < c; ++i) Traits<T>::put(v[i], buf.data() + i * Traits<T>::W); if (c) f.write(reinterpret_cast<const char*>(buf.data()), buf.size()); }
    static void rt(std::ifstream& f, std::vector<T>& v) { uint64_t c = 0; f.read(reinterpret_cast<char*>(&c), 8); std::vector<unsigned char> buf(c * Traits<T>::W); if (c) f.read(reinterpret_cast<char*>(buf.data()), buf.size()); v.resize(c); for (size_t i = 0; i < c; ++i) v[i] = Traits<T>::get(buf.data() + i * Traits<T>::W); }
    void save(const std::string& path) const {
        if (!compact) throw std::runtime_error("save needs the compact form (call compact_runs first)");
        std::ofstream f(path, std::ios::binary); if (!f) throw std::runtime_error("cannot write " + path);
        const char magic[8] = {'C','H','A','I','N','X','0','3'}; f.write(magic, 8);
        const uint32_t hdr[4] = {n, chains, static_cast<uint32_t>(max_size), static_cast<uint32_t>(Traits<T>::W)}; f.write(reinterpret_cast<const char*>(hdr), 16);
        const uint8_t flags = packed_tops ? 1 : 0; f.write(reinterpret_cast<const char*>(&flags), 1);
        wv(f, start_bits); wv(f, start_cum); wv(f, start_pos); wv(f, omega); wv(f, sigma); wv(f, traj_off); wv(f, traj_node); wv(f, res_w); wv(f, res_prefix); wv(f, residue_boff); wv(f, residue_bytes);
        for (int s = 2; s <= max_size; ++s) { const Layer& L = layers[s]; if (packed_tops) { f.write(reinterpret_cast<const char*>(&L.top_w), 1); wv(f, L.top_bytes); } else wt(f, L.top); wv(f, L.parent); wv(f, L.size); wv(f, L.entry); wv(f, L.runs); }
        if (!f) throw std::runtime_error("write failed " + path);
    }
    static ChainIndex load(const std::string& path) {
        std::ifstream f(path, std::ios::binary); if (!f) throw std::runtime_error("cannot read " + path);
        char magic[8]; f.read(magic, 8); if (std::memcmp(magic, "CHAINX03", 8) != 0) throw std::runtime_error("bad magic");
        uint32_t hdr[4]; f.read(reinterpret_cast<char*>(hdr), 16); ChainIndex ix; ix.n = hdr[0]; ix.chains = hdr[1]; ix.max_size = static_cast<int>(hdr[2]);
        if (hdr[3] != Traits<T>::W) throw std::runtime_error("count width mismatch");
        uint8_t flags = 0; f.read(reinterpret_cast<char*>(&flags), 1); ix.packed_tops = (flags & 1) != 0;
        rv(f, ix.start_bits); rv(f, ix.start_cum); rv(f, ix.start_pos); rv(f, ix.omega); rv(f, ix.sigma); rv(f, ix.traj_off); rv(f, ix.traj_node); rv(f, ix.res_w); rv(f, ix.res_prefix); rv(f, ix.residue_boff); rv(f, ix.residue_bytes);
        ix.layers.resize(ix.max_size + 1);
        for (int s = 2; s <= ix.max_size; ++s) { Layer& L = ix.layers[s]; if (ix.packed_tops) { f.read(reinterpret_cast<char*>(&L.top_w), 1); rv(f, L.top_bytes); } else rt(f, L.top); rv(f, L.parent); rv(f, L.size); rv(f, L.entry); rv(f, L.runs); }
        if (!f) throw std::runtime_error("read failed " + path);
        ix.compact = true; ix.finish(); return ix;
    }
};

} // namespace chainindex
