// ChainIndex: the final r = 1 all-size community/value index (CHAINS.md,
// RESULTS_CHAINS.md, RESULTS_FINAL.md).  Vertices carry aligned labels: every
// chain is one id range, so the vertex-to-chain map is a bitmap with rank.
// Flat arrays only; one file on disk; queries are constant-step walks plus
// output.  Header-only, templated on the count type T (W = sizeof(T)).
#pragma once
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

namespace chainindex {

static constexpr uint32_t kNone = UINT32_MAX;

template<class T> struct Traits {                       // raw little-endian bytes <-> T
    static constexpr size_t W = sizeof(T);
    static void put(const T& x, unsigned char* p) { std::memcpy(p, &x, W); }
    static T get(const unsigned char* p) { T x; std::memcpy(&x, p, W); return x; }
};
template<unsigned Bits> struct Traits<boost::multiprecision::number<boost::multiprecision::cpp_int_backend<Bits, Bits, boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>> {
    using T = boost::multiprecision::number<boost::multiprecision::cpp_int_backend<Bits, Bits, boost::multiprecision::unsigned_magnitude, boost::multiprecision::unchecked, void>>;
    static constexpr size_t W = Bits / 8;
    static void put(const T& x, unsigned char* p) { std::memset(p, 0, W); std::vector<unsigned char> tmp; boost::multiprecision::export_bits(x, std::back_inserter(tmp), 8, false); std::memcpy(p, tmp.data(), std::min(tmp.size(), W)); }
    static T get(const unsigned char* p) { T x; boost::multiprecision::import_bits(x, p, p + W, 8, false); return x; }
};

template<class T> struct ChainIndex {
    // header
    uint32_t n = 0, chains = 0; int max_size = 2;
    // block 1: vertex -> chain (bitmap of chain starts + rank directory), chain -> vertex range
    std::vector<uint64_t> start_bits; std::vector<uint32_t> start_cum; std::vector<uint32_t> start_pos;   // start_pos[chains] = n
    // block 2: per chain
    std::vector<uint8_t> omega, sigma; std::vector<uint32_t> traj_off, traj_node;   // traj_node[traj_off[c] + (s-2)] = own node of chain c at size s
    std::vector<uint32_t> residue_off; std::vector<T> residue;                     // residue values for 2 <= s < sigma(c)
    // block 3: per size s (index s, 0 and 1 unused)
    struct Layer { std::vector<T> top; std::vector<uint32_t> parent, size, jump, bucket, slice; };   // slice: chain ids in DFS order, own-first
    std::vector<Layer> layers;
    // binomials for the certified tail, C(a, b) for a <= max_size
    std::vector<std::vector<T>> binom;

    // ---- derived accessors
    uint32_t chain_of(uint32_t v) const { const uint32_t w = v >> 6; const uint64_t mask = (v & 63) == 63 ? ~0ull : ((2ull << (v & 63)) - 1);
        return start_cum[w] + static_cast<uint32_t>(__builtin_popcountll(start_bits[w] & mask)) - 1; }
    uint32_t own_node(uint32_t c, int s) const { return s > omega[c] || s < 2 ? kNone : traj_node[traj_off[c] + static_cast<size_t>(s - 2)]; }
    uint32_t climb(int s, uint32_t x, const T& k) const {                // highest ancestor with top >= k; tops decrease upward
        const Layer& L = layers[s];
        while (L.parent[x] != kNone && L.top[L.parent[x]] >= k) {
            const uint32_t j = L.jump[x];
            if (j != kNone && L.top[j] >= k) x = j; else x = L.parent[x];
        }
        return x;
    }
    // ---- queries (internal labels)
    T value(uint32_t v, int s) const {
        const uint32_t c = chain_of(v);
        if (s < 2 || s > omega[c]) return T{0};
        if (s >= sigma[c]) return binom[omega[c] - 1][s - 1];
        return residue[residue_off[c] + static_cast<size_t>(s - 2)];
    }
    // community of v at (s, k) as vertex-id ranges appended to `ranges` (pairs lo, hi); returns the node or kNone
    uint32_t community_ranges(uint32_t v, int s, const T& k, std::vector<uint32_t>& ranges) const {
        const uint32_t c = chain_of(v); const uint32_t x = own_node(c, s); if (x == kNone) return kNone;
        const Layer& L = layers[s]; const uint32_t node = climb(s, x, k);
        const uint32_t b = L.bucket[node], hi = node + L.size[node], e = hi < L.bucket.size() ? L.bucket[hi] : static_cast<uint32_t>(L.slice.size());
        for (uint32_t j = b; j < e; ++j) { const uint32_t cc = L.slice[j]; ranges.push_back(start_pos[cc]); ranges.push_back(start_pos[cc + 1]); }
        return node;
    }
    static void expand(const std::vector<uint32_t>& ranges, std::vector<uint32_t>& out) {
        size_t total = 0; for (size_t j = 0; j < ranges.size(); j += 2) total += ranges[j + 1] - ranges[j];
        out.resize(total); uint32_t* w = out.data();
        for (size_t j = 0; j < ranges.size(); j += 2) for (uint32_t x = ranges[j], hi = ranges[j + 1]; x < hi; ++x) *w++ = x;
    }
    bool member(uint32_t u, uint32_t v, int s, const T& k) const {
        const uint32_t xv = own_node(chain_of(v), s), xu = own_node(chain_of(u), s); if (xv == kNone || xu == kNone) return false;
        const uint32_t node = climb(s, xv, k); return xu >= node && xu < node + layers[s].size[node];
    }
    // ladder of v at size s: (level, community vertex count) for every nucleus containing v, from its own level upward
    void ladder(uint32_t v, int s, std::vector<std::pair<T, uint64_t>>& out) const {
        out.clear(); const uint32_t x = own_node(chain_of(v), s); if (x == kNone) return;
        const Layer& L = layers[s];
        for (uint32_t node = x; node != kNone; node = L.parent[node]) {
            const uint32_t b = L.bucket[node], hi = node + L.size[node], e = hi < L.bucket.size() ? L.bucket[hi] : static_cast<uint32_t>(L.slice.size());
            uint64_t count = 0; for (uint32_t j = b; j < e; ++j) { const uint32_t cc = L.slice[j]; count += start_pos[cc + 1] - start_pos[cc]; }
            out.emplace_back(L.top[node], count);
        }
    }
    // ---- bytes
    uint64_t bytes_map() const { return 8ull * start_bits.size() + 4ull * start_cum.size() + 4ull * start_pos.size(); }
    uint64_t bytes_chains() const { return omega.size() + sigma.size() + 4ull * traj_off.size() + 4ull * traj_node.size() + 4ull * residue_off.size() + Traits<T>::W * residue.size(); }
    uint64_t bytes_layers() const { uint64_t b = 0; for (const auto& L : layers) b += Traits<T>::W * L.top.size() + 4ull * (L.parent.size() + L.size.size() + L.jump.size() + L.bucket.size() + L.slice.size()); return b; }
    uint64_t bytes_total() const { return bytes_map() + bytes_chains() + bytes_layers(); }
    uint64_t node_count() const { uint64_t c = 0; for (const auto& L : layers) c += L.top.size(); return c; }

    // ---- jump pointers (Myers' single skip pointer per node): built from parent links, nodes in DFS preorder (parents before children)
    static void build_jumps(Layer& L) {
        const size_t N = L.top.size(); L.jump.assign(N, kNone); std::vector<uint32_t> depth(N, 0);
        for (uint32_t x = 0; x < N; ++x) {
            const uint32_t p = L.parent[x]; if (p == kNone) { depth[x] = 0; L.jump[x] = kNone; continue; }
            depth[x] = depth[p] + 1;
            const uint32_t jp = L.jump[p];
            if (jp != kNone && L.jump[jp] != kNone && depth[p] - depth[jp] == depth[jp] - depth[L.jump[jp]]) L.jump[x] = L.jump[jp];
            else L.jump[x] = p;
        }
    }
    void finish() {
        for (auto& L : layers) if (!L.top.empty()) build_jumps(L);
        binom.assign(max_size + 1, std::vector<T>(max_size + 1, T{0}));
        for (int a = 0; a <= max_size; ++a) { binom[a][0] = T{1}; for (int b = 1; b <= a; ++b) binom[a][b] = binom[a - 1][b - 1] + (b <= a - 1 ? binom[a - 1][b] : T{0}); }
    }

    // ---- disk format: magic, header, then arrays as (u64 count, raw bytes)
    template<class V> static void wv(std::ofstream& f, const std::vector<V>& v) { const uint64_t c = v.size(); f.write(reinterpret_cast<const char*>(&c), 8); if (c) f.write(reinterpret_cast<const char*>(v.data()), c * sizeof(V)); }
    template<class V> static void rv(std::ifstream& f, std::vector<V>& v) { uint64_t c = 0; f.read(reinterpret_cast<char*>(&c), 8); v.resize(c); if (c) f.read(reinterpret_cast<char*>(v.data()), c * sizeof(V)); }
    static void wt(std::ofstream& f, const std::vector<T>& v) { const uint64_t c = v.size(); f.write(reinterpret_cast<const char*>(&c), 8); std::vector<unsigned char> buf(c * Traits<T>::W); for (size_t i = 0; i < c; ++i) Traits<T>::put(v[i], buf.data() + i * Traits<T>::W); if (c) f.write(reinterpret_cast<const char*>(buf.data()), buf.size()); }
    static void rt(std::ifstream& f, std::vector<T>& v) { uint64_t c = 0; f.read(reinterpret_cast<char*>(&c), 8); std::vector<unsigned char> buf(c * Traits<T>::W); if (c) f.read(reinterpret_cast<char*>(buf.data()), buf.size()); v.resize(c); for (size_t i = 0; i < c; ++i) v[i] = Traits<T>::get(buf.data() + i * Traits<T>::W); }
    void save(const std::string& path) const {
        std::ofstream f(path, std::ios::binary); if (!f) throw std::runtime_error("cannot write " + path);
        const char magic[8] = {'C','H','A','I','N','X','0','1'}; f.write(magic, 8);
        const uint32_t hdr[4] = {n, chains, static_cast<uint32_t>(max_size), static_cast<uint32_t>(Traits<T>::W)}; f.write(reinterpret_cast<const char*>(hdr), 16);
        wv(f, start_bits); wv(f, start_cum); wv(f, start_pos); wv(f, omega); wv(f, sigma); wv(f, traj_off); wv(f, traj_node); wv(f, residue_off); wt(f, residue);
        for (int s = 2; s <= max_size; ++s) { const Layer& L = layers[s]; wt(f, L.top); wv(f, L.parent); wv(f, L.size); wv(f, L.bucket); wv(f, L.slice); }
        if (!f) throw std::runtime_error("write failed " + path);
    }
    static ChainIndex load(const std::string& path) {
        std::ifstream f(path, std::ios::binary); if (!f) throw std::runtime_error("cannot read " + path);
        char magic[8]; f.read(magic, 8); if (std::memcmp(magic, "CHAINX01", 8) != 0) throw std::runtime_error("bad magic");
        uint32_t hdr[4]; f.read(reinterpret_cast<char*>(hdr), 16); ChainIndex ix; ix.n = hdr[0]; ix.chains = hdr[1]; ix.max_size = static_cast<int>(hdr[2]);
        if (hdr[3] != Traits<T>::W) throw std::runtime_error("count width mismatch");
        rv(f, ix.start_bits); rv(f, ix.start_cum); rv(f, ix.start_pos); rv(f, ix.omega); rv(f, ix.sigma); rv(f, ix.traj_off); rv(f, ix.traj_node); rv(f, ix.residue_off); rt(f, ix.residue);
        ix.layers.resize(ix.max_size + 1);
        for (int s = 2; s <= ix.max_size; ++s) { Layer& L = ix.layers[s]; rt(f, L.top); rv(f, L.parent); rv(f, L.size); rv(f, L.bucket); rv(f, L.slice); }
        if (!f) throw std::runtime_error("read failed " + path);
        ix.finish(); return ix;
    }
};

} // namespace chainindex
