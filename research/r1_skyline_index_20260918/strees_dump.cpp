// Writes the per-vertex S trees of a graph (one tree and one DFS array per size, the storage baseline of the paper) as a
// flat binary file, materialized from a chain index file, so that a general-purpose compressor can be run on it.
// Byte accounting is the one of chain_index_tool (baseline_vertex_bytes) with 8-byte values: per size, every node has
// value (8), parent (4), first and last DFS position (4 + 4); every active (vertex, size) pair has its DFS-array entry (4)
// and its own-node pointer (4); per vertex, offsets (4 (n + 1)), omega and sigma (1 + 1), residue offsets (8 (n + 1)) and
// the residues (8 each).  Usage: strees_dump <index.cx> <out.bin> [<index.cx>.perm]   -> one JSON line with the byte counts.
// Without the permutation file the vertex ids are the index's aligned labels; with it (perm[degeneracy label] = aligned
// label, written by the tool) they are the labels of the decomposition's input, the degeneracy order of the graph.
#include "../../src-r1index/chain_index.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
using namespace chainindex;
using Index = ChainIndex<double>;

int main(int argc, char** argv) {
    if (argc != 3 && argc != 4) { std::cerr << "usage: strees_dump <index.cx> <out.bin> [<index.cx>.perm]\n"; return 2; }
    Index ix = Index::load(argv[1]); const uint32_t n = ix.n;
    std::vector<uint32_t> perm(n), inv(n); for (uint32_t v = 0; v < n; ++v) perm[v] = inv[v] = v;   // identity: aligned labels
    if (argc == 4) { std::ifstream pf(argv[3], std::ios::binary); pf.read(reinterpret_cast<char*>(perm.data()), static_cast<std::streamsize>(n) * 4);
        if (!pf) { std::cerr << "cannot read permutation\n"; return 2; } for (uint32_t v = 0; v < n; ++v) inv[perm[v]] = v; }
    std::ofstream out(argv[2], std::ios::binary); uint64_t bytes = 0, node_bytes = 0, pair_bytes = 0, vertex_bytes = 0, residue_bytes = 0;
    auto put = [&](const void* p, size_t k) { out.write(static_cast<const char*>(p), static_cast<std::streamsize>(k)); bytes += k; };
    std::vector<uint32_t> dfs, own;
    for (int s = 2; s <= ix.max_size; ++s) {
        const Index::Layer& L = ix.layers[s]; const uint32_t N = static_cast<uint32_t>(L.size.size()); const size_t R = L.runs.size() / 2;
        std::vector<uint64_t> pre(R + 1, 0); for (size_t j = 0; j < R; ++j) pre[j + 1] = pre[j] + (L.runs[2 * j + 1] - L.runs[2 * j]);
        auto pos = [&](uint32_t r, uint32_t v) { return r < R ? pre[r] + (v - L.runs[2 * static_cast<size_t>(r)]) : pre[R]; };
        for (uint32_t x = 0; x < N; ++x) {                                  // nodes in preorder: value, parent, [lo, hi) in the DFS array
            const double val = ix.top_at(L, x); const uint32_t par = L.parent[x];
            const size_t e = static_cast<size_t>(x) + L.size[x];
            const uint32_t lo = static_cast<uint32_t>(pos(L.entry[2 * static_cast<size_t>(x)], L.entry[2 * static_cast<size_t>(x) + 1]));
            const uint32_t hi = static_cast<uint32_t>(pos(L.entry[2 * e], L.entry[2 * e + 1]));
            put(&val, 8); put(&par, 4); put(&lo, 4); put(&hi, 4); node_bytes += 20;
        }
        dfs.clear(); for (size_t j = 0; j < R; ++j) for (uint32_t v = L.runs[2 * j]; v < L.runs[2 * j + 1]; ++v) dfs.push_back(inv[v]);   // the DFS array of size s, in the chosen labels
        own.resize(dfs.size()); for (size_t i = 0; i < dfs.size(); ++i) own[i] = ix.own_node(ix.chain_of(perm[dfs[i]]), s);              // own node of every active vertex
        put(dfs.data(), dfs.size() * 4); put(own.data(), own.size() * 4); pair_bytes += dfs.size() * 8;
    }
    std::vector<uint32_t> off(n + 1, 0); std::vector<uint8_t> om(n), sg(n); std::vector<uint64_t> roff(n + 1, 0); std::vector<double> res;
    for (uint32_t v = 0; v < n; ++v) { const uint32_t c = ix.chain_of(perm[v]); om[v] = static_cast<uint8_t>(ix.omega[c]); sg[v] = static_cast<uint8_t>(ix.sigma[c]);
        off[v + 1] = off[v] + (ix.omega[c] >= 2 ? ix.omega[c] - 1 : 0); roff[v + 1] = roff[v] + (ix.omega[c] >= 2 ? ix.sigma[c] - 2 : 0);
        for (int s = 2; s < ix.sigma[c] && s <= ix.omega[c]; ++s) res.push_back(ix.value(perm[v], s)); }
    put(off.data(), off.size() * 4); put(om.data(), n); put(sg.data(), n); vertex_bytes += off.size() * 4 + 2ull * n;
    put(roff.data(), roff.size() * 8); vertex_bytes += roff.size() * 8; put(res.data(), res.size() * 8); residue_bytes += res.size() * 8;
    out.close();
    std::cout << "{\"labels\":\"" << (argc == 4 ? "degeneracy" : "aligned") << "\",\"n\":" << n << ",\"s_max\":" << ix.max_size << ",\"strees_bytes\":" << bytes << ",\"node_bytes\":" << node_bytes << ",\"pair_bytes\":" << pair_bytes
              << ",\"vertex_bytes\":" << vertex_bytes << ",\"residue_bytes\":" << residue_bytes << ",\"index_bytes\":" << ix.bytes_total() << "}\n";
    return 0;
}
