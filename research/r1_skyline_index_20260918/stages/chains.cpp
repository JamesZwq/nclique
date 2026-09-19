// Hierarchy-equivalence chains: group vertices by the tuple of their own
// canonical nodes (X_2(v), ..., X_omega(v)).  Counts distinct chains and
// projects index bytes for S trees over chains.  Counting only.
#define main skyline_stage1_count_main
#include "../count.cpp"
#undef main
#include <unordered_map>

template<class T> static void chains(const Input& in, unsigned bits) {
    const Graph& g = in.graph; const Vertex n = g.n; const int S = std::max(2, static_cast<int>(in.d) + 1);
    Layout layout(g, S); layout.prepare(n); typename Kernel<T>::Combinations choose(in.d + 1, S);
    terminal::Index ti(S); terminal::build(g, ti, 0); ti.prepare(n);
    auto out = terminal::Solver<T>::solve(g, ti, choose, in.ordinary); const auto& core = out.common.data.core;
    std::vector<int> omega(n, 0), sigma(n, 0);
    for (Vertex v = 0; v < n; ++v) { for (int s = 2; s <= S; ++s) if (core[static_cast<size_t>(s) * n + v] > T{0}) omega[v] = s;
        sigma[v] = omega[v] + 1; for (int s = 2; s <= omega[v]; ++s) if (sigma[v] == omega[v] + 1 && cpp_int(core[static_cast<size_t>(s) * n + v]) == choose_int(omega[v] - 1, s - 1)) sigma[v] = s; }
    std::vector<std::vector<int>> tuple(n); uint64_t nodes = 0;
    for (int s = 2; s <= S; ++s) { auto tr = make_tree(g, ti, core, s); nodes += tr.nodes.size();
        for (Vertex v = 0; v < n; ++v) if (tr.leaf[v] >= 0) tuple[v].push_back(tr.leaf[v]); }
    std::map<std::vector<int>, uint32_t> ids; std::vector<uint32_t> chain(n); std::vector<uint64_t> members;
    uint64_t active = 0;
    for (Vertex v = 0; v < n; ++v) { if (omega[v] < 2) continue; ++active; auto [it, fresh] = ids.emplace(tuple[v], static_cast<uint32_t>(ids.size())); chain[v] = it->second; if (fresh) members.push_back(0); ++members[it->second]; }
    const uint64_t C = ids.size();
    uint64_t pairs_v = 0, pairs_c = 0, residue_v = 0, residue_c = 0, twin_classes = 0; uint64_t largest = 0;
    { std::vector<std::vector<Vertex>> groups; twins(g, groups); twin_classes = groups.size(); }
    std::vector<uint8_t> seen(C, 0);
    for (Vertex v = 0; v < n; ++v) if (omega[v] >= 2) { pairs_v += omega[v] - 1; residue_v += sigma[v] - 2;
        if (!seen[chain[v]]) { seen[chain[v]] = 1; pairs_c += omega[v] - 1; residue_c += sigma[v] - 2; } }
    for (auto m : members) largest = std::max(largest, m);
    // consistency: all vertices of a chain share kappa at every size (L1)
    { std::vector<Vertex> rep(C, absent); for (Vertex v = 0; v < n; ++v) if (omega[v] >= 2) { auto& r = rep[chain[v]]; if (r == absent) r = v; else { require(omega[v] == omega[r], "chain omega"); for (int s = 2; s <= omega[v]; ++s) require(core[static_cast<size_t>(s) * n + v] == core[static_cast<size_t>(s) * n + r], "chain kappa"); } } }
    const uint64_t W = bits / 8;
    const uint64_t node_bytes = nodes * (W + 12);
    const uint64_t strees_v = node_bytes + 8 * pairs_v + 4 * (n + 1) + (2 * n + 4 * (n + 1) + W * residue_v);            // per-vertex S trees + Block D
    const uint64_t chains_map = node_bytes + 8 * pairs_c + 4 * (C + 1) + 4 * n + 4 * n + 4 * (C + 1) + (2 * C + 4 * (C + 1) + W * residue_c); // chain id per vertex + chain->vertex CSR
    const uint64_t chains_bitmap = node_bytes + 8 * pairs_c + 4 * (C + 1) + (n + 7) / 8 + 4 * (C + 1) + (2 * C + 4 * (C + 1) + W * residue_c); // aligned labels: bitmap + rank
    std::cout << std::fixed << std::setprecision(3) << "{\"n\":" << n << ",\"active\":" << active << ",\"twin_classes\":" << twin_classes << ",\"chains\":" << C
        << ",\"canonical_nodes\":" << nodes << ",\"largest_chain\":" << largest << ",\"pairs_vertex\":" << pairs_v << ",\"pairs_chain\":" << pairs_c
        << ",\"residue_vertex\":" << residue_v << ",\"residue_chain\":" << residue_c << ",\"bytes_strees_vertex\":" << strees_v
        << ",\"bytes_chains_map\":" << chains_map << ",\"bytes_chains_bitmap\":" << chains_bitmap
        << ",\"ratio_map\":" << double(strees_v) / chains_map << ",\"ratio_bitmap\":" << double(strees_v) / chains_bitmap << "}\n";
}
int main(int argc, char** argv) {
    try { require(argc == 2, "usage: chains <graph>"); Input in = prepare(argv[1]); Layout l(in.graph, std::max(2, static_cast<int>(in.d) + 1)); l.prepare(in.graph.n);
        const unsigned w = width(count_bound(in.graph, l, in.d));
        if (w == 64) chains<uint64_t>(in, w); else if (w == 128) chains<unsigned __int128>(in, w); else if (w == 256) chains<boost::multiprecision::uint256_t>(in, w); else chains<boost::multiprecision::uint512_t>(in, w);
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
