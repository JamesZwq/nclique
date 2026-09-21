// Case study: which clique size is the right one?  Ground-truth communities (SNAP top-5000) against the (s, k)-nuclei of
// their members, for every size s and every level k of the chain index.  For every community C and every member v:
// at every size s in [2, omega(v)] the ladder of nuclei containing v (own node, then its ancestors) is scored by F1
// against C, where |A cap C| is counted through the members of C (own node of u at s inside A's subtree) and |A| comes
// from the run bounds, so no community is ever listed.  The climb stops when the bound 2|C| / (|A| + |C|) falls below
// the best F1 found for the query.  Also a k-core sanity check of the file-label permutation.
// Usage: case_groundtruth <g.edges> <g.cx> <g.cmty>   -> one JSON line.
#include "../../src-r1index/chain_index.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace chainindex;
using Index = ChainIndex<double>; using Clock = std::chrono::steady_clock;

static std::vector<uint32_t> core_numbers(const std::string& path, uint32_t& n) {   // k-core numbers of the file graph (bucket peel)
    std::ifstream f(path); uint64_t m; f >> n >> m; std::vector<std::pair<uint32_t, uint32_t>> e(m); for (auto& [a, b] : e) f >> a >> b;
    std::vector<uint32_t> deg(n, 0), off(n + 1, 0); for (auto& [a, b] : e) { ++deg[a]; ++deg[b]; }
    for (uint32_t v = 0; v < n; ++v) off[v + 1] = off[v] + deg[v]; std::vector<uint32_t> adj(off[n]), fill(off.begin(), off.end() - 1);
    for (auto& [a, b] : e) { adj[fill[a]++] = b; adj[fill[b]++] = a; }
    uint32_t md = 0; for (uint32_t v = 0; v < n; ++v) md = std::max(md, deg[v]);
    std::vector<uint32_t> bin(md + 2, 0), pos(n), vert(n), d = deg; for (uint32_t v = 0; v < n; ++v) ++bin[d[v]];
    uint32_t start = 0; for (uint32_t k = 0; k <= md; ++k) { const uint32_t c = bin[k]; bin[k] = start; start += c; }
    for (uint32_t v = 0; v < n; ++v) { pos[v] = bin[d[v]]; vert[pos[v]] = v; ++bin[d[v]]; }
    for (uint32_t k = md + 1; k-- > 1;) bin[k] = bin[k - 1]; bin[0] = 0;
    for (uint32_t i = 0; i < n; ++i) { const uint32_t v = vert[i];
        for (uint32_t j = off[v]; j < off[v + 1]; ++j) { const uint32_t u = adj[j]; if (d[u] > d[v]) { const uint32_t du = d[u], pu = pos[u], pw = bin[du], w = vert[pw];
            if (u != w) { pos[u] = pw; vert[pu] = w; pos[w] = pu; vert[pw] = u; } ++bin[du]; --d[u]; } } }
    return d;
}

int main(int argc, char** argv) {
    if (argc != 4) { std::cerr << "usage: case_groundtruth <g.edges> <g.cx> <g.cmty>\n"; return 2; }
    const Index ix = Index::load(argv[2]); const uint32_t n = ix.n;
    std::vector<uint32_t> perm(n); { std::ifstream pf(std::string(argv[2]) + ".perm", std::ios::binary); pf.read(reinterpret_cast<char*>(perm.data()), static_cast<std::streamsize>(n) * 4); if (!pf) { std::cerr << "no perm\n"; return 2; } }
    { uint32_t nf; const auto core = core_numbers(argv[1], nf); if (nf != n) { std::cerr << "n mismatch\n"; return 2; }   // sanity: value(., 2) is the k-core number
      for (uint32_t v = 0; v < n; ++v) if (ix.value(perm[v], 2) != double(core[v])) { std::cerr << "k-core mismatch at " << v << "\n"; return 2; } }
    std::vector<std::vector<uint32_t>> cmty; { std::ifstream cf(argv[3]); std::string line; while (std::getline(cf, line)) { std::istringstream ss(line); std::vector<uint32_t> c; uint32_t x; while (ss >> x) c.push_back(perm[x]); if (c.size() >= 2) cmty.push_back(c); } }
    const int S = ix.max_size; std::vector<uint32_t> own_u;   // own nodes of the members of C at the current size
    // aggregates: fixed-size own level, fixed-size best level, best over sizes (own), best over sizes and levels; histogram of the best size
    const int F = std::min(S, 12); std::vector<double> sum_own(F + 1, 0), sum_lad(F + 1, 0); std::vector<uint64_t> cnt_fixed(F + 1, 0);
    double sum_best_own = 0, sum_best_all = 0; uint64_t queries = 0, ladder_steps = 0; std::vector<uint64_t> hist_best(S + 1, 0), hist_best_own(S + 1, 0);
    uint64_t improved_over_core = 0, improved_over_core_ladder = 0;   // best (s,k) strictly better than the best k-core (s = 2, best k)
    const auto t0 = Clock::now();
    for (const auto& C : cmty) {
        const double csz = static_cast<double>(C.size());
        for (uint32_t v : C) {
            const uint32_t cv = ix.chain_of(v); const int ov = ix.omega[cv]; if (ov < 2) continue; ++queries;
            double best_all = 0, best_own = 0; int best_s = 2, best_s_own = 2; std::vector<double> lad_at(F + 1, 0), own_at(F + 1, 0);
            for (int s = 2; s <= ov; ++s) {
                const Index::Layer& L = ix.layers[s]; own_u.resize(C.size());
                for (size_t i = 0; i < C.size(); ++i) own_u[i] = ix.own_node(ix.chain_of(C[i]), s);
                uint32_t A = ix.own_node(cv, s); double best_here = 0; bool first = true;
                while (A != kNone) { Index::Runs r; Index::node_runs(L, A, r); const double sz = static_cast<double>(Index::total(r)); ++ladder_steps;
                    if (2 * csz / (sz + csz) < best_here) break;                     // no ancestor can beat the best of this size
                    const uint32_t hi = A + L.size[A]; uint32_t inter = 0; for (uint32_t x : own_u) inter += (x != kNone && x >= A && x < hi);
                    const double f1 = 2.0 * inter / (sz + csz);
                    if (first) { if (s <= F) own_at[s] = f1; if (f1 > best_own) { best_own = f1; best_s_own = s; } first = false; }
                    if (f1 > best_here) best_here = f1;
                    A = L.parent[A]; }
                if (s <= F) lad_at[s] = best_here;
                if (best_here > best_all) { best_all = best_here; best_s = s; }
            }
            for (int s = 2; s <= std::min(ov, F); ++s) { sum_own[s] += own_at[s]; sum_lad[s] += lad_at[s]; ++cnt_fixed[s]; }
            sum_best_own += best_own; sum_best_all += best_all; ++hist_best[best_s]; ++hist_best_own[best_s_own];
            improved_over_core += best_all > lad_at[2] + 1e-12; improved_over_core_ladder += best_own > own_at[2] + 1e-12;
        }
    }
    const double secs = std::chrono::duration<double>(Clock::now() - t0).count();
    std::cout.precision(4); std::cout << std::fixed << "{\"n\":" << n << ",\"s_max\":" << S << ",\"communities\":" << cmty.size() << ",\"queries\":" << queries << ",\"ladder_steps\":" << ladder_steps << ",\"seconds\":" << secs
              << ",\"mean_best_own\":" << sum_best_own / queries << ",\"mean_best_all\":" << sum_best_all / queries << ",\"better_than_best_core\":" << double(improved_over_core) / queries << ",\"own_better_than_core_own\":" << double(improved_over_core_ladder) / queries;
    std::cout << ",\"fixed\":{"; for (int s = 2; s <= F; ++s) std::cout << (s > 2 ? "," : "") << "\"" << s << "\":{\"queries\":" << cnt_fixed[s] << ",\"mean_own\":" << (cnt_fixed[s] ? sum_own[s] / cnt_fixed[s] : 0) << ",\"mean_ladder\":" << (cnt_fixed[s] ? sum_lad[s] / cnt_fixed[s] : 0) << ",\"mean_own_all\":" << sum_own[s] / queries << ",\"mean_ladder_all\":" << sum_lad[s] / queries << "}";
    std::cout << "},\"best_size_hist\":{"; bool first = true; for (int s = 2; s <= S; ++s) if (hist_best[s]) { std::cout << (first ? "" : ",") << "\"" << s << "\":" << hist_best[s]; first = false; }
    std::cout << "},\"best_size_hist_own\":{"; first = true; for (int s = 2; s <= S; ++s) if (hist_best_own[s]) { std::cout << (first ? "" : ",") << "\"" << s << "\":" << hist_best_own[s]; first = false; }
    std::cout << "}}\n";
    return 0;
}
