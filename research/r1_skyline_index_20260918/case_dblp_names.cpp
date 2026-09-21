// Case study: drill-down by clique size on the DBLP coauthorship graph with author names and venues.  For every anchor
// author and every size s, the own-level community (the (s, kappa_s(v))-nucleus of the anchor) is listed from the chain
// index: its size, the query time, its topical purity (share of members with a paper at the anchor's main venue, and the
// share of members whose own main venue is the community's most common one), and the member names when the community is
// small.  Usage: case_dblp_names <g.cx> <authors.tsv> <author_venues.tsv> <anchor name> [<anchor name> ...]  -> JSON.
#include "../../src-r1index/chain_index.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <cmath>
#include <tuple>
#include <unordered_map>
using namespace chainindex;
using Index = ChainIndex<double>; using Clock = std::chrono::steady_clock;

static std::string json(const std::string& s) { std::string o; for (char c : s) { if (c == '"' || c == '\\') o += '\\'; if (static_cast<unsigned char>(c) < 32) continue; o += c; } return o; }

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: case_dblp_names <g.cx> <authors.tsv> <author_venues.tsv> <anchor> ...\n"; return 2; }
    const Index ix = Index::load(argv[1]); const uint32_t n = ix.n;
    std::vector<uint32_t> perm(n), inv(n); { std::ifstream pf(std::string(argv[1]) + ".perm", std::ios::binary); pf.read(reinterpret_cast<char*>(perm.data()), static_cast<std::streamsize>(n) * 4); if (!pf) { std::cerr << "no perm\n"; return 2; } }
    for (uint32_t v = 0; v < n; ++v) inv[perm[v]] = v;
    std::vector<std::string> name(n); std::unordered_map<std::string, uint32_t> id_of; id_of.reserve(n * 2);
    { std::ifstream af(argv[2]); std::string line; while (std::getline(af, line)) { const auto t = line.find('\t'); const uint32_t id = static_cast<uint32_t>(std::stoul(line.substr(0, t))); if (id < n) { name[id] = line.substr(t + 1); id_of[name[id]] = id; } } }
    std::unordered_map<std::string, uint32_t> venue_id; std::vector<std::string> venue_name; std::vector<std::vector<std::pair<uint32_t, uint32_t>>> venues(n);   // per author: (venue, papers)
    { std::ifstream vf(argv[3]); std::string line; while (std::getline(vf, line)) { const auto t = line.find('\t'); if (t == std::string::npos) continue; auto it = id_of.find(line.substr(0, t)); if (it == id_of.end()) continue;
        size_t p = t + 1; while (p < line.size()) { size_t q = line.find('|', p); if (q == std::string::npos) q = line.size(); const std::string item = line.substr(p, q - p); const auto c = item.rfind(':');
            if (c != std::string::npos) { const std::string vn = item.substr(0, c); uint32_t vid; auto vt = venue_id.find(vn); if (vt == venue_id.end()) { vid = static_cast<uint32_t>(venue_name.size()); venue_id[vn] = vid; venue_name.push_back(vn); } else vid = vt->second;
                venues[it->second].emplace_back(vid, static_cast<uint32_t>(std::stoul(item.substr(c + 1)))); }
            p = q + 1; } } }
    uint32_t corr = UINT32_MAX; { auto it = venue_id.find("CoRR"); if (it != venue_id.end()) corr = it->second; }   // arXiv listings are not a venue
    auto main_venue = [&](uint32_t v) { uint32_t best = UINT32_MAX, cnt = 0; for (auto [vid, k] : venues[v]) if (vid != corr && k > cnt) { cnt = k; best = vid; } return best; };
    auto has_venue = [&](uint32_t v, uint32_t vid) { for (auto [w, k] : venues[v]) if (w == vid) return true; return false; };
    std::vector<uint32_t> ids(static_cast<size_t>(n) + Index::kSlack); std::vector<uint32_t> mv_count;
    const int sizes[] = {2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50};
    if (std::string(argv[4]) == "--scan") {   // sizes of the own-level community of a uniform sample of authors at s = 2 .. 10: how common is the drill-down?
        std::mt19937_64 rng(20260921); const int S[] = {2, 3, 4, 5, 6, 8, 10}; const uint64_t sample = 200000; uint64_t taken = 0;
        std::map<int, std::map<int, uint64_t>> hist;   // s -> log10 bucket of the size -> count
        uint64_t narrow = 0, wide2 = 0;                // authors active at s = 8 whose community goes from > 1000 at s = 2 to <= 50 at s = 8
        std::vector<std::tuple<uint64_t, uint64_t, uint32_t>> examples;   // (size at 8, size at 2, vertex)
        while (taken < sample) { const uint32_t vfile = static_cast<uint32_t>(rng() % n); const uint32_t vv = perm[vfile]; const uint32_t cc = ix.chain_of(vv); if (ix.omega[cc] < 8) continue; ++taken;
            uint64_t s2 = 0, s8 = 0;
            for (int s : S) { if (s > ix.omega[cc]) break; Index::Runs r; uint32_t nd; ix.community_runs(vv, s, ix.value(vv, s), r, nd); const uint64_t sz = Index::total(r); int b = 0; while (b < 7 && sz >= static_cast<uint64_t>(std::pow(10.0, b + 1))) ++b; ++hist[s][b]; if (s == 2) s2 = sz; if (s == 8) s8 = sz; }
            if (s2 > 1000) { ++wide2; if (s8 <= 50) { ++narrow; if (examples.size() < 4000) examples.emplace_back(s8, s2, vfile); } } }
        std::cout << "{\"sampled_active_at_8\":" << taken << ",\"wide_at_2\":" << wide2 << ",\"narrow_at_8\":" << narrow << ",\"hist\":{";
        bool f1 = true; for (auto& [s, h] : hist) { std::cout << (f1 ? "" : ",") << "\"" << s << "\":{"; bool f2 = true; for (auto& [b, c] : h) { std::cout << (f2 ? "" : ",") << "\"1e" << b << "\":" << c; f2 = false; } std::cout << "}"; f1 = false; }
        std::cout << "},\"examples\":["; std::sort(examples.begin(), examples.end()); size_t shown = 0;
        for (auto& [s8, s2, vfile] : examples) { if (shown >= 60) break; const uint32_t mvA = main_venue(vfile); if (mvA == UINT32_MAX) continue;
            std::cout << (shown ? "," : "") << "{\"name\":\"" << json(name[vfile]) << "\",\"main_venue\":\"" << json(venue_name[mvA]) << "\",\"size2\":" << s2 << ",\"size8\":" << s8 << "}"; ++shown; }
        std::cout << "]}\n"; return 0; }
    std::cout << "{\"n\":" << n << ",\"s_max\":" << ix.max_size << ",\"anchors\":[";
    for (int a = 4; a < argc; ++a) {
        auto it = id_of.find(argv[a]); if (it == id_of.end()) { std::cerr << "unknown author " << argv[a] << "\n"; continue; }
        const uint32_t vf = it->second, v = perm[vf]; const uint32_t c = ix.chain_of(v); const int om = ix.omega[c]; const uint32_t mvA = main_venue(vf);
        std::cout << (a > 4 ? "," : "") << "{\"name\":\"" << json(argv[a]) << "\",\"omega\":" << om << ",\"main_venue\":\"" << (mvA == UINT32_MAX ? "" : json(venue_name[mvA])) << "\",\"levels\":[";
        bool first = true;
        for (int s : sizes) { if (s > om) break;
            const double k = ix.value(v, s); Index::Runs r; uint32_t node = 0; std::array<double, 5> ts{};
            for (int pass = 0; pass < 6; ++pass) { const auto t0 = Clock::now(); ix.community_runs(v, s, k, r, node); Index::expand(r, ids.data()); const double el = std::chrono::duration<double, std::nano>(Clock::now() - t0).count(); if (pass) ts[pass - 1] = el; }
            std::sort(ts.begin(), ts.end()); const uint64_t sz = Index::total(r);
            uint64_t with_anchor_venue = 0; std::map<uint32_t, uint32_t> mv; 
            for (uint64_t i = 0; i < sz; ++i) { const uint32_t u = inv[ids[i]]; if (mvA != UINT32_MAX && has_venue(u, mvA)) ++with_anchor_venue; const uint32_t m = main_venue(u); if (m != UINT32_MAX) ++mv[m]; }
            uint32_t top_v = UINT32_MAX, top_c = 0; for (auto [vid, cnt] : mv) if (cnt > top_c) { top_c = cnt; top_v = vid; }
            std::cout << (first ? "" : ",") << "{\"s\":" << s << ",\"k\":" << k << ",\"size\":" << sz << ",\"ranges\":" << r.count() << ",\"query_ns\":" << ts[2]
                      << ",\"share_with_anchor_venue\":" << (sz ? double(with_anchor_venue) / sz : 0) << ",\"top_main_venue\":\"" << (top_v == UINT32_MAX ? "" : json(venue_name[top_v])) << "\",\"top_main_venue_share\":" << (sz ? double(top_c) / sz : 0);
            if (sz <= 40) { std::cout << ",\"members\":["; for (uint64_t i = 0; i < sz; ++i) std::cout << (i ? "," : "") << "\"" << json(name[inv[ids[i]]]) << "\""; std::cout << "]"; }
            std::cout << "}"; first = false; }
        std::cout << "]}";
    }
    std::cout << "]}\n";
    return 0;
}
