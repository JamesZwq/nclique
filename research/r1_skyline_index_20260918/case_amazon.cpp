// Case study on the Amazon co-purchase network (SNAP com-amazon) with product titles and categories (SNAP amazon-meta).
// A product's own-level community at size s is the (s, kappa_s(q))-nucleus containing it.  Purity of a community for a
// query q: the share of its members that share a leaf category with q (fine), and the share that share a subject with
// q (the third level of the category tree, e.g. "Religion & Spirituality", "Jazz", "Comedy").
// Usage: case_amazon <g.cx> <g.map> <amazon-meta.tsv> --scan            -> sample of products, purity and size per s
//        case_amazon <g.cx> <g.map> <amazon-meta.tsv> --query <Id> ...   -> drill-down of the given products (metadata Ids)
//        case_amazon <g.cx> <g.map> <amazon-meta.tsv> --find <substring> -> list products whose title contains it
#include "../../src-r1index/chain_index.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <unordered_map>
using namespace chainindex;
using Index = ChainIndex<double>; using Clock = std::chrono::steady_clock;

static std::string json(const std::string& s) { std::string o; for (char c : s) { if (c == '"' || c == '\\') o += '\\'; if (static_cast<unsigned char>(c) < 32) continue; o += c; } return o; }
struct Product { std::string title, group; std::vector<uint32_t> leaf, subject; };

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: see header\n"; return 2; }
    const Index ix = Index::load(argv[1]); const uint32_t n = ix.n;
    std::vector<uint32_t> perm(n), inv(n); { std::ifstream pf(std::string(argv[1]) + ".perm", std::ios::binary); pf.read(reinterpret_cast<char*>(perm.data()), static_cast<std::streamsize>(n) * 4); if (!pf) { std::cerr << "no perm\n"; return 2; } }
    for (uint32_t v = 0; v < n; ++v) inv[perm[v]] = v;
    std::vector<uint32_t> orig(n); std::unordered_map<uint32_t, uint32_t> label_of; { std::ifstream mf(argv[2]); for (uint32_t v = 0; v < n; ++v) { mf >> orig[v]; label_of[orig[v]] = v; } }
    std::vector<Product> prod(n); std::unordered_map<std::string, uint32_t> cat_id; std::vector<std::string> cat_name;
    auto intern = [&](const std::string& s) { auto it = cat_id.find(s); if (it != cat_id.end()) return it->second; const uint32_t id = static_cast<uint32_t>(cat_name.size()); cat_id[s] = id; cat_name.push_back(s); return id; };
    { std::ifstream mf(argv[3]); std::string line; uint32_t got = 0;
      while (std::getline(mf, line)) { std::vector<std::string> f; size_t p = 0; while (true) { size_t q = line.find('\t', p); f.push_back(line.substr(p, q == std::string::npos ? std::string::npos : q - p)); if (q == std::string::npos) break; p = q + 1; }
          if (f.size() < 5) continue; auto it = label_of.find(static_cast<uint32_t>(std::stoul(f[0]))); if (it == label_of.end()) continue; Product& pr = prod[it->second]; pr.group = f[1]; pr.title = f[3]; ++got;
          size_t a = 0; while (a < f[4].size()) { size_t b = f[4].find(';', a); if (b == std::string::npos) b = f[4].size(); const std::string path = f[4].substr(a, b - a);
              std::vector<std::string> comp; size_t c = 0; while (true) { size_t d = path.find('|', c); comp.push_back(path.substr(c, d == std::string::npos ? std::string::npos : d - c)); if (d == std::string::npos) break; c = d + 1; }
              if (comp.size() >= 2 && !comp.back().empty()) pr.leaf.push_back(intern(comp.back()));
              if (comp.size() >= 4) pr.subject.push_back(intern(comp[3])); else if (comp.size() >= 2) pr.subject.push_back(intern(comp.back()));
              a = b + 1; } }
      std::cerr << "products with metadata: " << got << " of " << n << "\n"; }
    std::vector<uint32_t> mark(cat_name.size() + 1, 0); uint32_t stamp = 0; std::vector<uint32_t> ids(static_cast<size_t>(n) + Index::kSlack);
    // purity of the member list ids[0, sz) for query label q (file labels after inv): linear in the members, array counters
    std::vector<uint32_t> qmark(cat_name.size() + 1, 0), cnt(cat_name.size() + 1, 0), touched; uint32_t qstamp = 0;
    auto purity = [&](uint32_t q, uint64_t sz, double& leaf_share, double& subject_share, std::string& top_subject, double& top_share) {
        ++qstamp; for (uint32_t c : prod[q].leaf) qmark[c] = qstamp; uint32_t ls = 0;
        for (uint64_t i = 0; i < sz; ++i) { const uint32_t u = inv[ids[i]]; for (uint32_t c : prod[u].leaf) if (qmark[c] == qstamp) { ++ls; break; } }
        ++qstamp; for (uint32_t c : prod[q].subject) qmark[c] = qstamp; uint32_t ss = 0; touched.clear();
        for (uint64_t i = 0; i < sz; ++i) { const uint32_t u = inv[ids[i]]; bool hit = false;
            for (uint32_t c : prod[u].subject) { if (qmark[c] == qstamp) hit = true; if (cnt[c] == 0) touched.push_back(c); if (cnt[c] < (1u << 31)) ++cnt[c]; }
            ss += hit; }
        uint32_t best = 0, bc = 0; for (uint32_t c : touched) { if (cnt[c] > bc) { bc = cnt[c]; best = c; } }
        for (uint32_t c : touched) cnt[c] = 0;   // a member with the same subject on two paths counts twice; rare and harmless for the top subject
        leaf_share = sz ? double(ls) / sz : 0; subject_share = sz ? double(ss) / sz : 0; top_subject = bc ? cat_name[best] : ""; top_share = sz ? std::min(1.0, double(bc) / sz) : 0; };
    const std::string mode = argv[4];
    if (mode == "--find") { for (uint32_t v = 0; v < n; ++v) if (prod[v].title.find(argv[5]) != std::string::npos) std::cout << orig[v] << "\t" << prod[v].group << "\t" << ix.omega[ix.chain_of(perm[v])] << "\t" << prod[v].title << "\n"; return 0; }
    if (mode == "--scan") {
        std::mt19937_64 rng(20260921); const int S[] = {2, 3, 4, 5, 6, 7}; const uint64_t sample = 5000; uint64_t taken = 0;
        std::map<int, std::array<double, 4>> acc; std::map<int, uint64_t> cnt, pure; std::map<int, std::map<int, uint64_t>> hist; std::map<int, uint64_t> best_s; std::map<int, std::vector<uint64_t>> sizes_at;
        while (taken < sample) { const uint32_t vf = static_cast<uint32_t>(rng() % n); const uint32_t v = perm[vf]; const uint32_t c = ix.chain_of(v); if (ix.omega[c] < 3 || prod[vf].leaf.empty()) continue; ++taken;
            double bestp = -1; int bs = 2;
            for (int s : S) { if (s > ix.omega[c]) break; Index::Runs r; uint32_t nd; ix.community_runs(v, s, ix.value(v, s), r, nd); const uint64_t sz = static_cast<uint64_t>(Index::expand(r, ids.data()) - ids.data());
                double ls, ss, ts; std::string tn; purity(vf, sz, ls, ss, tn, ts); auto& a = acc[s]; a[0] += sz; a[1] += ls; a[2] += ss; a[3] += ts; ++cnt[s]; pure[s] += ls >= 0.8; sizes_at[s].push_back(sz);
                int b = 0; while (b < 6 && sz >= static_cast<uint64_t>(std::pow(10.0, b + 1))) ++b; ++hist[s][b]; if (ls > bestp + 1e-12) { bestp = ls; bs = s; } }
            ++best_s[bs]; }
        std::cout << "{\"sampled\":" << taken << ",\"per_size\":{"; bool f = true;
        for (auto& [s, a] : acc) { auto& sv = sizes_at[s]; std::sort(sv.begin(), sv.end());
            std::cout << (f ? "" : ",") << "\"" << s << "\":{\"queries\":" << cnt[s] << ",\"median_size\":" << sv[sv.size() / 2] << ",\"mean_size\":" << a[0] / cnt[s] << ",\"mean_leaf_share\":" << a[1] / cnt[s] << ",\"mean_subject_share\":" << a[2] / cnt[s] << ",\"mean_top_subject_share\":" << a[3] / cnt[s] << ",\"pure_leaf_08\":" << double(pure[s]) / cnt[s] << ",\"size_hist\":{";
            bool g = true; for (auto& [b, k] : hist[s]) { std::cout << (g ? "" : ",") << "\"1e" << b << "\":" << k; g = false; } std::cout << "}}"; f = false; }
        std::cout << "},\"best_size_by_leaf_purity\":{"; f = true; for (auto& [s, k] : best_s) { std::cout << (f ? "" : ",") << "\"" << s << "\":" << k; f = false; } std::cout << "}}\n"; return 0; }
    if (mode == "--query") {
        std::cout << "{\"queries\":["; bool firstq = true;
        for (int a = 5; a < argc; ++a) { auto it = label_of.find(static_cast<uint32_t>(std::stoul(argv[a]))); if (it == label_of.end()) { std::cerr << "unknown id " << argv[a] << "\n"; continue; }
            const uint32_t vf = it->second, v = perm[vf]; const uint32_t c = ix.chain_of(v); const int om = ix.omega[c];
            std::cout << (firstq ? "" : ",") << "{\"id\":" << orig[vf] << ",\"title\":\"" << json(prod[vf].title) << "\",\"group\":\"" << prod[vf].group << "\",\"omega\":" << om << ",\"subjects\":["; bool fs = true; for (uint32_t s : prod[vf].subject) { std::cout << (fs ? "" : ",") << "\"" << json(cat_name[s]) << "\""; fs = false; } std::cout << "],\"levels\":["; firstq = false;
            for (int s = 2; s <= om; ++s) { const double k = ix.value(v, s); Index::Runs r; uint32_t nd; std::array<double, 5> ts{}; uint64_t sz = 0;
                for (int pass = 0; pass < 6; ++pass) { const auto t0 = Clock::now(); ix.community_runs(v, s, k, r, nd); sz = static_cast<uint64_t>(Index::expand(r, ids.data()) - ids.data()); const double el = std::chrono::duration<double, std::nano>(Clock::now() - t0).count(); if (pass) ts[pass - 1] = el; }
                std::sort(ts.begin(), ts.end()); double ls, ss, tsh; std::string tn; purity(vf, sz, ls, ss, tn, tsh);
                std::cout << (s > 2 ? "," : "") << "{\"s\":" << s << ",\"k\":" << k << ",\"size\":" << sz << ",\"ranges\":" << r.count() << ",\"query_ns\":" << ts[2] << ",\"leaf_share\":" << ls << ",\"subject_share\":" << ss << ",\"top_subject\":\"" << json(tn) << "\",\"top_subject_share\":" << tsh;
                if (sz <= 30) { std::cout << ",\"members\":["; for (uint64_t i = 0; i < sz; ++i) { const uint32_t u = inv[ids[i]]; std::cout << (i ? "," : "") << "\"" << json(prod[u].title) << " (" << prod[u].group << ")\""; } std::cout << "]"; }
                std::cout << "}"; }
            std::cout << "]}"; }
        std::cout << "]}\n"; return 0; }
    std::cerr << "unknown mode\n"; return 2;
}
