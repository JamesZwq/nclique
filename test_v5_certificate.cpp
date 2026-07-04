// v5 zero-change certificate helper (decisive experiment; docs v5 design).
//
// Loads base graph G-e, inserts (u,v), and evaluates the certificate:
//   for every seed x in {u,v} ∪ (N(u)∩N(v)):
//     rd'(x) = # s-cliques K ∋ x in G' whose OTHER members all have
//              base pop-rank > rank(x)      (capped at core(x)+1)
//     certificate requires rd'(x) <= core(x).
// Prints: CERT fires=<0/1> nseeds=<k> maxover=<max(rd'-core) clamped> cert_us=<t>
//
// Usage: v5cert <base.edges> <s> <core_rank.tsv> <u> <v>
//   core_rank.tsv: "<orig_id>\t<core>\t<pop_rank>" per line ('#' comments ok).
//
// Uncounted-vertex convention: a base vertex with core -1/absent has rank
// +inf-equivalent handling: it was never peeled (isolated wrt s-cliques), so
// it cannot be an "earlier" member of any surviving clique; give it rank = a
// value larger than all real ranks so it never blocks a seed's suffix.

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

static int S;
static std::vector<std::vector<uint32_t>> adj;
static std::vector<double> coreB;
static std::vector<int64_t> rankB;   // base pop-rank; large sentinel if absent
static std::vector<std::vector<double>> nCrT;

static double nCr(size_t n, int k) {
    if (k < 0 || (size_t)k > n) return 0.0;
    while (nCrT.size() <= n) {
        size_t r = nCrT.size();
        std::vector<double> row(S + 1, 0.0);
        row[0] = 1.0;
        for (int j = 1; j <= S && (size_t)j <= r; ++j)
            row[j] = nCrT[r - 1][j - 1] + (j <= S ? nCrT[r - 1][j] : 0.0);
        nCrT.push_back(std::move(row));
    }
    return nCrT[n][k];
}

static double edgesWithinCapped(const std::vector<uint32_t> &P, double cap) {
    if (P.size() < 2 || cap <= 0.0) return 0.0;
    double e = 0.0;
    for (uint32_t w : P) {
        const auto &nw = adj[w];
        size_t i = 0, j = 0;
        while (i < P.size() && j < nw.size()) {
            if (P[i] < nw[j]) ++i;
            else if (P[i] > nw[j]) ++j;
            else { if (nw[j] > w) { e += 1.0; if (e >= cap) return e; } ++i; ++j; }
        }
    }
    return e;
}

// capped s-clique count with `held` required members, candidate list P, piv
// accumulated pivots; returns v with (v>=cap) iff (true>=cap), exact below.
static double sctCap(std::vector<uint32_t> &P, int held, int piv, double cap) {
    if (held == S) return 1.0;
    if (held + piv + (int)P.size() < S) return 0.0;
    if (P.empty()) return nCr(piv, S - held);
    if (held == S - 1) return (double)piv + (double)P.size();
    if (held == S - 2) {
        double base = nCr(piv, 2) + (double)piv * (double)P.size();
        if (base >= cap) return base;
        return base + edgesWithinCapped(P, cap - base);
    }
    uint32_t p = P[0]; size_t best = 0;
    for (uint32_t c : P) {
        size_t cnt = 0; const auto &nc = adj[c];
        size_t i = 0, j = 0;
        while (i < P.size() && j < nc.size()) {
            if (P[i] < nc[j]) ++i;
            else if (P[i] > nc[j]) ++j;
            else { ++cnt; ++i; ++j; }
        }
        if (cnt > best || (cnt == best && c == P[0])) { best = cnt; p = c; }
    }
    auto inter = [](const std::vector<uint32_t> &a, const std::vector<uint32_t> &b,
                    std::vector<uint32_t> &o) {
        o.clear(); size_t i = 0, j = 0;
        while (i < a.size() && j < b.size()) {
            if (a[i] < b[j]) ++i;
            else if (a[i] > b[j]) ++j;
            else { o.push_back(a[i]); ++i; ++j; }
        }
    };
    double res = 0.0; std::vector<uint32_t> sub;
    inter(P, adj[p], sub);
    res += sctCap(sub, held, piv + 1, cap);
    if (res >= cap) return res;
    std::vector<uint32_t> pool = P;
    pool.erase(std::find(pool.begin(), pool.end(), p));
    std::vector<uint32_t> branch; const auto &np = adj[p];
    for (uint32_t v : pool)
        if (!std::binary_search(np.begin(), np.end(), v)) branch.push_back(v);
    for (uint32_t v : branch) {
        std::vector<uint32_t> pv; inter(pool, adj[v], pv);
        res += sctCap(pv, held + 1, piv, cap);
        if (res >= cap) return res;
        pool.erase(std::find(pool.begin(), pool.end(), v));
    }
    return res;
}

// rd'(x): s-cliques through x whose other members are all ranked AFTER x.
// DYN_PESS=1: pessimistic tie-break — a SAME-core member counts as "after"
// regardless of its (arbitrary intra-batch) rank, since some valid peel order
// could place it before x. This over-counts rd', making the certificate
// SOUND against any intra-batch ordering. Otherwise strict-rank (naive).
static int PESS = -1;
static double rdPrime(uint32_t x, double cap) {
    if (PESS < 0) PESS = std::getenv("DYN_PESS") ? 1 : 0;
    std::vector<uint32_t> P;
    P.reserve(adj[x].size());
    const int64_t rx = rankB[x];
    const double cx = coreB[x];
    for (uint32_t y : adj[x]) {
        bool after;
        if (PESS) after = (coreB[y] > cx) || (coreB[y] == cx);  // same-core = after
        else after = rankB[y] > rx;
        if (after) P.push_back(y);
    }
    return sctCap(P, 1, 0, cap);
}

int main(int argc, char **argv) {
    if (argc < 6) {
        std::fprintf(stderr, "usage: %s <base.edges> <s> <core_rank.tsv> <u> <v>\n", argv[0]);
        return 1;
    }
    S = std::atoi(argv[2]);
    const uint32_t U = (uint32_t)std::atoll(argv[4]);
    const uint32_t V = (uint32_t)std::atoll(argv[5]);
    FILE *f = std::fopen(argv[1], "r");
    if (!f) { std::perror("graph"); return 1; }
    size_t n, m;
    if (std::fscanf(f, "%zu %zu", &n, &m) != 2) return 1;
    adj.assign(n, {});
    for (size_t i = 0; i < m; ++i) {
        uint32_t a, b;
        if (std::fscanf(f, "%u %u", &a, &b) != 2) break;
        if (a == b) continue;
        adj[a].push_back(b); adj[b].push_back(a);
    }
    std::fclose(f);
    for (auto &l : adj) {
        std::sort(l.begin(), l.end());
        l.erase(std::unique(l.begin(), l.end()), l.end());
    }
    coreB.assign(n, 0.0);
    const int64_t BIG = (int64_t)4e18;
    rankB.assign(n, BIG);   // absent/uncounted => never blocks a suffix
    f = std::fopen(argv[3], "r");
    if (!f) { std::perror("corerank"); return 1; }
    char line[256];
    while (std::fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        unsigned long long id; double c; long long rk;
        if (std::sscanf(line, "%llu %lf %lld", &id, &c, &rk) == 3 && id < n) {
            coreB[id] = c; rankB[id] = rk;
        }
    }
    std::fclose(f);
    nCrT.reserve(256); nCr(0, 0);

    auto t0 = std::chrono::steady_clock::now();
    // seeds = {u,v} ∪ (N(u)∩N(v)), computed in the BASE graph, THEN insert e.
    std::vector<uint32_t> W;
    {
        const auto &nu = adj[U], &nv = adj[V];
        size_t i = 0, j = 0;
        while (i < nu.size() && j < nv.size()) {
            if (nu[i] < nv[j]) ++i;
            else if (nu[i] > nv[j]) ++j;
            else { W.push_back(nu[i]); ++i; ++j; }
        }
    }
    adj[U].insert(std::lower_bound(adj[U].begin(), adj[U].end(), V), V);
    adj[V].insert(std::lower_bound(adj[V].begin(), adj[V].end(), U), U);

    std::vector<uint32_t> seeds = {U, V};
    for (uint32_t w : W) seeds.push_back(w);

    bool fires = true;
    int64_t maxover = -(int64_t)4e18;
    for (uint32_t x : seeds) {
        const double cap = coreB[x] + 1.0;   // only the <= core(x) decision matters
        double rd = rdPrime(x, cap);
        int64_t over = (int64_t)rd - (int64_t)coreB[x];   // >0 => this seed fails
        if (over > maxover) maxover = over;
        if (rd > coreB[x]) fires = false;
    }
    auto t1 = std::chrono::steady_clock::now();
    double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
    std::printf("CERT fires=%d nseeds=%zu maxover=%lld cert_us=%.0f\n",
                fires ? 1 : 0, seeds.size(), (long long)maxover, us);
    return 0;
}
