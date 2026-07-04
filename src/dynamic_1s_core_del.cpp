// Dynamic (1,s)-core maintenance — DELETION core-update (symmetric to v6).
//
// Delete e=(u,v): removes exactly the s-cliques containing BOTH u,v. Cores only
// FALL. Unlike insertion, deletion needs NO optimism/ring handling: a vertex
// falls as soon as its ACTUAL support drops below its core, a monotone cascade.
// Grow the region from CONFIRMED FALLERS; the pinned peel (boundary frozen at
// old core = new core for non-fallers) yields the fallen cores exactly.
//
//   Phase 1: W = N(u)∩N(v); if G[W] has no (s-2)-clique the deleted edge is in
//            no s-clique -> no change, exit. Then REMOVE e.
//   Phase 2: region = seeds; repeat { pinned-peel region; grow by
//            clq-nbr(fallers) } until no new faller. Cap rounds -> FALLBACK.
//   Emit CHANGED (core decreased) from the final peel.
//
// Correctness prototype (fresh index per edge): exact keys via full recompute
// per removal; region tiny so cheap. Leaf surgery (§23) is a separate
// index-maintenance concern (cost measured by E2), not needed here.
//
// Usage: dynamic_1s_core_del <graph.edges> <s> <core.tsv> <u> <v>
//   graph.edges must CONTAIN edge (u,v); core.tsv = cores of that graph.

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

static int S;
static std::vector<std::vector<uint32_t>> adj;
static std::vector<double> coreBase;
static std::vector<uint8_t> alive;
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

// ==================== bitset SCT engine (copied verbatim from v4) ==========
static size_t Unb = 0, Uw = 0;
static std::vector<uint64_t> rowsBuf;
static std::vector<uint32_t> lid2g, g2l, g2lStamp;
static uint32_t stampCur = 0;
static inline uint64_t *rowOf(size_t lid) { return rowsBuf.data() + lid * Uw; }
static inline void bitSet(uint64_t *b, size_t i) { b[i >> 6] |= 1ull << (i & 63); }
static inline void bitClr(uint64_t *b, size_t i) { b[i >> 6] &= ~(1ull << (i & 63)); }
static inline size_t pcAll(const uint64_t *a) {
    size_t s = 0; for (size_t w = 0; w < Uw; ++w) s += (size_t)std::popcount(a[w]);
    return s;
}
static inline size_t pcAnd(const uint64_t *a, const uint64_t *b) {
    size_t s = 0; for (size_t w = 0; w < Uw; ++w) s += (size_t)std::popcount(a[w] & b[w]);
    return s;
}
static void buildUniverse(const uint32_t *verts, size_t cnt) {
    ++stampCur; Unb = cnt; Uw = (cnt + 63) / 64; if (Uw == 0) Uw = 1;
    lid2g.assign(verts, verts + cnt);
    for (size_t i = 0; i < cnt; ++i) { g2l[verts[i]] = (uint32_t)i; g2lStamp[verts[i]] = stampCur; }
    rowsBuf.assign(Unb * Uw, 0ull);
    for (size_t i = 0; i < cnt; ++i) {
        uint64_t *r = rowOf(i);
        for (uint32_t g : adj[verts[i]]) if (g2lStamp[g] == stampCur) bitSet(r, g2l[g]);
    }
}
static std::vector<std::vector<uint64_t>> scrT, scrPool, scrD;
static inline uint64_t *getScr(std::vector<std::vector<uint64_t>> &s, int depth) {
    if ((int)s.size() <= depth) s.resize(depth + 1);
    if (s[depth].size() < Uw) s[depth].assign(Uw, 0ull);
    return s[depth].data();
}
static double bitSct(const uint64_t *P, int held, int piv, double cap, int depth) {
    if (held == S) return 1.0;
    size_t cnt = pcAll(P);
    if (held + piv + (int)cnt < S) return 0.0;
    if (cnt == 0) return nCr(piv, S - held);
    if (held == S - 1) return (double)piv + (double)cnt;
    if (held == S - 2) {
        double base = nCr(piv, 2) + (double)piv * (double)cnt;
        if (base >= cap) return base;
        double need2 = 2.0 * (cap - base), sum2 = 0.0;
        for (size_t w = 0; w < Uw && sum2 < need2; ++w) {
            uint64_t m = P[w];
            while (m) {
                size_t lid = (w << 6) + (size_t)std::countr_zero(m); m &= m - 1;
                sum2 += (double)pcAnd(P, rowOf(lid));
                if (sum2 >= need2) break;
            }
        }
        return base + sum2 * 0.5;
    }
    size_t p = 0, best = 0; bool first = true;
    for (size_t w = 0; w < Uw; ++w) {
        uint64_t m = P[w];
        while (m) {
            size_t lid = (w << 6) + (size_t)std::countr_zero(m); m &= m - 1;
            size_t c = pcAnd(P, rowOf(lid));
            if (first || c > best) { best = c; p = lid; first = false; }
        }
    }
    uint64_t *T = getScr(scrT, depth);
    { const uint64_t *rp = rowOf(p); for (size_t w = 0; w < Uw; ++w) T[w] = P[w] & rp[w]; }
    double res = bitSct(T, held, piv + 1, cap, depth + 1);
    if (res >= cap) return res;
    uint64_t *pool = getScr(scrPool, depth); uint64_t *D = getScr(scrD, depth);
    { const uint64_t *rp = rowOf(p);
      for (size_t w = 0; w < Uw; ++w) { pool[w] = P[w]; D[w] = P[w] & ~rp[w]; }
      bitClr(pool, p); bitClr(D, p); }
    for (size_t w = 0; w < Uw; ++w) {
        uint64_t m = D[w];
        while (m) {
            size_t v = (w << 6) + (size_t)std::countr_zero(m); m &= m - 1;
            const uint64_t *rv = rowOf(v);
            for (size_t x = 0; x < Uw; ++x) T[x] = pool[x] & rv[x];
            res += bitSct(T, held + 1, piv, cap, depth + 1);
            if (res >= cap) return res;
            bitClr(pool, v);
        }
    }
    return res;
}
static double bitCountList(const std::vector<uint32_t> &P0, int held, double cap) {
    if (held + (int)P0.size() < S) return 0.0;
    buildUniverse(P0.data(), P0.size());
    static std::vector<uint64_t> full;
    full.assign(Uw, 0ull);
    for (size_t w = 0; w + 1 < Uw; ++w) full[w] = ~0ull;
    full[Uw - 1] = (Unb & 63) ? ((1ull << (Unb & 63)) - 1) : ~0ull;
    if (Unb == 0) full[0] = 0ull;
    return bitSct(full.data(), held, 0, cap, 0);
}
static const double INF = std::numeric_limits<double>::infinity();

// ==================== v6 ====================

// clique-neighbours of x in G' (current adj): y in N(x) sharing >= s-2 common
// neighbours with x (necessary to co-inhabit an s-clique). Over-approximation
// of "shares an s-clique with x" (sound for completeness, tight in practice).
static void cliqueNbrs(uint32_t x, std::vector<uint32_t> &out) {
    const int k = S - 2;
    const auto &ax = adj[x];
    for (uint32_t y : ax) {
        if (k <= 0) { out.push_back(y); continue; }
        const auto &ay = adj[y];
        size_t i = 0, j = 0; int c = 0; bool ok = false;
        while (i < ax.size() && j < ay.size()) {
            if (ax[i] < ay[j]) ++i;
            else if (ax[i] > ay[j]) ++j;
            else { if (++c >= k) { ok = true; break; } ++i; ++j; }
        }
        if (ok) out.push_back(y);
    }
}

// s-cliques through x within the current alive set (exact key).
static double supportOf(uint32_t x) {
    std::vector<uint32_t> P;
    P.reserve(adj[x].size());
    for (uint32_t y : adj[x]) if (alive[y]) P.push_back(y);
    return bitCountList(P, 1, INF);
}

int main(int argc, char **argv) {
    if (argc < 6) {
        std::fprintf(stderr, "usage: %s <base.edges> <s> <core.tsv> <u> <v>\n", argv[0]);
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
    for (auto &l : adj) { std::sort(l.begin(), l.end()); l.erase(std::unique(l.begin(), l.end()), l.end()); }
    coreBase.assign(n, 0.0);
    f = std::fopen(argv[3], "r");
    if (!f) { std::perror("cores"); return 1; }
    char line[256];
    while (std::fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        unsigned long long id; double c;
        if (std::sscanf(line, "%llu %lf", &id, &c) == 2 && id < n) coreBase[id] = c;
    }
    std::fclose(f);
    nCrT.reserve(256); nCr(0, 0);
    g2l.assign(n, 0); g2lStamp.assign(n, 0);
    alive.assign(n, 0);
    std::vector<uint8_t> inRegion(n, 0);
    const bool dbg = std::getenv("DYN1S_DEBUG") != nullptr;
    size_t CAP = std::max<size_t>(4096, n / 16);
    if (const char *ce = std::getenv("DYN1S_CAP")) CAP = (size_t)std::atoll(ce);

    auto t0 = std::chrono::steady_clock::now();

    // -------- Phase 1: W, new-clique check, insert e --------
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
    double newCliques = ((int)W.size() + 2 >= S) ? bitCountList(W, 2, 1.0) : 0.0;
    if (newCliques == 0.0) {
        double us = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - t0).count();
        std::printf("STATS region=0 pinned=0 rounds=0 pops=0 changed=0 "
                    "newcliques=0 insert_us=%.0f\n", us);
        return 0;
    }
    // REMOVE e (W was computed with e present; W = N(u)∩N(v) is unaffected by
    // e itself). All Phase-2 counting is now in G-e.
    adj[U].erase(std::lower_bound(adj[U].begin(), adj[U].end(), V));
    adj[V].erase(std::lower_bound(adj[V].begin(), adj[V].end(), U));

    // -------- Phase 2: region = seeds; grow from FALLERS (no ring/optimism) ---
    std::vector<uint32_t> region;
    auto addRegion = [&](uint32_t x) { if (!inRegion[x]) { inRegion[x] = 1; region.push_back(x); } };
    addRegion(U); addRegion(V);
    for (uint32_t w : W) addRegion(w);

    std::vector<double> newCore(n, -1.0);
    int rounds = 0;
    size_t peelPops = 0, pinnedSz = 0;
    bool fallback = false;

    while (true) {
        ++rounds;
        // ---- pinned peel on `region`, boundary = N(region)\region frozen ----
        std::vector<uint32_t> pinned;
        {
            std::vector<uint8_t> inP(n, 0);
            for (uint32_t x : region)
                for (uint32_t y : adj[x])
                    if (!inRegion[y] && !inP[y]) { inP[y] = 1; pinned.push_back(y); }
        }
        std::sort(pinned.begin(), pinned.end(),
                  [](uint32_t a, uint32_t b) { return coreBase[a] < coreBase[b]; });
        pinnedSz = pinned.size();
        // DELETION: cores FALL, so a region vertex can drop below any pinned
        // level and below the region's min old core. The L0 fast-start
        // (insertion-only, where cores rise) is therefore DISABLED: include all
        // pinned (j0=0) and start minCore at 0.
        size_t j0 = 0;

        for (uint32_t x : region) alive[x] = 1;
        for (size_t k = j0; k < pinned.size(); ++k) alive[pinned[k]] = 1;

        std::vector<uint32_t> pool = region;
        std::vector<double> key(n, 0.0);
        for (uint32_t x : region) key[x] = supportOf(x);

        auto removeAndRecount = [&](uint32_t z) {
            alive[z] = 0;
            for (uint32_t y : adj[z]) if (alive[y] && inRegion[y]) key[y] = supportOf(y);
        };

        size_t jc = j0;
        double minCore = 0.0;   // deletion: cores fall; no L0 floor
        for (uint32_t x : region) newCore[x] = -1.0;
        while (!pool.empty()) {
            size_t bi = 0;
            for (size_t i = 1; i < pool.size(); ++i) if (key[pool[i]] < key[pool[bi]]) bi = i;
            double rmin = key[pool[bi]];
            double pmin = (jc < pinned.size()) ? coreBase[pinned[jc]] : INF;
            if (pmin <= rmin) {
                uint32_t z = pinned[jc++];
                if (coreBase[z] > minCore) minCore = coreBase[z];
                removeAndRecount(z); ++peelPops;
            } else {
                uint32_t z = pool[bi]; pool[bi] = pool.back(); pool.pop_back();
                if (key[z] > minCore) minCore = key[z];
                newCore[z] = minCore;
                removeAndRecount(z); ++peelPops;
            }
        }
        // reset alive for next round (region grew; recompute fresh)
        for (uint32_t x : region) alive[x] = 0;
        for (size_t k = j0; k < pinned.size(); ++k) alive[pinned[k]] = 0;

        // ---- grow from confirmed FALLERS. Deletion needs NO ring/optimism:
        // a vertex falls iff its ACTUAL support (in G-e, boundary pinned)
        // dropped below its core — a monotone-down cascade. When x falls, its
        // clique-neighbours may lose x as a high-level supporter and fall too,
        // so grow to clq-nbr(fallers). Fall-chain connectivity (dual of Cor 3a)
        // ⟹ region reaches all fallers.
        std::vector<uint32_t> fallers;
        for (uint32_t x : region) if (newCore[x] < coreBase[x]) fallers.push_back(x);
        std::vector<uint32_t> newv, tmp;
        for (uint32_t x : fallers) {
            tmp.clear(); cliqueNbrs(x, tmp);
            for (uint32_t y : tmp) if (!inRegion[y]) newv.push_back(y);
        }
        std::sort(newv.begin(), newv.end());
        newv.erase(std::unique(newv.begin(), newv.end()), newv.end());
        if (dbg)
            std::fprintf(stderr, "[dbg] round=%d region=%zu pinned=%zu fallers=%zu new=%zu\n",
                         rounds, region.size(), pinnedSz, fallers.size(), newv.size());
        if (newv.empty()) break;
        if (region.size() + newv.size() > CAP || rounds >= 40) { fallback = true; break; }
        for (uint32_t y : newv) addRegion(y);
    }

    auto t1 = std::chrono::steady_clock::now();
    double us = std::chrono::duration<double, std::micro>(t1 - t0).count();

    if (fallback) {
        std::printf("FALLBACK\n");
        std::printf("STATS region=%zu pinned=%zu rounds=%d pops=%zu changed=0 "
                    "newcliques=%.0f insert_us=%.0f\n",
                    region.size(), pinnedSz, rounds, peelPops, newCliques, us);
        return 0;
    }
    size_t nChanged = 0;
    for (uint32_t x : region) {
        if (newCore[x] > coreBase[x])   // deletion: cores must not RISE
            std::fprintf(stderr, "WARN rise vertex=%u old=%.0f new=%.0f\n", x, coreBase[x], newCore[x]);
        if (newCore[x] != coreBase[x]) {
            std::printf("CHANGED %u %.0f %.0f\n", x, coreBase[x], newCore[x]);
            ++nChanged;
        }
    }
    std::printf("STATS region=%zu pinned=%zu rounds=%d pops=%zu changed=%zu "
                "newcliques=%.0f insert_us=%.0f\n",
                region.size(), pinnedSz, rounds, peelPops, nChanged, newCliques, us);
    return 0;
}
