// Dynamic (1,s)-core maintenance prototype — insert-only, v2 (M2 algorithm +
// M1 bitset counting engine).
//
// v2 (see docs/dynamic_v2_spec.md) replaces v1's expansion-round loop with a
// three-phase static-discovery pipeline. Persistent state between updates:
// core[] + adjacency ONLY (no CPI/SDCT, no hidden caches).
//
// Insert (u,v):
//   Phase 1 (§7): W = N(u)∩N(v); count new (s-2)-cliques in G[W] via one SCT.
//                 If 0 -> no core changes (Cor 3c), STATS, exit (fast path).
//                 Then insert e into adjacency (G' from here on).
//   Phase 2 (§6): trigger closure from {u,v} using the STATIC admission test
//                 PASS(y) [OS(y) >= c(y)+1 with the TS-disjunct]; then, STRICTLY
//                 AFTERWARD, the eviction cascade (EOS test). Respect the |C| cap
//                 (default max(4096,n/16)) with FALLBACK output.
//   Phase 3 (§8): ONE pinned peel on the final C (τ-view, L0 fast-start, exact
//                 delta subtraction, max-rule). No rounds.
//
// M1 engine (§11): every count runs on a LOCAL COMPACT UNIVERSE — global ids
// mapped to dense local ids, adjacency rows as bitsets (ceil(|U|/64) words),
// pivot selection by popcount, branch sets by word-AND, member predicates
// (alive / region / core>=τ) folded into row masks. Counts are exact integers
// (<2^53, stored in double), so all outputs are bit-identical to the M2
// vector engine: capped counts agree on every threshold decision, uncapped
// counts (peel keys) agree exactly, hence identical pops and CHANGED lines.
//
// Usage:
//   dynamic_1s_core <base.edges> <s> <core_base.tsv> <u> <v>
// core_base.tsv: "<orig_id>\t<core>" lines ('#' comments ok); missing id = 0.
// Output: "CHANGED x old new" lines + one "STATS ..." line (or "FALLBACK" + STATS).

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

static int S;                                  // clique size s
static std::vector<std::vector<uint32_t>> adj; // sorted adjacency
static std::vector<double> coreBase;           // maintained core[]
static std::vector<uint8_t> inRegion;          // region membership (C)
static std::vector<std::vector<double>> nCrT;  // Pascal, cols 0..S

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

// ==================== M1 bitset engine (§11) ====================
// One active local universe at a time. Universe = an ordered list of global
// vertex ids; rows = adjacency restricted to the universe, as bitsets.

static size_t Unb = 0, Uw = 0;             // universe size, words per row
static std::vector<uint64_t> rowsBuf;      // Unb rows × Uw words
static std::vector<uint32_t> lid2g;        // local id -> global id
static std::vector<uint32_t> g2l;          // global -> local (valid iff stamped)
static std::vector<uint32_t> g2lStamp;
static uint32_t stampCur = 0;

static inline uint64_t *rowOf(size_t lid) { return rowsBuf.data() + lid * Uw; }
static inline void bitSet(uint64_t *b, size_t i) { b[i >> 6] |= 1ull << (i & 63); }
static inline void bitClr(uint64_t *b, size_t i) { b[i >> 6] &= ~(1ull << (i & 63)); }
static inline size_t pcAll(const uint64_t *a) {
    size_t s = 0;
    for (size_t w = 0; w < Uw; ++w) s += (size_t)std::popcount(a[w]);
    return s;
}
static inline size_t pcAnd(const uint64_t *a, const uint64_t *b) {
    size_t s = 0;
    for (size_t w = 0; w < Uw; ++w) s += (size_t)std::popcount(a[w] & b[w]);
    return s;
}
// zero bits [from, to) of b
static inline void clearRange(uint64_t *b, size_t from, size_t to) {
    if (from >= to) return;
    size_t wf = from >> 6, wt = (to - 1) >> 6;
    uint64_t mf = ~0ull << (from & 63);            // bits >= from in word wf
    uint64_t mt = (to & 63) ? ((1ull << (to & 63)) - 1) : ~0ull; // bits < to in wt
    if (wf == wt) { b[wf] &= ~(mf & mt); return; }
    b[wf] &= ~mf;
    for (size_t w = wf + 1; w < wt; ++w) b[w] = 0;
    b[wt] &= ~mt;
}

static void buildUniverse(const uint32_t *verts, size_t cnt) {
    ++stampCur;
    Unb = cnt;
    Uw = (cnt + 63) / 64;
    if (Uw == 0) Uw = 1;
    lid2g.assign(verts, verts + cnt);
    for (size_t i = 0; i < cnt; ++i) { g2l[verts[i]] = (uint32_t)i; g2lStamp[verts[i]] = stampCur; }
    rowsBuf.assign(Unb * Uw, 0ull);
    for (size_t i = 0; i < cnt; ++i) {
        uint64_t *r = rowOf(i);
        for (uint32_t g : adj[verts[i]])
            if (g2lStamp[g] == stampCur) bitSet(r, g2l[g]);
    }
}

// per-depth recursion scratch (grown on demand; depth is bounded by the pool
// shrinking by >=1 per level, in practice ~max-clique-size)
static std::vector<std::vector<uint64_t>> scrT, scrPool, scrD;
static inline uint64_t *getScr(std::vector<std::vector<uint64_t>> &s, int depth) {
    if ((int)s.size() <= depth) s.resize(depth + 1);
    if (s[depth].size() < Uw) s[depth].assign(Uw, 0ull);
    return s[depth].data();
}

// SCT recursion on the current universe. P = candidate bitset (all adjacent to
// every held and every pivot vertex, maintained by construction), held = #
// required members, piv = accumulated pivots. Capped: returns v with
// (v >= cap) IFF (true count >= cap), and v == true count when v < cap; pass
// cap = +inf for exact counts (Phase-3 keys, §16a). P is never modified.
static double bitSct(const uint64_t *P, int held, int piv, double cap, int depth) {
    if (held == S) return 1.0;
    size_t cnt = pcAll(P);
    if (held + piv + (int)cnt < S) return 0.0;
    if (cnt == 0) return nCr(piv, S - held);
    // closed form one level above the leaves
    if (held == S - 1) return (double)piv + (double)cnt;
    // two levels above the leaves: C(piv,2) + piv*|P| + E(G[P])
    if (held == S - 2) {
        double base = nCr(piv, 2) + (double)piv * (double)cnt;
        if (base >= cap) return base;
        double need2 = 2.0 * (cap - base);  // early-exit on DOUBLE-counted sum
        double sum2 = 0.0;
        for (size_t w = 0; w < Uw && sum2 < need2; ++w) {
            uint64_t m = P[w];
            while (m) {
                size_t lid = (w << 6) + (size_t)std::countr_zero(m);
                m &= m - 1;
                sum2 += (double)pcAnd(P, rowOf(lid));
                if (sum2 >= need2) break;
            }
        }
        return base + sum2 * 0.5;   // exact when below cap (sum2 even)
    }
    // pivot p = argmax |P ∩ N(p)| by popcount
    size_t p = 0, best = 0;
    bool first = true;
    for (size_t w = 0; w < Uw; ++w) {
        uint64_t m = P[w];
        while (m) {
            size_t lid = (w << 6) + (size_t)std::countr_zero(m);
            m &= m - 1;
            size_t c = pcAnd(P, rowOf(lid));
            if (first || c > best) { best = c; p = lid; first = false; }
        }
    }
    // pivot branch: P ∩ N(p)
    uint64_t *T = getScr(scrT, depth);
    {
        const uint64_t *rp = rowOf(p);
        for (size_t w = 0; w < Uw; ++w) T[w] = P[w] & rp[w];
    }
    double res = bitSct(T, held, piv + 1, cap, depth + 1);
    if (res >= cap) return res;
    // hold branches: v ∈ (P \ {p}) \ N(p), snapshot list, sequentially-shrinking pool
    uint64_t *pool = getScr(scrPool, depth);
    uint64_t *D = getScr(scrD, depth);
    {
        const uint64_t *rp = rowOf(p);
        for (size_t w = 0; w < Uw; ++w) { pool[w] = P[w]; D[w] = P[w] & ~rp[w]; }
        bitClr(pool, p);
        bitClr(D, p);
    }
    for (size_t w = 0; w < Uw; ++w) {
        uint64_t m = D[w];
        while (m) {
            size_t v = (w << 6) + (size_t)std::countr_zero(m);
            m &= m - 1;
            const uint64_t *rv = rowOf(v);
            for (size_t x = 0; x < Uw; ++x) T[x] = pool[x] & rv[x];
            res += bitSct(T, held + 1, piv, cap, depth + 1);
            if (res >= cap) return res;
            bitClr(pool, v);
        }
    }
    return res;
}

// Count s-cliques through `held` implicit held vertices whose candidate list is
// P0 (global ids; the member predicate is already applied by the caller when
// constructing P0 — never inside the recursion). Builds a local universe over
// P0 and runs the bitset SCT.
static double bitCountList(const std::vector<uint32_t> &P0, int held, double cap) {
    if (held + (int)P0.size() < S) return 0.0;  // == recursion's first check (piv=0)
    buildUniverse(P0.data(), P0.size());
    static std::vector<uint64_t> full;
    full.assign(Uw, 0ull);
    for (size_t w = 0; w + 1 < Uw; ++w) full[w] = ~0ull;
    full[Uw - 1] = (Unb & 63) ? ((1ull << (Unb & 63)) - 1) : ~0ull;
    if (Unb == 0) full[0] = 0ull;
    return bitSct(full.data(), held, 0, cap, 0);
}

static const double INF = std::numeric_limits<double>::infinity();

// |N(x) ∩ N(y)| >= k ? (early-exit merge on global adjacency). Used by the
// OPTIONAL clique-sharing trigger (§6.2 Remark): chain-adjacent partners share
// an s-clique, hence >= s-2 common neighbors, so restricting the closure
// trigger to |N(x)∩N(y)| >= s-2 is COMPLETE (still R* ⊆ C) with fewer tests.
static bool commonAtLeast(uint32_t x, uint32_t y, int k) {
    if (k <= 0) return true;
    const auto &ax = adj[x], &ay = adj[y];
    size_t i = 0, j = 0;
    int c = 0;
    while (i < ax.size() && j < ay.size()) {
        if (ax[i] < ay[j]) ++i;
        else if (ax[i] > ay[j]) ++j;
        else { if (++c >= k) return true; ++i; ++j; }
    }
    return false;
}

int main(int argc, char **argv) {
    if (argc < 6) {
        std::fprintf(stderr,
                     "usage: %s <base.edges> <s> <core_base.tsv> <u> <v>\n", argv[0]);
        return 1;
    }
    S = std::atoi(argv[2]);
    const uint32_t U = (uint32_t)std::atoll(argv[4]);
    const uint32_t V = (uint32_t)std::atoll(argv[5]);

    // ---- load graph (original ids) ----
    FILE *f = std::fopen(argv[1], "r");
    if (!f) { std::perror("graph"); return 1; }
    size_t n, m;
    if (std::fscanf(f, "%zu %zu", &n, &m) != 2) return 1;
    adj.assign(n, {});
    for (size_t i = 0; i < m; ++i) {
        uint32_t a, b;
        if (std::fscanf(f, "%u %u", &a, &b) != 2) break;
        if (a == b) continue;
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
    std::fclose(f);
    for (auto &l : adj) {
        std::sort(l.begin(), l.end());
        l.erase(std::unique(l.begin(), l.end()), l.end());
    }

    // ---- load maintained cores ----
    coreBase.assign(n, 0.0);
    f = std::fopen(argv[3], "r");
    if (!f) { std::perror("cores"); return 1; }
    char line[256];
    while (std::fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        unsigned long long id;
        double c;
        if (std::sscanf(line, "%llu %lf", &id, &c) == 2 && id < n) coreBase[id] = c;
    }
    std::fclose(f);

    nCrT.reserve(1024);
    nCr(0, 0);
    // engine infrastructure (persistent buffers in a maintained process;
    // allocated at load time, outside the timed update)
    g2l.assign(n, 0);
    g2lStamp.assign(n, 0);

    const bool dbg = std::getenv("DYN1S_DEBUG") != nullptr;
    // OPTIONAL clique-sharing trigger (§6.2 Remark). Default OFF = faithful M2
    // (plain N(x) trigger). ON: test y only if |N(x)∩N(y)| >= s-2 (complete,
    // fewer tests). Cores stay bit-identical (C still ⊇ R*); only the candidate
    // stats change. This is the measured, spec-sanctioned admission tightening.
    const bool ctrig = std::getenv("DYN1S_CTRIG") != nullptr;
    size_t CAP = std::max<size_t>(4096, n / 16);
    if (const char *ce = std::getenv("DYN1S_CAP")) CAP = (size_t)std::atoll(ce);

    auto t0 = std::chrono::steady_clock::now();

    // ================= Phase 1 (§7): seeds + early exit =================
    // W = N_G(u) ∩ N_G(v)  (computed in G, before inserting e).
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
    // Count new (s-2)-cliques in G[W]: holding {u,v} (adjacent to all of W and,
    // after insertion, to each other), each (s-2)-clique in G[W] completes into
    // one new s-clique through e. Only the ==0 test is load-bearing (Cor 3c), so
    // cap at 1; the stat is emitted as a 0/1 "creates-new-clique" indicator.
    double newCliques = bitCountList(W, 2, 1.0) > 0.0 ? 1.0 : 0.0;

    auto tp1 = std::chrono::steady_clock::now();
    double p1_us = std::chrono::duration<double, std::micro>(tp1 - t0).count();

    if (newCliques == 0.0) {
        // Corollary 3c: no s-clique of G' contains both u,v => no core changes.
        double us = p1_us;
        if (dbg)
            std::fprintf(stderr,
                         "[dbg] EARLY-EXIT newcliques=0 |W|=%zu p1_us=%.0f\n",
                         W.size(), p1_us);
        std::printf("STATS region=0 pinned=0 tested=0 evicted=0 rounds=0 pops=0 "
                    "changed=0 l0=0 pinned_skipped=0 newcliques=0 "
                    "lambda_intervals=0 lambda_us=0 lambda_empty=0 "
                    "p1_us=%.0f p2_us=0 p3_us=0 insert_us=%.0f\n",
                    p1_us, us);
        return 0;
    }

    // Insert e NOW: all Phase-2/3 counting is in G'.
    adj[U].insert(std::lower_bound(adj[U].begin(), adj[U].end(), V), V);
    adj[V].insert(std::lower_bound(adj[V].begin(), adj[V].end(), U), U);
    tp1 = std::chrono::steady_clock::now();
    p1_us = std::chrono::duration<double, std::micro>(tp1 - t0).count();

    // ===== Phase 1.5 (§18): active-level filter Λ̂ via the seed mini-peel =====
    // Seeds S = {u,v} ∪ W. Pinned peel on (S, b ≡ +∞): boundary never exits, so
    // key(x) = # s-cliques through x in G' avoiding DEAD seeds, and the peel
    // computes UB(x) = ccore_{S,∞}(x) >= c'(x) (Lemma 8). Keys are capped at
    // c(x)+256; a seed that saturates gets UB = ∞ and stays PERMANENTLY alive
    // (sound by Lemma 7); deltas are never applied — unsaturated keys are
    // recomputed from scratch (capped) after each pop, and only for pool
    // members ADJACENT to the popped seed (exact: removing z changes only
    // counts of cliques containing z, all of which pass through neighbors of z).
    // Λ̂ = ∪_x (c(x), UB(x)] covers every riser's c+1 (Lemma 9).
    std::vector<uint32_t> seeds;
    seeds.reserve(W.size() + 2);
    seeds.push_back(U);
    seeds.push_back(V);
    for (uint32_t w : W) seeds.push_back(w);
    const size_t NS = seeds.size();
    std::vector<uint8_t> isDeadSeed(n, 0);
    std::vector<double> sKey(NS, 0.0), sUB(NS, 0.0);
    std::vector<uint8_t> sState(NS, 0);   // 0 alive-unsat, 1 saturated(∞), 2 dead
    auto seedKey = [&](uint32_t x) -> double {
        const double kcap = coreBase[x] + 256.0;
        std::vector<uint32_t> P;
        P.reserve(adj[x].size());
        for (uint32_t y : adj[x])
            if (!isDeadSeed[y]) P.push_back(y);
        return bitCountList(P, 1, kcap);
    };
    std::vector<uint32_t> mpool;          // indices into seeds, alive unsaturated
    for (size_t i = 0; i < NS; ++i) {
        sKey[i] = seedKey(seeds[i]);
        if (sKey[i] >= coreBase[seeds[i]] + 256.0) sState[i] = 1;  // saturated
        else mpool.push_back((uint32_t)i);
    }
    double mMini = 0.0;
    while (!mpool.empty()) {
        size_t bi = 0;
        for (size_t i = 1; i < mpool.size(); ++i)
            if (sKey[mpool[i]] < sKey[mpool[bi]]) bi = i;
        uint32_t ip = mpool[bi];
        mpool[bi] = mpool.back();
        mpool.pop_back();
        if (sKey[ip] > mMini) mMini = sKey[ip];
        sUB[ip] = mMini;
        sState[ip] = 2;
        isDeadSeed[seeds[ip]] = 1;
        const auto &az = adj[seeds[ip]];
        for (uint32_t j : mpool)
            if (std::binary_search(az.begin(), az.end(), seeds[j]))
                sKey[j] = seedKey(seeds[j]);
    }
    // assemble sorted disjoint intervals (lo, hi]
    std::vector<std::pair<double, double>> lambdaIv;
    {
        std::vector<std::pair<double, double>> iv;
        for (size_t i = 0; i < NS; ++i) {
            double lo = coreBase[seeds[i]];
            double hi = (sState[i] == 1) ? INF : sUB[i];
            if (hi > lo) iv.push_back({lo, hi});
        }
        std::sort(iv.begin(), iv.end());
        for (auto &p : iv) {
            if (!lambdaIv.empty() && p.first <= lambdaIv.back().second)
                lambdaIv.back().second = std::max(lambdaIv.back().second, p.second);
            else
                lambdaIv.push_back(p);
        }
    }
    // ℓ ∈ Λ̂ ⟺ some interval has lo < ℓ <= hi  (O(log) membership)
    auto inLambda = [&](double l) -> bool {
        size_t lo = 0, hi = lambdaIv.size();
        while (lo < hi) {              // first interval with .first >= l
            size_t mid = (lo + hi) / 2;
            if (lambdaIv[mid].first < l) lo = mid + 1;
            else hi = mid;
        }
        return lo > 0 && l <= lambdaIv[lo - 1].second;
    };

    auto tl1 = std::chrono::steady_clock::now();
    double lambda_us = std::chrono::duration<double, std::micro>(tl1 - tp1).count();

    if (lambdaIv.empty()) {
        // Corollary 9a: no active level anywhere => R* = ∅, no core changes.
        auto t1 = std::chrono::steady_clock::now();
        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        if (dbg)
            std::fprintf(stderr,
                         "[dbg] LAMBDA-EMPTY seeds=%zu lambda_us=%.0f p1_us=%.0f\n",
                         NS, lambda_us, p1_us);
        std::printf("STATS region=0 pinned=0 tested=0 evicted=0 rounds=0 pops=0 "
                    "changed=0 l0=0 pinned_skipped=0 newcliques=1 "
                    "lambda_intervals=0 lambda_us=%.0f lambda_empty=1 "
                    "p1_us=%.0f p2_us=0 p3_us=0 insert_us=%.0f\n",
                    lambda_us, p1_us, us);
        return 0;
    }

    // ================= Phase 2 (§6): discovery =================
    inRegion.assign(n, 0);        // C membership
    std::vector<uint8_t> tested(n, 0);
    // TS memo: tsVal[z] holds either the EXACT TS(z) (tsExact[z]=1) or a
    // CONFIRMED lower bound TS(z) >= tsVal[z] (tsExact[z]=0). -1 = unset. This
    // caps each TS computation at the actual threshold it is compared against
    // (usually tiny), not a global max.
    std::vector<double> tsVal(n, -1.0);
    std::vector<uint8_t> tsExact(n, 0);
    std::vector<uint32_t> Cvec;             // admitted vertices (may later evict)

    // Exact decision "TS(z) >= thr" (TS = total s-clique support of z in G').
    auto TS_ge = [&](uint32_t z, double thr) -> bool {
        if (tsVal[z] >= 0.0) {
            if (tsExact[z]) return tsVal[z] >= thr;   // exact known
            if (tsVal[z] >= thr) return true;         // confirmed >= tsVal >= thr
            // confirmed only up to tsVal < thr: recompute with a higher cap
        }
        double v = bitCountList(adj[z], 1, thr);
        if (v >= thr) { tsVal[z] = thr; tsExact[z] = 0; return true; }
        tsVal[z] = v; tsExact[z] = 1; return false;    // exact TS(z) = v < thr
    };

    // Static admission test PASS(y): OS(y) >= c(y)+1, where OS counts s-cliques
    // through y whose other members each satisfy c(z) >= ℓ_y OR the optimistic
    // branch. §19 hook 2: the optimistic branch additionally requires
    // c(z)+1 ∈ Λ̂ (a member needed under optimism is a co-riser, in-band by
    // Lemma 9), besides TS(z) >= ℓ_y. TS is lazy+memoized; the OS count is
    // capped at ℓ_y — only "OS >= ℓ_y" matters (exact decision).
    auto PASS = [&](uint32_t y) -> bool {
        const double ly = coreBase[y] + 1.0;
        std::vector<uint32_t> P;
        P.reserve(adj[y].size());
        for (uint32_t z : adj[y])
            if (coreBase[z] >= ly ||
                (inLambda(coreBase[z] + 1.0) && TS_ge(z, ly)))
                P.push_back(z);
        return bitCountList(P, 1, ly) >= ly;
    };

    // seeds
    inRegion[U] = inRegion[V] = 1;
    tested[U] = tested[V] = 1;
    Cvec.push_back(U); Cvec.push_back(V);
    std::vector<uint32_t> frontier = {U, V};
    size_t testedCount = 0;   // vertices tested by PASS (excludes seeds)
    bool fallback = false;

    for (size_t fh = 0; fh < frontier.size(); ++fh) {
        uint32_t x = frontier[fh];
        for (uint32_t y : adj[x]) {
            if (tested[y]) continue;
            // §19 hook 1 (trigger filter): a riser y has c(y)+1 ∈ Λ̂ (Lemma 9);
            // out-of-band y cannot rise — skip AND mark tested (the test is
            // static, independent of the triggering frontier vertex).
            if (!inLambda(coreBase[y] + 1.0)) { tested[y] = 1; continue; }
            // clique-sharing trigger: skip (WITHOUT marking tested, so a better
            // frontier neighbor may still reach y) unless y shares >= s-2 common
            // neighbors with x (necessary to co-inhabit an s-clique with x).
            if (ctrig && !commonAtLeast(x, y, S - 2)) continue;
            tested[y] = 1;
            ++testedCount;
            if (PASS(y)) {
                inRegion[y] = 1;
                Cvec.push_back(y);
                frontier.push_back(y);
                if (Cvec.size() > CAP) { fallback = true; break; }
            }
        }
        if (fallback) break;
    }

    auto tp2a = std::chrono::steady_clock::now();
    double p2test_us = std::chrono::duration<double, std::micro>(tp2a - tl1).count();

    if (fallback) {
        auto t1 = std::chrono::steady_clock::now();
        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        double p2_us = std::chrono::duration<double, std::micro>(t1 - tl1).count();
        if (dbg)
            std::fprintf(stderr,
                         "[dbg] FALLBACK |C|=%zu > CAP=%zu tested=%zu "
                         "p1_us=%.0f p2_us=%.0f\n",
                         Cvec.size(), CAP, testedCount, p1_us, p2_us);
        std::printf("FALLBACK\n");
        std::printf("STATS region=%zu pinned=0 tested=%zu evicted=0 rounds=1 "
                    "pops=0 changed=0 l0=0 pinned_skipped=0 newcliques=%.0f "
                    "lambda_intervals=%zu lambda_us=%.0f lambda_empty=0 "
                    "p1_us=%.0f p2_us=%.0f p3_us=0 insert_us=%.0f\n",
                    Cvec.size(), testedCount, newCliques, lambdaIv.size(),
                    lambda_us, p1_us, p2_us, us);
        return 0;
    }

    // ---- eviction cascade (§6.3), STRICTLY after the closure saturates ----
    // EOS_C(x): s-cliques through x whose other members each satisfy
    // c(z) >= c(x)+1 OR the optimistic branch, which per §19 hook 3 requires
    // z ∈ C AND c(z)+1 ∈ Λ̂ AND TS(z) >= c(x)+1.  Prune x when < c(x)+1.
    size_t evicted = 0;
    {
        auto EOScheck = [&](uint32_t x) -> double {
            const double thr = coreBase[x] + 1.0;
            std::vector<uint32_t> P;
            P.reserve(adj[x].size());
            for (uint32_t z : adj[x])
                if (coreBase[z] >= thr ||
                    (inRegion[z] && inLambda(coreBase[z] + 1.0) && TS_ge(z, thr)))
                    P.push_back(z);
            return bitCountList(P, 1, thr);   // only "EOS >= thr" matters
        };
        std::vector<uint8_t> inQ(n, 0);
        std::vector<uint32_t> wl;
        for (uint32_t x : Cvec)
            if (x != U && x != V) { wl.push_back(x); inQ[x] = 1; }
        for (size_t wh = 0; wh < wl.size(); ++wh) {
            uint32_t x = wl[wh];
            inQ[x] = 0;
            if (!inRegion[x] || x == U || x == V) continue;
            if (EOScheck(x) < coreBase[x] + 1.0) {
                inRegion[x] = 0;   // evict
                ++evicted;
                // only neighbors still in C can lose EOS value -> re-check them
                for (uint32_t y : adj[x])
                    if (inRegion[y] && y != U && y != V && !inQ[y]) {
                        wl.push_back(y);
                        inQ[y] = 1;
                    }
            }
        }
    }

    // final region
    std::vector<uint32_t> region;
    region.reserve(Cvec.size());
    for (uint32_t x : Cvec)
        if (inRegion[x]) region.push_back(x);

    auto tp2 = std::chrono::steady_clock::now();
    double p2evict_us = std::chrono::duration<double, std::micro>(tp2 - tp2a).count();
    double p2_us = std::chrono::duration<double, std::micro>(tp2 - tl1).count();

    // ================= Phase 3 (§8): one pinned peel =================
    // pinned boundary = N(region) \ region, frozen key coreBase.
    std::vector<uint32_t> pinned;
    {
        std::vector<uint8_t> inP(n, 0);
        for (uint32_t x : region)
            for (uint32_t y : adj[x])
                if (!inRegion[y] && !inP[y]) { inP[y] = 1; pinned.push_back(y); }
    }
    std::sort(pinned.begin(), pinned.end(),
              [](uint32_t a, uint32_t b) { return coreBase[a] < coreBase[b]; });
    size_t pinnedSz = pinned.size();

    // L0 fast-start.
    double L0 = INF;
    for (uint32_t x : region) L0 = std::min(L0, coreBase[x]);
    size_t j0 = std::lower_bound(pinned.begin(), pinned.end(), L0,
                                 [](uint32_t a, double v) { return coreBase[a] < v; })
                - pinned.begin();
    size_t pinnedSkipped = j0;

    // Local universe for the peel: region [0, R), then alive pinned (ascending
    // coreBase) [R, Unb). The τ-view "pinned with core >= τ" is then the id
    // range [cut(τ), Unb): a single clearRange folds the predicate into masks.
    const size_t R = region.size();
    {
        std::vector<uint32_t> univ;
        univ.reserve(R + pinned.size() - j0);
        for (uint32_t x : region) univ.push_back(x);
        for (size_t k = j0; k < pinned.size(); ++k) univ.push_back(pinned[k]);
        buildUniverse(univ.data(), univ.size());
    }
    std::vector<double> pinCore;             // ascending, local ids R+i
    pinCore.reserve(pinned.size() - j0);
    for (size_t k = j0; k < pinned.size(); ++k) pinCore.push_back(coreBase[pinned[k]]);
    auto viewCut = [&](double tau) -> size_t {
        return R + (size_t)(std::lower_bound(pinCore.begin(), pinCore.end(), tau) -
                            pinCore.begin());
    };

    std::vector<uint64_t> aliveM(Uw, 0ull);
    for (size_t i = 0; i < Unb; ++i) bitSet(aliveM.data(), i);

    std::vector<uint64_t> Pbuf(Uw), Dbuf(Uw);   // top-level P for key counts

    // key(x) = # s-cliques through x within alive ∩ view(τ_x)
    auto supportOfL = [&](size_t lx) -> double {
        const uint64_t *rx = rowOf(lx);
        for (size_t w = 0; w < Uw; ++w) Pbuf[w] = rx[w] & aliveM[w];
        clearRange(Pbuf.data(), R, viewCut(coreBase[lid2g[lx]]));
        return bitSct(Pbuf.data(), 1, 0, INF, 0);
    };
    // exact removal delta: cliques through {x,z} with other members in x's view;
    // computed BEFORE z leaves aliveM
    auto supportDeltaL = [&](size_t lx, size_t lz) -> double {
        const uint64_t *rx = rowOf(lx), *rz = rowOf(lz);
        for (size_t w = 0; w < Uw; ++w) Dbuf[w] = rx[w] & rz[w] & aliveM[w];
        clearRange(Dbuf.data(), R, viewCut(coreBase[lid2g[lx]]));
        return bitSct(Dbuf.data(), 2, 0, INF, 0);
    };

    std::vector<double> key(R, 0.0), newCoreL(R, -1.0);
    std::vector<uint32_t> pool;                 // local region ids
    pool.reserve(R);
    for (size_t i = 0; i < R; ++i) { key[i] = supportOfL(i); pool.push_back((uint32_t)i); }

    size_t peelPops = 0;

    // Remove z, keeping every alive region neighbor's key exact via delta
    // subtraction (computed BEFORE clearing z from aliveM). τ-view skip: a
    // pinned z below y's level was never counted in y's view.
    auto removeAndUpdate = [&](size_t lz) {
        const bool zPinned = lz >= R;
        const double zCore = zPinned ? pinCore[lz - R] : 0.0;
        const uint64_t *rz = rowOf(lz);
        size_t regWords = (R + 63) / 64;
        for (size_t w = 0; w < regWords; ++w) {
            uint64_t m = rz[w] & aliveM[w];
            if (w == regWords - 1 && (R & 63))
                m &= (1ull << (R & 63)) - 1;    // region ids only
            while (m) {
                size_t ly = (w << 6) + (size_t)std::countr_zero(m);
                m &= m - 1;
                if (zPinned && zCore < coreBase[lid2g[ly]]) continue;
                key[ly] -= supportDeltaL(ly, lz);
                if (key[ly] < 0) {
                    std::fprintf(stderr, "WARN negkey vertex=%u key=%.0f\n",
                                 lid2g[ly], key[ly]);
                    key[ly] = 0;
                }
            }
        }
        bitClr(aliveM.data(), lz);
    };

    size_t jc = 0;        // cursor into pinCore (ascending, past-skipped only)
    double minCore = L0;  // no region vertex can pop below L0
    while (!pool.empty()) {
        size_t bi = 0;
        for (size_t i = 1; i < pool.size(); ++i)
            if (key[pool[i]] < key[pool[bi]]) bi = i;
        double rmin = key[pool[bi]];
        double pmin = (jc < pinCore.size()) ? pinCore[jc] : INF;
        if (pmin <= rmin) {
            size_t lz = R + jc;
            ++jc;
            if (pinCore[lz - R] > minCore) minCore = pinCore[lz - R];
            removeAndUpdate(lz);
            ++peelPops;
        } else {
            uint32_t lz = pool[bi];
            pool[bi] = pool.back();
            pool.pop_back();
            if (key[lz] > minCore) minCore = key[lz];
            newCoreL[lz] = minCore;
            removeAndUpdate(lz);
            ++peelPops;
        }
    }

    auto tp3 = std::chrono::steady_clock::now();
    double p3_us = std::chrono::duration<double, std::micro>(tp3 - tp2).count();
    double us = std::chrono::duration<double, std::micro>(tp3 - t0).count();

    if (dbg)
        std::fprintf(stderr,
                     "[dbg] |C|=%zu region=%zu pinned=%zu tested=%zu evicted=%zu "
                     "skipped=%zu pops=%zu l0=%.0f newcliques=%.0f "
                     "lambda_iv=%zu lambda_us=%.0f "
                     "p1_us=%.0f p2_us=%.0f (test=%.0f evict=%.0f) p3_us=%.0f "
                     "insert_us=%.0f\n",
                     Cvec.size(), region.size(), pinnedSz, testedCount, evicted,
                     pinnedSkipped, peelPops, L0, newCliques, lambdaIv.size(),
                     lambda_us, p1_us, p2_us, p2test_us, p2evict_us, p3_us, us);

    // ---- emit CHANGED + STATS ----
    size_t nChanged = 0;
    for (size_t i = 0; i < R; ++i) {
        uint32_t x = region[i];
        if (newCoreL[i] < coreBase[x])
            std::fprintf(stderr, "WARN drop vertex=%u old=%.0f new=%.0f\n",
                         x, coreBase[x], newCoreL[i]);
        if (newCoreL[i] != coreBase[x]) {
            std::printf("CHANGED %u %.0f %.0f\n", x, coreBase[x], newCoreL[i]);
            ++nChanged;
        }
    }
    std::printf("STATS region=%zu pinned=%zu tested=%zu evicted=%zu rounds=1 "
                "pops=%zu changed=%zu l0=%.0f pinned_skipped=%zu newcliques=%.0f "
                "lambda_intervals=%zu lambda_us=%.0f lambda_empty=0 "
                "p1_us=%.0f p2_us=%.0f p3_us=%.0f insert_us=%.0f\n",
                region.size(), pinnedSz, testedCount, evicted, peelPops, nChanged,
                L0, pinnedSkipped, newCliques, lambdaIv.size(), lambda_us,
                p1_us, p2_us, p3_us, us);
    return 0;
}
