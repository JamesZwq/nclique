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
#include <sys/resource.h>
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

// ==================== v4 persistent clique index (§22-§25) ====================
// Level-s clique forest: flat leaves (H, Π), Π SORTED BY CORE VALUE (ascending,
// ties by id); per-vertex incidence lists (leafId<<1 | isPivot); maintained
// per-vertex total s-clique support[] for the CURRENT graph. Only leaves with
// |H| <= s <= |H|+|Π| are stored (others encode no s-clique). Every clique of
// the graph is encoded at exactly one leaf as H ∪ σ, σ ⊆ Π; each leaf member
// is adjacent to every other member (H clique, Π pairwise adjacent, Π ~ H).
static std::vector<uint32_t> leafData;   // per leaf: H members then Π members
static std::vector<uint64_t> leafStart;  // CSR into leafData (nLeaves+1)
static std::vector<uint8_t> leafHLen;    // |H| (<= s <= 16)
static std::vector<std::vector<uint32_t>> vLeaves; // vertex -> (leafId<<1)|isPiv
static std::vector<double> supportG;     // total s-clique support, current graph

static size_t idxLeaves() { return leafHLen.size(); }

// Append one leaf (global ids). Sorts Π by (core, id); updates incidences and
// support[] attribution (held: C(|Π|, s-|H|); pivot: C(|Π|-1, s-|H|-1)).
static void idxEmitLeaf(const std::vector<uint32_t> &H, std::vector<uint32_t> &Pi) {
    const int hl = (int)H.size(), pl = (int)Pi.size();
    if (hl > S || hl + pl < S) return;               // s-filter
    std::sort(Pi.begin(), Pi.end(), [](uint32_t a, uint32_t b) {
        return coreBase[a] != coreBase[b] ? coreBase[a] < coreBase[b] : a < b;
    });
    uint32_t li = (uint32_t)idxLeaves();
    leafHLen.push_back((uint8_t)hl);
    for (uint32_t w : H) leafData.push_back(w);
    for (uint32_t w : Pi) leafData.push_back(w);
    leafStart.push_back(leafData.size());
    const double aH = nCr(pl, S - hl), aP = nCr(pl - 1, S - hl - 1);
    for (uint32_t w : H) { vLeaves[w].push_back(li << 1); supportG[w] += aH; }
    for (uint32_t w : Pi) { vLeaves[w].push_back((li << 1) | 1u); supportG[w] += aP; }
}

// Emitting SCT recursion on the CURRENT bitset universe (no closed forms —
// leaves must be materialized). H carries GLOBAL ids; PiL carries local ids
// (mapped at emit). Prunes: |H| > s (all deeper leaves have bigger H);
// |H|+|Π|+|P| < s (cannot reach an s-clique).
static void idxEmitSct(const uint64_t *P, std::vector<uint32_t> &H,
                       std::vector<uint32_t> &PiL, int depth) {
    if ((int)H.size() > S) return;
    size_t cnt = pcAll(P);
    if ((int)(H.size() + PiL.size() + cnt) < S) return;
    if (cnt == 0) {
        std::vector<uint32_t> Pi;
        Pi.reserve(PiL.size());
        for (uint32_t l : PiL) Pi.push_back(lid2g[l]);
        idxEmitLeaf(H, Pi);
        return;
    }
    // pivot p = argmax |P ∩ N(p)|
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
    uint64_t *T = getScr(scrT, depth);
    {
        const uint64_t *rp = rowOf(p);
        for (size_t w = 0; w < Uw; ++w) T[w] = P[w] & rp[w];
    }
    PiL.push_back((uint32_t)p);
    idxEmitSct(T, H, PiL, depth + 1);
    PiL.pop_back();
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
            H.push_back(lid2g[v]);
            idxEmitSct(T, H, PiL, depth + 1);
            H.pop_back();
            bitClr(pool, v);
        }
    }
}

// V4.1: build the forest for the CURRENT adjacency (base graph at load).
// Degeneracy order; per-vertex subproblem on later neighbors (small universes,
// bitset engine). Cliques partition by their earliest-ordered member.
static void idxBuild(size_t n) {
    vLeaves.assign(n, {});
    supportG.assign(n, 0.0);
    leafData.clear();
    leafHLen.clear();
    leafStart.assign(1, 0);
    std::vector<uint32_t> deg(n), order;
    order.reserve(n);
    std::vector<int64_t> rank(n, -1);
    {
        uint32_t maxd = 0;
        for (size_t i = 0; i < n; ++i) {
            deg[i] = (uint32_t)adj[i].size();
            maxd = std::max(maxd, deg[i]);
        }
        std::vector<std::vector<uint32_t>> bucket(maxd + 1);
        for (size_t i = 0; i < n; ++i) bucket[deg[i]].push_back((uint32_t)i);
        std::vector<uint8_t> done(n, 0);
        uint32_t cur = 0;
        while (order.size() < n) {
            while (cur <= maxd && bucket[cur].empty()) ++cur;
            if (cur > maxd) break;
            uint32_t v = bucket[cur].back();
            bucket[cur].pop_back();
            if (done[v] || deg[v] != cur) continue;   // stale entry
            done[v] = 1;
            rank[v] = (int64_t)order.size();
            order.push_back(v);
            for (uint32_t w : adj[v])
                if (!done[w]) {
                    --deg[w];
                    bucket[deg[w]].push_back(w);
                    if (deg[w] < cur) cur = deg[w];
                }
        }
    }
    std::vector<uint32_t> H, PiL, P0;
    for (uint32_t v : order) {
        P0.clear();
        for (uint32_t w : adj[v])
            if (rank[w] > rank[v]) P0.push_back(w);
        if (1 + (int)P0.size() < S) continue;         // cannot reach s
        buildUniverse(P0.data(), P0.size());
        static std::vector<uint64_t> full;
        full.assign(Uw, 0ull);
        for (size_t w = 0; w + 1 < Uw; ++w) full[w] = ~0ull;
        full[Uw - 1] = (Unb & 63) ? ((1ull << (Unb & 63)) - 1) : ~0ull;
        H.assign(1, v);
        PiL.clear();
        idxEmitSct(full.data(), H, PiL, 0);
    }
}

// Filtered s-clique count through x, capped, member predicate
// pass(w) = coreBase[w] >= thr OR extra(w). Exact leaf-sum (§25): per leaf,
// H∖{x} must all pass; qualifying Π counted via binary search on the
// core-sorted Π (suffix passes by first disjunct) + extra() scan below thr.
template <class Extra>
static double leafCount(uint32_t x, double thr, double cap, Extra &&extra) {
    double res = 0.0;
    for (uint32_t enc : vLeaves[x]) {
        const uint32_t li = enc >> 1;
        const bool isPiv = enc & 1u;
        const uint32_t *d = leafData.data() + leafStart[li];
        const int hl = leafHLen[li];
        const int pl = (int)(leafStart[li + 1] - leafStart[li]) - hl;
        bool ok = true;
        for (int i = 0; i < hl; ++i) {
            uint32_t w = d[i];
            if (w == x) continue;
            if (!(coreBase[w] >= thr || extra(w))) { ok = false; break; }
        }
        if (!ok) continue;
        const uint32_t *pi = d + hl;
        int lo = 0, hi = pl;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (coreBase[pi[mid]] < thr) lo = mid + 1;
            else hi = mid;
        }
        int q = pl - lo;
        for (int i = 0; i < lo; ++i)
            if (extra(pi[i])) ++q;
        if (!isPiv) {
            res += nCr(q, S - hl);
        } else {
            const bool xq = (coreBase[x] >= thr || extra(x));
            res += nCr(q - (xq ? 1 : 0), S - hl - 1);
        }
        if (res >= cap) return res;
    }
    return res;
}

// ---- V4.3 index-native peel scratch (epoch-stamped, grow-only) ----
static std::vector<uint32_t> leafEpochV;   // per-leaf activation epoch
static std::vector<uint32_t> leafActId;    // per-leaf active id (valid iff epoch)
static uint32_t leafEpoch = 0;

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
    // Streaming mode (V4.4 Tier-2 gate): <base.edges> <s> <core.tsv> --edges <file>
    // processes the file's "u v" lines sequentially, REUSING the maintained
    // index/support/cores across updates (Lemma 10 end-to-end); final output is
    // the accumulated per-vertex diff vs the loaded cores.
    const bool streaming = (std::strcmp(argv[4], "--edges") == 0);
    std::vector<std::pair<uint32_t, uint32_t>> edgeList;
    if (streaming) {
        FILE *ef = std::fopen(argv[5], "r");
        if (!ef) { std::perror("edges"); return 1; }
        unsigned long long a, b;
        while (std::fscanf(ef, "%llu %llu", &a, &b) == 2)
            edgeList.push_back({(uint32_t)a, (uint32_t)b});
        std::fclose(ef);
    } else {
        edgeList.push_back({(uint32_t)std::atoll(argv[4]),
                            (uint32_t)std::atoll(argv[5])});
    }

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

    // ===== V4.1: build the persistent clique index (untimed, like graph load)
    {
        auto tb0 = std::chrono::steady_clock::now();
        idxBuild(n);
        auto tb1 = std::chrono::steady_clock::now();
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
        double rssMB = ru.ru_maxrss / (1024.0 * 1024.0);
#else
        double rssMB = ru.ru_maxrss / 1024.0;
#endif
        std::printf("INDEX leaves=%zu entries=%zu build_ms=%.0f rss_mb=%.0f\n",
                    idxLeaves(), leafData.size(),
                    std::chrono::duration<double, std::milli>(tb1 - tb0).count(),
                    rssMB);
    }

    const bool dbg = std::getenv("DYN1S_DEBUG") != nullptr;
    // OPTIONAL clique-sharing trigger (§6.2 Remark). Default OFF = faithful M2
    // (plain N(x) trigger). ON: test y only if |N(x)∩N(y)| >= s-2 (complete,
    // fewer tests). Cores stay bit-identical (C still ⊇ R*); only the candidate
    // stats change. This is the measured, spec-sanctioned admission tightening.
    const bool ctrig = std::getenv("DYN1S_CTRIG") != nullptr;
    size_t CAP = std::max<size_t>(4096, n / 16);
    if (const char *ce = std::getenv("DYN1S_CAP")) CAP = (size_t)std::atoll(ce);

    // Snapshot for streaming-mode final diff.
    std::vector<double> coreOrig = coreBase;

    for (size_t eIdx = 0; eIdx < edgeList.size(); ++eIdx) {
    const uint32_t U = edgeList[eIdx].first;
    const uint32_t V = edgeList[eIdx].second;
    if (streaming) std::printf("EDGE %u %u\n", U, V);

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
    // Build + APPEND the EdgeTree T_e (Lemma 10): SCT on G[W] with held {u,v}
    // prepended. Emitting it here (Phase 1) makes the index valid for G'
    // BEFORE discovery, so TS(z) = support[z] is O(1) with no side structures.
    // T_e s-filtered leaves exist ⟺ e creates a new s-clique (Cor 3c check).
    size_t leavesBefore = idxLeaves();
    if ((int)W.size() + 2 >= S) {
        buildUniverse(W.data(), W.size());
        static std::vector<uint64_t> fullW;
        fullW.assign(Uw, 0ull);
        for (size_t w = 0; w + 1 < Uw; ++w) fullW[w] = ~0ull;
        fullW[Uw - 1] = (Unb & 63) ? ((1ull << (Unb & 63)) - 1) : ~0ull;
        std::vector<uint32_t> He = {U, V}, PiLe;
        idxEmitSct(fullW.data(), He, PiLe, 0);
    }
    double newCliques = (idxLeaves() > leavesBefore) ? 1.0 : 0.0;

    // Insert e into adjacency NOW (in every case — the graph is G' from here;
    // in streaming mode the next update needs current adjacency even when e
    // changes no core).
    adj[U].insert(std::lower_bound(adj[U].begin(), adj[U].end(), V), V);
    adj[V].insert(std::lower_bound(adj[V].begin(), adj[V].end(), U), U);

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
        continue;
    }

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
    // Index-backed seed key: filtered count over leaves(x), pred = not-dead-seed
    // (no core-threshold shortcut: thr = INF forces the extra() scan on every
    // Π member; #leaves(x) <= support[x] bounds unsaturated seeds, and the cap
    // early-exits saturated ones). Counts in G' (T_e appended in Phase 1).
    auto seedKey = [&](uint32_t x) -> double {
        const double kcap = coreBase[x] + 256.0;
        return leafCount(x, INF, kcap,
                         [&](uint32_t w) { return !isDeadSeed[w]; });
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
        continue;
    }

    // ================= Phase 2 (§6): discovery =================
    inRegion.assign(n, 0);        // C membership
    std::vector<uint8_t> tested(n, 0);
    std::vector<uint32_t> Cvec;             // admitted vertices (may later evict)

    // V4.2: TS(z) in G' = supportG[z], O(1) — the index is G'-valid since the
    // Phase-1 T_e append. Replaces the per-call capped SCT + memo entirely.
    auto TS_ge = [&](uint32_t z, double thr) -> bool {
        return supportG[z] >= thr;
    };

    // Static admission test PASS(y): OS(y) >= c(y)+1, where OS counts s-cliques
    // through y whose other members each satisfy c(z) >= ℓ_y OR the optimistic
    // branch (§19 hook 2: c(z)+1 ∈ Λ̂ AND TS(z) >= ℓ_y). V4.2: evaluated as an
    // exact leaf-sum over leaves(y), capped at ℓ_y (early exit) — the
    // core-threshold disjunct is a suffix of each core-sorted Π.
    auto PASS = [&](uint32_t y) -> bool {
        const double ly = coreBase[y] + 1.0;
        return leafCount(y, ly, ly, [&](uint32_t w) {
                   return inLambda(coreBase[w] + 1.0) && supportG[w] >= ly;
               }) >= ly;
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
        if (streaming) { std::fprintf(stderr, "FATAL streaming fallback at edge %u %u\n", U, V); return 2; }
        return 0;
    }

    // ---- eviction cascade (§6.3 + §20 DELTA-MAINTAINED), STRICTLY after the
    // closure saturates. EOS_C(x): s-cliques through x whose other members each
    // satisfy c(z) >= c(x)+1 OR (z ∈ C AND c(z)+1 ∈ Λ̂ AND TS(z) >= c(x)+1,
    // §19 hook 3). Prune x when EOS < c(x)+1.
    //
    // §20: ev[x] is EXACT below cap_x = c(x)+1+MARGIN (evExact=1) or a
    // confirmed ">= cap_x at last compute" (saturated, evExact=0). EOS_C(x)
    // changes ONLY when a NEIGHBOR z of x is evicted (all members of counted
    // cliques are neighbors), and only if z's qualification was C-DEPENDENT:
    // z qualified via the optimistic branch and NOT via c(z) >= thr_x (a
    // first-disjunct member survives the eviction as a member — subtracting
    // for it would leave ev stale-LOW and wrongly evict). On such a touch:
    //   exact ev[y]: subtract the pair count through {y,z} under y's CURRENT
    //     predicate, computed BEFORE dropping z (cap = ev[y]+1: a delta never
    //     exceeds ev[y], so the capped count is always exact);
    //   saturated ev[y]: recompute (capped) AFTER dropping z.
    // All decisions (ev[x] < c(x)+1) are made on exact-at-decision-time
    // values. The fixpoint is order-independent (EOS is monotone in C: an
    // evictable vertex stays evictable), so the final region equals the
    // recompute-cascade's.
    size_t evicted = 0;
    size_t evDeltas = 0, evRecomputes = 0;
    {
        double MARGIN = 8.0;   // §20 exact-window width (DYN1S_EVMARGIN to tune)
        if (const char *me = std::getenv("DYN1S_EVMARGIN")) MARGIN = std::atof(me);
        // Cost-gated hybrid EOS (both branches EXACT ⟹ identical decisions):
        // leaf-sum evaluation must scan x's whole incidence when the test
        // FAILS (fixed mass), which regresses on dense-cluster hubs whose
        // filtered bitset universe shrinks as C collapses. Use leaves when the
        // incidence is small, bitset otherwise (measured crossover ~1k).
        size_t EOSLEAF_MAX = 1024;
        if (const char *le = std::getenv("DYN1S_EOSLEAF_MAX"))
            EOSLEAF_MAX = (size_t)std::atoll(le);
        auto EOScnt = [&](uint32_t x) -> double {   // capped at c(x)+1+MARGIN
            const double thr = coreBase[x] + 1.0;
            auto extra = [&](uint32_t w) {
                return inRegion[w] && inLambda(coreBase[w] + 1.0) &&
                       supportG[w] >= thr;
            };
            if (vLeaves[x].size() <= EOSLEAF_MAX)
                return leafCount(x, thr, thr + MARGIN, extra);
            std::vector<uint32_t> P;
            P.reserve(adj[x].size());
            for (uint32_t z : adj[x])
                if (coreBase[z] >= thr || extra(z)) P.push_back(z);
            return bitCountList(P, 1, thr + MARGIN);
        };
        std::vector<double> ev(n, 0.0);
        std::vector<uint8_t> evExact(n, 0), inQ(n, 0), evDirty(n, 0);
        std::vector<uint32_t> Q;
        for (uint32_t x : Cvec) {
            if (x == U || x == V) continue;
            ev[x] = EOScnt(x);
            evExact[x] = (ev[x] < coreBase[x] + 1.0 + MARGIN) ? 1 : 0;
            if (evExact[x] && ev[x] < coreBase[x] + 1.0) { Q.push_back(x); inQ[x] = 1; }
        }
        for (size_t qh = 0; qh < Q.size(); ++qh) {
            uint32_t x = Q[qh];
            inQ[x] = 0;
            if (!inRegion[x] || x == U || x == V) continue;
            // LAZY batched recompute: a saturated entry touched while queued is
            // recomputed ONCE here, covering every touch accumulated in the
            // meantime (the FIFO between enqueue and pop is the batching
            // window). Decisions below stay exact-at-decision-time.
            if (evDirty[x]) {
                ev[x] = EOScnt(x);
                ++evRecomputes;
                evExact[x] = (ev[x] < coreBase[x] + 1.0 + MARGIN) ? 1 : 0;
                evDirty[x] = 0;
            }
            if (!(evExact[x] && ev[x] < coreBase[x] + 1.0)) continue;  // exact at decision time
            ++evicted;
            for (uint32_t y : adj[x]) {
                if (!inRegion[y] || y == U || y == V) continue;
                const double thry = coreBase[y] + 1.0;
                if (coreBase[x] >= thry) continue;   // x stays a valid member
                if (!(inLambda(coreBase[x] + 1.0) && TS_ge(x, thry)))
                    continue;                        // x was never counted for y
                // V4.2: recompute-only lazy cascade. Stage-B measured the exact
                // pair-delta path firing negligibly (34 deltas vs 10k
                // recomputes on the worst edge) and the leaf-based EOScnt makes
                // recomputes cheap — so every qualifying touch just defers ONE
                // batched recompute to pop time (exact-at-decision unchanged).
                evDirty[y] = 1;
                if (!inQ[y]) { Q.push_back(y); inQ[y] = 1; }
            }
            inRegion[x] = 0;                          // drop x's C-membership
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

    const size_t R = region.size();
    std::vector<double> pinCore;             // ascending, index i = rank j0+i
    pinCore.reserve(pinned.size() - j0);
    for (size_t k = j0; k < pinned.size(); ++k) pinCore.push_back(coreBase[pinned[k]]);

    std::vector<double> key(R, 0.0), newCoreL(R, -1.0);
    std::vector<uint32_t> pool;                 // region local ids
    pool.reserve(R);
    for (size_t i = 0; i < R; ++i) pool.push_back((uint32_t)i);
    size_t peelPops = 0;
    size_t jc = 0;        // cursor into pinCore (rank = j0 + jc)
    double minCore = L0;  // no region vertex can pop below L0

    const char *peelEnv = std::getenv("DYN1S_PEEL");
    const bool bitsetPeel = peelEnv && std::strcmp(peelEnv, "bitset") == 0;

    if (bitsetPeel) {
        // ---- M1/Stage-B bitset peel (kept for A/B; DYN1S_PEEL=bitset) ----
        // Local universe: region [0, R), then alive pinned (ascending core)
        // [R, Unb). τ-view = clearRange over the pinned id range.
        {
            std::vector<uint32_t> univ;
            univ.reserve(R + pinned.size() - j0);
            for (uint32_t x : region) univ.push_back(x);
            for (size_t k = j0; k < pinned.size(); ++k) univ.push_back(pinned[k]);
            buildUniverse(univ.data(), univ.size());
        }
        auto viewCut = [&](double tau) -> size_t {
            return R + (size_t)(std::lower_bound(pinCore.begin(), pinCore.end(), tau) -
                                pinCore.begin());
        };
        std::vector<uint64_t> aliveM(Uw, 0ull);
        for (size_t i = 0; i < Unb; ++i) bitSet(aliveM.data(), i);
        std::vector<uint64_t> Pbuf(Uw), Dbuf(Uw);
        auto supportOfL = [&](size_t lx) -> double {
            const uint64_t *rx = rowOf(lx);
            for (size_t w = 0; w < Uw; ++w) Pbuf[w] = rx[w] & aliveM[w];
            clearRange(Pbuf.data(), R, viewCut(coreBase[lid2g[lx]]));
            return bitSct(Pbuf.data(), 1, 0, INF, 0);
        };
        auto supportDeltaL = [&](size_t lx, size_t lz) -> double {
            const uint64_t *rx = rowOf(lx), *rz = rowOf(lz);
            for (size_t w = 0; w < Uw; ++w) Dbuf[w] = rx[w] & rz[w] & aliveM[w];
            clearRange(Dbuf.data(), R, viewCut(coreBase[lid2g[lx]]));
            return bitSct(Dbuf.data(), 2, 0, INF, 0);
        };
        for (size_t i = 0; i < R; ++i) key[i] = supportOfL(i);
        auto removeAndUpdate = [&](size_t lz) {
            const bool zPinned = lz >= R;
            const double zCore = zPinned ? pinCore[lz - R] : 0.0;
            const uint64_t *rz = rowOf(lz);
            size_t regWords = (R + 63) / 64;
            for (size_t w = 0; w < regWords; ++w) {
                uint64_t m = rz[w] & aliveM[w];
                if (w == regWords - 1 && (R & 63))
                    m &= (1ull << (R & 63)) - 1;
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
    } else {
        // ================= V4.3 index-native peel =================
        // Pinned deaths happen in RANK order (rank = index in the core-sorted
        // pinned array; unique per vertex, core-monotone). At cursor rank cr,
        // pinned w is alive ⟺ rank(w) >= cr, and visible to region y ⟺
        // rank(w) >= rankTau[y] (:= first rank with core >= τ_y). So each
        // leaf's visible pinned-pivot count is ONE suffix count over its
        // rank-sorted pinned-pivot list at cutoff m = max(cr, rankTau[y]), and
        // its pinned-H members are all visible ⟺ minPHR >= m. Rank cutoffs
        // move by +1 per pinned death and affect ONLY leaves containing that
        // rank — so pinned deaths need no counter mutation, just key deltas
        // over the dying vertex's active leaves. Region deaths update the
        // leaf's dH/aR counters. Keys stay exact at every pop (§16a).
        ++stampCur;   // vertex -> (region lid | pinned rank) map, via g2l stamps
        for (size_t i = 0; i < R; ++i) {
            g2l[region[i]] = (uint32_t)(i << 1);            // even: region lid
            g2lStamp[region[i]] = stampCur;
        }
        for (size_t r = 0; r < pinned.size(); ++r) {
            g2l[pinned[r]] = (uint32_t)((r << 1) | 1u);     // odd: pinned rank
            g2lStamp[pinned[r]] = stampCur;
        }
        std::vector<uint32_t> rankTau(R);
        for (size_t i = 0; i < R; ++i)
            rankTau[i] = (uint32_t)(std::lower_bound(
                             pinned.begin(), pinned.end(), coreBase[region[i]],
                             [](uint32_t a, double v) { return coreBase[a] < v; }) -
                         pinned.begin());
        if (leafEpochV.size() < idxLeaves()) {
            leafEpochV.resize(idxLeaves(), 0);
            leafActId.resize(idxLeaves(), 0);
        }
        ++leafEpoch;
        struct AL {
            uint32_t k, dH, aR, nPP;
            uint32_t minPHR;
            uint32_t membOff, membCnt, ppOff;
        };
        std::vector<AL> als;
        std::vector<uint32_t> aMemb;    // (region lid << 1) | isH
        std::vector<uint32_t> aPP;      // pinned ranks, sorted asc per leaf
        std::vector<std::vector<uint32_t>> incR(R);   // (actId << 1) | zIsH
        std::vector<std::vector<uint32_t>> incP(pinned.size() - j0); // actId
        const uint32_t NORANK = 0xFFFFFFFFu;
        for (uint32_t x : region) {
            for (uint32_t enc : vLeaves[x]) {
                const uint32_t li = enc >> 1;
                if (leafEpochV[li] == leafEpoch) continue;
                leafEpochV[li] = leafEpoch;
                leafActId[li] = 0xFFFFFFFFu;   // default: not activated
                const uint32_t *d = leafData.data() + leafStart[li];
                const int hl = leafHLen[li];
                const int tot = (int)(leafStart[li + 1] - leafStart[li]);
                AL a;
                a.k = (uint32_t)(S - hl);
                a.dH = 0;
                a.membOff = (uint32_t)aMemb.size();
                a.ppOff = (uint32_t)aPP.size();
                uint32_t minPHR = NORANK, aR = 0;
                bool valid = true;
                for (int i = 0; i < tot && valid; ++i) {
                    uint32_t w = d[i];
                    if (g2lStamp[w] != stampCur) { valid = false; break; }
                    uint32_t e2 = g2l[w];
                    if (i < hl) {                        // H member
                        if (e2 & 1u) minPHR = std::min(minPHR, e2 >> 1);
                        else aMemb.push_back(((e2 >> 1) << 1) | 1u);
                    } else {                             // Π member
                        if (e2 & 1u) aPP.push_back(e2 >> 1);
                        else { aMemb.push_back(((e2 >> 1) << 1)); ++aR; }
                    }
                }
                if (!valid) {   // member outside U: impossible by N[x] ⊆ U; guard
                    std::fprintf(stderr, "WARN leaf member outside U leaf=%u\n", li);
                    aMemb.resize(a.membOff);
                    aPP.resize(a.ppOff);
                    continue;
                }
                a.membCnt = (uint32_t)aMemb.size() - a.membOff;
                a.nPP = (uint32_t)aPP.size() - a.ppOff;
                a.minPHR = minPHR;
                a.aR = aR;
                std::sort(aPP.begin() + a.ppOff, aPP.end());
                uint32_t actId = (uint32_t)als.size();
                leafActId[li] = actId;
                als.push_back(a);
                for (uint32_t idx = a.membOff; idx < a.membOff + a.membCnt; ++idx)
                    incR[aMemb[idx] >> 1].push_back((actId << 1) | (aMemb[idx] & 1u));
                for (uint32_t idx = a.ppOff; idx < a.ppOff + a.nPP; ++idx)
                    if (aPP[idx] >= j0) incP[aPP[idx] - j0].push_back(actId);
                if (minPHR != NORANK && minPHR >= j0)
                    incP[minPHR - j0].push_back(actId);   // H-death event carrier
            }
        }
        // contribution of leaf a to region y's key at cursor rank cr;
        // aRd = temporary delta on aR (for before/after of a region-Π death)
        auto contrib = [&](const AL &a, uint32_t yl, bool yIsH, uint32_t cr,
                           int aRd) -> double {
            if (a.dH) return 0.0;
            uint32_t m = std::max(cr, rankTau[yl]);
            if (a.minPHR != NORANK && a.minPHR < m) return 0.0;
            const uint32_t *pp = aPP.data() + a.ppOff;
            size_t lo = std::lower_bound(pp, pp + a.nPP, m) - pp;
            int nvis = (int)a.aR + aRd + (int)(a.nPP - lo);
            return yIsH ? nCr(nvis, a.k) : nCr(nvis - 1, a.k - 1);
        };
        std::vector<uint8_t> aliveR(R, 1);
        const uint32_t rank0 = (uint32_t)j0;
        for (size_t i = 0; i < R; ++i)
            for (uint32_t e : incR[i])
                key[i] += contrib(als[e >> 1], (uint32_t)i, e & 1u, rank0, 0);
        auto guardKey = [&](uint32_t yl) {
            if (key[yl] < 0) {
                std::fprintf(stderr, "WARN negkey vertex=%u key=%.0f\n",
                             region[yl], key[yl]);
                key[yl] = 0;
            }
        };
        while (!pool.empty()) {
            size_t bi = 0;
            for (size_t i = 1; i < pool.size(); ++i)
                if (key[pool[i]] < key[pool[bi]]) bi = i;
            double rmin = key[pool[bi]];
            double pmin = (jc < pinCore.size()) ? pinCore[jc] : INF;
            const uint32_t cr = rank0 + (uint32_t)jc;
            if (pmin <= rmin) {
                // pinned death at rank cr: keys of region members sharing its
                // leaves lose exactly the m: cr -> cr+1 suffix change.
                if (pinCore[jc] > minCore) minCore = pinCore[jc];
                for (uint32_t actId : incP[jc]) {
                    const AL &a = als[actId];
                    for (uint32_t idx = a.membOff; idx < a.membOff + a.membCnt; ++idx) {
                        uint32_t yl = aMemb[idx] >> 1;
                        bool yIsH = aMemb[idx] & 1u;
                        if (!aliveR[yl]) continue;
                        if (rankTau[yl] > cr) continue;   // τ-view skip
                        key[yl] -= contrib(a, yl, yIsH, cr, 0) -
                                   contrib(a, yl, yIsH, cr + 1, 0);
                        guardKey(yl);
                    }
                }
                ++jc;
                ++peelPops;
            } else {
                uint32_t lz = pool[bi];
                pool[bi] = pool.back();
                pool.pop_back();
                if (key[lz] > minCore) minCore = key[lz];
                newCoreL[lz] = minCore;
                aliveR[lz] = 0;
                for (uint32_t e : incR[lz]) {
                    uint32_t actId = e >> 1;
                    bool zIsH = e & 1u;
                    AL &a = als[actId];
                    for (uint32_t idx = a.membOff; idx < a.membOff + a.membCnt; ++idx) {
                        uint32_t yl = aMemb[idx] >> 1;
                        bool yIsH = aMemb[idx] & 1u;
                        if (yl == lz || !aliveR[yl]) continue;
                        double before = contrib(a, yl, yIsH, cr, 0);
                        double after = zIsH ? 0.0 : contrib(a, yl, yIsH, cr, -1);
                        key[yl] -= before - after;
                        guardKey(yl);
                    }
                    if (zIsH) a.dH += 1;
                    else a.aR -= 1;
                }
                ++peelPops;
            }
        }
    }

    auto tp3 = std::chrono::steady_clock::now();
    double p3_us = std::chrono::duration<double, std::micro>(tp3 - tp2).count();
    double us = std::chrono::duration<double, std::micro>(tp3 - t0).count();

    if (dbg)
        std::fprintf(stderr,
                     "[dbg] |C|=%zu region=%zu pinned=%zu tested=%zu evicted=%zu "
                     "skipped=%zu pops=%zu l0=%.0f newcliques=%.0f "
                     "lambda_iv=%zu lambda_us=%.0f ev_deltas=%zu ev_recomputes=%zu "
                     "p1_us=%.0f p2_us=%.0f (test=%.0f evict=%.0f) p3_us=%.0f "
                     "insert_us=%.0f\n",
                     Cvec.size(), region.size(), pinnedSz, testedCount, evicted,
                     pinnedSkipped, peelPops, L0, newCliques, lambdaIv.size(),
                     lambda_us, evDeltas, evRecomputes,
                     p1_us, p2_us, p2test_us, p2evict_us, p3_us, us);

    // ---- emit CHANGED + apply maintained state ----
    size_t nChanged = 0;
    std::vector<uint32_t> changedNow;
    for (size_t i = 0; i < R; ++i) {
        uint32_t x = region[i];
        if (newCoreL[i] < coreBase[x])
            std::fprintf(stderr, "WARN drop vertex=%u old=%.0f new=%.0f\n",
                         x, coreBase[x], newCoreL[i]);
        if (newCoreL[i] != coreBase[x]) {
            if (!streaming)
                std::printf("CHANGED %u %.0f %.0f\n", x, coreBase[x], newCoreL[i]);
            changedNow.push_back(x);
            ++nChanged;
        }
    }
    // Apply new cores to the maintained state, then restore the index
    // invariant "leaf Π sorted by (core, id)" for every leaf containing a
    // changed vertex (V4.2's per-leaf binary searches depend on it; changed
    // sets are tiny so this is O(few leaf sorts) per update).
    for (size_t i = 0; i < R; ++i) {
        uint32_t x = region[i];
        if (newCoreL[i] >= 0.0) coreBase[x] = newCoreL[i];
    }
    for (uint32_t x : changedNow) {
        for (uint32_t enc : vLeaves[x]) {
            uint32_t li = enc >> 1;
            uint32_t *b = leafData.data() + leafStart[li] + leafHLen[li];
            uint32_t *e2 = leafData.data() + leafStart[li + 1];
            std::sort(b, e2, [](uint32_t a, uint32_t c) {
                return coreBase[a] < coreBase[c] ||
                       (coreBase[a] == coreBase[c] && a < c);
            });
        }
    }
    std::printf("STATS region=%zu pinned=%zu tested=%zu evicted=%zu rounds=1 "
                "pops=%zu changed=%zu l0=%.0f pinned_skipped=%zu newcliques=%.0f "
                "lambda_intervals=%zu lambda_us=%.0f lambda_empty=0 "
                "ev_deltas=%zu ev_recomputes=%zu "
                "p1_us=%.0f p2_us=%.0f p2test_us=%.0f p2evict_us=%.0f "
                "p3_us=%.0f insert_us=%.0f\n",
                region.size(), pinnedSz, testedCount, evicted, peelPops, nChanged,
                L0, pinnedSkipped, newCliques, lambdaIv.size(), lambda_us,
                evDeltas, evRecomputes,
                p1_us, p2_us, p2test_us, p2evict_us, p3_us, us);
    }   // end per-edge update loop

    if (streaming) {
        size_t totalChanged = 0;
        for (size_t x = 0; x < n; ++x)
            if (coreBase[x] != coreOrig[x]) {
                std::printf("CHANGED %zu %.0f %.0f\n", x, coreOrig[x], coreBase[x]);
                ++totalChanged;
            }
        std::printf("STREAM_DONE edges=%zu changed=%zu\n",
                    edgeList.size(), totalChanged);
    }
    return 0;
}
