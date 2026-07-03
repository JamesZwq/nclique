// E1 — CPI additivity (Lemma 10) bit-exact check.
//
// Verifies docs/dynamic_v2_spec.md §22 Lemma 10 as a per-vertex DELTA:
//   count_G[x] == count_base[x] + Delta_e[x]   (bit-exact, every vertex x)
// where
//   count_G[x]    = # s-cliques of G containing x,
//   base          = G - e,   e = (u,v) sampled from E(G),
//   count_base[x] = # s-cliques of base containing x,
//   Delta_e[x]    = per-vertex attribution of the edge-local SCT built on
//                   base[W] (W = common nbrs of u,v in base) with held pair
//                   {u,v} prepended, i.e. # s-cliques of G containing x,u,v.
//
// All three quantities are computed by the SAME pivoter/SCT recursion with the
// per-leaf attribution of §11 (held member: C(|Pi|, s-|H|); pivot member:
// C(|Pi|-1, s-|H|-1)), Pascal-triangle doubles throughout. Counts are exact
// integers stored in double (<2^53), so the equality is bit-exact.
//
// Self-contained (no COMMON_SOURCES): own degeneracy order + forward-neighbor
// pivoter + bitset pool. Kept out of the src/ GLOB (lives at repo root).
//
// Usage: test_cpi_additivity <graph.edges> <s> <numEdges> [seed=42]

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>

static int S;                                  // clique size s
static size_t N;                               // vertex count
static std::vector<std::vector<uint32_t>> adj; // sorted adjacency (current graph)

// ---- Pascal-triangle nCr in double (grown on demand), same convention as repo.
static std::vector<std::vector<double>> nCrT;
static double nCr(long long n, long long k) {
    if (k < 0 || n < 0 || k > n) return 0.0;
    while ((long long)nCrT.size() <= n) {
        long long r = (long long)nCrT.size();
        std::vector<double> row(S + 1, 0.0);
        row[0] = 1.0;
        for (int j = 1; j <= S && j <= r; ++j)
            row[j] = nCrT[r - 1][j - 1] + (j <= S ? nCrT[r - 1][j] : 0.0);
        nCrT.push_back(std::move(row));
    }
    return nCrT[n][k];
}

// ======================= bitset local universe =======================
// One active local universe at a time (per top-level vertex). Universe = an
// ordered list of GLOBAL ids; rows[i] = adjacency of universe[i] to the rest of
// the universe, as a bitset (Uw words). lid2g maps local->global.
static size_t Uw = 0;
static std::vector<uint64_t> rowsBuf;    // Unb rows x Uw words
static std::vector<uint32_t> lid2g;      // local id -> global id
static std::vector<uint32_t> g2l;        // global -> local (valid iff stamped)
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

// Build a universe over `verts` (global ids). Rows carry FULL mutual adjacency
// among the universe (symmetric). `fwd` is the forward-neighbor list used to
// discover edges cheaply: every intra-universe edge (f,g) is found once via the
// lower-ranked endpoint's fwd list; we set both bits.
static void buildUniverse(const uint32_t *verts, size_t cnt,
                          const std::vector<std::vector<uint32_t>> &fwd) {
    ++stampCur;
    Uw = (cnt + 63) / 64;
    if (Uw == 0) Uw = 1;
    if (lid2g.size() < cnt) lid2g.resize(cnt);
    for (size_t i = 0; i < cnt; ++i) { lid2g[i] = verts[i]; g2l[verts[i]] = (uint32_t)i; g2lStamp[verts[i]] = stampCur; }
    rowsBuf.assign(cnt * Uw, 0ull);
    for (size_t i = 0; i < cnt; ++i) {
        for (uint32_t g : fwd[verts[i]]) {
            if (g2lStamp[g] == stampCur) {
                size_t j = g2l[g];
                bitSet(rowOf(i), j);
                bitSet(rowOf(j), i);
            }
        }
    }
}

// per-depth recursion scratch
static std::vector<std::vector<uint64_t>> scrT, scrPool, scrD;
static inline uint64_t *getScr(std::vector<std::vector<uint64_t>> &s, int depth) {
    if ((int)s.size() <= depth) s.resize(depth + 1);
    if (s[depth].size() < Uw) s[depth].assign(Uw, 0ull);
    return s[depth].data();
}

// Per-vertex s-clique attribution recursion (full pivoter, to P-empty leaves).
// P = candidate bitset (local ids). Hg = held GLOBAL ids, Pig = pivot GLOBAL
// ids. count[] indexed by global id. Also accumulates leaf stats when statOn.
static double *g_count = nullptr;
static bool g_statOn = false;
static long long g_leaves = 0, g_leafSizeSum = 0;

static void attribSct(const uint64_t *P, std::vector<uint32_t> &Hg,
                      std::vector<uint32_t> &Pig, int depth) {
    const int h = (int)Hg.size();
    if (h > S) return;                        // >s required members: no s-clique
    const int pv = (int)Pig.size();
    const size_t cnt = pcAll(P);
    if (h + pv + (int)cnt < S) return;        // max reachable clique < s
    if (cnt == 0) {                           // leaf (H, Pi)
        if (g_statOn) { ++g_leaves; g_leafSizeSum += (h + pv); }
        const double ch = nCr(pv, S - h);         // held vertex share
        if (ch != 0.0) for (uint32_t g : Hg) g_count[g] += ch;
        const double cp = nCr(pv - 1, S - h - 1);  // pivot vertex share
        if (cp != 0.0) for (uint32_t g : Pig) g_count[g] += cp;
        return;
    }
    // pivot p = argmax |P & N(p)|
    size_t p = 0, best = 0; bool first = true;
    for (size_t w = 0; w < Uw; ++w) {
        uint64_t m = P[w];
        while (m) {
            size_t lid = (w << 6) + (size_t)std::countr_zero(m);
            m &= m - 1;
            size_t c = pcAnd(P, rowOf(lid));
            if (first || c > best) { best = c; p = lid; first = false; }
        }
    }
    // pivot branch: recurse on P & N(p), p added to Pi
    uint64_t *T = getScr(scrT, depth);
    const uint64_t *rp = rowOf(p);
    for (size_t w = 0; w < Uw; ++w) T[w] = P[w] & rp[w];
    Pig.push_back(lid2g[p]);
    attribSct(T, Hg, Pig, depth + 1);
    Pig.pop_back();
    // hold branches: v in (P \ {p}) \ N(p), sequentially-shrinking pool
    uint64_t *pool = getScr(scrPool, depth);
    uint64_t *D = getScr(scrD, depth);
    for (size_t w = 0; w < Uw; ++w) { pool[w] = P[w]; D[w] = P[w] & ~rp[w]; }
    bitClr(pool, p);
    bitClr(D, p);
    for (size_t w = 0; w < Uw; ++w) {
        uint64_t m = D[w];
        while (m) {
            size_t v = (w << 6) + (size_t)std::countr_zero(m);
            m &= m - 1;
            const uint64_t *rv = rowOf(v);
            uint64_t *T2 = getScr(scrT, depth);   // reuse T at this depth
            for (size_t x = 0; x < Uw; ++x) T2[x] = pool[x] & rv[x];
            Hg.push_back(lid2g[v]);
            attribSct(T2, Hg, Pig, depth + 1);
            Hg.pop_back();
            bitClr(pool, v);
        }
    }
}

// ---- degeneracy order (min-degree elimination) + forward-neighbor lists ----
static void degeneracyForward(std::vector<std::vector<uint32_t>> &fwd) {
    std::vector<uint32_t> deg(N), rank(N, 0);
    size_t maxd = 0;
    for (size_t v = 0; v < N; ++v) { deg[v] = (uint32_t)adj[v].size(); maxd = std::max(maxd, (size_t)deg[v]); }
    std::vector<std::vector<uint32_t>> bucket(maxd + 1);
    std::vector<size_t> pos(N);
    for (size_t v = 0; v < N; ++v) { pos[v] = bucket[deg[v]].size(); bucket[deg[v]].push_back((uint32_t)v); }
    std::vector<uint8_t> removed(N, 0);
    size_t cur = 0;
    for (size_t iter = 0; iter < N; ++iter) {
        while (cur <= maxd && bucket[cur].empty()) ++cur;
        // pick a min-degree vertex
        uint32_t x = bucket[cur].back();
        bucket[cur].pop_back();
        rank[x] = (uint32_t)iter;
        removed[x] = 1;
        for (uint32_t y : adj[x]) {
            if (removed[y]) continue;
            uint32_t d = deg[y];
            // remove y from bucket[d]
            size_t py = pos[y];
            uint32_t last = bucket[d].back();
            bucket[d][py] = last; pos[last] = py; bucket[d].pop_back();
            deg[y] = d - 1;
            pos[y] = bucket[d - 1].size(); bucket[d - 1].push_back(y);
            if (d - 1 < cur) cur = d - 1;
        }
    }
    fwd.assign(N, {});
    for (size_t v = 0; v < N; ++v)
        for (uint32_t y : adj[v])
            if (rank[v] < rank[y]) fwd[v].push_back(y);
}

// Full per-vertex s-clique count of the current graph into count[] (size N).
static void fullCount(std::vector<double> &count) {
    std::fill(count.begin(), count.end(), 0.0);
    std::vector<std::vector<uint32_t>> fwd;
    degeneracyForward(fwd);
    g_count = count.data();
    g_statOn = false;
    std::vector<uint32_t> Hg, Pig;
    std::vector<uint64_t> full;
    for (size_t x = 0; x < N; ++x) {
        const auto &fx = fwd[x];
        if ((int)fx.size() + 1 < S) continue;       // can't reach an s-clique
        buildUniverse(fx.data(), fx.size(), fwd);
        full.assign(Uw, 0ull);
        for (size_t i = 0; i < fx.size(); ++i) bitSet(full.data(), i);
        Hg.clear(); Pig.clear();
        Hg.push_back((uint32_t)x);                  // x held (min-rank of its cliques)
        attribSct(full.data(), Hg, Pig, 0);
    }
}

// Delta_e[x]: per-vertex attribution of s-cliques through both u,v.
// SCT on base[W] (= G[W]) with held pair {u,v}. Records leaf stats.
static void deltaCount(uint32_t u, uint32_t v, const std::vector<uint32_t> &W,
                       std::vector<double> &delta, long long &leaves,
                       long long &leafSizeSum) {
    std::fill(delta.begin(), delta.end(), 0.0);
    leaves = 0; leafSizeSum = 0;
    if ((int)W.size() + 2 < S) return;              // need >= s-2 in W
    // Forward lists restricted to W (any order works; use ascending id).
    // Build a symmetric universe over W directly from adjacency intersected W.
    ++stampCur;
    for (uint32_t w : W) { g2lStamp[w] = stampCur; g2l[w] = 0; }
    // assign local ids
    for (size_t i = 0; i < W.size(); ++i) g2l[W[i]] = (uint32_t)i;
    Uw = (W.size() + 63) / 64; if (Uw == 0) Uw = 1;
    if (lid2g.size() < W.size()) lid2g.resize(W.size());
    for (size_t i = 0; i < W.size(); ++i) lid2g[i] = W[i];
    rowsBuf.assign(W.size() * Uw, 0ull);
    for (size_t i = 0; i < W.size(); ++i)
        for (uint32_t g : adj[W[i]])
            if (g2lStamp[g] == stampCur) bitSet(rowOf(i), g2l[g]);
    std::vector<uint64_t> full(Uw, 0ull);
    for (size_t i = 0; i < W.size(); ++i) bitSet(full.data(), i);
    g_count = delta.data();
    g_statOn = true; g_leaves = 0; g_leafSizeSum = 0;
    std::vector<uint32_t> Hg = {u, v}, Pig;
    attribSct(full.data(), Hg, Pig, 0);
    leaves = g_leaves; leafSizeSum = g_leafSizeSum;
    g_statOn = false;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: %s <graph.edges> <s> <numEdges> [seed=42]\n", argv[0]);
        return 1;
    }
    const char *fpath = argv[1];
    S = std::atoi(argv[2]);
    const int numEdges = std::atoi(argv[3]);
    const uint64_t seed = (argc >= 5) ? (uint64_t)std::atoll(argv[4]) : 42ull;

    // ---- load graph ----
    FILE *f = std::fopen(fpath, "r");
    if (!f) { std::perror("graph"); return 1; }
    size_t n, m;
    if (std::fscanf(f, "%zu %zu", &n, &m) != 2) return 1;
    N = n;
    adj.assign(n, {});
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    edges.reserve(m);
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
    // unique undirected edge list (a<b)
    for (uint32_t a = 0; a < n; ++a)
        for (uint32_t b : adj[a])
            if (b > a) edges.push_back({a, b});

    nCrT.reserve(1024);
    nCr(0, 0);
    g2l.assign(n, 0);
    g2lStamp.assign(n, 0);

    // ---- reference full count on G ----
    std::vector<double> countG(n, 0.0), countBase(n, 0.0), delta(n, 0.0);
    fullCount(countG);
    double totG = 0.0; for (double c : countG) totG += c;

    // ---- sample edges (seed) ----
    std::vector<uint32_t> idx(edges.size());
    for (uint32_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::mt19937_64 rng(seed);
    std::shuffle(idx.begin(), idx.end(), rng);
    int K = std::min<int>(numEdges, (int)edges.size());

    std::printf("== E1 %s  s=%d  seed=%llu ==\n", fpath, S, (unsigned long long)seed);
    std::printf("n=%zu m_unique=%zu  totalG(sum count_G)=%.0f  (#s-cliques=%.0f)\n",
                n, edges.size(), totG, totG / (double)S);

    long long totMismatch = 0, edgesTested = 0;
    long long sumLeaves = 0, sumLeafSize = 0;
    long long minLeaves = -1, maxLeaves = 0, minLS = -1, maxLS = 0;
    long long firstFailEdges = 0;
    for (int t = 0; t < K; ++t) {
        uint32_t u = edges[idx[t]].first, v = edges[idx[t]].second;
        // base = G - e
        adj[u].erase(std::lower_bound(adj[u].begin(), adj[u].end(), v));
        adj[v].erase(std::lower_bound(adj[v].begin(), adj[v].end(), u));
        fullCount(countBase);
        // W = common nbrs in base (== common nbrs in G, since e joins u,v only)
        std::vector<uint32_t> W;
        {
            const auto &nu = adj[u], &nv = adj[v];
            size_t i = 0, j = 0;
            while (i < nu.size() && j < nv.size()) {
                if (nu[i] < nv[j]) ++i;
                else if (nu[i] > nv[j]) ++j;
                else { W.push_back(nu[i]); ++i; ++j; }
            }
        }
        long long lv, ls;
        deltaCount(u, v, W, delta, lv, ls);
        // restore e
        adj[u].insert(std::lower_bound(adj[u].begin(), adj[u].end(), v), v);
        adj[v].insert(std::lower_bound(adj[v].begin(), adj[v].end(), u), u);

        // per-vertex bit-exact check
        long long mism = 0;
        int shown = 0;
        for (size_t x = 0; x < n; ++x) {
            if (countBase[x] + delta[x] != countG[x]) {
                ++mism;
                if (shown < 5) {
                    std::printf("  MISMATCH edge=(%u,%u) x=%zu base=%.0f delta=%.0f "
                                "base+delta=%.0f G=%.0f\n",
                                u, v, x, countBase[x], delta[x],
                                countBase[x] + delta[x], countG[x]);
                    ++shown;
                }
            }
        }
        totMismatch += mism;
        ++edgesTested;
        if (mism) ++firstFailEdges;
        sumLeaves += lv; sumLeafSize += ls;
        if (minLeaves < 0 || lv < minLeaves) minLeaves = lv;
        if (lv > maxLeaves) maxLeaves = lv;
        if (minLS < 0 || ls < minLS) minLS = ls;
        if (ls > maxLS) maxLS = ls;
    }

    std::printf("edges_tested=%lld  mismatched_vertices_total=%lld  failing_edges=%lld\n",
                edgesTested, totMismatch, firstFailEdges);
    std::printf("edge-tree stats: leaves[avg=%.2f min=%lld max=%lld]  "
                "sum_leaf_size[avg=%.2f min=%lld max=%lld]\n",
                edgesTested ? (double)sumLeaves / (double)edgesTested : 0.0,
                minLeaves < 0 ? 0 : minLeaves, maxLeaves,
                edgesTested ? (double)sumLeafSize / (double)edgesTested : 0.0,
                minLS < 0 ? 0 : minLS, maxLS);
    std::printf("VERDICT: %s\n", totMismatch == 0 ? "PASS (0 mismatches)" : "FAIL");
    return totMismatch == 0 ? 0 : 2;
}
