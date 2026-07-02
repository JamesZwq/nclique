// Dynamic (1,s)-core maintenance prototype — insert-only v0.
//
// Persistent state between updates: core[] + adjacency ONLY (no CPI/SDCT).
// Insert (u,v):
//   1. seed region A = {u,v} ∪ (N(u)∩N(v))
//   2. scoped re-peel on G[A ∪ N(A)] with pinned boundary: boundary vertex w
//      keeps fixed key core_old(w) (removed when the peel level passes it,
//      never assigned); region supports are recounted per-vertex with a
//      pivot (SCT) recursion over N(x)∩alive.
//   3. verify-and-expand: riser x (c_x -> c'_x) may lift an outside
//      neighbor y iff c_x <= core_old(y) < c'_x (alive-at-level flip);
//      add such y to A and redo until fixpoint.
//
// Usage:
//   dynamic_1s_core <base.edges> <s> <core_base.tsv> <u> <v>
// core_base.tsv: "<orig_id>\t<core>" lines ('#' comments ok); missing id = 0.
// Output: "CHANGED x old new" lines + one "STATS ..." line.

#include <algorithm>
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
static std::vector<uint8_t> alive;             // peel-alive flag
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

// ---- per-vertex support: count s-cliques containing x within alive set ----
// Pivoter recursion; P sorted; held = |H| (x included), piv = |Π|.
static double sctCount(std::vector<uint32_t> &P, int held, int piv) {
    if (held == S) return 1.0;
    if (held + piv + (int)P.size() < S) return 0.0;
    if (P.empty()) return nCr(piv, S - held);
    // Exact closed form one level above the leaves: with S-held == 1, each
    // pivot and each P member individually completes the held set into one
    // s-clique (all are adjacent to every held vertex by construction), so
    // the count is piv + |P|. Provable by induction on the recursion below:
    // pivot branch gives (piv+1) + |P∩N(p)|, hold branches give one clique
    // per v in P\N(p)\{p}, and |P∩N(p)| + |P\N(p)| = |P| with p in P\N(p).
    if (held == S - 1) return (double)piv + (double)P.size();
    // pivot p = argmax |P ∩ N(p)|
    uint32_t p = P[0];
    size_t best = 0;
    std::vector<uint32_t> tmp;
    for (uint32_t c : P) {
        size_t cnt = 0;
        const auto &nc = adj[c];
        size_t i = 0, j = 0;
        while (i < P.size() && j < nc.size()) {
            if (P[i] < nc[j]) ++i;
            else if (P[i] > nc[j]) ++j;
            else { ++cnt; ++i; ++j; }
        }
        if (cnt > best || (cnt == best && c == P[0])) { best = cnt; p = c; }
    }
    auto inter = [](const std::vector<uint32_t> &a, const std::vector<uint32_t> &b,
                    std::vector<uint32_t> &out) {
        out.clear();
        size_t i = 0, j = 0;
        while (i < a.size() && j < b.size()) {
            if (a[i] < b[j]) ++i;
            else if (a[i] > b[j]) ++j;
            else { out.push_back(a[i]); ++i; ++j; }
        }
    };
    double res = 0.0;
    std::vector<uint32_t> sub;
    inter(P, adj[p], sub);
    res += sctCount(sub, held, piv + 1);
    // hold branches: v in P \ N(p) \ {p}, sequentially removed pool
    std::vector<uint32_t> pool = P;
    pool.erase(std::find(pool.begin(), pool.end(), p));
    std::vector<uint32_t> branch;
    {
        const auto &np = adj[p];
        for (uint32_t v : pool) {
            if (!std::binary_search(np.begin(), np.end(), v)) branch.push_back(v);
        }
    }
    for (uint32_t v : branch) {
        std::vector<uint32_t> pv;
        inter(pool, adj[v], pv);
        res += sctCount(pv, held + 1, piv);
        pool.erase(std::find(pool.begin(), pool.end(), v));
    }
    return res;
}

static double supportOf(uint32_t x) {
    std::vector<uint32_t> P;
    P.reserve(adj[x].size());
    for (uint32_t y : adj[x])
        if (alive[y]) P.push_back(y);
    return sctCount(P, 1, 0);
}

// Exact per-removal delta: # s-cliques containing BOTH x and z whose other
// members are all alive. Must be called BEFORE setting alive[z]=0 (the other
// S-2 members must be alive); z itself is excluded automatically (z ∉ N(z)),
// as is x (x ∉ N(x)). held = {x,z}, pool P = N(x) ∩ N(z) ∩ alive.
static double supportDelta(uint32_t x, uint32_t z) {
    std::vector<uint32_t> P;
    const auto &ax = adj[x], &az = adj[z];
    size_t i = 0, j = 0;
    while (i < ax.size() && j < az.size()) {
        if (ax[i] < az[j]) ++i;
        else if (ax[i] > az[j]) ++j;
        else { if (alive[ax[i]]) P.push_back(ax[i]); ++i; ++j; }
    }
    return sctCount(P, 2, 0);
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

    auto t0 = std::chrono::steady_clock::now();

    // ---- insert edge ----
    adj[U].insert(std::lower_bound(adj[U].begin(), adj[U].end(), V), V);
    adj[V].insert(std::lower_bound(adj[V].begin(), adj[V].end(), U), U);

    // ---- seed region: {u,v} ∪ (N(u)∩N(v)) ----
    std::vector<uint8_t> inA(n, 0);
    std::vector<uint32_t> A;
    auto addA = [&](uint32_t x) { if (!inA[x]) { inA[x] = 1; A.push_back(x); } };
    addA(U); addA(V);
    {
        const auto &nu = adj[U], &nv = adj[V];
        size_t i = 0, j = 0;
        while (i < nu.size() && j < nv.size()) {
            if (nu[i] < nv[j]) ++i;
            else if (nu[i] > nv[j]) ++j;
            else { addA(nu[i]); ++i; ++j; }
        }
    }

    std::vector<double> newCore(n, -1.0); // computed cores for region members
    int rounds = 0;
    size_t pinnedSz = 0, peelPops = 0;

    while (true) {
        ++rounds;
        // pinned boundary = N(A) \ A
        std::vector<uint32_t> pinned;
        std::vector<uint8_t> inP(n, 0);
        for (uint32_t x : A)
            for (uint32_t y : adj[x])
                if (!inA[y] && !inP[y]) { inP[y] = 1; pinned.push_back(y); }
        pinnedSz = pinned.size();
        // sort pinned ascending by frozen core_old so the peel can batch-drain
        // same-level runs of them instead of scanning them in the min-pool.
        std::sort(pinned.begin(), pinned.end(),
                  [](uint32_t a, uint32_t b) { return coreBase[a] < coreBase[b]; });

        alive.assign(n, 0);
        for (uint32_t x : A) alive[x] = 1;
        for (uint32_t w : pinned) alive[w] = 1;

        // pool = REGION vertices only; keys = live support. Initial keys are
        // one full supportOf per region vertex; afterwards keys stay exact
        // via per-removal delta subtraction. Pinned vertices never enter the
        // scan-min pool: they are consumed via the sorted `pinned`/`j`
        // cursor, one per iteration, at their frozen coreBase level.
        std::vector<uint32_t> pool;
        std::vector<double> key(n, 0.0);
        for (uint32_t x : A) { key[x] = supportOf(x); pool.push_back(x); }

        size_t j = 0; // cursor into sorted `pinned`

        // Remove z from the alive set, keeping every alive REGION neighbor's
        // key exact via delta subtraction: x loses exactly the s-cliques that
        // contain both x and z (all other members alive). Deltas are computed
        // BEFORE alive[z]=0, so the counted cliques' other members are alive;
        // z is excluded from its own count automatically (z ∉ N(z)). Exact at
        // every step — no verification or rollback needed, and each delta is
        // only a common-neighborhood-sized pivoter call.
        auto removeAndUpdate = [&](uint32_t z) {
            for (uint32_t y : adj[z]) {
                if (alive[y] && inA[y]) {
                    key[y] -= supportDelta(y, z);
                    if (key[y] < 0) {
                        std::fprintf(stderr, "WARN negkey vertex=%u key=%.0f\n",
                                     y, key[y]);
                        key[y] = 0;
                    }
                }
            }
            alive[z] = 0;
        };

        double minCore = 0.0;
        while (!pool.empty()) {
            size_t bi = 0;
            for (size_t i = 1; i < pool.size(); ++i)
                if (key[pool[i]] < key[pool[bi]]) bi = i;
            double rmin = key[pool[bi]];
            double pmin = (j < pinned.size()) ? coreBase[pinned[j]]
                                               : std::numeric_limits<double>::infinity();

            if (pmin <= rmin) {
                // remove ONE pinned vertex at its frozen level
                uint32_t z = pinned[j++];
                if (coreBase[z] > minCore) minCore = coreBase[z];
                removeAndUpdate(z);
                ++peelPops;
            } else {
                // pop the region argmin
                uint32_t z = pool[bi];
                pool[bi] = pool.back();
                pool.pop_back();
                if (key[z] > minCore) minCore = key[z];
                newCore[z] = minCore;
                removeAndUpdate(z);
                ++peelPops;
            }
        }

        // verify-and-expand: riser x may lift outside y iff
        // core_old(x) <= core_old(y) <= core_new(x)
        std::vector<uint32_t> added;
        for (uint32_t x : A) {
            if (newCore[x] <= coreBase[x]) continue;
            for (uint32_t y : adj[x])
                if (!inA[y] && coreBase[x] <= coreBase[y] && coreBase[y] <= newCore[x])
                    added.push_back(y);
        }

        // Gated same-level optimistic promotion (zero-net-rise deadlock
        // breaker): a non-riser region vertex x can be a SILENT riser whose
        // true gain was cancelled by same-level pinned neighbors dying at
        // their frozen coreBase == coreBase[x]; since x never rises, the
        // riser-window rule above never promotes them. Test x under
        // optimism: count s-cliques through x restricted to neighbors that
        // could plausibly be alive at level c+1 — pinned survivors
        // (coreBase >= c), region members that ended >= c+1, and stuck
        // same-level region peers (optimistically assumed to rise together).
        // If even that optimistic support reaches c+1, x's same-level pinned
        // neighbors are the only possible blockers: promote them into A and
        // redo the round so they get peeled dynamically instead of frozen.
        for (uint32_t x : A) {
            if (newCore[x] != coreBase[x]) continue;
            double c = coreBase[x];
            // short-circuit: no same-level pinned neighbor => the promotion
            // set is empty whatever hyp says; skip the expensive count.
            bool hasCandidate = false;
            for (uint32_t y : adj[x])
                if (!inA[y] && coreBase[y] == c) { hasCandidate = true; break; }
            if (!hasCandidate) continue;
            std::vector<uint32_t> P;
            P.reserve(adj[x].size());
            for (uint32_t y : adj[x]) {
                if (inA[y]) {
                    // any region peer that REACHED level c is optimistically
                    // alive at c+1 (covers full risers, stuck same-level
                    // peers, and partial risers capped at exactly c that may
                    // rise further in the mutual-rise group).
                    if (newCore[y] >= c) P.push_back(y);
                } else if (coreBase[y] >= c) {
                    P.push_back(y);
                }
            }
            double hyp = sctCount(P, 1, 0);
            if (hyp >= c + 1) {
                for (uint32_t y : adj[x])
                    if (!inA[y] && coreBase[y] == c)
                        added.push_back(y);
            }
        }

        if (added.empty()) break;
        if (rounds >= 40) {
            // safety cap: bail out; the driver treats FALLBACK edges as
            // needing a full recompute. No CHANGED lines in this case.
            auto t1f = std::chrono::steady_clock::now();
            double usf =
                std::chrono::duration<double, std::micro>(t1f - t0).count();
            std::printf("FALLBACK\n");
            std::printf("STATS region=%zu pinned=%zu rounds=%d pops=%zu "
                        "changed=0 insert_us=%.0f\n",
                        A.size(), pinnedSz, rounds, peelPops, usf);
            return 0;
        }
        for (uint32_t y : added) addA(y);
    }

    auto t1 = std::chrono::steady_clock::now();
    double us = std::chrono::duration<double, std::micro>(t1 - t0).count();

    size_t nChanged = 0;
    for (uint32_t x : A) {
        if (newCore[x] < coreBase[x])
            std::fprintf(stderr, "WARN drop vertex=%u old=%.0f new=%.0f\n",
                         x, coreBase[x], newCore[x]);
        if (newCore[x] != coreBase[x]) {
            std::printf("CHANGED %u %.0f %.0f\n", x, coreBase[x], newCore[x]);
            ++nChanged;
        }
    }
    std::printf("STATS region=%zu pinned=%zu rounds=%d pops=%zu changed=%zu insert_us=%.0f\n",
                A.size(), pinnedSz, rounds, peelPops, nChanged, us);
    return 0;
}
