// E2 — deletion surgery cost distribution (docs/dynamic_v2_spec.md §23, §26).
//
// For 300 edges (u,v) sampled (seed 42) from E(G), record over the level-s CPI
// leaves that contain BOTH u and v (the leaves §23's surgery must touch):
//   leafCount(u,v) = |leaves(u) ∩ leaves(v)|
//   sigmaL(u,v)    = Σ_{L : u,v ∈ L} |L|   (|L| = |H|+|Π|)
// via ONE SDCT walk per (graph,s) using the real callback-based walker
// (SDCT_Augmented_NoTree, bkRecurse_NoTree onLeaf), so the leaf set is exactly
// the index v4 would maintain. Then report median/p90/p99/max of both.
//
// A leaf "contains" vertex x iff x ∈ keepV(H) ∪ dropV(Π). Edges are sampled in
// the internal (degeneracy-relabeled) id space — an edge is an edge, the cost
// distribution is invariant to relabeling. Kept out of the src/ GLOB.
//
// Usage: test_e2_surgery <graph.edges> <s> [numPairs=300] [seed=42]

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#include <utility>
#include <iostream>

#include "src/graph/Graph.h"
#include "src/degeneracy_algorithm_cliques_V.h"
#include "src/misc.h"
#include "src/Global/Global.h"
#include "src/SDCT_Augmented.h"

extern double nCr[1001][401];

static double pct(std::vector<long long> v, double p) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    if (p <= 0) return (double)v.front();
    if (p >= 1) return (double)v.back();
    double r = p * (double)(v.size() - 1);
    size_t lo = (size_t)r;
    double frac = r - (double)lo;
    if (lo + 1 >= v.size()) return (double)v[lo];
    return (double)v[lo] * (1.0 - frac) + (double)v[lo + 1] * frac;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <graph.edges> <s> [numPairs=300] [seed=42]\n", argv[0]);
        return 1;
    }
    const char *fpath = argv[1];
    const int s = std::atoi(argv[2]);
    const int numPairs = (argc >= 4) ? std::atoi(argv[3]) : 300;
    const uint64_t seed = (argc >= 5) ? (uint64_t)std::atoll(argv[4]) : 42ull;

    Graph g(fpath);
    populate_nCr();
    const size_t n = g.getGraphNodeSize();
    daf::vListMap.resize((daf::Size)n + 1);
    std::memset(daf::vListMap.data(), -1, ((size_t)n + 1) * sizeof(daf::Size));
    g.sortByDegeneracyOrder();

    // ---- internal-id unique edge list (a<b) ----
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    for (uint32_t v = 0; v < n; ++v)
        for (daf::Size w : g.getNbrVec(v))
            if (w > v) edges.push_back({v, (uint32_t)w});

    // ---- sample pairs (seed) ----
    std::vector<uint32_t> idx(edges.size());
    for (uint32_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::mt19937_64 rng(seed);
    std::shuffle(idx.begin(), idx.end(), rng);
    int K = std::min<int>(numPairs, (int)edges.size());

    // ---- pair indices ----
    std::vector<uint8_t> isQueried(n, 0);
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> pairsOf(n); // under smaller endpoint: (other, pid)
    for (int p = 0; p < K; ++p) {
        uint32_t a = edges[idx[p]].first, b = edges[idx[p]].second; // a<b
        isQueried[a] = 1; isQueried[b] = 1;
        pairsOf[a].push_back({b, (uint32_t)p});
    }

    std::vector<long long> leafCount(K, 0), sigmaL(K, 0);

    // ---- single SDCT walk (level-s CPI leaves) ----
    std::vector<uint32_t> stampArr(n, 0);
    uint32_t cur = 0;
    std::vector<uint32_t> present;
    present.reserve(64);
    long long totalLeaves = 0;

    auto onLeaf = [&](daf::Size /*leafId*/,
                      const daf::StaticVector<int> &keepV,
                      const daf::StaticVector<int> &dropV) {
        ++totalLeaves;
        ++cur;
        present.clear();
        const long long leafSize = (long long)keepV.size() + (long long)dropV.size();
        for (daf::Size i = 0; i < keepV.size(); ++i) {
            uint32_t x = (uint32_t)keepV[i];
            stampArr[x] = cur;
            if (isQueried[x]) present.push_back(x);
        }
        for (daf::Size i = 0; i < dropV.size(); ++i) {
            uint32_t x = (uint32_t)dropV[i];
            stampArr[x] = cur;
            if (isQueried[x]) present.push_back(x);
        }
        if (present.size() < 2) return;
        for (uint32_t qv : present) {
            for (auto &pr : pairsOf[qv]) {
                if (stampArr[pr.first] == cur) {      // both endpoints in this leaf
                    ++leafCount[pr.second];
                    sigmaL[pr.second] += leafSize;
                }
            }
        }
    };

    // forwarding-reference template deduces OnLeafFn from an rvalue; move the
    // named lambda in (its captures outlive the call in this scope).
    size_t emitted = SDCT_Augmented_NoTree(g, s, s, std::move(onLeaf));

    // ---- report ----
    long long zeros = 0;
    for (int p = 0; p < K; ++p) if (leafCount[p] == 0) ++zeros;

    std::printf("== E2 %s  s=%d  pairs=%d  seed=%llu ==\n",
                fpath, s, K, (unsigned long long)seed);
    std::printf("n=%zu  m_unique=%zu  total_level-s_leaves(emitted)=%zu  walk_leaves=%lld\n",
                n, edges.size(), emitted, totalLeaves);
    std::printf("pairs_with_zero_shared_leaves=%lld (%.1f%%)\n",
                zeros, 100.0 * (double)zeros / (double)K);
    std::printf("|leaves(u)∩leaves(v)| : median=%.1f p90=%.1f p99=%.1f max=%.0f\n",
                pct(leafCount, 0.5), pct(leafCount, 0.90), pct(leafCount, 0.99),
                pct(leafCount, 1.0));
    std::printf("Σ|L| (touched leaves) : median=%.1f p90=%.1f p99=%.1f max=%.0f\n",
                pct(sigmaL, 0.5), pct(sigmaL, 0.90), pct(sigmaL, 0.99),
                pct(sigmaL, 1.0));
    return 0;
}
