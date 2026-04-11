#pragma once
//
// MaxCliqEnum: lightweight maximal clique enumerator
// Degeneracy-ordered BK with P∪X pivoting.
// Returns vector<vector<daf::Size>> of maximal cliques.
//

#include "Global/Global.h"
#include <algorithm>
#include <vector>

namespace maxcliq {

// Sorted intersection of two sorted vectors
inline std::vector<int> intersect(const std::vector<int> &a, const std::vector<int> &b) {
    std::vector<int> res;
    auto i = a.begin(), j = b.begin();
    while (i != a.end() && j != b.end()) {
        if (*i < *j) ++i;
        else if (*i > *j) ++j;
        else { res.push_back(*i); ++i; ++j; }
    }
    return res;
}

// Get sorted neighbor list for vertex v
inline std::vector<int> getNeighbors(int v, const Graph &g) {
    auto [nb, ne] = g.getNbr(v);
    std::vector<int> res(ne - nb);
    for (int i = nb; i < ne; ++i) res[i - nb] = g.adj_list[i];
    std::sort(res.begin(), res.end());
    return res;
}

// Count |N(u) ∩ P| where P is sorted
inline int countInP(int u, const Graph &g, const std::vector<int> &P) {
    auto [nb, ne] = g.getNbr(u);
    int count = 0;
    auto pi = P.begin();
    for (int i = nb; i < ne && pi != P.end(); ++i) {
        int w = g.adj_list[i];
        while (pi != P.end() && *pi < w) ++pi;
        if (pi != P.end() && *pi == w) { count++; ++pi; }
    }
    return count;
}

// BK recursion with P∪X pivoting
// clique: current R (being built)
// P, X: sorted candidate and exclusion sets
static void bkRecurse(
    std::vector<daf::Size> &clique,
    std::vector<int> P, std::vector<int> X,
    const Graph &g,
    std::vector<std::vector<daf::Size>> &output)
{
    if (P.empty()) {
        if (X.empty() && !clique.empty()) {
            output.push_back(clique);
        }
        return;
    }

    // Choose pivot from P∪X maximizing |N(u) ∩ P|
    int pivot = P[0], bestCnt = -1;
    for (int u : P) {
        int c = countInP(u, g, P);
        if (c > bestCnt) { bestCnt = c; pivot = u; }
    }
    for (int u : X) {
        int c = countInP(u, g, P);
        if (c > bestCnt) { bestCnt = c; pivot = u; }
    }

    // Candidates = P \ N(pivot) (sorted)
    auto pivotNbrs = getNeighbors(pivot, g);
    std::vector<int> candidates;
    {
        auto pi = pivotNbrs.begin();
        for (int v : P) {
            while (pi != pivotNbrs.end() && *pi < v) ++pi;
            if (pi == pivotNbrs.end() || *pi != v) {
                candidates.push_back(v); // v is NOT a neighbor of pivot
            }
        }
    }

    for (int v : candidates) {
        clique.push_back(v);

        auto vNbrs = getNeighbors(v, g);
        auto newP = intersect(P, vNbrs);
        auto newX = intersect(X, vNbrs);

        bkRecurse(clique, std::move(newP), std::move(newX), g, output);

        clique.pop_back();

        // Move v from P to X
        P.erase(std::lower_bound(P.begin(), P.end(), v));
        X.insert(std::lower_bound(X.begin(), X.end(), v), v);
    }
}

} // namespace maxcliq

// Enumerate all maximal cliques in graph g (degeneracy-ordered)
// g must be degeneracy-sorted (adj_list ordered by degeneracy)
static std::vector<std::vector<daf::Size>> enumerateMaximalCliques(const Graph &g) {
    int n = g.getGraphNodeSize();
    std::vector<std::vector<daf::Size>> result;
    std::vector<bool> processed(n, false);

    for (int v = 0; v < n; ++v) {
        // P = later neighbors (not yet processed)
        // X = earlier neighbors (already processed)
        std::vector<int> P, X;
        auto [nb, ne] = g.getNbr(v);
        for (int i = nb; i < ne; ++i) {
            int u = g.adj_list[i];
            if (processed[u]) X.push_back(u);
            else P.push_back(u);
        }
        std::sort(P.begin(), P.end());
        std::sort(X.begin(), X.end());

        std::vector<daf::Size> clique = {(daf::Size)v};
        maxcliq::bkRecurse(clique, std::move(P), std::move(X), g, result);

        processed[v] = true;
    }

    // Sort each clique's vertex list
    for (auto &c : result) std::sort(c.begin(), c.end());

    return result;
}
