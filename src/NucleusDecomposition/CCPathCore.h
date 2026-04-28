// CCPath core algorithms — self-contained, no SDCT dependencies.
// Faithful C++ port of the algorithms in solved.py:
//   * insert_antichain
//   * tuple_to_threshold, impossible, covers_whole_path
//   * count_with_extra_lower
//   * inclusion_exclusion_terms
//   * support_count
//   * lazy_delete_tuple, first_failing_split_by_vector
//   * choose_split_vector, controlled_split
//
// Verified by 2000 random tests in solved.py against brute-force
// enumeration. This header is portable: only depends on <vector>,
// <algorithm>, <utility>, <cstdint>, <cstddef>.

#ifndef CCPATH_CORE_H
#define CCPATH_CORE_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace ccpath {

using Vec = std::vector<int16_t>;

// All entry points needing binomial coefficients take a `double nCr_fn(int,int)`
// function pointer (or template equivalent). Callers provide a wrapper around
// whatever nCr table they have. Keeps this header free of global state.

// Component-wise <=.
inline bool leq(const Vec &a, const Vec &b) {
    const int m = (int)a.size();
    for (int i = 0; i < m; ++i) if (a[i] > b[i]) return false;
    return true;
}

// Component-wise max into out.
inline void component_max(const Vec &a, const Vec &b, Vec &out) {
    const int m = (int)a.size();
    out.resize(m);
    for (int i = 0; i < m; ++i) out[i] = a[i] > b[i] ? a[i] : b[i];
}

inline Vec zeros_vec(int m) { return Vec((size_t)m, 0); }

// Insert a forbidden vector into an antichain, keeping only Pareto-minimal
// vectors. If some existing p <= a, a is dominated and dropped. Otherwise
// any q with a <= q is removed before appending a.
inline void insert_antichain(std::vector<Vec> &A, const Vec &a) {
    for (const auto &p : A) {
        if (leq(p, a)) return;
    }
    auto it = std::remove_if(A.begin(), A.end(),
                             [&a](const Vec &q) { return leq(a, q); });
    A.erase(it, A.end());
    A.push_back(a);
}

struct CCPath {
    Vec h;                       // hold count per class
    Vec n;                       // pivot count per class
    int T = 0;                   // sum y[c] target
    Vec ell;                     // lower bound per class
    Vec u;                       // upper bound per class
    std::vector<Vec> forbidden;  // antichain of forbidden thresholds

    // Optional metadata — not used by the core algorithms but preserved
    // through copies/moves so callers can carry application-specific
    // payloads (e.g., which tuples are still represented by this path,
    // which graph-level class each path-class corresponds to).
    std::vector<uint32_t> tupleIdxs;
    std::vector<int32_t> classIds;

    int m() const { return (int)n.size(); }

    bool quick_feasible() const {
        const int M = m();
        int sumL = 0, sumU = 0;
        for (int i = 0; i < M; ++i) {
            if (ell[i] > u[i]) return false;
            sumL += ell[i];
            sumU += u[i];
        }
        return sumL <= T && T <= sumU;
    }

    static CCPath initial(const Vec &h, const Vec &n, int s) {
        CCPath p;
        p.h = h;
        p.n = n;
        int sumH = 0;
        for (auto x : h) sumH += (int)x;
        p.T = s - sumH;
        p.ell.assign(n.size(), 0);
        p.u = n;
        return p;
    }
};

// a[c] = max(0, jvec[c] - h[c])
inline Vec tuple_to_threshold(const CCPath &p, const Vec &jvec) {
    const int M = p.m();
    Vec a(M);
    for (int c = 0; c < M; ++c) {
        int v = (int)jvec[c] - (int)p.h[c];
        a[c] = v > 0 ? (int16_t)v : (int16_t)0;
    }
    return a;
}

inline bool impossible(const CCPath &p, const Vec &a) {
    const int M = p.m();
    for (int i = 0; i < M; ++i) if (p.u[i] < a[i]) return true;
    return false;
}

inline bool covers_whole_path(const CCPath &p, const Vec &a) {
    const int M = p.m();
    for (int i = 0; i < M; ++i) if (p.ell[i] < a[i]) return false;
    return true;
}

// count_with_extra_lower — uses path.n[c] in the binomial weight.
// nCr_fn(n_minus_b, y_minus_b) returns C(n - b, y - b).
template <typename NCr>
inline double count_with_extra_lower_impl(const CCPath &p, const Vec &b,
                                          const Vec &extra_lower, NCr nCr_fn) {
    const int M = p.m();
    if ((int)b.size() != M) return 0.0;
    if ((int)extra_lower.size() != M) return 0.0;

    std::vector<int> L(M), U(M);
    for (int c = 0; c < M; ++c) {
        if (b[c] < 0 || (int)b[c] > (int)p.n[c]) return 0.0;
        int Lc = p.ell[c];
        if ((int)b[c] > Lc) Lc = (int)b[c];
        if ((int)extra_lower[c] > Lc) Lc = (int)extra_lower[c];
        L[c] = Lc;
        U[c] = (int)p.u[c];
        if (L[c] > U[c]) return 0.0;
    }
    int sumL = 0, sumU = 0;
    for (int c = 0; c < M; ++c) { sumL += L[c]; sumU += U[c]; }
    if (sumL > p.T || sumU < p.T) return 0.0;

    std::vector<double> dp((size_t)p.T + 1, 0.0);
    std::vector<double> ndp((size_t)p.T + 1, 0.0);
    dp[0] = 1.0;

    for (int c = 0; c < M; ++c) {
        std::fill(ndp.begin(), ndp.end(), 0.0);
        const int bc = (int)b[c];
        const int nc = (int)p.n[c];
        for (int total = 0; total <= p.T; ++total) {
            double w = dp[(size_t)total];
            if (w == 0.0) continue;
            int max_y = U[c];
            if (p.T - total < max_y) max_y = p.T - total;
            for (int y = L[c]; y <= max_y; ++y) {
                ndp[(size_t)(total + y)] += w * nCr_fn(nc - bc, y - bc);
            }
        }
        dp.swap(ndp);
    }
    return dp[(size_t)p.T];
}

// inclusion-exclusion terms over forbidden antichain. Returns (lower,coeff)
// pairs where lower is component-wise max over a subset of forbidden and
// coeff is the signed contribution.
inline std::vector<std::pair<Vec, int>>
inclusion_exclusion_terms(const std::vector<Vec> &A, int m_dim) {
    std::vector<std::pair<Vec, int>> terms;
    terms.emplace_back(zeros_vec(m_dim), 1);

    Vec mlower((size_t)m_dim);
    for (const auto &a : A) {
        // Snapshot existing terms.
        std::vector<std::pair<Vec, int>> snapshot = terms;
        for (auto &kv : snapshot) {
            component_max(kv.first, a, mlower);
            // Find existing entry with this lower.
            bool merged = false;
            for (auto &t : terms) {
                if (t.first == mlower) {
                    t.second -= kv.second;
                    merged = true;
                    break;
                }
            }
            if (!merged) terms.emplace_back(mlower, -kv.second);
        }
        // Remove zero-coefficient entries.
        terms.erase(std::remove_if(terms.begin(), terms.end(),
                                   [](const std::pair<Vec, int> &kv) {
                                       return kv.second == 0;
                                   }),
                    terms.end());
    }
    return terms;
}

// support_count — uses given nCr function for the weight binomial.
template <typename NCr>
inline double support_count_impl(const CCPath &p, const Vec &b, NCr nCr_fn) {
    for (const auto &a : p.forbidden) {
        if (leq(a, b)) return 0.0;
    }
    auto terms = inclusion_exclusion_terms(p.forbidden, p.m());
    double total = 0.0;
    for (auto &kv : terms) {
        total += (double)kv.second
               * count_with_extra_lower_impl(p, b, kv.first, nCr_fn);
    }
    return total;
}

inline std::vector<CCPath> lazy_delete_tuple(const CCPath &path, const Vec &jvec) {
    Vec a = tuple_to_threshold(path, jvec);
    if (impossible(path, a)) return { path };
    if (covers_whole_path(path, a)) return {};
    CCPath out = path;
    insert_antichain(out.forbidden, a);
    return { out };
}

inline std::vector<CCPath>
first_failing_split_by_vector(const CCPath &path, const Vec &a) {
    if (impossible(path, a)) return { path };
    if (covers_whole_path(path, a)) return {};

    const int M = path.m();
    std::vector<int> D;
    D.reserve((size_t)M);
    for (int c = 0; c < M; ++c) if (path.ell[c] < a[c]) D.push_back(c);

    std::vector<CCPath> children;
    children.reserve(D.size());

    for (size_t i = 0; i < D.size(); ++i) {
        int d = D[i];
        CCPath child = path;
        child.forbidden.clear();
        for (size_t k = 0; k < i; ++k) {
            int prev = D[k];
            if (a[prev] > child.ell[prev]) child.ell[prev] = a[prev];
        }
        int cap = (int)a[d] - 1;
        if (cap < (int)child.u[d]) child.u[d] = (int16_t)cap;
        if (child.quick_feasible()) children.push_back(std::move(child));
    }
    return children;
}

inline double split_score(const CCPath &path, const Vec &a) {
    double score = 1.0;
    bool active = false;
    const int M = path.m();
    for (int c = 0; c < M; ++c) {
        int ac = (int)a[c];
        int lc = (int)path.ell[c];
        int uc = (int)path.u[c];
        if (ac > lc) {
            if (ac > uc) return 0.0;
            active = true;
            int denom = uc - lc + 1;
            if (denom < 1) denom = 1;
            score *= (double)(uc - ac + 1) / (double)denom;
        }
    }
    return active ? score : 1.0;
}

inline Vec choose_split_vector(const CCPath &path) {
    size_t bestIdx = 0;
    double bestScore = -1.0;
    for (size_t i = 0; i < path.forbidden.size(); ++i) {
        double s = split_score(path, path.forbidden[i]);
        if (s > bestScore) { bestScore = s; bestIdx = i; }
    }
    return path.forbidden[bestIdx];
}

inline std::vector<CCPath> controlled_split(const CCPath &path, int kmax) {
    if ((int)path.forbidden.size() <= kmax) return { path };

    Vec a_star = choose_split_vector(path);
    auto raw = first_failing_split_by_vector(path, a_star);
    std::vector<CCPath> result;
    result.reserve(raw.size() * 2);

    for (auto &child : raw) {
        bool discard = false;
        for (const auto &q : path.forbidden) {
            if (q == a_star) continue;
            if (impossible(child, q)) continue;
            if (covers_whole_path(child, q)) { discard = true; break; }
            insert_antichain(child.forbidden, q);
        }
        if (discard) continue;
        if (!child.quick_feasible()) continue;
        auto sub = controlled_split(child, kmax);
        for (auto &c : sub) result.push_back(std::move(c));
    }
    return result;
}

// Convenience overloads using std::function-style nCr.
// Caller can pass a lambda computing C(n, k).
using NCrFn = double (*)(int, int);

inline double count_with_extra_lower(const CCPath &p, const Vec &b,
                                     const Vec &extra_lower, NCrFn nCr_fn) {
    return count_with_extra_lower_impl(p, b, extra_lower, nCr_fn);
}

inline double support_count(const CCPath &p, const Vec &b, NCrFn nCr_fn) {
    return support_count_impl(p, b, nCr_fn);
}

}  // namespace ccpath

#endif  // CCPATH_CORE_H
