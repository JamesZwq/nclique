// test_region_forbidden.cpp — M1 adversarial test for the region-native
// PEEL foundation: a region (a maximal clique) maps to a degenerate
// ccpath::CCPath (no pivot holds: h=0, ell=0, u=n, T=s), and the s-clique
// support of a tuple AFTER some tuples are peeled is
//   ccpath::support_count(region_path, b)   with the peeled tuples inserted
//   into the path's forbidden antichain.
// We validate this against VERTEX-LEVEL brute force using SINGLETON classes
// (each vertex its own class), so a composition vector is exactly a vertex
// set and the answer is unambiguous. This pins the mapping + the
// forbidden-antichain inclusion-exclusion before any peel is built on top.
//
// Build: g++ -O3 -std=c++17 -I../src/NucleusDecomposition \
//          -o test_region_forbidden test_region_forbidden.cpp
#include "CCPathCore.h"
#include <cstdio>
#include <random>
#include <vector>
#include <algorithm>
using namespace std;
using ccpath::Vec;

// Pascal table C(n,k) as double, n up to 64.
static double PASCAL[65][65];
static void initPascal() {
    for (int n = 0; n <= 64; n++) {
        PASCAL[n][0] = 1.0;
        for (int k = 1; k <= n; k++)
            PASCAL[n][k] = PASCAL[n-1][k-1] + (k <= n-1 ? PASCAL[n-1][k] : 0.0);
    }
}
static double nCr_fn(int n, int k) {
    if (k < 0 || n < 0 || k > n) return 0.0;
    return PASCAL[n][k];
}

// VERTEX-LEVEL brute force: region = complete clique on N vertices
// (singleton classes). bpos = vertex set of the tuple. forb = list of
// peeled vertex sets. Count s-subsets Y of [N] with bpos subset Y and,
// for every f in forb, NOT (f subset Y).
static long long brute(int N, const vector<int>& bpos,
                       const vector<vector<int>>& forb, int s) {
    vector<int> idx(s);
    long long cnt = 0;
    // iterate all s-subsets via combination indices
    vector<int> comb(s);
    for (int i = 0; i < s; i++) comb[i] = i;
    auto contains = [](const vector<int>& Y, const vector<int>& sub) {
        for (int x : sub) if (!binary_search(Y.begin(), Y.end(), x)) return false;
        return true;
    };
    if (s > N) return 0;
    while (true) {
        // comb is sorted ascending, the current s-subset
        if (contains(comb, bpos)) {
            bool dead = false;
            for (auto& f : forb) if (contains(comb, f)) { dead = true; break; }
            if (!dead) cnt++;
        }
        // next combination
        int i = s - 1;
        while (i >= 0 && comb[i] == N - s + i) i--;
        if (i < 0) break;
        comb[i]++;
        for (int j = i + 1; j < s; j++) comb[j] = comb[j-1] + 1;
    }
    return cnt;
}

// ENGINE: region as CCPath (singleton classes => n = all ones, M = N),
// tuple b and forbidden as 0/1 composition vectors over the N classes.
static double engine(int N, const vector<int>& bpos,
                     const vector<vector<int>>& forb, int s) {
    Vec n((size_t)N, 1), zeros((size_t)N, 0);
    ccpath::CCPath p = ccpath::CCPath::initial(zeros, n, s);  // h=0,ell=0,u=n,T=s
    Vec b((size_t)N, 0);
    for (int x : bpos) b[(size_t)x] = 1;
    for (auto& f : forb) {
        Vec fv((size_t)N, 0);
        for (int x : f) fv[(size_t)x] = 1;
        ccpath::insert_antichain(p.forbidden, fv);
    }
    return ccpath::support_count(p, b, nCr_fn);
}

int main() {
    initPascal();
    // ---- fixed hand cases ----
    // K6, r=2 s=4, tuple {0,1}, peel {1,2}: expect 6 - 3 = 3.
    {
        double e = engine(6, {0,1}, {{1,2}}, 4);
        long long b = brute(6, {0,1}, {{1,2}}, 4);
        printf("hand K6 r2 s4 peel{1,2}: engine=%.0f brute=%lld %s\n",
               e, b, (e==b?"OK":"FAIL"));
    }
    // K6, r=2 s=4, tuple {0,1}, no peel: expect C(4,2)=6.
    {
        double e = engine(6, {0,1}, {}, 4);
        long long b = brute(6, {0,1}, {}, 4);
        printf("hand K6 r2 s4 nopeel:   engine=%.0f brute=%lld %s\n",
               e, b, (e==b?"OK":"FAIL"));
    }
    // ---- randomized adversarial trials ----
    mt19937 rng(20260618);
    int trials = 4000, fails = 0;
    for (int t = 0; t < trials; t++) {
        int N = 5 + (int)(rng() % 6);          // 5..10
        int r = 1 + (int)(rng() % 3);          // 1..3
        int s = r + 1 + (int)(rng() % 3);      // r+1..r+3
        if (s > N) continue;
        // random tuple b: r distinct vertices
        vector<int> perm(N); for (int i = 0; i < N; i++) perm[i] = i;
        shuffle(perm.begin(), perm.end(), rng);
        vector<int> bpos(perm.begin(), perm.begin() + r);
        sort(bpos.begin(), bpos.end());
        // 0..3 random forbidden tuples (each r distinct vertices)
        int nf = (int)(rng() % 4);
        vector<vector<int>> forb;
        for (int k = 0; k < nf; k++) {
            shuffle(perm.begin(), perm.end(), rng);
            vector<int> f(perm.begin(), perm.begin() + r);
            sort(f.begin(), f.end());
            forb.push_back(f);
        }
        double e = engine(N, bpos, forb, s);
        long long b = brute(N, bpos, forb, s);
        if ((long long)(e + 0.5) != b) {
            fails++;
            if (fails <= 8)
                printf("FAIL N=%d r=%d s=%d engine=%.1f brute=%lld\n",
                       N, r, s, e, b);
        }
    }
    printf("random trials: %d, fails=%d  %s\n",
           trials, fails, (fails==0?"[ALL EXACT]":"[MISMATCH]"));
    return fails == 0 ? 0 : 1;
}
