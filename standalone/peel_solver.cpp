// Standalone correct solver for the abstract peeling problem.
//
// Input file format:
//   Line 1: numClasses T  (T = s - h = total pivots needed per s-clique)
//   Line 2: np_1 np_2 ... np_m  (pivot count per class)
//   Line 3: numTuples
//   For each tuple (next numTuples lines):
//     r j_1 j_2 ... j_m  (r = tuple size, j_c = count from class c, most j_c = 0 omitted)
//     Format: r followed by pairs (classIdx, count), ended by -1
//     Example: 3 0 2 1 3  means r=3, j_{c0}=2, j_{c1}=1, j_{c2}=3... no that's bad.
//     Simpler format: r numEntries (classIdx count) ...
//
// Actually, let me use the simplest possible format:
//
//   Line 1: m T   (m = numClasses, T = target total pivots = s - h)
//   Line 2: np[0] np[1] ... np[m-1]
//   Line 3: K     (K = number of tuples)
//   Next K lines: each tuple as: r c1 c2 ... cr  (class indices, sorted, with repeats)
//     Example: "3 0 0 1" means r=3, tuple = (class0, class0, class1)
//
// Output: for each tuple, its core value after peeling.
//
// Algorithm: brute-force correct peeling.
//   1. Enumerate ALL s-cliques (all Q_s with sum = T, 0 <= n_c <= np[c])
//   2. For each tuple: compute support = number of s-cliques containing its r-clique
//   3. Peel in order of support (standard peeling)
//
// This is O(product(np[c]) * K) — only feasible for small inputs.
// But it serves as the ground truth.

#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <functional>

using namespace std;

// nCr table
double nCr[501][501];
void buildNCr() {
    for (int i = 0; i <= 500; i++) {
        nCr[i][0] = 1;
        for (int j = 1; j <= i; j++)
            nCr[i][j] = nCr[i-1][j-1] + nCr[i-1][j];
    }
}

struct Tuple {
    int r;
    vector<int> key;            // sorted class indices (with repeats)
    map<int, int> classCounts;  // class -> count
    long long mult;             // number of r-cliques represented
    double support;
    int coreValue;
    bool peeled;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Usage: %s input.txt\n", argv[0]);
        return 1;
    }
    buildNCr();

    FILE* f = fopen(argv[1], "r");
    if (!f) { printf("Cannot open %s\n", argv[1]); return 1; }

    int m, T;
    fscanf(f, "%d %d", &m, &T);

    vector<int> np(m);
    for (int i = 0; i < m; i++) fscanf(f, "%d", &np[i]);

    int K;
    fscanf(f, "%d", &K);

    vector<Tuple> tuples(K);
    for (int ti = 0; ti < K; ti++) {
        fscanf(f, "%d", &tuples[ti].r);
        tuples[ti].key.resize(tuples[ti].r);
        for (int j = 0; j < tuples[ti].r; j++) {
            fscanf(f, "%d", &tuples[ti].key[j]);
        }
        sort(tuples[ti].key.begin(), tuples[ti].key.end());
        // Build classCounts
        for (int c : tuples[ti].key) tuples[ti].classCounts[c]++;
        // Compute mult = product of C(np[c], j_c)
        tuples[ti].mult = 1;
        for (auto& [c, jc] : tuples[ti].classCounts) {
            tuples[ti].mult *= (long long)llround(nCr[np[c]][jc]);
        }
        tuples[ti].peeled = false;
        tuples[ti].coreValue = -1;
    }
    fclose(f);

    printf("Classes: %d, Target T: %d, Tuples: %d\n", m, T, K);

    // === Enumerate ALL s-cliques ===
    // s-clique = vector of (n_c for each class), sum = T, 0 <= n_c <= np[c]
    // For each s-clique: determine which tuples' r-cliques it contains.

    // An s-clique Q contains r-clique from tuple τ iff:
    //   for each class c: n_Q[c] >= j_c (enough vertices from each class)
    //   AND the r-clique vertices are a subset of Q's vertices.
    //
    // Number of r-cliques from τ contained in Q:
    //   product over c of C(n_Q[c], j_c)
    //
    // For support: support(τ) = Σ_Q (1 if Q contains at least one r-clique from τ)
    //            = Σ_Q (1 if ∀c: n_Q[c] >= j_c)
    //
    // Wait, that counts s-cliques that COULD contain an r-clique from τ.
    // But support(τ) should be per-specific-r-clique.
    // support(τ) = per-C_r' count = (Σ_Q #r-cliques-from-τ-in-Q) / mult(τ)

    // Enumerate s-cliques by recursion over classes
    struct SClique {
        vector<int> nQ; // nQ[c] = pivots from class c
    };

    vector<SClique> allSCliques;
    {
        SClique cur;
        cur.nQ.resize(m, 0);
        function<void(int, int)> enumerate = [&](int classIdx, int remaining) {
            if (classIdx == m) {
                if (remaining == 0) allSCliques.push_back(cur);
                return;
            }
            for (int nc = 0; nc <= min(np[classIdx], remaining); nc++) {
                cur.nQ[classIdx] = nc;
                enumerate(classIdx + 1, remaining - nc);
            }
        };
        enumerate(0, T);
    }

    printf("Total s-cliques: %zu\n", allSCliques.size());

    // === Compute initial support ===
    // For each tuple τ: support = Σ_Q [Π_c C(nQ[c], j_c)] / mult(τ)
    for (int ti = 0; ti < K; ti++) {
        double total = 0;
        for (auto& Q : allSCliques) {
            double ways = 1;
            for (auto& [c, jc] : tuples[ti].classCounts) {
                if (Q.nQ[c] < jc) { ways = 0; break; }
                ways *= nCr[Q.nQ[c]][jc];
            }
            total += ways;
        }
        tuples[ti].support = total / tuples[ti].mult;
    }

    printf("\nInitial support:\n");
    for (int ti = 0; ti < K; ti++) {
        printf("  τ%d: support=%.0f mult=%lld key=(",
               ti, tuples[ti].support, tuples[ti].mult);
        for (int j = 0; j < (int)tuples[ti].key.size(); j++)
            printf("%s%d", j?",":"", tuples[ti].key[j]);
        printf(")\n");
    }

    // === Track alive s-cliques ===
    vector<bool> alive(allSCliques.size(), true);

    // === Peeling ===
    int coreLevel = 0;
    map<int, long long> coreDist;

    for (int step = 0; step < K; step++) {
        // Find min-support alive tuple
        int minIdx = -1;
        double minSup = 1e18;
        for (int ti = 0; ti < K; ti++) {
            if (!tuples[ti].peeled && tuples[ti].support < minSup) {
                minSup = tuples[ti].support;
                minIdx = ti;
            }
        }
        if (minIdx == -1) break;

        // Peel it
        tuples[minIdx].peeled = true;
        int level = (int)llround(tuples[minIdx].support);
        coreLevel = max(coreLevel, level);
        tuples[minIdx].coreValue = coreLevel;
        coreDist[coreLevel] += tuples[minIdx].mult;

        // Kill s-cliques containing this tuple's r-cliques
        // An s-clique Q contains τ's r-clique iff ∀c: nQ[c] >= m_c
        auto& tau = tuples[minIdx];
        for (int qi = 0; qi < (int)allSCliques.size(); qi++) {
            if (!alive[qi]) continue;
            bool contains = true;
            for (auto& [c, jc] : tau.classCounts) {
                if (allSCliques[qi].nQ[c] < jc) { contains = false; break; }
            }
            if (!contains) continue;

            // This s-clique dies
            alive[qi] = false;

            // Update support of all alive tuples
            for (int ti = 0; ti < K; ti++) {
                if (tuples[ti].peeled) continue;
                // How many r-cliques from τ_ti are in this s-clique?
                double ways = 1;
                for (auto& [c2, jc2] : tuples[ti].classCounts) {
                    if (allSCliques[qi].nQ[c2] < jc2) { ways = 0; break; }
                    ways *= nCr[allSCliques[qi].nQ[c2]][jc2];
                }
                tuples[ti].support -= ways / tuples[ti].mult;
            }
        }
    }

    // === Output ===
    printf("\n=== Core Decomposition ===\n");
    for (int ti = 0; ti < K; ti++) {
        printf("  τ%d: core=%d  mult=%lld  key=(",
               ti, tuples[ti].coreValue, tuples[ti].mult);
        for (int j = 0; j < (int)tuples[ti].key.size(); j++)
            printf("%s%d", j?",":"", tuples[ti].key[j]);
        printf(")\n");
    }

    printf("\n=== Core Distribution ===\n");
    for (auto& [core, cnt] : coreDist) {
        printf("  core=%d count=%lld\n", core, cnt);
    }

    return 0;
}
