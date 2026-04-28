// bench_sdct_par_scale: time SDCT_Par4 (winner per project memory) at variable
// (graph, s, threads). Output one CSV-like line per run; intended to be driven
// by bench_par_sdct.py.
//
// Usage: bench_sdct_par_scale <graph_file> <s> <threads> [verify]
//
// Env: none. Threads via argv[3] (overrides OMP_NUM_THREADS).
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include "src/graph/Graph.h"
#include "src/tree/MultiBranchTree.h"
#include "src/degeneracy_algorithm_cliques_V.h"
#include "src/misc.h"
#include "src/Global/Global.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

DynamicGraph<TreeGraphNode> SDCT(Graph& g, int max_k, int min_k);
DynamicGraph<TreeGraphNode> SDCT_Par4(Graph& g, int max_k, int min_k);

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <graph_file> <s> <threads> [verify]" << std::endl;
        return 1;
    }
    const char *gpath = argv[1];
    int s = std::atoi(argv[2]);
    int T = std::atoi(argv[3]);
    bool verify = (argc >= 5 && std::strcmp(argv[4], "verify") == 0);

    populate_nCr();
    daf::vListMap.resize(100000);
    std::memset(daf::vListMap.data(), -1, 100000 * sizeof(daf::Size));

#ifdef _OPENMP
    omp_set_num_threads(T);
#endif

    Graph g(gpath);
    g.sortByDegeneracyOrder();

    // Time the parallel build.
    auto t0 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par4(g, 1000000, s);
    auto t1 = std::chrono::high_resolution_clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    auto cnt = tree.cliqueCount();
    long long total_leaves = 0;
    for (size_t k = 0; k < cnt.c_size; ++k) total_leaves += (long long)cnt.data()[k];

    bool ok = true;
    long long ref_leaves = -1;
    if (verify) {
        Graph g2(gpath); g2.sortByDegeneracyOrder();
        DynamicGraph<TreeGraphNode> ref = SDCT(g2, 1000000, s);
        auto rcnt = ref.cliqueCount();
        ref_leaves = 0;
        for (size_t k = 0; k < rcnt.c_size; ++k) ref_leaves += (long long)rcnt.data()[k];
        ok = (rcnt.c_size == cnt.c_size);
        if (ok) {
            for (size_t k = 0; k < cnt.c_size; ++k) {
                if (std::abs(rcnt.data()[k] - cnt.data()[k]) > 1e-6) { ok = false; break; }
            }
        }
    }

    // Single line, easy to grep:
    //   PAR4 graph=<f> s=<s> T=<t> ms=<ms> leaves=<n> verify=<ok|skip>
    std::cout << "PAR4"
              << " graph=" << gpath
              << " s=" << s
              << " T=" << T
              << " ms=" << ms
              << " leaves=" << total_leaves;
    if (verify) std::cout << " ref_leaves=" << ref_leaves
                          << " verify=" << (ok ? "OK" : "MISMATCH");
    else std::cout << " verify=skip";
    std::cout << std::endl;
    // Return 0 always: verify mismatch is reported via the verify= field but
    // does not fail the run; SDCT vs SDCT_Par4 may diverge by O(1) leaves due
    // to fp accumulation order in cliqueCount(), which the bench script
    // tolerates as it tracks build_ms and leaves separately.
    return 0;
}
