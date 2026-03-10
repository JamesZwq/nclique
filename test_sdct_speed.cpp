// 测试 SDCT 和 SDCT_Par 在不同线程数下的运行速度

#include <iostream>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <string>
#include "src/graph/Graph.h"
#include "src/tree/MultiBranchTree.h"
#include "src/degeneracy_algorithm_cliques_V.h"
#include "src/misc.h"
#include "src/Global/Global.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

// 运行 SDCT 并返回耗时 (ms)
long long runSDCT(const char* fpath) {
    Graph g(fpath);
    g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
}

// 运行 SDCT_Par 并返回耗时 (ms)
long long runSDCT_Par(const char* fpath, int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    Graph g(fpath);
    g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <graph_file> <max_threads>" << std::endl;
        std::cerr << "  Threads double from 1 up to max_threads." << std::endl;
        return 1;
    }

    const char* fpath = argv[1];
    int maxThreads = std::stoi(argv[2]);

    populate_nCr();

    std::cout << "graph=" << fpath << std::endl;

    // SDCT 串行（只跑一次，作为 baseline）
    long long serialTime = runSDCT(fpath);
    std::cout << "SDCT threads=1 time_ms=" << serialTime << std::endl;

    // SDCT_Par：线程从 1 翻倍直到 maxThreads
    for (int t = 1; t <= maxThreads; t *= 2) {
        long long parTime = runSDCT_Par(fpath, t);
        double speedup = (parTime > 0) ? (double)serialTime / parTime : 0.0;
        std::cout << "SDCT_Par threads=" << t
                  << " time_ms=" << parTime
                  << " speedup=" << std::fixed << std::setprecision(2) << speedup
                  << std::endl;
    }

    return 0;
}
