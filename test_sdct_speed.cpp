// 测试 SDCT, SDCT_Par, SDCT_Par2 在不同线程数下的运行速度

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

long long runSDCT(const char* fpath) {
    Graph g(fpath);
    g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
}

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

long long runSDCT_Par2(const char* fpath, int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    Graph g(fpath);
    g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par2(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
}

long long runSDCT_Par3(const char* fpath, int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    Graph g(fpath);
    g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par3(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
}

long long runSDCT_Par4(const char* fpath, int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    Graph g(fpath); g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par4(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <graph_file> <max_threads>" << std::endl;
        return 1;
    }

    const char* fpath = argv[1];
    int maxThreads = std::stoi(argv[2]);

    populate_nCr();

    std::cout << "graph=" << fpath << std::endl;

    // SDCT baseline (serial)
    long long serialTime = runSDCT(fpath);
    std::cout << "SDCT threads=1 time_ms=" << serialTime << std::endl;

    // SDCT_Par, SDCT_Par2, SDCT_Par3: threads double from 1 to maxThreads
    for (int t = 1; t <= maxThreads; t *= 2) {
        long long par1Time = runSDCT_Par(fpath, t);
        double speedup1 = (par1Time > 0) ? (double)serialTime / par1Time : 0.0;
        std::cout << "SDCT_Par threads=" << t
                  << " time_ms=" << par1Time
                  << " speedup=" << std::fixed << std::setprecision(2) << speedup1
                  << std::endl;

        long long par2Time = runSDCT_Par2(fpath, t);
        double speedup2 = (par2Time > 0) ? (double)serialTime / par2Time : 0.0;
        std::cout << "SDCT_Par2 threads=" << t
                  << " time_ms=" << par2Time
                  << " speedup=" << std::fixed << std::setprecision(2) << speedup2
                  << std::endl;

        long long par3Time = runSDCT_Par3(fpath, t);
        double speedup3 = (par3Time > 0) ? (double)serialTime / par3Time : 0.0;
        std::cout << "SDCT_Par3 threads=" << t
                  << " time_ms=" << par3Time
                  << " speedup=" << std::fixed << std::setprecision(2) << speedup3
                  << std::endl;

        long long par4Time = runSDCT_Par4(fpath, t);
        double speedup4 = (par4Time > 0) ? (double)serialTime / par4Time : 0.0;
        std::cout << "SDCT_Par4 threads=" << t
                  << " time_ms=" << par4Time
                  << " speedup=" << std::fixed << std::setprecision(2) << speedup4
                  << std::endl;
    }

    return 0;
}
