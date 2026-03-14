// 测试 SDCT, SDCT_Par, SDCT_Par2 在不同线程数下的运行速度和正确性

#include <iostream>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
#include <cstring>
#include "src/graph/Graph.h"
#include "src/tree/MultiBranchTree.h"
#include "src/degeneracy_algorithm_cliques_V.h"
#include "src/misc.h"
#include "src/Global/Global.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

struct TestResult {
    std::string name;
    long long timeMs;
    std::vector<double> cliqueCount;
    bool isCorrect;
};

std::pair<long long, DynamicGraph<TreeGraphNode>> runSDCT(const char* fpath) {
    Graph g(fpath);
    g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    return {ms, tree};
}

std::pair<long long, DynamicGraph<TreeGraphNode>> runSDCT_Par(const char* fpath, int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    Graph g(fpath);
    g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    return {ms, tree};
}

std::pair<long long, DynamicGraph<TreeGraphNode>> runSDCT_Par2(const char* fpath, int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    Graph g(fpath);
    g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par2(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    return {ms, tree};
}

std::pair<long long, DynamicGraph<TreeGraphNode>> runSDCT_Par3(const char* fpath, int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    Graph g(fpath);
    g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par3(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    return {ms, tree};
}

std::pair<long long, DynamicGraph<TreeGraphNode>> runSDCT_Par4(const char* fpath, int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    Graph g(fpath); g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par4(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    return {ms, tree};
}

std::pair<long long, DynamicGraph<TreeGraphNode>> runSDCT_Par5(const char* fpath, int threads) {
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    Graph g(fpath); g.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> tree = SDCT_Par5(g, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    long long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    return {ms, tree};
}

bool compareCliqueCount(const std::vector<double>& ref, const std::vector<double>& test) {
    if (ref.size() != test.size()) return false;
    for (size_t i = 0; i < ref.size(); i++) {
        if (std::abs(ref[i] - test[i]) > 1e-6) return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <graph_file> <max_threads>" << std::endl;
        return 1;
    }

    const char* fpath = argv[1];
    int maxThreads = std::stoi(argv[2]);

    populate_nCr();
    daf::vListMap.resize(100000);
    std::memset(daf::vListMap.data(), -1, 100000 * sizeof(daf::Size));

    std::cout << "\n========================================" << std::endl;
    std::cout << "SDCT Correctness and Performance Test" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Graph: " << fpath << std::endl;
    std::cout << "Max threads: " << maxThreads << std::endl << std::endl;

    // Run SDCT baseline (serial) and get reference clique count
    std::cout << "Running SDCT (reference)..." << std::endl;
    auto [serialTime, refTree] = runSDCT(fpath);
    auto refCount = refTree.cliqueCount();
    std::vector<double> refVec(refCount.data(), refCount.data() + refCount.c_size);
    
    std::cout << "SDCT: time=" << serialTime << "ms" << std::endl;
    std::cout << "Reference clique counts (first 15): ";
    for (size_t i = 0; i < std::min((size_t)15, refVec.size()); i++) {
        if (refVec[i] > 0) std::cout << "k" << i << "=" << (long long)refVec[i] << " ";
    }
    std::cout << std::endl << std::endl;

    std::vector<TestResult> allResults;
    allResults.push_back({"SDCT", serialTime, refVec, true});

    // Test all parallel versions with different thread counts
    std::vector<std::string> versions = {"SDCT_Par", "SDCT_Par2", "SDCT_Par3", "SDCT_Par4", "SDCT_Par5"};
    
    for (int t = 1; t <= maxThreads; t *= 2) {
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Testing with " << t << " thread(s)" << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        
        for (const auto& version : versions) {
            std::cout << version << " (threads=" << t << ")... " << std::flush;
            
            long long time;
            DynamicGraph<TreeGraphNode> tree;
            
            if (version == "SDCT_Par") {
                auto result = runSDCT_Par(fpath, t);
                time = result.first;
                tree = result.second;
            } else if (version == "SDCT_Par2") {
                auto result = runSDCT_Par2(fpath, t);
                time = result.first;
                tree = result.second;
            } else if (version == "SDCT_Par3") {
                auto result = runSDCT_Par3(fpath, t);
                time = result.first;
                tree = result.second;
            } else if (version == "SDCT_Par4") {
                auto result = runSDCT_Par4(fpath, t);
                time = result.first;
                tree = result.second;
            } else if (version == "SDCT_Par5") {
                auto result = runSDCT_Par5(fpath, t);
                time = result.first;
                tree = result.second;
            }
            
            auto count = tree.cliqueCount();
            std::vector<double> countVec(count.data(), count.data() + count.c_size);
            bool correct = compareCliqueCount(refVec, countVec);
            
            double speedup = (time > 0) ? (double)serialTime / time : 0.0;
            
            std::cout << "time=" << time << "ms, speedup=" << std::fixed << std::setprecision(2) 
                      << speedup << "x, ";
            
            if (correct) {
                std::cout << "✓ CORRECT" << std::endl;
            } else {
                std::cout << "✗ INCORRECT" << std::endl;
                // Print first mismatch
                for (size_t i = 0; i < std::min(refVec.size(), countVec.size()); i++) {
                    if (std::abs(refVec[i] - countVec[i]) > 1e-6) {
                        std::cout << "  First mismatch at k=" << i << ": expected=" 
                                  << (long long)refVec[i] << ", got=" << (long long)countVec[i] << std::endl;
                        break;
                    }
                }
            }
            
            allResults.push_back({version + "_t" + std::to_string(t), time, countVec, correct});
        }
        std::cout << std::endl;
    }

    // Summary
    std::cout << "\n========================================" << std::endl;
    std::cout << "SUMMARY" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::left << std::setw(20) << "Implementation" 
              << std::right << std::setw(12) << "Time (ms)" 
              << std::setw(12) << "Speedup" 
              << std::setw(15) << "Status" << std::endl;
    std::cout << std::string(59, '-') << std::endl;
    
    for (const auto& result : allResults) {
        double speedup = (result.timeMs > 0) ? (double)serialTime / result.timeMs : 0.0;
        std::cout << std::left << std::setw(20) << result.name
                  << std::right << std::setw(12) << result.timeMs
                  << std::setw(11) << std::fixed << std::setprecision(2) << speedup << "x"
                  << std::setw(15) << (result.isCorrect ? "✓ PASS" : "✗ FAIL")
                  << std::endl;
    }
    std::cout << std::string(59, '=') << std::endl;

    // Check if all are correct
    bool allCorrect = true;
    for (const auto& result : allResults) {
        if (!result.isCorrect) {
            allCorrect = false;
            break;
        }
    }

    if (allCorrect) {
        std::cout << "\n✓ ALL TESTS PASSED - All implementations produce identical cliqueCount!" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ SOME TESTS FAILED - Some implementations produce different results!" << std::endl;
        return 1;
    }
}
