// 测试高级并行算法
// 目标：在8线程上达到3倍以上加速

#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include "src/graph/Graph.h"
#include "src/tree/MultiBranchTree.h"
#include "src/NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "src/degeneracy_algorithm_cliques_V.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// 验证两个结果是否相同
bool verifyResults(
    const std::vector<std::pair<std::vector<daf::Size>, int>> &result1,
    const std::vector<std::pair<std::vector<daf::Size>, int>> &result2) {
    
    if (result1.size() != result2.size()) {
        std::cout << "ERROR: Result sizes differ: " << result1.size() 
                  << " vs " << result2.size() << std::endl;
        return false;
    }
    
    // 创建映射：clique -> core value
    std::map<std::vector<daf::Size>, int> map1, map2;
    
    for (const auto &[clique, core] : result1) {
        auto sorted = clique;
        std::sort(sorted.begin(), sorted.end());
        map1[sorted] = core;
    }
    
    for (const auto &[clique, core] : result2) {
        auto sorted = clique;
        std::sort(sorted.begin(), sorted.end());
        map2[sorted] = core;
    }
    
    // 比较
    int mismatches = 0;
    for (const auto &[clique, core1] : map1) {
        auto it = map2.find(clique);
        if (it == map2.end()) {
            std::cout << "ERROR: Clique not found in result2" << std::endl;
            mismatches++;
            if (mismatches >= 10) break;
        } else if (it->second != core1) {
            std::cout << "ERROR: Core value mismatch for clique: " 
                      << core1 << " vs " << it->second << std::endl;
            mismatches++;
            if (mismatches >= 10) break;
        }
    }
    
    if (mismatches > 0) {
        std::cout << "Total mismatches: " << mismatches << std::endl;
        return false;
    }
    
    std::cout << "✓ Results match perfectly!" << std::endl;
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <graph_file> <r> <s> [num_threads]" << std::endl;
        std::cout << "Example: " << argv[0] << " data/com-dblp.ungraph.txt 3 4 8" << std::endl;
        return 1;
    }
    
    const char* graphFile = argv[1];
    int r = std::stoi(argv[2]);
    int s = std::stoi(argv[3]);
    int numThreads = (argc >= 5) ? std::stoi(argv[4]) : 8;
    
#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif
    
    std::cout << "========================================" << std::endl;
    std::cout << "Advanced Parallel Nucleus Decomposition Test" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Graph: " << graphFile << std::endl;
    std::cout << "Parameters: r=" << r << ", s=" << s << std::endl;
    std::cout << "Threads: " << numThreads << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // 加载图
    std::cout << "Loading graph..." << std::endl;
    Graph edgeGraph(graphFile);
    edgeGraph.printGraphInfo();
    edgeGraph.sortByDegeneracyOrder();
    
    // 构建tree
    std::cout << "\nBuilding clique tree..." << std::endl;
    auto treeGraph = SDCT_Par(edgeGraph, 1000000, 0);
    std::cout << "Tree leaves: " << treeGraph.adj_list.size() << std::endl;
    
    // 构建treeGraphV
    std::cout << "Building treeGraphV..." << std::endl;
    DynamicGraphSet<TreeGraphNode> treeGraphV(treeGraph, edgeGraph.getGraphNodeSize(), s);
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test 1: Reference Algorithm (Single-threaded)" << std::endl;
    std::cout << "========================================" << std::endl;
    
    auto refTree = treeGraph.clone();
    auto refTreeGraphV = treeGraphV.clone();
    
    auto start = std::chrono::high_resolution_clock::now();
    auto refResult = NucleusCoreDecompositionRClique(refTree, edgeGraph, refTreeGraphV, r, s);
    auto end = std::chrono::high_resolution_clock::now();
    auto refTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    std::cout << "\nReference time: " << refTime << " ms" << std::endl;
    std::cout << "Result size: " << refResult.size() << std::endl;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test 2: Advanced Parallel Algorithm" << std::endl;
    std::cout << "========================================" << std::endl;
    
    auto testTree = treeGraph.clone();
    auto testTreeGraphV = treeGraphV.clone();
    
    start = std::chrono::high_resolution_clock::now();
    auto testResult = AdvancedParallel::NucleusCoreDecompositionAdvancedParallel(
        testTree, edgeGraph, testTreeGraphV, r, s);
    end = std::chrono::high_resolution_clock::now();
    auto testTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    std::cout << "\nAdvanced parallel time: " << testTime << " ms" << std::endl;
    std::cout << "Result size: " << testResult.size() << std::endl;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Performance Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Reference (1 thread):  " << refTime << " ms" << std::endl;
    std::cout << "Advanced (" << numThreads << " threads): " << testTime << " ms" << std::endl;
    std::cout << "Speedup: " << (double)refTime / testTime << "x" << std::endl;
    
    if ((double)refTime / testTime >= 3.0) {
        std::cout << "✓ TARGET ACHIEVED: 3x+ speedup!" << std::endl;
    } else {
        std::cout << "✗ Target not met (need 3x+)" << std::endl;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Correctness Verification" << std::endl;
    std::cout << "========================================" << std::endl;
    
    bool correct = verifyResults(refResult, testResult);
    
    if (correct) {
        std::cout << "\n✓ ALL TESTS PASSED!" << std::endl;
    } else {
        std::cout << "\n✗ CORRECTNESS CHECK FAILED!" << std::endl;
        return 1;
    }
    
    return 0;
}
