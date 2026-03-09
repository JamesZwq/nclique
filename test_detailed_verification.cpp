// 增强的测试程序 - 详细验证正确性

#include <iostream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <map>
#include "src/graph/Graph.h"
#include "src/tree/MultiBranchTree.h"
#include "src/NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "src/degeneracy_algorithm_cliques_V.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void printSeparator() {
    std::cout << "\n========================================" << std::endl;
}

bool verifyResults(
    const std::vector<std::pair<std::vector<daf::Size>, int>>& ref,
    const std::vector<std::pair<std::vector<daf::Size>, int>>& test,
    bool verbose = false) {
    
    if (ref.size() != test.size()) {
        std::cout << "✗ Size mismatch: ref=" << ref.size() << ", test=" << test.size() << std::endl;
        return false;
    }
    
    // 创建映射：clique -> core value
    std::map<std::vector<daf::Size>, int> refMap, testMap;
    
    for (const auto& [clique, core] : ref) {
        auto sortedClique = clique;
        std::sort(sortedClique.begin(), sortedClique.end());
        refMap[sortedClique] = core;
    }
    
    for (const auto& [clique, core] : test) {
        auto sortedClique = clique;
        std::sort(sortedClique.begin(), sortedClique.end());
        testMap[sortedClique] = core;
    }
    
    // 检查每个clique的core值
    int mismatches = 0;
    int checked = 0;
    
    for (const auto& [clique, refCore] : refMap) {
        checked++;
        
        auto it = testMap.find(clique);
        if (it == testMap.end()) {
            if (verbose) {
                std::cout << "✗ Clique not found in test results: ";
                for (auto v : clique) std::cout << v << " ";
                std::cout << std::endl;
            }
            mismatches++;
            if (mismatches > 10) break;
        } else if (it->second != refCore) {
            if (verbose) {
                std::cout << "✗ Core value mismatch for clique: ";
                for (auto v : clique) std::cout << v << " ";
                std::cout << " (ref=" << refCore << ", test=" << it->second << ")" << std::endl;
            }
            mismatches++;
            if (mismatches > 10) break;
        }
    }
    
    if (mismatches > 0) {
        std::cout << "✗ Found " << mismatches << " mismatches out of " << checked << " cliques" << std::endl;
        return false;
    }
    
    std::cout << "✓ All " << checked << " cliques verified correctly" << std::endl;
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <graph_file> <r> <s> [num_threads]" << std::endl;
        return 1;
    }
    
    const char* fpath = argv[1];
    int r = std::stoi(argv[2]);
    int s = std::stoi(argv[3]);
    
#ifdef _OPENMP
    int numThreads = (argc >= 5) ? std::stoi(argv[4]) : omp_get_max_threads();
    omp_set_num_threads(numThreads);
#else
    int numThreads = 1;
#endif
    
    printSeparator();
    std::cout << "Detailed Correctness Verification" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Graph: " << fpath << std::endl;
    std::cout << "Parameters: r=" << r << ", s=" << s << std::endl;
    std::cout << "Threads: " << numThreads << std::endl;
    printSeparator();
    
    // 加载图
    std::cout << "\n[1/4] Loading graph..." << std::endl;
    Graph edgeGraph(fpath);
    edgeGraph.printGraphInfo();
    edgeGraph.sortByDegeneracyOrder();
    
    // 构建 tree
    std::cout << "\n[2/4] Building tree structure..." << std::endl;
    auto treeGraph = SDCT_Par(edgeGraph, 1000000, 0);
    std::cout << "Tree leaves: " << treeGraph.adj_list.size() << std::endl;
    
    // 构建 treeGraphV
    std::cout << "\n[3/4] Building treeGraphV..." << std::endl;
    DynamicGraphSet<TreeGraphNode> treeGraphV(treeGraph, edgeGraph.getGraphNodeSize(), s);
    
    // 测试Reference算法
    printSeparator();
    std::cout << "[4/4] Running Reference Algorithm" << std::endl;
    printSeparator();
    
    auto refTree = treeGraph.clone();
    auto refTreeGraphV = treeGraphV.clone();
    
    auto refStart = std::chrono::high_resolution_clock::now();
    auto refResult = NucleusCoreDecompositionRCliqueRef(
        refTree, edgeGraph, refTreeGraphV, r, s);
    auto refEnd = std::chrono::high_resolution_clock::now();
    auto refTime = std::chrono::duration_cast<std::chrono::milliseconds>(refEnd - refStart).count();
    
    std::cout << "\nReference: " << refTime << " ms" << std::endl;
    std::cout << "Result size: " << refResult.size() << std::endl;
    
    // 测试Ultra Parallel算法
    printSeparator();
    std::cout << "Running Ultra Parallel Algorithm" << std::endl;
    printSeparator();
    
    auto ultraTree = treeGraph.clone();
    auto ultraTreeGraphV = treeGraphV.clone();
    
    auto ultraStart = std::chrono::high_resolution_clock::now();
    auto ultraResult = UltraParallel::NucleusCoreDecompositionUltraParallel(
        ultraTree, edgeGraph, ultraTreeGraphV, r, s);
    auto ultraEnd = std::chrono::high_resolution_clock::now();
    auto ultraTime = std::chrono::duration_cast<std::chrono::milliseconds>(ultraEnd - ultraStart).count();
    
    std::cout << "\nUltra Parallel: " << ultraTime << " ms" << std::endl;
    std::cout << "Result size: " << ultraResult.size() << std::endl;
    
    // 详细验证
    printSeparator();
    std::cout << "DETAILED CORRECTNESS VERIFICATION" << std::endl;
    printSeparator();
    
    std::cout << "\nVerifying all cliques and core values..." << std::endl;
    bool correct = verifyResults(refResult, ultraResult, true);
    
    printSeparator();
    if (correct) {
        std::cout << "✓✓✓ VERIFICATION PASSED ✓✓✓" << std::endl;
        std::cout << "All " << refResult.size() << " cliques have correct core values" << std::endl;
    } else {
        std::cout << "✗✗✗ VERIFICATION FAILED ✗✗✗" << std::endl;
        std::cout << "Some cliques have incorrect core values!" << std::endl;
        return 1;
    }
    
    // 性能比较
    std::cout << "\nPerformance:" << std::endl;
    std::cout << "  Reference:     " << refTime << " ms" << std::endl;
    std::cout << "  Ultra Parallel: " << ultraTime << " ms" << std::endl;
    
    double speedup = (double)refTime / ultraTime;
    std::cout << "  Speedup:       " << std::fixed << std::setprecision(2) << speedup << "x" << std::endl;
    
    printSeparator();
    
    return 0;
}
