// 测试Ultra Parallel算法的性能

#include <iostream>
#include <chrono>
#include <iomanip>
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
    std::cout << "Ultra-Parallel Nucleus Decomposition Test" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Graph: " << fpath << std::endl;
    std::cout << "Parameters: r=" << r << ", s=" << s << std::endl;
    std::cout << "Threads: " << numThreads << std::endl;
    printSeparator();
    
    // 加载图
    std::cout << "\n[1/5] Loading graph..." << std::endl;
    Graph edgeGraph(fpath);
    edgeGraph.printGraphInfo();
    edgeGraph.sortByDegeneracyOrder();
    
    // 构建 tree
    std::cout << "\n[2/5] Building tree structure..." << std::endl;
    auto treeGraph = SDCT_Par(edgeGraph, 1000000, 0);
    std::cout << "Tree leaves: " << treeGraph.adj_list.size() << std::endl;
    
    // 构建 treeGraphV
    std::cout << "\n[3/5] Building treeGraphV..." << std::endl;
    DynamicGraphSet<TreeGraphNode> treeGraphV(treeGraph, edgeGraph.getGraphNodeSize(), s);
    
    // 测试Ultra Parallel算法
    printSeparator();
    std::cout << "[4/5] Running ULTRA PARALLEL Algorithm" << std::endl;
    printSeparator();
    
    auto ultraTree = treeGraph.clone();
    auto ultraTreeGraphV = treeGraphV.clone();
    
    auto ultraStart = std::chrono::high_resolution_clock::now();
    auto ultraResult = UltraParallel::NucleusCoreDecompositionUltraParallel(
        ultraTree, edgeGraph, ultraTreeGraphV, r, s);
    auto ultraEnd = std::chrono::high_resolution_clock::now();
    auto ultraTime = std::chrono::duration_cast<std::chrono::milliseconds>(ultraEnd - ultraStart).count();
    
    // 测试参考算法（单线程）
    printSeparator();
    std::cout << "[5/5] Running REFERENCE Algorithm (for comparison)" << std::endl;
    printSeparator();
    
#ifdef _OPENMP
    omp_set_num_threads(1);  // 强制单线程
#endif
    
    auto refTree = treeGraph.clone();
    auto refTreeGraphV = treeGraphV.clone();
    
    auto refStart = std::chrono::high_resolution_clock::now();
    auto refResult = NucleusCoreDecompositionRCliqueRef(
        refTree, edgeGraph, refTreeGraphV, r, s);
    auto refEnd = std::chrono::high_resolution_clock::now();
    auto refTime = std::chrono::duration_cast<std::chrono::milliseconds>(refEnd - refStart).count();
    
    // 结果比较
    printSeparator();
    std::cout << "PERFORMANCE COMPARISON" << std::endl;
    printSeparator();
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\nReference (1 thread):  " << refTime << " ms" << std::endl;
    std::cout << "Ultra Parallel (" << numThreads << " threads): " << ultraTime << " ms" << std::endl;
    
    double speedup = (double)refTime / ultraTime;
    std::cout << "\nSpeedup: " << speedup << "x" << std::endl;
    
    if (speedup >= 3.0) {
        std::cout << "✓ TARGET ACHIEVED! (>= 3x speedup)" << std::endl;
    } else {
        std::cout << "✗ Target not met (need >= 3x speedup)" << std::endl;
    }
    
    // 验证正确性
    std::cout << "\nCorrectness Check:" << std::endl;
    std::cout << "Reference result size: " << refResult.size() << std::endl;
    std::cout << "Ultra result size:     " << ultraResult.size() << std::endl;
    
    if (refResult.size() == ultraResult.size()) {
        std::cout << "✓ Result sizes match" << std::endl;
        
        // 简单验证：检查前几个结果
        bool match = true;
        int checkCount = std::min(10, (int)refResult.size());
        for (int i = 0; i < checkCount; ++i) {
            if (refResult[i].second != ultraResult[i].second) {
                match = false;
                break;
            }
        }
        
        if (match) {
            std::cout << "✓ Sample results match (first " << checkCount << " entries)" << std::endl;
        } else {
            std::cout << "✗ WARNING: Results differ!" << std::endl;
        }
    } else {
        std::cout << "✗ WARNING: Result sizes differ!" << std::endl;
    }
    
    printSeparator();
    
    return 0;
}
