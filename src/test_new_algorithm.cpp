// 测试新算法的程序

#include <iostream>
#include <chrono>
#include "graph/Graph.h"
#include "tree/MultiBranchTree.h"
#include "NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "degeneracy_algorithm_cliques_V.h"

// 声明新算法
namespace OptimizedNucleus {
std::vector<std::pair<std::vector<daf::Size>, double>> 
NucleusCoreDecompositionOptimized(
    DynamicGraph<TreeGraphNode>& tree,
    const Graph& edgeGraph,
    DynamicGraphSet<TreeGraphNode>& treeGraphV,
    daf::CliqueSize r,
    daf::CliqueSize s,
    int numThreads);
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <graph_file> <r> <s> [num_threads]" << std::endl;
        return 1;
    }
    
    const char* fpath = argv[1];
    int r = std::stoi(argv[2]);
    int s = std::stoi(argv[3]);
    int numThreads = (argc >= 5) ? std::stoi(argv[4]) : 16;
    
    std::cout << "========================================" << std::endl;
    std::cout << "Testing New Nucleus Decomposition Algorithm" << std::endl;
    std::cout << "Graph: " << fpath << std::endl;
    std::cout << "r=" << r << ", s=" << s << std::endl;
    std::cout << "Threads: " << numThreads << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // 加载图
    Graph edgeGraph(fpath);
    edgeGraph.printGraphInfo();
    edgeGraph.sortByDegeneracyOrder();
    
    // 构建 tree
    auto treeGraph = SDCT_Par(edgeGraph, 1000000, 0);
    std::cout << "Tree leaves: " << treeGraph.adj_list.size() << std::endl;
    
    // 构建 treeGraphV
    DynamicGraphSet<TreeGraphNode> treeGraphV(treeGraph, edgeGraph.getGraphNodeSize(), s);
    
    // 测试新算法
    std::cout << "\n=== Testing NEW Algorithm ===" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    auto newResult = OptimizedNucleus::NucleusCoreDecompositionOptimized(
        treeGraph, edgeGraph, treeGraphV, r, s, numThreads);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto newTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    std::cout << "\nNEW Algorithm: " << newTime << " ms" << std::endl;
    std::cout << "Result size: " << newResult.size() << std::endl;
    
    // 测试旧算法（参考）
    std::cout << "\n=== Testing OLD Algorithm (Reference) ===" << std::endl;
    auto refTree = treeGraph.clone();
    auto refTreeGraphV = treeGraphV.clone();
    
    start = std::chrono::high_resolution_clock::now();
    
    auto oldResult = NucleusCoreDecompositionRCliqueRef(
        refTree, edgeGraph, refTreeGraphV, r, s);
    
    end = std::chrono::high_resolution_clock::now();
    auto oldTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    std::cout << "\nOLD Algorithm: " << oldTime << " ms" << std::endl;
    std::cout << "Result size: " << oldResult.size() << std::endl;
    
    // 比较结果
    std::cout << "\n=== Comparison ===" << std::endl;
    std::cout << "Speedup: " << (double)oldTime / newTime << "x" << std::endl;
    
    if (newResult.size() != oldResult.size()) {
        std::cout << "WARNING: Result sizes differ!" << std::endl;
    } else {
        std::cout << "Result sizes match!" << std::endl;
    }
    
    return 0;
}
