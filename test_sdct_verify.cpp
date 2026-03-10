// 测试 SDCT_Par 和 SDCT 的 cliqueCount 输出是否一致
// double值不能直接比较，使用相对误差

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

bool approxEqual(double a, double b, double relTol = 1e-9, double absTol = 1e-9) {
    if (a == b) return true;
    double diff = std::abs(a - b);
    double maxVal = std::max(std::abs(a), std::abs(b));
    if (diff <= absTol) return true;
    return diff <= relTol * maxVal;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <graph_file> [threads]" << std::endl;
        return 1;
    }

    const char* fpath = argv[1];
    int numThreads = (argc >= 3) ? std::stoi(argv[2]) : 4;

    populate_nCr();

    std::cout << "========================================" << std::endl;
    std::cout << "SDCT vs SDCT_Par cliqueCount Verification" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Graph: " << fpath << std::endl;
    std::cout << "Threads for SDCT_Par: " << numThreads << std::endl;
    std::cout << std::endl;

    // 运行 SDCT（串行版本）
    std::cout << "[1/2] Running SDCT (serial)..." << std::endl;
    std::cout.flush();
    Graph edgeGraphSerial(fpath);
    edgeGraphSerial.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> treeSerial = SDCT(edgeGraphSerial, 1000000, 0);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto serialTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "  SDCT done: " << serialTime << " ms"
              << ", leaves: " << treeSerial.adj_list.size() << std::endl;

    // 运行 SDCT_Par（并行版本）
    std::cout << "[2/2] Running SDCT_Par (" << numThreads << " threads)..." << std::endl;
    std::cout.flush();
#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif
    Graph edgeGraphPar(fpath);
    edgeGraphPar.sortByDegeneracyOrder();
    auto t3 = std::chrono::high_resolution_clock::now();
    DynamicGraph<TreeGraphNode> treePar = SDCT_Par(edgeGraphPar, 1000000, 0);
    auto t4 = std::chrono::high_resolution_clock::now();
    auto parTime = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count();
    std::cout << "  SDCT_Par done: " << parTime << " ms"
              << ", leaves: " << treePar.adj_list.size() << std::endl;
    std::cout << std::endl;

    // 计算 cliqueCount()
    auto countsSerial = treeSerial.cliqueCount();
    auto countsPar    = treePar.cliqueCount();

    std::cout << "========================================" << std::endl;
    std::cout << "cliqueCount() comparison" << std::endl;
    std::cout << "  SDCT    size: " << countsSerial.size() << std::endl;
    std::cout << "  SDCT_Par size: " << countsPar.size() << std::endl;
    std::cout << std::endl;

    std::cout << std::left
              << std::setw(5)  << "k"
              << std::setw(20) << "SDCT"
              << std::setw(20) << "SDCT_Par"
              << "Match" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    bool allMatch = true;
    daf::Size maxK = std::max(countsSerial.size(), countsPar.size());

    for (daf::Size k = 0; k < maxK; ++k) {
        double sVal = (k < countsSerial.size()) ? countsSerial[k] : 0.0;
        double pVal = (k < countsPar.size())    ? countsPar[k]    : 0.0;
        if (sVal == 0.0 && pVal == 0.0) continue;

        bool match = approxEqual(sVal, pVal);
        if (!match) allMatch = false;

        std::cout << std::left
                  << std::setw(5)  << k
                  << std::setw(20) << std::fixed << std::setprecision(1) << sVal
                  << std::setw(20) << pVal
                  << (match ? "OK" : "MISMATCH !!")
                  << std::endl;
    }

    std::cout << std::string(60, '-') << std::endl;
    std::cout << std::endl;

    // 总结
    std::cout << "========================================" << std::endl;
    if (allMatch) {
        std::cout << "PASS: SDCT and SDCT_Par produce identical cliqueCount" << std::endl;
    } else {
        std::cout << "FAIL: MISMATCH detected!" << std::endl;
    }
    std::cout << "========================================" << std::endl;

    return allMatch ? 0 : 1;
}
