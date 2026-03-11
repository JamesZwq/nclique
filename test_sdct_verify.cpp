// 测试 SDCT, SDCT_Par, SDCT_Par2 的 cliqueCount 输出是否一致

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
        std::cerr << "Usage: " << argv[0] << " <graph_file> [threads]" << std::endl;
        return 1;
    }

    const char* fpath = argv[1];
    int numThreads = (argc >= 3) ? std::stoi(argv[2]) : 4;

    populate_nCr();

    std::cout << "========================================" << std::endl;
    std::cout << "SDCT vs SDCT_Par vs SDCT_Par2 Verification" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Graph: " << fpath << ", Threads: " << numThreads << std::endl << std::endl;

    // Run SDCT
    std::cout << "[1/3] SDCT (serial)..." << std::endl;
    Graph g1(fpath); g1.sortByDegeneracyOrder();
    auto t1 = std::chrono::high_resolution_clock::now();
    auto treeSerial = SDCT(g1, 1000000, 0);
    auto serialMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - t1).count();
    std::cout << "  done: " << serialMs << " ms, leaves: " << treeSerial.adj_list.size() << std::endl;

    // Run SDCT_Par
#ifdef _OPENMP
    omp_set_num_threads(numThreads);
#endif
    std::cout << "[2/4] SDCT_Par (" << numThreads << " threads)..." << std::endl;
    Graph g2(fpath); g2.sortByDegeneracyOrder();
    auto t2 = std::chrono::high_resolution_clock::now();
    auto treePar = SDCT_Par(g2, 1000000, 0);
    auto parMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - t2).count();
    std::cout << "  done: " << parMs << " ms, leaves: " << treePar.adj_list.size() << std::endl;

    // Run SDCT_Par2
    std::cout << "[3/4] SDCT_Par2 (" << numThreads << " threads)..." << std::endl;
    Graph g3(fpath); g3.sortByDegeneracyOrder();
    auto t3 = std::chrono::high_resolution_clock::now();
    auto treePar2 = SDCT_Par2(g3, 1000000, 0);
    auto par2Ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - t3).count();
    std::cout << "  done: " << par2Ms << " ms, leaves: " << treePar2.adj_list.size() << std::endl;

    // Run SDCT_Par4
    std::cout << "[4/4] SDCT_Par4 (" << numThreads << " threads)..." << std::endl;
    Graph g4(fpath); g4.sortByDegeneracyOrder();
    auto t4 = std::chrono::high_resolution_clock::now();
    auto treePar4 = SDCT_Par4(g4, 1000000, 0);
    auto par4Ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - t4).count();
    std::cout << "  done: " << par4Ms << " ms, leaves: " << treePar4.adj_list.size() << std::endl;
    std::cout << std::endl;

    auto cSerial = treeSerial.cliqueCount();
    auto cPar    = treePar.cliqueCount();
    auto cPar2   = treePar2.cliqueCount();
    auto cPar4   = treePar4.cliqueCount();

    std::cout << "k    SDCT             SDCT_Par         SDCT_Par2        SDCT_Par4        Par Match  Par2 Match Par4 Match" << std::endl;
    std::cout << std::string(100, '-') << std::endl;

    bool allMatch = true;
    daf::Size maxK = std::max({cSerial.size(), cPar.size(), cPar2.size(), cPar4.size()});
    for (daf::Size k = 0; k < maxK; k++) {
        double sv  = (k < cSerial.size()) ? cSerial[k] : 0.0;
        double pv  = (k < cPar.size())    ? cPar[k]    : 0.0;
        double p2v = (k < cPar2.size())   ? cPar2[k]   : 0.0;
        double p4v = (k < cPar4.size())   ? cPar4[k]   : 0.0;
        if (sv == 0.0 && pv == 0.0 && p2v == 0.0 && p4v == 0.0) continue;

        bool m1 = approxEqual(sv, pv);
        bool m2 = approxEqual(sv, p2v);
        bool m4 = approxEqual(sv, p4v);
        if (!m1 || !m2 || !m4) allMatch = false;

        std::cout << std::left
                  << std::setw(5)  << k
                  << std::setw(17) << std::fixed << std::setprecision(1) << sv
                  << std::setw(17) << pv
                  << std::setw(17) << p2v
                  << std::setw(17) << p4v
                  << std::setw(11) << (m1 ? "OK" : "MISMATCH")
                  << std::setw(11) << (m2 ? "OK" : "MISMATCH")
                  << (m4 ? "OK" : "MISMATCH")
                  << std::endl;
    }

    std::cout << std::string(100, '-') << std::endl;
    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    if (allMatch) {
        std::cout << "PASS: All produce identical cliqueCount" << std::endl;
    } else {
        std::cout << "FAIL: MISMATCH detected!" << std::endl;
    }
    std::cout << "Speedup SDCT_Par:  " << std::fixed << std::setprecision(2)
              << (parMs > 0 ? (double)serialMs/parMs : 0.0) << "x" << std::endl;
    std::cout << "Speedup SDCT_Par2: " << (par2Ms > 0 ? (double)serialMs/par2Ms : 0.0) << "x" << std::endl;
    std::cout << "Speedup SDCT_Par4: " << (par4Ms > 0 ? (double)serialMs/par4Ms : 0.0) << "x" << std::endl;
    std::cout << "========================================" << std::endl;

    return allMatch ? 0 : 1;
}
