#include <iostream>
#include <vector>
#include <cstring>
#include "src/graph/Graph.h"
#include "src/degeneracy_algorithm_cliques_V.h"
#include "src/misc.h"

extern double nCr[1001][401];

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s <graphFile>\n", argv[0]);
        return 1;
    }

    const char *fpath = argv[1];
    printf("Loading graph from %s...\n", fpath);
    Graph edgeGraph(fpath);
    edgeGraph.printGraphInfo();
    
    populate_nCr();
    daf::vListMap.resize(edgeGraph.n + 1);
    memset(daf::vListMap.data(), -1, edgeGraph.n * sizeof(daf::Size));

    edgeGraph.sortByDegeneracyOrder(false);
    printf("Graph sorted by degeneracy order\n\n");

    // Test SDCT
    printf("Testing SDCT...\n");
    auto tree_sdct = SDCT(edgeGraph, 1000000, 0);
    auto count_sdct = tree_sdct.cliqueCount();
    printf("SDCT clique counts: ");
    for (size_t i = 0; i < count_sdct.c_size && i < 20; i++) {
        if (count_sdct[i] > 0) printf("k%zu=%.0f ", i, count_sdct[i]);
    }
    printf("\n\n");

    // Test SDCT_Par2
    printf("Testing SDCT_Par2...\n");
    auto tree_par2 = SDCT_Par2(edgeGraph, 1000000, 0);
    auto count_par2 = tree_par2.cliqueCount();
    printf("SDCT_Par2 clique counts: ");
    for (size_t i = 0; i < count_par2.c_size && i < 20; i++) {
        if (count_par2[i] > 0) printf("k%zu=%.0f ", i, count_par2[i]);
    }
    printf("\n");
    bool par2_ok = (count_sdct.c_size == count_par2.c_size);
    if (par2_ok) {
        for (size_t i = 0; i < count_sdct.c_size; i++) {
            if (count_sdct[i] > 0) {
                // 使用相对误差而不是绝对误差，因为这些是非常大的数字
                double relError = std::abs(count_sdct[i] - count_par2[i]) / count_sdct[i];
                if (relError > 1e-10) {  // 允许 0.00000001% 的相对误差
                    par2_ok = false;
                    printf("MISMATCH at k=%zu: SDCT=%.0f, Par2=%.0f, relError=%.2e\n", i, count_sdct[i], count_par2[i], relError);
                    break;
                }
            }
        }
    }
    printf("Result: %s\n\n", par2_ok ? "✓ CORRECT" : "✗ INCORRECT");

    // Test SDCT_Par3
    printf("Testing SDCT_Par3...\n");
    auto tree_par3 = SDCT_Par3(edgeGraph, 1000000, 0);
    auto count_par3 = tree_par3.cliqueCount();
    printf("SDCT_Par3 clique counts: ");
    for (size_t i = 0; i < count_par3.c_size && i < 20; i++) {
        if (count_par3[i] > 0) printf("k%zu=%.0f ", i, count_par3[i]);
    }
    printf("\n");
    bool par3_ok = (count_sdct.c_size == count_par3.c_size);
    if (par3_ok) {
        for (size_t i = 0; i < count_sdct.c_size; i++) {
            if (count_sdct[i] > 0) {
                double relError = std::abs(count_sdct[i] - count_par3[i]) / count_sdct[i];
                if (relError > 1e-10) {
                    par3_ok = false;
                    printf("MISMATCH at k=%zu: SDCT=%.0f, Par3=%.0f, relError=%.2e\n", i, count_sdct[i], count_par3[i], relError);
                    break;
                }
            }
        }
    }
    printf("Result: %s\n\n", par3_ok ? "✓ CORRECT" : "✗ INCORRECT");

    // Test SDCT_Par4
    printf("Testing SDCT_Par4...\n");
    auto tree_par4 = SDCT_Par4(edgeGraph, 1000000, 0);
    auto count_par4 = tree_par4.cliqueCount();
    printf("SDCT_Par4 clique counts: ");
    for (size_t i = 0; i < count_par4.c_size && i < 20; i++) {
        if (count_par4[i] > 0) printf("k%zu=%.0f ", i, count_par4[i]);
    }
    printf("\n");
    bool par4_ok = (count_sdct.c_size == count_par4.c_size);
    if (par4_ok) {
        for (size_t i = 0; i < count_sdct.c_size; i++) {
            if (count_sdct[i] > 0) {
                double relError = std::abs(count_sdct[i] - count_par4[i]) / count_sdct[i];
                if (relError > 1e-10) {
                    par4_ok = false;
                    printf("MISMATCH at k=%zu: SDCT=%.0f, Par4=%.0f, relError=%.2e\n", i, count_sdct[i], count_par4[i], relError);
                    break;
                }
            }
        }
    }
    printf("Result: %s\n\n", par4_ok ? "✓ CORRECT" : "✗ INCORRECT");

    // Test SDCT_Par5
    printf("Testing SDCT_Par5...\n");
    auto tree_par5 = SDCT_Par5(edgeGraph, 1000000, 0);
    auto count_par5 = tree_par5.cliqueCount();
    printf("SDCT_Par5 clique counts: ");
    for (size_t i = 0; i < count_par5.c_size && i < 20; i++) {
        if (count_par5[i] > 0) printf("k%zu=%.0f ", i, count_par5[i]);
    }
    printf("\n");
    bool par5_ok = (count_sdct.c_size == count_par5.c_size);
    if (par5_ok) {
        for (size_t i = 0; i < count_sdct.c_size; i++) {
            if (count_sdct[i] > 0) {
                double relError = std::abs(count_sdct[i] - count_par5[i]) / count_sdct[i];
                if (relError > 1e-10) {
                    par5_ok = false;
                    printf("MISMATCH at k=%zu: SDCT=%.0f, Par5=%.0f, relError=%.2e\n", i, count_sdct[i], count_par5[i], relError);
                    break;
                }
            }
        }
    }
    printf("Result: %s\n\n", par5_ok ? "✓ CORRECT" : "✗ INCORRECT");

    // Summary
    printf("\n========================================\n");
    printf("SUMMARY\n");
    printf("========================================\n");
    printf("SDCT:     ✓ Reference\n");
    printf("SDCT_Par2: %s\n", par2_ok ? "✓ CORRECT" : "✗ INCORRECT");
    printf("SDCT_Par3: %s\n", par3_ok ? "✓ CORRECT" : "✗ INCORRECT");
    printf("SDCT_Par4: %s\n", par4_ok ? "✓ CORRECT" : "✗ INCORRECT");
    printf("SDCT_Par5: %s\n", par5_ok ? "✓ CORRECT" : "✗ INCORRECT");
    printf("========================================\n");

    if (par2_ok && par3_ok && par4_ok && par5_ok) {
        printf("\n✓ ALL TESTS PASSED!\n");
        return 0;
    } else {
        printf("\n✗ SOME TESTS FAILED!\n");
        return 1;
    }
}
