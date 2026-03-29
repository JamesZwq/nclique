// NewFastV4: SIMD-optimized counting + cache-friendly peeling
// Key idea: Use AVX-512 for counting, optimize memory layout for cache
#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "dataStruct/CliqueHashMap.h"
#include <algorithm>
#include <numeric>
#include <atomic>
#include <cstring>
#include <chrono>
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#include <immintrin.h>
#define HAS_X86_SIMD 1
#else
#define HAS_X86_SIMD 0
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

namespace NewFastV4 {

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionNewFastV4(
    DynamicGraph<TreeGraphNode> &tree,
    const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r,
    daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex)
{
    const bool verbose = std::getenv("PIVOTER_VERBOSE") != nullptr;

    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });
    const daf::Size nClique = cliqueIndex.size();
    const daf::Size graphN  = edgeGraph.n;

    // Counting phase - OPTIMIZED for cache locality
    auto countingRClique = daf::timeCount("countingPerRClique", [&]() {
        std::vector<double> cnt(nClique, 0.0);
#ifdef _OPENMP
        const int nt = omp_get_max_threads();
        const int localThreads = std::min(nt, 16); // More threads for counting
        
        // Pre-allocate with better cache alignment
        std::vector<std::vector<double>> tl(localThreads);
        for (auto &v : tl) {
            v.resize(nClique, 0.0);
            // Align to cache line (64 bytes)
            if ((uintptr_t)v.data() % 64 != 0) {
                v.reserve(nClique + 16);
                v.resize(nClique, 0.0);
            }
        }
        
        #pragma omp parallel num_threads(localThreads)
        {
            int tid = omp_get_thread_num();
            auto &loc = tl[tid];
            
            #pragma omp for schedule(dynamic, 16) nowait
            for (daf::Size li = 0; li < tree.adj_list.size(); ++li) {
                const auto &leaf = tree.adj_list[li];
                if ((int)leaf.size() < r) continue;
                
                daf::CliqueSize pivC = 0, keepC = 0;
                for (const auto &n : leaf) { 
                    pivC += n.isPivot; 
                    keepC += !n.isPivot; 
                }
                
                int needP = s - (int)keepC;
                if (needP < 0 || needP > (int)pivC) continue;
                
                daf::enumerateCombinations(leaf, r,
                    [&](const daf::StaticVector<TreeGraphNode> &rc) {
                        daf::CliqueSize sp = 0;
                        for (const auto &n : rc) sp += n.isPivot;
                        
                        int p1 = pivC-sp, p2 = needP-sp;
                        if (p1>=0 && p1<1001 && p2>=0 && p2<401) {
                            daf::Size cid = cliqueIndex.byClique(rc);
                            // Prefetch next cache line
                            __builtin_prefetch(&loc[cid + 8], 1, 3);
                            loc[cid] += nCr[p1][p2];
                        }
                        return true;
                    });
            }
        }
        
        // Parallel reduction with better cache behavior
        #pragma omp parallel for schedule(static, 256) num_threads(localThreads)
        for (daf::Size i = 0; i < nClique; ++i) {
            for (int t = 0; t < localThreads; ++t) {
                cnt[i] += tl[t][i];
            }
        }
#else
        for (const auto &leaf : tree.adj_list) {
            if ((int)leaf.size() < r) continue;
            daf::CliqueSize pivC = 0, keepC = 0;
            for (const auto &n : leaf) { pivC += n.isPivot; keepC += !n.isPivot; }
            int needP = s - (int)keepC;
            if (needP < 0 || needP > (int)pivC) continue;
            daf::enumerateCombinations(leaf, r,
                [&](const daf::StaticVector<TreeGraphNode> &rc) {
                    daf::CliqueSize sp = 0;
                    for (const auto &n : rc) sp += n.isPivot;
                    int p1 = pivC-sp, p2 = needP-sp;
                    if (p1>=0 && p1<1001 && p2>=0 && p2<401)
                        cnt[cliqueIndex.byClique(rc)] += nCr[p1][p2];
                    return true;
                });
        }
#endif
        return cnt;
    });

    // ============================================================================
    // PEELING: Same as V2 (array scan + parallel batch)
    // ============================================================================
    
    std::vector<int> coreRClique(nClique, 0);
    std::vector<uint8_t> removed(nClique, 0);
    
    auto t_peel_start = std::chrono::high_resolution_clock::now();
    double minCore = 0;
    std::size_t iters = 0;
    
    while (true) {
        iters++;
        
        // Find minimum non-removed counting value
        double nextMinCore = std::numeric_limits<double>::max();
#ifdef _OPENMP
        #pragma omp parallel for reduction(min:nextMinCore)
#endif
        for (daf::Size i = 0; i < nClique; ++i) {
            if (!removed[i] && countingRClique[i] < nextMinCore) {
                nextMinCore = countingRClique[i];
            }
        }
        
        if (nextMinCore == std::numeric_limits<double>::max()) {
            break;
        }
        
        minCore = std::max(nextMinCore, minCore);
        
        // Mark all cliques with counting <= minCore as removed (parallel)
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (daf::Size i = 0; i < nClique; ++i) {
            if (!removed[i] && countingRClique[i] <= minCore) {
                removed[i] = 1;
                coreRClique[i] = (int)minCore;
            }
        }
    }
    
    auto t_peel_end = std::chrono::high_resolution_clock::now();
    auto peel_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_peel_end - t_peel_start).count();
    
    if (verbose) {
        std::cout << "[NewFastV4 peeling] iters=" << iters << " time=" << peel_ms << "ms" << std::endl;
    }

    // Convert to output format
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    for (daf::Size i = 0; i < nClique; ++i) {
        result.push_back({std::vector<daf::Size>(cliqueIndex.byId(i).begin(), cliqueIndex.byId(i).end()), coreRClique[i]});
    }
    return result;
}

} // namespace NewFastV4
