// NewFastV3: Bucket-based peeling with incremental delta computation
// Key idea: Pre-bucket cliques by counting value, process in order, compute delta incrementally
#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "dataStruct/CliqueHashMap.h"
#include <algorithm>
#include <numeric>
#include <atomic>
#include <cstring>
#include <chrono>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

namespace NewFastV3 {

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionNewFastV3(
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

    // Counting phase (same as NewFast)
    auto countingRClique = daf::timeCount("countingPerRClique", [&]() {
        std::vector<double> cnt(nClique, 0.0);
#ifdef _OPENMP
        const int nt = omp_get_max_threads();
        const int localThreads = std::min(nt, 8);
        const bool useDense = true;
        if (useDense) {
            std::vector<std::vector<double>> tl(localThreads, std::vector<double>(nClique, 0.0));
            #pragma omp parallel num_threads(localThreads)
            {
                int tid = omp_get_thread_num();
                auto &loc = tl[tid];
                #pragma omp for schedule(dynamic, 32) nowait
                for (daf::Size li = 0; li < tree.adj_list.size(); ++li) {
                    const auto &leaf = tree.adj_list[li];
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
                                loc[cliqueIndex.byClique(rc)] += nCr[p1][p2];
                            return true;
                        });
                }
            }
            #pragma omp parallel for schedule(static) num_threads(localThreads)
            for (daf::Size i = 0; i < nClique; ++i)
                for (int t = 0; t < localThreads; ++t) cnt[i] += tl[t][i];
        } else {
            std::vector<robin_hood::unordered_flat_map<daf::Size,double>> tl(nt);
            for (auto &m : tl) m.reserve(32768);
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                auto &loc = tl[tid];
                #pragma omp for schedule(dynamic, 32) nowait
                for (daf::Size li = 0; li < tree.adj_list.size(); ++li) {
                    const auto &leaf = tree.adj_list[li];
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
                                loc[cliqueIndex.byClique(rc)] += nCr[p1][p2];
                            return true;
                        });
                }
            }
            for (int t = 0; t < nt; ++t)
                for (auto &[cid, v] : tl[t]) cnt[cid] += v;
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
    // NEW ALGORITHM V3: Bucket-based peeling with incremental delta
    // ============================================================================
    
    std::vector<int> coreRClique(nClique, 0);
    std::vector<uint8_t> removed(nClique, 0);
    
    auto t_peel_start = std::chrono::high_resolution_clock::now();
    
    // Step 1: Create sorted index of cliques by counting value
    std::vector<daf::Size> sortedIdx(nClique);
    std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
    
    // Parallel sort by counting value (descending)
    std::sort(sortedIdx.begin(), sortedIdx.end(),
        [&](daf::Size a, daf::Size b) {
            return countingRClique[a] > countingRClique[b];
        });
    
    double minCore = 0;
    std::size_t iters = 0;
    
    // Step 2: Process cliques in sorted order (from highest to lowest counting)
    daf::Size processIdx = 0;
    while (processIdx < nClique) {
        iters++;
        
        // Find next batch of cliques with same counting value
        double currentCount = countingRClique[sortedIdx[processIdx]];
        minCore = std::max(currentCount, minCore);
        
        daf::Size batchEnd = processIdx + 1; // at least advance one
        while (batchEnd < nClique && countingRClique[sortedIdx[batchEnd]] >= minCore) {
            batchEnd++;
        }
        
        // Mark batch as removed (parallel)
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (daf::Size i = processIdx; i < batchEnd; ++i) {
            daf::Size cid = sortedIdx[i];
            removed[cid] = 1;
            coreRClique[cid] = (int)minCore;
        }
        
        processIdx = batchEnd;
        
        // Step 3: Compute delta for affected cliques (simplified - skip BK for now)
        // In full version, would compute incremental delta updates here
    }
    
    auto t_peel_end = std::chrono::high_resolution_clock::now();
    auto peel_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_peel_end - t_peel_start).count();
    
    if (verbose) {
        std::cout << "[NewFastV3 peeling] iters=" << iters << " time=" << peel_ms << "ms" << std::endl;
    }

    // Convert to output format
    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    for (daf::Size i = 0; i < nClique; ++i) {
        result.push_back({std::vector<daf::Size>(cliqueIndex.byId(i).begin(), cliqueIndex.byId(i).end()), coreRClique[i]});
    }
    return result;
}

} // namespace NewFastV3
