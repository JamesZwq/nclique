// NewFastHybrid: Simplified fusion - adaptive counting + V2 peeling
#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include <algorithm>
#include <numeric>
#include <atomic>
#include <cstring>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

namespace NewFastHybrid {

std::vector<std::pair<std::vector<daf::Size>, int>>
NucleusCoreDecompositionNewFastHybrid(
    DynamicGraph<TreeGraphNode> &tree,
    const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r,
    daf::CliqueSize s)
{
    const bool verbose = std::getenv("PIVOTER_VERBOSE") != nullptr;

    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });
    const daf::Size nClique = cliqueIndex.size();
    const daf::Size graphN  = edgeGraph.n;

    // Adaptive counting based on nClique size
    auto countingRClique = daf::timeCount("countingPerRClique", [&]() {
        std::vector<double> cnt(nClique, 0.0);
#ifdef _OPENMP
        const int nt = omp_get_max_threads();
        const int localThreads = std::min(nt, 8);
        const bool useDense = true; // Always use dense for simplicity
        
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
                        if (p1>=0 && p1<1001 && p2>=0 && p2<401) {
                            loc[cliqueIndex.byClique(rc)] += nCr[p1][p2];
                        }
                        return true;
                    });
            }
        }
        
        #pragma omp parallel for schedule(static) num_threads(localThreads)
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

    // V2-style peeling: array scan + parallel batch
    std::vector<int> coreRClique(nClique, 0);
    std::vector<uint8_t> removed(nClique, 0);
    
    double minCore = 0;
    std::size_t iters = 0;
    
    while (true) {
        iters++;
        
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

    // Convert to output format
    std::vector<std::pair<std::vector<daf::Size>, int>> result;
    for (daf::Size i = 0; i < nClique; ++i) {
        result.push_back({std::vector<daf::Size>(cliqueIndex.byId(i).begin(), cliqueIndex.byId(i).end()), coreRClique[i]});
    }
    return result;
}

} // namespace NewFastHybrid
