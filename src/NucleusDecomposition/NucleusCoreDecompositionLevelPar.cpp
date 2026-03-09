// Level-based Parallel Nucleus Decomposition
// 核心思想：一次性并行处理所有相同 support 值的 r-cliques

#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "BK/BronKerboschRmRClique.hpp"
#include <vector>
#include <algorithm>
#include <atomic>
#include <unordered_map>
#include <ranges>
#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

namespace LevelParallel {

std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRCliqueLevelPar(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    
    auto time_start = std::chrono::high_resolution_clock::now();
    
    // Step 1: Build clique index
    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });
    
    const daf::Size nClique = cliqueIndex.size();
    const daf::Size graphN = edgeGraph.n;
    
    // Step 2: 计算初始 counting（并行优化版本）
    auto countingRClique = daf::timeCount("countingPerRClique", [&]() {
        std::vector<double> rCliqueSCounting(nClique, 0.0);
#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
        std::vector<std::vector<double>> thread_locals(nthreads, std::vector<double>(nClique, 0.0));
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto &local = thread_locals[tid];
#pragma omp for schedule(dynamic, 64)
            for (daf::Size leafIdx = 0; leafIdx < tree.adj_list.size(); ++leafIdx) {
                const auto &leaf = tree.adj_list[leafIdx];
                if (leaf.size() < r) continue;
                daf::CliqueSize pivotC = 0, keepC = 0;
                for (const auto &i : leaf) {
                    if (i.isPivot) pivotC++; else keepC++;
                }
                int needPivot = s - static_cast<int>(keepC);
                daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                    daf::CliqueSize subNumPovit = 0;
                    for (const auto &node : rClique) if (node.isPivot) subNumPovit++;
                    double ncrValue = nCr[pivotC - subNumPovit][needPivot - subNumPovit];
                    auto id = cliqueIndex.byClique(rClique);
                    if (id < nClique) local[id] += ncrValue;
                    return true;
                });
            }
        }
#pragma omp parallel for schedule(static)
        for (daf::Size i = 0; i < nClique; ++i) {
            for (int t = 0; t < nthreads; ++t) {
                rCliqueSCounting[i] += thread_locals[t][i];
            }
        }
#else
        for (const auto &leaf : tree.adj_list) {
            if (leaf.size() < r) continue;
            daf::CliqueSize pivotC = 0, keepC = 0;
            for (const auto &i : leaf) { if (i.isPivot) pivotC++; else keepC++; }
            int needPivot = s - static_cast<int>(keepC);
            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                daf::CliqueSize subNumPovit = 0;
                for (const auto &node : rClique) if (node.isPivot) subNumPovit++;
                rCliqueSCounting[cliqueIndex.byClique(rClique)] += nCr[pivotC - subNumPovit][needPivot - subNumPovit];
                return true;
            });
        }
#endif
        return rCliqueSCounting;
    });
    
    daf::log_memory("Other index(incloud head)");
    
    // Step 3: Level-based Parallel Peeling
    std::vector<int> coreRClique(nClique, 0);
    std::vector<bool> removed(nClique, false);
    
    // 使用 atomic 的 counting 数组
    std::vector<std::atomic<double>> atomicCounting(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        atomicCounting[i].store(countingRClique[i], std::memory_order_relaxed);
    }
    
    // 用于标记受影响的 leaves
    std::vector<daf::Size> changedLeafIndex(graphN, std::numeric_limits<daf::Size>::max());
    
    std::cout << "=========================begin=========================" << std::endl;
    
    int iterationCount = 0;
    size_t totalRemoved = 0;
    
    while (totalRemoved < nClique) {
        // Step 3.1: 找到当前最小 support 值
        double minSupport = std::numeric_limits<double>::max();
        for (daf::Size i = 0; i < nClique; ++i) {
            if (!removed[i]) {
                double val = atomicCounting[i].load(std::memory_order_relaxed);
                if (val < minSupport) {
                    minSupport = val;
                }
            }
        }
        
        if (minSupport == std::numeric_limits<double>::max()) break;
        
        // Step 3.2: 收集所有 support == minSupport 的 r-cliques
        std::vector<daf::Size> currentLevel;
        for (daf::Size i = 0; i < nClique; ++i) {
            if (!removed[i]) {
                double val = atomicCounting[i].load(std::memory_order_relaxed);
                if (std::abs(val - minSupport) < 0.001) {  // 浮点数容差
                    currentLevel.push_back(i);
                }
            }
        }
        
        if (currentLevel.empty()) break;
        
        // Step 3.3: 标记这些 r-cliques 为已移除
        for (auto cliqueId : currentLevel) {
            coreRClique[cliqueId] = static_cast<int>(minSupport);
            removed[cliqueId] = true;
        }
        totalRemoved += currentLevel.size();
        
        std::cout << "minCore: " << static_cast<int>(minSupport) 
                  << " removed: " << currentLevel.size() 
                  << " total: " << totalRemoved << "/" << nClique << std::endl;
        
        // Step 3.4: 并行收集所有受影响的 leaves
        std::vector<daf::Size> changedLeaf;
        std::vector<std::vector<daf::Size>> removedRCliqueIdForLeaf;
        
#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
        std::vector<std::vector<std::pair<daf::Size, daf::Size>>> threadLeafPairs(nthreads);
        
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto &localPairs = threadLeafPairs[tid];
            
#pragma omp for schedule(dynamic, 4)
            for (size_t idx = 0; idx < currentLevel.size(); ++idx) {
                auto rmCliqueId = currentLevel[idx];
                auto rClique = cliqueIndex.byId(rmCliqueId);
                daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                    [&](const TreeGraphNode &uClique) {
                        localPairs.emplace_back(idx, uClique.v);
                    });
            }
        }
        
        // 合并并排序
        std::vector<std::pair<daf::Size, daf::Size>> allPairs;
        size_t totalSize = 0;
        for (const auto &tp : threadLeafPairs) totalSize += tp.size();
        allPairs.reserve(totalSize);
        for (const auto &tp : threadLeafPairs) {
            allPairs.insert(allPairs.end(), tp.begin(), tp.end());
        }
        std::sort(allPairs.begin(), allPairs.end(), 
                  [](const auto &a, const auto &b) { return a.second < b.second; });
        
        for (const auto &p : allPairs) {
            daf::Size origIdx = p.first;
            daf::Size leafV = p.second;
            auto rmCliqueId = currentLevel[origIdx];
            auto &leafId = changedLeafIndex[leafV];
            if (leafId == std::numeric_limits<daf::Size>::max()) {
                leafId = removedRCliqueIdForLeaf.size();
                removedRCliqueIdForLeaf.emplace_back();
                changedLeaf.push_back(leafV);
                removedRCliqueIdForLeaf.back().reserve(10);
            }
            removedRCliqueIdForLeaf[leafId].emplace_back(rmCliqueId);
        }
#else
        for (auto rmCliqueId : currentLevel) {
            auto rClique = cliqueIndex.byId(rmCliqueId);
            daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                [&](const TreeGraphNode &uClique) {
                    auto &leafId = changedLeafIndex[uClique.v];
                    if (leafId == std::numeric_limits<daf::Size>::max()) {
                        leafId = removedRCliqueIdForLeaf.size();
                        removedRCliqueIdForLeaf.emplace_back();
                        changedLeaf.push_back(uClique.v);
                        removedRCliqueIdForLeaf.back().reserve(10);
                    }
                    removedRCliqueIdForLeaf[leafId].emplace_back(rmCliqueId);
                });
        }
#endif
        
        // Step 3.5: 并行处理每个受影响的 leaf
#ifdef _OPENMP
        struct LeafUpdate {
            std::vector<std::pair<daf::Size, double>> increments;
            std::vector<std::pair<daf::Size, double>> decrements;
        };
        std::vector<LeafUpdate> leafUpdates(changedLeaf.size());
        
#pragma omp parallel
        {
            static thread_local daf::StaticVector<daf::Size> threadMap;
            if (threadMap.maxSize < graphN + 1) threadMap.resize(graphN + 1);
            
#pragma omp for schedule(dynamic, 16)
            for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
                auto leafId = changedLeaf[idx];
                const auto leaf = tree.adj_list[leafId];
                auto leafIndex = changedLeafIndex[leafId];
                auto removedR = removedRCliqueIdForLeaf[leafIndex];
                
                // 计算 decrements
                for (auto rmId : removedR) {
                    auto rClique = cliqueIndex.byId(rmId);
                    daf::CliqueSize pivotC = 0, keepC = 0;
                    for (const auto &node : leaf) {
                        if (node.isPivot) pivotC++; else keepC++;
                    }
                    daf::CliqueSize subNumPovit = 0;
                    for (const auto &node : rClique) {
                        // rClique 是 span<const daf::Size>，需要在 leaf 中查找对应的 TreeGraphNode
                        for (const auto &leafNode : leaf) {
                            if (leafNode.v == node) {
                                if (leafNode.isPivot) subNumPovit++;
                                break;
                            }
                        }
                    }
                    int needPivot = s - static_cast<int>(keepC);
                    if (pivotC >= subNumPovit && needPivot >= subNumPovit &&
                        pivotC - subNumPovit < 1001 && needPivot - subNumPovit < 401) {
                        double ncrValue = nCr[pivotC - subNumPovit][needPivot - subNumPovit];
                        leafUpdates[idx].decrements.emplace_back(rmId, ncrValue);
                    }
                }
                
                // 使用 Bron-Kerbosch 生成新 leaves 并计算 increments
                auto mapped = removedR | std::views::transform(
                    [&](const daf::Size id) { return cliqueIndex.byId(id); });
                
                std::vector<std::vector<TreeGraphNode>> newLeaves;
                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                        auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                        if (newLeaf.size() >= s) {
                            newLeaves.push_back(std::move(newLeaf));
                        }
                    },
                    &threadMap);
                
                for (const auto &newLeaf : newLeaves) {
                    daf::CliqueSize newPivotC = 0, newKeepC = 0;
                    for (const auto &j : newLeaf) {
                        if (j.isPivot) newPivotC++; else newKeepC++;
                    }
                    daf::Size needPivot = s - newKeepC;
                    daf::enumerateCombinations(newLeaf, r, [&](const daf::StaticVector<TreeGraphNode> &rclique) {
                        daf::CliqueSize subNumPovit = 0;
                        for (const auto &node : rclique) if (node.isPivot) subNumPovit++;
                        if (subNumPovit <= needPivot && newPivotC >= subNumPovit &&
                            newPivotC - subNumPovit < 1001 && needPivot >= subNumPovit && needPivot - subNumPovit < 401) {
                            double ncrValue = nCr[newPivotC - subNumPovit][needPivot - subNumPovit];
                            auto cliqueId = cliqueIndex.byClique(rclique);
                            leafUpdates[idx].increments.emplace_back(cliqueId, ncrValue);
                        }
                        return true;
                    });
                }
            }
        }
        
        // Step 3.6: 并行应用所有更新（使用原子操作）
#pragma omp parallel for schedule(static)
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            for (const auto &[cliqueId, value] : leafUpdates[idx].decrements) {
                double current = atomicCounting[cliqueId].load(std::memory_order_relaxed);
                while (!atomicCounting[cliqueId].compare_exchange_weak(
                    current, current - value, std::memory_order_relaxed)) {}
            }
            for (const auto &[cliqueId, value] : leafUpdates[idx].increments) {
                double current = atomicCounting[cliqueId].load(std::memory_order_relaxed);
                while (!atomicCounting[cliqueId].compare_exchange_weak(
                    current, current + value, std::memory_order_relaxed)) {}
            }
        }
#else
        // 串行版本的更新逻辑
        // ... (省略，与 Ref 版本类似)
#endif
        
        // 清理
        for (auto leafV : changedLeaf) {
            changedLeafIndex[leafV] = std::numeric_limits<daf::Size>::max();
        }
        
        iterationCount++;
    }
    
    std::cout << "Total iterations: " << iterationCount << std::endl;
    
    // Step 4: 构造返回结果
    std::vector<std::pair<std::vector<daf::Size>, int>> result;
    result.reserve(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueVec(clique.begin(), clique.end());
        result.emplace_back(std::move(cliqueVec), coreRClique[i]);
    }
    
    auto time_end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    std::cout << "Level-based parallel peeling took: " << elapsed << " ms" << std::endl;
    
    return result;
}

} // namespace LevelParallel
