// Ultra-Parallel Nucleus Decomposition
// 目标：8线程达到3倍以上加速
// 
// 核心创新策略：
// 1. Lock-Free Batch Processing：批量处理相同support的r-cliques，完全无锁
// 2. Hierarchical Parallelism：三层并行（batch级、leaf级、computation级）
// 3. Memory-Efficient Design：使用线程局部缓冲区避免同步
// 4. Smart Load Balancing：动态任务分配 + work stealing
// 5. Optimized Data Structures：使用flat arrays替代复杂结构

#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "BK/BronKerboschRmRClique.hpp"
#include <boost/heap/d_ary_heap.hpp>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <chrono>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

namespace UltraParallel {

// ============================================================================
// 高性能数据结构
// ============================================================================

// 线程局部工作缓冲区（避免动态分配）
struct alignas(64) ThreadLocalWorkspace {
    // 用于BK算法的工作空间
    daf::StaticVector<daf::Size> workMap;
    
    // 用于收集更新的缓冲区
    std::vector<std::pair<daf::Size, double>> increments;
    std::vector<std::pair<daf::Size, double>> decrements;
    std::vector<std::vector<TreeGraphNode>> newLeaves;
    
    // 用于leaf处理的临时数据
    std::vector<daf::Size> tempCliqueIds;
    
    explicit ThreadLocalWorkspace(size_t graphSize) {
        workMap.resize(graphSize + 1);
        increments.reserve(5000);
        decrements.reserve(5000);
        newLeaves.reserve(50);
        tempCliqueIds.reserve(100);
    }
    
    void clear() {
        increments.clear();
        decrements.clear();
        newLeaves.clear();
        tempCliqueIds.clear();
    }
};

// Leaf处理任务
struct LeafTask {
    daf::Size leafId;
    daf::Size removedStartIdx;  // 在全局removedCliques数组中的起始位置
    daf::Size removedCount;     // 该leaf需要处理的removed cliques数量
};

// ============================================================================
// 核心计算函数
// ============================================================================

// 并行计算初始support（高度优化版本）
std::vector<double> computeInitialSupportOptimized(
    const DynamicGraph<TreeGraphNode> &tree,
    const StaticCliqueIndex &cliqueIndex,
    daf::CliqueSize r,
    daf::CliqueSize s) {
    
    const daf::Size nClique = cliqueIndex.size();
    std::vector<double> support(nClique, 0.0);
    
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
    std::vector<std::vector<double>> threadLocal(nthreads, std::vector<double>(nClique, 0.0));
    
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto &local = threadLocal[tid];
        
        // 使用dynamic调度以平衡负载
        #pragma omp for schedule(dynamic, 32) nowait
        for (daf::Size leafIdx = 0; leafIdx < tree.adj_list.size(); ++leafIdx) {
            const auto &leaf = tree.adj_list[leafIdx];
            const auto leafSize = leaf.size();
            if (leafSize < r) continue;
            
            // 快速计算pivot和keep数量
            daf::CliqueSize pivotC = 0, keepC = 0;
            for (const auto &node : leaf) {
                pivotC += node.isPivot;
                keepC += !node.isPivot;
            }
            
            const int needPivot = s - static_cast<int>(keepC);
            if (needPivot < 0 || needPivot > pivotC) continue;
            
            // 枚举所有r-cliques
            daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                daf::CliqueSize subNumPivot = 0;
                for (const auto &node : rClique) {
                    subNumPivot += node.isPivot;
                }
                
                const int p1 = pivotC - subNumPivot;
                const int p2 = needPivot - subNumPivot;
                
                // 边界检查
                if (p1 >= 0 && p1 < 1001 && p2 >= 0 && p2 < 401) {
                    const double ncrValue = nCr[p1][p2];
                    const auto id = cliqueIndex.byClique(rClique);
                    if (id < nClique) {
                        local[id] += ncrValue;
                    }
                }
                return true;
            });
        }
    }
    
    // 并行归约（使用SIMD友好的访问模式）
    #pragma omp parallel for schedule(static, 256)
    for (daf::Size i = 0; i < nClique; ++i) {
        double sum = 0.0;
        for (int t = 0; t < nthreads; ++t) {
            sum += threadLocal[t][i];
        }
        support[i] = sum;
    }
#else
    // 串行版本
    for (const auto &leaf : tree.adj_list) {
        if (leaf.size() < r) continue;
        
        daf::CliqueSize pivotC = 0, keepC = 0;
        for (const auto &node : leaf) {
            pivotC += node.isPivot;
            keepC += !node.isPivot;
        }
        
        const int needPivot = s - static_cast<int>(keepC);
        if (needPivot < 0 || needPivot > pivotC) continue;
        
        daf::enumerateCombinations(leaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            daf::CliqueSize subNumPivot = 0;
            for (const auto &node : rClique) {
                subNumPivot += node.isPivot;
            }
            
            const int p1 = pivotC - subNumPivot;
            const int p2 = needPivot - subNumPivot;
            
            if (p1 >= 0 && p1 < 1001 && p2 >= 0 && p2 < 401) {
                support[cliqueIndex.byClique(rClique)] += nCr[p1][p2];
            }
            return true;
        });
    }
#endif
    
    return support;
}

// 并行收集受影响的leaves（优化版本）
void collectAffectedLeavesOptimized(
    const std::vector<daf::Size> &removedCliques,
    const StaticCliqueIndex &cliqueIndex,
    const DynamicGraphSet<TreeGraphNode> &treeGraphV,
    std::vector<LeafTask> &tasks,
    std::vector<daf::Size> &globalRemovedCliques,
    std::vector<daf::Size> &leafIndexMap) {
    
    tasks.clear();
    globalRemovedCliques.clear();
    
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
    
    // 第一阶段：并行收集所有(cliqueIdx, leafId)对
    std::vector<std::vector<std::pair<daf::Size, daf::Size>>> threadPairs(nthreads);
    
    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto &localPairs = threadPairs[tid];
        localPairs.reserve(removedCliques.size() * 8);
        
        #pragma omp for schedule(dynamic, 4) nowait
        for (size_t idx = 0; idx < removedCliques.size(); ++idx) {
            const auto cliqueId = removedCliques[idx];
            auto rClique = cliqueIndex.byId(cliqueId);
            
            // 找到所有包含这个r-clique的leaves
            daf::intersect_dense_sets_multi(rClique, const_cast<DynamicGraphSet<TreeGraphNode>&>(treeGraphV).adj_list,
                [&](const TreeGraphNode &node) {
                    localPairs.emplace_back(idx, node.v);
                });
        }
    }
    
    // 第二阶段：合并并排序
    size_t totalSize = 0;
    for (const auto &tp : threadPairs) totalSize += tp.size();
    
    std::vector<std::pair<daf::Size, daf::Size>> allPairs;
    allPairs.reserve(totalSize);
    for (const auto &tp : threadPairs) {
        allPairs.insert(allPairs.end(), tp.begin(), tp.end());
    }
    
    // 按leafId排序以便分组
    std::sort(allPairs.begin(), allPairs.end(),
              [](const auto &a, const auto &b) { 
                  return a.second < b.second || (a.second == b.second && a.first < b.first);
              });
    
    // 第三阶段：构建任务列表
    if (allPairs.empty()) return;
    
    daf::Size currentLeaf = allPairs[0].second;
    daf::Size startIdx = 0;
    
    for (size_t i = 0; i < allPairs.size(); ++i) {
        const auto [cliqueIdx, leafId] = allPairs[i];
        
        if (leafId != currentLeaf) {
            // 创建任务
            LeafTask task;
            task.leafId = currentLeaf;
            task.removedStartIdx = startIdx;
            task.removedCount = globalRemovedCliques.size() - startIdx;
            tasks.push_back(task);
            
            currentLeaf = leafId;
            startIdx = globalRemovedCliques.size();
        }
        
        globalRemovedCliques.push_back(removedCliques[cliqueIdx]);
    }
    
    // 最后一个任务
    LeafTask task;
    task.leafId = currentLeaf;
    task.removedStartIdx = startIdx;
    task.removedCount = globalRemovedCliques.size() - startIdx;
    tasks.push_back(task);
    
#else
    // 串行版本
    std::unordered_map<daf::Size, std::vector<daf::Size>> leafToCliques;
    
    for (auto cliqueId : removedCliques) {
        auto rClique = cliqueIndex.byId(cliqueId);
        
        daf::intersect_dense_sets_multi(rClique, const_cast<DynamicGraphSet<TreeGraphNode>&>(treeGraphV).adj_list,
            [&](const TreeGraphNode &node) {
                leafToCliques[node.v].push_back(cliqueId);
            });
    }
    
    for (const auto &[leafId, cliques] : leafToCliques) {
        LeafTask task;
        task.leafId = leafId;
        task.removedStartIdx = globalRemovedCliques.size();
        task.removedCount = cliques.size();
        tasks.push_back(task);
        
        globalRemovedCliques.insert(globalRemovedCliques.end(), cliques.begin(), cliques.end());
    }
#endif
}

// 处理单个leaf的更新（核心计算函数）
void processLeafUpdateOptimized(
    const LeafTask &task,
    const std::vector<daf::Size> &globalRemovedCliques,
    const DynamicGraph<TreeGraphNode> &tree,
    const StaticCliqueIndex &cliqueIndex,
    daf::CliqueSize r,
    daf::CliqueSize s,
    ThreadLocalWorkspace &workspace) {
    
    const auto &leaf = tree.adj_list[task.leafId];
    if (leaf.empty()) return;
    
    // 计算当前leaf的pivot和keep数量
    daf::CliqueSize pivotC = 0, keepC = 0;
    for (const auto &node : leaf) {
        pivotC += node.isPivot;
        keepC += !node.isPivot;
    }
    const int needPivot = s - static_cast<int>(keepC);
    
    // 收集该leaf需要处理的removed cliques
    workspace.tempCliqueIds.clear();
    for (daf::Size i = 0; i < task.removedCount; ++i) {
        workspace.tempCliqueIds.push_back(globalRemovedCliques[task.removedStartIdx + i]);
    }
    
    // 计算decrements（移除的r-cliques对support的影响）
    for (auto rmCliqueId : workspace.tempCliqueIds) {
        const auto rClique = cliqueIndex.byId(rmCliqueId);
        
        daf::CliqueSize subNumPivot = 0;
        for (const auto &v : rClique) {
            for (const auto &node : leaf) {
                if (node.v == v && node.isPivot) {
                    subNumPivot++;
                    break;
                }
            }
        }
        
        const int p1 = pivotC - subNumPivot;
        const int p2 = needPivot - subNumPivot;
        
        if (p1 >= 0 && p1 < 1001 && p2 >= 0 && p2 < 401) {
            const double ncrValue = nCr[p1][p2];
            workspace.decrements.emplace_back(rmCliqueId, ncrValue);
        }
    }
    
    // 使用Bron-Kerbosch生成新的leaves
    auto mapped = workspace.tempCliqueIds | std::views::transform(
        [&](daf::Size id) { return cliqueIndex.byId(id); });
    
    bkRmClique::removeRClique(leaf, mapped, r, s,
        [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
            auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
            if (newLeaf.size() >= s) {
                workspace.newLeaves.push_back(std::move(newLeaf));
            }
        },
        &workspace.workMap);
    
    // 计算increments（新leaves对support的贡献）
    for (const auto &newLeaf : workspace.newLeaves) {
        daf::CliqueSize newPivotC = 0, newKeepC = 0;
        for (const auto &node : newLeaf) {
            newPivotC += node.isPivot;
            newKeepC += !node.isPivot;
        }
        
        const int newNeedPivot = s - static_cast<int>(newKeepC);
        if (newNeedPivot < 0 || newNeedPivot > newPivotC) continue;
        
        daf::enumerateCombinations(newLeaf, r, [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            daf::CliqueSize subNumPivot = 0;
            for (const auto &node : rClique) {
                subNumPivot += node.isPivot;
            }
            
            const int p1 = newPivotC - subNumPivot;
            const int p2 = newNeedPivot - subNumPivot;
            
            if (p1 >= 0 && p1 < 1001 && p2 >= 0 && p2 < 401) {
                const double ncrValue = nCr[p1][p2];
                const auto cliqueId = cliqueIndex.byClique(rClique);
                workspace.increments.emplace_back(cliqueId, ncrValue);
            }
            return true;
        });
    }
}

// 使用lock-free atomic操作批量应用support更新
void applySupportUpdatesLockFree(
    const std::vector<ThreadLocalWorkspace> &workspaces,
    std::vector<std::atomic<double>> &atomicSupport,
    const std::vector<bool> &removed) {
    
#ifdef _OPENMP
    // 收集所有需要更新的clique IDs
    std::unordered_set<daf::Size> affectedCliques;
    for (const auto &ws : workspaces) {
        for (const auto &[id, _] : ws.increments) {
            if (!removed[id]) affectedCliques.insert(id);
        }
        for (const auto &[id, _] : ws.decrements) {
            if (!removed[id]) affectedCliques.insert(id);
        }
    }
    
    std::vector<daf::Size> affectedVec(affectedCliques.begin(), affectedCliques.end());
    
    // 并行处理每个受影响的clique
    #pragma omp parallel for schedule(dynamic, 64)
    for (size_t i = 0; i < affectedVec.size(); ++i) {
        const daf::Size cliqueId = affectedVec[i];
        double delta = 0.0;
        
        // 收集该clique的所有增量
        for (const auto &ws : workspaces) {
            for (const auto &[id, val] : ws.increments) {
                if (id == cliqueId) delta += val;
            }
            for (const auto &[id, val] : ws.decrements) {
                if (id == cliqueId) delta -= val;
            }
        }
        
        // 原子更新
        if (delta != 0.0) {
            double current = atomicSupport[cliqueId].load(std::memory_order_relaxed);
            double newVal;
            do {
                newVal = current + delta;
            } while (!atomicSupport[cliqueId].compare_exchange_weak(
                current, newVal, std::memory_order_relaxed));
        }
    }
#else
    // 串行版本
    for (const auto &ws : workspaces) {
        for (const auto &[id, val] : ws.increments) {
            if (!removed[id]) {
                double current = atomicSupport[id].load(std::memory_order_relaxed);
                atomicSupport[id].store(current + val, std::memory_order_relaxed);
            }
        }
        for (const auto &[id, val] : ws.decrements) {
            if (!removed[id]) {
                double current = atomicSupport[id].load(std::memory_order_relaxed);
                atomicSupport[id].store(current - val, std::memory_order_relaxed);
            }
        }
    }
#endif
}

// ============================================================================
// 主算法
// ============================================================================

std::vector<std::pair<std::vector<daf::Size>, int>> 
NucleusCoreDecompositionUltraParallel(
    DynamicGraph<TreeGraphNode> &tree,
    const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r,
    daf::CliqueSize s) {
    
    auto timeStart = std::chrono::high_resolution_clock::now();
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Ultra-Parallel Nucleus Decomposition" << std::endl;
    std::cout << "Target: 3x+ speedup on 8 threads" << std::endl;
#ifdef _OPENMP
    std::cout << "OpenMP threads: " << omp_get_max_threads() << std::endl;
#else
    std::cout << "OpenMP: DISABLED (serial mode)" << std::endl;
#endif
    std::cout << "========================================\n" << std::endl;
    
    // Step 1: 构建clique索引
    StaticCliqueIndex cliqueIndex(r);
    auto t1 = std::chrono::high_resolution_clock::now();
    cliqueIndex.build(tree, edgeGraph.adj_list.size());
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Clique index built: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() 
              << " ms" << std::endl;
    
    const daf::Size nClique = cliqueIndex.size();
    const daf::Size graphN = edgeGraph.n;
    
    std::cout << "Total r-cliques: " << nClique << std::endl;
    std::cout << "Graph vertices: " << graphN << std::endl;
    std::cout << "Tree leaves: " << tree.adj_list.size() << std::endl;
    
    // Step 2: 并行计算初始support
    t1 = std::chrono::high_resolution_clock::now();
    auto support = computeInitialSupportOptimized(tree, cliqueIndex, r, s);
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Initial support computed: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() 
              << " ms" << std::endl;
    
    // Step 3: 初始化数据结构
    std::vector<int> coreValues(nClique, 0);
    std::vector<bool> removed(nClique, false);
    
    // 使用atomic的support数组
    std::vector<std::atomic<double>> atomicSupport(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        atomicSupport[i].store(support[i], std::memory_order_relaxed);
    }
    
    // 初始化线程局部工作空间
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif
    std::vector<ThreadLocalWorkspace> workspaces;
    workspaces.reserve(nthreads);
    for (int i = 0; i < nthreads; ++i) {
        workspaces.emplace_back(graphN);
    }
    
    std::cout << "\n=== Starting Peeling Process ===" << std::endl;
    
    int iteration = 0;
    size_t totalRemoved = 0;
    
    std::vector<daf::Size> leafIndexMap(graphN, std::numeric_limits<daf::Size>::max());
    std::vector<LeafTask> tasks;
    std::vector<daf::Size> globalRemovedCliques;
    
    // Step 4: 批量peeling循环
    while (totalRemoved < nClique) {
        auto iterStart = std::chrono::high_resolution_clock::now();
        
        // 4.1 找到最小support值
        double minSupport = std::numeric_limits<double>::max();
        for (daf::Size i = 0; i < nClique; ++i) {
            if (!removed[i]) {
                double val = atomicSupport[i].load(std::memory_order_relaxed);
                if (val < minSupport) {
                    minSupport = val;
                }
            }
        }
        
        if (minSupport == std::numeric_limits<double>::max()) break;
        
        // 4.2 收集所有support == minSupport的r-cliques（批量处理）
        std::vector<daf::Size> currentBatch;
        currentBatch.reserve(nClique / 100);
        
        for (daf::Size i = 0; i < nClique; ++i) {
            if (!removed[i]) {
                double val = atomicSupport[i].load(std::memory_order_relaxed);
                if (std::abs(val - minSupport) < 0.001) {
                    currentBatch.push_back(i);
                    removed[i] = true;
                    coreValues[i] = static_cast<int>(minSupport);
                }
            }
        }
        
        if (currentBatch.empty()) break;
        
        totalRemoved += currentBatch.size();
        
        // 4.3 并行收集受影响的leaves
        auto t3 = std::chrono::high_resolution_clock::now();
        collectAffectedLeavesOptimized(currentBatch, cliqueIndex, treeGraphV,
                                      tasks, globalRemovedCliques, leafIndexMap);
        auto t4 = std::chrono::high_resolution_clock::now();
        
        if (tasks.empty()) {
            if (iteration % 10 == 0) {
                std::cout << "Iter " << iteration 
                          << " | minCore=" << static_cast<int>(minSupport)
                          << " | batch=" << currentBatch.size()
                          << " | NO AFFECTED LEAVES"
                          << " | removed=" << totalRemoved << "/" << nClique
                          << std::endl;
            }
            iteration++;
            continue;
        }
        
        // 4.4 并行处理所有受影响的leaves
        auto t5 = std::chrono::high_resolution_clock::now();
        
#ifdef _OPENMP
        // 清空所有工作空间
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < nthreads; ++i) {
            workspaces[i].clear();
        }
        
        // 并行处理所有leaf任务
        #pragma omp parallel
        {
            const int tid = omp_get_thread_num();
            auto &workspace = workspaces[tid];
            
            #pragma omp for schedule(dynamic, 8) nowait
            for (size_t i = 0; i < tasks.size(); ++i) {
                processLeafUpdateOptimized(tasks[i], globalRemovedCliques,
                                          tree, cliqueIndex, r, s, workspace);
            }
        }
#else
        // 串行版本
        workspaces[0].clear();
        for (const auto &task : tasks) {
            processLeafUpdateOptimized(task, globalRemovedCliques,
                                      tree, cliqueIndex, r, s, workspaces[0]);
        }
#endif
        
        auto t6 = std::chrono::high_resolution_clock::now();
        
        // 4.5 批量应用support更新（lock-free）
        applySupportUpdatesLockFree(workspaces, atomicSupport, removed);
        
        auto t7 = std::chrono::high_resolution_clock::now();
        
        auto iterEnd = std::chrono::high_resolution_clock::now();
        
        if (iteration % 10 == 0 || currentBatch.size() > 100) {
            auto collectTime = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count();
            auto processTime = std::chrono::duration_cast<std::chrono::microseconds>(t6 - t5).count();
            auto updateTime = std::chrono::duration_cast<std::chrono::microseconds>(t7 - t6).count();
            auto totalIterTime = std::chrono::duration_cast<std::chrono::milliseconds>(iterEnd - iterStart).count();
            
            std::cout << "Iter " << iteration 
                      << " | minCore=" << static_cast<int>(minSupport)
                      << " | batch=" << currentBatch.size()
                      << " | leaves=" << tasks.size()
                      << " | removed=" << totalRemoved << "/" << nClique
                      << " | time=" << totalIterTime << "ms"
                      << " (collect=" << collectTime << "us, process=" << processTime << "us, update=" << updateTime << "us)"
                      << std::endl;
        }
        
        iteration++;
    }
    
    auto timeEnd = std::chrono::high_resolution_clock::now();
    auto totalTime = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
    
    std::cout << "\n=== Peeling Complete ===" << std::endl;
    std::cout << "Total iterations: " << iteration << std::endl;
    std::cout << "Total time: " << totalTime << " ms" << std::endl;
    
    // Step 5: 构造返回结果
    std::vector<std::pair<std::vector<daf::Size>, int>> result;
    result.reserve(nClique);
    
    for (daf::Size i = 0; i < nClique; ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueVec(clique.begin(), clique.end());
        result.emplace_back(std::move(cliqueVec), coreValues[i]);
    }
    
    return result;
}

} // namespace UltraParallel
