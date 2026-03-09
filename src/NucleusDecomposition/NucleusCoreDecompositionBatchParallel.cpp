// 全新设计：批量并行 Nucleus Decomposition
// 核心思想：批量收集 + 完全并行处理 + 延迟更新

#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "BK/BronKerboschRmRClique.hpp"
#include <boost/heap/d_ary_heap.hpp>
#include <atomic>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

extern double nCr[1001][401];

namespace BatchParallel {

// 批量更新结构
struct BatchUpdate {
    std::vector<daf::Size> leavesToRemove;
    std::vector<std::vector<TreeGraphNode>> leavesToAdd;
    std::vector<std::pair<daf::Size, double>> supportIncr;
    std::vector<std::pair<daf::Size, double>> supportDecr;
};

std::vector<std::pair<std::vector<daf::Size>, int>> NucleusCoreDecompositionRCliqueBatchParallel(
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
    
    // Step 3: 使用 atomic counting
    std::vector<std::atomic<double>> atomicCounting(nClique);
#pragma omp parallel for schedule(static)
    for (daf::Size i = 0; i < nClique; ++i) {
        atomicCounting[i].store(countingRClique[i], std::memory_order_relaxed);
    }
    
    std::vector<int> coreRClique(nClique, 0);
    std::vector<std::atomic<bool>> rCliqueInHeap(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        rCliqueInHeap[i].store(true, std::memory_order_relaxed);
    }
    
    // Heap 结构
    using HeapNode = std::pair<double, daf::Size>;
    using Heap = boost::heap::d_ary_heap<HeapNode, boost::heap::arity<4>, 
                                         boost::heap::mutable_<true>,
                                         boost::heap::compare<std::greater<HeapNode>>>;
    Heap heap;
    std::vector<typename Heap::handle_type> heapHandles(nClique);
    
    for (daf::Size i = 0; i < nClique; ++i) {
        heapHandles[i] = heap.push({countingRClique[i], i});
    }
    
    std::cout << "=========================begin=========================" << std::endl;
    std::cout << "Using Batch Parallel Algorithm (Aggressive Optimization)" << std::endl;
    
    std::vector<daf::Size> changedLeafIndex(graphN, std::numeric_limits<daf::Size>::max());
    int iterationCount = 0;
    
    while (!heap.empty()) {
        auto [minCore, minId] = heap.top();
        heap.pop();
        
        if (!rCliqueInHeap[minId].load(std::memory_order_relaxed)) continue;
        
        coreRClique[minId] = std::max(static_cast<int>(minCore), 
                                      static_cast<int>(atomicCounting[minId].load(std::memory_order_relaxed)));
        rCliqueInHeap[minId].store(false, std::memory_order_relaxed);
        
        // 批量收集相同 core 值的 r-cliques
        std::vector<daf::Size> currentRemoveRcliqueIds;
        currentRemoveRcliqueIds.push_back(minId);
        
        while (!heap.empty() && std::abs(heap.top().first - minCore) < 0.001) {
            auto [_, id] = heap.top();
            heap.pop();
            if (rCliqueInHeap[id].load(std::memory_order_relaxed)) {
                coreRClique[id] = std::max(static_cast<int>(minCore),
                                          static_cast<int>(atomicCounting[id].load(std::memory_order_relaxed)));
                rCliqueInHeap[id].store(false, std::memory_order_relaxed);
                currentRemoveRcliqueIds.push_back(id);
            }
        }
        
        if (iterationCount % 100 == 0) {
            std::cout << "minCore: " << static_cast<int>(minCore) 
                      << " heap size: " << heap.size() 
                      << " batch size: " << currentRemoveRcliqueIds.size() << std::endl;
        }
        iterationCount++;
        
        // 并行收集受影响的 leaves
        std::vector<daf::Size> changedLeaf;
        std::vector<std::vector<daf::Size>> removedRCliqueIdForLeaf;
        
#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
        std::vector<std::vector<std::pair<daf::Size, daf::Size>>> threadPairs(nthreads);
        
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            auto &localPairs = threadPairs[tid];
            
#pragma omp for schedule(dynamic, 4)
            for (size_t idx = 0; idx < currentRemoveRcliqueIds.size(); ++idx) {
                auto rmCliqueId = currentRemoveRcliqueIds[idx];
                auto rClique = cliqueIndex.byId(rmCliqueId);
                daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
                    [&](const TreeGraphNode &uClique) {
                        localPairs.emplace_back(idx, uClique.v);
                    });
            }
        }
        
        std::vector<std::pair<daf::Size, daf::Size>> allPairs;
        size_t totalSize = 0;
        for (const auto &tp : threadPairs) totalSize += tp.size();
        allPairs.reserve(totalSize);
        for (const auto &tp : threadPairs) {
            allPairs.insert(allPairs.end(), tp.begin(), tp.end());
        }
        std::sort(allPairs.begin(), allPairs.end(), 
                  [](const auto &a, const auto &b) { return a.second < b.second; });
        
        for (const auto &p : allPairs) {
            daf::Size origIdx = p.first;
            daf::Size leafV = p.second;
            auto rmCliqueId = currentRemoveRcliqueIds[origIdx];
            auto &leafId = changedLeafIndex[leafV];
            if (leafId == std::numeric_limits<daf::Size>::max()) {
                leafId = removedRCliqueIdForLeaf.size();
                removedRCliqueIdForLeaf.emplace_back();
                changedLeaf.push_back(leafV);
            }
            removedRCliqueIdForLeaf[leafId].emplace_back(rmCliqueId);
        }
#endif
        
        // 完全并行处理所有 leaves
        struct LeafResult {
            std::vector<std::vector<TreeGraphNode>> newLeaves;
            std::vector<std::pair<daf::Size, double>> incr;
            std::vector<std::pair<daf::Size, double>> decr;
        };
        std::vector<LeafResult> leafResults(changedLeaf.size());
        
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
                daf::CliqueSize pivotC = 0, keepC = 0;
                for (const auto &node : leaf) {
                    if (node.isPivot) pivotC++; else keepC++;
                }
                
                for (auto rmId : removedR) {
                    auto rClique = cliqueIndex.byId(rmId);
                    daf::CliqueSize subNumPovit = 0;
                    for (const auto &v : rClique) {
                        for (const auto &node : leaf) {
                            if (node.v == v && node.isPivot) {
                                subNumPovit++;
                                break;
                            }
                        }
                    }
                    int needPivot = s - static_cast<int>(keepC);
                    if (pivotC >= subNumPovit && needPivot >= subNumPovit &&
                        pivotC - subNumPovit < 1001 && needPivot - subNumPovit < 401) {
                        double ncrValue = nCr[pivotC - subNumPovit][needPivot - subNumPovit];
                        leafResults[idx].decr.emplace_back(rmId, ncrValue);
                    }
                }
                
                // 使用 Bron-Kerbosch 生成新 leaves
                auto mapped = removedR | std::views::transform(
                    [&](const daf::Size id) { return cliqueIndex.byId(id); });
                
                bkRmClique::removeRClique(leaf, mapped, r, s,
                    [&](const bkRmClique::Bitset &c, const bkRmClique::Bitset &pivots) {
                        auto newLeaf = bkRmClique::coverToVertex(c, pivots, leaf);
                        if (newLeaf.size() >= s) {
                            leafResults[idx].newLeaves.push_back(std::move(newLeaf));
                        }
                    },
                    &threadMap);
                
                // 计算 increments
                for (const auto &newLeaf : leafResults[idx].newLeaves) {
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
                            leafResults[idx].incr.emplace_back(cliqueId, ncrValue);
                        }
                        return true;
                    });
                }
            }
        }
        
        // 并行更新 support（使用 atomic）
#pragma omp parallel for schedule(dynamic, 16)
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            for (const auto &p : leafResults[idx].incr) {
                union { double d; uint64_t u; } old_val, new_val;
                std::atomic<uint64_t> *atomic_ptr = 
                    reinterpret_cast<std::atomic<uint64_t>*>(&atomicCounting[p.first]);
                old_val.u = atomic_ptr->load(std::memory_order_relaxed);
                do {
                    new_val.d = old_val.d + p.second;
                } while (!atomic_ptr->compare_exchange_weak(
                    old_val.u, new_val.u, std::memory_order_relaxed));
            }
            for (const auto &p : leafResults[idx].decr) {
                union { double d; uint64_t u; } old_val, new_val;
                std::atomic<uint64_t> *atomic_ptr = 
                    reinterpret_cast<std::atomic<uint64_t>*>(&atomicCounting[p.first]);
                old_val.u = atomic_ptr->load(std::memory_order_relaxed);
                do {
                    new_val.d = old_val.d - p.second;
                } while (!atomic_ptr->compare_exchange_weak(
                    old_val.u, new_val.u, std::memory_order_relaxed));
            }
        }
        
        // 串行更新 tree/treeGraphV（保持安全）
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            auto leafId = changedLeaf[idx];
            const auto leaf = tree.adj_list[leafId];
            
            for (auto leafV : leaf) {
                if (leafV.isPivot) treeGraphV.removeNbr(leafV.v, {leafId, true});
                else treeGraphV.removeNbr(leafV.v, {leafId, false});
            }
            
            for (auto &newLeaf : leafResults[idx].newLeaves) {
                auto newId = tree.addNode(newLeaf);
                const auto &stored = tree.adj_list[newId];
                for (auto i : stored) {
                    if (i.isPivot) treeGraphV.addNbr(i.v, {newId, true});
                    else treeGraphV.addNbr(i.v, {newId, false});
                }
                if (newId >= changedLeafIndex.size())
                    changedLeafIndex.resize(newId * 2, std::numeric_limits<daf::Size>::max());
            }
            
            tree.removeNode(leafId);
        }
        
        // 更新 heap
        std::unordered_set<daf::Size> needUpdate;
        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            for (const auto &p : leafResults[idx].decr) {
                needUpdate.insert(p.first);
            }
        }
        for (auto cliqueId : needUpdate) {
            if (rCliqueInHeap[cliqueId].load(std::memory_order_relaxed)) {
                countingRClique[cliqueId] = atomicCounting[cliqueId].load(std::memory_order_relaxed);
                heap.update(heapHandles[cliqueId], {countingRClique[cliqueId], cliqueId});
            }
        }
        
        // 清理
        for (auto leafV : changedLeaf) {
            changedLeafIndex[leafV] = std::numeric_limits<daf::Size>::max();
        }
    }
    
    std::cout << "Total iterations: " << iterationCount << std::endl;
    
    // 构造返回结果
    std::vector<std::pair<std::vector<daf::Size>, int>> result;
    result.reserve(nClique);
    for (daf::Size i = 0; i < nClique; ++i) {
        auto clique = cliqueIndex.byId(i);
        std::vector<daf::Size> cliqueVec(clique.begin(), clique.end());
        result.emplace_back(std::move(cliqueVec), coreRClique[i]);
    }
    
    auto time_end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
    std::cout << "Batch parallel algorithm took: " << elapsed << " ms" << std::endl;
    
    return result;
}

} // namespace BatchParallel
