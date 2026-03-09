/*
 * 全新的 Nucleus Decomposition 算法
 * 
 * 核心思想：
 * 1. Index-based: 使用整数 ID 代替 clique 内容比较
 * 2. Batch processing: 批量处理同一 core 值的所有 r-cliques
 * 3. Lock-free parallel: 使用分区和无锁数据结构
 * 4. Cache-friendly: 优化内存访问模式
 * 
 * 与原算法的区别：
 * - 原算法：逐个处理 r-clique，每次都要查找和更新
 * - 新算法：批量处理，预先分区，减少同步开销
 */

#include "../NucleusDecomposition/NCliqueCoreDecomposition.h"
#include <vector>
#include <algorithm>
#include <atomic>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

extern double nCr[1001][401];

namespace NucleusNew {

// ============================================================================
// 数据结构：紧凑的 r-clique 表示
// ============================================================================

struct RCliqueCompact {
    uint32_t id;           // r-clique ID
    uint32_t leafId;       // 所属的 leaf ID
    float support;         // support 值（使用 float 节省内存）
    uint16_t core;         // core 值
    uint8_t size;          // clique 大小
    uint8_t flags;         // 标志位
    
    // 顶点列表存储在外部数组中，通过 offset 访问
    uint32_t vertexOffset; // 顶点列表的起始位置
};

// ============================================================================
// 核心数据结构：分区的 r-clique 存储
// ============================================================================

class PartitionedRCliqueStore {
public:
    static constexpr int NUM_PARTITIONS = 64; // 分区数量
    
    struct Partition {
        std::vector<RCliqueCompact> cliques;
        std::vector<uint32_t> vertices;  // 所有 cliques 的顶点
        std::atomic<int> lock{0};        // 细粒度锁
        
        void reserve(size_t n) {
            cliques.reserve(n);
            vertices.reserve(n * 10); // 估计每个 clique 10 个顶点
        }
    };
    
    std::vector<Partition> partitions;
    
    PartitionedRCliqueStore() : partitions(NUM_PARTITIONS) {}
    
    // 根据 clique ID 获取分区
    int getPartition(uint32_t cliqueId) const {
        return cliqueId % NUM_PARTITIONS;
    }
    
    // 添加 r-clique（线程安全）
    void addClique(uint32_t id, uint32_t leafId, 
                   const std::vector<uint32_t>& verts, float support) {
        int p = getPartition(id);
        auto& part = partitions[p];
        
        // 简单的自旋锁
        while (part.lock.exchange(1, std::memory_order_acquire) == 1) {
            // spin
        }
        
        RCliqueCompact clique;
        clique.id = id;
        clique.leafId = leafId;
        clique.support = support;
        clique.core = 0;
        clique.size = verts.size();
        clique.flags = 0;
        clique.vertexOffset = part.vertices.size();
        
        part.cliques.push_back(clique);
        part.vertices.insert(part.vertices.end(), verts.begin(), verts.end());
        
        part.lock.store(0, std::memory_order_release);
    }
    
    // 获取 clique 的顶点
    std::vector<uint32_t> getVertices(int partition, size_t cliqueIdx) const {
        const auto& part = partitions[partition];
        const auto& clique = part.cliques[cliqueIdx];
        
        auto begin = part.vertices.begin() + clique.vertexOffset;
        auto end = begin + clique.size;
        
        return std::vector<uint32_t>(begin, end);
    }
};

// ============================================================================
// 核心算法：批量 peeling
// ============================================================================

class BatchPeelingAlgorithm {
public:
    BatchPeelingAlgorithm(int numThreads = 16) : numThreads_(numThreads) {}
    
    // 主函数：计算 nucleus decomposition
    std::vector<uint16_t> compute(
        const DynamicGraph<TreeGraphNode>& tree,
        const Graph& edgeGraph,
        int r, int s) {
        
        std::cout << "Starting new batch peeling algorithm..." << std::endl;
        std::cout << "Threads: " << numThreads_ << std::endl;
        
        // 步骤 1: 构建紧凑的 r-clique 存储
        auto store = buildCompactStore(tree, edgeGraph, r, s);
        
        // 步骤 2: 计算初始 support
        computeInitialSupport(store, tree, r, s);
        
        // 步骤 3: 批量 peeling
        auto cores = batchPeeling(store, tree, edgeGraph, r, s);
        
        return cores;
    }
    
private:
    int numThreads_;
    
    // 构建紧凑存储
    PartitionedRCliqueStore buildCompactStore(
        const DynamicGraph<TreeGraphNode>& tree,
        const Graph& edgeGraph,
        int r, int s) {
        
        std::cout << "Building compact store..." << std::endl;
        
        PartitionedRCliqueStore store;
        
        // 预分配空间
        size_t totalLeaves = tree.adj_list.size();
        for (auto& part : store.partitions) {
            part.reserve(totalLeaves / store.NUM_PARTITIONS + 1000);
        }
        
        // 并行构建
        std::atomic<uint32_t> nextId{0};
        
        #pragma omp parallel for schedule(dynamic, 64) num_threads(numThreads_)
        for (size_t leafIdx = 0; leafIdx < totalLeaves; ++leafIdx) {
            const auto& leaf = tree.adj_list[leafIdx];
            
            if (leaf.size() < r) continue;
            
            // 枚举所有 r-cliques
            enumerateRCliques(leaf, r, [&](const std::vector<uint32_t>& clique) {
                uint32_t id = nextId.fetch_add(1, std::memory_order_relaxed);
                store.addClique(id, leafIdx, clique, 0.0f);
            });
        }
        
        std::cout << "Total r-cliques: " << nextId.load() << std::endl;
        
        return store;
    }
    
    // 枚举 r-cliques
    template<typename Func>
    void enumerateRCliques(const std::vector<TreeGraphNode>& leaf, 
                          int r, Func&& callback) {
        std::vector<uint32_t> current;
        current.reserve(r);
        
        std::function<void(size_t, int)> enumerate = 
            [&](size_t start, int remaining) {
            if (remaining == 0) {
                callback(current);
                return;
            }
            
            for (size_t i = start; i <= leaf.size() - remaining; ++i) {
                current.push_back(leaf[i].v);
                enumerate(i + 1, remaining - 1);
                current.pop_back();
            }
        };
        
        enumerate(0, r);
    }
    
    // 计算初始 support
    void computeInitialSupport(
        PartitionedRCliqueStore& store,
        const DynamicGraph<TreeGraphNode>& tree,
        int r, int s) {
        
        std::cout << "Computing initial support..." << std::endl;
        
        // 并行计算每个分区
        #pragma omp parallel for schedule(static) num_threads(numThreads_)
        for (int p = 0; p < store.NUM_PARTITIONS; ++p) {
            auto& part = store.partitions[p];
            
            for (auto& clique : part.cliques) {
                // 计算 support（简化版本）
                const auto& leaf = tree.adj_list[clique.leafId];
                
                int pivotC = 0, keepC = 0;
                for (const auto& node : leaf) {
                    if (node.isPivot) pivotC++;
                    else keepC++;
                }
                
                int needPivot = s - keepC;
                
                // 计算这个 r-clique 的 pivot 数量
                auto verts = store.getVertices(p, &clique - &part.cliques[0]);
                int subNumPivot = 0;
                for (uint32_t v : verts) {
                    // 检查是否是 pivot（需要从 leaf 中查找）
                    for (const auto& node : leaf) {
                        if (node.v == v && node.isPivot) {
                            subNumPivot++;
                            break;
                        }
                    }
                }
                
                if (subNumPivot <= needPivot && 
                    pivotC - subNumPivot >= 0 &&
                    pivotC - subNumPivot < 1001 &&
                    needPivot - subNumPivot >= 0 &&
                    needPivot - subNumPivot < 401) {
                    clique.support = nCr[pivotC - subNumPivot][needPivot - subNumPivot];
                } else {
                    clique.support = 0.0f;
                }
            }
        }
    }
    
    // 批量 peeling
    std::vector<uint16_t> batchPeeling(
        PartitionedRCliqueStore& store,
        const DynamicGraph<TreeGraphNode>& tree,
        const Graph& edgeGraph,
        int r, int s) {
        
        std::cout << "Starting batch peeling..." << std::endl;
        
        // 收集所有 cliques 到一个数组
        std::vector<std::pair<float, uint32_t>> supportAndId;
        
        for (int p = 0; p < store.NUM_PARTITIONS; ++p) {
            const auto& part = store.partitions[p];
            for (size_t i = 0; i < part.cliques.size(); ++i) {
                supportAndId.emplace_back(part.cliques[i].support, part.cliques[i].id);
            }
        }
        
        // 排序
        std::sort(supportAndId.begin(), supportAndId.end());
        
        std::cout << "Total cliques to process: " << supportAndId.size() << std::endl;
        
        // 批量处理
        std::vector<uint16_t> cores(supportAndId.size());
        uint16_t currentCore = 0;
        
        size_t processed = 0;
        while (processed < supportAndId.size()) {
            float minSupport = supportAndId[processed].first;
            currentCore = std::max(currentCore, (uint16_t)minSupport);
            
            // 找到所有 support <= currentCore 的 cliques
            size_t batchEnd = processed;
            while (batchEnd < supportAndId.size() && 
                   supportAndId[batchEnd].first <= currentCore) {
                batchEnd++;
            }
            
            size_t batchSize = batchEnd - processed;
            std::cout << "Processing batch: " << batchSize << " cliques at core " 
                      << currentCore << std::endl;
            
            // 批量处理这些 cliques
            #pragma omp parallel for schedule(dynamic, 64) num_threads(numThreads_)
            for (size_t i = processed; i < batchEnd; ++i) {
                uint32_t cliqueId = supportAndId[i].second;
                cores[cliqueId] = currentCore;
            }
            
            processed = batchEnd;
            
            // TODO: 更新剩余 cliques 的 support
            // 这部分需要实现 BK 算法来重新计算
        }
        
        return cores;
    }
};

} // namespace NucleusNew

// ============================================================================
// 对外接口：与原算法兼容
// ============================================================================

std::vector<std::pair<std::vector<daf::Size>, int>> 
NucleusCoreDecompositionNew(
    DynamicGraph<TreeGraphNode>& tree,
    const Graph& edgeGraph,
    DynamicGraphSet<TreeGraphNode>& treeGraphV,
    daf::CliqueSize r,
    daf::CliqueSize s) {
    
    // 使用新算法
    NucleusNew::BatchPeelingAlgorithm algo(16); // 16 threads
    auto cores = algo.compute(tree, edgeGraph, r, s);
    
    // 转换输出格式（与原算法兼容）
    std::vector<std::pair<std::vector<daf::Size>, int>> result;
    
    // TODO: 转换格式
    
    return result;
}
