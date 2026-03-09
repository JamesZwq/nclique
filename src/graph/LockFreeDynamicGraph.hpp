// Lock-Free Dynamic Graph for Parallel Nucleus Decomposition
// 使用无锁数据结构实现完全并行的 tree 操作

#pragma once

#include <atomic>
#include <vector>
#include <memory>
#include <cstdint>
#include "../Global/Global.h"

namespace LockFree {

// 无锁的节点结构
template<typename NodeType>
struct LockFreeNode {
    std::vector<NodeType> data;
    std::atomic<bool> valid;  // 节点是否有效
    std::atomic<uint64_t> version;  // 版本号，用于检测并发修改
    
    LockFreeNode() : valid(true), version(0) {}
    
    LockFreeNode(const std::vector<NodeType>& d) 
        : data(d), valid(true), version(0) {}
};

// 无锁的动态图
template<typename NodeType>
class LockFreeDynamicGraph {
private:
    std::vector<std::unique_ptr<LockFreeNode<NodeType>>> nodes;
    std::atomic<daf::Size> nodeCount;
    std::atomic<daf::Size> capacity;
    
    // 预分配的节点池，避免频繁分配
    static constexpr daf::Size INITIAL_CAPACITY = 1000000;
    
public:
    LockFreeDynamicGraph() : nodeCount(0), capacity(INITIAL_CAPACITY) {
        nodes.resize(INITIAL_CAPACITY);
        for (daf::Size i = 0; i < INITIAL_CAPACITY; ++i) {
            nodes[i] = std::make_unique<LockFreeNode<NodeType>>();
            nodes[i]->valid.store(false, std::memory_order_relaxed);
        }
    }
    
    // 并行安全的添加节点
    daf::Size addNode(const std::vector<NodeType>& data) {
        daf::Size id = nodeCount.fetch_add(1, std::memory_order_relaxed);
        
        // 如果超出容量，扩展（使用 critical section，但很少发生）
        if (id >= capacity.load(std::memory_order_relaxed)) {
            #pragma omp critical(expand_capacity)
            {
                daf::Size currentCap = capacity.load(std::memory_order_relaxed);
                if (id >= currentCap) {
                    daf::Size newCap = currentCap * 2;
                    nodes.resize(newCap);
                    for (daf::Size i = currentCap; i < newCap; ++i) {
                        nodes[i] = std::make_unique<LockFreeNode<NodeType>>();
                        nodes[i]->valid.store(false, std::memory_order_relaxed);
                    }
                    capacity.store(newCap, std::memory_order_release);
                }
            }
        }
        
        // 设置节点数据
        nodes[id]->data = data;
        nodes[id]->valid.store(true, std::memory_order_release);
        nodes[id]->version.fetch_add(1, std::memory_order_relaxed);
        
        return id;
    }
    
    // 并行安全的移除节点（标记为无效）
    void removeNode(daf::Size id) {
        if (id < capacity.load(std::memory_order_relaxed)) {
            nodes[id]->valid.store(false, std::memory_order_release);
            nodes[id]->version.fetch_add(1, std::memory_order_relaxed);
        }
    }
    
    // 并行安全的读取节点
    bool getNode(daf::Size id, std::vector<NodeType>& out) const {
        if (id >= capacity.load(std::memory_order_acquire)) {
            return false;
        }
        
        uint64_t v1 = nodes[id]->version.load(std::memory_order_acquire);
        bool valid = nodes[id]->valid.load(std::memory_order_acquire);
        
        if (!valid) return false;
        
        out = nodes[id]->data;
        
        uint64_t v2 = nodes[id]->version.load(std::memory_order_acquire);
        
        // 如果版本号改变，说明有并发修改，重试
        if (v1 != v2) {
            return getNode(id, out);
        }
        
        return true;
    }
    
    // 获取节点数量
    daf::Size size() const {
        return nodeCount.load(std::memory_order_relaxed);
    }
    
    // 克隆（用于测试）
    LockFreeDynamicGraph clone() const {
        LockFreeDynamicGraph result;
        daf::Size count = nodeCount.load(std::memory_order_relaxed);
        
        for (daf::Size i = 0; i < count; ++i) {
            if (nodes[i]->valid.load(std::memory_order_relaxed)) {
                result.addNode(nodes[i]->data);
            }
        }
        
        return result;
    }
};

} // namespace LockFree
