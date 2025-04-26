//
// Created by 张文谦 on 25-4-8.
//

#ifndef CLIQUEHASHMAP_H
#define CLIQUEHASHMAP_H
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <graph/Graph.h>

#include "../Global/Global.h"


struct Clique {
    // std::vector<daf::Size> vertices;
    // daf::Size *vertices;
    daf::StaticVector<daf::Size> vertices;
    std::size_t hashValue;
    // 构造函数接收一个顶点列表，并排序以保证顺序一致性
    explicit Clique(const daf::StaticVector<daf::Size> &vs) : vertices(vs) {
        std::ranges::sort(vertices);
        hashValue = compute_hash();
    }

    // 重载 == 运算符，便于 unordered_map 判断两个 clique 是否相等
    bool operator==(const Clique &other) const {
        return vertices == other.vertices;
    }

private:
    // 使用与之前相同的混合技术计算哈希值
    [[nodiscard]] std::size_t compute_hash() const {
        std::size_t seed = 0;
        for (const daf::Size &v: vertices) {
            seed ^= std::hash<daf::Size>()(v) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


template<typename T>
class EdgeHashMap {
public:
    // 用来把 (u,v) 编码成一个 64 位整数：高 32 位放小的节点，低 32 位放大的节点
    static inline std::uint64_t encode(std::size_t u, std::size_t v) noexcept {
        if (u > v) std::swap(u, v);
        return (std::uint64_t(u) << 32) | std::uint64_t(v);
    }

    static inline std::pair<std::size_t, std::size_t> decode(std::uint64_t key) noexcept {
        std::size_t u = key >> 32;
        std::size_t v = key & 0xFFFFFFFFULL;
        return {u, v};
    }

    // 恒等哈希：key 本身就是 64 位整数
    struct IdentityHash {
        std::size_t operator()(std::uint64_t x) const noexcept {
            return std::size_t(x);
        }
    };

    // 完全去掉 edgeCliqueCount，只保留 map
    explicit EdgeHashMap(std::size_t reserve_size = 0) {
        map.reserve(reserve_size);
    }

    // operator[]：直接编码，并返回 T&，不存在时会默认构造
    T& operator[](const std::pair<std::size_t, std::size_t>& edge) {
        auto key = encode(edge.first, edge.second);
        return map[key];
    }

    const T& operator[](const std::pair<std::size_t, std::size_t>& edge) const {
        auto key = encode(edge.first, edge.second);
        auto it = map.find(key);
        if (it == map.end()) {
            throw std::out_of_range("Edge not found");
        }
        return it->second;
    }

    T& get(const std::size_t u, const std::size_t v) {
        auto key = encode(u, v);
        return map[key];
    }

    // 如果你想只查不插，可以加个 find 接口
    T* find(const std::pair<std::size_t, std::size_t>& edge) {
        auto key = encode(edge.first, edge.second);
        auto it = map.find(key);
        return it == map.end() ? nullptr : &it->second;
    }

    friend std::ostream& operator<<(std::ostream& os, const EdgeHashMap& ehm) {
        os << "{";
        bool first = true;
        for (const auto& pair : ehm.map) {
            if (!first) {
                os << ", ";
            }
            os << decode(pair.first) << ": " << pair.second;
            first = false;
        }
        os << "}";
        return os;
    }

    void print() const {
        for (const auto& pair : map) {
            std::cout << decode(pair.first) << " : " << pair.second << std::endl;
        }
    }

private:
    std::unordered_map<std::uint64_t, T, IdentityHash> map;
};

// 自定义哈希函数，采用类似 boost::hash_combine 的方法
template<>
struct std::hash<Clique> {
    std::size_t operator()(const Clique &c) const noexcept {
        return c.hashValue;
    }
};
#endif //CLIQUEHASHMAP_H
