//
// Created by 张文谦 on 25-5-2.
//

#ifndef EDGESET_H
#define EDGESET_H
/***********************************************************
 *  EdgeSet.h  ——  轻量级无向边集合
 *  C++17 单头文件；直接 #include "EdgeSet.h" 即可
 ***********************************************************/
#pragma once
#include <set>
#include <string>
#include <fstream>
#include <stdexcept>
#include <utility>

class EdgeSet {
    using Edge = std::pair<int,int>;
    std::set<Edge> edges_;

    static Edge normalize(int u, int v) noexcept {
        return (u < v) ? Edge(u,v) : Edge(v,u);
    }

public:
    /* ---------- 构造与加载 ---------- */


    /// 直接从文件加载
    explicit EdgeSet(const std::string& path  = "/Users/zhangwenqian/UNSW/pivoter/a.txt")            { load(path); }

    /// 重新加载（会清空原内容）
    void load(const std::string& path) {
        std::ifstream fin(path);
        if (!fin) throw std::runtime_error("Cannot open " + path);

        edges_.clear();

        int u, v;
        while (fin >> u >> v) {
            if (u == v) continue; // 自环
            if (u > v) std::swap(u, v);
            insert(u, v);
        }
    }

    /* ---------- 基本操作 ---------- */
    bool insert(int u, int v) {
        return edges_.insert(normalize(u, v)).second;
    }
    bool erase(int u, int v) {
        return edges_.erase(normalize(u, v)) > 0;
    }
    bool contains(int u, int v) const {
        return edges_.count(normalize(u, v)) != 0;
    }

    bool check(int u, int v) const {
        if (u > v) std::swap(u, v);
        if (!contains(u, v)) return false;
        if (edges_.count(normalize(u, v)) != 0) {
            std::cout << "Edge (" << u << ", " << v << ") exists." << std::endl;
        }
        return true;
    }

    std::size_t size()   const noexcept { return edges_.size(); }
    bool         empty() const noexcept { return edges_.empty(); }

    /* ---------- 迭代访问 ---------- */
    auto begin() const noexcept { return edges_.begin(); }
    auto end()   const noexcept { return edges_.end();   }

    /* ---------- 导出 ---------- */
    void save(const std::string& path, long nVert = 0) const {
        std::ofstream fout(path);
        if (!fout) throw std::runtime_error("Cannot open " + path);

        fout << nVert << ' ' << edges_.size() << '\n';
        for (auto [u, v] : edges_) fout << u << ' ' << v << '\n';
    }
};

#endif //EDGESET_H
