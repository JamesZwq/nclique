//
// Created by _ on 25-5-2.
//

#ifndef EDGESET_H
#define EDGESET_H
/***********************************************************
 *  EdgeSet.h  ——  
 *  C++17 ； #include "EdgeSet.h" 
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
    /* ----------  ---------- */


    /// 
    explicit EdgeSet(const std::string& path  = "~/_/pivoter/a.txt")            { load(path); }

    /// （）
    void load(const std::string& path) {
        std::ifstream fin(path);
        if (!fin) throw std::runtime_error("Cannot open " + path);

        edges_.clear();

        int u, v;
        while (fin >> u >> v) {
            if (u == v) continue; // 
            if (u > v) std::swap(u, v);
            insert(u, v);
        }
    }

    /* ----------  ---------- */
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

    /* ----------  ---------- */
    auto begin() const noexcept { return edges_.begin(); }
    auto end()   const noexcept { return edges_.end();   }

    /* ----------  ---------- */
    void save(const std::string& path, long nVert = 0) const {
        std::ofstream fout(path);
        if (!fout) throw std::runtime_error("Cannot open " + path);

        fout << nVert << ' ' << edges_.size() << '\n';
        for (auto [u, v] : edges_) fout << u << ' ' << v << '\n';
    }
};

#endif //EDGESET_H
