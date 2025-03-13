//
// Created by 张文谦 on 25-3-4.
//

#ifndef NCLIQUECOREDECOMPOSITION_H
#define NCLIQUECOREDECOMPOSITION_H
#include "MultiBranchTree.h"
#include <tbb/spin_mutex.h>
#include <map>

struct ThreadSafeMap {
    std::map<daf::Size, daf::Size, std::greater<>> map;
    tbb::spin_mutex mutex;

    void insert(daf::Size key, daf::Size value) {
        tbb::spin_mutex::scoped_lock lock(mutex);
        map[key] = value;
    }

    daf::Size add(daf::Size key, daf::Size value) {
        tbb::spin_mutex::scoped_lock lock(mutex);
        map[key] += value;
        return map[key];
    }

    daf::Size get(daf::Size key) {
        tbb::spin_mutex::scoped_lock lock(mutex);
        return map[key];
    }

    friend std::ostream &operator<<(std::ostream &os, const ThreadSafeMap &map) {
        os << "{";
        bool first = true;
        for (const auto &pair : map.map) {
            if (!first) {
                os << ", ";
            }
            os << pair.first << ": " << pair.second;
            first = false;
        }
        os << "}";
        return os;
    }
};

void baseNucleusCoreDecomposition(const MultiBranchTree &tree, daf::CliqueSize k);
void baseNucleusCoreDecompositionPar(const MultiBranchTree &tree, daf::CliqueSize k);
void baseNucleusCoreDecompositionParHash(const MultiBranchTree &tree, daf::CliqueSize k);


#endif //NCLIQUECOREDECOMPOSITION_H
