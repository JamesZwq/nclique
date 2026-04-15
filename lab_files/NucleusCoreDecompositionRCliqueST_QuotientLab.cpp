//
// Quotient Lab for r>=3:
// Measure clean-leaf quotient compression potential, then run exact V20.
//

#include "NCliqueCoreDecomposition.h"
#include "../BK/BronKerboschRmRClique.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

extern double nCr[1001][401];

namespace {

struct ExplicitLeafEntryEstimate {
    daf::Size cliqueId;
    double ncrValue;
    uint8_t positions[8];
};

struct QuotientCleanStateEntry {
    uint8_t subNumPivot;
    uint64_t multiplicity;
    double contribution;
};

struct QuotientDeltaStateEntry {
    uint8_t removedPivotCount;
    uint8_t survivorPivotCount;
    uint8_t overlapPivotCount;
    uint64_t multiplicity;
    double delta;
};

struct QuotientLeafStateIndex {
    uint64_t cleanBegin{};
    uint32_t cleanCount{};
    uint64_t deltaBegin{};
    uint32_t deltaCount{};
};

struct QuotientLeafStat {
    daf::Size leafId{};
    int leafSize{};
    int keepC{};
    int pivotC{};
    uint64_t explicitEntries{};
    uint64_t quotientClasses{};
    uint64_t worstRefinedClasses{};
    uint64_t worstDeltaClasses{};
    uint64_t worstDeltaAffected{};
    uint64_t maxCleanMultiplicity{};
    uint64_t maxDeltaMultiplicity{};
    long double ratio{};
    long double refinedRatio{};
    long double deltaRatio{};
    long double deltaTouchFraction{};
};

struct QuotientSurveySummary {
    bool valid{};
    daf::CliqueSize r{};
    daf::CliqueSize s{};
    uint64_t totalLeaves{};
    uint64_t activeLeaves{};
    int maxLeafSize{};
    uint64_t totalExplicitEntries{};
    uint64_t totalQuotientClasses{};
};

struct DeltaStateSummary {
    uint64_t classCount{};
    uint64_t affectedCount{};
    uint64_t maxMultiplicity{};
};

static QuotientSurveySummary gQuotientSurveySummary{};

struct SparseCliqueKey {
    uint8_t len{};
    std::array<daf::Size, 8> verts{};

    bool operator==(const SparseCliqueKey &other) const noexcept {
        if (len != other.len) return false;
        for (uint8_t i = 0; i < len; ++i) {
            if (verts[i] != other.verts[i]) return false;
        }
        return true;
    }
};

struct SparseCliqueKeyHash {
    size_t operator()(const SparseCliqueKey &key) const noexcept {
        size_t h = static_cast<size_t>(key.len);
        for (uint8_t i = 0; i < key.len; ++i) {
            h ^= std::hash<daf::Size>()(key.verts[i]) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

static bool hasPositiveContribution(int pivotC, int keepC, int subNumPivot,
                                    daf::CliqueSize r, daf::CliqueSize s);

static uint64_t packLeafPosKey(const uint8_t *positions, daf::CliqueSize r) {
    uint64_t key = static_cast<uint64_t>(r);
    for (daf::CliqueSize i = 0; i < r; ++i) {
        key = (key << 8) | static_cast<uint64_t>(positions[i]);
    }
    return key;
}

static uint64_t packLeafPosKey(const std::vector<uint8_t> &positions) {
    uint64_t key = static_cast<uint64_t>(positions.size());
    for (auto pos : positions) {
        key = (key << 8) | static_cast<uint64_t>(pos);
    }
    return key;
}

static SparseCliqueKey makeSparseCliqueKey(const std::vector<TreeGraphNode> &leaf,
                                           const std::vector<uint8_t> &positions,
                                           daf::CliqueSize r) {
    SparseCliqueKey key;
    key.len = static_cast<uint8_t>(r);
    for (daf::CliqueSize j = 0; j < r; ++j) {
        key.verts[j] = leaf[positions[j]].v;
    }
    std::sort(key.verts.begin(), key.verts.begin() + r);
    return key;
}

template <class Map>
static void accumulateLeafSparseContribution(const std::vector<TreeGraphNode> &leaf,
                                             daf::CliqueSize r, daf::CliqueSize s,
                                             Map &target, double scale,
                                             uint64_t *entryCount = nullptr) {
    if (leaf.size() < r) return;
    int keepC = 0, pivotC = 0;
    for (const auto &node : leaf) {
        if (node.isPivot) pivotC++;
        else keepC++;
    }
    daf::enumerateCombinations(leaf, r,
        [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            int subNumPivot = 0;
            SparseCliqueKey key;
            key.len = static_cast<uint8_t>(r);
            for (daf::CliqueSize j = 0; j < r; ++j) {
                if (rClique[j].isPivot) subNumPivot++;
                key.verts[j] = rClique[j].v;
            }
            if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return true;
            const int row = pivotC - subNumPivot;
            const int col = static_cast<int>(s) - keepC - subNumPivot;
            if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return true;
            const double contrib = nCr[row][col];
            if (contrib <= 0.0) return true;
            std::sort(key.verts.begin(), key.verts.begin() + r);
            target[key] += scale * contrib;
            if (entryCount) (*entryCount)++;
            return true;
        });
}

template <class Map>
static void accumulateLeafSparseIncidence(const std::vector<TreeGraphNode> &leaf,
                                          daf::CliqueSize r, daf::CliqueSize s,
                                          Map &target, uint64_t *entryCount = nullptr) {
    if (leaf.size() < r) return;
    int keepC = 0, pivotC = 0;
    for (const auto &node : leaf) {
        if (node.isPivot) pivotC++;
        else keepC++;
    }
    daf::enumerateCombinations(leaf, r,
        [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            int subNumPivot = 0;
            SparseCliqueKey key;
            key.len = static_cast<uint8_t>(r);
            for (daf::CliqueSize j = 0; j < r; ++j) {
                if (rClique[j].isPivot) subNumPivot++;
                key.verts[j] = rClique[j].v;
            }
            if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return true;
            const int row = pivotC - subNumPivot;
            const int col = static_cast<int>(s) - keepC - subNumPivot;
            if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return true;
            const double contrib = nCr[row][col];
            if (contrib <= 0.0) return true;
            std::sort(key.verts.begin(), key.verts.begin() + r);
            target[key] += 1;
            if (entryCount) (*entryCount)++;
            return true;
        });
}

template <class Set>
static void collectLeafSparseKeys(const std::vector<TreeGraphNode> &leaf,
                                  daf::CliqueSize r, daf::CliqueSize s,
                                  Set &target, uint64_t *entryCount = nullptr) {
    if (leaf.size() < r) return;
    int keepC = 0, pivotC = 0;
    for (const auto &node : leaf) {
        if (node.isPivot) pivotC++;
        else keepC++;
    }
    daf::enumerateCombinations(leaf, r,
        [&](const daf::StaticVector<TreeGraphNode> &rClique) {
            int subNumPivot = 0;
            SparseCliqueKey key;
            key.len = static_cast<uint8_t>(r);
            for (daf::CliqueSize j = 0; j < r; ++j) {
                if (rClique[j].isPivot) subNumPivot++;
                key.verts[j] = rClique[j].v;
            }
            if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return true;
            const int row = pivotC - subNumPivot;
            const int col = static_cast<int>(s) - keepC - subNumPivot;
            if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return true;
            const double contrib = nCr[row][col];
            if (contrib <= 0.0) return true;
            std::sort(key.verts.begin(), key.verts.begin() + r);
            target.insert(key);
            if (entryCount) (*entryCount)++;
            return true;
        });
}

static double leafContributionForKey(const std::vector<TreeGraphNode> &leaf,
                                     const SparseCliqueKey &key,
                                     daf::CliqueSize r, daf::CliqueSize s) {
    if (leaf.size() < r || key.len != static_cast<uint8_t>(r)) return 0.0;
    int keepC = 0, pivotC = 0;
    for (daf::Size i = 0; i < leaf.size(); ++i) {
        daf::vListMap[leaf[i].v] = i;
        if (leaf[i].isPivot) pivotC++;
        else keepC++;
    }
    int subNumPivot = 0;
    for (daf::CliqueSize j = 0; j < r; ++j) {
        const daf::Size pos = daf::vListMap[key.verts[j]];
        if (pos >= leaf.size() || leaf[pos].v != key.verts[j]) return 0.0;
        if (leaf[pos].isPivot) subNumPivot++;
    }
    if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return 0.0;
    const int row = pivotC - subNumPivot;
    const int col = static_cast<int>(s) - keepC - subNumPivot;
    if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return 0.0;
    return nCr[row][col];
}

struct LeafKeyRole {
    bool valid = false;
    int keepC = 0;
    int pivotC = 0;
    int subNumPivot = 0;
    double contribution = 0.0;
};

static LeafKeyRole analyzeLeafKeyRole(const std::vector<TreeGraphNode> &leaf,
                                      const SparseCliqueKey &key,
                                      daf::CliqueSize r, daf::CliqueSize s) {
    LeafKeyRole role;
    if (leaf.size() < r || key.len != static_cast<uint8_t>(r)) return role;
    for (daf::Size i = 0; i < leaf.size(); ++i) {
        daf::vListMap[leaf[i].v] = i;
        if (leaf[i].isPivot) role.pivotC++;
        else role.keepC++;
    }
    for (daf::CliqueSize j = 0; j < r; ++j) {
        const daf::Size pos = daf::vListMap[key.verts[j]];
        if (pos >= leaf.size() || leaf[pos].v != key.verts[j]) return role;
        if (leaf[pos].isPivot) role.subNumPivot++;
    }
    if (!hasPositiveContribution(role.pivotC, role.keepC, role.subNumPivot, r, s)) return role;
    const int row = role.pivotC - role.subNumPivot;
    const int col = static_cast<int>(s) - role.keepC - role.subNumPivot;
    if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return role;
    role.contribution = nCr[row][col];
    role.valid = role.contribution > 0.0;
    return role;
}

static uint32_t activeLeafIncidenceForKey(
    const SparseCliqueKey &key, daf::CliqueSize r,
    const std::vector<std::vector<TreeGraphNode>> &activeLeaves,
    const std::vector<uint8_t> &leafAlive,
    std::vector<robin_hood::unordered_flat_set<daf::Size>> &activeLeafByVertex) {
    daf::StaticVector<daf::Size> verts;
    verts.resize(r);
    for (daf::CliqueSize j = 0; j < r; ++j) verts[j] = key.verts[j];
    uint32_t count = 0;
    daf::intersect_dense_sets_multi(verts, activeLeafByVertex,
        [&](const daf::Size &leafId) {
            if (leafId < activeLeaves.size() && leafAlive[leafId]) {
                count++;
            }
        });
    verts.free();
    return count;
}

static void accumulateExactLeafSparseDeltaBK(
    const std::vector<TreeGraphNode> &leaf,
    const std::vector<SparseCliqueKey> &removedKeys,
    daf::CliqueSize r, daf::CliqueSize s,
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> &deltaMap,
    uint64_t &oldEntries,
    uint64_t &newEntries,
    uint64_t &subleafCount) {
    accumulateLeafSparseContribution(leaf, r, s, deltaMap, -1.0, &oldEntries);

    daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
    for (daf::Size i = 0; i < leaf.size(); ++i) {
        mapRef[leaf[i].v] = i;
    }

    std::vector<std::vector<daf::Size>> conflictSets;
    conflictSets.reserve(removedKeys.size());
    for (const auto &rmKey : removedKeys) {
        std::vector<daf::Size> one;
        one.reserve(r);
        bool valid = true;
        for (daf::CliqueSize j = 0; j < r; ++j) {
            daf::Size pos = mapRef[rmKey.verts[j]];
            if (pos >= leaf.size() || leaf[pos].v != rmKey.verts[j]) {
                valid = false;
                break;
            }
            one.push_back(rmKey.verts[j]);
        }
        if (valid) conflictSets.push_back(std::move(one));
    }
    if (conflictSets.empty()) return;

    std::vector<TreeGraphNode> leafCopy(leaf.begin(), leaf.end());
    bkRmClique::removeRClique(leafCopy, conflictSets, r, s,
        [&](const bkRmClique::Bitset &cover, const bkRmClique::Bitset &pivots) {
            auto newLeaf = bkRmClique::coverToVertex(cover, pivots, leaf);
            subleafCount++;
            accumulateLeafSparseContribution(newLeaf, r, s, deltaMap, 1.0, &newEntries);
        });
}

static std::vector<std::vector<TreeGraphNode>> splitLeafByRemovedKeys(
    const std::vector<TreeGraphNode> &leaf,
    const std::vector<SparseCliqueKey> &removedKeys,
    daf::CliqueSize r, daf::CliqueSize s,
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> &deltaMap,
    uint64_t *oldEntries = nullptr,
    uint64_t *newEntries = nullptr) {
    if (oldEntries) {
        accumulateLeafSparseContribution(leaf, r, s, deltaMap, -1.0, oldEntries);
    } else {
        accumulateLeafSparseContribution(leaf, r, s, deltaMap, -1.0);
    }

    daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
    for (daf::Size i = 0; i < leaf.size(); ++i) {
        mapRef[leaf[i].v] = i;
    }

    std::vector<std::vector<daf::Size>> conflictSets;
    conflictSets.reserve(removedKeys.size());
    for (const auto &rmKey : removedKeys) {
        std::vector<daf::Size> one;
        one.reserve(r);
        bool valid = true;
        for (daf::CliqueSize j = 0; j < r; ++j) {
            daf::Size pos = mapRef[rmKey.verts[j]];
            if (pos >= leaf.size() || leaf[pos].v != rmKey.verts[j]) {
                valid = false;
                break;
            }
            one.push_back(rmKey.verts[j]);
        }
        if (valid) conflictSets.push_back(std::move(one));
    }
    if (conflictSets.empty()) return {};

    std::vector<std::vector<TreeGraphNode>> newLeaves;
    std::vector<TreeGraphNode> leafCopy(leaf.begin(), leaf.end());
    bkRmClique::removeRClique(leafCopy, conflictSets, r, s,
        [&](const bkRmClique::Bitset &cover, const bkRmClique::Bitset &pivots) {
            auto newLeaf = bkRmClique::coverToVertex(cover, pivots, leaf);
            if (newEntries) {
                accumulateLeafSparseContribution(newLeaf, r, s, deltaMap, 1.0, newEntries);
            } else {
                accumulateLeafSparseContribution(newLeaf, r, s, deltaMap, 1.0);
            }
            newLeaves.push_back(std::move(newLeaf));
        });
    return newLeaves;
}

static void addLeafToVertexIndex(
    const std::vector<TreeGraphNode> &leaf, daf::Size leafId,
    std::vector<robin_hood::unordered_flat_set<daf::Size>> &activeLeafByVertex,
    daf::CliqueSize r) {
    if (leaf.size() < r) return;
    for (const auto &node : leaf) {
        activeLeafByVertex[node.v].insert(leafId);
    }
}

static void removeLeafFromVertexIndex(
    const std::vector<TreeGraphNode> &leaf, daf::Size leafId,
    std::vector<robin_hood::unordered_flat_set<daf::Size>> &activeLeafByVertex,
    daf::CliqueSize r) {
    if (leaf.size() < r) return;
    for (const auto &node : leaf) {
        activeLeafByVertex[node.v].erase(leafId);
    }
}

static void addLeafToVertexAndPivotIndex(
    const std::vector<TreeGraphNode> &leaf, daf::Size leafId,
    std::vector<robin_hood::unordered_flat_set<daf::Size>> &activeLeafByVertex,
    std::vector<robin_hood::unordered_flat_set<daf::Size>> &activePivotLeafByVertex,
    std::vector<uint16_t> &activeKeepCount,
    std::vector<uint16_t> &activePivotCount,
    daf::CliqueSize r) {
    if (leafId >= activeKeepCount.size()) activeKeepCount.resize(leafId + 1, 0);
    if (leafId >= activePivotCount.size()) activePivotCount.resize(leafId + 1, 0);
    activeKeepCount[leafId] = 0;
    activePivotCount[leafId] = 0;
    if (leaf.size() < r) return;
    for (const auto &node : leaf) {
        activeLeafByVertex[node.v].insert(leafId);
        if (node.isPivot) {
            activePivotLeafByVertex[node.v].insert(leafId);
            activePivotCount[leafId]++;
        } else {
            activeKeepCount[leafId]++;
        }
    }
}

static void removeLeafFromVertexAndPivotIndex(
    const std::vector<TreeGraphNode> &leaf, daf::Size leafId,
    std::vector<robin_hood::unordered_flat_set<daf::Size>> &activeLeafByVertex,
    std::vector<robin_hood::unordered_flat_set<daf::Size>> &activePivotLeafByVertex,
    daf::CliqueSize r) {
    if (leaf.size() < r) return;
    for (const auto &node : leaf) {
        activeLeafByVertex[node.v].erase(leafId);
        if (node.isPivot) activePivotLeafByVertex[node.v].erase(leafId);
    }
}

static double leafContributionForKeyIndexed(
    daf::Size leafId, const SparseCliqueKey &key,
    const std::vector<uint16_t> &activeKeepCount,
    const std::vector<uint16_t> &activePivotCount,
    std::vector<robin_hood::unordered_flat_set<daf::Size>> &activePivotLeafByVertex,
    daf::CliqueSize r, daf::CliqueSize s) {
    if (leafId >= activeKeepCount.size() || key.len != static_cast<uint8_t>(r)) return 0.0;
    const int keepC = static_cast<int>(activeKeepCount[leafId]);
    const int pivotC = static_cast<int>(activePivotCount[leafId]);
    int subNumPivot = 0;
    for (daf::CliqueSize j = 0; j < r; ++j) {
        if (activePivotLeafByVertex[key.verts[j]].find(leafId) != activePivotLeafByVertex[key.verts[j]].end()) {
            subNumPivot++;
        }
    }
    if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return 0.0;
    const int row = pivotC - subNumPivot;
    const int col = static_cast<int>(s) - keepC - subNumPivot;
    if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return 0.0;
    return nCr[row][col];
}

static uint64_t combSmall(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    uint64_t res = 1;
    for (int i = 1; i <= k; ++i) {
        res = (res * static_cast<uint64_t>(n - k + i)) / static_cast<uint64_t>(i);
    }
    return res;
}

template <class Fn>
static void forEachChoose(const std::vector<uint8_t> &src, int need,
                          std::vector<uint8_t> &picked, Fn &&fn, int start = 0) {
    if (need == 0) {
        fn();
        return;
    }
    for (int i = start; i <= static_cast<int>(src.size()) - need; ++i) {
        picked.push_back(src[i]);
        forEachChoose(src, need - 1, picked, fn, i + 1);
        picked.pop_back();
    }
}

template <class Fn>
static void forEachExpandedCleanClass(const std::vector<TreeGraphNode> &leaf,
                                      daf::CliqueSize r, int qP, Fn &&fn) {
    std::vector<uint8_t> keepPos;
    std::vector<uint8_t> pivotPos;
    for (uint8_t i = 0; i < leaf.size(); ++i) {
        if (leaf[i].isPivot) pivotPos.push_back(i);
        else keepPos.push_back(i);
    }
    const int needKeep = static_cast<int>(r) - qP;
    if (qP < 0 || needKeep < 0) return;
    if (qP > static_cast<int>(pivotPos.size()) || needKeep > static_cast<int>(keepPos.size())) return;

    std::vector<uint8_t> pickedPivot, pickedKeep, merged;
    merged.reserve(r);
    forEachChoose(pivotPos, qP, pickedPivot, [&]() {
        forEachChoose(keepPos, needKeep, pickedKeep, [&]() {
            merged.clear();
            std::merge(pickedPivot.begin(), pickedPivot.end(),
                       pickedKeep.begin(), pickedKeep.end(),
                       std::back_inserter(merged));
            fn(merged);
        });
    });
}

template <class Fn>
static void forEachExpandedCanonicalDeltaClass(const std::vector<TreeGraphNode> &leaf,
                                               daf::CliqueSize r, int remP,
                                               int qP, int overlapP, Fn &&fn) {
    std::vector<uint8_t> keepPos;
    std::vector<uint8_t> pivotPos;
    for (uint8_t i = 0; i < leaf.size(); ++i) {
        if (leaf[i].isPivot) pivotPos.push_back(i);
        else keepPos.push_back(i);
    }

    const int remK = static_cast<int>(r) - remP;
    const int needKeep = static_cast<int>(r) - qP;
    if (remP < 0 || qP < 0 || overlapP < 0 || remK < 0 || needKeep < 0) return;
    if (remP > static_cast<int>(pivotPos.size()) || remK > static_cast<int>(keepPos.size())) return;
    if (qP > static_cast<int>(pivotPos.size()) || needKeep > static_cast<int>(keepPos.size())) return;
    if (overlapP > remP || overlapP > qP) return;

    std::vector<uint8_t> removedPivotSet(pivotPos.begin(), pivotPos.begin() + remP);
    std::vector<uint8_t> removedKeepSet(keepPos.begin(), keepPos.begin() + remK);
    std::vector<uint8_t> freePivotSet(pivotPos.begin() + remP, pivotPos.end());
    const int freeNeed = qP - overlapP;
    if (freeNeed < 0 || freeNeed > static_cast<int>(freePivotSet.size())) return;

    std::vector<uint8_t> pickedOverlap, pickedFree, pickedKeep;
    std::vector<uint8_t> pivMerged, merged;
    pivMerged.reserve(qP);
    merged.reserve(r);
    forEachChoose(removedPivotSet, overlapP, pickedOverlap, [&]() {
        forEachChoose(freePivotSet, freeNeed, pickedFree, [&]() {
            forEachChoose(keepPos, needKeep, pickedKeep, [&]() {
                if (pickedOverlap == removedPivotSet &&
                    pickedFree.empty() &&
                    pickedKeep == removedKeepSet) {
                    return;
                }
                pivMerged.clear();
                merged.clear();
                std::merge(pickedOverlap.begin(), pickedOverlap.end(),
                           pickedFree.begin(), pickedFree.end(),
                           std::back_inserter(pivMerged));
                std::merge(pivMerged.begin(), pivMerged.end(),
                           pickedKeep.begin(), pickedKeep.end(),
                           std::back_inserter(merged));
                fn(merged);
            });
        });
    });
}

static bool hasPositiveContribution(int pivotC, int keepC, int subNumPivot,
                                    daf::CliqueSize r, daf::CliqueSize s) {
    const int needPivot = static_cast<int>(s) - keepC;
    if (subNumPivot > needPivot) return false;
    const int row = pivotC - subNumPivot;
    const int col = needPivot - subNumPivot;
    if (row < 0 || row >= 1001 || col < 0 || col >= 401) return false;
    if (col > row) return false;
    return nCr[row][col] > 0.0;
}

static uint64_t countBoundedClassStates(const std::vector<int> &caps, int target) {
    uint64_t cnt = 0;
    const int g = static_cast<int>(caps.size());
    if (g == 0) return 0;
    std::vector<int> pick(g, 0);
    std::function<void(int, int)> dfs = [&](int idx, int left) {
        if (idx + 1 == g) {
            if (left >= 0 && left <= caps[idx]) cnt++;
            return;
        }
        const int hi = std::min(caps[idx], left);
        for (int x = 0; x <= hi; ++x) dfs(idx + 1, left - x);
    };
    dfs(0, target);
    return cnt;
}

static uint64_t worstOneStepRefinedClasses(int keepC, int pivotC,
                                           daf::CliqueSize r, daf::CliqueSize s) {
    const int low = std::max(0, static_cast<int>(r) - keepC);
    const int high = std::min(static_cast<int>(r), pivotC);
    uint64_t worst = 0;
    for (int remP = low; remP <= high; ++remP) {
        if (!hasPositiveContribution(pivotC, keepC, remP, r, s)) continue;
        const int remK = static_cast<int>(r) - remP;
        if (remK < 0 || remK > keepC || remP > pivotC) continue;
        std::vector<int> caps = {
            remK,
            keepC - remK,
            remP,
            pivotC - remP
        };
        uint64_t cls = countBoundedClassStates(caps, static_cast<int>(r));
        if (cls > 0) cls -= 1; // the removed r-clique itself disappears
        worst = std::max(worst, cls);
    }
    return worst;
}

static DeltaStateSummary worstOneRemovedDeltaState(
    int keepC, int pivotC, daf::CliqueSize r, daf::CliqueSize s) {
    const int needPivot = static_cast<int>(s) - keepC;
    const int low = std::max(0, static_cast<int>(r) - keepC);
    const int high = std::min(static_cast<int>(r), pivotC);

    uint64_t bestClasses = 0;
    uint64_t bestAffected = 0;
    uint64_t bestMaxMultiplicity = 0;
    long double bestRatio = std::numeric_limits<long double>::infinity();

    for (int remP = low; remP <= high; ++remP) {
        if (!hasPositiveContribution(pivotC, keepC, remP, r, s)) continue;

        uint64_t classCnt = 0;
        uint64_t affectedCnt = 0;
        uint64_t maxMultiplicity = 0;

        for (int qP = low; qP <= high; ++qP) {
            if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;

            const int xLow = std::max(0, qP - (pivotC - remP));
            const int xHigh = std::min(qP, remP);
            for (int x = xLow; x <= xHigh; ++x) {
                uint64_t mult = combSmall(remP, x) *
                                combSmall(pivotC - remP, qP - x) *
                                combSmall(keepC, static_cast<int>(r) - qP);
                if (qP == remP && x == remP && mult > 0) {
                    mult -= 1;
                }
                if (mult == 0) continue;

                const int unionP = remP + qP - x;
                if (unionP > needPivot) continue;
                const int row = pivotC - unionP;
                const int col = needPivot - unionP;
                if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) continue;
                if (nCr[row][col] <= 0.0) continue;

                classCnt++;
                affectedCnt += mult;
                maxMultiplicity = std::max(maxMultiplicity, mult);
            }
        }

        if (classCnt == 0 || affectedCnt == 0) continue;
        long double ratio = static_cast<long double>(affectedCnt) /
                            static_cast<long double>(classCnt);
        if (ratio < bestRatio) {
            bestRatio = ratio;
            bestClasses = classCnt;
            bestAffected = affectedCnt;
            bestMaxMultiplicity = maxMultiplicity;
        }
    }

    if (bestClasses == 0) return {};
    return {bestClasses, bestAffected, bestMaxMultiplicity};
}

static bool chooseContainsAll(const std::vector<uint8_t> &chosen, const std::vector<uint8_t> &need) {
    size_t i = 0, j = 0;
    while (i < chosen.size() && j < need.size()) {
        if (chosen[i] == need[j]) { ++i; ++j; }
        else if (chosen[i] < need[j]) ++i;
        else return false;
    }
    return j == need.size();
}

static uint64_t countExpandedCleanClass(const std::vector<TreeGraphNode> &leaf,
                                        daf::CliqueSize r, int qP) {
    uint64_t cnt = 0;
    forEachExpandedCleanClass(leaf, r, qP, [&](const std::vector<uint8_t> &) { ++cnt; });
    return cnt;
}

static uint64_t countExpandedCanonicalDeltaClass(const std::vector<TreeGraphNode> &leaf,
                                                 daf::CliqueSize r, int remP,
                                                 int qP, int overlapP) {
    uint64_t cnt = 0;
    forEachExpandedCanonicalDeltaClass(leaf, r, remP, qP, overlapP,
                                       [&](const std::vector<uint8_t> &) { ++cnt; });
    return cnt;
}

static void maybeVerifyOneRemovedFormula(const DynamicGraph<TreeGraphNode> &tree,
                                         daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_VERIFY")) return;

    int checkedLeaves = 0;
    for (daf::Size leafId = 0; leafId < tree.adj_list.size() && checkedLeaves < 3; ++leafId) {
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.size() < r || leaf.size() > 18) continue;

        int keepC = 0, pivotC = 0;
        std::vector<uint8_t> pivotPos;
        for (uint8_t i = 0; i < leaf.size(); ++i) {
            if (leaf[i].isPivot) {
                pivotPos.push_back(i);
                pivotC++;
            } else {
                keepC++;
            }
        }
        const int needPivot = static_cast<int>(s) - keepC;
        if (needPivot < 0 || needPivot > pivotC) continue;

        struct RCliqueInfo {
            std::vector<uint8_t> pos;
            std::vector<uint8_t> pivs;
            int qP{};
        };
        std::vector<RCliqueInfo> rCliques;
        daf::enumerateCombinationsWithIdx(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique, const size_t* idx) {
                int qP = 0;
                std::vector<uint8_t> pos;
                std::vector<uint8_t> pivs;
                for (daf::Size j = 0; j < r; ++j) {
                    uint8_t p = static_cast<uint8_t>(idx[j]);
                    pos.push_back(p);
                    if (rClique[j].isPivot) {
                        pivs.push_back(p);
                        qP++;
                    }
                }
                if (hasPositiveContribution(pivotC, keepC, qP, r, s)) {
                    rCliques.push_back({std::move(pos), std::move(pivs), qP});
                }
                return true;
            });
        if (rCliques.empty()) continue;

        std::vector<std::vector<uint8_t>> sCliques;
        std::vector<uint8_t> chosen;
        std::function<void(int, int)> dfs = [&](int idx, int left) {
            if (left == 0) {
                sCliques.push_back(chosen);
                return;
            }
            for (int i = idx; i <= pivotC - left; ++i) {
                chosen.push_back(pivotPos[i]);
                dfs(i + 1, left - 1);
                chosen.pop_back();
            }
        };
        dfs(0, needPivot);
        if (sCliques.empty()) continue;

        bool cleanOk = true;
        std::map<int, uint64_t> actualCleanMultiplicity;
        for (const auto &Q : rCliques) actualCleanMultiplicity[Q.qP]++;
        const int low = std::max(0, static_cast<int>(r) - keepC);
        const int high = std::min(static_cast<int>(r), pivotC);
        for (int qP = low; qP <= high; ++qP) {
            if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
            uint64_t expectMult =
                combSmall(pivotC, qP) * combSmall(keepC, static_cast<int>(r) - qP);
            if (actualCleanMultiplicity[qP] != expectMult) {
                cleanOk = false;
                break;
            }
            const int row = pivotC - qP;
            const int col = needPivot - qP;
            if (!(row >= 0 && row < 1001 && col >= 0 && col < 401 && col <= row &&
                  nCr[row][col] > 0.0)) {
                cleanOk = false;
                break;
            }
        }

        std::map<int, size_t> remSeen;
        bool deltaOk = true;
        bool expandOk = true;
        for (const auto &R : rCliques) {
            if (remSeen[R.qP] > 0) continue;
            remSeen[R.qP]++;

            std::map<std::pair<int, int>, std::pair<uint64_t, int>> groupStat;
            for (const auto &Q : rCliques) {
                bool same = (Q.pos == R.pos);
                if (same) continue;

                int overlapP = 0;
                size_t i = 0, j = 0;
                while (i < Q.pivs.size() && j < R.pivs.size()) {
                    if (Q.pivs[i] == R.pivs[j]) { overlapP++; ++i; ++j; }
                    else if (Q.pivs[i] < R.pivs[j]) ++i;
                    else ++j;
                }

                std::vector<uint8_t> unionPivs;
                std::set_union(Q.pivs.begin(), Q.pivs.end(),
                               R.pivs.begin(), R.pivs.end(),
                               std::back_inserter(unionPivs));
                int actualDelta = 0;
                for (const auto &S : sCliques) {
                    if (chooseContainsAll(S, unionPivs)) actualDelta++;
                }

                auto key = std::make_pair(Q.qP, overlapP);
                auto &slot = groupStat[key];
                slot.first++;
                if (slot.second == 0) slot.second = actualDelta + 1;
                else if (slot.second != actualDelta + 1) deltaOk = false;
            }

            for (const auto &[key, val] : groupStat) {
                const int qP = key.first;
                const int overlapP = key.second;
                const uint64_t actualMult = val.first;
                const int actualDelta = val.second - 1;

                uint64_t expectMult = combSmall(R.qP, overlapP) *
                                      combSmall(pivotC - R.qP, qP - overlapP) *
                                      combSmall(keepC, static_cast<int>(r) - qP);
                if (qP == R.qP && overlapP == R.qP && expectMult > 0) expectMult -= 1;
                uint64_t expandMult = countExpandedCanonicalDeltaClass(
                    leaf, r, R.qP, qP, overlapP);

                const int unionP = R.qP + qP - overlapP;
                int expectDelta = 0;
                if (unionP <= needPivot) {
                    expectDelta = static_cast<int>(std::llround(
                        nCr[pivotC - unionP][needPivot - unionP]));
                }

                if (actualMult != expectMult || actualDelta != expectDelta) {
                    deltaOk = false;
                    break;
                }
                if (expandMult != expectMult) {
                    expandOk = false;
                    break;
                }
            }
            if (!deltaOk || !expandOk) break;
        }

        for (int qP = low; qP <= high && expandOk; ++qP) {
            if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
            uint64_t expectMult =
                combSmall(pivotC, qP) * combSmall(keepC, static_cast<int>(r) - qP);
            uint64_t expandMult = countExpandedCleanClass(leaf, r, qP);
            if (expandMult != expectMult) {
                expandOk = false;
            }
        }

        bool updateOk = true;
        for (int remP = low; remP <= high && updateOk; ++remP) {
            if (!hasPositiveContribution(pivotC, keepC, remP, r, s)) continue;
            const int remK = static_cast<int>(r) - remP;
            if (remK < 0 || remK > keepC) {
                updateOk = false;
                break;
            }

            std::vector<uint8_t> removedPivotSet(pivotPos.begin(), pivotPos.begin() + remP);
            std::vector<uint8_t> keepPos;
            for (uint8_t i = 0; i < leaf.size(); ++i) {
                if (!leaf[i].isPivot) keepPos.push_back(i);
            }
            std::vector<uint8_t> removedKeepSet(keepPos.begin(), keepPos.begin() + remK);
            std::vector<uint8_t> removedPos;
            std::merge(removedPivotSet.begin(), removedPivotSet.end(),
                       removedKeepSet.begin(), removedKeepSet.end(),
                       std::back_inserter(removedPos));

            std::map<std::vector<uint8_t>, int> actualAfter;
            std::map<std::vector<uint8_t>, int> expectAfter;
            for (const auto &Q : rCliques) {
                if (Q.pos == removedPos) continue;

                int actualNew = 0;
                for (const auto &S : sCliques) {
                    if (!chooseContainsAll(S, Q.pivs)) continue;
                    if (chooseContainsAll(S, removedPivotSet)) continue;
                    actualNew++;
                }
                actualAfter[Q.pos] = actualNew;

                int row = pivotC - Q.qP;
                int col = needPivot - Q.qP;
                expectAfter[Q.pos] = static_cast<int>(std::llround(nCr[row][col]));
            }

            for (int qP = low; qP <= high && updateOk; ++qP) {
                if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
                int base = static_cast<int>(std::llround(nCr[pivotC - qP][needPivot - qP]));
                const int xLow = std::max(0, qP - (pivotC - remP));
                const int xHigh = std::min(qP, remP);
                for (int x = xLow; x <= xHigh && updateOk; ++x) {
                    const int unionP = remP + qP - x;
                    if (unionP > needPivot) continue;
                    int delta = static_cast<int>(std::llround(nCr[pivotC - unionP][needPivot - unionP]));
                    forEachExpandedCanonicalDeltaClass(leaf, r, remP, qP, x,
                        [&](const std::vector<uint8_t> &pos) {
                            auto it = expectAfter.find(pos);
                            if (it == expectAfter.end()) {
                                updateOk = false;
                                return;
                            }
                            it->second = base - delta;
                        });
                }
            }

            if (!updateOk || actualAfter.size() != expectAfter.size()) {
                updateOk = false;
                break;
            }
            for (const auto &[pos, expected] : expectAfter) {
                auto it = actualAfter.find(pos);
                if (it == actualAfter.end() || it->second != expected) {
                    updateOk = false;
                    break;
                }
            }
        }

        std::cout << "Quotient one-removed formula verify: leaf=" << leafId
                  << " size=" << leaf.size()
                  << " clean=" << (cleanOk ? "OK" : "FAIL")
                  << " delta=" << (deltaOk ? "OK" : "FAIL")
                  << " expand=" << (expandOk ? "OK" : "FAIL")
                  << " update=" << (updateOk ? "OK" : "FAIL") << std::endl;
        checkedLeaves++;
    }
}

static void printQuotientStats(const DynamicGraph<TreeGraphNode> &tree,
                               daf::CliqueSize r, daf::CliqueSize s) {
    maybeVerifyOneRemovedFormula(tree, r, s);
    gQuotientSurveySummary = {};
    uint64_t totalLeaves = 0;
    uint64_t activeLeaves = 0;
    uint64_t totalExplicitEntries = 0;
    uint64_t totalQuotientClasses = 0;
    uint64_t totalWorstRefinedClasses = 0;
    uint64_t totalWorstDeltaClasses = 0;
    uint64_t totalWorstDeltaAffected = 0;
    uint64_t maxExplicitEntries = 0;
    uint64_t maxCleanMultiplicity = 0;
    uint64_t maxDeltaMultiplicity = 0;
    uint64_t keepGtRLeaves = 0;
    uint64_t keepGtRExplicitEntries = 0;
    uint64_t keepGtRQuotientClasses = 0;
    int maxLeafSize = 0;
    std::vector<long double> ratios;
    std::vector<long double> refinedRatios;
    std::vector<long double> deltaRatios;
    std::vector<long double> cleanExpansionPeaks;
    std::vector<long double> deltaExpansionPeaks;
    std::vector<long double> deltaTouchFractions;
    std::vector<QuotientLeafStat> topLeaves;
    topLeaves.reserve(8);

    for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.size() < r) continue;
        totalLeaves++;

        int keepC = 0, pivotC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }

        uint64_t explicitEntries = 0;
        uint64_t quotientClasses = 0;
        uint64_t leafMaxCleanMultiplicity = 0;
        const int low = std::max(0, static_cast<int>(r) - keepC);
        const int high = std::min(static_cast<int>(r), pivotC);
        for (int subNumPivot = low; subNumPivot <= high; ++subNumPivot) {
            if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) continue;
            uint64_t mult = combSmall(pivotC, subNumPivot) *
                            combSmall(keepC, static_cast<int>(r) - subNumPivot);
            explicitEntries += mult;
            leafMaxCleanMultiplicity = std::max(leafMaxCleanMultiplicity, mult);
            quotientClasses++;
        }

        if (explicitEntries == 0 || quotientClasses == 0) continue;
        if (keepC > static_cast<int>(r)) {
            keepGtRLeaves++;
            keepGtRExplicitEntries += explicitEntries;
            keepGtRQuotientClasses += quotientClasses;
        }
        uint64_t worstRefinedClasses = worstOneStepRefinedClasses(keepC, pivotC, r, s);
        if (worstRefinedClasses == 0) worstRefinedClasses = quotientClasses;
        auto deltaSummary = worstOneRemovedDeltaState(keepC, pivotC, r, s);
        uint64_t worstDeltaClasses = deltaSummary.classCount;
        uint64_t worstDeltaAffected = deltaSummary.affectedCount;
        uint64_t leafMaxDeltaMultiplicity = deltaSummary.maxMultiplicity;
        if (worstDeltaClasses == 0) {
            worstDeltaClasses = quotientClasses;
            worstDeltaAffected = explicitEntries;
            leafMaxDeltaMultiplicity = leafMaxCleanMultiplicity;
        }

        activeLeaves++;
        totalExplicitEntries += explicitEntries;
        totalQuotientClasses += quotientClasses;
        totalWorstRefinedClasses += worstRefinedClasses;
        totalWorstDeltaClasses += worstDeltaClasses;
        totalWorstDeltaAffected += worstDeltaAffected;
        maxExplicitEntries = std::max(maxExplicitEntries, explicitEntries);
        maxCleanMultiplicity = std::max(maxCleanMultiplicity, leafMaxCleanMultiplicity);
        maxDeltaMultiplicity = std::max(maxDeltaMultiplicity, leafMaxDeltaMultiplicity);
        maxLeafSize = std::max(maxLeafSize, static_cast<int>(leaf.size()));

        long double ratio = static_cast<long double>(explicitEntries) /
                            static_cast<long double>(quotientClasses);
        long double refinedRatio = static_cast<long double>(explicitEntries) /
                                   static_cast<long double>(worstRefinedClasses);
        long double deltaRatio = static_cast<long double>(worstDeltaAffected) /
                                 static_cast<long double>(worstDeltaClasses);
        long double deltaTouchFraction = static_cast<long double>(worstDeltaAffected) /
                                         static_cast<long double>(explicitEntries);
        ratios.push_back(ratio);
        refinedRatios.push_back(refinedRatio);
        deltaRatios.push_back(deltaRatio);
        cleanExpansionPeaks.push_back(static_cast<long double>(leafMaxCleanMultiplicity));
        deltaExpansionPeaks.push_back(static_cast<long double>(leafMaxDeltaMultiplicity));
        deltaTouchFractions.push_back(deltaTouchFraction);

        QuotientLeafStat stat{
            leafId,
            static_cast<int>(leaf.size()),
            keepC,
            pivotC,
            explicitEntries,
            quotientClasses,
            worstRefinedClasses,
            worstDeltaClasses,
            worstDeltaAffected,
            leafMaxCleanMultiplicity,
            leafMaxDeltaMultiplicity,
            ratio,
            refinedRatio,
            deltaRatio,
            deltaTouchFraction
        };
        topLeaves.push_back(stat);
    }

    std::sort(ratios.begin(), ratios.end());
    std::sort(refinedRatios.begin(), refinedRatios.end());
    std::sort(deltaRatios.begin(), deltaRatios.end());
    std::sort(cleanExpansionPeaks.begin(), cleanExpansionPeaks.end());
    std::sort(deltaExpansionPeaks.begin(), deltaExpansionPeaks.end());
    std::sort(deltaTouchFractions.begin(), deltaTouchFractions.end());
    std::sort(topLeaves.begin(), topLeaves.end(),
              [](const QuotientLeafStat &a, const QuotientLeafStat &b) {
                  if (a.ratio != b.ratio) return a.ratio > b.ratio;
                  return a.explicitEntries > b.explicitEntries;
              });

    std::cout << "=================== QuotientLab (clean leaf) ===================" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Leaves (size>=r): " << totalLeaves << std::endl;
    std::cout << "  Active leaves:    " << activeLeaves << std::endl;
    std::cout << "  Max leaf size:    " << maxLeafSize << std::endl;
    std::cout << "  Max explicit C(|L|,r)-like entries: " << maxExplicitEntries << std::endl;
    std::cout << "  Total explicit entries: " << totalExplicitEntries << std::endl;
    std::cout << "  Total quotient classes: " << totalQuotientClasses << std::endl;
    std::cout << "  Total one-step refined classes (worst case): "
              << totalWorstRefinedClasses << std::endl;
    std::cout << "  Total one-removed delta classes (worst case): "
              << totalWorstDeltaClasses << std::endl;
    std::cout << "  Total one-removed affected entries (worst case): "
              << totalWorstDeltaAffected << std::endl;
    std::cout << "  Leaves with keep>r: " << keepGtRLeaves << std::endl;
    std::cout << "  Explicit entries from keep>r leaves: "
              << keepGtRExplicitEntries << std::endl;
    std::cout << "  Quotient classes from keep>r leaves: "
              << keepGtRQuotientClasses << std::endl;
    std::cout << "  Largest clean class multiplicity: " << maxCleanMultiplicity << std::endl;
    std::cout << "  Largest one-removed delta class multiplicity: "
              << maxDeltaMultiplicity << std::endl;
    std::cout << "  Estimated state memory (explicit clean): "
              << (totalExplicitEntries * sizeof(ExplicitLeafEntryEstimate)) / 1024 / 1024
              << " MB" << std::endl;
    std::cout << "  Estimated state memory (quotient clean): "
              << (totalQuotientClasses * sizeof(QuotientCleanStateEntry)) / 1024 / 1024
              << " MB" << std::endl;
    std::cout << "  Estimated state memory (one-step refined): "
              << (totalWorstRefinedClasses * sizeof(QuotientCleanStateEntry)) / 1024 / 1024
              << " MB" << std::endl;
    std::cout << "  Estimated state memory (one-removed delta): "
              << (totalWorstDeltaClasses * sizeof(QuotientDeltaStateEntry)) / 1024 / 1024
              << " MB" << std::endl;

    if (totalQuotientClasses > 0) {
        long double totalRatio =
            static_cast<long double>(totalExplicitEntries) /
            static_cast<long double>(totalQuotientClasses);
        long double totalRefinedRatio =
            static_cast<long double>(totalExplicitEntries) /
            static_cast<long double>(totalWorstRefinedClasses);
        long double totalDeltaRatio =
            static_cast<long double>(totalWorstDeltaAffected) /
            static_cast<long double>(totalWorstDeltaClasses);
        long double totalDeltaTouchFraction =
            static_cast<long double>(totalWorstDeltaAffected) /
            static_cast<long double>(totalExplicitEntries);
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "  Total clean quotient compression: " << totalRatio << "x" << std::endl;
        std::cout << "  Total one-step refined compression: " << totalRefinedRatio << "x" << std::endl;
        std::cout << "  Total one-removed delta compression: " << totalDeltaRatio << "x" << std::endl;
        std::cout << std::setprecision(6)
                  << "  Total one-removed touch fraction: "
                  << (100.0L * totalDeltaTouchFraction) << "%" << std::endl
                  << std::setprecision(2);
        if (!ratios.empty()) {
            auto pick = [&](double q) {
                size_t idx = std::min(
                    ratios.size() - 1,
                    static_cast<size_t>(q * static_cast<double>(ratios.size() - 1)));
                return ratios[idx];
            };
            auto pickRefined = [&](double q) {
                size_t idx = std::min(
                    refinedRatios.size() - 1,
                    static_cast<size_t>(q * static_cast<double>(refinedRatios.size() - 1)));
                return refinedRatios[idx];
            };
            auto pickDelta = [&](double q) {
                size_t idx = std::min(
                    deltaRatios.size() - 1,
                    static_cast<size_t>(q * static_cast<double>(deltaRatios.size() - 1)));
                return deltaRatios[idx];
            };
            std::cout << "  Per-leaf compression median/p90/max: "
                      << pick(0.50) << "x / "
                      << pick(0.90) << "x / "
                      << ratios.back() << "x" << std::endl;
            std::cout << "  Per-leaf one-step refined median/p90/max: "
                      << pickRefined(0.50) << "x / "
                      << pickRefined(0.90) << "x / "
                      << refinedRatios.back() << "x" << std::endl;
            std::cout << "  Per-leaf one-removed delta median/p90/max: "
                      << pickDelta(0.50) << "x / "
                      << pickDelta(0.90) << "x / "
                      << deltaRatios.back() << "x" << std::endl;
            auto pickCleanPeak = [&](double q) {
                size_t idx = std::min(
                    cleanExpansionPeaks.size() - 1,
                    static_cast<size_t>(q * static_cast<double>(cleanExpansionPeaks.size() - 1)));
                return cleanExpansionPeaks[idx];
            };
            auto pickDeltaPeak = [&](double q) {
                size_t idx = std::min(
                    deltaExpansionPeaks.size() - 1,
                    static_cast<size_t>(q * static_cast<double>(deltaExpansionPeaks.size() - 1)));
                return deltaExpansionPeaks[idx];
            };
            auto pickTouch = [&](double q) {
                size_t idx = std::min(
                    deltaTouchFractions.size() - 1,
                    static_cast<size_t>(q * static_cast<double>(deltaTouchFractions.size() - 1)));
                return deltaTouchFractions[idx];
            };
            std::cout << "  Largest clean class multiplicity median/p90/max: "
                      << pickCleanPeak(0.50) << " / "
                      << pickCleanPeak(0.90) << " / "
                      << cleanExpansionPeaks.back() << std::endl;
            std::cout << "  Largest delta class multiplicity median/p90/max: "
                      << pickDeltaPeak(0.50) << " / "
                      << pickDeltaPeak(0.90) << " / "
                      << deltaExpansionPeaks.back() << std::endl;
            std::cout << std::setprecision(4)
                      << "  One-removed touch fraction median/p90/max: "
                      << (100.0L * pickTouch(0.50)) << "% / "
                      << (100.0L * pickTouch(0.90)) << "% / "
                      << (100.0L * deltaTouchFractions.back()) << "%" << std::endl
                      << std::setprecision(2);
        }
        std::cout.unsetf(std::ios::floatfield);
    }

    std::cout << "  Top quotient leaves:" << std::endl;
    for (size_t i = 0; i < std::min<size_t>(5, topLeaves.size()); ++i) {
        const auto &st = topLeaves[i];
        std::cout << "    leaf=" << st.leafId
                  << " size=" << st.leafSize
                  << " keep=" << st.keepC
                  << " pivot=" << st.pivotC
                  << " explicit=" << st.explicitEntries
                  << " classes=" << st.quotientClasses
                  << " refined=" << st.worstRefinedClasses
                  << " deltaCls=" << st.worstDeltaClasses
                  << " deltaAff=" << st.worstDeltaAffected
                  << " maxCleanMult=" << st.maxCleanMultiplicity
                  << " maxDeltaMult=" << st.maxDeltaMultiplicity
                  << " ratio=" << std::fixed << std::setprecision(2) << st.ratio << "x"
                  << " refinedRatio=" << st.refinedRatio << "x"
                  << " deltaRatio=" << st.deltaRatio << "x"
                  << std::setprecision(6)
                  << " touch=" << (100.0L * st.deltaTouchFraction) << "%"
                  << std::setprecision(2)
                  << std::endl;
    }
    std::cout.unsetf(std::ios::floatfield);
    std::cout << "===============================================================" << std::endl;

    gQuotientSurveySummary.valid = true;
    gQuotientSurveySummary.r = r;
    gQuotientSurveySummary.s = s;
    gQuotientSurveySummary.totalLeaves = totalLeaves;
    gQuotientSurveySummary.activeLeaves = activeLeaves;
    gQuotientSurveySummary.maxLeafSize = maxLeafSize;
    gQuotientSurveySummary.totalExplicitEntries = totalExplicitEntries;
    gQuotientSurveySummary.totalQuotientClasses = totalQuotientClasses;
}

static void maybeCompareIndexCoverage(const DynamicGraph<TreeGraphNode> &tree,
                                      const Graph &edgeGraph,
                                      daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_COMPARE_INDEX")) return;

    StaticCliqueIndex cliqueIndex(r);
    daf::timeCount("clique Index build (Quotient compare)", [&]() {
        cliqueIndex.build(tree, edgeGraph.adj_list.size());
    });

    robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> seenKeepGtR;
    uint64_t uniqueKeepGtR = 0;
    uint64_t coveredKeepGtR = 0;
    uint64_t missedKeepGtR = 0;

    for (const auto &leaf : tree.adj_list) {
        if (leaf.size() < r) continue;
        int keepC = 0, pivotC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }
        if (keepC <= static_cast<int>(r)) continue;

        daf::enumerateCombinations(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                int subNumPivot = 0;
                SparseCliqueKey key;
                key.len = static_cast<uint8_t>(r);
                for (daf::CliqueSize j = 0; j < r; ++j) {
                    if (rClique[j].isPivot) subNumPivot++;
                    key.verts[j] = rClique[j].v;
                }
                if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return true;
                std::sort(key.verts.begin(), key.verts.begin() + r);
                if (!seenKeepGtR.emplace(key).second) return true;
                uniqueKeepGtR++;
                auto cid = cliqueIndex.tryByClique(std::span<const daf::Size>(key.verts.data(), r));
                if (cid < cliqueIndex.size()) coveredKeepGtR++;
                else missedKeepGtR++;
                return true;
            });
    }

    std::cout << "================ Quotient Index Coverage =================" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Unique keep>r positive cliques: " << uniqueKeepGtR << std::endl;
    std::cout << "  Covered by StaticCliqueIndex:   " << coveredKeepGtR << std::endl;
    std::cout << "  Missed by StaticCliqueIndex:    " << missedKeepGtR << std::endl;
    std::cout << "==========================================================" << std::endl;
}

static void maybeBuildPrototypeState(const DynamicGraph<TreeGraphNode> &tree,
                                     daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_BUILD_STATE")) return;

    std::vector<QuotientLeafStateIndex> leafIndex(tree.adj_list.size());
    std::vector<QuotientCleanStateEntry> cleanEntries;
    std::vector<QuotientDeltaStateEntry> deltaEntries;
    cleanEntries.reserve(tree.adj_list.size() * 2);
    deltaEntries.reserve(tree.adj_list.size() * 6);

    uint64_t activeLeaves = 0;
    uint32_t maxCleanPerLeaf = 0;
    uint32_t maxDeltaPerLeaf = 0;

    for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.size() < r) continue;

        int keepC = 0, pivotC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }
        const int needPivot = static_cast<int>(s) - keepC;
        const int low = std::max(0, static_cast<int>(r) - keepC);
        const int high = std::min(static_cast<int>(r), pivotC);

        auto &idx = leafIndex[leafId];
        idx.cleanBegin = cleanEntries.size();
        idx.deltaBegin = deltaEntries.size();

        for (int qP = low; qP <= high; ++qP) {
            if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
            int row = pivotC - qP;
            int col = needPivot - qP;
            cleanEntries.push_back({
                static_cast<uint8_t>(qP),
                combSmall(pivotC, qP) * combSmall(keepC, static_cast<int>(r) - qP),
                nCr[row][col]
            });
            idx.cleanCount++;
        }

        for (int remP = low; remP <= high; ++remP) {
            if (!hasPositiveContribution(pivotC, keepC, remP, r, s)) continue;
            for (int qP = low; qP <= high; ++qP) {
                if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
                const int xLow = std::max(0, qP - (pivotC - remP));
                const int xHigh = std::min(qP, remP);
                for (int x = xLow; x <= xHigh; ++x) {
                    uint64_t mult = combSmall(remP, x) *
                                    combSmall(pivotC - remP, qP - x) *
                                    combSmall(keepC, static_cast<int>(r) - qP);
                    if (mult == 0) continue;
                    const int unionP = remP + qP - x;
                    if (unionP > needPivot) continue;
                    const int row = pivotC - unionP;
                    const int col = needPivot - unionP;
                    if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) continue;
                    double delta = nCr[row][col];
                    if (delta <= 0.0) continue;
                    deltaEntries.push_back({
                        static_cast<uint8_t>(remP),
                        static_cast<uint8_t>(qP),
                        static_cast<uint8_t>(x),
                        mult,
                        delta
                    });
                    idx.deltaCount++;
                }
            }
        }

        if (idx.cleanCount > 0 || idx.deltaCount > 0) {
            activeLeaves++;
            maxCleanPerLeaf = std::max(maxCleanPerLeaf, idx.cleanCount);
            maxDeltaPerLeaf = std::max(maxDeltaPerLeaf, idx.deltaCount);
        }
    }

    uint64_t totalBytes =
        leafIndex.size() * sizeof(QuotientLeafStateIndex) +
        cleanEntries.size() * sizeof(QuotientCleanStateEntry) +
        deltaEntries.size() * sizeof(QuotientDeltaStateEntry);

    std::cout << "================ Quotient Prototype State =================" << std::endl;
    std::cout << "  Active leaves with state: " << activeLeaves << std::endl;
    std::cout << "  Clean state entries:      " << cleanEntries.size() << std::endl;
    std::cout << "  Delta state entries:      " << deltaEntries.size() << std::endl;
    std::cout << "  Max clean entries/leaf:   " << maxCleanPerLeaf << std::endl;
    std::cout << "  Max delta entries/leaf:   " << maxDeltaPerLeaf << std::endl;
    if (activeLeaves > 0) {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "  Avg clean entries/leaf:   "
                  << static_cast<long double>(cleanEntries.size()) / activeLeaves << std::endl;
        std::cout << "  Avg delta entries/leaf:   "
                  << static_cast<long double>(deltaEntries.size()) / activeLeaves << std::endl;
        std::cout.unsetf(std::ios::floatfield);
    }
    std::cout << "  Prototype state memory:   "
              << totalBytes / 1024 / 1024 << " MB" << std::endl;
    std::cout << "===========================================================" << std::endl;
}

static DeltaStateSummary deltaStateForRemovedPivotCount(
    int keepC, int pivotC, daf::CliqueSize r, daf::CliqueSize s, int remP) {
    const int needPivot = static_cast<int>(s) - keepC;
    const int low = std::max(0, static_cast<int>(r) - keepC);
    const int high = std::min(static_cast<int>(r), pivotC);
    if (remP < low || remP > high) return {};
    if (!hasPositiveContribution(pivotC, keepC, remP, r, s)) return {};

    uint64_t classCnt = 0;
    uint64_t affectedCnt = 0;
    uint64_t maxMultiplicity = 0;

    for (int qP = low; qP <= high; ++qP) {
        if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
        const int xLow = std::max(0, qP - (pivotC - remP));
        const int xHigh = std::min(qP, remP);
        for (int x = xLow; x <= xHigh; ++x) {
            uint64_t mult = combSmall(remP, x) *
                            combSmall(pivotC - remP, qP - x) *
                            combSmall(keepC, static_cast<int>(r) - qP);
            if (qP == remP && x == remP && mult > 0) mult -= 1;
            if (mult == 0) continue;
            const int unionP = remP + qP - x;
            if (unionP > needPivot) continue;
            const int row = pivotC - unionP;
            const int col = needPivot - unionP;
            if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) continue;
            if (nCr[row][col] <= 0.0) continue;
            classCnt++;
            affectedCnt += mult;
            maxMultiplicity = std::max(maxMultiplicity, mult);
        }
    }
    return {classCnt, affectedCnt, maxMultiplicity};
}

static void maybeRunOneRoundPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {
    if (!std::getenv("PIVOTER_QUOTIENT_ONE_ROUND")) return;

    try {
    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "================ Quotient One-Round Prototype ================" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip fused global clique-id build" << std::endl;
        std::cout << "  Reason: raw clique layer is too large for this prototype" << std::endl;
        std::cout << "===============================================================" << std::endl;
        return;
    }

    (void)prebuiltIndex;
    StaticCliqueIndex cliqueIndex(r);
    std::vector<int> leafPivotC(tree.adj_list.size(), 0);
    std::vector<int> leafNeedPivot(tree.adj_list.size(), 0);
    for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
        const auto &leaf = tree.adj_list[leafId];
        int pivotC = 0, keepC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }
        leafPivotC[leafId] = pivotC;
        leafNeedPivot[leafId] = static_cast<int>(s) - keepC;
    }

    std::vector<double> countingRClique;
    std::vector<std::unordered_map<uint64_t, daf::Size>> leafLocalCliqueId(tree.adj_list.size());
    daf::timeCount("fused build+counting (Quotient one-round)", [&]() {
        cliqueIndex.buildWithFullEnum(tree, edgeGraph.adj_list.size(),
            [&](daf::Size leafIdx, StaticCliqueIndex::Id cliqueId, daf::CliqueSize subNumPivot,
                const uint8_t *positions) {
                if (cliqueId >= countingRClique.size()) {
                    countingRClique.resize(cliqueId + 1, 0.0);
                }
                leafLocalCliqueId[leafIdx].emplace(packLeafPosKey(positions, r), cliqueId);
                if (!hasPositiveContribution(
                        leafPivotC[leafIdx],
                        static_cast<int>(s) - leafNeedPivot[leafIdx],
                        subNumPivot, r, s)) {
                    return;
                }
                const int row = leafPivotC[leafIdx] - static_cast<int>(subNumPivot);
                const int col = leafNeedPivot[leafIdx] - static_cast<int>(subNumPivot);
                if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return;
                countingRClique[cliqueId] += nCr[row][col];
            });
    });
    std::cout << "Quotient one-round: counting ready" << std::endl;

    if (countingRClique.empty()) return;
    double minCore = *std::min_element(countingRClique.begin(), countingRClique.end());
    std::vector<daf::Size> currentRemoveRcliqueIds;
    currentRemoveRcliqueIds.reserve(countingRClique.size() / 16 + 1);
    for (daf::Size i = 0; i < countingRClique.size(); ++i) {
        if (countingRClique[i] <= minCore) currentRemoveRcliqueIds.push_back(i);
    }
    std::cout << "Quotient one-round: first frontier ready" << std::endl;

    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<std::vector<daf::Size>> removedRCliqueIdForLeaf;
    std::vector<daf::Size> changedLeaf;
    removedRCliqueIdForLeaf.reserve(tree.adj_list.size() / 16 + 1);
    changedLeaf.reserve(tree.adj_list.size() / 16 + 1);

    for (auto rmRCliqueId : currentRemoveRcliqueIds) {
        auto rClique = cliqueIndex.byId(rmRCliqueId);
        daf::intersect_dense_sets_multi(rClique, treeGraphV.adj_list,
            [&](const TreeGraphNode &uClique) {
                auto &leafIdx = changedLeafIndex[uClique.v];
                if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                    leafIdx = removedRCliqueIdForLeaf.size();
                    removedRCliqueIdForLeaf.emplace_back();
                    changedLeaf.push_back(uClique.v);
                }
                removedRCliqueIdForLeaf[leafIdx].push_back(rmRCliqueId);
            });
    }
    std::cout << "Quotient one-round: changed leaves ready" << std::endl;

    uint64_t handledLeaves = 0;
    uint64_t singleRemovedLeaves = 0;
    uint64_t multiRemovedLeaves = 0;
    uint64_t handledDeltaClasses = 0;
    uint64_t handledTouchedCliques = 0;
    uint64_t handledExplicitEntries = 0;
    uint64_t emitMiss = 0;
    std::vector<long double> handledTouchFractions;
    handledTouchFractions.reserve(changedLeaf.size());
    std::unordered_map<daf::Size, double> emittedDelta;
    emittedDelta.reserve(1 << 20);
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> sparseTupleDelta;
    sparseTupleDelta.reserve(1 << 20);

    for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
        daf::Size leafId = changedLeaf[idx];
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.size() < r) continue;
        const auto &removedR = removedRCliqueIdForLeaf[changedLeafIndex[leafId]];
        if (removedR.size() != 1) {
            if (!removedR.empty()) multiRemovedLeaves++;
            continue;
        }
        singleRemovedLeaves++;

        int keepC = 0, pivotC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }
        const int low = std::max(0, static_cast<int>(r) - keepC);
        const int high = std::min(static_cast<int>(r), pivotC);

        auto rmClique = cliqueIndex.byId(removedR[0]);
        daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
        for (int i = 0; i < static_cast<int>(leaf.size()); ++i) mapRef[leaf[i].v] = static_cast<daf::Size>(i);

        int remP = 0;
        for (auto v : rmClique) {
            daf::Size pos = mapRef[v];
            if (pos >= leaf.size()) {
                remP = -1;
                break;
            }
            if (leaf[pos].isPivot) remP++;
        }
        if (remP < low || remP > high) continue;

        uint64_t explicitEntries = 0;
        for (int qP = low; qP <= high; ++qP) {
            if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
            explicitEntries += combSmall(pivotC, qP) *
                               combSmall(keepC, static_cast<int>(r) - qP);
        }
        if (explicitEntries == 0) continue;

        auto deltaSummary = deltaStateForRemovedPivotCount(keepC, pivotC, r, s, remP);
        if (deltaSummary.classCount == 0 || deltaSummary.affectedCount == 0) continue;

        handledLeaves++;
        handledDeltaClasses += deltaSummary.classCount;
        handledTouchedCliques += deltaSummary.affectedCount;
        handledExplicitEntries += explicitEntries;
        handledTouchFractions.push_back(
            static_cast<long double>(deltaSummary.affectedCount) /
            static_cast<long double>(explicitEntries));

        const int needPivot = static_cast<int>(s) - keepC;
        for (int qP = low; qP <= high; ++qP) {
            if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
            const int xLow = std::max(0, qP - (pivotC - remP));
            const int xHigh = std::min(qP, remP);
            for (int x = xLow; x <= xHigh; ++x) {
                const int unionP = remP + qP - x;
                if (unionP > needPivot) continue;
                const int row = pivotC - unionP;
                const int col = needPivot - unionP;
                if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) continue;
                double delta = nCr[row][col];
                if (delta <= 0.0) continue;

                forEachExpandedCanonicalDeltaClass(leaf, r, remP, qP, x,
                    [&](const std::vector<uint8_t> &pos) {
                        auto it = leafLocalCliqueId[leafId].find(packLeafPosKey(pos));
                        if (it == leafLocalCliqueId[leafId].end()) {
                            emitMiss++;
                            return;
                        }
                        sparseTupleDelta[makeSparseCliqueKey(leaf, pos, r)] += delta;
                        emittedDelta[it->second] += delta;
                    });
            }
        }
    }
    std::cout << "Quotient one-round: delta emission ready" << std::endl;

    std::sort(handledTouchFractions.begin(), handledTouchFractions.end());
    auto pickTouch = [&](double q) -> long double {
        if (handledTouchFractions.empty()) return 0.0;
        size_t idx = std::min(
            handledTouchFractions.size() - 1,
            static_cast<size_t>(q * static_cast<double>(handledTouchFractions.size() - 1)));
        return handledTouchFractions[idx];
    };

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "================ Quotient One-Round Prototype ================" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  First-round min core:      " << minCore << std::endl;
    std::cout << "  Removed r-cliques:         " << currentRemoveRcliqueIds.size() << std::endl;
    std::cout << "  Changed leaves:            " << changedLeaf.size() << std::endl;
    std::cout << "  Single-removed leaves:     " << singleRemovedLeaves << std::endl;
    std::cout << "  Multi-removed leaves:      " << multiRemovedLeaves << std::endl;
    std::cout << "  Quotient-handled leaves:   " << handledLeaves << std::endl;
    std::cout << "  Handled delta classes:     " << handledDeltaClasses << std::endl;
    std::cout << "  Handled touched cliques:   " << handledTouchedCliques << std::endl;
    std::cout << "  Handled explicit entries:  " << handledExplicitEntries << std::endl;
    std::cout << "  Unique clique deltas:      " << emittedDelta.size() << std::endl;
    std::cout << "  Unique sparse deltas:      " << sparseTupleDelta.size() << std::endl;
    std::cout << "  Counting lookup misses:    0" << std::endl;
    std::cout << "  Delta lookup misses:       " << emitMiss << std::endl;
    std::cout << "  Sparse/global delta agree: "
              << (sparseTupleDelta.size() == emittedDelta.size() ? "YES" : "NO")
              << std::endl;
    if (handledExplicitEntries > 0) {
        std::cout << std::fixed << std::setprecision(6)
                  << "  Handled touch fraction:    "
                  << (100.0L * static_cast<long double>(handledTouchedCliques) /
                      static_cast<long double>(handledExplicitEntries))
                  << "%" << std::endl
                  << std::setprecision(4)
                  << "  Per-leaf touch median/p90: "
                  << (100.0L * pickTouch(0.50)) << "% / "
                  << (100.0L * pickTouch(0.90)) << "%" << std::endl
                  << std::setprecision(2);
    }
    std::cout << "  Prototype time:            " << elapsedMs << " ms" << std::endl;
    std::cout << "===============================================================" << std::endl;
    } catch (const std::exception &e) {
        std::cout << "Quotient one-round prototype exception: "
                  << e.what() << std::endl;
    }
}

static void maybeRunOneRoundSparseSupportPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_SPARSE_SUPPORT")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "============ Quotient One-Round Sparse Support ==============" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip sparse support build" << std::endl;
        std::cout << "  Reason: global sparse support map would still be too large" << std::endl;
        std::cout << "=============================================================" << std::endl;
        return;
    }

    std::vector<int> leafPivotC(tree.adj_list.size(), 0);
    std::vector<int> leafKeepC(tree.adj_list.size(), 0);
    std::vector<int> leafNeedPivot(tree.adj_list.size(), 0);
    for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
        const auto &leaf = tree.adj_list[leafId];
        int pivotC = 0, keepC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }
        leafPivotC[leafId] = pivotC;
        leafKeepC[leafId] = keepC;
        leafNeedPivot[leafId] = static_cast<int>(s) - keepC;
    }

    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> supportByKey;
    supportByKey.reserve(static_cast<size_t>(std::min<uint64_t>(rawCliqueCount, 10000000ULL)));
    daf::timeCount("sparse support build (Quotient one-round)", [&]() {
        for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
            const auto &leaf = tree.adj_list[leafId];
            if (leaf.size() < r) continue;
            daf::enumerateCombinations(leaf, r,
                [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                    int subNumPivot = 0;
                    SparseCliqueKey key;
                    key.len = static_cast<uint8_t>(r);
                    for (daf::CliqueSize j = 0; j < r; ++j) {
                        if (rClique[j].isPivot) subNumPivot++;
                        key.verts[j] = rClique[j].v;
                    }
                    if (!hasPositiveContribution(leafPivotC[leafId], leafKeepC[leafId], subNumPivot, r, s)) {
                        return true;
                    }
                    std::sort(key.verts.begin(), key.verts.begin() + r);
                    const int row = leafPivotC[leafId] - subNumPivot;
                    const int col = leafNeedPivot[leafId] - subNumPivot;
                    if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return true;
                    supportByKey[key] += nCr[row][col];
                    return true;
                });
        }
    });

    if (supportByKey.empty()) return;

    double minCore = std::numeric_limits<double>::max();
    for (const auto &[key, value] : supportByKey) {
        (void)key;
        minCore = std::min(minCore, value);
    }

    std::vector<SparseCliqueKey> currentRemoveKeys;
    currentRemoveKeys.reserve(supportByKey.size() / 16 + 1);
    for (const auto &[key, value] : supportByKey) {
        if (value <= minCore) currentRemoveKeys.push_back(key);
    }

    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
    std::vector<daf::Size> changedLeaf;
    removedKeyForLeaf.reserve(tree.adj_list.size() / 16 + 1);
    changedLeaf.reserve(tree.adj_list.size() / 16 + 1);

    daf::StaticVector<daf::Size> rmVerts;
    rmVerts.resize(r);
    for (const auto &rmKey : currentRemoveKeys) {
        for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
        daf::intersect_dense_sets_multi(rmVerts, treeGraphV.adj_list,
            [&](const TreeGraphNode &uClique) {
                auto &leafIdx = changedLeafIndex[uClique.v];
                if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                    leafIdx = removedKeyForLeaf.size();
                    removedKeyForLeaf.emplace_back();
                    changedLeaf.push_back(uClique.v);
                }
                removedKeyForLeaf[leafIdx].push_back(rmKey);
            });
    }

    uint64_t handledLeaves = 0;
    uint64_t singleRemovedLeaves = 0;
    uint64_t multiRemovedLeaves = 0;
    uint64_t handledDeltaClasses = 0;
    uint64_t handledTouchedCliques = 0;
    uint64_t handledExplicitEntries = 0;
    std::vector<long double> handledTouchFractions;
    handledTouchFractions.reserve(changedLeaf.size());
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> sparseDelta;
    sparseDelta.reserve(1 << 20);

    for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
        const daf::Size leafId = changedLeaf[idx];
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.size() < r) continue;
        const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
        if (removed.size() != 1) {
            if (!removed.empty()) multiRemovedLeaves++;
            continue;
        }
        singleRemovedLeaves++;

        const int keepC = leafKeepC[leafId];
        const int pivotC = leafPivotC[leafId];
        const int low = std::max(0, static_cast<int>(r) - keepC);
        const int high = std::min(static_cast<int>(r), pivotC);

        daf::StaticVector<daf::Size> &mapRef = daf::vListMap;
        for (int i = 0; i < static_cast<int>(leaf.size()); ++i) {
            mapRef[leaf[i].v] = static_cast<daf::Size>(i);
        }

        int remP = 0;
        bool validRemoved = true;
        for (daf::CliqueSize j = 0; j < r; ++j) {
            daf::Size pos = mapRef[removed[0].verts[j]];
            if (pos >= leaf.size()) {
                validRemoved = false;
                break;
            }
            if (leaf[pos].isPivot) remP++;
        }
        if (!validRemoved || remP < low || remP > high) continue;

        uint64_t explicitEntries = 0;
        for (int qP = low; qP <= high; ++qP) {
            if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
            explicitEntries += combSmall(pivotC, qP) *
                               combSmall(keepC, static_cast<int>(r) - qP);
        }
        if (explicitEntries == 0) continue;

        auto deltaSummary = deltaStateForRemovedPivotCount(keepC, pivotC, r, s, remP);
        if (deltaSummary.classCount == 0 || deltaSummary.affectedCount == 0) continue;

        handledLeaves++;
        handledDeltaClasses += deltaSummary.classCount;
        handledTouchedCliques += deltaSummary.affectedCount;
        handledExplicitEntries += explicitEntries;
        handledTouchFractions.push_back(
            static_cast<long double>(deltaSummary.affectedCount) /
            static_cast<long double>(explicitEntries));

        const int needPivot = static_cast<int>(s) - keepC;
        for (int qP = low; qP <= high; ++qP) {
            if (!hasPositiveContribution(pivotC, keepC, qP, r, s)) continue;
            const int xLow = std::max(0, qP - (pivotC - remP));
            const int xHigh = std::min(qP, remP);
            for (int x = xLow; x <= xHigh; ++x) {
                const int unionP = remP + qP - x;
                if (unionP > needPivot) continue;
                const int row = pivotC - unionP;
                const int col = needPivot - unionP;
                if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) continue;
                const double delta = nCr[row][col];
                if (delta <= 0.0) continue;

                forEachExpandedCanonicalDeltaClass(leaf, r, remP, qP, x,
                    [&](const std::vector<uint8_t> &pos) {
                        sparseDelta[makeSparseCliqueKey(leaf, pos, r)] += delta;
                    });
            }
        }
    }
    rmVerts.free();

    std::sort(handledTouchFractions.begin(), handledTouchFractions.end());
    auto pickTouch = [&](double q) -> long double {
        if (handledTouchFractions.empty()) return 0.0;
        size_t idx = std::min(
            handledTouchFractions.size() - 1,
            static_cast<size_t>(q * static_cast<double>(handledTouchFractions.size() - 1)));
        return handledTouchFractions[idx];
    };

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "============ Quotient One-Round Sparse Support ==============" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Sparse support states:     " << supportByKey.size() << std::endl;
    std::cout << "  First-round min core:      " << minCore << std::endl;
    std::cout << "  Removed sparse cliques:    " << currentRemoveKeys.size() << std::endl;
    std::cout << "  Changed leaves:            " << changedLeaf.size() << std::endl;
    std::cout << "  Single-removed leaves:     " << singleRemovedLeaves << std::endl;
    std::cout << "  Multi-removed leaves:      " << multiRemovedLeaves << std::endl;
    std::cout << "  Quotient-handled leaves:   " << handledLeaves << std::endl;
    std::cout << "  Handled delta classes:     " << handledDeltaClasses << std::endl;
    std::cout << "  Handled touched cliques:   " << handledTouchedCliques << std::endl;
    std::cout << "  Handled explicit entries:  " << handledExplicitEntries << std::endl;
    std::cout << "  Unique sparse deltas:      " << sparseDelta.size() << std::endl;
    if (handledExplicitEntries > 0) {
        std::cout << std::fixed << std::setprecision(6)
                  << "  Handled touch fraction:    "
                  << (100.0L * static_cast<long double>(handledTouchedCliques) /
                      static_cast<long double>(handledExplicitEntries))
                  << "%" << std::endl
                  << std::setprecision(4)
                  << "  Per-leaf touch median/p90: "
                  << (100.0L * pickTouch(0.50)) << "% / "
                  << (100.0L * pickTouch(0.90)) << "%" << std::endl
                  << std::setprecision(2);
    }
    std::cout << "  Prototype time:            " << elapsedMs << " ms" << std::endl;
    std::cout << "=============================================================" << std::endl;
}

#include "NucleusCoreDecompositionRCliqueST_QuotientLab_FrontierAnalysis.inc"

struct LowSupportBuildStats {
    uint64_t keptUpdates = 0;
    uint64_t evictedUpdates = 0;
    uint64_t skippedUpdates = 0;
};

struct LowSupportBuildResult {
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> lowSupport;
    robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> overTau;
    LowSupportBuildStats stats;
};

struct AdaptiveLowPhaseStats {
    double tau = 0.0;
    double capTau = 0.0;
    double windowCapTau = 0.0;
    int lookahead = 0;
    bool fullBand = false;
    bool twoPhaseBuild = false;
    bool cacheRebandBuild = false;
    bool singleOverMapBuild = false;
    bool cappedLowBuild = false;
    uint64_t rebuildMs = 0;
    uint64_t initialLowKeys = 0;
    uint64_t initialBufferedKeys = 0;
    uint64_t initialOverExactKeys = 0;
    uint64_t initialOverTauKeys = 0;
    uint64_t initialTrackedOverKeys = 0;
    LowSupportBuildStats buildStats;
    int rounds = 0;
    uint64_t frontierKeys = 0;
    uint64_t exactLeaves = 0;
    uint64_t candidateKeys = 0;
    uint64_t spawnedSubleaves = 0;
    double maxRoundMin = 0.0;
    uint64_t remainingLowKeys = 0;
    uint64_t phaseMs = 0;
    uint64_t leafUpdateMs = 0;
    uint64_t recomputeMs = 0;
    uint64_t bucketUpdateMs = 0;
    uint64_t deltaFastKeys = 0;
    uint64_t deltaPrunedKeys = 0;
    uint64_t deltaStateKeys = 0;
    uint64_t deltaOverKeys = 0;
    uint64_t deltaLeafBoundKeys = 0;
    uint64_t deltaOverLowerKeys = 0;
    uint64_t fullRecomputeKeys = 0;
    uint64_t overCacheHits = 0;
    uint64_t overCacheStores = 0;
    uint64_t zeroDeltaKeys = 0;
};

struct BandedLowSupportBuildResult {
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> lowSupport;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> bufferedSupport;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> overExactSupport;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> overExactSeed;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> overLowerBound;
    robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> overWindow;
    LowSupportBuildStats stats;
    bool usedSingleOverMap = false;
    bool usedCappedLowState = false;
};

struct BucketedBandedLowState {
    std::unordered_map<SparseCliqueKey, uint16_t, SparseCliqueKeyHash> lowValue;
    std::unordered_map<SparseCliqueKey, uint16_t, SparseCliqueKeyHash> bufferedValue;
    std::unordered_map<SparseCliqueKey, uint16_t, SparseCliqueKeyHash> overValue;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> overLowerBound;
    std::vector<robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash>> lowBuckets;
    std::vector<robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash>> bufferedBuckets;
    robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> overWindow;
};

static LowSupportBuildResult buildLowSupportLayer(
    const DynamicGraph<TreeGraphNode> &tree, daf::CliqueSize r, daf::CliqueSize s,
    double tau, uint64_t rawCliqueCount) {
    LowSupportBuildResult result;
    result.lowSupport.reserve(static_cast<size_t>(
        std::min<uint64_t>(rawCliqueCount / 8 + 1024, 5000000ULL)));
    result.overTau.reserve(static_cast<size_t>(
        std::min<uint64_t>(rawCliqueCount / 4 + 1024, 10000000ULL)));

    for (const auto &leaf : tree.adj_list) {
        if (leaf.size() < r) continue;
        int keepC = 0, pivotC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }
        daf::enumerateCombinations(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                int subNumPivot = 0;
                SparseCliqueKey key;
                key.len = static_cast<uint8_t>(r);
                for (daf::CliqueSize j = 0; j < r; ++j) {
                    if (rClique[j].isPivot) subNumPivot++;
                    key.verts[j] = rClique[j].v;
                }
                if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return true;
                const int row = pivotC - subNumPivot;
                const int col = static_cast<int>(s) - keepC - subNumPivot;
                if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return true;
                const double contrib = nCr[row][col];
                if (contrib <= 0.0) return true;
                std::sort(key.verts.begin(), key.verts.begin() + r);
                if (result.overTau.find(key) != result.overTau.end()) {
                    result.stats.skippedUpdates++;
                    return true;
                }
                auto it = result.lowSupport.find(key);
                if (it == result.lowSupport.end()) {
                    if (contrib <= tau) {
                        result.lowSupport.emplace(key, contrib);
                        result.stats.keptUpdates++;
                    } else {
                        result.overTau.insert(key);
                        result.stats.evictedUpdates++;
                    }
                    return true;
                }
                const double newValue = it->second + contrib;
                if (newValue <= tau) {
                    it->second = newValue;
                    result.stats.keptUpdates++;
                } else {
                    result.lowSupport.erase(it);
                    result.overTau.insert(key);
                    result.stats.evictedUpdates++;
                }
                return true;
            });
    }
    return result;
}

static LowSupportBuildResult buildLowSupportLayerActive(
    const std::vector<std::vector<TreeGraphNode>> &activeLeaves,
    const std::vector<uint8_t> &leafAlive,
    daf::CliqueSize r, daf::CliqueSize s,
    double tau, uint64_t reserveHint) {
    LowSupportBuildResult result;
    result.lowSupport.reserve(static_cast<size_t>(
        std::min<uint64_t>(reserveHint / 8 + 1024, 5000000ULL)));
    result.overTau.reserve(static_cast<size_t>(
        std::min<uint64_t>(reserveHint / 4 + 1024, 10000000ULL)));

    for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
        if (!leafAlive[leafId]) continue;
        const auto &leaf = activeLeaves[leafId];
        if (leaf.size() < r) continue;
        int keepC = 0, pivotC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }
        daf::enumerateCombinations(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                int subNumPivot = 0;
                SparseCliqueKey key;
                key.len = static_cast<uint8_t>(r);
                for (daf::CliqueSize j = 0; j < r; ++j) {
                    if (rClique[j].isPivot) subNumPivot++;
                    key.verts[j] = rClique[j].v;
                }
                if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return true;
                const int row = pivotC - subNumPivot;
                const int col = static_cast<int>(s) - keepC - subNumPivot;
                if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return true;
                const double contrib = nCr[row][col];
                if (contrib <= 0.0) return true;
                std::sort(key.verts.begin(), key.verts.begin() + r);
                if (result.overTau.find(key) != result.overTau.end()) {
                    result.stats.skippedUpdates++;
                    return true;
                }
                auto it = result.lowSupport.find(key);
                if (it == result.lowSupport.end()) {
                    if (contrib <= tau) {
                        result.lowSupport.emplace(key, contrib);
                        result.stats.keptUpdates++;
                    } else {
                        result.overTau.insert(key);
                        result.stats.evictedUpdates++;
                    }
                    return true;
                }
                const double newValue = it->second + contrib;
                if (newValue <= tau) {
                    it->second = newValue;
                    result.stats.keptUpdates++;
                } else {
                    result.lowSupport.erase(it);
                    result.overTau.insert(key);
                    result.stats.evictedUpdates++;
                }
                return true;
            });
    }
    return result;
}

// Variant that keeps exact support for over-tau keys (eliminates first-time full recompute)
struct ExactCacheBuildResult {
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> lowSupport;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> overTauExact;
    LowSupportBuildStats stats;
};

static ExactCacheBuildResult buildLowSupportLayerActiveWithExactCache(
    const std::vector<std::vector<TreeGraphNode>> &activeLeaves,
    const std::vector<uint8_t> &leafAlive,
    daf::CliqueSize r, daf::CliqueSize s,
    double tau, uint64_t reserveHint) {
    ExactCacheBuildResult result;
    result.lowSupport.reserve(static_cast<size_t>(
        std::min<uint64_t>(reserveHint / 8 + 1024, 5000000ULL)));
    result.overTauExact.reserve(static_cast<size_t>(
        std::min<uint64_t>(reserveHint / 4 + 1024, 10000000ULL)));

    for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
        if (!leafAlive[leafId]) continue;
        const auto &leaf = activeLeaves[leafId];
        if (leaf.size() < r) continue;
        int keepC = 0, pivotC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }
        daf::enumerateCombinations(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                int subNumPivot = 0;
                SparseCliqueKey key;
                key.len = static_cast<uint8_t>(r);
                for (daf::CliqueSize j = 0; j < r; ++j) {
                    if (rClique[j].isPivot) subNumPivot++;
                    key.verts[j] = rClique[j].v;
                }
                if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return true;
                const int row = pivotC - subNumPivot;
                const int col = static_cast<int>(s) - keepC - subNumPivot;
                if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return true;
                const double contrib = nCr[row][col];
                if (contrib <= 0.0) return true;
                std::sort(key.verts.begin(), key.verts.begin() + r);
                // Keep accumulating for over-tau keys (don't skip)
                auto overIt = result.overTauExact.find(key);
                if (overIt != result.overTauExact.end()) {
                    overIt->second += contrib;
                    result.stats.keptUpdates++;
                    return true;
                }
                auto it = result.lowSupport.find(key);
                if (it == result.lowSupport.end()) {
                    if (contrib <= tau) {
                        result.lowSupport.emplace(key, contrib);
                        result.stats.keptUpdates++;
                    } else {
                        result.overTauExact.emplace(key, contrib);
                        result.stats.evictedUpdates++;
                    }
                    return true;
                }
                const double newValue = it->second + contrib;
                if (newValue <= tau) {
                    it->second = newValue;
                    result.stats.keptUpdates++;
                } else {
                    result.lowSupport.erase(it);
                    result.overTauExact.emplace(key, newValue);
                    result.stats.evictedUpdates++;
                }
                return true;
            });
    }
    return result;
}

// Two-phase low-support build for large graphs.
//
// Problem: standard build enumerates C(|L|,r) per leaf. On web-it-2004
// this is 28.8B combinations — infeasible even with skip optimization
// because the combination GENERATION itself (not hash ops) dominates.
//
// Phase A (leaf-centric): process a small number of leaves to build a
//   large overTau set. After a few leaves, most keys are in overTau.
//   Stop when new key rate drops below a threshold.
//
// Phase B (key-centric): for each remaining lowSupport candidate key,
//   verify exact support via vertex-to-leaf intersection. This is
//   O(|lowSupport| × intersection_cost), which is fast when lowSupport
//   is small.
//
// On web-it-2004: Phase A processes ~100 leaves (6.6s), Phase B verifies
// ~200K keys (2s). Total ~9s vs infeasible 28,000s for full enumeration.
static LowSupportBuildResult buildLowSupportTwoPhase(
    const std::vector<std::vector<TreeGraphNode>> &activeLeaves,
    const std::vector<uint8_t> &leafAlive,
    std::vector<robin_hood::unordered_flat_set<daf::Size>> &activeLeafByVertex,
    daf::CliqueSize r, daf::CliqueSize s,
    double tau, uint64_t reserveHint) {
    LowSupportBuildResult result;
    result.lowSupport.reserve(1 << 20);
    result.overTau.reserve(static_cast<size_t>(
        std::min<uint64_t>(reserveHint / 4 + 1024, 10000000ULL)));

    // Collect active leaf indices sorted by size descending (process densest first)
    std::vector<daf::Size> leafOrder;
    leafOrder.reserve(activeLeaves.size());
    for (daf::Size i = 0; i < activeLeaves.size(); ++i) {
        if (leafAlive[i] && activeLeaves[i].size() >= r) leafOrder.push_back(i);
    }
    std::sort(leafOrder.begin(), leafOrder.end(), [&](daf::Size a, daf::Size b) {
        return activeLeaves[a].size() > activeLeaves[b].size();
    });

    auto tPhaseA = std::chrono::high_resolution_clock::now();

    // Phase A: leaf-centric enumeration with class-level skip + overTau skip
    uint64_t phaseALeaves = 0;
    uint64_t phaseAEnumerated = 0;
    uint64_t phaseAClassSkipped = 0;
    uint64_t prevLowSize = 0;
    const uint64_t minPhaseALeaves = 10;
    const uint64_t checkInterval = 10;

    for (daf::Size li = 0; li < leafOrder.size(); ++li) {
        const daf::Size leafId = leafOrder[li];
        const auto &leaf = activeLeaves[leafId];
        int keepC = 0, pivotC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }

        // Precompute per-class contribution values
        std::vector<double> classContrib(static_cast<size_t>(r) + 1, 0.0);
        bool anyClassNeeded = false;
        for (int p = 0; p <= static_cast<int>(r); ++p) {
            const int needKeep = static_cast<int>(r) - p;
            if (p > pivotC || needKeep > keepC) continue;
            if (!hasPositiveContribution(pivotC, keepC, p, r, s)) continue;
            const int row = pivotC - p;
            const int col = static_cast<int>(s) - keepC - p;
            if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) continue;
            const double c = nCr[row][col];
            if (c <= 0.0) continue;
            if (c > tau) {
                phaseAClassSkipped += static_cast<uint64_t>(combSmall(pivotC, p))
                                   * combSmall(keepC, needKeep);
                continue;
            }
            classContrib[p] = c;
            anyClassNeeded = true;
        }
        if (!anyClassNeeded) { phaseALeaves++; continue; }

        // Enumerate once per leaf, filtering by class contribution
        daf::enumerateCombinations(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                int snp = 0;
                for (daf::CliqueSize j = 0; j < r; ++j) {
                    if (rClique[j].isPivot) snp++;
                }
                const double contrib = classContrib[snp];
                if (contrib <= 0.0) return true; // class skipped
                SparseCliqueKey key;
                key.len = static_cast<uint8_t>(r);
                for (daf::CliqueSize j = 0; j < r; ++j) key.verts[j] = rClique[j].v;
                std::sort(key.verts.begin(), key.verts.begin() + r);
                phaseAEnumerated++;

                if (result.overTau.find(key) != result.overTau.end()) {
                    result.stats.skippedUpdates++;
                    return true;
                }
                auto it = result.lowSupport.find(key);
                if (it == result.lowSupport.end()) {
                    result.lowSupport.emplace(key, contrib);
                    result.stats.keptUpdates++;
                } else {
                    const double nv = it->second + contrib;
                    if (nv <= tau) {
                        it->second = nv;
                        result.stats.keptUpdates++;
                    } else {
                        result.lowSupport.erase(it);
                        result.overTau.insert(key);
                        result.stats.evictedUpdates++;
                    }
                }
                return true;
            });

        phaseALeaves++;

        // Check convergence: if lowSupport stopped changing, switch to Phase B
        if (phaseALeaves >= minPhaseALeaves && phaseALeaves % checkInterval == 0) {
            const uint64_t curLow = result.lowSupport.size();
            const double changeRate = (prevLowSize > 0)
                ? std::abs(static_cast<double>(curLow) - static_cast<double>(prevLowSize))
                  / static_cast<double>(prevLowSize)
                : 1.0;
            prevLowSize = curLow;
            if (changeRate < 0.001 && phaseALeaves >= 50) break; // stable enough
        }
    }

    auto tPhaseAEnd = std::chrono::high_resolution_clock::now();
    auto phaseAMs = std::chrono::duration_cast<std::chrono::milliseconds>(tPhaseAEnd - tPhaseA).count();

    // Phase B: key-centric verification
    // For each remaining lowSupport key, compute exact support via intersection
    uint64_t phaseBVerified = 0;
    uint64_t phaseBEvicted = 0;
    std::vector<SparseCliqueKey> toEvict;
    for (auto &[key, currentSupport] : result.lowSupport) {
        // Recompute exact support from ALL active leaves
        double exactSupport = 0.0;
        daf::StaticVector<daf::Size> verts;
        verts.resize(r);
        for (daf::CliqueSize j = 0; j < r; ++j) verts[j] = key.verts[j];
        daf::intersect_dense_sets_multi(verts, activeLeafByVertex,
            [&](const daf::Size &leafId) {
                if (leafId < activeLeaves.size() && leafAlive[leafId]) {
                    exactSupport += leafContributionForKey(activeLeaves[leafId], key, r, s);
                }
            });
        verts.free();
        phaseBVerified++;
        if (exactSupport > tau) {
            toEvict.push_back(key);
            phaseBEvicted++;
        } else {
            currentSupport = exactSupport;
        }
    }
    for (const auto &key : toEvict) {
        result.lowSupport.erase(key);
        result.overTau.insert(key);
    }

    auto tPhaseBEnd = std::chrono::high_resolution_clock::now();
    auto phaseBMs = std::chrono::duration_cast<std::chrono::milliseconds>(tPhaseBEnd - tPhaseAEnd).count();

    std::cout << "  TwoPhase build:"
              << " phaseA_leaves=" << phaseALeaves
              << " phaseA_enum=" << phaseAEnumerated
              << " phaseA_classSkip=" << phaseAClassSkipped
              << " phaseA_ms=" << phaseAMs
              << " phaseB_verified=" << phaseBVerified
              << " phaseB_evicted=" << phaseBEvicted
              << " phaseB_ms=" << phaseBMs
              << " lowKeys=" << result.lowSupport.size()
              << " overKeys=" << result.overTau.size() << std::endl;

    return result;
}

static BandedLowSupportBuildResult buildBandedLowSupportLayerActive(
    const std::vector<std::vector<TreeGraphNode>> &activeLeaves,
    const std::vector<uint8_t> &leafAlive,
    daf::CliqueSize r, daf::CliqueSize s,
    double activeCapTau, double windowCapTau, uint64_t reserveHint,
    bool keepOverLowerBound,
    double exactOverCapTau = -1.0,
    bool keepExactOverSupport = false,
    double trackedOverCapTau = -1.0,
    bool preferSingleOverMap = false,
    bool preferCappedLowState = false) {
    BandedLowSupportBuildResult result;
    if (exactOverCapTau < windowCapTau) exactOverCapTau = windowCapTau;
    if (trackedOverCapTau < exactOverCapTau) {
        trackedOverCapTau = keepExactOverSupport && !keepOverLowerBound
            ? std::numeric_limits<double>::infinity()
            : exactOverCapTau;
    }
    const bool splitExactSeed =
        keepOverLowerBound && keepExactOverSupport &&
        std::isfinite(trackedOverCapTau) && trackedOverCapTau > exactOverCapTau + 1e-9;
    const double lowerBoundCap =
        splitExactSeed ? trackedOverCapTau + 1.0 : std::numeric_limits<double>::infinity();
    const bool largeReserve =
        std::getenv("PIVOTER_QUOTIENT_BUILD_LARGE_RESERVE") != nullptr;
    const bool requestCappedLowState =
        std::getenv("PIVOTER_QUOTIENT_BUILD_CAPPED_LOW") != nullptr;
    const bool requestSingleOverMap =
        std::getenv("PIVOTER_QUOTIENT_BUILD_SINGLE_OVER_MAP") != nullptr;
    const bool useCappedLowState =
        (requestCappedLowState || preferCappedLowState) &&
        exactOverCapTau <= 65535.0 &&
        std::abs(exactOverCapTau - std::round(exactOverCapTau)) <= 1e-9;
    const bool useSingleOverMap =
        (requestSingleOverMap || preferSingleOverMap) &&
        keepExactOverSupport &&
        !splitExactSeed &&
        trackedOverCapTau <= exactOverCapTau + 1e-9;
    result.usedSingleOverMap = useSingleOverMap;
    result.usedCappedLowState = useCappedLowState;
    const uint16_t exactCapInt = useCappedLowState
        ? static_cast<uint16_t>(std::llround(exactOverCapTau))
        : 0;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> exactWindow;
    robin_hood::unordered_flat_map<SparseCliqueKey, uint16_t, SparseCliqueKeyHash> cappedWindow;
    if (useCappedLowState) {
        cappedWindow.reserve(static_cast<size_t>(
            std::min<uint64_t>(reserveHint / 8 + 1024, largeReserve ? 20000000ULL : 5000000ULL)));
    } else {
        exactWindow.reserve(static_cast<size_t>(
            std::min<uint64_t>(reserveHint / 8 + 1024, largeReserve ? 20000000ULL : 5000000ULL)));
    }
    if (splitExactSeed) {
        result.overExactSeed.reserve(static_cast<size_t>(
            std::min<uint64_t>(reserveHint / 8 + 1024, largeReserve ? 20000000ULL : 5000000ULL)));
    }
    if (useSingleOverMap) {
        result.overLowerBound.reserve(static_cast<size_t>(
            std::min<uint64_t>(reserveHint / 4 + 1024, largeReserve ? 40000000ULL : 20000000ULL)));
    } else if (keepOverLowerBound) {
        result.overLowerBound.reserve(static_cast<size_t>(
            std::min<uint64_t>(reserveHint / 8 + 1024, largeReserve ? 20000000ULL : 5000000ULL)));
    }
    if (keepExactOverSupport && !keepOverLowerBound) {
        result.overLowerBound.reserve(static_cast<size_t>(
            std::min<uint64_t>(reserveHint / 4 + 1024, largeReserve ? 30000000ULL : 10000000ULL)));
    }
    if (!useSingleOverMap) {
        result.overWindow.reserve(static_cast<size_t>(
            std::min<uint64_t>(reserveHint / 4 + 1024, largeReserve ? 30000000ULL : 10000000ULL)));
    }
    const bool sortLeavesDesc =
        std::getenv("PIVOTER_QUOTIENT_BUILD_SORT_LEAVES") != nullptr;
    std::vector<daf::Size> leafOrder;
    if (sortLeavesDesc) {
        leafOrder.reserve(activeLeaves.size());
        for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
            if (leafAlive[leafId] && activeLeaves[leafId].size() >= r) {
                leafOrder.push_back(leafId);
            }
        }
        std::sort(leafOrder.begin(), leafOrder.end(),
                  [&](daf::Size a, daf::Size b) {
                      return activeLeaves[a].size() > activeLeaves[b].size();
                  });
    }

    auto processLeaf = [&](daf::Size leafId) {
        const auto &leaf = activeLeaves[leafId];
        if (leaf.size() < r) return;
        int keepC = 0, pivotC = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivotC++;
            else keepC++;
        }
        daf::enumerateCombinations(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                int subNumPivot = 0;
                SparseCliqueKey key;
                key.len = static_cast<uint8_t>(r);
                for (daf::CliqueSize j = 0; j < r; ++j) {
                    if (rClique[j].isPivot) subNumPivot++;
                    key.verts[j] = rClique[j].v;
                }
                if (!hasPositiveContribution(pivotC, keepC, subNumPivot, r, s)) return true;
                const int row = pivotC - subNumPivot;
                const int col = static_cast<int>(s) - keepC - subNumPivot;
                if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return true;
                const double contrib = nCr[row][col];
                if (contrib <= 0.0) return true;
                std::sort(key.verts.begin(), key.verts.begin() + r);
                if (useSingleOverMap) {
                    auto overIt = result.overLowerBound.find(key);
                    if (overIt != result.overLowerBound.end()) {
                        overIt->second += contrib;
                        result.stats.skippedUpdates++;
                        return true;
                    }
                } else if (result.overWindow.find(key) != result.overWindow.end()) {
                    if (keepOverLowerBound) {
                        auto &lower = result.overLowerBound[key];
                        lower = std::min(lowerBoundCap, lower + contrib);
                    }
                    if (splitExactSeed) {
                        auto exactIt = result.overExactSeed.find(key);
                        if (exactIt != result.overExactSeed.end()) {
                            const double updated = exactIt->second + contrib;
                            if (updated <= trackedOverCapTau) {
                                exactIt->second = updated;
                            } else {
                                result.overExactSeed.erase(exactIt);
                            }
                        }
                    } else if (keepExactOverSupport && !keepOverLowerBound) {
                        auto overIt = result.overLowerBound.find(key);
                        if (overIt != result.overLowerBound.end()) {
                            const double updated = overIt->second + contrib;
                            if (updated <= trackedOverCapTau) {
                                overIt->second = updated;
                            } else {
                                result.overLowerBound.erase(overIt);
                            }
                        }
                    }
                    result.stats.skippedUpdates++;
                    return true;
                }
                if (useCappedLowState) {
                    const uint16_t contribInt = static_cast<uint16_t>(std::llround(contrib));
                    auto it = cappedWindow.find(key);
                    if (it == cappedWindow.end()) {
                        if (contribInt <= exactCapInt) {
                            cappedWindow.emplace(key, contribInt);
                            result.stats.keptUpdates++;
                        } else {
                            if (useSingleOverMap) {
                                result.overLowerBound.emplace(key, static_cast<double>(contribInt));
                            } else if (keepOverLowerBound) {
                                result.overWindow.insert(key);
                                result.overLowerBound.emplace(
                                    key, std::min(lowerBoundCap, static_cast<double>(contribInt)));
                            } else {
                                result.overWindow.insert(key);
                            }
                            if (!useSingleOverMap &&
                                splitExactSeed &&
                                static_cast<double>(contribInt) <= trackedOverCapTau) {
                                result.overExactSeed.emplace(key, static_cast<double>(contribInt));
                            } else if (!useSingleOverMap &&
                                       !keepOverLowerBound &&
                                       keepExactOverSupport &&
                                       static_cast<double>(contribInt) <= trackedOverCapTau) {
                                result.overLowerBound.emplace(key, static_cast<double>(contribInt));
                            }
                            result.stats.evictedUpdates++;
                        }
                        return true;
                    }
                    const uint32_t newValueInt =
                        static_cast<uint32_t>(it->second) + static_cast<uint32_t>(contribInt);
                    if (newValueInt <= static_cast<uint32_t>(exactCapInt)) {
                        it->second = static_cast<uint16_t>(newValueInt);
                        result.stats.keptUpdates++;
                    } else {
                        cappedWindow.erase(it);
                        const double newValue = static_cast<double>(newValueInt);
                        if (useSingleOverMap) {
                            result.overLowerBound[key] = newValue;
                        } else if (keepOverLowerBound) {
                            result.overWindow.insert(key);
                            result.overLowerBound[key] = std::min(lowerBoundCap, newValue);
                        } else if (keepExactOverSupport && newValue <= trackedOverCapTau) {
                            result.overWindow.insert(key);
                            result.overLowerBound[key] = newValue;
                        } else {
                            result.overWindow.insert(key);
                            result.overLowerBound.erase(key);
                        }
                        if (!useSingleOverMap &&
                            splitExactSeed && newValue <= trackedOverCapTau) {
                            result.overExactSeed[key] = newValue;
                        }
                        result.stats.evictedUpdates++;
                    }
                } else {
                    auto it = exactWindow.find(key);
                    if (it == exactWindow.end()) {
                        if (contrib <= exactOverCapTau) {
                            exactWindow.emplace(key, contrib);
                            result.stats.keptUpdates++;
                        } else {
                            if (useSingleOverMap) {
                                result.overLowerBound.emplace(key, contrib);
                            } else if (keepOverLowerBound) {
                                result.overWindow.insert(key);
                                result.overLowerBound.emplace(key, std::min(lowerBoundCap, contrib));
                            } else {
                                result.overWindow.insert(key);
                            }
                            if (!useSingleOverMap &&
                                splitExactSeed && contrib <= trackedOverCapTau) {
                                result.overExactSeed.emplace(key, contrib);
                            } else if (!useSingleOverMap &&
                                       !keepOverLowerBound &&
                                       keepExactOverSupport && contrib <= trackedOverCapTau) {
                                result.overLowerBound.emplace(key, contrib);
                            }
                            result.stats.evictedUpdates++;
                        }
                        return true;
                    }
                    const double newValue = it->second + contrib;
                    if (newValue <= exactOverCapTau) {
                        it->second = newValue;
                        result.stats.keptUpdates++;
                    } else {
                        exactWindow.erase(it);
                        if (useSingleOverMap) {
                            result.overLowerBound[key] = newValue;
                        } else if (keepOverLowerBound) {
                            result.overWindow.insert(key);
                            result.overLowerBound[key] = std::min(lowerBoundCap, newValue);
                        } else if (keepExactOverSupport && newValue <= trackedOverCapTau) {
                            result.overWindow.insert(key);
                            result.overLowerBound[key] = newValue;
                        } else {
                            result.overWindow.insert(key);
                            result.overLowerBound.erase(key);
                        }
                        if (!useSingleOverMap &&
                            splitExactSeed && newValue <= trackedOverCapTau) {
                            result.overExactSeed[key] = newValue;
                        }
                        result.stats.evictedUpdates++;
                    }
                }
                return true;
            });
    };

    if (sortLeavesDesc) {
        for (const daf::Size leafId : leafOrder) {
            processLeaf(leafId);
        }
    } else {
        for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
            if (!leafAlive[leafId]) continue;
            processLeaf(leafId);
        }
    }

    const size_t exactSize = useCappedLowState ? cappedWindow.size() : exactWindow.size();
    result.lowSupport.reserve(exactSize / 2 + 1);
    result.bufferedSupport.reserve(exactSize / 2 + 1);
    result.overExactSupport.reserve(exactSize / 2 + 1);
    if (useCappedLowState) {
        for (auto &entry : cappedWindow) {
            const double value = static_cast<double>(entry.second);
            if (value <= activeCapTau) {
                result.lowSupport.emplace(entry.first, value);
            } else if (value <= windowCapTau) {
                result.bufferedSupport.emplace(entry.first, value);
            } else {
                result.overExactSupport.emplace(entry.first, value);
            }
        }
    } else {
        for (auto &entry : exactWindow) {
            if (entry.second <= activeCapTau) {
                result.lowSupport.emplace(entry.first, entry.second);
            } else if (entry.second <= windowCapTau) {
                result.bufferedSupport.emplace(entry.first, entry.second);
            } else {
                result.overExactSupport.emplace(entry.first, entry.second);
            }
        }
    }
    return result;
}

static int supportBucketIndex(double value, double capTau) {
    if (!(value > 0.0) || value > capTau + 1e-9) return -1;
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1e-9) return -1;
    const int bucket = static_cast<int>(rounded);
    if (bucket <= 0) return -1;
    return bucket;
}

static BucketedBandedLowState makeBucketedBandedLowState(
    BandedLowSupportBuildResult &&build, double activeCapTau, double windowCapTau,
    double exactOverCapTau = -1.0) {
    BucketedBandedLowState state;
    if (exactOverCapTau < windowCapTau) exactOverCapTau = windowCapTau;
    const int activeCap = std::max(1, static_cast<int>(std::ceil(activeCapTau)));
    const int windowCap = std::max(activeCap, static_cast<int>(std::ceil(windowCapTau)));
    const int overCap = std::max(windowCap, static_cast<int>(std::ceil(exactOverCapTau)));
    state.lowBuckets.resize(static_cast<size_t>(activeCap + 1));
    state.bufferedBuckets.resize(static_cast<size_t>(windowCap + 1));
    state.lowValue.reserve(build.lowSupport.size() * 2 + 1);
    state.bufferedValue.reserve(build.bufferedSupport.size() * 2 + 1);
    state.overValue.reserve(build.overExactSupport.size() * 2 + 1);
    state.overLowerBound = std::move(build.overLowerBound);

    for (auto &[key, value] : build.lowSupport) {
        const int bucket = supportBucketIndex(value, activeCapTau);
        if (bucket < 0 || bucket > activeCap) continue;
        state.lowValue.emplace(key, static_cast<uint16_t>(bucket));
        state.lowBuckets[static_cast<size_t>(bucket)].insert(key);
    }
    for (auto &[key, value] : build.bufferedSupport) {
        const int bucket = supportBucketIndex(value, windowCapTau);
        if (bucket < 0 || bucket > windowCap) continue;
        state.bufferedValue.emplace(key, static_cast<uint16_t>(bucket));
        state.bufferedBuckets[static_cast<size_t>(bucket)].insert(key);
    }
    for (auto &[key, value] : build.overExactSupport) {
        const int bucket = supportBucketIndex(value, exactOverCapTau);
        if (bucket <= windowCap || bucket > overCap) continue;
        state.overValue.emplace(key, static_cast<uint16_t>(bucket));
    }
    state.overWindow = std::move(build.overWindow);
    return state;
}

static int firstNonEmptyBucket(
    const std::vector<robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash>> &buckets) {
    for (size_t i = 1; i < buckets.size(); ++i) {
        if (!buckets[i].empty()) return static_cast<int>(i);
    }
    return -1;
}

static int chooseLookaheadFromOverSupportCache(
    const std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> &overSupportCache,
    double nextTau, double tauMax, uint64_t targetKeys) {
    const int nextBucket = std::max(1, static_cast<int>(std::ceil(nextTau - 1e-9)));
    const int maxBucket = std::max(nextBucket, static_cast<int>(std::floor(tauMax + 1e-9)));
    if (nextBucket > maxBucket || overSupportCache.empty()) return 0;

    std::vector<uint64_t> bucketCount(static_cast<size_t>(maxBucket + 1), 0);
    for (const auto &[key, support] : overSupportCache) {
        (void)key;
        const int bucket = supportBucketIndex(support, tauMax);
        if (bucket >= nextBucket && bucket <= maxBucket) {
            bucketCount[static_cast<size_t>(bucket)]++;
        }
    }

    uint64_t accum = 0;
    for (int bucket = nextBucket; bucket <= maxBucket; ++bucket) {
        accum += bucketCount[static_cast<size_t>(bucket)];
        if (accum >= targetKeys) {
            return bucket - nextBucket;
        }
    }
    return maxBucket - nextBucket;
}

static void insertBucketedKey(
    const SparseCliqueKey &key, int bucket,
    std::unordered_map<SparseCliqueKey, uint16_t, SparseCliqueKeyHash> &valueMap,
    std::vector<robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash>> &buckets);

static bool rebandBucketedStateFromExactCache(
    const robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> &persistentOverWindow,
    const std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> &overSupportCache,
    double activeCapTau, double windowCapTau, double exactOverCapTau,
    BucketedBandedLowState &newState) {
    if (overSupportCache.empty() && persistentOverWindow.empty()) return false;
    const int activeCap = std::max(1, static_cast<int>(std::ceil(activeCapTau)));
    const int windowCap = std::max(activeCap, static_cast<int>(std::ceil(windowCapTau)));
    newState = BucketedBandedLowState{};
    newState.lowBuckets.resize(static_cast<size_t>(activeCap + 1));
    newState.bufferedBuckets.resize(static_cast<size_t>(windowCap + 1));
    newState.lowValue.reserve(overSupportCache.size() / 8 + 1);
    newState.bufferedValue.reserve(overSupportCache.size() / 4 + 1);
    newState.overValue.reserve(overSupportCache.size() / 2 + 1);
    newState.overWindow.reserve(persistentOverWindow.size() + overSupportCache.size() * 2 + 1);
    for (const auto &key : persistentOverWindow) {
        newState.overWindow.insert(key);
    }

    auto placeExactSupport =
        [&](const SparseCliqueKey &key, double support,
            std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> *nextCache) -> bool {
            if (!(support > 0.0)) return true;
            if (support > exactOverCapTau) {
                newState.overWindow.insert(key);
                if (nextCache) {
                    (*nextCache)[key] = support;
                }
                return true;
            }
            if (support > windowCapTau) {
                const int bucket = supportBucketIndex(support, exactOverCapTau);
                if (bucket <= 0) return false;
                newState.overValue.emplace(key, static_cast<uint16_t>(bucket));
                return true;
            }
            const int bucket = supportBucketIndex(support, windowCapTau);
            if (bucket <= 0) return false;
            if (support <= activeCapTau) {
                insertBucketedKey(key, bucket, newState.lowValue, newState.lowBuckets);
            } else {
                insertBucketedKey(key, bucket, newState.bufferedValue, newState.bufferedBuckets);
            }
            return true;
        };

    for (const auto &[key, support] : overSupportCache) {
        if (!placeExactSupport(key, support, nullptr)) {
            return false;
        }
    }
    return true;
}

static void eraseBucketedKey(
    const SparseCliqueKey &key,
    std::unordered_map<SparseCliqueKey, uint16_t, SparseCliqueKeyHash> &valueMap,
    std::vector<robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash>> &buckets) {
    auto it = valueMap.find(key);
    if (it == valueMap.end()) return;
    const auto bucket = static_cast<size_t>(it->second);
    if (bucket < buckets.size()) buckets[bucket].erase(key);
    valueMap.erase(it);
}

static void insertBucketedKey(
    const SparseCliqueKey &key, int bucket,
    std::unordered_map<SparseCliqueKey, uint16_t, SparseCliqueKeyHash> &valueMap,
    std::vector<robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash>> &buckets) {
    if (bucket <= 0 || static_cast<size_t>(bucket) >= buckets.size()) return;
    valueMap[key] = static_cast<uint16_t>(bucket);
    buckets[static_cast<size_t>(bucket)].insert(key);
}

static void maybeRunLowSupportPrototype(
    const DynamicGraph<TreeGraphNode> &tree, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_LOW_SUPPORT")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "============= Quotient Low-Support Prototype =============" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip low-support prototype" << std::endl;
        std::cout << "  Reason: raw clique layer is too large for this prototype" << std::endl;
        std::cout << "==========================================================" << std::endl;
        return;
    }

    double tau = 1.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_LOW_SUPPORT_TAU")) {
        tau = std::max(0.0, std::atof(env));
    }

    LowSupportBuildResult build;

    daf::timeCount("low-support build (Quotient frontier)", [&]() {
        build = buildLowSupportLayer(tree, r, s, tau, rawCliqueCount);
    });

    bool compared = false;
    bool exactAgree = false;
    uint64_t exactFrontier = 0;
    uint64_t exactLowSupport = 0;
    double exactMinCore = std::numeric_limits<double>::quiet_NaN();
    if (std::getenv("PIVOTER_QUOTIENT_LOW_SUPPORT_COMPARE")) {
        std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> exactSupport;
        exactSupport.reserve(static_cast<size_t>(std::min<uint64_t>(rawCliqueCount, 10000000ULL)));
        for (const auto &leaf : tree.adj_list) {
            accumulateLeafSparseContribution(leaf, r, s, exactSupport, 1.0);
        }
        if (!exactSupport.empty()) {
            compared = true;
            exactMinCore = std::numeric_limits<double>::max();
            for (const auto &[key, value] : exactSupport) {
                (void)key;
                exactMinCore = std::min(exactMinCore, value);
                if (value <= tau) exactLowSupport++;
            }
            for (const auto &[key, value] : exactSupport) {
                if (value <= exactMinCore) exactFrontier++;
            }
            exactAgree = true;
            if (build.lowSupport.size() != exactLowSupport) {
                exactAgree = false;
            } else {
                for (const auto &[key, value] : exactSupport) {
                    if (value <= tau) {
                        auto it = build.lowSupport.find(key);
                        if (it == build.lowSupport.end() || std::abs(it->second - value) > 1e-9) {
                            exactAgree = false;
                            break;
                        }
                    }
                }
            }
        }
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "============= Quotient Low-Support Prototype =============" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Threshold tau:            " << tau << std::endl;
    std::cout << "  Low-support keys kept:    " << build.lowSupport.size() << std::endl;
    std::cout << "  Over-tau keys tracked:    " << build.overTau.size() << std::endl;
    std::cout << "  Kept updates:             " << build.stats.keptUpdates << std::endl;
    std::cout << "  Evicted updates:          " << build.stats.evictedUpdates << std::endl;
    std::cout << "  Skipped updates:          " << build.stats.skippedUpdates << std::endl;
    if (compared) {
        std::cout << "  Exact min core:           " << exactMinCore << std::endl;
        std::cout << "  Exact frontier size:      " << exactFrontier << std::endl;
        std::cout << "  Exact low-support keys:   " << exactLowSupport << std::endl;
        std::cout << "  Low-support exact agree:  " << (exactAgree ? "YES" : "NO") << std::endl;
    }
    std::cout << "  Prototype time:           " << elapsedMs << " ms" << std::endl;
    std::cout << "==========================================================" << std::endl;
}

static void maybeRunFirstFrontierLowSupportPrototype(
    const DynamicGraph<TreeGraphNode> &tree,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_FIRST_FRONTIER_LOW")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "========== Quotient First-Frontier Low-Support ==========" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip first-frontier low-support prototype" << std::endl;
        std::cout << "  Reason: raw clique layer is too large for this prototype" << std::endl;
        std::cout << "=========================================================" << std::endl;
        return;
    }

    LowSupportBuildResult build;
    daf::timeCount("first-frontier low-support build", [&]() {
        build = buildLowSupportLayer(tree, r, s, 1.0, rawCliqueCount);
    });
    if (build.lowSupport.empty()) return;

    std::vector<SparseCliqueKey> currentRemoveKeys;
    currentRemoveKeys.reserve(build.lowSupport.size());
    for (const auto &[key, value] : build.lowSupport) {
        if (value <= 1.0) currentRemoveKeys.push_back(key);
    }
    if (currentRemoveKeys.empty()) return;

    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
    std::vector<daf::Size> changedLeaf;
    removedKeyForLeaf.reserve(tree.adj_list.size() / 16 + 1);
    changedLeaf.reserve(tree.adj_list.size() / 16 + 1);

    daf::StaticVector<daf::Size> rmVerts;
    rmVerts.resize(r);
    for (const auto &rmKey : currentRemoveKeys) {
        for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
        daf::intersect_dense_sets_multi(rmVerts, treeGraphV.adj_list,
            [&](const TreeGraphNode &uClique) {
                auto &leafIdx = changedLeafIndex[uClique.v];
                if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                    leafIdx = removedKeyForLeaf.size();
                    removedKeyForLeaf.emplace_back();
                    changedLeaf.push_back(uClique.v);
                }
                removedKeyForLeaf[leafIdx].push_back(rmKey);
            });
    }
    rmVerts.free();

    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> sparseDelta;
    sparseDelta.reserve(1 << 20);
    uint64_t oldEntries = 0;
    uint64_t newEntries = 0;
    uint64_t generatedSubleaves = 0;
    uint64_t exactLeaves = 0;
    uint64_t singleRemovedLeaves = 0;
    uint64_t multiRemovedLeaves = 0;

    for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
        const daf::Size leafId = changedLeaf[idx];
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.size() < r) continue;
        const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
        if (removed.empty()) continue;
        if (removed.size() == 1) singleRemovedLeaves++;
        else multiRemovedLeaves++;
        exactLeaves++;
        accumulateExactLeafSparseDeltaBK(
            leaf, removed, r, s, sparseDelta, oldEntries, newEntries, generatedSubleaves);
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "========== Quotient First-Frontier Low-Support ==========" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Frontier keys:             " << currentRemoveKeys.size() << std::endl;
    std::cout << "  Changed leaves:            " << changedLeaf.size() << std::endl;
    std::cout << "  Exact-handled leaves:      " << exactLeaves << std::endl;
    std::cout << "  Single-removed leaves:     " << singleRemovedLeaves << std::endl;
    std::cout << "  Multi-removed leaves:      " << multiRemovedLeaves << std::endl;
    std::cout << "  Generated subleaves:       " << generatedSubleaves << std::endl;
    std::cout << "  Old local entries:         " << oldEntries << std::endl;
    std::cout << "  New local entries:         " << newEntries << std::endl;
    std::cout << "  Unique sparse deltas:      " << sparseDelta.size() << std::endl;
    std::cout << "  Prototype time:            " << elapsedMs << " ms" << std::endl;
    std::cout << "=========================================================" << std::endl;
}

static void maybeRunFirstFrontierNextPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_FIRST_FRONTIER_NEXT")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "========= Quotient First-Frontier Next Prototype =========" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip first-frontier next prototype" << std::endl;
        std::cout << "  Reason: raw clique layer is too large for this prototype" << std::endl;
        std::cout << "==========================================================" << std::endl;
        return;
    }

    double tau = 2.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_FIRST_FRONTIER_NEXT_TAU")) {
        tau = std::max(1.0, std::atof(env));
    }

    LowSupportBuildResult build;
    daf::timeCount("first-frontier next low-support build", [&]() {
        build = buildLowSupportLayer(tree, r, s, tau, rawCliqueCount);
    });
    if (build.lowSupport.empty()) return;

    std::vector<SparseCliqueKey> frontierKeys;
    frontierKeys.reserve(build.lowSupport.size());
    for (const auto &[key, value] : build.lowSupport) {
        if (value <= 1.0) frontierKeys.push_back(key);
    }
    if (frontierKeys.empty()) return;

    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
    std::vector<daf::Size> changedLeaf;
    removedKeyForLeaf.reserve(tree.adj_list.size() / 16 + 1);
    changedLeaf.reserve(tree.adj_list.size() / 16 + 1);

    daf::StaticVector<daf::Size> rmVerts;
    rmVerts.resize(r);
    for (const auto &rmKey : frontierKeys) {
        for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
        daf::intersect_dense_sets_multi(rmVerts, treeGraphV.adj_list,
            [&](const TreeGraphNode &uClique) {
                auto &leafIdx = changedLeafIndex[uClique.v];
                if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                    leafIdx = removedKeyForLeaf.size();
                    removedKeyForLeaf.emplace_back();
                    changedLeaf.push_back(uClique.v);
                }
                removedKeyForLeaf[leafIdx].push_back(rmKey);
            });
    }
    rmVerts.free();

    std::vector<std::vector<TreeGraphNode>> activeLeaves;
    activeLeaves.reserve(tree.adj_list.size() * 2 + 16);
    for (const auto &leaf : tree.adj_list) {
        activeLeaves.emplace_back(leaf.begin(), leaf.end());
    }
    std::vector<uint8_t> leafAlive(activeLeaves.size(), 1);
    std::vector<robin_hood::unordered_flat_set<daf::Size>> activeLeafByVertex(edgeGraph.n);
    std::vector<robin_hood::unordered_flat_set<daf::Size>> activePivotLeafByVertex(edgeGraph.n);
    std::vector<uint16_t> activeKeepCount(activeLeaves.size(), 0);
    std::vector<uint16_t> activePivotCount(activeLeaves.size(), 0);
    for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
        addLeafToVertexAndPivotIndex(activeLeaves[leafId], leafId,
            activeLeafByVertex, activePivotLeafByVertex, activeKeepCount, activePivotCount, r);
    }

    robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> frontierSet;
    frontierSet.reserve(frontierKeys.size() * 2 + 1);
    for (const auto &key : frontierKeys) frontierSet.insert(key);

    robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> candidateKeys;
    candidateKeys.reserve(frontierKeys.size() * 8 + 1024);
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> dummyDelta;
    dummyDelta.reserve(1 << 20);
    uint64_t oldEntries = 0;
    uint64_t newEntries = 0;
    uint64_t generatedSubleaves = 0;
    uint64_t exactLeaves = 0;
    uint64_t singleRemovedLeaves = 0;
    uint64_t multiRemovedLeaves = 0;

    for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
        const daf::Size leafId = changedLeaf[idx];
        auto &leaf = activeLeaves[leafId];
        if (!leafAlive[leafId] || leaf.size() < r) continue;
        const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
        if (removed.empty()) continue;
        if (removed.size() == 1) singleRemovedLeaves++;
        else multiRemovedLeaves++;
        exactLeaves++;

        collectLeafSparseKeys(leaf, r, s, candidateKeys);
        auto newLeaves = splitLeafByRemovedKeys(
            leaf, removed, r, s, dummyDelta, &oldEntries, &newEntries);

        removeLeafFromVertexIndex(leaf, leafId, activeLeafByVertex, r);
        leafAlive[leafId] = 0;
        generatedSubleaves += newLeaves.size();

        for (auto &newLeaf : newLeaves) {
            collectLeafSparseKeys(newLeaf, r, s, candidateKeys);
            daf::Size newLeafId = activeLeaves.size();
            activeLeaves.push_back(std::move(newLeaf));
            leafAlive.push_back(1);
            addLeafToVertexIndex(activeLeaves.back(), newLeafId, activeLeafByVertex, r);
        }
    }

    for (const auto &key : frontierKeys) {
        candidateKeys.erase(key);
    }

    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> recomputedTouched;
    recomputedTouched.reserve(candidateKeys.size() / 2 + 1);
    double touchedMin = std::numeric_limits<double>::infinity();
    for (const auto &key : candidateKeys) {
        daf::StaticVector<daf::Size> verts;
        verts.resize(r);
        for (daf::CliqueSize j = 0; j < r; ++j) verts[j] = key.verts[j];
        double support = 0.0;
        daf::intersect_dense_sets_multi(verts, activeLeafByVertex,
            [&](const daf::Size &leafId) {
                if (leafId < activeLeaves.size() && leafAlive[leafId]) {
                    support += leafContributionForKey(activeLeaves[leafId], key, r, s);
                }
            });
        verts.free();
        if (support > 0.0) {
            recomputedTouched.emplace(key, support);
            touchedMin = std::min(touchedMin, support);
        }
    }

    double untouchedMin = std::numeric_limits<double>::infinity();
    uint64_t untouchedLowKeys = 0;
    for (const auto &[key, value] : build.lowSupport) {
        if (value <= 1.0) continue;
        if (candidateKeys.find(key) != candidateKeys.end()) continue;
        untouchedLowKeys++;
        untouchedMin = std::min(untouchedMin, value);
    }

    double nextMin = std::min(untouchedMin, touchedMin);
    uint64_t nextFrontier = 0;
    if (std::isfinite(nextMin)) {
        for (const auto &[key, value] : build.lowSupport) {
            if (value <= 1.0) continue;
            if (candidateKeys.find(key) != candidateKeys.end()) continue;
            if (std::abs(value - nextMin) <= 1e-9) nextFrontier++;
        }
        for (const auto &[key, value] : recomputedTouched) {
            (void)key;
            if (std::abs(value - nextMin) <= 1e-9) nextFrontier++;
        }
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "========= Quotient First-Frontier Next Prototype =========" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Threshold tau:            " << tau << std::endl;
    std::cout << "  Frontier keys:            " << frontierKeys.size() << std::endl;
    std::cout << "  Changed leaves:           " << changedLeaf.size() << std::endl;
    std::cout << "  Exact-handled leaves:     " << exactLeaves << std::endl;
    std::cout << "  Single-removed leaves:    " << singleRemovedLeaves << std::endl;
    std::cout << "  Multi-removed leaves:     " << multiRemovedLeaves << std::endl;
    std::cout << "  Generated subleaves:      " << generatedSubleaves << std::endl;
    std::cout << "  Candidate touched keys:   " << candidateKeys.size() << std::endl;
    std::cout << "  Recomputed touched keys:  " << recomputedTouched.size() << std::endl;
    std::cout << "  Untouched low keys:       " << untouchedLowKeys << std::endl;
    std::cout << "  Old local entries:        " << oldEntries << std::endl;
    std::cout << "  New local entries:        " << newEntries << std::endl;
    std::cout << "  Touched min:              "
              << (std::isfinite(touchedMin) ? touchedMin : -1.0) << std::endl;
    std::cout << "  Untouched low min:        "
              << (std::isfinite(untouchedMin) ? untouchedMin : -1.0) << std::endl;
    std::cout << "  Next frontier min:        "
              << (std::isfinite(nextMin) ? nextMin : -1.0) << std::endl;
    std::cout << "  Next frontier size:       " << nextFrontier << std::endl;
    std::cout << "  Prototype time:           " << elapsedMs << " ms" << std::endl;
    std::cout << "==========================================================" << std::endl;
}

static void maybeRunMultiRoundLowSupportPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_LOW")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "========== Quotient Multi-Round Low-Support ==========" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip multi-round low-support prototype" << std::endl;
        std::cout << "  Reason: raw clique layer is too large for this prototype" << std::endl;
        std::cout << "=======================================================" << std::endl;
        return;
    }

    double tau = 2.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_LOW_TAU")) {
        tau = std::max(1.0, std::atof(env));
    }
    int maxRounds = 5;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_MAX")) {
        maxRounds = std::max(1, std::atoi(env));
    }

    LowSupportBuildResult build;
    daf::timeCount("multi-round low-support build", [&]() {
        build = buildLowSupportLayer(tree, r, s, tau, rawCliqueCount);
    });
    if (build.lowSupport.empty()) return;

    std::vector<std::vector<TreeGraphNode>> activeLeaves;
    activeLeaves.reserve(tree.adj_list.size() * 2 + 16);
    for (const auto &leaf : tree.adj_list) {
        activeLeaves.emplace_back(leaf.begin(), leaf.end());
    }
    std::vector<uint8_t> leafAlive(activeLeaves.size(), 1);
    std::vector<robin_hood::unordered_flat_set<daf::Size>> activeLeafByVertex(edgeGraph.n);
    for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
        addLeafToVertexIndex(activeLeaves[leafId], leafId, activeLeafByVertex, r);
    }

    auto currentLow = std::move(build.lowSupport);
    int completedRounds = 0;
    uint64_t totalGeneratedSubleaves = 0;
    uint64_t totalExactLeaves = 0;
    uint64_t totalCandidateKeys = 0;
    double maxObservedMin = 0.0;

    for (int round = 1; round <= maxRounds; ++round) {
        if (currentLow.empty()) break;

        double roundMin = std::numeric_limits<double>::max();
        for (const auto &[key, value] : currentLow) {
            (void)key;
            roundMin = std::min(roundMin, value);
        }
        if (!std::isfinite(roundMin) || roundMin > tau) break;
        maxObservedMin = std::max(maxObservedMin, roundMin);

        std::vector<SparseCliqueKey> frontierKeys;
        frontierKeys.reserve(currentLow.size() / 16 + 1);
        for (const auto &[key, value] : currentLow) {
            if (std::abs(value - roundMin) <= 1e-9) frontierKeys.push_back(key);
        }
        if (frontierKeys.empty()) break;

        std::vector<daf::Size> changedLeafIndex(activeLeaves.size(), std::numeric_limits<daf::Size>::max());
        std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
        std::vector<daf::Size> changedLeaf;
        removedKeyForLeaf.reserve(activeLeaves.size() / 16 + 1);
        changedLeaf.reserve(activeLeaves.size() / 16 + 1);

        daf::StaticVector<daf::Size> rmVerts;
        rmVerts.resize(r);
        for (const auto &rmKey : frontierKeys) {
            for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
            daf::intersect_dense_sets_multi(rmVerts, activeLeafByVertex,
                [&](const daf::Size &leafId) {
                    if (leafId >= activeLeaves.size() || !leafAlive[leafId]) return;
                    auto &leafIdx = changedLeafIndex[leafId];
                    if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                        leafIdx = removedKeyForLeaf.size();
                        removedKeyForLeaf.emplace_back();
                        changedLeaf.push_back(leafId);
                    }
                    removedKeyForLeaf[leafIdx].push_back(rmKey);
                });
        }
        rmVerts.free();

        robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> frontierSet;
        frontierSet.reserve(frontierKeys.size() * 2 + 1);
        for (const auto &key : frontierKeys) frontierSet.insert(key);

        robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> candidateKeys;
        candidateKeys.reserve(frontierKeys.size() * 8 + 1024);
        std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> dummyDelta;
        dummyDelta.reserve(1 << 20);
        uint64_t roundOldEntries = 0;
        uint64_t roundNewEntries = 0;
        uint64_t roundGeneratedSubleaves = 0;
        uint64_t roundExactLeaves = 0;
        uint64_t roundSingleRemovedLeaves = 0;
        uint64_t roundMultiRemovedLeaves = 0;

        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            const daf::Size leafId = changedLeaf[idx];
            auto &leaf = activeLeaves[leafId];
            if (!leafAlive[leafId] || leaf.size() < r) continue;
            const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
            if (removed.empty()) continue;
            if (removed.size() == 1) roundSingleRemovedLeaves++;
            else roundMultiRemovedLeaves++;
            roundExactLeaves++;

            collectLeafSparseKeys(leaf, r, s, candidateKeys);
            auto newLeaves = splitLeafByRemovedKeys(
                leaf, removed, r, s, dummyDelta, &roundOldEntries, &roundNewEntries);

            removeLeafFromVertexIndex(leaf, leafId, activeLeafByVertex, r);
            leafAlive[leafId] = 0;
            roundGeneratedSubleaves += newLeaves.size();
            totalGeneratedSubleaves += newLeaves.size();

            for (auto &newLeaf : newLeaves) {
                collectLeafSparseKeys(newLeaf, r, s, candidateKeys);
                daf::Size newLeafId = activeLeaves.size();
                activeLeaves.push_back(std::move(newLeaf));
                leafAlive.push_back(1);
                addLeafToVertexIndex(activeLeaves.back(), newLeafId, activeLeafByVertex, r);
            }
        }

        for (const auto &key : frontierKeys) {
            candidateKeys.erase(key);
        }

        std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> recomputedTouched;
        recomputedTouched.reserve(candidateKeys.size() / 2 + 1);
        double touchedMin = std::numeric_limits<double>::infinity();
        for (const auto &key : candidateKeys) {
            daf::StaticVector<daf::Size> verts;
            verts.resize(r);
            for (daf::CliqueSize j = 0; j < r; ++j) verts[j] = key.verts[j];
            double support = 0.0;
            daf::intersect_dense_sets_multi(verts, activeLeafByVertex,
                [&](const daf::Size &leafId) {
                    if (leafId < activeLeaves.size() && leafAlive[leafId]) {
                        support += leafContributionForKey(activeLeaves[leafId], key, r, s);
                    }
                });
            verts.free();
            if (support > 0.0) {
                recomputedTouched.emplace(key, support);
                touchedMin = std::min(touchedMin, support);
            }
        }

        std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> nextLow;
        nextLow.reserve(currentLow.size() + recomputedTouched.size());
        double untouchedMin = std::numeric_limits<double>::infinity();
        uint64_t untouchedLowKeys = 0;
        for (const auto &[key, value] : currentLow) {
            if (frontierSet.find(key) != frontierSet.end()) continue;
            if (candidateKeys.find(key) != candidateKeys.end()) continue;
            nextLow.emplace(key, value);
            untouchedLowKeys++;
            untouchedMin = std::min(untouchedMin, value);
        }
        for (const auto &[key, value] : recomputedTouched) {
            if (value <= tau) nextLow[key] = value;
        }

        double nextMin = std::numeric_limits<double>::infinity();
        for (const auto &[key, value] : nextLow) {
            (void)key;
            nextMin = std::min(nextMin, value);
        }
        uint64_t nextFrontier = 0;
        if (std::isfinite(nextMin)) {
            for (const auto &[key, value] : nextLow) {
                (void)key;
                if (std::abs(value - nextMin) <= 1e-9) nextFrontier++;
            }
        }

        totalExactLeaves += roundExactLeaves;
        totalCandidateKeys += candidateKeys.size();
        completedRounds = round;

        std::cout << "---- Quotient low-support round " << round << " ----" << std::endl;
        std::cout << "  Round min:                " << roundMin << std::endl;
        std::cout << "  Frontier keys:            " << frontierKeys.size() << std::endl;
        std::cout << "  Changed leaves:           " << changedLeaf.size() << std::endl;
        std::cout << "  Exact-handled leaves:     " << roundExactLeaves << std::endl;
        std::cout << "  Single-removed leaves:    " << roundSingleRemovedLeaves << std::endl;
        std::cout << "  Multi-removed leaves:     " << roundMultiRemovedLeaves << std::endl;
        std::cout << "  Generated subleaves:      " << roundGeneratedSubleaves << std::endl;
        std::cout << "  Candidate touched keys:   " << candidateKeys.size() << std::endl;
        std::cout << "  Recomputed touched keys:  " << recomputedTouched.size() << std::endl;
        std::cout << "  Untouched low keys:       " << untouchedLowKeys << std::endl;
        std::cout << "  Old local entries:        " << roundOldEntries << std::endl;
        std::cout << "  New local entries:        " << roundNewEntries << std::endl;
        std::cout << "  Touched min:              "
                  << (std::isfinite(touchedMin) ? touchedMin : -1.0) << std::endl;
        std::cout << "  Untouched low min:        "
                  << (std::isfinite(untouchedMin) ? untouchedMin : -1.0) << std::endl;
        std::cout << "  Next low keys:            " << nextLow.size() << std::endl;
        std::cout << "  Next frontier min:        "
                  << (std::isfinite(nextMin) ? nextMin : -1.0) << std::endl;
        std::cout << "  Next frontier size:       " << nextFrontier << std::endl;

        currentLow.swap(nextLow);
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();
    std::cout << "========== Quotient Multi-Round Low-Support ==========" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Threshold tau:            " << tau << std::endl;
    std::cout << "  Completed rounds:         " << completedRounds << std::endl;
    std::cout << "  Total exact leaves:       " << totalExactLeaves << std::endl;
    std::cout << "  Total candidate keys:     " << totalCandidateKeys << std::endl;
    std::cout << "  Total spawned subleaves:  " << totalGeneratedSubleaves << std::endl;
    std::cout << "  Max observed round min:   " << maxObservedMin << std::endl;
    std::cout << "  Remaining low keys:       " << currentLow.size() << std::endl;
    std::cout << "  Prototype time:           " << elapsedMs << " ms" << std::endl;
    std::cout << "=======================================================" << std::endl;
}

static void maybeRunAdaptiveLowSupportPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_ADAPTIVE_LOW")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "========== Quotient Adaptive Low-Support ==========" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip adaptive low-support prototype" << std::endl;
        std::cout << "  Reason: raw clique layer is too large for this prototype" << std::endl;
        std::cout << "===================================================" << std::endl;
        return;
    }

    double tauStart = 2.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ADAPTIVE_LOW_TAU_START")) {
        tauStart = std::max(1.0, std::atof(env));
    }
    double tauMax = 8.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ADAPTIVE_LOW_TAU_MAX")) {
        tauMax = std::max(tauStart, std::atof(env));
    }
    int maxRounds = 200;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_MAX")) {
        maxRounds = std::max(1, std::atoi(env));
    }
    std::vector<std::vector<TreeGraphNode>> activeLeaves;
    activeLeaves.reserve(tree.adj_list.size() * 2 + 16);
    for (const auto &leaf : tree.adj_list) {
        activeLeaves.emplace_back(leaf.begin(), leaf.end());
    }
    std::vector<uint8_t> leafAlive(activeLeaves.size(), 1);
    std::vector<robin_hood::unordered_flat_set<daf::Size>> activeLeafByVertex(edgeGraph.n);
    for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
        addLeafToVertexIndex(activeLeaves[leafId], leafId, activeLeafByVertex, r);
    }

    int completedRounds = 0;
    int tauPhaseCount = 0;
    uint64_t totalExactLeaves = 0;
    uint64_t totalCandidateKeys = 0;
    uint64_t totalSpawnedSubleaves = 0;
    double maxObservedMin = 0.0;
    double lastTau = tauStart;
    uint64_t totalRebuildMs = 0;
    std::vector<AdaptiveLowPhaseStats> phaseStats;

    for (double tau = tauStart; tau <= tauMax && completedRounds < maxRounds; tau += 1.0) {
        lastTau = tau;
        tauPhaseCount++;
        AdaptiveLowPhaseStats phase;
        phase.tau = tau;
        phase.capTau = tau;
        LowSupportBuildResult build;
        auto rebuildStart = std::chrono::high_resolution_clock::now();
        daf::timeCount("adaptive low-support rebuild", [&]() {
            build = buildLowSupportLayerActive(activeLeaves, leafAlive, r, s, tau, rawCliqueCount);
        });
        phase.rebuildMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - rebuildStart).count());
        totalRebuildMs += phase.rebuildMs;
        phase.initialLowKeys = build.lowSupport.size();
        phase.initialOverTauKeys = build.overTau.size();
        phase.buildStats = build.stats;
        auto currentLow = std::move(build.lowSupport);
        if (currentLow.empty()) {
            phaseStats.push_back(phase);
            continue;
        }

        while (!currentLow.empty() && completedRounds < maxRounds) {
            double roundMin = std::numeric_limits<double>::max();
            for (const auto &[key, value] : currentLow) {
                (void)key;
                roundMin = std::min(roundMin, value);
            }
            if (!std::isfinite(roundMin) || roundMin > tau) break;
            maxObservedMin = std::max(maxObservedMin, roundMin);
            phase.maxRoundMin = std::max(phase.maxRoundMin, roundMin);

            std::vector<SparseCliqueKey> frontierKeys;
            frontierKeys.reserve(currentLow.size() / 16 + 1);
            for (const auto &[key, value] : currentLow) {
                if (std::abs(value - roundMin) <= 1e-9) frontierKeys.push_back(key);
            }
            if (frontierKeys.empty()) break;
            phase.frontierKeys += frontierKeys.size();

            std::vector<daf::Size> changedLeafIndex(activeLeaves.size(), std::numeric_limits<daf::Size>::max());
            std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
            std::vector<daf::Size> changedLeaf;
            removedKeyForLeaf.reserve(activeLeaves.size() / 16 + 1);
            changedLeaf.reserve(activeLeaves.size() / 16 + 1);

            daf::StaticVector<daf::Size> rmVerts;
            rmVerts.resize(r);
            for (const auto &rmKey : frontierKeys) {
                for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
                daf::intersect_dense_sets_multi(rmVerts, activeLeafByVertex,
                    [&](const daf::Size &leafId) {
                        if (leafId >= activeLeaves.size() || !leafAlive[leafId]) return;
                        auto &leafIdx = changedLeafIndex[leafId];
                        if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                            leafIdx = removedKeyForLeaf.size();
                            removedKeyForLeaf.emplace_back();
                            changedLeaf.push_back(leafId);
                        }
                        removedKeyForLeaf[leafIdx].push_back(rmKey);
                    });
            }
            rmVerts.free();

            robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> frontierSet;
            frontierSet.reserve(frontierKeys.size() * 2 + 1);
            for (const auto &key : frontierKeys) frontierSet.insert(key);

            robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> candidateKeys;
            candidateKeys.reserve(frontierKeys.size() * 8 + 1024);
            std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> dummyDelta;
            dummyDelta.reserve(1 << 20);
            uint64_t roundOldEntries = 0;
            uint64_t roundNewEntries = 0;
            uint64_t roundGeneratedSubleaves = 0;
            uint64_t roundExactLeaves = 0;

            for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
                const daf::Size leafId = changedLeaf[idx];
                auto &leaf = activeLeaves[leafId];
                if (!leafAlive[leafId] || leaf.size() < r) continue;
                const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
                if (removed.empty()) continue;
                roundExactLeaves++;
                phase.exactLeaves++;

                collectLeafSparseKeys(leaf, r, s, candidateKeys);
                auto newLeaves = splitLeafByRemovedKeys(
                    leaf, removed, r, s, dummyDelta, &roundOldEntries, &roundNewEntries);

                removeLeafFromVertexIndex(leaf, leafId, activeLeafByVertex, r);
                leafAlive[leafId] = 0;
                roundGeneratedSubleaves += newLeaves.size();
                totalSpawnedSubleaves += newLeaves.size();
                phase.spawnedSubleaves += newLeaves.size();

                for (auto &newLeaf : newLeaves) {
                    collectLeafSparseKeys(newLeaf, r, s, candidateKeys);
                    daf::Size newLeafId = activeLeaves.size();
                    activeLeaves.push_back(std::move(newLeaf));
                    leafAlive.push_back(1);
                    addLeafToVertexIndex(activeLeaves.back(), newLeafId, activeLeafByVertex, r);
                }
            }

            for (const auto &key : frontierKeys) {
                candidateKeys.erase(key);
            }

            std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> recomputedTouched;
            recomputedTouched.reserve(candidateKeys.size() / 2 + 1);
            for (const auto &key : candidateKeys) {
                daf::StaticVector<daf::Size> verts;
                verts.resize(r);
                for (daf::CliqueSize j = 0; j < r; ++j) verts[j] = key.verts[j];
                double support = 0.0;
                daf::intersect_dense_sets_multi(verts, activeLeafByVertex,
                    [&](const daf::Size &leafId) {
                        if (leafId < activeLeaves.size() && leafAlive[leafId]) {
                            support += leafContributionForKey(activeLeaves[leafId], key, r, s);
                        }
                    });
                verts.free();
                if (support > 0.0) recomputedTouched.emplace(key, support);
            }

            std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> nextLow;
            nextLow.reserve(currentLow.size() + recomputedTouched.size());
            double nextMin = std::numeric_limits<double>::infinity();
            uint64_t nextFrontier = 0;
            for (const auto &[key, value] : currentLow) {
                if (frontierSet.find(key) != frontierSet.end()) continue;
                if (candidateKeys.find(key) != candidateKeys.end()) continue;
                nextLow.emplace(key, value);
                nextMin = std::min(nextMin, value);
            }
            for (const auto &[key, value] : recomputedTouched) {
                if (value <= tau) {
                    nextLow[key] = value;
                    nextMin = std::min(nextMin, value);
                }
            }
            if (std::isfinite(nextMin)) {
                for (const auto &[key, value] : nextLow) {
                    (void)key;
                    if (std::abs(value - nextMin) <= 1e-9) nextFrontier++;
                }
            }

            completedRounds++;
            phase.rounds++;
            totalExactLeaves += roundExactLeaves;
            totalCandidateKeys += candidateKeys.size();
            phase.candidateKeys += candidateKeys.size();

            std::cout << "---- Quotient adaptive-low round " << completedRounds << " ----" << std::endl;
            std::cout << "  Tau:                      " << tau << std::endl;
            std::cout << "  Round min:                " << roundMin << std::endl;
            std::cout << "  Frontier keys:            " << frontierKeys.size() << std::endl;
            std::cout << "  Changed leaves:           " << changedLeaf.size() << std::endl;
            std::cout << "  Exact-handled leaves:     " << roundExactLeaves << std::endl;
            std::cout << "  Generated subleaves:      " << roundGeneratedSubleaves << std::endl;
            std::cout << "  Candidate touched keys:   " << candidateKeys.size() << std::endl;
            std::cout << "  Recomputed touched keys:  " << recomputedTouched.size() << std::endl;
            std::cout << "  Old local entries:        " << roundOldEntries << std::endl;
            std::cout << "  New local entries:        " << roundNewEntries << std::endl;
            std::cout << "  Next low keys:            " << nextLow.size() << std::endl;
            std::cout << "  Next frontier min:        "
                      << (std::isfinite(nextMin) ? nextMin : -1.0) << std::endl;
            std::cout << "  Next frontier size:       " << nextFrontier << std::endl;

            currentLow.swap(nextLow);
        }

        phase.remainingLowKeys = currentLow.size();
        phaseStats.push_back(phase);
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();
    std::cout << "========== Quotient Adaptive Low-Support ==========" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Tau start/max:            " << tauStart << " / " << tauMax << std::endl;
    std::cout << "  Tau phases used:          " << tauPhaseCount << std::endl;
    std::cout << "  Completed rounds:         " << completedRounds << std::endl;
    std::cout << "  Total exact leaves:       " << totalExactLeaves << std::endl;
    std::cout << "  Total candidate keys:     " << totalCandidateKeys << std::endl;
    std::cout << "  Total spawned subleaves:  " << totalSpawnedSubleaves << std::endl;
    std::cout << "  Max observed round min:   " << maxObservedMin << std::endl;
    std::cout << "  Total rebuild time:       " << totalRebuildMs << " ms" << std::endl;
    std::cout << "  Final tau:                " << lastTau << std::endl;
    std::cout << "  Prototype time:           " << elapsedMs << " ms" << std::endl;
    for (const auto &phase : phaseStats) {
        std::cout << "  Phase tau=" << phase.tau
                  << " cap=" << phase.capTau
                  << " rebuild_ms=" << phase.rebuildMs
                  << " init_low=" << phase.initialLowKeys
                  << " init_over=" << phase.initialOverTauKeys
                  << " rounds=" << phase.rounds
                  << " frontier_keys=" << phase.frontierKeys
                  << " exact_leaves=" << phase.exactLeaves
                  << " cand_keys=" << phase.candidateKeys
                  << " spawned=" << phase.spawnedSubleaves
                  << " max_round_min=" << phase.maxRoundMin
                  << " remain_low=" << phase.remainingLowKeys
                  << " kept=" << phase.buildStats.keptUpdates
                  << " evicted=" << phase.buildStats.evictedUpdates
                  << " skipped=" << phase.buildStats.skippedUpdates
                  << std::endl;
    }
    std::cout << "===================================================" << std::endl;
}

static void maybeRunBufferedAdaptiveLowSupportPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "========== Quotient Buffered Adaptive Low-Support ==========" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip buffered adaptive prototype" << std::endl;
        std::cout << "  Reason: raw clique layer is too large for this prototype" << std::endl;
        std::cout << "============================================================" << std::endl;
        return;
    }

    double tauStart = 1.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_START")) {
        tauStart = std::max(1.0, std::atof(env));
    }
    double tauMax = 8.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_MAX")) {
        tauMax = std::max(tauStart, std::atof(env));
    }
    int overExactBand = 0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_OVER_BAND")) {
        overExactBand = std::max(0, std::atoi(env));
    }
    int lookahead = 1;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_LOOKAHEAD")) {
        lookahead = std::max(0, std::atoi(env));
    }
    const bool autoLookahead = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO") != nullptr;
    uint64_t autoThresholdKeys = 300000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_THRESHOLD_KEYS")) {
        autoThresholdKeys = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoRawCliqueThreshold = 1000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_RAW_THRESHOLD")) {
        autoRawCliqueThreshold = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoWideRawCliqueThreshold = 20000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_WIDE_RAW_THRESHOLD")) {
        autoWideRawCliqueThreshold = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoWideExplicitMin = 10000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_WIDE_EXPLICIT_MIN")) {
        autoWideExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoWideExplicitMax = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_WIDE_EXPLICIT_MAX")) {
        autoWideExplicitMax = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    int autoCompactMaxLeafSize = 20;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_MAX_LEAF")) {
        autoCompactMaxLeafSize = std::max(1, std::atoi(env));
    }
    uint64_t autoCompactExplicitMax = 20000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_EXPLICIT_MAX")) {
        autoCompactExplicitMax = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    double autoCompactCompressionMax = 12.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_COMPRESSION_MAX")) {
        autoCompactCompressionMax = std::max(1.0, std::atof(env));
    }
    int autoCompactLookahead = 4;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_LOOKAHEAD")) {
        autoCompactLookahead = std::max(1, std::atoi(env));
    }
    int autoCacheRebandMinLeafSize = 40;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_MIN_LEAF")) {
        autoCacheRebandMinLeafSize = std::max(1, std::atoi(env));
    }
    uint64_t autoCacheRebandExplicitMin = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_EXPLICIT_MIN")) {
        autoCacheRebandExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    double autoCacheRebandCompressionMin = 50.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_COMPRESSION_MIN")) {
        autoCacheRebandCompressionMin = std::max(1.0, std::atof(env));
    }
    int autoCacheRebandLookahead = 3;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_LOOKAHEAD")) {
        autoCacheRebandLookahead = std::max(1, std::atoi(env));
    }
    const bool enableAutoCacheRebandHeavy =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_HEAVY") != nullptr;
    int autoTwoPhaseMinLeafSize = 40;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE_MIN_LEAF")) {
        autoTwoPhaseMinLeafSize = std::max(1, std::atoi(env));
    }
    uint64_t autoTwoPhaseExplicitMin = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE_EXPLICIT_MIN")) {
        autoTwoPhaseExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool enableAutoTwoPhase =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE") != nullptr;
    const bool enableCacheReband =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND") != nullptr;
    double autoFullBandOverRatio = 8.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_FULL_OVER_RATIO")) {
        autoFullBandOverRatio = std::max(1.0, std::atof(env));
    }
    double autoFullBandRebuildFactor = 2.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_FULL_REBUILD_FACTOR")) {
        autoFullBandRebuildFactor = std::max(1.0, std::atof(env));
    }
    int autoFullBandMinRemainingTau = 2;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_FULL_MIN_REMAINING")) {
        autoFullBandMinRemainingTau = std::max(1, std::atoi(env));
    }
    int maxRounds = 200;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_MAX")) {
        maxRounds = std::max(1, std::atoi(env));
    }
    const bool enableDeltaFast = std::getenv("PIVOTER_QUOTIENT_BUCKETED_DELTA_FAST") != nullptr;

    std::vector<std::vector<TreeGraphNode>> activeLeaves;
    activeLeaves.reserve(tree.adj_list.size() * 2 + 16);
    for (const auto &leaf : tree.adj_list) {
        activeLeaves.emplace_back(leaf.begin(), leaf.end());
    }
    std::vector<uint8_t> leafAlive(activeLeaves.size(), 1);
    std::vector<robin_hood::unordered_flat_set<daf::Size>> activeLeafByVertex(edgeGraph.n);
    for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
        addLeafToVertexIndex(activeLeaves[leafId], leafId, activeLeafByVertex, r);
    }

    int completedRounds = 0;
    int rebuildCount = 0;
    uint64_t totalExactLeaves = 0;
    uint64_t totalCandidateKeys = 0;
    uint64_t totalSpawnedSubleaves = 0;
    double maxObservedMin = 0.0;
    double lastTau = tauStart;
    uint64_t totalRebuildMs = 0;
    std::vector<AdaptiveLowPhaseStats> phaseStats;

    double currentTau = tauStart;
    double currentCapTau = -1.0;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> currentLow;

    while (currentTau <= tauMax && completedRounds < maxRounds) {
        if (currentLow.empty()) {
            int phaseLookahead = lookahead;
            if (autoLookahead) {
                phaseLookahead = 1;
                if (phaseStats.empty()) {
                    if (static_cast<uint64_t>(rawCliqueCount) <= autoRawCliqueThreshold) {
                        phaseLookahead = lookahead;
                    }
                } else {
                    const auto &prev = phaseStats.back();
                    if (prev.initialLowKeys <= autoThresholdKeys) phaseLookahead = lookahead;
                }
            }
            currentCapTau = std::min(tauMax, currentTau + static_cast<double>(phaseLookahead));
            AdaptiveLowPhaseStats phase;
            phase.tau = currentTau;
            phase.capTau = currentCapTau;
            phase.lookahead = phaseLookahead;
            LowSupportBuildResult build;
            auto phaseStart = std::chrono::high_resolution_clock::now();
            auto rebuildStart = std::chrono::high_resolution_clock::now();
            daf::timeCount("buffered adaptive low-support rebuild", [&]() {
                build = buildLowSupportLayerActive(activeLeaves, leafAlive, r, s, currentCapTau, rawCliqueCount);
            });
            phase.rebuildMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - rebuildStart).count());
            totalRebuildMs += phase.rebuildMs;
            rebuildCount++;
            phase.initialLowKeys = build.lowSupport.size();
            phase.initialOverTauKeys = build.overTau.size();
            phase.buildStats = build.stats;
            currentLow = std::move(build.lowSupport);
            if (currentLow.empty()) {
                phase.phaseMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - phaseStart).count());
                phaseStats.push_back(phase);
                currentTau = currentCapTau + 1.0;
                lastTau = currentTau;
                continue;
            }

            while (!currentLow.empty() && completedRounds < maxRounds) {
                double roundMin = std::numeric_limits<double>::max();
                for (const auto &[key, value] : currentLow) {
                    (void)key;
                    roundMin = std::min(roundMin, value);
                }
                if (!std::isfinite(roundMin)) {
                    currentLow.clear();
                    break;
                }
                if (roundMin > currentCapTau) break;
                if (roundMin > currentTau) {
                    currentTau = roundMin;
                    lastTau = currentTau;
                    if (currentTau > tauMax) break;
                    if (currentTau > currentCapTau) break;
                    continue;
                }

                maxObservedMin = std::max(maxObservedMin, roundMin);
                phase.maxRoundMin = std::max(phase.maxRoundMin, roundMin);

                std::vector<SparseCliqueKey> frontierKeys;
                frontierKeys.reserve(currentLow.size() / 16 + 1);
                for (const auto &[key, value] : currentLow) {
                    if (std::abs(value - roundMin) <= 1e-9) frontierKeys.push_back(key);
                }
                if (frontierKeys.empty()) break;
                phase.frontierKeys += frontierKeys.size();

                std::vector<daf::Size> changedLeafIndex(activeLeaves.size(), std::numeric_limits<daf::Size>::max());
                std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
                std::vector<daf::Size> changedLeaf;
                removedKeyForLeaf.reserve(activeLeaves.size() / 16 + 1);
                changedLeaf.reserve(activeLeaves.size() / 16 + 1);

                daf::StaticVector<daf::Size> rmVerts;
                rmVerts.resize(r);
                for (const auto &rmKey : frontierKeys) {
                    for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
                    daf::intersect_dense_sets_multi(rmVerts, activeLeafByVertex,
                        [&](const daf::Size &leafId) {
                            if (leafId >= activeLeaves.size() || !leafAlive[leafId]) return;
                            auto &leafIdx = changedLeafIndex[leafId];
                            if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                                leafIdx = removedKeyForLeaf.size();
                                removedKeyForLeaf.emplace_back();
                                changedLeaf.push_back(leafId);
                            }
                            removedKeyForLeaf[leafIdx].push_back(rmKey);
                        });
                }
                rmVerts.free();

                robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> frontierSet;
                frontierSet.reserve(frontierKeys.size() * 2 + 1);
                for (const auto &key : frontierKeys) frontierSet.insert(key);

                robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> candidateKeys;
                candidateKeys.reserve(frontierKeys.size() * 8 + 1024);
                std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> localDelta;
                localDelta.reserve(1 << 20);
                uint64_t roundOldEntries = 0;
                uint64_t roundNewEntries = 0;
                uint64_t roundGeneratedSubleaves = 0;
                uint64_t roundExactLeaves = 0;
                const auto leafUpdateStart = std::chrono::high_resolution_clock::now();

                for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
                    const daf::Size leafId = changedLeaf[idx];
                    auto &leaf = activeLeaves[leafId];
                    if (!leafAlive[leafId] || leaf.size() < r) continue;
                    const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
                    if (removed.empty()) continue;
                    roundExactLeaves++;
                    phase.exactLeaves++;

                    collectLeafSparseKeys(leaf, r, s, candidateKeys);
                    collectLeafSparseKeys(leaf, r, s, candidateKeys);
                    auto newLeaves = splitLeafByRemovedKeys(
                        leaf, removed, r, s, localDelta, &roundOldEntries, &roundNewEntries);

                    removeLeafFromVertexIndex(leaf, leafId, activeLeafByVertex, r);
                    leafAlive[leafId] = 0;
                    roundGeneratedSubleaves += newLeaves.size();
                    totalSpawnedSubleaves += newLeaves.size();
                    phase.spawnedSubleaves += newLeaves.size();

                    for (auto &newLeaf : newLeaves) {
                        collectLeafSparseKeys(newLeaf, r, s, candidateKeys);
                        daf::Size newLeafId = activeLeaves.size();
                        activeLeaves.push_back(std::move(newLeaf));
                        leafAlive.push_back(1);
                        addLeafToVertexIndex(activeLeaves.back(), newLeafId, activeLeafByVertex, r);
                    }
                }

                for (const auto &key : frontierKeys) {
                    candidateKeys.erase(key);
                }

                std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> recomputedTouched;
                recomputedTouched.reserve(candidateKeys.size() / 2 + 1);
                for (const auto &key : candidateKeys) {
                    daf::StaticVector<daf::Size> verts;
                    verts.resize(r);
                    for (daf::CliqueSize j = 0; j < r; ++j) verts[j] = key.verts[j];
                    double support = 0.0;
                    daf::intersect_dense_sets_multi(verts, activeLeafByVertex,
                        [&](const daf::Size &leafId) {
                            if (leafId < activeLeaves.size() && leafAlive[leafId]) {
                                support += leafContributionForKey(activeLeaves[leafId], key, r, s);
                            }
                        });
                    verts.free();
                    if (support > 0.0) recomputedTouched.emplace(key, support);
                }

                std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> nextLow;
                nextLow.reserve(currentLow.size() + recomputedTouched.size());
                double nextMin = std::numeric_limits<double>::infinity();
                uint64_t nextFrontier = 0;
                for (const auto &[key, value] : currentLow) {
                    if (frontierSet.find(key) != frontierSet.end()) continue;
                    if (candidateKeys.find(key) != candidateKeys.end()) continue;
                    nextLow.emplace(key, value);
                    nextMin = std::min(nextMin, value);
                }
                for (const auto &[key, value] : recomputedTouched) {
                    if (value <= currentCapTau) {
                        nextLow[key] = value;
                        nextMin = std::min(nextMin, value);
                    }
                }
                if (std::isfinite(nextMin)) {
                    for (const auto &[key, value] : nextLow) {
                        (void)key;
                        if (std::abs(value - nextMin) <= 1e-9) nextFrontier++;
                    }
                }

                completedRounds++;
                phase.rounds++;
                totalExactLeaves += roundExactLeaves;
                totalCandidateKeys += candidateKeys.size();
                phase.candidateKeys += candidateKeys.size();

                std::cout << "---- Quotient buffered-low round " << completedRounds << " ----" << std::endl;
                std::cout << "  Tau/cap:                  " << currentTau << " / " << currentCapTau << std::endl;
                std::cout << "  Round min:                " << roundMin << std::endl;
                std::cout << "  Frontier keys:            " << frontierKeys.size() << std::endl;
                std::cout << "  Changed leaves:           " << changedLeaf.size() << std::endl;
                std::cout << "  Exact-handled leaves:     " << roundExactLeaves << std::endl;
                std::cout << "  Generated subleaves:      " << roundGeneratedSubleaves << std::endl;
                std::cout << "  Candidate touched keys:   " << candidateKeys.size() << std::endl;
                std::cout << "  Recomputed touched keys:  " << recomputedTouched.size() << std::endl;
                std::cout << "  Old local entries:        " << roundOldEntries << std::endl;
                std::cout << "  New local entries:        " << roundNewEntries << std::endl;
                std::cout << "  Next low keys:            " << nextLow.size() << std::endl;
                std::cout << "  Next frontier min:        "
                          << (std::isfinite(nextMin) ? nextMin : -1.0) << std::endl;
                std::cout << "  Next frontier size:       " << nextFrontier << std::endl;

                currentLow.swap(nextLow);
                lastTau = currentTau;
            }

            phase.remainingLowKeys = currentLow.size();
            phase.phaseMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - phaseStart).count());
            phaseStats.push_back(phase);
        }

        currentTau = currentCapTau + 1.0;
        lastTau = currentTau;
        currentLow.clear();
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();
    std::cout << "========== Quotient Buffered Adaptive Low-Support ==========" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Tau start/max:            " << tauStart << " / " << tauMax << std::endl;
    std::cout << "  Exact over band:          " << overExactBand << std::endl;
    std::cout << "  Lookahead:                " << lookahead << std::endl;
    std::cout << "  Auto lookahead:           " << (autoLookahead ? "YES" : "NO") << std::endl;
    if (autoLookahead) {
        std::cout << "  Auto threshold keys:      " << autoThresholdKeys << std::endl;
        std::cout << "  Auto raw threshold:       " << autoRawCliqueThreshold << std::endl;
    }
    std::cout << "  Rebuild count:            " << rebuildCount << std::endl;
    std::cout << "  Completed rounds:         " << completedRounds << std::endl;
    std::cout << "  Total exact leaves:       " << totalExactLeaves << std::endl;
    std::cout << "  Total candidate keys:     " << totalCandidateKeys << std::endl;
    std::cout << "  Total spawned subleaves:  " << totalSpawnedSubleaves << std::endl;
    std::cout << "  Max observed round min:   " << maxObservedMin << std::endl;
    std::cout << "  Total rebuild time:       " << totalRebuildMs << " ms" << std::endl;
    std::cout << "  Final tau:                " << lastTau << std::endl;
    std::cout << "  Prototype time:           " << elapsedMs << " ms" << std::endl;
    for (const auto &phase : phaseStats) {
        std::cout << "  Phase tau=" << phase.tau
                  << " cap=" << phase.capTau
                  << " lookahead=" << phase.lookahead
                  << " rebuild_ms=" << phase.rebuildMs
                  << " phase_ms=" << phase.phaseMs
                  << " init_low=" << phase.initialLowKeys
                  << " init_over=" << phase.initialOverTauKeys
                  << " rounds=" << phase.rounds
                  << " frontier_keys=" << phase.frontierKeys
                  << " exact_leaves=" << phase.exactLeaves
                  << " cand_keys=" << phase.candidateKeys
                  << " spawned=" << phase.spawnedSubleaves
                  << " max_round_min=" << phase.maxRoundMin
                  << " remain_low=" << phase.remainingLowKeys
                  << " kept=" << phase.buildStats.keptUpdates
                  << " evicted=" << phase.buildStats.evictedUpdates
                  << " skipped=" << phase.buildStats.skippedUpdates
                  << std::endl;
    }
    std::cout << "============================================================" << std::endl;
}

static void maybeRunBandedBufferedAdaptiveLowSupportPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_BUFFERED_BANDED_LOW")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "======= Quotient Banded Buffered Adaptive Low-Support =======" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip banded buffered adaptive prototype" << std::endl;
        std::cout << "  Reason: raw clique layer is too large for this prototype" << std::endl;
        std::cout << "=============================================================" << std::endl;
        return;
    }

    double tauStart = 1.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_START")) {
        tauStart = std::max(1.0, std::atof(env));
    }
    double tauMax = 8.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_MAX")) {
        tauMax = std::max(tauStart, std::atof(env));
    }
    int overExactBand = 0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_OVER_BAND")) {
        overExactBand = std::max(0, std::atoi(env));
    }
    int lookahead = 1;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_LOOKAHEAD")) {
        lookahead = std::max(1, std::atoi(env));
    }
    const bool autoLookahead = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO") != nullptr;
    uint64_t autoThresholdKeys = 300000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_THRESHOLD_KEYS")) {
        autoThresholdKeys = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoRawCliqueThreshold = 1000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_RAW_THRESHOLD")) {
        autoRawCliqueThreshold = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoWideRawCliqueThreshold = 20000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_WIDE_RAW_THRESHOLD")) {
        autoWideRawCliqueThreshold = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoWideExplicitMin = 10000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_WIDE_EXPLICIT_MIN")) {
        autoWideExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoWideExplicitMax = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_WIDE_EXPLICIT_MAX")) {
        autoWideExplicitMax = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    int autoCompactMaxLeafSize = 20;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_MAX_LEAF")) {
        autoCompactMaxLeafSize = std::max(1, std::atoi(env));
    }
    uint64_t autoCompactExplicitMax = 20000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_EXPLICIT_MAX")) {
        autoCompactExplicitMax = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    double autoCompactCompressionMax = 12.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_COMPRESSION_MAX")) {
        autoCompactCompressionMax = std::max(1.0, std::atof(env));
    }
    int autoCompactLookahead = 4;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_LOOKAHEAD")) {
        autoCompactLookahead = std::max(1, std::atoi(env));
    }
    int autoCacheRebandMinLeafSize = 40;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_MIN_LEAF")) {
        autoCacheRebandMinLeafSize = std::max(1, std::atoi(env));
    }
    uint64_t autoCacheRebandExplicitMin = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_EXPLICIT_MIN")) {
        autoCacheRebandExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    double autoCacheRebandCompressionMin = 50.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_COMPRESSION_MIN")) {
        autoCacheRebandCompressionMin = std::max(1.0, std::atof(env));
    }
    int autoCacheRebandLookahead = 3;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_LOOKAHEAD")) {
        autoCacheRebandLookahead = std::max(1, std::atoi(env));
    }
    const bool enableAutoCacheRebandHeavy =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_HEAVY") != nullptr;
    int autoTwoPhaseMinLeafSize = 40;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE_MIN_LEAF")) {
        autoTwoPhaseMinLeafSize = std::max(1, std::atoi(env));
    }
    uint64_t autoTwoPhaseExplicitMin = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE_EXPLICIT_MIN")) {
        autoTwoPhaseExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool enableAutoTwoPhase =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE") != nullptr;
    double autoFullBandOverRatio = 8.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_FULL_OVER_RATIO")) {
        autoFullBandOverRatio = std::max(1.0, std::atof(env));
    }
    double autoFullBandRebuildFactor = 2.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_FULL_REBUILD_FACTOR")) {
        autoFullBandRebuildFactor = std::max(1.0, std::atof(env));
    }
    int autoFullBandMinRemainingTau = 2;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_FULL_MIN_REMAINING")) {
        autoFullBandMinRemainingTau = std::max(1, std::atoi(env));
    }
    int maxRounds = 200;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_MAX")) {
        maxRounds = std::max(1, std::atoi(env));
    }
    const bool enableDeltaFast = std::getenv("PIVOTER_QUOTIENT_BUCKETED_DELTA_FAST") != nullptr;

    std::vector<std::vector<TreeGraphNode>> activeLeaves;
    activeLeaves.reserve(tree.adj_list.size() * 2 + 16);
    for (const auto &leaf : tree.adj_list) {
        activeLeaves.emplace_back(leaf.begin(), leaf.end());
    }
    std::vector<uint8_t> leafAlive(activeLeaves.size(), 1);
    std::vector<robin_hood::unordered_flat_set<daf::Size>> activeLeafByVertex(edgeGraph.n);
    std::vector<robin_hood::unordered_flat_set<daf::Size>> activePivotLeafByVertex(edgeGraph.n);
    std::vector<uint16_t> activeKeepCount(activeLeaves.size(), 0);
    std::vector<uint16_t> activePivotCount(activeLeaves.size(), 0);
    for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
        addLeafToVertexAndPivotIndex(activeLeaves[leafId], leafId,
                                     activeLeafByVertex, activePivotLeafByVertex,
                                     activeKeepCount, activePivotCount, r);
    }

    int completedRounds = 0;
    int rebuildCount = 0;
    uint64_t totalExactLeaves = 0;
    uint64_t totalCandidateKeys = 0;
    uint64_t totalSpawnedSubleaves = 0;
    double maxObservedMin = 0.0;
    double lastTau = tauStart;
    uint64_t totalRebuildMs = 0;
    std::vector<AdaptiveLowPhaseStats> phaseStats;

    double currentTau = tauStart;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> currentLow;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> bufferedLow;
    double currentCapTau = -1.0;
    double currentWindowCapTau = -1.0;

    while (currentTau <= tauMax && completedRounds < maxRounds) {
        if (currentLow.empty() && !bufferedLow.empty()) {
            currentTau = std::max(currentTau, currentCapTau + 1.0);
            if (currentTau > tauMax) break;
            currentLow.swap(bufferedLow);
            currentCapTau = currentWindowCapTau;
            lastTau = currentTau;
        }

        if (currentLow.empty()) {
            int phaseLookahead = lookahead;
            if (autoLookahead) {
                phaseLookahead = 1;
                if (phaseStats.empty()) {
                    if (static_cast<uint64_t>(rawCliqueCount) <= autoRawCliqueThreshold) {
                        phaseLookahead = lookahead;
                    }
                } else {
                    const auto &prev = phaseStats.back();
                    if (prev.initialLowKeys <= autoThresholdKeys) phaseLookahead = lookahead;
                }
            }
            const int activeExtra = std::max(1, phaseLookahead - 1);
            currentCapTau = std::min(tauMax, currentTau + static_cast<double>(activeExtra));
            currentWindowCapTau = std::min(tauMax, currentTau + static_cast<double>(phaseLookahead));

            AdaptiveLowPhaseStats phase;
            phase.tau = currentTau;
            phase.capTau = currentCapTau;
            phase.windowCapTau = currentWindowCapTau;
            phase.lookahead = phaseLookahead;

            BandedLowSupportBuildResult build;
            auto phaseStart = std::chrono::high_resolution_clock::now();
            auto rebuildStart = std::chrono::high_resolution_clock::now();
            daf::timeCount("banded buffered adaptive low-support rebuild", [&]() {
                build = buildBandedLowSupportLayerActive(
                    activeLeaves, leafAlive, r, s,
                    currentCapTau, currentWindowCapTau, rawCliqueCount,
                    false,
                    std::min(tauMax, currentWindowCapTau + static_cast<double>(overExactBand)));
            });
            phase.rebuildMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - rebuildStart).count());
            totalRebuildMs += phase.rebuildMs;
            rebuildCount++;
            phase.initialLowKeys = build.lowSupport.size();
            phase.initialBufferedKeys = build.bufferedSupport.size();
            phase.initialOverTauKeys = build.overWindow.size();
            phase.buildStats = build.stats;
            currentLow = std::move(build.lowSupport);
            bufferedLow = std::move(build.bufferedSupport);

            if (currentLow.empty() && bufferedLow.empty()) {
                phase.phaseMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - phaseStart).count());
                phaseStats.push_back(phase);
                currentTau = currentWindowCapTau + 1.0;
                lastTau = currentTau;
                continue;
            }

            while (completedRounds < maxRounds) {
                if (currentLow.empty()) {
                    if (bufferedLow.empty()) break;
                    currentTau = std::max(currentTau, currentCapTau + 1.0);
                    if (currentTau > tauMax) break;
                    currentLow.swap(bufferedLow);
                    currentCapTau = currentWindowCapTau;
                    lastTau = currentTau;
                }

                double roundMin = std::numeric_limits<double>::max();
                for (const auto &[key, value] : currentLow) {
                    (void)key;
                    roundMin = std::min(roundMin, value);
                }
                if (!std::isfinite(roundMin)) {
                    currentLow.clear();
                    bufferedLow.clear();
                    break;
                }
                if (roundMin > currentWindowCapTau) break;
                if (roundMin > currentTau) {
                    currentTau = roundMin;
                    lastTau = currentTau;
                    if (currentTau > tauMax) break;
                    if (currentTau > currentCapTau && !bufferedLow.empty()) {
                        currentLow.swap(bufferedLow);
                        currentCapTau = currentWindowCapTau;
                    }
                    if (currentTau > currentWindowCapTau) break;
                    continue;
                }

                maxObservedMin = std::max(maxObservedMin, roundMin);
                phase.maxRoundMin = std::max(phase.maxRoundMin, roundMin);

                std::vector<SparseCliqueKey> frontierKeys;
                frontierKeys.reserve(currentLow.size() / 16 + 1);
                for (const auto &[key, value] : currentLow) {
                    if (std::abs(value - roundMin) <= 1e-9) frontierKeys.push_back(key);
                }
                if (frontierKeys.empty()) break;
                phase.frontierKeys += frontierKeys.size();

                std::vector<daf::Size> changedLeafIndex(activeLeaves.size(), std::numeric_limits<daf::Size>::max());
                std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
                std::vector<daf::Size> changedLeaf;
                removedKeyForLeaf.reserve(activeLeaves.size() / 16 + 1);
                changedLeaf.reserve(activeLeaves.size() / 16 + 1);

                daf::StaticVector<daf::Size> rmVerts;
                rmVerts.resize(r);
                for (const auto &rmKey : frontierKeys) {
                    for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
                    daf::intersect_dense_sets_multi(rmVerts, activeLeafByVertex,
                        [&](const daf::Size &leafId) {
                            if (leafId >= activeLeaves.size() || !leafAlive[leafId]) return;
                            auto &leafIdx = changedLeafIndex[leafId];
                            if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                                leafIdx = removedKeyForLeaf.size();
                                removedKeyForLeaf.emplace_back();
                                changedLeaf.push_back(leafId);
                            }
                            removedKeyForLeaf[leafIdx].push_back(rmKey);
                        });
                }
                rmVerts.free();

                robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> frontierSet;
                frontierSet.reserve(frontierKeys.size() * 2 + 1);
                for (const auto &key : frontierKeys) frontierSet.insert(key);

                robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> candidateKeys;
                candidateKeys.reserve(frontierKeys.size() * 8 + 1024);
                std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> localDelta;
                localDelta.reserve(1 << 20);
                uint64_t roundOldEntries = 0;
                uint64_t roundNewEntries = 0;
                uint64_t roundGeneratedSubleaves = 0;
                uint64_t roundExactLeaves = 0;
                const auto leafUpdateStart = std::chrono::high_resolution_clock::now();

                for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
                    const daf::Size leafId = changedLeaf[idx];
                    auto &leaf = activeLeaves[leafId];
                    if (!leafAlive[leafId] || leaf.size() < r) continue;
                    const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
                    if (removed.empty()) continue;
                    roundExactLeaves++;
                    phase.exactLeaves++;

                    collectLeafSparseKeys(leaf, r, s, candidateKeys);
                    auto newLeaves = splitLeafByRemovedKeys(
                        leaf, removed, r, s, localDelta, &roundOldEntries, &roundNewEntries);

                    removeLeafFromVertexIndex(leaf, leafId, activeLeafByVertex, r);
                    leafAlive[leafId] = 0;
                    roundGeneratedSubleaves += newLeaves.size();
                    totalSpawnedSubleaves += newLeaves.size();
                    phase.spawnedSubleaves += newLeaves.size();

                    for (auto &newLeaf : newLeaves) {
                        collectLeafSparseKeys(newLeaf, r, s, candidateKeys);
                        daf::Size newLeafId = activeLeaves.size();
                        activeLeaves.push_back(std::move(newLeaf));
                        leafAlive.push_back(1);
                        addLeafToVertexIndex(activeLeaves.back(), newLeafId, activeLeafByVertex, r);
                    }
                }

                for (const auto &key : frontierKeys) {
                    candidateKeys.erase(key);
                }

                std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> recomputedTouched;
                recomputedTouched.reserve(candidateKeys.size() / 2 + 1);
                for (const auto &key : candidateKeys) {
                    daf::StaticVector<daf::Size> verts;
                    verts.resize(r);
                    for (daf::CliqueSize j = 0; j < r; ++j) verts[j] = key.verts[j];
                    double support = 0.0;
                    daf::intersect_dense_sets_multi(verts, activeLeafByVertex,
                        [&](const daf::Size &leafId) {
                            if (leafId < activeLeaves.size() && leafAlive[leafId]) {
                                support += leafContributionForKey(activeLeaves[leafId], key, r, s);
                            }
                        });
                    verts.free();
                    if (support > 0.0) recomputedTouched.emplace(key, support);
                }

                std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> nextLow;
                std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> nextBuffered;
                nextLow.reserve(currentLow.size() + recomputedTouched.size() / 2 + 1);
                nextBuffered.reserve(bufferedLow.size() + recomputedTouched.size() / 2 + 1);
                double nextMin = std::numeric_limits<double>::infinity();
                uint64_t nextFrontier = 0;

                for (const auto &[key, value] : currentLow) {
                    if (frontierSet.find(key) != frontierSet.end()) continue;
                    if (candidateKeys.find(key) != candidateKeys.end()) continue;
                    nextLow.emplace(key, value);
                    nextMin = std::min(nextMin, value);
                }
                for (const auto &[key, value] : bufferedLow) {
                    if (candidateKeys.find(key) != candidateKeys.end()) continue;
                    nextBuffered.emplace(key, value);
                }
                for (const auto &[key, value] : recomputedTouched) {
                    if (value <= currentCapTau) {
                        nextLow[key] = value;
                        nextMin = std::min(nextMin, value);
                    } else if (value <= currentWindowCapTau) {
                        nextBuffered[key] = value;
                    }
                }
                if (std::isfinite(nextMin)) {
                    for (const auto &[key, value] : nextLow) {
                        (void)key;
                        if (std::abs(value - nextMin) <= 1e-9) nextFrontier++;
                    }
                }

                completedRounds++;
                phase.rounds++;
                totalExactLeaves += roundExactLeaves;
                totalCandidateKeys += candidateKeys.size();
                phase.candidateKeys += candidateKeys.size();

                std::cout << "---- Quotient banded-buffered-low round " << completedRounds << " ----" << std::endl;
                std::cout << "  Tau/cap/window:           " << currentTau << " / "
                          << currentCapTau << " / " << currentWindowCapTau << std::endl;
                std::cout << "  Round min:                " << roundMin << std::endl;
                std::cout << "  Frontier keys:            " << frontierKeys.size() << std::endl;
                std::cout << "  Changed leaves:           " << changedLeaf.size() << std::endl;
                std::cout << "  Exact-handled leaves:     " << roundExactLeaves << std::endl;
                std::cout << "  Generated subleaves:      " << roundGeneratedSubleaves << std::endl;
                std::cout << "  Candidate touched keys:   " << candidateKeys.size() << std::endl;
                std::cout << "  Recomputed touched keys:  " << recomputedTouched.size() << std::endl;
                std::cout << "  Old local entries:        " << roundOldEntries << std::endl;
                std::cout << "  New local entries:        " << roundNewEntries << std::endl;
                std::cout << "  Next low keys:            " << nextLow.size() << std::endl;
                std::cout << "  Next buffered keys:       " << nextBuffered.size() << std::endl;
                std::cout << "  Next frontier min:        "
                          << (std::isfinite(nextMin) ? nextMin : -1.0) << std::endl;
                std::cout << "  Next frontier size:       " << nextFrontier << std::endl;

                currentLow.swap(nextLow);
                bufferedLow.swap(nextBuffered);
                lastTau = currentTau;
            }

            phase.remainingLowKeys = currentLow.size() + bufferedLow.size();
            phase.phaseMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - phaseStart).count());
            phaseStats.push_back(phase);
        }

        currentTau = currentWindowCapTau + 1.0;
        lastTau = currentTau;
        currentLow.clear();
        bufferedLow.clear();
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();
    std::cout << "======= Quotient Banded Buffered Adaptive Low-Support =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Tau start/max:            " << tauStart << " / " << tauMax << std::endl;
    std::cout << "  Lookahead:                " << lookahead << std::endl;
    std::cout << "  Auto lookahead:           " << (autoLookahead ? "YES" : "NO") << std::endl;
    if (autoLookahead) {
        std::cout << "  Auto threshold keys:      " << autoThresholdKeys << std::endl;
        std::cout << "  Auto raw threshold:       " << autoRawCliqueThreshold << std::endl;
    }
    std::cout << "  Rebuild count:            " << rebuildCount << std::endl;
    std::cout << "  Completed rounds:         " << completedRounds << std::endl;
    std::cout << "  Total exact leaves:       " << totalExactLeaves << std::endl;
    std::cout << "  Total candidate keys:     " << totalCandidateKeys << std::endl;
    std::cout << "  Total spawned subleaves:  " << totalSpawnedSubleaves << std::endl;
    std::cout << "  Max observed round min:   " << maxObservedMin << std::endl;
    std::cout << "  Total rebuild time:       " << totalRebuildMs << " ms" << std::endl;
    std::cout << "  Final tau:                " << lastTau << std::endl;
    std::cout << "  Prototype time:           " << elapsedMs << " ms" << std::endl;
    for (const auto &phase : phaseStats) {
        std::cout << "  Phase tau=" << phase.tau
                  << " cap=" << phase.capTau
                  << " window=" << phase.windowCapTau
                  << " lookahead=" << phase.lookahead
                  << " rebuild_ms=" << phase.rebuildMs
                  << " phase_ms=" << phase.phaseMs
                  << " init_low=" << phase.initialLowKeys
                  << " init_buffered=" << phase.initialBufferedKeys
                  << " init_over_exact=" << phase.initialOverExactKeys
                  << " init_over=" << phase.initialOverTauKeys
                  << " rounds=" << phase.rounds
                  << " frontier_keys=" << phase.frontierKeys
                  << " exact_leaves=" << phase.exactLeaves
                  << " cand_keys=" << phase.candidateKeys
                  << " spawned=" << phase.spawnedSubleaves
                  << " max_round_min=" << phase.maxRoundMin
                  << " remain_low=" << phase.remainingLowKeys
                  << " kept=" << phase.buildStats.keptUpdates
                  << " evicted=" << phase.buildStats.evictedUpdates
                  << " skipped=" << phase.buildStats.skippedUpdates
                  << std::endl;
    }
    std::cout << "=============================================================" << std::endl;
}

static void maybeRunBucketedBandedBufferedAdaptiveLowSupportPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_BUFFERED_BUCKETED_LOW")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "==== Quotient Bucketed Banded Buffered Adaptive Low-Support ====" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip bucketed banded buffered adaptive prototype" << std::endl;
        std::cout << "  Reason: raw clique layer is too large for this prototype" << std::endl;
        std::cout << "===============================================================" << std::endl;
        return;
    }

    double tauStart = 1.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_START")) {
        tauStart = std::max(1.0, std::atof(env));
    }
    double tauMax = 8.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_MAX")) {
        tauMax = std::max(tauStart, std::atof(env));
    }
    int overExactBand = 0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_OVER_BAND")) {
        overExactBand = std::max(0, std::atoi(env));
    }
    int lookahead = 1;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_LOOKAHEAD")) {
        lookahead = std::max(1, std::atoi(env));
    }
    const bool autoLookahead = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO") != nullptr;
    uint64_t autoThresholdKeys = 300000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_THRESHOLD_KEYS")) {
        autoThresholdKeys = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoRawCliqueThreshold = 1000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_RAW_THRESHOLD")) {
        autoRawCliqueThreshold = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoWideRawCliqueThreshold = 20000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_WIDE_RAW_THRESHOLD")) {
        autoWideRawCliqueThreshold = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoWideExplicitMin = 10000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_WIDE_EXPLICIT_MIN")) {
        autoWideExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoWideExplicitMax = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_WIDE_EXPLICIT_MAX")) {
        autoWideExplicitMax = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    int autoCompactMaxLeafSize = 20;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_MAX_LEAF")) {
        autoCompactMaxLeafSize = std::max(1, std::atoi(env));
    }
    uint64_t autoCompactExplicitMax = 20000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_EXPLICIT_MAX")) {
        autoCompactExplicitMax = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    double autoCompactCompressionMax = 12.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_COMPRESSION_MAX")) {
        autoCompactCompressionMax = std::max(1.0, std::atof(env));
    }
    int autoCompactLookahead = 4;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_COMPACT_LOOKAHEAD")) {
        autoCompactLookahead = std::max(1, std::atoi(env));
    }
    int autoCacheRebandMinLeafSize = 40;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_MIN_LEAF")) {
        autoCacheRebandMinLeafSize = std::max(1, std::atoi(env));
    }
    uint64_t autoCacheRebandExplicitMin = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_EXPLICIT_MIN")) {
        autoCacheRebandExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    double autoCacheRebandCompressionMin = 50.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_COMPRESSION_MIN")) {
        autoCacheRebandCompressionMin = std::max(1.0, std::atof(env));
    }
    int autoCacheRebandLookahead = 3;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_LOOKAHEAD")) {
        autoCacheRebandLookahead = std::max(1, std::atoi(env));
    }
    int autoTwoPhaseMinLeafSize = 40;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE_MIN_LEAF")) {
        autoTwoPhaseMinLeafSize = std::max(1, std::atoi(env));
    }
    uint64_t autoTwoPhaseExplicitMin = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE_EXPLICIT_MIN")) {
        autoTwoPhaseExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool enableAutoTwoPhase =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE") != nullptr;
    const bool enableAutoCacheRebandHeavy =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_HEAVY") != nullptr;
    const bool forceSingleOverMapHeavy =
        std::getenv("PIVOTER_QUOTIENT_BUILD_SINGLE_OVER_MAP_HEAVY") != nullptr;
    const bool disableAutoSingleOverMapHeavy =
        std::getenv("PIVOTER_QUOTIENT_BUILD_SINGLE_OVER_MAP_HEAVY_OFF") != nullptr;
    const bool forceOverSetHeavy =
        std::getenv("PIVOTER_QUOTIENT_BUILD_OVER_SET_HEAVY") != nullptr;
    double autoFullBandOverRatio = 8.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_FULL_OVER_RATIO")) {
        autoFullBandOverRatio = std::max(1.0, std::atof(env));
    }
    double autoFullBandRebuildFactor = 2.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_FULL_REBUILD_FACTOR")) {
        autoFullBandRebuildFactor = std::max(1.0, std::atof(env));
    }
    int autoFullBandMinRemainingTau = 2;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_FULL_MIN_REMAINING")) {
        autoFullBandMinRemainingTau = std::max(1, std::atoi(env));
    }
    int maxRounds = 200;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_MAX")) {
        maxRounds = std::max(1, std::atoi(env));
    }
    const bool enableCacheReband =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND") != nullptr;
    const bool enableLimitedCacheReband =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND_LIMITED") != nullptr;
    const bool enableAutoLimitedCacheReband =
        std::getenv("PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND_LIMITED_AUTO") != nullptr;
    uint64_t autoLimitedExplicitMin = 5000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND_LIMITED_AUTO_EXPLICIT_MIN")) {
        autoLimitedExplicitMin = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    uint64_t autoLimitedExplicitMax = 50000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND_LIMITED_AUTO_EXPLICIT_MAX")) {
        autoLimitedExplicitMax = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    double autoLimitedCompressionMax = 20.0;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND_LIMITED_AUTO_COMPRESSION_MAX")) {
        autoLimitedCompressionMax = std::max(1.0, std::atof(env));
    }
    const bool forceOverLowerBound =
        std::getenv("PIVOTER_QUOTIENT_BUCKETED_DELTA_OVER_LOWER") != nullptr;
    const bool disableAutoOverLowerBound =
        std::getenv("PIVOTER_QUOTIENT_BUCKETED_DELTA_OVER_LOWER_OFF") != nullptr;
    const bool useOverLowerBound =
        forceOverLowerBound || (enableCacheReband && !disableAutoOverLowerBound);

    std::vector<std::vector<TreeGraphNode>> activeLeaves;
    activeLeaves.reserve(tree.adj_list.size() * 2 + 16);
    for (const auto &leaf : tree.adj_list) {
        activeLeaves.emplace_back(leaf.begin(), leaf.end());
    }
    std::vector<uint8_t> leafAlive(activeLeaves.size(), 1);
    std::vector<robin_hood::unordered_flat_set<daf::Size>> activeLeafByVertex(edgeGraph.n);
    for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
        addLeafToVertexIndex(activeLeaves[leafId], leafId, activeLeafByVertex, r);
    }

    int completedRounds = 0;
    int rebuildCount = 0;
    uint64_t totalExactLeaves = 0;
    uint64_t totalCandidateKeys = 0;
    uint64_t totalSpawnedSubleaves = 0;
    double maxObservedMin = 0.0;
    double lastTau = tauStart;
    uint64_t totalRebuildMs = 0;
    std::vector<AdaptiveLowPhaseStats> phaseStats;
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> overSupportCache;
    overSupportCache.reserve(1 << 20);
    robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> persistentOverWindow;
    bool autoForceFullBandNext = false;
    int autoSuggestedLookaheadNext = 0;
    const bool hasSurveySummary = gQuotientSurveySummary.valid &&
                                  gQuotientSurveySummary.r == r &&
                                  gQuotientSurveySummary.s == s;
    const uint64_t surveyExplicitEntries =
        hasSurveySummary ? gQuotientSurveySummary.totalExplicitEntries : 0ULL;
    const int surveyMaxLeafSize =
        hasSurveySummary ? gQuotientSurveySummary.maxLeafSize : 0;
    const double surveyCompression =
        (hasSurveySummary && gQuotientSurveySummary.totalQuotientClasses > 0)
            ? static_cast<double>(gQuotientSurveySummary.totalExplicitEntries) /
                  static_cast<double>(gQuotientSurveySummary.totalQuotientClasses)
            : std::numeric_limits<double>::infinity();

    double currentTau = tauStart;
    BucketedBandedLowState state;
    double currentCapTau = -1.0;
    double currentWindowCapTau = -1.0;
    double currentOverExactCapTau = -1.0;

    while (currentTau <= tauMax && completedRounds < maxRounds) {
        if (state.lowValue.empty() && !state.bufferedValue.empty()) {
            currentTau = std::max(currentTau, currentCapTau + 1.0);
            if (currentTau > tauMax) break;
            state.lowValue.swap(state.bufferedValue);
            state.lowBuckets.swap(state.bufferedBuckets);
            currentCapTau = currentWindowCapTau;
            lastTau = currentTau;
        }

        if (state.lowValue.empty()) {
            int phaseLookahead = lookahead;
            bool phaseFullBand = false;
            if (autoLookahead) {
                phaseLookahead = 1;
                if (phaseStats.empty()) {
                    if (static_cast<uint64_t>(rawCliqueCount) <= autoRawCliqueThreshold) {
                        phaseLookahead = std::max(
                            1, static_cast<int>(std::ceil(tauMax - currentTau + 1.0)));
                        phaseFullBand = true;
                    } else if (enableAutoCacheRebandHeavy &&
                               enableCacheReband &&
                               hasSurveySummary &&
                               surveyMaxLeafSize >= autoCacheRebandMinLeafSize &&
                               surveyExplicitEntries >= autoCacheRebandExplicitMin &&
                               surveyCompression >= autoCacheRebandCompressionMin) {
                        phaseLookahead = std::max(phaseLookahead, autoCacheRebandLookahead);
                    } else if (hasSurveySummary &&
                               surveyMaxLeafSize <= autoCompactMaxLeafSize &&
                               surveyExplicitEntries <= autoCompactExplicitMax &&
                               surveyCompression <= autoCompactCompressionMax) {
                        phaseLookahead = std::max(phaseLookahead, autoCompactLookahead);
                    } else if (hasSurveySummary &&
                               static_cast<uint64_t>(rawCliqueCount) <= autoWideRawCliqueThreshold &&
                               surveyExplicitEntries >= autoWideExplicitMin &&
                               surveyExplicitEntries <= autoWideExplicitMax) {
                        phaseLookahead = std::max(phaseLookahead, lookahead + 1);
                    }
                } else if (autoForceFullBandNext) {
                    phaseLookahead = std::max(
                        1, static_cast<int>(std::ceil(tauMax - currentTau + 1.0)));
                    phaseFullBand = true;
                    autoForceFullBandNext = false;
                    autoSuggestedLookaheadNext = 0;
                } else if (autoSuggestedLookaheadNext > 0) {
                    phaseLookahead = std::max(1, autoSuggestedLookaheadNext);
                    if (currentTau + static_cast<double>(phaseLookahead) >= tauMax - 1e-9) {
                        phaseFullBand = true;
                    }
                    autoSuggestedLookaheadNext = 0;
                } else {
                    const auto &prev = phaseStats.back();
                    if (prev.initialLowKeys <= autoThresholdKeys) phaseLookahead = lookahead;
                }
            }
            const int activeExtra = std::max(1, phaseLookahead - 1);
            currentCapTau = std::min(tauMax, currentTau + static_cast<double>(activeExtra));
            currentWindowCapTau = std::min(tauMax, currentTau + static_cast<double>(phaseLookahead));
            currentOverExactCapTau = std::min(
                tauMax, currentWindowCapTau + static_cast<double>(overExactBand));

            AdaptiveLowPhaseStats phase;
            phase.tau = currentTau;
            phase.capTau = currentCapTau;
            phase.windowCapTau = currentWindowCapTau;
            phase.lookahead = phaseLookahead;
            phase.fullBand = phaseFullBand;
            phase.twoPhaseBuild = false;
            phase.cacheRebandBuild = false;

            BandedLowSupportBuildResult build;
            auto phaseStart = std::chrono::high_resolution_clock::now();
            auto rebuildStart = std::chrono::high_resolution_clock::now();
            const bool useOverCache = std::getenv("PIVOTER_QUOTIENT_BUCKETED_OVER_CACHE_OFF") == nullptr;
            const bool seedCacheFromBuild = useOverCache &&
                std::getenv("PIVOTER_QUOTIENT_BUCKETED_SEED_CACHE_OFF") == nullptr;
            bool usedCacheReband = false;
            double trackedOverCapTau = -1.0;
            const bool useLimitedCacheReband =
                (enableLimitedCacheReband ||
                 (enableAutoLimitedCacheReband &&
                  hasSurveySummary &&
                  surveyExplicitEntries >= autoLimitedExplicitMin &&
                  surveyExplicitEntries <= autoLimitedExplicitMax &&
                  surveyCompression <= autoLimitedCompressionMax)) &&
                enableCacheReband &&
                phaseStats.empty() &&
                !phaseFullBand;
            if (useLimitedCacheReband) {
                trackedOverCapTau = tauMax;
            }
            const bool useTwoPhaseFirstBuild =
                enableAutoTwoPhase &&
                phaseStats.empty() &&
                !phaseFullBand &&
                phaseLookahead == 1 &&
                std::abs(currentTau - 1.0) <= 1e-9 &&
                std::abs(currentCapTau - currentWindowCapTau) <= 1e-9 &&
                hasSurveySummary &&
                surveyMaxLeafSize >= autoTwoPhaseMinLeafSize &&
                surveyExplicitEntries >= autoTwoPhaseExplicitMin;
            const bool useSingleOverMapHeavy =
                (forceSingleOverMapHeavy ||
                 (enableCacheReband && !disableAutoSingleOverMapHeavy)) &&
                enableCacheReband &&
                phaseStats.empty() &&
                !phaseFullBand &&
                phaseLookahead == 1 &&
                std::abs(currentTau - 1.0) <= 1e-9 &&
                std::abs(currentCapTau - currentWindowCapTau) <= 1e-9 &&
                hasSurveySummary &&
                surveyMaxLeafSize >= autoCacheRebandMinLeafSize &&
                surveyExplicitEntries >= autoCacheRebandExplicitMin &&
                surveyCompression >= autoCacheRebandCompressionMin;
            const bool useOverSetHeavy =
                forceOverSetHeavy &&
                phaseStats.empty() &&
                !phaseFullBand &&
                phaseLookahead == 1 &&
                std::abs(currentTau - 1.0) <= 1e-9 &&
                std::abs(currentCapTau - currentWindowCapTau) <= 1e-9;
            daf::timeCount("bucketed banded buffered adaptive low-support rebuild", [&]() {
                if (enableCacheReband &&
                    currentOverExactCapTau <= currentWindowCapTau + 1e-9 &&
                    rebandBucketedStateFromExactCache(
                        persistentOverWindow, overSupportCache,
                        currentCapTau, currentWindowCapTau, currentOverExactCapTau, state)) {
                    usedCacheReband = true;
                } else if (useTwoPhaseFirstBuild) {
                    auto twoPhase = buildLowSupportTwoPhase(
                        activeLeaves, leafAlive, activeLeafByVertex,
                        r, s, currentCapTau, rawCliqueCount);
                    build.lowSupport = std::move(twoPhase.lowSupport);
                    build.overWindow = std::move(twoPhase.overTau);
                    build.stats = twoPhase.stats;
                } else {
                    const bool buildUseOverLowerBound =
                        useOverSetHeavy ? false : useOverLowerBound;
                    const bool buildSeedCacheFromBuild =
                        useOverSetHeavy ? false : seedCacheFromBuild;
                    build = buildBandedLowSupportLayerActive(
                        activeLeaves, leafAlive, r, s,
                        currentCapTau, currentWindowCapTau, rawCliqueCount,
                        buildUseOverLowerBound, currentOverExactCapTau, buildSeedCacheFromBuild,
                        trackedOverCapTau,
                        useSingleOverMapHeavy,
                        useSingleOverMapHeavy);
                }
            });
            phase.twoPhaseBuild = useTwoPhaseFirstBuild;
            phase.cacheRebandBuild = usedCacheReband;
            phase.rebuildMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - rebuildStart).count());
            totalRebuildMs += phase.rebuildMs;
            rebuildCount++;
            if (usedCacheReband) {
                phase.initialLowKeys = state.lowValue.size();
                phase.initialBufferedKeys = state.bufferedValue.size();
                phase.initialOverExactKeys = state.overValue.size();
                phase.initialOverTauKeys = state.overWindow.size();
                phase.initialTrackedOverKeys = std::max<uint64_t>(
                    phase.initialOverTauKeys,
                    static_cast<uint64_t>(overSupportCache.size()));
            } else {
                phase.singleOverMapBuild = build.usedSingleOverMap;
                phase.cappedLowBuild = build.usedCappedLowState;
                phase.initialLowKeys = build.lowSupport.size();
                phase.initialBufferedKeys = build.bufferedSupport.size();
                phase.initialOverExactKeys = build.overExactSupport.size();
                phase.initialOverTauKeys = build.overWindow.size();
                phase.buildStats = build.stats;
                // Seed overSupportCache from exact over-window values computed during build
                overSupportCache.clear();
                if (seedCacheFromBuild) {
                    if (!build.overExactSeed.empty()) {
                        overSupportCache = std::move(build.overExactSeed);
                    } else if (!build.overLowerBound.empty()) {
                        overSupportCache = std::move(build.overLowerBound);
                        build.overLowerBound.clear();
                    }
                }
                phase.initialTrackedOverKeys = std::max<uint64_t>(
                    phase.initialOverTauKeys,
                    static_cast<uint64_t>(overSupportCache.size()));
                state = makeBucketedBandedLowState(
                    std::move(build), currentCapTau, currentWindowCapTau, currentOverExactCapTau);
            }

            if (state.lowValue.empty() && state.bufferedValue.empty()) {
                persistentOverWindow = std::move(state.overWindow);
                phase.phaseMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - phaseStart).count());
                phaseStats.push_back(phase);
                currentTau = currentWindowCapTau + 1.0;
                lastTau = currentTau;
                state = BucketedBandedLowState{};
                continue;
            }

            while (completedRounds < maxRounds) {
                if (state.lowValue.empty()) {
                    if (state.bufferedValue.empty()) break;
                    currentTau = std::max(currentTau, currentCapTau + 1.0);
                    if (currentTau > tauMax) break;
                    state.lowValue.swap(state.bufferedValue);
                    state.lowBuckets.swap(state.bufferedBuckets);
                    currentCapTau = currentWindowCapTau;
                    lastTau = currentTau;
                }

                const int roundMinBucket = firstNonEmptyBucket(state.lowBuckets);
                if (roundMinBucket < 0) {
                    state.lowValue.clear();
                    state.bufferedValue.clear();
                    break;
                }
                const double roundMin = static_cast<double>(roundMinBucket);
                if (roundMin > currentWindowCapTau) break;
                if (roundMin > currentTau) {
                    currentTau = roundMin;
                    lastTau = currentTau;
                    if (currentTau > tauMax) break;
                    if (currentTau > currentCapTau && !state.bufferedValue.empty()) {
                        state.lowValue.swap(state.bufferedValue);
                        state.lowBuckets.swap(state.bufferedBuckets);
                        currentCapTau = currentWindowCapTau;
                    }
                    if (currentTau > currentWindowCapTau) break;
                    continue;
                }

                maxObservedMin = std::max(maxObservedMin, roundMin);
                phase.maxRoundMin = std::max(phase.maxRoundMin, roundMin);

                std::vector<SparseCliqueKey> frontierKeys;
                frontierKeys.reserve(state.lowBuckets[static_cast<size_t>(roundMinBucket)].size());
                for (const auto &key : state.lowBuckets[static_cast<size_t>(roundMinBucket)]) {
                    frontierKeys.push_back(key);
                }
                if (frontierKeys.empty()) break;
                phase.frontierKeys += frontierKeys.size();
                for (const auto &key : frontierKeys) {
                    state.lowValue.erase(key);
                }
                state.lowBuckets[static_cast<size_t>(roundMinBucket)].clear();

                std::vector<daf::Size> changedLeafIndex(activeLeaves.size(), std::numeric_limits<daf::Size>::max());
                std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
                std::vector<daf::Size> changedLeaf;
                removedKeyForLeaf.reserve(activeLeaves.size() / 16 + 1);
                changedLeaf.reserve(activeLeaves.size() / 16 + 1);

                daf::StaticVector<daf::Size> rmVerts;
                rmVerts.resize(r);
                for (const auto &rmKey : frontierKeys) {
                    for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
                    daf::intersect_dense_sets_multi(rmVerts, activeLeafByVertex,
                        [&](const daf::Size &leafId) {
                            if (leafId >= activeLeaves.size() || !leafAlive[leafId]) return;
                            auto &leafIdx = changedLeafIndex[leafId];
                            if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                                leafIdx = removedKeyForLeaf.size();
                                removedKeyForLeaf.emplace_back();
                                changedLeaf.push_back(leafId);
                            }
                            removedKeyForLeaf[leafIdx].push_back(rmKey);
                        });
                }
                rmVerts.free();

                std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> localDelta;
                localDelta.reserve(1 << 20);
                uint64_t roundOldEntries = 0;
                uint64_t roundNewEntries = 0;
                uint64_t roundGeneratedSubleaves = 0;
                uint64_t roundExactLeaves = 0;
                uint64_t roundCandidateKeys = 0;
                const auto leafUpdateStart = std::chrono::high_resolution_clock::now();

                for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
                    const daf::Size leafId = changedLeaf[idx];
                    auto &leaf = activeLeaves[leafId];
                    if (!leafAlive[leafId] || leaf.size() < r) continue;
                    const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
                    if (removed.empty()) continue;
                    roundExactLeaves++;
                    phase.exactLeaves++;

                    auto newLeaves = splitLeafByRemovedKeys(
                        leaf, removed, r, s, localDelta, &roundOldEntries, &roundNewEntries);

                    removeLeafFromVertexIndex(leaf, leafId, activeLeafByVertex, r);
                    leafAlive[leafId] = 0;
                    roundGeneratedSubleaves += newLeaves.size();
                    totalSpawnedSubleaves += newLeaves.size();
                    phase.spawnedSubleaves += newLeaves.size();

                    for (auto &newLeaf : newLeaves) {
                        daf::Size newLeafId = activeLeaves.size();
                        activeLeaves.push_back(std::move(newLeaf));
                        leafAlive.push_back(1);
                        addLeafToVertexIndex(activeLeaves.back(), newLeafId, activeLeafByVertex, r);
                    }
                }
                phase.leafUpdateMs += static_cast<uint64_t>(
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - leafUpdateStart).count());

                robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> frontierSet;
                frontierSet.reserve(frontierKeys.size() * 2 + 1);
                for (const auto &key : frontierKeys) frontierSet.insert(key);
                std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> recomputedTouched;
                recomputedTouched.reserve(localDelta.size() / 2 + 1);
                const auto recomputeStart = std::chrono::high_resolution_clock::now();
                const bool useDeltaFast = std::getenv("PIVOTER_QUOTIENT_BUCKETED_DELTA_FAST") != nullptr;
                const bool useDeltaPrune = std::getenv("PIVOTER_QUOTIENT_BUCKETED_DELTA_PRUNE") != nullptr;
                const bool useLeafBound = std::getenv("PIVOTER_QUOTIENT_BUCKETED_DELTA_LEAF_BOUND") != nullptr;
                const bool useOverCache = std::getenv("PIVOTER_QUOTIENT_BUCKETED_OVER_CACHE_OFF") == nullptr;
                for (const auto &[key, delta] : localDelta) {
                    if (frontierSet.find(key) != frontierSet.end()) continue;
                    roundCandidateKeys++;
                    if (std::abs(delta) <= 1e-12) {
                        phase.zeroDeltaKeys++;
                        continue;
                    }
                    if (useOverCache) {
                        auto cacheIt = overSupportCache.find(key);
                        if (cacheIt != overSupportCache.end()) {
                            const double newSupport = cacheIt->second + delta;
                            overSupportCache.erase(cacheIt);
                            phase.overCacheHits++;
                            if (newSupport > currentOverExactCapTau) {
                                state.overWindow.insert(key);
                                if (useOverLowerBound) state.overLowerBound[key] = newSupport;
                                overSupportCache.emplace(key, newSupport);
                                phase.overCacheStores++;
                                continue;
                            }
                            if (newSupport > currentWindowCapTau) {
                                const int bucket = supportBucketIndex(newSupport, currentOverExactCapTau);
                                if (bucket > 0) {
                                    state.overValue.emplace(key, static_cast<uint16_t>(bucket));
                                }
                                continue;
                            }
                            if (newSupport > 0.0) {
                                recomputedTouched.emplace(key, newSupport);
                                continue;
                            }
                        }
                    }
                    auto oldOverExactIt = state.overValue.find(key);
                    if (oldOverExactIt != state.overValue.end()) {
                        const double newSupport = static_cast<double>(oldOverExactIt->second) + delta;
                        state.overValue.erase(oldOverExactIt);
                        phase.deltaStateKeys++;
                        if (newSupport > currentOverExactCapTau) {
                            state.overWindow.insert(key);
                            if (useOverLowerBound) state.overLowerBound[key] = newSupport;
                            if (useOverCache && newSupport > currentWindowCapTau) {
                                overSupportCache[key] = newSupport;
                                phase.overCacheStores++;
                            }
                            continue;
                        }
                        if (newSupport > currentWindowCapTau) {
                            const int bucket = supportBucketIndex(newSupport, currentOverExactCapTau);
                            if (bucket > 0) {
                                state.overValue.emplace(key, static_cast<uint16_t>(bucket));
                            }
                            continue;
                        }
                        if (newSupport > 0.0) {
                            recomputedTouched.emplace(key, newSupport);
                        }
                        continue;
                    }
                    if (useOverLowerBound) {
                        auto oldOverLowerIt = state.overLowerBound.find(key);
                        if (oldOverLowerIt != state.overLowerBound.end()) {
                            const double newLowerBound = oldOverLowerIt->second + delta;
                            if (newLowerBound > currentWindowCapTau) {
                                oldOverLowerIt->second = newLowerBound;
                                phase.deltaOverLowerKeys++;
                                continue;
                            }
                            state.overLowerBound.erase(oldOverLowerIt);
                        }
                    }
                    auto overIt = state.overWindow.find(key);
                    if (overIt != state.overWindow.end() && delta > 0.0) {
                        phase.deltaOverKeys++;
                        continue;
                    }
                    if (overIt != state.overWindow.end()) {
                        if (useLeafBound) {
                            const double leafLowerBound = static_cast<double>(activeLeafIncidenceForKey(
                                key, r, activeLeaves, leafAlive, activeLeafByVertex));
                            if (leafLowerBound + delta > currentWindowCapTau) {
                                phase.deltaLeafBoundKeys++;
                                continue;
                            }
                        }
                        state.overWindow.erase(overIt);
                    }
                    int oldBucket = -1;
                    auto oldLowIt = state.lowValue.find(key);
                    if (oldLowIt != state.lowValue.end()) {
                        oldBucket = static_cast<int>(oldLowIt->second);
                    } else {
                        auto oldBufIt = state.bufferedValue.find(key);
                        if (oldBufIt != state.bufferedValue.end()) {
                            oldBucket = static_cast<int>(oldBufIt->second);
                        }
                    }
                    eraseBucketedKey(key, state.lowValue, state.lowBuckets);
                    eraseBucketedKey(key, state.bufferedValue, state.bufferedBuckets);
                    if (oldBucket > 0) {
                        const double newSupport = static_cast<double>(oldBucket) + delta;
                        phase.deltaStateKeys++;
                        if (newSupport > currentOverExactCapTau) {
                            state.overWindow.insert(key);
                            if (useOverLowerBound) state.overLowerBound[key] = newSupport;
                            if (useOverCache) {
                                overSupportCache[key] = newSupport;
                                phase.overCacheStores++;
                            }
                            continue;
                        }
                        if (newSupport > currentWindowCapTau) {
                            const int bucket = supportBucketIndex(newSupport, currentOverExactCapTau);
                            if (bucket > 0) {
                                state.overValue.emplace(key, static_cast<uint16_t>(bucket));
                            }
                            continue;
                        }
                        if (newSupport > 0.0) {
                            recomputedTouched.emplace(key, newSupport);
                        }
                        continue;
                    }
                    if (useDeltaFast || useDeltaPrune) {
                        if (useDeltaFast && oldBucket > 0) {
                            const double newSupport = static_cast<double>(oldBucket) + delta;
                            if (newSupport > 0.0) {
                                const int bucket = supportBucketIndex(newSupport, currentWindowCapTau);
                                if (bucket > 0) {
                                    recomputedTouched.emplace(key, newSupport);
                                    phase.deltaFastKeys++;
                                    continue;
                                }
                            }
                        }
                        if (useDeltaPrune && oldBucket < 0 && delta >= -1e-9) {
                            phase.deltaPrunedKeys++;
                            continue;
                        }
                    }
                    // Key not in any cache or state. With build-time cache seeding,
                    // this means the key had zero support at build time (filtered
                    // by hasPositiveContribution). Its support hasn't changed since
                    // then (no cache/state entry means no prior touch). So:
                    //   new_support = 0 + delta = delta
                    phase.fullRecomputeKeys++;
                    const double support = delta;
                    if (support > 0.0) {
                        if (useOverCache && support > currentWindowCapTau) {
                            if (support > currentOverExactCapTau) {
                                state.overWindow.insert(key);
                                if (useOverLowerBound) state.overLowerBound[key] = support;
                                overSupportCache[key] = support;
                                phase.overCacheStores++;
                            } else {
                                const int bucket = supportBucketIndex(support, currentOverExactCapTau);
                                if (bucket > 0) {
                                    state.overValue.emplace(key, static_cast<uint16_t>(bucket));
                                }
                            }
                        } else if (support > currentWindowCapTau) {
                            const int bucket = supportBucketIndex(support, currentOverExactCapTau);
                            if (bucket > 0) {
                                state.overValue.emplace(key, static_cast<uint16_t>(bucket));
                            } else {
                                state.overWindow.insert(key);
                                if (useOverLowerBound) state.overLowerBound[key] = support;
                            }
                        } else {
                            state.overWindow.erase(key);
                            state.overValue.erase(key);
                            if (useOverLowerBound) state.overLowerBound.erase(key);
                            recomputedTouched.emplace(key, support);
                        }
                    }
                }
                phase.recomputeMs += static_cast<uint64_t>(
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - recomputeStart).count());

                const auto bucketUpdateStart = std::chrono::high_resolution_clock::now();
                for (const auto &[key, value] : recomputedTouched) {
                    const int bucket = supportBucketIndex(value, currentWindowCapTau);
                    if (bucket < 0) continue;
                    if (value <= currentCapTau) {
                        insertBucketedKey(key, bucket, state.lowValue, state.lowBuckets);
                    } else {
                        insertBucketedKey(key, bucket, state.bufferedValue, state.bufferedBuckets);
                    }
                }

                const int nextMinBucket = firstNonEmptyBucket(state.lowBuckets);
                const double nextMin = (nextMinBucket < 0)
                    ? std::numeric_limits<double>::infinity()
                    : static_cast<double>(nextMinBucket);
                const uint64_t nextFrontier = (nextMinBucket < 0)
                    ? 0
                    : state.lowBuckets[static_cast<size_t>(nextMinBucket)].size();
                phase.bucketUpdateMs += static_cast<uint64_t>(
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::high_resolution_clock::now() - bucketUpdateStart).count());

                completedRounds++;
                phase.rounds++;
                totalExactLeaves += roundExactLeaves;
                totalCandidateKeys += roundCandidateKeys;
                phase.candidateKeys += roundCandidateKeys;

                std::cout << "---- Quotient bucketed-banded-low round " << completedRounds << " ----" << std::endl;
                std::cout << "  Tau/cap/window:           " << currentTau << " / "
                          << currentCapTau << " / " << currentWindowCapTau << std::endl;
                std::cout << "  Round min:                " << roundMin << std::endl;
                std::cout << "  Frontier keys:            " << frontierKeys.size() << std::endl;
                std::cout << "  Changed leaves:           " << changedLeaf.size() << std::endl;
                std::cout << "  Exact-handled leaves:     " << roundExactLeaves << std::endl;
                std::cout << "  Generated subleaves:      " << roundGeneratedSubleaves << std::endl;
                std::cout << "  Candidate touched keys:   " << roundCandidateKeys << std::endl;
                std::cout << "  Recomputed touched keys:  " << recomputedTouched.size() << std::endl;
                std::cout << "  Old local entries:        " << roundOldEntries << std::endl;
                std::cout << "  New local entries:        " << roundNewEntries << std::endl;
                std::cout << "  Next low keys:            " << state.lowValue.size() << std::endl;
                std::cout << "  Next buffered keys:       " << state.bufferedValue.size() << std::endl;
                std::cout << "  Next frontier min:        "
                          << (std::isfinite(nextMin) ? nextMin : -1.0) << std::endl;
                std::cout << "  Next frontier size:       " << nextFrontier << std::endl;

                lastTau = currentTau;
            }

            phase.remainingLowKeys = state.lowValue.size() + state.bufferedValue.size();
            phase.phaseMs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - phaseStart).count());
            if (autoLookahead && !phase.fullBand && phase.lookahead == 1) {
                const uint64_t phaseWorkMs =
                    phase.leafUpdateMs + phase.recomputeMs + phase.bucketUpdateMs;
                const double overRatio = phase.initialLowKeys == 0
                    ? (phase.initialTrackedOverKeys > 0 ? std::numeric_limits<double>::infinity() : 0.0)
                    : static_cast<double>(phase.initialTrackedOverKeys) /
                        static_cast<double>(phase.initialLowKeys);
                const bool hasRemainingTau =
                    currentWindowCapTau + static_cast<double>(autoFullBandMinRemainingTau) <= tauMax + 1e-9;
                if (hasRemainingTau &&
                    overRatio >= autoFullBandOverRatio &&
                    phase.rebuildMs >= static_cast<uint64_t>(
                        autoFullBandRebuildFactor * std::max<uint64_t>(1, phaseWorkMs))) {
                    const double nextTau = currentWindowCapTau + 1.0;
                    const int remainingLookahead = std::max(
                        1, static_cast<int>(std::ceil(tauMax - nextTau + 1.0)));
                    const int suggestedLookahead = chooseLookaheadFromOverSupportCache(
                        overSupportCache, nextTau, tauMax, autoThresholdKeys);
                    if (suggestedLookahead > 1 && suggestedLookahead < remainingLookahead) {
                        autoSuggestedLookaheadNext = suggestedLookahead;
                        autoForceFullBandNext = false;
                    } else {
                        autoSuggestedLookaheadNext = 0;
                        autoForceFullBandNext = true;
                    }
                }
            }
            phaseStats.push_back(phase);
        }

        currentTau = currentWindowCapTau + 1.0;
        lastTau = currentTau;
        persistentOverWindow = std::move(state.overWindow);
        state = BucketedBandedLowState{};
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();
    std::cout << "==== Quotient Bucketed Banded Buffered Adaptive Low-Support ====" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Tau start/max:            " << tauStart << " / " << tauMax << std::endl;
    std::cout << "  Lookahead:                " << lookahead << std::endl;
        std::cout << "  Auto lookahead:           " << (autoLookahead ? "YES" : "NO") << std::endl;
        std::cout << "  Over-support cache:       "
                  << (std::getenv("PIVOTER_QUOTIENT_BUCKETED_OVER_CACHE_OFF") == nullptr ? "YES" : "NO")
                  << std::endl;
    if (autoLookahead) {
        std::cout << "  Auto threshold keys:      " << autoThresholdKeys << std::endl;
        std::cout << "  Auto raw threshold:       " << autoRawCliqueThreshold << std::endl;
        std::cout << "  Auto wide raw threshold:  " << autoWideRawCliqueThreshold << std::endl;
        std::cout << "  Auto wide explicit min:   " << autoWideExplicitMin << std::endl;
        std::cout << "  Auto wide explicit max:   " << autoWideExplicitMax << std::endl;
        std::cout << "  Survey explicit entries:  " << surveyExplicitEntries << std::endl;
        std::cout << "  Auto full over ratio:     " << autoFullBandOverRatio << std::endl;
        std::cout << "  Auto full rebuild factor: " << autoFullBandRebuildFactor << std::endl;
        std::cout << "  Auto full min remain:     " << autoFullBandMinRemainingTau << std::endl;
    }
    std::cout << "  Rebuild count:            " << rebuildCount << std::endl;
    std::cout << "  Completed rounds:         " << completedRounds << std::endl;
    std::cout << "  Total exact leaves:       " << totalExactLeaves << std::endl;
    std::cout << "  Total candidate keys:     " << totalCandidateKeys << std::endl;
    std::cout << "  Total spawned subleaves:  " << totalSpawnedSubleaves << std::endl;
    std::cout << "  Max observed round min:   " << maxObservedMin << std::endl;
    std::cout << "  Total rebuild time:       " << totalRebuildMs << " ms" << std::endl;
    std::cout << "  Final tau:                " << lastTau << std::endl;
    std::cout << "  Prototype time:           " << elapsedMs << " ms" << std::endl;
    for (const auto &phase : phaseStats) {
        std::cout << "  Phase tau=" << phase.tau
                  << " cap=" << phase.capTau
                  << " window=" << phase.windowCapTau
                  << " lookahead=" << phase.lookahead
                  << " full_band=" << (phase.fullBand ? "yes" : "no")
                  << " two_phase=" << (phase.twoPhaseBuild ? "yes" : "no")
                  << " cache_reband=" << (phase.cacheRebandBuild ? "yes" : "no")
                  << " single_over_map=" << (phase.singleOverMapBuild ? "yes" : "no")
                  << " capped_low=" << (phase.cappedLowBuild ? "yes" : "no")
                  << " rebuild_ms=" << phase.rebuildMs
                  << " phase_ms=" << phase.phaseMs
                  << " leaf_ms=" << phase.leafUpdateMs
                  << " recomp_ms=" << phase.recomputeMs
                  << " bucket_ms=" << phase.bucketUpdateMs
                  << " delta_fast=" << phase.deltaFastKeys
                  << " delta_pruned=" << phase.deltaPrunedKeys
                  << " delta_state=" << phase.deltaStateKeys
                  << " delta_over=" << phase.deltaOverKeys
                  << " delta_leaf_bound=" << phase.deltaLeafBoundKeys
                  << " delta_over_lb=" << phase.deltaOverLowerKeys
                  << " full_recomp=" << phase.fullRecomputeKeys
                  << " over_cache_hit=" << phase.overCacheHits
                  << " over_cache_store=" << phase.overCacheStores
                  << " zero_delta=" << phase.zeroDeltaKeys
                  << " init_low=" << phase.initialLowKeys
                  << " init_buffered=" << phase.initialBufferedKeys
                  << " init_over=" << phase.initialOverTauKeys
                  << " init_tracked_over=" << phase.initialTrackedOverKeys
                  << " rounds=" << phase.rounds
                  << " frontier_keys=" << phase.frontierKeys
                  << " exact_leaves=" << phase.exactLeaves
                  << " cand_keys=" << phase.candidateKeys
                  << " spawned=" << phase.spawnedSubleaves
                  << " max_round_min=" << phase.maxRoundMin
                  << " remain_low=" << phase.remainingLowKeys
                  << " kept=" << phase.buildStats.keptUpdates
                  << " evicted=" << phase.buildStats.evictedUpdates
                  << " skipped=" << phase.buildStats.skippedUpdates
                  << std::endl;
    }
    std::cout << "===============================================================" << std::endl;
}

static void maybeRunOneRoundSparseExactPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_SPARSE_EXACT")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "=========== Quotient One-Round Sparse Exact ============" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip sparse exact first round" << std::endl;
        std::cout << "  Reason: global sparse support map would be too large" << std::endl;
        std::cout << "=========================================================" << std::endl;
        return;
    }

    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> supportByKey;
    supportByKey.reserve(static_cast<size_t>(std::min<uint64_t>(rawCliqueCount, 10000000ULL)));
    daf::timeCount("sparse support build (Quotient exact)", [&]() {
        for (const auto &leaf : tree.adj_list) {
            accumulateLeafSparseContribution(leaf, r, s, supportByKey, 1.0);
        }
    });
    if (supportByKey.empty()) return;

    double minCore = std::numeric_limits<double>::max();
    for (const auto &[key, value] : supportByKey) {
        (void)key;
        minCore = std::min(minCore, value);
    }

    std::vector<SparseCliqueKey> currentRemoveKeys;
    currentRemoveKeys.reserve(supportByKey.size() / 16 + 1);
    for (const auto &[key, value] : supportByKey) {
        if (value <= minCore) currentRemoveKeys.push_back(key);
    }

    std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
    std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
    std::vector<daf::Size> changedLeaf;
    removedKeyForLeaf.reserve(tree.adj_list.size() / 16 + 1);
    changedLeaf.reserve(tree.adj_list.size() / 16 + 1);

    daf::StaticVector<daf::Size> rmVerts;
    rmVerts.resize(r);
    for (const auto &rmKey : currentRemoveKeys) {
        for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
        daf::intersect_dense_sets_multi(rmVerts, treeGraphV.adj_list,
            [&](const TreeGraphNode &uClique) {
                auto &leafIdx = changedLeafIndex[uClique.v];
                if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                    leafIdx = removedKeyForLeaf.size();
                    removedKeyForLeaf.emplace_back();
                    changedLeaf.push_back(uClique.v);
                }
                removedKeyForLeaf[leafIdx].push_back(rmKey);
            });
    }
    rmVerts.free();

    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> sparseDelta;
    sparseDelta.reserve(1 << 20);
    uint64_t oldEntries = 0;
    uint64_t newEntries = 0;
    uint64_t generatedSubleaves = 0;
    uint64_t exactLeaves = 0;
    uint64_t singleRemovedLeaves = 0;
    uint64_t multiRemovedLeaves = 0;

    for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
        const daf::Size leafId = changedLeaf[idx];
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.size() < r) continue;
        const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
        if (removed.empty()) continue;
        if (removed.size() == 1) singleRemovedLeaves++;
        else multiRemovedLeaves++;
        exactLeaves++;
        accumulateExactLeafSparseDeltaBK(
            leaf, removed, r, s, sparseDelta, oldEntries, newEntries, generatedSubleaves);
    }

    robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> removedSet;
    removedSet.reserve(currentRemoveKeys.size() * 2 + 1);
    for (const auto &key : currentRemoveKeys) removedSet.insert(key);

    for (const auto &[key, delta] : sparseDelta) {
        if (removedSet.find(key) != removedSet.end()) continue;
        auto &slot = supportByKey[key];
        slot += delta;
    }
    for (auto it = supportByKey.begin(); it != supportByKey.end();) {
        if (removedSet.find(it->first) != removedSet.end() || it->second <= 0.0) {
            it = supportByKey.erase(it);
        } else {
            ++it;
        }
    }

    double nextMin = std::numeric_limits<double>::max();
    for (const auto &[key, value] : supportByKey) {
        (void)key;
        nextMin = std::min(nextMin, value);
    }
    uint64_t nextFrontier = 0;
    if (!supportByKey.empty()) {
        for (const auto &[key, value] : supportByKey) {
            (void)key;
            if (value <= nextMin) nextFrontier++;
        }
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "=========== Quotient One-Round Sparse Exact ============" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Sparse support states:     " << supportByKey.size() << std::endl;
    std::cout << "  First-round min core:      " << minCore << std::endl;
    std::cout << "  Removed sparse cliques:    " << currentRemoveKeys.size() << std::endl;
    std::cout << "  Changed leaves:            " << changedLeaf.size() << std::endl;
    std::cout << "  Exact-handled leaves:      " << exactLeaves << std::endl;
    std::cout << "  Single-removed leaves:     " << singleRemovedLeaves << std::endl;
    std::cout << "  Multi-removed leaves:      " << multiRemovedLeaves << std::endl;
    std::cout << "  Generated subleaves:       " << generatedSubleaves << std::endl;
    std::cout << "  Old local entries:         " << oldEntries << std::endl;
    std::cout << "  New local entries:         " << newEntries << std::endl;
    std::cout << "  Unique sparse deltas:      " << sparseDelta.size() << std::endl;
    if (!supportByKey.empty()) {
        std::cout << "  Surviving sparse states:   " << supportByKey.size() << std::endl;
        std::cout << "  Next frontier min core:    " << nextMin << std::endl;
        std::cout << "  Next frontier size:        " << nextFrontier << std::endl;
    } else {
        std::cout << "  Surviving sparse states:   0" << std::endl;
    }
    std::cout << "  Prototype time:            " << elapsedMs << " ms" << std::endl;
    std::cout << "=========================================================" << std::endl;
}

static void maybeRunMultiRoundSparseDeathOnlyPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_SPARSE")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "========== Quotient Multi-Round Sparse Prototype ==========" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip multi-round sparse prototype" << std::endl;
        std::cout << "  Reason: global sparse support map would be too large" << std::endl;
        std::cout << "===========================================================" << std::endl;
        return;
    }

    int maxRounds = 5;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_MAX")) {
        maxRounds = std::max(1, std::atoi(env));
    }

    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> supportByKey;
    supportByKey.reserve(static_cast<size_t>(std::min<uint64_t>(rawCliqueCount, 10000000ULL)));
    daf::timeCount("sparse support build (Quotient multi-round)", [&]() {
        for (const auto &leaf : tree.adj_list) {
            accumulateLeafSparseContribution(leaf, r, s, supportByKey, 1.0);
        }
    });
    if (supportByKey.empty()) return;

    std::vector<uint8_t> deadLeaf(tree.adj_list.size(), 0);
    int completedRounds = 0;
    bool splitEncountered = false;
    uint64_t totalDeadLeaves = 0;

    for (int round = 1; round <= maxRounds; ++round) {
        double minCore = std::numeric_limits<double>::max();
        for (const auto &[key, value] : supportByKey) {
            (void)key;
            minCore = std::min(minCore, value);
        }
        if (!std::isfinite(minCore)) break;

        std::vector<SparseCliqueKey> currentRemoveKeys;
        currentRemoveKeys.reserve(supportByKey.size() / 16 + 1);
        for (const auto &[key, value] : supportByKey) {
            if (value <= minCore) currentRemoveKeys.push_back(key);
        }
        if (currentRemoveKeys.empty()) break;

        std::vector<daf::Size> changedLeafIndex(tree.adj_list.size(), std::numeric_limits<daf::Size>::max());
        std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
        std::vector<daf::Size> changedLeaf;
        removedKeyForLeaf.reserve(tree.adj_list.size() / 16 + 1);
        changedLeaf.reserve(tree.adj_list.size() / 16 + 1);

        daf::StaticVector<daf::Size> rmVerts;
        rmVerts.resize(r);
        for (const auto &rmKey : currentRemoveKeys) {
            for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
            daf::intersect_dense_sets_multi(rmVerts, treeGraphV.adj_list,
                [&](const TreeGraphNode &uClique) {
                    if (deadLeaf[uClique.v]) return;
                    auto &leafIdx = changedLeafIndex[uClique.v];
                    if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                        leafIdx = removedKeyForLeaf.size();
                        removedKeyForLeaf.emplace_back();
                        changedLeaf.push_back(uClique.v);
                    }
                    removedKeyForLeaf[leafIdx].push_back(rmKey);
                });
        }
        rmVerts.free();

        std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> sparseDelta;
        sparseDelta.reserve(1 << 20);
        uint64_t roundOldEntries = 0;
        uint64_t roundNewEntries = 0;
        uint64_t roundSubleaves = 0;
        uint64_t roundSingleRemovedLeaves = 0;
        uint64_t roundMultiRemovedLeaves = 0;

        for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
            const daf::Size leafId = changedLeaf[idx];
            const auto &leaf = tree.adj_list[leafId];
            if (leaf.size() < r) continue;
            const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
            if (removed.empty()) continue;
            if (removed.size() == 1) roundSingleRemovedLeaves++;
            else roundMultiRemovedLeaves++;
            accumulateExactLeafSparseDeltaBK(
                leaf, removed, r, s, sparseDelta,
                roundOldEntries, roundNewEntries, roundSubleaves);
        }

        std::cout << "---- Quotient sparse round " << round << " ----" << std::endl;
        std::cout << "  Min core:                 " << minCore << std::endl;
        std::cout << "  Removed sparse cliques:   " << currentRemoveKeys.size() << std::endl;
        std::cout << "  Changed leaves:           " << changedLeaf.size() << std::endl;
        std::cout << "  Single-removed leaves:    " << roundSingleRemovedLeaves << std::endl;
        std::cout << "  Multi-removed leaves:     " << roundMultiRemovedLeaves << std::endl;
        std::cout << "  Generated subleaves:      " << roundSubleaves << std::endl;
        std::cout << "  Old local entries:        " << roundOldEntries << std::endl;
        std::cout << "  New local entries:        " << roundNewEntries << std::endl;
        std::cout << "  Unique sparse deltas:     " << sparseDelta.size() << std::endl;

        if (roundSubleaves > 0) {
            splitEncountered = true;
            break;
        }

        robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> removedSet;
        removedSet.reserve(currentRemoveKeys.size() * 2 + 1);
        for (const auto &key : currentRemoveKeys) removedSet.insert(key);

        for (const auto &[key, delta] : sparseDelta) {
            if (removedSet.find(key) != removedSet.end()) continue;
            auto &slot = supportByKey[key];
            slot += delta;
        }
        for (auto it = supportByKey.begin(); it != supportByKey.end();) {
            if (removedSet.find(it->first) != removedSet.end() || it->second <= 0.0) {
                it = supportByKey.erase(it);
            } else {
                ++it;
            }
        }
        for (auto leafId : changedLeaf) {
            if (!deadLeaf[leafId]) {
                deadLeaf[leafId] = 1;
                totalDeadLeaves++;
            }
        }

        completedRounds = round;
        if (supportByKey.empty()) break;
        double nextMin = std::numeric_limits<double>::max();
        for (const auto &[key, value] : supportByKey) {
            (void)key;
            nextMin = std::min(nextMin, value);
        }
        uint64_t nextFrontier = 0;
        for (const auto &[key, value] : supportByKey) {
            (void)key;
            if (value <= nextMin) nextFrontier++;
        }
        std::cout << "  Surviving sparse states:  " << supportByKey.size() << std::endl;
        std::cout << "  Next frontier min core:   " << nextMin << std::endl;
        std::cout << "  Next frontier size:       " << nextFrontier << std::endl;
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();
    std::cout << "========== Quotient Multi-Round Sparse Prototype ==========" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Completed rounds:         " << completedRounds << std::endl;
    std::cout << "  Split encountered:        " << (splitEncountered ? "YES" : "NO") << std::endl;
    std::cout << "  Dead leaves so far:       " << totalDeadLeaves << std::endl;
    std::cout << "  Remaining sparse states:  " << supportByKey.size() << std::endl;
    std::cout << "  Prototype time:           " << elapsedMs << " ms" << std::endl;
    std::cout << "===========================================================" << std::endl;
}

static void maybeRunMultiRoundSparseDynamicPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_DYNAMIC")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "========= Quotient Multi-Round Dynamic Prototype =========" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip dynamic sparse prototype" << std::endl;
        std::cout << "  Reason: global sparse support map would be too large" << std::endl;
        std::cout << "==========================================================" << std::endl;
        return;
    }

    int maxRounds = 5;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_MAX")) {
        maxRounds = std::max(1, std::atoi(env));
    }

    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> supportByKey;
    supportByKey.reserve(static_cast<size_t>(std::min<uint64_t>(rawCliqueCount, 10000000ULL)));
    daf::timeCount("sparse support build (Quotient dynamic)", [&]() {
        for (const auto &leaf : tree.adj_list) {
            accumulateLeafSparseContribution(leaf, r, s, supportByKey, 1.0);
        }
    });
    if (supportByKey.empty()) return;

    std::vector<std::vector<TreeGraphNode>> activeLeaves;
    activeLeaves.reserve(tree.adj_list.size() * 2 + 16);
    for (const auto &leaf : tree.adj_list) {
        activeLeaves.emplace_back(leaf.begin(), leaf.end());
    }
    std::vector<uint8_t> leafAlive(activeLeaves.size(), 1);
    std::vector<robin_hood::unordered_flat_set<daf::Size>> activeLeafByVertex(edgeGraph.n);
    for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
        addLeafToVertexIndex(activeLeaves[leafId], leafId, activeLeafByVertex, r);
    }

    int completedRounds = 0;
    uint64_t totalDeadLeaves = 0;
    uint64_t totalSpawnedLeaves = 0;
    uint64_t totalSplitLeaves = 0;
    uint64_t totalSingleChildLeaves = 0;
    uint64_t totalMultiChildLeaves = 0;
    int firstCoreGTOneRound = 0;
    int firstSplitRound = 0;
    int roundsWithSplit = 0;
    uint64_t maxSpawnedLeavesInRound = 0;
    double maxObservedCore = 0.0;

    for (int round = 1; round <= maxRounds; ++round) {
        double minCore = std::numeric_limits<double>::max();
        for (const auto &[key, value] : supportByKey) {
            (void)key;
            minCore = std::min(minCore, value);
        }
        if (!std::isfinite(minCore)) break;
        maxObservedCore = std::max(maxObservedCore, minCore);
        if (minCore > 1.0 && firstCoreGTOneRound == 0) firstCoreGTOneRound = round;

        std::vector<SparseCliqueKey> currentRemoveKeys;
        currentRemoveKeys.reserve(supportByKey.size() / 16 + 1);
        for (const auto &[key, value] : supportByKey) {
            if (value <= minCore) currentRemoveKeys.push_back(key);
        }
        if (currentRemoveKeys.empty()) break;

        std::vector<daf::Size> changedLeafIndex(activeLeaves.size(), std::numeric_limits<daf::Size>::max());
        std::vector<std::vector<SparseCliqueKey>> removedKeyForLeaf;
        std::vector<daf::Size> changedLeaf;
        removedKeyForLeaf.reserve(activeLeaves.size() / 16 + 1);
        changedLeaf.reserve(activeLeaves.size() / 16 + 1);

        daf::StaticVector<daf::Size> rmVerts;
        rmVerts.resize(r);
        for (const auto &rmKey : currentRemoveKeys) {
            for (daf::CliqueSize j = 0; j < r; ++j) rmVerts[j] = rmKey.verts[j];
            daf::intersect_dense_sets_multi(rmVerts, activeLeafByVertex,
                [&](const daf::Size &leafId) {
                    if (leafId >= activeLeaves.size() || !leafAlive[leafId]) return;
                    auto &leafIdx = changedLeafIndex[leafId];
                    if (leafIdx == std::numeric_limits<daf::Size>::max()) {
                        leafIdx = removedKeyForLeaf.size();
                        removedKeyForLeaf.emplace_back();
                        changedLeaf.push_back(leafId);
                    }
                    removedKeyForLeaf[leafIdx].push_back(rmKey);
                });
        }
        rmVerts.free();

        std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> sparseDelta;
        sparseDelta.reserve(1 << 20);
        uint64_t roundOldEntries = 0;
        uint64_t roundNewEntries = 0;
        uint64_t roundSpawnedLeaves = 0;
        uint64_t roundSingleRemovedLeaves = 0;
        uint64_t roundMultiRemovedLeaves = 0;
        uint64_t roundChangedLeaves = changedLeaf.size();
        uint64_t roundSplitLeaves = 0;
        uint64_t roundSingleChildLeaves = 0;
        uint64_t roundMultiChildLeaves = 0;

        for (auto leafId : changedLeaf) {
            if (leafId >= activeLeaves.size() || !leafAlive[leafId]) continue;
            const auto &leaf = activeLeaves[leafId];
            const auto &removed = removedKeyForLeaf[changedLeafIndex[leafId]];
            if (removed.empty()) continue;
            if (removed.size() == 1) roundSingleRemovedLeaves++;
            else roundMultiRemovedLeaves++;

            auto newLeaves = splitLeafByRemovedKeys(
                leaf, removed, r, s, sparseDelta, &roundOldEntries, &roundNewEntries);
            if (!newLeaves.empty()) {
                roundSplitLeaves++;
                totalSplitLeaves++;
                if (newLeaves.size() == 1) {
                    roundSingleChildLeaves++;
                    totalSingleChildLeaves++;
                } else {
                    roundMultiChildLeaves++;
                    totalMultiChildLeaves++;
                }
            }

            removeLeafFromVertexIndex(leaf, leafId, activeLeafByVertex, r);
            leafAlive[leafId] = 0;
            totalDeadLeaves++;

            for (auto &newLeaf : newLeaves) {
                daf::Size newLeafId = activeLeaves.size();
                activeLeaves.push_back(std::move(newLeaf));
                leafAlive.push_back(1);
                addLeafToVertexIndex(activeLeaves.back(), newLeafId, activeLeafByVertex, r);
                roundSpawnedLeaves++;
                totalSpawnedLeaves++;
            }
        }
        if (roundSplitLeaves > 0) {
            roundsWithSplit++;
            if (firstSplitRound == 0) firstSplitRound = round;
        }
        maxSpawnedLeavesInRound = std::max(maxSpawnedLeavesInRound, roundSpawnedLeaves);

        robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> removedSet;
        removedSet.reserve(currentRemoveKeys.size() * 2 + 1);
        for (const auto &key : currentRemoveKeys) removedSet.insert(key);

        for (const auto &[key, delta] : sparseDelta) {
            if (removedSet.find(key) != removedSet.end()) continue;
            auto &slot = supportByKey[key];
            slot += delta;
        }
        for (auto it = supportByKey.begin(); it != supportByKey.end();) {
            if (removedSet.find(it->first) != removedSet.end() || it->second <= 0.0) {
                it = supportByKey.erase(it);
            } else {
                ++it;
            }
        }

        completedRounds = round;

        std::cout << "---- Quotient dynamic round " << round << " ----" << std::endl;
        std::cout << "  Min core:                 " << minCore << std::endl;
        std::cout << "  Removed sparse cliques:   " << currentRemoveKeys.size() << std::endl;
        std::cout << "  Changed leaves:           " << roundChangedLeaves << std::endl;
        std::cout << "  Single-removed leaves:    " << roundSingleRemovedLeaves << std::endl;
        std::cout << "  Multi-removed leaves:     " << roundMultiRemovedLeaves << std::endl;
        std::cout << "  Split leaves:             " << roundSplitLeaves << std::endl;
        std::cout << "  Single-child leaves:      " << roundSingleChildLeaves << std::endl;
        std::cout << "  Multi-child leaves:       " << roundMultiChildLeaves << std::endl;
        std::cout << "  Old local entries:        " << roundOldEntries << std::endl;
        std::cout << "  New local entries:        " << roundNewEntries << std::endl;
        std::cout << "  Unique sparse deltas:     " << sparseDelta.size() << std::endl;
        std::cout << "  Spawned subleaves:        " << roundSpawnedLeaves << std::endl;
        std::cout << "  Surviving sparse states:  " << supportByKey.size() << std::endl;

        if (supportByKey.empty()) break;
        double nextMin = std::numeric_limits<double>::max();
        for (const auto &[key, value] : supportByKey) {
            (void)key;
            nextMin = std::min(nextMin, value);
        }
        uint64_t nextFrontier = 0;
        for (const auto &[key, value] : supportByKey) {
            (void)key;
            if (value <= nextMin) nextFrontier++;
        }
        std::cout << "  Next frontier min core:   " << nextMin << std::endl;
        std::cout << "  Next frontier size:       " << nextFrontier << std::endl;
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();
    std::cout << "========= Quotient Multi-Round Dynamic Prototype =========" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Completed rounds:         " << completedRounds << std::endl;
    std::cout << "  First core>1 round:       " << firstCoreGTOneRound << std::endl;
    std::cout << "  First split round:        " << firstSplitRound << std::endl;
    std::cout << "  Rounds with split:        " << roundsWithSplit << std::endl;
    std::cout << "  Dead leaves so far:       " << totalDeadLeaves << std::endl;
    std::cout << "  Split leaves so far:      " << totalSplitLeaves << std::endl;
    std::cout << "  Single-child leaves:      " << totalSingleChildLeaves << std::endl;
    std::cout << "  Multi-child leaves:       " << totalMultiChildLeaves << std::endl;
    std::cout << "  Spawned leaves so far:    " << totalSpawnedLeaves << std::endl;
    std::cout << "  Max spawned in one round: " << maxSpawnedLeavesInRound << std::endl;
    std::cout << "  Max observed core:        " << maxObservedCore << std::endl;
    std::cout << "  Remaining sparse states:  " << supportByKey.size() << std::endl;
    std::cout << "  Prototype time:           " << elapsedMs << " ms" << std::endl;
    std::cout << "==========================================================" << std::endl;
}

// ================================================================
// Class-level peeling prototype for r=3.
//
// Instead of 338M individual keys, tracks ~848K quotient classes.
#include "NucleusCoreDecompositionRCliqueST_QuotientLab_ExperimentalPeeling.inc"

} // namespace

std::vector<std::pair<std::vector<daf::Size>, double>> NucleusCoreDecompositionRClique_ST_QuotientLab(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    if (!std::getenv("PIVOTER_QUOTIENT_SKIP_SURVEY")) {
        printQuotientStats(tree, r, s);
    }
    if (std::getenv("PIVOTER_QUOTIENT_VERTEX_ANALYSIS")) {
        std::vector<uint32_t> vtxLeafCount(edgeGraph.n, 0);
        std::vector<uint32_t> vtxPivotCount(edgeGraph.n, 0);
        std::vector<uint32_t> vtxKeepCount(edgeGraph.n, 0);
        for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
            if (tree.adj_list[leafId].size() < r) continue;
            for (const auto &node : tree.adj_list[leafId]) {
                if (node.v < edgeGraph.n) {
                    vtxLeafCount[node.v]++;
                    if (node.isPivot) vtxPivotCount[node.v]++;
                    else vtxKeepCount[node.v]++;
                }
            }
        }
        std::map<uint32_t, uint32_t> pivotDist;
        uint32_t pivotRare1 = 0, pivotRare5 = 0, totalPivotActive = 0;
        for (daf::Size v = 0; v < edgeGraph.n; ++v) {
            if (vtxPivotCount[v] > 0) {
                pivotDist[vtxPivotCount[v]]++;
                totalPivotActive++;
                if (vtxPivotCount[v] == 1) pivotRare1++;
                if (vtxPivotCount[v] <= 5) pivotRare5++;
            }
        }
        std::cout << "=== Vertex-to-leaf analysis ===" << std::endl;
        std::cout << "  Pivot-active vertices: " << totalPivotActive << std::endl;
        std::cout << "  pivotLeafCount=1: " << pivotRare1 << std::endl;
        std::cout << "  pivotLeafCount<=5: " << pivotRare5 << std::endl;
        int shown = 0;
        for (auto &[c, n] : pivotDist) {
            if (shown < 10) std::cout << "  pivotCount=" << c << " -> " << n << " vertices" << std::endl;
            shown++;
        }
        std::cout << "===============================" << std::endl;
    }
    if (std::getenv("PIVOTER_QUOTIENT_VERTEX_SUPPORT")) {
        // Vertex-level support analysis for r=3.
        // For a triangle K=(va,vb,vc) with degeneracy positions pos(va)<pos(vb)<pos(vc):
        //   support(K) = nCr[pivotC_a - 2][s-3] + Σ_{u: pos(u)<pos(va), u adj all 3} nCr[pivotC_u - 3][s-4]
        // For s=4: support = (pivotC_a - 2) + |{u: pos(u)<pos(va), u adj all 3}|
        // This is computable from the graph structure without leaf enumeration.
        auto tStart = std::chrono::high_resolution_clock::now();

        // Build per-vertex info from tree
        struct VtxInfo { daf::Size pivotC = 0; daf::Size keepLeafId = 0; bool hasKeepLeaf = false; };
        std::vector<VtxInfo> vtxInfo(edgeGraph.n);
        for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
            const auto &leaf = tree.adj_list[leafId];
            if (leaf.size() < r) continue;
            daf::Size pivotC = 0;
            daf::Size keepV = 0;
            for (const auto &node : leaf) {
                if (node.isPivot) pivotC++;
                else keepV = node.v;
            }
            for (const auto &node : leaf) {
                if (!node.isPivot && node.v < edgeGraph.n) {
                    vtxInfo[node.v].pivotC = pivotC;
                    vtxInfo[node.v].keepLeafId = leafId;
                    vtxInfo[node.v].hasKeepLeaf = true;
                }
            }
        }

        // Find vertices with smallest pivotC (smallest keep-leaf → lowest support triangles)
        struct Candidate { daf::Size v; daf::Size pivotC; daf::Size keepLeafId; };
        std::vector<Candidate> candidates;
        for (daf::Size v = 0; v < edgeGraph.n; ++v) {
            if (vtxInfo[v].hasKeepLeaf && vtxInfo[v].pivotC >= r) {
                candidates.push_back({v, vtxInfo[v].pivotC, vtxInfo[v].keepLeafId});
            }
        }
        std::sort(candidates.begin(), candidates.end(), [](const Candidate &a, const Candidate &b) {
            return a.pivotC < b.pivotC;
        });

        // For s=4: keep-leaf contribution for class p=2 is pivotC-2
        // For general s: nCr[pivotC-2][s-3]
        double minSupport = std::numeric_limits<double>::max();
        daf::Size minSupportVertex = 0;
        uint64_t minSupportTriangles = 0;
        int checked = 0;
        for (const auto &cand : candidates) {
            double keepContrib = nCr[cand.pivotC - 2][static_cast<int>(s) - 3];
            if (keepContrib >= minSupport) break; // sorted by pivotC, so no better candidate exists
            // Count triangles from this vertex's keep-leaf with min n_commonEarlier
            // The minimum n_commonEarlier = 0 (triangles with no common earlier neighbor)
            // So min support from this vertex = keepContrib + 0 = keepContrib
            if (keepContrib < minSupport) {
                minSupport = keepContrib;
                minSupportVertex = cand.v;
                // Count how many triangles have n_commonEarlier=0
                const auto &leaf = tree.adj_list[cand.keepLeafId];
                // All C(pivotC, 2) triangles in the keep-leaf have support ≥ keepContrib
                // Those with support = keepContrib have n_commonEarlier = 0
                minSupportTriangles = static_cast<uint64_t>(cand.pivotC) * (cand.pivotC - 1) / 2;
            }
            checked++;
            if (checked >= 100) break; // only check first 100 candidates
        }

        auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - tStart).count();

        std::cout << "=== Vertex-level support analysis ===" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Total vertices with keep-leaf: " << candidates.size() << std::endl;
        std::cout << "  Smallest pivotC: " << (candidates.empty() ? 0 : candidates[0].pivotC) << std::endl;
        std::cout << "  Largest pivotC: " << (candidates.empty() ? 0 : candidates.back().pivotC) << std::endl;
        std::cout << "  Estimated min support: " << minSupport << std::endl;
        std::cout << "  Min support vertex: " << minSupportVertex << std::endl;
        std::cout << "  Triangles at min support (upper bound): " << minSupportTriangles << std::endl;
        std::cout << "  Analysis time: " << elapsedMs << " ms" << std::endl;
        // Show distribution of pivotC
        int shown = 0;
        for (const auto &c : candidates) {
            if (shown < 10) std::cout << "  pivotC=" << c.pivotC << " vertex=" << c.v << std::endl;
            shown++;
        }
        std::cout << "=====================================" << std::endl;
    }
    if (std::getenv("PIVOTER_QUOTIENT_TWO_PHASE_BUILD")) {
        const auto rawCliqueCount = tree.cliqueCount(r);
        std::vector<std::vector<TreeGraphNode>> activeLeaves;
        activeLeaves.reserve(tree.adj_list.size());
        for (const auto &leaf : tree.adj_list) {
            activeLeaves.emplace_back(leaf.begin(), leaf.end());
        }
        std::vector<uint8_t> leafAlive(activeLeaves.size(), 1);
        std::vector<robin_hood::unordered_flat_set<daf::Size>> activeLeafByVertex(edgeGraph.n);
        for (daf::Size leafId = 0; leafId < activeLeaves.size(); ++leafId) {
            addLeafToVertexIndex(activeLeaves[leafId], leafId, activeLeafByVertex, r);
        }
        double tau = 1.0;
        if (const char *env = std::getenv("PIVOTER_QUOTIENT_TWO_PHASE_TAU"))
            tau = std::max(1.0, std::atof(env));
        auto buildResult = daf::timeCount("TwoPhase low-support build", [&]() {
            return buildLowSupportTwoPhase(activeLeaves, leafAlive, activeLeafByVertex,
                                           r, s, tau, rawCliqueCount);
        });
        std::cout << "  TwoPhase result: lowKeys=" << buildResult.lowSupport.size()
                  << " overKeys=" << buildResult.overTau.size() << std::endl;
    }
    maybeCompareIndexCoverage(tree, edgeGraph, r, s);
    maybeBuildPrototypeState(tree, r, s);
    maybeRunOneRoundPrototype(tree, edgeGraph, treeGraphV, r, s, prebuiltIndex);
    maybeRunSparseSupportSharingSurvey(tree, r, s);
    maybeRunFirstFrontierOwnershipAnalysis(tree, treeGraphV, r, s);
    maybeRunPartialClassAnalysis(tree, treeGraphV, r, s);
    maybeRunLowSupportPrototype(tree, r, s);
    maybeRunFirstFrontierLowSupportPrototype(tree, treeGraphV, r, s);
    maybeRunFirstFrontierNextPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunMultiRoundLowSupportPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunAdaptiveLowSupportPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunBufferedAdaptiveLowSupportPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunBandedBufferedAdaptiveLowSupportPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunBucketedBandedBufferedAdaptiveLowSupportPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunOneRoundSparseSupportPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunOneRoundSparseExactPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunMultiRoundSparseDeathOnlyPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunMultiRoundSparseDynamicPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunClassLevelPeelingPrototype(tree, edgeGraph, treeGraphV, r, s);
    maybeRunVertexLevelPeelingPrototype(tree, edgeGraph, treeGraphV, r, s);
    if (const char *only = std::getenv("PIVOTER_QUOTIENT_LAB_ONLY");
        only != nullptr && std::string(only) == "1") {
        std::cout << "QuotientLab-only mode: skip exact decomposition." << std::endl;
        return {};
    }
    return NucleusCoreDecompositionRClique_ST_V20(
        tree, edgeGraph, treeGraphV, r, s, prebuiltIndex);
}
