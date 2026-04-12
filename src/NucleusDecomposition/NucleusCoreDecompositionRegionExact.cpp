#include "NCliqueCoreDecomposition.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

extern double nCr[1001][401];
extern std::vector<bool> g_maxCliqueTags;

namespace {

using TupleKey = std::vector<daf::Size>;

struct TupleHash {
    size_t operator()(const TupleKey &t) const noexcept {
        size_t h = t.size();
        for (auto x : t) {
            h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

struct ClassInfo {
    std::vector<daf::Size> regionIds;
    std::vector<daf::Size> vertices;
    daf::Size size = 0;
};

struct RTupleInfo {
    TupleKey key;
    std::vector<std::pair<daf::Size, int>> counts;
    long double multiplicity = 0.0L;
    long double support = 0.0L;
    std::vector<std::pair<daf::Size, long double>> incidentSTuples;
};

struct STupleInfo {
    TupleKey key;
    std::vector<std::pair<daf::Size, int>> counts;
    long double multiplicity = 0.0L;
    bool alive = true;
    std::vector<std::pair<daf::Size, long double>> incidentRTuples;
};

static long double chooseCount(daf::Size n, daf::Size k) {
    if (k > n) return 0.0L;
    if (k == 0 || k == n) return 1.0L;
    if (n <= 1000 && k <= 400) return static_cast<long double>(nCr[n][k]);
    return static_cast<long double>(binom_u128(n, k));
}

static std::vector<std::pair<daf::Size, int>> compressCounts(const TupleKey &key) {
    std::vector<std::pair<daf::Size, int>> counts;
    for (daf::Size cid : key) {
        if (counts.empty() || counts.back().first != cid) {
            counts.emplace_back(cid, 1);
        } else {
            counts.back().second++;
        }
    }
    return counts;
}

static long double tupleMultiplicity(
    const std::vector<std::pair<daf::Size, int>> &counts,
    const std::vector<ClassInfo> &classes) {

    long double mult = 1.0L;
    for (const auto &[cid, cnt] : counts) {
        mult *= chooseCount(classes[cid].size, static_cast<daf::Size>(cnt));
    }
    return mult;
}

static long double extensionCount(
    const std::vector<std::pair<daf::Size, int>> &hCounts,
    const std::vector<std::pair<daf::Size, int>> &tCounts,
    const std::vector<ClassInfo> &classes) {

    long double ways = 1.0L;
    size_t i = 0;
    size_t j = 0;
    while (i < hCounts.size()) {
        const auto &[hCid, hCnt] = hCounts[i];
        int tCnt = 0;
        if (j < tCounts.size() && tCounts[j].first == hCid) {
            tCnt = tCounts[j].second;
            ++j;
        }
        if (tCnt > hCnt) return 0.0L;
        ways *= chooseCount(classes[hCid].size - static_cast<daf::Size>(tCnt),
                            static_cast<daf::Size>(hCnt - tCnt));
        ++i;
    }
    if (j != tCounts.size()) return 0.0L;
    return ways;
}

template<typename EmitFunc>
static void enumerateRegionTuples(
    const std::vector<daf::Size> &classIds,
    const std::vector<ClassInfo> &classes,
    daf::CliqueSize k,
    EmitFunc &&emit) {

    if (classIds.empty()) return;
    std::vector<int> localCounts(classIds.size(), 0);
    TupleKey current;
    current.reserve(k);

    std::function<void(int, int)> dfs = [&](int start, int remaining) {
        if (remaining == 0) {
            emit(current);
            return;
        }
        for (int i = start; i < static_cast<int>(classIds.size()); ++i) {
            const daf::Size cid = classIds[i];
            if (localCounts[i] >= static_cast<int>(classes[cid].size)) continue;
            current.push_back(cid);
            localCounts[i]++;
            dfs(i, remaining - 1);
            localCounts[i]--;
            current.pop_back();
        }
    };

    dfs(0, static_cast<int>(k));
}

static std::vector<TupleKey> enumerateUniqueSubTuples(const TupleKey &key, daf::CliqueSize r) {
    std::unordered_set<TupleKey, TupleHash> subSet;
    TupleKey current;
    current.reserve(r);

    std::function<void(int, int)> dfs = [&](int start, int remaining) {
        if (remaining == 0) {
            subSet.insert(current);
            return;
        }
        for (int i = start; i <= static_cast<int>(key.size()) - remaining; ++i) {
            current.push_back(key[i]);
            dfs(i + 1, remaining - 1);
            current.pop_back();
        }
    };

    dfs(0, static_cast<int>(r));
    std::vector<TupleKey> result;
    result.reserve(subSet.size());
    for (const auto &sub : subSet) result.push_back(sub);
    return result;
}

static double normalizeCore(long double value) {
    long double rounded = std::round(value);
    if (std::fabsl(value - rounded) < 1e-9L) {
        return static_cast<double>(rounded);
    }
    return static_cast<double>(value);
}

template<typename EmitFunc>
static void enumerateTupleInstances(
    const std::vector<std::pair<daf::Size, int>> &counts,
    const std::vector<ClassInfo> &classes,
    EmitFunc &&emit) {

    std::vector<daf::Size> current;
    size_t totalSize = 0;
    for (const auto &[_, cnt] : counts) totalSize += static_cast<size_t>(cnt);
    current.reserve(totalSize);

    std::function<void(size_t)> dfsClass = [&](size_t idx) {
        if (idx == counts.size()) {
            std::vector<daf::Size> clique = current;
            std::sort(clique.begin(), clique.end());
            emit(clique);
            return;
        }

        const auto &[cid, need] = counts[idx];
        const auto &verts = classes[cid].vertices;
        std::vector<daf::Size> picked;
        picked.reserve(static_cast<size_t>(need));

        std::function<void(size_t, int)> dfsPick = [&](size_t start, int remain) {
            if (remain == 0) {
                const size_t oldSize = current.size();
                current.insert(current.end(), picked.begin(), picked.end());
                dfsClass(idx + 1);
                current.resize(oldSize);
                return;
            }
            if (start >= verts.size()) return;
            for (size_t i = start; i <= verts.size() - static_cast<size_t>(remain); ++i) {
                picked.push_back(verts[i]);
                dfsPick(i + 1, remain - 1);
                picked.pop_back();
            }
        };

        dfsPick(0, need);
    };

    dfsClass(0);
}

} // namespace

std::vector<std::pair<std::vector<daf::Size>, double>>
NucleusCoreDecompositionRClique_RegionExact(
    DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s,
    StaticCliqueIndex *prebuiltIndex) {

    (void)treeGraphV;
    auto timeStart = std::chrono::high_resolution_clock::now();
    auto step0 = timeStart;
    const daf::Size numVertices = edgeGraph.getGraphNodeSize();
    const daf::Size numPaths = tree.adj_list.size();
    const daf::Size INVALID = static_cast<daf::Size>(-1);

    std::cout << "======= Region Exact Nucleus Decomposition =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;

    // Step 1: maximal-clique regions.
    std::vector<bool> pathValid(numPaths, false);
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        pathValid[pid] = tree.adj_list[pid].size() >= s;
    }

    std::vector<daf::Size> regionOf(numPaths, INVALID);
    std::vector<daf::Size> maximalPathIds;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (!pathValid[pid]) continue;
        const bool isMax = (pid < g_maxCliqueTags.size()) ? g_maxCliqueTags[pid] : true;
        if (isMax) {
            regionOf[pid] = maximalPathIds.size();
            maximalPathIds.push_back(pid);
        }
    }

    std::vector<std::vector<daf::Size>> regionVerts(maximalPathIds.size());
    std::vector<std::vector<daf::Size>> vertexRegions(numVertices);
    for (daf::Size rid = 0; rid < maximalPathIds.size(); ++rid) {
        const auto &leaf = tree.adj_list[maximalPathIds[rid]];
        auto &verts = regionVerts[rid];
        verts.reserve(leaf.size());
        for (const auto &node : leaf) verts.push_back(node.v);
        std::sort(verts.begin(), verts.end());
        for (daf::Size v : verts) {
            vertexRegions[v].push_back(rid);
        }
    }

    daf::Size orphanRegions = 0;
    for (daf::Size pid = 0; pid < numPaths; ++pid) {
        if (!pathValid[pid] || regionOf[pid] != INVALID) continue;

        std::vector<daf::Size> pathVerts;
        pathVerts.reserve(tree.adj_list[pid].size());
        for (const auto &node : tree.adj_list[pid]) pathVerts.push_back(node.v);
        std::sort(pathVerts.begin(), pathVerts.end());

        daf::Size rareV = pathVerts.front();
        size_t rareCount = vertexRegions[rareV].size();
        for (daf::Size v : pathVerts) {
            if (vertexRegions[v].size() < rareCount) {
                rareV = v;
                rareCount = vertexRegions[v].size();
            }
        }

        for (daf::Size rid : vertexRegions[rareV]) {
            const auto &cand = regionVerts[rid];
            if (cand.size() < pathVerts.size()) continue;
            if (std::includes(cand.begin(), cand.end(), pathVerts.begin(), pathVerts.end())) {
                regionOf[pid] = rid;
                break;
            }
        }

        if (regionOf[pid] == INVALID) {
            regionOf[pid] = regionVerts.size();
            regionVerts.push_back(pathVerts);
            for (daf::Size v : pathVerts) vertexRegions[v].push_back(regionOf[pid]);
            orphanRegions++;
        }
    }
    auto step1 = std::chrono::high_resolution_clock::now();

    // Step 2: overlap classes.
    struct ProfileHash {
        size_t operator()(const std::vector<daf::Size> &p) const noexcept {
            size_t h = p.size();
            for (auto x : p) {
                h ^= std::hash<daf::Size>()(x) + 0x9e3779b9ULL + (h << 6) + (h >> 2);
            }
            return h;
        }
    };

    std::vector<ClassInfo> classes;
    std::vector<daf::Size> classOf(numVertices, INVALID);
    std::unordered_map<std::vector<daf::Size>, daf::Size, ProfileHash> profileToClass;
    for (daf::Size v = 0; v < numVertices; ++v) {
        auto &profile = vertexRegions[v];
        if (profile.empty()) continue;
        std::sort(profile.begin(), profile.end());
        auto it = profileToClass.find(profile);
        if (it == profileToClass.end()) {
            const daf::Size cid = classes.size();
            profileToClass.emplace(profile, cid);
            classes.push_back({profile, {v}, 1});
            classOf[v] = cid;
        } else {
            classes[it->second].size++;
            classes[it->second].vertices.push_back(v);
            classOf[v] = it->second;
        }
    }

    std::vector<std::vector<daf::Size>> classesInRegion(regionVerts.size());
    for (daf::Size cid = 0; cid < classes.size(); ++cid) {
        for (daf::Size rid : classes[cid].regionIds) {
            classesInRegion[rid].push_back(cid);
        }
    }
    for (auto &cids : classesInRegion) {
        std::sort(cids.begin(), cids.end());
        cids.erase(std::unique(cids.begin(), cids.end()), cids.end());
    }
    auto step2 = std::chrono::high_resolution_clock::now();

    // Step 3: exact r-class / s-class tuple sets.
    std::unordered_set<TupleKey, TupleHash> rTupleSet;
    std::unordered_set<TupleKey, TupleHash> sTupleSet;
    daf::Size maxClassesPerRegion = 0;
    for (daf::Size rid = 0; rid < classesInRegion.size(); ++rid) {
        const auto &cids = classesInRegion[rid];
        maxClassesPerRegion = std::max(maxClassesPerRegion, static_cast<daf::Size>(cids.size()));
        enumerateRegionTuples(cids, classes, r, [&](const TupleKey &key) {
            rTupleSet.insert(key);
        });
        enumerateRegionTuples(cids, classes, s, [&](const TupleKey &key) {
            sTupleSet.insert(key);
        });
    }

    std::vector<RTupleInfo> rTuples;
    rTuples.reserve(rTupleSet.size());
    std::unordered_map<TupleKey, daf::Size, TupleHash> rTupleIndex;
    for (const auto &key : rTupleSet) {
        RTupleInfo info;
        info.key = key;
        info.counts = compressCounts(key);
        info.multiplicity = tupleMultiplicity(info.counts, classes);
        const daf::Size tid = rTuples.size();
        rTupleIndex.emplace(info.key, tid);
        rTuples.push_back(std::move(info));
    }

    std::vector<STupleInfo> sTuples;
    sTuples.reserve(sTupleSet.size());
    for (const auto &key : sTupleSet) {
        STupleInfo info;
        info.key = key;
        info.counts = compressCounts(key);
        info.multiplicity = tupleMultiplicity(info.counts, classes);
        sTuples.push_back(std::move(info));
    }
    auto step3 = std::chrono::high_resolution_clock::now();

    // Step 4: exact incidence and initial support.
    long double totalIncidenceWeight = 0.0L;
    size_t totalIncidenceEdges = 0;
    for (daf::Size sid = 0; sid < sTuples.size(); ++sid) {
        auto &sInfo = sTuples[sid];
        const auto subTuples = enumerateUniqueSubTuples(sInfo.key, r);
        for (const auto &subKey : subTuples) {
            auto it = rTupleIndex.find(subKey);
            if (it == rTupleIndex.end()) continue;
            const daf::Size tid = it->second;
            const long double ext = extensionCount(
                sInfo.counts, rTuples[tid].counts, classes);
            if (ext <= 0.0L) continue;
            sInfo.incidentRTuples.emplace_back(tid, ext);
            rTuples[tid].incidentSTuples.emplace_back(sid, ext);
            rTuples[tid].support += ext;
            totalIncidenceWeight += ext;
            totalIncidenceEdges++;
        }
    }
    auto step4 = std::chrono::high_resolution_clock::now();

    long double totalSupportByClasses = 0.0L;
    long double totalRCliques = 0.0L;
    for (const auto &tInfo : rTuples) {
        totalSupportByClasses += tInfo.multiplicity * tInfo.support;
        totalRCliques += tInfo.multiplicity;
    }

    long double totalSupportCPI = 0.0L;
    for (const auto &leaf : tree.adj_list) {
        daf::Size keep = 0;
        daf::Size pivot = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pivot++;
            else keep++;
        }
        if (keep + pivot < s || pivot < s - keep) continue;
        totalSupportCPI += chooseCount(pivot, s - keep) * chooseCount(s, r);
    }

    // Step 5: exact compressed peeling.
    using HeapEntry = std::pair<long double, daf::Size>;
    std::priority_queue<HeapEntry, std::vector<HeapEntry>, std::greater<>> pq;
    std::vector<long double> coreOfRTuple(rTuples.size(), 0.0L);
    std::vector<long double> supportOfRTuple(rTuples.size(), 0.0L);
    std::vector<bool> peeledRTuple(rTuples.size(), false);

    for (daf::Size tid = 0; tid < rTuples.size(); ++tid) {
        supportOfRTuple[tid] = rTuples[tid].support;
        pq.emplace(rTuples[tid].support, tid);
    }

    daf::Size peeledCount = 0;
    while (!pq.empty()) {
        const auto [curSupport, tid] = pq.top();
        pq.pop();
        if (peeledRTuple[tid]) continue;
        if (std::fabsl(curSupport - supportOfRTuple[tid]) > 1e-12L) continue;

        peeledRTuple[tid] = true;
        coreOfRTuple[tid] = supportOfRTuple[tid];
        peeledCount++;

        for (const auto &[sid, _] : rTuples[tid].incidentSTuples) {
            auto &sInfo = sTuples[sid];
            if (!sInfo.alive) continue;
            sInfo.alive = false;
            for (const auto &[otherTid, ext] : sInfo.incidentRTuples) {
                if (peeledRTuple[otherTid]) continue;
                supportOfRTuple[otherTid] =
                    std::max(coreOfRTuple[tid], supportOfRTuple[otherTid] - ext);
                pq.emplace(supportOfRTuple[otherTid], otherTid);
            }
        }
    }
    auto step5 = std::chrono::high_resolution_clock::now();

    std::vector<std::pair<std::vector<daf::Size>, double>> result;
    if (totalRCliques <= static_cast<long double>(std::numeric_limits<size_t>::max())) {
        result.reserve(static_cast<size_t>(totalRCliques));
    }
    for (daf::Size tid = 0; tid < rTuples.size(); ++tid) {
        const double coreValue = normalizeCore(coreOfRTuple[tid]);
        enumerateTupleInstances(rTuples[tid].counts, classes,
                                [&](const std::vector<daf::Size> &clique) {
                                    result.emplace_back(clique, coreValue);
                                });
    }
    auto step6 = std::chrono::high_resolution_clock::now();

    auto timeEnd = std::chrono::high_resolution_clock::now();
    const auto totalMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
    const auto step1Ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(step1 - step0).count();
    const auto step2Ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(step2 - step1).count();
    const auto step3Ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(step3 - step2).count();
    const auto step4Ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(step4 - step3).count();
    const auto step5Ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(step5 - step4).count();
    const auto step6Ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(step6 - step5).count();

    std::cout << "  Regions: " << regionVerts.size()
              << " (orphans=" << orphanRegions << ")" << std::endl;
    std::cout << "  Overlap classes: " << classes.size()
              << ", max classes/region: " << maxClassesPerRegion << std::endl;
    std::cout << "  r-class tuples: " << rTuples.size()
              << ", s-class tuples: " << sTuples.size() << std::endl;
    std::cout << "  Incidence edges: " << totalIncidenceEdges
              << ", incidence weight sum: " << totalIncidenceWeight << std::endl;
    std::cout << "  Total r-cliques (compressed): " << std::fixed << std::setprecision(0)
              << totalRCliques << std::endl;
    std::cout << "  Total support sum (compressed): " << totalSupportByClasses << std::endl;
    std::cout << "  Total support sum (CPI):        " << totalSupportCPI << std::endl;
    std::cout << "  Peeled r-class tuples: " << peeledCount << " / " << rTuples.size() << std::endl;
    std::cout << "  Step ms [regions/classes/tuples/incidence/peel/output]: "
              << step1Ms << " / " << step2Ms << " / " << step3Ms
              << " / " << step4Ms << " / " << step5Ms << " / " << step6Ms << std::endl;
    std::cout << "  Total time: " << totalMs << " ms" << std::endl;
    std::cout << "==================================================" << std::endl;

    return result;
}
