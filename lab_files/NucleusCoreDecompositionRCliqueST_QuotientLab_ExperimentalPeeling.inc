// ================================================================
// Class-level peeling prototype for r=3.
// ================================================================
static void maybeRunClassLevelPeelingPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_CLASS_PEEL")) return;
    if (r != 3) { std::cout << "Class peel: only r=3\n"; return; }

    auto tStart = std::chrono::high_resolution_clock::now();

    (void)edgeGraph;
    (void)treeGraphV;

    struct ClassEntry {
        daf::Size leafId;
        int keepC;
        int pivotC;
        int subNumPivot;
        double contribution;
        uint64_t multiplicity;
    };
    std::vector<ClassEntry> allClasses;
    allClasses.reserve(tree.adj_list.size() * 2);
    uint64_t totalIncidences = 0;

    for (daf::Size lid = 0; lid < tree.adj_list.size(); ++lid) {
        const auto &leaf = tree.adj_list[lid];
        if (leaf.size() < r) continue;
        daf::Size kc = 0, pc = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pc++; else kc++;
        }
        for (int q = 0; q <= static_cast<int>(r); ++q) {
            int needKeep = static_cast<int>(r) - q;
            if (q > static_cast<int>(pc) || needKeep > static_cast<int>(kc)) continue;
            int row = static_cast<int>(pc) - q;
            int col = static_cast<int>(s) - static_cast<int>(kc) - q;
            if (row < 0 || col < 0 || col > row || row >= 1001 || col >= 401) continue;
            double contrib = nCr[row][col];
            if (contrib <= 0.0) continue;
            uint64_t mult = combSmall(pc, q) * combSmall(kc, needKeep);
            if (mult == 0) continue;
            allClasses.push_back({lid, static_cast<int>(kc), static_cast<int>(pc), q, contrib, mult});
            totalIncidences += mult;
        }
    }

    std::sort(allClasses.begin(), allClasses.end(),
        [](const ClassEntry &a, const ClassEntry &b) {
            return a.contribution < b.contribution;
        });

    std::cout << "======= Class-Level Peeling (Generalized) =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Total class entries: " << allClasses.size() << std::endl;
    std::cout << "  Total incidences (with dup): " << totalIncidences << std::endl;

    if (!allClasses.empty()) {
        const double minClassContrib = allClasses.front().contribution;
        struct Phase1ClassSummary {
            uint64_t classCount = 0;
            uint64_t multiplicity = 0;
        };
        std::map<std::pair<int, int>, Phase1ClassSummary> phase1ByKeepAndPivotRole;
        uint64_t phase1Classes = 0;
        uint64_t phase1Multiplicity = 0;
        for (const auto &ce : allClasses) {
            if (std::abs(ce.contribution - minClassContrib) > 1e-9) break;
            phase1Classes++;
            phase1Multiplicity += ce.multiplicity;
            auto &slot = phase1ByKeepAndPivotRole[{ce.keepC, ce.subNumPivot}];
            slot.classCount++;
            slot.multiplicity += ce.multiplicity;
        }
        std::cout << "  Phase-1 class contrib:     " << minClassContrib << std::endl;
        std::cout << "  Phase-1 class entries:     " << phase1Classes << std::endl;
        std::cout << "  Phase-1 incidences:        " << phase1Multiplicity << std::endl;
        std::cout << "  Phase-1 by (keepC, q):" << std::endl;
        for (const auto &[shape, stat] : phase1ByKeepAndPivotRole) {
            std::cout << "    keepC=" << shape.first
                      << " q=" << shape.second
                      << " classes=" << stat.classCount
                      << " mult=" << stat.multiplicity << std::endl;
        }
    }

    uint64_t totalPeeled = 0;
    double maxCore = 0;
    int rounds = 0;
    std::map<int, uint64_t> coreDist;

    for (const auto &ce : allClasses) {
        double core = std::max(maxCore, ce.contribution);
        maxCore = core;
        totalPeeled += ce.multiplicity;
        int coreInt = static_cast<int>(std::round(core));
        coreDist[coreInt] += ce.multiplicity;
        rounds++;
        if (rounds <= 10 || rounds % 10000 == 0) {
            std::cout << "  Peel round " << rounds
                      << ": leaf=" << ce.leafId
                      << " p=" << ce.subNumPivot
                      << " contrib=" << ce.contribution
                      << " mult=" << ce.multiplicity
                      << " core=" << core
                      << " total=" << totalPeeled << std::endl;
        }
    }

    std::cout << "  Core distribution (class-level, with duplicates):" << std::endl;
    int shown = 0;
    for (auto &[c, n] : coreDist) {
        if (shown++ < 10 || c == coreDist.rbegin()->first)
            std::cout << "    core=" << c << " count=" << n << std::endl;
    }
    if (coreDist.size() > 11) std::cout << "    ... (" << coreDist.size() << " levels total)" << std::endl;

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "  Total rounds:     " << rounds << std::endl;
    std::cout << "  Total peeled:     " << totalPeeled << " (p=2 class only)" << std::endl;
    std::cout << "  Max core:         " << maxCore << std::endl;
    std::cout << "  Time:             " << elapsedMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;
}

// ================================================================
// Vertex-level peeling prototype for r=3.
// ================================================================
static void maybeRunVertexLevelPeelingPrototype(
    const DynamicGraph<TreeGraphNode> &tree, const Graph &edgeGraph,
    DynamicGraphSet<TreeGraphNode> &treeGraphV, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_VERTEX_PEEL")) return;
    if (r != 3) {
        std::cout << "Vertex-level peeling: only r=3 supported, got r=" << r << std::endl;
        return;
    }

    auto tStart = std::chrono::high_resolution_clock::now();
    int maxRounds = 200;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_MULTI_ROUND_MAX")) {
        maxRounds = std::max(1, std::atoi(env));
    }

    struct LeafInfo {
        daf::Size keepVertex;
        daf::Size keepC, pivotC;
        bool alive;
    };
    std::vector<LeafInfo> leafInfo(tree.adj_list.size());
    for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.size() < r) { leafInfo[leafId] = {0, 0, 0, false}; continue; }
        daf::Size kc = 0, pc = 0, kv = 0;
        for (const auto &node : leaf) {
            if (node.isPivot) pc++;
            else { kc++; kv = node.v; }
        }
        leafInfo[leafId] = {kv, kc, pc, true};
    }

    struct VtxLeafEntry { daf::Size leafId; bool isPivot; };
    std::vector<std::vector<VtxLeafEntry>> vtxLeaves(edgeGraph.n);
    for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
        if (!leafInfo[leafId].alive) continue;
        for (const auto &node : tree.adj_list[leafId]) {
            if (node.v < edgeGraph.n) {
                vtxLeaves[node.v].push_back({leafId, node.isPivot});
            }
        }
    }

    robin_hood::unordered_flat_set<SparseCliqueKey, SparseCliqueKeyHash> peeled;
    auto makeKey3 = [](daf::Size a, daf::Size b, daf::Size c) -> SparseCliqueKey {
        SparseCliqueKey k;
        k.len = 3;
        daf::Size v[3] = {a, b, c};
        std::sort(v, v + 3);
        k.verts[0] = v[0];
        k.verts[1] = v[1];
        k.verts[2] = v[2];
        return k;
    };

    auto computeTriangleSupport = [&](daf::Size v0, daf::Size v1, daf::Size v2) -> double {
        double support = 0.0;
        for (const auto &entry0 : vtxLeaves[v0]) {
            if (!leafInfo[entry0.leafId].alive) continue;
            const auto &leaf = tree.adj_list[entry0.leafId];
            bool found1 = false, found2 = false;
            for (const auto &node : leaf) {
                if (node.v == v1) found1 = true;
                if (node.v == v2) found2 = true;
            }
            if (!found1 || !found2) continue;
            for (const auto &node : leaf) {
                daf::Size d = node.v;
                if (d == v0 || d == v1 || d == v2) continue;
                bool poisoned = false;
                if (!peeled.empty()) {
                    if (peeled.find(makeKey3(v0, v1, d)) != peeled.end()) poisoned = true;
                    if (!poisoned && peeled.find(makeKey3(v0, v2, d)) != peeled.end()) poisoned = true;
                    if (!poisoned && peeled.find(makeKey3(v1, v2, d)) != peeled.end()) poisoned = true;
                }
                if (!poisoned) support += 1.0;
            }
        }
        return support;
    };

    std::vector<daf::Size> vertexOrder;
    vertexOrder.reserve(edgeGraph.n);
    for (daf::Size v = 0; v < edgeGraph.n; ++v) {
        bool hasKeepLeaf = false;
        for (const auto &entry : vtxLeaves[v]) {
            if (!entry.isPivot && leafInfo[entry.leafId].alive) {
                hasKeepLeaf = true;
                break;
            }
        }
        if (hasKeepLeaf) vertexOrder.push_back(v);
    }

    auto getKeepLeafPivotC = [&](daf::Size v) -> daf::Size {
        for (const auto &entry : vtxLeaves[v]) {
            if (!entry.isPivot && leafInfo[entry.leafId].alive) {
                return leafInfo[entry.leafId].pivotC;
            }
        }
        return 999999;
    };
    std::sort(vertexOrder.begin(), vertexOrder.end(), [&](daf::Size a, daf::Size b) {
        return getKeepLeafPivotC(a) < getKeepLeafPivotC(b);
    });

    int completedRounds = 0;
    uint64_t totalPeeled = 0;
    double maxCoreValue = 0.0;
    uint64_t totalTrianglesChecked = 0;

    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> supportMap;
    std::map<double, std::vector<SparseCliqueKey>> buckets;

    auto ensureInMap = [&](daf::Size v0, daf::Size v1, daf::Size v2) {
        SparseCliqueKey key;
        key.len = 3;
        daf::Size vs[3] = {v0, v1, v2};
        std::sort(vs, vs + 3);
        key.verts[0] = vs[0];
        key.verts[1] = vs[1];
        key.verts[2] = vs[2];
        if (peeled.find(key) != peeled.end()) return;
        if (supportMap.find(key) != supportMap.end()) return;
        double sup = computeTriangleSupport(vs[0], vs[1], vs[2]);
        totalTrianglesChecked++;
        if (sup > 0.0) {
            supportMap[key] = sup;
            buckets[sup].push_back(key);
        }
    };

    daf::Size nextLevelIdx = 0;
    auto enumerateUpToSupport = [&](double maxSup) {
        while (nextLevelIdx < vertexOrder.size()) {
            daf::Size v = vertexOrder[nextLevelIdx];
            daf::Size pivotC = getKeepLeafPivotC(v);
            double keepContrib = (pivotC >= 2) ? nCr[pivotC - 2][static_cast<int>(s) - 3] : 0.0;
            if (keepContrib > maxSup) break;

            daf::Size keepLeafId = 0;
            bool found = false;
            for (const auto &entry : vtxLeaves[v]) {
                if (!entry.isPivot && leafInfo[entry.leafId].alive) {
                    keepLeafId = entry.leafId;
                    found = true;
                    break;
                }
            }
            if (!found) {
                nextLevelIdx++;
                continue;
            }

            const auto &leaf = tree.adj_list[keepLeafId];
            std::vector<daf::Size> pivots;
            for (const auto &node : leaf) {
                if (node.isPivot) pivots.push_back(node.v);
            }
            for (daf::Size i = 0; i < pivots.size(); ++i) {
                for (daf::Size j = i + 1; j < pivots.size(); ++j) {
                    ensureInMap(v, pivots[i], pivots[j]);
                }
            }
            for (daf::Size i = 0; i < pivots.size(); ++i) {
                for (daf::Size j = i + 1; j < pivots.size(); ++j) {
                    for (daf::Size k = j + 1; k < pivots.size(); ++k) {
                        ensureInMap(pivots[i], pivots[j], pivots[k]);
                    }
                }
            }
            nextLevelIdx++;
        }
    };

    (void)enumerateUpToSupport;
    (void)treeGraphV;

    for (daf::Size leafId = 0; leafId < tree.adj_list.size(); ++leafId) {
        if (!leafInfo[leafId].alive) continue;
        const auto &leaf = tree.adj_list[leafId];
        if (leaf.size() < r) continue;
        const auto &li = leafInfo[leafId];
        daf::enumerateCombinations(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                int snp = 0;
                SparseCliqueKey key;
                key.len = 3;
                for (daf::CliqueSize j = 0; j < r; ++j) {
                    if (rClique[j].isPivot) snp++;
                    key.verts[j] = rClique[j].v;
                }
                std::sort(key.verts.begin(), key.verts.begin() + r);
                if (!hasPositiveContribution(li.pivotC, li.keepC, snp, r, s)) return true;
                const int row = li.pivotC - snp;
                const int col = static_cast<int>(s) - li.keepC - snp;
                if (row < 0 || row >= 1001 || col < 0 || col >= 401 || col > row) return true;
                const double contrib = nCr[row][col];
                if (contrib <= 0.0) return true;
                supportMap[key] += contrib;
                totalTrianglesChecked++;
                return true;
            });
    }
    for (auto &[key, sup] : supportMap) {
        buckets[sup].push_back(key);
    }

    std::map<int, int> supportDist;
    for (auto &[k, v] : supportMap) supportDist[static_cast<int>(std::round(v))]++;
    std::cout << "  Initial enumeration: " << supportMap.size() << " unique keys" << std::endl;
    std::cout << "  Support distribution (first 20):" << std::endl;
    int shown = 0;
    for (auto &[sv, c] : supportDist) {
        if (shown++ < 20) {
            std::cout << "    support=" << sv << " count=" << c << std::endl;
        }
    }
    std::cout << "    max support=" << supportDist.rbegin()->first << std::endl;

    while (completedRounds < maxRounds) {
        double minBucket = -1.0;
        while (!buckets.empty()) {
            auto it = buckets.begin();
            auto &vec = it->second;
            vec.erase(std::remove_if(vec.begin(), vec.end(), [&](const SparseCliqueKey &k) {
                if (peeled.find(k) != peeled.end()) return true;
                auto sit = supportMap.find(k);
                if (sit == supportMap.end()) return true;
                return std::abs(sit->second - it->first) > 1e-9;
            }), vec.end());
            if (vec.empty()) {
                buckets.erase(it);
                continue;
            }
            minBucket = it->first;
            break;
        }

        if (minBucket < 0) break;

        double coreValue = std::max(maxCoreValue, static_cast<double>(minBucket));
        maxCoreValue = coreValue;

        std::vector<SparseCliqueKey> frontierKeys;
        while (!buckets.empty() && buckets.begin()->first <= coreValue + 1e-9) {
            auto bit = buckets.begin();
            for (auto &k : bit->second) {
                if (peeled.find(k) == peeled.end()) {
                    auto sit = supportMap.find(k);
                    if (sit != supportMap.end() && sit->second <= coreValue + 1e-9) {
                        frontierKeys.push_back(k);
                    }
                }
            }
            buckets.erase(bit);
        }
        uint64_t roundPeeled = 0;

        daf::Size fi = 0;
        while (fi < frontierKeys.size()) {
            const auto &key = frontierKeys[fi++];
            if (peeled.find(key) != peeled.end()) continue;
            auto sit = supportMap.find(key);
            if (sit == supportMap.end()) continue;
            if (sit->second > coreValue + 1e-9) continue;

            peeled.insert(key);
            supportMap.erase(sit);
            roundPeeled++;
            totalPeeled++;

            daf::Size kv[3] = {key.verts[0], key.verts[1], key.verts[2]};
            for (const auto &entry0 : vtxLeaves[kv[0]]) {
                if (!leafInfo[entry0.leafId].alive) continue;
                const auto &leaf = tree.adj_list[entry0.leafId];
                bool f1 = false, f2 = false;
                for (const auto &node : leaf) {
                    if (node.v == kv[1]) f1 = true;
                    if (node.v == kv[2]) f2 = true;
                }
                if (!f1 || !f2) continue;

                bool kIncludesKeep = false;
                for (int vi = 0; vi < 3; ++vi) {
                    for (const auto &node : leaf) {
                        if (!node.isPivot && node.v == kv[vi]) {
                            kIncludesKeep = true;
                            break;
                        }
                    }
                    if (kIncludesKeep) break;
                }

                auto updateNeighbor = [&](daf::Size d) {
                    for (int rep = 0; rep < 3; ++rep) {
                        daf::Size nvs[3];
                        for (int vi = 0, ni = 0; vi < 3; ++vi) {
                            if (vi == rep) nvs[ni++] = d;
                            else nvs[ni++] = kv[vi];
                        }
                        std::sort(nvs, nvs + 3);
                        SparseCliqueKey nkey;
                        nkey.len = 3;
                        nkey.verts[0] = nvs[0];
                        nkey.verts[1] = nvs[1];
                        nkey.verts[2] = nvs[2];
                        if (peeled.find(nkey) != peeled.end()) return;
                        auto nit = supportMap.find(nkey);
                        if (nit == supportMap.end()) return;
                        nit->second -= 1.0;
                        if (nit->second < 0.0) nit->second = 0.0;
                        if (nit->second <= coreValue + 1e-9) {
                            frontierKeys.push_back(nkey);
                        } else {
                            buckets[nit->second].push_back(nkey);
                        }
                    }
                };

                auto isAlreadyDestroyed = [&](daf::Size d) -> bool {
                    for (int rep = 0; rep < 3; ++rep) {
                        daf::Size nvs[3];
                        for (int vi = 0, ni = 0; vi < 3; ++vi) {
                            if (vi == rep) nvs[ni++] = d;
                            else nvs[ni++] = kv[vi];
                        }
                        std::sort(nvs, nvs + 3);
                        SparseCliqueKey nk;
                        nk.len = 3;
                        nk.verts[0] = nvs[0];
                        nk.verts[1] = nvs[1];
                        nk.verts[2] = nvs[2];
                        if (peeled.find(nk) != peeled.end()) return true;
                    }
                    return false;
                };

                if (kIncludesKeep) {
                    for (const auto &node : leaf) {
                        if (node.v == kv[0] || node.v == kv[1] || node.v == kv[2]) continue;
                        if (!isAlreadyDestroyed(node.v)) updateNeighbor(node.v);
                    }
                } else {
                    for (const auto &node : leaf) {
                        if (!node.isPivot) {
                            if (!isAlreadyDestroyed(node.v)) updateNeighbor(node.v);
                            break;
                        }
                    }
                }
            }
        }

        completedRounds++;
        if (completedRounds <= 20 || completedRounds % 100 == 0) {
            std::cout << "  VtxPeel round " << completedRounds
                      << ": coreValue=" << coreValue
                      << " peeled=" << roundPeeled
                      << " total=" << totalPeeled
                      << " supportMap=" << supportMap.size()
                      << " buckets=" << buckets.size() << std::endl;
        }
    }

    if (!supportMap.empty()) {
        double minRemaining = 1e18, maxRemaining = 0;
        int inBuckets = 0;
        for (auto &[k, v] : supportMap) {
            minRemaining = std::min(minRemaining, v);
            maxRemaining = std::max(maxRemaining, v);
            if (buckets.count(v) && !buckets[v].empty()) inBuckets++;
        }
        std::cout << "  Remaining in supportMap: " << supportMap.size()
                  << " minSupport=" << minRemaining << " maxSupport=" << maxRemaining
                  << " inBuckets=" << inBuckets
                  << " bucketCount=" << buckets.size() << std::endl;
    }

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "======= Vertex-Level Peeling Prototype =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Completed rounds:        " << completedRounds << std::endl;
    std::cout << "  Total peeled triangles:  " << totalPeeled << std::endl;
    std::cout << "  Total checked triangles: " << totalTrianglesChecked << std::endl;
    std::cout << "  Max core value:          " << maxCoreValue << std::endl;
    std::cout << "  Prototype time:          " << elapsedMs << " ms" << std::endl;
    std::cout << "==============================================" << std::endl;
}
