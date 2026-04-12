static void maybeRunSparseSupportSharingSurvey(
    const DynamicGraph<TreeGraphNode> &tree, daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_SUPPORT_SHARING")) return;

    auto tStart = std::chrono::high_resolution_clock::now();
    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "============ Quotient Sparse Support Sharing =============" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip support-sharing survey" << std::endl;
        std::cout << "  Reason: global sparse key map would be too large" << std::endl;
        std::cout << "==========================================================" << std::endl;
        return;
    }

    robin_hood::unordered_flat_map<SparseCliqueKey, uint32_t, SparseCliqueKeyHash> leafIncidence;
    leafIncidence.reserve(static_cast<size_t>(std::min<uint64_t>(rawCliqueCount, 10000000ULL)));
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> supportByKey;
    supportByKey.reserve(static_cast<size_t>(std::min<uint64_t>(rawCliqueCount, 10000000ULL)));
    uint64_t totalLeafKeyIncidence = 0;
    daf::timeCount("sparse key incidence build (Quotient sharing)", [&]() {
        for (const auto &leaf : tree.adj_list) {
            accumulateLeafSparseIncidence(leaf, r, s, leafIncidence, &totalLeafKeyIncidence);
            accumulateLeafSparseContribution(leaf, r, s, supportByKey, 1.0);
        }
    });
    if (leafIncidence.empty() || supportByKey.empty()) return;

    uint64_t uniqueKeys = 0;
    uint64_t sharedKeys = 0;
    uint64_t shared2Keys = 0;
    uint64_t shared3pKeys = 0;
    uint64_t uniqueIncidence = 0;
    uint64_t sharedIncidence = 0;
    uint32_t maxIncidence = 0;
    std::vector<uint32_t> incidences;
    incidences.reserve(leafIncidence.size());
    for (const auto &[key, cnt] : leafIncidence) {
        (void)key;
        incidences.push_back(cnt);
        maxIncidence = std::max(maxIncidence, cnt);
        if (cnt == 1) {
            uniqueKeys++;
            uniqueIncidence += cnt;
        } else {
            sharedKeys++;
            sharedIncidence += cnt;
            if (cnt == 2) shared2Keys++;
            else shared3pKeys++;
        }
    }
    std::sort(incidences.begin(), incidences.end());
    auto pickInc = [&](double q) -> uint32_t {
        if (incidences.empty()) return 0;
        size_t idx = std::min(
            incidences.size() - 1,
            static_cast<size_t>(q * static_cast<double>(incidences.size() - 1)));
        return incidences[idx];
    };

    double minCore = std::numeric_limits<double>::max();
    for (const auto &[key, value] : supportByKey) {
        (void)key;
        minCore = std::min(minCore, value);
    }
    uint64_t frontierKeys = 0;
    uint64_t frontierUniqueKeys = 0;
    uint64_t frontierSharedKeys = 0;
    uint64_t frontierIncidence = 0;
    uint64_t frontierUniqueIncidence = 0;
    uint64_t frontierSharedIncidence = 0;
    uint32_t frontierMaxIncidence = 0;
    std::vector<uint32_t> frontierIncidences;
    frontierIncidences.reserve(supportByKey.size() / 64 + 1);
    for (const auto &[key, value] : supportByKey) {
        if (value > minCore) continue;
        auto it = leafIncidence.find(key);
        if (it == leafIncidence.end()) continue;
        const uint32_t cnt = it->second;
        frontierKeys++;
        frontierIncidence += cnt;
        frontierMaxIncidence = std::max(frontierMaxIncidence, cnt);
        frontierIncidences.push_back(cnt);
        if (cnt == 1) {
            frontierUniqueKeys++;
            frontierUniqueIncidence += cnt;
        } else {
            frontierSharedKeys++;
            frontierSharedIncidence += cnt;
        }
    }
    std::sort(frontierIncidences.begin(), frontierIncidences.end());
    auto pickFrontierInc = [&](double q) -> uint32_t {
        if (frontierIncidences.empty()) return 0;
        size_t idx = std::min(
            frontierIncidences.size() - 1,
            static_cast<size_t>(q * static_cast<double>(frontierIncidences.size() - 1)));
        return frontierIncidences[idx];
    };

    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "============ Quotient Sparse Support Sharing =============" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Sparse keys:               " << leafIncidence.size() << std::endl;
    std::cout << "  Total leaf-key incidence:  " << totalLeafKeyIncidence << std::endl;
    std::cout << std::fixed << std::setprecision(4)
              << "  Mean leaf incidence:       "
              << (static_cast<long double>(totalLeafKeyIncidence) /
                  static_cast<long double>(leafIncidence.size())) << std::endl
              << "  Unique-key ratio:          "
              << (100.0L * static_cast<long double>(uniqueKeys) /
                  static_cast<long double>(leafIncidence.size())) << "%" << std::endl
              << "  Shared-key ratio:          "
              << (100.0L * static_cast<long double>(sharedKeys) /
                  static_cast<long double>(leafIncidence.size())) << "%" << std::endl
              << "  Unique-incidence ratio:    "
              << (100.0L * static_cast<long double>(uniqueIncidence) /
                  static_cast<long double>(totalLeafKeyIncidence)) << "%" << std::endl
              << "  Shared-incidence ratio:    "
              << (100.0L * static_cast<long double>(sharedIncidence) /
                  static_cast<long double>(totalLeafKeyIncidence)) << "%" << std::endl
              << std::setprecision(2);
    std::cout << "  Keys with incidence=2:     " << shared2Keys << std::endl;
    std::cout << "  Keys with incidence>=3:    " << shared3pKeys << std::endl;
    std::cout << "  Incidence median/p90/max:  "
              << pickInc(0.50) << " / " << pickInc(0.90) << " / " << maxIncidence << std::endl;
    std::cout << "  Frontier min core:         " << minCore << std::endl;
    std::cout << "  Frontier keys:             " << frontierKeys << std::endl;
    if (frontierKeys > 0) {
        std::cout << std::fixed << std::setprecision(4)
                  << "  Frontier mean incidence:   "
                  << (static_cast<long double>(frontierIncidence) /
                      static_cast<long double>(frontierKeys)) << std::endl
                  << "  Frontier unique-key ratio: "
                  << (100.0L * static_cast<long double>(frontierUniqueKeys) /
                      static_cast<long double>(frontierKeys)) << "%" << std::endl
                  << "  Frontier shared-key ratio: "
                  << (100.0L * static_cast<long double>(frontierSharedKeys) /
                      static_cast<long double>(frontierKeys)) << "%" << std::endl
                  << "  Frontier unique-inc ratio: "
                  << (100.0L * static_cast<long double>(frontierUniqueIncidence) /
                      static_cast<long double>(frontierIncidence)) << "%" << std::endl
                  << "  Frontier shared-inc ratio: "
                  << (100.0L * static_cast<long double>(frontierSharedIncidence) /
                      static_cast<long double>(frontierIncidence)) << "%" << std::endl
                  << std::setprecision(2);
        std::cout << "  Frontier incidence med/p90/max: "
                  << pickFrontierInc(0.50) << " / " << pickFrontierInc(0.90)
                  << " / " << frontierMaxIncidence << std::endl;
    }
    std::cout << "  Survey time:               " << elapsedMs << " ms" << std::endl;
    std::cout << "==========================================================" << std::endl;
}

static void maybeRunFirstFrontierOwnershipAnalysis(
    const DynamicGraph<TreeGraphNode> &tree,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_FRONTIER_OWNERSHIP")) return;
    if (r != 3) {
        std::cout << "Frontier ownership: only r=3 supported, got r=" << r << std::endl;
        return;
    }

    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "======= Quotient First-Frontier Ownership =======" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip ownership analysis" << std::endl;
        std::cout << "  Reason: exact sparse support map would be too large" << std::endl;
        std::cout << "==================================================" << std::endl;
        return;
    }

    auto tStart = std::chrono::high_resolution_clock::now();
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> supportByKey;
    supportByKey.reserve(static_cast<size_t>(std::min<uint64_t>(rawCliqueCount, 10000000ULL)));
    daf::timeCount("frontier ownership support build", [&]() {
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

    std::vector<SparseCliqueKey> frontierKeys;
    frontierKeys.reserve(supportByKey.size() / 16 + 1);
    for (const auto &[key, value] : supportByKey) {
        if (std::abs(value - minCore) <= 1e-9) frontierKeys.push_back(key);
    }

    uint64_t onlyP2Keep1 = 0;
    uint64_t onlyP3 = 0;
    uint64_t mixedP2P3 = 0;
    uint64_t hasKeepGt1 = 0;
    uint64_t hasOtherClass = 0;
    uint64_t exactOneP2Witness = 0;
    uint64_t exactOneWitness = 0;
    uint64_t totalWitnessLeaves = 0;
    uint64_t maxWitnessLeaves = 0;
    std::unordered_map<uint64_t, uint64_t> frontierByClass;
    std::unordered_map<uint64_t, uint64_t> multiplicityByClass;
    std::vector<long double> classCoverage;
    std::map<std::pair<int, int>, uint64_t> frontierByKeepAndPivotRole;
    std::unordered_map<uint64_t, std::pair<int, int>> classShapeByKey;

    daf::StaticVector<daf::Size> verts;
    verts.resize(r);
    for (const auto &key : frontierKeys) {
        for (daf::CliqueSize j = 0; j < r; ++j) verts[j] = key.verts[j];

        bool hasP2Keep1Witness = false;
        bool hasP3Witness = false;
        bool hasKeepGt1Witness = false;
        bool hasOtherWitness = false;
        uint64_t witnessLeaves = 0;
        uint64_t p2Keep1Witnesses = 0;
        daf::Size uniqueLeafId = std::numeric_limits<daf::Size>::max();
        int uniqueSubNumPivot = -1;
        int uniqueKeepC = -1;
        int uniquePivotC = -1;

        daf::intersect_dense_sets_multi(verts, treeGraphV.adj_list,
            [&](const TreeGraphNode &uClique) {
                const auto &leaf = tree.adj_list[uClique.v];
                const auto role = analyzeLeafKeyRole(leaf, key, r, s);
                if (!role.valid) return;
                witnessLeaves++;
                uniqueLeafId = uClique.v;
                uniqueSubNumPivot = role.subNumPivot;
                uniqueKeepC = role.keepC;
                uniquePivotC = role.pivotC;
                if (role.keepC > 1) hasKeepGt1Witness = true;
                if (role.keepC == 1 && role.subNumPivot == 2) {
                    hasP2Keep1Witness = true;
                    p2Keep1Witnesses++;
                } else if (role.subNumPivot == 3) {
                    hasP3Witness = true;
                } else {
                    hasOtherWitness = true;
                }
            });

        totalWitnessLeaves += witnessLeaves;
        maxWitnessLeaves = std::max(maxWitnessLeaves, witnessLeaves);
        if (witnessLeaves == 1) exactOneWitness++;
        if (p2Keep1Witnesses == 1 && witnessLeaves == 1) exactOneP2Witness++;
        if (witnessLeaves == 1 &&
            uniqueLeafId != std::numeric_limits<daf::Size>::max() &&
            uniqueSubNumPivot >= 0) {
            const uint64_t classKey =
                (static_cast<uint64_t>(uniqueLeafId) << 8) |
                static_cast<uint64_t>(uniqueSubNumPivot);
            frontierByClass[classKey]++;
            frontierByKeepAndPivotRole[{uniqueKeepC, uniqueSubNumPivot}]++;
            classShapeByKey.emplace(classKey, std::make_pair(uniqueKeepC, uniqueSubNumPivot));
            if (multiplicityByClass.find(classKey) == multiplicityByClass.end()) {
                const int needKeep = static_cast<int>(r) - uniqueSubNumPivot;
                multiplicityByClass.emplace(
                    classKey,
                    combSmall(uniquePivotC, uniqueSubNumPivot) *
                    combSmall(uniqueKeepC, needKeep));
            }
        }

        if (hasKeepGt1Witness) {
            hasKeepGt1++;
        } else if (hasOtherWitness) {
            hasOtherClass++;
        } else if (hasP2Keep1Witness && hasP3Witness) {
            mixedP2P3++;
        } else if (hasP2Keep1Witness) {
            onlyP2Keep1++;
        } else if (hasP3Witness) {
            onlyP3++;
        } else {
            hasOtherClass++;
        }
    }
    verts.free();

    uint64_t totalFrontierClassMultiplicity = 0;
    uint64_t fullCoveredClasses = 0;
    uint64_t frontierKeysInFullCoveredClasses = 0;
    uint64_t frontierKeysInPartialClasses = 0;
    uint64_t classMultiplicityInFullCoveredClasses = 0;
    uint64_t classMultiplicityInPartialClasses = 0;
    std::map<std::pair<int, int>, uint64_t> fullCoveredFrontierByShape;
    std::map<std::pair<int, int>, uint64_t> partialFrontierByShape;
    for (const auto &[classKey, frontierCount] : frontierByClass) {
        const auto multIt = multiplicityByClass.find(classKey);
        if (multIt == multiplicityByClass.end() || multIt->second == 0) continue;
        const uint64_t mult = multIt->second;
        const auto shapeIt = classShapeByKey.find(classKey);
        totalFrontierClassMultiplicity += mult;
        const long double cover =
            static_cast<long double>(frontierCount) /
            static_cast<long double>(mult);
        classCoverage.push_back(cover);
        if (frontierCount == mult) {
            fullCoveredClasses++;
            frontierKeysInFullCoveredClasses += frontierCount;
            classMultiplicityInFullCoveredClasses += mult;
            if (shapeIt != classShapeByKey.end()) {
                fullCoveredFrontierByShape[shapeIt->second] += frontierCount;
            }
        } else {
            frontierKeysInPartialClasses += frontierCount;
            classMultiplicityInPartialClasses += mult;
            if (shapeIt != classShapeByKey.end()) {
                partialFrontierByShape[shapeIt->second] += frontierCount;
            }
        }
    }

    const auto fullCoveredCountOf = [&](int keepC, int q) -> uint64_t {
        auto it = fullCoveredFrontierByShape.find({keepC, q});
        return it == fullCoveredFrontierByShape.end() ? 0ULL : it->second;
    };
    const uint64_t policyCoreSafe =
        fullCoveredCountOf(2, 1) + fullCoveredCountOf(3, 0);
    const uint64_t policyCorePlusP2 =
        policyCoreSafe + fullCoveredCountOf(1, 2);
    const uint64_t policyAllFullCovered = frontierKeysInFullCoveredClasses;

    const auto ratioToFrontier = [&](uint64_t count) -> double {
        if (frontierKeys.empty()) return 0.0;
        return static_cast<double>(count) /
               static_cast<double>(frontierKeys.size());
    };
    std::sort(classCoverage.begin(), classCoverage.end());
    auto pickCoverage = [&](double q) -> long double {
        if (classCoverage.empty()) return 0.0L;
        size_t idx = static_cast<size_t>(
            q * static_cast<double>(classCoverage.size() - 1));
        return classCoverage[idx];
    };

    const auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "======= Quotient First-Frontier Ownership =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
    std::cout << "  First frontier core:        " << minCore << std::endl;
    std::cout << "  First frontier keys:        " << frontierKeys.size() << std::endl;
    std::cout << "  only p2 keepC=1:           " << onlyP2Keep1 << std::endl;
    std::cout << "  only p3:                   " << onlyP3 << std::endl;
    std::cout << "  mixed p2+p3:               " << mixedP2P3 << std::endl;
    std::cout << "  has keepC>1 witness:       " << hasKeepGt1 << std::endl;
    std::cout << "  has other witness:         " << hasOtherClass << std::endl;
    std::cout << "  exact one witness:         " << exactOneWitness << std::endl;
    std::cout << "  exact one p2 witness:      " << exactOneP2Witness << std::endl;
    std::cout << "  frontier witness classes:  " << frontierByClass.size() << std::endl;
    std::cout << "  class mult total:          " << totalFrontierClassMultiplicity << std::endl;
    std::cout << "  full-covered classes:      " << fullCoveredClasses << std::endl;
    std::cout << "  frontier in full classes:  " << frontierKeysInFullCoveredClasses << std::endl;
    std::cout << "  frontier in part classes:  " << frontierKeysInPartialClasses << std::endl;
    std::cout << "  mult in full classes:      " << classMultiplicityInFullCoveredClasses << std::endl;
    std::cout << "  mult in part classes:      " << classMultiplicityInPartialClasses << std::endl;
    std::cout << "  mean class coverage:       "
              << std::fixed << std::setprecision(6)
              << ((totalFrontierClassMultiplicity == 0)
                  ? 0.0
                  : static_cast<double>(frontierKeys.size()) /
                        static_cast<double>(totalFrontierClassMultiplicity))
              << std::endl;
    std::cout << "  class coverage med/p90:    "
              << pickCoverage(0.50) << " / " << pickCoverage(0.90) << std::endl;
    std::cout << "  exact frontier by (keepC, q):" << std::endl;
    for (const auto &[shape, count] : frontierByKeepAndPivotRole) {
        std::cout << "    keepC=" << shape.first
                  << " q=" << shape.second
                  << " count=" << count << std::endl;
    }
    std::cout << "  full-covered frontier by (keepC, q):" << std::endl;
    for (const auto &[shape, count] : fullCoveredFrontierByShape) {
        std::cout << "    keepC=" << shape.first
                  << " q=" << shape.second
                  << " count=" << count << std::endl;
    }
    std::cout << "  partial frontier by (keepC, q):" << std::endl;
    for (const auto &[shape, count] : partialFrontierByShape) {
        std::cout << "    keepC=" << shape.first
                  << " q=" << shape.second
                  << " count=" << count << std::endl;
    }
    std::cout << "  policy core-safe:          " << policyCoreSafe
              << " (" << std::fixed << std::setprecision(6)
              << ratioToFrontier(policyCoreSafe) << ")" << std::endl;
    std::cout << "  policy core+p2:            " << policyCorePlusP2
              << " (" << ratioToFrontier(policyCorePlusP2) << ")" << std::endl;
    std::cout << "  policy all full-covered:   " << policyAllFullCovered
              << " (" << ratioToFrontier(policyAllFullCovered) << ")" << std::endl;
    std::cout << "  mean witness leaves:       "
              << std::fixed << std::setprecision(4)
              << (frontierKeys.empty()
                  ? 0.0
                  : static_cast<double>(totalWitnessLeaves) /
                        static_cast<double>(frontierKeys.size()))
              << std::endl
              << std::setprecision(2);
    std::cout << "  max witness leaves:        " << maxWitnessLeaves << std::endl;
    std::cout << "  Analysis time:             " << elapsedMs << " ms" << std::endl;
    std::cout << "==================================================" << std::endl;
}

static void maybeRunPartialClassAnalysis(
    const DynamicGraph<TreeGraphNode> &tree,
    DynamicGraphSet<TreeGraphNode> &treeGraphV,
    daf::CliqueSize r, daf::CliqueSize s) {
    if (!std::getenv("PIVOTER_QUOTIENT_PARTIAL_CLASS_ANALYSIS")) return;
    if (r != 3) {
        std::cout << "Partial class analysis: only r=3 supported, got r=" << r << std::endl;
        return;
    }

    const auto rawCliqueCount = tree.cliqueCount(r);
    uint64_t rawCliqueGuard = 100000000ULL;
    if (const char *env = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_MAX_RAW")) {
        rawCliqueGuard = static_cast<uint64_t>(std::strtoull(env, nullptr, 10));
    }
    const bool forceRun = std::getenv("PIVOTER_QUOTIENT_ONE_ROUND_FORCE") != nullptr;
    if (!forceRun && static_cast<uint64_t>(rawCliqueCount) > rawCliqueGuard) {
        std::cout << "======= Quotient Partial-Class Analysis =======" << std::endl;
        std::cout << "  r=" << r << " s=" << s << std::endl;
        std::cout << "  Raw leaf-enum clique count: " << rawCliqueCount << std::endl;
        std::cout << "  Raw clique guard:           " << rawCliqueGuard << std::endl;
        std::cout << "  Status: skip partial-class analysis" << std::endl;
        std::cout << "  Reason: exact sparse support map would be too large" << std::endl;
        std::cout << "================================================" << std::endl;
        return;
    }

    auto tStart = std::chrono::high_resolution_clock::now();
    std::unordered_map<SparseCliqueKey, double, SparseCliqueKeyHash> supportByKey;
    supportByKey.reserve(static_cast<size_t>(std::min<uint64_t>(rawCliqueCount, 10000000ULL)));
    daf::timeCount("partial class support build", [&]() {
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

    std::unordered_map<uint64_t, uint64_t> frontierByClass;
    std::unordered_map<uint64_t, uint64_t> multiplicityByClass;
    std::unordered_map<uint64_t, std::pair<int, int>> classShapeByKey;
    std::unordered_map<uint64_t, daf::Size> leafIdByClass;

    std::vector<SparseCliqueKey> frontierKeys;
    frontierKeys.reserve(supportByKey.size() / 16 + 1);
    for (const auto &[key, value] : supportByKey) {
        if (std::abs(value - minCore) > 1e-9) continue;
        frontierKeys.push_back(key);

        daf::StaticVector<daf::Size> verts;
        verts.resize(r);
        for (daf::CliqueSize j = 0; j < r; ++j) verts[j] = key.verts[j];
        daf::intersect_dense_sets_multi(verts, treeGraphV.adj_list,
            [&](const TreeGraphNode &uClique) {
                const auto &leaf = tree.adj_list[uClique.v];
                const auto role = analyzeLeafKeyRole(leaf, key, r, s);
                if (!role.valid) return;
                const uint64_t classKey =
                    (static_cast<uint64_t>(uClique.v) << 8) |
                    static_cast<uint64_t>(role.subNumPivot);
                frontierByClass[classKey]++;
                classShapeByKey.emplace(classKey, std::make_pair(role.keepC, role.subNumPivot));
                leafIdByClass.emplace(classKey, uClique.v);
                if (multiplicityByClass.find(classKey) == multiplicityByClass.end()) {
                    const int needKeep = static_cast<int>(r) - role.subNumPivot;
                    multiplicityByClass.emplace(
                        classKey,
                        combSmall(role.pivotC, role.subNumPivot) *
                        combSmall(role.keepC, needKeep));
                }
            });
        verts.free();
    }

    struct PartialShapeSummary {
        uint64_t classCount = 0;
        uint64_t frontierMembers = 0;
        uint64_t nonFrontierMembers = 0;
        uint64_t plus1Members = 0;
        uint64_t plus2Members = 0;
        uint64_t plus3pMembers = 0;
    };

    std::map<std::pair<int, int>, PartialShapeSummary> partialByShape;
    uint64_t partialClassCount = 0;

    for (const auto &[classKey, frontierCount] : frontierByClass) {
        const auto multIt = multiplicityByClass.find(classKey);
        const auto shapeIt = classShapeByKey.find(classKey);
        const auto leafIt = leafIdByClass.find(classKey);
        if (multIt == multiplicityByClass.end() || shapeIt == classShapeByKey.end() ||
            leafIt == leafIdByClass.end()) {
            continue;
        }
        if (frontierCount >= multIt->second) continue;
        partialClassCount++;

        const auto shape = shapeIt->second;
        auto &stat = partialByShape[shape];
        stat.classCount++;
        stat.frontierMembers += frontierCount;

        const auto &leaf = tree.adj_list[leafIt->second];
        daf::enumerateCombinations(leaf, r,
            [&](const daf::StaticVector<TreeGraphNode> &rClique) {
                int subNumPivot = 0;
                SparseCliqueKey key;
                key.len = static_cast<uint8_t>(r);
                for (daf::CliqueSize j = 0; j < r; ++j) {
                    if (rClique[j].isPivot) subNumPivot++;
                    key.verts[j] = rClique[j].v;
                }
                if (subNumPivot != shape.second) return true;
                std::sort(key.verts.begin(), key.verts.begin() + r);
                auto supIt = supportByKey.find(key);
                if (supIt == supportByKey.end()) return true;
                const double support = supIt->second;
                if (std::abs(support - minCore) <= 1e-9) {
                    return true;
                }
                stat.nonFrontierMembers++;
                const double gap = support - minCore;
                if (std::abs(gap - 1.0) <= 1e-9) stat.plus1Members++;
                else if (std::abs(gap - 2.0) <= 1e-9) stat.plus2Members++;
                else if (gap >= 3.0 - 1e-9) stat.plus3pMembers++;
                return true;
            });
    }

    const auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tStart).count();

    std::cout << "======= Quotient Partial-Class Analysis =======" << std::endl;
    std::cout << "  r=" << r << " s=" << s << std::endl;
    std::cout << "  First frontier core:        " << minCore << std::endl;
    std::cout << "  First frontier keys:        " << frontierKeys.size() << std::endl;
    std::cout << "  Partial classes:            " << partialClassCount << std::endl;
    std::cout << "  Partial shape breakdown:" << std::endl;
    for (const auto &[shape, stat] : partialByShape) {
        std::cout << "    keepC=" << shape.first
                  << " q=" << shape.second
                  << " classes=" << stat.classCount
                  << " frontier=" << stat.frontierMembers
                  << " non_frontier=" << stat.nonFrontierMembers
                  << " plus1=" << stat.plus1Members
                  << " plus2=" << stat.plus2Members
                  << " plus3p=" << stat.plus3pMembers << std::endl;
    }
    std::cout << "  Analysis time:             " << elapsedMs << " ms" << std::endl;
    std::cout << "================================================" << std::endl;
}
