#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

namespace {

using Clock = chrono::steady_clock;

struct BigUInt {
    static constexpr uint32_t kBase = 1000000000U;
    vector<uint32_t> digits;

    BigUInt(uint64_t value = 0) {
        if (value == 0) {
            return;
        }
        while (value > 0) {
            digits.push_back(static_cast<uint32_t>(value % kBase));
            value /= kBase;
        }
    }

    bool isZero() const {
        return digits.empty();
    }

    void trim() {
        while (!digits.empty() && digits.back() == 0) {
            digits.pop_back();
        }
    }

    friend bool operator==(const BigUInt& lhs, const BigUInt& rhs) {
        return lhs.digits == rhs.digits;
    }

    friend bool operator<(const BigUInt& lhs, const BigUInt& rhs) {
        if (lhs.digits.size() != rhs.digits.size()) {
            return lhs.digits.size() < rhs.digits.size();
        }
        for (int i = static_cast<int>(lhs.digits.size()) - 1; i >= 0; --i) {
            if (lhs.digits[i] != rhs.digits[i]) {
                return lhs.digits[i] < rhs.digits[i];
            }
        }
        return false;
    }

    BigUInt& operator+=(const BigUInt& other) {
        const size_t n = max(digits.size(), other.digits.size());
        digits.resize(n, 0);
        uint64_t carry = 0;
        for (size_t i = 0; i < n; ++i) {
            uint64_t sum = carry + digits[i];
            if (i < other.digits.size()) {
                sum += other.digits[i];
            }
            digits[i] = static_cast<uint32_t>(sum % kBase);
            carry = sum / kBase;
        }
        if (carry > 0) {
            digits.push_back(static_cast<uint32_t>(carry));
        }
        return *this;
    }

    BigUInt& operator-=(const BigUInt& other) {
        uint64_t borrow = 0;
        for (size_t i = 0; i < digits.size(); ++i) {
            uint64_t lhsValue = digits[i];
            uint64_t rhsValue = borrow + (i < other.digits.size() ? other.digits[i] : 0U);
            if (lhsValue < rhsValue) {
                digits[i] = static_cast<uint32_t>(lhsValue + kBase - rhsValue);
                borrow = 1;
            } else {
                digits[i] = static_cast<uint32_t>(lhsValue - rhsValue);
                borrow = 0;
            }
        }
        trim();
        return *this;
    }

    BigUInt mulUint64(uint64_t factor) const {
        if (isZero() || factor == 0) {
            return 0;
        }
        BigUInt result;
        result.digits.resize(digits.size(), 0);
        uint64_t carry = 0;
        for (size_t i = 0; i < digits.size(); ++i) {
            unsigned __int128 cur =
                static_cast<unsigned __int128>(digits[i]) * factor + carry;
            result.digits[i] = static_cast<uint32_t>(cur % kBase);
            carry = static_cast<uint64_t>(cur / kBase);
        }
        while (carry > 0) {
            result.digits.push_back(static_cast<uint32_t>(carry % kBase));
            carry /= kBase;
        }
        result.trim();
        return result;
    }
};

ostream& operator<<(ostream& out, const BigUInt& value) {
    if (value.digits.empty()) {
        out << '0';
        return out;
    }
    out << value.digits.back();
    for (int i = static_cast<int>(value.digits.size()) - 2; i >= 0; --i) {
        out << setw(9) << setfill('0') << value.digits[i];
    }
    out << setfill(' ');
    return out;
}

struct Instance {
    int m = 0;
    int T = 0;
    vector<int> upper;
    vector<int> base;
    vector<vector<int>> boxes;
};

struct SolveStats {
    long long recursiveCalls = 0;
    long long memoHits = 0;
    long long feasibleMemoHits = 0;
    long long ieSubsets = 0;
    long long cellStates = 0;
    long long coveredCells = 0;
};

struct MethodResult {
    string name;
    bool supported = true;
    string note;
    BigUInt value = 0;
    double elapsedMs = 0.0;
    SolveStats stats;
};

string bigIntToString(const BigUInt& value) {
    ostringstream oss;
    oss << value;
    return oss.str();
}

string encodeVector(const vector<int>& values) {
    string key;
    key.reserve(values.size() * 4 + 1);
    for (int value : values) {
        key.append(to_string(value));
        key.push_back(',');
    }
    return key;
}

string encodeState(const vector<int>& lower,
                   const vector<int>& upper,
                   const vector<vector<int>>& boxes) {
    string key = encodeVector(lower);
    key.push_back('|');
    key += encodeVector(upper);
    key.push_back('|');
    for (const auto& box : boxes) {
        key += encodeVector(box);
        key.push_back(';');
    }
    return key;
}

bool dominates(const vector<int>& lhs, const vector<int>& rhs) {
    for (int i = 0; i < (int)lhs.size(); ++i) {
        if (lhs[i] > rhs[i]) {
            return false;
        }
    }
    return true;
}

int sumVector(const vector<int>& values) {
    return accumulate(values.begin(), values.end(), 0);
}

bool isRegionFeasible(const vector<int>& lower,
                      const vector<int>& upper,
                      int total) {
    int minSum = 0;
    int maxSum = 0;
    for (int i = 0; i < (int)lower.size(); ++i) {
        if (lower[i] > upper[i]) {
            return false;
        }
        minSum += lower[i];
        maxSum += upper[i];
    }
    return minSum <= total && total <= maxSum;
}

struct NormalizedBoxes {
    bool fullCover = false;
    vector<vector<int>> boxes;
};

NormalizedBoxes normalizeBoxes(const vector<int>& lower,
                               const vector<int>& upper,
                               int total,
                               const vector<vector<int>>& boxes,
                               bool pruneDominance) {
    vector<vector<int>> effective;
    effective.reserve(boxes.size());

    for (const auto& box : boxes) {
        vector<int> current(box.size());
        bool impossible = false;
        int lowerSum = 0;
        for (int i = 0; i < (int)box.size(); ++i) {
            current[i] = max(lower[i], box[i]);
            if (current[i] > upper[i]) {
                impossible = true;
                break;
            }
            lowerSum += current[i];
        }
        if (impossible || lowerSum > total) {
            continue;
        }
        if (current == lower) {
            return {true, {}};
        }
        effective.push_back(std::move(current));
    }

    sort(effective.begin(), effective.end());
    effective.erase(unique(effective.begin(), effective.end()), effective.end());

    if (!pruneDominance) {
        return {false, effective};
    }

    sort(effective.begin(), effective.end(), [](const vector<int>& lhs,
                                                const vector<int>& rhs) {
        int lhsSum = sumVector(lhs);
        int rhsSum = sumVector(rhs);
        if (lhsSum != rhsSum) {
            return lhsSum < rhsSum;
        }
        return lhs < rhs;
    });

    vector<vector<int>> minimal;
    for (const auto& box : effective) {
        bool contained = false;
        for (const auto& kept : minimal) {
            if (dominates(kept, box)) {
                contained = true;
                break;
            }
        }
        if (!contained) {
            minimal.push_back(box);
        }
    }
    return {false, minimal};
}

class UnionCounter {
public:
    explicit UnionCounter(Instance instance)
        : instance_(std::move(instance)) {}

    BigUInt countTotalFeasible() {
        return countFeasible(instance_.base, instance_.upper);
    }

    BigUInt countNaive() {
        vector<int> current(instance_.m, 0);
        return countNaiveRec(0, instance_.T, current);
    }

    BigUInt countBranchBasic(SolveStats* stats = nullptr) {
        activeStats_ = stats;
        BigUInt value = countUnionRec(instance_.base, instance_.upper, instance_.boxes,
                                      false, false, stats);
        activeStats_ = nullptr;
        return value;
    }

    BigUInt countBranchFast(SolveStats* stats = nullptr) {
        activeStats_ = stats;
        NormalizedBoxes raw =
            normalizeBoxes(instance_.base, instance_.upper, instance_.T,
                           instance_.boxes, false);
        NormalizedBoxes pruned =
            normalizeBoxes(instance_.base, instance_.upper, instance_.T,
                           instance_.boxes, true);
        bool usePruning =
            pruned.fullCover || pruned.boxes.size() < raw.boxes.size();
        BigUInt value = countUnionRec(instance_.base, instance_.upper, instance_.boxes,
                                      usePruning, false, stats);
        activeStats_ = nullptr;
        return value;
    }

    BigUInt countInclusionExclusion(SolveStats* stats = nullptr,
                                    int maxBoxes = 24) {
        activeStats_ = stats;
        NormalizedBoxes normalized =
            normalizeBoxes(instance_.base, instance_.upper, instance_.T,
                           instance_.boxes, true);
        if (normalized.fullCover) {
            BigUInt value = countFeasible(instance_.base, instance_.upper);
            activeStats_ = nullptr;
            return value;
        }
        if (normalized.boxes.size() > static_cast<size_t>(maxBoxes)) {
            activeStats_ = nullptr;
            throw runtime_error("inclusion-exclusion disabled: too many minimal boxes");
        }
        BigUInt positive = 0;
        BigUInt negative = 0;
        vector<int> current = instance_.base;
        countIeRec(0, 0, current, normalized.boxes, &positive, &negative, stats);
        activeStats_ = nullptr;
        if (positive < negative) {
            throw runtime_error("inclusion-exclusion internal error: negative total");
        }
        positive -= negative;
        return positive;
    }

    BigUInt countThresholdLattice(SolveStats* stats = nullptr) {
        activeStats_ = stats;
        NormalizedBoxes normalized =
            normalizeBoxes(instance_.base, instance_.upper, instance_.T,
                           instance_.boxes, true);
        if (normalized.fullCover) {
            BigUInt value = countFeasible(instance_.base, instance_.upper);
            activeStats_ = nullptr;
            return value;
        }
        if (normalized.boxes.empty()) {
            activeStats_ = nullptr;
            return 0;
        }

        vector<int> activeDims;
        vector<vector<int>> levels;
        for (int c = 0; c < instance_.m; ++c) {
            vector<int> values = {instance_.base[c], instance_.upper[c] + 1};
            for (const auto& box : normalized.boxes) {
                if (box[c] > instance_.base[c]) {
                    values.push_back(box[c]);
                }
            }
            sort(values.begin(), values.end());
            values.erase(unique(values.begin(), values.end()), values.end());
            if (values.size() > 2) {
                activeDims.push_back(c);
                levels.push_back(std::move(values));
            }
        }

        if (activeDims.empty()) {
            activeStats_ = nullptr;
            return 0;
        }

        vector<int> lower = instance_.base;
        vector<int> upper = instance_.upper;
        BigUInt total = 0;
        enumerateCellsRec(0, activeDims, levels, normalized.boxes,
                          &lower, &upper, &total, stats);
        activeStats_ = nullptr;
        return total;
    }

private:
    int choosePivotBox(const vector<int>& lower,
                       const vector<vector<int>>& boxes) const {
        int bestIndex = 0;
        int bestActiveDims = numeric_limits<int>::max();
        int bestLowerSum = -1;
        for (int i = 0; i < (int)boxes.size(); ++i) {
            int activeDims = 0;
            int lowerSum = 0;
            for (int c = 0; c < instance_.m; ++c) {
                if (boxes[i][c] > lower[c]) {
                    ++activeDims;
                }
                lowerSum += boxes[i][c];
            }
            if (activeDims < bestActiveDims ||
                (activeDims == bestActiveDims && lowerSum > bestLowerSum)) {
                bestIndex = i;
                bestActiveDims = activeDims;
                bestLowerSum = lowerSum;
            }
        }
        return bestIndex;
    }

    BigUInt countNaiveRec(int index, int remaining, vector<int>& current) {
        if (index == instance_.m) {
            if (remaining != 0) {
                return 0;
            }
            for (const auto& box : instance_.boxes) {
                bool inside = true;
                for (int c = 0; c < instance_.m; ++c) {
                    if (current[c] < box[c]) {
                        inside = false;
                        break;
                    }
                }
                if (inside) {
                    return 1;
                }
            }
            return 0;
        }

        BigUInt total = 0;
        int minValue = instance_.base[index];
        int maxValue = min(instance_.upper[index], remaining);
        int suffixUpper = 0;
        int suffixLower = 0;
        for (int c = index + 1; c < instance_.m; ++c) {
            suffixUpper += instance_.upper[c];
            suffixLower += instance_.base[c];
        }

        for (int value = minValue; value <= maxValue; ++value) {
            int nextRemaining = remaining - value;
            if (nextRemaining < suffixLower || nextRemaining > suffixUpper) {
                continue;
            }
            current[index] = value;
            total += countNaiveRec(index + 1, nextRemaining, current);
        }
        return total;
    }

    void countIeRec(int index,
                    int pickedCount,
                    vector<int>& currentLower,
                    const vector<vector<int>>& boxes,
                    BigUInt* positive,
                    BigUInt* negative,
                    SolveStats* stats) {
        if (index == static_cast<int>(boxes.size())) {
            if (pickedCount == 0) {
                return;
            }
            if (stats) {
                ++stats->ieSubsets;
            }
            BigUInt value = countFeasible(currentLower, instance_.upper);
            if ((pickedCount & 1) == 1) {
                *positive += value;
            } else {
                *negative += value;
            }
            return;
        }

        countIeRec(index + 1, pickedCount, currentLower, boxes,
                   positive, negative, stats);

        vector<int> saved = currentLower;
        bool impossible = false;
        int lowerSum = 0;
        for (int c = 0; c < instance_.m; ++c) {
            currentLower[c] = max(currentLower[c], boxes[index][c]);
            if (currentLower[c] > instance_.upper[c]) {
                impossible = true;
                break;
            }
            lowerSum += currentLower[c];
        }
        if (!impossible && lowerSum <= instance_.T) {
            countIeRec(index + 1, pickedCount + 1, currentLower, boxes,
                       positive, negative, stats);
        }
        currentLower.swap(saved);
    }

    void enumerateCellsRec(int depth,
                           const vector<int>& activeDims,
                           const vector<vector<int>>& levels,
                           const vector<vector<int>>& boxes,
                           vector<int>* lower,
                           vector<int>* upper,
                           BigUInt* total,
                           SolveStats* stats) {
        if (depth == static_cast<int>(activeDims.size())) {
            if (stats) {
                ++stats->cellStates;
            }
            bool covered = false;
            for (const auto& box : boxes) {
                bool inside = true;
                for (int c = 0; c < instance_.m; ++c) {
                    if ((*lower)[c] < box[c]) {
                        inside = false;
                        break;
                    }
                }
                if (inside) {
                    covered = true;
                    break;
                }
            }
            if (!covered) {
                return;
            }
            if (stats) {
                ++stats->coveredCells;
            }
            *total += countFeasible(*lower, *upper);
            return;
        }

        int c = activeDims[depth];
        int oldLower = (*lower)[c];
        int oldUpper = (*upper)[c];
        for (int i = 0; i + 1 < static_cast<int>(levels[depth].size()); ++i) {
            (*lower)[c] = levels[depth][i];
            (*upper)[c] = levels[depth][i + 1] - 1;
            enumerateCellsRec(depth + 1, activeDims, levels, boxes,
                              lower, upper, total, stats);
        }
        (*lower)[c] = oldLower;
        (*upper)[c] = oldUpper;
    }

    BigUInt countFeasible(const vector<int>& lower,
                          const vector<int>& upper) {
        if (!isRegionFeasible(lower, upper, instance_.T)) {
            return 0;
        }

        string key = encodeVector(lower);
        key.push_back('|');
        key += encodeVector(upper);
        auto it = feasibleMemo_.find(key);
        if (it != feasibleMemo_.end()) {
            if (activeStats_ != nullptr) {
                ++activeStats_->feasibleMemoHits;
            }
            return it->second;
        }

        int minSum = sumVector(lower);
        int target = instance_.T - minSum;
        vector<int> caps(instance_.m, 0);
        int capSum = 0;
        for (int i = 0; i < instance_.m; ++i) {
            caps[i] = upper[i] - lower[i];
            capSum += caps[i];
        }
        if (target < 0 || target > capSum) {
            feasibleMemo_[key] = 0;
            return 0;
        }

        vector<BigUInt> dp(target + 1);
        dp[0] = 1;
        for (int cap : caps) {
            vector<BigUInt> next(target + 1);
            BigUInt window = 0;
            for (int s = 0; s <= target; ++s) {
                window += dp[s];
                if (s - cap - 1 >= 0) {
                    window -= dp[s - cap - 1];
                }
                next[s] = window;
            }
            dp.swap(next);
        }

        feasibleMemo_[key] = dp[target];
        return dp[target];
    }

    BigUInt countUnionRec(const vector<int>& lower,
                         const vector<int>& upper,
                         const vector<vector<int>>& boxes,
                         bool pruneDominance,
                         bool useMemo,
                         SolveStats* stats) {
        if (stats) {
            ++stats->recursiveCalls;
        }
        if (!isRegionFeasible(lower, upper, instance_.T)) {
            return 0;
        }

        NormalizedBoxes normalized =
            normalizeBoxes(lower, upper, instance_.T, boxes, pruneDominance);
        if (normalized.fullCover) {
            return countFeasible(lower, upper);
        }
        if (normalized.boxes.empty()) {
            return 0;
        }

        if (useMemo) {
            string stateKey = encodeState(lower, upper, normalized.boxes);
            auto it = unionMemo_.find(stateKey);
            if (it != unionMemo_.end()) {
                if (stats) {
                    ++stats->memoHits;
                }
                return it->second;
            }

            BigUInt value =
                countUnionCore(lower, upper, normalized.boxes, pruneDominance,
                               useMemo, stats);
            unionMemo_[stateKey] = value;
            return value;
        }

        return countUnionCore(lower, upper, normalized.boxes, pruneDominance,
                              useMemo, stats);
    }

    BigUInt countUnionCore(const vector<int>& lower,
                          const vector<int>& upper,
                          const vector<vector<int>>& boxes,
                          bool pruneDominance,
                          bool useMemo,
                          SolveStats* stats) {
        int pivotIndex = choosePivotBox(lower, boxes);
        vector<int> pivot = boxes[pivotIndex];

        BigUInt total = countFeasible(pivot, upper);

        vector<vector<int>> remaining = boxes;
        remaining.erase(remaining.begin() + pivotIndex);

        for (int splitDim = 0; splitDim < instance_.m; ++splitDim) {
            if (pivot[splitDim] <= lower[splitDim]) {
                continue;
            }

            vector<int> nextLower = lower;
            vector<int> nextUpper = upper;
            for (int earlier = 0; earlier < splitDim; ++earlier) {
                nextLower[earlier] = max(nextLower[earlier], pivot[earlier]);
            }
            nextUpper[splitDim] = min(nextUpper[splitDim], pivot[splitDim] - 1);

            total += countUnionRec(nextLower, nextUpper, remaining,
                                   pruneDominance, useMemo, stats);
        }
        return total;
    }

    Instance instance_;
    SolveStats* activeStats_ = nullptr;
    unordered_map<string, BigUInt> feasibleMemo_;
    unordered_map<string, BigUInt> unionMemo_;
};

Instance readInstance(const string& path) {
    ifstream in(path);
    if (!in) {
        throw runtime_error("Cannot open input file: " + path);
    }

    Instance instance;
    in >> instance.m >> instance.T;
    instance.upper.resize(instance.m);
    instance.base.resize(instance.m);
    for (int i = 0; i < instance.m; ++i) {
        in >> instance.upper[i];
    }
    for (int i = 0; i < instance.m; ++i) {
        in >> instance.base[i];
    }
    int k = 0;
    in >> k;
    instance.boxes.assign(k, vector<int>(instance.m, 0));
    for (int i = 0; i < k; ++i) {
        for (int c = 0; c < instance.m; ++c) {
            in >> instance.boxes[i][c];
        }
    }
    if (!in) {
        throw runtime_error("Malformed input file: " + path);
    }
    return instance;
}

Instance randomInstance(mt19937& rng,
                        int m,
                        int total,
                        int k,
                        int upperMin,
                        int upperMax,
                        int witnessWidth) {
    Instance instance;
    instance.m = m;
    instance.T = total;
    instance.upper.resize(m);

    uniform_int_distribution<int> upperDist(upperMin, upperMax);
    int upperSum = 0;
    for (int c = 0; c < m; ++c) {
        instance.upper[c] = upperDist(rng);
        upperSum += instance.upper[c];
    }
    if (upperSum < total) {
        int deficit = total - upperSum;
        for (int c = 0; c < m && deficit > 0; ++c) {
            int add = min(deficit, upperMax);
            instance.upper[c] += add;
            deficit -= add;
        }
    }

    instance.base.assign(m, 0);
    int baseBudget = total / 3;
    uniform_int_distribution<int> classDist(0, m - 1);
    while (baseBudget > 0) {
        int c = classDist(rng);
        if (instance.base[c] < instance.upper[c]) {
            ++instance.base[c];
            --baseBudget;
        }
    }

    instance.boxes.assign(k, vector<int>(m, 0));
    uniform_int_distribution<int> widthDist(1, max(1, witnessWidth));
    for (int i = 0; i < k; ++i) {
        int width = min(m, widthDist(rng));
        vector<int> dims(m);
        iota(dims.begin(), dims.end(), 0);
        shuffle(dims.begin(), dims.end(), rng);
        dims.resize(width);

        int remaining = total;
        for (int c = 0; c < m; ++c) {
            instance.boxes[i][c] = 0;
        }
        for (int idx = 0; idx < width; ++idx) {
            int c = dims[idx];
            int cap = min(instance.upper[c], remaining);
            if (cap == 0) {
                continue;
            }
            uniform_int_distribution<int> reqDist(1, cap);
            instance.boxes[i][c] = reqDist(rng);
            remaining -= instance.boxes[i][c] / 2;
            remaining = max(0, remaining);
        }
    }

    return instance;
}

Instance makeCyclicOverlapInstance(int m,
                                   int total,
                                   int upperCap,
                                   int rotation,
                                   bool anchoredBase) {
    Instance instance;
    instance.m = m;
    instance.T = total;
    instance.upper.assign(m, upperCap);
    instance.base.assign(m, 0);
    if (anchoredBase) {
        for (int c = 0; c < m; c += 2) {
            instance.base[c] = 1;
        }
    }

    for (int i = 0; i < m; ++i) {
        vector<int> box(m, 0);
        int a = (i + rotation) % m;
        int b = (i + 1 + rotation) % m;
        box[a] = anchoredBase ? 3 : 4;
        box[b] = anchoredBase ? 3 : 4;
        if (anchoredBase) {
            for (int c = 0; c < m; c += 2) {
                box[c] = max(box[c], 1);
            }
        }
        instance.boxes.push_back(box);
    }

    for (int i = 0; i < m / 2; ++i) {
        vector<int> box(m, 0);
        int a = (2 * i + rotation) % m;
        int b = (2 * i + 2 + rotation) % m;
        int c = (2 * i + 4 + rotation) % m;
        box[a] = 2;
        box[b] = 2;
        box[c] = 2;
        if (anchoredBase) {
            for (int d = 0; d < m; d += 2) {
                box[d] = max(box[d], 1);
            }
        }
        instance.boxes.push_back(box);
    }
    return instance;
}

Instance makeDominatedRingInstance(int m,
                                   int total,
                                   int upperCap,
                                   int rotation) {
    Instance instance;
    instance.m = m;
    instance.T = total;
    instance.upper.assign(m, upperCap);
    instance.base.assign(m, 0);

    for (int i = 0; i < m; ++i) {
        int a = (i + rotation) % m;
        int b = (i + 1 + rotation) % m;

        for (const auto& [ra, rb] : vector<pair<int, int>>{
                 {3, 3}, {4, 3}, {3, 4}, {4, 4}}) {
            vector<int> box(m, 0);
            box[a] = ra;
            box[b] = rb;
            instance.boxes.push_back(box);
        }
    }
    return instance;
}

Instance makeFewBoxesManyLevelsInstance(int rotation) {
    Instance instance;
    instance.m = 8;
    instance.T = 20;
    instance.upper.assign(instance.m, 8);
    instance.base.assign(instance.m, 0);

    vector<pair<int, int>> pairs = {
        {0, 1}, {1, 2}, {2, 3}, {3, 4},
        {4, 5}, {5, 6}, {6, 7}, {7, 0}
    };
    vector<pair<int, int>> levels = {
        {2, 5}, {3, 6}, {4, 2}, {5, 3},
        {6, 4}, {2, 6}, {3, 4}, {5, 2}
    };

    for (int i = 0; i < static_cast<int>(pairs.size()); ++i) {
        vector<int> box(instance.m, 0);
        int a = (pairs[i].first + rotation) % instance.m;
        int b = (pairs[i].second + rotation) % instance.m;
        box[a] = levels[i].first;
        box[b] = levels[i].second;
        instance.boxes.push_back(box);
    }
    return instance;
}

Instance makeAllPairsLowLevelsInstance(int m,
                                       int total,
                                       int upperCap,
                                       int threshold) {
    Instance instance;
    instance.m = m;
    instance.T = total;
    instance.upper.assign(m, upperCap);
    instance.base.assign(m, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            vector<int> box(m, 0);
            box[i] = threshold;
            box[j] = threshold;
            instance.boxes.push_back(box);
        }
    }
    return instance;
}

void buildSubsetBoxesRec(int nextDim,
                         int need,
                         int threshold,
                         vector<int>* current,
                         vector<vector<int>>* boxes) {
    if (need == 0) {
        boxes->push_back(*current);
        return;
    }
    if (nextDim + need > static_cast<int>(current->size())) {
        return;
    }
    for (int c = nextDim; c < static_cast<int>(current->size()); ++c) {
        (*current)[c] = threshold;
        buildSubsetBoxesRec(c + 1, need - 1, threshold, current, boxes);
        (*current)[c] = 0;
    }
}

Instance makeAllSubsetsLowLevelsInstance(int m,
                                         int total,
                                         int upperCap,
                                         int threshold,
                                         int subsetSize) {
    Instance instance;
    instance.m = m;
    instance.T = total;
    instance.upper.assign(m, upperCap);
    instance.base.assign(m, 0);
    vector<int> current(m, 0);
    buildSubsetBoxesRec(0, subsetSize, threshold, &current, &instance.boxes);
    return instance;
}

bool verifyInstance(const Instance& instance,
                    bool printFailure) {
    UnionCounter counter(instance);
    BigUInt naive = counter.countNaive();

    SolveStats basicStats;
    UnionCounter counterBasic(instance);
    BigUInt basic = counterBasic.countBranchBasic(&basicStats);

    SolveStats fastStats;
    UnionCounter counterFast(instance);
    BigUInt fast = counterFast.countBranchFast(&fastStats);

    SolveStats ieStats;
    UnionCounter counterIe(instance);
    BigUInt ie = counterIe.countInclusionExclusion(&ieStats);

    SolveStats latticeStats;
    UnionCounter counterLattice(instance);
    BigUInt lattice = counterLattice.countThresholdLattice(&latticeStats);

    bool ok = (naive == basic) && (naive == fast) &&
              (naive == ie) && (naive == lattice);
    if (!ok && printFailure) {
        cerr << "Verification failed.\n";
        cerr << "naive=" << bigIntToString(naive)
             << " basic=" << bigIntToString(basic)
             << " fast=" << bigIntToString(fast)
             << " ie=" << bigIntToString(ie)
             << " lattice=" << bigIntToString(lattice) << "\n";
        cerr << "m=" << instance.m << " T=" << instance.T << "\nupper:";
        for (int x : instance.upper) cerr << ' ' << x;
        cerr << "\nbase:";
        for (int x : instance.base) cerr << ' ' << x;
        cerr << "\nboxes:\n";
        for (const auto& box : instance.boxes) {
            for (int x : box) cerr << x << ' ';
            cerr << "\n";
        }
    }
    return ok;
}

template <typename Fn>
pair<BigUInt, double> timedRun(Fn&& fn) {
    auto start = Clock::now();
    BigUInt value = fn();
    auto end = Clock::now();
    chrono::duration<double, milli> elapsed = end - start;
    return {value, elapsed.count()};
}

template <typename Fn>
MethodResult runMethod(const string& name, Fn&& fn) {
    MethodResult result;
    result.name = name;
    try {
        auto [value, elapsedMs] = timedRun(fn);
        result.value = std::move(value);
        result.elapsedMs = elapsedMs;
    } catch (const exception& ex) {
        string message = ex.what();
        if (message.find("disabled") != string::npos) {
            result.supported = false;
            result.note = std::move(message);
        } else {
            throw;
        }
    }
    return result;
}

void printSolveReport(const Instance& instance,
                      const BigUInt& totalFeasible,
                      const vector<MethodResult>& methods) {
    cout << "Problem summary\n";
    cout << "  classes=" << instance.m
         << " total=" << instance.T
         << " boxes=" << instance.boxes.size() << "\n";
    cout << "  upper:";
    for (int value : instance.upper) cout << ' ' << value;
    cout << "\n  base:";
    for (int value : instance.base) cout << ' ' << value;
    cout << "\n";

    cout << "\nCounts\n";
    cout << "  feasible_containing_tau'=" << bigIntToString(totalFeasible) << "\n";
    for (const auto& method : methods) {
        if (method.supported) {
            cout << "  dead_" << method.name << "="
                 << bigIntToString(method.value) << "\n";
        } else {
            cout << "  dead_" << method.name << "=UNSUPPORTED"
                 << " (" << method.note << ")\n";
        }
    }

    cout << "\nTiming (ms)\n";
    for (const auto& method : methods) {
        if (!method.supported) {
            cout << "  " << method.name << "=SKIPPED\n";
            continue;
        }
        cout << "  " << method.name << "=" << method.elapsedMs;
        if (method.stats.recursiveCalls > 0) {
            cout << " recursive_calls=" << method.stats.recursiveCalls;
        }
        if (method.stats.memoHits > 0 || method.stats.feasibleMemoHits > 0) {
            cout << " memo_hits=" << method.stats.memoHits
                 << " feasible_memo_hits=" << method.stats.feasibleMemoHits;
        }
        if (method.stats.ieSubsets > 0) {
            cout << " ie_subsets=" << method.stats.ieSubsets;
        }
        if (method.stats.cellStates > 0) {
            cout << " cell_states=" << method.stats.cellStates
                 << " covered_cells=" << method.stats.coveredCells;
        }
        cout << "\n";
    }
}

void runSelfTest(int cases, uint32_t seed) {
    mt19937 rng(seed);
    for (int i = 0; i < cases; ++i) {
        int m = uniform_int_distribution<int>(3, 8)(rng);
        int total = uniform_int_distribution<int>(4, 16)(rng);
        int k = uniform_int_distribution<int>(1, 12)(rng);
        int width = uniform_int_distribution<int>(1, min(4, m))(rng);
        Instance instance = randomInstance(rng, m, total, k, 1, 6, width);
        if (!verifyInstance(instance, true)) {
            throw runtime_error("self-test failed");
        }
    }
    cout << "Self-test passed: " << cases << " random instances, seed=" << seed << "\n";
}

void runBenchmark(int casesPerProfile, uint32_t seed) {
    struct Profile {
        string name;
        function<Instance(int)> generator;
    };

    vector<Profile> profiles = {
        {"few_boxes_many_levels", [](int idx) {
             return makeFewBoxesManyLevelsInstance(idx % 8);
         }},
        {"all_pairs_low_levels", [](int) {
             return makeAllPairsLowLevelsInstance(10, 16, 4, 2);
         }},
        {"all_quadruples_low_levels", [](int) {
             return makeAllSubsetsLowLevelsInstance(10, 16, 4, 2, 4);
         }},
        {"cyclic_overlap", [](int idx) {
             return makeCyclicOverlapInstance(8, 16, 6, idx % 8, false);
         }},
        {"anchored_overlap", [](int idx) {
             return makeCyclicOverlapInstance(8, 16, 6, idx % 8, true);
         }},
        {"dominated_ring", [](int idx) {
             return makeDominatedRingInstance(10, 20, 6, idx % 10);
         }},
    };

    cout << "profile,cases,naive_ms_avg,ie_ms_avg,lattice_ms_avg,"
            "branch_basic_ms_avg,branch_fast_ms_avg,"
            "naive_over_ie,naive_over_lattice,naive_over_branch_fast\n";

    for (const auto& profile : profiles) {
        double naiveMs = 0.0;
        double ieMs = 0.0;
        int ieRuns = 0;
        double latticeMs = 0.0;
        double basicMs = 0.0;
        double fastMs = 0.0;

        for (int i = 0; i < casesPerProfile; ++i) {
            Instance instance = profile.generator(i + seed);

            MethodResult naive = runMethod("naive", [&] {
                UnionCounter counter(instance);
                return counter.countNaive();
            });
            SolveStats ieStats;
            MethodResult ie = runMethod("ie", [&] {
                UnionCounter counter(instance);
                return counter.countInclusionExclusion(&ieStats);
            });
            ie.stats = ieStats;
            SolveStats latticeStats;
            MethodResult lattice = runMethod("lattice", [&] {
                UnionCounter counter(instance);
                return counter.countThresholdLattice(&latticeStats);
            });
            lattice.stats = latticeStats;
            SolveStats basicStats;
            MethodResult basic = runMethod("branch_basic", [&] {
                UnionCounter counter(instance);
                return counter.countBranchBasic(&basicStats);
            });
            basic.stats = basicStats;
            SolveStats fastStats;
            MethodResult fast = runMethod("branch_fast", [&] {
                UnionCounter counter(instance);
                return counter.countBranchFast(&fastStats);
            });
            fast.stats = fastStats;

            vector<MethodResult> methods = {naive, ie, lattice, basic, fast};
            for (const auto& method : methods) {
                if (!method.supported) {
                    continue;
                }
                if (!(method.value == naive.value)) {
                    cerr << "Benchmark correctness mismatch in profile "
                         << profile.name << " method=" << method.name << "\n";
                    throw runtime_error("benchmark failed");
                }
            }

            naiveMs += naive.elapsedMs;
            if (ie.supported) {
                ieMs += ie.elapsedMs;
                ++ieRuns;
            }
            if (lattice.supported) {
                latticeMs += lattice.elapsedMs;
            }
            basicMs += basic.elapsedMs;
            fastMs += fast.elapsedMs;
        }

        naiveMs /= casesPerProfile;
        latticeMs /= casesPerProfile;
        basicMs /= casesPerProfile;
        fastMs /= casesPerProfile;
        double ieAvg = (ieRuns > 0) ? ieMs / ieRuns : -1.0;

        cout << profile.name << ','
             << casesPerProfile << ','
             << naiveMs << ',';
        if (ieRuns > 0) {
            cout << ieAvg;
        } else {
            cout << "NA";
        }
        cout << ','
             << latticeMs << ','
             << basicMs << ','
             << fastMs << ',';
        if (ieRuns > 0 && ieAvg > 0.0) {
            cout << (naiveMs / ieAvg);
        } else {
            cout << "NA";
        }
        cout << ','
             << (latticeMs > 0.0 ? naiveMs / latticeMs : 0.0) << ','
             << (fastMs > 0.0 ? naiveMs / fastMs : 0.0) << "\n";
    }
}

void printUsage(const char* program) {
    cerr << "Usage:\n";
    cerr << "  " << program << " solve <input.txt>\n";
    cerr << "  " << program
         << " method <input.txt> <total|naive|ie|lattice|branch_basic|branch_fast>\n";
    cerr << "  " << program << " selftest [cases] [seed]\n";
    cerr << "  " << program << " bench [cases_per_profile] [seed]\n";
    cerr << "\nInput format for solve:\n";
    cerr << "  line1: m T\n";
    cerr << "  line2: upper[0..m-1]\n";
    cerr << "  line3: base[0..m-1]\n";
    cerr << "  line4: k\n";
    cerr << "  next k lines: box_i[0..m-1]\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            printUsage(argv[0]);
            return 1;
        }

        string mode = argv[1];
        if (mode == "solve") {
            if (argc != 3) {
                printUsage(argv[0]);
                return 1;
            }
            Instance instance = readInstance(argv[2]);

            UnionCounter totalCounter(instance);
            BigUInt totalFeasible = totalCounter.countTotalFeasible();

            vector<MethodResult> methods;
            methods.push_back(runMethod("naive", [&] {
                UnionCounter counter(instance);
                return counter.countNaive();
            }));

            SolveStats ieStats;
            MethodResult ie = runMethod("ie", [&] {
                UnionCounter counter(instance);
                return counter.countInclusionExclusion(&ieStats);
            });
            ie.stats = ieStats;
            methods.push_back(ie);

            SolveStats latticeStats;
            MethodResult lattice = runMethod("lattice", [&] {
                UnionCounter counter(instance);
                return counter.countThresholdLattice(&latticeStats);
            });
            lattice.stats = latticeStats;
            methods.push_back(lattice);

            SolveStats basicStats;
            MethodResult basic = runMethod("branch_basic", [&] {
                UnionCounter counter(instance);
                return counter.countBranchBasic(&basicStats);
            });
            basic.stats = basicStats;
            methods.push_back(basic);

            SolveStats fastStats;
            MethodResult fast = runMethod("branch_fast", [&] {
                UnionCounter counter(instance);
                return counter.countBranchFast(&fastStats);
            });
            fast.stats = fastStats;
            methods.push_back(fast);

            const BigUInt& reference = methods.front().value;
            for (const auto& method : methods) {
                if (!method.supported) {
                    continue;
                }
                if (!(method.value == reference)) {
                    cerr << "Mismatch detected: reference="
                         << bigIntToString(reference)
                         << " method=" << method.name
                         << " value=" << bigIntToString(method.value) << "\n";
                    return 2;
                }
            }

            printSolveReport(instance, totalFeasible, methods);
            return 0;
        }

        if (mode == "selftest") {
            int cases = (argc >= 3) ? stoi(argv[2]) : 200;
            uint32_t seed = (argc >= 4) ? static_cast<uint32_t>(stoul(argv[3])) : 42u;
            runSelfTest(cases, seed);
            return 0;
        }

        if (mode == "method") {
            if (argc != 4) {
                printUsage(argv[0]);
                return 1;
            }
            Instance instance = readInstance(argv[2]);
            string method = argv[3];

            SolveStats stats;
            auto [value, elapsedMs] = timedRun([&] {
                UnionCounter counter(instance);
                if (method == "total") {
                    return counter.countTotalFeasible();
                }
                if (method == "naive") {
                    return counter.countNaive();
                }
                if (method == "ie") {
                    return counter.countInclusionExclusion(&stats);
                }
                if (method == "lattice") {
                    return counter.countThresholdLattice(&stats);
                }
                if (method == "branch_basic") {
                    return counter.countBranchBasic(&stats);
                }
                if (method == "branch_fast") {
                    return counter.countBranchFast(&stats);
                }
                throw runtime_error("Unknown method: " + method);
            });

            cout << "method=" << method << "\n";
            cout << "count=" << bigIntToString(value) << "\n";
            cout << "time_ms=" << elapsedMs << "\n";
            if (method == "branch_basic" || method == "branch_fast") {
                cout << "recursive_calls=" << stats.recursiveCalls << "\n";
                cout << "memo_hits=" << stats.memoHits << "\n";
                cout << "feasible_memo_hits=" << stats.feasibleMemoHits << "\n";
            }
            if (method == "ie") {
                cout << "ie_subsets=" << stats.ieSubsets << "\n";
            }
            if (method == "lattice") {
                cout << "cell_states=" << stats.cellStates << "\n";
                cout << "covered_cells=" << stats.coveredCells << "\n";
            }
            return 0;
        }

        if (mode == "bench") {
            int cases = (argc >= 3) ? stoi(argv[2]) : 20;
            uint32_t seed = (argc >= 4) ? static_cast<uint32_t>(stoul(argv[3])) : 123u;
            runBenchmark(cases, seed);
            return 0;
        }

        printUsage(argv[0]);
        return 1;
    } catch (const exception& ex) {
        cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}
