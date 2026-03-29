#pragma once
#include <vector>
#include <set>
#include <algorithm>
#include <cstdint>

// Hybrid Bucket + std::set priority queue for peeling.
// Bucket array for support ≤ THRESHOLD (O(1) ops).
// Overflow set for support > THRESHOLD (O(log n) ops).
// Handles R≥3 large s where nCr values can reach 10^308.
struct HybridBucketPQ {
    static constexpr double THRESHOLD = 5000000.0;

    std::vector<std::vector<uint32_t>> buckets;
    std::set<std::pair<double, uint32_t>> overflowSet;
    std::vector<int> bucket_of;        // -1 = in overflow
    std::vector<uint32_t> pos_in_bucket;
    std::vector<double> overflowStoredVal;
    int curBucket = 0;

    void init(uint32_t n, const std::vector<double> &count) {
        double rawMax = 0;
        for (uint32_t i = 0; i < n; ++i)
            rawMax = std::max(rawMax, count[i]);
        int maxB = (int)std::min(rawMax, THRESHOLD);
        buckets.assign(maxB + 2, {});
        overflowSet.clear();
        bucket_of.assign(n, -1);
        pos_in_bucket.resize(n);
        overflowStoredVal.assign(n, -1);

        for (uint32_t i = 0; i < n; ++i) {
            if (count[i] <= THRESHOLD) {
                int b = (int)count[i];
                bucket_of[i] = b;
                pos_in_bucket[i] = (uint32_t)buckets[b].size();
                buckets[b].push_back(i);
            } else {
                overflowSet.insert({count[i], i});
                overflowStoredVal[i] = count[i];
            }
        }
        curBucket = 0;
    }

    // Move id after its count changed. val = new count value.
    void move(uint32_t id, double val) {
        val = std::max(0.0, val);
        int oldB = bucket_of[id];

        if (oldB == -1)
            overflowSet.erase({overflowStoredVal[id], id});

        if (val <= THRESHOLD) {
            int newB = (int)val;
            if (oldB >= 0 && newB == oldB) return;
            if (oldB >= 0) {
                auto &v = buckets[oldB];
                uint32_t p = pos_in_bucket[id];
                if (p < (uint32_t)v.size() - 1) {
                    uint32_t last = v.back();
                    v[p] = last;
                    pos_in_bucket[last] = p;
                }
                v.pop_back();
            }
            bucket_of[id] = newB;
            pos_in_bucket[id] = (uint32_t)buckets[newB].size();
            buckets[newB].push_back(id);
            if (newB < curBucket) curBucket = newB;
        } else {
            if (oldB >= 0) {
                auto &v = buckets[oldB];
                uint32_t p = pos_in_bucket[id];
                if (p < (uint32_t)v.size()) {
                    if (p < (uint32_t)v.size() - 1) {
                        uint32_t last = v.back();
                        v[p] = last;
                        pos_in_bucket[last] = p;
                    }
                    v.pop_back();
                }
            }
            overflowSet.insert({val, id});
            overflowStoredVal[id] = val;
            bucket_of[id] = -1;
        }
    }

    // Drain overflow items below threshold into buckets.
    template<typename InHeapFn>
    void drainOverflow(const std::vector<double> &count, InHeapFn &&inHeap) {
        while (!overflowSet.empty()) {
            uint32_t id = overflowSet.begin()->second;
            if (!inHeap(id)) { overflowSet.erase(overflowSet.begin()); continue; }
            if (count[id] <= THRESHOLD) {
                overflowSet.erase(overflowSet.begin());
                int b = (int)count[id];
                bucket_of[id] = b;
                pos_in_bucket[id] = (uint32_t)buckets[b].size();
                buckets[b].push_back(id);
            } else break;
        }
    }

    // Pop min-support items. Returns false if empty.
    template<typename InHeapFn>
    bool popMin(double &minCore,
                std::vector<uint32_t> &out,
                const std::vector<double> &count,
                InHeapFn &&inHeap,
                auto &&markPopped, // void(uint32_t id, double core)
                uint32_t &remaining) {
        drainOverflow(count, inHeap);

        while (curBucket < (int)buckets.size() && buckets[curBucket].empty())
            curBucket++;

        if (curBucket >= (int)buckets.size()) {
            while (!overflowSet.empty()) {
                uint32_t id = overflowSet.begin()->second;
                overflowSet.erase(overflowSet.begin());
                if (!inHeap(id)) continue;
                minCore = std::max(count[id], minCore);
                markPopped(id, minCore);
                out.push_back(id);
                remaining--;
                while (!overflowSet.empty()) {
                    uint32_t next = overflowSet.begin()->second;
                    if (!inHeap(next)) { overflowSet.erase(overflowSet.begin()); continue; }
                    if (count[next] <= minCore) {
                        overflowSet.erase(overflowSet.begin());
                        markPopped(next, minCore);
                        out.push_back(next);
                        remaining--;
                    } else break;
                }
                return !out.empty();
            }
            return false;
        }

        minCore = std::max((double)curBucket, minCore);
        while (curBucket < (int)buckets.size() && !buckets[curBucket].empty()
               && curBucket <= (int)minCore) {
            while (!buckets[curBucket].empty()) {
                auto id = buckets[curBucket].back();
                buckets[curBucket].pop_back();
                markPopped(id, minCore);
                out.push_back(id);
                remaining--;
            }
            if (curBucket + 1 < (int)buckets.size() && !buckets[curBucket + 1].empty()
                && (curBucket + 1) <= (int)minCore)
                curBucket++;
            else break;
        }
        return !out.empty();
    }
};
