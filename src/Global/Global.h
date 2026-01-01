//
// Created by _ on 24-7-29.
//
// #pragma once
#ifndef SUBGRAPHMATCHING_GLOBAL_H
#define SUBGRAPHMATCHING_GLOBAL_H
#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>
#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <ranges>
#include <chrono>
#include <map>
#include <span>
#include <tbb/spin_mutex.h>
#include <csignal>
#include <memory>
// #include <google/dense_hash_set>

#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>

#include <iomanip>
#include <random>
#if defined(__APPLE__)
#include <mach/mach.h>
#endif
#include "dataStruct/robin_hood.h"
#include "dataStruct/Bitset.hpp"
#include "Global/ostreamOverload.hpp"


#define MAX_CSIZE 400

namespace daf {
    using Size = uint32_t;
    using CliqueSize = uint16_t;
    static constexpr Size UNDEFINED = std::numeric_limits<Size>::max();

    void coverGraphToAdjFileDisk(const std::string &intput);

    void coverGraphToAdjFileMemory(const std::string &intput);

    constexpr Size INVALID_SIZE = std::numeric_limits<Size>::max();

    static constexpr Size accuarcy = 100000;

    template<typename T>
    void printArray(const T *array, Size size) {
        for (Size i = 0; i < size; ++i) {
            std::cout << "{" << i << ": " << array[i] << "} ";
        }
        std::cout << std::endl;
    }


    inline void log_memory(const std::string& label) {
        long rss_kb = 0;

#if defined(__linux__)
        // ==========================================
        // Linux 实现: 读取 /proc/self/status
        // ==========================================
        std::ifstream status_file("/proc/self/status");
        std::string line;
        while (std::getline(status_file, line)) {
            if (line.substr(0, 6) == "VmRSS:") {
                // line 格式通常是 "VmRSS:    1234 kB"
                // 我们提取数字部分
                std::string val_str = line.substr(7);
                // 去掉尾部的 " kB" (简单处理，直接打印整行也可以，这里为了统一定义个变量)
                // 实际上直接打印 line 也是可以的
                std::cout << "[Memory-Linux] " << label << ": " << line.substr(6) << std::endl;
                return;
            }
        }

#elif defined(__APPLE__)
        // ==========================================
        // macOS (Darwin) 实现: 使用 task_info
        // ==========================================
        struct task_basic_info t_info;
        mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

        if (task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count) == KERN_SUCCESS) {
            // macOS 返回的是字节 (Bytes)，我们需要除以 1024 换算成 KB，以便与 Linux 统一
            rss_kb = t_info.resident_size / 1024;
            std::cout << "[Memory-Mac] " << label << ": " << rss_kb << " kB" << std::endl;
        } else {
            std::cerr << "[Memory-Mac] Failed to get memory info" << std::endl;
        }

#else
        // ==========================================
        // 其他系统 (如 Windows)
        // ==========================================
        std::cout << "[Memory] " << label << ": OS not supported for memory tracking" << std::endl;
#endif
    }

    // 用法:
    // load_graph(G);
    // log_memory("After Graph Load");
    // build_tree(T);
    // log_memory("After Tree Build");

    template<typename T>
    class StaticVector {
        static_assert(std::is_copy_constructible_v<T>,
                      "StaticVector requires T to be copy-constructible");

    public:
        using iterator = T *;
        using const_iterator = const T *;
        using value_type = T;
        using Size = daf::Size;

        //  / 
        Size c_size = 0;
        Size maxSize = 0;

    private:
        // ，/；
        std::shared_ptr<T[]> data_;

    public:
        explicit StaticVector(Size cap = MAX_CSIZE)
            : c_size(0), maxSize(cap),
              data_(cap
                        ? std::shared_ptr<T[]>(
                            new T[cap], std::default_delete<T[]>()) // C++17 
                        : std::shared_ptr<T[]>()) {
        }

        // /：（），O(1)
        StaticVector(const StaticVector &) = default;

        StaticVector &operator=(const StaticVector &) = default;

        StaticVector(StaticVector &&) noexcept = default;

        StaticVector &operator=(StaticVector &&) noexcept = default;

        // （shared_ptr ）
        ~StaticVector() = default;

        // 
        explicit operator T *() { return data(); }
        T *getData() const { return data(); }
        T *data() const { return data_.get(); }

        void push_back_with_check(const T &v) {
            if (c_size >= maxSize) reAllocate(grow_capacity());
            data()[c_size++] = v;
        }

        void push_back(const T &v) {
            if (c_size >= maxSize) reAllocate(grow_capacity());
            data()[c_size++] = v;
        }

        template<typename... Args>
        void emplace_back(Args &&... args) {
            if (c_size >= maxSize) reAllocate(grow_capacity());
            new(data() + c_size) T(std::forward<Args>(args)...);
            ++c_size;
        }

        void pop_back() { --c_size; }

        // O(1) ：，
        void remove(Size index) {
#ifndef NDEBUG
            if (index >= c_size) throw std::out_of_range("Index out of range.");
#endif
            data()[index] = std::move(data()[--c_size]);
        }

        // （：，）
        void extend(const T *begin, Size sz) {
            std::memcpy(data() + c_size, begin, sz * sizeof(T));
            c_size += sz;
        }

        template<typename BinaryPred = std::equal_to<T> >
        void unique(BinaryPred comp = BinaryPred()) noexcept {
            if (c_size < 2) return;
            Size write = 1;
            for (Size read = 1; read < c_size; ++read) {
                if (!comp(data()[read], data()[write - 1])) {
                    data()[write++] = data()[read];
                }
            }
            c_size = write;
        }

        [[nodiscard]] bool empty() const { return c_size == 0; }
        void clear() { c_size = 0; }

        [[nodiscard]] StaticVector deepCopy() const {
            StaticVector copy(maxSize);
            copy.c_size = c_size;
            if constexpr (std::is_trivially_copyable_v<T>) {
                std::memcpy(copy.data(), data(), c_size * sizeof(T));
            } else {
                for (Size i = 0; i < c_size; ++i) copy.data()[i] = data()[i];
            }
            return copy;
        }

        void resize(Size newSize) {
            reAllocate(newSize);
            c_size = newSize;
        }

        void reserve(Size newCap) {
            if (newCap > maxSize) reAllocate(newCap);
        }

        void swap(StaticVector &other) {
            std::swap(c_size, other.c_size);
            std::swap(maxSize, other.maxSize);
            data_.swap(other.data_);
        }

        // ：（）
        void free() {
            data_.reset();
            c_size = 0;
            maxSize = 0;
        }

        T &front() { return data()[0]; }
        const T &front() const { return data()[0]; }
        T &back() { return data()[c_size - 1]; }
        const T &back() const { return data()[c_size - 1]; }

        T &operator[](Size i) { return data()[i]; }
        const T &operator[](Size i) const { return data()[i]; }

        [[nodiscard]] Size size() const { return c_size; }

        // （）
        friend std::ostream &operator<<(std::ostream &os, const StaticVector &v) {
            if constexpr (std::is_same_v<T, bool>) os << std::boolalpha;
            os << "[";
            for (Size i = 0; i < v.c_size; ++i) {
                os << v.data()[i];
                if (i + 1 != v.c_size) os << ", ";
            }
            os << "]";
            return os;
        }

        bool operator==(const StaticVector &o) const {
            if (c_size != o.c_size || maxSize != o.maxSize) return false;
            if constexpr (std::is_trivially_copyable_v<T>) {
                return std::memcmp(data(), o.data(), c_size * sizeof(T)) == 0;
            } else {
                for (Size i = 0; i < c_size; ++i)
                    if (!(data()[i] == o.data()[i])) return false;
                return true;
            }
        }

        void print(const std::string &name = "") const {
            if constexpr (std::is_same_v<T, bool>) std::cout << std::boolalpha;
            if (!name.empty()) std::cout << name << ": ";
            std::cout << *this << '\n';
        }

        iterator begin() { return data(); }
        iterator end() { return data() + c_size; }
        const_iterator begin() const { return data(); }
        const_iterator end() const { return data() + c_size; }
        const_iterator cbegin() const { return data(); }
        const_iterator cend() const { return data() + c_size; }

    private:
        // 1.5x （）
        Size grow_capacity() const {
            return maxSize ? (maxSize + (maxSize >> 1) + 1) : Size(1);
        }

        void reAllocate(Size newCap) {
            if (newCap <= maxSize) return;
            std::shared_ptr<T[]> newData(new T[newCap], std::default_delete<T[]>());
            if constexpr (std::is_trivially_copyable_v<T>) {
                std::memcpy(newData.get(), data(), c_size * sizeof(T));
            } else {
                for (Size i = 0; i < c_size; ++i) newData[i] = std::move(data()[i]);
            }
            data_.swap(newData); // 
            maxSize = newCap;
        }
    };

    template<typename T>
    class MutexStaticVector {
    public:
        T *data;
        tbb::spin_mutex *mutexes;

        explicit MutexStaticVector(Size size) : data(new T[size]), mutexes(new tbb::spin_mutex[size]) {
        }

        void add(daf::Size index, T value) {
            tbb::spin_mutex::scoped_lock lock(mutexes[index]);
            data[index] += value;
        }

        void set(daf::Size index, T value) {
            tbb::spin_mutex::scoped_lock lock(mutexes[index]);
            data[index] = value;
        }

        bool setIfBigger(daf::Size index, T value) {
            tbb::spin_mutex::scoped_lock lock(mutexes[index]);
            if (data[index] < value) {
                data[index] = value;
                return true;
            }
            return false;
        }
    };

    template<typename Func, typename... Args>
    auto timeCount(const std::string &name, Func &&func, Args &&... args)
        -> decltype(std::forward<Func>(func)(std::forward<Args>(args)...)) {
        timespec start{}, end{};
        clock_gettime(CLOCK_MONOTONIC, &start);

        if constexpr (std::is_void_v<decltype(std::forward<Func>(func)(std::forward<Args>(args)...))>) {
            std::forward<Func>(func)(std::forward<Args>(args)...);
            clock_gettime(CLOCK_MONOTONIC, &end);

            double elapsed = (static_cast<double>(end.tv_sec) - static_cast<double>(start.tv_sec)) * 1e3
                             + (static_cast<double>(end.tv_nsec) - static_cast<double>(start.tv_nsec)) / 1e6;
            std::cout << name << " took: " << elapsed << " ms" << std::endl;
            return;
        } else {
            auto result = std::forward<Func>(func)(std::forward<Args>(args)...);
            clock_gettime(CLOCK_MONOTONIC, &end);

            double elapsed = (static_cast<double>(end.tv_sec) - static_cast<double>(start.tv_sec)) * 1e3
                             + (static_cast<double>(end.tv_nsec) - static_cast<double>(start.tv_nsec)) / 1e6;
            std::cout << name << " took: " << elapsed << " ms" << std::endl;

            return result;
        }
    }

    template<typename Func, typename... Args>
    auto timeCount(Func &&func, Args &&... args)
        -> decltype(timeCount(std::string(), std::forward<Func>(func), std::forward<Args>(args)...)) {
        return timeCount(std::string(), std::forward<Func>(func), std::forward<Args>(args)...);
    }

    template<typename T>
    bool printCandidate(StaticVector<T> &keep, StaticVector<T> &combination) {
        for (int x: keep) std::cout << x << " ";
        for (int x: combination) std::cout << x << " ";
        std::cout << "\n";
        return true;
    }

    template<typename T>
    concept ValidCallback = requires(std::vector<T> vec, std::function<bool(StaticVector<T> &, StaticVector<T> &)> f)
    {
        { f(vec) } -> std::same_as<bool>;
    };

    // input a function, the format of the function is:
    // bool function(StaticVector<T> &keep, StaticVector<T> &combination)
    // return true if continue to enumerate, false if stop to enumerate
    // template<
    //     typename Container,
    //     typename Callback,
    //     typename T = typename Container::value_type>
    // void enumerateCombinations(const Container &items,
    //                            size_t r,
    //                            Callback cb) {
    //     size_t n = items.size();
    //     if (r > n) return; //  r 
    //
    //     std::vector<T> combination;
    //     combination.reserve(r);
    //
    //     // ： items[start..n)  k 
    //     std::function<bool(size_t, size_t)> backtrack =
    //             [&](size_t start, size_t k) -> bool {
    //         if (k == 0) {
    //             // ，
    //             return cb(combination);
    //         }
    //         // ： k 
    //         if (n - start < k) {
    //             return true; // ，
    //         }
    //         // ： [start .. n-k] 
    //         for (size_t i = start; i + k <= n; ++i) {
    //             combination.push_back(items[i]);
    //             //  k-1 ， i+1 
    //             if (!backtrack(i + 1, k - 1)) {
    //                 return false; // 
    //             }
    //             combination.pop_back();
    //         }
    //         return true; // ，
    //     };
    //
    //     // 
    //     backtrack(0, r);
    // }

    template<class T, class U>
    concept EqualityComparable = requires(const T &a, const U &b)
    {
        { a == b } -> std::convertible_to<bool>;
    };

    template<
        typename S,
        typename R,
        typename Fn
    >
        requires EqualityComparable<typename S::value_type, typename R::value_type>

    auto enumerateCombinations(const S &keep,
                               const R &drop,
                               size_t r,
                               Fn cb) -> bool {
        using T = typename S::value_type;

        size_t keepCount = keep.size();
        size_t dropCount = drop.size();
        if (keepCount > r || keepCount + dropCount < r) {
            return true; // ，
        }

        size_t needDrop = r - keepCount;
        daf::StaticVector<T> comb(needDrop);
        comb.c_size = needDrop; // 
        // ， vector

        std::function<bool(size_t, size_t)> dfs =
                [&](size_t start, size_t choose) -> bool {
            if (choose == 0) {
                return cb(keep, comb);
            }
            // ：
            if (dropCount - start < choose) {
                return true;
            }
            for (size_t i = start; i + choose <= dropCount; ++i) {
                comb[needDrop - choose] = drop[i];
                if (!dfs(i + 1, choose - 1))
                    return false;
            }
            return true;
        };

        if (needDrop == 0) {
            cb(keep, comb);
        } else {
            dfs(0, needDrop);
        }
        comb.free(); //  dfs 
        return true;
    }

    template<typename Container, typename Fn>
    bool enumerateCombinations__(const Container &C, size_t size, Fn cb) {
        size_t n = C.size();
        if (size > n) return true;
        daf::StaticVector<typename Container::value_type> comb(size);
        comb.c_size = size; // 
        // dfs(pos, left):  C[pos..]  left 
        std::function<bool(size_t, size_t)> dfs = [&](size_t pos, size_t left) -> bool {
            if (left == 0) {
                if constexpr (std::is_void_v<
                    std::invoke_result_t<Fn, decltype(comb)>
                >) {
                    cb(comb);
                    return true;
                } else {
                    return cb(comb);
                }
            }
            // ：
            if (n - pos < left) return true;
            // ： [pos .. n-left] 
            for (size_t i = pos; i + left <= n; ++i) {
                comb[size - left] = C[i];
                if (!dfs(i + 1, left - 1)) return false;
            }
            return true;
        };

        bool cont = dfs(0, size);
        comb.free(); // ** dfs **
        return cont;
    }

    // template<class Container, class Fn>
    // void enumerateCombinations(const Container &C, std::size_t k, Fn &&cb) {
    //     const std::size_t n = C.size();
    //     if (k == 0 || k > n) return;
    //
    //     // index buffer
    //     daf::StaticVector<std::size_t> idx(k);
    //     // use shot vector to avoid heap allocation
    //     idx.c_size = k;
    //     for (std::size_t i = 0; i < k; ++i) idx[i] = i;
    //
    //     daf::StaticVector<typename Container::value_type> comb(k);
    //     comb.c_size = k;
    //
    //     auto emit = [&]() {
    //         for (std::size_t t = 0; t < k; ++t) comb[t] = C[idx[t]];
    //         if constexpr (std::is_same_v<void,
    //             decltype(std::forward<Fn>(cb)(comb))>)
    //             cb(comb);
    //         else
    //             if (!cb(comb)) return false;
    //         return true;
    //     };
    //     if (!emit()) return;
    //
    //     while (true) {
    //         std::size_t i = k;
    //         // find right-most idx[i] that can ++
    //         while (i > 0 && idx[i - 1] == n - k + i - 1) --i;
    //         if (i == 0) break; // finished
    //
    //         --i;
    //         ++idx[i];
    //         for (std::size_t j = i + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
    //         if (!emit()) break;
    //     }
    // }


    template<class Container, class Fn, size_t K>
    inline void enumerateCombinations_fixedK(const Container &C, Fn &&cb) {
        static_assert(K >= 1, "K must be >= 1");
        const size_t n = C.size();
        if (K > n) return;

        // 1)  idx、 comb（）
        std::array<size_t, K> idx;
        for (size_t i = 0; i < K; ++i) idx[i] = i;

        daf::StaticVector<typename Container::value_type> comb(K);
        comb.c_size = K; //  StaticVector ，

        auto emit = [&]() -> bool {
            //  comb
            for (size_t t = 0; t < K; ++t) comb[t] = C[idx[t]];
            if constexpr (std::is_same_v<void, decltype(std::forward<Fn>(cb)(comb))>) {
                cb(comb);
                return true;
            } else {
                return cb(comb);
            }
        };

        if (!emit()) return;

        while (true) {
            // 2) “”（ K ，）
            size_t i = K;
            while (i > 0 && idx[i - 1] == n - K + i - 1) --i;
            if (i == 0) break;
            --i;
            ++idx[i];
            for (size_t j = i + 1; j < K; ++j) idx[j] = idx[j - 1] + 1;
            if (!emit()) break;
        }
    }

    // ：
    template<class Container, class Fn>
    void enumerateCombinations(const Container &C, size_t k, Fn &&cb) {
        switch (k) {
            case 1: return enumerateCombinations_fixedK<Container, Fn, 1>(C, std::forward<Fn>(cb));
            case 2: return enumerateCombinations_fixedK<Container, Fn, 2>(C, std::forward<Fn>(cb));
            case 3: return enumerateCombinations_fixedK<Container, Fn, 3>(C, std::forward<Fn>(cb));
            case 4: return enumerateCombinations_fixedK<Container, Fn, 4>(C, std::forward<Fn>(cb));
            case 5: return enumerateCombinations_fixedK<Container, Fn, 5>(C, std::forward<Fn>(cb));
            case 6: return enumerateCombinations_fixedK<Container, Fn, 6>(C, std::forward<Fn>(cb));
            case 7: return enumerateCombinations_fixedK<Container, Fn, 7>(C, std::forward<Fn>(cb));
            case 8: return enumerateCombinations_fixedK<Container, Fn, 8>(C, std::forward<Fn>(cb));
            default: break; // 
        }

        // ===== （）=====
        const std::size_t n = C.size();
        if (k == 0 || k > n) return;

        daf::StaticVector<std::size_t> idx(k);
        idx.c_size = k;
        for (std::size_t i = 0; i < k; ++i) idx[i] = i;

        daf::StaticVector<typename Container::value_type> comb(k);
        comb.c_size = k;

        auto emit = [&]() {
            for (std::size_t t = 0; t < k; ++t) comb[t] = C[idx[t]];
            if constexpr (std::is_same_v<void, decltype(std::forward<Fn>(cb)(comb))>)
                cb(comb);
            else
                if (!cb(comb)) return false;
            return true;
        };
        if (!emit()) return;

        while (true) {
            std::size_t i = k;
            while (i > 0 && idx[i - 1] == n - k + i - 1) --i;
            if (i == 0) break;
            --i;
            ++idx[i];
            for (std::size_t j = i + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
            if (!emit()) break;
        }
    }


    template<typename Container, typename T, typename Fn>
    bool enumerateCombinations(const Container &C,
                               size_t size,
                               T *buf,
                               Fn cb) {
        size_t n = C.size();
        if (size > n) return true;

        std::function<bool(size_t, size_t)> dfs = [&](size_t pos, size_t left) -> bool {
            if (left == 0) {
                // return cb(buf, size);
                if constexpr (std::is_void_v<
                    std::invoke_result_t<Fn,
                        decltype(buf), // buf ， T*
                        decltype(size) // size_t
                    >
                >) {
                    cb(buf, size);
                    return true;
                } else {
                    return cb(buf, size);
                }
            }
            if (n - pos < left) return true;
            for (size_t i = pos; i + left <= n; ++i) {
                buf[size - left] = C[i];
                if (!dfs(i + 1, left - 1)) return false;
            }
            return true;
        };
        return dfs(0, size);
    }

    // -----------------------------------------------------------------------------
    //  S ∪ R  k  1  S 
    // S, R: ，S::value_type != R::value_type 
    // Fn:  bool(const Sval* bufS, size_t sCount, const Rval* bufR, size_t rCount)
    // -----------------------------------------------------------------------------
    template<typename S, typename R, typename Fn>
    bool enumAtLeastOneFromTwo(const S &s,
                               const R &r,
                               size_t k,
                               Fn cb) {
        size_t sz = s.size(), zr = r.size();
        assert(k > 0 && k <= sz + zr);

        using Sval = typename S::value_type;
        using Rval = typename R::value_type;

        // 
        std::vector<Sval> bufS(k);
        std::vector<Rval> bufR(k);

        // t = 1..min(k, sz)
        for (size_t t = 1; t <= std::min(k, sz); ++t) {
            bool contS = enumerateCombinations(s, t, bufS.data(), [&](Sval *sdata, size_t scount) {
                //  R  k-t
                return enumerateCombinations(r, k - t, bufR.data(), [&](Rval *rdata, size_t rcount) {
                    //  S  R 
                    return cb(sdata, scount, rdata, rcount);
                });
            });
            if (!contS) return false;
        }
        return true;
    }


    template<typename It1, typename It2, typename Func>
    inline void intersect_with_callback(It1 first1, It1 last1,
                                        It2 first2, It2 last2,
                                        Func &&f) noexcept {
        while (first1 != last1 && first2 != last2) {
            if (*first1 < *first2) {
                ++first1;
            } else if (*first2 < *first1) {
                ++first2;
            } else {
                f(*first1);
                ++first1;
                ++first2;
            }
        }
    }

    template<
        std::ranges::forward_range R1,
        std::ranges::forward_range R2,
        typename Func
    >
        requires std::equality_comparable_with<
                     std::ranges::range_value_t<R1>,
                     std::ranges::range_value_t<R2>
                 > &&
                 std::invocable<Func, std::ranges::range_value_t<R1> >
    inline void intersect_with_callback(const R1 &c1,
                                        const R2 &c2,
                                        Func &&f) noexcept {
        intersect_with_callback(
            std::begin(c1), std::end(c1),
            std::begin(c2), std::end(c2),
            std::forward<Func>(f)
        );
    }

    // ：
    // template<typename It1, typename It2, typename Func>
    // inline void intersect_with_callback_pair(It1 first1, It1 last1,
    //                                          It2 first2, It2 last2,
    //                                          Func &&f) noexcept {
    //     while (first1 != last1 && first2 != last2) {
    //         if (*first1 < *first2) {
    //             ++first1;
    //         } else if (*first2 < *first1) {
    //             ++first2;
    //         } else {
    //             // “” f
    //             f(*first1, *first2);
    //             ++first1;
    //             ++first2;
    //         }
    //     }
    // }

    // 
    template<
        std::ranges::forward_range R1,
        std::ranges::forward_range R2,
        typename Func
    >
        requires
        // ／
        std::totally_ordered_with<
            std::ranges::range_reference_t<R1>,
            std::ranges::range_reference_t<R2>
        > &&
        //  Func ，
        std::invocable<
            Func,
            std::ranges::range_reference_t<R1>,
            std::ranges::range_reference_t<R2>
        >
    inline void intersect_with_callback_pair(const R1 &c1,
                                             const R2 &c2,
                                             Func &&f) noexcept {
        auto first1 = std::begin(c1), last1 = std::end(c1);
        auto first2 = std::begin(c2), last2 = std::end(c2);
        while (first1 != last1 && first2 != last2) {
            if (*first1 < *first2) {
                ++first1;
            } else if (*first2 < *first1) {
                ++first2;
            } else {
                //  f  TreeGraphNode& 
                f(*first1, *first2);
                ++first1;
                ++first2;
            }
        }
    }

    inline void printProgress(std::size_t curr, std::size_t total) {
        constexpr int barWidth = 50;
        double ratio = static_cast<double>(curr) / total;
        int filled = static_cast<int>(ratio * barWidth);

        std::cout << '\r' // ，
                << '[';
        for (int i = 0; i < barWidth; ++i)
            std::cout << (i < filled ? '=' : ' ');
        std::cout << "] "
                << std::setw(6) << std::fixed << std::setprecision(2)
                << (ratio * 100.0) << "%  "
                << '(' << curr << '/' << total << ')'
                << std::flush; // 
    }

    template<typename T, typename Func>
    void intersect_dense_sets(const robin_hood::unordered_flat_set<T> &A,
                              const robin_hood::unordered_flat_set<T> &B,
                              Func callback) noexcept {
        // ， find 
        const auto *small = &A, *large = &B;
        bool swapped = false;
        if (B.size() < A.size()) {
            std::swap(small, large);
            swapped = true;
        }
        for (const auto &x: *small) {
            auto it = large->find(x);
            if (it != large->end()) {
                // x:  small； *it:  large
                // callback(x, *it);
                if (swapped) {
                    callback(*it, x); //  swapped， large 
                } else {
                    callback(x, *it); //  small 
                }
            }
        }
    }

    template<class T,
        class IndexRange, // ：std::vector<size_t>、std::span<size_t> 
        class Func>
    void intersect_dense_sets_multi(IndexRange &indices, // 
                                    std::vector<robin_hood::unordered_flat_set<T> > &adj_list,
                                    Func &&callback) // 
        noexcept {
        /* ---------- 0.  ---------- */
        if (indices.empty()) return;

        /* ---------- 1. “” ---------- */
        std::size_t best_idx = std::numeric_limits<std::size_t>::max();
        std::size_t best_sz = std::numeric_limits<std::size_t>::max();

        for (auto id: indices) {
            const auto &S = adj_list[id];
            if (S.size() < best_sz) {
                best_sz = S.size();
                best_idx = id;
            }
        }
        if (best_sz == 0) return; //  ⇒ 

        const auto &base = adj_list[best_idx];

        /* ---------- 2.  ---------- */
        std::vector<const robin_hood::unordered_flat_set<T> *> others;
        others.reserve(indices.size() - 1);
        for (auto id: indices)
            if (id != best_idx) others.push_back(&adj_list[id]);

        /* ---------- 3.  base， ---------- */
        for (const auto &x: base) {
            bool ok = true;
            for (const auto *pSet: others) {
                if (pSet->find(x) == pSet->end()) {
                    ok = false;
                    break;
                }
            }
            if (ok) callback(x); // x 
        }
    }

    extern daf::StaticVector<daf::Size> vListMap;
    extern daf::StaticVector<daf::Size> globalCSR;
}

struct TreeGraphNode {
    daf::Size v: 63;
    daf::Size isPivot: 1;

    constexpr TreeGraphNode() = default;

    constexpr TreeGraphNode(daf::Size v, bool isPivot) : v(v), isPivot(isPivot) {
    }

    constexpr explicit TreeGraphNode(daf::Size v) : v(v), isPivot(false) {
    }

    constexpr TreeGraphNode(daf::Size v, daf::Size isPivot) : v(v), isPivot(isPivot) {
    }

    // <<
    friend std::ostream &operator<<(std::ostream &os, const TreeGraphNode &node) {
        os << "(" << node.v << ", ";
        if (node.isPivot) {
            os << "\033[31mDrop\033[0m";
        } else {
            os << "\033[32mKeep\033[0m";
        }
        os << ")";
        return os;
    }

    operator daf::Size() const noexcept { return v; }
    // == operater with unit64
    constexpr bool operator==(const TreeGraphNode &other) const {
        return v == other.v;
    }

    bool operator!=(const TreeGraphNode &other) const {
        return !(*this == other);
    }

    bool operator<(const TreeGraphNode &other) const {
        return v < other.v;
    }

    bool operator>(const TreeGraphNode &other) const {
        return v > other.v;
    }


    constexpr bool operator==(const daf::Size other) const {
        return v == other;
    }


    static const TreeGraphNode EMPTYKEY;
    static const TreeGraphNode DELETEDKEY;
};


inline const TreeGraphNode TreeGraphNode::EMPTYKEY{
    (daf::Size(1) << 63) - 1, false
};
inline const TreeGraphNode TreeGraphNode::DELETEDKEY{
    (daf::Size(1) << 63) - 2, false
};

template<>
struct std::hash<TreeGraphNode> {
    size_t operator()(const TreeGraphNode &node) const noexcept {
        return std::hash<daf::Size>()(node.v);
    }
};

template<>
struct std::hash<daf::StaticVector<daf::Size> > {
    std::size_t operator()(const daf::StaticVector<daf::Size> &c) const noexcept {
        std::size_t seed = 0;
        for (const daf::Size &v: c) {
            seed ^= std::hash<daf::Size>()(v) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

#endif //SUBGRAPHMATCHING_GLOBAL_H