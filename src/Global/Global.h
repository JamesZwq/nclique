//
// Created by 张文谦 on 24-7-29.
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
// #include <google/dense_hash_set>

#include <iomanip>
#include <random>

#include "dataStruct/robin_hood.h"


#define MAX_CSIZE 400

template<typename U, typename V>
std::ostream &operator<<(std::ostream &os, const std::pair<U, V> &p);

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec);

template<typename Key, typename Value, typename Compare>
std::ostream &operator<<(std::ostream &os, const std::map<Key, Value, Compare> &m);


template<typename Iterator>
std::ostream &operator<<(std::ostream &os, const std::ranges::subrange<Iterator> &subrange) {
    os << "[";
    bool first = true;
    for (auto it = subrange.begin(); it != subrange.end(); ++it) {
        if (!first) {
            os << ", ";
        }
        os << *it;
        first = false;
    }
    os << "]";
    return os;
}

template<typename U, typename V>
std::ostream &operator<<(std::ostream &os, const std::pair<U, V> &p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}


template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
    os << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i != vec.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

template<typename T>
std::ostream &operator<<(std::ostream & lhs, const std::span<T> & rhs) {
    lhs << "[";
    for (size_t i = 0; i < rhs.size(); ++i) {
        lhs << rhs[i];
        if (i != rhs.size() - 1) {
            lhs << ", ";
        }
    }
    lhs << "]";
    return lhs;
}

template<typename Key, typename Value, typename Compare>
std::ostream &operator<<(std::ostream &os, const std::map<Key, Value, Compare> &m) {
    os << "{";
    bool first = true;
    for (const auto &pair: m) {
        if (!first) {
            os << ", ";
        }
        os << pair.first << ": " << pair.second;
        first = false;
    }
    os << "}";
    return os;
}

// 全局，放一个公共 .h/.hpp 里
// template<typename T>
// std::ostream &operator<<(std::ostream &os, const google::dense_hash_set<T> &s) {
//     os << "{";
//     bool first = true;
//     for (const auto &x: s) {
//         if (!first) os << ", ";
//         first = false;
//         os << x;
//     }
//     os << "}";
//     return os;
// }

namespace daf {
    using Size = uint64_t;
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


    template<typename T>
    class StaticVector {
        static_assert(std::is_copy_constructible_v<T>, "StaticVector requires T to be trivially copyable");

    public:
        typedef T *iterator;
        typedef const T *const_iterator;
        using value_type = T;

        daf::Size c_size;
        daf::Size maxSize;
        T *data;


        explicit StaticVector(const daf::Size maxSize = MAX_CSIZE) : c_size(0), maxSize(maxSize + 1),
                                                                     data(new T[maxSize]) {
        }

        // copy operator
        StaticVector(const StaticVector &other) : c_size(other.c_size), maxSize(other.maxSize), data(other.data) {
        }

        StaticVector &operator=(const StaticVector &other) {
            if (this == &other) {
                return *this;
            }
            c_size = other.c_size;
            maxSize = other.maxSize;
            data = other.data;
            // std::cout << "copy constructor" << std::endl;
            return *this;
        }

        bool operator==(const StaticVector &vertices) const {
            if (c_size != vertices.c_size || maxSize != vertices.maxSize) {
                return false;
            }
            // for (daf::Size i = 0; i < c_size; i++) {
            //     if (data[i] != vertices.data[i]) {
            //         return false;
            //     }
            // }
            return std::memcmp(data, vertices.data, c_size * sizeof(T)) == 0;
        }


        StaticVector(StaticVector &&other) noexcept : c_size(other.c_size), data(other.data), maxSize(other.maxSize) {
            other.data = nullptr;
            other.c_size = 0;
        }

        explicit operator T *() {
            return data;
        }


        void push_back_with_check(const T &value) {
            if (c_size >= maxSize) {
                reAllocate(maxSize * 1.5);
                // std::cout << "Reallocating memory for StaticVector, current size: " << c_size
                //           << ", new max size: " << maxSize << std::endl;
            }
            data[c_size++] = value;
        }


        void push_back(const T &value) {
            if (c_size >= maxSize) {
                reAllocate(maxSize * 1.5);
                // std::cout << "Reallocating memory for StaticVector, current size: " << c_size
                //           << ", new max size: " << maxSize << std::endl;
            }
            data[c_size++] = value;
        }

        template<typename... Args>
        void emplace_back(Args &&... args) {
            // 对于原始内存，可用 placement new；如果 data 已经是 T 类型数组直接赋值也行
            if (c_size >= maxSize) {
                reAllocate(maxSize * 1.5);
            }
            new(&data[c_size]) T(std::forward<Args>(args)...);
            ++c_size;
        }

        void pop_back() {
            c_size--;
        }

        // erase
        void remove(daf::Size index) {
#ifndef NDEBUG
            if (index >= c_size) {
                throw std::out_of_range("Index out of range.");
            }
#endif
            data[index] = data[--c_size];
        }

        void extend(const T *begin, daf::Size size) {
            std::memcpy(data + c_size, begin, size * sizeof(T));
            c_size += size;
        }

        template<typename BinaryPred = std::equal_to<T> >
        void unique(BinaryPred comp = BinaryPred()) noexcept {
            if (c_size < 2) return;
            daf::Size write = 1;
            for (daf::Size read = 1; read < c_size; ++read) {
                if (!comp(data[read], data[write - 1])) {
                    data[write++] = data[read];
                }
            }
            c_size = write;
        }

        [[nodiscard]] bool empty() const {
            return c_size == 0;
        }

        void clear() {
            c_size = 0;
        }

        [[nodiscard]] StaticVector deepCopy() const {
            StaticVector copy(maxSize);
            copy.c_size = c_size;
            std::memcpy(copy.data, data, c_size * sizeof(T));
            return copy;
        }

        void resize(Size newSize) {
            reAllocate(newSize);
            c_size = newSize;
        }

        // reserve
        void reserve(Size newSize) {
            if (newSize > maxSize) {
                T *newData = new T[newSize];
                std::memcpy(newData, data, c_size * sizeof(T));
                delete[] data;
                data = newData;
                maxSize = newSize;
            }
        }


        void swap(StaticVector &other) {
            std::swap(c_size, other.c_size);
            std::swap(maxSize, other.maxSize);
            std::swap(data, other.data);
        }

        void free() const {
            delete[] data;
        }

        T &front() {
            return data[0];
        }

        [[nodiscard]] const T &front() const {
            return data[0];
        }

        T &back() {
            return data[c_size - 1];
        }

        [[nodiscard]] const T &back() const {
            return data[c_size - 1];
        }

        T &operator[](daf::Size index) {
#if DEBUG
            if (index >= c_size) {
                throw std::out_of_range("Index out of range.");
            }
#endif
            return data[index];
        }

        const T &operator[](daf::Size index) const {
#if DEBUG
            if (index >= c_size) {
                throw std::out_of_range("Index out of range.");
            }
#endif
            return data[index];
        }

        [[nodiscard]] daf::Size size() const {
            return c_size;
        }

        friend std::ostream &operator<<(std::ostream &os, const StaticVector &vector) {
            // if type is bool
            if constexpr (std::is_same_v<T, bool>) {
                os << std::boolalpha;
            }
            os << "[";
            for (daf::Size i = 0; i < vector.c_size; i++) {
                os << vector.data[i];
                if (i != vector.c_size - 1) {
                    os << ", ";
                }
            }
            os << "]";
            return os;
        }

        void print(const std::string &name = "") const {
            if constexpr (std::is_same_v<T, bool>) {
                std::cout << std::boolalpha;
            }
            std::cout << name << ": [";
            for (daf::Size i = 0; i < c_size; i++) {
                std::cout << data[i];
                if (i != c_size - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
        }

        [[nodiscard]] T *getData() const {
            return data;
        }

        // 添加迭代器接口，支持范围 for 循环
        iterator begin() {
            return data;
        }

        iterator end() {
            return data + c_size;
        }

        [[nodiscard]] const_iterator begin() const {
            return data;
        }

        [[nodiscard]] const_iterator end() const {
            return data + c_size;
        }

        [[nodiscard]] const_iterator cbegin() const {
            return data;
        }

        [[nodiscard]] const_iterator cend() const {
            return data + c_size;
        }

    private:
        // void resize(Size newSize) {
        //     if (newSize <= maxSize) {
        //         // 如果只是缩小，直接改变当前大小即可
        //         c_size = newSize;
        //         return;
        //     }
        //     // 扩容：申请新内存，注意先申请再 delete[]，保证异常安全
        //     T* newData = new T[newSize];
        //     // 拷贝旧数据——如果 T 是 POD 类型，可换成 std::memcpy；否则用 std::move 或 std::copy
        //     if constexpr(std::is_trivially_copyable<T>::value) {
        //         std::memcpy(newData, data, sizeof(T) * c_size);
        //     } else {
        //         for (Size i = 0; i < c_size; ++i) {
        //             newData[i] = std::move(data[i]);
        //         }
        //     }
        //     // 释放旧内存
        //     delete[] data;
        //     // 指针、容量、当前大小更新
        //     data    = newData;
        //     maxSize = newSize;
        //     c_size  = newSize;
        // }
        void reAllocate(Size newSize) {
            if (newSize > maxSize) {
                T *newData = new T[newSize];
                if constexpr (std::is_trivially_copyable<T>::value) {
                    std::memcpy(newData, data, sizeof(T) * c_size);
                } else {
                    for (Size i = 0; i < c_size; ++i) {
                        newData[i] = std::move(data[i]);
                    }
                }
                delete[] data;
                data = newData;
                maxSize = newSize;
            }
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
    //     if (r > n) return; // 不足以选出 r 个
    //
    //     std::vector<T> combination;
    //     combination.reserve(r);
    //
    //     // 用递归回溯：从 items[start..n) 中选 k 个
    //     std::function<bool(size_t, size_t)> backtrack =
    //             [&](size_t start, size_t k) -> bool {
    //         if (k == 0) {
    //             // 选够了，调用回调
    //             return cb(combination);
    //         }
    //         // 剪枝：剩余元素不足 k 个
    //         if (n - start < k) {
    //             return true; // 这一支不够，回到上一层继续
    //         }
    //         // 枚举：在 [start .. n-k] 范围内选一个
    //         for (size_t i = start; i + k <= n; ++i) {
    //             combination.push_back(items[i]);
    //             // 还需选 k-1 个，从 i+1 开始
    //             if (!backtrack(i + 1, k - 1)) {
    //                 return false; // 回调请求提前终止
    //             }
    //             combination.pop_back();
    //         }
    //         return true; // 本层枚举完毕，继续上一层
    //     };
    //
    //     // 触发回溯
    //     backtrack(0, r);
    // }

    template<class T, class U>
concept EqualityComparable = requires(const T& a, const U& b) {
        { a == b } -> std::convertible_to<bool>;
};

    template<
        typename S,
        typename R,
        typename Fn
    >
        requires EqualityComparable<typename S::value_type, typename R::value_type>

    auto enumerateCombinations(const S& keep,
                               const R& drop,
                               size_t r,
                               Fn cb) -> bool {
        using T = typename S::value_type;

        size_t keepCount = keep.size();
        size_t dropCount = drop.size();
        if (keepCount > r || keepCount + dropCount < r) {
            return true; // 条件不满足，直接跳过
        }

        size_t needDrop = r - keepCount;
        daf::StaticVector<T> comb(needDrop);
        comb.c_size = needDrop; // 设置当前大小
        // 为了支持任意可迭代容器，先拷贝至 vector

        std::function<bool(size_t, size_t)> dfs =
            [&](size_t start, size_t choose) -> bool {
                if (choose == 0) {
                    return cb(keep, comb);
                }
                // 剪枝：剩余元素不足
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
        comb.free(); // 仅在 dfs 完全结束后释放
        return true;
    }

    template<typename Container, typename Fn>
    bool enumerateCombinations(const Container &C, size_t size, Fn cb) {
        size_t n = C.size();
        if (size > n) return true;
        daf::StaticVector<typename Container::value_type> comb(size);
        comb.c_size = size; // 设置当前大小
        // dfs(pos, left): 从 C[pos..] 还要选 left 个
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
            // 剪枝：剩余元素不足
            if (n - pos < left) return true;
            // 枚举：在 [pos .. n-left] 中取一个
            for (size_t i = pos; i + left <= n; ++i) {
                comb[size - left] = C[i];
                if (!dfs(i + 1, left - 1)) return false;
            }
            return true;
        };

        bool cont = dfs(0, size);
        comb.free();          // **仅在 dfs 完全结束后释放**
        return cont;
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
                if constexpr ( std::is_void_v<
                    std::invoke_result_t<Fn,
                        decltype(buf),    // buf 的指针类型，比如 T*
                        decltype(size)    // size_t
                    >
                > ) {
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
    // 从 S ∪ R 中枚举所有大小为 k 且至少包含 1 个来自 S 的组合
    // S, R: 随机访问容器，S::value_type != R::value_type 也可
    // Fn: 回调 bool(const Sval* bufS, size_t sCount, const Rval* bufR, size_t rCount)
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

        // 预分配组合缓冲区
        std::vector<Sval> bufS(k);
        std::vector<Rval> bufR(k);

        // t = 1..min(k, sz)
        for (size_t t = 1; t <= std::min(k, sz); ++t) {
            bool contS = enumerateCombinations(s, t, bufS.data(), [&](Sval *sdata, size_t scount) {
                // 再从 R 中选 k-t
                return enumerateCombinations(r, k - t, bufR.data(), [&](Rval *rdata, size_t rcount) {
                    // 回调分别传入 S 部分和 R 部分
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

    // 迭代器版：二参回调
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
    //             // 这里把两个“相等”但来自不同容器的引用都传给 f
    //             f(*first1, *first2);
    //             ++first1;
    //             ++first2;
    //         }
    //     }
    // }

    // 容器版转发
    template<
        std::ranges::forward_range R1,
        std::ranges::forward_range R2,
        typename Func
    >
        requires
        // 确保两个引用能比较大小／相等
        std::totally_ordered_with<
            std::ranges::range_reference_t<R1>,
            std::ranges::range_reference_t<R2>
        > &&
        // 确保 Func 可以被调用，且参数是两个引用
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
                // 这里 f 拿到的就是 TreeGraphNode& 而非拷贝
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

        std::cout << '\r' // 回到行首，覆盖上一帧
                << '[';
        for (int i = 0; i < barWidth; ++i)
            std::cout << (i < filled ? '=' : ' ');
        std::cout << "] "
                << std::setw(6) << std::fixed << std::setprecision(2)
                << (ratio * 100.0) << "%  "
                << '(' << curr << '/' << total << ')'
                << std::flush; // 立即刷新
    }

    template<typename T, typename Func>
    void intersect_dense_sets(const robin_hood::unordered_flat_set<T> &A,
                              const robin_hood::unordered_flat_set<T> &B,
                              Func callback) noexcept {
        // 先挑小的遍历，保证 find 次数最少
        const auto *small = &A, *large = &B;
        bool swapped = false;
        if (B.size() < A.size()) {
            std::swap(small, large);
            swapped = true;
        }
        for (const auto &x: *small) {
            auto it = large->find(x);
            if (it != large->end()) {
                // x: 来自 small； *it: 来自 large
                // callback(x, *it);
                if (swapped) {
                    callback(*it, x); // 如果 swapped，先是 large 的元素
                } else {
                    callback(x, *it); // 否则先是 small 的元素
                }
            }
        }
    }

    extern daf::StaticVector<daf::Size> vListMap;


    inline void StaticVectorTest() {
        /* ===== 参数 ===== */
        constexpr std::size_t N = 5'000'000; // 若内存足够可调大
        constexpr int ROUNDS = 3;

        /* ===== 生成一个随机索引序列，用于“随机写入” ===== */
        std::vector<std::size_t> idx(N);
        std::iota(idx.begin(), idx.end(), 0);
        std::mt19937_64 rng(42);
        std::shuffle(idx.begin(), idx.end(), rng);

        auto ms = [](auto ns) { return ns.count() * 1e-6; };
        using clock = std::chrono::high_resolution_clock;

        auto bench = [&](auto &&container, const char *name) {
            std::chrono::nanoseconds best_seq_write{};
            std::chrono::nanoseconds best_rand_write{};
            std::chrono::nanoseconds best_iter{};

            best_seq_write = best_rand_write = best_iter = std::chrono::nanoseconds::max();

            for (int r = 0; r < ROUNDS; ++r) {
                /* ---------- 顺序写 ---------- */
                auto t0 = clock::now();
                for (std::size_t i = 0; i < N; ++i) container[i] = static_cast<int>(i);
                auto t1 = clock::now();
                best_seq_write = std::min(best_seq_write, t1 - t0);

                /* ---------- 随机写 ---------- */
                auto t2 = clock::now();
                for (std::size_t i = 0; i < N; ++i) container[idx[i]] = static_cast<int>(i);
                auto t3 = clock::now();
                best_rand_write = std::min(best_rand_write, t3 - t2);

                /* ---------- 顺序遍历 ---------- */
                volatile long long sink = 0;
                auto t4 = clock::now();
                for (std::size_t i = 0; i < N; ++i) sink += container[i];
                auto t5 = clock::now();
                best_iter = std::min(best_iter, t5 - t4);
            }

            std::cout << name
                    << "  seq-write: " << ms(best_seq_write)
                    << " ms  | rand-write: " << ms(best_rand_write)
                    << " ms  | iterate: " << ms(best_iter) << " ms\n";
        };

        /* ===== 1) std::vector ===== */
        std::vector<int> stdVec(N);
        bench(stdVec, "std::vector");

        /* ===== 2) StaticVector ===== */
        StaticVector<int> staticVec(N);
        staticVec.resize(N); // 保证 c_size==N，直接用 operator[]
        bench(staticVec, "StaticVector");

        /* ===== 3) C-array ===== */
        auto *rawArr = new int[N];
        bench(rawArr, "C-array");
        delete[] rawArr;

        /* ===== 4) std::array ===== */
        static std::array<int, N> stdArr; // 静态放 data 段；或改成 thread_local
        bench(stdArr, "std::array");
    }
}

struct TreeGraphNode {
    uint64_t v: 63;
    uint64_t isPivot: 1;

    constexpr TreeGraphNode() = default;

    constexpr TreeGraphNode(uint64_t v, bool isPivot) : v(v), isPivot(isPivot) {
    }

    constexpr explicit TreeGraphNode(uint64_t v) : v(v), isPivot(false) {
    }

    constexpr TreeGraphNode(uint64_t v, uint64_t isPivot) : v(v), isPivot(isPivot) {
    }

    // <<
    friend std::ostream &operator<<(std::ostream &os, const TreeGraphNode &node) {
        // os << "(" << node.v << ", " << node.isPivot << ")";
        if (node.isPivot) {
            os << "(" << node.v << ", Drop)";
        } else {
            os << "(" << node.v << ", Keep)";
        }
        return os;
    }

    operator uint64_t() const noexcept { return v; }
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
    (uint64_t(1) << 63) - 1, false
};
inline const TreeGraphNode TreeGraphNode::DELETEDKEY{
    (uint64_t(1) << 63) - 2, false
};

template<>
struct std::hash<TreeGraphNode> {
    size_t operator()(const TreeGraphNode &node) const noexcept {
        return std::hash<uint64_t>()(node.v);
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
