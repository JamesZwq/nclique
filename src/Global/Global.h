//
// Created by 张文谦 on 24-7-29.
//
#pragma once
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
#include <tbb/spin_mutex.h>


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


        explicit StaticVector(const daf::Size maxSize = MAX_CSIZE) : c_size(0), maxSize(maxSize), data(new T[maxSize]) {
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


        void push_back(const T &value) {
            data[c_size++] = value;
        }

        template <typename... Args>
        void emplace_back(Args&&... args) {
            // 对于原始内存，可用 placement new；如果 data 已经是 T 类型数组直接赋值也行
            new (&data[c_size]) T(std::forward<Args>(args)...);
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

        template<typename BinaryPred = std::equal_to<T>>
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
            delete[] data;
            data = new T[newSize];
            maxSize = newSize;
            c_size = newSize;
        }
        void free() const {
            delete[] data;
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
            std::cout << name <<": [";
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
    template<typename T,
        typename Callback = bool(*)(StaticVector<T> &, StaticVector<T> &)>
    void enumerateCombinations(StaticVector<T> &keep, StaticVector<T> &drop,
                               CliqueSize r,
                               Callback cb = printCandidate<T>) {
        // 检查输入条件是否满足
        if (keep.size() > r) {
            return;
        }
        if (keep.size() + drop.size() < r) {
            return;
        }

        // 需要从 drop 数组中选择的额外元素数量
        const CliqueSize needDrop = r - keep.size();
        StaticVector<T> combination(needDrop); // 用于存储当前组合

        // 定义回溯函数：start 表示从 drop 的哪个位置开始选，choose 表示还需要选几个元素
        std::function<bool(int, int)> backtrack = [&](int start, int choose) {
            if (choose == 0) {
                return cb(keep, combination);
            }
            // 剪枝：保证剩余元素足够选出需要的个数
            for (int i = start; i <= static_cast<int>(drop.size()) - choose; ++i) {
                combination.push_back(drop[i]);
                if (backtrack(i + 1, choose - 1)) {
                    combination.pop_back(); // 回溯：撤销选择
                } else {
                    combination.pop_back();
                    return false;
                }
            }
            return true;
        };

        // 如果无需从 drop 中选择任何元素，则直接调用回调函数处理 keep 数组
        if (needDrop == 0) {
            cb(keep, combination);
        } else {
            backtrack(0, needDrop);
        }

        combination.free();
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



    extern daf::StaticVector<daf::Size> vListMap;
}

template<>
struct std::hash<daf::StaticVector<daf::Size>> {
    std::size_t operator()(const daf::StaticVector<daf::Size> &c) const noexcept {
        std::size_t seed = 0;
        for (const daf::Size &v: c) {
            seed ^= std::hash<daf::Size>()(v) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

#endif //SUBGRAPHMATCHING_GLOBAL_H
