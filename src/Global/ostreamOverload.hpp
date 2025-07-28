//
// Created by 张文谦 on 25-7-27.
//

#ifndef OSTREAMOVERLOAD_HPP
#define OSTREAMOVERLOAD_HPP
#include <vector>
#include <iostream>
#include <numeric>
#include <ranges>
#include <map>
#include <span>

#include <random>
#include <type_traits>

#include "dataStruct/robin_hood.h"


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


// template<typename T>
// std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
//     os << "vec size: "<< vec.size() << " [";
//     for (size_t i = 0; i < vec.size(); ++i) {
//         os << vec[i];
//         if (i != vec.size() - 1) {
//             os << ", ";
//         }
//     }
//     os << "]";
//     return os;
// }
template<class T>
concept ostreamable =
    requires(std::ostream& os, T const& v) { os << v; };

// ─── 通用 Range 可打印概念 ────────────────────────────────
template<class R>
concept printable_range =
       std::ranges::input_range<R>                                               // 可遍历
    && requires(std::ostream& os, std::ranges::range_reference_t<R> ref) {       // 每个元素能被 << 输出
           os << ref;
       }
    // 排除字符串 / 字符数组 / C‑string
    && !std::same_as<std::remove_cvref_t<R>, std::string>
    && !std::same_as<std::remove_cvref_t<R>, std::string_view>
    && !(std::is_pointer_v<std::remove_cvref_t<R>> &&
         std::is_same_v<std::remove_cvref_t<std::remove_pointer_t<R>>, char>)
    && !(std::is_array_v<std::remove_cvref_t<R>> &&
         std::is_same_v<std::remove_all_extents_t<std::remove_cvref_t<R>>, char>);

// ─── 真正的通用输出重载 ───────────────────────────────────────────
template<printable_range R>
std::ostream& operator<<(std::ostream& os, R&& r)
{
    using std::ranges::begin, std::ranges::end, std::ranges::distance;

    bool isClique = false;
    // 检查范围内的元素是否为整数类型，并且范围长度大于1
    if constexpr (std::is_integral_v<std::ranges::range_value_t<R>>) {
        auto n = distance(r);
        if (n > 1) {
            isClique = true;
        }
    }

    if (isClique) {
        os << "[";
        bool first = true;
        for (auto&& elem : r) {
            if (!first) os << ", ";
            first = false;
            os << elem;
        }
        return os << "]";
    } else {
        os << "vec size: " << distance(r) << " [";
        bool first = true;
        for (auto&& elem : r) {
            if (!first) os << ", ";
            first = false;
            os << elem;
        }
        return os << "]";
    }
}
// Make the generic << visible to std::ranges views (transform_view, filter_view …)
namespace std::ranges {
    using ::operator<<;
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


inline void debug_break_impl() noexcept
{
#if defined(_MSC_VER)            // MSVC / clang-cl
    __debugbreak();              // 直接调用 IDE 自带的软中断
#elif defined(__clang__) || defined(__GNUC__)
#   if __has_builtin(__builtin_debugtrap)
    __builtin_debugtrap();   // clang / gcc 10+
#   elif (defined(__i386__) || defined(__x86_64__))
    __asm__ volatile("int3"); // x86 上最简单
#   else
    std::raise(SIGTRAP);     // 其它架构退而求其次
#   endif
#else
    std::raise(SIGTRAP);         // C 标准库方案，几乎处处可用
#endif
}

/* 只要把表达式写进去，满足时就停：
 *   DEBUG_BREAK_IF(ptr == nullptr);
 */
#define DEBUG_BREAK_IF(expr)     \
do {                         \
if (__builtin_expect(!!(expr), 0)) \
debug_break_impl();  \
} while (false)

/* 像 assert 一样用，但失败时不是 abort 而是断在原地 */
#define ASSERT_BREAK(expr)       \
DEBUG_BREAK_IF(!(expr))
#endif //OSTREAMOVERLOAD_HPP
