//
// Created by _ on 25-7-27.
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

// ostreamable concept must come before any usage in overloads
template<class T>
concept ostreamable =
    requires(std::ostream& os, T const& v) { os << v; };




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

// ───  Range  ────────────────────────────────
template<typename R>
auto operator<<(std::ostream& os, R const& r)
    -> std::enable_if_t<
           std::ranges::input_range<R> &&
           !std::same_as<std::remove_cvref_t<R>, std::string> &&
           !std::same_as<std::remove_cvref_t<R>, std::string_view> &&
           !(std::is_pointer_v<std::remove_cvref_t<R>> &&
             std::is_same_v<std::remove_pointer_t<std::remove_cvref_t<R>>, char>) &&
           !(std::is_array_v<std::remove_cvref_t<R>> &&
             std::is_same_v<std::remove_all_extents_t<std::remove_cvref_t<R>>, char>),
           std::ostream&
       >
{
    using std::ranges::begin;
    using std::ranges::end;
    using std::ranges::size;

    std::size_t n;
    if constexpr (std::ranges::sized_range<R>) {
        n = size(r);
    } else {
        n = 0;
        for (auto it = begin(r); it != end(r); ++it) ++n;
    }

    constexpr bool elem_is_int =
        std::integral<std::remove_cvref_t<std::ranges::range_value_t<R>>>;
    bool isClique = elem_is_int && n > 1;

    if (isClique) {
        os << "[";
    } else {
        os << "vec size: " << n << " [";
    }

    bool first = true;
    for (auto&& e : r) {
        if (!first) os << ", ";
        first = false;
        os << e;
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


// ， .h/.hpp 
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
    __debugbreak();              //  IDE 
#elif defined(__clang__) || defined(__GNUC__)
#   if __has_builtin(__builtin_debugtrap)
    __builtin_debugtrap();   // clang / gcc 10+
#   elif (defined(__i386__) || defined(__x86_64__))
    __asm__ volatile("int3"); // x86 
#   else
    std::raise(SIGTRAP);     // 
#   endif
#else
    std::raise(SIGTRAP);         // C ，
#endif
}

/* ，：
 *   DEBUG_BREAK_IF(ptr == nullptr);
 */
#define DEBUG_BREAK_IF(expr)     \
do {                         \
if (__builtin_expect(!!(expr), 0)) \
debug_break_impl();  \
} while (false)

/*  assert ， abort  */
#define ASSERT_BREAK(expr)       \
DEBUG_BREAK_IF(!(expr))
#endif //OSTREAMOVERLOAD_HPP
