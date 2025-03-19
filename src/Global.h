//
// Created by 张文谦 on 24-7-29.
//

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

#define MAX_CSIZE 400

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::pair<T, T> &pair);

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec);

template<typename Key, typename Value, typename Compare>
std::ostream &operator<<(std::ostream &os, const std::map<Key, Value, Compare> &m);


template<typename T>
std::ostream &operator<<(std::ostream &os, const std::pair<T, T> &pair) {
    os << "(" << pair.first << ", " << pair.second << ")";
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
    for (const auto &pair : m) {
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

    // class Fraction {
    //     void normalize() {
    //         if (denominator == 0) {
    //             throw std::invalid_argument("Denominator cannot be zero.");
    //         }
    //         int gcd = std::gcd(numerator, denominator);
    //         numerator /= gcd;
    //         denominator /= gcd;
    //         if (denominator < 0) {
    //             // 确保分母为正
    //             numerator = -numerator;
    //             denominator = -denominator;
    //         }
    //     }
    //
    // public:
    //     Size numerator; // 分子
    //     Size denominator; // 分母
    //     // int f_value;
    //     // 构造函数
    //     explicit Fraction(const Size num = 0, const Size denom = 1) : numerator(num), denominator(denom) {
    //     }
    //
    //     // copy constructor
    //     Fraction(const Fraction &other) = default;
    //
    //
    //     [[nodiscard]] Size get_min_Numer(const Size newDeno) const {
    //         return static_cast<Size>(std::ceil(static_cast<double>(numerator) / denominator * newDeno - 1e-9));
    //     }
    //
    //     // 获取分子和分母
    //     [[nodiscard]] Size getNumerator() const { return numerator; }
    //     [[nodiscard]] Size getDenominator() const { return denominator; }
    //
    //     // 转换为浮点数
    //     [[nodiscard]] double toDouble() const {
    //         return denominator == 0 ? 0.0 : static_cast<double>(numerator) / denominator;
    //     }
    //
    //     // 运算符重载
    //     Fraction operator+(const Fraction &other) const {
    //         const Size num = numerator * other.denominator + other.numerator * denominator;
    //         const Size denom = denominator * other.denominator;
    //         return Fraction(num, denom);
    //     }
    //
    //     Fraction operator-(const Fraction &other) const {
    //         const Size num = numerator * other.denominator - other.numerator * denominator;
    //         const Size denom = denominator * other.denominator;
    //         return Fraction(num, denom);
    //     }
    //
    //     Fraction operator*(const Fraction &other) const {
    //         return Fraction(numerator * other.numerator, denominator * other.denominator);
    //     }
    //
    //     Fraction operator/(const Fraction &other) const {
    //         if (other.numerator == 0) {
    //             throw std::invalid_argument("Division by zero.");
    //         }
    //         return Fraction(numerator * other.denominator, denominator * other.numerator);
    //     }
    //
    //     // 比较运算符
    //     bool operator==(const Fraction &other) const {
    //         return numerator * other.denominator == other.numerator * denominator;
    //     }
    //
    //     bool operator!=(const Fraction &other) const {
    //         return !(*this == other);
    //     }
    //
    //     bool operator<(const Fraction &other) const {
    //         return numerator * other.denominator < other.numerator * denominator;
    //     }
    //
    //     bool operator<=(const Fraction &other) const {
    //         return numerator * other.denominator <= other.numerator * denominator;
    //     }
    //
    //     bool operator>(const Fraction &other) const {
    //         return numerator * other.denominator > other.numerator * denominator;
    //     }
    //
    //     bool operator>=(const Fraction &other) const {
    //         return numerator * other.denominator >= other.numerator * denominator;
    //     }
    //
    //     // 输出流重载
    //     friend std::ostream &operator<<(std::ostream &os, Fraction &frac) {
    //         // os << frac.numerator;
    //         // if (frac.denominator != 1) {
    //         //     os << "/" << frac.denominator;
    //         // }
    //         // os << (frac.numerator * accuarcy) / frac.denominator;
    //         frac.normalize();
    //         os << frac.numerator << "/" << frac.denominator << " | " << (frac.numerator * accuarcy) / frac.denominator;
    //         return os;
    //     }
    //
    //     Fraction &operator=(const Fraction *fraction) {
    //         this->numerator = fraction->numerator;
    //         this->denominator = fraction->denominator;
    //         return *this;
    //     }
    //
    //     bool operator==(const int i) const {
    //         return denominator == 1 && numerator == i;
    //     }
    //
    //     operator Size() const {
    //         return numerator * accuarcy / denominator; // 示例转换逻辑
    //     }
    // };

    inline Size divide(const Size numerator, const Size denominator) {
        if (denominator == 0) return 0;
        return (numerator * accuarcy) / denominator;
    }

}


#endif //SUBGRAPHMATCHING_GLOBAL_H
