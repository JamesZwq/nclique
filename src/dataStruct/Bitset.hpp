//
// Created by 张文谦 on 25-8-4.
//

#ifndef BITSET_HPP
#define BITSET_HPP

#include <cstdint>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <limits>
#include <ostream>

#if defined(__x86_64__) || defined(_M_X64)
  #include <immintrin.h>
#endif
#if defined(__ARM_NEON)
  #include <arm_neon.h>
#endif

struct DynBitset {
    /* -------- 类型与常量 -------- */
    using word_t = std::uint64_t;
    static constexpr std::size_t W = 64;
    static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

    /* -------- 内存布局 -------- */
    std::size_t nbits = 0; ///< 真实位数
    std::vector<word_t> data; ///< ⌈nbits / 64⌉ 个 64-bit 块

    explicit DynBitset(std::size_t n = 0)
        : nbits(n), data((n + W - 1) / W, 0) {
    }

    /* -------- 容量 -------- */
    std::size_t size() const noexcept { return nbits; }

    /* -------- 单点操作 -------- */
    inline void set(std::size_t i) { data[i / W] |= word_t(1) << (i & 63); }
    inline void reset(std::size_t i) { data[i / W] &= ~(word_t(1) << (i & 63)); }
    inline bool test(std::size_t i) const { return (data[i / W] >> (i & 63)) & 1; }

    /* --- boost‑style no‑arg variants --- */
    inline void set() { set_all(); }
    inline void reset() { reset_all(); }

    /* -------- 整体操作 -------- */
    inline void set_all() {
        std::fill(data.begin(), data.end(), ~word_t(0));
        trim_tail();
    }

    inline void reset_all() { std::fill(data.begin(), data.end(), 0); }

    /* -------- 逻辑运算 (原地) -------- */
    inline void assign_and(const DynBitset &o) {
        // 支持不同长度时取交集；自赋值时仅需掩去尾部
        if (&o == this) { trim_tail(); return; }
        const std::size_t N = std::min(data.size(), o.data.size());

    #if defined(__AVX512F__)
        // 512-bit: 每次 8×u64
        std::size_t i = 0, L = N & ~std::size_t(7);
        for (; i < L; i += 8) {
            __m512i a = _mm512_loadu_si512(reinterpret_cast<const void*>(&data[i]));
            __m512i b = _mm512_loadu_si512(reinterpret_cast<const void*>(&o.data[i]));
            a = _mm512_and_si512(a, b);
            _mm512_storeu_si512(reinterpret_cast<void*>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] &= o.data[i];
    #elif defined(__AVX2__)
        // 256-bit: 每次 4×u64
        std::size_t i = 0, L = N & ~std::size_t(3);
        for (; i < L; i += 4) {
            __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&data[i]));
            __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&o.data[i]));
            a = _mm256_and_si256(a, b);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] &= o.data[i];
    #elif defined(__SSE2__)
        // 128-bit: 每次 2×u64
        std::size_t i = 0, L = N & ~std::size_t(1);
        for (; i < L; i += 2) {
            __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&data[i]));
            __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&o.data[i]));
            a = _mm_and_si128(a, b);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(&data[i]), a);
        }
        if (i < N) data[i] &= o.data[i];
    #elif defined(__ARM_NEON)
        // NEON: 每次 2×u64
        std::size_t i = 0, L = N & ~std::size_t(1);
        for (; i < L; i += 2) {
            uint64x2_t a = vld1q_u64(&data[i]);
            uint64x2_t b = vld1q_u64(&o.data[i]);
            a = vandq_u64(a, b);
            vst1q_u64(&data[i], a);
        }
        if (i < N) data[i] &= o.data[i];
    #else
        for (std::size_t i = 0; i < N; ++i) data[i] &= o.data[i];
    #endif
        // AND 理论上不会点亮尾部，但为了严谨性，仍然掩掉无效尾位
        trim_tail();
    }

    inline void assign_or(const DynBitset &o) {
        if (&o == this) { trim_tail(); return; }
        const std::size_t N = std::min(data.size(), o.data.size());

    #if defined(__AVX512F__)
        std::size_t i = 0, L = N & ~std::size_t(7);
        for (; i < L; i += 8) {
            __m512i a = _mm512_loadu_si512(reinterpret_cast<const void*>(&data[i]));
            __m512i b = _mm512_loadu_si512(reinterpret_cast<const void*>(&o.data[i]));
            a = _mm512_or_si512(a, b);
            _mm512_storeu_si512(reinterpret_cast<void*>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] |= o.data[i];
    #elif defined(__AVX2__)
        std::size_t i = 0, L = N & ~std::size_t(3);
        for (; i < L; i += 4) {
            __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&data[i]));
            __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&o.data[i]));
            a = _mm256_or_si256(a, b);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] |= o.data[i];
    #elif defined(__SSE2__)
        std::size_t i = 0, L = N & ~std::size_t(1);
        for (; i < L; i += 2) {
            __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&data[i]));
            __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&o.data[i]));
            a = _mm_or_si128(a, b);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(&data[i]), a);
        }
        if (i < N) data[i] |= o.data[i];
    #elif defined(__ARM_NEON)
        std::size_t i = 0, L = N & ~std::size_t(1);
        for (; i < L; i += 2) {
            uint64x2_t a = vld1q_u64(&data[i]);
            uint64x2_t b = vld1q_u64(&o.data[i]);
            a = vorrq_u64(a, b);
            vst1q_u64(&data[i], a);
        }
        if (i < N) data[i] |= o.data[i];
    #else
        for (std::size_t i = 0; i < N; ++i) data[i] |= o.data[i];
    #endif
        // OR 可能点亮无效尾位，必须掩掉
        trim_tail();
    }

    inline void assign_xor(const DynBitset &o) {
        if (&o == this) { reset_all(); return; }
        const std::size_t N = std::min(data.size(), o.data.size());

    #if defined(__AVX512F__)
        std::size_t i = 0, L = N & ~std::size_t(7);
        for (; i < L; i += 8) {
            __m512i a = _mm512_loadu_si512(reinterpret_cast<const void*>(&data[i]));
            __m512i b = _mm512_loadu_si512(reinterpret_cast<const void*>(&o.data[i]));
            a = _mm512_xor_si512(a, b);
            _mm512_storeu_si512(reinterpret_cast<void*>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] ^= o.data[i];
    #elif defined(__AVX2__)
        std::size_t i = 0, L = N & ~std::size_t(3);
        for (; i < L; i += 4) {
            __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&data[i]));
            __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&o.data[i]));
            a = _mm256_xor_si256(a, b);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] ^= o.data[i];
    #elif defined(__SSE2__)
        std::size_t i = 0, L = N & ~std::size_t(1);
        for (; i < L; i += 2) {
            __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&data[i]));
            __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&o.data[i]));
            a = _mm_xor_si128(a, b);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(&data[i]), a);
        }
        if (i < N) data[i] ^= o.data[i];
    #elif defined(__ARM_NEON)
        std::size_t i = 0, L = N & ~std::size_t(1);
        for (; i < L; i += 2) {
            uint64x2_t a = vld1q_u64(&data[i]);
            uint64x2_t b = vld1q_u64(&o.data[i]);
            a = veorq_u64(a, b);
            vst1q_u64(&data[i], a);
        }
        if (i < N) data[i] ^= o.data[i];
    #else
        for (std::size_t i = 0; i < N; ++i) data[i] ^= o.data[i];
    #endif
        trim_tail();
    }

    // 统计 (*this & o) 中 1 的个数；不产生临时对象
    inline std::size_t count_and(const DynBitset& o) const noexcept {
        const std::size_t N = std::min(data.size(), o.data.size());
        std::size_t c = 0;
        for (std::size_t i = 0; i < N; ++i) {
            c += __builtin_popcountll(data[i] & o.data[i]);
        }
        return c;
    }

    /* -------- 计数 / 判空 -------- */
    inline std::size_t count() const {
        // 同 boost 的 count()
        std::size_t c = 0;
        for (auto w: data) c += __builtin_popcountll(w);
        return c;
    }

    inline std::size_t popcount() const { return count(); } // 兼容老名字
    inline bool any() const {
        for (auto w: data) if (w) return true;
        return false;
    }

    inline bool none() const { return !any(); }

    /* -------- 查找 -------- */
    inline std::size_t find_first() const {
        for (std::size_t blk = 0, N = data.size(); blk < N; ++blk) {
            word_t w = data[blk];
            if (w) return blk * W + __builtin_ctzll(w);
        }
        return npos;
    }

    inline std::size_t find_next(std::size_t prev) const {
        ++prev;
        if (prev >= nbits) return npos;
        std::size_t blk = prev / W;
        word_t w = data[blk] & (~word_t(0) << (prev & 63));
        if (w) return blk * W + __builtin_ctzll(w);
        for (++blk; blk < data.size(); ++blk) {
            w = data[blk];
            if (w) return blk * W + __builtin_ctzll(w);
        }
        return npos;
    }

    /* -------- 遍历 -------- */
    template<class Fn>
    inline void for_each_bit(Fn &&fn) const {
        for (std::size_t blk = 0, N = data.size(); blk < N; ++blk) {
            word_t w = data[blk];
            while (w) {
                std::size_t b = __builtin_ctzll(w);
                fn(blk * W + b);
                w &= w - 1;
            }
        }
    }

private:
    /* -------- 抹掉末尾多余位 -------- */
    inline void trim_tail() {
        const unsigned extra = nbits & 63; // 0~63
        if (extra) {
            data.back() &= ((word_t(1) << extra) - 1);
        }
    }

    /* =====================  bitwise operators  ===================== */
public:
    friend inline DynBitset operator~(DynBitset a) {
        // unary NOT
        for (auto &w: a.data) w = ~w;
        a.trim_tail();
        return a;
    }

    friend inline DynBitset operator&(DynBitset a, const DynBitset &b) {
        a.assign_and(b);
        return a;
    }

    friend inline DynBitset operator|(DynBitset a, const DynBitset &b) {
        a.assign_or(b);
        return a;
    }

    friend inline DynBitset operator^(DynBitset a, const DynBitset &b) {
        a.assign_xor(b);
        return a;
    }

    inline DynBitset &operator&=(const DynBitset &o) {
        assign_and(o);
        return *this;
    }

    inline DynBitset &operator|=(const DynBitset &o) {
        assign_or(o);
        return *this;
    }

    inline DynBitset &operator^=(const DynBitset &o) {
        assign_xor(o);
        return *this;
    }

    // <<, output to ostream
    friend std::ostream &operator<<(std::ostream &os, const DynBitset &bs) {
        os << "DynBitset(" << bs.nbits << "): ";
        for (std::size_t i = 0; i < bs.nbits; ++i) {
            if (bs.test(i)) os << '1';
            else os << '0';
        }
        return os;
    }
    // >>, input from istream
};

inline void subtractVertexFromCounts(
    size_t v,
    std::vector<uint16_t> &cnt,
    const std::vector<uint16_t> &lim,
    DynBitset &violMask,
    const std::vector<DynBitset> &vtxMask) {
    vtxMask[v].for_each_bit([&](size_t c) {
        if (cnt[c] && --cnt[c] < lim[c]) violMask.reset(c);
    });
}

inline void addVertexToCounts(
    size_t v,
    std::vector<uint16_t> &cnt,
    const std::vector<uint16_t> &lim,
    DynBitset &violMask,
    const std::vector<DynBitset> &vtxMask) {
    vtxMask[v].for_each_bit([&](size_t c) {
        if (++cnt[c] >= lim[c]) violMask.set(c);
    });
}
#endif //BITSET_HPP
