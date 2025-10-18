//
// Created by _ on 25-8-4.
//

#ifndef BITSET_HPP
#define BITSET_HPP

#include <cstdint>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <limits>
#include <ostream>
#include <array>
#include <cassert>

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif
#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif

struct DynBitset {
    /* --------  -------- */
    using word_t = std::uint64_t;
    static constexpr std::size_t W = 64;
    static constexpr std::size_t MAX_BITS = 400;                 // 
    static constexpr std::size_t NWORDS  = (MAX_BITS + W - 1) / W;
    static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

    /* --------  --------
     * “ + ”：
     * - data  400 （7×u64），/
     * - nbits （<= MAX_BITS）
     */
    std::size_t nbits = MAX_BITS;                 ///< （）
    alignas(64) std::array<word_t, NWORDS> data{};///<  400 ，

    // ： MAX_BITS（400），
    explicit DynBitset() noexcept = default;

    // ；
    inline void setSize(std::size_t n) noexcept {
        if (this->nbits == n) return;
        const std::size_t oldWN = word_count();
        nbits = n;
        const std::size_t newWN = word_count();
        // ， word 
        for (std::size_t i = oldWN; i < newWN; ++i) data[i] = 0;
        trim_tail();
    }

    //  u64 
    inline std::size_t word_count() const noexcept { return (nbits + W - 1) / W; }

    // （ nbits ）
    static inline word_t tail_mask_from_bits(std::size_t n) noexcept {
        const unsigned e = static_cast<unsigned>(n & 63);
        return e ? ((word_t(1) << e) - 1) : ~word_t(0);
    }
    inline word_t tail_mask() const noexcept { return tail_mask_from_bits(nbits); }

    /* --------  -------- */
    std::size_t size() const noexcept { return nbits; }

    /* --------  -------- */
    inline void set(std::size_t i) {
        assert(i < nbits && "DynBitset::set out of range");
        data[i / W] |= word_t(1) << (i & 63);
    }
    inline void reset(std::size_t i) {
        assert(i < nbits && "DynBitset::reset out of range");
        data[i / W] &= ~(word_t(1) << (i & 63));
    }
    inline bool test(std::size_t i) const {
        assert(i < nbits && "DynBitset::test out of range");
        return (data[i / W] >> (i & 63)) & 1;
    }

    /* --- boost‑style no‑arg variants --- */
    inline void set() { set_all(); }
    inline void reset() { reset_all(); }

    /* --------  -------- */
    inline void set_all() {
        const std::size_t WN = word_count();
        if (WN == 0) return;
        //  WN-1  word  1
        for (std::size_t i = 0; i + 1 < WN; ++i) data[i] = ~word_t(0);
        //  word  nbits 
        data[WN - 1] = tail_mask();
        //  word  0（，）
        for (std::size_t i = WN; i < data.size(); ++i) data[i] = 0;
    }

    inline void reset_all() {
        const std::size_t WN = word_count();
        for (std::size_t i = 0; i < WN; ++i) data[i] = 0;
        // ：，
        for (std::size_t i = WN; i < data.size(); ++i) data[i] = 0;
    }

    /* --------  () -------- */
    inline void assign_and(const DynBitset &o) {
        if (&o == this) { trim_tail(); return; }
        const std::size_t WN  = word_count();
        const std::size_t WNo = o.word_count();
        const std::size_t N   = std::min(WN, WNo);

#if defined(__AVX512F__)
        std::size_t i = 0, L = N & ~std::size_t(7);
        for (; i < L; i += 8) {
            __m512i a = _mm512_loadu_si512(reinterpret_cast<const void *>(&data[i]));
            __m512i b = _mm512_loadu_si512(reinterpret_cast<const void *>(&o.data[i]));
            a = _mm512_and_si512(a, b);
            _mm512_storeu_si512(reinterpret_cast<void *>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] &= o.data[i];
#elif defined(__AVX2__)
        std::size_t i = 0, L = N & ~std::size_t(3);
        for (; i < L; i += 4) {
            __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&data[i]));
            __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&o.data[i]));
            a = _mm256_and_si256(a, b);
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] &= o.data[i];
#elif defined(__SSE2__)
        std::size_t i = 0, L = N & ~std::size_t(1);
        for (; i < L; i += 2) {
            __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i *>(&data[i]));
            __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i *>(&o.data[i]));
            a = _mm_and_si128(a, b);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(&data[i]), a);
        }
        if (i < N) data[i] &= o.data[i];
#elif defined(__ARM_NEON)
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
        // ，AND  0 ，
        for (std::size_t i = N; i < WN; ++i) data[i] = 0;
        trim_tail();
    }

    inline void assign_or(const DynBitset &o) {
        if (&o == this) { trim_tail(); return; }
        const std::size_t WN  = word_count();
        const std::size_t WNo = o.word_count();
        const std::size_t N   = std::min(WN, WNo);

#if defined(__AVX512F__)
        std::size_t i = 0, L = N & ~std::size_t(7);
        for (; i < L; i += 8) {
            __m512i a = _mm512_loadu_si512(reinterpret_cast<const void *>(&data[i]));
            __m512i b = _mm512_loadu_si512(reinterpret_cast<const void *>(&o.data[i]));
            a = _mm512_or_si512(a, b);
            _mm512_storeu_si512(reinterpret_cast<void *>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] |= o.data[i];
#elif defined(__AVX2__)
        std::size_t i = 0, L = N & ~std::size_t(3);
        for (; i < L; i += 4) {
            __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&data[i]));
            __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&o.data[i]));
            a = _mm256_or_si256(a, b);
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] |= o.data[i];
#elif defined(__SSE2__)
        std::size_t i = 0, L = N & ~std::size_t(1);
        for (; i < L; i += 2) {
            __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i *>(&data[i]));
            __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i *>(&o.data[i]));
            a = _mm_or_si128(a, b);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(&data[i]), a);
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
        trim_tail();
    }

    inline void assign_xor(const DynBitset &o) {
        if (&o == this) { reset_all(); return; }
        const std::size_t WN  = word_count();
        const std::size_t WNo = o.word_count();
        const std::size_t N   = std::min(WN, WNo);

#if defined(__AVX512F__)
        std::size_t i = 0, L = N & ~std::size_t(7);
        for (; i < L; i += 8) {
            __m512i a = _mm512_loadu_si512(reinterpret_cast<const void *>(&data[i]));
            __m512i b = _mm512_loadu_si512(reinterpret_cast<const void *>(&o.data[i]));
            a = _mm512_xor_si512(a, b);
            _mm512_storeu_si512(reinterpret_cast<void *>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] ^= o.data[i];
#elif defined(__AVX2__)
        std::size_t i = 0, L = N & ~std::size_t(3);
        for (; i < L; i += 4) {
            __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&data[i]));
            __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&o.data[i]));
            a = _mm256_xor_si256(a, b);
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(&data[i]), a);
        }
        for (; i < N; ++i) data[i] ^= o.data[i];
#elif defined(__SSE2__)
        std::size_t i = 0, L = N & ~std::size_t(1);
        for (; i < L; i += 2) {
            __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i *>(&data[i]));
            __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i *>(&o.data[i]));
            a = _mm_xor_si128(a, b);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(&data[i]), a);
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

    //  (*this & o)  1 ；
    inline std::size_t count_and(const DynBitset &o) const noexcept {
        const std::size_t N = std::min(word_count(), o.word_count());
        std::size_t c = 0;
        for (std::size_t i = 0; i < N; ++i) {
            c += __builtin_popcountll(data[i] & o.data[i]);
        }
        return c;
    }

    /* --------  /  -------- */
    inline std::size_t count() const {
        //  boost  count()
        std::size_t c = 0;
        const std::size_t WN = word_count();
        for (std::size_t i = 0; i < WN; ++i) c += __builtin_popcountll(data[i]);
        return c;
    }

    inline std::size_t popcount() const { return count(); } // 
    inline bool any() const {
        const std::size_t WN = word_count();
        for (std::size_t i = 0; i < WN; ++i) if (data[i]) return true;
        return false;
    }

    inline bool none() const { return !any(); }

    /* --------  -------- */
    inline std::size_t find_first() const {
        const std::size_t WN = word_count();
        if (WN == 0) return npos;
        for (std::size_t blk = 0; blk < WN; ++blk) {
            word_t w = data[blk];
            if (blk + 1 == WN) w &= tail_mask();
            if (w) return blk * W + __builtin_ctzll(w);
        }
        return npos;
    }

    inline std::size_t find_next(std::size_t prev) const {
        ++prev;
        if (prev >= nbits) return npos;
        const std::size_t WN = word_count();
        std::size_t blk = prev / W;
        word_t w = data[blk] & (~word_t(0) << (prev & 63));
        if (blk + 1 == WN) w &= tail_mask();
        if (w) return blk * W + __builtin_ctzll(w);
        for (++blk; blk < WN; ++blk) {
            w = data[blk];
            if (blk + 1 == WN) w &= tail_mask();
            if (w) return blk * W + __builtin_ctzll(w);
        }
        return npos;
    }

    /* --------  -------- */
    template<class Fn>
    inline void for_each_bit(Fn &&fn) const {
        const std::size_t WN = word_count();
        for (std::size_t blk = 0; blk < WN; ++blk) {
            word_t w = data[blk];
            if (blk + 1 == WN) w &= tail_mask();
            while (w) {
                std::size_t b = __builtin_ctzll(w);
                fn(blk * W + b);
                w &= w - 1;
            }
        }
    }

private:
    /* --------  -------- */
    inline void trim_tail() {
        const std::size_t WN = word_count();
        if (WN == 0) return;
        const unsigned extra = static_cast<unsigned>(nbits & 63); // 0~63
        if (extra) {
            data[WN - 1] &= ((word_t(1) << extra) - 1);
        }
    }

    /* =====================  bitwise operators  ===================== */
public:
    friend inline DynBitset operator~(DynBitset a) {
        const std::size_t WN = a.word_count();
        for (std::size_t i = 0; i < WN; ++i) a.data[i] = ~a.data[i];
        //  0（）
        for (std::size_t i = WN; i < a.data.size(); ++i) a.data[i] = 0;
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