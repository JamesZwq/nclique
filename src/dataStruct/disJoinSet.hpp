#pragma once
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cassert>

class BipartiteDSU {
public:

    /**
     *  @param aSize        A （）
     *  @param bInit         B 
     *  @param reserveHint   B （ 0）
     */
    BipartiteDSU(std::size_t aSize,
                 std::size_t bInit,
                 std::size_t reserveHint = 0)
        : aSize_(aSize), bCnt_(bInit)
    {
        parent_.reserve(aSize_ + bInit + reserveHint + 1);     // +1  0
        parent_.assign(aSize_ + bInit + 1, -1);                // 
    }

    /*---------  ： B   ---------*/
    //  newCnt ≤  bCnt_ ； [bCnt_+1, newCnt]  parent_  -1
    inline void expandB(std::size_t newCnt) {
        if (newCnt <= bCnt_) return;
        parent_.resize(aSize_ + newCnt + 1, -1);
        bCnt_ = newCnt;
    }

    /*---------  /（） ---------*/
    inline void uniteAA(std::size_t a1, std::size_t a2) { unite(idxA(a1), idxA(a2)); }
    inline void uniteAB(std::size_t a , std::size_t b ) { unite(idxA(a ), idxB(b )); }
    inline void uniteBB(std::size_t b1, std::size_t b2) { unite(idxB(b1), idxB(b2)); }

    inline bool sameAA(std::size_t a1, std::size_t a2) { return find(idxA(a1)) == find(idxA(a2)); }
    inline bool sameAB(std::size_t a , std::size_t b ) { return find(idxA(a )) == find(idxB(b )); }
    inline bool sameBB(std::size_t b1, std::size_t b2) { return find(idxB(b1)) == find(idxB(b2)); }

    inline std::size_t setSizeA(std::size_t a) { return setSizeRaw(idxA(a)); }
    inline std::size_t setSizeB(std::size_t b) { return setSizeRaw(idxB(b)); }

    friend std::ostream& operator<<(std::ostream &os, const BipartiteDSU &dsu) {
        os << "BipartiteDSU: aSize=" << dsu.aSize_ << ", bCnt=" << dsu.bCnt_ << "\n";
        os << "Parent: ";
        for (std::size_t i = 0; i < dsu.parent_.size(); ++i) {
            os << dsu.parent_[i] << " ";
        }
        return os;
    }

    void print() const {
        std::cout << *this << std::endl;
    }
private:
    /*---------   ---------*/
    std::size_t           aSize_;
    std::size_t           bCnt_;
    std::vector<int32_t>  parent_;          // parent_[x] < 0 ⇒ ，|value|=size

    // ：A  1…aSize，B  aSize+1…aSize+bCnt_
    // 0-based 
    inline std::size_t idxA(std::size_t a) {              // A: 0 … aSize_-1
        assert(a < aSize_);
        return a;
    }
    inline std::size_t idxB(std::size_t b) const {               // B: aSize_ … aSize_+bCnt_-1
        assert(b < bCnt_);
        return aSize_ + b;
    }

    inline std::size_t find(std::size_t x) {
        return parent_[x] < 0 ? x : parent_[x] = static_cast<int32_t>(find(parent_[x]));
    }
    inline void unite(std::size_t x, std::size_t y) {
        x = find(x); y = find(y);
        if (x == y) return;
        if (parent_[x] > parent_[y]) std::swap(x, y);   // x （-size ）
        parent_[x] += parent_[y];
        parent_[y]  = static_cast<int32_t>(x);
    }
    inline std::size_t setSizeRaw(std::size_t idx) {
        return static_cast<std::size_t>(-parent_[find(idx)]);
    }
};