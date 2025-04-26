// bench.cpp  ——  push_back + pop_back 基准
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

using clk  = std::chrono::steady_clock;
using Size = uint32_t;                   // == daf::Size，够用到 4 G
inline void compiler_barrier() { asm volatile("" ::: "memory"); }
/* ---------- 1. StaticVector：固定容量、仅移动 ---------- */
template<typename T>
class StaticVector {
    static_assert(std::is_trivially_copyable_v<T>);
public:
    explicit StaticVector(Size cap)
        : cap_(cap), sz_(0), buf_(new T[cap]) {}
    ~StaticVector() { delete[] buf_; }

    StaticVector(const StaticVector&)            = delete;
    StaticVector& operator=(const StaticVector&) = delete;

    StaticVector(StaticVector&& o) noexcept
        : cap_(o.cap_), sz_(o.sz_), buf_(o.buf_) { o.buf_ = nullptr; }
    StaticVector& operator=(StaticVector&& o) noexcept {
        if (this != &o) {
            delete[] buf_;
            cap_ = o.cap_;  sz_ = o.sz_;  buf_ = o.buf_;
            o.buf_ = nullptr;
        }
        return *this;
    }

    void push_back(const T& v) noexcept { buf_[sz_++] = v; compiler_barrier();}   // 调用方保证不溢出
    void pop_back()        noexcept { --sz_; compiler_barrier();}
    T&       operator[](Size i)       noexcept { return buf_[i]; }
    const T& operator[](Size i) const noexcept { return buf_[i]; }
    Size size() const noexcept { return sz_; }
    T* begin() noexcept { return buf_; }
    T* end()   noexcept { return buf_ + sz_; }

private:
    Size  cap_, sz_;
    T*    buf_;
};

/* ---------- 2. FixedVector：外部分配 + 零分支 push/pop ---------- */
template<typename T>
class FixedVector {
public:
    FixedVector(T* ptr, Size cap) : ptr_(ptr), cap_(cap), sz_(0) {}
    void push_back(const T& v) noexcept { ptr_[sz_++] = v; compiler_barrier();}
    void pop_back()            noexcept { --sz_; compiler_barrier();}
    T&       operator[](Size i)       noexcept { return ptr_[i]; }
    const T& operator[](Size i) const noexcept { return ptr_[i]; }
    Size size() const noexcept { return sz_; }
    T* begin() noexcept { return ptr_; }
    T* end()   noexcept { return ptr_ + sz_; }

private:
    T*   ptr_;
    Size cap_, sz_;
};

/* ---------- 3. 通用基准模板：交替 push 32 次 / pop 32 次 ---------- */
template<class Vec>
double bench_push_pop(Size ops, const std::string& tag)
{
    // 1. 至少能放下 ops 个元素
    Vec v(ops);

    /* ---------- 构建阶段 ---------- */
    auto t0 = clk::now();

    // push ops 元素
    for (Size i = 0; i < ops; ++i) v.push_back(i);

    // pop 一半，让容器里还剩 ops/2 个元素
    for (Size i = 0; i < ops / 2; ++i) v.pop_back();

    auto t1 = clk::now();
    /* ---------- 计时到此 ---------- */

    /* ---------- 强制读取数据，阻止 DCE ---------- */
    Size sum = 0;
//    for (auto x : v) sum += x;
    for (Size i = 0; i < v.size(); ++i) sum += v[i];
    volatile Size sink = sum; (void)sink;
    /* ------------------------------------------- */

    double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cout << std::setw(9) << tag << " | push "
              << ops << " / pop " << ops/2
              << " : " << std::fixed << std::setprecision(6)
              << secs << " s\n";
    return secs;
}

/* ---------- 4. 三种容器包装 ---------- */
// 4.1 StdVec：std::vector + reserve
struct StdVec {
    explicit StdVec(Size n) { v_.reserve(n); }
    void push_back(Size x) noexcept { v_.push_back(x); }
    void pop_back()         noexcept { v_.pop_back();  }
    Size& operator[](Size i) noexcept { return v_[i]; }
    const Size& operator[](Size i) const noexcept { return v_[i]; }
    Size size() const noexcept { return v_.size(); }
private:
    std::vector<Size> v_;
};

// 4.2 StaticVector 包装
struct StaticHolder {
    explicit StaticHolder(Size n) : sv_(n) {}
    void push_back(Size x) noexcept { sv_.push_back(x); }
    void pop_back()         noexcept { sv_.pop_back();  }
    Size& operator[](Size i) noexcept { return sv_[i]; }
    const Size& operator[](Size i) const noexcept { return sv_[i]; }
    Size size() const noexcept { return sv_.size(); }
private:
    StaticVector<Size> sv_;
};

// 4.3 FixedVector 包装
struct FixedHolder {
    explicit FixedHolder(Size n)
        : buf_(std::make_unique<Size[]>(n)), fv_(buf_.get(), n) {}
    void push_back(Size x) noexcept { fv_.push_back(x); }
    void pop_back()         noexcept { fv_.pop_back();  }
    Size& operator[](Size i) noexcept { return fv_[i]; }
    const Size& operator[](Size i) const noexcept { return fv_[i]; }
    Size size() const noexcept { return fv_.size(); }
private:
    std::unique_ptr<Size[]> buf_;
    FixedVector<Size>       fv_;
};

/* ---------- 5. main ---------- */
int main(int argc, char** argv)
{
    Size ops = (argc > 1 ? std::strtoull(argv[1], nullptr, 10) : 50'000'000);
    std::cout << "Total push+pop operations: " << ops << "\n";

    bench_push_pop<FixedHolder >(ops, "FixedVec");
    bench_push_pop<StaticHolder>(ops, "StaticVec");
    bench_push_pop<StdVec      >(ops, "StdVec");

    return 0;
}