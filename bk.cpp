// file: bitset_bench.cpp
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <chrono>
#include <iostream>
#include <iomanip>   // for fixed, setprecision, setw
#include <random>    // for std::mt19937_64, std::uniform_int_distribution
#include <vector>
#include <chrono>    // for timing
#include <boost/dynamic_bitset.hpp>
using namespace std;
using clk = chrono::high_resolution_clock;

template<typename Fn>
double timeit(Fn&& fn, unsigned rounds = 1)
{
    auto st = clk::now();
    while (rounds--) fn();
    return chrono::duration<double, std::nano>(clk::now() - st).count();
}

struct ManualBits
{
    vector<uint64_t> w;
    explicit ManualBits(size_t nBits) : w((nBits + 63) >> 6) {}
    void randomFill() { for(auto& x: w) x = rng64(); }
    static uint64_t rng64()
    {
        static thread_local std::mt19937_64 gen{1234567};
        return gen();
    }
    // 四种集合运算：this = this (op) other
    void op_and(const ManualBits& o) { for(size_t i=0;i<w.size();++i) w[i] &= o.w[i]; }
    void op_or (const ManualBits& o) { for(size_t i=0;i<w.size();++i) w[i] |= o.w[i]; }
    void op_xor(const ManualBits& o) { for(size_t i=0;i<w.size();++i) w[i] ^= o.w[i]; }
    void op_sub(const ManualBits& o) { for(size_t i=0;i<w.size();++i) w[i] &= ~o.w[i]; }
};

int main(int argc,char* argv[])
{
    if (argc!=3){ cerr<<"Usage: "<<argv[0]<<" Nbits Rounds\n"; return 1; }
    const size_t N = stoull(argv[1]);
    const unsigned R = stoul(argv[2]);      // 每个操作重复 R 次取平均
    cout<<fixed<<setprecision(3);

    /* ---------- boost::dynamic_bitset ---------- */
    boost::dynamic_bitset<> b1(N), b2(N);
    b1.set();  b2.reset();  // 先随便变一下以避免 lazy-all-zero
    // 填随机
    for(size_t i=0;i<N;++i) if(rand()%2) b1.set(i);
    for(size_t i=0;i<N;++i) if(rand()%2) b2.set(i);

    auto t_and_boost = timeit([&]{
        auto tmp = b1;
        tmp &= b2;
    }, R)/R;
    auto t_or_boost  = timeit([&]{
        auto tmp = b1;
        tmp |= b2;
    }, R)/R;
    auto t_xor_boost = timeit([&]{
        auto tmp = b1;
        tmp ^= b2;
    }, R)/R;
    auto t_sub_boost = timeit([&]{
        auto tmp = b1;
        tmp &= ~b2;
    }, R)/R;

    /* ---------- vector<bool> ---------- */
    vector<bool> vb1(N), vb2(N);
    generate(vb1.begin(), vb1.end(), []{ return rand()&1; });
    generate(vb2.begin(), vb2.end(), []{ return rand()&1; });

    auto t_and_vbool = timeit([&]{
        vector<bool> tmp = vb1;
        for(size_t i=0;i<N;++i) tmp[i] = tmp[i] & vb2[i];
    }, R)/R;
    auto t_or_vbool  = timeit([&]{
        vector<bool> tmp = vb1;
        for(size_t i=0;i<N;++i) tmp[i] = tmp[i] | vb2[i];
    }, R)/R;
    auto t_xor_vbool = timeit([&]{
        vector<bool> tmp = vb1;
        for(size_t i=0;i<N;++i) tmp[i] = tmp[i] ^ vb2[i];
    }, R)/R;
    auto t_sub_vbool = timeit([&]{
        vector<bool> tmp = vb1;
        for(size_t i=0;i<N;++i) tmp[i] = tmp[i] & !vb2[i];
    }, R)/R;

    /* ---------- 手写 uint64_t 数组 ---------- */
    ManualBits m1(N), m2(N); m1.randomFill(); m2.randomFill();

    auto t_and_man = timeit([&]{
        auto tmp = m1; tmp.op_and(m2);
    }, R)/R;
    auto t_or_man  = timeit([&]{
        auto tmp = m1; tmp.op_or(m2);
    }, R)/R;
    auto t_xor_man = timeit([&]{
        auto tmp = m1; tmp.op_xor(m2);
    }, R)/R;
    auto t_sub_man = timeit([&]{
        auto tmp = m1; tmp.op_sub(m2);
    }, R)/R;

    /* ---------- 报告 ---------- */
    cout<< "N bits = "<<N<<", R = "<<R <<"\n";
    cout<< left << setw(18) << "Operation"
        << setw(14) << "boost"
        << setw(14) << "vector<bool>"
        << setw(14) << "manual64"
        << "\n";

    auto line=[&](string op,double a,double b,double c){
        cout<< setw(18)<<op<< setw(14)<<a<< setw(14)<<b<< setw(14)<<c<<"\n";
    };
    line("AND (∩)",  t_and_boost, t_and_vbool, t_and_man);
    line("OR  (∪)",  t_or_boost , t_or_vbool , t_or_man);
    line("XOR (△)",  t_xor_boost, t_xor_vbool, t_xor_man);
    line("SUB (\\)", t_sub_boost, t_sub_vbool, t_sub_man);
}