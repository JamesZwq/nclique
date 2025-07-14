/*********************************************************************
 * split_bench.cpp  ──  Baseline-MU  VS  Forward-Checking-Prop
 *   • 流式打印每个极大合法子集       • O(1) 记内存：只保留计数+64bit 哈希
 *   • 校验：两实现 (cnt, hash) 完全一致 ⇒ 结果相同
 *
 *  编译:  clang++/g++  -std=c++20 -O3 -march=native split_bench.cpp -o bench
 *  运行:  ./bench  n  m  forbiddenLen  keepProb  k
 *          (若无参数默认 n=250 m=60 len=6 keep=0.1 k=80)
 *
 *  2025-06-21  Wenqian Zhang
 *********************************************************************/
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>
using namespace std;
using U64 = unsigned long long;
constexpr int MAXN = 400;

/*------------------------------------------------------------
 *  1. 简易 400-bit 定长 bitset
 *-----------------------------------------------------------*/
struct B400 {
    array<U64,(MAXN+63)/64> w{};               // 零初始化
    inline void set (int i)       { w[i>>6] |= 1ULL << (i & 63); }
    inline bool test(int i) const { return (w[i>>6] >> (i & 63)) & 1ULL; }
};

/*------------------------------------------------------------
 *  2. Streamer：流式输出 + 计数 + 64bit FNV-1a 哈希
 *-----------------------------------------------------------*/
struct Streamer {
    size_t   cnt   = 0;
    uint64_t hash  = 0;

    template<class VecInt>
    static uint64_t fastHash(const VecInt& v){       // FNV-1a 64
        uint64_t h = 146527ull + v.size();
        for(int x: v)
            h = (h ^ (uint64_t(x) + 0x9e3779b97f4a7c15ULL))
                * 11400714819323198485ULL;
        return h;
    }

    void emit(const std::vector<int>& idxSet, const std::vector<int>& mapNew2Old){
        /* 1) 打印（原元素索引顺序） */
        for (size_t i = 0; i < idxSet.size(); ++i) {
            int orig = mapNew2Old[idxSet[i]];
            std::printf("%d%c", orig, (i + 1 == idxSet.size() ? '\n' : ' '));
        }
        /* 2) 滚动哈希 + 计数 */
        hash ^= fastHash(idxSet);
        ++cnt;
    }
};

/*------------------------------------------------------------
 *  3. Baseline：Shrunk + Murakami-Uno (最小截集)  改正为补集输出
 *-----------------------------------------------------------*/
struct MU {
    int n, kNeed;
    vector<B400> edge;
    vector<vector<int>> adj;
    vector<int> C, U;                      // 当前 hitting-set C  与  未覆盖边 U
    Streamer* S; const vector<int>* map;   // map: new→old

    MU(int n_, int kNeed_, const vector<vector<int>>& forb,
       const vector<int>& map_, Streamer& st)
        : n(n_), kNeed(kNeed_), edge(forb.size()), adj(n_), S(&st), map(&map_)
    {
        for (int e = 0; e < (int)forb.size(); ++e) {
            for (int v : forb[e]) { edge[e].set(v); adj[v].push_back(e); }
        }
    }
    /* ---- 返回 A'\C (= 极大合法子集在缩减空间中的索引表示) ---- */
    vector<int> complement() const {
        vector<int> T; T.reserve(n - C.size());
        for (int v = 0; v < n; ++v)
            if (find(C.begin(), C.end(), v) == C.end())
                T.push_back(v);
        return T;
    }
    /* ---- 深搜 ---- */
    void dfs() {
        if ((int)C.size() > n - kNeed) return;           // k 裁剪
        if (U.empty()) {                                // C 已击中全部禁用边
            S->emit(complement(), *map);                // **补集输出**
            return;
        }
        /* pivot = freq 最大的元素 */
        vector<int> freq(n, 0);
        for (int e : U)
            for (int v = 0; v < n; ++v)
                if (edge[e].test(v)) ++freq[v];
        int pv = int(max_element(freq.begin(), freq.end()) - freq.begin());

        /* --- branch 1: 选 pv 进入 C --- */
        C.push_back(pv);
        vector<int> U2; U2.reserve(U.size());
        for (int e : U) if (!edge[e].test(pv)) U2.push_back(e);
        auto U_bak = std::move(U); U = std::move(U2); dfs(); U = std::move(U_bak);
        C.pop_back();

        /* --- branch 2: ban pv (必须在含 pv 的第一条边选别人) --- */
        if (adj[pv].empty()) return;
        int e0 = adj[pv][0];
        for (int v = 0; v < n; ++v) if (edge[e0].test(v) && v != pv) {
            if (find(C.begin(), C.end(), v) != C.end()) continue;
            C.push_back(v);
            vector<int> U3; U3.reserve(U.size());
            for (int e : U) if (!edge[e].test(v)) U3.push_back(e);
            auto U_bak2 = std::move(U); U = std::move(U3); dfs(); U = std::move(U_bak2);
            C.pop_back();
        }
    }
    void run() { C.clear(); U.resize(edge.size()); iota(U.begin(), U.end(), 0); dfs(); }
};

/*------------------------------------------------------------
 *  4. Forward-Checking + 单洞 Propagation 版本
 *-----------------------------------------------------------*/
struct FC {
    int n, kNeed;
    vector<vector<int>> E, adj;
    vector<int> rem;                 // edge 剩余可选元素数
    vector<uint8_t> stt;             // 0 free | 1 sel | 2 ban
    vector<int> sel, queue;
    Streamer* S; const vector<int>* map;

    FC(int n_, int kNeed_, const vector<vector<int>>& forb,
       const vector<int>& map_, Streamer& st)
        : n(n_), kNeed(kNeed_), E(forb), adj(n_), rem(forb.size()),
          stt(n_, 0), S(&st), map(&map_)
    {
        for (int e = 0; e < (int)E.size(); ++e) {
            rem[e] = E[e].size();
            for (int v : E[e]) adj[v].push_back(e);
        }
    }
    bool propagate(vector<int>& changedE, vector<int>& changedV) {
        for (int e : changedE) if (rem[e] == 1) queue.push_back(e);
        while (!queue.empty()) {
            int e = queue.back(); queue.pop_back();
            int lone = -1;
            for (int v : E[e]) if (stt[v] == 0) { lone = v; break; }
            if (lone == -1) return false;        // 矛盾
            stt[lone] = 2; changedV.push_back(lone);   // ban
            for (int ed : adj[lone]) {
                if (--rem[ed] == 0) return false;
                changedE.push_back(ed);
                if (rem[ed] == 1) queue.push_back(ed);
            }
        }
        return true;
    }
    void dfs(int freeCnt) {
        if ((int)sel.size() + freeCnt < kNeed) return;
        if (freeCnt == 0) { S->emit(sel, *map); return; }

        /* pivot = free & max degree */
        int pv = -1, best = -1;
        for (int v = 0; v < n; ++v) if (stt[v] == 0) {
            int d = 0; for (int e : adj[v]) if (rem[e]) ++d;
            if (d > best) { best = d; pv = v; }
        }

        /* ------ branch A: select pv ------ */
        {
            vector<int> chE, chV;
            stt[pv] = 1; sel.push_back(pv);       chV.push_back(pv);
            for (int e : adj[pv]) if (rem[e]) { rem[e] = 0; chE.push_back(e); }
            if (propagate(chE, chV))
                dfs(freeCnt - int(chV.size()) + 1);
            /* undo */
            for (int v : chV) stt[v] = 0;
            for (int e = 0; e < (int)E.size(); ++e) rem[e] = E[e].size();
            sel.pop_back();
        }
        /* ------ branch B: ban pv ------ */
        {
            vector<int> chE, chV;
            stt[pv] = 2; chV.push_back(pv);
            for (int e : adj[pv]) { --rem[e]; chE.push_back(e); }
            if (propagate(chE, chV))
                dfs(freeCnt - int(chV.size()) + 1);
            /* undo */
            for (int v : chV) stt[v] = 0;
            for (int e = 0; e < (int)E.size(); ++e) rem[e] = E[e].size();
        }
    }
    void run() { dfs(n); }
};

/*------------------------------------------------------------
 *  5. 随机用例生成
 *-----------------------------------------------------------*/
struct Case {
    vector<int> A;
    vector<vector<int>> S;
    vector<char> K;   // 0/1
    int k;
};
Case genRnd(int n, int m, int len, double keepProb, int k) {
    std::mt19937_64 rng(12345);
    Case c;
    c.A.resize(n); std::iota(c.A.begin(), c.A.end(), 0);
    c.K.assign(n, 0);
    for (int i = 0; i < n; ++i)
        if (std::uniform_real_distribution<>(0,1)(rng) < keepProb) c.K[i] = 1;

    std::uniform_int_distribution<int> pos(0, n-1);
    for (int i = 0; i < m; ++i) {
        std::unordered_set<int> st;
        while ((int)st.size() < len) st.insert(pos(rng));
        c.S.emplace_back(st.begin(), st.end());
    }
    c.k = k;
    return c;
}

/*------------------------------------------------------------
 *  6. 计时辅助
 *-----------------------------------------------------------*/
template<class F> long long ms(F&& fn){
    auto t0 = std::chrono::high_resolution_clock::now();
    fn();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
}

/*------------------------------------------------------------
 *  7. main
 *-----------------------------------------------------------*/
int main(int argc, char* argv[]) {
    int n = 250, m = 60, len = 6, k = 80; double pKeep = 0.1;
    if (argc >= 6) {
        n      = std::atoi(argv[1]);
        m      = std::atoi(argv[2]);
        len    = std::atoi(argv[3]);
        pKeep  = std::atof (argv[4]);
        k      = std::atoi(argv[5]);
    }
    Case cs = genRnd(n, m, len, pKeep, k);

    /* --- 预处理必选元素 K --- */
    vector<int> old2new(n, -1), new2old;
    for (int i = 0; i < n; ++i)
        if (!cs.K[i]) { old2new[i] = new2old.size(); new2old.push_back(i); }
    int kPrime = std::max(0, k - int(std::count(cs.K.begin(), cs.K.end(), 1)));

    vector<vector<int>> forb2;
    bool infeasible = false;
    for (auto& F : cs.S) {
        vector<int> tmp;
        for (int x : F) if (!cs.K[x]) tmp.push_back(old2new[x]);
        if (tmp.empty()) { infeasible = true; break; }
        forb2.push_back(std::move(tmp));
    }
    if (infeasible) { std::puts("UNSAT: K itself violates a forbidden set"); return 0; }

    std::printf("n=%d  |K|=%zu  n'=%zu  k=%d → k'=%d  |S|=%zu\n",
        n, std::count(cs.K.begin(), cs.K.end(), 1),
        new2old.size(), k, kPrime, cs.S.size());

    /* --- Baseline-MU ---------------------------------------------------- */
    Streamer stMU;
    long long tMU = ms([&]{
        MU mu(new2old.size(), kPrime, forb2, new2old, stMU);
        mu.run();
    });

    /* --- FC-Prop -------------------------------------------------------- */
    Streamer stFC;
    long long tFC = ms([&]{
        FC fc(new2old.size(), kPrime, forb2, new2old, stFC);
        fc.run();
    });

    std::printf("\nBaseline-MU : %lld ms   cnt=%zu   hash=0x%016llx\n",
                tMU, stMU.cnt, (unsigned long long)stMU.hash);
    std::printf(  "FC-Prop     : %lld ms   cnt=%zu   hash=0x%016llx\n",
                tFC, stFC.cnt, (unsigned long long)stFC.hash);
    std::puts( (stMU.cnt==stFC.cnt && stMU.hash==stFC.hash)?
               "SAME RESULT" : "DIFF RESULT" );
    return 0;
}