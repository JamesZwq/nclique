#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <set>  // 新增

// Boost 头文件
#include <boost/heap/d_ary_heap.hpp>
#include <boost/heap/pairing_heap.hpp>
#include <boost/heap/fibonacci_heap.hpp>

struct Node {
    int key;
};

// 用于 heap 和 multiset 的比较器（最小堆 / 最小有序集）
struct Compare {
    bool operator()(Node const& a, Node const& b) const {
        return a.key > b.key;
    }
};

// 通用的 Boost.Heap 基准函数
template<typename Heap>
void benchmark_heap(const std::string& name, int N, int M, std::mt19937_64& gen) {
    using Handle = typename Heap::handle_type;
    std::uniform_int_distribution<int> dist(0, N * 10);

    Heap heap;
    std::vector<Handle> handles;
    handles.reserve(N);

    // 插入
    auto t1 = std::chrono::steady_clock::now();
    for (int i = 0; i < N; ++i) {
        handles.push_back(heap.push(Node{ dist(gen) }));
    }
    auto t2 = std::chrono::steady_clock::now();

    // decrease-key
    auto t3 = std::chrono::steady_clock::now();
    for (int i = 0; i < M; ++i) {
        int idx = gen() % N;
        int new_key = dist(gen);
        Node &node = *handles[idx];
        if (new_key < node.key) {
            heap.update(handles[idx], Node{ new_key });
        }
    }
    auto t4 = std::chrono::steady_clock::now();

    // pop-all
    while (!heap.empty()) heap.pop();
    auto t5 = std::chrono::steady_clock::now();

    auto to_ms = [&](auto d){
        return std::chrono::duration_cast<std::chrono::milliseconds>(d).count();
    };
    std::cout << name
              << " | insert: "      << to_ms(t2 - t1) << " ms"
              << " | decrease-key: " << to_ms(t4 - t3) << " ms"
              << " | pop-all: "      << to_ms(t5 - t4) << " ms"
              << " | total: "      << to_ms(t5 - t1) << " ms"
              << "\n";
}

// 专用于 std::multiset 的基准函数
void benchmark_set(const std::string& name, int N, int M, std::mt19937_64& gen) {
    std::uniform_int_distribution<int> dist(0, N * 10);
    std::multiset<Node, std::function<bool(Node const&, Node const&)>> s(
        [](Node const& a, Node const& b){ return a.key < b.key; }
    );
    std::vector<std::multiset<Node, std::function<bool(Node const&, Node const&)>>::iterator> handles;
    handles.reserve(N);

    // 插入
    auto t1 = std::chrono::steady_clock::now();
    for (int i = 0; i < N; ++i) {
        handles.push_back(s.insert(Node{ dist(gen) }));
    }
    auto t2 = std::chrono::steady_clock::now();

    // “减少键” —— erase + insert
    auto t3 = std::chrono::steady_clock::now();
    for (int i = 0; i < M; ++i) {
        int idx = gen() % N;
        int new_key = dist(gen);
        auto it = handles[idx];
        if (new_key < it->key) {
            s.erase(it);
            handles[idx] = s.insert(Node{ new_key });
        }
    }
    auto t4 = std::chrono::steady_clock::now();

    // pop-all
    while (!s.empty()) s.erase(s.begin());
    auto t5 = std::chrono::steady_clock::now();

    auto to_ms = [&](auto d){
        return std::chrono::duration_cast<std::chrono::milliseconds>(d).count();
    };
    std::cout << name
              << " | insert: "      << to_ms(t2 - t1) << " ms"
              << " | decrease-key: " << to_ms(t4 - t3) << " ms"
              << " | pop-all: "      << to_ms(t5 - t4) << " ms"
              << " | total: "      << to_ms(t5 - t1) << " ms"
              << "\n";
}

int main() {
    const int N = 500000;   // 插入元素个数
    const int M = 2000000;   // decrease-key 操作次数
    std::mt19937_64 gen(123456);
    std::cout << "N = " << N << ", M = " << M << "\n";
    // d-ary heap (d=4)
    using DHeap = boost::heap::d_ary_heap<
        Node,
        boost::heap::arity<4>,
        boost::heap::mutable_<true>,
        boost::heap::compare<Compare>
    >;
    benchmark_heap<DHeap>("d-ary-heap (d=4)", N, M, gen);
    using DHeap8 = boost::heap::d_ary_heap<
        Node,
        boost::heap::arity<8>,
        boost::heap::mutable_<true>,
        boost::heap::compare<Compare>
    >;
    benchmark_heap<DHeap8>("d-ary-heap (d=8)", N, M, gen);
    // pairing heap
    using PHeap = boost::heap::pairing_heap<
        Node,
        boost::heap::mutable_<true>,
        boost::heap::compare<Compare>
    >;
    benchmark_heap<PHeap>("pairing-heap",           N, M, gen);

    // fibonacci heap
    using FHeap = boost::heap::fibonacci_heap<
        Node,
        boost::heap::mutable_<true>,
        boost::heap::compare<Compare>
    >;
    benchmark_heap<FHeap>("fibonacci-heap",         N, M, gen);

    // C++ 标准库 multiset 对比
    benchmark_set("std::multiset (erase+insert)",  N, M, gen);

    return 0;
}