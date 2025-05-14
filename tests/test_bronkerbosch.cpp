// test_bronkerbosch.cpp
#include "tree/BronKerbosch.h"

#include <cassert>
#include <iostream>
#include <set>
#include <algorithm>
#include <vector>
#include <Global/Global.h>

// --------- 帮助函数 ---------
template<class Bitset>
std::vector<int> bitset_to_vec(const Bitset& bs) {
    std::vector<int> res;
    bk::for_each_bit(bs, (int)bs.size(), [&](int v){
        res.push_back(v);
        return true;
    });
    std::sort(res.begin(), res.end());
    return res;
}

bool eq_multisets(std::vector<std::vector<int>> a,
                  std::vector<std::vector<int>> b) {
    for (auto& v: a) std::sort(v.begin(), v.end());
    for (auto& v: b) std::sort(v.begin(), v.end());
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    return a == b;
}

// --------- 运行单个案例 ---------
void run_case(
        const std::string& name,
        std::vector<TreeGraphNode> vList,
        daf::StaticVector<std::pair<daf::Size, daf::Size>> removeEdgeList,
        std::vector<std::vector<int>> expect)
{
    std::vector<std::vector<int>> got;
    bk::bronKerbosch(
        vList, removeEdgeList, /*minK=*/1,
        [&](const bk::Bitset& R, const bk::Bitset& /*pivots*/){
            got.emplace_back(bitset_to_vec(R));
        });

    assert(eq_multisets(got, expect) && "cliques mismatch");
    std::cout << "[PASS] " << name << '\n';
}

// --------- main ---------
int main() {
    // --------- 案例 1: K4 ---------
    {
        std::vector<TreeGraphNode> vList;
        for (int i = 0; i < 4; ++i) vList.emplace_back(i, false);
        daf::StaticVector<std::pair<daf::Size, daf::Size>> rm;
        run_case("K4", std::move(vList), rm, {{0,1,2,3}});
    }

    // --------- 案例 2: K4 - edge{0,1} ---------
    {
        std::vector<TreeGraphNode> vList;
        for (int i = 0; i < 4; ++i) vList.emplace_back(i, false);
        daf::StaticVector<std::pair<daf::Size, daf::Size>> rm;
        rm.emplace_back(0,1);
        run_case("K4 minus (0,1)", std::move(vList), rm,
                 {{0,2,3}, {1,2,3}});
    }

    // --------- 案例 3: triangle + isolated ---------
    {
        std::vector<TreeGraphNode> vList;
        for (int i = 0; i < 4; ++i) vList.emplace_back(i, false);
        daf::StaticVector<std::pair<daf::Size, daf::Size>> rm;
        rm.emplace_back(0,3);
        rm.emplace_back(1,3);
        rm.emplace_back(2,3);
        run_case("triangle + isolated", std::move(vList), rm,
                 {{0,1,2}, {3}});
    }

    std::cout << "All tests passed! 🎉\n";
    return 0;
}