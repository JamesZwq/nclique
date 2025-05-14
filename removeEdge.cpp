#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <set>
#include <algorithm>
#include <fstream>
using namespace std;

// 按 (u,v) 字典序排序
static bool cmpEdge(const pair<int,int>& a, const pair<int,int>& b) {
    return a.first < b.first || (a.first == b.first && a.second < b.second);
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        cerr << "用法: " << argv[0] << " 原图文件   待删边文件   输出文件 topNin待删边";
        return 1;
    }
    const char* origFile = argv[1];
    const char* delFile  = argv[2];
    const char* outFile  = argv[3];

    int N;
    long long E;
    vector<pair<int,int>> edges;
    vector<pair<int,int>> delEdges;

    // —— 1. 读入原图
    {
        ifstream fin(origFile);
        fin >> N >> E;
        edges.reserve(E);
        int u, v;
        while (fin >> u >> v) {
            edges.emplace_back(u, v);
        }
    }

    // —— 2. 读入待删边
    {
        // argv[4];
        int num = stoi(argv[4]);
        ifstream fin(delFile);
        // 不需要第一行头；直接读取所有 (u,v)
        int u, v, k;
        while (fin >> u >> v >> k) {
            delEdges.emplace_back(u, v);
            if (num-- <= 0) break;
        }
    }

    // —— 3. 在内存中对两组边分别排序
    sort(edges.begin(),    edges.end(),    cmpEdge);
    sort(delEdges.begin(), delEdges.end(), cmpEdge);

    // —— 4. 双指针归并：生成保留边
    vector<pair<int,int>> kept;
    kept.reserve(edges.size());
    size_t i = 0, j = 0;
    while (i < edges.size() && j < delEdges.size()) {
        const auto &ei = edges[i];
        const auto &dj = delEdges[j];
        if (ei < dj) {
            kept.push_back(ei);
            ++i;
        }
        else if (dj < ei) {
            ++j;
        }
        else {
            // 相等 => 删除
            ++i; ++j;
        }
    }
    // 剩余没有被删除的部分
    while (i < edges.size()) {
        kept.push_back(edges[i++]);
    }

    // —— 5. 写出结果
    ofstream fout(outFile);
    fout << N << " " << kept.size() << "\n";
    for (auto &e : kept) {
        fout << e.first << " " << e.second << "\n";
    }

    return 0;
}