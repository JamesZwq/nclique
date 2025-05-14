#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import networkx as nx
from typing import List, Tuple

level = 0
def enumerate_max_cliques(vertices: List[int],
                          deleted_edges: List[Tuple[int,int]]) -> List[List[int]]:
    """
    从完全图中删去 deleted_edges 后，枚举所有极大团。
    采用“覆盖/不覆盖”回溯，参见上文示例。
    """
    V = sorted(vertices)
    # 只保留真正要删掉的边
    adjList = []
    for i in range(len(V)):
        adjList.append([])
    # sV = set(V)

    isInDeleted = [False] * len(V)
    for u,v in deleted_edges:
        # if u in sV and v in sV:
        if u<v: u,v = v,u
        adjList[u].append(v)
        isInDeleted[u] = True
        isInDeleted[v] = True
    print("removed AdjList:", adjList)
    startR = []
    startP = []

    for v, deleted in enumerate(isInDeleted):
        # 只保留真正要删掉的边
        if deleted:
            startP.append(v)
        else :
            startR.append(v)
    out = []
    print("startR:", startR)
    print("startP:", startP)

    def dfs(R, P):
        global level
        """
        递归回溯函数。
        R: 当前团
        P: 剩余顶点
        """
        if P == set():
            out.append(R)
            level -= 1
            return

        # 选择最小的顶点作为 pivot

        # 遍历剩余顶点
        did = set()
        for v in P:
            level += 1
            did.add(v)
            dfs(R + [v], P - set(adjList[v]) - did)

        level -= 1



    # 从最小顶点开始回溯
    dfs(startR, set(startP))
    print(out)
    return out

def verify(vertices: List[int],
           deleted_edges: List[Tuple[int,int]]) -> bool:
    """
    对比我们算法和 NetworkX 在删边之后的极大团结果。
    返回 True 表示一致，False 表示不一致并打印差异。
    """
    our = enumerate_max_cliques(vertices, deleted_edges)

    # 用 NetworkX 构造剩余子图
    G = nx.Graph()
    G.add_nodes_from(vertices)
    for i in vertices:
        for j in vertices:
            if i < j and (i,j) not in deleted_edges and (j,i) not in deleted_edges:
                G.add_edge(i, j)

    nx_cliques = sorted([sorted(c) for c in nx.find_cliques(G)])

    our_set = {tuple(c) for c in our}
    nx_set  = {tuple(c) for c in nx_cliques}

    if our_set != nx_set:
        print("❌ 差异检测！")
        print("输入顶点：", vertices)
        print("删除的边：", deleted_edges)
        print("\n—— 我们算法的结果 ——")
        for c in our:
            print("  ", c)
        print("\n—— NetworkX 的结果 ——")
        for c in nx_cliques:
            print("  ", c)
        print("\n  仅在我们算法中：", our_set - nx_set)
        print("  仅在 NetworkX 中：", nx_set - our_set)
        return False

    return True

if __name__ == "__main__":
    # 先执行几个固定删除顺序的示例验证
    clique = [0,1,2,3,4,5,6,7]
    removes = [(2, 4), (3, 4)]
    print("clique:", clique)
    print("removes:", removes)

    enumerate_max_cliques(clique, removes)


    # 1 5 6 7
    # 1 5 6 7 2 3
    # 1 5 6 7 4
    # deleted = []
    # for e in removes:
    #     deleted.append(e)
    #     print(f"\n=== 验证：删除边 {e} 后 ===")
    #     ok = verify(clique, deleted)
    #     if not ok:
    #         exit(1)
    # print("\n固定案例全部通过！")

    # # --- 随机测试 1 000 次 ---
    # print("\n开始 1000 次随机测试...")
    # random.seed(42)
    # for t in range(1, 1001):
    #     k = random.randint(2, 10)           # 随机 clique 大小
    #     verts = list(range(1, k+1))         # 顶点编号 1..k
    #     # clique 完全图上的所有边
    #     all_edges = [(i,j) for i in verts for j in verts if i<j]
    #     # 随机删掉 0..len(all_edges) 条
    #     m = random.randint(0, len(all_edges))
    #     deleted = random.sample(all_edges, m)
    #
    #     if not verify(verts, deleted):
    #         print(f"\n在第 {t} 次随机测试时失败，见上面输出。")
    #         exit(1)
    #     if t % 100 == 0:
    #         print(f"  已完成 {t} 次……")
    # print("\n🎉 随机测试 1000 次全部通过！")