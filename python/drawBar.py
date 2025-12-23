# -*- coding: utf-8 -*-
from typing import List, Set, Tuple, Optional

# ---------------- I/O：从 .edge 文件读取 ----------------

def read_undirected_edge_file(path: str) -> Tuple[int, List[Set[int]]]:
    """
    读取无向图：第一行 'n m'，后续每行 'u v'（0-based）。
    返回：n, adj（长度 n 的 list[set[int]]）
    """
    with open(path, 'r') as f:
        first = f.readline()
        if not first:
            raise ValueError("Empty file")
        parts = first.strip().split()
        if len(parts) < 1:
            raise ValueError("First line must contain n [m]")
        n = int(parts[0])
        # m 可忽略；文件里行数不一定严格等于 m
        adj: List[Set[int]] = [set() for _ in range(n)]

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            sp = line.split()
            if len(sp) < 2:
                continue
            u = int(sp[0]); v = int(sp[1])
            if u == v:
                continue
            # 简单的范围检查（可注释掉以追求速度）
            if not (0 <= u < n and 0 <= v < n):
                raise ValueError(f"vertex out of range: ({u},{v}) with n={n}")
            adj[u].add(v)
            adj[v].add(u)
    return n, adj

# ---------------- 退化序 / 前向邻域 ----------------

def degeneracy_ordering(adj: List[Set[int]]) -> Tuple[List[int], List[int]]:
    """
    退化序（数组实现）。返回：
    - ordering: 按退化序移除的顶点序列（长度 n）
    - rank: 顶点 -> 在 ordering 中的位置（0..n-1）
    """
    n = len(adj)
    deg = [len(adj[v]) for v in range(n)]
    maxd = max(deg) if n > 0 else 0
    buckets: List[Set[int]] = [set() for _ in range(maxd + 1)]
    for v in range(n):
        buckets[deg[v]].add(v)

    ordering: List[int] = []
    removed = [False] * n
    cur = 0
    remaining = n

    while remaining > 0:
        while cur <= maxd and not buckets[cur]:
            cur += 1
        if cur > maxd:  # 图可能有孤立点/空桶，这里兜底
            # 直接把未移除的顶点加入（度可能已变负）
            for v in range(n):
                if not removed[v]:
                    ordering.append(v)
                    removed[v] = True
                    remaining -= 1
            break
        v = buckets[cur].pop()
        if removed[v]:
            continue
        ordering.append(v)
        removed[v] = True
        remaining -= 1
        # 降低邻居度
        for nb in list(adj[v]):
            if not removed[nb]:
                d = deg[nb]
                if nb in buckets[d]:
                    buckets[d].remove(nb)
                deg[nb] = d - 1
                buckets[d - 1].add(nb)

    rank = [0] * n
    ordering = range(n)
    for i, v in enumerate(ordering):
        rank[i] = i
    return ordering, rank

def bron_kerbosch_noX(adj: List[Set[int]]) -> List[List[int]]:
    """
    Bron–Kerbosch (No-X) with degeneracy ordering.
    输入:
        adj: 长度为 n 的无向邻接表, 顶点为 0..n-1, 且 adj[v] 不包含 v 自身
    输出:
        所有在该搜索策略下得到的团（按伪代码行为实现；注意这版无 X）
    """
    n = len(adj)

    # ---- 退化序 (ordering) 与 rank[v] ----
    deg = [len(adj[v]) for v in range(n)]
    maxd = max(deg) if n > 0 else 0
    buckets: List[Set[int]] = [set() for _ in range(maxd + 1)]
    for v in range(n):
        buckets[deg[v]].add(v)

    ordering: List[int] = []           # 退化序 v1..vn
    removed = [False] * n
    cur = 0
    remaining = n
    while remaining > 0:
        while cur <= maxd and not buckets[cur]:
            cur += 1
        if cur > maxd:
            # 兜底：如果出现空桶问题，直接把未移除的都推入
            for v in range(n):
                if not removed[v]:
                    ordering.append(v)
                    removed[v] = True
                    remaining -= 1
            break
        v = buckets[cur].pop()
        if removed[v]:
            continue
        ordering.append(v)
        removed[v] = True
        remaining -= 1
        # 降低邻居度
        for nb in list(adj[v]):
            if not removed[nb]:
                d = deg[nb]
                if nb in buckets[d]:
                    buckets[d].remove(nb)
                deg[nb] = d - 1
                buckets[d - 1].add(nb)

    rank = [0] * n
    for i, v in enumerate(ordering):
        rank[v] = i

    # ---- 根层前向邻域 N^+(v) ----
    N_plus: List[Set[int]] = [set() for _ in range(n)]
    for v in range(n):
        rv = rank[v]
        N_plus[v] = {w for w in adj[v] if rank[w] > rv}

    # ---- 递归枚举 ----
    cliques: List[List[int]] = []

    def BK(R: Set[int], P: Set[int]):
        # 与伪代码一致：P 为空则输出 R
        if not P:
            cliques.append(sorted(R))
            return
        # 枢轴: u_piv = argmax_{x in P} |P ∩ N(x)|  （注意这里用无向 N）
        u_piv = max(P, key=lambda x: len(P & adj[x]))
        # 仅对 P \ N(u_piv) 分支
        for v in list(P - adj[u_piv]):
            if v not in P:   # 被同层更新移除了就跳过
                continue
            P.remove(v)
            BK(R | {v}, P & adj[v])

    # 根循环：按退化序，每个 v 作为种子，P = N^+(v)
    for v in ordering:
        BK({v}, set(N_plus[v]))

    return cliques

def bron_kerbosch_with_X(adj: List[Set[int]]) -> List[List[int]]:
    """
    Bron–Kerbosch with X-set, seeded per vertex using degeneracy order.
    与论文伪代码一致：对每个根 v 按 σ 递增顺序，P=N^+(v), X=N(v)\\N^+(v)，递归中 pivot 从 P∪X 选。
    返回所有极大团（每个一次）。
    """
    n = len(adj)
    ordering, rank = degeneracy_ordering(adj)

    # 前向邻域 N^+(v)
    N_plus: List[Set[int]] = [set() for _ in range(n)]
    for v in range(n):
        rv = rank[v]
        N_plus[v] = {w for w in adj[v] if rank[w] > rv}

    cliques: List[List[int]] = []

    def BK(R: Set[int], P: Set[int], X: Set[int]):
        if not P and not X:
            cliques.append(sorted(R))
            return
        # Tomita pivot: u ∈ P∪X maximizing |P ∩ N(u)|
        if P or X:
            u = max(P | X, key=lambda x: len(P & adj[x]))
            branch = list(P - adj[u])
        else:
            branch = list(P)
        for v in branch:
            BK(R | {v}, P & adj[v], X & adj[v])
            P.discard(v)
            X.add(v)

    # one root per vertex, in degeneracy order
    for v in ordering:
        R = {v}
        P = set(N_plus[v])
        X = set(adj[v]) - P  # N^-(v)
        BK(R, P, X)

    return cliques

def forward_neighbors(adj: List[Set[int]], rank: List[int]) -> List[Set[int]]:
    """
    前向邻域 N^+(v) = { w in N(v) | rank[w] > rank[v] } （数组实现）
    """
    n = len(adj)
    N_plus: List[Set[int]] = [set() for _ in range(n)]
    for v in range(n):
        rv = rank[v]
        # 只保留 rank 更大的邻居
        N_plus[v] = {w for w in adj[v] if rank[w] > rv}
    return N_plus

# ---------------- Build–SCT 主体 ----------------

def build_sct_from_adj(adj: List[Set[int]], k_max: Optional[int] = None):
    """
    构建 SCT。返回所有根-叶路径：list of dicts:
        {"hold": set[int], "pivot": set[int]}
    备注：k_max 预留，用于未来剪枝（此实现未使用）。
    """
    n = len(adj)
    ordering, rank = degeneracy_ordering(adj)
    N_plus = forward_neighbors(adj, rank)
    N = adj  # 无向邻居

    paths = []

    def expand(P: Set[int], V_hold: Set[int], V_pivot: Set[int]):
        if not P:
            paths.append({"hold": set(V_hold), "pivot": set(V_pivot)})
            return
        # pivot: argmax_{x in P} |P ∩ N^+(x)|
        u_piv = max(P, key=lambda x: len(P & N_plus[x]))
        # 非邻分支集合：P \ N(u_piv)  （注意：N 为无向邻域）
        branch_vertices = list(P - N[u_piv])
        for v in branch_vertices:
            if v not in P:
                continue
            P.remove(v)  # 本层永久移除 v
            P_next = P & N[v]
            if v == u_piv:
                # 将 pivot 延期：加入 V_pivot（不进 hold）
                expand(P_next, V_hold, V_pivot | {v})
            else:
                # 非邻分支：v 加入 hold
                expand(P_next, V_hold | {v}, V_pivot)

    # 根：按退化序，每个 v 作为根；候选集为 N^+(v)
    for v in ordering:
        V_hold = {v}
        V_pivot = set()
        P0 = set(N_plus[v])
        expand(P0, V_hold, V_pivot)

    return paths


def build_sct_from_file(path: str, k_max: Optional[int] = None):
    n, adj = read_undirected_edge_file(path)
    return build_sct_from_adj(adj, k_max=k_max)


# ---------------- 简单命令行 ----------------
if __name__ == "__main__":
    import sys
    # if len(sys.argv) >= 2:
    #     for fp in sys.argv[1:]:
    #         sct = build_sct_from_file("/Users/zhangwenqian/UNSW/pivoter/a.edge")
    #         print(f"{fp}: #paths={len(sct)}")
    # else:
    #     print("Usage: python build_sct.py <graph.edge> [...]", file=sys.stderr)
    sct = build_sct_from_file("/Users/zhangwenqian/UNSW/pivoter/a.edge")
    bk = bron_kerbosch_with_X(read_undirected_edge_file("/Users/zhangwenqian/UNSW/pivoter/a.edge")[1])
    print(f"#paths={len(sct)}")
    for i, path in enumerate(sct):
        print(f"Path {i}: hold={sorted(path['hold'])}, pivot={sorted(path['pivot'])}")

    print(f"#bron_kerbosch_noX cliques={len(bk)}")
    for i, c in enumerate(bk):
        print(f"Clique {i}: {c}")