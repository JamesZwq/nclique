from typing import List, Tuple, Set

def enumerate_valid_maximal_sets(v_list, conflict_sets):
    keep   = {v for v, isP in v_list if not isP}
    pivots = [v for v, isP in v_list if  isP]
    involved = {v: [i for i,F in enumerate(conflict_sets) if v in F] for v in pivots}

    F_size  = [len(F) for F in conflict_sets]
    counter = [len(F & keep) for F in conflict_sets]

    res, cur, cur_pivs = [], set(keep), set()

    def can_add(v):
        """ 把 v 放进来会不会装满某个 F？ """
        return all(counter[i] + 1 < F_size[i] for i in involved[v])

    def is_maximal():
        """ 没有任何 pivot 再加进来仍合法 => 当前 clique 是极大 """
        return all(not can_add(v) for v in pivots if v not in cur_pivs)

    def dfs(idx):
        # 若已冲突，剪枝
        if any(counter[i] == F_size[i] for i in range(len(conflict_sets))):
            return
        if idx == len(pivots):
            if is_maximal():
                res.append((cur.copy(), cur_pivs.copy()))
            return

        v = pivots[idx]

        # 选 v
        if can_add(v):
            for i in involved[v]: counter[i] += 1
            cur.add(v); cur_pivs.add(v)
            dfs(idx+1)
            cur.remove(v); cur_pivs.remove(v)
            for i in involved[v]: counter[i] -= 1

        # 不选 v
        dfs(idx+1)

    pivots.sort(key=lambda v: -len(involved[v]))
    dfs(0)
    return res

# --------------------------- Demo ------------------------------------------
if __name__ == "__main__":
    v_list = [
        (1, False),
        (2, False), (4, True), (5, True), (6, True), (9, True)   # pivots
    ]
    conflict_sets = [
        {5, 6, 9},
    ]

    valid = enumerate_valid_maximal_sets(v_list, conflict_sets)

    print(f"总计 {len(valid)} 个合法 (Clique, PivotSubset)：\n")
    for clique, pivs in valid:
        holds = clique - pivs
        print(f"Clique: {sorted(clique)}  Pivots: {sorted(pivs)}  Holds: {sorted(holds)}")