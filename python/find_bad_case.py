#!/usr/bin/env python3
# chmod +x find_bad_case.py
import random, subprocess, tempfile, shutil, os, pathlib, sys

ORIG   = "/Users/zhangwenqian/UNSW/pivoter/dblp_removed.edge"
BIN1   = "/Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/degeneracy_cliques"
BIN2   = "/Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/main"
BIN3   = "/Users/zhangwenqian/UNSW/nucleus/nd/nucleus"

def read_graph(path):
    with open(path) as f:
        n, _ = map(int, f.readline().split())
        edges = [tuple(map(int, line.split())) for line in f]
    return n, edges

def write_graph(path, n, edges):
    with open(path, "w") as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")

def make_case(src, keep_ratio=1):
    """
    keep_ratio=0.1 ⇒ 仅保留 10 % 的现有边，其余全部删除
    """
    n, old_edges = read_graph(src)
    k = max(1, int(len(old_edges) * keep_ratio))       # 至少保留 1 条，避免空图
    new_edges = random.sample(old_edges, k)

    tmp = "/Users/zhangwenqian/UNSW/pivoter/graphs/case.txt"
    write_graph(tmp, n, new_edges)
    return tmp

def pipeline(gfile):
    cmd = (
        f"{BIN1} -i {gfile} -t V -d 1 -k 4 && "
        f"{BIN2} {gfile}.tree 2 4 {gfile} > out && "
        f"{BIN3} {gfile} 24 no && diff /Users/zhangwenqian/UNSW/pivoter/a /Users/zhangwenqian/UNSW/pivoter/a.tmp"
    )
    return subprocess.run(cmd, shell=True).returncode


pipeline(ORIG)
def main(max_try=1000000, keep_ratio=1):
    for i in range(1, max_try + 1):
        print(f"[+] 试第 {i} 次，保留 {keep_ratio*100:.1f}% 边")
        case = make_case(ORIG, keep_ratio)
        if pipeline(case) != 0:                      # 触发 diff
            dst = f"interesting_case_{i}.txt"
            shutil.copy(case, dst)
            print(f"[+] 找到失败用例（保留 {keep_ratio*100:.1f}% 边）→ {dst}")
            break
        else:
            os.remove(case)
    else:
        print("[-] 试完也没触发 diff，可调低 keep_ratio 或增大 max_try")

if __name__ == "__main__":
    # 用法:  python3 find_bad_case.py [max_try] [keep_ratio]
    main(*(map(float, sys.argv[1:])) if len(sys.argv) > 1 else [])