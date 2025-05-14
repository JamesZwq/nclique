#!/usr/bin/env python3
# chmod +x shrink_until_min_fail.py
import subprocess, tempfile, os, pathlib, shutil

ORIG = "/Users/zhangwenqian/UNSW/pivoter/small_garph.txt_13_K"
BIN1 = "/Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/degeneracy_cliques"
BIN2 = "/Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/main"
BIN3 = "/Users/zhangwenqian/UNSW/nucleus/nd/nucleus"

# ---------- 基础 I/O ------------- #
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

# ---------- 运行三段式管道 -------- #
def pipeline(gfile) -> int:
    """返回码：0=OK(输出相同)，1=diff 不同，其它>1=程序崩溃"""
    cmd = (
        f"{BIN1} -i {gfile} -t V -d 1 -k 0 && "
        f"{BIN2} {gfile}.tree 2 4 {gfile} > out && "
        f"{BIN3} {gfile} 24 no && diff a b"
    )
    return subprocess.run(cmd, shell=True).returncode

# ---------- 主过程 --------------- #
def main():
    n, edges = read_graph(ORIG)
    print(f"起始图：{n} 点, {len(edges)} 条边")

    # 先确保原始图确实 FAIL；若没有 diff 就没必要缩
    tmp0 = tempfile.mktemp(suffix=".txt", dir=str(pathlib.Path(ORIG).parent))
    write_graph(tmp0, n, edges)
    if pipeline(tmp0) == 0:
        print("❗️ 原始图已经 PASS（diff=0），没有要缩的内容")
        os.remove(tmp0)
        return
    os.remove(tmp0)

    round_no = 0
    while True:
        round_no += 1
        changed = False
        print(f"\n=== 第 {round_no} 轮扫描（当前 {len(edges)} 边） ===")
        for idx in range(len(edges)):
            # 每次都新写文件（简化处理 .tree）
            new_edges = edges[:idx] + edges[idx+1:]
            gfile = tempfile.mktemp(suffix=".txt", dir=str(pathlib.Path(ORIG).parent))
            write_graph(gfile, n, new_edges)

            rc = pipeline(gfile)

            # 清理生成的 .tree、临时图文件
            # for ext in ("", ".tree"):
            #     try:
            #         os.remove(gfile + ext)
            #     except FileNotFoundError:
            #         pass
            subprocess.run(f"rm -f {gfile} {gfile}.tree {gfile}.bin {gfile}._24_K", shell=True)


            if rc != 0:               # 删掉后仍 FAIL → 真能删
                    print(f"  [-] 删除边 {edges[idx]} ⇒ 仍 FAIL, 保留删除")
                    # writh remove edge to file
                    with open("/Users/zhangwenqian/UNSW/pivoter/remove.txt", "a") as f:
                        f.write(f"{edges[idx]}\n")
                    edges = new_edges     # 永久更新
                    changed = True
                    break                 # 边集变了，重新开始新一轮

        if not changed:
            print("\n✅ 找不到可进一步删除的边，缩减完成")
            break

    out_path = "minimal_fail.txt"
    write_graph(out_path, n, edges)
    print(f"最小失败用例写入 → {out_path}  (共 {len(edges)} 条边)")

if __name__ == "__main__":
    main()