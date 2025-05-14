#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

import subprocess
outputFile = Path("/Users/zhangwenqian/UNSW/pivoter/b.tmp")

def main():
    # 读取输入（命令行优先，其次标准输入）
    # if len(sys.argv) > 1:
    #     n = int(sys.argv[1])
    # else:
    #     n = int(sys.stdin.readline().strip())
    n = 40
    if not (1 <= n < 400):
        sys.exit("n 必须是 1‒399 之间的整数")

    m = n * (n - 1) // 2  # 边数

    orgEdges = []
    # 逐行写入文件
    # with outputFile.open("w") as f:
    for u in range(n - 1):
        for v in range(u + 1, n):
            # f.write(f"{u} {v}\n")  # 每条边一行
            orgEdges.append((u, v))

    # random remove 10% edges, write removed edge to file /Users/zhangwenqian/UNSW/pivoter/b
    # 10% of m
    import random
# /Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/BCtest && /Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/degeneracy_cliques -i /Users/zhangwenqian/UNSW/pivoter/b.tmp -t V -d 1
#     random.seed(0)
    while True:
        num_edges_to_remove = m // 10
        removed_edges = random.sample(orgEdges, num_edges_to_remove)
        edges = [edge for edge in orgEdges if edge not in removed_edges]
        # write removed edges to file
        with open("/Users/zhangwenqian/UNSW/pivoter/b", "w") as f:
            f.write(f"{n}\n")
            for edge in removed_edges:
                f.write(f"{edge[0]} {edge[1]}\n")
            # write remaining edges to file

        with outputFile.open("w") as f:
            f.write(f"{n} {len(edges)}\n")
            for edge in edges:
                f.write(f"{edge[0]} {edge[1]}\n")
        # 运行命令

        cmds = [
            "/Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/BCtest",
            "/Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/degeneracy_cliques -i /Users/zhangwenqian/UNSW/pivoter/b.tmp -t V -d 1 -k 0",
            "diff outA.txt outB.txt",   # 注意：diff 返回码 1 仅表示文件不同，不一定是错误
        ]

        for cmd in cmds:
            print(f"\n→ 正在执行: {cmd}")
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            # 单步输出
            if result.stdout:
                print("STDOUT:")
                print(result.stdout)
            if result.stderr:
                print("STDERR:")
                print(result.stderr)

            # 判断是否真正出错
            if cmd.startswith("diff"):
                # diff：0=相同，1=不同，2=有问题
                if result.returncode == 0:
                    print("diff 结果：文件相同")
                elif result.returncode == 1:
                    print("diff 结果：文件不同")
                    sys.exit(1)
                else:  # 2
                    print("diff 发生错误，终止脚本")
                    sys.exit(1)
            else:
                if result.returncode != 0:
                    print(f"命令失败 (returncode={result.returncode})，终止脚本")
                    sys.exit(1)


if __name__ == "__main__":
    main()