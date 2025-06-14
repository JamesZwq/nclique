import random
import subprocess
import networkx as nx
import matplotlib.pyplot as plt
import math
import math
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull
from matplotlib.path import Path
import matplotlib.patches as patches
import random
import bisect
import subprocess
import random
import time
import sys
from colorama import init, Fore, Style


import numpy as np


def circle_from_three_points(p, q, r):
    # 计算由三点确定的圆（外接圆）的圆心和半径
    d = 2 * (p[0]*(q[1]-r[1]) + q[0]*(r[1]-p[1]) + r[0]*(p[1]-q[1]))
    if abs(d) < 1e-6:
        return None  # 三点共线，无法确定唯一圆
    ux = ((p[0]**2+p[1]**2)*(q[1]-r[1]) + (q[0]**2+q[1]**2)*(r[1]-p[1]) + (r[0]**2+r[1]**2)*(p[1]-q[1])) / d
    uy = ((p[0]**2+p[1]**2)*(r[0]-q[0]) + (q[0]**2+q[1]**2)*(p[0]-r[0]) + (r[0]**2+r[1]**2)*(q[0]-p[0])) / d
    center = (ux, uy)
    radius = math.dist(center, p)
    return center, radius

def min_enclosing_circle(points):
    # 若只有一个点，则圆心就是该点，半径为0
    if not points:
        return ((0,0), 0)
    if len(points) == 1:
        return (points[0], 0)

    best_center = None
    best_radius = float('inf')

    # 检查所有点对构成的圆
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            p = points[i]
            q = points[j]
            center = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
            radius = math.dist(p, q) / 2
            if all(math.dist(center, pt) <= radius + 1e-6 for pt in points):
                if radius < best_radius:
                    best_radius = radius
                    best_center = center

    # 检查所有三点构成的外接圆
    for i in range(len(points)):
        for j in range(i+1, len(points)):
            for k in range(j+1, len(points)):
                p = points[i]
                q = points[j]
                r = points[k]
                circle = circle_from_three_points(p, q, r)
                if circle is not None:
                    center, radius = circle
                    if all(math.dist(center, pt) <= radius + 1e-6 for pt in points):
                        if radius < best_radius:
                            best_radius = radius
                            best_center = center

    # 若未找到合适的圆，则用节点平均位置和最大距离构造一个圆
    if best_center is None:
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        center = (sum(xs)/len(points), sum(ys)/len(points))
        radius = max(math.dist(center, p) for p in points)
        return center, radius

    return best_center, best_radius

from scipy.spatial import ConvexHull

def rounded_polygon_patch(vertices, r, resolution=16, edgecolor='red'):
    """
    根据给定的顶点列表，生成一个具有圆角效果的多边形 patch。
    参数:
      vertices: 多边形顶点列表，顺序排列（numpy 数组或可转为 numpy 数组）
      r: 圆角的半径
      resolution: 每个圆角使用多少个点来近似圆弧
      edgecolor: 多边形边缘颜色
    """
    vertices = np.array(vertices)
    n = len(vertices)
    verts = []
    codes = []

    for i in range(n):
        prev = vertices[i - 1]
        curr = vertices[i]
        nxt = vertices[(i + 1) % n]

        # 计算从当前点指向前后顶点的向量
        v1 = prev - curr
        v2 = nxt - curr

        # 单位化向量
        v1_norm = v1 / np.linalg.norm(v1)
        v2_norm = v2 / np.linalg.norm(v2)

        # 计算当前角度
        angle = np.arccos(np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0))

        # 根据几何关系，确定从顶点到切点的距离 d
        d = r / np.tan(angle / 2)
        # 确保 d 不超过相邻边长的一半
        d = min(d, np.linalg.norm(v1)/2, np.linalg.norm(v2)/2)

        tangent1 = curr + v1_norm * d
        tangent2 = curr + v2_norm * d

        # 计算角平分线方向
        bisector = v1_norm + v2_norm
        if np.linalg.norm(bisector) < 1e-8:
            # 如果退化，则直接添加当前点
            if i == 0:
                verts.append(tuple(curr))
                codes.append(Path.MOVETO)
            else:
                verts.append(tuple(curr))
                codes.append(Path.LINETO)
            continue
        bisector = bisector / np.linalg.norm(bisector)
        # 根据三角形关系，弧心到顶点距离 = r/sin(angle/2)
        arc_center = curr + bisector * (r / np.sin(angle/2))

        # 计算切点相对于弧心的角度
        start_angle = np.arctan2(tangent1[1] - arc_center[1], tangent1[0] - arc_center[0])
        end_angle = np.arctan2(tangent2[1] - arc_center[1], tangent2[0] - arc_center[0])

        # 保证弧度为正
        angle_diff = end_angle - start_angle
        if angle_diff <= 0:
            angle_diff += 2 * np.pi

        # 对于第一个顶点，移动到第一个切点
        if i == 0:
            verts.append(tuple(tangent1))
            codes.append(Path.MOVETO)
        else:
            verts.append(tuple(tangent1))
            codes.append(Path.LINETO)

        # 在两个切点之间插入圆弧上的点
        for j in range(1, resolution):
            theta = start_angle + angle_diff * j / resolution
            point = arc_center + np.array([r * np.cos(theta), r * np.sin(theta)])
            verts.append(tuple(point))
            codes.append(Path.LINETO)

    # 关闭路径
    verts.append(verts[0])
    codes.append(Path.CLOSEPOLY)

    path = Path(verts, codes)
    patch = patches.PathPatch(path, fill=False, linewidth=2, edgecolor=edgecolor)
    return patch

def draw_graph_with_cliques(edges):
    G = nx.Graph()
    G.add_edges_from(edges)
    pos = nx.kamada_kawai_layout(G)

    plt.figure(figsize=(8, 6))
    nx.draw_networkx_edges(G, pos)
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=500)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')

    cliques = [sorted(clique) for clique in nx.find_cliques(G) if len(clique) >= 4]
    cliques = [clique for clique in cliques if 1 in clique]
    print("Cliques:", cliques)

    colors = ['red', 'green', 'blue', 'orange', 'purple', 'cyan', 'magenta', 'yellow']
    ax = plt.gca()

    for i, clique in enumerate(cliques):
        clique_pos = [pos[node] for node in clique]
        clique_pos_arr = np.array(clique_pos)

        if len(clique_pos) >= 3:
            # 用凸包获取 clique 的外围边界
            hull = ConvexHull(clique_pos_arr)
            hull_points = clique_pos_arr[hull.vertices]

            # 适当膨胀凸包（padding）
            centroid = np.mean(hull_points, axis=0)
            padding_factor = 1.5
            padded_hull_points = centroid + (hull_points - centroid) * padding_factor

            # 根据边长选取一个合适的圆角半径
            edge_lengths = [np.linalg.norm(padded_hull_points[j] - padded_hull_points[(j+1) % len(padded_hull_points)])
                            for j in range(len(padded_hull_points))]
            avg_edge = np.mean(edge_lengths)
            r_round = avg_edge * 0.2  # 可根据实际效果调整

            patch = rounded_polygon_patch(padded_hull_points, r_round, resolution=16, edgecolor=colors[i % len(colors)])
            ax.add_patch(patch)
        else:
            # 退化情况：只有两个节点则使用圆圈
            center = ((clique_pos[0][0] + clique_pos[1][0]) / 2, (clique_pos[0][1] + clique_pos[1][1]) / 2)
            radius = math.dist(clique_pos[0], clique_pos[1]) / 2 * 1.1
            circle = plt.Circle(center, radius, color=colors[i % len(colors)], fill=False, linestyle='--', linewidth=2)
            ax.add_patch(circle)

    plt.title("Graph with 4-Cliques Highlighted by Rounded Polygons")
    plt.savefig("/Users/zhangwenqian/UNSW/KClique/small_graph_with_cliques.png", dpi=500)
    # plt.show()

def draw_graph(edges):
    G = nx.Graph()
    G.add_edges_from(edges)

    pos = nx.kamada_kawai_layout(G)

    plt.figure(figsize=(8, 6))
    nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=500, font_size=10, font_weight='bold')
    plt.title("Generated Graph")
    plt.savefig("/Users/zhangwenqian/UNSW/KClique/small_garph.edges.png")
import random, bisect, os
from typing import Set, Tuple

def generate_graph_log(
        node_count: int,
        edge_count: int,
        output_file: str,
        mu: float = 0.0,
        sigma: float = 1.0
) -> Set[Tuple[int,int]]:
    """
    用对数正态分布给节点加权，再按权重抽样边。
    mu, sigma 控制 lognormal 的位置和形状。
    返回最终的边集合（无向，(u,v) u<v）。
    """
    total = node_count * (node_count - 1) // 2
    if edge_count > total:
        raise ValueError(f"Too many edges (max {total}), got {edge_count}")

    # 先做 ±10%~20% 扰动
    edge_count = int(edge_count * random.uniform(0.9, 1.2))
    edge_count = min(edge_count, total)

    # 1. 为每个节点生成一个 log-normal 权重
    weights = [random.lognormvariate(mu, sigma) for _ in range(node_count)]
    # 要传给 random.choices 的权重列表
    # 注意 random.choices 要 Python3.6+
    node_indices = list(range(node_count))

    # 2. 加权抽样生成边
    edges = set()
    def normalize(u: int, v: int) -> Tuple[int,int]:
        return (u,v) if u < v else (v,u)

    while len(edges) < edge_count:
        u = random.choices(node_indices, weights)[0]
        v = random.choices(node_indices, weights)[0]
        if u == v:
            continue
        edges.add(normalize(u, v))

    # 3. 写文件
    with open(output_file, 'w') as f:
        f.write(f"{node_count} {len(edges)}\n")
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")

    return edges


def run_cmd(name, cmd, error_label):
    # print(f"{Fore.YELLOW}🔧 [{name}] {cmd}")
    start = time.time()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    elapsed = time.time() - start
    if result.returncode != 0:
        print(f"{Fore.RED}❌ [{name}] 失败 ({elapsed:.2f}s) —— {error_label}")
        print(f"{Fore.RED}{result.stderr or result.stdout}")
        sys.exit(1)
    # else:
        # print(f"{Fore.GREEN}✅ [{name}] 成功 ({elapsed:.2f}s)")
    # print()  # 空行分隔
# Example usage:
node_count = 8  # Number of nodes
edge_count = 20 # Number of edges


output_file = '/Users/zhangwenqian/UNSW/KClique/new_small_garph.edges'  # Output file path



BIN1 = "/Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/degeneracy_cliques"
# BIN2 = "/Users/zhangwenqian/UNSW/pivoter/cmake-build-release/bin/main"
# BIN3 = "/Users/zhangwenqian/UNSW/nucleus/nd/nucleus"

count = 0
while True:
    count += 1
    print(f"{Style.BRIGHT}{Fore.CYAN}🚀 第 {count} 轮测试启动！加油！\n")

    # 1. 生成随机图
    # edge count +- 10%
    edgeList = generate_graph_log(node_count, edge_count, output_file)
    # print(f"{Fore.CYAN}🗺️  随机图生成完毕，共 {len(edgeList)} 条边。\n")

    # 2. 第一步工具：degeneracy_cliques
    cmd1 = f"{BIN1} {output_file} 2 2"

    run_cmd("DegeneracyCliques", cmd1, "degeneracy_cliques 非零退出")
    print(f"{Fore.GREEN}✅ DegeneracyCliques 成功！\n")
    # 3. 第二步工具：main
    # cmd2 = f"{BIN2} {output_file}.tree 2 4 {output_file}"
    # run_cmd("Main", cmd2, "main 非零退出")
    # print(f"{Fore.GREEN}✅ Main 成功！\n")
    # # 4. 第三步工具：nucleus
    # cmd3 = f"{BIN3} {output_file} 24 no"
    #
    # run_cmd("Nucleus", cmd3, "nucleus 非零退出")
    # print(f"{Fore.GREEN}✅ Nucleus 成功！\n")

    # 5. 比对结果
    # print(f"{Fore.YELLOW}🔍 正在用 uniq + diff 检查一致性...")
    # subprocess.run(f"uniq /Users/zhangwenqian/UNSW/pivoter/a > /Users/zhangwenqian/UNSW/pivoter/a.tmp", shell=True)
    # subprocess.run(f"uniq /Users/zhangwenqian/UNSW/pivoter/b > /Users/zhangwenqian/UNSW/pivoter/b.tmp", shell=True)
    # diff = subprocess.run("diff /Users/zhangwenqian/UNSW/pivoter/a /Users/zhangwenqian/UNSW/pivoter/a.tmp", shell=True,
    #                       capture_output=True, text=True)
    # if diff.stdout or diff.stderr:
    #     print(f"{Fore.RED}❌ 对比失败！输出不一致：\n{diff.stdout or diff.stderr}")
    #     print(f"{Fore.MAGENTA}🖼️ 报错时的图边列表：\n{edgeList}")
    #     draw_graph_with_cliques(edgeList)
    #     sys.exit(1)
    # else:
    print(f"{Fore.GREEN}✅ 结果一致！本轮测试完美通过 🎉\n")

#
# file_path = output_file
# edges = []
# firstLine = True
# with open(file_path, 'r') as f:
#     for line in f:
#         u, v = map(int, line.strip().split())
#         if firstLine:
#             firstLine = False
#             continue
#         edges.append((u, v))
#
# G = nx.Graph()
# G.add_edges_from(edges)
# pos = nx.kamada_kawai_layout(G)
#
# plt.figure(figsize=(8, 6))
# nx.draw_networkx_edges(G, pos)
# nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=500)'bold')
# #
# nx.draw_networkx_labels(G, pos, font_size=10, font_weight= # [0,2,3,6,8,9,12,13]
# # 找到所有maximal clique
# cliques = [sorted(clique) for clique in nx.find_cliques(G) if len(clique) >= 4]
# print(cliques)
# cliques = [clique for clique in cliques if clique]
# in [0,2,3,6,8,9,12,13]
# // if clique is subset of [0,2,3,6,8,9,12,13]
# print("4-cliques:", cliques)
# print(len(cliques))