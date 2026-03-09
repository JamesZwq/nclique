#!/usr/bin/env python3
"""生成测试图用于多线程性能测试"""

import random
import sys

def generate_erdos_renyi(n, p, filename):
    """生成Erdos-Renyi随机图"""
    edges = set()
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < p:
                edges.add((i, j))
    
    with open(filename, 'w') as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")
    
    print(f"Generated graph: {n} nodes, {len(edges)} edges")
    print(f"Saved to: {filename}")

def generate_barabasi_albert(n, m, filename):
    """生成Barabasi-Albert无标度网络"""
    edges = set()
    
    # 初始完全图
    for i in range(m):
        for j in range(i+1, m):
            edges.add((i, j))
    
    # 优先连接
    degrees = [m-1] * m
    
    for new_node in range(m, n):
        # 选择m个节点连接
        total_degree = sum(degrees)
        targets = set()
        
        while len(targets) < m:
            r = random.random() * total_degree
            cumsum = 0
            for node, deg in enumerate(degrees):
                cumsum += deg
                if cumsum >= r:
                    targets.add(node)
                    break
        
        for target in targets:
            edges.add((target, new_node))
            degrees[target] += 1
        
        degrees.append(len(targets))
    
    with open(filename, 'w') as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")
    
    print(f"Generated BA graph: {n} nodes, {len(edges)} edges")
    print(f"Saved to: {filename}")

if __name__ == "__main__":
    random.seed(42)
    
    # 生成不同大小的测试图
    print("Generating test graphs...")
    print()
    
    # 小图 - 快速测试
    generate_erdos_renyi(50, 0.3, "test_small.edges")
    print()
    
    # 中图 - 性能测试
    generate_erdos_renyi(200, 0.2, "test_medium.edges")
    print()
    
    # 大图 - 压力测试
    generate_barabasi_albert(500, 5, "test_large.edges")
    print()
    
    print("All test graphs generated!")
