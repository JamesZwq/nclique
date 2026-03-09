#!/usr/bin/env python3
"""生成更大的测试图用于多线程性能测试"""

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
    return len(edges)

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
        
        attempts = 0
        while len(targets) < m and attempts < 1000:
            r = random.random() * total_degree
            cumsum = 0
            for node, deg in enumerate(degrees):
                cumsum += deg
                if cumsum >= r:
                    targets.add(node)
                    break
            attempts += 1
        
        for target in targets:
            edges.add((min(target, new_node), max(target, new_node)))
            degrees[target] += 1
        
        degrees.append(len(targets))
    
    with open(filename, 'w') as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")
    
    print(f"Generated BA graph: {n} nodes, {len(edges)} edges")
    print(f"Saved to: {filename}")
    return len(edges)

if __name__ == "__main__":
    random.seed(42)
    
    print("Generating larger test graphs for performance testing...")
    print()
    
    # 中等图 - 用于初步性能测试
    print("1. Medium graph (500 nodes, dense)...")
    generate_erdos_renyi(500, 0.15, "test_medium_500.edges")
    print()
    
    # 大图 - 用于性能测试
    print("2. Large graph (1000 nodes)...")
    generate_erdos_renyi(1000, 0.08, "test_large_1000.edges")
    print()
    
    # 超大图 - 用于压力测试
    print("3. Extra large graph (2000 nodes)...")
    generate_barabasi_albert(2000, 8, "test_xlarge_2000.edges")
    print()
    
    # 巨大图 - 用于最终测试
    print("4. Huge graph (3000 nodes)...")
    generate_barabasi_albert(3000, 10, "test_huge_3000.edges")
    print()
    
    print("All test graphs generated!")
    print("\nRecommended test parameters:")
    print("  - test_medium_500.edges: r=3, s=4")
    print("  - test_large_1000.edges: r=3, s=4")
    print("  - test_xlarge_2000.edges: r=3, s=5")
    print("  - test_huge_3000.edges: r=3, s=5")
