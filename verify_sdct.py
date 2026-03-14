#!/usr/bin/env python3
"""
验证所有 SDCT 版本的正确性
"""
import subprocess
import sys
import os
import tempfile

def run_test(binary_path, graph_file):
    """运行测试并返回结果"""
    try:
        result = subprocess.run(
            [binary_path, graph_file, "2", "3"],
            capture_output=True,
            text=True,
            timeout=300
        )
        return result.stdout + result.stderr
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    except Exception as e:
        return f"ERROR: {str(e)}"

def extract_clique_count(output):
    """从输出中提取 clique count"""
    for line in output.split('\n'):
        if 'TreeGraph Clique Count:' in line:
            # 下一行应该是数字
            continue
        if line.strip() and line[0].isdigit():
            try:
                return float(line.strip())
            except:
                pass
    return None

def main():
    if len(sys.argv) < 2:
        print("Usage: verify_sdct.py <graph_file>")
        sys.exit(1)
    
    graph_file = sys.argv[1]
    bin_dir = "/Users/zhangwenqian/UNSW/pivoter/build/bin"
    
    if not os.path.exists(graph_file):
        print(f"Error: Graph file not found: {graph_file}")
        sys.exit(1)
    
    print("=" * 60)
    print("SDCT Correctness Verification")
    print("=" * 60)
    print(f"Graph: {graph_file}\n")
    
    # 运行 degeneracy_cliques（使用 SDCT_Par）
    print("Running degeneracy_cliques (uses SDCT_Par)...")
    output = run_test(f"{bin_dir}/degeneracy_cliques", graph_file)
    
    # 查找 TreeGraph Clique Count
    lines = output.split('\n')
    for i, line in enumerate(lines):
        if 'TreeGraph Clique Count:' in line:
            print(f"Found clique count output at line {i}")
            # 打印接下来的几行
            for j in range(i, min(i+10, len(lines))):
                print(lines[j])
            break
    
    print("\n" + "=" * 60)
    print("Full output:")
    print("=" * 60)
    print(output[:2000])  # 打印前 2000 字符

if __name__ == "__main__":
    main()
