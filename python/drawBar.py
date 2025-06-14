from decimal import Decimal, getcontext
import math
import sys

# 提高小数运算的精度
getcontext().prec = 50

def collision_probability(m):
    """
    计算在 128-bit 哈希空间中，m 个元素的碰撞概率
    近似公式：p ≈ 1 - exp(- m*(m-1)/(2 * 2^128))
    对于 m 非常小的情况，可以近似 p ≈ m*(m-1)/(2 * 2^128)
    """
    m = Decimal(m)
    two_128 = Decimal(2) ** 128
    pair_count = m * (m - 1) / 2
    x = pair_count / two_128  # 这是 -log(1-p) 的近似值
    # 对于极小 x，可直接用 x 近似 p
    p_approx = x
    # 计算精确值 1 - exp(-x)，对于 x 极小也近似为 x
    p_exact = Decimal(1) - Decimal(math.exp(-float(x)))
    return p_approx, p_exact

def main():
    if len(sys.argv) != 4:
        print("用法: python collision.py <n_max> <clique_count> <k>")
        print("  <n_max> : 顶点范围最大值")
        print("  <clique_count> : 实际 clique 数量")
        print("  <k> : 每个 clique 的大小 (用于计算可能的总组合，非必需)")
        sys.exit(1)

    # 解析命令行参数
    try:
        n_max = int(sys.argv[1])
        clique_count = int(sys.argv[2])
        k = int(sys.argv[3])
    except ValueError:
        print("请确保所有输入都是整数。")
        sys.exit(1)

    # 计算可能的总 clique 数 (C(n_max, k))，如果 k > n_max，则忽略
    try:
        total_possible = math.comb(n_max, k)
    except ValueError:
        total_possible = None

    # 计算碰撞概率
    p_approx, p_exact = collision_probability(clique_count)

    # 输出结果
    print(f"输入参数: n_max = {n_max}, clique_count = {clique_count}, k = {k}")
    if total_possible is not None:
        print(f"可能的总 clique 数 = C({n_max}, {k}) = {total_possible}")
    else:
        print("无法计算总组合数 (可能 k > n_max)。")

    print(f"\n使用 {clique_count} 个元素计算 128 位哈希碰撞概率:")
    print(f"  近似碰撞概率 p ≈ {p_approx}")
    print(f"  精确计算 1 - exp(-x) ≈ {p_exact}")

if __name__ == "__main__":
    main()