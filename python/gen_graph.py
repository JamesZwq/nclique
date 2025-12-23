import networkx as nx
import os
import csv


def generate_formatted_stress_test():
    # ================= 配置区域 =================
    N = 1000  # 固定节点数

    # 生成 10 个密度点：从 0.1 (10%) 到 1.0 (100%)
    densities = [round(0.1 * i, 1) for i in range(1, 11)]

    output_dir = "stress_test_datasets"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    info_file = os.path.join(output_dir, "dataset_info.csv")
    # ===========================================

    print(f"Start generating {len(densities)} graphs with N={N}...")

    with open(info_file, 'w', newline='') as csvfile:
        fieldnames = ['filename', 'nodes', 'edges', 'density_target']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for p in densities:
            # 文件名使用 .edges 后缀，显得更专业
            filename = f"er_n{N}_p{p:.1f}.edges"
            filepath = os.path.join(output_dir, filename)

            print(f"Generating p={p:.1f} ... ", end="")

            # 1. 生成随机图 (固定种子 seed=42 保证结果可复现)
            G = nx.fast_gnp_random_graph(N, p, seed=42, directed=False)

            # 获取节点数和边数
            n_nodes = G.number_of_nodes()
            n_edges = G.number_of_edges()

            # 2. 写入文件，严格按照您要求的格式
            with open(filepath, 'w') as f:
                # [关键修改] 第一行写入：节点数 边数
                f.write(f"{n_nodes} {n_edges}\n")

                # 后续行写入：u v
                for u, v in G.edges():
                    f.write(f"{u} {v}\n")

            # 3. 记录元数据
            writer.writerow({
                'filename': filename,
                'nodes': n_nodes,
                'edges': n_edges,
                'density_target': p
            })

            print(f"Done. Header written: '{n_nodes} {n_edges}'")

    print(f"\nAll datasets generated in: ./{output_dir}")


if __name__ == "__main__":
    try:
        import networkx

        generate_formatted_stress_test()
    except ImportError:
        print("Error: NetworkX is missing. Run: pip install networkx")