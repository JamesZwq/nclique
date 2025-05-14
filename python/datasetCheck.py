def remap_and_check_edges(input_path, output_path=None):
    node_map = dict()   # 原点ID -> 新编号
    next_id = 0
    edge_set = set()

    has_self_loop = False
    has_duplicate = False

    new_edges = []

    with open(input_path, 'r') as f:
        first_line = f.readline().strip()
        print(f"跳过第一行: {first_line}")

        for line_num, line in enumerate(f, 2):
            u_str, v_str = line.strip().split()
            u, v = int(u_str), int(v_str)

            # 重新编号
            if u not in node_map:
                node_map[u] = next_id
                next_id += 1
            if v not in node_map:
                node_map[v] = next_id
                next_id += 1

            u_new, v_new = node_map[u], node_map[v]

            # 检查自环
            if u_new == v_new:
                print(f"自环发现于第 {line_num} 行: ({u}, {v}) -> ({u_new}, {v_new})")
                has_self_loop = True

            # 检查重复边（无向边）
            edge = (min(u_new, v_new), max(u_new, v_new))
            if edge in edge_set:
                print(f"重复边发现于第 {line_num} 行: ({u}, {v}) -> ({u_new}, {v_new})")
                has_duplicate = True
            else:
                edge_set.add(edge)
                new_edges.append(edge)

    print(f"\n重新编号节点总数: {len(node_map)}")
    if not has_self_loop and not has_duplicate:
        print("✅ 没有发现自环或重复边")

    if output_path:
        with open(output_path, 'w') as out:
            out.write(f"{len(node_map)} {len(new_edges)}\n")
            for u, v in new_edges:
                out.write(f"{u} {v}\n")
        print(f"✅ 新图数据写入至 {output_path}")

    return node_map  # 如果你还想要反向映射，可以改为返回两者


remap_and_check_edges(
    input_path="/Users/zhangwenqian/UNSW/pivoter/dblp_removed.edge",
    output_path="/Users/zhangwenqian/UNSW/pivoter/dblp_removed.edge"
)