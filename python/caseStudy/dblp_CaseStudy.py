import json
import os
import pandas as pd
from collections import defaultdict

# ================= 配置 =================
MAP_FILE = "/Users/zhangwenqian/UNSW/pivoter/python/caseStudy/han_ego_names.map"  # 你的名字映射文件
DATA_DIR = "/Users/zhangwenqian/UNSW/pivoter/python/caseStudy/vis_data/case_study_output"  # 数据文件夹
R_VAL = 3  # 固定 r=3
S_RANGE = range(4, 11)  # s 从 4 跑到 10


def load_name_map(map_file):
    id_to_name = {}
    if not os.path.exists(map_file):
        print(f"Warning: {map_file} not found. Showing IDs only.")
        return {}
    with open(map_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split(',', 1)
            if len(parts) == 2:
                id_to_name[int(parts[0])] = parts[1]
    return id_to_name


def analyze_batch():
    id_to_name = load_name_map(MAP_FILE)

    # 用于存储汇总数据的列表
    summary_data = []

    print(
        f"{'s':<5} | {'MaxCore':<10} | {'#Cliques':<10} | {'#Unique People (MaxCore)':<25} | {'Top Members in MaxCore'}")
    print("-" * 120)

    for s in S_RANGE:
        filename = f"metadata_r{R_VAL}_s{s}.json"
        filepath = os.path.join(DATA_DIR, filename)

        if not os.path.exists(filepath):
            print(f"Skipping {filename} (not found)")
            continue

        with open(filepath, 'r') as f:
            cliques = json.load(f)

        df = pd.DataFrame(cliques)

        # 1. 找到当前 s 下的最大 Core 值
        if df.empty:
            continue
        max_core_val = df['core'].max()

        # 2. 提取 Max-Core 的所有团
        max_core_cliques = df[df['core'] == max_core_val]
        num_cliques = len(max_core_cliques)

        # 3. 统计 Max-Core 涉及的唯一人数
        unique_people_ids = set()
        person_freq = defaultdict(int)

        for nodes in max_core_cliques['nodes']:
            for uid in nodes:
                unique_people_ids.add(uid)
                person_freq[uid] += 1

        num_people = len(unique_people_ids)

        # 4. 找出最核心的成员 (频率最高的前5名)
        sorted_people = sorted(person_freq.items(), key=lambda x: x[1], reverse=True)
        top_names = []
        for uid, count in sorted_people[:5]:
            name = id_to_name.get(uid, str(uid))
            # 简写名字以节省显示空间
            if len(name) > 15: name = name.split()[0] + "..."
            top_names.append(f"{name}({count})")

        top_members_str = ", ".join(top_names)

        print(f"{s:<5} | {max_core_val:<10.1f} | {num_cliques:<10} | {num_people:<25} | {top_members_str}")

        summary_data.append({
            "s": s,
            "max_core": max_core_val,
            "people_count": num_people,
            "top_members": [id_to_name.get(uid, str(uid)) for uid, _ in sorted_people[:10]]
        })

    return summary_data


if __name__ == "__main__":
    print(f"=== Batch Analysis for r={R_VAL}, s=[{min(S_RANGE)}..{max(S_RANGE)}] ===\n")
    data = analyze_batch()

    # 可以在这里加代码把 data 存成 CSV 方便画图