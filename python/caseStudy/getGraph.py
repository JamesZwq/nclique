import xml.sax
import gzip
import os
import sys
from collections import defaultdict

# ================= 配置区 =================
INPUT_XML = "dblp.xml.gz"  # 请确保文件存在
OUTPUT_EDGE = "han_ego_2hop.edges"
OUTPUT_MAP = "han_ego_names.map"

# 目标中心人物 (DBLP中名字通常是唯一的，有时会有 ' 0001' 后缀，这里用包含匹配)
TARGET_NAME = "Jiawei Han"

# ================= 全局数据结构 (内存优化) =================
# 1. 名字 -> ID 映射 (String -> Int)
name_to_id = {}
id_to_name = []  # (Int -> String)
next_id = 0

# 2. 论文存储 (List of Lists of Ints)
#    all_papers[i] = [author_id_1, author_id_2, ...]
all_papers = []

# 3. 倒排索引: 作者ID -> 参与的论文索引列表
#    author_to_paper_indices[author_id] = [paper_idx_1, paper_idx_2, ...]
author_to_paper_indices = defaultdict(list)


# ================= SAX 解析器 =================
class DBLPLoader(xml.sax.ContentHandler):
    def __init__(self):
        self.current_tag = ""
        self.authors = []
        self.paper_count = 0
        self.in_paper = False

        # 引用全局变量
        global next_id

    def startElement(self, tag, attributes):
        self.current_tag = tag
        if tag == "article" or tag == "inproceedings":
            self.in_paper = True
            self.authors = []

    def characters(self, content):
        if self.in_paper and self.current_tag == "author":
            content = content.strip()
            if content:
                self.authors.append(content)

    def endElement(self, tag):
        if tag == "article" or tag == "inproceedings":
            self.in_paper = False

            # 过滤掉作者过多的论文 (比如大型综述或物理实验，通常不适合做Co-author分析)
            # Jiawei Han的圈子通常 2-20 人比较合理
            if 1 < len(self.authors) < 100:
                global next_id

                # 将名字转为 ID
                current_paper_author_ids = []
                for name in self.authors:
                    if name not in name_to_id:
                        name_to_id[name] = next_id
                        id_to_name.append(name)  # 保持 ID 和 index 一致
                        next_id += 1
                    current_paper_author_ids.append(name_to_id[name])

                # 存入全局数据
                paper_idx = len(all_papers)
                all_papers.append(current_paper_author_ids)

                # 更新倒排索引
                for uid in current_paper_author_ids:
                    author_to_paper_indices[uid].append(paper_idx)

                self.paper_count += 1
                if self.paper_count % 100000 == 0:
                    print(f"Loaded {self.paper_count} papers... (Total Authors: {len(name_to_id)})")


def main():
    if not os.path.exists(INPUT_XML):
        print(f"Error: {INPUT_XML} not found.")
        return

    # 1. 加载数据 (这一步最耗时，但为了构建完美的 Induced Subgraph 是必须的)
    print(f"Phase 1: Loading DBLP data into memory from {INPUT_XML}...")
    handler = DBLPLoader()
    parser = xml.sax.make_parser()
    parser.setFeature(xml.sax.handler.feature_namespaces, 0)
    parser.setContentHandler(handler)

    try:
        with gzip.open(INPUT_XML, 'rt', encoding='utf-8') as f:
            parser.parse(f)
    except Exception as e:
        print(f"Parsing Error: {e}")
        return

    print(f"Data Loaded. Total Papers: {len(all_papers)}, Total Authors: {len(name_to_id)}")

    # 2. 寻找中心节点 ID
    print(f"Phase 2: Searching for target '{TARGET_NAME}'...")
    target_id = -1
    # 简单的模糊匹配，防止 "Jiawei Han 0001" 这种情况
    for name, uid in name_to_id.items():
        if name == TARGET_NAME or (TARGET_NAME in name and len(name) < len(TARGET_NAME) + 5):
            print(f"Found Target: {name} (ID: {uid})")
            target_id = uid
            break

    if target_id == -1:
        print(f"Error: Target '{TARGET_NAME}' not found in dataset.")
        return

    # 3. BFS 查找 2-hop 节点 (Core Decomposition 的关键是节点集合)
    print("Phase 3: Performing BFS to find 2-hop neighbors...")

    # Set of IDs
    nodes_0hop = {target_id}
    nodes_1hop = set()
    nodes_2hop = set()

    # Find 1-hop
    print("  - Expanding 1-hop...")
    for p_idx in author_to_paper_indices[target_id]:
        for neighbor in all_papers[p_idx]:
            if neighbor != target_id:
                nodes_1hop.add(neighbor)

    print(f"    Found {len(nodes_1hop)} 1-hop neighbors.")

    # Find 2-hop
    print("  - Expanding 2-hop (this might take a moment)...")
    for u in nodes_1hop:
        # u 的所有 paper
        for p_idx in author_to_paper_indices[u]:
            for neighbor in all_papers[p_idx]:
                if neighbor not in nodes_0hop and neighbor not in nodes_1hop:
                    nodes_2hop.add(neighbor)

    print(f"    Found {len(nodes_2hop)} 2-hop neighbors.")

    # 合并所有有效节点
    valid_nodes = nodes_0hop | nodes_1hop | nodes_2hop
    print(f"Total Nodes in 2-hop Ego Network: {len(valid_nodes)}")

    # 4. 构建 Induced Subgraph (保留这些节点之间的所有连边)
    print("Phase 4: Constructing Induced Subgraph...")

    # 为了让输出 ID 连续 (0 ~ N-1)，我们需要重新映射 ID
    # Old_ID -> New_Continuous_ID
    old_to_new = {old_id: i for i, old_id in enumerate(valid_nodes)}

    final_edges = set()

    # 策略：遍历有效节点参与的所有 Paper。
    # 只有当一条边的两端都在 valid_nodes 里时，才保留这条边。
    # 为了效率，我们只遍历 valid_nodes 涉及到的 papers

    checked_papers = set()

    processed_count = 0
    total_valid_list = list(valid_nodes)

    # 这种遍历方式：遍历所有在子图中的人，看他们的论文
    # 如果论文里的其他人也在子图中，就加边
    for u_old in valid_nodes:
        processed_count += 1
        if processed_count % 10000 == 0:
            sys.stdout.write(f"\r  Processing node {processed_count}/{len(valid_nodes)}")
            sys.stdout.flush()

        for p_idx in author_to_paper_indices[u_old]:
            if p_idx in checked_papers:
                continue
            checked_papers.add(p_idx)

            # 检查这篇 paper 里的所有 clique 边
            paper_authors = all_papers[p_idx]

            # 只有当两个作者都在 valid_nodes 里时，才添加边
            # 这是一个 O(k^2) 操作，但 k 很小 (author list size)

            # 优化：先筛选出这篇 paper 里有多少个有效节点
            valid_authors_in_paper = [a for a in paper_authors if a in valid_nodes]

            if len(valid_authors_in_paper) < 2:
                continue

            # 构建团 (Clique)
            valid_authors_in_paper.sort()
            for i in range(len(valid_authors_in_paper)):
                for j in range(i + 1, len(valid_authors_in_paper)):
                    u = valid_authors_in_paper[i]
                    v = valid_authors_in_paper[j]

                    # 转为新 ID 存储
                    new_u = old_to_new[u]
                    new_v = old_to_new[v]

                    if new_u < new_v:
                        final_edges.add((new_u, new_v))
                    else:
                        final_edges.add((new_v, new_u))

    print(f"\nGraph Construction Complete. Edges: {len(final_edges)}")

    # 5. 输出文件
    print(f"Phase 5: Writing to {OUTPUT_EDGE}...")

    sorted_edges = sorted(list(final_edges))

    with open(OUTPUT_EDGE, 'w') as f:
        # Header: N M
        f.write(f"{len(valid_nodes)} {len(sorted_edges)}\n")
        for u, v in sorted_edges:
            f.write(f"{u} {v}\n")

    print(f"Writing name map to {OUTPUT_MAP}...")
    with open(OUTPUT_MAP, 'w', encoding='utf-8') as f:
        # New_ID, Old_ID, Name
        # 按照新 ID 排序输出
        # old_to_new 是 {old: new} -> 反转一下
        new_to_old = {v: k for k, v in old_to_new.items()}
        for i in range(len(valid_nodes)):
            original_id = new_to_old[i]
            original_name = id_to_name[original_id]
            f.write(f"{i},{original_name}\n")

    print("All Done. Use this dataset for your Core Decomposition!")


if __name__ == "__main__":
    main()