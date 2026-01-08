import xml.sax
import gzip
import os
import sys

# ================= 配置区 =================
INPUT_XML = "dblp.xml.gz"  # 请确保下载了此文件
OUTPUT_EDGE = "dblp_coauthor.edges"
OUTPUT_MAP = "dblp_names.map"

# 【过滤器】只保留 Data Mining / Database 顶级会议
# 这样跑出来的 Case Study 才有明显的社区结构 (DB圈 vs DM圈)
TARGET_CONFERENCES = {
    "kdd", "icdm", "sdm", "pkdd", "wsdm",  # Data Mining
    "sigmod", "vldb", "icde", "edbt", "pods",  # Database
    "neurips", "icml", "aaai", "ijcai"  # AI (可选, 数据量大可注释掉)
}


# ================= SAX 解析器 =================
class DBLPHandler(xml.sax.ContentHandler):
    def __init__(self):
        self.current_tag = ""
        self.authors = []
        self.year = ""
        self.key = ""
        self.is_target_conf = False

        # 映射: Name -> ID
        self.name_to_id = {}
        self.next_id = 0  # ID 从 0 开始

        # 使用 set 存储边，自动去重 (u, v) 且保证 u < v
        self.edges = set()
        self.paper_count = 0

    def startElement(self, tag, attributes):
        self.current_tag = tag
        if tag == "article" or tag == "inproceedings":
            self.authors = []
            self.key = attributes.get("key", "")
            # key 格式通常为: "conf/kdd/HanY19"
            parts = self.key.split('/')
            if len(parts) > 1:
                conf_name = parts[1].lower()
                # 检查是否在目标会议列表中
                if not TARGET_CONFERENCES or conf_name in TARGET_CONFERENCES:
                    self.is_target_conf = True
                else:
                    self.is_target_conf = False
            else:
                self.is_target_conf = False

    def characters(self, content):
        if self.is_target_conf and self.current_tag == "author":
            content = content.strip()
            if content:
                self.authors.append(content)

    def endElement(self, tag):
        if not self.is_target_conf: return

        if tag == "article" or tag == "inproceedings":
            # 如果作者数 > 1 且 < 10，构建全连接团 (忽略超大作者列表)
            if 1 < len(self.authors) < 20:
                # 1. 获取/分配 ID
                current_ids = []
                for name in self.authors:
                    if name not in self.name_to_id:
                        self.name_to_id[name] = self.next_id
                        self.next_id += 1
                    current_ids.append(self.name_to_id[name])

                # 2. 生成边 (u, v) 其中 u < v
                # 排序 ID 列表方便两两组合
                current_ids.sort()
                for i in range(len(current_ids)):
                    for j in range(i + 1, len(current_ids)):
                        u, v = current_ids[i], current_ids[j]
                        self.edges.add((u, v))

                self.paper_count += 1
                if self.paper_count % 10000 == 0:
                    print(
                        f"Parsed {self.paper_count} papers... (Nodes: {len(self.name_to_id)}, Unique Edges: {len(self.edges)})")


def main():
    if not os.path.exists(INPUT_XML):
        print(f"Error: {INPUT_XML} not found. Please download from https://dblp.org/xml/")
        return

    print("Parsing DBLP XML (Filter: {} conferences)...".format(len(TARGET_CONFERENCES)))

    handler = DBLPHandler()
    parser = xml.sax.make_parser()
    parser.setFeature(xml.sax.handler.feature_namespaces, 0)
    parser.setContentHandler(handler)

    try:
        with gzip.open(INPUT_XML, 'rt', encoding='utf-8') as f:
            parser.parse(f)
    except Exception as e:
        print(f"Parsing Error (check if file is corrupted): {e}")
        return

    num_nodes = len(handler.name_to_id)
    raw_edges = list(handler.edges)  # 转为列表准备排序
    num_edges = len(raw_edges)

    print(f"\nParsing Complete!")
    print(f"Nodes: {num_nodes}")
    print(f"Edges: {num_edges}")

    # ================= 关键步骤：排序与输出 =================

    print("Sorting edges...")
    # Python 的默认元组排序就是先按第一个元素，再按第二个元素排序
    # (0, 2) 会排在 (0, 5) 前面，(0, 5) 会排在 (1, 2) 前面
    raw_edges.sort()

    print(f"Writing sorted edge list to {OUTPUT_EDGE}...")
    with open(OUTPUT_EDGE, 'w') as f:
        # 【要求3】第一行写入 点数 边数
        f.write(f"{num_nodes} {num_edges}\n")

        # 写入排序后的边
        for u, v in raw_edges:
            f.write(f"{u} {v}\n")

    print(f"Writing name map to {OUTPUT_MAP}...")
    with open(OUTPUT_MAP, 'w', encoding='utf-8') as f:
        # ID,Name
        # 反转字典 ID -> Name
        id_to_name = {v: k for k, v in handler.name_to_id.items()}
        # 按 ID 顺序写入
        for i in range(num_nodes):
            f.write(f"{i},{id_to_name[i]}\n")

    print("All Done.")


if __name__ == "__main__":
    main()