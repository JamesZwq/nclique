import networkx as nx

def read_graph_from_file(filename):
    """
    从文件中读取图信息，返回一个无向图对象
    """
    with open(filename, 'r', encoding='utf-8') as f:
        # 第一行包含顶点数 n 和边数 m
        n, m = map(int, f.readline().split())

        # 初始化无向图
        G = nx.Graph()
        # 先添加 n 个顶点（可选，后面直接加边也会自动创建顶点）
        G.add_nodes_from(range(n))

        # 读取 m 条边并加入到图中
        for _ in range(m):
            u, v = map(int, f.readline().split())
            G.add_edge(u, v)

    return G

def find_maximal_cliques(G):
    """
    使用 NetworkX 的 find_cliques 函数获取所有极大团
    返回结果是一个列表，其中每个元素是一个极大团的顶点集合
    """
    cliques = list(nx.find_cliques(G))
    #  sort
    for clique in cliques:
        clique.sort()
    # cliques sort by elements
    cliques.sort(key=lambda x: (len(x), x))
    return cliques

if __name__ == "__main__":
    # 文件名根据需要自行修改
    filename = "/Users/zhangwenqian/UNSW/KClique/newSmallGraph.txt"

    # 读取图
    graph = read_graph_from_file(filename)

    # 查找所有极大团
    maximal_cliques = find_maximal_cliques(graph)

    # 输出所有极大团
    print("所有极大团如下：")
    for clique in maximal_cliques:
        print(clique)