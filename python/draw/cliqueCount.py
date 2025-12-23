import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from pathlib import Path

# — Academic style settings —
rcParams['font.family']      = 'Times New Roman'
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm']      = 'Times New Roman'
rcParams['mathtext.it']      = 'Times New Roman Italic'
rcParams['mathtext.bf']      = 'Times New Roman Bold'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

fontsize = 28
BAR_WIDTH = 2
FIGSIZE = (6, 4)            # 固定物理尺寸（英寸）
AX_RECT = [0.25, 0.24, 0.74, 0.70]
# [left, bottom, width, height]  20% 左边距、24% 下边距，宽高相应缩一点

def read_clique_file(path: str):
    ks, counts = [], []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or not line[0].isdigit():
                continue
            k_str, v_str = line.split(",", 1)
            ks.append(int(k_str.strip()))
            counts.append(float(v_str.strip()))
    return np.array(ks), np.array(counts, dtype=float)

def plot_one(data_path: str, out_path: str, title: str = None, use_log=True):
    ks, counts = read_clique_file(data_path)

    # 固定图尺寸 + 固定坐标轴矩形，确保不同数据集图的“高度”和“坐标轴高度”完全一致
    fig = plt.figure(figsize=FIGSIZE)
    ax  = fig.add_axes(AX_RECT)

    ax.bar(ks, counts, width=BAR_WIDTH, facecolor='gray')
    if use_log:
        ax.set_yscale('log')

    ax.set_xlabel('Clique Size $(k)$', fontsize=fontsize, labelpad=3)
    ax.set_ylabel('#Cliques',          fontsize=fontsize, labelpad=3)
    # if title:
    #     ax.set_title(title, fontsize=fontsize*0.8)

    ax.tick_params(axis='x', rotation=0, labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)

    # 不用 tight/ bbox_inches='tight'，避免因标签长度不同而改变导出尺寸
    fig.savefig(out_path, dpi=300)
    plt.close(fig)

# ====== 分别导出两张，高度完全一致，y 轴各自独立 ======
plot_one(
    data_path='/Users/zhangwenqian/UNSW/pivoter/dblpClique.txt',
    out_path ='/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/clique_distribution_dblp.eps',
    title='dblp', use_log=True
)

plot_one(
    data_path='/Users/zhangwenqian/UNSW/pivoter/StanfordClique.txt',
    out_path ='/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/clique_distribution_stanford.eps',
    title='web-Stanford', use_log=True
)