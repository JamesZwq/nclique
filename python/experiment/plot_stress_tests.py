import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.ticker import ScalarFormatter

# ================= 配置 =================
INPUT_FILES = {
    'CBS_New': "data/stress_test_resultsCBSNew.txt",
    'CBS_Old': "data/stress_test_resultsCBS.txt",
    'SOTA': "data/stress_test_SOTA.txt"
}
OUTPUT_DIR = "plots_real_dataset"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# ===== 学术绘图风格 =====
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'legend.fontsize': 13,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'savefig.dpi': 300,
    'lines.linewidth': 2,
    'lines.markersize': 8,
    'figure.figsize': (6, 4.5)
})

COLOR_MAP = {
    'CBS': '#000000',  # Black
    'SOTA': '#D62728',  # Red
    'TreeBuild': '#1f77b4',  # Blue
    'Indexing': '#ff7f0e',  # Orange
    'Peeling': '#2ca02c',  # Green
}

MARKER_MAP = {'CBS': 'o', 'SOTA': 'x'}
LINE_MAP = {'CBS': '-', 'SOTA': '--'}


# ================= 解析逻辑 =================

def parse_logs(filepath, algo_name):
    data = []
    if not os.path.exists(filepath):
        print(f"[WARN] File not found: {filepath}")
        return pd.DataFrame()

    with open(filepath, 'r') as f:
        content = f.read()

    # 按运行块分割
    runs = content.split("================== RUN ==================")

    for run in runs:
        if not run.strip(): continue

        entry = {'Algorithm': algo_name}
        dataset_name = ""

        for line in run.strip().split('\n'):
            line = line.strip()

            # 基础信息
            if line.startswith("dataset"):
                dataset_name = line.split(":", 1)[1].strip()
                entry['dataset'] = dataset_name
            elif line.startswith("s ") or line.startswith("s\t"):
                entry['s'] = int(line.split(":")[1].strip())
            elif line.startswith("r ") or line.startswith("r\t"):
                entry['r'] = int(line.split(":")[1].strip())
            elif "User time (seconds):" in line:
                entry['time'] = float(line.split(":")[1].strip())
            elif "Maximum resident set size (kbytes):" in line:
                entry['memory_kb'] = int(line.split(":")[1].strip())
            elif "Exit status:" in line:
                entry['exit_status'] = int(line.split(":")[1].strip())

            # SOTA 特有
            elif "Data Structure Size:" in line:
                entry['index_size'] = int(line.split(":")[1].strip())

            # CBS 特有 (Breakdown)
            elif "Tree Build took:" in line:
                entry['t_tree'] = float(re.search(r"took: (\d+(\.\d+)?)", line).group(1)) / 1000.0
            elif "clique Index build took:" in line:
                entry['t_index'] = float(re.search(r"took: (\d+(\.\d+)?)", line).group(1)) / 1000.0
            elif "countingPerRClique took:" in line:
                entry['t_count'] = float(re.search(r"took: (\d+(\.\d+)?)", line).group(1)) / 1000.0
            elif "NucleusCoreDecomposition took:" in line:  # Peeling time
                entry['t_peel'] = float(re.search(r"took: (\d+(\.\d+)?)", line).group(1)) / 1000.0
            elif "Clique Index:" in line:
                match = re.search(r"Clique Index: (\d+) cliques", line)
                if match: entry['index_size'] = int(match.group(1))

        # 只处理 "real" 数据集
        if 'dataset' in entry and 'real' in entry['dataset']:
            # 提取密度 p
            p_match = re.search(r"_p(\d+\.\d+)", entry['dataset'])
            if p_match:
                entry['density'] = float(p_match.group(1))
                entry['mem_gb'] = entry.get('memory_kb', 0) / (1024 * 1024)
                data.append(entry)

    return pd.DataFrame(data)


# ================= 绘图逻辑 =================

def plot_comparison(df, r, s, metric, ylabel, title, filename, log_scale=False):
    """绘制 CBS vs SOTA 对比图"""
    plt.figure()
    ax = plt.gca()

    # 筛选数据
    sub = df[(df['r'] == r) & (df['s'] == s)].copy()

    has_data = False
    for algo in ['CBS', 'SOTA']:
        algo_data = sub[sub['Algorithm'] == algo].sort_values('density')
        if not algo_data.empty:
            # 过滤掉失败的运行 (exit_status != 0)
            algo_data = algo_data[algo_data['exit_status'] == 0]

            plt.plot(algo_data['density'], algo_data[metric],
                     label=algo, color=COLOR_MAP[algo],
                     marker=MARKER_MAP[algo], linestyle=LINE_MAP[algo])
            has_data = True

    if not has_data:
        plt.close()
        return

    plt.xlabel(r'Density ($p$)')
    plt.ylabel(ylabel)
    plt.title(title)

    if log_scale:
        plt.yscale('log')

    # 美化
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend(frameon=False)
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()

    save_path = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(save_path)
    print(f"Generated: {save_path}")
    plt.close()


def plot_breakdown(df, r, s):
    """绘制 CBS 时间分解堆叠图"""
    cbs_data = df[(df['Algorithm'] == 'CBS') & (df['r'] == r) & (df['s'] == s) & (df['exit_status'] == 0)].sort_values(
        'density')

    if cbs_data.empty: return

    # 检查是否有 breakdown 数据
    required_cols = ['t_tree', 't_index', 't_count', 't_peel']
    if not all(col in cbs_data.columns for col in required_cols):
        return

    # 填充 NaN 为 0
    cbs_data[required_cols] = cbs_data[required_cols].fillna(0)

    densities = cbs_data['density'].astype(str)
    x = np.arange(len(densities))
    width = 0.6

    plt.figure()
    ax = plt.gca()

    # 绘制堆叠柱状图
    p1 = plt.bar(x, cbs_data['t_tree'] + cbs_data['t_index'], width, label='Index Const.', color=COLOR_MAP['TreeBuild'])
    p2 = plt.bar(x, cbs_data['t_count'], width, bottom=cbs_data['t_tree'] + cbs_data['t_index'], label='Counting',
                 color=COLOR_MAP['Indexing'])
    # Counting 底部是 Index Const.
    bottom_peel = cbs_data['t_tree'] + cbs_data['t_index'] + cbs_data['t_count']
    p3 = plt.bar(x, cbs_data['t_peel'], width, bottom=bottom_peel, label='Peeling', color=COLOR_MAP['Peeling'])

    plt.xlabel(r'Density ($p$)')
    plt.ylabel('Time (seconds)')
    plt.title(f'Runtime Breakdown (CBS)\nReal Dataset ($r={r}, s={s}$)')
    plt.xticks(x, densities)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend(frameon=False)
    plt.tight_layout()

    save_path = os.path.join(OUTPUT_DIR, f"breakdown_r{r}_s{s}.png")
    plt.savefig(save_path)
    print(f"Generated: {save_path}")
    plt.close()


# ================= 主程序 =================

def main():
    print("Parsing Log Files...")
    df_cbs_new = parse_logs(INPUT_FILES['CBS_New'], 'CBS')
    df_cbs_old = parse_logs(INPUT_FILES['CBS_Old'], 'CBS')  # 用作补充
    df_sota = parse_logs(INPUT_FILES['SOTA'], 'SOTA')

    # 合并数据 (优先使用 New)
    df_cbs = pd.concat([df_cbs_new, df_cbs_old]).drop_duplicates(subset=['dataset', 'r', 's'], keep='first')
    df_all = pd.concat([df_cbs, df_sota], ignore_index=True)

    if df_all.empty:
        print("No 'real' dataset data found!")
        return

    # 获取所有参数组合
    params = df_all[['r', 's']].drop_duplicates().sort_values(['r', 's'])

    print(f"Found experiments for: {params.values.tolist()}")

    for _, row in params.iterrows():
        r, s = int(row['r']), int(row['s'])

        # 1. 运行时间对比 (Log Scale)
        plot_comparison(df_all, r, s, 'time', 'Runtime (s)',
                        f'Runtime Comparison ($r={r}, s={s}$)',
                        f'runtime_r{r}_s{s}.png', log_scale=True)

        # 2. 内存对比
        plot_comparison(df_all, r, s, 'mem_gb', 'Peak Memory (GB)',
                        f'Memory Usage ($r={r}, s={s}$)',
                        f'memory_r{r}_s{s}.png', log_scale=False)

        # 3. 索引大小对比 (Log Scale)
        plot_comparison(df_all, r, s, 'index_size', 'Index Elements (Count)',
                        f'Index Compression ($r={r}, s={s}$)',
                        f'index_size_r{r}_s{s}.png', log_scale=True)

        # 4. CBS 时间分解图
        plot_breakdown(df_cbs, r, s)

    print("\nAnalysis Complete. Check 'plots_real_dataset' folder.")


if __name__ == "__main__":
    main()