import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import os
import numpy as np

# ================= 配置区 =================
INPUT_FILES = {
    'CBS_Old': "data/stress_test_resultsCBS.txt",
    'SOTA': "data/stress_test_SOTA.txt",
    'CBS_Memory': "data/stress_test_resultsCBSMemory.txt"
}

OUTPUT_DIR = "plots_paper_clean"  # 新的输出文件夹
EPS_OUTPUT_DIR = "/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure"

if not os.path.exists(OUTPUT_DIR): os.makedirs(OUTPUT_DIR)
if not os.path.exists(EPS_OUTPUT_DIR):
    try:
        os.makedirs(EPS_OUTPUT_DIR)
    except:
        pass

# ===== 【核心修改】论文发表级绘图风格 =====
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],

    # 1. 字体极大化 (应对缩小排版)
    'font.size': 22,
    'axes.labelsize': 24,
    'axes.titlesize': 24,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 16, # Legend font can be slightly smaller

    # 2. 线条与边框加粗
    'lines.linewidth': 3.0,
    'lines.markersize': 10,
    'axes.linewidth': 1.5,  # 坐标轴线加粗
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,

    # 3. 纹理加粗
    'hatch.linewidth': 1.2,

    # 4. 去除网格 (User Request)
    'axes.grid': False,

    # 5. 画布尺寸 (4:3 比例，适合紧凑排版)
    'figure.figsize': (5, 3.8),

    # 6. 边距极小
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05
})

# ===== 样式映射 (高对比度黑白灰) =====
STYLE_MAP = {
    # --- 折线图 (Runtime & Memory) ---
    'CBS': {'color': 'black', 'linestyle': '-', 'marker': 'o', 'label': 'CND'},
    'SOTA': {'color': 'black', 'linestyle': '--', 'marker': 'x', 'label': 'SOTA'},

    # --- 柱状图 (Time Breakdown) ---
    'TreeBuild': {
        'color': '#E0E0E0', 'edgecolor': 'black', 'hatch': '', 'label': 'Index Const.'
    },
    'Indexing': {
        'color': '#808080', 'edgecolor': 'black', 'hatch': '//', 'label': 'Counting'
    },
    'Peeling': {
        'color': 'white', 'edgecolor': 'black', 'hatch': 'XX', 'label': 'Peeling'
    },

    # --- 柱状图 (Memory Breakdown) ---
    'Graph': {'color': '#404040', 'edgecolor': 'black', 'hatch': '', 'label': 'Graph'},
    'Tree': {'color': '#808080', 'edgecolor': 'black', 'hatch': '//', 'label': 'Tree'},
    'Auxiliary Index': {'color': '#C0C0C0', 'edgecolor': 'black', 'hatch': '..', 'label': 'Auxiliary Index'},
}


# ================= 数据解析 (保持逻辑不变) =================
def parse_logs(filepath, algo_name):
    data = []
    columns = ['Algorithm', 'dataset', 's', 'r', 'time', 'memory_kb', 'exit_status', 'index_size', 't_tree', 't_index',
               't_count', 't_peel', 'density', 'mem_gb', 'num_leaves', 'num_cliques']
    if not os.path.exists(filepath): return pd.DataFrame(columns=columns)
    with open(filepath, 'r') as f:
        content = f.read()
    runs = content.split("================== RUN ==================")
    for run in runs:
        if not run.strip(): continue
        entry = {'Algorithm': algo_name}
        t_peel_inner, t_peel_outer = 0.0, 0.0
        for line in run.strip().split('\n'):
            line = line.strip()
            if line.startswith("dataset"):
                try:
                    entry['dataset'] = line.split(":", 1)[1].strip()
                except:
                    pass
            elif line.startswith("s ") or line.startswith("s\t"):
                try:
                    entry['s'] = int(line.split(":")[1].strip())
                except:
                    pass
            elif line.startswith("r ") or line.startswith("r\t"):
                try:
                    entry['r'] = int(line.split(":")[1].strip())
                except:
                    pass
            elif "Exit status:" in line:
                try:
                    entry['exit_status'] = int(line.split(":")[1].strip())
                except:
                    pass
            elif "User time (seconds):" in line:
                try:
                    entry['time'] = float(line.split(":")[1].strip())
                except:
                    pass
            elif "Maximum resident set size (kbytes):" in line:
                try:
                    entry['memory_kb'] = int(line.split(":")[1].strip())
                except:
                    pass
            elif "Data Structure Size:" in line:
                try:
                    entry['index_size'] = int(line.split(":")[1].strip())
                except:
                    pass
            elif "Tree Build took:" in line:
                match = re.search(r"took: (\d+(\.\d*)?([eE][+-]?\d+)?)", line)
                if match: entry['t_tree'] = float(match.group(1)) / 1000.0
            elif "clique Index build took:" in line:
                match = re.search(r"took: (\d+(\.\d*)?([eE][+-]?\d+)?)", line)
                if match: entry['t_index'] = float(match.group(1)) / 1000.0
            elif "countingPerRClique took:" in line:
                match = re.search(r"took: (\d+(\.\d*)?([eE][+-]?\d+)?)", line)
                if match: entry['t_count'] = float(match.group(1)) / 1000.0
            elif line.startswith("time:") and "ms" in line:
                match = re.search(r"time:\s*(\d+(\.\d*)?([eE][+-]?\d+)?)", line)
                if match: t_peel_inner = float(match.group(1)) / 1000.0
            elif "NucleusCoreDecomposition took:" in line:
                match = re.search(r"took: (\d+(\.\d*)?([eE][+-]?\d+)?)", line)
                if match: t_peel_outer = float(match.group(1)) / 1000.0
            elif "Clique Index:" in line:
                match = re.search(r"Clique Index: (\d+) cliques", line)
                if match: entry['index_size'] = int(match.group(1))
            elif "nun Leaf:" in line:
                match = re.search(r"nun Leaf: (\d+)", line)
                if match: entry['num_leaves'] = int(match.group(1))
            elif "Clique count:" in line:
                match = re.search(r"Clique count: ([\d\.e\+]+)", line)
                if match: entry['num_cliques'] = float(match.group(1))

        entry['t_peel'] = t_peel_inner if t_peel_inner > 0 else t_peel_outer
        if 'dataset' in entry and 'inc' in entry['dataset']:
            p_match = re.search(r"_p(\d+\.\d+)", entry['dataset'])
            if p_match:
                entry['density'] = float(p_match.group(1))
                entry['mem_gb'] = entry.get('memory_kb', 0) / (1024 * 1024)
                data.append(entry)
    return pd.DataFrame(data) if data else pd.DataFrame(columns=columns)


def parse_memory_detailed_logs(filepath):
    static_cols = ['dataset', 'density', 'r', 's', 'Graph', 'Tree', 'Auxiliary Index']
    dynamic_cols = ['dataset', 'density', 'r', 's', 'minCore', 'heap_size', 'memory_kb', 'memory_mb']
    if not os.path.exists(filepath): return pd.DataFrame(columns=static_cols), pd.DataFrame(columns=dynamic_cols)
    data_static, data_dynamic = [], []
    with open(filepath, 'r') as f:
        content = f.read()
    runs = content.split("================== RUN ==================")
    for run in runs:
        if not run.strip(): continue
        dataset, r, s, density = "", 0, 0, 0.0
        mem_graph, mem_tree, mem_total = 0, 0, 0
        dynamic_series = []
        cur_min_core, cur_heap, exit_status = 0, 0, -1
        for line in run.strip().split('\n'):
            line = line.strip()
            if line.startswith("dataset"):
                try:
                    dataset = line.split(":", 1)[1].strip()
                except:
                    pass
                p_match = re.search(r"_p(\d+\.\d+)", dataset)
                if p_match: density = float(p_match.group(1))
            elif line.startswith("s ") or line.startswith("s\t"):
                try:
                    s = int(line.split(":")[1].strip())
                except:
                    pass
            elif line.startswith("r ") or line.startswith("r\t"):
                try:
                    r = int(line.split(":")[1].strip())
                except:
                    pass
            elif "Exit status:" in line:
                try:
                    exit_status = int(line.split(":")[1].strip())
                except:
                    pass
            if "[Memory-Linux] Graph Memory:" in line:
                match = re.search(r"(\d+)\s*kB", line)
                if match: mem_graph = int(match.group(1))
            elif "[Memory-Linux] Tree Memory:" in line:
                match = re.search(r"(\d+)\s*kB", line)
                if match: mem_tree = int(match.group(1))
            elif "[Memory-Linux] Other Index Memory:" in line:
                match = re.search(r"(\d+)\s*kB", line)
                if match: mem_total = int(match.group(1)) # This is the total memory
            if line.startswith("minCore:"):
                c_match = re.search(r"minCore:\s*(\d+(\.\d+)?)", line)
                h_match = re.search(r"heap size:\s*(\d+)", line)
                if c_match: cur_min_core = float(c_match.group(1))
                if h_match: cur_heap = int(h_match.group(1))
            elif "[Memory-Linux] while loop top:" in line:
                match = re.search(r"(\d+)\s*kB", line)
                if match:
                    dynamic_series.append(
                        {'minCore': cur_min_core, 'heap_size': cur_heap, 'memory_kb': int(match.group(1))})
        if dataset and exit_status == 0:
            # Correct calculation: Tree_exclusive = Tree_total - Graph, Aux = Total - Tree_total
            mem_tree_exclusive = mem_tree - mem_graph
            mem_aux = mem_total - mem_tree
            
            data_static.append({'dataset': dataset, 'density': density, 'r': r, 's': s, 
                                'Graph': mem_graph / 1024,
                                'Tree': mem_tree_exclusive / 1024, 
                                'Auxiliary Index': mem_aux / 1024})
            for entry in dynamic_series:
                entry.update(
                    {'dataset': dataset, 'density': density, 'r': r, 's': s, 'memory_mb': entry['memory_kb'] / 1024})
                data_dynamic.append(entry)
    return pd.DataFrame(data_static), pd.DataFrame(data_dynamic)


# ================= 绘图函数 (无网格，高对比) =================
def save_plot(filename):
    png_path = os.path.join(OUTPUT_DIR, filename)
    plt.savefig(png_path)
    print(f"Generated: {png_path}")
    if os.path.exists(EPS_OUTPUT_DIR):
        eps_path = os.path.join(EPS_OUTPUT_DIR, filename.replace('.png', '.eps'))
        plt.savefig(eps_path, format='eps')
    plt.close()


def format_xaxis(ax, x_vals):
    # 字体已经很大，旋转防止重叠
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    # 稀疏显示：如果标签多，只显示头尾和中间
    n = len(x_vals)
    if n > 6:
        step = max(2, n // 3)
        for i, label in enumerate(ax.xaxis.get_ticklabels()):
            if i % step != 0 and i != n - 1:
                label.set_visible(False)


def plot_comparison(df, r, s, metric, ylabel, filename, log_scale=False):
    if df.empty: return
    sub = df[(df['r'] == r) & (df['s'] == s)].copy()
    if 'exit_status' in sub.columns: sub = sub[sub['exit_status'] == 0]
    if sub.empty: return

    plt.figure()
    ax = plt.gca()
    # 去除顶部和右侧的框线，更加学术
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for algo in ['CBS', 'SOTA']:
        dat = sub[sub['Algorithm'] == algo].sort_values('density')
        if not dat.empty:
            st = STYLE_MAP[algo]
            plt.plot(dat['density'], dat[metric], color=st['color'], ls=st['linestyle'], marker=st['marker'],
                     label=st['label'])

    # 移除轴标签文字（因为LaTeX标题会说明），只保留刻度
    # plt.ylabel(ylabel)
    if log_scale: plt.yscale('log')

    save_plot(filename)


def plot_breakdown(df, r, s):
    cbs = df[(df['Algorithm'] == 'CBS') & (df['r'] == r) & (df['s'] == s)].copy()
    if 'exit_status' in cbs.columns: cbs = cbs[cbs['exit_status'] == 0]
    cols = ['t_tree', 't_index', 't_count', 't_peel']
    for col in cols:
        if col not in cbs.columns: cbs[col] = 0.0
    cbs = cbs.fillna(0)
    cbs['total_time'] = cbs[cols].sum(axis=1)
    cbs = cbs[cbs['total_time'] > 0].sort_values('density')
    if cbs.empty: return

    densities = cbs['density'].apply(lambda x: f"{x:.2f}").tolist()
    x = np.arange(len(densities))
    width = 0.85

    plt.figure()
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    vals1 = (cbs['t_tree'] + cbs['t_index']).values
    st1 = STYLE_MAP['TreeBuild']
    plt.bar(x, vals1, width, color=st1['color'], edgecolor=st1['edgecolor'], hatch=st1['hatch'], linewidth=1.5)

    vals2 = cbs['t_count'].values
    st2 = STYLE_MAP['Indexing']
    plt.bar(x, vals2, width, bottom=vals1, color=st2['color'], edgecolor=st2['edgecolor'], hatch=st2['hatch'],
            linewidth=1.5)

    vals3 = cbs['t_peel'].values
    st3 = STYLE_MAP['Peeling']
    plt.bar(x, vals3, width, bottom=vals1 + vals2, color=st3['color'], edgecolor=st3['edgecolor'], hatch=st3['hatch'],
            linewidth=1.5)

    plt.xticks(x, densities)
    format_xaxis(ax, densities)
    save_plot(f"breakdown_r{r}_s{s}.png")


def plot_memory_breakdown(df, r, s):
    if df.empty or 'r' not in df.columns: return
    sub = df[(df['r'] == r) & (df['s'] == s)].sort_values('density')
    if sub.empty: return

    densities = sub['density'].apply(lambda x: f"{x:.2f}").tolist()
    x = np.arange(len(densities))
    width = 0.85

    plt.figure()
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    bottom = np.zeros(len(sub))
    comps = ['Graph', 'Tree', 'Auxiliary Index']
    for comp in comps:
        if comp in sub.columns:
            vals = sub[comp].fillna(0).values
            st = STYLE_MAP[comp]
            plt.bar(x, vals, width, bottom=bottom, color=st['color'], edgecolor=st['edgecolor'], hatch=st['hatch'],
                    linewidth=1.5)
            bottom += vals

    plt.xticks(x, densities)
    format_xaxis(ax, densities)
    save_plot(f"memory_breakdown_r{r}_s{s}.png")


def plot_heap_evolution(df, r, s):
    if df.empty or 'r' not in df.columns: return
    sub = df[(df['r'] == r) & (df['s'] == s)].copy()
    if sub.empty: return

    plt.figure()
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    densities = sorted(sub['density'].unique())
    
    # Define styles for different density lines for grayscale
    line_colors = ['black', 'dimgray', 'gray', 'darkgray', 'silver']
    line_styles = ['-', '--', ':', '-.']

    for i, density in enumerate(densities):
        plot_data = sub[sub['density'] == density].sort_values('minCore', ascending=False)
        if not plot_data.empty:
            label = f'density={density:.2f}'
            # Cycle through styles to make lines distinguishable
            color = line_colors[i % len(line_colors)]
            linestyle = line_styles[i % len(line_styles)]
            
            # The Y-axis is 'memory_mb' which is the total memory usage
            plt.plot(plot_data['minCore'], plot_data['memory_mb'], 
                     color=color, ls=linestyle, label=label)

    # Only add legend if there are multiple lines to distinguish
    if len(densities) > 1:
        plt.legend()

    # X-axis is coreness k, reversed to show peeling from high k to low k
    ax.invert_xaxis()
    plt.xlabel("k (coreness)")
    plt.ylabel("Memory (MB)")

    save_plot(f"heap_evolution_r{r}_s{s}.png")


def plot_leaf_vs_cliques(df, r, s):
    if df.empty: return
    sub = df[(df['Algorithm'] == 'CBS') & (df['r'] == r) & (df['s'] == s)].copy()
    if 'exit_status' in sub.columns: sub = sub[sub['exit_status'] == 0]
    sub = sub.dropna(subset=['num_leaves', 'num_cliques'])
    if sub.empty: return

    sub = sub.sort_values('density')
    
    plt.figure()
    ax1 = plt.gca()
    ax1.spines['top'].set_visible(False)
    
    # Plot cliques on the primary y-axis (ax1)
    line1, = ax1.plot(sub['density'], sub['num_cliques'], color='black', marker='o', linestyle='-', label='Clique Count')
    ax1.tick_params(axis='y')
    ax1.set_yscale('log')
    ax1.set_xlabel('Density') # Add x-label
    ax1.set_ylabel('Clique Count (log scale)') # Add y-label for ax1

    # Create a second y-axis for the number of leaves
    ax2 = ax1.twinx()
    ax2.spines['top'].set_visible(False)
    
    # Use a dashed line for the secondary metric
    line2, = ax2.plot(sub['density'], sub['num_leaves'], color='black', marker='x', linestyle='--', label='Leaf Count')
    ax2.tick_params(axis='y')
    ax2.set_ylabel('Leaf Count') # Add y-label for ax2

    # Combine legends from both axes
    lines = [line1, line2]
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper left', frameon=False) # Add legend
    
    save_plot(f"leaf_clique_r{r}_s{s}.png")


def export_legends():
    print("Generating Standalone Legends...")
    # Comparison
    fig_leg = plt.figure(figsize=(8, 0.8))
    handles = [
        mlines.Line2D([], [], color=STYLE_MAP['CBS']['color'], marker=STYLE_MAP['CBS']['marker'],
                      ls=STYLE_MAP['CBS']['linestyle'], label=STYLE_MAP['CBS']['label'], linewidth=3, markersize=12),
        mlines.Line2D([], [], color=STYLE_MAP['SOTA']['color'], marker=STYLE_MAP['SOTA']['marker'],
                      ls=STYLE_MAP['SOTA']['linestyle'], label=STYLE_MAP['SOTA']['label'], linewidth=3, markersize=12)
    ]
    # rename CBS to CND in legend, rename SOTA to ARB in legend
    for handle in handles:
        if handle.get_label() == 'CND':
            handle.set_label('CND')
        elif handle.get_label() == 'SOTA':
            handle.set_label('ARB')
    fig_leg.legend(handles=handles, loc='center', ncol=2, frameon=False, fontsize=24)

    plt.axis('off')
    save_plot("legend_comparison.png")

    # Time Breakdown
    fig_leg = plt.figure(figsize=(10, 0.8))
    handles = []
    for k in ['TreeBuild', 'Indexing', 'Peeling']:
        st = STYLE_MAP[k]
        handles.append(
            mpatches.Patch(facecolor=st['color'], edgecolor=st['edgecolor'], hatch=st['hatch'], label=st['label'],
                           linewidth=1.5))
    fig_leg.legend(handles=handles, loc='center', ncol=3, frameon=False, fontsize=24)
    plt.axis('off')
    save_plot("legend_time_breakdown.png")

    # Memory Breakdown
    fig_leg = plt.figure(figsize=(12, 0.8))
    handles = []
    order = ['Graph', 'Tree', 'Auxiliary Index']
    for k in order:
        st = STYLE_MAP[k]
        handles.append(
            mpatches.Patch(facecolor=st['color'], edgecolor=st['edgecolor'], hatch=st['hatch'], label=st['label'],
                           linewidth=1.5))
    fig_leg.legend(handles=handles, loc='center', ncol=3, frameon=False, fontsize=22)
    plt.axis('off')
    save_plot("legend_memory_breakdown.png")


def main():
    print("Parsing Logs...")
    df_cbs_old = parse_logs(INPUT_FILES['CBS_Old'], 'CBS')
    df_sota = parse_logs(INPUT_FILES['SOTA'], 'SOTA')
    df_mem_static, df_mem_dynamic = parse_memory_detailed_logs(INPUT_FILES['CBS_Memory'])
    df_cbs = pd.concat([df_cbs_old], ignore_index=True)
    if not df_cbs.empty: df_cbs = df_cbs.drop_duplicates(subset=['dataset', 'r', 's'], keep='first')
    df_all = pd.concat([df_cbs, df_sota], ignore_index=True)
    if df_all.empty: return

    export_legends()

    params = df_all[['r', 's']].drop_duplicates().sort_values(['r', 's'])
    for _, row in params.iterrows():
        r, s = int(row['r']), int(row['s'])
        print(f"Plotting r={r}, s={s}...")
        plot_comparison(df_all, r, s, 'time', None, f'runtime_r{r}_s{s}.png', log_scale=True)
        plot_comparison(df_all, r, s, 'mem_gb', None, f'memory_r{r}_s{s}.png', log_scale=False)
        plot_breakdown(df_cbs, r, s)
        plot_memory_breakdown(df_mem_static, r, s)
        plot_heap_evolution(df_mem_dynamic, r, s)
        plot_leaf_vs_cliques(df_cbs, r, s)
    print("Done.")


if __name__ == "__main__":
    main()
