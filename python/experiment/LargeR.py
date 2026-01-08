import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.ticker import NullLocator
from matplotlib.lines import Line2D

# ================= 配置区域 =================
DATA_DIR = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/"
# 输出路径
OUTPUT_DIR = "/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/"
# 本地预览
PNG_DIR = "plots_row_bw_png/"

# 目标: 2个数据集 x 2个R值
TARGET_DATASETS = ["com-youtube", "soc-pokec-relationships"]
TARGET_R_VALUES = [5, 6]

# ================= 样式常量 (黑白 & 紧凑) =================
EPS_T = 1e-1
TIME_YMAX_SEC = 3600
HOLLOW_MARKERS = True
MARKER_EDGE_WIDTH = 1.5

# 纯黑白算法配置
ALGO_CONFIG = {
    'Ours': {
        'sources': ['CBS-noHi'],
        'label': 'CBS-noHi (Ours)',
        'color': 'black',  # 纯黑
        'fmt': 'o-',  # 实线 + 圆点
        'lw': 2.5, 'ms': 9  # 线条稍细一点点以适应黑白
    },
    'Baseline': {
        'sources': ['ARB', 'ARB-noHi'],
        'label': 'ARB-noHi (Baseline)',
        'color': 'black',  # 纯黑
        'fmt': 'x--',  # 虚线 + 叉号
        'lw': 2.5, 'ms': 9
    }
}

# 绘图全局设置 (字体保持大号，但图变矮了)
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 20,
    'axes.labelsize': 22,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'lines.linewidth': 2.5,
    'axes.linewidth': 1.5,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'mathtext.fontset': 'stix',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})


# ================= 数据加载 (保持不变) =================

def load_data(data_dir):
    dfs = []
    file_list = ["data1.csv", "data4.csv", "dataARBnoHi.csv", "data3.csv", "dataAlls.csv"]

    print(f"Loading files from {data_dir}...")
    for fname in file_list:
        path = os.path.join(data_dir, fname)
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                if 'dataset_name' in df.columns:
                    df['dataset_name'] = df['dataset_name'].apply(
                        lambda x: os.path.basename(str(x)).replace('.edges', '').replace('.adj', '').replace('.txt',
                                                                                                             ''))
                if fname == "dataAlls.csv" and 'dataset_name' in df.columns:
                    df = df[~df["dataset_name"].str.contains(r"\.p", na=False)]
                dfs.append(df)
            except:
                pass

    if not dfs: return pd.DataFrame()
    combined = pd.concat(dfs, ignore_index=True)

    if 'source' not in combined.columns: combined['source'] = 'CBS'
    if 'r' in combined.columns:
        mask_low_r = (combined['source'] == 'CBS') & (combined['r'].isin([1, 2]))
        combined.loc[mask_low_r, 'source'] = 'CBS-noHi'

    for col in ['r', 's', 'total_sec']:
        combined[col] = pd.to_numeric(combined[col], errors='coerce')

    if 'exit_status' in combined.columns:
        combined = combined[combined['exit_status'] == 0]
    # rename CBS to CND

    return combined


# ================= 辅助函数 =================

def _get_smart_ticks(labels, max_ticks=4):
    n = len(labels)
    if n <= max_ticks: return np.arange(n), labels
    indices = np.linspace(0, n - 1, num=max_ticks)
    indices = np.unique(np.round(indices).astype(int))
    return indices, [labels[i] for i in indices]


def _style_axes(ax):
    ax.tick_params(axis='both', which='major', direction='out', length=6, width=1.5)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    for spine in ['top', 'right']: ax.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']: ax.spines[spine].set_linewidth(1.5)


def save_plot(fig, filename_base):
    eps_path = os.path.join(OUTPUT_DIR, f"{filename_base}.eps")
    # 不使用透明度，确保纯黑白兼容性
    fig.savefig(eps_path, format='eps', bbox_inches='tight', pad_inches=0.02)

    png_path = os.path.join(PNG_DIR, f"{filename_base}.png")
    fig.savefig(png_path, format='png', dpi=300, bbox_inches='tight', pad_inches=0.02)
    print(f"Saved: {eps_path}")


def plot_standalone_legend():
    """生成纯黑白长条图例"""
    fig_leg = plt.figure(figsize=(10, 0.8))
    ax_leg = fig_leg.add_subplot(111)
    handles, labels = [], []

    for key in ['Ours', 'Baseline']:
        cfg = ALGO_CONFIG[key]
        marker = cfg['fmt'][0]
        ls = cfg['fmt'][1:]

        line = Line2D([0], [0], label=cfg['label'], color=cfg['color'],
                      linestyle=ls, linewidth=cfg['lw'], marker=marker, markersize=12,
                      markerfacecolor='none' if HOLLOW_MARKERS else 'auto',
                      markeredgecolor=cfg['color'], markeredgewidth=2.0)
        handles.append(line)

        # fig.legend(handles, [lab.replace("CBS", "CND") for lab in labels], loc="center",
                   # ncol=(ncol or min(4, len(labels))), frameon=True, framealpha=1.0, edgecolor="black")
        labels.append(cfg['label'].replace("CBS", "CND").replace("-noHi", ""))

    ax_leg.legend(handles, labels, loc='center', ncol=2, frameon=False, fontsize=20)
    ax_leg.axis('off')
    save_plot(fig_leg, "legend_row_bw")
    plt.close(fig_leg)


# ================= 主绘图逻辑 =================

def plot_specific_figures(df):
    plot_idx = 0

    for ds_name in TARGET_DATASETS:
        for r_val in TARGET_R_VALUES:

            df_r = df[(df['dataset_name'] == ds_name) & (df['r'] == r_val)]
            max_s = 20 if ds_name == 'com-youtube' else 35
            df_r = df_r[df_r["s"] <= max_s]
            df_r = df_r[df_r["r"] != df_r["s"]]

            s_vals = sorted(df_r["s"].unique())
            if not s_vals: continue

            # 【修改点】更矮的画布，压扁高度
            fig, ax = plt.subplots(figsize=(6.0, 3.5))
            has_data = False

            for algo_key in ['Ours', 'Baseline']:
                cfg = ALGO_CONFIG[algo_key]
                target_sources = cfg['sources']

                sub_df = df_r[df_r['source'].isin(target_sources)].copy()
                if sub_df.empty: continue

                agg_df = sub_df.groupby('s')['total_sec'].min().reset_index().sort_values('s')
                valid_data = agg_df[agg_df['s'].isin(s_vals)]
                if valid_data.empty: continue

                x_indices = [s_vals.index(s) for s in valid_data['s']]
                ys = valid_data['total_sec'].tolist()

                marker = cfg['fmt'][0]
                ls = cfg['fmt'][1:]

                ax.plot(x_indices, ys,
                        label=cfg['label'], color=cfg['color'], linestyle=ls,
                        linewidth=cfg['lw'], marker=marker, markersize=cfg['ms'],
                        markerfacecolor='none' if HOLLOW_MARKERS else 'auto',
                        markeredgecolor=cfg['color'],
                        markeredgewidth=1.8)
                has_data = True

            if not has_data:
                plt.close(fig)
                continue

            tick_positions, tick_labels = _get_smart_ticks(s_vals, max_ticks=4)
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels)
            ax.set_xlabel("$s$")

            ax.set_yscale('log')
            ax.set_ylim(bottom=0.1, top=3600)

            if plot_idx == 0:
                ax.set_ylabel("Time (s)")
            else:
                ax.set_ylabel("")

            _style_axes(ax)

            # 【修改点】彻底移除网格
            ax.grid(False)

            # 使用 _bw 后缀
            fname = f"{ds_name}_r{r_val}_final"
            save_plot(fig, fname)
            plt.close(fig)

            plot_idx += 1


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(PNG_DIR, exist_ok=True)

    df = load_data(DATA_DIR)
    if df.empty:
        print("No data loaded.")
        return

    df = df[~df['dataset_name'].str.contains("inc", case=False, na=False)]

    print("Generating 4 BW plots (r=5, 6)...")
    plot_specific_figures(df)
    plot_standalone_legend()
    print("\nDone! All files saved to:", OUTPUT_DIR)


if __name__ == "__main__":
    main()