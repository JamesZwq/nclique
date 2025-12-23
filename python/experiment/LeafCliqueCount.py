# LeafCliqueCount.py — 学术风格版（黑/灰、无网格、同轴；单一通用图例）
import re
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

# ========= 路径 =========
IN_PATH = Path("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdataLeafCount")
OUT_DIR = Path("/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/")
CSV_DIR = Path("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data")  # CSV 辅助输出
assert IN_PATH.exists(), f"Input file not found: {IN_PATH}"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CSV_DIR.mkdir(parents=True, exist_ok=True)

# ========= 风格（学术/黑白友好/无透明）=========
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 15,
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": False,
    "savefig.transparent": False,  # EPS 友好
})

LINEWIDTH = 1.2
MARKERSIZE = 3.6

# ========= 可配置项 =========
USE_DUAL_Y = False              # 单轴（叶与 s-clique 共用同一 y 轴）
LOG_Y = True                    # 单轴下是否用 log（需要线性就改 False）
SAVE_PDF_ALSO = False           # 若需要同时导出 PDF
FILENAME_STYLE_SR = True        # True: leafCount_<dataset>_<s><r>.eps；False: leafCount_<dataset>_s<s>_r<r>.eps
SAVE_SINGLE_LEGEND = True       # 仅生成一张通用图例 legend_leafCount_generic.(eps|pdf)

# ========= 读入原始文本 =========
text = IN_PATH.read_text()

# 按 RUN block 切分
blocks = re.split(r"=+\s*RUN\s*=+", text)

def parse_header(blk: str):
    """解析 dataset/s/r（按行锚定，鲁棒些）"""
    m_ds = re.search(r"^dataset\s*:\s*(.+)$", blk, flags=re.MULTILINE)
    m_s  = re.search(r"^s\s*:\s*(\d+)\s*$",   blk, flags=re.MULTILINE)
    m_r  = re.search(r"^r\s*:\s*(\d+)\s*$",   blk, flags=re.MULTILINE)
    if not (m_ds and m_s and m_r):
        return None
    dataset_path = m_ds.group(1).strip()
    s = int(m_s.group(1))
    r = int(m_r.group(1))
    ds_name = Path(dataset_path).stem
    return ds_name, s, r

def parse_rows(blk: str):
    """
    提取每行：minCore/heap/num_leaf/<k>-Clique count
    注意：r=10/11 时是 '10-Clique count'/'11-Clique count'，因此用 \d+-Clique。
    """
    pat = r"minCore:\s*(\d+)\s+heap size:\s*(\d+)\s+num Leaf:\s*(\d+)\s+\d+-Clique count:\s*([\deE\+\.-]+)"
    rows = re.findall(pat, blk)
    if not rows:
        return None
    df = pd.DataFrame(rows, columns=["minCore", "heap_size", "num_leaf", "clique_count"]).astype(float)
    return df

def save_csv(df: pd.DataFrame, ds_name: str, s: int, r: int):
    csv_path = CSV_DIR / f"leafCount_{ds_name}_{s}{r}.csv"
    df.to_csv(csv_path, index=False)

def make_filename(ds_name: str, s: int, r: int, ext: str = "eps") -> Path:
    if FILENAME_STYLE_SR:
        # 保持 `_sr` 习惯（例如 s=3,r=10 -> _310）
        sr_str = f"{s}{r}"
        fname = f"leafCount_{ds_name}_{sr_str}.{ext}"
    else:
        fname = f"leafCount_{ds_name}_s{s}_r{r}.{ext}"
    return OUT_DIR / fname

def make_generic_legend_filename(ext: str = "eps") -> Path:
    return OUT_DIR / f"legend_leafCount_generic.{ext}"

def save_generic_legend():
    """
    生成单一通用图例：黑色(leaf count) + 灰色(s-clique)
    与主图风格完全一致（marker、线宽、颜色）。
    """
    handles = [
        Line2D([0], [0], color="black", marker="o", linewidth=LINEWIDTH, markersize=MARKERSIZE, label="#Paths"),
        Line2D([0], [0], color="0.5",   marker="s", linewidth=LINEWIDTH, markersize=MARKERSIZE, label="#s-clique"),
    ]
    labels = [h.get_label() for h in handles]
    fig_leg = plt.figure(figsize=(2.0, 0.35))
    fig_leg.legend(handles, labels, loc="center", ncols=len(labels), frameon=False)
    fig_leg.tight_layout(pad=0.1)
    out_eps = make_generic_legend_filename("eps")
    fig_leg.savefig(out_eps, format="eps", bbox_inches="tight", pad_inches=0.05)
    if SAVE_PDF_ALSO:
        out_pdf = make_generic_legend_filename("pdf")
        fig_leg.savefig(out_pdf, format="pdf", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig_leg)
    print(f"[Saved legend] {out_eps}")

def plot_run(df: pd.DataFrame, ds_name: str, s: int, r: int):
    x = df.index

    # 单 y 轴：黑/灰色、无网格，同一 y 轴
    fig, ax = plt.subplots(figsize=(3, 2))

    # 黑色：leaf count；灰色：s-clique
    ax.plot(x, df["num_leaf"], marker="o", linewidth=LINEWIDTH, markersize=MARKERSIZE,
            label="leaf count", color="black")
    ax.plot(x, df["clique_count"], marker="s", linewidth=LINEWIDTH, markersize=MARKERSIZE,
            label="s-clique", color="0.5")  # 这里图中不显示图例

    ax.set_xlabel("#Iterations")
    # ax.set_ylabel("count")

    # 单轴是否 log
    if LOG_Y and ((df["num_leaf"] > 0) | (df["clique_count"] > 0)).any():
        ax.set_yscale("log")

    # 背景不要格子；不显示图例
    ax.grid(False)

    # ax.set_title(f"{ds_name} (s={s}, r={r})")
    fig.tight_layout(pad=1.0)

    # 保存主图
    out_eps = make_filename(ds_name, s, r, "eps")
    fig.savefig(out_eps, format="eps")
    if SAVE_PDF_ALSO:
        out_pdf = make_filename(ds_name, s, r, "pdf")
        fig.savefig(out_pdf, format="pdf")
    plt.close(fig)

    return [str(out_eps)]

def main():
    runs = []
    for blk in blocks:
        header = parse_header(blk)
        if not header:
            continue
        ds_name, s, r = header

        df = parse_rows(blk)
        if df is None or df.empty:
            # 没有匹配到 minCore 行，跳过（避免画空图）
            continue

        # 导出 CSV（便于复现/调试）
        save_csv(df, ds_name, s, r)

        # 作图（主图中不显示图例）
        outs = plot_run(df, ds_name, s, r)
        runs.append({"dataset": ds_name, "s": s, "r": r, "rows": len(df), "outputs": outs})
        print(f"[Saved] {', '.join(outs)} ({len(df)} rows)")

    # 只在最后生成一张通用图例
    if SAVE_SINGLE_LEGEND:
        save_generic_legend()

    if not runs:
        print("No RUN blocks with parsed rows were found.")
    else:
        print("\nSummary:")
        for rinfo in runs:
            print(f" - {rinfo['dataset']} (s={rinfo['s']}, r={rinfo['r']}): rows={rinfo['rows']}")

if __name__ == "__main__":
    main()