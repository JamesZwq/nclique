"""Hierarchy bench figure: 2 rows x 4 cols.

  Row 1 = wall time (ms, log y)
  Row 2 = peak memory (MB, log y)

Each panel: one graph, four lines (algos).  Reads paper_data/bench_r1_hierarchy.csv.
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

OUT      = Path(__file__).parent
DATA_CSV = Path("/Users/zhangwenqian/UNSW/pivoter/paper_data/bench_r1_hierarchy.csv")

GRAPHS = [
    ("com-youtube",              "com-youtube"),
    ("web-Google",               "web-Google"),
    ("soc-pokec-relationships",  "soc-pokec"),
    ("wiki-Talk",                "wiki-Talk"),
]

# (csv-key, paper-label, colour, marker)
ALGOS = [
    ("ST_V3",       r"SPIN$^{\star}$",                                 "#1f78b4", "o"),
    ("ST_V3_HIER",  r"SPIN$^{\star}$+BuildHier",                       "#33a02c", "s"),
    ("REF_R1",      r"CND",                                            "#e31a1c", "^"),
    ("REF_R1_HIER", r"CND+Hier",                                       "#ff7f00", "D"),
]


def load():
    rows = defaultdict(dict)  # (graph, s) -> {algo: row}
    with open(DATA_CSV) as f:
        for r in csv.DictReader(f):
            try: s = int(r["s"])
            except ValueError: continue
            rows[(r["graph"], s)][r["algo"]] = r
    return rows


def series(rows, graph, algo, metric):
    """Return (xs, ys) of valid OK cells for this algo on this graph.

    For metric == 'postbuild_ms' we report the algorithm's true post-CPI
    cost (excluding graph load and CPI build, which are common to all
    four algos and would otherwise dilute the comparison):

      ST_V3        -> peel_ms
      ST_V3_HIER   -> peel_ms + hier_ms
      REF_R1       -> total_ms (REF's "time: X ms" print of its own
                                function wall, post-build)
      REF_R1_HIER  -> total_ms (same; CND's union work is fused in)

    A few REF rows on wiki-Talk have a bogus total_ms (~1 ms while wall
    is several minutes) due to a stdout-flush race during the giant
    print_to_file dump; we fall back to wall_ms - build_ms there.
    """
    xs, ys = [], []
    for (g, s), cells in sorted(rows.items()):
        if g != graph: continue
        r = cells.get(algo)
        if not r or r.get("status") != "OK": continue
        try:
            if metric == "postbuild_ms":
                wall  = float(r.get("wall_ms",  -1) or -1)
                build = float(r.get("build_ms", -1) or -1)
                peel  = float(r.get("peel_ms",  -1) or -1)
                hier  = float(r.get("hier_ms",  -1) or -1)
                total = float(r.get("total_ms", -1) or -1)

                if algo.startswith("ST_V3"):
                    if peel < 0:
                        continue
                    y = peel + (hier if (algo.endswith("_HIER") and hier > 0) else 0.0)
                else:  # REF_R1 / REF_R1_HIER
                    if total >= 1.0 and (wall < 0 or total >= 0.005 * wall):
                        y = total
                    elif wall >= 0 and build >= 0:
                        y = max(wall - build, 0.0)
                    else:
                        continue
            else:
                y = float(r[metric])
        except (TypeError, ValueError): continue
        if y < 0: continue
        xs.append(s)
        ys.append(y)
    return xs, ys


def plot(out_pdf):
    rows = load()
    fig, axes = plt.subplots(2, len(GRAPHS), figsize=(14.5, 3.8), sharex=True)

    for col, (stem, lbl) in enumerate(GRAPHS):
        for row, (metric, ylabel) in enumerate([
            ("postbuild_ms", "peel + hier time (ms)"),
            ("mem_kB",       "memory (MB)"),
        ]):
            ax = axes[row][col]
            for algo, plbl, color, marker in ALGOS:
                xs, ys = series(rows, stem, algo, metric)
                if not xs: continue
                ys_plot = [y / 1024.0 for y in ys] if metric == "mem_kB" else ys
                ax.plot(xs, ys_plot, color=color, marker=marker,
                        markersize=4.5, linewidth=1.4, label=plbl)
            ax.set_yscale("log")
            ax.grid(True, which="major", alpha=0.25, linestyle=":")
            ax.tick_params(axis="both", which="major", labelsize=8.5)
            ax.tick_params(axis="both", which="minor", labelsize=0)
            for sp in ("top", "right"): ax.spines[sp].set_visible(False)
            if row == 0:
                ax.set_title(lbl, fontsize=10, fontweight="bold")
            if row == 1:
                ax.set_xlabel(r"$\boldsymbol{s}$", fontsize=10, fontweight="bold")
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=9.5)
            ax.set_xticks([3, 5, 7, 9, 11, 13, 15])

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=10,
               frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot(OUT / "fig_hier.pdf")
