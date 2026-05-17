"""com-friendster (1.806 B edges) bench figure: 2 panels (peel time, memory).

Each panel groups bars by s value with SPIN-star and CND side by side.
CND completes only at s=2 (its bar appears there); for s>=3 CND is
marked OOM with an x annotation at the failure RSS (524 GB).

Data is hardcoded from paper_data/friendster_billion/bench_billion_ref.csv
and the SPIN-star V3 logs in the same directory.
"""
from pathlib import Path

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).parent

# s, SPIN-star peel(s), SPIN-star mem(GB, 10^6 base), CND peel(s)/'oom', CND mem(GB)
# Source: paper_data/friendster_billion/bench_billion.csv (V3) and
# bench_billion_ref.csv (REF). peel_ms / 1000 -> seconds; time_max_rss_kB
# * 1000 / 1024 / 10^6 (kB-to-GB conversion: time-v outputs in kB=1000B by GNU).
# To avoid mixing GiB/GB conventions we use 10^6 GB throughout.
DATA = [
    ( 2,  344,  160, 1774.84, 400),
    ( 3,  562,  340, "oom",   523),
    ( 4,  676,  472, "oom",   524),
    ( 5,  624,  467, "oom",   524),
    ( 6,  651,  451, "oom",   524),
    (11,  167,  312, None,    None),
]

SPIN_COLOR = "#1f78b4"
CND_COLOR  = "#e31a1c"
OOM_COLOR  = "#888888"


def plot(out_pdf):
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 2.6))

    xs = np.arange(len(DATA))
    width = 0.36

    spin_peel = [r[1] for r in DATA]
    spin_mem  = [r[2] for r in DATA]
    cnd_peel  = [r[3] for r in DATA]
    cnd_mem   = [r[4] for r in DATA]
    s_labels  = [str(r[0]) for r in DATA]

    # ---- panel 1: peel time ----
    ax = axes[0]
    ax.bar(xs - width/2, spin_peel, width, color=SPIN_COLOR,
           edgecolor="white", linewidth=0.4, label=r"SPIN$^{\star}$")
    for i, v in enumerate(cnd_peel):
        if isinstance(v, (int, float)):
            ax.bar(xs[i] + width/2, v, width, color=CND_COLOR,
                   edgecolor="white", linewidth=0.4,
                   label="CND" if i == 0 else None)
        elif v == "oom":
            ax.text(xs[i] + width/2, 80, "OOM", ha="center", va="bottom",
                    fontsize=8, color=OOM_COLOR, rotation=90)
    ax.set_yscale("log")
    ax.set_xticks(xs); ax.set_xticklabels(s_labels)
    ax.set_xlabel(r"$\boldsymbol{s}$", fontsize=10, fontweight="bold")
    ax.set_ylabel("peel time (s)", fontsize=10)
    ax.set_title("peel time", fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=8.5)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.grid(True, axis="y", which="major", alpha=0.25, linestyle=":")
    ax.legend(fontsize=9, frameon=False, loc="upper right")

    # ---- panel 2: memory ----
    ax = axes[1]
    ax.bar(xs - width/2, spin_mem, width, color=SPIN_COLOR,
           edgecolor="white", linewidth=0.4, label=r"SPIN$^{\star}$")
    for i, v in enumerate(cnd_mem):
        if v is None: continue
        if cnd_peel[i] == "oom":
            ax.bar(xs[i] + width/2, v, width, color=OOM_COLOR,
                   edgecolor="white", linewidth=0.4, alpha=0.7,
                   hatch="//")
            ax.text(xs[i] + width/2, v, "OOM", ha="center", va="bottom",
                    fontsize=7.5, color=OOM_COLOR)
        else:
            ax.bar(xs[i] + width/2, v, width, color=CND_COLOR,
                   edgecolor="white", linewidth=0.4,
                   label="CND" if i == 0 else None)
    # 503 GB budget line
    ax.axhline(503, color="black", linestyle=":", linewidth=0.9)
    ax.text(len(xs) - 0.55, 488, "503 GB budget",
            ha="right", va="top", fontsize=7.5, color="black",
            bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
    ax.set_xticks(xs); ax.set_xticklabels(s_labels)
    ax.set_xlabel(r"$\boldsymbol{s}$", fontsize=10, fontweight="bold")
    ax.set_ylabel("peak memory (GB)", fontsize=10)
    ax.set_title("peak memory", fontsize=10.5, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=8.5)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.grid(True, axis="y", which="major", alpha=0.25, linestyle=":")
    ax.set_ylim(0, 600)
    ax.legend(fontsize=9, frameon=False, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    plot(OUT / "fig_friendster.pdf")
