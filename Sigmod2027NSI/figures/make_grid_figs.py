#!/usr/bin/env python3
# Exp-1 figures: build time and memory, ours vs serial CND, per cell.
# Data transcribed VERBATIM from the authoritative scout logs on tods2:
#   /data/wenqianz/grid_scout_batch1.out  (HEAD 9bfe8d6, 2026-07-22, serial CND)
#   /data/wenqianz/grid_scout_fem.out     (HEAD 320cf03, 2026-07-23, serial CND)
# Protocol: CND runs per cell (1800 s budget, 300 GB cap); ours runs one sweep per
# row and a cell is charged its marginal time; ours memory is the row's peak.
# Failure glyphs (teacher convention): red x = over time budget, red - = memory kill.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CELLS = ["3,4","3,5","3,6","4,5","4,6","4,7","5,6","5,7","5,8"]

# per graph: cnd_s (None+kind for fail), ours_marg_s, ours_row_mb (per row), cnd_mb
G = {
 "WI": dict(name="web-it-2004",
   cnd=[("x",None)]*3+[("m",None)]+[("x",None)]*2+[("x",None)]+[("x",None)]*2,
   cnd_s=[None]*9,
   ours=[0.73,0.05,0.03, 1.45,0.08,0.07, 2.94,0.12,0.12],
   cnd_mb=[None]*9, ours_mb=[512,512,512, 783,783,783, 1199,1199,1199]),
 "WG": dict(name="web-Google",
   cnd_s=[25.77,24.91,24.81, 59.72,57.56,56.66, 215.95,171.74,181.40],
   cnd=[(None,0)]*9,
   ours=[18.02,2.10,4.91, 35.74,1.31,1.10, 69.21,1.41,1.38],
   cnd_mb=[5191,5227,5265, 12513,12629,12722, 30736,31229,31504],
   ours_mb=[6937]*3+[11803]*3+[19153]*3),
 "DB": dict(name="com-dblp",
   cnd_s=[4.11,4.52,4.53, 43.61,43.31,43.81, 1075.23,1032.74,1045.16],
   cnd=[(None,0)]*9,
   ours=[1.40,0.07,0.05, 3.53,0.08,0.07, 10.76,0.18,0.18],
   cnd_mb=[996,1000,1008, 5701,5694,5703, 93860,93875,93902],
   ours_mb=[532]*3+[742]*3+[1776]*3),
 "AZ": dict(name="com-amazon",
   cnd_s=[2.58,2.43,2.12, 1.64,1.62,1.63, 1.20,1.20,1.15],
   cnd=[(None,0)]*9,
   ours=[0.47,0.05,0.03, 0.11,0.02,0.01, 0.01,0.01,0.01],  # 0.00 drawn at harness resolution
   cnd_mb=[519,520,528, 362,362,365, 242,243,240],
   ours_mb=[388]*3+[122]*3+[36]*3),
 "NS": dict(name="nasasrb",
   cnd_s=[60.28,52.60,50.40, 714.53,614.69,532.54, None,None,None],
   cnd=[(None,0)]*6+[("x",None)]*3,
   ours=[9.24,0.32,0.41, 22.27,0.37,0.37, 43.41,0.80,0.79],
   cnd_mb=[4450,4574,4606, 19256,19383,19443, None,None,None],
   ours_mb=[3448]*3+[7704]*3+[14520]*3),
 "PK1": dict(name="pkustk11",
   cnd_s=[142.97,140.65,139.04, None,None,None, None,None,None],
   cnd=[(None,0)]*3+[("x",None)]*3+[("x",None)]*3,
   ours=[0.72,0.09,0.05, 1.27,0.10,0.08, 1.97,0.13,0.12],
   cnd_mb=[17047,17921,18088, None,None,None, None,None,None],
   ours_mb=[809]*3+[1260]*3+[1788]*3),
 "PK3": dict(name="pkustk13",
   cnd_s=[325.91,366.66,376.47, None,None,None, None,None,None],
   cnd=[(None,0)]*3+[("x",None),("x",None),("x",None)]+[("m",None),("x",None),("x",None)],
   ours=[20.58,0.27,0.26, 64.68,0.70,0.69, 164.54,1.90,1.89],
   cnd_mb=[38786,46491,50483, None,None,None, None,None,None],
   ours_mb=[8237]*3+[23507]*3+[57243]*3),
 "PW": dict(name="pwtk",
   cnd_s=[290.36,283.73,271.89, None,None,None, None,None,None],
   cnd=[(None,0)]*3+[("x",None),("x",None),("x",None)]+[("m",None),("x",None),("x",None)],
   ours=[3.79,0.52,0.71, 7.98,0.41,0.33, 13.18,0.52,0.33],
   cnd_mb=[34406,36566,37125, None,None,None, None,None,None],
   ours_mb=[2441]*3+[4192]*3+[6516]*3),
}
# failure kinds per cell, aligned with the logs:
G["WI"]["cnd"] = [("x",None),("x",None),("x",None), ("m",None),("x",None),("x",None),
                  ("x",None),("x",None),("x",None)]
G["PK1"]["cnd"]= [(None,0)]*3+[("x",None)]*6
G["NS"]["cnd"] = [(None,0)]*6+[("x",None)]*3

ORDER = ["WI","WG","DB","AZ","NS","PK1","PK3","PW"]
BUDGET_S = 1800.0
BUDGET_MB = 300*1024.0

def panelize(metric):
    # figure* renders at \textwidth (~7.0 in); draw at that size so every glyph
    # keeps its true point size (skill rule: no element below 7 pt after rendering).
    fig, axes = plt.subplots(2, 4, figsize=(7.05, 3.0), sharey=False)
    x = range(9); w = 0.38
    for ax, key in zip(axes.flat, ORDER):
        g = G[key]
        for i in x:
            kind = g["cnd"][i][0]
            if metric == "time":
                cv, ov = g["cnd_s"][i], max(g["ours"][i], 0.01)
                top, ylab = BUDGET_S, "time (s)"
            else:
                cv, ov = g["cnd_mb"][i], g["ours_mb"][i]
                top, ylab = BUDGET_MB, "memory (MB)"
            if cv is not None:
                ax.bar(i - w/2, cv, w, color="white", edgecolor="black", linewidth=0.6, zorder=2)
            else:
                glyph = r"$\times$" if kind == "x" else r"$-$"
                ax.text(i - w/2, top*0.72, glyph, color="red", ha="center",
                        va="center", fontsize=9, fontweight="bold", zorder=4)
            ax.bar(i + w/2, ov, w, color="black", edgecolor="black", linewidth=0.4, zorder=3)
        ax.set_yscale("log")
        if metric == "time":
            ax.set_ylim(0.008, 4000); ax.axhline(BUDGET_S, lw=0.5, ls=":", color="gray")
        else:
            ax.set_ylim(20, 500000); ax.axhline(BUDGET_MB, lw=0.5, ls=":", color="gray")
        # compact 2-digit cell labels ("34" = cell (3,4)) fit unrotated at 7 pt
        ax.set_xticks(list(x)); ax.set_xticklabels([c.replace(",","") for c in CELLS], fontsize=7)
        ax.set_title(key, fontsize=8.5, fontweight="bold", pad=2)
        ax.tick_params(axis="y", labelsize=7); ax.spines[["top","right"]].set_visible(False)
    for ax in axes[:,0]:
        ax.set_ylabel("time (s)" if metric=="time" else "memory (MB)", fontsize=8)
    # shared legend strip
    import matplotlib.patches as mp
    handles=[mp.Patch(facecolor="white",edgecolor="black",label="CND (per cell)"),
             mp.Patch(facecolor="black",label="NSI build (marginal per cell)")]
    fig.legend(handles=handles, ncol=2, loc="upper center", fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout(rect=(0,0,1,0.94), h_pad=0.8, w_pad=0.6)
    out = f"fig_grid_{metric}.pdf"
    fig.savefig(out); print("wrote", out)

panelize("time")
panelize("memory")
