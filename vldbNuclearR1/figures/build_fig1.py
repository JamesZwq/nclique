"""
Figure 1: running example overview.
Panels: (a) graph G with colored blocks and core labels;
        (b) the 6 CPI leaves with V_h | V_p boxes;
        (c) per-vertex (1,3)-core values as a horizontal bar chart.
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType — VLDB/ACM forbid Type 3
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['text.usetex']  = False
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import networkx as nx
from pathlib import Path

OUT = Path(__file__).parent

# ------- Panel (a): graph -------
G = nx.Graph()
edges = [
    (1,2),(1,3),(1,4),(2,3),(2,4),(3,4),           # Block A
    (3,5),(3,6),(4,5),(4,6),(5,6),                  # Block B
    (6,7),(6,8),(7,8),                              # T1
    (2,9),(2,10),(9,10),                            # T2
]
G.add_edges_from(edges)
pos = {
    1: (-2.4, 1.2), 2: (-2.4, -1.2),
    3: (-0.8, 0.6), 4: (-0.8, -0.6),
    5: (0.8, 1.2), 6: (0.8, -1.2),
    7: (2.4, -0.6), 8: (2.4, -1.8),
    9: (-3.8, -1.8), 10: (-3.8, 0.0),
}
core = {1:3, 2:3, 3:3, 4:3, 5:3, 6:3, 7:1, 8:1, 9:1, 10:1}
# Monochrome palette to match the paper's figure-style rules: black for
# the "winner" (high-core) class, light gray for the "secondary" class.
core_color = {3: "black", 1: "#cccccc"}

# ------- Panel (b): CPI leaves -------
leaves = [
    ("L_1", [7], [8, 6], 1),
    ("L_2", [9], [2, 10], 1),
    ("L_3", [1], [2, 3, 4], 3),
    ("L_4", [2], [3, 4], 1),
    ("L_5", [3], [4, 5, 6], 3),
    ("L_6", [4], [5, 6], 1),
]

# ------- Figure layout: 3 panels side by side with width ratios -------
fig = plt.figure(figsize=(15, 4.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1.15, 0.75], wspace=0.22)

# -------- Panel (a) --------
ax_a = fig.add_subplot(gs[0, 0])
block_a = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
block_b = [(3,5),(3,6),(4,5),(4,6),(5,6)]
t1 = [(6,7),(6,8),(7,8)]
t2 = [(2,9),(2,10),(9,10)]
shared_edge = [(3,4)]
# Distinguish blocks by line style + width rather than color (paper rule:
# monochrome plus shape/style, not hue).
for es, style, lw in [(block_a, "solid",  1.6),
                       (block_b, "dashed", 1.6),
                       (t1,      "dotted", 1.4),
                       (t2,      "dotted", 1.4)]:
    nx.draw_networkx_edges(G, pos, edgelist=es, edge_color="#555555",
                            style=style, width=lw, alpha=0.9, ax=ax_a)
nx.draw_networkx_edges(G, pos, edgelist=shared_edge, edge_color="black",
                        width=2.6, ax=ax_a)

node_colors = [core_color[core[v]] for v in G.nodes()]
label_colors = {v: ("white" if core_color[core[v]] == "black" else "black")
                for v in G.nodes()}
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=720,
                       edgecolors="black", linewidths=1.1, ax=ax_a)
for v, (x, y) in pos.items():
    ax_a.text(x, y, f"$v_{{{v}}}$", ha="center", va="center",
              fontsize=10.5, color=label_colors[v])
ax_a.set_title(r"(a) Input graph $G$: $|V|{=}10,\ |E|{=}17,\ 10$ triangles at $s{=}3$",
               fontsize=11)
ax_a.legend(handles=[
    mpatches.Patch(facecolor="black",   edgecolor="black", label=r"$(1,3)$-core $= 3$"),
    mpatches.Patch(facecolor="#cccccc", edgecolor="black", label=r"$(1,3)$-core $= 1$"),
], loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2, fontsize=9, frameon=False)
ax_a.axis("off")

# -------- Panel (b): CPI leaves --------
ax_b = fig.add_subplot(gs[0, 1])
ax_b.set_xlim(0, 10); ax_b.set_ylim(0, 10)
ax_b.axis("off")
ax_b.set_title(r"(b) Clique Path Index: 6 leaves encoding all 10 triangles", fontsize=11)

# Draw each leaf as a rounded box
# Layout: 2 columns × 3 rows
layout = [
    (1.0, 7.5), (5.5, 7.5),
    (1.0, 4.5), (5.5, 4.5),
    (1.0, 1.5), (5.5, 1.5),
]
for (name, h, p, ncliq), (x, y) in zip(leaves, layout):
    w, hh = 4.0, 2.3
    box = FancyBboxPatch((x, y), w, hh, boxstyle="round,pad=0.08",
                         linewidth=1.2, edgecolor="#333", facecolor="#f5f5f5")
    ax_b.add_patch(box)
    # Name on top-left
    ax_b.text(x + 0.25, y + hh - 0.38, f"${name}$", fontsize=11, fontweight="bold")
    # V_h | V_p content
    h_str = ",".join(f"v_{{{v}}}" for v in h)
    p_str = ",".join(f"v_{{{v}}}" for v in p)
    ax_b.text(x + 0.25, y + hh - 1.05,
              f"$V_h{{=}}\\{{{h_str}\\}}$", fontsize=10)
    ax_b.text(x + 0.25, y + hh - 1.65,
              f"$V_p{{=}}\\{{{p_str}\\}}$", fontsize=10)
    # Cliques count (annotation in black, bold; no green accent)
    eta = 3 - len(h)
    ax_b.text(x + w - 0.25, y + 0.3,
              f"$\\binom{{{len(p)}}}{{{eta}}}{{=}}{ncliq}$ cliques",
              fontsize=9, ha="right", color="black", fontweight="bold")

# -------- Panel (c): core values as horizontal bar chart --------
ax_c = fig.add_subplot(gs[0, 2])
vs = list(range(1, 11))
cs = [core[v] for v in vs]
colors = [core_color[c] for c in cs]
bars = ax_c.barh([f"$v_{{{v}}}$" for v in vs], cs, color=colors,
                 edgecolor="black", linewidth=0.6)
ax_c.invert_yaxis()
ax_c.set_xlim(0, 3.5)
ax_c.set_xticks([0, 1, 2, 3])
ax_c.set_xlabel(r"(1,3)-core value", fontsize=10)
ax_c.set_title(r"(c) Peeling output", fontsize=11)
for bar, c in zip(bars, cs):
    ax_c.text(bar.get_width() + 0.08, bar.get_y() + bar.get_height() / 2,
              f"{c}", va="center", fontsize=9)
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
ax_c.grid(axis="x", alpha=0.2, linestyle=":")

fig.suptitle(r"Figure 1: Running example $G$, its CPI at $s{=}3$, and (1,3)-core values.",
             fontsize=11.5, y=1.02)
fig.tight_layout()
fig.savefig(OUT / "fig1_overview.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig1_overview.png", dpi=170, bbox_inches="tight")
print(f"Saved fig1_overview to {OUT}")
