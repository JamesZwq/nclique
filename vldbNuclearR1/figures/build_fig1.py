"""
Figure 1: running example — graph topology only.
The 6 CPI leaves and the (1,3)-core values are deliberately not drawn:
they live in Table tab:cpi-leaves and inline prose, respectively.
Per the paper figure-necessity rule, only the spatial/structural panel
earns figure space.
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType — VLDB/ACM forbid Type 3
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['text.usetex']  = False
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
from pathlib import Path

OUT = Path(__file__).parent

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
core_color = {3: "black", 1: "#cccccc"}

block_a = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
block_b = [(3,5),(3,6),(4,5),(4,6),(5,6)]
t1 = [(6,7),(6,8),(7,8)]
t2 = [(2,9),(2,10),(9,10)]
shared_edge = [(3,4)]

# Single panel: half-column width, tight layout.
fig, ax = plt.subplots(figsize=(3.4, 2.6))
for es, style, lw in [(block_a, "solid",  1.5),
                       (block_b, "dashed", 1.5),
                       (t1,      "dotted", 1.3),
                       (t2,      "dotted", 1.3)]:
    nx.draw_networkx_edges(G, pos, edgelist=es, edge_color="#555555",
                            style=style, width=lw, alpha=0.9, ax=ax)
nx.draw_networkx_edges(G, pos, edgelist=shared_edge, edge_color="black",
                        width=2.4, ax=ax)

node_colors = [core_color[core[v]] for v in G.nodes()]
label_colors = {v: ("white" if core_color[core[v]] == "black" else "black")
                for v in G.nodes()}
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=560,
                       edgecolors="black", linewidths=1.0, ax=ax)
for v, (x, y) in pos.items():
    ax.text(x, y, f"$v_{{{v}}}$", ha="center", va="center",
            fontsize=9, color=label_colors[v])

ax.legend(handles=[
    mpatches.Patch(facecolor="black",   edgecolor="black", label=r"$(1,3)$-core $= 3$"),
    mpatches.Patch(facecolor="#cccccc", edgecolor="black", label=r"$(1,3)$-core $= 1$"),
], loc="lower center", bbox_to_anchor=(0.5, -0.10), ncol=2, fontsize=8, frameon=False)
ax.axis("off")
ax.set_xlim(-4.4, 3.0); ax.set_ylim(-2.4, 2.0)

fig.savefig(OUT / "fig1_overview.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig1_overview.png", dpi=180, bbox_inches="tight")
print(f"Saved fig1_overview to {OUT}")
