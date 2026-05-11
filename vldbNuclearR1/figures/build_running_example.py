"""Draw the running-example graph for the paper (10 vertices)."""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType — VLDB/ACM forbid Type 3
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['text.usetex']  = False
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path

OUT = Path(__file__).parent

G = nx.Graph()
edges = [
    # Block A: K4 on {1,2,3,4}
    (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4),
    # Block B: K4 on {3,4,5,6}
    (3, 5), (3, 6), (4, 5), (4, 6), (5, 6),
    # Pendant triangle {6, 7, 8}
    (6, 7), (6, 8), (7, 8),
    # Pendant triangle {2, 9, 10}
    (2, 9), (2, 10), (9, 10),
]
G.add_edges_from(edges)

# Custom positions so K4 blocks are visually obvious
pos = {
    # Block A on the left
    1: (-2.4, 1.2), 2: (-2.4, -1.2),
    3: (-0.8, 0.6), 4: (-0.8, -0.6),
    # Block B on the right of Block A
    5: (0.8, 1.2), 6: (0.8, -1.2),
    # Pendant triangle T1 far right
    7: (2.4, -0.6), 8: (2.4, -1.8),
    # Pendant triangle T2 far lower-left
    9: (-3.8, -1.8), 10: (-3.8, 0.0),
}

# Core values from peeling trace
core = {1: 3, 2: 3, 3: 3, 4: 3, 5: 3, 6: 3,
        7: 1, 8: 1, 9: 1, 10: 1}
color_by_core = {1: "#e8a44c", 3: "#4c8be8"}

fig, ax = plt.subplots(figsize=(8.5, 4.6))
# Draw edges by block for color
block_a = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
block_b = [(3,5),(3,6),(4,5),(4,6),(5,6)]
t1 = [(6,7),(6,8),(7,8)]
t2 = [(2,9),(2,10),(9,10)]
for es, col, lw in [(block_a, "#3a5fcd", 1.7),
                     (block_b, "#1f8a3a", 1.7),
                     (t1, "#b05020", 1.3),
                     (t2, "#b05020", 1.3)]:
    nx.draw_networkx_edges(G, pos, edgelist=es, edge_color=col, width=lw, alpha=0.85, ax=ax)

# Highlight shared edge (3,4) in bold
nx.draw_networkx_edges(G, pos, edgelist=[(3,4)], edge_color="#7a2bb0", width=3.2, ax=ax)

node_colors = [color_by_core[core[v]] for v in G.nodes()]
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=780,
                       edgecolors="black", linewidths=1.2, ax=ax)
nx.draw_networkx_labels(G, pos, labels={v: f"$v_{{{v}}}$" for v in G.nodes()},
                        font_size=11, font_color="black", ax=ax)

# Legend
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
ax.legend(handles=[
    Patch(facecolor="#4c8be8", edgecolor="black", label=r"(1,3)-core = 3"),
    Patch(facecolor="#e8a44c", edgecolor="black", label=r"(1,3)-core = 1"),
    Line2D([0],[0], color="#3a5fcd", lw=1.7, label=r"Block A: $K_4(v_1,\ldots,v_4)$"),
    Line2D([0],[0], color="#1f8a3a", lw=1.7, label=r"Block B: $K_4(v_3,\ldots,v_6)$"),
    Line2D([0],[0], color="#7a2bb0", lw=3.0, label=r"shared edge $(v_3,v_4)$"),
    Line2D([0],[0], color="#b05020", lw=1.3, label=r"pendant triangles $T_1,T_2$"),
], loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=8, frameon=False)

ax.set_title(
    r"Running example $G$: two $K_4$ blocks sharing edge $(v_3,v_4)$ + two pendant triangles; "
    r"10 triangles at $s{=}3$",
    fontsize=11
)
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT / "running_example_graph.png", dpi=150, bbox_inches="tight")
fig.savefig(OUT / "running_example_graph.pdf", bbox_inches="tight")
print(f"Saved to {OUT}")
