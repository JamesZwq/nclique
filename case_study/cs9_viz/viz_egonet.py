"""
CS9: Case Study II analog — visualize (1,s)-core subgraphs as s varies.

Following Nuclear_CD_TODS Case Study II (Jiawei Han's ego network), we take
an "anchor" vertex (center of the densest region) and visualize its k-hop
neighborhood. For each s, color vertices by (1,s)-connected component.

Metrics reported per figure:
  - # surviving vertices
  - # connected components
  - size-weighted average separability m_in / m_cut
  - conductance m_cut / (2 m_in + m_cut)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams["pdf.fonttype"]=42  # TrueType, ACM/VLDB requirement
matplotlib.rcParams["ps.fonttype"]=42
import matplotlib.pyplot as plt
import networkx as nx
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).parent
GT_ROOT = ROOT.parent / "cs3_groundtruth"
CS6_ROOT = ROOT.parent / "cs6_grid"
N_SNAP = 425957

# ---- Load graph ----
print("Loading graph...")
adj = [set() for _ in range(N_SNAP)]
deg = np.zeros(N_SNAP, dtype=np.int32)
with open(GT_ROOT / "com-dblp.ungraph.txt") as f:
    for line in f:
        if line.startswith("#"): continue
        u, v = line.split(); u, v = int(u), int(v)
        if u != v:
            adj[u].add(v); adj[v].add(u)
for i in range(N_SNAP):
    deg[i] = len(adj[i])

# ---- Load cores ----
def load_mapping(path):
    df = pd.read_csv(path, sep="\t", comment="#", names=["new", "orig"], dtype=np.int64)
    m = np.full(int(df["new"].max()) + 1, -1, dtype=np.int64)
    m[df["new"].values] = df["orig"].values
    return m

def load_r1(path, mapping):
    df = pd.read_csv(path, sep="\t", comment="#", names=["vid", "core"],
                     dtype={"vid": np.int64, "core": np.float64})
    cv = np.zeros(N_SNAP)
    orig = mapping[df["vid"].values]
    mask = orig >= 0
    cv[orig[mask]] = df["core"].values[mask]
    return cv

print("Loading (1,s) cores for visualization s values...")
VIZ_S = [3, 5, 10]   # small s ("loose") → large s ("tight")
cores = {}
for s in VIZ_S:
    mp = load_mapping(CS6_ROOT / f"r1_s{s}_map.tsv")
    cores[s] = load_r1(CS6_ROOT / f"r1_s{s}.tsv", mp)

# ---- Find anchor: high-degree "hub" author whose 2-hop neighborhood is diverse ----
# Pick an author with DEGREE in top-100 but (1,10)-core NOT in the 113-clique
# (i.e., a hub with varied collaborations, not just one tight group)
candidates = np.where((deg > 50) & (deg < 300) & (cores[10] > 100) & (cores[10] < 1e8))[0]
# Prefer one whose 2-hop contains enough structural diversity
best_anchor = None
best_size = 0
for v in candidates[:500]:  # limit search
    hop2 = {int(v)}
    for u in adj[int(v)]:
        hop2.add(u)
        for w in adj[u]:
            hop2.add(w)
    if 400 <= len(hop2) <= 1500:  # target size range like the paper's case study
        if len(hop2) > best_size:
            best_anchor = int(v); best_size = len(hop2)
            if len(hop2) > 1000:
                break

if best_anchor is None:
    best_anchor = int(candidates[0]) if len(candidates) else int(np.argmax(cores[10]))
# Override: use the anchor identified by find_anchor.py (fragments cleanly at s=10)
anchor = 52213
print(f"\nAnchor vertex: {anchor} (degree={deg[anchor]}, (1,10)-core={cores[10][anchor]:.2e})")

# 2-hop ego net
hop1 = adj[anchor] | {anchor}
hop2 = set(hop1)
for u in hop1:
    hop2 |= adj[u]
egonet = hop2
print(f"2-hop ego net size: {len(egonet)}")
anchor_score = cores[10]

# ---- For each s, restrict cores to egonet and compute CCs ----
def components_induced(vertices, positive_score_vertices):
    """BFS-based CCs on induced subgraph of vertices ∩ positive_score_vertices."""
    survive = vertices & positive_score_vertices
    seen = set()
    comps = []
    for start in survive:
        if start in seen: continue
        comp = set()
        stack = [start]; seen.add(start); comp.add(start)
        while stack:
            v = stack.pop()
            for u in adj[v]:
                if u in survive and u not in seen:
                    seen.add(u); comp.add(u); stack.append(u)
        comps.append(comp)
    comps.sort(key=len, reverse=True)
    return comps

def measure(comps, egonet):
    """
    Per-component metrics:
      m_in  = internal edges of component
      m_cut = edges from component to NON-survivors (peeled vertices within egonet)
      sep   = m_in / m_cut
      cond  = m_cut / (2 m_in + m_cut)
    """
    total_v = sum(len(c) for c in comps)
    if not comps: return dict(n_comp=0, n_v=0, avg_sep=0, avg_cond=0, intra_frac=0)
    stats = []
    total_in = 0; total_cut = 0
    survive_all = set()
    for c in comps: survive_all |= c
    peeled = egonet - survive_all       # vertices that didn't survive
    for c in comps:
        m_in = 0; m_cut = 0
        for v in c:
            for u in adj[v]:
                if u == v: continue
                if u in c and u > v:
                    m_in += 1
                elif u in peeled:
                    m_cut += 1
        sep = m_in / m_cut if m_cut > 0 else (float("inf") if m_in > 0 else 0)
        cond = m_cut / (2 * m_in + m_cut) if (2 * m_in + m_cut) > 0 else 0
        stats.append({"size": len(c), "m_in": m_in, "m_cut": m_cut, "sep": sep, "cond": cond})
        total_in += m_in; total_cut += m_cut
    size_sum = sum(s["size"] for s in stats)
    sep_capped = [min(s["sep"], 100) for s in stats]
    wavg_sep = sum(s_val * st["size"] for s_val, st in zip(sep_capped, stats)) / size_sum
    wavg_cond = sum(s["cond"] * s["size"] for s in stats) / size_sum
    intra_frac = total_in / (total_in + total_cut) if (total_in + total_cut) > 0 else 0
    return dict(n_comp=len(comps), n_v=total_v, avg_sep=wavg_sep,
                avg_cond=wavg_cond, intra_frac=intra_frac, stats=stats)

# ---- Visualize 3 subplots: s=3, 5, 10 ----
print("\n=== Anchor-based ego network analysis ===")
results = {}
for s in VIZ_S:
    cv = cores[s]
    pos_vs = set(int(i) for i in np.where(cv > 0)[0])
    comps = components_induced(egonet, pos_vs)
    m = measure(comps, egonet)
    results[s] = {"comps": comps, "metrics": m}
    print(f"s={s:2d}: V={m['n_v']:5d}  CCs={m['n_comp']:3d}  "
          f"intra_frac={m['intra_frac']:.3f}  "
          f"avg_sep={m['avg_sep']:.3f}  avg_cond={m['avg_cond']:.3f}")
    # Top-5 components
    for i, st in enumerate(m["stats"][:5]):
        print(f"  CC{i+1}: size={st['size']:3d} m_in={st['m_in']:4d} m_cut={st['m_cut']:4d} sep={st['sep']:.2f}")

# ---- Build ONE shared layout from the largest egonet (s=3 survives all) ----
# All three panels reuse this layout, so peeling is visually a fade-out of
# vertices in place rather than three random force-graphs.
union_verts = set()
for s in VIZ_S:
    for c in results[s]["comps"]: union_verts |= c
union_verts.add(anchor)
G_all = nx.Graph()
for v in union_verts:
    G_all.add_node(v)
    for u in adj[v]:
        if u in union_verts and u != v:
            G_all.add_edge(u, v)
pos = nx.spring_layout(G_all, seed=42, iterations=120,
                      k=1.4 / max(np.sqrt(len(G_all)), 1))
# Pin anchor near the centre so each panel looks like the same ego-net.
ax_pos = pos[anchor]
shift  = np.array([0.0, 0.0]) - np.array(ax_pos)
for v in pos: pos[v] = (pos[v][0] + shift[0], pos[v][1] + shift[1])

# ---- Reference component coloring from s=10 (finest split) ----
ref_s = max(VIZ_S)
ref_comps = results[ref_s]["comps"]
K_cc      = min(8, len(ref_comps))
palette   = plt.cm.tab10(np.linspace(0, 1, K_cc))
comp_color = {}                  # vertex -> color, by its s=10 component
for i, c in enumerate(ref_comps[:K_cc]):
    for v in c: comp_color[v] = palette[i]

fig, axes = plt.subplots(1, len(VIZ_S), figsize=(4.2 * len(VIZ_S), 4.2))

for ax_i, s in enumerate(VIZ_S):
    ax = axes[ax_i]
    comps = results[s]["comps"]
    survive = set()
    for c in comps: survive |= c
    survive.add(anchor)

    # Edges: only those whose both endpoints survive at this s.
    edge_list = [(u, v) for u, v in G_all.edges()
                 if u in survive and v in survive]
    nx.draw_networkx_edges(G_all, pos, edgelist=edge_list, ax=ax,
                           alpha=0.30, width=0.5, edge_color="#888")

    # Peeled vertices: tiny gray dots (in place, so the reader sees the
    # neighborhood "shrinking" rather than a new graph each panel).
    peeled = union_verts - survive - {anchor}
    if peeled:
        px = [pos[v][0] for v in peeled]
        py = [pos[v][1] for v in peeled]
        ax.scatter(px, py, s=4, c="#dddddd", linewidths=0,
                   zorder=1)

    # Survivors coloured by s=10 reference component (consistent across
    # panels). Vertices not in any top-8 ref component use light gray.
    sx, sy, sc = [], [], []
    for v in survive:
        if v == anchor: continue
        sx.append(pos[v][0]); sy.append(pos[v][1])
        sc.append(comp_color.get(v, "#bbbbbb"))
    ax.scatter(sx, sy, s=22, c=sc, edgecolors="white",
               linewidths=0.3, alpha=0.95, zorder=2)

    # Anchor: always red star, labelled in the first panel only.
    ax.scatter([pos[anchor][0]], [pos[anchor][1]], marker="*",
               s=220, c="#d62728", edgecolors="white",
               linewidths=0.8, zorder=3)
    if ax_i == 0:
        ax.annotate("anchor", pos[anchor],
                    xytext=(10, 10), textcoords="offset points",
                    fontsize=9, color="#d62728", weight="bold")

    ax.set_title(f"$s={s}$", fontsize=11, fontweight="bold")
    ax.axis("off")
    ax.set_aspect("equal")

fig.tight_layout()
fig.savefig(ROOT / "cs9_egonet.png", dpi=150, bbox_inches="tight")
fig.savefig(ROOT / "cs9_egonet.pdf", bbox_inches="tight")
print(f"\nFigures saved to {ROOT}")
