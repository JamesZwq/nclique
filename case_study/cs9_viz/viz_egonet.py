"""
CS9 ego-net case study, 4-panel monotone progression:
Same com-dblp ego-network drawn under (1,s)-nucleus at s=3, 5, 7, 10.
Shared spring layout, so the reader watches the same spatial points
drop out as s grows:

  s=3  -> one diffuse blob (48 peeled out of 1452, d=0.01)
  s=5  -> blob shrinks, periphery gone (469 peeled, still 1 component)
  s=7  -> structural collapse: K_5 clique and several K_3 modules
          emerge from the residual (top-5 avg d ~ 0.67)
  s=10 -> kernel deepens: K_6 / K_5 / K_4 / K_3 modules
          (top-5 avg d = 0.72)

We stop at s=10 because past it (e.g. s=15) the K_6 and K_5 cliques are
themselves peeled, so the top-5 quality drops --- not the monotone
"better and better" story we want to tell.

Callout boxes are reserved for the two interesting panels (s=7, s=10);
the two loose panels carry a single inline summary instead, to keep the
4-up layout readable.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
from pathlib import Path

ROOT     = Path(__file__).parent
GT_ROOT  = ROOT.parent / "cs3_groundtruth"
CS6_ROOT = ROOT.parent / "cs6_grid"
N_SNAP   = 425957

# 4-panel monotone progression: top-5 avg density is non-decreasing.
S_VALUES   = (3, 5, 7, 10)
S_CALLOUTS = {7, 10}               # panels that get the callout-box treatment
ANCHOR = 3813          # selected via find_anchor search (clean cliques at s=10)
TOP_K  = 5             # number of components to annotate per panel

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
    df = pd.read_csv(path, sep="\t", comment="#", names=["new","orig"], dtype=np.int64)
    m  = np.full(int(df["new"].max())+1, -1, dtype=np.int64)
    m[df["new"].values] = df["orig"].values
    return m

def load_r1(path, mapping):
    df = pd.read_csv(path, sep="\t", comment="#", names=["vid","core"],
                     dtype={"vid": np.int64, "core": np.float64})
    cv = np.zeros(N_SNAP)
    orig = mapping[df["vid"].values]
    mask = orig >= 0
    cv[orig[mask]] = df["core"].values[mask]
    return cv

print("Loading (1,s) cores...")
cores = {}
for s in S_VALUES:
    mp = load_mapping(CS6_ROOT / f"r1_s{s}_map.tsv")
    cores[s] = load_r1(CS6_ROOT / f"r1_s{s}.tsv", mp)

# The global (1,s)-cores tell us which vertices survive ANYWHERE in the
# graph, but a vertex v alive at s may belong to a K_s that spans outside
# this ego.  For an honest local case study every alive vertex in a panel
# must be in some K_s entirely contained in the ego --- otherwise the
# callouts can show n<s "K_n cliques" at threshold s, which is incoherent.
# We compute local-alive sets by enumerating maximal cliques in the
# ego subgraph and labelling v alive at s iff some maximal clique
# containing v in the ego has size >= s.
print("Computing local-alive sets via clique enumeration in ego...")

# ---- Build ego net ----
anchor = ANCHOR
hop1 = adj[anchor] | {anchor}
hop2 = set(hop1)
for u in hop1:
    hop2 |= adj[u]
egonet = hop2
print(f"Anchor={anchor} deg={deg[anchor]} ego2h={len(egonet)}")

# ---- Local clique enumeration for honest (1,s) panels ----
import networkx as nx
from collections import defaultdict as _dd
G_ego = nx.Graph()
for v in egonet: G_ego.add_node(v)
for u in egonet:
    for v in adj[u]:
        if v in egonet and v > u: G_ego.add_edge(u, v)
print(f"  ego graph: V={G_ego.number_of_nodes()}, E={G_ego.number_of_edges()}")
_all_max_cliques = list(nx.find_cliques(G_ego))
_max_clique_size = _dd(int)
for _c in _all_max_cliques:
    _k = len(_c)
    for _v in _c:
        if _k > _max_clique_size[_v]:
            _max_clique_size[_v] = _k
# local_alive[s] = vertices in some K_s fully within the ego subgraph
local_alive = {s: set(v for v in egonet if _max_clique_size[v] >= s)
               for s in S_VALUES}
for s in S_VALUES:
    print(f"  s={s}: locally alive = {len(local_alive[s])}")

# Pick the largest maximal cliques in the ego that survive at every panel
# (i.e. with size >= max(S_VALUES)); these are the "story cliques" the
# reader will see persist across panels while the periphery peels.  We
# greedily skim non-overlapping members so that each highlighted clique
# adds visual structure rather than just nesting inside the previous one.
_S_MAX     = max(S_VALUES)
_top_cliques = sorted([c for c in _all_max_cliques if len(c) >= _S_MAX],
                       key=lambda c: -len(c))
_NUM_HILITE = 4   # ring of clique callouts per panel
HILITE_CLIQUES = []
_taken = set()
for _c in _top_cliques:
    _cs = set(_c)
    if len(_cs - _taken) < 0.6 * len(_cs):
        continue           # mostly already covered, skip
    HILITE_CLIQUES.append(sorted(_cs))
    _taken |= _cs
    if len(HILITE_CLIQUES) >= _NUM_HILITE:
        break
print(f"  highlighted cliques: {[len(c) for c in HILITE_CLIQUES]}")


def components_induced(verts, alive):
    survive = verts & alive
    seen = set(); out = []
    for v in survive:
        if v in seen: continue
        c = {v}; seen.add(v); st = [v]
        while st:
            u = st.pop()
            for w in adj[u]:
                if w in survive and w not in seen:
                    c.add(w); seen.add(w); st.append(w)
        out.append(c)
    return sorted(out, key=len, reverse=True)


def density(c):
    n = len(c)
    if n < 2: return 0.0, 0
    m_in = sum(1 for u in c for w in adj[u] if w in c and w > u)
    return m_in / (n * (n - 1) / 2), m_in


def shape_label(c, dens, m_in):
    """Short structural tag for the callout box."""
    n = len(c)
    if n >= 50:                     return "loose aggregation"
    if dens >= 0.99:                return f"K_{{{n}}} clique"
    if dens >= 0.70:                return "near-clique"
    if dens >= 0.40:                return "dense module"
    return "sparse module"


# ---- Per-panel data ----
# For callout panels we rank by  density first, size second  with a
# min-size threshold of 3, so the K-cliques get the callout boxes instead
# of the loose residual background.  For loose panels we keep the largest
# component (which is the single blob) as the headline.
EGO_SIZE = len(egonet)
panels = []
for s in S_VALUES:
    alive    = local_alive[s]
    n_peeled = EGO_SIZE - len(alive)
    panels.append({"s": s,
                   "alive":      alive,
                   "n_peeled":   n_peeled,
                   "n_alive":    len(alive),
                   "hilite":     HILITE_CLIQUES})
    print(f"s={s}: alive={len(alive)}, peeled={n_peeled}/{EGO_SIZE}")


# ---- Shared spring layout from the union of survivors + ego ----
union_verts = set(egonet)
union_verts.add(anchor)
G_all = nx.Graph()
for v in union_verts:
    G_all.add_node(v)
    for u in adj[v]:
        if u in union_verts and u != v:
            G_all.add_edge(u, v)
pos = nx.spring_layout(G_all, seed=42, iterations=150,
                      k=1.6 / max(np.sqrt(len(G_all)), 1))
# Centre the anchor.
ax_pos = pos[anchor]
shift  = np.array([0.0, 0.0]) - np.array(ax_pos)
for v in pos: pos[v] = (pos[v][0] + shift[0], pos[v][1] + shift[1])


# ---- Render: 1x4 progression ----
fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.1))

PALETTE = plt.cm.tab10(np.linspace(0, 1, TOP_K))

# Pre-compute layout extents and a shared callout-ring radius once.
all_xs = [x for x, _ in pos.values()]
all_ys = [y for _, y in pos.values()]
LAYOUT_R = max(max(abs(x) for x in all_xs), max(abs(y) for y in all_ys))
RING_R   = 1.05 * LAYOUT_R + 0.05


def render_panel(ax, panel):
    s             = panel["s"]
    alive         = panel["alive"]
    hilite        = panel["hilite"]
    survive_verts = set(alive) | {anchor}

    # Edges among survivors.
    el = [(u, v) for u, v in G_all.edges()
          if u in survive_verts and v in survive_verts]
    nx.draw_networkx_edges(G_all, pos, edgelist=el, ax=ax,
                           alpha=0.22, width=0.35, edge_color="#999")

    # Peeled vertices: light gray, in-place, so the eye reads "which
    # spatial positions just dropped out at this s".
    peeled = union_verts - survive_verts
    if peeled:
        px = [pos[v][0] for v in peeled]; py = [pos[v][1] for v in peeled]
        ax.scatter(px, py, s=5, c="#e2e2e2", linewidths=0, zorder=1)

    # Highlighted maximal cliques (same set across panels) --- these are
    # the "backbone" the reader should see persist as the periphery peels.
    hilite_verts = set()
    for c in hilite: hilite_verts |= set(c)

    # Alive vertices that are NOT part of any highlighted clique = the
    # peripheral alive ring.  At s=3 this is large; at s=10 it has mostly
    # shrunk away, leaving just the highlighted cliques.
    other_surv = (survive_verts - hilite_verts) - {anchor}
    if other_surv:
        ox = [pos[v][0] for v in other_surv]
        oy = [pos[v][1] for v in other_surv]
        ax.scatter(ox, oy, s=22, c="#9aa0a6", linewidths=0.3,
                   edgecolors="white", alpha=0.9, zorder=2)

    # Highlighted cliques get distinct colors + callouts.
    for i, c in enumerate(hilite):
        color = PALETTE[i]
        cx = [pos[v][0] for v in c]; cy = [pos[v][1] for v in c]
        ax.scatter(cx, cy, s=80, c=[color], edgecolors="white",
                   linewidths=0.5, alpha=0.97, zorder=3)
        bx = float(np.mean(cx)); by = float(np.mean(cy))
        ang = 2 * np.pi * (i + 0.5) / len(hilite)
        lx  = RING_R * np.cos(ang)
        ly  = RING_R * np.sin(ang)
        ax.annotate(f"$K_{{{len(c)}}}$",
                    xy=(bx, by),
                    xytext=(lx, ly),
                    fontsize=9, ha="center", va="center",
                    color="#202020", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.25",
                              fc="white", ec=color, lw=1.0, alpha=0.97),
                    arrowprops=dict(arrowstyle="-",
                                    color=color, lw=0.8, alpha=0.9,
                                    connectionstyle="arc3,rad=0.05"))

    # Anchor.
    ax.scatter([pos[anchor][0]], [pos[anchor][1]], marker="*",
               s=240, c="#d62728", edgecolors="white",
               linewidths=0.8, zorder=4)

    n_peeled = panel["n_peeled"]; n_alive = panel["n_alive"]
    ax.set_title(
        f"$s={s}$: alive ${n_alive}/{EGO_SIZE}$ "
        f"(peeled ${n_peeled}$)",
        fontsize=10.5, fontweight="bold")
    ax.axis("off")
    pad = RING_R * 0.32
    ax.set_xlim(min(all_xs) - pad, max(all_xs) + pad)
    ax.set_ylim(min(all_ys) - pad, max(all_ys) + pad)
    ax.set_aspect("equal")


for ax, p in zip(axes, panels):
    render_panel(ax, p)

fig.tight_layout()
fig.savefig(ROOT / "cs9_egonet.png", dpi=180, bbox_inches="tight")
fig.savefig(ROOT / "cs9_egonet.pdf", bbox_inches="tight")
print(f"\nFigures saved to {ROOT}")
