"""
CS9 (Nuclear_CD-style ego-net case study):
Two-panel comparison of a com-dblp ego-network under a loose vs a tight
(1,s)-nucleus.  Reader sees:

  (a) small s -> one diffuse component, low edge density.
  (b) large s -> several NEAR-CLIQUE modules, high edge density.

Each panel annotates its top-5 components with labelled callout boxes
showing (N members, density), in the style of Nuclear_CD's Jiawei-Han
figure (which used (% Active) extracted from DBLP paper data).
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

# Two-panel layout: s_LO vs s_HI.
S_LO = 3
S_HI = 10
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
for s in (S_LO, S_HI):
    mp = load_mapping(CS6_ROOT / f"r1_s{s}_map.tsv")
    cores[s] = load_r1(CS6_ROOT / f"r1_s{s}.tsv", mp)

# ---- Build ego net ----
anchor = ANCHOR
hop1 = adj[anchor] | {anchor}
hop2 = set(hop1)
for u in hop1:
    hop2 |= adj[u]
egonet = hop2
print(f"Anchor={anchor} deg={deg[anchor]} ego2h={len(egonet)}")


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
# For panel labels we rank by  size * density  (the "quality x scale" score)
# rather than by raw size, because at large s the largest component is the
# loose residual background; the *story* components are the small cliques.
panels = []
for s in (S_LO, S_HI):
    alive = set(int(i) for i in np.where(cores[s] > 0)[0])
    comps = components_induced(egonet, alive)
    info  = []
    for c in comps:
        d, m_in = density(c)
        info.append({"comp": c, "n": len(c), "m_in": m_in, "dens": d,
                     "tag":  shape_label(c, d, m_in),
                     "qscore": d * len(c) if len(c) >= 3 else 0})
    # rank for annotation: panel (a) shows the one big loose component, so
    # the "top by quality" reduces to that. Panel (b) reranks by qscore so
    # the dense cliques get the callouts.
    info_ranked = sorted(info, key=lambda x: -x["qscore"])
    top = info_ranked[:TOP_K]
    avgd = sum(it["dens"] for it in top) / max(len(top), 1)
    panels.append({"s": s, "comps": info, "comps_ranked": info_ranked,
                   "alive": alive, "avgd_top": avgd})
    print(f"s={s}: {len(comps)} components, top-{TOP_K} avg density = {avgd:.2f}")
    for it in top:
        print(f"    n={it['n']:3d}  d={it['dens']:.2f}  {it['tag']}")


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


# ---- Render ----
fig, axes = plt.subplots(1, 2, figsize=(13, 6.2))

PALETTE = plt.cm.tab10(np.linspace(0, 1, TOP_K))

def render_panel(ax, panel):
    s     = panel["s"]
    comps = panel["comps"]
    survive_verts = set()
    for it in comps: survive_verts |= it["comp"]
    survive_verts.add(anchor)

    # Edges among survivors.
    el = [(u, v) for u, v in G_all.edges()
          if u in survive_verts and v in survive_verts]
    nx.draw_networkx_edges(G_all, pos, edgelist=el, ax=ax,
                           alpha=0.25, width=0.5, edge_color="#888")

    # Peeled vertices: gray dots in place.
    peeled = union_verts - survive_verts - {anchor}
    if peeled:
        px = [pos[v][0] for v in peeled]; py = [pos[v][1] for v in peeled]
        ax.scatter(px, py, s=4, c="#dddddd", linewidths=0, zorder=1)

    # Survivors not in top-K-by-quality: light gray.
    topK = panel["comps_ranked"][:TOP_K]
    topK_verts = set()
    for it in topK: topK_verts |= it["comp"]
    other_surv = survive_verts - topK_verts - {anchor}
    if other_surv:
        ox = [pos[v][0] for v in other_surv]; oy = [pos[v][1] for v in other_surv]
        ax.scatter(ox, oy, s=14, c="#bbbbbb", linewidths=0.3,
                   edgecolors="white", alpha=0.85, zorder=2)

    # Top-K components: colored. Callouts placed on a ring around the panel
    # centre to prevent overlap, with a leader line back to each centroid.
    n_top   = len(topK)
    ring_r  = 1.05 * max(abs(x) for x, _ in pos.values()) + 0.05
    for i, it in enumerate(topK):
        c     = it["comp"]
        color = PALETTE[i]
        cx = [pos[v][0] for v in c]; cy = [pos[v][1] for v in c]
        ax.scatter(cx, cy, s=46, c=[color], edgecolors="white",
                   linewidths=0.5, alpha=0.95, zorder=3)
        bx = float(np.mean(cx)); by = float(np.mean(cy))
        # Place the box on a ring; angle distributes evenly so boxes don't pile.
        ang = 2 * np.pi * (i + 0.5) / n_top
        lx  = ring_r * np.cos(ang)
        ly  = ring_r * np.sin(ang)
        label = f"{it['tag']}\n({it['n']} members, $d{{=}}${it['dens']:.2f})"
        ax.annotate(label,
                    xy=(bx, by),
                    xytext=(lx, ly),
                    fontsize=9, ha="center", va="center",
                    color="#202020",
                    bbox=dict(boxstyle="round,pad=0.35",
                              fc="white", ec=color, lw=1.3, alpha=0.97),
                    arrowprops=dict(arrowstyle="-",
                                    color=color, lw=1.0, alpha=0.9,
                                    connectionstyle="arc3,rad=0.05"))

    # Anchor.
    ax.scatter([pos[anchor][0]], [pos[anchor][1]], marker="*",
               s=360, c="#d62728", edgecolors="white",
               linewidths=1.0, zorder=4)
    ax.annotate("anchor", pos[anchor],
                xytext=(8, 8), textcoords="offset points",
                fontsize=9.5, color="#d62728", weight="bold")

    quality_word = "loose"    if panel["avgd_top"] < 0.30 else \
                   "moderate" if panel["avgd_top"] < 0.60 else "dense"
    ax.set_title(
        f"$s={s}$: {len(comps)} components, top-{TOP_K} avg density "
        f"$= {panel['avgd_top']:.2f}$  ({quality_word})",
        fontsize=11.5, fontweight="bold")
    ax.axis("off")
    # Leave room for the ring of callout boxes.
    pad = ring_r * 0.35
    xs = [x for x, _ in pos.values()]; ys = [y for _, y in pos.values()]
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect("equal")


for ax, p in zip(axes, panels):
    render_panel(ax, p)

fig.tight_layout()
fig.savefig(ROOT / "cs9_egonet.png", dpi=180, bbox_inches="tight")
fig.savefig(ROOT / "cs9_egonet.pdf", bbox_inches="tight")
print(f"\nFigures saved to {ROOT}")
