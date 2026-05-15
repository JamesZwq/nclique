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
# For callout panels we rank by  density first, size second  with a
# min-size threshold of 3, so the K-cliques get the callout boxes instead
# of the loose residual background.  For loose panels we keep the largest
# component (which is the single blob) as the headline.
EGO_SIZE = len(egonet)
MIN_CALLOUT_N = 3
panels = []
for s in S_VALUES:
    alive = set(int(i) for i in np.where(cores[s] > 0)[0])
    comps = components_induced(egonet, alive)
    info  = []
    for c in comps:
        d, m_in = density(c)
        info.append({"comp": c, "n": len(c), "m_in": m_in, "dens": d,
                     "tag":  shape_label(c, d, m_in)})
    # Largest clique (perfect-density component) actually present at this s.
    max_clique = max((it["n"] for it in info
                      if it["dens"] >= 0.999 and it["n"] >= 2), default=1)
    # Vertices peeled out of the original ego-net.
    n_peeled = EGO_SIZE - sum(it["n"] for it in info)
    # Rank callouts by (density desc, size desc), restricted to dense
    # modules (d>=0.5, n>=3) so the loose residual blob never steals a
    # callout slot from a near-clique.  The residual is drawn in the
    # background gray; the callouts are reserved for "story" cliques.
    callout_pool = [it for it in info
                    if it["n"] >= MIN_CALLOUT_N and it["dens"] >= 0.5]
    callout_pool = sorted(callout_pool,
                          key=lambda x: (-x["dens"], -x["n"]))
    info_ranked  = callout_pool[:TOP_K]
    panels.append({"s": s, "comps": info, "comps_ranked": info_ranked,
                   "alive": alive,
                   "max_clique": max_clique,
                   "n_peeled":   n_peeled,
                   "n_alive":    EGO_SIZE - n_peeled})
    print(f"s={s}: {len(comps)} comps, peeled={n_peeled}, "
          f"max_clique=K_{max_clique}")
    for it in info_ranked:
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


# ---- Render: 1x4 progression ----
fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.1))

PALETTE = plt.cm.tab10(np.linspace(0, 1, TOP_K))

# Pre-compute layout extents and a shared callout-ring radius once.
all_xs = [x for x, _ in pos.values()]
all_ys = [y for _, y in pos.values()]
LAYOUT_R = max(max(abs(x) for x in all_xs), max(abs(y) for y in all_ys))
RING_R   = 1.05 * LAYOUT_R + 0.05


def render_panel(ax, panel):
    s        = panel["s"]
    comps    = panel["comps"]
    callouts = (s in S_CALLOUTS)
    survive_verts = set()
    for it in comps: survive_verts |= it["comp"]
    survive_verts.add(anchor)

    # Edges among survivors.
    el = [(u, v) for u, v in G_all.edges()
          if u in survive_verts and v in survive_verts]
    nx.draw_networkx_edges(G_all, pos, edgelist=el, ax=ax,
                           alpha=0.22, width=0.35, edge_color="#999")

    # Peeled vertices: light gray, in-place, so the eye reads "which
    # spatial positions just dropped out at this s".
    peeled = union_verts - survive_verts - {anchor}
    if peeled:
        px = [pos[v][0] for v in peeled]; py = [pos[v][1] for v in peeled]
        ax.scatter(px, py, s=5, c="#e2e2e2", linewidths=0, zorder=1)

    topK       = panel["comps_ranked"][:TOP_K]
    topK_verts = set()
    for it in topK: topK_verts |= it["comp"]

    # In the two loose panels (s=3, 5) every survivor goes the same blue:
    # there are no clique modules yet to extract, so the panel just shows
    # "still one big diffuse blob".  In the two callout panels the
    # survivors outside the highlighted cliques go medium gray (alive but
    # not a story component).
    if not callouts:
        loose_surv = survive_verts - {anchor}
        if loose_surv:
            ox = [pos[v][0] for v in loose_surv]
            oy = [pos[v][1] for v in loose_surv]
            ax.scatter(ox, oy, s=22, c="#4a90d9", edgecolors="white",
                       linewidths=0.3, alpha=0.85, zorder=3)
    else:
        other_surv = survive_verts - topK_verts - {anchor}
        if other_surv:
            ox = [pos[v][0] for v in other_surv]
            oy = [pos[v][1] for v in other_surv]
            ax.scatter(ox, oy, s=22, c="#9aa0a6", linewidths=0.3,
                       edgecolors="white", alpha=0.9, zorder=2)
        for i, it in enumerate(topK):
            c     = it["comp"]
            color = PALETTE[i]
            cx = [pos[v][0] for v in c]; cy = [pos[v][1] for v in c]
            ax.scatter(cx, cy, s=95, c=[color], edgecolors="white",
                       linewidths=0.6, alpha=0.97, zorder=3)
            bx = float(np.mean(cx)); by = float(np.mean(cy))
            ang = 2 * np.pi * (i + 0.5) / len(topK)
            lx  = RING_R * np.cos(ang)
            ly  = RING_R * np.sin(ang)
            label = f"{it['tag']}\n($n{{=}}{it['n']}$)"
            ax.annotate(label,
                        xy=(bx, by),
                        xytext=(lx, ly),
                        fontsize=7.5, ha="center", va="center",
                        color="#202020",
                        bbox=dict(boxstyle="round,pad=0.25",
                                  fc="white", ec=color, lw=1.0, alpha=0.97),
                        arrowprops=dict(arrowstyle="-",
                                        color=color, lw=0.8, alpha=0.9,
                                        connectionstyle="arc3,rad=0.05"))

    # Anchor.
    ax.scatter([pos[anchor][0]], [pos[anchor][1]], marker="*",
               s=240, c="#d62728", edgecolors="white",
               linewidths=0.8, zorder=4)

    # For the two loose panels, label the largest (single blob) component
    # in the corner so the reader doesn't need a ring of identical callouts.
    if not callouts:
        biggest = max(comps, key=lambda x: x["n"])
        ax.text(0.04, 0.96,
                f"loose aggregation\n($n{{=}}{biggest['n']}$)",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=8.5, color="#202020",
                bbox=dict(boxstyle="round,pad=0.3", fc="white",
                          ec="#4a90d9", lw=0.9, alpha=0.92))

    # Monotone headline numbers: peeled count (non-decreasing) and
    # largest emerging clique (non-decreasing).
    n_peeled = panel["n_peeled"]; n_alive = panel["n_alive"]
    kmax     = panel["max_clique"]
    clique_tag = (f"max clique $K_{{{kmax}}}$" if kmax >= 3
                  else "no $K_3$ yet")
    ax.set_title(
        f"$s={s}$: peeled ${n_peeled}/{EGO_SIZE}$, {clique_tag}",
        fontsize=10, fontweight="bold")
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
