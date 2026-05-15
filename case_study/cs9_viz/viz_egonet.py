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


# ---- Render: treemap, one rectangle per component, area proportional
#     to member count, fill colour encoding internal edge density.
import squarify
from matplotlib.colors import LinearSegmentedColormap

# Loose -> dense colour gradient.  Light blue at d=0, orange-red at d=1.
DENSITY_CMAP = LinearSegmentedColormap.from_list(
    "loose_to_dense", ["#cfe5f5", "#f6c277", "#ef6548", "#a50f15"])

# Show at most this many ranked components per panel; the remainder
# collapse into one "other" rectangle so the visual stays compact.
MAX_RECTS = 6

fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))


def render_panel(ax, panel):
    s     = panel["s"]
    comps = panel["comps"]
    ranked = panel["comps_ranked"]

    # Pick the top-MAX_RECTS rectangles; collapse the rest into "other".
    head = ranked[:MAX_RECTS]
    tail = ranked[MAX_RECTS:]
    sizes = [it["n"] for it in head]
    if tail:
        sizes.append(sum(it["n"] for it in tail))
    # Squarify needs sizes in non-increasing order; head is already by
    # quality but treemap looks neater sorted by raw size.
    order = sorted(range(len(sizes)), key=lambda i: -sizes[i])
    sizes = [sizes[i] for i in order]
    labels_head = list(head) + ([{"n": sum(it["n"] for it in tail),
                                  "dens": 0.0,
                                  "tag":  "other"}] if tail else [])
    labels_ordered = [labels_head[i] for i in order]

    # Layout the treemap inside a unit square.
    norm_sizes = squarify.normalize_sizes(sizes, 1.0, 1.0)
    rects      = squarify.squarify(norm_sizes, 0, 0, 1.0, 1.0)

    for it, r in zip(labels_ordered, rects):
        x, y, w, h = r["x"], r["y"], r["dx"], r["dy"]
        # "other" rectangle stays neutral grey to keep visual focus on
        # the top components.
        if it["tag"] == "other":
            face = "#e7e7e7"
        else:
            face = DENSITY_CMAP(min(it["dens"], 1.0))
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=face,
                                   edgecolor="white", linewidth=1.2))
        # In-rectangle label: only if the rectangle is large enough to
        # fit text.  Use white text on dark fills, dark on light.
        txt_color = "#101010" if (face == "#e7e7e7" or it["dens"] < 0.35) else "white"
        if w > 0.18 and h > 0.10:
            primary = f"$n{{=}}{it['n']}$"
            if it["tag"] != "other":
                secondary = f"$d{{=}}${it['dens']:.2f}"
            else:
                secondary = ""
            ax.text(x + w/2, y + h/2 + 0.012, primary,
                    ha="center", va="center", fontsize=10.5,
                    color=txt_color, fontweight="bold")
            if secondary:
                ax.text(x + w/2, y + h/2 - 0.040, secondary,
                        ha="center", va="center", fontsize=8.5,
                        color=txt_color)
        elif w > 0.08 and h > 0.05:
            # Tight rectangle: single-line "n=N (d=X)".
            tag = f"$n{{=}}{it['n']}$"
            if it["tag"] != "other":
                tag += f", $d{{=}}${it['dens']:.2f}"
            ax.text(x + w/2, y + h/2, tag,
                    ha="center", va="center", fontsize=8,
                    color=txt_color, fontweight="bold")
        # Smaller rectangles get no label - the colour already carries
        # the density signal and the legend below decodes it.

    quality_word = "loose"    if panel["avgd_top"] < 0.30 else \
                   "moderate" if panel["avgd_top"] < 0.60 else "dense"
    ax.set_title(
        f"$s={s}$: {len(comps)} components, top-{TOP_K} avg density "
        f"$= {panel['avgd_top']:.2f}$  ({quality_word})",
        fontsize=11.5, fontweight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)


for ax, p in zip(axes, panels):
    render_panel(ax, p)

# Density legend bar at the bottom of the figure.
import matplotlib as mpl
cax = fig.add_axes([0.20, 0.04, 0.60, 0.025])
sm  = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1), cmap=DENSITY_CMAP)
cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
cbar.set_label("internal edge density $d$ (loose $\\rightarrow$ near-clique)",
               fontsize=9)
cbar.ax.tick_params(labelsize=8)

fig.tight_layout(rect=[0, 0.10, 1, 1])
fig.savefig(ROOT / "cs9_egonet.png", dpi=180, bbox_inches="tight")
fig.savefig(ROOT / "cs9_egonet.pdf", bbox_inches="tight")
print(f"\nFigures saved to {ROOT}")
