"""
CS10 viz: per-graph join tree at s=5, with SHARED log-y axis = core value.

The shared y-axis makes the "deep / shallow" contrast across domains
visible at a glance:
  - com-dblp     -> nodes spread over many decades of k     (deep, branchy)
  - com-youtube  -> nodes pile up at low k                  (shallow)
  - web-Stanford -> intermediate                            (medium)

Data:
    cs10_hierarchy_tree.csv with columns
        graph,s,node_id,k_birth,k_death,size_birth,size_death,parent_id

Outputs:
    cs10_hierarchy.{pdf,png}
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ROOT     = Path(__file__).parent
TREE_CSV = ROOT / "cs10_hierarchy_tree.csv"
GRAPHS   = [
    ("com-dblp",     "#1f78b4"),
    ("com-youtube",  "#e31a1c"),
    ("web-Stanford", "#33a02c"),
]
S_FOR_TREE = 5
TOP_K      = 80


def load_tree():
    by_key = defaultdict(list)
    with open(TREE_CSV) as f:
        for r in csv.DictReader(f):
            by_key[(r["graph"], int(r["s"]))].append({
                "id":          int(r["node_id"]),
                "k_birth":     float(r["k_birth"]),
                "k_death":     float(r["k_death"]),
                "size_birth":  int(r["size_birth"]),
                "size_death":  int(r["size_death"]),
                "parent":      int(r["parent_id"]),
                "persistence": float(r["k_birth"]) - float(r["k_death"]),
            })
    return by_key


def filter_topk_with_ancestors(nodes, k):
    by_id    = {n["id"]: n for n in nodes}
    ranked   = sorted(nodes, key=lambda n: -n["persistence"])
    seed_ids = {n["id"] for n in ranked[:k]}
    kept     = set(seed_ids)
    for nid in list(seed_ids):
        p = by_id[nid]["parent"]
        while p != -1 and p not in kept:
            kept.add(p)
            p = by_id[p]["parent"]
    return [by_id[nid] for nid in kept]


def horizontal_layout(kept_nodes):
    """Recursively place children left-to-right; return x coord per node id."""
    by_id    = {n["id"]: n for n in kept_nodes}
    children = defaultdict(list)
    for n in kept_nodes:
        if n["parent"] in by_id:
            children[n["parent"]].append(n["id"])
    roots = [n["id"] for n in kept_nodes if n["parent"] not in by_id]

    xs     = {}
    cursor = [0.0]

    def place(nid):
        kids = sorted(children[nid], key=lambda i: -by_id[i]["persistence"])
        if not kids:
            x = cursor[0]; cursor[0] += 1.0
            xs[nid] = x; return x
        cx = [place(c) for c in kids]
        x  = sum(cx) / len(cx)
        xs[nid] = x; return x

    for r in roots:
        place(r); cursor[0] += 1.5
    return xs, by_id, children


def render_panel(ax, graph_label, color, kept, x_of, by_id, children):
    if not kept:
        ax.text(0.5, 0.5, "no branches", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="#999")
        ax.set_title(graph_label, fontsize=10, fontweight="bold")
        return

    smin = max(min(n["size_death"] for n in kept), 1)
    smax = max(n["size_death"] for n in kept)
    def marker_area(sz):
        if smax == smin: return 24.0
        t = (math.log(max(sz, 1)) - math.log(smin)) / (math.log(smax) - math.log(smin))
        return 12.0 + 90.0 * t

    segs = []
    for n in kept:
        for c_id in children[n["id"]]:
            x0 = x_of[n["id"]]; y0 = n["k_birth"] + 1
            x1 = x_of[c_id];    y1 = by_id[c_id]["k_birth"] + 1
            segs.append([(x0, y0), (x1, y1)])
    if segs:
        ax.add_collection(LineCollection(segs, colors="0.7",
                                         linewidths=0.6, zorder=1))

    xs = [x_of[n["id"]]       for n in kept]
    ys = [n["k_birth"] + 1    for n in kept]
    ar = [marker_area(n["size_death"]) for n in kept]
    ax.scatter(xs, ys, s=ar, c=color, edgecolors="white",
               linewidths=0.4, alpha=0.9, zorder=2)

    ax.set_title(graph_label, fontsize=10, fontweight="bold")
    ax.set_yscale("log")
    ax.set_xticks([])
    for sp in ("top", "right", "bottom"): ax.spines[sp].set_visible(False)
    ax.grid(True, axis="y", which="major", alpha=0.25, linestyle=":")


def main():
    by_key = load_tree()

    panels = []
    y_lo, y_hi = float("inf"), 0
    for gname, color in GRAPHS:
        nodes = by_key.get((gname, S_FOR_TREE), [])
        if not nodes:
            panels.append((gname, color, [], None, None, None)); continue
        kept = filter_topk_with_ancestors(nodes, TOP_K)
        x_of, byid, ch = horizontal_layout(kept)
        panels.append((gname, color, kept, x_of, byid, ch))
        for n in kept:
            y = n["k_birth"] + 1
            if y < y_lo and y > 0: y_lo = y
            if y > y_hi:           y_hi = y

    y_lo = max(y_lo / 2, 0.5)
    y_hi = y_hi * 2

    fig, axes = plt.subplots(1, len(GRAPHS), figsize=(8.5, 3.2), sharey=True)
    if len(GRAPHS) == 1: axes = [axes]
    for ax, (gname, color, kept, x_of, byid, ch) in zip(axes, panels):
        render_panel(ax, gname, color, kept, x_of, byid, ch)
        ax.set_ylim(y_lo, y_hi)

    axes[0].set_ylabel(r"core value $k$  (log)", fontsize=9.5)

    fig.tight_layout()
    pdf = ROOT / "cs10_hierarchy.pdf"
    png = ROOT / "cs10_hierarchy.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=200, bbox_inches="tight")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
