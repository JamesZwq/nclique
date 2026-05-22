#!/usr/bin/env python3
"""Case Study 4: Jiawei Han 1-hop ego, component decomposition under r=1 nucleus.

Mirrors NuclearCD CS-II's framing:
  - For each s in {LOW, HIGH}: take alive subgraph (core_{1,s}(v) >= 1) inside Han's 1-hop ego.
  - Find connected components, sort by size.
  - Report aggregate stats: |alive|, #comps, intra-component edge fraction,
    size-weighted avg separability m_in/m_cut, size-weighted avg conductance phi(C)=m_cut/(2 m_in + m_cut).
  - Top-K components: per-component sep, cond, active ratio (fraction of members in
    at least one DBLP paper containing >=3 within-component co-authors), top-3 representative names.
  - Saves cs9_components.json with all stats, and renders cs9_egonet.pdf as a 2-panel
    figure colored by component (top-K distinct colors, others gray).
"""
import json, time, sys
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx

DIR     = Path("/Users/zhangwenqian/UNSW/pivoter/case_study/cs9_dblp_real")
ANCHOR  = 1373214
ANCHOR_NAME = "Jiawei Han"
S_LOW, S_HIGH = 3, 10
PANELS  = (S_LOW, S_HIGH)
TOPK    = {S_LOW: 5, S_HIGH: 8}    # NuclearCD: top-5 at low s, top-8 at high s
MIN_COMP_SIZE = 4                  # ignore tiny components

# -------- load --------
t0 = time.time()
id2name = {}
with (DIR/"authors.tsv").open() as f:
    for line in f:
        i, n = line.rstrip("\n").split("\t", 1)
        id2name[int(i)] = n
print(f"authors loaded ({len(id2name)}) [{int(time.time()-t0)}s]")

adj = defaultdict(set)
with (DIR/"edges.txt").open() as f:
    for line in f:
        u, v = line.split()
        u, v = int(u), int(v)
        adj[u].add(v); adj[v].add(u)
print(f"adj loaded [{int(time.time()-t0)}s]")

hop1 = adj[ANCHOR] | {ANCHOR}
print(f"|1-hop ego| = {len(hop1)}")

core = {s: {} for s in PANELS}
for s in PANELS:
    with (DIR/f"r1_s{s}.tsv").open() as f:
        for line in f:
            if line.startswith("#"): continue
            i, c = line.rstrip().split()
            i = int(i)
            if i in hop1: core[s][i] = int(c)

# NuclearCD definition: "surviving" = participates in at least one (r,s)-nucleus.
# That is, core_{r,s}(v) >= 1 (the unit of the nucleus core number).
alive = {s: {v for v in hop1 if core[s].get(v, 0) >= 1} for s in PANELS}
print({s: len(alive[s]) for s in PANELS})

adj_hop1 = {v: (adj[v] & hop1) for v in hop1}

# -------- papers restricted to hop1 --------
hop1_names = {id2name[v]: v for v in hop1 if v in id2name}
papers = []
with (DIR/"papers.tsv").open() as f:
    for line in f:
        names = [a for a in line.rstrip("\n").split("\t") if a]
        mem = {hop1_names[a] for a in names if a in hop1_names}
        if len(mem) >= 3:
            papers.append(mem)
print(f"papers >=3 hop1-authors: {len(papers)} [{int(time.time()-t0)}s]")

# pre-bucket papers by author for quick lookup
papers_of = defaultdict(list)
for pi, pa in enumerate(papers):
    for v in pa:
        papers_of[v].append(pi)

# -------- per-s analysis --------
results = {}
panel_data = {}     # s -> (comps_sorted, comp_id_map, sizes, sep, cond, active)
for s in PANELS:
    A = alive[s]
    Gs = nx.Graph(); Gs.add_nodes_from(A)
    for u in A:
        for v in (adj_hop1[u] & A):
            if u < v: Gs.add_edge(u, v)

    comps = sorted(nx.connected_components(Gs), key=len, reverse=True)
    comp_id = {}
    for ci, C in enumerate(comps):
        for v in C: comp_id[v] = ci

    sizes = [len(C) for C in comps]
    m_in = [0]*len(comps); m_cut = [0]*len(comps)
    for (u, v) in Gs.edges():
        if comp_id[u] == comp_id[v]:
            m_in[comp_id[u]] += 1
        else:
            m_cut[comp_id[u]] += 1
            m_cut[comp_id[v]] += 1
    total_intra = sum(m_in)
    total_inter = sum(m_cut) // 2     # double-counted above
    intra_frac  = total_intra / max(total_intra + total_inter, 1)

    sep  = [(m_in[i] / m_cut[i]) if m_cut[i] > 0 else float("inf") for i in range(len(comps))]
    cond = [(m_cut[i] / (2*m_in[i] + m_cut[i])) if (2*m_in[i] + m_cut[i]) > 0 else 0.0
            for i in range(len(comps))]

    SEP_CAP = 20.0    # for size-weighted average: cap isolated-component sep
    sep_capped = [min(x, SEP_CAP) for x in sep]
    total_sz   = sum(sizes)
    avg_sep    = sum(sep_capped[i]*sizes[i] for i in range(len(comps))) / total_sz
    avg_cond   = sum(cond[i]*sizes[i]       for i in range(len(comps))) / total_sz

    # active ratio per comp under within-component collaboration threshold = 3
    active = []
    n_collab_papers = []   # how many papers had >=3 from comp
    for ci, C in enumerate(comps):
        active_members = set()
        cnt = 0
        candidate_papers = set()
        for v in C:
            for pi in papers_of.get(v, []):
                candidate_papers.add(pi)
        for pi in candidate_papers:
            inter = papers[pi] & C
            if len(inter) >= 3:
                cnt += 1
                active_members |= inter
        active.append(len(active_members) / len(C) if C else 0.0)
        n_collab_papers.append(cnt)

    # top-K
    k = TOPK[s]
    eligible = [i for i in range(len(comps)) if sizes[i] >= MIN_COMP_SIZE]
    top_idx = sorted(eligible, key=lambda i: -sizes[i])[:k]
    top_avg_active = (sum(active[i]*sizes[i] for i in top_idx) /
                      sum(sizes[i] for i in top_idx)) if top_idx else 0.0
    top_avg_active_unw = sum(active[i] for i in top_idx) / max(len(top_idx), 1)

    # representative names per top component: highest within-comp degree
    def reps(C, n=3):
        deg = [(v, len(adj_hop1[v] & C)) for v in C]
        deg.sort(key=lambda kv: -kv[1])
        return [id2name.get(v, f"id={v}") for v, _ in deg[:n]]

    top_info = []
    for i in top_idx:
        top_info.append({
            "comp_id":   i,
            "size":      sizes[i],
            "m_in":      m_in[i],
            "m_cut":     m_cut[i],
            "sep":       sep[i] if sep[i] != float("inf") else None,
            "cond":      cond[i],
            "active":    active[i],
            "collab_papers": n_collab_papers[i],
            "reps":      reps(comps[i], 4),
            "has_anchor": ANCHOR in comps[i],
        })

    results[s] = {
        "n_alive": len(A),
        "n_comps": len(comps),
        "n_comps_ge4": sum(1 for sz in sizes if sz >= MIN_COMP_SIZE),
        "total_intra_edges": total_intra,
        "total_inter_edges": total_inter,
        "intra_frac":  intra_frac,
        "avg_sep":     avg_sep,
        "avg_cond":    avg_cond,
        "topK_avg_active":     top_avg_active,
        "topK_avg_active_unw": top_avg_active_unw,
        "top": top_info,
    }
    panel_data[s] = (comps, comp_id, sizes, sep, cond, active, top_idx)

    print(f"\n=== s={s} ===")
    print(f"  alive={len(A)}  comps={len(comps)} (>=4: {sum(1 for sz in sizes if sz>=MIN_COMP_SIZE)})")
    print(f"  intra_frac={intra_frac:.3f}  avg_sep={avg_sep:.3f}  avg_cond={avg_cond:.3f}")
    print(f"  top-{k} avg active ratio (size-weighted) = {top_avg_active:.3f}")
    print(f"  top-{k} avg active ratio (unweighted)    = {top_avg_active_unw:.3f}")
    for info in top_info:
        anc = " [anchor]" if info["has_anchor"] else ""
        sep_str = f"{info['sep']:.2f}" if info["sep"] is not None else "inf"
        print(f"   C{info['comp_id']:>3d}: |C|={info['size']:>4d} "
              f"sep={sep_str:>6s} cond={info['cond']:.2f} "
              f"act={info['active']:.2f} #collabP={info['collab_papers']:>3d}{anc}")
        print(f"       reps: {', '.join(info['reps'])}")

(DIR/"cs9_components.json").write_text(json.dumps(results, indent=2, default=str))
print(f"\nwrote cs9_components.json [{int(time.time()-t0)}s]")

# -------- render figure --------
print("building global layout...")
keep = hop1
G = nx.Graph(); G.add_nodes_from(keep)
for u in keep:
    for v in adj_hop1[u]:
        if u < v: G.add_edge(u, v)
pos = nx.spring_layout(G, seed=42, iterations=120, k=1.6/(len(keep)**0.5))
print("  layout done")

# Color palette for top components (qualitative, colorblind-friendlyish)
PALETTE = [
    "#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a",
    "#b15928", "#a6cee3", "#fb9a99", "#b2df8a", "#fdbf6f",
]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
for ax, s in zip(axes, PANELS):
    comps, comp_id, sizes, sep, cond, active, top_idx = panel_data[s]
    A = alive[s]
    peeled = keep - A

    # peeled (gray)
    nx.draw_networkx_nodes(G, pos, nodelist=peeled,
        node_color="#e6e6e6", node_size=3, alpha=0.55, ax=ax, linewidths=0)

    # alive but not in a top component: light gray
    top_set = set(top_idx)
    in_top = set().union(*(comps[i] for i in top_idx)) if top_idx else set()
    other_alive = (A - in_top) - {ANCHOR}
    nx.draw_networkx_nodes(G, pos, nodelist=other_alive,
        node_color="#9a9a9a", node_size=8, alpha=0.7, ax=ax, linewidths=0)

    # each top-K component: distinct color
    for k_idx, ci in enumerate(top_idx):
        c = PALETTE[k_idx % len(PALETTE)]
        nodes = list(comps[ci] - {ANCHOR})
        nx.draw_networkx_nodes(G, pos, nodelist=nodes,
            node_color=c, node_size=22, alpha=0.95, ax=ax,
            linewidths=0.3, edgecolors="white")

    # anchor on top
    if ANCHOR in A:
        ax.scatter([pos[ANCHOR][0]], [pos[ANCHOR][1]], marker="*",
                   s=180, c="black", zorder=6, linewidths=0)

    n_comps_show = len(top_idx)
    ax.set_title(f"$s={s}$:  {len(A)} alive,  {results[s]['n_comps_ge4']} components ($|C|\\ge 4$),  "
                 f"top-{n_comps_show} colored",
                 fontsize=10)
    ax.axis("off")

# build a single legend that shows top-K rank colors
max_k = max(TOPK.values())
leg_handles = [mpatches.Patch(color=PALETTE[i % len(PALETTE)],
                              label=f"rank {i+1}") for i in range(max_k)]
leg_handles += [mpatches.Patch(color="#9a9a9a", label="other alive"),
                mpatches.Patch(color="#e6e6e6", label="peeled")]
fig.legend(handles=leg_handles, loc="lower center", ncol=max_k+2,
           fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig(DIR/"cs9_egonet.pdf", bbox_inches="tight")
fig.savefig(DIR/"cs9_egonet.png", bbox_inches="tight", dpi=180)
print("wrote cs9_egonet.pdf / .png")
