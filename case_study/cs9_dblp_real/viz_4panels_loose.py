#!/usr/bin/env python3
"""Four panels on Han's 1-hop ego with LOOSE module extraction.

For each (r,s), keep the full (r,s)-alive subgraph, then extract
modules with looser density (0.5) and larger min_size (20). At small s
the alive set is big and includes peripheral collaborators, so some
modules pick up loose teams with <100% active. At large s the alive
set is already restricted to vertices in s-cliques, so any extracted
module is inherently tighter -> ~100% active.

r=1 vs r=3 at the same s should give similar active%, showing that
moving from r=1 to r=3 doesn't sharpen the community quality.
"""
import json, time
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt
import networkx as nx

DIR = Path("/Users/zhangwenqian/UNSW/pivoter/case_study/cs9_dblp_real")
ANCHOR = 1373214
MIN_MOD, DENS = 25, 0.25
# "active" = member in >=1 paper with >=COLLAB_K within-module coauthors
COLLAB_K = 5
ACT_PAPERS_REQ = 1

t0 = time.time()
def stamp(m): print(f"[{int(time.time()-t0):>3d}s] {m}", flush=True)

# ---- ids ----
id2name = {}
with (DIR/"authors.tsv").open() as f:
    for line in f:
        i,n = line.rstrip("\n").split("\t",1); id2name[int(i)]=n

new2orig = {}; orig2new = {}
with (DIR/"han_ego.map").open() as f:
    for line in f:
        if line.startswith("#"): continue
        i,v = line.rstrip().split("\t"); new2orig[int(i)] = int(v); orig2new[int(v)] = int(i)

# ---- adjacency on hop1 ----
adj = defaultdict(set)
with (DIR/"edges.txt").open() as f:
    for line in f:
        u,v=line.split(); u,v=int(u),int(v); adj[u].add(v); adj[v].add(u)
hop1 = adj[ANCHOR] | {ANCHOR}
hop1_set = set(hop1)
adj_hop1 = {v: adj[v]&hop1_set for v in hop1}
stamp(f"hop1 = {len(hop1)}")

# ---- papers ----
hop1_names = {id2name[v]:v for v in hop1 if v in id2name}
papers = []; venues = []
with (DIR/"papers_venue.tsv").open() as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2: continue
        v = parts[0]
        mem = {hop1_names[a] for a in parts[1:] if a in hop1_names}
        if len(mem) >= 3:
            papers.append(mem); venues.append(v)
papers_of = defaultdict(list)
for pi, pa in enumerate(papers):
    for v in pa: papers_of[v].append(pi)
stamp(f"papers >=3 hop1: {len(papers)}")

# ---- labels ----
VENUE_LABEL = [
    (("KDD", "ICDM", "WSDM", "DMKD"),               "Data Mining"),
    (("SIGMOD", "VLDB", "ICDE", "PODS", "EDBT",
      "TKDE", "Inf. Syst.", "Data Knowl. Eng."),    "DB Lab"),
    (("ACL", "NAACL", "EMNLP", "COLING", "TAC",
      "Comput. Linguistics"),                        "NLP Group"),
    (("CVPR", "ICCV", "ECCV", "BMVC", "IJCV", "TPAMI"),
                                                     "Vision Group"),
    (("SIGIR", "WWW", "CIKM", "RecSys"),            "IR / Web"),
    (("NeurIPS", "NIPS", "ICML", "AAAI", "IJCAI",
      "ICLR", "UAI", "JMLR"),                        "AI / ML"),
    (("S&P", "USENIX Security", "CCS", "NDSS"),     "Security"),
    (("BIBM", "Bioinformatics", "Nucleic Acids",
      "BMC Bioinform.", "PLOS Comput. Biol."),       "Bioinformatics"),
]
MANUAL_LABEL = {
    "Manling Li":            "Multimodal NLP",
    "Haoran Zhang":          "NLP (newer)",
    "Zhikun Zhang":          "ML Safety",
    "Yu-Che Tsai":           "Multimodal Web",
    "Vidhisha Balachandran": "NLP Group",
}
def infer_label(C):
    cnt = Counter()
    for pi in {p for v in C for p in papers_of.get(v,[])}:
        if len(papers[pi]&C)>=3 and venues[pi]!="CoRR":
            cnt[venues[pi]] += 1
    score = defaultdict(int)
    for ven,c in cnt.most_common(25):
        for keys,label in VENUE_LABEL:
            if any(k in ven for k in keys): score[label] += c; break
    if score: return max(score.items(), key=lambda x:x[1])[0]
    return cnt.most_common(1)[0][0][:18] if cnt else None
def label_for(C):
    reps = [id2name.get(v,"?") for v,_ in sorted(
        [(v,len(adj_hop1[v]&C)) for v in C], key=lambda x:-x[1])[:4]]
    for r in reps:
        if r in MANUAL_LABEL: return MANUAL_LABEL[r]
    return infer_label(C) or "Research Cluster"

def comp_metrics(C, A):
    m_in=m_cut=0
    for u in C:
        for v in adj_hop1[u]:
            if v not in A or u==v: continue
            if v in C:
                if u<v: m_in += 1
            else: m_cut += 1
    cands = set()
    for v in C: cands.update(papers_of.get(v,[]))
    paper_count = defaultdict(int)   # per-member: number of K-coauthor collab papers
    nP = 0
    for pi in cands:
        inter = papers[pi] & C
        if len(inter) >= COLLAB_K:
            nP += 1
            for v in inter:
                paper_count[v] += 1
    actives = {v for v, k in paper_count.items() if k >= ACT_PAPERS_REQ}
    return m_in, m_cut, (len(actives)/len(C) if C else 0.0), nP

# ---- per-vertex core ----
def core_r1(s):
    core = {v:0 for v in hop1}
    with (DIR/f"han_r1_s{s}.tsv").open() as f:
        for line in f:
            if line.startswith("#"): continue
            i,c = line.rstrip().split(); core[new2orig[int(i)]] = int(c)
    return core

def core_r3(s):
    core = defaultdict(int)
    with (DIR/f"han_r3_s{s}.tsv").open() as f:
        for line in f:
            if line.startswith("#"): continue
            ids,c = line.rstrip().split("\t"); c = int(c)
            for x in ids.split():
                v = new2orig[int(x)]
                if c > core[v]: core[v] = c
    return core

# ---- looser modules ----
def find_loose_modules(pool):
    pool = set(pool); pool.discard(ANCHOR); found = []
    while len(pool) >= MIN_MOD and len(found) < 300:
        deg = {x:len(adj_hop1[x]&pool) for x in pool}
        if not deg: break
        seed = max(deg, key=deg.get)
        if deg[seed] < MIN_MOD-1: pool.discard(seed); continue
        cl = {seed}
        cand = sorted(adj_hop1[seed]&pool, key=lambda x:-len(adj_hop1[x]&pool))
        for c in cand:
            if c in cl: continue
            # accept if at least DENS fraction of edges to current cluster exist
            if len(adj_hop1[c]&cl) < DENS*len(cl): continue
            cl.add(c)
        if len(cl) >= MIN_MOD:
            e = sum(1 for u in cl for w in cl if u<w and w in adj_hop1[u])
            p = len(cl)*(len(cl)-1)//2
            if p and e/p >= DENS: found.append(cl); pool -= cl
            else: pool.discard(seed)
        else: pool.discard(seed)
    return found

# ---- run all 4 panels ----
PANELS_SPEC = [(1,4), (3,4), (1,10), (3,10)]
panel = {}
for r, s in PANELS_SPEC:
    stamp(f"computing (r={r},s={s})...")
    core = core_r1(s) if r == 1 else core_r3(s)
    alive_set = {v for v in hop1 if core.get(v,0) >= 1}
    # loose modules on FULL alive subgraph
    mods = sorted(find_loose_modules(alive_set), key=len, reverse=True)
    TOPK = 5
    top = mods[:TOPK]
    info = []
    for C in top:
        m_in, m_cut, act, nP = comp_metrics(C, alive_set)
        density = m_in / max(len(C)*(len(C)-1)//2, 1)
        reps = [id2name.get(v,"?") for v,_ in sorted(
            [(v,len(adj_hop1[v]&C)) for v in C], key=lambda x:-x[1])[:4]]
        info.append({"size":len(C), "m_in":m_in, "m_cut":m_cut,
                     "density":round(density,2),
                     "active":round(act,2), "papers":nP,
                     "label":label_for(C), "reps":reps})
    panel[(r,s)] = {"alive_size": len(alive_set), "n_mods": len(mods),
                    "alive_set": alive_set, "top": top, "info": info}
    stamp(f"  (r={r},s={s}) alive={len(alive_set)} #mods={len(mods)} sizes={[len(c) for c in mods[:8]]}")
    for d in info:
        print(f"    |C|={d['size']:>3d} dens={d['density']:.2f} P={d['papers']:>3d} act={d['active']:.2f}  "
              f"{d['label']}  reps={', '.join(d['reps'][:3])}")

# ---- layout (more separation for clearer cluster visuals) ----
stamp("layout...")
G = nx.Graph(); G.add_nodes_from(hop1)
for u in hop1:
    for v in adj_hop1[u]:
        if u<v: G.add_edge(u,v)
pos = nx.spring_layout(G, seed=42, iterations=200, k=2.8/(len(hop1)**0.5))
stamp("layout done")

# Color palette tuned for visibility (colorblind-friendly, saturated)
PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
           "#8c564b", "#17becf", "#e377c2"]

# Mapping module label -> color, kept CONSISTENT across panels so the
# same working group appears in the same colour everywhere.
LABEL_COLOR = {}
def color_for_label(lbl, k_idx):
    if lbl not in LABEL_COLOR:
        LABEL_COLOR[lbl] = PALETTE[len(LABEL_COLOR) % len(PALETTE)]
    return LABEL_COLOR[lbl]

def render_panel(ax, r, s, data):
    A = data["alive_set"]; top = data["top"]; info = data["info"]
    peeled = hop1_set - A
    # very faint peeled nodes
    nx.draw_networkx_nodes(G, pos, nodelist=peeled, ax=ax,
        node_color="#d8d8d8", node_size=3, alpha=0.30, linewidths=0)
    # alive-but-not-in-top
    top_set = set().union(*top) if top else set()
    others = (A - top_set) - {ANCHOR}
    nx.draw_networkx_nodes(G, pos, nodelist=others, ax=ax,
        node_color="#9a9a9a", node_size=8, alpha=0.55, linewidths=0)
    # top-K modules with larger markers + consistent colour by label
    centroids = []
    for k_idx, C in enumerate(top):
        color = color_for_label(info[k_idx]["label"], k_idx)
        nodes = [v for v in C if v != ANCHOR]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax,
            node_color=color, node_size=42, alpha=0.95,
            edgecolors="white", linewidths=0.7)
        xs = np.array([pos[v][0] for v in nodes])
        ys = np.array([pos[v][1] for v in nodes])
        centroids.append((xs.mean(), ys.mean()))
    # anchor: clean red ring + italic label below
    if ANCHOR in A:
        ax.scatter([pos[ANCHOR][0]], [pos[ANCHOR][1]], marker="o",
                   s=220, facecolor="white", edgecolors="#cc0000",
                   linewidths=2.0, zorder=8)
        ax.text(pos[ANCHOR][0], pos[ANCHOR][1] - 0.055,
                "Jiawei Han 0001",
                ha="center", va="top", fontsize=8.5, style="italic",
                color="#cc0000", zorder=9)
    # inline labels: place each label NEAR its cluster centroid
    # offset radially outward from the global figure center so labels
    # don't sit inside the cluster.
    if centroids:
        xs_a = np.array([pos[v][0] for v in hop1])
        ys_a = np.array([pos[v][1] for v in hop1])
        cx_g = (xs_a.min()+xs_a.max())/2; cy_g = (ys_a.min()+ys_a.max())/2
        span_x = xs_a.max()-xs_a.min(); span_y = ys_a.max()-ys_a.min()
        for k_idx, (cx, cy) in enumerate(centroids):
            d = info[k_idx]; color = color_for_label(d["label"], k_idx)
            # radial offset outward from global centre
            dx, dy = cx - cx_g, cy - cy_g
            norm = (dx**2 + dy**2)**0.5 + 1e-9
            offx = (dx/norm) * span_x * 0.18
            offy = (dy/norm) * span_y * 0.18
            lx, ly = cx + offx, cy + offy
            text = (f"{d['label']}\n"
                    f"{d['size']} m, {int(d['active']*100)}% active")
            ax.text(lx, ly, text, ha="center", va="center", fontsize=8.5,
                    family="serif", linespacing=1.2, fontweight="bold",
                    color=color, zorder=10,
                    bbox=dict(boxstyle="round,pad=0.30",
                              fc="white", ec=color, lw=1.3, alpha=0.95))
    avg_act = (sum(d['active']*d['size'] for d in info) / max(sum(d['size'] for d in info), 1)) if info else 0
    method = "ours" if r == 1 else "NuclearCD"
    ax.set_title(f"$(r,s)\\!=\\!({r},{s})$ ({method}):  "
                 f"{len(A)} alive,  {data['n_mods']} modules,  "
                 f"avg active = {int(avg_act*100)}\\%",
                 fontsize=11.5, pad=6)
    ax.axis("off")

# 2x2 grid with breathing room
fig, axes = plt.subplots(2, 2, figsize=(14.5, 12),
                         gridspec_kw={"hspace":0.16, "wspace":0.05})
xs_all = np.array([pos[v][0] for v in hop1])
ys_all = np.array([pos[v][1] for v in hop1])
xr = (xs_all.min()-0.06, xs_all.max()+0.06)
yr = (ys_all.min()-0.06, ys_all.max()+0.06)
for ax_row in axes:
    for ax in ax_row:
        ax.set_xlim(*xr); ax.set_ylim(*yr)
for ax, (r,s) in zip(axes.flatten(), PANELS_SPEC):
    render_panel(ax, r, s, panel[(r,s)])

fig.subplots_adjust(left=0.02, right=0.99, top=0.96, bottom=0.02)
fig.savefig(DIR/"cs9_egonet.pdf", bbox_inches="tight")
fig.savefig(DIR/"cs9_egonet.png", bbox_inches="tight", dpi=180)
stamp("wrote cs9_egonet.pdf")

summary = {f"r{r}_s{s}": {
    "alive_size": d["alive_size"], "n_mods": d["n_mods"], "top": d["info"],
} for (r,s), d in panel.items()}
(DIR/"cs9_4panels.json").write_text(json.dumps(summary, indent=2))
stamp("wrote cs9_4panels.json")
