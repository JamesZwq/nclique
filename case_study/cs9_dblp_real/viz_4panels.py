#!/usr/bin/env python3
"""Four-panel comparison on Han's 1-hop ego: (1,4), (3,4), (1,10), (3,10)."""
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
MIN_MOD, DENS = 10, 0.85

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

# ---- load full hop1 adj ----
adj = defaultdict(set)
with (DIR/"edges.txt").open() as f:
    for line in f:
        u,v = line.split(); u,v = int(u), int(v); adj[u].add(v); adj[v].add(u)
hop1 = adj[ANCHOR] | {ANCHOR}
hop1_set = set(hop1)
adj_hop1 = {v: adj[v] & hop1_set for v in hop1}
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

# ---- semantic labels ----
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
    actives = set(); nP = 0
    for pi in cands:
        inter = papers[pi] & C
        if len(inter) >= 3: nP += 1; actives |= inter
    return m_in, m_cut, (len(actives)/len(C) if C else 0.0), nP

# ---- r=1 helpers ----
def load_r1(s):
    core = {}
    with (DIR/f"han_r1_s{s}.tsv").open() as f:
        for line in f:
            if line.startswith("#"): continue
            i,c = line.rstrip().split(); core[new2orig[int(i)]] = int(c)
    return {v for v in hop1 if core.get(v,0) >= 1}

def find_modules(pool):
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
            if len(adj_hop1[c]&cl) < DENS*len(cl): continue
            cl.add(c)
        if len(cl) >= MIN_MOD:
            e = sum(1 for u in cl for w in cl if u<w and w in adj_hop1[u])
            p = len(cl)*(len(cl)-1)//2
            if p and e/p >= DENS: found.append(cl); pool -= cl
            else: pool.discard(seed)
        else: pool.discard(seed)
    return found

# ---- r=3 helpers ----
class UF:
    def __init__(self,n): self.p = list(range(n))
    def find(self,x):
        r = x
        while self.p[r]!=r: r = self.p[r]
        while self.p[x]!=r: self.p[x], x = r, self.p[x]
        return r
    def union(self,a,b):
        ra,rb = self.find(a),self.find(b)
        if ra!=rb: self.p[ra]=rb

def load_r3_tri_components(s):
    """Triangle-connected components of the (3,s)-alive subgraph (core>=1)."""
    tris = []
    with (DIR/f"han_r3_s{s}.tsv").open() as f:
        for line in f:
            if line.startswith("#"): continue
            ids, _ = line.rstrip().split("\t")
            a,b,cc = ids.split()
            tris.append((new2orig[int(a)], new2orig[int(b)], new2orig[int(cc)]))
    uf = UF(len(tris))
    edge_to = defaultdict(list)
    for i, (a,b,c) in enumerate(tris):
        for e in [(a,b),(a,c),(b,c)]: edge_to[e].append(i)
    for tl in edge_to.values():
        for j in range(1, len(tl)): uf.union(tl[0], tl[j])
    groups = defaultdict(set)
    for i, (a,b,c) in enumerate(tris):
        r = uf.find(i); groups[r].update((a,b,c))
    comps = sorted(groups.values(), key=len, reverse=True)
    alive = set().union(*comps) if comps else set()
    return comps, alive

# ---- compute all 4 panels ----
PANELS_SPEC = [(1,4), (3,4), (1,10), (3,10)]
panel_results = {}
for r, s in PANELS_SPEC:
    stamp(f"computing (r={r},s={s})...")
    if r == 1:
        # r=1: dense modules (greedy near-cliques) on the (1,s)-alive subgraph
        A = load_r1(s)
        comps = sorted(find_modules(A), key=len, reverse=True)
        method = "dense modules"
    else:
        # r=3: triangle-connected nuclei (NuclearCD's native partition) on (3,s)-alive
        all_comps, A = load_r3_tri_components(s)
        comps = [c for c in all_comps if len(c) >= 4]
        method = "tri.-conn. nuclei"
    K_used = None
    TOPK = 5
    top = comps[:TOPK]
    info = []
    for C in top:
        m_in, m_cut, act, nP = comp_metrics(C, A)
        reps = [id2name.get(v,"?") for v,_ in sorted(
            [(v,len(adj_hop1[v]&C)) for v in C], key=lambda x:-x[1])[:4]]
        info.append({"size":len(C), "m_in":m_in, "m_cut":m_cut,
                     "active":round(act,2), "papers":nP,
                     "label":label_for(C), "reps":reps})
    panel_results[(r,s)] = {
        "alive_size": len(A),
        "n_comps": len(comps),
        "K_used": K_used,
        "method": method,
        "alive": A, "top": top, "info": info,
    }
    stamp(f"  (r={r},s={s}) alive={len(A)} comps={len(comps)} K_used={K_used} sizes={[len(c) for c in comps[:6]]}")
    for d in info:
        print(f"    |C|={d['size']:>3d} P={d['papers']:>3d} act={d['active']:.2f}  {d['label']}  reps={', '.join(d['reps'][:3])}")

# ---- spring layout ----
stamp("layout...")
G = nx.Graph(); G.add_nodes_from(hop1)
for u in hop1:
    for v in adj_hop1[u]:
        if u < v: G.add_edge(u,v)
pos = nx.spring_layout(G, seed=42, iterations=120, k=1.6/(len(hop1)**0.5))
stamp("layout done")

PALETTE = ["#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a",
           "#b15928", "#a6cee3", "#fb9a99"]

def render_panel(ax, r, s, data):
    A = data["alive"]; top = data["top"]; info = data["info"]
    peeled = hop1_set - A
    nx.draw_networkx_nodes(G, pos, nodelist=peeled, ax=ax,
        node_color="#ededed", node_size=3, alpha=0.45, linewidths=0)
    top_set = set().union(*top) if top else set()
    others = (A - top_set) - {ANCHOR}
    nx.draw_networkx_nodes(G, pos, nodelist=others, ax=ax,
        node_color="#a8a8a8", node_size=6, alpha=0.6, linewidths=0)
    centroids = []
    for k_idx, C in enumerate(top):
        color = PALETTE[k_idx % len(PALETTE)]
        nodes = [v for v in C if v != ANCHOR]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax,
            node_color=color, node_size=14, alpha=0.92,
            edgecolors="white", linewidths=0.3)
        xs = np.array([pos[v][0] for v in nodes])
        ys = np.array([pos[v][1] for v in nodes])
        centroids.append((xs.mean(), ys.mean()))
    if ANCHOR in A:
        ax.scatter([pos[ANCHOR][0]], [pos[ANCHOR][1]], marker="o",
                   s=140, facecolor="#e6e6e6", edgecolors="#cc0000",
                   linewidths=1.5, zorder=8)
        ax.text(pos[ANCHOR][0], pos[ANCHOR][1] - 0.04, "Jiawei Han 0001",
                ha="center", va="top", fontsize=6.5, style="italic",
                color="#cc0000", zorder=9,
                bbox=dict(boxstyle="round,pad=0.14", fc="white",
                          ec="#cc0000", lw=0.6))
    n = len(top)
    if n > 0:
        angles = np.linspace(40, 400, n+1)[:-1] * np.pi/180
        xs_a = np.array([pos[v][0] for v in hop1])
        ys_a = np.array([pos[v][1] for v in hop1])
        cx = (xs_a.min()+xs_a.max())/2; cy = (ys_a.min()+ys_a.max())/2
        rx = (xs_a.max()-xs_a.min())*0.68
        ry = (ys_a.max()-ys_a.min())*0.68
        for k_idx, d in enumerate(info):
            a = angles[k_idx]
            bx, by = cx + rx*np.cos(a), cy + ry*np.sin(a)
            color = PALETTE[k_idx % len(PALETTE)]
            text = (f"{d['label']}\n({d['size']} M)\n"
                    f"({d['papers']} P, {int(d['active']*100)}%)")
            ax.text(bx, by, text, ha="center", va="center", fontsize=6.6,
                    family="serif", linespacing=1.1,
                    bbox=dict(boxstyle="round,pad=0.22",
                              fc="white", ec=color, lw=1.0), zorder=10)
            ctx, cty = centroids[k_idx]
            ax.annotate("", xy=(ctx, cty), xytext=(bx, by),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.6, alpha=0.7),
                zorder=4)
    ax.set_title(f"$(r,s)=({r},{s})$:  {len(A)} alive,  {data['n_comps']} {data['method']}",
                 fontsize=10)
    ax.axis("off")

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
xs_all = np.array([pos[v][0] for v in hop1])
ys_all = np.array([pos[v][1] for v in hop1])
xr = (xs_all.min()-0.04, xs_all.max()+0.04)
yr = (ys_all.min()-0.04, ys_all.max()+0.04)
for ax_row in axes:
    for ax in ax_row:
        ax.set_xlim(*xr); ax.set_ylim(*yr)
for ax, (r,s) in zip(axes.flatten(), PANELS_SPEC):
    render_panel(ax, r, s, panel_results[(r,s)])

fig.tight_layout()
fig.savefig(DIR/"cs9_egonet.pdf", bbox_inches="tight")
fig.savefig(DIR/"cs9_egonet.png", bbox_inches="tight", dpi=170)
stamp("wrote cs9_egonet.pdf")

# JSON summary (drop alive/top sets — keep top info dicts only)
summary = {}
for k, v in panel_results.items():
    summary[f"r{k[0]}_s{k[1]}"] = {
        "alive_size": v["alive_size"], "n_comps": v["n_comps"],
        "K_used": v["K_used"], "method": v["method"], "top": v["info"],
    }
(DIR/"cs9_4panels.json").write_text(json.dumps(summary, indent=2))
stamp("wrote cs9_4panels.json")
