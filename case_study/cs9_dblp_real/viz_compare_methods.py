#!/usr/bin/env python3
"""Side-by-side: our (1,3)-nucleus vs NuclearCD (3,4)-nucleus on Han's 1-hop ego.

Left  panel: our (1,3) — Han ego analysed with r=1, s=3 (our method)
Right panel: NuclearCD (3,4) — same ego analysed with r=3, s=4 (reference)

For (3,4) we project the triangle-level core onto vertices
(vertex_core(v) = max_{t ni v} triangle_core(t)) and partition the alive
subgraph into connected components (which fragments naturally at r=3).

NuclearCD CS-II style annotations on both panels: top-K modules carry a
callout box with semantic label, member count, paper count, % active.
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
from matplotlib.patches import FancyBboxPatch
import networkx as nx

DIR = Path("/Users/zhangwenqian/UNSW/pivoter/case_study/cs9_dblp_real")
ANCHOR = 1373214
MIN_MOD, DENS = 10, 0.85

# ---- load original ids + adjacency + papers ----
t0 = time.time()
id2name = {}
with (DIR/"authors.tsv").open() as f:
    for line in f:
        i,n=line.rstrip("\n").split("\t",1); id2name[int(i)]=n

adj = defaultdict(set)
with (DIR/"edges.txt").open() as f:
    for line in f:
        u,v=line.split(); u,v=int(u),int(v); adj[u].add(v); adj[v].add(u)
print(f"adj loaded [{int(time.time()-t0)}s]")

hop1 = sorted(adj[ANCHOR] | {ANCHOR})
hop1_set = set(hop1)
adj_hop1 = {v: adj[v]&hop1_set for v in hop1}

# new_id <-> orig_id for Han ego (from han_ego.map)
new2orig = {}; orig2new = {}
with (DIR/"han_ego.map").open() as f:
    for line in f:
        if line.startswith("#"): continue
        i,v = line.rstrip().split("\t"); i,v = int(i),int(v)
        new2orig[i]=v; orig2new[v]=i
ANCHOR_NEW = orig2new[ANCHOR]
print(f"ego size = {len(hop1)}, Han new_id = {ANCHOR_NEW}")

# papers
hop1_names = {id2name[v]:v for v in hop1 if v in id2name}
papers = []
papers_venue = []
with (DIR/"papers_venue.tsv").open() as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts)<2: continue
        venue = parts[0]
        mem = {hop1_names[a] for a in parts[1:] if a in hop1_names}
        if len(mem)>=3:
            papers.append(mem); papers_venue.append(venue)
papers_of = defaultdict(list)
for pi,pa in enumerate(papers):
    for v in pa: papers_of[v].append(pi)
print(f"papers >=3 hop1: {len(papers)} [{int(time.time()-t0)}s]")

# ---- our (1,3) on full graph, restricted to hop1 ----
core_15 = {v:0 for v in hop1}
with (DIR/"r1_s3.tsv").open() as f:
    for line in f:
        if line.startswith("#"): continue
        i,c = line.rstrip().split(); i=int(i)
        if i in hop1_set: core_15[i]=int(c)
alive_13 = {v for v in hop1 if core_15[v] >= 1}
print(f"(1,3) alive: {len(alive_13)}")

# ---- NuclearCD (3,4) on Han ego: project triangle core -> vertex core + triangle list ----
vertex_core_34 = defaultdict(int)
triangles_34 = []   # list of (a,b,c) orig-id triples that are alive (core >= 1)
with (DIR/"han_r3_s4.tsv").open() as f:
    for line in f:
        if line.startswith("#"): continue
        ids, c = line.rstrip().split("\t")
        c = int(c)
        ids_orig = tuple(sorted(new2orig[int(x)] for x in ids.split()))
        if c >= 1:
            triangles_34.append(ids_orig)
        for v in ids_orig:
            if c > vertex_core_34[v]: vertex_core_34[v] = c
alive_34 = {v for v in hop1 if vertex_core_34.get(v,0) >= 1}
print(f"(3,4) alive: {len(alive_34)},  alive triangles: {len(triangles_34)}")

# Triangle-connected components: union-find on triangles sharing an edge.
# Build edge -> list of triangle indices, then union triangles that share edges.
class UF:
    def __init__(self, n): self.p = list(range(n))
    def find(self, x):
        while self.p[x]!=x:
            self.p[x] = self.p[self.p[x]]; x = self.p[x]
        return x
    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb: self.p[ra] = rb

uf = UF(len(triangles_34))
edge_to_tris = defaultdict(list)
for ti, (a,b,c) in enumerate(triangles_34):
    for e in [(a,b),(a,c),(b,c)]:
        edge_to_tris[e].append(ti)
for tris in edge_to_tris.values():
    for j in range(1, len(tris)):
        uf.union(tris[0], tris[j])

# Group triangles by root, then collect vertex sets
group_vs = defaultdict(set)
for ti, (a,b,c) in enumerate(triangles_34):
    r = uf.find(ti)
    group_vs[r].update((a,b,c))
comps_34_tri = sorted(group_vs.values(), key=len, reverse=True)
print(f"(3,4) triangle-connected nuclei: {len(comps_34_tri)} (>=4: {sum(1 for c in comps_34_tri if len(c)>=4)})")
print(f"  top sizes: {[len(c) for c in comps_34_tri[:8]]}")

# ---- partition definitions ----
def find_modules_dense(pool, anchor_excluded=True):
    pool = set(pool)
    if anchor_excluded: pool.discard(ANCHOR)
    found = []
    while len(pool)>=MIN_MOD and len(found)<200:
        deg = {x:len(adj_hop1[x]&pool) for x in pool}
        if not deg: break
        seed = max(deg, key=deg.get)
        if deg[seed]<MIN_MOD-1: pool.discard(seed); continue
        cl={seed}
        cand=sorted(adj_hop1[seed]&pool, key=lambda x:-len(adj_hop1[x]&pool))
        for c in cand:
            if c in cl: continue
            if len(adj_hop1[c]&cl)<DENS*len(cl): continue
            cl.add(c)
        if len(cl)>=MIN_MOD:
            e=sum(1 for u in cl for w in cl if u<w and w in adj_hop1[u])
            p=len(cl)*(len(cl)-1)//2
            if p and e/p>=DENS:
                found.append(cl); pool-=cl
            else: pool.discard(seed)
        else: pool.discard(seed)
    return found

def connected_components_of_alive(A):
    """Components of the alive subgraph (no anchor exclusion)."""
    G = nx.Graph(); G.add_nodes_from(A)
    for u in A:
        for v in adj_hop1[u]&A:
            if u<v: G.add_edge(u,v)
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    return [set(c) for c in comps]

# ---- semantic label ----
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
        if len(papers[pi]&C)>=3 and papers_venue[pi]!="CoRR":
            cnt[papers_venue[pi]] += 1
    if not cnt: return None
    score = defaultdict(int)
    for ven,c in cnt.most_common(25):
        for keys,label in VENUE_LABEL:
            if any(k in ven for k in keys):
                score[label] += c; break
    return max(score.items(), key=lambda x:x[1])[0] if score else cnt.most_common(1)[0][0][:18]

def label_for(C):
    reps = [id2name.get(v,"?") for v,_ in sorted(
        [(v,len(adj_hop1[v]&C)) for v in C], key=lambda x:-x[1])[:5]]
    for r in reps:
        if r in MANUAL_LABEL: return MANUAL_LABEL[r]
    return infer_label(C) or "Research Cluster"

def comp_metrics(C, A):
    """m_in, m_cut (to other-alive), active%, # collab papers (>=3 within-C)."""
    m_in=m_cut=0
    for u in C:
        for v in adj_hop1[u]:
            if v not in A or u==v: continue
            if v in C:
                if u<v: m_in += 1
            else:
                m_cut += 1
    cands=set()
    for v in C: cands.update(papers_of.get(v,[]))
    actives=set(); nP=0
    for pi in cands:
        inter = papers[pi]&C
        if len(inter)>=3: nP+=1; actives|=inter
    act = len(actives)/len(C) if C else 0.0
    return m_in, m_cut, act, nP

# ---- panels ----
# LEFT: our (1,3) — greedy dense modules (since r=1 alive is one giant comp)
mods_13 = find_modules_dense(alive_13)
mods_13 = sorted(mods_13, key=len, reverse=True)
TOPK_13 = 5
top_13 = mods_13[:TOPK_13]

# RIGHT: NuclearCD (3,4) at level K>=10: triangle-connected nuclei from triangles with core>=10.
# At K=1 on this 1-hop ego the (3,4)-nucleus is one giant blob (Han is central),
# so we report the K=10 shell of the (3,4) hierarchy where 5 sub-nuclei separate.
K_RIGHT = 10
class UF2:
    def __init__(self,n): self.p=list(range(n))
    def find(self,x):
        while self.p[x]!=x: self.p[x]=self.p[self.p[x]]; x=self.p[x]
        return x
    def union(self,a,b):
        ra,rb=self.find(a),self.find(b)
        if ra!=rb: self.p[ra]=rb
keep_tris = []
with (DIR/"han_r3_s4.tsv").open() as f:
    for line in f:
        if line.startswith("#"): continue
        ids,c = line.rstrip().split("\t"); c=int(c)
        if c>=K_RIGHT:
            keep_tris.append(tuple(sorted(new2orig[int(x)] for x in ids.split())))
uf2 = UF2(len(keep_tris))
edge_to2 = defaultdict(list)
for i,(a,b,cc) in enumerate(keep_tris):
    for e in [(a,b),(a,cc),(b,cc)]: edge_to2[e].append(i)
for tl in edge_to2.values():
    for j in range(1,len(tl)): uf2.union(tl[0], tl[j])
groups = defaultdict(set)
for i,(a,b,cc) in enumerate(keep_tris):
    r = uf2.find(i); groups[r].update((a,b,cc))
comps_34 = sorted(groups.values(), key=len, reverse=True)
alive_34_K = set().union(*comps_34) if comps_34 else set()
big_34 = [c for c in comps_34 if len(c)>=4]
TOPK_34 = min(5, len(big_34))
top_34 = big_34[:TOPK_34]
print(f"\n(1,3) modules: {len(mods_13)} top sizes={[len(c) for c in mods_13[:8]]}")
print(f"(3,4) nuclei at shell K>={K_RIGHT}: {len(comps_34)} (>=4: {len(big_34)}), alive={len(alive_34_K)}")
print(f"  top sizes: {[len(c) for c in comps_34[:10]]}")

# ---- summarise per panel ----
def panel_summary(top, A, name):
    info=[]
    for C in top:
        m_in,m_cut,act,nP = comp_metrics(C, A)
        info.append({
            "size":len(C), "m_in":m_in, "m_cut":m_cut,
            "active":round(act,2), "papers":nP,
            "label":label_for(C),
            "reps":[id2name.get(v,"?") for v,_ in sorted(
                [(v,len(adj_hop1[v]&C)) for v in C], key=lambda x:-x[1])[:4]],
        })
    return info

info_13 = panel_summary(top_13, alive_13, "(1,3)")
info_34 = panel_summary(top_34, alive_34_K, "(3,4)")

print("\n=== (1,3) top modules ===")
for d in info_13:
    print(f"  |C|={d['size']:>3d} m_in={d['m_in']:>3d} m_cut={d['m_cut']:>4d} "
          f"act={d['active']:.2f} P={d['papers']:>3d}  label={d['label']}")
    print(f"     reps: {', '.join(d['reps'])}")
print("\n=== (3,4) top components ===")
for d in info_34:
    print(f"  |C|={d['size']:>3d} m_in={d['m_in']:>3d} m_cut={d['m_cut']:>4d} "
          f"act={d['active']:.2f} P={d['papers']:>3d}  label={d['label']}")
    print(f"     reps: {', '.join(d['reps'])}")

# save summary
summary = {
    "ours_1_3":   {"n_alive":len(alive_13), "n_modules":len(mods_13), "top":info_13},
    "nuclearcd_3_4": {"n_alive":len(alive_34), "n_comps":len(comps_34),
                     "n_comps_ge4":len(big_34), "top":info_34},
}
(DIR/"cs9_compare.json").write_text(json.dumps(summary, indent=2))

# ---- spring layout (shared) ----
print(f"\nlayout [{int(time.time()-t0)}s]...")
G = nx.Graph(); G.add_nodes_from(hop1)
for u in hop1:
    for v in adj_hop1[u]:
        if u<v: G.add_edge(u,v)
pos = nx.spring_layout(G, seed=42, iterations=140, k=1.6/(len(hop1)**0.5))
print(f"done [{int(time.time()-t0)}s]")

PALETTE = ["#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a",
           "#b15928", "#a6cee3", "#fb9a99"]

def render_panel(ax, A, modules, info, title):
    peeled = set(hop1) - A
    # peeled
    nx.draw_networkx_nodes(G, pos, nodelist=peeled, ax=ax,
        node_color="#e9e9e9", node_size=4, alpha=0.5, linewidths=0)
    # other alive
    top_set = set().union(*modules) if modules else set()
    others = (A - top_set) - {ANCHOR}
    nx.draw_networkx_nodes(G, pos, nodelist=others, ax=ax,
        node_color="#b8b8b8", node_size=10, alpha=0.6, linewidths=0)
    # top modules
    centroids = []
    for k_idx, C in enumerate(modules):
        color = PALETTE[k_idx % len(PALETTE)]
        nodes = [v for v in C if v != ANCHOR]
        sizes = [25 + 6*len(adj_hop1[v]&C) for v in nodes]
        sizes = [min(s_, 220) for s_ in sizes]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax,
            node_color=color, node_size=sizes, alpha=0.92,
            edgecolors="white", linewidths=0.6)
        xs = np.array([pos[v][0] for v in nodes])
        ys = np.array([pos[v][1] for v in nodes])
        centroids.append((xs.mean(), ys.mean()))
    # anchor
    if ANCHOR in A:
        ax.scatter([pos[ANCHOR][0]], [pos[ANCHOR][1]], marker="o",
                   s=240, facecolor="#e6e6e6", edgecolors="#cc0000",
                   linewidths=1.8, zorder=8)
        ax.text(pos[ANCHOR][0], pos[ANCHOR][1] - 0.045, "Jiawei Han 0001",
                ha="center", va="top", fontsize=8, style="italic",
                color="#cc0000", zorder=9,
                bbox=dict(boxstyle="round,pad=0.18", fc="white",
                          ec="#cc0000", lw=0.8))
    # annotation boxes
    n = len(modules)
    if n>0:
        angles = np.linspace(35, 395, n+1)[:-1] * np.pi/180
        xs_a = np.array([pos[v][0] for v in hop1])
        ys_a = np.array([pos[v][1] for v in hop1])
        cx = (xs_a.min()+xs_a.max())/2; cy=(ys_a.min()+ys_a.max())/2
        rx = (xs_a.max()-xs_a.min())*0.62
        ry = (ys_a.max()-ys_a.min())*0.62
        for k_idx, d in enumerate(info):
            a = angles[k_idx]
            bx, by = cx + rx*np.cos(a), cy + ry*np.sin(a)
            color = PALETTE[k_idx % len(PALETTE)]
            text = (f"{d['label']}\n"
                    f"({d['size']} Members)\n"
                    f"({d['papers']} Papers, {int(d['active']*100)}% Active)")
            ax.text(bx, by, text, ha="center", va="center", fontsize=8.5,
                    family="serif", linespacing=1.18,
                    bbox=dict(boxstyle="round,pad=0.32",
                              fc="white", ec=color, lw=1.4), zorder=10)
            ctx, cty = centroids[k_idx]
            ax.annotate("", xy=(ctx, cty), xytext=(bx, by),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.8, alpha=0.75),
                zorder=4)
    ax.set_title(title, fontsize=11)
    ax.axis("off")

fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4))
xs_all = np.array([pos[v][0] for v in hop1])
ys_all = np.array([pos[v][1] for v in hop1])
xr = (xs_all.min()-0.05, xs_all.max()+0.05)
yr = (ys_all.min()-0.05, ys_all.max()+0.05)
for ax in axes: ax.set_xlim(*xr); ax.set_ylim(*yr)
render_panel(axes[0], alive_13, top_13, info_13,
             f"Ours $(r,s)=(1,3)$:  {len(alive_13)} alive,  "
             f"{len(mods_13)} dense modules")
render_panel(axes[1], alive_34_K, top_34, info_34,
             f"NuclearCD $(r,s)=(3,4)$ shell-{K_RIGHT}:  {len(alive_34_K)} alive,  "
             f"{len(big_34)} triangle-connected nuclei")

fig.tight_layout()
fig.savefig(DIR/"cs9_egonet.pdf", bbox_inches="tight")
fig.savefig(DIR/"cs9_egonet.png", bbox_inches="tight", dpi=180)
print(f"wrote cs9_egonet.pdf / .png [{int(time.time()-t0)}s]")
