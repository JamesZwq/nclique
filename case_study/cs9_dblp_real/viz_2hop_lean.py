#!/usr/bin/env python3
"""(1,3) vs NuclearCD (3,4) on Han's 2-hop coauthor ego — memory-efficient.

Stream the 491MB triangle file in one pass to compute per-vertex max core.
"alive at (3,4)" = vertex with at least one alive triangle.
For visualisation we restrict to alive vertices only; layout is done on
that much smaller subgraph.
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
MIN_MOD, DENS = 10, 0.85

t0 = time.time()
def stamp(msg): print(f"[{int(time.time()-t0):>4d}s] {msg}", flush=True)

# ---- ids ----
id2name = {}
with (DIR/"authors.tsv").open() as f:
    for line in f:
        i,n=line.rstrip("\n").split("\t",1); id2name[int(i)]=n
stamp(f"authors {len(id2name)}")

new2orig = {}; orig2new = {}
with (DIR/"han_ego_2hop.map").open() as f:
    for line in f:
        if line.startswith("#"): continue
        i,v = line.rstrip().split("\t"); i,v=int(i),int(v)
        new2orig[i]=v; orig2new[v]=i
N = max(new2orig)+1
ANCHOR_NEW = orig2new[ANCHOR]
stamp(f"ego N={N}  Han new_id={ANCHOR_NEW}")

# ---- (1,3) per-vertex core ----
core_13 = [0]*N
with (DIR/"han_2hop_r1_s3.tsv").open() as f:
    for line in f:
        if line.startswith("#"): continue
        i,c=line.rstrip().split(); i,c=int(i),int(c)
        if i<N: core_13[i]=c
alive_13_new = {i for i in range(N) if core_13[i]>=1}
stamp(f"(1,3) alive = {len(alive_13_new)}")

# ---- (3,4) per-vertex max core via streaming ----
max_core_34 = [0]*N
n_tris = 0
with (DIR/"han_2hop_r3_s4.tsv").open() as f:
    for line in f:
        if line.startswith("#"): continue
        # format: "a b c\tcore"
        tab = line.rfind("\t")
        c = int(line[tab+1:])
        if c < 1: continue
        ids = line[:tab].split()
        a,b,cc = int(ids[0]),int(ids[1]),int(ids[2])
        if c > max_core_34[a]: max_core_34[a]=c
        if c > max_core_34[b]: max_core_34[b]=c
        if c > max_core_34[cc]: max_core_34[cc]=c
        n_tris += 1
alive_34_new = {i for i in range(N) if max_core_34[i]>=1}
stamp(f"(3,4) alive triangles {n_tris}, alive vertices = {len(alive_34_new)}")

# convert to orig ids
alive_13 = {new2orig[i] for i in alive_13_new}
alive_34 = {new2orig[i] for i in alive_34_new}
view = alive_13 | alive_34
stamp(f"alive_13 orig {len(alive_13)}, alive_34 orig {len(alive_34)}, view {len(view)}")

# ---- load adj restricted to view ----
adj_view = defaultdict(set)
with (DIR/"edges.txt").open() as f:
    for line in f:
        u,v=line.split(); u,v=int(u),int(v)
        if u in view and v in view:
            adj_view[u].add(v); adj_view[v].add(u)
stamp(f"adj_view |E|={sum(len(s) for s in adj_view.values())//2}")

# ---- (3,4) components: connected components of alive_34 in adj_view ----
G34 = nx.Graph(); G34.add_nodes_from(alive_34)
for u in alive_34:
    for v in adj_view[u]:
        if v in alive_34 and u<v: G34.add_edge(u,v)
comps_34 = sorted([set(c) for c in nx.connected_components(G34)], key=len, reverse=True)
stamp(f"(3,4) components: {len(comps_34)} top sizes={[len(c) for c in comps_34[:8]]}")

# If still one giant: threshold by max_core_34 quantile
giant = comps_34[0] if comps_34 else set()
if len(giant) > 0.6*len(alive_34) and len(comps_34) <= 3:
    # threshold by upper quartile
    cores_alive = sorted([max_core_34[orig2new[v]] for v in alive_34], reverse=True)
    Q = cores_alive[len(cores_alive)//4]   # 75th percentile
    stamp(f"  too-giant: try threshold core >= {Q} (75th pct)")
    inner_34 = {v for v in alive_34 if max_core_34[orig2new[v]] >= Q}
    G34i = nx.Graph(); G34i.add_nodes_from(inner_34)
    for u in inner_34:
        for v in adj_view[u]:
            if v in inner_34 and u<v: G34i.add_edge(u,v)
    comps_34 = sorted([set(c) for c in nx.connected_components(G34i)], key=len, reverse=True)
    alive_34_eff = inner_34
    LABEL_34 = f"(3,4) inner shell core$\\geq${Q}"
    stamp(f"  inner comps {len(comps_34)} top sizes={[len(c) for c in comps_34[:8]]}")
else:
    alive_34_eff = alive_34
    LABEL_34 = "(3,4) alive"

# ---- (1,3) dense modules in alive_13 (anchor excluded) ----
def find_modules(pool):
    pool=set(pool); pool.discard(ANCHOR); found=[]
    while len(pool)>=MIN_MOD and len(found)<300:
        deg={x:len(adj_view[x]&pool) for x in pool}
        if not deg: break
        seed=max(deg,key=deg.get)
        if deg[seed]<MIN_MOD-1: pool.discard(seed); continue
        cl={seed}
        cand=sorted(adj_view[seed]&pool, key=lambda x:-len(adj_view[x]&pool))
        for c in cand:
            if c in cl: continue
            if len(adj_view[c]&cl)<DENS*len(cl): continue
            cl.add(c)
        if len(cl)>=MIN_MOD:
            e=sum(1 for u in cl for w in cl if u<w and w in adj_view[u])
            p=len(cl)*(len(cl)-1)//2
            if p and e/p>=DENS: found.append(cl); pool-=cl
            else: pool.discard(seed)
        else: pool.discard(seed)
    return found

mods_13 = sorted(find_modules(alive_13), key=len, reverse=True)
stamp(f"(1,3) modules {len(mods_13)} top sizes={[len(c) for c in mods_13[:8]]}")

# ---- papers in view ----
view_names = {id2name[v]:v for v in view if v in id2name}
papers=[]; venues=[]
with (DIR/"papers_venue.tsv").open() as f:
    for line in f:
        parts=line.rstrip("\n").split("\t")
        if len(parts)<2: continue
        venue=parts[0]
        mem={view_names[a] for a in parts[1:] if a in view_names}
        if len(mem)>=3:
            papers.append(mem); venues.append(venue)
papers_of=defaultdict(list)
for pi,pa in enumerate(papers):
    for v in pa: papers_of[v].append(pi)
stamp(f"papers >=3 in view {len(papers)}")

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
def infer_label(C):
    cnt=Counter()
    for pi in {p for v in C for p in papers_of.get(v,[])}:
        if len(papers[pi]&C)>=3 and venues[pi]!="CoRR":
            cnt[venues[pi]]+=1
    score=defaultdict(int)
    for ven,c in cnt.most_common(25):
        for keys,label in VENUE_LABEL:
            if any(k in ven for k in keys):
                score[label]+=c; break
    if score: return max(score.items(),key=lambda x:x[1])[0]
    return cnt.most_common(1)[0][0][:18] if cnt else None

def comp_metrics(C, A):
    m_in=m_cut=0
    for u in C:
        for v in adj_view[u]:
            if v not in A or u==v: continue
            if v in C:
                if u<v: m_in+=1
            else: m_cut+=1
    cands=set()
    for v in C: cands.update(papers_of.get(v,[]))
    actives=set(); nP=0
    for pi in cands:
        inter=papers[pi]&C
        if len(inter)>=3: nP+=1; actives|=inter
    return m_in,m_cut,(len(actives)/len(C) if C else 0),nP

def summarise(top, A):
    info=[]
    for C in top:
        m_in,m_cut,act,nP=comp_metrics(C,A)
        reps=[id2name.get(v,"?") for v,_ in sorted(
            [(v,len(adj_view[v]&C)) for v in C], key=lambda x:-x[1])[:5]]
        info.append({"size":len(C),"m_in":m_in,"m_cut":m_cut,
                     "active":round(act,2),"papers":nP,
                     "label":infer_label(C) or "Research Cluster",
                     "reps":reps})
    return info

TOPK=5
big_34=[c for c in comps_34 if len(c)>=10]
top_13=mods_13[:TOPK]
top_34=big_34[:TOPK]
info_13=summarise(top_13, alive_13)
info_34=summarise(top_34, alive_34_eff)

print("\n=== (1,3) top modules ===")
for d in info_13:
    print(f"  |C|={d['size']:>4d} P={d['papers']:>4d} act={d['active']:.2f}  "
          f"label={d['label']}  reps={', '.join(d['reps'][:3])}")
print(f"\n=== (3,4) {LABEL_34} top components ===")
print(f"  total components: {len(comps_34)} (>=10: {len(big_34)})")
for d in info_34:
    print(f"  |C|={d['size']:>4d} P={d['papers']:>4d} act={d['active']:.2f}  "
          f"label={d['label']}  reps={', '.join(d['reps'][:3])}")

(DIR/"cs9_2hop.json").write_text(json.dumps({
    "ours_1_3":{"alive":len(alive_13),"modules":len(mods_13),"top":info_13},
    "nuclearcd_3_4":{"alive":len(alive_34),"alive_eff":len(alive_34_eff),
                     "components":len(comps_34),"comps_ge10":len(big_34),
                     "label":LABEL_34,"top":info_34},
}, indent=2))
stamp("saved cs9_2hop.json")

# ---- layout on view (union of alive sets, can be 5-50k typically) ----
stamp(f"layout |V|={len(view)}...")
G = nx.Graph(); G.add_nodes_from(view)
for u in view:
    for v in adj_view[u]:
        if u<v: G.add_edge(u,v)
stamp(f"  |E|={G.number_of_edges()}")
pos = nx.spring_layout(G, seed=42, iterations=60, k=2.0/(len(view)**0.5))
stamp(f"layout done")

PALETTE = ["#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a",
           "#b15928", "#a6cee3", "#fb9a99"]

def render(ax, A, top, info, title):
    peeled = view - A
    nx.draw_networkx_nodes(G, pos, nodelist=peeled, ax=ax,
        node_color="#ededed", node_size=2, alpha=0.4, linewidths=0)
    top_set = set().union(*top) if top else set()
    others = (A - top_set) - {ANCHOR}
    nx.draw_networkx_nodes(G, pos, nodelist=others, ax=ax,
        node_color="#a6a6a6", node_size=4, alpha=0.55, linewidths=0)
    centroids=[]
    for k_idx, C in enumerate(top):
        color = PALETTE[k_idx % len(PALETTE)]
        nodes = [v for v in C if v != ANCHOR]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax,
            node_color=color, node_size=12, alpha=0.92,
            edgecolors="white", linewidths=0.3)
        xs = np.array([pos[v][0] for v in nodes])
        ys = np.array([pos[v][1] for v in nodes])
        centroids.append((xs.mean(), ys.mean()))
    if ANCHOR in A:
        ax.scatter([pos[ANCHOR][0]], [pos[ANCHOR][1]], marker="o",
                   s=200, facecolor="#e6e6e6", edgecolors="#cc0000",
                   linewidths=1.6, zorder=8)
        ax.text(pos[ANCHOR][0], pos[ANCHOR][1] - 0.04, "Jiawei Han 0001",
                ha="center", va="top", fontsize=8, style="italic",
                color="#cc0000", zorder=9,
                bbox=dict(boxstyle="round,pad=0.18", fc="white",
                          ec="#cc0000", lw=0.7))
    n = len(top)
    if n>0:
        angles = np.linspace(35, 395, n+1)[:-1] * np.pi/180
        xs_a = np.array([pos[v][0] for v in view])
        ys_a = np.array([pos[v][1] for v in view])
        cx = (xs_a.min()+xs_a.max())/2; cy=(ys_a.min()+ys_a.max())/2
        rx = (xs_a.max()-xs_a.min())*0.62
        ry = (ys_a.max()-ys_a.min())*0.62
        for k_idx, d in enumerate(info):
            a = angles[k_idx]
            bx, by = cx + rx*np.cos(a), cy + ry*np.sin(a)
            color = PALETTE[k_idx % len(PALETTE)]
            text = (f"{d['label']}\n({d['size']} Members)\n"
                    f"({d['papers']} Papers, {int(d['active']*100)}% Active)")
            ax.text(bx, by, text, ha="center", va="center", fontsize=8.0,
                    family="serif", linespacing=1.15,
                    bbox=dict(boxstyle="round,pad=0.28",
                              fc="white", ec=color, lw=1.3), zorder=10)
            ctx, cty = centroids[k_idx]
            ax.annotate("", xy=(ctx, cty), xytext=(bx, by),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.8, alpha=0.75),
                zorder=4)
    ax.set_title(title, fontsize=10.5)
    ax.axis("off")

fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6))
for ax in axes:
    xs_all = np.array([pos[v][0] for v in view])
    ys_all = np.array([pos[v][1] for v in view])
    ax.set_xlim(xs_all.min()-0.05, xs_all.max()+0.05)
    ax.set_ylim(ys_all.min()-0.05, ys_all.max()+0.05)
render(axes[0], alive_13, top_13, info_13,
       f"Ours $(r,s)=(1,3)$ on Han 2-hop:  {len(alive_13)} alive,  "
       f"{len(mods_13)} dense modules")
render(axes[1], alive_34_eff, top_34, info_34,
       f"NuclearCD $(r,s)=(3,4)$ on Han 2-hop:  {len(alive_34_eff)} alive in "
       f"{LABEL_34},  {len(big_34)} comps ($|C|\\geq 10$)")
fig.tight_layout()
fig.savefig(DIR/"cs9_egonet.pdf", bbox_inches="tight")
fig.savefig(DIR/"cs9_egonet.png", bbox_inches="tight", dpi=180)
stamp(f"wrote cs9_egonet.pdf")
