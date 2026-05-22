#!/usr/bin/env python3
"""NuclearCD CS-II style figure: 2 panels (s=3 vs s=10) with annotated dense modules.

Each top-K module gets:
  - distinct color
  - bigger node markers (members)
  - annotation box at panel margin with:
      "<Lab Name>"  /  "(N Members)"  /  "(N Papers, X% Active)"
  - colored arrow from box to module centroid
Anchor (Jiawei Han) labeled in red italic in center.
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
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import networkx as nx

DIR     = Path("/Users/zhangwenqian/UNSW/pivoter/case_study/cs9_dblp_real")
ANCHOR  = 1373214
PANELS  = (3, 10)
TOPK    = {3: 5, 10: 3}
MIN_MOD, DENS = 10, 0.85

# ---- load ----
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

hop1 = adj[ANCHOR] | {ANCHOR}
adj_hop1 = {v: adj[v]&hop1 for v in hop1}

core = {s:{} for s in PANELS}
for s in PANELS:
    with (DIR/f"r1_s{s}.tsv").open() as f:
        for line in f:
            if line.startswith("#"): continue
            i,c=line.rstrip().split(); i=int(i)
            if i in hop1: core[s][i]=int(c)
alive = {s: {v for v in hop1 if core[s].get(v,0)>=1} for s in PANELS}

# papers
hop1_names = {id2name[v]:v for v in hop1 if v in id2name}
papers = []
papers_venue = []
with (DIR/"papers_venue.tsv").open() as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2: continue
        venue = parts[0]
        mem = {hop1_names[a] for a in parts[1:] if a in hop1_names}
        if len(mem)>=3:
            papers.append(mem); papers_venue.append(venue)
papers_of = defaultdict(list)
for pi,pa in enumerate(papers):
    for v in pa: papers_of[v].append(pi)
print(f"papers >=3 hop1: {len(papers)} [{int(time.time()-t0)}s]")

# ---- module extraction (greedy near-clique, anchor excluded) ----
def find_modules(pool):
    pool=set(pool); pool.discard(ANCHOR)
    found=[]
    while len(pool)>=MIN_MOD and len(found)<200:
        deg={x:len(adj_hop1[x]&pool) for x in pool}
        if not deg: break
        seed=max(deg,key=deg.get)
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

# ---- semantic label: from top non-CoRR venue + manual override by representative author ----
VENUE_LABEL = [
    (("KDD", "ICDM", "WSDM", "DMKD"),                          "Data Mining"),
    (("SIGMOD", "VLDB", "ICDE", "PODS", "EDBT", "TKDE",
      "Inf. Syst.", "Data Knowl. Eng."),                       "DB Lab"),
    (("ACL", "NAACL", "EMNLP", "COLING", "TAC", "Comput. Linguistics"),
                                                                "NLP Group"),
    (("CVPR", "ICCV", "ECCV", "BMVC", "IJCV", "TPAMI"),        "Vision Group"),
    (("SIGIR", "WWW", "CIKM", "RecSys"),                       "IR / Web"),
    (("NeurIPS", "NIPS", "ICML", "AAAI", "IJCAI", "ICLR",
      "UAI", "JMLR"),                                          "AI / ML"),
    (("S&P", "USENIX Security", "CCS", "NDSS"),                "Security"),
]
# Manual overrides keyed by top representative member's name (after venue inference)
MANUAL_LABEL = {
    "Manling Li":       "Multimodal NLP",   # ACL/EMNLP heavy, multimodal event
    "Haoran Zhang":     "NLP (newer)",      # NAACL/ACL demos
    "Zhikun Zhang":     "ML Safety",        # privacy/security/ML (Suman Jana, etc.)
    "Yu-Che Tsai":      "Multimodal Web",   # WWW companion + multimodal authors
    "Vidhisha Balachandran": "NLP Group",   # Yulia Tsvetkov leads
}

def infer_label(C):
    cnt = Counter()
    for pi in {p for v in C for p in papers_of.get(v,[])}:
        if len(papers[pi] & C) >= 3 and papers_venue[pi] != "CoRR":
            cnt[papers_venue[pi]] += 1
    if not cnt: return None
    top = cnt.most_common(20)
    score = defaultdict(int)
    for ven, c in top:
        for keys, label in VENUE_LABEL:
            if any(k in ven for k in keys):
                score[label] += c; break
    if score:
        return max(score.items(), key=lambda x: x[1])[0]
    return top[0][0][:18]

def label_for(C):
    reps = [id2name.get(v,"?") for v,_ in sorted(
        [(v,len(adj_hop1[v]&C)) for v in C], key=lambda x:-x[1])[:4]]
    for r in reps:
        if r in MANUAL_LABEL:
            return MANUAL_LABEL[r]
    return infer_label(C) or "Research Cluster"

# ---- panel data ----
results = {}
panel = {}
for s in PANELS:
    A = alive[s]
    mods = find_modules(A)
    top_idx = sorted(range(len(mods)), key=lambda i: -len(mods[i]))[:TOPK[s]]
    edges_alive = sum(1 for u in A for v in adj_hop1[u]&A if u<v)
    mod_id={}
    for ci,M in enumerate(mods):
        for v in M: mod_id[v]=ci
    m_in=[0]*len(mods); m_cut=[0]*len(mods)
    for u in A:
        if u not in mod_id: continue
        ci=mod_id[u]
        for v in adj_hop1[u]:
            if v not in A or v==u: continue
            if v in mod_id and mod_id[v]==ci:
                if u<v: m_in[ci]+=1
            else:
                m_cut[ci]+=1
    sep = [(m_in[i]/m_cut[i]) if m_cut[i]>0 else float("inf") for i in range(len(mods))]
    cond= [(m_cut[i]/(2*m_in[i]+m_cut[i])) if (2*m_in[i]+m_cut[i])>0 else 0.0 for i in range(len(mods))]
    info=[]
    for ci in top_idx:
        C = mods[ci]
        cands=set()
        for v in C: cands.update(papers_of.get(v,[]))
        actives=set(); nP=0
        for pi in cands:
            inter = papers[pi]&C
            if len(inter)>=3:
                nP+=1; actives|=inter
        label = label_for(C)
        info.append({
            "ci": ci, "size": len(C), "m_in": m_in[ci], "m_cut": m_cut[ci],
            "sep": (None if sep[ci]==float("inf") else round(sep[ci],2)),
            "cond": round(cond[ci],2),
            "papers": nP, "active": round(len(actives)/len(C),2),
            "label": label,
            "reps": [id2name.get(v,"?") for v,_ in sorted(
                [(v,len(adj_hop1[v]&C)) for v in C], key=lambda x:-x[1])[:4]],
        })
    sizes=[len(c) for c in mods]; tot=sum(sizes)
    intra_frac = sum(m_in)/max(edges_alive,1)
    SEP_CAP=20.0
    avg_sep = sum(min(sep[i],SEP_CAP)*sizes[i] for i in range(len(mods)))/tot
    avg_cond= sum(cond[i]*sizes[i] for i in range(len(mods)))/tot
    results[s] = {
        "n_alive": len(A), "n_modules": len(mods),
        "intra_frac": round(intra_frac,3),
        "avg_sep": round(avg_sep,3),
        "avg_cond": round(avg_cond,3),
        "top": info,
    }
    panel[s] = (mods, top_idx, info)
    print(f"\ns={s}: alive={len(A)} mods={len(mods)} intra={intra_frac:.2f} sep={avg_sep:.2f} cond={avg_cond:.2f}")
    for d in info:
        print(f"  M{d['ci']:>2d} sz={d['size']:>3d}  label={d['label']:<14s} "
              f"papers={d['papers']:>3d} active={d['active']*100:.0f}%  reps={', '.join(d['reps'])}")

(DIR/"cs9_modules.json").write_text(json.dumps(results, indent=2))

# ---- shared spring layout ----
print(f"\nlaying out [{int(time.time()-t0)}s]...")
G = nx.Graph(); G.add_nodes_from(hop1)
for u in hop1:
    for v in adj_hop1[u]:
        if u<v: G.add_edge(u,v)
pos = nx.spring_layout(G, seed=42, iterations=140, k=1.6/(len(hop1)**0.5))
print(f"  done [{int(time.time()-t0)}s]")

# Per-module color palette
PALETTE = ["#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a",
           "#b15928", "#a6cee3", "#fb9a99"]

# ---- render ----
def annotate_modules(ax, s, panel_xy_range):
    mods, top_idx, info = panel[s]
    A = alive[s]
    peeled = hop1 - A

    # peeled vertices (faint)
    nx.draw_networkx_nodes(G, pos, nodelist=peeled, ax=ax,
        node_color="#e9e9e9", node_size=4, alpha=0.5, linewidths=0)

    top_set = set().union(*(mods[i] for i in top_idx)) if top_idx else set()
    others = (A - top_set) - {ANCHOR}
    nx.draw_networkx_nodes(G, pos, nodelist=others, ax=ax,
        node_color="#b8b8b8", node_size=10, alpha=0.6, linewidths=0)

    # top-K modules
    centroids = {}
    for k_idx, ci in enumerate(top_idx):
        color = PALETTE[k_idx % len(PALETTE)]
        nodes = [v for v in mods[ci] if v != ANCHOR]
        # node size scaled by within-module degree
        sizes = [25 + 6*len(adj_hop1[v] & mods[ci]) for v in nodes]
        sizes = [min(s_, 220) for s_ in sizes]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax,
            node_color=color, node_size=sizes, alpha=0.92,
            edgecolors="white", linewidths=0.6)
        xs = np.array([pos[v][0] for v in nodes])
        ys = np.array([pos[v][1] for v in nodes])
        centroids[ci] = (xs.mean(), ys.mean())

    # anchor on top
    if ANCHOR in A:
        ax.scatter([pos[ANCHOR][0]], [pos[ANCHOR][1]], marker="o",
                   s=240, facecolor="#e6e6e6", edgecolors="#cc0000",
                   linewidths=1.8, zorder=8)
        ax.text(pos[ANCHOR][0], pos[ANCHOR][1] - 0.045, "Jiawei Han 0001",
                ha="center", va="top", fontsize=7.5, style="italic",
                color="#cc0000", zorder=9,
                bbox=dict(boxstyle="round,pad=0.18", fc="white",
                          ec="#cc0000", lw=0.8))

    # annotation boxes around the perimeter
    n = len(top_idx)
    if n == 0: return
    # box positions: distribute around the panel in clockwise order starting top
    if n <= 3:
        angles = np.linspace(45, 405, n+1)[:-1] * np.pi/180   # NE clockwise to NW
    else:
        angles = np.linspace(20, 380, n+1)[:-1] * np.pi/180
    xmin, xmax, ymin, ymax = panel_xy_range
    cx, cy = (xmin+xmax)/2, (ymin+ymax)/2
    rx, ry = (xmax-xmin)*0.55, (ymax-ymin)*0.55
    for k_idx, ci in enumerate(top_idx):
        a = angles[k_idx]
        # box anchor
        bx, by = cx + rx*np.cos(a), cy + ry*np.sin(a)
        d = info[k_idx]
        color = PALETTE[k_idx % len(PALETTE)]
        text = (f"{d['label']}\n"
                f"({d['size']} Members)\n"
                f"({d['papers']} Papers, {int(d['active']*100)}% Active)")
        ha = "center"
        ax.text(bx, by, text, ha=ha, va="center", fontsize=8.5,
                family="serif", linespacing=1.18,
                bbox=dict(boxstyle="round,pad=0.32",
                          fc="white", ec=color, lw=1.4),
                zorder=10)
        # arrow from box to centroid
        ctx, cty = centroids[ci]
        ax.annotate("", xy=(ctx, cty), xytext=(bx, by),
                    arrowprops=dict(arrowstyle="-",
                        color=color, lw=0.8, alpha=0.75,
                        connectionstyle="arc3,rad=0.0"),
                    zorder=4)

fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4))
xs_all = np.array([pos[v][0] for v in hop1])
ys_all = np.array([pos[v][1] for v in hop1])
xr = (xs_all.min()-0.05, xs_all.max()+0.05)
yr = (ys_all.min()-0.05, ys_all.max()+0.05)
panel_xy_range = (xr[0], xr[1], yr[0], yr[1])
for ax, s in zip(axes, PANELS):
    ax.set_xlim(*xr); ax.set_ylim(*yr)
    annotate_modules(ax, s, panel_xy_range)
    sep = results[s]["avg_sep"]; cond = results[s]["avg_cond"]
    ax.set_title(f"$s={s}$:   {results[s]['n_alive']} alive,   "
                 f"{results[s]['n_modules']} modules    "
                 f"(avg $m_{{\\rm in}}/m_{{\\rm cut}}={sep:.2f}$,   $\\phi={cond:.2f}$)",
                 fontsize=11)
    ax.axis("off")

fig.tight_layout()
fig.savefig(DIR/"cs9_egonet.pdf", bbox_inches="tight")
fig.savefig(DIR/"cs9_egonet.png", bbox_inches="tight", dpi=180)
print(f"\nwrote cs9_egonet.pdf / .png [{int(time.time()-t0)}s]")
