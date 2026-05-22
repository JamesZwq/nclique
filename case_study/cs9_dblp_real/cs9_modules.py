#!/usr/bin/env python3
"""Case Study 4: Jiawei Han 1-hop ego, nucleus-module decomposition under r=1.

Mirrors NuclearCD CS-II's framing, adapted for r=1:
  At r=1 the (1,s)-nucleus alive subgraph is naturally one giant connected
  component (every coauthor of Han is reachable through a common collaborator).
  We therefore decompose the alive subgraph into disjoint dense modules
  via greedy near-clique extraction (size >= MIN_MOD, density >= DENS) on
  the anchor-removed alive subgraph; these modules play the role of
  "nucleus components" in NuclearCD CS-II.

For each s in {3, 10}, we report:
  - |alive|, # modules
  - intra-module edge fraction = edges within modules / edges in alive subgraph
  - size-weighted avg separability m_in/m_cut  (m_cut = module->other-alive edges)
  - size-weighted avg conductance phi(C) = m_cut / (2 m_in + m_cut)
  - top-K modules: per-module sep, cond, active ratio
       active ratio = fraction of members in >=1 DBLP paper with >=3
       within-module coauthors
  - top-K avg active ratio (size-weighted and unweighted)

Saves cs9_modules.json + cs9_egonet.pdf (2 panels coloured by module).
"""
import json, time
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
TOPK    = {S_LOW: 5, S_HIGH: 8}
MIN_MOD = 10           # min near-clique size
DENS    = 0.85         # min internal edge density
MAX_MOD = 200          # safety cap

# ---- load ----
t0 = time.time()
id2name = {}
with (DIR/"authors.tsv").open() as f:
    for line in f:
        i, n = line.rstrip("\n").split("\t", 1); id2name[int(i)] = n
print(f"authors loaded ({len(id2name)}) [{int(time.time()-t0)}s]")

adj = defaultdict(set)
with (DIR/"edges.txt").open() as f:
    for line in f:
        u, v = line.split(); u, v = int(u), int(v); adj[u].add(v); adj[v].add(u)
print(f"adj loaded [{int(time.time()-t0)}s]")

hop1 = adj[ANCHOR] | {ANCHOR}
print(f"|1-hop ego| = {len(hop1)}")

core = {s: {} for s in PANELS}
for s in PANELS:
    with (DIR/f"r1_s{s}.tsv").open() as f:
        for line in f:
            if line.startswith("#"): continue
            i, c = line.rstrip().split(); i = int(i)
            if i in hop1: core[s][i] = int(c)

# "alive" = participates in at least one (1,s)-nucleus (core >= 1)
alive = {s: {v for v in hop1 if core[s].get(v, 0) >= 1} for s in PANELS}
print({s: len(alive[s]) for s in PANELS})

adj_hop1 = {v: (adj[v] & hop1) for v in hop1}

# ---- papers ----
hop1_names = {id2name[v]: v for v in hop1 if v in id2name}
papers = []
with (DIR/"papers.tsv").open() as f:
    for line in f:
        names = [a for a in line.rstrip("\n").split("\t") if a]
        mem = {hop1_names[a] for a in names if a in hop1_names}
        if len(mem) >= 3:
            papers.append(mem)
print(f"papers >=3 hop1-authors: {len(papers)} [{int(time.time()-t0)}s]")
papers_of = defaultdict(list)
for pi, pa in enumerate(papers):
    for v in pa: papers_of[v].append(pi)

# ---- greedy near-clique extractor ----
def find_modules(pool, anchor_excluded=True, density=DENS, min_size=MIN_MOD,
                 adj_local=None, max_n=MAX_MOD):
    pool = set(pool)
    if anchor_excluded: pool.discard(ANCHOR)
    if adj_local is None: adj_local = adj_hop1
    found = []
    while len(pool) >= min_size and len(found) < max_n:
        deg = {x: len(adj_local[x] & pool) for x in pool}
        if not deg: break
        seed = max(deg, key=deg.get)
        if deg[seed] < min_size - 1: pool.discard(seed); continue
        clique = {seed}
        cand = sorted(adj_local[seed] & pool, key=lambda x: -len(adj_local[x] & pool))
        for c in cand:
            if c in clique: continue
            if len(adj_local[c] & clique) < density * len(clique): continue
            clique.add(c)
        if len(clique) >= min_size:
            e = sum(1 for u in clique for w in clique
                    if u < w and w in adj_local[u])
            p = len(clique)*(len(clique)-1)//2
            if p and e/p >= density:
                found.append(clique); pool -= clique
            else: pool.discard(seed)
        else: pool.discard(seed)
    return found

# ---- analysis per s ----
results = {}
panel_data = {}
for s in PANELS:
    A = alive[s]
    mods = find_modules(A)
    print(f"\ns={s}: alive={len(A)} modules={len(mods)} sizes={[len(c) for c in mods]}")

    # global edges in alive subgraph
    edges_alive = 0
    for u in A:
        edges_alive += sum(1 for v in adj_hop1[u] if v in A and u < v)

    # module assignment
    mod_id = {}
    for ci, M in enumerate(mods):
        for v in M: mod_id[v] = ci

    sizes = [len(M) for M in mods]
    m_in  = [0]*len(mods)
    m_cut = [0]*len(mods)    # edges from module to other-alive (incl. other modules + non-module alive + anchor)
    for u in A:
        if u not in mod_id: continue
        ci = mod_id[u]
        for v in adj_hop1[u]:
            if v not in A or v == u: continue
            if v in mod_id and mod_id[v] == ci:
                if u < v: m_in[ci] += 1
            else:
                m_cut[ci] += 1
    total_intra = sum(m_in)
    intra_frac  = total_intra / max(edges_alive, 1)

    sep  = [(m_in[i] / m_cut[i]) if m_cut[i] > 0 else float("inf") for i in range(len(mods))]
    cond = [(m_cut[i] / (2*m_in[i] + m_cut[i])) if (2*m_in[i] + m_cut[i]) > 0 else 0.0
            for i in range(len(mods))]
    SEP_CAP = 20.0
    sep_cap = [min(x, SEP_CAP) for x in sep]
    total_sz = sum(sizes) if sizes else 1
    avg_sep  = (sum(sep_cap[i]*sizes[i] for i in range(len(mods)))/total_sz) if mods else 0.0
    avg_cond = (sum(cond[i]*sizes[i]    for i in range(len(mods)))/total_sz) if mods else 0.0

    # active ratio per module: members in any >=3-within-module paper
    active = []
    n_collab = []
    for ci, M in enumerate(mods):
        cands = set()
        for v in M: cands.update(papers_of.get(v, []))
        actives = set(); cnt = 0
        for pi in cands:
            inter = papers[pi] & M
            if len(inter) >= 3:
                cnt += 1; actives |= inter
        active.append(len(actives)/len(M) if M else 0.0)
        n_collab.append(cnt)

    k = TOPK[s]
    top_idx = sorted(range(len(mods)), key=lambda i: -sizes[i])[:k]
    if top_idx:
        topk_act_w   = sum(active[i]*sizes[i] for i in top_idx) / sum(sizes[i] for i in top_idx)
        topk_act_unw = sum(active[i] for i in top_idx) / len(top_idx)
    else:
        topk_act_w = topk_act_unw = 0.0

    def reps(M, n=4):
        deg = [(v, len(adj_hop1[v] & M)) for v in M]
        deg.sort(key=lambda kv: -kv[1])
        return [id2name.get(v, f"id={v}") for v, _ in deg[:n]]

    top_info = []
    for i in top_idx:
        top_info.append({
            "module_id": i, "size": sizes[i],
            "m_in": m_in[i], "m_cut": m_cut[i],
            "sep":  None if sep[i] == float("inf") else round(sep[i], 3),
            "cond": round(cond[i], 3),
            "active": round(active[i], 3),
            "collab_papers": n_collab[i],
            "reps": reps(mods[i], 5),
        })

    results[s] = {
        "n_alive": len(A),
        "n_modules": len(mods),
        "edges_alive": edges_alive,
        "intra_frac": round(intra_frac, 3),
        "avg_sep":  round(avg_sep, 3),
        "avg_cond": round(avg_cond, 3),
        "topK_avg_active":     round(topk_act_w, 3),
        "topK_avg_active_unw": round(topk_act_unw, 3),
        "top": top_info,
    }
    panel_data[s] = (mods, sizes, sep, cond, active, top_idx, mod_id)

    print(f"  alive={len(A)}  modules={len(mods)}  edges_alive={edges_alive}")
    print(f"  intra_frac={intra_frac:.3f}  avg_sep={avg_sep:.3f}  avg_cond={avg_cond:.3f}")
    print(f"  top-{k} avg active (size-w) = {topk_act_w:.3f}   (unw) = {topk_act_unw:.3f}")
    for info in top_info:
        sep_s = "inf" if info["sep"] is None else f"{info['sep']:.2f}"
        print(f"   M{info['module_id']:>2d}: |C|={info['size']:>4d} m_in={info['m_in']:>4d} "
              f"m_cut={info['m_cut']:>4d} sep={sep_s:>6s} cond={info['cond']:.2f} "
              f"act={info['active']:.2f} collabP={info['collab_papers']:>3d}")
        print(f"        reps: {', '.join(info['reps'])}")

(DIR/"cs9_modules.json").write_text(json.dumps(results, indent=2))
print(f"\nwrote cs9_modules.json [{int(time.time()-t0)}s]")

# ---- render figure: 2 panels ----
print("building layout...")
G = nx.Graph(); G.add_nodes_from(hop1)
for u in hop1:
    for v in adj_hop1[u]:
        if u < v: G.add_edge(u, v)
pos = nx.spring_layout(G, seed=42, iterations=120, k=1.6/(len(hop1)**0.5))
print(" done")

PALETTE = ["#1f78b4", "#e31a1c", "#33a02c", "#ff7f00", "#6a3d9a",
           "#b15928", "#a6cee3", "#fb9a99"]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
for ax, s in zip(axes, PANELS):
    mods, sizes, sep, cond, active, top_idx, mod_id = panel_data[s]
    A = alive[s]
    peeled = hop1 - A

    nx.draw_networkx_nodes(G, pos, nodelist=peeled,
        node_color="#e6e6e6", node_size=3, alpha=0.55, ax=ax, linewidths=0)
    top_set = set().union(*(mods[i] for i in top_idx)) if top_idx else set()
    others = (A - top_set) - {ANCHOR}
    nx.draw_networkx_nodes(G, pos, nodelist=others,
        node_color="#9a9a9a", node_size=8, alpha=0.7, ax=ax, linewidths=0)
    for k_idx, ci in enumerate(top_idx):
        color = PALETTE[k_idx % len(PALETTE)]
        nodes = list(mods[ci] - {ANCHOR})
        nx.draw_networkx_nodes(G, pos, nodelist=nodes,
            node_color=color, node_size=24, alpha=0.95, ax=ax,
            linewidths=0.3, edgecolors="white")
    if ANCHOR in A:
        ax.scatter([pos[ANCHOR][0]], [pos[ANCHOR][1]], marker="*",
                   s=180, c="black", zorder=6, linewidths=0)

    ax.set_title(f"$s={s}$:  {len(A)} alive,  {results[s]['n_modules']} modules,  "
                 f"top-{len(top_idx)} coloured",
                 fontsize=10)
    ax.axis("off")

max_k = max(TOPK.values())
leg_handles = [mpatches.Patch(color=PALETTE[i % len(PALETTE)], label=f"rank {i+1}")
               for i in range(max_k)]
leg_handles += [mpatches.Patch(color="#9a9a9a", label="other alive"),
                mpatches.Patch(color="#e6e6e6", label="peeled")]
fig.legend(handles=leg_handles, loc="lower center", ncol=max_k+2,
           fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig(DIR/"cs9_egonet.pdf", bbox_inches="tight")
fig.savefig(DIR/"cs9_egonet.png", bbox_inches="tight", dpi=180)
print("wrote cs9_egonet.pdf / .png")
