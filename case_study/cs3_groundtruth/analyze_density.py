"""
CS3 density analysis: induced subgraph density on top-K vertices.
Compares (1,s)-core vs k-core on intrinsic denseness (edge + clique density),
independent of ground-truth community labels.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams["pdf.fonttype"]=42
matplotlib.rcParams["ps.fonttype"]=42
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent
S_VALUES = [2, 3, 5, 10, 15, 20, 30]
N = 425957

print("Loading graph adj...")
adj = [set() for _ in range(N)]
with open(ROOT / "com-dblp.ungraph.txt") as f:
    for line in f:
        if line.startswith("#"): continue
        u, v = line.split()
        u, v = int(u), int(v)
        if u != v:
            adj[u].add(v); adj[v].add(u)

cores = {}
for s in S_VALUES:
    df = pd.read_csv(ROOT / f"com-dblp-snap_def_s{s}.tsv", sep="\t", comment="#",
                     names=["vid", "core"], dtype={"vid": np.int64, "core": np.float64})
    cv = np.zeros(N)
    cv[df["vid"].values] = df["core"].values
    cores[s] = cv
kc = np.load(ROOT / "kcore.npy")

def induced_metrics(top_set):
    top_set = set(top_set)
    edges = 0
    for v in top_set:
        edges += len(adj[v] & top_set)
    edges //= 2
    n = len(top_set)
    density = edges / (n * (n - 1) / 2) if n > 1 else 0
    triangles = 0
    for v in top_set:
        nbrs_in = adj[v] & top_set
        for u in nbrs_in:
            if u > v:
                triangles += len(adj[u] & nbrs_in)
    triangles //= 3
    return edges, density, triangles

Ks = [100, 500, 1000, 2000, 5000, 10000]
results = []
for K in Ks:
    for s in S_VALUES:
        cv = cores[s]
        nz = int((cv > 0).sum())
        k_eff = min(K, nz)
        top = [int(i) for i in np.argsort(cv)[-k_eff:]]
        e, d, t = induced_metrics(top)
        results.append({"method": f"(1,{s})-core", "K": K, "top_size": k_eff,
                        "edges": e, "density": d, "triangles": t,
                        "tri_per_v": t / k_eff if k_eff else 0})
    top_kc = [int(i) for i in np.argsort(kc)[-K:]]
    e, d, t = induced_metrics(top_kc)
    results.append({"method": "k-core", "K": K, "top_size": K,
                    "edges": e, "density": d, "triangles": t,
                    "tri_per_v": t / K})

df = pd.DataFrame(results)
df.to_csv(ROOT / "cs3_density.csv", index=False)
print("\n=== Induced subgraph density ===")
for K in Ks:
    print(f"\n--- K = {K} ---")
    sub = df[df["K"] == K]
    print(sub.to_string(index=False))

# Plot: density vs K for each method
fig, axes = plt.subplots(1, 2, figsize=(13, 4.2))
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(S_VALUES)))
for ax, col, ylabel, title in [
    (axes[0], "density", "Edge density", "CS3(c): Induced edge density vs K"),
    (axes[1], "tri_per_v", "Triangles per vertex", "CS3(d): Triangle density vs K"),
]:
    for s, c in zip(S_VALUES, colors):
        sub = df[df["method"] == f"(1,{s})-core"].sort_values("K")
        ax.plot(sub["K"], sub[col], marker="o", color=c, label=f"s={s}", markersize=4, linewidth=1.2)
    sub_kc = df[df["method"] == "k-core"].sort_values("K")
    ax.plot(sub_kc["K"], sub_kc[col], marker="s", color="red", label="k-core", markersize=6, linewidth=1.8, linestyle="--")
    ax.set_xscale("log")
    if col == "tri_per_v":
        ax.set_yscale("log")
    ax.set_xlabel("K (top ranking)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3, linestyle=":")
    ax.legend(ncol=2, fontsize=8, loc="best")
fig.tight_layout()
fig.savefig(ROOT / "cs3_density.png", dpi=150)
fig.savefig(ROOT / "cs3_density.pdf")
print(f"\nSaved plots to {ROOT}")
