"""
Figure 2: SPIN* peel trace on the running example as a vertex h-value
trajectory across the ten deletion iterations.

Layout: 10 rows (one per vertex v_1..v_10), columns for iter 0
(initial sdeg) plus iters 1..10 (state after each pop). Each cell:
  - if v is still alive at iter t: residual support (numeric), shaded
    by value;
  - the cell where v is popped: highlighted, labelled "kappa=K";
  - cells after v is popped: blank.

Header row gives the popped vertex and the assigned kappa at each iter,
so the figure also shows the peel order.

Parallel in layout to the SPIN figure for Example 5.2: same 10 rows,
same cell style, just more columns because SPIN* finalises one vertex
per iteration.

The trajectory data is regenerated from compute_cpi.simulate_peeling
on every run, so it stays in sync with the running example.
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['text.usetex']  = False
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from math import comb
from pathlib import Path
from compute_cpi import build_graph, sdct_cpi, simulate_peeling

OUT = Path(__file__).parent

adj = build_graph()
s = 3
leaves, _ = sdct_cpi(adj, s)
trace, _ = simulate_peeling(adj, leaves, s, verbose=False)

N = 10
T = len(trace)            # 10 iterations

# Initial sdeg
h0 = [0] * N
for (H, P) in leaves:
    eta = s - len(H)
    for u in H:
        h0[u-1] += comb(len(P), eta)
    for u in P:
        if eta >= 1:
            h0[u-1] += comb(len(P)-1, eta-1)

# Build per-iter snapshot:
#   snap[t][i] is either an integer (still alive), the string
#   'pop:K' for the iteration where v_{i+1} is popped, or None
#   for iterations after the pop.
popped_at = {}   # vertex -> (iter_round, kappa)
for row in trace:
    popped_at[row['victim']] = (row['round'], row['core'])

snap = [[None] * N for _ in range(T + 1)]   # +1 for iter 0
snap[0] = h0[:]

for row in trace:
    t = row['round']
    s_after = row['supp_after']
    for v in range(1, N + 1):
        if v in popped_at and popped_at[v][0] <= t:
            if popped_at[v][0] == t:
                snap[t][v-1] = f"pop:{popped_at[v][1]}"
            else:
                snap[t][v-1] = None
        else:
            snap[t][v-1] = int(s_after.get(v, 0))

# Color ramp (matches SPIN figure)
def color_for(v):
    if v is None or v == 0:
        return '#fafafa'
    if v == 1:    return '#dde8f5'
    if v <= 3:    return '#a9c4e6'
    if v <= 5:    return '#6d96cf'
    return '#3a6ebb'

def text_color_for(v):
    return '#1a1a1a' if v <= 3 else 'white'

POP_FACE = '#ffe9b3'
POP_EDGE = '#1f3a8a'

cell_w = 0.78
cell_h = 0.55
gap_y  = 0.10
left   = 1.55

fig, ax = plt.subplots(figsize=(8.4, 3.9))

TOP = 0.95
total_h = TOP + N * (cell_h + gap_y) + 0.25
ax.set_xlim(-0.3, left + (T + 1) * cell_w + 0.4)
ax.set_ylim(total_h, -0.7)
ax.axis('off')

# Header: vertex column label + iter columns
ax.text(left / 2 - 0.05, -0.30, 'vertex',
        ha='center', va='center', fontsize=9.0, fontweight='bold')

# Per-iter header info: pop / kappa
# iter 0 has no pop event
ax.text(left + 0.5 * cell_w, -0.55, "iter 0",
        ha='center', va='center', fontsize=7.5, color='#666')
ax.text(left + 0.5 * cell_w,  0.10, "(init)",
        ha='center', va='center', fontsize=7.5, color='#666',
        fontstyle='italic')

for row in trace:
    t = row['round']
    x = left + t * cell_w + cell_w / 2
    ax.text(x, -0.55, f"iter {t}",
            ha='center', va='center', fontsize=7.5, color='#666')
    ax.text(x, -0.30, f"pop $v_{{{row['victim']}}}$",
            ha='center', va='center', fontsize=8.4, fontweight='bold')
    ax.text(x,  0.10, rf"$\kappa\!=\!{row['core']}$",
            ha='center', va='center', fontsize=8.0, color='#1f3a8a')

# Divider under header
ax.plot([0.0, left + (T + 1) * cell_w],
        [TOP - 0.10, TOP - 0.10],
        linewidth=0.6, color='#999')

# Rows
for i in range(N):
    y0 = TOP + i * (cell_h + gap_y)
    y_mid = y0 + cell_h / 2

    ax.text(left / 2 - 0.05, y_mid, f"$v_{{{i+1}}}$",
            ha='center', va='center', fontsize=9.5)

    for t in range(T + 1):
        x0 = left + t * cell_w
        cell = snap[t][i]
        if cell is None:
            # blank cell after pop
            continue
        if isinstance(cell, str) and cell.startswith("pop:"):
            kappa = int(cell.split(":")[1])
            ax.add_patch(mpatches.Rectangle(
                (x0, y0), cell_w * 0.92, cell_h,
                facecolor=POP_FACE, edgecolor=POP_EDGE, linewidth=1.2))
            ax.text(x0 + cell_w * 0.46, y_mid, rf"$\kappa\!=\!{kappa}$",
                    ha='center', va='center',
                    fontsize=7.8, fontweight='bold', color=POP_EDGE)
        else:
            v = int(cell)
            ax.add_patch(mpatches.Rectangle(
                (x0, y0), cell_w * 0.92, cell_h,
                facecolor=color_for(v), edgecolor='#bbb', linewidth=0.4))
            ax.text(x0 + cell_w * 0.46, y_mid, f"{v}",
                    ha='center', va='center', fontsize=8.6,
                    fontweight='bold', color=text_color_for(v))

fig.savefig(OUT / "fig2_peeling.pdf", bbox_inches='tight')
fig.savefig(OUT / "fig2_peeling.png", dpi=180, bbox_inches='tight')
print(f"Saved fig2_peeling to {OUT}")
