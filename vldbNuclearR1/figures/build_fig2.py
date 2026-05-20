"""
Figure 2: SPIN* peel trace on the running example, as leaf lifelines.

Each leaf gets one horizontal lifeline showing its state across the 10
iterations. Iteration columns at the top carry the popped vertex and
the kappa value. Each cell in a lifeline is either:
  - blank (the iteration's popped vertex is not on this leaf)
  - a coloured tick marking which case fired:
       C  (pivot shrink, alive)        -> green
       C' (pivot -> death)             -> orange
       A  (hold removed -> path dies)  -> red
  - a grey strip after the leaf has died

The leaf's pivot-count trajectory is annotated to the right.

Trace is the deterministic output of compute_cpi.simulate_peeling on
the running example.
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType -- VLDB/ACM forbid Type 3
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['text.usetex']  = False
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from compute_cpi import build_graph, sdct_cpi, simulate_peeling

OUT = Path(__file__).parent

adj = build_graph()
s = 3
leaves, _ = sdct_cpi(adj, s)
trace, core = simulate_peeling(adj, leaves, s, verbose=False)

L = len(leaves)
R = len(trace)

# case_grid[i][r] = case label on leaf i at iteration r, or None
# 'A' / 'C' / "C'" / 'skip'
case_grid = [[None] * R for _ in range(L)]
for r_idx, row in enumerate(trace):
    for e in row['events']:
        i = e['leaf'] - 1
        if e.get('status') == 'skip (dead)':
            case_grid[i][r_idx] = 'skip'
        else:
            tag = e['case'].split(' ')[0]   # "A" / "C" / "C'"
            case_grid[i][r_idx] = tag

# Lifeline metadata: pivot trajectory (initial pivot count -> per-step)
init_p = [len(p) for h, p in leaves]
eta    = [s - len(h) for h, p in leaves]

# ---- Layout ----
# 6 lifelines, each with R=10 cells of width 1.
# Header height: 2 rows (pop, kappa); top of plot.
fig, ax = plt.subplots(figsize=(7.0, 2.6))

cell_w = 1.0
cell_h = 0.62
gap_y  = 0.25

# Colours
RED    = '#d94545'   # Case A
GREEN  = '#3a8a3a'   # Case C (pivot shrink)
ORANGE = '#e58a2a'   # Case C' (pivot -> death)
GREY   = '#d0d0d0'   # already-dead strip
LIVE   = '#fafafa'   # alive but no event this iter
HOLDV  = '#fff3d6'   # hold-vertex band (subtle)

ax.set_xlim(-2.4, R + 0.4)
# total height: header (~1.4) + 6 lifelines (each cell_h + gap_y)
TOP = 1.4
total_h = TOP + L * (cell_h + gap_y)
ax.set_ylim(total_h, -0.2)
ax.axis('off')

# --- Header: iteration columns ---
for r_idx, row in enumerate(trace):
    x = r_idx + 0.5
    # iter number
    ax.text(x, 0.05, f"{row['round']}", ha='center', va='center',
            fontsize=7.5, color='#666')
    # popped vertex
    ax.text(x, 0.45, f"$v_{{{row['victim']}}}$", ha='center', va='center',
            fontsize=8.5, fontweight='bold')
    # kappa
    ax.text(x, 0.90, f"{row['core']}", ha='center', va='center',
            fontsize=8.0, color='#1f3a8a')

# Header row labels
ax.text(-0.05, 0.05, "iter", ha='right', va='center',
        fontsize=7.5, fontweight='bold', color='#666')
ax.text(-0.05, 0.45, "pop",  ha='right', va='center',
        fontsize=8.0, fontweight='bold')
ax.text(-0.05, 0.90, r"$\kappa_s$", ha='right', va='center',
        fontsize=8.0, fontweight='bold', color='#1f3a8a')

# Divider line under header
ax.plot([-2.2, R + 0.2], [TOP - 0.15, TOP - 0.15],
        linewidth=0.6, color='#999')

# --- Leaf lifelines ---
for i in range(L):
    h, p = leaves[i]
    y0 = TOP + i * (cell_h + gap_y)
    y_mid = y0 + cell_h / 2

    # Left label: leaf id + composition
    hold_str  = ','.join(f"v_{{{v}}}" for v in h)
    pivot_str = ','.join(f"v_{{{v}}}" for v in p)
    label = (rf"$L_{i+1}$: $V_h\!=\!\{{{hold_str}\}}$, "
             rf"$V_p\!=\!\{{{pivot_str}\}}$, $\eta\!=\!{eta[i]}$")
    ax.text(-0.05, y_mid, label, ha='right', va='center', fontsize=7.4)

    # Has the leaf died yet?
    died_at = None
    for r_idx in range(R):
        c = case_grid[i][r_idx]
        if c in ('A', "C'"):
            died_at = r_idx
            break

    # Background strip: alive (light) until death, grey after
    if died_at is None:
        ax.add_patch(mpatches.Rectangle((0, y0), R, cell_h,
                     facecolor=LIVE, edgecolor='#bbb', linewidth=0.3))
    else:
        ax.add_patch(mpatches.Rectangle((0, y0), died_at + 1, cell_h,
                     facecolor=LIVE, edgecolor='#bbb', linewidth=0.3))
        if died_at + 1 < R:
            ax.add_patch(mpatches.Rectangle((died_at + 1, y0),
                         R - died_at - 1, cell_h,
                         facecolor=GREY, edgecolor='#bbb', linewidth=0.3,
                         hatch='///', alpha=0.55))

    # Event markers
    for r_idx in range(R):
        c = case_grid[i][r_idx]
        x = r_idx + 0.5
        if c == 'C':
            ax.add_patch(mpatches.Circle((x, y_mid), 0.20,
                         facecolor=GREEN, edgecolor='white', linewidth=0.6))
            ax.text(x, y_mid, 'C', ha='center', va='center',
                    fontsize=7.5, fontweight='bold', color='white')
        elif c == "C'":
            ax.add_patch(mpatches.FancyBboxPatch(
                (x - 0.28, y_mid - 0.22), 0.56, 0.44,
                boxstyle="round,pad=0.01",
                facecolor=ORANGE, edgecolor='white', linewidth=0.6))
            ax.text(x, y_mid, "C$'$", ha='center', va='center',
                    fontsize=7.5, fontweight='bold', color='white')
        elif c == 'A':
            ax.add_patch(mpatches.FancyBboxPatch(
                (x - 0.28, y_mid - 0.22), 0.56, 0.44,
                boxstyle="round,pad=0.01",
                facecolor=RED, edgecolor='white', linewidth=0.6))
            ax.text(x, y_mid, 'A', ha='center', va='center',
                    fontsize=7.5, fontweight='bold', color='white')
        # 'skip' cells render as part of the grey post-death strip; no marker
        # None (vertex not on this leaf) renders as background

# --- Legend ---
legend_handles = [
    mpatches.Patch(facecolor=GREEN,  edgecolor='white',
                   label="C: pivot shrinks (leaf alive)"),
    mpatches.Patch(facecolor=ORANGE, edgecolor='white',
                   label=r"C$'$: pivot $\to$ death"),
    mpatches.Patch(facecolor=RED,    edgecolor='white',
                   label='A: hold removed (death)'),
    mpatches.Patch(facecolor=GREY,   edgecolor='#bbb',
                   hatch='///', alpha=0.55,
                   label='leaf already dead'),
]
ax.legend(handles=legend_handles, loc='lower center',
          bbox_to_anchor=(0.42, -0.14), ncol=4, fontsize=7.2,
          frameon=False, handlelength=1.2, columnspacing=1.0)

fig.savefig(OUT / "fig2_peeling.pdf", bbox_inches='tight')
fig.savefig(OUT / "fig2_peeling.png", dpi=180, bbox_inches='tight')
print(f"Saved fig2_peeling to {OUT}")
