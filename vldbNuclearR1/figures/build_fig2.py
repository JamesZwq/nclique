"""
Figure 2: Immutable-CPI peeling trace (state matrix only).

The 10 rounds * 6 leaves matrix is the only genuinely visual content:
each cell shows the case (A / C / C' / skip) at that (round, leaf), and
the per-row colour pattern lets the reader scan how cases distribute
over time.  The Round-5 nCr delta worked example (formerly panel b)
lives inline in the body prose where it belongs.

Trace is the deterministic output of compute_cpi.simulate_peeling on the
running example.
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

# ---- Replay peeling trace ----
adj = build_graph()
s = 3
leaves, _ = sdct_cpi(adj, s)
trace, core = simulate_peeling(adj, leaves, s, verbose=False)

L = len(leaves)
R = len(trace)

# Per-(round, leaf) state matrix
case_grid     = [[None] * L for _ in range(R)]
delta_grid    = [[None] * L for _ in range(R)]
pivot_grid    = [[None] * L for _ in range(R)]

# Track pivot count over time
rp_now = [len(leaves[i][1]) for i in range(L)]
alive  = [True] * L

for r_idx, row in enumerate(trace):
    # Apply events to update rp_now
    for e in row['events']:
        i = e['leaf'] - 1
        if e.get('status') == 'skip (dead)':
            case_grid[r_idx][i] = 'skip'
        else:
            case_grid[r_idx][i]  = e['case'][0:2]  # 'A ', 'C ', "C'"
            delta_grid[r_idx][i] = (e.get('Dh', 0), e.get('Dp', 0))
            if e['case'].startswith('C') and 'p' in e:
                # parse '3->2'
                new_p = int(e['p'].split('\u2192')[1])
                rp_now[i] = new_p
            if e.get('path_dies'):
                alive[i] = False
    # snapshot pivot counts at end of round
    for i in range(L):
        pivot_grid[r_idx][i] = rp_now[i] if alive[i] or case_grid[r_idx][i] else rp_now[i]

# ---- Figure layout: state matrix only, half-column width ----
fig, ax_a = plt.subplots(figsize=(3.6, 3.2))
ax_a.set_xlim(-0.5, L + 0.5); ax_a.set_ylim(R + 0.5, -1.5)
ax_a.set_xticks([]); ax_a.set_yticks([])
ax_a.axis('off')

# Color palette
COL = {
    'A':    '#d94545',  # hold removed -> path dies
    'C ':   '#3a8a3a',  # pivot shrink, alive
    "C'":   '#e58a2a',  # pivot -> death
    'skip': '#cfcfcf',  # already dead
    None:   '#fafafa',  # leaf not touched
}

# Header row: round | popped | core | L1..L6
hdr_y = -0.6
ax_a.text(-0.45, hdr_y, "round", fontsize=9, fontweight='bold')
ax_a.text(0.55,  hdr_y, "pop",   fontsize=9, fontweight='bold')
ax_a.text(1.45,  hdr_y, r"$\kappa$", fontsize=9, fontweight='bold')
for i in range(L):
    ax_a.text(2.2 + i, hdr_y, f"$L_{i+1}$", fontsize=9, fontweight='bold',
              ha='center')

cell_w, cell_h = 0.92, 0.78
for r_idx, row in enumerate(trace):
    y = r_idx
    ax_a.text(-0.45, y, f"{row['round']}", fontsize=8.5)
    ax_a.text( 0.55, y, f"$v_{{{row['victim']}}}$", fontsize=8.5)
    ax_a.text( 1.45, y, f"{row['core']}", fontsize=8.5, ha='left')
    for i in range(L):
        c = case_grid[r_idx][i]
        if c is None:
            face = COL[None]
            text = ''
        else:
            face = COL.get(c, '#ffffff')
            if c == 'skip':
                text = '\u2014'
            elif c == 'A':
                text = 'A'
            elif c == 'C ':
                text = 'C'
            elif c == "C'":
                text = "C\u2032"
            else:
                text = c
        rect = mpatches.FancyBboxPatch(
            (2.2 + i - cell_w/2, y - cell_h/2 + 0.05),
            cell_w, cell_h, boxstyle="round,pad=0.02",
            linewidth=0.5, edgecolor='#888', facecolor=face)
        ax_a.add_patch(rect)
        # case letter (top-left)
        ax_a.text(2.2 + i - 0.32, y + 0.10, text, fontsize=8.5,
                  fontweight='bold', color='white' if c not in (None, 'skip') else '#555')
        # delta numbers (right, only for live cases)
        d = delta_grid[r_idx][i]
        if d is not None and c in ('A', 'C ', "C'"):
            ax_a.text(2.2 + i + 0.30, y + 0.10,
                      f"$\\Delta_h{{=}}{d[0]}$", fontsize=7.0, color='white')
            ax_a.text(2.2 + i + 0.30, y - 0.18,
                      f"$\\Delta_p{{=}}{d[1]}$", fontsize=7.0, color='white')

# Legend
legend_handles = [
    mpatches.Patch(facecolor=COL['A'],   edgecolor='#888', label='Case A: hold removed (dies)'),
    mpatches.Patch(facecolor=COL['C '],  edgecolor='#888', label='Case C: pivot shrinks (alive)'),
    mpatches.Patch(facecolor=COL["C'"],  edgecolor='#888', label="Case C$'$: pivot $\\to$ death"),
    mpatches.Patch(facecolor=COL['skip'],edgecolor='#888', label='already dead (skip)'),
]
ax_a.legend(handles=legend_handles, loc='lower center',
            bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=7.5, frameon=False)

fig.savefig(OUT / "fig2_peeling.pdf", bbox_inches='tight')
fig.savefig(OUT / "fig2_peeling.png", dpi=180, bbox_inches='tight')
print(f"Saved fig2_peeling to {OUT}")
