"""
Figure 2: Immutable-CPI peeling trace and per-leaf integer state.

Panel (a): Peeling state matrix --- 10 rounds (rows) by 6 leaves (cols).
           Each cell shows the case (A/C/C'/skip) and the leaf's pivot
           count p_P after the round. Color-coded by case.
Panel (b): nCr delta zoom for Round 5 (peel v2 from L3, Case C):
           |V_p| 3 -> 2, supports of v1, v3, v4 fall by Dh = 2, Dp = 1.

Trace is the deterministic output of compute_cpi.simulate_peeling on the
running example. We reuse compute_cpi to keep the figure in lockstep with
the algorithm description in the paper.
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

# ---- Figure layout ----
fig = plt.figure(figsize=(15, 4.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.4, 1.0], wspace=0.20)

# ============ Panel (a): State matrix ============
ax_a = fig.add_subplot(gs[0, 0])
ax_a.set_xlim(-0.5, L + 0.5); ax_a.set_ylim(R + 0.5, -1.5)
ax_a.set_xticks([]); ax_a.set_yticks([])
ax_a.set_title(r"(a) Peeling trace: 10 rounds $\times$ 6 \cpi leaves "
               r"(running example, $s{=}3$)", fontsize=11)
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
            bbox_to_anchor=(0.5, -0.12), ncol=4, fontsize=8.5, frameon=False)

# ============ Panel (b): Delta zoom (Round 5) ============
ax_b = fig.add_subplot(gs[0, 1])
ax_b.set_xlim(0, 10); ax_b.set_ylim(0, 10)
ax_b.axis('off')
ax_b.set_title(r"(b) Delta zoom on Round 5: peel $v_2$ from $L_3$",
               fontsize=11)

# Top: leaf state before
ax_b.text(0.3, 9.0, r"Before:", fontsize=10, fontweight='bold')
ax_b.text(0.3, 8.3,
          r"$L_3 = (V_h{=}\{v_1\},\ V_p{=}\{v_2,v_3,v_4\}),\ \eta{=}2,\ p{=}3$",
          fontsize=9.5)
ax_b.text(0.3, 7.6,
          r"encodes $\binom{3}{2}{=}3$ triangles "
          r"$\{v_1{,}v_2{,}v_3\}, \{v_1{,}v_2{,}v_4\}, \{v_1{,}v_3{,}v_4\}$",
          fontsize=9)

# Arrow
ax_b.annotate('', xy=(5.0, 6.6), xytext=(5.0, 7.2),
              arrowprops=dict(arrowstyle='->', lw=1.4, color='#3a3a3a'))
ax_b.text(5.2, 6.85, r"peel $v_2$ (a pivot of $L_3$, Case C)",
          fontsize=9, color='#3a3a3a')

# After
ax_b.text(0.3, 6.0, r"After:", fontsize=10, fontweight='bold')
ax_b.text(0.3, 5.3,
          r"$L_3$ \emph{unchanged in memory}; counter $p{:}3{\to}2$, "
          r"liveness still alive",
          fontsize=9.5)
ax_b.text(0.3, 4.6,
          r"effective encoding shrinks to $\binom{2}{2}{=}1$ triangle "
          r"$\{v_1,v_3,v_4\}$",
          fontsize=9)

# Delta box
ax_b.add_patch(mpatches.FancyBboxPatch(
    (0.3, 1.6), 9.4, 2.6, boxstyle="round,pad=0.10",
    linewidth=1.0, edgecolor='#3a8a3a', facecolor='#eaf6ea'))
ax_b.text(0.6, 3.6, r"\textbf{Counter-based delta} (Lemma~2):",
          fontsize=10, fontweight='bold', color='#1f5a1f')
ax_b.text(0.6, 2.9,
          r"$\Delta_{\mathrm{hold}} = \binom{p}{\eta} - \binom{p{-}d}{\eta}$"
          r"$= \binom{3}{2}-\binom{2}{2} = 3-1 = 2$",
          fontsize=10)
ax_b.text(0.6, 2.2,
          r"$\Delta_{\mathrm{pivot}} = \binom{p{-}1}{\eta{-}1} - \binom{p{-}d{-}1}{\eta{-}1}$"
          r"$= \binom{2}{1}-\binom{1}{1} = 2-1 = 1$",
          fontsize=10)

# Support updates
ax_b.text(0.3, 1.0, r"\textbf{Support updates:} "
          r"$v_1$ (hold) $-2$;\ \ $v_3,v_4$ (pivots) $-1$ each",
          fontsize=9.5)
ax_b.text(0.3, 0.3,
          r"\emph{No tree mutation, no hash op, $O(|V_h|+|V_p|{-}1)$ work}.",
          fontsize=9, style='italic', color='#1f5a1f')

fig.suptitle(r"Figure 2: Peeling trace and counter-based delta on the running example.",
             fontsize=11.5, y=1.02)
fig.tight_layout()
fig.savefig(OUT / "fig2_peeling.pdf", bbox_inches='tight')
fig.savefig(OUT / "fig2_peeling.png", dpi=170, bbox_inches='tight')
print(f"Saved fig2_peeling to {OUT}")
