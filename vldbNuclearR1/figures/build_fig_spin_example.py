"""
Figure for Example 5.2 (SPIN on the running example): vertex h-value
trajectory across full-pass iterations.

Layout: 10 rows (one per vertex v_1..v_10), columns for passes
0 (initial), 1 (after one full pass), 2 (= pass 1, certifies the
fixed point). Each cell shows h^(t)[v]; the converged value equals
kappa_s(v).

This parallels Figure 2 for SPIN*: same rows-are-vertices layout, just
fewer columns because SPIN converges in one full pass on this
instance.

Data:
  h^(0) = sdeg = [3, 4, 6, 6, 3, 4, 1, 1, 1, 1]
  h^(1) = h^(2) = [3, 3, 3, 3, 3, 3, 1, 1, 1, 1]   (= kappa_3)
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['text.usetex']  = False
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

OUT = Path(__file__).parent

# Per-pass h trajectory, verified against compute_cpi.
H = [
    [3, 4, 6, 6, 3, 4, 1, 1, 1, 1],   # pass 0
    [3, 3, 3, 3, 3, 3, 1, 1, 1, 1],   # pass 1
    [3, 3, 3, 3, 3, 3, 1, 1, 1, 1],   # pass 2 -- equal to pass 1, fixed point
]
PASSES = len(H)
N = 10

# Final kappa = converged values
KAPPA = H[-1]

# Color ramp by value (light = small, dark = large)
def color_for(v):
    # blue gradient
    if v <= 0:    return '#fafafa'
    if v == 1:    return '#dde8f5'
    if v <= 3:    return '#a9c4e6'
    if v <= 5:    return '#6d96cf'
    return '#3a6ebb'

def text_color_for(v):
    return '#1a1a1a' if v <= 3 else 'white'

cell_w = 0.95
cell_h = 0.55
gap_y  = 0.10
left   = 1.8     # space for vertex labels

fig, ax = plt.subplots(figsize=(3.4, 3.9))

TOP = 0.6
total_h = TOP + N * (cell_h + gap_y) + 0.25
ax.set_xlim(-0.2, left + PASSES * cell_w + 0.2)
ax.set_ylim(total_h, -0.4)
ax.axis('off')

# Header: pass labels
header_y = 0.10
ax.text(left / 2 - 0.05, header_y - 0.3, 'vertex',
        ha='center', va='center', fontsize=9.0, fontweight='bold')
ax.text(left / 2 - 0.05, header_y, r"$h^{(t)}[v]$",
        ha='center', va='center', fontsize=8.5, color='#666')
labels = [r'$h^{(0)}$' + '\n' + r'$=\!\sdeg$',
          r'$h^{(1)}$',
          r'$h^{(2)}\!=\!h^{(1)}$' + '\n(converged)']
labels_short = [r'$h^{(0)}\!=\!\sdeg$', r'$h^{(1)}$',
                r'$h^{(2)}\!=\!h^{(1)}$']
# matplotlib mathtext does not support \sdeg; use \mathrm{sdeg}
labels_short = [
    r"$h^{(0)}\!=\!\mathrm{sdeg}$",
    r"$h^{(1)}$",
    r"$h^{(2)}\!=\!h^{(1)}$",
]
for t in range(PASSES):
    x = left + t * cell_w + cell_w / 2
    ax.text(x, header_y - 0.30, f"pass {t}",
            ha='center', va='center', fontsize=8.0, color='#666')
    ax.text(x, header_y, labels_short[t],
            ha='center', va='center', fontsize=8.2,
            fontweight='bold' if t == PASSES - 1 else 'normal',
            color='#1f3a8a' if t == PASSES - 1 else '#222')

# Divider line
ax.plot([0.0, left + PASSES * cell_w],
        [TOP - 0.10, TOP - 0.10],
        linewidth=0.6, color='#999')

# Rows: one per vertex
for i in range(N):
    y0 = TOP + i * (cell_h + gap_y)
    y_mid = y0 + cell_h / 2

    # vertex label
    ax.text(left / 2 - 0.05, y_mid, f"$v_{{{i+1}}}$",
            ha='center', va='center', fontsize=9.5)

    for t in range(PASSES):
        x0 = left + t * cell_w
        v = H[t][i]
        is_final = (t == PASSES - 1)
        face = color_for(v)
        edge = '#1f3a8a' if is_final else '#bbb'
        lw = 1.0 if is_final else 0.4
        ax.add_patch(mpatches.Rectangle(
            (x0, y0), cell_w * 0.92, cell_h,
            facecolor=face, edgecolor=edge, linewidth=lw))
        ax.text(x0 + cell_w * 0.46, y_mid, f"{v}",
                ha='center', va='center',
                fontsize=9.0,
                fontweight='bold',
                color=text_color_for(v))

# Annotation: "= kappa_s" under the last column
note_y = TOP + N * (cell_h + gap_y) + 0.05
ax.text(left + (PASSES - 0.5) * cell_w, note_y,
        r"$=\kappa_3(v)$",
        ha='center', va='top', fontsize=8.5,
        color='#1f3a8a', fontweight='bold')

fig.savefig(OUT / "fig_spin_example.pdf", bbox_inches='tight')
fig.savefig(OUT / "fig_spin_example.png", dpi=180, bbox_inches='tight')
print(f"Saved fig_spin_example to {OUT}")
