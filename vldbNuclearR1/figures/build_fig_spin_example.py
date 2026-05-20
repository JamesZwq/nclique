"""
Figure for Example 5.2 (SPIN on the running example): visualise
vertexSupport(v_4, h) as four path-lifelines, parallel in layout to
Figure 2 for SPIN*.

Each path through v_4 (L_3, L_4, L_5, L_6) gets one horizontal row.
Columns are thresholds k = 1, ..., h[v_4] = 6.  Each cell shows
W_{v_4}(P, h, k) = number of v_4-incident encoded s-cliques whose other
vertices all have h >= k.  Once W drops to 0 the rest of the row is a
hatched grey strip (parallel to "leaf already dead" in Figure 2).

A bottom aggregate row gives the cumulative sum sum_P W and highlights
the largest k satisfying sum >= k, which is Phi_h(v_4) = kappa_3(v_4).

Numbers below are computed by hand from Table 4.1 and the initial
supports h[v_i] = sdeg(v_i) of the running example.
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['text.usetex']  = False
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

OUT = Path(__file__).parent

# Initial threshold witness counts W_{v_4}(P, h, k) for k = 1..6.
# Derivation (h initialised to sdeg = initial s-clique support):
#   h[v_1]=3, h[v_2]=4, h[v_3]=6, h[v_4]=6, h[v_5]=3, h[v_6]=4,
#   h[v_7]=h[v_8]=h[v_9]=h[v_10]=1.
#
#  L_3 = ({v_1}, {v_2,v_3,v_4}), v_4 is pivot, eta=2.
#    hold gate: h[v_1] >= k iff k <= 3.
#    Q_k = pivots from {v_2,v_3} with h >= k: 2 for k<=4, 1 for k=5,... 0 for k>6.
#    W = C(Q_k, eta-1) = C(Q_k, 1) when gate holds, else 0.
#       k=1..3: gate holds, Q_k=2, W=2
#       k=4:    gate fails, W=0
#       k=5,6:  W=0
LEAVES = [
    ("L_3", r"$v_4$ pivot, $\eta\!=\!2$", [2, 2, 2, 0, 0, 0]),
    # L_4 = ({v_2}, {v_3,v_4}), v_4 pivot.
    #   gate h[v_2]>=k iff k<=4. Q_k = {v_3} with h=6 -> 1 for k<=6.
    #   W = C(Q_k, 1) = 1 when gate holds, else 0.
    ("L_4", r"$v_4$ pivot, $\eta\!=\!2$", [1, 1, 1, 1, 0, 0]),
    # L_5 = ({v_3}, {v_4,v_5,v_6}), v_4 pivot.
    #   gate h[v_3]>=k iff k<=6. Q_k = pivots {v_5,v_6} with h>=k.
    #     h[v_5]=3, h[v_6]=4. Q_k = 2 for k<=3, 1 for k=4, 0 for k>=5.
    #   W = C(Q_k, 1).  k=1..3:2, k=4:1, k=5,6:0.
    ("L_5", r"$v_4$ pivot, $\eta\!=\!2$", [2, 2, 2, 1, 0, 0]),
    # L_6 = ({v_4}, {v_5,v_6}), v_4 hold.
    #   No other hold so gate trivially holds. Q_k = pivots {v_5,v_6} with h>=k.
    #     k<=3: 2; k=4: 1; k>=5: 0.
    #   W = C(Q_k, eta) = C(Q_k, 2). k<=3:1, k=4:0, k>=5:0.
    ("L_6", r"$v_4$ hold,  $\eta\!=\!2$",  [1, 1, 1, 0, 0, 0]),
]

K_MAX = 6
totals = [sum(L[2][k-1] for L in LEAVES) for k in range(1, K_MAX + 1)]
# largest k with total >= k
phi = max((k for k in range(1, K_MAX + 1) if totals[k-1] >= k), default=0)

fig, ax = plt.subplots(figsize=(6.0, 2.5))

GREEN  = '#3a8a3a'
LIGHT  = '#eaf4ea'
LIVE   = '#fafafa'
GREY   = '#d0d0d0'
HIGHL  = '#1f3a8a'
RED    = '#c0392b'

cell_w = 1.0
cell_h = 0.60
gap_y  = 0.20
rows   = len(LEAVES) + 1   # +1 for the aggregate row
TOP    = 0.85
total_h = TOP + rows * (cell_h + gap_y) + 0.35

ax.set_xlim(-2.6, K_MAX + 0.5)
ax.set_ylim(total_h, -0.2)
ax.axis('off')

# --- Header: threshold k columns ---
ax.text(-0.05, 0.20, r"threshold $k$:", ha='right', va='center',
        fontsize=8.5, fontweight='bold', color='#444')
for k in range(1, K_MAX + 1):
    x = k - 0.5
    ax.text(x, 0.20, f"{k}", ha='center', va='center',
            fontsize=9.0, fontweight='bold',
            color=HIGHL if k == phi else '#444')

# divider under header
ax.plot([-2.4, K_MAX + 0.3], [0.55, 0.55], linewidth=0.6, color='#999')

# --- Path rows ---
for i, (Lname, role, vals) in enumerate(LEAVES):
    y0 = TOP + i * (cell_h + gap_y)
    y_mid = y0 + cell_h / 2

    label = rf"${Lname}$ ({role})"
    ax.text(-0.05, y_mid, label, ha='right', va='center', fontsize=8.0)

    # find death k (first 0)
    death_at = next((k for k in range(1, K_MAX + 1) if vals[k-1] == 0), K_MAX + 1)

    # alive background up to death_at - 1, grey afterward
    if death_at > 1:
        ax.add_patch(mpatches.Rectangle(
            (0, y0), death_at - 1, cell_h,
            facecolor=LIGHT, edgecolor='#bbb', linewidth=0.3))
    if death_at <= K_MAX:
        ax.add_patch(mpatches.Rectangle(
            (death_at - 1, y0), K_MAX - death_at + 1, cell_h,
            facecolor=GREY, edgecolor='#bbb', linewidth=0.3,
            hatch='///', alpha=0.55))

    # Witness values
    for k in range(1, K_MAX + 1):
        x = k - 0.5
        w = vals[k-1]
        if w > 0:
            ax.text(x, y_mid, f"{w}", ha='center', va='center',
                    fontsize=8.2, fontweight='bold', color=GREEN)
        # cells with w == 0 sit in the hatched dead strip; no text

# --- Aggregate row: sum_P W ---
agg_y0 = TOP + len(LEAVES) * (cell_h + gap_y) + 0.10
agg_mid = agg_y0 + cell_h / 2

ax.text(-0.05, agg_mid, r"$\sum_{P\ni v_4} W_{v_4}(P,h,k)$",
        ha='right', va='center', fontsize=8.5, fontweight='bold')

# separator line above aggregate
ax.plot([0, K_MAX], [agg_y0 - 0.06, agg_y0 - 0.06],
        linewidth=0.6, color='#999')

for k in range(1, K_MAX + 1):
    x = k - 0.5
    s = totals[k-1]
    cond = s >= k
    is_phi = (k == phi)
    face = LIGHT if cond else '#fafafa'
    if is_phi:
        face = '#ffe9b3'
    ax.add_patch(mpatches.Rectangle(
        (k - 1, agg_y0), 1.0, cell_h,
        facecolor=face,
        edgecolor=HIGHL if is_phi else '#bbb',
        linewidth=1.2 if is_phi else 0.3))
    txt_color = HIGHL if is_phi else ('#1a1a1a' if cond else '#999')
    ax.text(x, agg_mid, f"{s}", ha='center', va='center',
            fontsize=9.0,
            fontweight='bold' if is_phi else 'normal',
            color=txt_color)

# Annotate Phi
ax.text(K_MAX + 0.6, agg_mid,
        rf"$\Phi_h(v_4)\!=\!{phi}\!=\!\kappa_3(v_4)$",
        ha='left', va='center', fontsize=8.5, color=HIGHL,
        fontweight='bold')

# Legend
ax.text(-2.5, total_h - 0.10,
        r"Green digit: $W_{v_4}(P,h,k)\!>\!0$ (path qualifies at threshold $k$). "
        r"Hatched grey: $W\!=\!0$ (no qualifying clique through $v_4$).",
        ha='left', va='center', fontsize=7.0, color='#555')

fig.savefig(OUT / "fig_spin_example.pdf", bbox_inches='tight')
fig.savefig(OUT / "fig_spin_example.png", dpi=180, bbox_inches='tight')
print(f"Saved fig_spin_example to {OUT}")
