#!/usr/bin/env python3
# Scalability figures for the NSI paper (§149 data). House style: monochrome, solid black =
# primary, dashed gray = secondary; no top/right spines; wide-and-flat; ~9pt fonts.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 200})

# Panel A: whole-spectrum time vs order r (log y). webit + dblp.
rA = [3,4,5,6,7]
webit_A = [0.81,1.52,2.90,5.29,14.82]
dblp_A  = [2.92,5.94,18.14,47.08,97.01]
fig, ax = plt.subplots(figsize=(3.3,2.0))
ax.plot(rA, webit_A, "-o", color="black", markerfacecolor="black", label="web-it")
ax.plot(rA, dblp_A, "--s", color="0.45", markerfacecolor="white", label="com-dblp")
ax.set_yscale("log"); ax.set_xlabel("order $r$"); ax.set_ylabel("time (s)")
ax.set_xticks(rA); ax.legend(frameon=False, loc="upper left", fontsize=8)
fig.tight_layout(pad=0.3); fig.savefig("scal_r.pdf"); plt.close(fig)

# Panel B: whole-spectrum time vs range width smax (linear y). Near-flat = the headline.
smaxB = [6,7,8,9,10,11,12]
webit_B = [1.41,1.47,1.56,1.58,1.62,1.70,1.73]
dblp_B  = [5.27,5.34,5.46,5.47,5.58,6.00,6.09]
fig, ax = plt.subplots(figsize=(3.3,2.0))
ax.plot(smaxB, webit_B, "-o", color="black", markerfacecolor="black", label="web-it")
ax.plot(smaxB, dblp_B, "--s", color="0.45", markerfacecolor="white", label="com-dblp")
ax.set_xlabel("range width $s_{\\max}$"); ax.set_ylabel("time (s)")
ax.set_ylim(0, 7); ax.set_xticks(smaxB); ax.legend(frameon=False, loc="center left", fontsize=8)
fig.tight_layout(pad=0.3); fig.savefig("scal_s.pdf"); plt.close(fig)
print("wrote scal_r.pdf, scal_s.pdf")
