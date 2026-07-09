#!/usr/bin/env python3
# Figure 1: the running example. 9 vertices, three maximal cliques M1,M2,M3, with the (3,4)-nucleus
# answer overlaid as two nested tinted contours (teacher style: the only colourful figure).
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams.update({"font.size": 10})

pos = {  # v1,v2 M1-only | v3,v4,v5 M1cap M2 | v6,v7 M2cap M3 | v8,v9 M3-only
 "v_1":(0.0,1.7),"v_2":(0.0,0.3),
 "v_3":(1.6,2.15),"v_4":(1.6,1.0),"v_5":(1.6,-0.15),
 "v_6":(3.2,1.7),"v_7":(3.2,0.3),
 "v_8":(4.8,1.7),"v_9":(4.8,0.3)}
M1=["v_1","v_2","v_3","v_4","v_5"]; M2=["v_3","v_4","v_5","v_6","v_7"]; M3=["v_6","v_7","v_8","v_9"]

fig, ax = plt.subplots(figsize=(3.4,2.3))

# nested nucleus contours (drawn first, behind everything)
def blob(verts, pad, color, z):
    xs=[pos[v][0] for v in verts]; ys=[pos[v][1] for v in verts]
    x0,x1,y0,y1=min(xs)-pad,max(xs)+pad,min(ys)-pad,max(ys)+pad
    ax.add_patch(FancyBboxPatch((x0,y0),x1-x0,y1-y0,
        boxstyle="round,pad=0.02,rounding_size=0.35",
        fc=color, ec="none", zorder=z))
# 1-(3,4)-nucleus = all 9 (outer, lightest); 2-(3,4)-nucleus = M1 u M2 = v1..v7 (inner)
blob(list(pos), 0.62, "#dfe8f3", 0)
blob(M1+M2,      0.42, "#b9cde8", 1)
ax.text(4.85,-0.75,"1-(3,4)-nucleus",fontsize=7.5,color="#4a6da3",ha="right")
ax.text(0.05,2.72,"2-(3,4)-nucleus",fontsize=7.5,color="#2f4f7f")

# edges: all pairs within each clique
seen=set()
for clq in (M1,M2,M3):
    for i in range(len(clq)):
        for j in range(i+1,len(clq)):
            e=tuple(sorted((clq[i],clq[j])))
            if e in seen: continue
            seen.add(e)
            (x0,y0),(x1,y1)=pos[e[0]],pos[e[1]]
            ax.plot([x0,x1],[y0,y1],"-",color="0.55",lw=0.8,zorder=3)

# nodes
for v,(x,y) in pos.items():
    ax.plot(x,y,"o",ms=15,mfc="white",mec="black",mew=1.1,zorder=4)
    ax.text(x,y,"$"+v+"$",fontsize=7.5,ha="center",va="center",zorder=5)

ax.set_xlim(-0.9,5.7); ax.set_ylim(-0.95,2.95); ax.axis("off")
fig.tight_layout(pad=0.1); fig.savefig("fig1_running.pdf"); print("wrote fig1_running.pdf")
