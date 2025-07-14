import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

# — Academic style settings —
rcParams['font.family']      = 'Times New Roman'
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm']      = 'Times New Roman'
rcParams['mathtext.it']      = 'Times New Roman Italic'
rcParams['mathtext.bf']      = 'Times New Roman Bold'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

# — Path to your data file —
DATA_PATH = '/Users/zhangwenqian/UNSW/pivoter/dblpClique.txt'
# DATA_PATH = '/Users/zhangwenqian/UNSW/pivoter/StanfordClique.txt'
# — Read and parse only the lines that start with a clique‐size integer —
ks = []
counts = []
with open(DATA_PATH, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or not line[0].isdigit():
            continue
        parts = line.split(',')
        # parts[0] = k, parts[1] = count
        ks.append(int(parts[0].strip()))
        counts.append(float(parts[1].strip()))

ks = np.array(ks)
counts = np.array(counts)

# — Plotting —
fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(ks, counts,
       width=0.5,
       facecolor='black',
       edgecolor='black')

ax.set_yscale('log')
ax.set_xlabel('Clique Size $(k)$', fontsize=18, labelpad=8)
ax.set_ylabel('Number of Cliques', fontsize=18, labelpad=8)

ax.tick_params(axis='x', rotation=-40, labelsize=14)
ax.tick_params(axis='y', labelsize=14)

plt.tight_layout()
plt.savefig('/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/clique_distribution_dblp.eps', dpi=300)
# plt.savefig('/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/clique_distribution_stanford.eps', dpi=300)
# plt.show()
