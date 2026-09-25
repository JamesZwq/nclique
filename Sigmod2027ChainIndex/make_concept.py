#!/usr/bin/env python3
"""The conceptual figures of the paper, drawn on the running example (Figure 1) as hand-authored SVG and
converted to vector PDF with Chrome (the paper's font, Linux Libertine, is embedded as TrueType).

  figures/fig_running.pdf   Figure 1: the graph, and the community of x at sizes 2, 3 and 4 (shaded)
  figures/fig_pair.pdf      Figure 2 (Section 3): Baseline beside ChainIndex, one row per size: the same bars (nodes)
                            over the vertex arrays (44 positions) and over the chains with their runs, the label array
                            with the chain records, and the query of x at (3, 2) in orange
  figures/fig_build.pdf     Section 7: order replay at size 2 (two orders) and the refinement of the chains in the trie

Every object is recomputed here by brute force from the definitions (nuclei, own nodes, chains, ranks, depth-first
arrays, runs, entries, forward counts); the assertions pin the numbers the text prints.
Usage: python3 make_concept.py [fig_name ...]   (default: all four PDFs; add --png for previews in figures/preview/)
"""
import itertools, os, subprocess, sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / 'figures'
OUT.mkdir(exist_ok=True)
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
BLUE, ORANGE, GREY, SHADE, FAINT = "#1f5fa8", "#c2531a", "#8a8a8a", "#dfe9f5", "#cfcfcf"
ORANGE_SHADE = "#f6e3d7"

# ---------------------------------------------------------------- the running example, by brute force ----
names = ['b1', 'b2', 'b3', 'a1', 'a2', 'a3', 'a4', 'a5', 'u1', 'u2', 'u3', 'w1', 'w2', 'w3', 'x']
A = ['a1', 'a2', 'a3', 'a4', 'a5']; B = ['x', 'b1', 'b2', 'b3']; U = ['u1', 'u2', 'u3']; W = ['w1', 'w2', 'w3']
edges = set()
for S in (A, B):
    for p, q in itertools.combinations(S, 2): edges.add(frozenset((p, q)))
for u in U:
    for w in W: edges.add(frozenset((u, w)))
for y in U + W: edges.add(frozenset(('x', y)))
edges.add(frozenset(('a5', 'b3')))

def cliques(S, s):
    return [frozenset(c) for c in itertools.combinations(sorted(S), s)
            if all(frozenset((p, q)) in edges for p, q in itertools.combinations(c, 2))]

def nuclei(s, k):
    S = set(names)
    while True:
        cl = cliques(S, s); cnt = defaultdict(int)
        for c in cl:
            for v in c: cnt[v] += 1
        bad = {v for v in S if cnt[v] < k}
        if not bad: break
        S -= bad
    cl = cliques(S, s); parent = {v: v for v in S}
    def find(v):
        while parent[v] != v: parent[v] = parent[parent[v]]; v = parent[v]
        return v
    for c in cl:
        c = sorted(c)
        for v in c[1:]: parent[find(v)] = find(c[0])
    comps = defaultdict(set)
    for v in S:
        if cnt[v] >= k: comps[find(v)].add(v)
    return [frozenset(c) for c in comps.values()]

omega = {v: max(len(c) for k in range(1, 6) for c in cliques(set(names), k) if v in c) for v in names}
kappa = {v: {} for v in names}
for s in range(2, 6):
    k = 1
    while True:
        ns = nuclei(s, k)
        if not ns: break
        for N in ns:
            for v in N: kappa[v][s] = k
        k += 1
SIZES = [2, 3, 4, 5]
nodes = {}; own = {}
for s in SIZES:
    seen = {}
    for k in range(1, 8):
        for N in nuclei(s, k): seen.setdefault(N, []).append(k)
    nodes[s] = {N: max(ks) for N, ks in seen.items()}          # node -> top
    for v in names:
        if s in kappa[v]: own[(v, s)] = next(N for N in seen if v in N and max(seen[N]) == kappa[v][s])
NAME = {frozenset(A): 'A', frozenset(B): 'B', frozenset(B + U + W): 'N'}   # the sets the text names

def children(s):
    """parent map of the forest T_s (largest strictly containing node)."""
    par = {}
    for N in nodes[s]:
        sup = [M for M in nodes[s] if N < M]
        par[N] = min(sup, key=len) if sup else None
    return par

# chains = classes of equal trajectories (vertices in no edge: none here)
traj = {v: tuple(own[(v, s)] for s in SIZES if s in kappa[v]) for v in names}
chain_of_traj = {}
for v in names: chain_of_traj.setdefault(traj[v], []).append(v)
assert len(chain_of_traj) == 4

# preorder of each forest = the traversal of the layout pass: roots and children by the smallest chain rank they
# contain.  Ranks are the lexicographic order of the keys (preorder numbers of own nodes size by size); the two
# are consistent because a subtree's smallest rank is its first chain in preorder.  Computed by fixed point.
def preorder(s, rank):
    par = children(s); kids = defaultdict(list)
    for N, P in par.items():
        if P is not None: kids[P].append(N)
    def smallest_rank(N):
        return min(rank[chain_of_traj_key[v]] for v in N)
    order = []
    def visit(N):
        order.append(N)
        for C in sorted(kids[N], key=smallest_rank): visit(C)
    for R in sorted([N for N, P in par.items() if P is None], key=smallest_rank): visit(R)
    return {N: i for i, N in enumerate(order)}, kids

chain_of_traj_key = {v: traj[v] for v in names}
rank = {t: i for i, t in enumerate(chain_of_traj)}                # provisional
for _ in range(4):
    pre = {s: preorder(s, rank)[0] for s in SIZES}
    keys = {t: tuple(pre[SIZES[i]][t[i]] for i in range(len(t))) for t in chain_of_traj}
    rank = {t: i for i, t in enumerate(sorted(chain_of_traj, key=lambda t: keys[t]))}
pre = {s: preorder(s, rank)[0] for s in SIZES}; kids = {s: preorder(s, rank)[1] for s in SIZES}
chains = sorted(chain_of_traj, key=lambda t: rank[t])           # trajectories in rank order
members = {rank[t]: chain_of_traj[t] for t in chains}
assert [sorted(members[i]) for i in range(4)] == [['b1', 'b2', 'b3'], ['a1', 'a2', 'a3', 'a4', 'a5'],
                                                  ['u1', 'u2', 'u3', 'w1', 'w2', 'w3'], ['x']]
label_order = [v for i in range(4) for v in members[i]]          # aligned labels 0..14
label = {v: i for i, v in enumerate(label_order)}
start = {}; pos = 0
for i in range(4): start[i] = pos; pos += len(members[i])
start[4] = pos
assert [start[i] for i in range(5)] == [0, 3, 8, 14, 15]
chain_of = {v: i for i in range(4) for v in members[i]}
omega_c = {i: omega[members[i][0]] for i in range(4)}
kappa_c = {i: kappa[members[i][0]] for i in range(4)}
from math import comb
sigma_c = {}
for i in range(4):
    om = omega_c[i]; sig = om + 1
    for s in range(2, om + 1):
        if kappa_c[i][s] == comb(om - 1, s - 1): sig = s; break
    sigma_c[i] = sig
assert [omega_c[i] for i in range(4)] == [4, 5, 3, 4] and [sigma_c[i] for i in range(4)] == [2, 2, 4, 3]
residues = {i: [(s, kappa_c[i][s]) for s in range(2, sigma_c[i])] for i in range(4)}
assert residues == {0: [], 1: [], 2: [(2, 4), (3, 3)], 3: [(2, 4)]}

# per-size layers: depth-first array over chains, runs, entries
layers = {}
for s in SIZES:
    active = [i for i in range(4) if s in kappa_c[i]]
    ownc = defaultdict(list)                                       # node -> own chains
    for i in active: ownc[own[(members[i][0], s)]].append(i)
    arr = []; seg = {}
    def visit(N):
        seg[N] = [len(arr), None]
        items = [('c', i) for i in ownc[N]] + [('n', C) for C in kids[s][N]]
        def key(it): return it[1] if it[0] == 'c' else min(rank[traj[v]] for v in it[1])
        for it in sorted(items, key=key):
            if it[0] == 'c': arr.append(it[1])
            else: visit(it[1])
        seg[N][1] = len(arr)
    roots = sorted([N for N in nodes[s] if all(not (N < M) for M in nodes[s])], key=lambda N: min(rank[traj[v]] for v in N))
    for R in roots: visit(R)
    runs = []                                                       # (lo, hi) label ranges, maximal consecutive ranks
    for c in arr:
        if runs and runs[-1][2] == c - 1: runs[-1] = (runs[-1][0], start[c + 1], c)
        else: runs.append((start[c], start[c + 1], c))
    runs = [(lo, hi) for lo, hi, _ in runs]
    # entry of a node: run index holding the first label of its segment, and that label
    run_of_pos = []
    for c in arr:
        run_of_pos.append(next(r for r, (lo, hi) in enumerate(runs) if lo <= start[c] < hi))
    entry = {N: (run_of_pos[seg[N][0]], start[arr[seg[N][0]]]) for N in nodes[s]}
    sentinel = (len(runs), runs[-1][1])
    layers[s] = dict(arr=arr, runs=runs, entry=entry, sentinel=sentinel, ownc=ownc, roots=roots, seg=seg)
assert layers[3]['runs'] == [(0, 3), (8, 15), (3, 8)] and layers[3]['sentinel'] == (3, 8)
assert layers[2]['runs'] == [(0, 15)] and layers[4]['runs'] == [(0, 3), (14, 15), (3, 8)] and layers[5]['runs'] == [(3, 8)]
N3 = frozenset(B + U + W); assert layers[3]['entry'][N3] == (0, 0) and layers[3]['entry'][frozenset(A)] == (2, 3)

# order replay at size 2 on the two orders of Example 7.2
def replay(order, s=2):
    f = {}
    for i, v in enumerate(order):
        suffix = set(order[i:]); f[v] = sum(1 for c in cliques(suffix, s) if v in c)
    Umax = {}; m = 0
    for v in order: m = max(m, f[v]); Umax[v] = m
    inside = {}
    for v in order:
        Wv = {u for u in names if Umax[u] >= Umax[v]}; inside[v] = sum(1 for c in cliques(Wv, s) if v in c)
    ok = all(inside[v] >= Umax[v] for v in order)
    return f, Umax, inside, ok
ORDER1 = ['b1', 'b2', 'b3', 'a1', 'a2', 'a3', 'a4', 'a5', 'u1', 'u2', 'u3', 'w1', 'w2', 'w3', 'x']
ORDER2 = ['x', 'a1', 'a2', 'a3', 'a4', 'a5', 'b1', 'b2', 'b3', 'u1', 'u2', 'u3', 'w1', 'w2', 'w3']
f1, U1, in1, ok1 = replay(ORDER1); f2, U2, in2, ok2 = replay(ORDER2)
assert ok1 and all(U1[v] == kappa[v][2] for v in names) and not ok2
assert [f1[v] for v in ORDER1] == [3, 2, 2, 4, 3, 2, 1, 0, 4, 4, 4, 1, 1, 1, 0]

# ---------------------------------------------------------------- SVG helpers ----
def sub(name):
    """b1 -> b<sub>1</sub> as tspans (italic vertex names, upright subscripts)."""
    if name == 'x': return 'x'
    return f'{name[0]}<tspan font-size="70%" dy="1.6" font-style="normal">{name[1:]}</tspan><tspan dy="-1.6"> </tspan>'

def setname(vs):
    """a comma-free short name for a vertex set: A, B, N, or the list of members."""
    return NAME.get(vs, None)

class SVG:
    def __init__(self, w, h):
        self.w, self.h, self.parts = w, h, []
    def text(self, x, y, s, size=7.0, anchor='start', italic=False, fill='#000', weight=None, extra=''):
        st = f'font-size="{size}px" text-anchor="{anchor}" fill="{fill}"'
        if italic: st += ' font-style="italic"'
        if weight: st += f' font-weight="{weight}"'
        self.parts.append(f'<text x="{x:.2f}" y="{y:.2f}" {st} {extra}>{s}</text>')
    def line(self, x1, y1, x2, y2, color='#000', w=0.6, dash=None, cap='butt'):
        d = f' stroke-dasharray="{dash}"' if dash else ''
        self.parts.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{color}" stroke-width="{w}" stroke-linecap="{cap}"{d}/>')
    def rect(self, x, y, w, h, fill='none', stroke='#000', sw=0.6, rx=0, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ''
        self.parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')
    def circle(self, x, y, r, fill='#fff', stroke='#000', sw=0.6):
        self.parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
    def path(self, d, stroke='#000', w=0.6, fill='none', dash=None, cap='round'):
        dd = f' stroke-dasharray="{dash}"' if dash else ''
        self.parts.append(f'<path d="{d}" stroke="{stroke}" stroke-width="{w}" fill="{fill}" stroke-linecap="{cap}" stroke-linejoin="round"{dd}/>')
    def arrow(self, x1, y1, x2, y2, color='#000', w=0.6):
        import math
        self.line(x1, y1, x2, y2, color, w)
        ang = math.atan2(y2 - y1, x2 - x1); L = 2.6
        for da in (2.6, -2.6):
            self.line(x2, y2, x2 - L * math.cos(ang + da), y2 - L * math.sin(ang + da), color, w, cap='round')
    def bracket(self, x1, x2, y, tick=2.2, color='#000', w=0.6):
        """a horizontal bracket under [x1, x2] with end ticks pointing up."""
        self.line(x1, y, x2, y, color, w); self.line(x1, y, x1, y - tick, color, w); self.line(x2, y, x2, y - tick, color, w)
    def write(self, name):
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w:.1f} {self.h:.1f}" '
               f'font-family="Linux Libertine, Libertine, Times New Roman, serif">'
               f'<rect width="{self.w:.1f}" height="{self.h:.1f}" fill="#fff"/>' + '\n'.join(self.parts) + '</svg>')
        (OUT / f'{name}.svg').write_text(svg)
        html = OUT / f'.{name}.print.html'
        html.write_text(f"<!doctype html><meta charset='utf-8'><style>@page{{size:{self.w}pt {self.h}pt;margin:0}}"
                        f"html,body{{margin:0;padding:0}}svg{{display:block;width:{self.w}pt;height:{self.h}pt}}</style>" + svg)
        pdf = OUT / f'{name}.pdf'
        subprocess.run([CHROME, "--headless", "--disable-gpu", "--no-pdf-header-footer", f"--print-to-pdf={pdf}", f"file://{html}"], capture_output=True)
        html.unlink()
        if '--png' in sys.argv:
            (OUT / 'preview').mkdir(exist_ok=True)
            subprocess.run(['pdftoppm', '-r', '200', '-png', '-singlefile', str(pdf), str(OUT / 'preview' / name)])
        print('wrote', pdf.relative_to(HERE), f'{self.w}x{self.h}pt')

# ---------------------------------------------------------------- Figure 1: the graph ----
# layout (pt): A as a pentagon on the left, B as a diamond, x, and the K_{3,3} as two rows to the right of x
import math
POS = {}
for i in range(5):
    ang = math.radians(90 + 72 * i); POS[f'a{i+1}'] = (34 + 25 * math.cos(ang), 60 - 25 * math.sin(ang))
POS.update({'b3': (80, 60), 'b1': (106, 40), 'b2': (106, 80), 'x': (136, 60)})
for i in range(3):
    POS[f'u{i+1}'] = (170 + 31 * i, 24); POS[f'w{i+1}'] = (170 + 31 * i, 96)

def draw_graph(sv, ox, oy, scale, shaded=frozenset(), labels=True, r=6.2, lw=0.6):
    P = {v: (ox + x * scale, oy + y * scale) for v, (x, y) in POS.items()}
    for e in sorted(edges, key=lambda e: sorted(e)):
        p, q = sorted(e); dashed = e == frozenset(('a5', 'b3'))
        sv.line(*P[p], *P[q], color='#555' if not dashed else '#555', w=lw * scale if scale < 1 else lw, dash='1.6,1.4' if dashed else None)
    for v in names:
        x, y = P[v]; on = v in shaded
        sv.circle(x, y, r * scale, fill=('#8fb0d8' if on else '#fff'), stroke=(BLUE if on else '#000'), sw=(0.9 if on else 0.6) * max(scale, 0.7))
        if labels: sv.text(x, y + 2.3, sub(v), size=6.6, anchor='middle', italic=True)
    return P

def fig_running():
    PW, PH = 241.0, 158.0                      # 2026-09-22: the blank band under (a) removed; the panel letter stays in the caption
    sv = SVG(PW, PH)
    P = draw_graph(sv, 0, -6, 1.0)
    # the named sets, as light labels beside them
    sv.text(34, 92, 'A', size=7.5, anchor='middle', italic=True)
    sv.text(93, 94, 'B', size=7.5, anchor='middle', italic=True)
    sv.text(4, 8, '(a)', size=7.2)
    # (b)-(d): the community of x at its own level of sizes 2, 3, 4
    coms = [(2, frozenset(['x'] + U + W)), (3, frozenset(B + U + W)), (4, frozenset(B))]
    for j, (s, com) in enumerate(coms):
        assert com == own[('x', s)], (s, sorted(com))
        sc = 0.32; ox = 1 + j * 80; oy = 106
        draw_graph(sv, ox, oy, sc, shaded=com, labels=False, r=6.2, lw=0.5)
        k = kappa['x'][s]
        sv.text(ox + 118 * sc, oy + 118 * sc + 9, f'({"bcd"[j]}) s = {s}: level {k}, {len(com)} vertices', size=6.6, anchor='middle')
    sv.write('fig_running')

# ---------------------------------------------------------------- Section 3: one tree and one array per size ----
def strees_array(s):
    """the depth-first array of T_s: vertices in label order inside each node, nodes as the layout traversal."""
    L = layers[s]; out = []
    for c in L['arr']: out += members[c]
    return out

# ---------------------------------------------------------------- Section 7: replay and refinement ----
def fig_build(trie_only=True):
    # 2026-09-22: the paper shows the trie alone at column width (the replay tables are Example 7.2); trie_only=False draws both panels
    PW, PH = (241.0, 136.0) if trie_only else (506.0, 150.0)
    sv = SVG(PW, PH)
    # (a) order replay at size 2, two orders
    if not trie_only: sv.text(2, 9, '(a) order replay at s = 2', size=7.2)
    CW, CH = 14.2, 10.0; X0 = 34.0
    def block(y, order, f, Umax, inside, ok, title):
        sv.text(X0, y - 3, title, size=6.2, fill='#555')
        rows = [('π', [sub(v) for v in order], True),
                ('f<tspan font-size="70%" dy="1.6">2</tspan><tspan dy="-1.6"> </tspan>', [str(f[v]) for v in order], False),
                ('U<tspan font-size="70%" dy="1.6">2</tspan><tspan dy="-1.6"> </tspan>', [str(Umax[v]) for v in order], False),
                ('inside', [str(inside[v]) for v in order], False)]
        for r, (lab, cells, it) in enumerate(rows):
            yy = y + r * CH
            sv.text(X0 - 3, yy + 8.2, lab, size=6.4, anchor='end', italic=(r < 3))
            for i, cell in enumerate(cells):
                x = X0 + i * CW
                bad = (r == 3 and inside[order[i]] < Umax[order[i]])
                sv.rect(x, yy, CW, CH, fill=(ORANGE_SHADE if bad else 'none'), stroke='#000', sw=0.4)
                sv.text(x + CW / 2, yy + 8.2, cell, size=6.3, anchor='middle', italic=it, fill=(ORANGE if bad else '#000'))
        sv.text(X0 + 15 * CW, y + 4 * CH + 8, ('the certificate holds: U<tspan font-size="70%" dy="1.6">2</tspan><tspan dy="-1.6"> = κ</tspan><tspan font-size="70%" dy="1.6">2</tspan>' if ok
                                                else 'the certificate fails: inside &lt; U<tspan font-size="70%" dy="1.6">2</tspan><tspan dy="-1.6"> at </tspan><tspan font-style="italic">a</tspan><tspan font-size="70%" dy="1.6">1</tspan>'),
                size=6.2, anchor='end', fill=(BLUE if ok else ORANGE))
    if not trie_only:
        block(22, ORDER1, f1, U1, in1, ok1, 'the degeneracy order')
        block(84, ORDER2, f2, U2, in2, ok2, 'x first')
        sv.text(2, PH - 4, 'inside: the edges of v inside {u : U<tspan font-size="70%" dy="1.6">2</tspan><tspan dy="-1.6">(u) ≥ U</tspan><tspan font-size="70%" dy="1.6">2</tspan><tspan dy="-1.6">(v)}</tspan>', size=5.8, fill='#555')
    # (b) the refinement of the chains in the prefix trie
    TX = 30.0 if trie_only else 322.0
    if not trie_only: sv.text(TX - 4, 9, '(b) the chains, refined size by size in the trie', size=7.2)
    trie = {}          # path of own-node ids -> vertices that pass through
    stop = defaultdict(list)
    for v in names:
        path = ()
        for s in SIZES:
            if s not in kappa[v]: break
            path = path + (pre[s][own[(v, s)]],); trie.setdefault(path, set()).add(v)
        stop[path].append(v)
    paths = sorted(trie)
    xpos = {}; col = [0]; COLW = 48.0
    def place(p):
        kidsp = [q for q in paths if len(q) == len(p) + 1 and q[:len(p)] == p]
        if stop.get(p):
            xpos[p] = TX + 20 + col[0] * COLW; col[0] += 1
            for q in kidsp: place(q)
        else:
            xs = []
            for q in kidsp: place(q); xs.append(xpos[q])
            xpos[p] = sum(xs) / len(xs)
    place(())
    Y = {0: 8, 1: 25, 2: 43, 3: 61, 4: 79, 5: 97} if trie_only else {0: 18, 1: 36, 2: 56, 3: 76, 4: 96, 5: 116}
    if trie_only: sv.h = 124.0
    NW, NH = 34, 11
    for p in [()] + paths:
        d = len(p); x = xpos[p]; y = Y[d]
        if d == 0:
            sv.circle(x, y, 2.0, fill='#000')
        else:
            s = SIZES[d - 1]; N = next(N for N in nodes[s] if pre[s][N] == p[-1]); nm = setname(N)
            lab = f'{p[-1]}' + (f' = {nm}' if nm else '')
            sv.rect(x - NW / 2, y - NH / 2, NW, NH, rx=2, sw=0.6)
            sv.text(x, y + 2.4, lab, size=6.0, anchor='middle')
            parent = p[:-1]; pxx, pyy = xpos[parent], Y[d - 1]
            sv.line(pxx, pyy + (2.0 if d == 1 else NH / 2), x, y - NH / 2, w=0.5)
        if stop.get(p):
            # the vertices that stop here hang below the node as a leaf in the node's own column
            vs = sorted(stop[p], key=lambda v: label[v]); ly = Y[d + 1]
            txt = ', '.join(sub(v) for v in vs) if len(vs) <= 3 else f'{sub(vs[0])}, …, {sub(vs[-1])}'
            sv.line(x, y + (NH / 2 if d else 2.0), x, ly - NH / 2, w=0.5, color=BLUE)
            sv.rect(x - 19, ly - NH / 2, 38, NH, rx=2, sw=0.6, stroke=BLUE, dash='1.6,1.2')
            sv.text(x, ly + 2.4, txt, size=5.8, anchor='middle', italic=True, fill=BLUE)
            c = chain_of[vs[0]]
            sv.text(x, ly + 11.5, f'chain {c}', size=5.6, anchor='middle', fill='#555')
            sv.text(x, ly + 18.5, f'labels [{start[c]}, {start[c+1]})', size=5.6, anchor='middle', fill='#555')
    for d, s in enumerate(SIZES):
        sv.text(TX - 4, Y[d + 1] + 2.4, f's = {s}', size=6.0, anchor='end', fill='#777')
    sv.write('fig_build')

# ---------------------------------------------------------------- Figure 2: the baseline beside ChainIndex ----
# 2026-09-24 cold read: tints separated in lightness for grayscale print (L* 88, 95, 80, 72; adjacent pairs differ by >= 7),
# chroma kept low (10-21) so that no tint dominates
CHAIN_TINT = {0: '#d0ddee', 1: '#f7f0e1', 2: '#b6cfad', 3: '#b7a9cf'}

def kappa_txt(s, k):
    return f'<tspan font-style="italic">κ</tspan><tspan font-size="75%" dy="1.7">{s}</tspan><tspan dy="-1.7"> = {k}</tspan>'

def fig_pair():
    """2026-09-24 (user): one figure, (a) Baseline beside (b) ChainIndex, one row per size and the rows aligned.
    Both panels draw the tree of a size the same way: bars stacked by depth, each over the members of its node and
    labelled with its name and top.  (a) puts them over the depth-first array of the vertices of the size (44 positions
    over the four sizes); (b) over the chains in depth-first order, whose runs are the brackets labelled with their
    label ranges, under the one label array with the four chains and their records.  Cells are tinted by chain.
    Orange: the community query of x at (3, 2), one slice in (a), two runs in (b), two label ranges on top.
    Cold read 2026-09-24: larger type (names 7 pt, bars 6.6 pt, records and runs 6.4 pt, gray numbers 5.8 pt), the
    traced bar outlined only, the rows of (b) labelled, the records tagged, the record of chain 3 in black."""
    PW = 506.0
    CW, CH = 13.8, 12.0                # vertex and label cells
    CC = 30.0                          # chain cells of (b)
    BH, BG = 10.5, 1.6                 # bar height, gap between bar rows
    XA, XB = 26.0, 262.0               # left edges of the cells of (a) and of (b)
    GRAY = '#555'
    sv = SVG(PW, 1)
    hot_s, hot_N = 3, own[('x', 3)]
    assert hot_N == frozenset(B + U + W)
    sv.text(2, 9.0, '(a) Baseline', size=7.8)
    sv.text(XB - 20, 9.0, '(b) ChainIndex', size=7.8)
    # (b) the label array, the chains, their records: one column at the right edge, beside the rows
    # 2026-09-25 (user: the figure takes too much space): the column replaces the header band over (b), so the rows
    # of both panels start at the top.
    LX, LH, LY = 424.0, 11.0, 19.0     # left edge of the label cells, cell height, top of the column
    sv.text(LX + CW / 2, LY - 4.0, 'label', size=5.8, anchor='middle', fill=GRAY)
    for i, v in enumerate(label_order):
        yy = LY + i * LH
        sv.rect(LX, yy, CW, LH, fill=CHAIN_TINT[chain_of[v]], stroke='#000', sw=0.5)
        sv.text(LX - 3.0, yy + 7.9, str(i), size=5.8, anchor='end', fill=GRAY)
        sv.text(LX + CW / 2, yy + 8.2, sub(v), size=7.0, anchor='middle', italic=True)
    for lo, hi in ((0, 3), (8, 15)):                                   # the answer of the traced query
        sv.line(LX - 11.5, LY + lo * LH + 0.6, LX - 11.5, LY + hi * LH - 0.6, color=ORANGE, w=1.2)
    sv.rect(LX, LY + 14 * LH, CW, LH, stroke=ORANGE, sw=1.2)           # the query vertex x
    RX = LX + CW + 4.0
    sv.text(RX + 3.5, LY - 4.0, 'record', size=5.8, fill=GRAY)
    for c in range(4):
        lo, hi = LY + start[c] * LH, LY + start[c + 1] * LH
        sv.line(RX, lo + 1, RX, hi - 1, w=0.6); sv.line(RX, lo + 1, RX - 2.4, lo + 1, w=0.6); sv.line(RX, hi - 1, RX - 2.4, hi - 1, w=0.6)
        lines = [f'{c}: <tspan font-style="italic">ω</tspan> = {omega_c[c]}, <tspan font-style="italic">σ</tspan> = {sigma_c[c]}']
        if residues[c]: lines.append('\u2002\u2002' + ', '.join(kappa_txt(s, k) for s, k in residues[c]))
        ty = lo + 8.0 if hi - lo < 2 * LH else (lo + hi) / 2 - 4.0 * (len(lines) - 1) + 2.2
        for t, line in enumerate(lines):
            sv.text(RX + 3.5, ty + t * 8.0, line, size=6.4)
    y = LY
    total_a = 0; total_runs = 0
    for j, s in enumerate(SIZES):
        L = layers[s]; par = children(s)
        def dep(N): return 0 if par[N] is None else 1 + dep(par[N])
        maxd = max(dep(N) for N in nodes[s])
        # numbering the bars of (b) by their left ends (outer first on a tie) is the preorder the text uses
        lr = sorted(nodes[s], key=lambda N: (L['seg'][N][0], dep(N)))
        assert [pre[s][N] for N in lr] == list(range(len(lr))), s
        assert len({L['seg'][N][0] for N in nodes[s]}) == len(nodes[s]), s   # no two bars share a left end
        if j: sv.line(2, y - 3.2, LX - 16.0, y - 3.2, color='#dddddd', w=0.4)
        cy = y + (maxd + 1) * (BH + BG)
        sv.text(XA - 4, cy + 8.6, f's = {s}', size=7.0, anchor='end')
        sv.text(XB - 4, cy + 8.6, f's = {s}', size=6.4, anchor='end', fill=GRAY)
        # (a): the depth-first array of the vertices, the bars over it
        arr_v = strees_array(s); total_a += len(arr_v)
        pa = {v: XA + i * CW for i, v in enumerate(arr_v)}
        for v in arr_v:
            hot = (s == hot_s and v == 'x')
            sv.rect(pa[v], cy, CW, CH, fill=CHAIN_TINT[chain_of[v]], stroke=(ORANGE if hot else '#000'), sw=(1.2 if hot else 0.5))
            sv.text(pa[v] + CW / 2, cy + 8.8, sub(v), size=7.0, anchor='middle', italic=True)
        # (b): the chains in depth-first order
        for i, c in enumerate(L['arr']):
            hot = (s == hot_s and c == 3)
            sv.rect(XB + i * CC, cy, CC, CH, fill=CHAIN_TINT[c], stroke=(ORANGE if hot else '#000'), sw=(1.2 if hot else 0.5))
            sv.text(XB + i * CC + CC / 2, cy + 8.8, str(c), size=7.0, anchor='middle')
        # the same bars in both panels; the traced bar outlined in orange, no fill
        for N in nodes[s]:
            nm = setname(N); lab = (f'{nm}, top {nodes[s][N]}' if nm else f'top {nodes[s][N]}')
            hot = (s == hot_s and N == hot_N); by = y + dep(N) * (BH + BG)
            a0, a1 = L['seg'][N]
            vlo = members[L['arr'][a0]][0]; vhi = members[L['arr'][a1 - 1]][-1]
            for lo, hi in ((pa[vlo], pa[vhi] + CW), (XB + a0 * CC, XB + a1 * CC)):
                sv.rect(lo + 0.9, by, hi - lo - 1.8, BH, rx=1.5, fill='#fff',
                        stroke=(ORANGE if hot else '#000'), sw=(1.2 if hot else 0.55))
                sv.text((lo + hi) / 2, by + 7.6, lab, size=6.6, anchor='middle')
        # (b): the runs under the chains
        ry = cy + CH + 4.5
        for r, (lo_l, hi_l) in enumerate(L['runs']):
            idx = [i for i, c in enumerate(L['arr']) if lo_l <= start[c] < hi_l]
            assert idx == list(range(idx[0], idx[-1] + 1))
            hot = (s == hot_s and (lo_l, hi_l) in ((0, 3), (8, 15)))
            x0, x1 = XB + idx[0] * CC, XB + (idx[-1] + 1) * CC
            sv.bracket(x0 + 1.5, x1 - 1.5, ry, tick=2.4, w=(0.9 if hot else 0.6), color=(ORANGE if hot else '#000'))
            sv.text((x0 + x1) / 2, ry + 7.8, f'[{lo_l}, {hi_l})', size=6.4, anchor='middle', fill=(ORANGE if hot else '#000'))
            total_runs += 1
        y = ry + 7.8 + 6.5
    assert total_a == 44 and total_runs == 8
    sv.h = max(y - 1.0, LY + 15 * LH + 18.0)
    sv.write('fig_pair')

FIGS = {'fig_running': fig_running, 'fig_pair': fig_pair, 'fig_build': fig_build}

if __name__ == '__main__':
    only = [a for a in sys.argv[1:] if a in FIGS]
    for name in (only or list(FIGS)): FIGS[name]()
