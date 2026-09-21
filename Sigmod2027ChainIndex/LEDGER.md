# Ledger (frozen 2026-09-21, paper-architect Phase 0)

Every unit of prose reads this file. Nothing here is invented; every number has its source.

## Story line (one sentence)

A community query over clique-core hierarchies names the clique size as well as the level, and
serving every size has meant one hierarchy per size; we show that the vertices no query can tell
apart are exactly those with the same own node at every size (chains), that this partition is the
coarsest an exact index can use, that three layout decisions turn the partition into ranges of
memory answered in constant time, that one clique tree builds every size at once, and that the
index is 1.28x to 48.7x smaller than one tree per size while listing at the speed of a memory copy.

## Slot ledger (teacher pattern -> this paper)

| slot | filling |
|---|---|
| model / object | $(1,s)$-nucleus decomposition; the community of $v$ at (size $s$, level $k$) |
| task | index every size at once: value $(v,s)$ and community $(v,s,k)$ |
| state of the art | one hierarchy per size, `\strees` (tree + DFS array per size); per-size construction `\cnd` |
| the one central limitation (bold) | **one hierarchy per size: a vertex is stored once for every size it lives in** |
| research question | How to index every clique-core hierarchy of a graph at once, storing and retrieving each vertex once? |
| the one concept | *chains* (hierarchy equivalence classes = same own node at every size) |
| challenge nouns (3) | repetition across sizes; contiguity of communities; construction of every size |
| mirrored limitations (numbers) | (1) `\strees` stores a vertex once per size: 190 MB on web-uk-2005 for 129,632 vertices (181 arrays per vertex on average); twins buy only 1.03 to 1.28x (design study). (2) deduplicated layouts (SGL-style skyline over chains) list 2 to 3 times slower than a tree layout (design study). (3) one decomposition per size: 15.7 h of build and peel on web-uk-2005 (499 sizes) with `\cnd` |
| our idea features (3, echo the challenges) | Chains (coarsest lossless partition); Aligned labels, chain order and run arrays (constant-time ranges); One clique tree for every size (order replay) |
| contributions (4) | Chains (Theorem chains, Proposition coarsest); ChainIndex (Theorems onerange, retrieval; Theorem tail for values); Construction (Theorem replay); Comprehensive experimental evaluation |
| headline numbers (verbatim everywhere) | 1.28x to 48.7x smaller than STrees, median 7.0x; locate 6 to 39 ns; listing within 0.2x to 2.8x of a memory copy (median 1.02); build 2.8x to 814x faster than one decomposition per size (median 16x); 29 graphs |
| algorithm names | `\query` (Algorithm Query), `\build` (Algorithm Build) |
| Exp list | Exp-1 size vs STrees; Exp-2 size by part; Exp-3 design study; Exp-4 query latency (+levels figure); Exp-5 latency by clique size and answer size; Exp-6 construction; Exp-7 scalability; Exp-8 vs CND; Exp-9 case study I (ground truth); Exp-10 case study II (Amazon zoom) |
| case studies | I: SNAP ground truth, best (s,k) per query (dblp 3.6x k-core; best s spread 2..8); II: Amazon with titles/categories, median community 157,331 -> 38 -> 11 -> 8 -> 7, purity 0.10 -> 0.44 -> 0.70 -> 0.79 -> 0.83 |

## Notation (one symbol per quantity; first defined where stated)

| symbol | meaning | defined |
|---|---|---|
| $G=(V,E)$, $n$ | graph, number of vertices | Prelim |
| $s$ | clique size (query parameter) | Prelim |
| $k$ | level | Prelim |
| $\omg{v}$ | size of a largest clique containing $v$ | Prelim |
| $\core{s}{v}$ | core value of $v$ at size $s$ | Prelim (Def nucleus) |
| $\Tree{s}$ | forest of canonical nodes of size $s$ | Sec hierarchy (Def node) |
| $\own{s}{v}$ | own node of $v$ at size $s$ | Def node |
| top($X$) | largest level of node $X$ = smallest core value of its vertices | Def node |
| $\traj{v}$, $\chain{v}$ | trajectory (tuple of own nodes) and chain (set of vertices with equal trajectories; isolated vertices form one last chain) | Def chain |
| $f_s(v)=\binom{\omg{v}-1}{s-1}$ | floor | Theorem tail |
| $\sig{v}$ | certification point | after Theorem tail |
| $C, N, P, R$ | chains; canonical nodes over all sizes; (chain,size) pairs; runs over all sizes | Sec size |
| $(lo,hi)$, $r_X$, $v_X$ | run endpoints; entry of node $X$: run index and label; sentinel $(R, hi(R-1))$; answer = first range, whole runs, last range | Sec runs |
| $M$ | canonical nodes over all sizes (size accounting) | Sec size |
| $d_s(W,u)$, $f_s$, $U_s$ | s-cliques of $u$ in $G[W]$; forward count; prefix max | Sec construction |

## Running example (Figure 1) -- computed by tools/example_check.py (2026-09-21)

Vertices: $A=\{a_1..a_5\}$ ($K_5$), $B=\{x,b_1,b_2,b_3\}$ ($K_4$), $U=\{u_1,u_2,u_3\}$, $W=\{w_1,w_2,w_3\}$ with all $u$-$w$ edges ($K_{3,3}$) and $x$ adjacent to all six ($W\cup U$ written "the six vertices of $W$" in the paper: $W$ denotes the six), edge $a_5b_3$ in no triangle. 15 vertices, 32 edges.

omega: $a_i$: 5; $b_i$, $x$: 4; $u_i$, $w_i$: 3.

core values ($\kappa_2,\kappa_3,\kappa_4,\kappa_5$): $a_i$: 4, 6, 4, 1; $b_i$: 3, 3, 1, 0; $x$: 4, 3, 1, 0; $u_i,w_i$: 4, 3, 0, 0.

canonical nodes (vertex set, level interval):
- size 2: root = all 15 vertices (1..3); $W\cup\{x\}$ = 7 vertices (4..4); $A$ (4..4). Own nodes: $a_i \to A$; $b_i \to$ root; $x,u_i,w_i \to W\cup\{x\}$.
- size 3: $N=B\cup W$ = 10 vertices (1..3); $A$ (1..6). Own: $b_i,x,u_i,w_i \to N$; $a_i \to A$.
- size 4: $A$ (1..4); $B$ (1..1). Own: $a_i\to A$; $b_i,x\to B$.
- size 5: $A$ (1..1). Own: $a_i \to A$.

chains (4): $\{b_1,b_2,b_3\}$ (root, $N$, $B$); $\{a_1..a_5\}$ ($A,A,A,A$); $\{u_1..w_3\}$ ($W\cup\{x\}$, $N$); $\{x\}$ ($W\cup\{x\}$, $N$, $B$).

certified tail: $a_i$: $\sigma=2$ ($\kappa_2=4=\binom41$); $b_i$: $\sigma=2$ ($\kappa_2=3=\binom31$); $x$: $\sigma=3$ ($\kappa_2=4>\binom31=3$, $\kappa_3=3=\binom32$); $u_i,w_i$: $\sigma=4=\omega+1$ ($\kappa_2=4>2$, $\kappa_3=3>1$). Residues stored: chain $\{x\}$: $\kappa_2=4$; chain $\{u,w\}$: $\kappa_2=4,\kappa_3=3$. Three numbers.

chain ranks (lexicographic by preorder of own nodes with $\Tree{2}$ preorder root=0, $A$=1, $W\cup\{x\}$=2): $\{b_i\}<\{a_i\}<\{u,w\}<\{x\}$; labels $[0,3),[3,8),[8,14),[14,15)$; vertex order $b_1b_2b_3\,a_1..a_5\,u_1u_2u_3w_1w_2w_3\,x$.

DFS arrays, runs, entries:
- size 2: array = rank order $\{b\},\{a\},\{u,w\},\{x\}$ (root's own chain $\{b\}$ first, then children by smallest rank: $A$ then $W\cup\{x\}$); one run $[0,15)$; entries root $\to$ (0, 0), $A\to$ (0, 3), $W\cup\{x\}\to$ (0, 8), sentinel (1, 15). Community of $x$ at $(2,4)$ = $[8,15)$; at $(2,3)$ = $[0,15)$.
- size 3: preorder $N$ (0), $A$ (1); array $\{b\},\{u,w\},\{x\},\{a\}$; runs $[0,3),[8,15),[3,8)$; entries $N\to$ (0, 0), $A\to$ (2, 3), sentinel (3, 8). Community of $b_1$ at $(3,3)$: head $[0,3)$, whole run $[8,15)$, no tail.
- size 4: preorder $A$ (0) [smallest rank 1], $B$ (1) [smallest rank 0]? -- roots ordered by smallest rank: $B$ contains rank 0 ($\{b\}$), $A$ contains rank 1: order $B$ then $A$; array $\{b\},\{x\},\{a\}$; runs $[0,3),[14,15),[3,8)$; entries $B\to(0,0)$, $A\to(2,3)$, sentinel (3, 8).
- size 5: $A$ alone; array $\{a\}$; run $[3,8)$; entry $A\to(0,3)$, sentinel (1, 8).

query walks: value$(u_1,3)$: label 8, chain 2, $\omega=3,\sigma=4$, residue at size 3 = 3. value$(a_2,4)$: chain 1, $\sigma=2\le4$, $\binom43=4$. community$(x,3,2)$: chain 3, own node $N$, no parent, entries $N$=(0,0), $A$=(2,3): head $[0,3)$, run $[8,15)$: 10 vertices, 2 ranges.

## Claim -> evidence

| claim | evidence |
|---|---|
| chains = hierarchy equivalence | Theorem chains (proof) |
| chains are the coarsest lossless partition | Proposition coarsest |
| size-2 communities are one range | Theorem onerange |
| any community = head + whole runs + tail, fewest ranges, O(1) after climb | Theorem retrieval |
| values below sigma only | Theorem tail |
| replay certificate exact | Theorem replay |
| 1.28x-48.7x smaller, median 7.0x | Exp size (tables/size, size_stats) |
| locate 6-39 ns; list 0.2x-2.8x of memcpy | Exp latency (profile records) |
| build 2.8x-814x faster than per-size CND | Exp CND (prior records) |
| best s varies per query; Amazon zoom | Exp case studies (case/*.json) |
| correctness | selftest 34,075 graphs (four forms) |

## Definition dependency

nucleus -> core value -> community -> (laminar) canonical node -> own node -> chain (hierarchy equivalence)
own node -> trajectory -> chain order (preorder keys) -> aligned labels -> runs -> entries -> retrieval
omega -> floor -> certified tail -> sigma -> residue
clique tree paths -> forward counts -> replay certificate

## Section skeleton (house order)

1 Introduction (Variant B): P1 model; P2 task; P3 prose definition; Figure 1 + Example 1; Applications (3 sstitle); Existing Methods and Key Limitation; Research Question; Challenges (3); Our Idea (3); Contributions (4).
2 Preliminaries: nucleus (Def + Example), core value, community, queries, Problem Statement.
3 One Hierarchy per Size (the state of the art, fair): canonical node / own node (Def + Example), STrees layout, cost, per-size construction, closing remark with the three limitations.
4 Chains: certified tail? (no: tail belongs to values, Sec 5) -- Theorem chains (+Example, Remark), Lemma union, Lemma twins, Proposition coarsest, how many chains (numbers).
5 The Chain Index: overview (Figure layout); labels; order (Theorem onerange, Corollary kcore); runs and entries (Theorem retrieval + Example); values (Theorem tail + Example); size.
6 Queries: Algorithm Query (float + walk + example); complexity.
7 Construction: Figure pipeline; one clique tree; Theorem replay; Algorithm Build (float + walk + example); complexity + memory.
8 Experiments: Hardware / Algorithms / Datasets / Metrics; Exp-1..8; Exp-9, Exp-10 case studies.
9 Related Work. 10 Conclusion.

## Theory only (user rule, 2026-09-21, "任何实现上面的东西都不要在 paper 里面写")
Nothing implementation-level anywhere in the paper: no bitmap / rank directory / population count, no byte widths or
doubles, no cache, SIMD, prefix sums, sanitizers, compiler flags, "in-process", "loaded".  The index is defined by what
it stores (labels, chain starts as a rank structure with n + o(n) bits, records, layers with runs and entries), its size
is stated in words and bits, and the experiments keep only the measurement protocol and the numbers.  Gate: the
implementation-term grep of paper-architect/theory-not-implementation.md must hit only measurement units.

## Prose pass (2026-09-22, writing-paper-prose + five cold readers)
- Terms fixed at first use: chain(v) = number of chain starts at or before the label, minus one (a rank query); layer =
  tree of a size with its runs; traversal order = the depth-first order with children by smallest rank (a preorder of
  T_s, distinct from pre_s, used by Theorem retrieval and Section 6); "node nearest the root with top >= k on the path
  from X_s(v)" replaces "highest ancestor-or-self"; active chain = size <= clique number; K, K' for nuclei in proofs
  (N is the example's community); machines are server 1 (tods1) and server 2 (tods2).
- Corrections found by the cold reads: the listing-time ratio is STrees/ChainIndex (0.2-2.8, median 1.02); the intro
  and abstract now quote the inverse (ChainIndex takes 0.35-5.3x the copy's time, median 0.98, "about as fast");
  the skyline layout of SGL over CHAINS is 1.01-1.12x larger and 2.3-6.3x slower (the 0.99-2.28x figure was over twin
  classes); Table 5's best fixed size is now at its best level (dblp s=5: 0.147, so per-query is 1.2x, not 1.3x);
  amazon k-core F1 0.474; web-BerkStan CND build+peel is 2.0 h (2.25 h end-to-end).
- Algorithm 2 now resets active[] per size, joins only live paths with a newly active vertex, and records the value of
  each trie node per size (value[node][s]); the walk says where the residues live between sizes.
- Gates: median 17 words, over-25 15.7%, >=35 0.9%, which-tails 2.4%, so-rate 7.2%, proofs 16.3%. Page count 17.

## Symbol inventory (2026-09-22 audit; counts = uses in $...$ over sections/*.tex)
Universal, free everywhere: G, V, n, s, k, v, u (vertices), O(.).
Model (Sections 2-3, used throughout): kappa_s(v) core value (79), omega(v) clique number (34), (s,k)-nucleus (20),
T_s forest (13), X_s(v) own node (20); bare omega, sigma = a chain's clique number / certification point (records,
Algorithm 1).  Section 4: traj(v) (5); K, K' = nuclei inside proofs (Sections 3, 4, 5.4); P = a lossless partition.
Running example: A, B, N, a_i, b_i, u_i, w_j, x (only names allowed in the introduction besides s, k).
Index (Section 5): chain(v), C (number of chains), pre_s, X/Y (nodes), r_X, v_X (entry), lo(r), hi(r), R (runs),
sigma(v), c (a chain), t and m and g_s inside the Kruskal-Katona lemma, L = a largest clique (proof of 5.7a).
Queries (Section 6): the arrays of Algorithm 1 (parent, top, jump, size, entry, trajectory, residue), d = climb depth.
Construction (Section 7): T clique tree, H, Q, Z path sets, W / W_k vertex sets, d_s(W,u), pi, f_s forward count,
U_s running maximum, ell = a top in Section 3 only, the arrays of Algorithm 2 (active, node, value).
Experiments: s_max, p (sample fraction), q (query product), D (ground-truth community), F1.
Resolved clashes: M (clique vs node count), Q (nucleus vs optional set), f_s (floor vs forward count), P (partition vs
pair count), C (chains vs community), t (top vs real), k (level vs set count), N (example vs proof nucleus).
Inlined (used <= 3 times): the floor symbol, P and M of the size accounting, key(c), chain(u)=chain(w), E.

## Model name (user decision 2026-09-22): clique core, not nucleus
The paper's object is the (s,k)-core, "the k-core with cliques of s vertices in place of edges" (Definition 2.1,
macro \Nuc typesets (s,k)-core); the per-size hierarchy is the clique-core hierarchy (title, intro, research question,
contributions, conclusion).  The word nucleus appears only where the origin is credited: one sentence in the intro,
one after Definition 2.1 ("an (s,k)-core is a k-(1,s)-nucleus in the terminology of Sariyuce et al."), the CND
description (its name), and the Related Work paragraph (the (r,s) lineage and the edge/larger-clique extension).
Never present the work as "the r = 1 special case"; the (r,s) family is background, not the frame.


## Filler pass (2026-09-22, user: "写论文就是简单明了地写清楚你要做的东西就可以了")

Rule applied: a paragraph keeps only what the reader must know next; a sentence that restates an earlier block,
previews a later section, closes a paragraph by summing it up, or explains a query the index does not define is
removed. Definitions, theorems, proofs, algorithms, figures, tables and every number are unchanged (commit a174f3e,
92 sentences fewer, 17 -> 16 body+refs pages before the sentence splits, 15 body pages).

Removed or merged, by section:
- Abstract: the four mechanism sentences (labels, order, per-size tree, value storage) -> one claim sentence; construction -> one sentence.
- Introduction: "the size is the parameter" kept once (Applications); Example 1.1 states the three communities of x and points to
  Figure 1 (the counts are Example 2.2); "Its limitation is repetition" label sentence; Our idea 24 -> 16 sentences (entries,
  runs as maximal stretches, the climb, "read off the trees", "tree built as values arrive" left to Sections 5 and 7);
  Contributions: chain-index and construction items one or two sentences each with their theorems; Organization two sentences.
- Preliminaries: the three sentences on membership and level-list queries ("Neither is a separate query").
- Section 3: the preview of Sections 4-7; the top facts keep one reason each (4 sentences fewer); "The answer is one memory
  copy"; "Every size has its own tree"; Remark 3.3 (re-listed the three difficulties) -> one bridge sentence.
- Section 4: one preview sentence instead of two; the two commentary sentences after Definition 4.1; Remark 4.5 and Remark 4.8;
  "One kind of vertices ... needs no theorem"; "How many chains" merged into the Twins block (the 0.26-5.67 ratio kept);
  "Relation to other groupings" (SGL, EquiTruss) removed, Section 9 carries both.
- Section 5: the subsection roadmap and "Example 5.6 walks the same query"; the two id/label sentences merged; the k-core index
  comparison after Corollary 5.3 (Section 9 has it); the two-sentence preview of Section 6 after Example 5.6; Remark 5.10 (dblp
  1,246,933 pairs / 84,010 / 9,737 residues) moved to Exp-2 as one sentence; "Sizes with small values cost little per value".
- Section 6: the opener no longer redefines chain(v) and the entry; "Reporting the ranges costs one step per range / Expanding
  ... one step per vertex" (the Complexity block says it).
- Section 7: opener 9 -> 5 sentences ("No stage holds all values" stays in Cost only); "The removal order of that peel is the order
  tried at the next size" (said twice already); "These recorded values are the values of the chains"; the last Cost sentence.
- Experiments: "which ask for a core value"; the Exp-5 closer (constant per range vs per label, said in Exp-4); com-youtube two
  sentences merged.
- Related work: the theorem restatement after "Chains are the same idea taken across sizes".
- Conclusion: the label sentence naming the three challenges.

Kept on purpose (judged not filler): Applications (three sstitles, the paper's motivation), the three Challenges, Exp-2's
share explanation (lowest on large cliques / highest on sparse graphs), Case Study protocol sentences (how the scan reads
the index), "No step of either procedure touches the graph".
