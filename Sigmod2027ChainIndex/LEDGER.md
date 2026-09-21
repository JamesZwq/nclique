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
| Exp list | Exp-1 size vs STrees; Exp-2 where the bytes go; Exp-3 design study; Exp-4 query latency (+levels figure); Exp-5 latency by clique size and answer size; Exp-6 construction; Exp-7 scalability; Exp-8 vs CND; Exp-9 case study I (ground truth); Exp-10 case study II (Amazon zoom) |
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
