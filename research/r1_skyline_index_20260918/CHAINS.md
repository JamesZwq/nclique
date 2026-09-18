# Hierarchy-Equivalence Chains: The Atomic Unit For r = 1 All-Size Storage

Date: 2026-09-19. Sections 1-5 are the theory and counts written before
any chain index existed (Section 4 bytes are projections); Section 6 points
to the measured S trees over chains, Section 9 to the final layout
([RESULTS_FINAL.md](RESULTS_FINAL.md)).

## 1. Motivation

SGL's gain (Zong et al., SIGMOD 2026) comes from atomic units that are
neither nested nor overlapping, so that a query touches each vertex once.
In the bipartite setting a vertex needs one unit per skyline pair (mu
units), because non-dominated cores overlap. In the (1,s)-nucleus setting
the nuclei of one size are laminar (F4) and sizes are related by the
shifted refinement (F3), so the atomic units can be a PARTITION of the
vertices: storage redundancy and retrieval redundancy vanish together.

## 2. Definition And Lemmas

Notation as in THEORY.md: X_s(v) is the own canonical node of v at size s
(the (s, kappa_s(v))-nucleus containing v), for 2 <= s <= omega(v).

Definition (chain). For omega(v) >= 2, chain(v) = (X_2(v), ..., X_omega(v)(v)),
a tuple of node ids, one per size. The chain class of v is the set of
vertices with the same tuple.

Lemma C1 (chains are exactly hierarchy equivalence). For u, v with
omega >= 2 the following are equivalent:
(a) chain(u) = chain(v);
(b) kappa_s(u) = kappa_s(v) for every s, and for every s and every k with
1 <= k <= kappa_s(u), u and v lie in the same (s, k)-nucleus.
Proof. (a) implies (b): X_s(u) = X_s(v) gives kappa_s(u) = k_hi(X_s(u)) =
kappa_s(v) by L1, and for k <= kappa_s(u) the (s, k)-nucleus containing u
is the unique ancestor-or-self of X_s(u) whose interval contains k, the
same node for v. (b) implies (a): with kappa_s(u) = kappa_s(v) = k, the
(s, k)-nucleus containing u is X_s(u) by definition, contains v, and is
therefore the (s, kappa_s(v))-nucleus containing v, that is X_s(v).

Lemma C2 (every nucleus is a union of chain classes). If u lies in the
(s, k)-nucleus N, k >= 1, and chain(w) = chain(u), then w lies in N.
Proof. kappa_s(w) = kappa_s(u) >= k, and X_s(w) = X_s(u) is a subset of N
by F4 because its top level kappa_s(u) is at least k.

Lemma C3 (chains coarsen twins). N[u] = N[w] implies chain(u) = chain(w).
Proof. F5: twins share every core value and every nucleus; apply C1.

Consequently a canonical node's own classes, its DFS slices and the
values of Block D can all be taken per chain class: every quantity the
index stores per vertex is constant on a chain class (C1), and every
query answer is a union of chain classes (C2).

## 3. Counts On The Five Local Inputs

`chains.cpp` (builds the trees of stage 1, groups vertices by their node
tuple, checks that kappa agrees on every size inside every class):

| Graph | n | twin classes | chains | canonical nodes | largest chain | (v,s) pairs -> (chain,s) pairs | residue cells -> per chain |
|---|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 5,242 | 4,105 | 715 | 1,517 | 796 | 20,045 -> 2,079 | 1,058 -> 176 |
| ca-HepPh | 12,008 | 9,559 | 1,135 | 3,538 | 1,699 | 172,963 -> 8,468 | 6,004 -> 1,731 |
| com-dblp | 317,080 | 256,983 | 13,459 | 33,979 | 60,398 | 1,246,933 -> 70,208 | 84,010 -> 9,737 |
| web-Stanford | 281,903 | 270,303 | 17,963 | 53,239 | 25,928 | 1,607,989 -> 151,364 | 872,260 -> 113,583 |
| amazon0302 | 262,111 | 253,766 | 40,867 | 65,134 | 12,148 | 824,991 -> 143,567 | 355,255 -> 57,614 |

The vertex axis collapses 6.4x (amazon) to 23.6x (dblp); twin classes
gave 1.03x to 1.28x. On all five inputs the number of chains is below
the number of canonical nodes; this is an observation, not a theorem.

## 4. Projected Bytes For S Trees Over Chains

Layout: the stage-2 baseline with chain classes in place of twin classes.
Per canonical node W + 12 bytes; per (chain, s) pair 8 bytes (node id and
DFS array entry); per chain 2 bytes (omega, sigma) and W bytes per residue
size; per-chain offsets 4 bytes; the vertex-to-chain map either as one
chain id per vertex plus a chain-to-vertex CSR (4 + 4 bytes per vertex),
or, with vertex labels aligned so that every chain is one id range, as a
bitmap with rank (1 bit per vertex) and one range per chain.

| Graph | S trees per vertex, with values | chains + id map | ratio | chains + bitmap | ratio |
|---|---:|---:|---:|---:|---:|
| ca-GrQc | 251,592 | 100,338 | 2.51 | 59,058 | 4.26 |
| ca-HepPh | 1,851,592 | 390,774 | 4.74 | 296,211 | 6.25 |
| com-dblp | 15,441,844 | 4,393,946 | 3.51 | 1,896,941 | 8.14 |
| web-Stanford | 23,725,810 | 5,691,074 | 4.17 | 3,471,088 | 6.84 |
| amazon0302 | 13,365,766 | 5,581,166 | 2.39 | 3,517,042 | 3.80 |

Against the naive representation (per-vertex S trees plus one explicit
value per active (v, s) pair) the aligned layout is 17x smaller on dblp
(32.9 MB to 1.9 MB) and 23x on HepPh (6.8 MB to 0.3 MB).

Queries keep the S-tree form: a community is one contiguous slice of
chain ids (with aligned labels, a list of vertex-id ranges), membership
is one rank plus one interval test, values are one lookup per chain.

## 5. Caveats

- Projections from counts; no chain index has been built or timed.
- The bitmap variant needs the index to use aligned vertex labels; the
  permutation from input labels costs 4 bytes per vertex unless the graph
  is stored in that order.
- chains <= canonical nodes holds on the five inputs without proof.
- The skyline dedup of stage 2 can still be applied on top of chains,
  but the (chain, s) pairs are now a minor block; it is not needed.

## 6. Measured (2026-09-19)

Done: [RESULTS_CHAINS.md](RESULTS_CHAINS.md). S trees over chains are
2.5x to 4.5x smaller than S trees over twins as measured (3.9x to 8.3x
projected with aligned labels) and 1.7x to 8.6x faster on community
listing; the skyline over chains is dominated.

Final module (later the same day): [RESULTS_FINAL.md](RESULTS_FINAL.md).
Aligned labels, lexicographic ranks and run arrays (Section 9): 3.7x to
7.8x fewer bytes than per-vertex S trees, communities located in O(1)
(3-40 ns) and listed at 0.08-0.19 ns per vertex.

## 7. How Many Chains Can There Be

Upper bounds that hold: chains <= n, chains <= twin classes (C3), and
chains <= the number of distinct tuples, which is at most the product
over sizes of |T_s|. The observation chains <= sum_s |T_s| (the total
number of canonical nodes) is NOT a theorem. Recombination defeats it:
take two vertices u, u' with the same own node at size 3 but different
own nodes at size 2 (equal kappa_3 in one 3-nucleus, different kappa_2),
and two vertices w, w' with a second common size-3 node and the same two
size-2 nodes; all four share one size-4 node. That is four chains on five
nodes, and stacking the construction over t sizes gives 2^t chains on
2t + 1 nodes. Nothing in F3 forbids it: the size-2 own node only has to
lie inside the container of the size-3 own node, and containers nest.

What the counts say is that real graphs recombine little: on the five
inputs the ratio chains / nodes is 0.47, 0.32, 0.40, 0.34 and 0.63. A
vertex's core values across sizes are coherent, so the tuples cluster.
The index size is therefore bounded by n bits plus data proportional to
the number of chains, which is empirically below the number of canonical
nodes, not by a theorem in terms of the hierarchy alone.

## 8. Toward r >= 2

For r >= 2 the items are r-cliques and the nuclei partition the r-cliques
of core value at least k into s-connected components. The chain of an
r-clique is again its tuple of own nodes, and Lemmas C1 and C2 go through
verbatim (per-size laminarity and the own-node definition are all that is
used). Two instances of one pattern (the class multiset of the NSI paper)
lie in exactly the same maximal cliques, so whenever some maximal clique
of at least s vertices contains the pattern they share an s-clique and
hence a nucleus at every level, and otherwise both have value zero; so
patterns refine chains, and the r-clique-to-chain map factors through the
paper's pattern lookup. For certified patterns the paper's forest already
answers communities without per-pattern storage; chains would replace the
per-cell residue nodes. This is a theory step, not started.

## 9. The Final Layout: Aligned Labels, Lexicographic Ranks, Run Arrays

Implemented in `chain_index.hpp` (2026-09-19); measured in
[RESULTS_FINAL.md](RESULTS_FINAL.md). Everything below is about labels
and arrays; the hierarchy itself is the canonical-node forest of THEORY.md.

Definitions. Fix for every size s a preorder pre_s of T_s (any child
order). The key of a chain c is the tuple key(c) = (pre_2(X_2(c)),
pre_3(X_3(c)), ..., pre_omega(X_omega(c))). Ranks are assigned in the
lexicographic order of keys; the vertices with omega < 2 form one chain
ranked last. Labels are aligned to ranks: chain r is the label interval
[start_pos[r], start_pos[r + 1]).

Lemma C4 (size 2 is one range). For every k >= 1 the labels of a
(2, k)-nucleus N form one interval.
Proof. N is the vertex set of the subtree of T_2 rooted at some node X.
A chain c lies in N iff X_2(c) is a descendant-or-self of X (X_2(c) is
the (2, kappa_2(c))-nucleus containing c; it is inside N iff c is in N,
by laminarity). The preorder ids of that subtree form an interval, and
the chains whose first key lies in an interval are contiguous in the
lexicographic order, whatever their later keys; the inactive chain is
outside. Contiguous ranks are contiguous labels.

Lemma C5 (DFS arrays). Let A_s be the sequence of chains with omega >= s
produced by a DFS of T_s that visits, at every node, its own chains and
its child subtrees in ascending order of the smallest rank they contain
(roots in the same order). Then (i) the chains of any node's subtree are
a contiguous segment of A_s; (ii) the first chain of node x's segment is
the smallest rank in x's subtree; (iii) A_2 is the rank order.
Proof. (i) is the DFS. (ii) by induction: the first item visited at x is
the one with the smallest key, either x's smallest own rank or the child
whose subtree minimum is smallest, and the first chain emitted inside
that child is its minimum. (iii) at s = 2 the own chains of x have first
key pre_2(x), smaller than every descendant's, and the minima of the
children are ordered as their preorder ids because ranks are ordered by
the first key; so the DFS visits nodes in preorder and, within a node,
chains by rank, which is the rank order. All active chains have ranks
0 .. C - 2, so A_2 lists consecutive ranks: one run of labels.

Run array. Positions j, j + 1 of A_s lie in one run iff the two chains
have consecutive ranks (then their label intervals abut). Store the
maximal runs as (lo, hi) label pairs, and for every node x its entry:
the index r_x of the run containing the first position b_x of its
segment and the label v_x = start_pos[A_s[b_x]]; add a sentinel entry
(M, .) after the last node, M the number of runs.

Lemma C6 (retrieval). Let y be the node following x's subtree in
preorder (the sentinel if none). The labels of x's community are
[v_x, v_y) if r_x = r_y, and otherwise [v_x, hi(r_x)), the whole runs
r_x + 1 .. r_y - 1, and [lo(r_y), v_y) (empty when b_y starts a run).
The number of ranges reported is the number of maximal runs of the
segment, one at s = 2. After the climb, locating the answer is O(1).
Proof. The segment is the position interval [b_x, b_y). Runs are maximal
in the whole array, so inside the segment they can only be cut at the
two boundaries: the run holding b_x contributes from v_x on, the run
holding b_y - 1 contributes up to v_y, and the runs strictly between are
whole. Each piece is a maximal run of the segment. At s = 2 the whole
array is one run (C5 iii), so every community is one range (C4 again).

Bytes. Per node: W (top) + 4 (parent) + 4 (subtree size) + 8 (entry); the
jump pointer (4) is derived at load time and not stored. Per run: 8. The
build form stores instead 4 per node (bucket) and 4 per (chain, s) pair;
the compact form is smaller exactly when 8 runs + 4 nodes < 4 pairs.
Both forms are exact; the compact form is the file format.
