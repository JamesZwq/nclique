# Hierarchy-Equivalence Chains: The Atomic Unit For r = 1 All-Size Storage

Date: 2026-09-19. Theory plus counts. No index over chains implemented yet;
the byte figures in Section 4 are projections from measured counts.

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

## 6. Next Step

Replace the twin classes by chain classes in `index.cpp` (the baseline
code does not depend on what a class is), keep the brute-force selftest,
and measure bytes and latency of S trees over chains on the five inputs.
