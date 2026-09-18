# Skyline Community Index For All-Size Clique Cores (r = 1)

Date: 2026-09-18. Theory only. No C++, no measurement, no production, paper
or default change. The derivation follows `research/PROOF_PROTOCOL.md`.

Provenance. The index adapts the SGL index of Zong, Meng, Yang, Wang, Qiu
and Zhu, "Zero-Redundancy Search for Bi-Components in Bipartite Graphs",
Proc. ACM Manag. Data 4(3), Article 252 (SIGMOD 2026). Their predecessor
value, gamma-groups, summary graph and zero-redundancy search are theirs.
What this note supplies for (1,s)-nuclei is the dominance order through
the Kruskal-Katona shadow, the level-exactness and chain lemmas that
collapse their summary graph into per-size merge trees with one
cross-size pointer per node, the certified-tail collapse, and the twin
quotient. No literature-wide novelty is claimed beyond that mapping.

## 1. Problem, Symbols, Queries

G = (V, E) is a finite simple graph, n = |V|. For s >= 2 and W subset V,
d_s(W, v) is the number of s-cliques of G[W] containing v. The core value
is kappa_s(v) = max over W containing v of min over u in W of d_s(W, u).

Core. C_s^k = {v : kappa_s(v) >= k}. It is the largest W with minimum
d_s(W, .) at least k: the union of two such sets is again such a set,
because degrees only grow under union, so the largest one exists and every
vertex with kappa_s(v) >= k lies in it.

Nucleus. Two vertices of W are s-connected in G[W] when a chain of
s-cliques of G[W] joins them, consecutive cliques sharing at least one
vertex. An (s, k)-nucleus is a maximal s-connected subset of C_s^k in
G[C_s^k]. For k >= 1 every vertex of C_s^k lies in an s-clique of
G[C_s^k], so the (s, k)-nuclei partition C_s^k. An s-clique of G[C_s^k]
that meets a nucleus lies inside it. T_s is the family of all (s, k)-nuclei
over k >= 1.

Largest clique. omega(v) is the size of a largest clique containing v.
kappa_s(v) >= 1 iff v lies in an s-clique (the clique itself is a witness),
so {s : kappa_s(v) > 0} = [2, omega(v)] when omega(v) >= 2 and is empty
otherwise.

Binomials. C(a, b) = 0 unless 0 <= b <= a.

Cascade and shadow. For k >= 1 and s >= 1 there is a unique s-cascade
k = C(a_s, s) + C(a_{s-1}, s-1) + ... + C(a_t, t) with
a_s > a_{s-1} > ... > a_t >= t >= 1. Define
sigma_s(k) = C(a_s, s-1) + C(a_{s-1}, s-2) + ... + C(a_t, t-1) and
sigma_s(0) = 0. Kruskal-Katona theorem (external): any family of k distinct
s-element sets has at least sigma_s(k) distinct (s-1)-element subsets.
sigma_s is nondecreasing in k. The index s of sigma_s is the size of the
sets whose count is k, NOT the clique size of the layer; see F2.

Queries.
- Q1 value: (v, s) -> kappa_s(v).
- Q2 community: (v, s, k) with 1 <= k <= kappa_s(v) -> the vertex set of
  the (s, k)-nucleus containing v.
- Q3 largest size: (v, k) -> max {s : kappa_s(v) >= k}.
- Q4 membership: (u, v, s, k) -> whether u lies in the (s, k)-nucleus of v.

All answers are exact. Values use the existing W-byte exact integers.

## 2. Structural Facts

F1 (support nesting). kappa_{s+1}(v) >= 1 implies kappa_s(v) >= s.
Proof. Let W witness kappa_{s+1}(v) >= 1. Every u in W lies in an
(s+1)-clique of G[W]; dropping any one of its other s vertices gives s
distinct s-cliques of G[W] through u. So min_u d_s(W, u) >= s.

F2 (shadow bound). kappa_s(v) >= sigma_s(kappa_{s+1}(v)), and
C_{s+1}^k subset C_s^{sigma_s(k)}.
Proof. Let W witness kappa_{s+1}(v) = k >= 1. Fix u in W. The (s+1)-cliques
of G[W] through u are at least k; removing u from each gives at least k
distinct s-element subsets of N(u) meet W (the links). By Kruskal-Katona
the links have at least sigma_s(k) distinct (s-1)-element subsets, and each
such subset together with u is an s-clique of G[W] through u. Hence
d_s(W, u) >= sigma_s(k) for every u in W, so kappa_s(v) >= sigma_s(k).
For the inclusion take W = C_{s+1}^k. The bound is tight on K_m:
kappa_{s+1} = C(m-1, s) is a single-term s-cascade whose shadow is
C(m-1, s-1) = kappa_s. (This is the shifted containment law verified
earlier by brute force and on real graphs with zero violations.)

F3 (refinement). Every (s+1, k)-nucleus N' lies inside exactly one
(s, sigma_s(k))-nucleus.
Proof. Vertex containment is F2. Take u, w in N' and a chain of
(s+1)-cliques of G[C_{s+1}^k] joining them, consecutive cliques sharing a
vertex x_i. Replace each clique by an s-element subset that keeps the
shared vertices with both neighbours in the chain (an (s+1)-clique has
room for two prescribed vertices when s >= 2) and, for the end cliques,
keeps u or w. Consecutive subsets still share x_i, all lie in
G[C_s^{sigma_s(k)}], so u and w are s-connected there. Uniqueness: the
(s, sigma_s(k))-nuclei partition C_s^{sigma_s(k)}.

F4 (laminarity within one size). For k' >= k every (s, k')-nucleus lies
inside exactly one (s, k)-nucleus. Proof: C_s^{k'} subset C_s^k, and a
chain inside the smaller induced subgraph is a chain inside the larger.

F5 (true twins). If N[u] = N[w] then kappa_s(u) = kappa_s(w) for every s,
and for every k <= kappa_s(u) the vertices u and w lie in the same
(s, k)-nucleus. N[u] = N[w] holds iff u and w lie in exactly the same
maximal cliques.
Proof. The transposition of u and w is an automorphism, so the values
agree. For k >= 1 take an s-clique K of G[C_s^k] through u; w is adjacent
to every vertex of K other than itself, so K together with w is a clique,
and an s-element subset containing u and w is an s-clique of G[C_s^k]
because kappa_s(w) = kappa_s(u) >= k. Equivalence of the two twin
definitions: a maximal clique through u lies in N[u] = N[w], so adding w
keeps it a clique and maximality forces w in it; conversely every
neighbour y of u lies with u in some maximal clique, which then contains
w, so y is adjacent to w or equals w, and symmetrically.

F6 (dominance order). Write sigma^{s' -> s} = sigma_s o sigma_{s+1} o ...
o sigma_{s'-1} for s' > s and the identity for s' = s. Define
(s', k') >= (s, k) iff s' >= s and sigma^{s' -> s}(k') >= k.
Then >= is a partial order on pairs with k, k' >= 1, and
(s', k') >= (s, k) implies C_{s'}^{k'} subset C_s^k and that every
(s', k')-nucleus lies inside exactly one (s, k)-nucleus.
Proof. Reflexivity is trivial. Transitivity: each sigma_t is nondecreasing,
so sigma^{s'' -> s}(k'') = sigma^{s' -> s}(sigma^{s'' -> s'}(k'')) >=
sigma^{s' -> s}(k') >= k. Antisymmetry: s' >= s and s >= s' force s = s',
then k' >= k >= k'. The containments follow by iterating F2 and F3 down
from s' to s, then F4 from level sigma^{s' -> s}(k') down to k.

F7 (trajectory and skyline). Fix v with omega(v) >= 2 and put
P_s(v) = (s, kappa_s(v)) for 2 <= s <= omega(v). Define
delta_s(v) = kappa_s(v) - sigma_s(kappa_{s+1}(v)) >= 0 (F2) for
2 <= s < omega(v), and the skyline
Sky(v) = {s in [2, omega(v)] : s = omega(v) or delta_s(v) > 0}.
For s <= omega(v) let nxt(v, s) = min {t in Sky(v) : t >= s} (it exists
because omega(v) is in Sky(v)) and prev(v, s') = max {t in Sky(v) : t < s'}
for s' in Sky(v), with prev = 0 when no such t exists.
(a) P_{s+1}(v) >= P_s(v) iff delta_s(v) = 0.
(b) For every s <= omega(v) with s' = nxt(v, s):
kappa_s(v) = sigma^{s' -> s}(kappa_{s'}(v)) and P_{s'}(v) >= P_s(v).
(c) The maximal elements of {P_s(v)} under >= are exactly the P_s(v) with
s in Sky(v).
Proof. (a) sigma_s(kappa_{s+1}) >= kappa_s iff equality holds, since the
reverse inequality is F2. (b) delta_t(v) = 0 for s <= t < s', so
kappa_t = sigma_t(kappa_{t+1}) for those t, and composition gives the
identity; dominance follows. (c) A pair at s in Sky(v) with s < omega(v)
has delta_s > 0; any s'' > s gives
sigma^{s'' -> s}(kappa_{s''}) = sigma_s(sigma^{s'' -> s+1}(kappa_{s''}))
<= sigma_s(kappa_{s+1}) < kappa_s, using iterated F2 for the inner
inequality and monotonicity of sigma_s, so no later pair dominates it;
no pair has larger s than omega(v). A pair with delta_s = 0 is dominated
by P_{s+1}(v) by (a).

F8 (certified tail). Let f_s(v) = C(omega(v)-1, s-1); then
kappa_s(v) >= f_s(v) (a largest clique is a witness). If
kappa_t(v) = f_t(v) for some t <= omega(v), then kappa_s(v) = f_s(v) for
every s in [t, omega(v)], and delta_s(v) = 0 for t <= s < omega(v).
Let sigma(v) be the least such t, and sigma(v) = omega(v) + 1 if none
exists. Consequently Sky(v) meets [sigma(v), omega(v)] only in omega(v).
Proof. Strictness lemma: for s >= 2, a >= s and k > C(a, s) one has
sigma_s(k) > C(a, s-1). Indeed the leading cascade term of k is
C(a_s, s) with a_s >= a; if a_s > a then sigma_s(k) >= C(a_s, s-1) >
C(a, s-1) because s-1 >= 1; if a_s = a then k > C(a, s) forces further
terms, each contributing at least 1 to the shadow. Now suppose
kappa_{s+1}(v) > f_{s+1}(v) = C(omega-1, s) for some s+1 in
[t+1, omega(v)]. Then kappa_s(v) >= sigma_s(kappa_{s+1}(v)) >
C(omega-1, s-1) = f_s(v). Descending from any s+1 > t with a strict value
reaches a strict value at t, contradicting kappa_t = f_t. The delta
statement is sigma_s(C(omega-1, s)) = C(omega-1, s-1).

Counterexamples to tempting variants. Same-k containment across sizes is
false: in K_8, kappa_2 = 7 < kappa_3 = 21, so at k = 10 the size-3 core is
the whole graph while the size-2 core is empty. Same-k components can
cross. Values are not monotone or unimodal in s. Every cross-size
statement in this note therefore goes through sigma, never through a
fixed k.

## 3. Canonical Nodes And Four Lemmas

Definition (canonical node). A vertex set N that is an (s, k)-nucleus for
at least one k is a canonical node of T_s. Its level set
I(N) = {k : N is an (s, k)-nucleus} is an integer interval
[k_lo(N), k_hi(N)] with k_hi(N) = min over v in N of kappa_s(v).
Proof. If N is a nucleus at k_1 < k_2 and k_1 < k < k_2, F4 puts N (as a
k_2-nucleus) inside one k-nucleus N_k and N_k inside the k_1-nucleus that
contains it, which is N; so N_k = N. Every vertex of a k-nucleus has
kappa_s >= k, so k_hi <= min. At level m = min, N lies in C_s^m; chains
joining vertices of N never leave N (an s-clique of G[C_s^k] meeting N
lies in N), so N is s-connected in G[C_s^m], and a larger s-connected set
at level m would contradict maximality at the lower level k_lo. So N is
a nucleus at m and k_hi = m.

Parent. For k_lo(N) > 1 the parent of N is the (s, k_lo(N)-1)-nucleus P
containing N. P is a canonical node with k_hi(P) = k_lo(N) - 1: if P were
still a nucleus at level k_lo(N), then P and N would be two nuclei at
that level with N inside P, forcing P = N and N a nucleus at
k_lo(N) - 1, contradicting the definition of k_lo. Parents make T_s a
forest whose roots are the nuclei at level 1. Descendants of N have
intervals strictly above I(N); ancestors strictly below. Two canonical
nodes of T_s are nested or disjoint (F4).

Definition (own node). For v with kappa_s(v) >= 1, X_s(v) is the canonical
node of T_s that contains v and has kappa_s(v) in its interval; it is the
(s, kappa_s(v))-nucleus containing v.

Lemma L1 (level exactness). k_hi(X_s(v)) = kappa_s(v).
Proof. v lies in X_s(v), so the minimum k_hi is at most kappa_s(v);
kappa_s(v) lies in the interval, so kappa_s(v) <= k_hi.

Definition (container). For a canonical node M of T_{s+1} let
A_hi(M) be the canonical node of T_s that contains M and has
sigma_s(k_hi(M)) in its interval. By F3 applied at level k_hi(M), at
which M is a nucleus, A_hi(M) exists and is unique.

Lemma L2 (chain). If delta_s(v) = 0 then X_s(v) = A_hi(X_{s+1}(v)).
Proof. By L1, k_hi(X_{s+1}(v)) = kappa_{s+1}(v), so A_hi(X_{s+1}(v)) is
the (s, sigma_s(kappa_{s+1}(v)))-nucleus containing X_{s+1}(v), hence
containing v. With delta_s(v) = 0 that level is kappa_s(v), so this
nucleus is the (s, kappa_s(v))-nucleus containing v, which is X_s(v).

Lemma L3 (membership). Let N be a canonical node of T_s and k in I(N).
A vertex u with kappa_s(u) >= 1 lies in the (s, k)-nucleus N iff X_s(u)
is N or a descendant of N.
Proof. If X_s(u) is N or a descendant, X_s(u) is a subset of N, so u is
in N. Conversely, if u is in N then kappa_s(u) >= k, and the
(s, kappa_s(u))-nucleus containing u lies inside the (s, k)-nucleus
containing u (F4), which is N; a canonical node inside N with a higher
interval is N itself or a descendant.

Lemma L4 (enumeration). Fix s, k, and a canonical node N of T_s with k in
I(N). Put Adm_s = {N} together with all descendants of N, and for t >= s
put Adm_{t+1} = {M in T_{t+1} : A_hi(M) in Adm_t}. For a vertex u with
s <= omega(u) and s' = nxt(u, s): u lies in N iff X_{s'}(u) in Adm_{s'}.
Proof. delta_t(u) = 0 for s <= t < s', so L2 iterated gives
X_s(u) = A_hi^{(s'-s)}(X_{s'}(u)). By L3, u in N iff X_s(u) in Adm_s.
By the definition of the Adm sets, X_{s'}(u) in Adm_{s'} iff
A_hi(X_{s'}(u)) in Adm_{s'-1} iff ... iff A_hi^{(s'-s)}(X_{s'}(u))
in Adm_s.

Remark. L1 is what makes one cross-size pointer per canonical node
sufficient. The SGL paper needs level-exact nodes and explicit co-nesting
edges because its nodes are keyed by the vertex pair; here every vertex
attached to a canonical node sits exactly at that node's top level, and
the container of the node at its top level is the container of each of
those vertices at their own level (L2).

## 4. The Index

Classes. Vertices are grouped into true-twin classes (F5). All per-vertex
data below is stored once per class; class(v) is one label per vertex and
the class-to-vertex lists form one CSR over V. Sky, nxt, prev, omega,
sigma and X_s are class properties.

Block A, skeleton, one per size s in [2, s_max]. The canonical nodes of
T_s numbered in DFS preorder of the forest, so that a node and its
descendants occupy the id range [id, id + size). Per node: k_hi (W bytes),
parent id, subtree size, A_hi (id in T_{s-1}; absent for s = 2), and the
offset of its bucket in Block B. Per size s >= 2 also the reverse
cross-size lists: for every node Z of T_s, the ids of the nodes M of
T_{s+1} with A_hi(M) = Z, stored as one CSR ordered by the id of Z, so
that the lists of an id range of T_s form one contiguous slice.

Block B, entries. For every class c and every s in Sky(c) one entry
(c, gamma) with gamma = prev(c, s), placed in the bucket of X_s(c) and
sorted within the bucket by gamma ascending. A bucket may be empty.

Block C, location. For every class c the ids of X_s(c) for s in Sky(c),
ascending in s; the size of a node is known from the layer it belongs
to, so s is not stored.

Block D, values. Per class: omega(c) and sigma(c) (small integers).
Variant A: kappa_s(c) for 2 <= s < sigma(c), explicit, W bytes each.
Variant B: kappa_s(c) only for s in Sky(c) with s < sigma(c).
For s >= sigma(c) the value is C(omega(c)-1, s-1) by F8 and is not
stored; for s > omega(c) it is zero.

Size. Let E = sum over classes of |Sky(c)|, N_T = sum over s of |T_s|,
and n_cls the number of classes. Blocks A to C take
N_T (W + 12) + N_T 4 + 5 E + 4 E + 8 n_cls bytes, the class maps take
4 n + 4 n_cls bytes, and Block D takes 4 n_cls + W R bytes with
R = sum over classes of (sigma(c) - 2) in variant A or the number of
residue skyline entries in variant B.

Baseline for comparison ("S trees"). Per size s: the same canonical nodes
with their intervals, plus the DFS position of every class active at s.
Community = one contiguous position interval, membership = one
comparison. Its per-class cost is sum over s of n_cls_active(s) positions,
against E entries here; the node terms coincide except for the A_hi and
reverse-list words, which only this index carries.

## 5. Queries

Q2, community (v, s, k).
1. c = class(v). Scan Location(c) for s' = nxt(c, s), the first stored
   node whose layer is at least s. If none, s > omega(c) and the answer is
   empty (kappa_s(v) = 0, so k <= kappa_s(v) fails).
2. X = X_{s'}(c). Apply X = A_hi(X) exactly s' - s times. By F7(b) every
   delta_t(c) with s <= t < s' is zero, so by L2 iterated the result is
   X_s(c), whose k_hi equals kappa_s(v) by L1.
3. N = X_s(c). While N has a parent and k <= k_hi(parent(N)), replace N by
   its parent. On exit k lies in I(N): initially k <= kappa_s(v) = k_hi(N);
   moving up keeps k <= k_hi and stops when k > k_hi(parent) = k_lo(N)-1.
4. Adm_s is the id range [id(N), id(N) + size(N)). Output every entry
   (c', gamma) in the buckets of that range. Then take the reverse lists of
   the range (one CSR slice) as Adm_{s+1}; output every entry of their
   buckets with gamma < s; take the union of the reverse lists of the nodes
   of Adm_{s+1} as Adm_{s+2}; output with the same filter; continue until
   the admissible set is empty or the last layer is reached.
5. Expand each output class through the class-to-vertex CSR.

Theorem Q2. The output is exactly the vertex set of the (s, k)-nucleus N
containing v, each vertex once.
Proof. Let u be in N and c' = class(u). Then kappa_s(u) >= k >= 1, so
s <= omega(c') and s' = nxt(c', s) exists, with prev(c', s') < s <= s'.
By L4, X_{s'}(c') lies in Adm_{s'}, so the entry (c', prev) of c' at
layer s' is in a visited bucket and passes the filter gamma < s. It is
the only entry of c' that passes: an entry at layer t >= s with
prev(c', t) < s has no skyline between s and t, so t = nxt(c', s).
Conversely, take an output entry (c', gamma) at layer t >= s with
gamma < s. Then t in Sky(c'), no skyline of c' lies in [s, t), so
t = nxt(c', s); its node lies in Adm_t, so by L4 every vertex of c' lies
in N. At layer s every entry of a bucket in the range passes the filter
automatically, since prev(c', s) < s for all entries at s. Finally F5
puts all vertices of a class into the same nucleus.

Cost. Step 1 reads at most |Sky(c)| ids; step 2 makes s' - s pointer
steps; step 3 at most the depth of T_s; step 4 visits the admissible
nodes and reads the output entries. The admissible nodes of one layer t
are disjoint subsets of N, so their number is at most |N meet C_t^1|, and
the total number of node visits is at most sum over u in N of
(omega(u) - s + 1). No value arithmetic occurs: k is compared with node
tops only. No graph access occurs.

Q4, membership (u, v, s, k). Locate N for v as in steps 1 to 3. Compute
X_s(class(u)) by steps 1 and 2 for u; if s > omega(class(u)) answer no.
Answer yes iff id(X_s(class(u))) lies in [id(N), id(N) + size(N)) (L3).
Cost: two location scans, two chains, one tree walk.

Q1, value (v, s). c = class(v). If s > omega(c): 0. If s >= sigma(c):
C(omega(c)-1, s-1) (F8), one binomial. Otherwise, variant A returns the
stored residue value. Variant B finds s' = nxt(c, s), converts the stored
kappa_{s'}(c) to its (s'-1)-cascade (O(len) binary searches, each a
W-byte binomial evaluation), applies s' - s shadow steps, each a shift of
the coefficient sequence followed by re-canonicalisation of a possible
trailing C(a_1, 0) = 1 term, and evaluates the final cascade at size s.
By F7(b) the result is kappa_s(v). Variant B costs microseconds where
variant A costs nanoseconds; the cascade route removes the repeated
big-integer shadow evaluations that made the earlier anchor-plus-delta
prototype slow, but it does not reach constant time.

Q3, largest size (v, k). For s in [sigma(c), omega(c)] the value is
C(omega(c)-1, s-1), unimodal in s, so the largest s with value at least k
on that range is found by a scan from omega(c) downward or a binary
search on the decreasing side. For s < sigma(c) the values are not
monotone (F8 remark), so the residue is scanned. Cost O(omega(c)).

## 6. Construction

Input. The all-size engine (the terminal replay solver, or the sweep
engine) produces, for s = 2, 3, ..., the vector kappa_s, the peel order
of layer s (nondecreasing kappa_s), and the path index whose rows are
(H, Q, X) with the family H union T union (at most one x in X), T subset
Q, where H union Q union {x} is a clique for each x and X is independent
(the terminal representation; plain rows have X empty).

Step 1, connectivity of layer s. Process vertices in reverse peel order,
grouped by level k from the largest downward. Activating v: for every
row L containing v (reverse CSR), update the row's counters of active
holds, active optional members q' and active choices z'. Call the row
s-live when all its holds are active and h + q' + [z' >= 1] >= s. When a
row becomes s-live, union all its active members; afterwards union every
newly activated member with the row's representative.
Correctness. (i) Every union joins two vertices that lie in a common
s-clique of the active induced subgraph: H union Q' is a clique, every
x in X' is adjacent to all of it, and the live condition guarantees an
s-element clique through any prescribed member together with H, so all
members are chained through H. (ii) Every s-clique K of the active
subgraph is H union T union (at most one x) for exactly one row, all
holds of that row lie in K and are active, |K| = s makes the row s-live,
and the members of K are unioned. Hence after finishing level k the
union-find components are exactly the (s, k)-nuclei.
Canonical nodes. After finishing level k, every root whose component
changed during level k (a level-k vertex joined it, possibly merging
older components) receives a new node with k_hi = k, whose children are
the nodes of the merged components (their intervals start at k+1);
unchanged roots keep their node. The node of v's root after level
kappa_s(v) is X_s(v) (L1: its top is kappa_s(v)); record it as
leaf_s[v]. Number the nodes in DFS preorder and compute subtree sizes.
Cost O(sum over rows valid at s of |L| times alpha(n)) plus
O(n_active(s)).

Step 2, skyline of layer s, run with one layer of lag so that
kappa_{s+1} is known. For each class c active at s:
delta_s(c) = kappa_s(c) - sigma_s(kappa_{s+1}(c)). Shortcuts: if
s >= sigma(c) (certified, F8) then delta = 0 without arithmetic; if
kappa_{s+1}(c) = 0 then s = omega(c) and s is in Sky(c). Otherwise
compute the s-cascade of kappa_{s+1}(c) (O(s) binary searches on
binomials) and its shadow; cache by value, since many classes share
values. For s in Sky(c): append the entry (c, prev(c, s)) to the bucket
of leaf_s[c], append leaf_s[c] to Location(c), set prev(c, .) = s.
sigma(c) is detected as the first s with kappa_s(c) = C(omega(c)-1, s-1)
and stays by F8.

Step 3, cross-size links, after T_{s+1} is built and numbered. For each
node M of T_{s+1} take its creating vertex u (kappa_{s+1}(u) = k_hi(M),
u in M) and l = sigma_s(k_hi(M)) (one cascade per node). Start at
leaf_s[u] = X_s(u), whose top kappa_s(u) is at least l by F2, and walk up
while the parent's top is at least l. The node reached contains u, has l
in its interval, and by F3 is the (s, l)-nucleus containing M, that is
A_hi(M). If delta_s(u) = 0 no walk is needed (L2). Then build the
reverse lists of T_s from the A_hi values of T_{s+1}.

Working memory. leaf arrays for two consecutive layers, one union-find
over n, the row counters, and the engine's own state. The path index is
read once per layer.

Total cost. The all-size decomposition, plus connectivity
O(sum over s of sum over rows valid at s of |L| alpha(n)), plus skyline
detection O(sum over s of n_cls_active(s)) cascade computations minus the
certified shortcuts, plus cross links O(N_T (depth + cascade)). The
connectivity term has the same order as the peel's own target work.

## 7. Comparison, Break-Even, Build Test

Against S trees. Both keep the canonical nodes. S trees store one
position per active (class, size) pair; this index stores one entry per
skyline pair and adds two words per canonical node and one gamma byte
per entry. In words, the saving is
sum over s of n_cls_active(s) - E, the added cost is 2 N_T + E/4.
With per-class averages avg_traj (mean trajectory length) and mu (mean
skyline length), the membership part shrinks by the factor avg_traj / mu.
From the earlier saturation measurement (per vertex, not per class, and
excluding ca-HepPh whose high sizes overflowed there) the estimated
avg_traj - mu is 2.7 on GrQc, 3.2 on CondMat, 9.0 on AstroPh, 2.7 on
com-dblp, 1.4 on amazon0302 and 0.4 on cit-HepPh. N_T is unmeasured.

Build test. After the decomposition, E and N_T are counted in one pass.
Build this index only when
n_cls (avg_traj - mu) > 2 N_T + mu n_cls / 4;
otherwise keep S trees. The test is exact for the membership blocks and
ignores only the values, which are common to both designs.

Against the maximal-clique forest. The forest keys nodes by maximal
cliques and needs the class-to-clique lists, whose total length is
sum over v of |M(v)|; on cit-HepPh that is 100 per vertex and the forest
loses to S trees by 27x. This index never touches maximal cliques: in
the certified region a class has a single entry at its top size, the
canonical nodes there are the counterpart of the forest's clique
components, and witness joins need no rule because connectivity is
computed per size.

Query time. Community listing is output-bound in both S trees and this
index, but S trees read one contiguous interval while this index follows
pointers through admissible nodes and copies bucket prefixes. Expect a
constant-factor loss, not an asymptotic one. Membership is one comparison
in S trees and two chains plus a walk here.

## 8. Boundary Cases

- omega(v) <= 1: no entries, no location, kappa_s(v) = 0 for all s >= 2.
- omega(v) = 2: Sky(v) = {2}; certified iff kappa_2(v) = 1.
- s = 2: ordinary cores; T_2 is the k-core component hierarchy; sigma_2
  takes 2-cascades because kappa_3 counts triangles whose links are
  2-element sets.
- kappa_omega(v) >= 2: sigma(v) = omega(v) + 1, nothing is certified,
  omega(v) is still the top entry and its value is stored.
- Equal values: classes with equal kappa in one nucleus share the node;
  distinct classes at equal level in different components have distinct
  nodes.
- Level 1 has several nuclei: T_s is a forest; DFS numbering per tree.
- Empty buckets: canonical nodes without skyline entries are kept; they
  are containers on chains and connectors between siblings. Removing
  them requires sibling links (the co-nested edges of SGL); without those
  links a query whose target is the removed node stops at one child and
  loses the others. This is the one place where SGL's extra edge type is
  genuinely needed if node count must be reduced.
- Overflow: values are W-byte checked integers; cascades use small
  integers; sigma comparisons use W-byte binomials.
- Layers beyond the largest omega are empty; s_max = max omega.

## 9. Adverse Cases And Limits

- Low saturation (cit-HepPh, 13 percent): E approaches the number of
  active pairs, the index is larger than S trees by the added words, and
  the build test rejects it.
- Fragmented hierarchies at small s inflate N_T for both designs.
- Deep merge trees make the walk in Q2 step 3 and in construction step 3
  long; level-ancestor jump pointers cost one word per node if needed.
- Certified vertices queried at small s follow chains of length up to
  omega(v) - s; bounded by 239 on the local inputs.
- Variant B value queries are microseconds, not nanoseconds; graphs with
  many uncertified classes make variant A as large as the dense table.
- Nothing here reduces build time; the decomposition still computes
  every layer, and the certified tail only saves skyline arithmetic.

## 10. What Is New Here And What Is Not

From SGL: the predecessor value, the one-entry-per-vertex retrieval, the
summary-graph search. From earlier work in this repository: F2 (shifted
containment, verified), the saturation statistics, the observation that
anchor-plus-delta value reconstruction is slow. New in this note, for
this problem: the partial order F6 built from iterated shadows; L1 and L2,
which reduce SGL's summary graph to per-size merge trees with one
container pointer per node and no co-nesting edges; L4 and Theorem Q2,
which make community search free of value arithmetic; the certified-tail
collapse to a single entry (F8 with F7); the twin quotient (F5); and the
exact break-even against S trees. Not claimed: any speedup, any measured
size, literature novelty of the retrieval idea, or applicability beyond
r = 1.

## 11. Measurement Plan Before Any Prototype

On the five local full-range inputs, using the terminal solver's core
matrix and path index, count: (1) E per class and per vertex, and
avg_traj - mu; (2) N_T, depths, and the number of empty canonical nodes;
(3) the share of residue cells (s < sigma) with delta = 0, which decides
variant A against B; (4) the number of active (class, size) pairs, the
S-tree cost. Then apply the build test of Section 7. These are counts,
not an index implementation.

## 12. Proof Gate

Dated 2026-09-18. Proved above: F1 to F8, the canonical-node interval
and parent facts, L1 to L4, Theorem Q2, the correctness of Q1, Q3, Q4,
and of construction steps 1 to 3. External dependencies: the
Kruskal-Katona theorem and monotonicity of its shadow function; the
partition property of the path index (every clique has exactly one row
with H inside it and the clique inside the row's family) and the clique
property of rows; the earlier verification of F2 on real graphs. No
performance or size claim is made; the benefit gate of the proof
protocol is open until the counts of Section 11 exist.
