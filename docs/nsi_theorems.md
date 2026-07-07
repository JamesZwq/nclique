# The Nucleus Spectrum Index: Theorem Set (working draft, feeds the paper's theory section)

Status 2026-07-07. Every statement here is either PROVEN below (self-contained proofs) or
explicitly marked as a sketch/conjecture. Empirical validation pointers reference SigmodPlus
sections. Notation: G a simple graph; an l-clique is a set of l pairwise-adjacent vertices.

## D1. Definitions

Fix r < s. For a family F of s-cliques of G and an r-clique R', the F-support is
sup_F(R') = |{S in F : R' subset of S}|; R' is covered by F if sup_F(R') >= 1.
The (r,s)-core number of an r-clique R is

  kappa_{r,s}(R) = max { k : there is a family F of s-cliques with R covered by F and
                          min over all r-cliques R' covered by F of sup_F(R') >= k },

and kappa_{r,s}(R) = 0 if R is in no s-clique. (Max-min form; P9 states equivalence with the
standard peeling computation.) c(R) denotes the size of the largest clique of G containing R;
equivalently the size of the largest maximal clique (region) containing R.
The spectrum of R is the table s -> kappa_{r,s}(R) (and across r, the (r,s) surface).

## T1. Clique lower bound

THEOREM. kappa_{r,s}(R) >= C(c(R) - r, s - r).

PROOF. If c(R) < s both sides handle trivially (right side is 0). Otherwise let K be a clique
containing R with |K| = c(R) and take F = all s-cliques inside K. Every r-clique R' covered by
F lies in K and has sup_F(R') = C(c(R) - r, s - r) exactly, and F contains an s-clique
containing R. The max-min definition gives the bound. QED

## T2. Kruskal-Katona shadow bound (the s-direction)

THEOREM. For s >= r + 2: kappa_{r,s-1}(R) >= g_a(kappa_{r,s}(R)), where a = s - r and g_a is
the Lovasz form of the Kruskal-Katona lower-shadow function: for m >= 0 write m = C(x, a) with
real x >= a; then g_a(m) = C(x, a - 1).

PROOF. Let k = kappa_{r,s}(R), realized by the family F. Let F' = all (s-1)-cliques contained
in members of F. Every r-clique R' covered by F' lies inside some member of F. Consider the
extension family E(R') = { S \ R' : S in F, S contains R' }, an (s-r)-uniform family with
|E(R')| = sup_F(R') >= k. For every T in the lower shadow of E(R') (T an (s-1-r)-subset of a
member of E(R')), R' union T is an (s-1)-clique inside a member of F, hence in F', so
sup_{F'}(R') >= |shadow(E(R'))| >= g_a(|E(R')|) >= g_a(k) by Kruskal-Katona (Lovasz form) and
monotonicity of g_a. F' contains an (s-1)-clique containing R (drop one vertex of S \ R for
any S in F containing R; s - r >= 2 makes this possible). The max-min definition gives the
claim. QED

REMARK (zero slack at binomials). If k = C(y, a) for an INTEGER y then g_a(k) = C(y, a - 1)
exactly; on a clique K_c the chain of clique values C(c-r, s-r) is g-tight at every s.

## T3. The s-chain certificate (the sweep engine's workhorse, SigmodPlus 130)

THEOREM. Let c = c(R) and s >= r + 2. If kappa_{r,s-1}(R) = C(c - r, s - 1 - r), then
kappa_{r,s}(R) = C(c - r, s - r).

PROOF. Lower bound: T1. Upper bound: if c < s then no s-clique contains R (else c >= s) and
kappa_{r,s}(R) = 0 = C(c - r, s - r). Otherwise c >= s; suppose kappa_{r,s}(R) >= C(c-r, s-r) + 1.
By T2, kappa_{r,s-1}(R) >= g_a(C(c-r, a) + 1) with a = s - r >= 2. The Lovasz x(m) solving
C(x, a) = m is strictly increasing in m, and x -> C(x, a-1) is strictly increasing for
x > a - 1; at the integer point x = c - r (>= a since c >= s) we get
g_a(C(c-r, a) + 1) > g_a(C(c-r, a)) = C(c-r, a-1) = kappa_{r,s-1}(R), a contradiction. QED

COROLLARY (absorbing chain). If R is certified at cell s (kappa equals the clique value), the
hypothesis of T3 holds at s + 1 with the same c, so R stays certified at every larger s.
This is the one-build, ascending-s sweep: cold-peel the boundary cell, then each cell certifies
from its predecessor by one integer comparison per pattern; only the residue is peeled.
Validated bit-exact: SigmodPlus 130 (99.96-100 percent certified on clique-dominated graphs,
including the former loss cells; 47-51 percent on soc-Epinions).

## T4. Subclique bound (the r-direction; valid but not tight)

THEOREM. For an (r+1)-clique R+: kappa_{r+1,s}(R+) <= min over r-subcliques R of R+ of
kappa_{r,s}(R). More generally every r-clique covered by a family realizing
k = kappa_{r+1,s}(R+) has (r,s)-core at least k.

PROOF. Let F realize k for R+ at level (r+1, s). Let R' be any r-clique covered by F, say
R' subset of S in F. Pick u in S \ R'; then R' union {u} is an (r+1)-clique covered by F, so it
lies in at least k members of F, each of which also contains R'. Hence sup_F(R') >= k, and F
witnesses kappa_{r,s}(R') >= k in the max-min definition. Applied to the r-subcliques of R+
(covered because R+ is), this gives the claim. QED

REMARK (no squeeze). On a c-clique the bound reads C(c-r, s-r) against a true value of
C(c-r-1, s-r-1): strict slack, so T4 alone cannot certify in the r-direction (measured
consequence: SigmodPlus 131-132, the diagonal band).

## T5. The diagonal +1 theorem (NEW; generalizes the classical truss-core inequality)

THEOREM. For every (r+1)-clique R+:

  kappa_{r+1,r+2}(R+) <= ( min over r-subcliques R of R+ of kappa_{r,r+1}(R) ) - 1.

More generally, every r-clique covered by a family realizing k = kappa_{r+1,r+2}(R+) has
kappa_{r,r+1} >= k + 1. (The case r = 1 is the classical truss-core relation
"truss number <= core number - 1".)

PROOF. Let F be a family of (r+2)-cliques realizing k for R+ at level (r+1, r+2), and let
F' = all (r+1)-cliques contained in members of F. Let R' be any r-clique covered by F', and
pick W = R' union {u} in F' (u a single vertex since members of F' have r + 1 vertices). W is
an (r+1)-clique contained in a member of F, hence covered by F, so at least k distinct members
S_1, ..., S_k of F contain W; each S_i = W union {w_i} with the w_i pairwise distinct and not
in W. Now F' contains the k + 1 distinct (r+1)-cliques
  R' union {u} = W  and  R' union {w_i} (i = 1..k),
the latter because R' union {w_i} is an (r+1)-subset of the clique S_i in F. Hence
sup_{F'}(R') >= k + 1, and F' (which covers every r-subclique of R+ because R+ lies in a member
of F) witnesses kappa_{r,r+1}(R') >= k + 1. QED

REMARK (zero slack on cliques). On a c-clique: left side c - r - 1, right side (c - r) - 1:
equality. So T5, unlike T4, can SQUEEZE in the r-direction along the t = 1 diagonal:

COROLLARY (diagonal chain certificate). If min over r-subcliques R of R+ of kappa_{r,r+1}(R)
equals c(R+) - r, then kappa_{r+1,r+2}(R+) = c(R+) - r - 1 exactly (T1 floor meets the T5
ceiling). This is the diagonal analogue of T3 and the theoretical basis for a one-pass
diagonal sweep. MEASURED (SigmodPlus 133): zero violations across 156M+ 4-cliques on three
graphs, and the ceiling is an EQUALITY (kappa_{4,5} = min_sub kappa_{3,4} - 1) for 100.0
percent of 4-cliques on ca-GrQc (329k) and ca-HepPh (150.3M), 23.4 percent on soc-Epinions.

## T6. No-early-death and host-1 exactness (all s; used by the band engine at t = 1)

THEOREM. Let M be a maximal clique with |M| >= s and L_M = C(|M| - r, s - r). In the standard
peel (P9): (i) no r-clique contained in M is removed before level L_M; (ii) every r-clique R
whose ONLY maximal clique is M (a host-1 clique) has kappa_{r,s}(R) = L_M.

PROOF. (i) Induction on levels. While no r-clique inside M has been removed, every s-clique
inside M is intact, so every r-clique R' inside M has current support at least
C(|M| - r, s - r) = L_M; hence at any level lower than L_M none is removable. (ii) Every
s-clique containing a host-1 R lies in a maximal clique containing R, which must be M, so
R's initial support is exactly L_M; by (i) it is intact until level L_M, where R is removed
with support exactly L_M. QED

REMARK. (ii) is the engine's SKIP_H1 rule (previously verified empirically with 0 exceptions);
(i) at t = 1 is the no-early-death lemma that makes the band engine's region waves need no
already-dead correction (SigmodPlus 132).

## T7. Quotient invariance (patterns are core-classes)

Classes: vertices with identical maximal-clique membership. The pattern of an r-clique is the
multiset of classes of its vertices.

THEOREM. Two r-cliques with the same pattern have the same kappa_{r,s} for every s.

PROOF. It suffices to show the transposition (u v) of two same-class vertices maps cliques to
cliques (then it is an automorphism of the clique complex, hence of every s-clique hypergraph,
and same-pattern r-cliques are connected by such transpositions). Let K be a clique with
u in K, v not in K. K lies in some maximal clique M containing u; since v has the same
maximal-clique membership as u, v is in M as well, so (K \ {u}) union {v} is a subset of M
and hence a clique. Cliques containing both or neither of u, v map to themselves or
symmetrically. QED

REMARK. This is the exactness of the engine's pattern-level peel (previously brute-verified,
0 mismatches in 735k checks, SigmodPlus 110); the proof removes the empirical caveat.

## T8. Witness ownership split (t = 1; exactness of the band engine's replay)

THEOREM. Fix s = r + 1 and call (r+1)-cliques witnesses. (i) If a witness W contains an
r-subclique R' whose only maximal clique is M, then M is the only maximal clique containing W,
and W's first sub-removal happens exactly at level L_M = |M| - r; (ii) witnesses whose every
r-subclique is multi-host (contained in at least two maximal cliques) are exactly the
witnesses not of type (i), and their first sub-removal is an event of the multi-host peel.

PROOF. (i) Any maximal clique containing W contains R', so it equals M; thus W and all its
subcliques lie in M. By T6(i) nothing in M is removed before L_M, and by T6(ii) R' is removed
at L_M, so W's first sub-removal is at L_M exactly. (ii) is the complement by definition; the
first removal among multi-host subcliques is by definition a multi-host peel event, and by (i)
it cannot be preceded by a host-1 removal (any host-1 subclique would place W in type (i)). QED

REMARK. (i) is the region wave (closed-form credit, no enumeration of host-1 cliques);
(ii) is the deadW-deduplicated per-death enumeration. The two classes partition the witnesses,
so every witness death is credited exactly once: this is the exactness of diag_band.cpp,
validated bit-exact on 6/6 gates (SigmodPlus 132).

## P9. Peel = max-min, and the certification principle

PROPOSITION (standard; proof sketch). The peeling process (repeatedly remove an r-clique of
minimum current support; an s-clique survives while all its r-subcliques survive; the core of
R is the level at which R is removed, with the monotone clamp) computes exactly the max-min
quantity of D1. Sketch: (<=) when R is removed at level k, the surviving family at that moment
witnesses min-support >= k for everything still covered... (>=) any family F with min-support
k >= 1 survives intact below level k by induction, so its members are removed at levels >= k.
Consequence (certification principle): if L(R) <= kappa(R) <= U(R) are any valid bounds and
L(R) = U(R), then kappa(R) is known without peeling R; T3 and the T5 corollary are the two
zero-slack instances used by the engines.

## P10. No free spectrum (output-independence lower bound; sketch)

PROPOSITION (sketch, from the two independent derivations in SigmodPlus 128). For every omega
there are graphs on which the Theta(omega^2) cell values kappa-max over cells (r,s), r < s <=
omega, vary independently: for each cell one can splice a disjoint gadget that changes that
cell's value while leaving all other cells and all boundary supports unchanged (a graph rich
in (s-1)-cliques can be s-clique-free, so no cross-s certificate exists in general). Hence no
algorithm outputs the exact spectrum with o(omega^2) cell-distinguishing work, and the
certified fraction of T3/T5 is necessarily graph-dependent (measured: 99.96-100 percent on
clique-dominated graphs, about half on social graphs). STATUS: proof sketch; the gadget
construction should be written out before the paper claims it as a theorem.

## Summary table

  T1 clique floor            proven   every cell, every graph
  T2 KK shadow               proven   s-direction transfer, zero slack at binomials
  T3 s-chain certificate     proven   the sweep engine (SigmodPlus 130), absorbing
  T4 subclique bound         proven   r-direction, has slack, cannot certify alone
  T5 diagonal +1             proven   NEW; generalizes truss-core; zero slack on cliques
  T6 no-early-death/host-1   proven   band engine waves + SKIP_H1, now unconditional
  T7 quotient invariance     proven   pattern peel exactness, empirical caveat removed
  T8 ownership split         proven   band engine replay exactness
  P9 peel = max-min          standard sketch included
  P10 no free spectrum       sketch   write out the gadget before claiming
