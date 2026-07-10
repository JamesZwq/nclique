# Framing Recommendation: Batch Decomposition and Class-Compressed Indexing

## Executive verdict

The proposed two-pillar framing is the strongest honest framing, with one important sharpening:

> **One structural idea, two consequences.** Relevant maximal-clique profiles expose interchangeable vertex classes. The same class quotient (i) supports exact batch decomposition over a bounded \((r,s)\) plane and (ii) compresses the resulting plane into a queryable index.

This is stronger than presenting transfer theory, the CPI, replay, compression, and hierarchy as five coequal contributions. It gives the paper one intellectual center and two clearly different system consequences:

1. **Batch decomposition:** replace cell-at-a-time recomputation with family-at-a-time construction. One immutable, \(r\)-independent class-SCT is shared across the bounded plane; exact certificates infer some neighboring-cell values, and shell-order replay computes the uncertified residue.
2. **Class-compressed indexing:** replace value-at-a-time materialization with a structural representation. A pattern represents an orbit of interchangeable \(r\)-cliques, certified tails are reconstructed arithmetically, and only uncertified residue values are materialized.

The whole plane should be the paper's **headline scope**, because Universal CPI is a genuine theorem and the implemented multi-\(r\) engine is bit-exact. It must not be sold as a cheaper-than-\(r\)-by-\(r\) construction. The exact claim is that one universal structure is sufficient for the plane and yields one unified index; substantial pattern, map, boundary, and residue work remains specific to each \(r\).

Do not describe CND as the only baseline. Use two baselines, explicitly assigned to two different questions:

- **External effectiveness baseline:** CND run once per requested cell. This is the state-of-the-art comparison and supports the orders-of-magnitude result.
- **Internal sharing baseline:** the optimized fixed-\(r\) engine run once per \(r\). This isolates the value and cost of the multi-\(r\) design. The honest result is competitive construction time, not a speedup. Its benefit is consolidation into one universal, queryable plane index.

This distinction preempts the obvious reviewer charge that the paper chose a weak baseline to make cross-parameter reuse look stronger than it is.

## Recommended paper theorem and terminology

### One-sentence thesis

> **The bounded \((r,s)\)-nucleus plane can be treated as one exact indexable object: an \(r\)-independent class quotient supports every cell, cross-cell certification and shell-order replay batch their decomposition, and the same classes compress the outputs into patterns.**

A shorter sentence for the Introduction is:

> **NSI replaces cell-at-a-time decomposition with family-at-a-time construction, and value-at-a-time materialization with class-compressed indexing.**

### Terms that must be separated

Use the following terminology consistently:

- A **spectrum row** fixes \(r\) and varies \(s\).
- The **bounded nucleus plane** contains all admissible cells with \(r\in[r_{\min},r_{\max}]\) and \(s\in[\max\{r+1,S_{\min}\},S_{\max}]\).
- A **shared structural build** means maximal-clique processing, the relevant class quotient, and one immutable class-SCT are built once. It does not mean that all per-\(r\) work disappears.
- A **certified cell value** is inferred exactly. A cell is not “free” unless the measured residue is empty; even then, construction performs pattern checks.
- A pattern's multiplicity is
  \[
  m(P)=\prod_c \binom{w_c}{b_c^P}.
  \]
  The expression \(\binom{w}{r}\) is only the one-class special case.

The current fixed-\(r\) Problem Statement and the Related Work sentence that places cross-\(r\) indexing in future work must be removed. They directly contradict the new headline.

## Recommended title

> **NSI: Batch Decomposition and Class-Compressed Indexing of the \((r,s)\)-Nucleus Plane**

This title states the two contributions and the plane scope without promising a one-cell build. Delete the current phrase “at the Cost of One Cell.” It is false for the plane and only conditionally true for a row when certification leaves little or no residue.

## Paste-ready abstract

> The \((r,s)\)-nucleus exposes cohesive subgraphs at different clique orders and witness sizes, so exploratory analysis requires a family of cells rather than one decomposition. State-of-the-art algorithms, however, solve each cell from scratch, while explicitly materializing all cell values is prohibitive when a graph contains billions or trillions of \(r\)-cliques. We present NSI, an exact index that replaces cell-at-a-time processing with batch decomposition and class-compressed storage over a bounded \((r,s)\) plane. Our key observation is that relevant maximal-clique profiles partition vertices into interchangeable classes: we prove that one immutable, \(r\)-independent class-based clique tree exactly counts support for every cell in the plane, and that peeling weighted class patterns is equivalent to peeling individual \(r\)-cliques. Exact chain certificates, conditional diagonal certificates, and shell-order replay infer certified values and compute the uncertified residue without changing any from-scratch core value. The same classes compress the output: one pattern represents \(\prod_c\binom{w_c}{b_c}\) \(r\)-cliques, certified tails require only a boundary value and arithmetic, and only residue values are materialized. Against the state-of-the-art per-cell method CND, batch construction is up to \(1{,}600\times\) faster on clique-structured graphs and completes whole spectra on five graphs where CND cannot finish one cell within the memory budget; NSI stores the resulting family at \(0.001\)--\(2\) bytes per represented \(r\)-clique and answers in-memory queries in nanoseconds. Building across \(r\) is competitive with running our optimized fixed-\(r\) engine once per \(r\), rather than cheaper; its value is one exact, reusable index over the plane and the elimination of repeated per-cell decomposition within each indexed family.

Before submission, replace “competitive” in the last sentence with a measured range, for example “\([x,y]\times\) the summed fixed-\(r\) time,” once the final uncontended whole-plane benchmark is frozen. If the final range is noisy, retain “competitive” and give the exact per-dataset values in the evaluation.

## Paste-ready contribution list

Use four contributions, ordered from foundation to consequence to evidence. Do not make the hierarchy or case studies a fifth headline contribution; they are derived uses of the index.

### 1. Universal class representation

> **A universal class representation for the bounded nucleus plane.** We prove that relevant maximal-clique profiles induce exact symmetry classes and construct one immutable class-based succinct clique tree, independent of \(r\), whose disjoint leaves exactly count support for every admissible \((r,s)\) in a bounded plane. We further prove that weighted-pattern peeling is equivalent to individual \(r\)-clique peeling, so one pattern can safely represent \(\prod_c\binom{w_c}{b_c}\) interchangeable \(r\)-cliques.

### 2. Exact batch decomposition

> **An exact family-at-a-time decomposition algorithm.** We design a batch constructor that reuses the universal class tree across all indexed cells. Within each spectrum row, the chain certificate identifies values forced by the preceding cell; between adjacent row boundaries, the diagonal bound certifies a value only when it meets the clique floor; all remaining patterns are computed by shell-order replay. We prove replay exact even when a certified pattern is removed at its known core level while its current support is larger, closing the gap between certification and from-scratch peeling.

### 3. Class-compressed plane index

> **A class-compressed index for exact plane queries.** We introduce NSI, which stores the shared class structure once and stores, for each \(r\), pattern metadata, a boundary value, the first certified cell, and residue dictionaries. Certified values are reconstructed by binomial arithmetic rather than materialized. The resulting index supports exact point and spectrum-row queries without rerunning a decomposition and has size proportional to the shared structure, pattern columns, and uncertified residues rather than to the number of represented \(r\)-cliques times the number of cells.

### 4. Evaluation with both external and internal baselines

> **An evidence-separated experimental evaluation.** Against CND, the state-of-the-art cell-at-a-time method, our batch construction is up to \(1{,}600\times\) faster and up to \(140\times\) more memory-efficient on clique-structured graphs, and it completes whole spectra on five graphs where CND cannot finish one cell within the \(300\) GB budget. NSI occupies \(0.001\)--\(2\) bytes per represented \(r\)-clique and answers in-memory queries in nanoseconds. Against the stronger internal control of running our optimized fixed-\(r\) engine once per \(r\), whole-plane construction is competitive rather than faster, showing that the multi-\(r\) contribution is a unified shared index—not the elimination of inherently \(r\)-specific pattern and peeling work.

If the final paper has space, the first two items may be split into “Universal CPI” and “Exact certification and replay,” yielding five items. Four is cleaner because it preserves the two-pillar story.

## How to present the multi-\(r\) plane honestly

### The exact message

Use this paragraph in the Introduction immediately after presenting the two pillars:

> **Scope across the plane.** The class tree is universal but the entire algorithm is not \(r\)-free. We build the maximal-clique profile classes and one immutable class-SCT once for the bounded plane. For each \(r\), NSI still constructs an \(r\)-specific pattern column and may perform boundary and residue peeling. Consequently, our claim is not that plane construction is asymptotically cheaper than running an optimized fixed-\(r\) constructor once per \(r\). Our claim is that one exact counting structure suffices for all cells, that certified cells avoid from-scratch decomposition, and that the result is one compact index answering any indexed \((r,s)\) query.

Use this formal cost statement in the construction section:

> Let \(T_{\mathrm{shared}}\) be the cost of maximal-clique enumeration, profile-class construction, and the universal class-SCT. For each order \(r\), let \(T_{\mathrm{col}}(r)\) include mergeability classification, pattern and pattern-to-leaf maps, the boundary computation, \(O(|\mathcal P_r|)\) certification checks per higher cell, and exact residue replay. Plane construction costs
> \[
> T_{\mathrm{shared}}+\sum_{r=r_{\min}}^{r_{\max}}T_{\mathrm{col}}(r).
> \]
> If every higher-cell pattern certifies, that cell requires \(O(|\mathcal P_r|)\) checks; in the worst case, residue replay can approach a full peel. Universal CPI proves exact coverage, not a universal speedup.

Use this result sentence in the evaluation:

> The plane builder is competitive with the sum of independently optimized fixed-\(r\) sweeps, but it does not beat that sum consistently because each order retains its own pattern and peeling state. The shared design instead produces one universal structure and one queryable plane index, while avoiding CND's repeated cell-at-a-time computation.

### Where the plane should live

The plane belongs in the title, abstract, thesis, problem statement, contribution list, and conclusion. It should not first appear as a late extension.

Pedagogically, however, the algorithm should still be explained in two stages:

1. Explain a fixed-\(r\) spectrum row: boundary peel, chain certification, and residue replay.
2. Add a main-body subsection titled **“From One Spectrum Row to the Bounded Plane.”** State Universal CPI, show what is shared, show what remains per \(r\), explain the conditional diagonal certificate, and give the outer \(r\)-loop pseudocode.

This organization makes the plane a headline contribution without forcing readers to understand two-dimensional control flow before they understand one row.

The index section must then show a plane layout, not the current fixed-\(r\) NSI layout:

- shared class map, class weights, region profiles, region sizes, and one universal class-SCT;
- a directory of \(r\)-columns;
- per-column pattern dictionary, multiplicities, clique floors, boundary cores, first-certified cells, and residue dictionaries;
- a query algorithm whose input includes \(r\), \(s\), and the vertices of the queried \(r\)-clique.

If the implementation does not yet serialize and query this multi-\(r\) layout, either implement and benchmark it before submission or weaken “plane index” to “plane construction with fixed-\(r\) index columns.” A hostile reviewer will notice the difference.

### Required plane experiment

Add a dedicated research question: **“What does sharing across \(r\) save, and what remains \(r\)-specific?”** Report, on several easy and hard graphs:

| Comparison | Purpose | Required reporting |
|---|---|---|
| Plane NSI vs CND over every indexed cell | External practical value | total time, peak RSS, completed cells, identical budgets |
| Plane NSI vs sum of optimized fixed-\(r\) sweeps | Honest cross-\(r\) control | total time ratio and peak RSS; do not expect a speedup |
| Shared build vs per-\(r\) columns | Cost attribution | MCE/classes/tree, patterns, maps, boundaries, residue replay |
| With vs without diagonal certification | Cross-row ablation | certified boundary fraction, residue size, time |
| Plane index size and query time | Index claim | absolute MB, shared bytes counted once, per-column bytes, point and row-query latency |

Report at least one graph where plane construction loses to summed fixed-\(r\), if one exists. That boundary result will make the large CND wins more believable.

## Claim--evidence ledger

| Claim | Evidence that must appear in the paper |
|---|---|
| One class-SCT covers every bounded \((r,s)\) cell | Universal CPI theorem plus the disjoint exact-cover theorem for `Gen` |
| Pattern computation equals individual \(r\)-clique computation | Relevant-class symmetry and weighted-pattern peeling equivalence |
| Scheduled certified deaths preserve residue cores | Shell-order replay theorem and complete replay algorithm |
| Certification along \(s\) is exact and absorbing | Clique floor, Kruskal--Katona shadow bound, and chain certificate |
| A cross-\(r\) boundary value is exact | Only the conditional squeeze \(U=L\); the diagonal theorem alone is an upper bound |
| The index is compact | Formal size expression plus absolute sizes and bytes per represented clique |
| Queries take nanoseconds | Defined in-memory query protocol, latency distribution, and clear inclusion/exclusion of pattern formation and clique validation |
| Batch construction beats the state of the art | End-to-end same-machine CND comparison over all requested cells |
| Whole-plane sharing is worthwhile but not cheaper | Plane-vs-summed-fixed-\(r\) experiment and phase breakdown |
| Implementation is faithful | Bit-exact comparison with the trusted fixed-\(r\) engine; this validates implementation but does not replace the proofs |

## Hostile reviewer attacks and exact preemptions

### Attack 1: “The title says plane, but the formal problem and algorithm fix \(r\).”

**Preemption:** Redefine the problem over \([r_{\min},r_{\max}]\times[S_{\min},S_{\max}]\), with \(r<s\). Include Universal CPI, the outer \(r\)-loop, the plane index layout, and plane experiments in the main paper. Remove the Related Work claim that cross-\(r\) indexing is future work.

Exact sentence:

> Given \(G\), an order interval \([r_{\min},r_{\max}]\), and a witness-size cap \(S_{\max}\), we build an exact index for every admissible cell \((r,s)\) with \(r_{\min}\le r\le r_{\max}\) and \(r<s\le S_{\max}\).

### Attack 2: “One universal tree does not mean one cheap plane computation.”

**Preemption:** State the additive cost formula and list all per-\(r\) state. Put “competitive with summed fixed-\(r\)” in the abstract and evaluation.

Exact sentence:

> Universal CPI is a sufficiency theorem for counting, not a claim that all \(r\)-specific work vanishes.

### Attack 3: “CND is a straw baseline because the authors can repeat their own optimized fixed-\(r\) engine.”

**Preemption:** Show both comparisons and explain their roles before presenting results.

Exact sentence:

> CND measures improvement over the state of the art; repeated fixed-\(r\) construction measures the marginal effect of sharing across \(r\).

### Attack 4: “The paper claims every additional cell is free.”

**Preemption:** Replace “free” by a conditional cost statement.

Exact sentence:

> A higher cell costs \(O(|\mathcal P_r|)\) certification checks plus exact peeling of its uncertified residue; it is near-free only when the measured residue is empty or small.

### Attack 5: “The diagonal theorem is only an upper bound, not a cross-row recurrence.”

**Preemption:** State that it certifies only when the upper bound meets the clique floor; otherwise peel the boundary residue. Never schedule a death at the upper bound alone.

Exact sentence:

> The diagonal bound certifies \(Q\) only when \(U_r(Q)=L_r(Q)\); when \(U_r(Q)>L_r(Q)\), NSI makes no inference and invokes exact replay.

### Attack 6: “Pattern peeling silently changes support multiplicities.”

**Preemption:** Integrate the weighted-pattern equivalence theorem, not merely an appendix reference. Emphasize per-representative support and one-time witness death.

Exact sentence:

> Pattern multiplicity controls how many outputs the pattern represents; it is not a multiplier in support or decrement, because a witness containing several representatives dies once.

### Attack 7: “Replay schedules a certified item before ordinary peeling would remove it.”

**Preemption:** Replace the current Level Invariance argument with Shell-Order Replay. The present Level Invariance theorem does not prove this stronger operation.

Exact sentence:

> Shell-order replay may delete a certified core-\(k\) pattern at level \(k\) even when its current support exceeds \(k\), because the support protecting the \((k+1)\)-superlevel set is internal to that set.

### Attack 8: “Class symmetry is asserted too broadly.”

**Preemption:** Define profiles using maximal cliques of size at least \(S_{\min}\), treat empty profiles carefully, and retain the caveat that the induced permutations preserve relevant clique families, not necessarily edges found only in smaller maximal cliques.

### Attack 9: “Compression is dataset-specific and can be exponential in the worst case.”

**Preemption:** State the target regime and the worst case. Report pattern counts, residue counts, and coverage on social as well as clique-structured graphs.

Exact sentence:

> NSI exploits repeated maximal-clique profiles; when profiles are nearly unique or certification is sparse, its pattern and residue columns can be large, and no worst-case compression over explicit clique materialization is claimed.

### Attack 10: “Bytes per \(r\)-clique and nanosecond queries are under-specified.”

**Preemption:** Report absolute MB beside normalized size, define the denominator for a multi-\(r\) plane, count shared storage once, and disclose whether query timing includes vertex-to-class mapping, pattern formation, validation that the input is a clique, allocation, and I/O. Prefer median plus p95 over one average.

### Attack 11: “The index returns core numbers, not connected nuclei.”

**Preemption:** Preserve the corrected connectivity definition. State that NSI answers core-value queries; connected nuclei are obtained by an additional witness-connectivity/union-find pass. Do not let “indexing the decomposition” imply that component membership is returned unless it is actually stored or constructed.

Exact sentence:

> NSI indexes exact core values; connected \(k\)-nuclei are derived from the corresponding superlevel set using \(s\)-witness connectivity.

### Attack 12: “The hierarchy is claimed to be free.”

**Preemption:** Say “without rerunning decomposition,” not “free.” Give the union-find cost and measured overhead if hierarchy remains prominent. Treat the cohesion landscape as a derived capability, not a third pillar.

### Attack 13: “The case study varies \(s\) but calls it a support threshold.”

**Preemption:** Fix the case studies before using them as application evidence. \(s\) is the witness clique size, not the nucleus threshold \(k\). For every survival curve, specify the selected \(k\), or define a reproducible rule for choosing \(k\) at each cell.

### Attack 14: “Bit-exact experiments are presented as the proof of correctness.”

**Preemption:** State that theory proves correctness and differential tests validate the implementation.

Exact sentence:

> The theorems establish exactness; bit-exact agreement with independent fixed-\(r\) runs is an implementation check.

### Attack 15: “The largest results are censored speedups against failed runs.”

**Preemption:** Compute speedup only where both methods finish. For budget failures, state the budget and report the qualitative result separately: CND did not complete one cell, whereas NSI completed the requested family. Do not turn a timeout into a numeric speedup.

## Claims that must not appear

Delete or avoid the following claims in the title, abstract, contributions, captions, and conclusion:

1. **“The whole plane is built at the cost of one cell.”** False; per-\(r\) patterns, maps, boundaries, and residues remain.
2. **“The plane is cheaper than running the fixed-\(r\) engine once per \(r\).”** Not supported; the measured result is competitive, not cheaper.
3. **“One shared build computes every cell for free.”** “Shared build” covers the universal structure only. Certified cells still require checks, and residues may require peeling.
4. **“Every cell after the boundary is free.”** This is conditional on certification coverage and residue size.
5. **“Universal CPI proves the plane is fast.”** It proves exact counting coverage. Speed is empirical and data-dependent.
6. **“The algorithm is \(r\)-independent.”** The class-SCT is \(r\)-independent; pattern columns and peeling are not.
7. **“The diagonal theorem determines the next row.”** It gives an upper bound and a certificate only when that bound meets the floor.
8. **“One pattern represents \(\binom{w}{r}\) cliques” as the general rule.** The general multiplicity is a product of binomials across classes.
9. **“Deviations from the closed-form surface are always small.”** The experiments may show small residues on particular datasets; the theory does not.
10. **“NSI is always smaller/faster than per-cell methods.”** The advantage depends on class compression and certification coverage; current results include near-parity cases.
11. **“The method has a polynomial-size index.”** The class tree and pattern set can still be combinatorial in the worst case.
12. **“NSI never enumerates cliques” without qualification.** Say it does not materialize all individual \(r\)- or \(s\)-cliques during indexed construction; maximal-clique and pattern enumeration still occur.
13. **“Nanosecond query” as an end-to-end user latency.** It is an in-memory microbenchmark under a precisely stated input protocol.
14. **“The hierarchy is free.”** It avoids new decompositions but still requires connectivity processing.
15. **“CND cannot compute the spectrum.”** Say it did not complete under the stated time or memory budget on the tested machine.
16. **“First index over the plane” as an unqualified novelty claim.** Prefer “to our knowledge, the first exact index spanning both clique order \(r\) and witness size \(s\),” after a careful related-work check.
17. **“Bit-exactness proves correctness.”** It validates the implementation; the formal theorems prove correctness.

## Section-level prescription

1. **Title and Abstract:** replace the one-cell claim with the recommended two-pillar framing and the honest fixed-\(r\) comparison.
2. **Introduction:** make class symmetry the key observation; present batch computation and class compression as its two consequences. Define the two baselines before reporting results.
3. **Preliminaries:** define the bounded plane, relevant classes for \(S_{\min}\), patterns for each \(r\), and pattern multiplicity. Preserve nucleus connectivity.
4. **Straightforward Methods:** compare online per-cell CND, naive materialization, repeated optimized fixed-\(r\) columns, and NSI. The repeated fixed-\(r\) method must be visible in the comparison table.
5. **Theory:** keep the clique floor, shadow bound, and chain certificate; integrate Universal CPI and weighted-pattern equivalence. Describe the diagonal result as a conditional certificate, not an exact recurrence.
6. **Index:** formalize the plane layout, query algorithm, and size bound
   \[
   O\!\left(|V|+|\mathcal C|+|\mathcal L|+\sum_r\left(|\mathcal P_r|+\sum_s|D_{r,s}|\right)\right),
   \]
   plus explicitly stored region/profile metadata. State exactly what a point query includes.
7. **Construction:** present general product-of-binomials leaf counting first and \(\binom{w}{s-r}\) only as a special case. Replace Level Invariance as the replay justification with Shell-Order Replay and include complete pseudocode. Explain fixed-\(r\) first, then the outer plane loop.
8. **Hierarchy:** demote to a derived query/application. Replace “free” with “no additional decompositions,” and report construction overhead if retained.
9. **Experiments:** keep the strong fixed-row CND results, add the whole-plane experiment and the summed-fixed-\(r\) control, and separate measured facts from interpretation. State that SpecND is single-threaded and CND's available parallelism precisely.
10. **Case Studies:** repair the missing \(k\) semantics before claiming that an \(s\)-sweep changes the surviving set.
11. **Related Work:** position against single-cell decomposition, fixed-cell threshold indexes, and parameter-sweep indexes. Remove “cross-\(r\) is future work.”
12. **Conclusion:** end on the two supported consequences: family-at-a-time exact computation and class-compressed querying. Repeat that Universal CPI unifies the plane but does not remove per-\(r\) work.

## Final recommended framing in one paragraph

> Existing methods treat the \((r,s)\)-nucleus plane as unrelated cells: they recompute each cell and, if queries are needed later, would have to materialize every value. NSI shows that the plane has a shared class structure. One \(r\)-independent class-SCT exactly supports every bounded cell; certification and shell-order replay turn it into a family-at-a-time exact constructor, while the same classes collapse individual \(r\)-cliques into patterns and certified arithmetic tails. This yields orders-of-magnitude gains over the state-of-the-art per-cell baseline and a compact nanosecond-query index. The whole-plane build is deliberately claimed only to be competitive with repeated optimized fixed-\(r\) construction: the contribution is one exact reusable representation and the removal of repeated per-cell work, not the disappearance of inherently \(r\)-specific computation.
