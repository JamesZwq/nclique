# Section Rewrite Plan

## Introduction.tex

Current role:
Compact overview and contribution list.

Problem:
It introduces many module names early and reads as a component list.

Desired role:
Make the invariant feel inevitable from the tension between general CPI and vertex peeling.

Required rewrite:
Start with higher-order vertex cohesion.
Introduce general CPI as powerful but mutable.
Show why \(r=1\) changes the state model.
Derive LocalH, PeelH, ParaBuild, and BuildHier as consequences.

Dependencies:
Prelim definitions must not be required too early.

Key claims:
Static CPI is exact for \(r=1\); counter-only updates suffice; pipeline is faster and useful.

Evidence:
Theorem `thm:vertex-removal`, Lemma `lem:support-delta`, algorithms, experiments, case studies.

## Preliminaries.tex

Current role:
Definitions and baseline.

Problem:
Definition of nucleus uses \(k\)-admissible cliques before residual-subgraph support is fully intuitive.

Desired role:
Define the target object and CPI with minimum notation.

Required rewrite:
Define support in an induced subgraph where needed.
Define core number, CPI path, support formula, mutable baseline.

Dependencies:
None.

Key claims:
CPI path support can be computed by binomial counts.

Evidence:
Definition `def:path`, Lemma `lem:supp`, prior CPI citation.

## StaticCPI.tex

Current role:
Core theorem and delta.

Problem:
The "why specific to \(r=1\)" point is too short.

Desired role:
Conceptual core.

Required rewrite:
Separate static structure, dynamic counters, path death, and exact support equivalence.

Dependencies:
CPI path and support formula.

Key claims:
Vertex peeling preserves CPI structure; counters suffice.

Evidence:
Theorem `thm:vertex-removal`, Lemma `lem:support-delta`, figures.

## HIndex.tex

Current role:
H-operator, LocalH, PeelH, equivalence.

Problem:
Refresh primitive and \(O(1)\) incremental update are mixed in a way that may invite theory attacks.

Desired role:
Show two realizations of one fixed-point computation.

Required rewrite:
Add state-machine explanations before Refresh, LocalH, PeelH.
Make \(U\) explicit.

Dependencies:
Static CPI theorem and support delta.

Key claims:
LocalH computes fixed point; PeelH removes synchronous redundancy and computes same values.

Evidence:
Lemmas `lem:hindex-fixed-point`, `lem:equivalence`; algorithms.

## ParaBuild.tex

Current role:
Parallel construction.

Problem:
It sounds like an isolated engineering add-on.

Desired role:
Consequence of bottleneck shift.

Required rewrite:
Start from phase breakdown.
Explain seed ownership, static output, thread-local arenas, deterministic merge.

Dependencies:
CPI construction and static CPI consumer layout.

Key claims:
Construction can stream static paths without hot-path synchronization.

Evidence:
Algorithm `alg:parabuild`, Table `tab:par`, phase breakdown.

## BuildHier.tex

Current role:
Reverse scan algorithm.

Problem:
Assumes elder-rule terminology.

Desired role:
Explain hierarchy as a byproduct of sorted deletion.

Required rewrite:
Give intuition first, formal setup second.

Dependencies:
Core values and deletion log.

Key claims:
Reverse scan builds join tree in \(O((n+m)\alpha(n))\).

Evidence:
Algorithm `alg:buildhier`, complexity proof.

## Experimental.tex

Current role:
Large figure-by-figure report.

Problem:
Claims and figures are not organized as RQs.

Desired role:
Claim-driven evaluation.

Required rewrite:
Use RQ1 correctness, RQ2 end-to-end, RQ3 phase/ablation, RQ4 parallel construction, RQ5 LocalH, RQ6 stress.
Mark unclear source claims with TODO comments.

Dependencies:
Claim-evidence ledger.

Key claims:
Exactness, speed, memory, phase attribution, parallel scaling, boundary behavior.

Evidence:
Existing tables, figures, `benchmark_all_results.csv`.

## CaseStudies.tex

Current role:
Many case studies with practical guidelines.

Problem:
Long and sometimes decorative.

Desired role:
Defend why \(r=1\) matters.

Required rewrite:
Tie each case study to vertex-centric usefulness and cost against higher \(r\).

Dependencies:
Definitions of \((r,s)\) for comparison.

Key claims:
\((1,s)\) can match higher-\(r\) quality on evaluated tasks while being much cheaper.

Evidence:
Tables `tab:cs3`, `tab:cs4`, `tab:cs7`, figures.

## RelatedWork.tex

Current role:
Short list.

Problem:
Not enough positioning.

Desired role:
Category-based contrast.

Required rewrite:
Use enumeration-based nucleus decomposition, parallel/local peeling, clique counting/compact representations, hierarchy methods.

Dependencies:
Main method already known.

Key claims:
Prior work either handles general \(r\) with mutable state, counts cliques, or extracts hierarchies separately.

Evidence:
Citations already in draft.

## Conclusion.tex

Current role:
Dense summary.

Problem:
Too many modules in one paragraph.

Desired role:
Short return to thesis.

Required rewrite:
State state-model change, consequences, empirical implication, limitation.

Dependencies:
Whole paper.

Key claims:
Static symbolic CPI is sufficient for \(r=1\).

Evidence:
All prior sections.
