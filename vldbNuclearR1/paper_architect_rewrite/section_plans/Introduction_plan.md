# Introduction Plan

Section purpose:
Make the paper's dependency chain visible.

Reader knowledge before:
The reader knows core decomposition and may know nucleus decomposition, but not CPI internals.

Reader knowledge after:
The reader understands why general CPI is powerful but costly, why \(r=1\) exposes a static-path invariant, and why the rest of the pipeline follows.

Paragraph contracts:
1. Motivate higher-order vertex cohesion.
2. Introduce CPI as the current exact compact representation and identify mutable residual state as the bottleneck.
3. State the \(r=1\) invariant.
4. Derive counter-only peeling and the two \(h\)-index realizations.
5. Explain bottleneck shift to construction and ParaBuild.
6. Explain hierarchy as deletion-log byproduct.
7. List non-overlapping contributions.

Key definitions introduced:
\((1,s)\)-nucleus informally, CPI informally, hold/pivot behavior informally.

Claims made:
Static CPI is exact for \(r=1\); PeelH removes LocalH redundancy; ParaBuild addresses shifted bottleneck; BuildHier follows from sorted deletion; experiments support speed and usefulness.

Evidence used:
Theorem `thm:vertex-removal`, Lemma `lem:support-delta`, Algorithms `alg:localH`, `alg:peelH`, `alg:parabuild`, `alg:buildhier`, Experiments, CaseStudies.

Sentence-level risks:
Avoid too many module names before the reader knows why each exists.
Mark ten-graph headline as TODO where local evidence is unclear.
