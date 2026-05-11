# Rewrite Diff Summary

## Story Changes

The original draft was organized around named components.
The rewrite is organized around a dependency chain.
The central story is now:
general CPI is powerful but mutable, \(r=1\) preserves path structure, static paths reduce peeling to counters, counter-only peeling exposes construction, and the deletion log gives hierarchy.

## Contribution Changes

The contribution list now separates:
1. static-\cpi{} state theorem,
2. counter arithmetic,
3. fixed-point algorithms,
4. construction after bottleneck shift,
5. hierarchy from deletion log,
6. empirical and application evidence.

This avoids treating StaticCPI, LocalH, PeelH, ParaBuild, and BuildHier as unrelated modules.

## Section-Level Changes

Introduction:
Rewritten from scratch as a technical story.

Preliminaries:
Definitions are reordered and simplified.
Mutable baseline is described as the specific state model being removed.

StaticCPI:
Expanded into the conceptual core.
The static/dynamic split and \(r\ge2\) boundary are explicit.

HIndex:
LocalH and PeelH are presented as two realizations of one fixed-point computation.
Each algorithm has a state-machine explanation before pseudocode.

ParaBuild:
Reframed as a consequence of bottleneck shifting.

BuildHier:
Reframed as a reverse scan of a sorted deletion log, with elder-rule intuition first.

Experimental:
Rewritten around RQs.
The locally CSV-backed subset is separated from broader draft headline numbers.
TODO comments mark claims whose raw source is not in `figures/benchmark_all_results.csv`.

CaseStudies:
Reframed as the defense of practical relevance for \(r=1\).
The universal claim was weakened to a scoped empirical claim.

RelatedWork:
Rewritten by categories rather than citation list.

Conclusion:
Shortened and returned to the state-model thesis.

## Terminology Changes

"Static CPI" now means stored hold/pivot path sets remain unchanged while liveness and pivot counters change.
"Counter-only" is used for residual state updates, not for construction.
"Ours" is scoped to V3 in experiments.
"Higher-\(r\)" claims are scoped to evaluated tasks.

## Claims Strengthened

The hidden invariant is stated earlier and more explicitly.
The dependency between StaticCPI, PeelH, ParaBuild, and BuildHier is stronger.
The \(r=1\) boundary is clearer.

## Claims Weakened Or Marked TODO

The ten-graph \(762\)-configuration headline is marked TODO because the local `figures/benchmark_all_results.csv` supports a \(505\)-cell three-graph subset.
LocalH timing claims are marked TODO for raw-log verification.
Stress-test exact memory values are marked TODO for raw CSV verification.
Parallel-construction aggregate speedup is marked TODO for raw-source verification.

## Risks Remaining

The rewritten proof of PeelH equivalence should be checked carefully for equal-level bucket behavior.
The source formatting mostly follows sentence-line style, but theorem/proof environments may need a stricter percent-line pass.
The final compile succeeded from the repository root.
It still leaves one non-fatal overfull-vbox warning and inherited bibliography metadata warnings.
