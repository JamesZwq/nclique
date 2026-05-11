# Sentence Utility Report

This audit samples the rewritten sections and checks whether each paragraph has a job.
Full sentence-by-sentence auditing would be longer than the paper itself, so this report audits representative paragraphs and records the delete-test rule used for the whole rewrite.

## Introduction

Paragraph: opening motivation.
Sentence role: define why edge cores are insufficient and introduce \((1,s)\).
Reader uncertainty reduced: what problem is being solved and why \(s\)-cliques matter.
What breaks if deleted: the reader sees CPI before knowing why high-order vertex cohesion matters.
Decision: Keep.

Paragraph: CPI baseline.
Sentence role: explain the existing representation and its cost.
Reader uncertainty reduced: why the paper is not just another clique enumeration paper.
What breaks if deleted: the static invariant has no target to improve.
Decision: Keep.

Paragraph: \(r=1\) invariant.
Sentence role: state the hidden invariant in plain terms.
Reader uncertainty reduced: why a special case can change the state model.
What breaks if deleted: the paper becomes a list of components again.
Decision: Keep.

Paragraph: contributions.
Sentence role: convert modules into dependency chain.
Reader uncertainty reduced: how each contribution follows from the previous one.
What breaks if deleted: reviewers may see overlapping components.
Decision: Keep, with TODO retained for unverified headline.

## Preliminaries

Paragraph: graph notation.
Sentence role: define only symbols needed later.
Reader uncertainty reduced: base graph objects.
What breaks if deleted: later definitions become underspecified.
Decision: Keep.

Paragraph: mutable baseline.
Sentence role: identify the exact state removed by the new method.
Reader uncertainty reduced: what "mutable CPI" means concretely.
What breaks if deleted: the contribution loses contrast.
Decision: Keep.

## StaticCPI

Paragraph: static/dynamic split.
Sentence role: separate stored path structure from live counters.
Reader uncertainty reduced: what remains static and what changes.
What breaks if deleted: theorem reads like a trivial vertex deletion statement.
Decision: Keep.

Proof of Theorem 1.
Sentence role: handle hold, pivot, and unrelated cases.
Reader uncertainty reduced: why counter decrement is support-equivalent to set deletion.
What breaks if deleted: core correctness has no structural basis.
Decision: Keep.

Boundary paragraph.
Sentence role: prevent overclaiming to \(r\ge2\).
Reader uncertainty reduced: why the paper does not subsume general CPI.
What breaks if deleted: hostile prior-work reviewer can attack the scope.
Decision: Keep.

## HIndex

Opening paragraph.
Sentence role: connect static counters to core computation.
Reader uncertainty reduced: why \(h\)-index appears after StaticCPI.
What breaks if deleted: LocalH and PeelH feel unrelated.
Decision: Keep.

Refresh state-machine paragraph.
Sentence role: explain input state, mutable state, invariant, transition, output, cost driver before pseudocode.
Reader uncertainty reduced: how to read Algorithm 1.
What breaks if deleted: pseudocode has no state contract.
Decision: Keep.

PeelH pseudocode explanation.
Sentence role: distinguish semantic Refresh from incremental delta implementation.
Reader uncertainty reduced: why the algorithm remains equivalent while the implementation is faster.
What breaks if deleted: the \(O(1)\) delta claim looks unsupported.
Decision: Keep.

## ParaBuild

Opening paragraph.
Sentence role: show bottleneck shift.
Reader uncertainty reduced: why construction belongs in the main story.
What breaks if deleted: ParaBuild looks like an isolated engineering trick.
Decision: Keep.

Scope paragraph.
Sentence role: limit claims about parallelism.
Reader uncertainty reduced: what is and is not being claimed.
What breaks if deleted: systems reviewer may read ideal-scaling overclaim.
Decision: Keep.

## BuildHier

Opening paragraphs.
Sentence role: explain super-level filtration before elder-rule mechanics.
Reader uncertainty reduced: what hierarchy means.
What breaks if deleted: algorithm appears without motivation.
Decision: Keep.

Complexity paragraph.
Sentence role: identify edge scan and union-find cost.
Reader uncertainty reduced: why CPI is not consulted.
What breaks if deleted: hierarchy cost claim lacks support.
Decision: Keep.

## Experiments

RQ paragraph.
Sentence role: map evaluation to claims.
Reader uncertainty reduced: what each experiment proves.
What breaks if deleted: section returns to figure-by-figure reporting.
Decision: Keep.

CSV-backed subset paragraph.
Sentence role: separate locally verified evidence from broader draft claims.
Reader uncertainty reduced: which numbers are directly grounded in `benchmark_all_results.csv`.
What breaks if deleted: unsupported headline numbers may be overclaimed.
Decision: Keep.

Takeaway paragraph.
Sentence role: separate measured facts by evidence strength.
Reader uncertainty reduced: what is verified, partially verified, and pending.
What breaks if deleted: paper hides evidence gaps.
Decision: Keep.

## CaseStudies

Opening paragraph.
Sentence role: state practical relevance question and scope.
Reader uncertainty reduced: why case studies are not decoration.
What breaks if deleted: section becomes a collection of unrelated plots.
Decision: Keep.

Cross-\(r\) paragraphs.
Sentence role: show quality/cost tradeoff.
Reader uncertainty reduced: why optimizing \(r=1\) matters.
What breaks if deleted: hostile reviewer can call \(r=1\) unimportant.
Decision: Keep.

Takeaway paragraph.
Sentence role: state scoped conclusion and limitation.
Reader uncertainty reduced: avoids universal application claim.
What breaks if deleted: section overclaims practical generality.
Decision: Keep.

## RelatedWork

Category paragraphs.
Sentence role: place prior work by problem solved and limitation.
Reader uncertainty reduced: why this paper needs a different idea.
What breaks if deleted: related work becomes a citation list.
Decision: Keep.

## Conclusion

First paragraph.
Sentence role: return to the state-model change.
Reader uncertainty reduced: what the reader should remember.
What breaks if deleted: conclusion becomes a module list.
Decision: Keep.

Future-work sentence.
Sentence role: state boundary and next question.
Reader uncertainty reduced: what remains unresolved.
What breaks if deleted: no research outlook.
Decision: Keep.
