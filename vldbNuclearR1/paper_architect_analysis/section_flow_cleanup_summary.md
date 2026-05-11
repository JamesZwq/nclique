# Section Flow Cleanup Summary

## Subsections Merged Or Removed

- Section 4 was reduced from five subsections to three narrative subsections.
- `The H-Operator` and `The Refresh Primitive` were merged into `4.1 A fixed-point view over static CPI`.
- `From Synchronous Refresh to Event-Driven Peel` was removed as a standalone subsection and folded into `4.3 PeelH: the event-driven realization`.
- The remaining Section 4 subsections are:
  - `4.1 A fixed-point view over static CPI`
  - `4.2 LocalH: the direct synchronous solver`
  - `4.3 PeelH: the event-driven realization`
- Existing labels were preserved by attaching the old labels to the merged subsections: `sec:hop`, `sec:refresh`, `sec:localh`, `sec:redundancy`, and `sec:peelh`.

## Section Openings Rewritten

- `Preliminaries`: added a bridge from the informal thesis to the definitions needed for support, nuclei, CPI paths, and the mutable baseline.
- `StaticCPI`: added a bridge from the mutable residual-index baseline to the specific question of whether vertex peeling needs mutable CPI state.
- `HIndex`: rewrote the opening to connect static path counters to the remaining problem of computing core values.
- `ParaBuild`: rewrote the opening around bottleneck shift after cheap peeling.
- `BuildHier`: rewrote the opening around the deletion log as an additional output of PeelH and weakened the hierarchy framing to the edge-connected core-value filtration.
- `Experimental`: rewrote the opening so the RQs follow the dependency chain instead of reading as a checklist.
- `CaseStudies`: added a bridge from algorithmic efficiency to the practical question of whether the \(r=1\) target is useful.
- `RelatedWork`: added a framing paragraph that positions prior work around peeling, compact clique representations, local updates, construction, and hierarchy outputs.
- `Conclusion`: added a return to the paper's opening tension between compact clique encoding and mutable residual state.

## Transitions Added

- Added a transition before Theorem 1 in `StaticCPI` to state that the theorem formalizes path death versus pivot-counter reduction.
- Added a fixed-point bridge before the H-operator definition in Section 4.
- Added a bridge from Lemma `per-path-witness` to `refreshOp`, explaining why the refresh primitive is shared by LocalH and PeelH.
- Added a transition before `LocalH` that explains it as a literal synchronous fixed-point solver.
- Added a transition before the LocalH correctness lemma explaining why the proof makes LocalH a semantic reference.
- Added a transition before `PeelH` that explains the redundancy in LocalH and the local state needed to avoid global rounds.
- Rewrote the ParaBuild pre-algorithm prose around thread-local ownership and deterministic merge.
- Rewrote the BuildHier pre-algorithm prose around the reverse scan invariant.
- Reworked the experiment opening so RQ1--RQ6 are tied to the method's dependency chain.

## Remaining Abrupt Places

- The theorem and proof blocks are still dense, but their statements and labels were intentionally preserved.
- The experiment setup still uses compact paragraph headings for hardware, datasets, algorithms, and methodology; these are mechanical but useful for reviewer lookup.
- The case-study section still uses numbered paragraph headings because each case study asks a distinct empirical question.
- Some TODO comments remain where earlier passes marked unverified experimental claims; they do not enter the PDF.
