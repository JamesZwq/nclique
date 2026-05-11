# Final Acceptance Checklist

1. Can the paper be summarized in one sentence?
Yes.
For \((1,s)\), vertex peeling preserves CPI structure, so exact decomposition can run over static paths with counter updates.

2. Is the central tension visible by the end of the introduction?
Yes.
The introduction contrasts general mutable CPI with the \(r=1\) invariant.

3. Is every key term defined before technical use?
Mostly yes.
\((1,s)\), CPI path, hold set, pivot set, target, support delta, H-operator, deletion log, and hierarchy are defined before algorithmic use.

4. Are contributions non-overlapping?
Yes.
They are written as structural theorem, counter arithmetic, fixed-point algorithms, construction consequence, hierarchy consequence, and evidence.

5. Does every contribution have evidence?
Yes, with TODOs for some experimental source files.

6. Does every algorithm have an invariant?
Yes.
Refresh, LocalH, PeelH, ParaBuild, and BuildHier each have state-machine prose before pseudocode.

7. Does every figure support a claim?
Mostly yes.
Figure `fig:overview` supports path encoding.
Figure `fig:peel` supports counter events.
Experiment figures support runtime/memory and stress behavior.
Case-study figures support granularity, cross-\(r\) tradeoffs, ego-network resolution, and hierarchy.

8. Are limitations and boundary cases explicit?
Yes.
The rewrite repeatedly states that static CPI is for \(r=1\), not general \((r,s)\).

9. Can a hostile reviewer call this "just a special case"?
Yes.
That remains a likely attack.

10. Where does the paper rebut that?
The introduction identifies \(r=1\) as the vertex-centric case.
StaticCPI shows it has a different state model.
CaseStudies show practical quality/cost value for the evaluated tasks.

11. Does the abstract tell a story rather than list components?
Yes.
It moves from high-order decomposition to mutable CPI, then invariant, then static-counter pipeline, then evidence.

12. Does the introduction make the method feel inevitable?
Mostly yes.
The dependency chain now derives LocalH, PeelH, ParaBuild, and BuildHier from the static-CPI invariant.

Remaining acceptance risks:
The experimental headline numbers need raw-source verification.
The theorem/proof of PeelH equivalence should be reviewed by a theory coauthor for edge cases around simultaneous equal-level peeling.
Some source files may need stricter sentence-line `%` formatting inside theorem/proof environments before final submission.
