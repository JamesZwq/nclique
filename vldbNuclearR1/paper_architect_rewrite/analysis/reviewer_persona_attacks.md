# Reviewer Persona Attacks

## Algorithm/Theory Reviewer

Objection:
The static CPI theorem is obvious because vertices are removed, not edges.

Required revision:
Separate structural preservation from support accounting.
The theorem is not just that the clique remains a clique; it is that the encoded family after each event is support-equivalent to static sets plus counters, which preserves exact peeling.

Objection:
PeelH equivalence to LocalH is underspecified.

Required revision:
Give a state-machine explanation before pseudocode.
State the invariant over unfinalized vertices and the path-counter state.
Avoid suggesting a full Refresh can become \(O(1)\) unless the update is explicitly tied to one changed path and the affected refresh count \(U\).

Objection:
\(U\) could be large.

Required revision:
Define \(U\) as measured work and avoid hiding it.
Explain when the asymptotic win is from removing \(\log n\) and mutable hash/tree updates rather than from reducing every possible incidence.

## Systems/Performance Reviewer

Objection:
The headline speedup is a pipeline effect, not just PeelH.

Required revision:
Make that explicit.
Use phase breakdown and ablation to attribute speedup to static state, tree-free construction, event-driven peeling, and parallel construction.

Objection:
The ten-graph headline is not locally reproducible from `figures/benchmark_all_results.csv`.

Required revision:
Add TODO verification comments near claims whose source is outside the allowed local evidence.
Report the local 505-cell subset separately where possible.

Objection:
Parallel construction may be load-imbalanced.

Required revision:
Say contiguous seed chunks are the presented design, and use the measured scaling table instead of claiming ideal scaling.

## Graph Mining/Application Reviewer

Objection:
Why care about \(r=1\) if nucleus decomposition was designed for general \(r\)?

Required revision:
Position \(r=1\) as the vertex-ranking case.
Use case studies to show that tuning \(s\) can match higher-\(r\) quality for the reported tasks.
Avoid saying higher \(r\) is never useful.

Objection:
Case studies are mostly DBLP and YouTube.

Required revision:
Frame them as evidence of practical relevance, not universal validation.
Use hierarchy examples across DBLP, YouTube, and web-Stanford as qualitative support.

## Hostile Prior-Work Expert

Objection:
This is only a specialization of NuclearCD.

Required revision:
Acknowledge that the construction and support formula come from prior CPI work.
State the contribution as a new state model for the \(r=1\) specialization: static paths plus counters.

Objection:
The paper overclaims that mutable CPI is unnecessary.

Required revision:
Always qualify with \(r=1\).
Explicitly say the general \(r\ge2\) case can break path structure and is not subsumed.

Objection:
LocalH is a known \(h\)-index core algorithm.

Required revision:
Present LocalH as the known fixed-point view lifted to the implicit \(s\)-clique hypergraph encoded by CPI.
The new part is the static-CPI counter evaluation and the event-driven realization.
