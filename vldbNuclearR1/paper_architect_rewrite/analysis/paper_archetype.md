# Paper Archetype

Dominant archetype: hidden-invariant plus bottleneck-shifting paper.

The hidden invariant is that vertex peeling for the \((1,s)\)-nucleus removes vertices but does not delete edges.
This preserves the structural validity of every Clique Path Index path.
The paper should make this invariant visible early, because it explains why a general mutable CPI can become a static symbolic index in the special \(r=1\) case.

The bottleneck shift follows from the invariant.
Once peeling no longer mutates CPI paths and only updates counters, residual support maintenance becomes cheap.
The dominant cost then moves to CPI construction, which motivates ParaBuild as part of the same dependency chain rather than as a separate systems add-on.

Secondary archetype: exact-special-case paper.
The work does not claim to solve general \((r,s)\)-nucleus decomposition.
It instead argues that the important vertex-centric case has a stronger invariant than the general problem, and that exploiting this invariant produces an exact algorithmic and systems pipeline.

Reviewer risk: a hostile reader may call the work "just \(r=1\)."
The rebuttal must be that \(r=1\) is the vertex-ranking case, is empirically useful in the case studies, and admits a qualitatively different state model with static paths and counter-only updates.
