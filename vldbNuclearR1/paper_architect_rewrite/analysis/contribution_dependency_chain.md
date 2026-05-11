# Contribution Dependency Chain

## Chain Form

1. Structural theorem:
For \(r=1\), vertex peeling preserves the structural validity of each CPI path.
The only dynamic state needed per path is a liveness bit, the fixed target, and an effective pivot counter.

2. Counter arithmetic:
Because path state is static except for counters, every support loss has a closed-form binomial delta.
This turns residual-index maintenance into arithmetic over path counters.

3. Core computation:
The counter view supports the same \(h\)-index fixed point as classical core decomposition.
LocalH states the synchronous fixed-point computation.
PeelH realizes the same computation in deletion order without repeated global refreshes.

4. Construction consequence:
Once peeling is cheap, CPI construction becomes a visible or dominant cost.
Static CPI consumers need only a path list and incidence CSR, so construction can stream emitted paths into thread-local arenas and merge deterministically.

5. Hierarchy consequence:
PeelH already emits vertices in nondecreasing \(\kappa_s\).
The hierarchy is therefore a reverse scan of the deletion log rather than a separate sorted pass.

6. Empirical consequence:
The combined pipeline is exact against the mutable baseline and reduces measured runtime and memory on the reported benchmarks.

7. Use-case consequence:
Case studies show that the vertex-centric \(r=1\) setting is not merely a technical corner case; it can match higher-\(r\) quality metrics at much lower cost on the evaluated community tasks.

## Non-overlap Rule

Do not list StaticCPI, LocalH, PeelH, ParaBuild, and BuildHier as independent modules.
Each contribution should depend on the previous one.
The introduction should make the method feel inevitable from the invariant.
