# HIndex Plan

Section purpose:
Show that LocalH and PeelH are two realizations of the same static-CPI fixed point.

Reader knowledge before:
Static CPI paths have counters and support deltas.

Reader knowledge after:
The reader knows the H-operator view, the refresh primitive, LocalH, PeelH, their invariants, and their complexity.

Paragraph contracts:
1. State why the H-operator is needed after static counters exist.
2. Define \(\mathcal{H}\).
3. Define per-path witness count.
4. Give refresh state machine and algorithm.
5. Give LocalH state machine, algorithm, fixed-point proof, complexity.
6. Explain LocalH redundancy.
7. Give PeelH state machine, algorithm, equivalence proof, complexity.
8. Link to running example figure.

Key definitions introduced:
\(\mathcal{H}\), witness count \(W_v\), Refresh, LocalH, PeelH, \(U\).

Claims made:
LocalH converges to \(\kappa_s\); PeelH returns the same values with event-driven updates.

Evidence used:
Algorithms `alg:refresh`, `alg:localH`, `alg:peelH`; Lemmas `lem:per-path-witness`, `lem:hindex-fixed-point`, `lem:equivalence`; Theorems `thm:localh-complexity`, `thm:peelh-complexity`.

Sentence-level risks:
The \(O(1)\) incremental claim must be scoped to one changed path contribution and counted through \(U\).
