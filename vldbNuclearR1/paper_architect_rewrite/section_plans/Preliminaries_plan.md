# Preliminaries Plan

Section purpose:
Define the target problem and the CPI representation with enough precision for the static theorem.

Reader knowledge before:
Basic graph terminology.

Reader knowledge after:
The reader knows \((1,s)\)-core values, CPI paths, support formula, and the mutable baseline's state.

Paragraph contracts:
1. Define graph notation.
2. Define clique and \(s\)-clique support.
3. Define \((1,s)\)-nucleus and core number.
4. Define CPI path, hold set, pivot set, target, and index size.
5. State the support formula.
6. Describe mutable peeling baseline.
7. State the exact problem.

Key definitions introduced:
\(k\)-clique, \(\sdeg\), \(k\)-\((1,s)\)-nucleus, \(\kappa_s\), CPI path, \(\eta_P\), \(\Sig\).

Claims made:
CPI encodes every \(s\)-clique exactly once; mutable baseline maintains residual state.

Evidence used:
Definitions and Lemma `lem:supp`; citations `NuclearCD`, `BronKerbosch`, `eppstein`, `tomita`.

Sentence-level risks:
Do not overload \(k\) and \(s\).
Define each symbol before use.
