# StaticCPI Plan

Section purpose:
Prove the hidden invariant and make it the conceptual core of the paper.

Reader knowledge before:
CPI paths encode \(s\)-cliques through hold and pivot sets.

Reader knowledge after:
The reader knows which state remains static, which counters change, and why this is exact only for \(r=1\).

Paragraph contracts:
1. State the static/dynamic split.
2. Present vertex-removal theorem.
3. Prove hold removal, pivot removal, and unaffected cases.
4. Explain \(r\ge2\) boundary.
5. State support-delta lemma.
6. Link theorem to running example and figure.

Key definitions introduced:
Effective pivot counter \(p_{\TPath}\), liveness bit, support delta.

Claims made:
Path sets are never modified; hold peel kills path; pivot peel decrements counter; support losses are binomial differences.

Evidence used:
Theorem `thm:vertex-removal`, Lemma `lem:support-delta`, Figure `fig:overview`, Example `ex:peel`.

Sentence-level risks:
Do not say CPI is static for general \((r,s)\).
Do not confuse structural validity with liveness.
