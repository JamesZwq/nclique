# Sentence-Level Rules

1. Every sentence must answer one reader question, define one object, state one claim, give one piece of evidence, or create one necessary transition.

2. A sentence may not remain merely because it is true.
If deleting it does not change the reader's ability to follow the argument, delete it.

3. Introduce at most one new concept per sentence.
If a sentence contains two new terms, split it or delay one term.

4. Prefer concrete verbs.
Use count, store, remove, update, scan, encode, preserve, avoid, merge, emit, finalise.
Avoid leverage, utilize, facilitate, conduct, perform, and propose unless the sentence would otherwise become awkward.

5. Separate mechanism from evidence.
A mechanism sentence says why the algorithm works.
An evidence sentence says where a theorem, algorithm, table, figure, or experiment supports the claim.

6. Do not front-load internal names.
Introduce StaticCPI, LocalH, PeelH, ParaBuild, and BuildHier only after the reader knows the problem each name solves.

7. Keep \(r=1\) qualifications visible.
Never let a sentence imply that static CPI solves general \((r,s)\) decomposition.

8. Use sentence-line LaTeX formatting for prose.
Put one sentence on one source line.
Put a standalone `%` line between consecutive sentences in the same paragraph.
Use a blank line only for a real paragraph break.

9. Do not force sentence-line formatting inside equations, algorithms, tables, bibliography entries, macros, or figure code.

10. Use TODO comments for unsupported or unclear evidence.
Format:
`% TODO: verify this claim against [file/figure/table].`
