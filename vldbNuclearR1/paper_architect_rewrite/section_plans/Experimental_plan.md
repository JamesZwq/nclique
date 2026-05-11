# Experimental Plan

Section purpose:
Evaluate claims, not figures.

Reader knowledge before:
The full algorithmic pipeline.

Reader knowledge after:
The reader knows which claims are experimentally supported and where evidence is incomplete or scoped.

Paragraph contracts:
1. State RQs.
2. Hardware, datasets, algorithms, methodology.
3. RQ1 correctness.
4. RQ2 end-to-end time/memory.
5. RQ3 phase breakdown and ablation.
6. RQ4 parallel construction.
7. RQ5 LocalH behavior.
8. RQ6 dense stress behavior.
9. Per-graph summary with TODO where source is unclear.

Key definitions introduced:
REF, ST, V2, V3, phase names.

Claims made:
Exactness against REF; static state reduces runtime/memory; construction becomes bottleneck; ParaBuild scales; LocalH matches semantics but is slower; dense stress remains favorable where both finish.

Evidence used:
Tables `tab:datasets`, `tab:mem`, `tab:bd-time`, `tab:bd-mem`, `tab:par`, `tab:per-graph`; Figures `fig:exp-time`, `fig:exp-mem`, `fig:stress-time`, `fig:stress-mem`; CSV `figures/benchmark_all_results.csv`.

Sentence-level risks:
Do not invent numbers.
Use TODO comments where existing figures/tables are not locally backed by `benchmark_all_results.csv`.
