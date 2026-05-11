# ParaBuild Plan

Section purpose:
Frame parallel construction as the consequence of bottleneck shifting.

Reader knowledge before:
PeelH makes support maintenance cheap and consumes a static CPI layout.

Reader knowledge after:
The reader knows why seed-owned construction is parallel and why the output can be streamed.

Paragraph contracts:
1. Use phase breakdown to motivate construction.
2. Explain seed ownership of \(s\)-cliques.
3. Explain thread-local arenas and deterministic merge.
4. Present algorithm.
5. State correctness and limitations.

Key definitions introduced:
Seed ownership, thread-local arena, deterministic merge.

Claims made:
Construction has no hot-path synchronization; merge is deterministic; ParaBuild builds the same CPI semantics.

Evidence used:
Algorithm `alg:parabuild`, Table `tab:par`, prior BK/Tomita citations.

Sentence-level risks:
Do not overclaim ideal parallel scaling or claim PeelH parallelization is unnecessary in all settings.
