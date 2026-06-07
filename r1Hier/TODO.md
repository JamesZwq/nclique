# r1Hier — TODO / Next Steps

**Last updated:** 2026-06-08
See `HANDOFF.md` for full context. Ordered by value/priority.

## P0 — core gaps in the index story

- [ ] **Materialize the hierarchy T_s, not just κ_s.** Query currently returns
      core values. The user's vision is "return the full s-clique core hierarchy
      tree T_s". Needs, per level s: connected components of the s-clique-core
      induced subgraph, nested across core values. Approach: union-find over the
      level-s structure; decide storage (per-s component forest vs a single
      cross-s persistent tree). This is the missing headline capability.
- [ ] **Parallelize the peel.** `ST_V3_Peel` is single-thread and dominates
      `spectrum_ms` on dense graphs (e.g. AstroPh peel ~233ms of 300ms). The
      recompute is already OpenMP. A parallel peel would cut index-build time and
      shrink the `shared` baseline further (re-check whether dense still wins by
      the same margin — it should, but measure).

## P1 — stronger experiments / fairness

- [ ] **Parallelize `shared` baseline too**, then re-run the 3-tier query table.
      The user explicitly asked whether the fair baseline can be pushed lower.
- [ ] **bytes-per-vertex / Σ(clique incidences) vs ω & density plot.** Show when
      the CPI tree memory blows up (cit-HepPh 1339 B/vtx) vs stays small.
- [ ] **Larger graphs** (soc-pokec, com-lj, com-orkut, com-friendster). Check
      universal-CPI memory (Σ) scaling and whether build-once still wins. Watch
      the 1e18 clamp; report fuzzy fractions. Friendster needs the server.
- [ ] **Middle-tier as a first-class query path** is now in `bench_spectrum_index`
      (`shared_ms`); consider also reporting its *memory* (= CPI tree) as the
      "no-materialization" point on the size/query Pareto curve.

## P2 — compression (weak, only if pursuing the compression angle)

- [ ] **Trajectory equivalence classes** — the one lever that can beat the
      anchor's avgTraj/2 ceiling. Vertices in the same dense region often share an
      *identical* trajectory (all K_m vertices do). Group vertices by exact
      trajectory; store each distinct trajectory once + a per-vertex reference.
      Measure #distinct trajectories vs n_active and resulting compression. This
      is the deciding experiment for whether compression is worth a paper.
- [ ] Bit-level size accounting (s* in ⌈log₂ω⌉ bits, varint deltas) instead of
      counting "values" — only changes absolute numbers, not the ranking.

## P3 — robustness / polish

- [ ] **128-bit / exact big-int support value** to remove the 2⁵³ clamp on dense
      high-ω graphs (or document it as an accepted limit). Currently `fuzzy`.
- [ ] Serialize the index to disk + measure load time (so "index size" is a real
      on-disk artifact, not just resident bytes).
- [ ] Fold the index targets into `CMakeLists` test/verify harness.

## Done (for reference — see HANDOFF.md §4)

- [x] ω + saturation/compression measurement on 9 graphs (`bench_spectrum_sat`).
- [x] KK Shifted Containment Law + support nesting verified on real graphs (law_viol=0).
- [x] Build-once universal CPI (`SDCT_Augmented_NoTree_Universal`), 1.6×–233× faster build.
- [x] Materialized dense + anchor+delta index; per-vertex + all-s correctness (0 mismatches).
- [x] Query latency vs 2 fair baselines (dense 113×–88000× vs shared).
- [x] Memory: index vs CPI tree (dense 0.4–0.6×, down to 0.07×, of the tree).

## Watch-outs (so we don't lose work again)

- `docs/` and `bench_plots/` are **gitignored** — keep canonical copies of
  findings/CSVs in `r1Hier/artifacts/` and commit them.
- The repo `main` is shared with concurrent SIGMOD-paper work; `git log -N` may
  show *their* commits on top. Our spectrum commits: `1dd29bd`, `52e6e07`,
  `b70358c`, `3933a79` (all ancestors of HEAD, verified safe).
