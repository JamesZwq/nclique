# Nucleus Spectrum: ω and Saturation Measurements on Real Graphs

Measured by `bench_spectrum_sat` (source: `bench_spectrum_sat.cpp`). This answers
"how compressible is the r=1 vertex s-clique-core *spectrum* {κ_s(v) : s=2..ω}?"
using the existing SPIN★ build/peel as the engine.

## 1. Problem

For each vertex `v` the spectrum is its trajectory `κ_2(v), κ_3(v), …, κ_{s*}(v)`
where `s*(v)=max{s : κ_s(v)>0}`. We encode each trajectory as

* anchor `(s*, κ_{s*})`  +  sparse deltas `δ_s = κ_s − ∂_s(κ_{s+1}) ≥ 0`

and reconstruct downward via the **Kruskal–Katona lower shadow** `∂_s`. A vertex is
**saturated** iff all its `δ_s = 0` (the whole trajectory follows from the anchor).
Compression ratio = (naive: store every κ_s) / (anchor+delta: `2·n_active + #δ≠0`).

## 2. Method / reproducibility

* Engine: `NCliqueVertexCoreDecomposition_ST_V3_Build/_Peel` (= SPIN★, the trusted
  path used by `degeneracy_cliques` under `PIVOTER_RUN_ST_V3`).
* **The CPI build is NOT reusable across s.** `SDCT_Augmented_NoTree(g, max_k=s,
  min_k=1)` emits a leaf only when the maximal-clique size `cSize ≥ s`
  (`bkRecurse_NoTree`, line 51), with a per-s `(keepV,dropV)` decomposition. Hence
  we rebuild the CPI for every `s`; ω is detected as the first `s` with
  `numLeaves==0`. Streaming keeps memory at O(n) (two κ arrays + per-vertex δ count).
* `∂_s` ported and self-tested against hand-verified karate-club values (`--selftest`).
* Correctness gate: κ_s is **bit-identical** to the `degeneracy_cliques` ST_V3 dump
  (checked at s=3 and s=6 on `mini_diff_8v`), and the K₅/K₆ clique signatures
  `[4,6,4,1]` / `[5,10,10,5,1]` reproduce exactly.

Build & run:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j 12 \
      --target bench_spectrum_sat
./build/bin/bench_spectrum_sat <graph.edges> [--smax N] [--sort degen]
./run_spectrum_sat.sh <out_dir> graphs/*.edges     # batch -> CSV
```
Machine: Darwin arm64 (note: `long double == double` here, so κ ≥ 2⁵³ is "fuzzy").
Single-thread. Data: `bench_plots/spectrum_sat.csv`.

## 3. Results (single-thread, degeneracy order)

| graph | n | m | ω | avg traj | sat% | comp | law_viol | fuzzy |
|---|---|---|---|---|---|---|---|---|
| bio-celegans | 453 | 2.0K | 9 | 3.76 | 44% | 1.20× | 0 | 0 |
| ca-GrQc | 5.2K | 14K | 44 | 3.82 | **89%** | 1.80× | 0 | 0 |
| ca-HepTh | 9.9K | 26K | 32 | 2.75 | 77% | 1.19× | 0 | 0 |
| ca-CondMat | 23K | 93K | 26 | 4.49 | 79% | 1.94× | 0 | 0 |
| ca-AstroPh | 19K | 198K | 57 | 10.95 | 66% | **3.67×** | 0 | 0 |
| cit-HepPh | 35K | 421K | 19 | 4.91 | 13% | **0.90×** | 0 | 0 |
| amazon0302 | 262K | 900K | 7 | 3.15 | 31% | 1.13× | 0 | 0 |
| com-dblp | 317K | 1.05M | 114 | 3.93 | 86%¹ | 1.76×¹ | 0 | 19334 |
| ca-HepPh | 89K | 118K | 239 | 14.41 | 76%¹ | 2.85×¹ | 0 | 76069 |

¹ exact-regime (κ<2⁵³) figures; see §5. All-vertex conservative numbers in the CSV.

## 4. Findings (measured facts)

1. **The Shifted Containment Law holds exactly on real graphs.** `law_viol=0` and
   `nest_viol=0` on all 9 graphs (per-vertex `κ_s ≥ ∂_s(κ_{s+1})` and support
   nesting). Validation beyond the earlier toy/karate checks, up to 317K vertices.
2. **Compression follows `comp = avgTraj / (2 + avg δ≠0 per vertex)` — fits the
   measurements exactly.** The 2-value anchor is a fixed per-vertex cost, so even at
   100% saturation compression is capped at `avgTraj/2`.
3. **Saturation is graph-class-dependent (13%–89%).** Collaboration networks (`ca-*`)
   saturate well (66–89%); a citation network (`cit-HepPh`) saturates at only 13%.
4. **Compression is modest and not universal: 0.90×–3.67×.** Best when trajectories
   are long (`ca-AstroPh`, dense clique structure, avg traj 10.95 → 3.67×); on
   `cit-HepPh` the anchor overhead makes the encoding **expand** (0.90×).

## 5. Limitations

* **Numerical clamp.** On high-ω dense graphs the s-clique core numbers exceed 2⁵³
  and the SPIN★ peel clamps support at `kBucketKeyClamp=1e18`; the high-s spectrum is
  not exactly representable in double. Such (v,s) are flagged `fuzzy` and excluded
  from the trustworthy stats (`*_exact` columns); they are NOT counted as law
  violations. dblp: only 280/317080 vertices fuzzy. HepPh: 743/12006.
* The compression accounting counts *values*, not bits. A bit-level model (s* in
  ⌈log₂ω⌉ bits, varint deltas) would shift absolute numbers but not the qualitative
  ranking or the `avgTraj/2` ceiling.

## 6. Conclusion / next

The KK anchor+delta encoding is a **correct** but **modest** compressor of the
spectrum (1.2–2× typical, 3.7× best, can expand). The strongest result is the clean
empirical confirmation of the Shifted Containment Law on real graphs. The fixed
2-value anchor is the bottleneck; the promising next lever is **trajectory
equivalence classes** — vertices in the same dense region often share an *identical*
trajectory (all K_m vertices do), so storing each distinct trajectory once and
referencing it could break the `avgTraj/2` ceiling. That is the next experiment to
run before judging the index direction.
