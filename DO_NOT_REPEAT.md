# DO NOT REPEAT

Everything here was already tried, measured, or paid for. Re-doing any of it costs days and yields
what is already written down. Two lists: **dead ends** (ideas that were killed, with the measurement
that killed them) and **traps** (mistakes that produced wrong results, with what they looked like).

Rule of thumb that generated most of this file: **a cheap probe before an expensive build.** Almost
every entry below cost one counter or one run to settle, and would have cost days to discover by
building the thing first.

---

## PART 1 — DEAD ENDS (do not re-attempt without new evidence)

### D1. SLP step 2 "grouped walk" on the 71-163x redundancy figure — the figure was inflated
§203 green-lit it on `witnessProbes / newWitnesses` = 71x (ca-AstroPh) / 163x (web-Google).
That counter increments **before** the ell/u feasibility test, so it lumps bounds-rejects (which
grouping cannot remove) with hash-rejects (which it can). The clean counter already existed on the
fixed-r path (`ppYEnum`, incremented **after** the test) and reads **2.2-9.8x**.
**Partial reversal, do read this:** §214 later measured the per-cell *pop machinery* at **60-84%** of
row peel, which the witness figure never covered. So the **cohort/skeleton** form is green-lit;
only the mode-1 (truncated inclusion-exclusion) variant stays dead on the antichain-blowup risk.

### D2. Active-class list for addDelta's DFS (§208 "G3-lite") — nothing to skip
Hypothesis: the DFS scans all `Mloc` coordinates per node, so a precomputed active list would help.
Measured: **85-93% of scanned coordinates are already usable**, so the ceiling is 1.1-1.4x.
Cost of finding out: one counter pair. addDelta's cost is the enumeration itself, not scan overhead.

### D3. Row-wide polynomial supInit (§214 "layer 2") — 3-8% ceiling
Elegant (support is a coefficient of one generating polynomial, so a whole row could share it) and
it stays in the paper's algorithm statement. But certified patterns already skip supInit in a sweep,
so the measured share is **3-8% of row peel**. Parked; build nothing until something makes it matter.

### D4. Residue cross-cell delta coding — attacks 0.1% of the index
Ranked first on intuition, then measured: residue is **0.0-0.1%** of the index file; the pattern
table is 86-99% (§221). Revisit only if a low-certification graph shows a different anatomy.

### D5. `twinFrac` as the structural predictor of compression — refuted
Proposed as "fraction of vertices in a class of size >= 2 drives compression", supported on 6 graphs,
then **killed by web-it-2004**: twinFrac 0.422 with compression 515x, against ca-CondMat at twinFrac
0.504 with 1.62x. Mechanism: `mult(P) = prod_c C(n_c, b_c)` is convex and explosive in class SIZE,
so one class of 430 outweighs 200 classes of size 2, which twinFrac weights identically.
**Replacement (correct and exact): compression = E_P[prod_c C(n_c, b_c)] (§220).**

### D6. Hybrid guard "if W > tau, fall back to per-clique peeling" — ceiling is parity, and it is
### self-defeating
W is only known after the front end, so the fallback pays our whole front end (which CND never pays)
and then runs something CND-equivalent: **provably <= CND, never better**. Worse, W falls
monotonically in r while CND's #r-cliques explodes, so the guard fires on the small cells where both
finish in seconds and stays off where we actually win. It also contradicts the project's thesis
(`feedback_own_algorithm_not_graft`). Analysis in full at the §215/§216 discussion.

### D7. Precision guarding in the deconvolution supInit — bought nothing, cost speed
A per-coefficient cancellation test fired on **3-44% of LEGITIMATE exact zeros**, each trip costing a
full `supFast` recompute, to insure against a risk that 89M A/B comparisons had measured at **exactly
zero error**. Removed; 6/6 cells still bit-exact without it. Standing instruction: precision issues
that are storage artefacts (a count a `double` cannot hold but Python's unbounded int could) are NOT
logic errors and are to be ignored.

### D8. Reusing the E3 CND logs — they are parallel-CND
26 cell logs sit in `tods2:/home/wenqianz/nsi_e3/`. They look reusable and are not: webit (3,4) reads
316s / 307GB there (a 300GB prlimit abort, i.e. time-to-OOM) against **5581s / 101GB** for a serial
run that COMPLETES. Every CND number produced before 2026-07-21 is quarantined (§217).

### D9. Billion-scale web (it-2004 1.03B edges, uk-2005 783M) — MCE wall
Both are downloaded and converted (`/data/wenqianz/*.full.edges`, verified against their headers).
Both **exceed a 3-hour MCE budget** on tods2. Do not just raise the budget: CND cannot touch those
graphs either, so the outcome would be "neither finishes", which is not a result. If billion scale is
required, the open options are a k-core prefilter (needs justification in the paper) or a faster MCE.

### D10. Graphs excluded from the roster, each for a measured structural reason
| graph | why | number |
|---|---|---|
| com-youtube | real losing cells at r=3 | (3,5) 0.99x, (3,6) 0.77x |
| soc-pokec | twin-free: nothing to quotient | compression **1.001x** |
| com-lj / com-orkut | same social structure | deprioritised after pokec |
| dblp-coauthor | **pattern wall** (a third boundary type) | MCE finishes in ~1h, then incidences blow the 500M cap at both (3,4) and (4,6) |
| email-Eu-core | twin-free | compression 1.002x, 740 classes but no compression |
Social is absent **by structure**, and W says so before any build. That absence is a contribution
(the applicability characterization), not something to fix.

---

## PART 2 — TRAPS (these produced WRONG results, not just slow ones)

### T1. zsh `set -- $spec` does not word-split — silently empty loops
Cost two separate incidents. Once it produced a bit-exactness gate that reported **PASS on 0 lines of
output**. **Rule: loops go in a file and run under `bash`, never inline.** Every gate now checks the
reference output is non-empty before comparing.

### T2. Total-time comparison on MCE-dominated graphs — manufactures false wins AND false losses
Reported "the deconvolution loses on three graphs" from total time. Phase attribution refuted it: its
peel never regressed (0.998x-2.075x). com-amazon's "loss" came from a code path that reported
`deconv=0` and never executed; web-it's came entirely from MCE (5.55 -> 5.96s), which runs before any
deconvolution code; com-youtube's "1.027x win" was the same MCE noise pointing the other way.
**Rule: attribute the metric to the phase you changed, then verify the untouched phases are equal.**

### T3. A number measured under a profiling env var is not comparable to a production number
`PIVOTER_PEEL_PROFILE` hashes every incidence. The first deconvolution figures (125.43s -> 73.42s,
1.71x) were both inflated ~34%; clean is 93.53s -> 46.12s = 2.03x. The ratio survived, the absolutes
did not.

### T4. macOS `maximum resident set size` is noisier than the effect it was being used to measure
Reported a "+14% memory regression" that did not exist. The bisect: RSS differed by 1.83GB at a
checkpoint where **both runs execute identical code**. The allocator high-water mark differed by 1.5%,
and Linux measured the memory ratio at 1.000-1.005. **Use `peak memory footprint` on macOS, and
prefer the Linux numbers for any memory claim.**

### T5. Optimising by reading a fat struct in a hot loop — cache-hostile, and it inverted the result
The clamp pre-filter first read `pats[q].alive/.key` per bucket candidate: random access into ~150MB
of `Pat` structs, replacing a loop that was sequential over `leafFlat` and never touched `pats[]`
until it matched. Result: +13.5% on ca-AstroPh but **-8% on ca-HepPh**. Fixed with a compact
`patLive[q] = alive ? key : -1` array (9.4MB, cache-resident): 1.04-1.36x on all three graphs.
**Rule: when two graphs disagree in SIGN, that is a signal to diagnose, not to average or to pick one.**

### T6. A table sized from data that a new format no longer stores — silent zeros
NSI2 sizes its binomial table from the **stored** pattern `cP` values. NSI3 stores no patterns, so the
table stayed at `smax=8` and `C2(41,1)` silently returned 0: **every certified answer was 0 while the
index structure and its size looked perfect.** A size-only or structure-only check would have called
it a triumph. Caught only by comparing answers.

### T7. A verification workload that is not what it claims to be
`verify-cp` fed r=3 samples to the r=4 and r=5 columns as 4- and 5-tuples, which are not cliques, and
reported **1,577 bogus mismatches**. It now infers r from the workload, and `nsi_query sample R COUNT`
generates real r-cliques per r. Separately, `nsi_baseline sample` wrote its summary to stdout, so the
workload file began with text and `fscanf("%d")` produced an EMPTY workload that scored 0/0 on
everything.

### T8. The same change can flip SIGN between machines on the same graph
Sparse footprints speed up ca-AstroPh's peel 1.15-1.22x on local ARM and regressed it (~0.94x) on
tods2 Xeon. Cache and prefetch behaviour differ by architecture. **Paper numbers come from tods2;
local numbers are direction-finding only.**

### T9. Running two timing benchmarks concurrently
Recorded in the project long before this session (contention once inflated a peel 3.7x and produced a
wrong "does not scale" verdict) and hit again anyway. Timing runs get the machine to themselves.

### T10. "It ran without error" does not mean it measured anything
A dirty server tree makes `git pull` abort silently, so a run rebuilds from OLD source, exits 0 on
every graph, and produces zero instrumentation. **Gate long runs on HEAD hash + a grep for the symbol
you expect + build exit code, before waiting.** The scout scripts do this and it has already caught a
stale-hash launch (`FATAL: HEAD mismatch`).

### T11. Tooling failures that look like experimental results
Hit repeatedly: the server `awk` rejects ternaries inside `printf` args; macOS `bash` 3.2 has no
`declare -A`; the clang LSP reports ~10 bogus errors on the engine because it lacks
`-I../src/NucleusDecomposition`; LAW's own tarball URL 404s (use Maven Central) and `BVGraph` needs
`commons-math3` on the classpath. **Confirm a failure is real before interpreting it.**

### T12. Default budgets silently truncate big graphs
The engine's MCE budget defaults to **120s**, which aborts Queen_4147 (163M edges) and dblp-coauthor
before any measurement. `SCT_MAX_INC` similarly guards pattern explosion (exit 7). Both are correct
behaviours that look like failures; raise them explicitly for big graphs (`--mce-budget`, `SCT_MAX_INC`).

---

## PART 3 — THINGS THAT ARE TRUE AND CHEAP TO RE-CHECK

Read these before proposing anything, because they bound most ideas:
- `compression = E_P[prod_c C(n_c, b_c)]` and `W = hostSz_avg / compression` — both computable from
  the front end alone via `SCT_W_ONLY=1`, which also prints classes / clsMean / clsMax.
- `#classes <= n` always (classes partition the vertices), so classes are ALWAYS far fewer than
  r-cliques. That is **not** compression: email-Eu-core has 740 classes against 102,747 r-cliques and
  compression 1.002x.
- The class alphabet is **r-independent** (pkustk13 and Queen_4147 have identical class counts at
  r=3,4,5); only s enters, through the `|M| >= s` region filter.
- P10 lower-bounds the residue: no algorithm avoids per-cell work there. The index's shape is
  therefore optimal in kind, and that is the right way to present it.

---

## PART 4 — ADDED 2026-07-29 (§225-§229)

### D11. "Port the §210 optimization stack into the plane engine" (old debt #1) — RE-RANKED, measure first
It was the top build-side item for weeks. The plane build's phase split says the §210 stack aims at
the wrong phase: the peel is 13-50% of the plane build, the leaf/pattern maps were 28-52%, and the
MCE + class quotient is **0.3%** (pkustk13: 1.3s of 415.7s). That last number also retires the
standing suspicion that maximal-clique enumeration is the plane engine's bottleneck. §228 fixed the
map waste instead; re-measure the split before porting anything.

### D12. T6 host-1 certification above the first row — worth nothing
The diagonal certificate already certifies everything host-1 would, everywhere except r=rMin (which
has no previous row to transfer from). Measured on ca-GrQc: host-1-AND-still-uncertified is 44.5% at
r=3 and **0.00% at r=4 and r=5**. Ceiling on the graph where it matters most (pkustk13) is ~4% of the
build, against a soundness argument that has to be made at the universal-pattern level. §229.

### T13. Sorting a tuple IN PLACE inside a recursion that still owns its prefix
The archive expansion sorted `cur[0..r)` to emit a canonical key, which permutes the blocks the OUTER
recursion levels still own, so every tuple after the first came out corrupted. It presented as a 3.6%
"miss" rate on a workload that is 100% real cliques -- i.e. as a plausible-looking property of the
data, not as a crash. **Sort into a separate buffer.** Latent in `nsi_baseline.cpp` since it was
written; found only because a pattern-sampled workload MUST hit an archive built from those patterns.

### T14. A gate whose reference file is in the wrong FORMAT reports NO-REF, not FAIL
`rowfile R FILE` wants R vertices per line; the bench format is `R v1 ... vR`. Feeding one to the
other produces an empty reference, and a gate that only compares outputs would have called that a
pass. Every gate here checks the reference is non-empty FIRST (the same rule T1 produced), and it
fired immediately.

### Cheap probes that are now one command
- `nsi_query INDEX anatomy` — exact byte breakdown of a loaded index plus the price of a packed
  encoding, self-checked against the file size. Run this BEFORE proposing any format change.
- `nsi_query NSI2 archive --vs SLIM` — archive size, index size, and the ratio, no expansion.
- `[nsi-plane-col] ... residue=N leaf-maps=built|SKIPPED` and `[nsi-plane-probe §229]` — how much of
  a plane column is certified, and how much host-1 would add, printed by every build.
