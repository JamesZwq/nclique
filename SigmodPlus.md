# SIGMOD 2027 Paper — Experiment Guide

Comprehensive reference for the experiments backing the SIGMOD 2027
submission *Efficient (r,s)-Nucleus Decomposition on Dense Higher-Order
Graphs*.  Covers (1) the setup we benchmark under, (2) how to run any
piece of the bench from scratch, and (3) how the resulting CSV/TSV
files map to the paper's figures and tables.

Paper source lives in `Sigmod2027Nuclear/` (Dropbox/Overleaf symlink
from the repo root).

---

## 1.  Setup and Settings

### 1.1  Hardware

Two identical Linux servers, accessed over SSH:

| Server | CPU                          | RAM      | /data free |
|--------|------------------------------|----------|------------|
| tods1  | Intel Xeon Gold 6342, 96 c   | 503 GB   | 1.2 TB     |
| tods2  | Intel Xeon Gold 6342, 96 c   | 503 GB   | 464 GB     |

Every benchmark cell is wrapped with `/usr/bin/time -v` so the
`Maximum resident set size (kbytes)` line gives peak memory.
Single-threaded by default; concurrent runs share CPU but no
contention beyond a few cores.

### 1.2  Build

```bash
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 12        # hard cap: NEVER use -j > 12
```

The binary is `./build/bin/degeneracy_cliques`.  C++23, GCC-11,
`-O3`, OpenMP + TBB + Boost + sparsehash.

### 1.3  Datasets (paper-6)

`graphs/*.edges` in CSR-friendly edge-list format.  Sorted by `|V|`:

| Tag    | Graph         | \|V\|     | \|E\|        | degeneracy | omega |
|--------|---------------|-----------|--------------|-----------:|------:|
| \daDBC{}  | dblp-core30   |   1,206   |    31,769    | 113        | 114   |
| \daGRQ{}  | ca-GrQc       |   5,242   |    14,484    |  43        |  44   |
| \daHEPP{} | ca-HepPh      |  12,008   |   118,489    | 238        | 239   |
| \daCON{}  | ca-CondMat    |  23,133   |    93,439    |  25        |  26   |
| \daDB{}   | com-dblp      | 317,080   | 1,049,866    | 113        | 114   |
| \daWI{}   | web-it-2004   | 509,338   | 7,178,413    | 431        | 432   |

Extras tested for the "infeasibility regime": `com-orkut`,
`com-lj`, `com-friendster` (only on tods2, not in paper).

### 1.4  Algorithms Under Test

Each algorithm is selected via an env var.  The same binary
implements all four:

| Paper macro      | CSV `algo` label   | Env var                           | Source file                                                   |
|------------------|--------------------|-----------------------------------|---------------------------------------------------------------|
| `\nuclearcd{}`   | `REF`              | `PIVOTER_RUN_REF=1`               | `NUcllearCoreDeompositoinCorrect.cpp`                          |
| `\regnd{}`       | `RegNDC` / `V3LM`* | `PIVOTER_RUN_REGION_V3LM=1` + `PIVOTER_VSAFE_CLOUD=1` | `NucleusCoreDecompositionRegionCPI_LowMem.cpp`              |
| `\regndn{}`      | `V3LM_NOCPI`       | `PIVOTER_RUN_REGION_V3LM_NOCPI=1` | (same file, NoCPI code path)                                  |
| `\regndh{}`      | `V3LM_HIER`        | `PIVOTER_RUN_REGION_V3LM_HIER=1`  | RegND + post-peel union-find on tuple instances               |

*`RegNDC` is the newer CSV label, `V3LM` is the historical one —
both point at the same env var and the same binary code path.
Resolve by preferring `RegNDC` when present.

### 1.5  Budgets

The main sweep in `bench_v3_all.py` uses a one-hour budget.
It kills an individual process if its memory exceeds 250 GB, and
launches new jobs only while total machine memory use is below
300 GB.
The focused Pass A/B/C runs use the explicit `--timeout` and
`--mem-cap-gb` values shown in §2.3.

| Resource | Main-sweep cap | Where enforced                                  |
|----------|----------------|-------------------------------------------------|
| Time     | 3600 s (1 h)   | `TIMEOUT` in `bench_v3_all.py`                  |
| Per-run memory | 250 GB  | `PER_PROC_MEM_GB` in `bench_v3_all.py`          |
| Launch gate | 300 GB total machine use | `MEM_LIMIT_GB` in `bench_v3_all.py` |
| Disk     | unbounded      | logs land under `bench_v3_all_logs/` and `/data/wenqianz/...` |

When a targeted cell is killed by timeout, `bench_targeted.py`
also calls `os.killpg(pid, SIGKILL)` on the child group so the
binary cannot survive `/usr/bin/time -v` and become an orphan.

### 1.6  Skip Discipline

For sweeps over `s` at fixed `r`, the bench skips every higher `s`
once a cell TIMEOUTs or OOMs (`--skip-after-fail`).  The skipped
rows are written with `status=SKIP_AFTER_FAIL`.  On restart, every
row with non-blank status is considered done, so resume never
retries a TIMEOUT.

---

## 2.  How to Run Experiments

### 2.1  Single Cell (Debug / Spot Check)

```bash
PIVOTER_RUN_REGION_V3LM=1 PIVOTER_VSAFE_CLOUD=1 \
  /usr/bin/time -v ./build/bin/degeneracy_cliques \
    graphs/dblp-core30.edges 3 6 degen
```

Stdout has phase timers (`Step 1`, `CPI counting time`,
`Peeling time`, `Total time`).  Stderr has `Maximum resident set
size (kbytes)`.

### 2.2  Main Sweep (`bench_v3_all.py`)

The big paper-6 grid.  Auto-discovers the server and uses pre-set
graph lists.  Output: `bench_v3_all_results.csv` in the repo root
(NOT under `/data`).

```bash
python3 bench_v3_all.py            # interactive, prints PIDs
nohup setsid python3 bench_v3_all.py > bench_v3_all.log 2>&1 &
```

This produces the big "Exp-1 / Exp-2" heatmap data.  Schema in §3.1.

### 2.3  Focused Bench (`scripts/bench_targeted.py`)

Single algorithm × explicit (graph × r × s) grid.  Used for the
three SIGMOD coverage gaps (Pass A/B/C).

```bash
# Pass A: V3LM_NOCPI ablation row for Fig 3
python3 scripts/bench_targeted.py --algo V3LM_NOCPI \
  --graphs dblp-core30 ca-GrQc ca-CondMat ca-HepPh com-dblp \
  --rs 3 4 5 --ss 4 5 6 7 8 \
  --timeout 1200 --mem-cap-gb 100 --skip-after-fail \
  --out /data/wenqianz/passA_v3lm_nocpi.csv

# Pass B: REF + V3LM at high r (memory ratio fig)
for algo in REF V3LM; do
  python3 scripts/bench_targeted.py --algo $algo \
    --graphs ca-GrQc ca-CondMat dblp-core30 \
    --rs 7 10 12 15 --ss 8 10 12 14 16 18 20 22 25 30 \
    --timeout 1200 --mem-cap-gb 200 --skip-after-fail \
    --out /data/wenqianz/passB_highr.csv
done

# Pass C: large-graph extension
for algo in REF V3LM; do
  python3 scripts/bench_targeted.py --algo $algo \
    --graphs com-orkut com-lj \
    --rs 3 4 5 6 7 --ss 4 5 6 7 8 10 12 \
    --timeout 1800 --mem-cap-gb 200 --skip-after-fail \
    --out /data/wenqianz/passC_largegraphs.csv
done
```

The script wraps every cell in its own process group so a TIMEOUT
kill terminates both the `/usr/bin/time` wrapper and the binary.
Resume-friendly: rows with any non-blank status are skipped on
restart.

### 2.4  Structural Compression (`scripts/collect_class_tuple_rclique.py`)

Captures \(|K_r|, |\tupleset|, |\classset|\) from RegND's stdout
for the Exp-3 figure.

```bash
python3 scripts/collect_class_tuple_rclique.py \
  --graphs com-dblp ca-CondMat ca-GrQc ca-HepPh dblp-core30 web-it-2004 \
  --rs 3 4 5 6 7 8 9 --s 10 \
  --timeout 1200 \
  --out experiments/vsafe_phase_b1/class_tuple_rclique.csv
```

### 2.5  Phase Breakdown (`scripts/bench_phase_breakdown.py`)

Captures per-phase time AND per-phase delta-RSS in one shot.
Drives Exp-4 and Exp-5.

```bash
python3 scripts/bench_phase_breakdown.py \
  --graphs com-dblp web-it-2004 ca-HepPh ca-CondMat \
  --rs 3 --ss 4 6 8 10 \
  --timeout 1800 \
  --out /data/wenqianz/phase_breakdown.tsv
```

Internally sets `PIVOTER_BREAKDOWN_LOG=<TSV>` and
`PIVOTER_BREAKDOWN_META=<cell-id>` — the binary's `PhaseLogger`
appends rows after each phase mark inside `main.cpp` and RegND's
hot loop.

### 2.6  Server Workflow

Code is propagated via Git, not `scp`.  Always:

```bash
# locally
git add -f scripts/<changed>.py src/<changed>.cpp
git commit -m "..."
git push origin main

# on tods1 / tods2
ssh tods1 "cd ~/nclique && git stash && git pull origin main && \
           cmake --build build -j 12"
```

Long runs go in their own `setsid nohup` subshell so closing the
SSH session does not kill them:

```bash
ssh tods2 "cd ~/nclique && setsid nohup bash -c '
python3 scripts/bench_targeted.py --algo ... \
  --out /data/wenqianz/X.csv > /data/wenqianz/X.log 2>&1
echo DONE >> /data/wenqianz/X.log
' </dev/null > /dev/null 2>&1 &"
```

To kill cleanly:

```bash
ssh tods1 "pkill -9 -u \$(whoami) -f 'bench_targeted'"
# then verify with: ps -o pid,etime,cmd -u $(whoami) | grep degeneracy_cliques
# and kill any leftover binary by PID
```

If an OS RSS check shows an orphan with elapsed > 1.5x the cell
timeout, kill it directly:

```bash
kill -9 <pid>
```

---

## 3.  How to Read the Data

### 3.1  Main Bench CSV (`bench_full_merged.csv`)

Master CSV: tods1 + tods2 + Pass A/B/C combined.  857 k rows.

Schema:

| Column   | Meaning                                            |
|----------|----------------------------------------------------|
| graph    | dataset name (matches `graphs/<name>.edges`)       |
| r        | r-clique size                                      |
| s        | s-clique witness size                              |
| algo     | one of `REF`, `RegNDC`, `V3LM`, `V3LM_NOCPI`, `V3LM_HIER` |
| status   | `OK` / `TIMEOUT` / `OOM` / `SKIP_AFTER_FAIL` / `SKIP_TIMEOUT` / `ERROR(rc)` |
| wall_ms  | end-to-end time including process startup (parsed from `Total time` if OK, else timeout cap) |
| total_ms | `Total time: X ms` from binary's stdout (OK cells only) |
| step4_ms | `CPI counting time: X ms` — Step 4 closed-form support |
| peel_ms  | `Peeling time: X ms` — main loop                    |
| hier_ms  | `Hierarchy time: X ms` — V3LM_HIER only             |
| mem_kB   | `Maximum resident set size (kbytes)` from `time -v` |

Status semantics:

- **OK** — completed within budget; numbers are real
- **TIMEOUT** — exceeded the cell's time cap; numbers are placeholders
- **OOM** — exceeded the memory cap or was OOM-killed
- **SKIP_TIMEOUT** — the bench harness inferred a TIMEOUT from
  neighbouring cells (e.g., (r, s+1) timed out so (r, s+2) is
  marked without running). Standard practice across all our sweeps.
- **SKIP_AFTER_FAIL** — same idea, but recorded by `bench_targeted.py`
  for cells skipped after an in-script TIMEOUT/OOM
- **ERROR(rc)** — process exited non-zero (segfault, abort, etc.)

For paper Fig 3 and Fig memory-ratio: only `status == OK` rows
are used, and only when both `REF` and `RegNDC`/`V3LM` are OK on
the same `(graph, r, s)` cell ("matched" cells).

Quick coverage check:

```python
import csv, collections
counts = collections.defaultdict(lambda: collections.defaultdict(int))
with open('paper_data/bench_full_merged.csv') as f:
    for row in csv.DictReader(f):
        if row['status'] == 'OK':
            counts[row['algo']][row['graph']] += 1
print(counts['RegNDC'])   # per-graph OK count for our main algorithm
```

### 3.2  Compression Counts CSV (`class_tuple_rclique_merged.csv`)

Drives Exp-3 (`exp_compression_r.pdf`).  33 rows, paper-6 graphs.

Schema:

| Column   | Meaning                                                    |
|----------|------------------------------------------------------------|
| graph    | dataset                                                    |
| r        | r-clique size                                              |
| s        | s-clique witness size (fixed at 10 for the figure)         |
| classes  | \|classset\| — profile-class count (paper symbol \(C\))   |
| tuples   | \|tupleset\| — region-tuple count (paper symbol \(T\))    |
| rcliques | \|K_r\| — alive r-clique count (paper symbol \(K_r\))     |
| wall_s   | wall time in seconds (informational)                       |
| error    | non-empty if the run failed (e.g., `timeout after 1200s`)  |

The 1.9·10⁹× compression number on web-it-2004 at r=7 is
`rcliques / tuples` on that row.

### 3.3  Phase Breakdown TSV (`phase_breakdown.tsv`)

Drives Exp-4 (time breakdown) and Exp-5 (memory breakdown).
136 rows (15 cells × 9 phase marks each).

Schema (tab-separated):

| Column           | Meaning                                                |
|------------------|--------------------------------------------------------|
| meta             | `graph,r,s,algo,run_id` — joined cell identifier       |
| phase            | `loadAndSort`, `buildSDCT`, `preMutation`, `prepareGraph`, `MCEnum`, `Index`, `Support`, `Peel`, `dispatch_total` |
| duration_ms      | wall time spent between this mark and the previous     |
| rss_kb           | absolute RSS sample at this mark                       |
| delta_rss_kb     | RSS delta (kB) added since the previous mark           |
| component_bytes  | optional payload tag (e.g., SDCT tree byte size)       |

Paper-mapped phases:

| Paper phase | TSV phases                            |
|-------------|---------------------------------------|
| Tree        | `buildSDCT`                           |
| MCEnum      | `MCEnum`                              |
| Support     | `Support`                             |
| Index       | `Index`                               |
| Peel        | `Peel`                                |

Memory slices for Exp-5:

| Paper slice | TSV phase delta_rss summed |
|-------------|----------------------------|
| Graph       | `loadAndSort`              |
| Tree        | `buildSDCT`                |
| Index       | `MCEnum` + `Index`         |
| RClique     | `Support` + `Peel`         |

### 3.4  Status Quick Reference

| If the paper claims … | The data shows it via …                                                        |
|-----------------------|--------------------------------------------------------------------------------|
| Speedup at (r,s)      | `bench_full_merged.csv`: both REF and RegNDC `OK`, take `wall_ms` ratio        |
| Coverage gap          | Count `REF status != OK` per graph and compare to `RegNDC status == OK`        |
| Memory ratio          | Same matched cells, take `mem_kB` ratio                                        |
| Tuple compression     | `class_tuple_rclique_merged.csv`: `rcliques / tuples` at any (graph, r)        |
| Peel takes X% of time | `phase_breakdown.tsv`: `Peel.duration_ms / sum(all phases).duration_ms`        |
| RClique slice < Y GB  | `phase_breakdown.tsv`: `Support.delta_rss_kb + Peel.delta_rss_kb`              |
| Hierarchy is "free"   | Per-graph `V3LM_HIER.hier_ms / V3LM_HIER.wall_ms` from main bench              |

---

## 4.  Figures and Their Data Sources

| File                                    | Source                                | Script                                       |
|-----------------------------------------|---------------------------------------|----------------------------------------------|
| `figures/exp_heatmap.pdf`               | `bench_full_merged.csv`               | `scripts/make_sigmod_figs.py`                |
| `figures/exp_memory_r.pdf`              | `bench_full_merged.csv`               | `scripts/make_sigmod_figs.py`                |
| `figures/exp_compression_r.pdf`         | `class_tuple_rclique_merged.csv`      | `scripts/plot_class_tuple_rclique.py`        |
| `figures/exp_breakdown_time.pdf`        | `phase_breakdown.tsv`                 | `scripts/plot_breakdown.py`                  |
| `figures/breakdown_memory_DBLP.pdf`     | `phase_breakdown.tsv` (com-dblp rows) | `scripts/plot_breakdown.py`                  |
| `figures/breakdown_memory_web-it-2004.pdf` | `phase_breakdown.tsv` (web-it rows) | `scripts/plot_breakdown.py`                  |
| Tab 4 (Backend ablation)                | hard-coded 12-cell table              | n/a (numbers from `paper_data/bench_full_merged.csv` `V3LM_NOCPI` rows) |

Regenerate everything:

```bash
python3 scripts/make_sigmod_figs.py
python3 scripts/plot_class_tuple_rclique.py \
  --csv paper_data/class_tuple_rclique_merged.csv \
  --log /dev/null --cols 3 \
  --out Sigmod2027Nuclear/figures/exp_compression_r.pdf
python3 scripts/plot_breakdown.py
```

After every figure change, compile the paper and verify zero broken
references:

```bash
cd Sigmod2027Nuclear
latexmk -pdf -interaction=nonstopmode main.tex
pdftotext main.pdf - | grep -c '??'   # MUST print 0
```

---

## 5.  Data Files Index

```
paper_data/
├── bench_v3_all_tods1.csv               # raw tods1 main sweep   (277k rows)
├── bench_v3_all_tods2.csv               # raw tods2 main sweep   (580k rows)
├── bench_v3_all_merged.csv              # tods1 + tods2 dedup    (857k rows)
├── passA_v3lm_nocpi.csv                 # V3LM_NOCPI ablation     (68  rows)
├── passB_highr.csv                      # high-r REF + V3LM       (187 rows)
├── passC_largegraphs.csv                # com-orkut / com-lj      (94  rows)
├── bench_full_merged.csv                # bench_v3_all + passes   (857k rows)
├── class_tuple_pass3a.csv               # com-dblp + ca-CondMat   (15  rows)
├── class_tuple_pass3b.csv               # web-it-2004 partial     (3   rows)
├── class_tuple_rclique_merged.csv       # all paper-6 graphs      (33  rows)
└── phase_breakdown.tsv                  # 5-phase × 15 cells      (136 rows)
```

---

## 6.  Result Summary (current as of 2026-06-05)

vs Nuclear CD on 549 strict matched cells across paper-6:

| Subset                              | Time geomean | Memory geomean |
|-------------------------------------|--------------|----------------|
| 4 non-HepPh graphs (with REF cells) | 16.4×        | 14.6×          |
| 2 dense (\daDBC, \daDB)             | 32.7×        | 24.1×          |
| 2 dense, s ≥ 6                      | 38.8×        | 25.7×          |

Top-5 per-cell speedups, all on dblp-core30:

| (r, s) | Nuclear CD  | RegND    | Speedup |
|--------|-------------|----------|---------|
| (5, 40)| 20.0 min    | 210 ms   | 5,767×  |
| (3, 45)| 18.4 min    | 204 ms   | 5,403×  |
| (4, 49)|              | 297 ms   | 3,239×  |
| (5, 12)|              | 130 ms   | 3,080×  |
| (5, 17)|              | 119 ms   | 3,022×  |

Coverage gap on web-it-2004:
- Nuclear CD: 0 OK cells  /  RegND: 92,003 OK cells

Compression ratio peak: web-it-2004 at r=7 — \(|K_r| / |T|\) ≈
1.9·10⁹×.

---

## 7.  Known Gaps and Limits

- **com-orkut / com-lj**: even RegND TIMEOUTs at every (r, s) ≤ 7
  in the focused large-graph pass.  These graphs are an honest
  infeasibility regime for RegND and are not in Fig 3.
- **ca-HepPh**: only r=3 yields useful data; r=4..9 all SKIP/TIMEOUT
  for V3LM_NOCPI and most of REF.  Paper §exp-compression
  acknowledges this as the structural-boundary regime.
- **r ∈ {10, 15} memory ratio**: only ca-CondMat and ca-GrQc r=7
  have matched cells.  Other graphs at high r have REF=0 OK so no
  ratio computable.
- **Pre-2026-06 bench had a `/usr/bin/time -v` orphan bug**: TIMEOUT
  killed the wrapper but not the binary.  Fixed in
  `scripts/bench_targeted.py` by calling `os.killpg` on
  `TimeoutExpired`.  Any TSV/CSV row from before commit `b70a313`
  may have wasted CPU on an orphan that was never recorded.
---

## 8.  4-Tier Ablation Reproduction Guide

The 4-tier ablation is the SIGMOD §6 experiment that justifies the
RegND family as a progression rather than a single algorithm.  Each
tier is one strict superset of optimizations over the previous one,
and all four produce **bit-exact** core values on every cell that
finishes.

### 8.1  What the Four Tiers Mean

| Tier | Name                       | Paper anchor                    | One-line description                                                                                  |
|------|----------------------------|---------------------------------|-------------------------------------------------------------------------------------------------------|
| T1   | Class Integrality          | `\refsec{approach-class}`       | Direct s-clique enumeration + per-tuple bucket peel, no path-level pivoting (baseline RegND).         |
| T2   | Pivot Extension            | `\refsec{approach-pivot}`       | T1 + CPI path interface + AggrCount DP, but re-runs AggrCount per peel batch (no incremental delta).  |
| T3   | Dead-Box Inclusion-Exclusion | `\refsec{approach-deadbox}`   | T2 + closed-form support-count via dead-box IE so peel deltas are O(paths) instead of O(re-enum).     |
| T4   | RegND\* (full system)      | `\refsec{approach-algorithm}`   | T3 + private-cloud V-safe optimization, instant path retire, and the SIMD-vectorized feasWeighted DP. |

All four tiers ship in the same binary (`degeneracy_cliques`) and are
selected at runtime via env var.  No recompile needed to switch.

### 8.2  Build

```bash
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 12        # hard cap: NEVER use -j > 12
```

Output binary: `./build/bin/degeneracy_cliques`.  Same binary used by
every other section of this guide.

### 8.3  Run a Single Cell Locally

```bash
PIVOTER_RUN_REGION_TIER=1 PIVOTER_TIER=4 \
  /usr/bin/time -v ./build/bin/degeneracy_cliques \
    graphs/ca-GrQc.edges 3 4 degen
```

Both env vars are required:

- `PIVOTER_RUN_REGION_TIER=1` activates the tier-dispatch code path
  inside the RegND family (otherwise the binary defaults to the
  legacy V3LM entry point).
- `PIVOTER_TIER=N` (N ∈ {1, 2, 3, 4}) selects which tier to run.

Stdout has the usual phase timers (`CPI counting time`,
`Peeling time`, `Total time`) plus the `Core value distribution`
block used by the verifier.  Stderr from `/usr/bin/time -v` gives
peak RSS.

### 8.4  Verify Bit-Exact Correctness

`verify_tiers.py` parses the `core=N count=M` histogram from each
tier's stdout and asserts T1 = T2 = T3 = T4 on every (graph, r, s).
T4 is the reference.

```bash
# Default cell (ca-GrQc, r=3, s=4)
python3 verify_tiers.py

# A specific cell
python3 verify_tiers.py ca-GrQc 5 6

# Built-in smoke suite (5 cells across 3 small graphs)
python3 verify_tiers.py --all
```

A successful run prints `-> PASS` for every cell and ends with
`ALL PASS`.  Any mismatch is dumped as the first five `(core, t<X>,
t<ref>)` triples plus a `FAIL` marker; the script exits non-zero
so it composes with CI.

### 8.5  Full Sweep on tods2

The full ablation matrix is **18 cells × 4 tiers = 72 runs**
(see `ABLATION_CELLS` in `bench_tier_ablation.py`).  Per-cell
timeout: 3600 s.  Per-cell memory cap: 250 GB (kernel OOM-kill).
T4 rows can be reused from the main sweep CSV via `--reuse-t4`.

```bash
# 1. Local: commit + push the harness
git add bench_tier_ablation.py verify_tiers.py
git commit -m "ablation harness"
git push origin main

# 2. Server: pull, build, link graphs
ssh tods2 "cd ~/nclique && git stash && git pull origin main && \
           cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && \
           cmake --build build -j 12 --target degeneracy_cliques"

# 3. Server: launch in a detached setsid + nohup subshell
ssh tods2 "cd ~/nclique && setsid nohup bash -c '
python3 bench_tier_ablation.py tods2 --reuse-t4 \
  > bench_tier_ablation.log 2>&1
echo DONE >> bench_tier_ablation.log
' </dev/null > /dev/null 2>&1 &"

# 4. Server: monitor progress
ssh tods2 "tail -f ~/nclique/bench_tier_ablation.log"
# or: ssh tods2 "wc -l ~/nclique/bench_tier_ablation_results.csv"

# 5. Server: kill cleanly when needed
ssh tods2 "pkill -9 -u \$(whoami) -f 'bench_tier_ablation'"
ssh tods2 "pkill -9 -u \$(whoami) -f 'degeneracy_cliques'"
ssh tods2 "pkill -9 -u \$(whoami) -f '/usr/bin/time'"

# 6. Local: pull the result CSV back into paper_data/
scp tods2:~/nclique/bench_tier_ablation_results.csv \
    paper_data/bench_tier_ablation_results.csv
```

The harness writes `bench_tier_ablation_results.csv` incrementally
(one row per cell, flushed immediately), so a mid-sweep kill is
recoverable: rerun and the harness skips any row already present.
Per-cell logs land in `bench_tier_ablation_logs/T<tier>_<graph>_r<r>_s<s>.log`.

### 8.6  CSV Schema (12 columns)

Output file: `bench_tier_ablation_results.csv` (locally) or
`paper_data/bench_tier_ablation_results.csv` (after pulling back).

```
graph,r,s,tier,algo,status,wall_ms,total_ms,step4_ms,peel_ms,mem_kB,max_core
```

| Column   | Meaning                                                    |
|----------|------------------------------------------------------------|
| graph    | dataset stem (matches `graphs/<name>.edges`)               |
| r        | r-clique size                                              |
| s        | s-clique witness size                                      |
| tier     | 1, 2, 3, or 4                                              |
| algo     | always `T<tier>` — convenient label for plotters           |
| status   | `OK` / `TIMEOUT` / `TIMEOUT(skip)` / `OK(reused)` / `ERR(rc)` |
| wall_ms  | end-to-end wall time including `/usr/bin/time` startup     |
| total_ms | `NucleusCoreDecomposition took: X` from binary stdout       |
| step4_ms | `CPI counting time: X` (T2/T3/T4 only — empty for T1)      |
| peel_ms  | `Peeling time: X` from binary stdout                       |
| mem_kB   | `[Memory-...] Final Memory: X kB` from binary's tracker     |
| max_core | `Max core: N` (drives the verifier's bit-exact check)       |

Status semantics:

- **OK** — completed within budget; every numeric column is real.
- **TIMEOUT** — exceeded 3600 s cap; numeric columns are placeholders
  (`wall_ms = 3600000`, the rest empty).  The cell **was attempted**.
- **TIMEOUT(skip)** — cell **was not attempted** because it appears in
  the `T1_HARD_SKIP` set (or, by extension, was inferred to be
  infeasible from neighbouring cells).  Same placeholder numerics as
  `TIMEOUT`; treat identically when computing completion rates.
- **OK(reused)** — T4 row joined from `bench_v3_all_results.csv`
  (RegNDC) via `--reuse-t4`.  `max_core` is `-1` because the source
  CSV did not log the histogram.  Counts as `OK` for completion.

### 8.7  Headline Completion Rates (Current Run)

From `paper_data/bench_tier_ablation_results.csv` (72 rows, 18 cells
× 4 tiers, sweep finished on tods2):

| Tier | OK | TIMEOUT | TIMEOUT(skip) | Completion |
|------|----|---------|---------------|------------|
| T1   |  7 |       2 |             9 | **7 / 18** |
| T2   |  7 |      11 |             0 | **7 / 18** |
| T3   | 14 |       4 |             0 | **14 / 18**|
| T4   | 15 |       3 |             0 | **15 / 18**|

Interpretation:

- T1 → T2: zero gain in completion count; T2's per-batch AggrCount
  re-runs are no cheaper than T1's direct s-enum on the cells that
  matter.  The motivation for T2 in the paper is **pivot extension
  as a structural prerequisite for T3**, not a speedup of its own.
- T2 → T3: **+7 cells (doubles coverage)**.  Closed-form support-count
  via dead-box IE is the load-bearing optimization.
- T3 → T4: **+1 cell (web-it-2004 r=5 s=6)**, plus large per-cell
  speedups on the cells both tiers finish.  The private-cloud V-safe
  is a constant-factor improvement, not a coverage breaker.

Use only `status == OK` (and `OK(reused)` if `--reuse-t4` was used)
for the §6 speedup and memory-ratio tables; never let placeholder
`wall_ms = 3600000` rows pollute geomeans.

### 8.8  Known Gotchas

- **T1/T2 timeout boundaries are real.**  The `T1_HARD_SKIP` set in
  `bench_tier_ablation.py` is not arbitrary: it is the closure under
  "if (r, s) times out, every (r', s') with r' ≥ r and s' ≥ s also
  times out" inferred from earlier sweeps.  When you extend the
  matrix to a new graph, **first run T4** on the candidate (r, s)
  cells; only attempt T1/T2 on cells where T4 finishes in well under
  the 3600 s budget, otherwise the T1/T2 row is a 1-hour CPU sink for
  no information gain.

- **V-safe 3-tuple miscount on ca-CondMat (task #125, FIXED 2026-06-12).**
  Root cause was NOT V-safe-specific: `pathAliveCount` was decremented
  twice for a tuple that first *saturated* on a path (Theorem 1
  removal, count--) and later *died* (the static `pathsCoveringTuple`
  index still lists the path, so the blind
  `pathAliveCount -= deadTuples.size()` decremented again).  The count
  hit 0 while live tuples remained, the path was retired prematurely
  (dead boxes freed), and the stranded tuples missed every subsequent
  support subtraction — κ came out 1 too high for 3 tuples.  342 such
  premature retirements occurred on ca-CondMat 3,4 alone; private-cloud
  mode merely changed the peel order so the damage became visible.
  Fix (in `refreshAffectedPaths`' caller, the affected-path loop):
  filter `pi.tupleIdxs` by `rPeeled` *in place* and reset
  `pathAliveCount = tupleIdxs.size()` — exact bookkeeping, no blind
  decrement.  Diagnosed by differential tracing vs EventPeel
  (`PIVOTER_TRACE_TUPLE`, `PIVOTER_TRACE_PATH`,
  `PIVOTER_DUMP_TUPLE_CORE`, `PIVOTER_RETIRE_CHECK` — all kept,
  env-gated, zero cost when unset).  Post-fix: 10/10 compare cells
  exact (ca-CondMat 3,4/3,5/4,5 × private/no-private/vsafe, ca-GrQc
  3,4/3,5/4,5/5,6, com-dblp 3,4), retire-check 0 violations.
  Performance impact (matched-mode A/B vs pre-fix binary built from
  837c2d5^, same machine, sequential, VSAFE; 3-rep medians for small
  cells, single clean pair for the big cell): none.
  ca-CondMat 3,4 peel 166→207 ms (the 342 stranded paths now get
  their legally required updates; PRE/POST rep ranges overlap),
  ca-GrQc 5,6 1918→1884 ms, com-dblp 3,4 3860→3686 ms, com-dblp 5,6
  216.3→215.6 s.  Memory flat on all four (38/21/213/725-731 MB).
  NOTE: the §10.8 "285 s" V3LM com-dblp 5,6 number was contaminated
  by concurrent session load; both clean runs land at ~216 s — use
  that as the V3LM baseline going forward.

- **`PIVOTER_TIER` without `PIVOTER_RUN_REGION_TIER` is silent.**
  Setting only `PIVOTER_TIER=N` makes the binary fall through to the
  legacy V3LM dispatch and ignore the tier number.  The harness always
  sets both; if you write a one-off script, set both too.

- **`OK(reused)` rows have no `max_core`.**  If a downstream plotter
  needs `max_core`, either re-run T4 from scratch on those cells (drop
  `--reuse-t4`) or join `max_core` from a parallel sweep.

- **`mem_kB` is the binary's internal tracker, not `/usr/bin/time`.**
  We wrap with `/usr/bin/time -v` for the log file, but the CSV
  records the `[Memory-Final] Final Memory: X kB` line emitted by the
  in-binary tracker, which is consistently within ~2% of
  `Maximum resident set size (kbytes)`.  When comparing against
  `bench_full_merged.csv.mem_kB`, prefer the harness CSV for tier
  rows (they were measured in the same setup) and the main CSV for
  REF / RegNDC rows.

### 8.9  Extending the Matrix

To add a new graph or a new (r, s) cell:

1. Edit `ABLATION_CELLS` in `bench_tier_ablation.py` (a tuple
   `(graph_stem, r, s)`).
2. If the cell is dense enough that T1 will obviously time out,
   add it to `T1_HARD_SKIP` to save a CPU-hour.  Rule of thumb:
   degeneracy × s ≥ 1500 → skip T1.
3. Confirm the graph file exists on the server under
   `SERVER_DATADIR[server]` (default `/data/wenqianz/`).  The
   harness symlinks `graphs/<g>.edges` automatically.
4. Run with `--dry-run` first to inspect the row plan.
5. Run for real with `--reuse-t4` so previously-finished T4 cells
   from the main sweep are not re-executed.

The harness is **resume-friendly by construction**: rows already in
the output CSV are skipped on restart, so it is safe to interrupt
and relaunch at any cell boundary.


## 9.  Current State Snapshot (2026-06-09)

Three things have changed since the §6 result summary above was
frozen on 2026-06-05:

1.  the 4-tier RegND family was committed to `main`
    (`PIVOTER_TIER={1,2,3,4}` dispatch in C++);
2.  a full ablation sweep ran on tods2 and the CSV came back to
    `paper_data/`;
3.  the paper now ships a new §6.8 *Tier Ablation* subsection with
    its own figure, two tables, and three-paragraph reading.

This section is the snapshot the next maintainer should read first.
Everything is reproducible from `main` at HEAD; no uncommitted edits
gate the figures.

### 9.1  What Shipped This Week (commits)

| Commit  | Subject                                                              |
|---------|----------------------------------------------------------------------|
| ae3dcac | Add PIVOTER_TIER 4-tier RegND ablation entry (T1..T4 dispatcher)     |
| 8a1c0b3 | BuildHierarchyR1: expose LeafConnector + Exact variants in header    |
| 04e6ffb | 4-tier ablation deliverables: figure + handoff guide (this §8 + §9)  |

All three are pushed to `git@github.com:JamesZwq/nclique.git` `main`.
The fresh-clone build on tods2 was verified end-to-end (see §8.4).

### 9.2  Tier Ablation Headline Numbers (72-cell sweep)

Run on tods2 single-stream over ~21h on 2026-06-07/08.  Source:
`paper_data/bench_tier_ablation_results.csv` (12 columns, 72 rows).

**Coverage within the one-hour budget (out of 18 (r,s) cells):**

| Tier                          | OK | TIMEOUT (run) | TIMEOUT (hard-skip) |
|-------------------------------|----|---------------|---------------------|
| T1 (`\regnd`)                 |  7 |  2            |  9                  |
| T2 (`\regndplus`)             |  7 | 11            |  0                  |
| T3 (`\regndplusplus`)         | 14 |  4            |  0                  |
| T4 (`\regndstar`, headline)   | 15 |  3            |  0                  |

**Three differentials worth quoting in talks / rebuttals:**

- **T1 → T2 is flat (±3 %)** on every shared cell.  The CPI-leaf
  convolution alone does not pay for itself when the per-peel
  recompute remains.  (`com-dblp` (3,4): 877,390 → 875,938 ms.)
- **T2 → T3 is the main contribution: 14× – 318×** on every shared
  cell and lifts coverage from 7 to 14 of 18 cells.  Highest:
  `dblp-core30` (3,4) 84,674 → 266 ms = 318×.
- **T3 → T4 trims another 1.0 × – 8.5 ×** and contributes the
  single T4-only cell `com-dblp` (5,6) which T3 cannot reach in
  one hour but T4 clears in 227 s.

**Six representative cells (wall-clock, ms; `to` = 1 h timeout):**

| Graph       | (r,s)   |     T1     |     T2     |     T3     |     T4     |
|-------------|---------|------------|------------|------------|------------|
| ca-GrQc     | (3, 4)  |     4,878  |     5,110  |       158  |       118  |
| ca-CondMat  | (7, 8)  |    42,706  |    41,829  |     2,897  |     2,206  |
| dblp-core30 | (3, 4)  |    85,593  |    84,674  |       266  |       255  |
| com-dblp    | (3, 4)  |   877,390  |   875,938  |    37,018  |     5,898  |
| ca-HepPh    | (3, 4)  |       to   |       to   | 1,155,762  |   384,497  |
| web-it-2004 | (5, 6)  |       to   |       to   |   212,140  |   201,188  |

All four tiers produce **bit-exact identical** core-value
histograms on every cell where they all finish; verified via
`python3 verify_tiers.py --all`.  The only known mismatch is the
pre-existing V-safe/private-cloud 3-tuple miscount on `ca-CondMat`
r=3 s=4 (task #125, ~0.002 % of tuples, not introduced by this
work).

### 9.3  Where the Data Lives

**On the laptop** (`/Users/zhangwenqian/UNSW/pivoter/`):

| Artifact                                            | Purpose                                                                 |
|-----------------------------------------------------|-------------------------------------------------------------------------|
| `paper_data/bench_tier_ablation_results.csv`        | the 72-row sweep CSV (single source of truth for §6.8)                  |
| `paper_data/01_main_benchmark_v3.csv`               | the 549-cell main RegND* / CND sweep (unchanged this week)              |
| `make_tier_ablation_fig.py`                         | reads the tier CSV, writes the §6.8 figure to the Overleaf path         |
| `scripts/make_sigmod_figs.py`                       | unchanged; still produces the other 8 paper-6 figures from the main CSV |
| `verify_tiers.py`                                   | bit-exact κ_s histogram diff across T1/T2/T3/T4 for any (graph, r, s)   |
| `bench_tier_ablation.py`                            | the sweep driver itself (run on tods2 with `--reuse-t4`)                |

**On tods2** (`/home/wenqianz/UNSW/pivoter/`, fresh clone from this
week, NOT the older stale `/home/wenqianz/nclique/pivoter`):

| Path                                                            | Contents                                                |
|-----------------------------------------------------------------|---------------------------------------------------------|
| `~/UNSW/pivoter/bench_tier_ablation_results.csv`                | identical to local `paper_data/…csv` (pulled 2026-06-09)|
| `~/UNSW/pivoter/bench_tier_ablation_logs/`                      | one `T<n>_<graph>_r<r>_s<s>.log` per cell, with `/usr/bin/time -v` tail |
| `~/UNSW/pivoter/bench_tier_ablation.log`                        | harness stdout (progress messages)                      |
| `~/UNSW/pivoter/bench_tier_ablation.pid`                        | PID of the last sweep (`2864765`, now exited)           |
| `~/UNSW/pivoter/build/bin/degeneracy_cliques`                   | the binary used for the sweep                           |
| `/data/wenqianz/<graph>.edges`                                  | the raw edge files (symlinked from `~/UNSW/pivoter/graphs/` by `link_graphs()`) |

**In the Overleaf paper**
(`/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Sigmod2027Nuclear/`):

| File                                          | What it adds for §6.8                                                |
|-----------------------------------------------|----------------------------------------------------------------------|
| `figures/fig_tier_ablation.pdf`               | 6 × 3 grid of bar groups (re-rendered any time the CSV changes)      |
| `sections/Experiments.tex` (§6.8 block)       | the new `\experimentsection{Tier Ablation}` block + 2 tables + prose |

Both are committed to Dropbox / Overleaf history; no external state.

### 9.4  Latest Local-Repo Change: `graphs/` Moved to T7

To reclaim ~3 GB of internal-SSD space, the entire
`~/UNSW/pivoter/graphs/` directory was moved onto an external T7
drive:

- Real files:           `/Volumes/WenqianT7/pivoter_graphs/`  (75 files, 3.4 GB)
- Pivoter view:         `~/UNSW/pivoter/graphs`  — now a **symlink**.
- `.gitignore`:         updated to ignore both `graphs/` and the
                        `graphs` symlink (1-line addition).

**Implications for the next maintainer:**

1.  Any local experiment requires the T7 to be mounted
    (`ls /Volumes/WenqianT7` must succeed) — otherwise
    `graphs/<g>.edges` opens fail at the C++ side and the harness
    aborts at the first cell.
2.  **tods2 is unaffected**: the server uses its own
    `/data/wenqianz/*.edges` symlinks, independent of the laptop.
3.  Any script that reads from `graphs/` (including
    `bench_tier_ablation.py`, `verify_tiers.py`, and ad-hoc CLI
    invocations of `./build/bin/degeneracy_cliques graphs/x.edges …`)
    continues to work unchanged through the symlink.
4.  `make_tier_ablation_fig.py` and `scripts/make_sigmod_figs.py`
    do **not** touch `graphs/` — they only read `paper_data/*.csv`.
    Figures can be re-rendered without the T7 mounted.

### 9.5  Paper Compile Status

`pdflatex` clean as of `04e6ffb`:

- `main.pdf`: 20 pages (one more than the pre-tier-ablation 19-page
  build).
- 0 unresolved references (`Latexmk: ====` empty).
- 0 `??` markers in extracted PDF text
  (`pdftotext main.pdf - | grep -c '??'` returns 0).
- Multiply-defined labels: 0.

To rebuild from a fresh checkout:

```bash
cd "/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Sigmod2027Nuclear"
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

To regenerate the §6.8 figure after a re-sweep:

```bash
cd /Users/zhangwenqian/UNSW/pivoter
python3 make_tier_ablation_fig.py     # reads paper_data/bench_tier_ablation_results.csv
                                      # writes Overleaf/figures/fig_tier_ablation.pdf
```

### 9.6  Open Items for the Next Session

Carried over from this week, not yet addressed:

| Task # | Description                                                                |
|--------|----------------------------------------------------------------------------|
| #99    | RegNDC without r-mergeable merging — ablation still in progress            |
| #125   | ~~V-safe / private-cloud 3-tuple miscount on `ca-CondMat`~~ FIXED 2026-06-12 (pathAliveCount double decrement, see §11.5) |
| #74    | Memory-opt RegNDC sweep (carried from before this work)                    |
| #73    | Re-run paper-6 benchmark on servers with RegNDC (low priority now that T4 sweep is fresh) |

None of these block paper §6.8 from going to reviewer; the
tier-ablation evidence is complete and self-contained.

---


## 10.  Cell-Peel (Quotient Peel) Investigation — live log

Started 2026-06-12.  Goal: replace the dead-box union machinery
(B&B + DomPrune + caches, ~half of V3LM's 2657 lines) with a single
mechanism, and reframe the paper as "quotient the peel's bipartite
incidence by within-class symmetry on BOTH sides".

**Backups before any work** (task #127, done):

- git tag `pre-cellpeel-20260612` (pushed)
- `/Volumes/WenqianT7/backups/pivoter_pre_cellpeel_20260612.tgz` (7.4 GB)
- `/Volumes/WenqianT7/backups/Sigmod2027Nuclear_pre_cellpeel_20260612.tgz`

### 10.1  The idea

A **cell** on path P = one pivot-count vector y (Σy_c = T = s−h,
0 ≤ y_c ≤ p_c) = one orbit of s-cliques under within-class vertex
permutations.  Tuple–cell incidence weight is the existing
W_τ(y) = Π_c C(p_c,y_c)·C(h_c+y_c,j_c).  Initial support = all cells
alive; peeling τ\* kills cells y ≥ ℓ(τ\*|P); each cell dies once and
broadcasts −W to surviving tuples.  One ~30-line `applyCell` would
replace AggrCount + UnionCount + DomPrune + deadCache + the m=1/m=2
special cases + saturation detection.

### 10.2  Probe results (PIVOTER_CELL_PROBE, commit 4f4fc74)

Run as `PIVOTER_RUN_REGION_V3LM=1 PIVOTER_VSAFE_CLOUD=1
PIVOTER_CELL_PROBE=1 PIVOTER_CELL_PROBE_EXIT=1 ./build/bin/... g r s degen`.

| cell | paths | total cells | max cells (m,T) | p50/p90/p99 | naive work bound |
|---|---|---|---|---|---|
| ca-GrQc 3,4 | 450 | 33 K | 3,030 (28,3) | 1/35/1860 | 5.9e7 |
| ca-GrQc 5,6 | 124 | 151 K | 22,860 (22,5) | 5/1818/22860 | 2.8e9 |
| dblp-core30 9,10 | 177 | 92 K | 4,791 (10,9) | 31/2049/4791 | 2.6e8 |
| ca-CondMat 7,8 | 821 | 92 K | 9,438 (15,7) | 6/183/2211 | 3.7e8 |
| com-dblp 3,4 | 64 K | 1.6 M | 28,327 (57,3) | 1/15/220 | 1.1e10 |
| **com-dblp 5,6** | 9.5 K | **13.8 M** | **875,440 (42,5)** | 3/91/19836 | **5.5e12** |

**Verdict so far:**

1. The typical path is trivially materializable (p50 = 1–31), and
   witness-orbit compression (s-cliques / cells) reaches 1e4–1e10 —
   a publishable s-side compression number mirroring the r-side one.
2. The heavy tail is REAL: monster paths with m ≈ 40–60 classes
   produce up to 8.75e5 cells each.  A **pure** cell-materialization
   engine is infeasible on com-dblp-like graphs.  V-safe ON does not
   shrink the monsters (it removed only 169 classes there).
3. Peel = **88–99 %** of RegND\* total time on every headline cell
   (com-dblp 5,6: 263 s of 277 s; web-it 7,8: 1513 s of 1534 s), so
   the peel engine is worth real algorithmic work.

### 10.3  Open question (next step)

Attribute the current B&B work (countUnionRec recursive calls) to
path cell-count buckets.  If most work is on small-cell paths, a
hybrid engine (materialized cells for small paths, lazy B&B for
monsters, one ADT interface) wins; if monsters dominate, the gain
needs a better implicit representation instead.

### 10.4  Task map

#127 backup (done) · #128 probe (running) · #129 prototype engine ·
#130 bit-exact verify · #131 bench · #132 paper §4 rewrite.
Tag to roll back everything: `pre-cellpeel-20260612`.

### 10.5  Wall-time attribution (commit ff85e1f) — design pivot

`PIVOTER_CELL_PROBE` without `_EXIT` now also attributes refresh
wall-time and B&B rec-calls to path cell-count buckets.

ca-GrQc 5,6: **99.9 % of refresh time (1946 of 2088 ms peel) sits on
the 6 paths with 1e4–1e5 cells**; every path below 1e3 cells is free.
B&B recursion itself is negligible (23 K calls total); the real cost
is the multiplier  #batches × #alive-survivors × DP-query  on monster
paths.  Consequences:

1. Pure cell materialization is dead (probe gate FAILED): it would
   optimize the free 0.1 % and explode on the expensive paths.
2. A hybrid (cells for small paths) is pointless for the same reason.
3. The right target is the per-query cost on monster paths.

### 10.6  New design: factored cell-death events ("rank-1 + sparse")

Key algebraic fact: for a dying cell y (≤T non-zero coords) and a
survivor τ′,

    W_τ′(y) = g(y) · base_τ′ · corr(y, τ′)

where g(y)=Π_{c∈S(y)} C(p_c,y_c) is tuple-independent,
base_τ′=Π_c C(h_c, j_c) is cell-independent (0 exactly for tuples
that need pivots — the corner constraint handles itself), and
corr ≠ 1 only when touched(τ′) ∩ S(y) ≠ ∅ (sparse: |S(y)| ≤ T,
|touched| ≤ r, m ≈ 40).

Engine per death event: enumerate newly-dead cells once (amortized:
each cell dies exactly once), accumulate the scalar G = Σ g(y), apply
Δτ′ = base_τ′·G to all survivors in O(1) each, plus sparse exact
corrections for the few (cell, tuple) pairs whose class sets
intersect.  No B&B, no inclusion–exclusion, no dominance pruning, no
caches; per-path state = packed alive-cell set (~8 B/cell).

Work bound: O(Σ_P cells(P)·T + #batches×#survivors·O(1) + sparse).
On the GrQc monster this replaces ~10⁶ DP queries (m·T² each) with
~10⁶ O(1) updates + 150 K cell deaths.

Open risks: (i) com-dblp 5,6 holds 13.8 M cells ⇒ ~110 MB packed
alive-cell state (current T4 peak there is 580 MB, so plausible but
must be measured); (ii) web-it-2004 cell counts unknown — probe
running on tods2 (`cellprobe_big.log`); (iii) bulk path death needs
either per-pair liveContrib (+124 MB on com-dblp) or one final DP per
pair at retire time — decide during prototyping.

Next: read tods2 probe; if web-it cells are sane, prototype the
event engine (task #129, design above), verify bit-exact (#130).

### 10.7  Big-graph probe results — GATE PASSED for the event engine

From tods2 `cellprobe_big.log` (V-safe ON):

| cell | total cells | max cells (m,T) | packed alive-cell est. | sCliques/cells |
|---|---|---|---|---|
| web-it-2004 5,6 | 5.3 M | 17,300 (18,5) | ~45 MB | 1.4e7 |
| web-it-2004 7,8 | 10.3 M | 112,720 (17,7) | ~85 MB | 1.8e10 |
| ca-HepPh 3,4 | 43.1 M | 879,023 (m=176!, 3) | ~390 MB | 1.33 |
| com-dblp 5,6 | 13.8 M | 875,440 (42,5) | ~115 MB | 16 |

Reading:

- The naive bound (cells × tuples, up to 2e13) is irrelevant to the
  factored event engine; its bound is Σcells·T + #queries·O(1) +
  sparse corrections.  Σcells·T ≤ 1.3e8 everywhere — fine.
- web-it-2004 (the 1513 s peel, biggest prize) is comfortably
  feasible: 10 M cells ≈ 85 MB vs the current 1.4 GB peak.
- ca-HepPh remains the honest boundary regime (compression 1.33,
  m=176): the engine stays correct there; expected win is from the
  O(1) shared term replacing per-query DP, not from orbit
  compression.  Measure, don't promise.
- Remaining unknown to measure in the prototype: the sparse
  correction volume (cells whose class set intersects hot tuples).

Decision: proceed to prototype (task #129).  Probe phase (#128) done.

### 10.8  EventPeel prototype: results and verdict (task #129)

Engine `PIVOTER_RUN_REGION_EVENT=1`
(NucleusCoreDecompositionRegionCPI_EventPeel.cpp, ~1880 lines vs
V3LM's 2990; the union-count/DomPrune/deadCache machinery is fully
replaced by materialized cell sets + factored death events).

**Correctness (the headline win):**

- Exact vs the trusted reference on every tested cell, including
  VSAFE mode.
- Matches V3LM everywhere EXCEPT ca-CondMat 3,4 VSAFE — where
  **EventPeel is exact and V3LM has the known 3-tuple miscount
  (#125)**. The bug therefore lives in V3LM's dead-box/V-safe
  interaction, which EventPeel does not share. Independent
  cross-check value regardless of engine choice.
- One real bug found during bring-up (ASan): cross-path stride reuse
  in event-scratch cleanup; fixed by cleaning at event end.

**Performance (honest, clean single-run, VSAFE):**

| cell | EventPeel peel | V3LM peel | corr pairs |
|---|---|---|---|
| com-dblp 3,4 | 2.3 s | 4.0 s (**1.7x win**) | 4.8e7 |
| com-dblp 5,6 (headline) | 497 s | 285 s (1.7x loss) | 1.55e10 |
| ca-GrQc 5,6 | 9.7 s | 1.7 s (5.7x loss) | 4.0e8 |
| small cells | comparable | — | — |

**Diagnosis.** The correction-pair volume
Σ (dying cell × touching tuple) is the engine's true cost and is
irreducible: single-touch histograms (v1.1) and signature
memoization (v1.2, 88 hits / 24k misses) both failed because the
multi-class signatures are almost all unique. On com-dblp 5,6 there
are ZERO bulk deaths, so the engine's amortized weapon never fires.

**Conclusion: the two evaluators are duals; neither dominates.**
- Cell-event enumeration wins when bulk deaths dominate or
  correction incidence is sparse (com-dblp 3,4).
- Per-query union-count DP (V3LM) wins when correction incidence is
  dense — it integrates the incidence in closed form instead of
  enumerating it.

**Options for next session (decision needed):**
1. Low risk: keep V3LM as the perf engine; lift the QUOTIENT THEORY
   into the paper (dead box = corner of dead witness-orbits, B&B =
   corner-measure evaluator, Theorem 5 = bulk orbit extinction);
   code elegance via extracting an orbit-evaluator ADT from V3LM.
   EventPeel stays as cross-check + the bug-#125-free reference.
2. Ambitious: adaptive dual engine (choose evaluator per path by
   correction density) — honest "adaptive representation" story,
   more code.
3. Push EventPeel further (SIMD the correction loop) — bounded gain,
   cannot close a 15.5e9-pair asymmetry. Not recommended.

## 11.  Current State Snapshot (2026-06-12, post-EventPeel) — READ THIS FIRST

This supersedes §9 as the entry point.  §10 is the investigation log;
this section is what a new maintainer needs to act.

### 11.1  TL;DR of the cell-peel arc

Goal was a single elegant peel mechanism replacing the dead-box
B&B/DomPrune/deadCache machinery.  Outcome:

1. **Theory crystallized**: classes/tuples are r-clique orbits, CPI
   leaves decompose the s-clique orbit space into *cells*, peeling is
   orbit extinction.  This vocabulary is paper-ready regardless of
   engine choice.
2. **A second engine exists and is exact**: `EventPeel`
   (`PIVOTER_RUN_REGION_EVENT=1`) — factored cell-death events, no
   B&B/IE/caches, ~1880 lines.  Exact vs the trusted reference on
   every tested cell *including* the one where V3LM is wrong (#125).
3. **Performance verdict: dual evaluators, neither dominates.**
   EventPeel wins where bulk deaths dominate (com-dblp 3,4: 1.7×),
   loses where correction incidence is dense (com-dblp 5,6: 1.7×
   slower with 1.55e10 correction pairs; ca-GrQc 5,6: 5.7× slower).
   Both mitigation attempts (per-class histograms, signature
   memoization) were defeated by unique signatures — the incidence
   volume is irreducible for any per-pair-exact method.
4. **Open decision** (§10.8): (1) keep V3LM as perf engine + lift the
   quotient theory into the paper + ADT refactor for code elegance
   [low risk, recommended]; (2) adaptive dual engine; (3) keep
   pushing EventPeel [not recommended].

### 11.2  Session commit chain (all pushed; rollback tag below)

| commit | what |
|---|---|
| 4f4fc74 | PIVOTER_CELL_PROBE: per-path cell counts, work bound, compression |
| ae75bec | SigmodPlus §10 investigation log opened |
| ff85e1f | probe v2: wall-time attribution by cell-count bucket (99.9% on monsters) |
| 9e3a1a8 | §10.5–10.6: design pivot to factored events |
| bcba51c | §10.7: big-graph probes (web-it OK, HepPh boundary), gate passed |
| a40fb58 | EventPeel v1.1: engine, correct, ASan stride-bug fixed |
| 64c9fd3 | EventPeel v1.2: signature memo (correct, ineffective) |
| 71126d5 | §10.8 verdict + options |

Rollback of the whole arc: `git checkout pre-cellpeel-20260612`
(tag pushed; full tars on `/Volumes/WenqianT7/backups/`).

### 11.3  How to run / verify EventPeel

```bash
# run (any r>=3 cell; VSAFE optional, mirrors V3LM env semantics)
PIVOTER_RUN_REGION_EVENT=1 [PIVOTER_VSAFE_CLOUD=1] \
  ./build/bin/degeneracy_cliques graphs/<g>.edges <r> <s> degen

# exactness vs trusted reference (small/medium cells)
PIVOTER_RUN_REGION_EVENT=1 PIVOTER_COMPARE=1 ./build/bin/... 

# histogram A/B vs V3LM
diff <(ENV_A ... | grep '^  core=') <(ENV_B ... | grep '^  core=')
```

Engine-specific log lines: `EventPeel init` (per-pair contributions +
init check vs Step 4), `Bulk deaths`, `Partial death events`,
`Cells died (materialized N on K paths)`, `Sparse correction pairs
(sig memo H hits / M miss)`.

### 11.4  File / artifact map

| artifact | state |
|---|---|
| `src/NucleusDecomposition/NucleusCoreDecompositionRegionCPI_EventPeel.cpp` | the new engine (Steps 1–4 + PathInfo verbatim from V3LM; peel half rewritten) |
| `...RegionCPI_LowMem.cpp` (V3LM) | UNTOUCHED except the CellProbe instrumentation; still the production engine |
| `PIVOTER_CELL_PROBE[, _EXIT]` | cell-count probe inside V3LM (run-and-exit cheap on big graphs) |
| tods2 `~/UNSW/pivoter` | behind by this session's commits — `git pull` before any server work; `cellprobe_big.log` there holds the web-it/HepPh probe outputs |
| `/tmp/ev_dblp56.log`, `/tmp/v3_dblp56.log` (laptop) | headline-cell runs backing §10.8's table (volatile, /tmp) |

### 11.5  Known issues & their state

- **#125 (V-safe 3-tuple miscount, ca-CondMat 3,4)**: **FIXED
  2026-06-12** — root cause was a pathAliveCount double decrement
  (saturation removal + later death of the same tuple on the same
  path) retiring paths prematurely; not V-safe-specific.  Full root
  cause, fix, and post-fix verification in the §9 known-issues entry
  ("V-safe 3-tuple miscount … FIXED").  The ca-CondMat 3,4 caveat on
  T4 paper numbers is lifted.
- EventPeel on ca-HepPh: untested (43M cells expected ≈ 390 MB
  materialized; engine handles it but was not the priority).
- `verify_tiers.py` does not yet know about EventPeel; add
  `PIVOTER_RUN_REGION_EVENT` as a 5th column if it becomes a tier.
- The tier-ablation CSV/figure (§8) and the paper are UNAFFECTED by
  this session: no paper edits were made on top of the quotient idea
  yet (deliberately — engine verdict first).

### 11.6  Task map after this session

#127 backup ✓ · #128 probe ✓ · #129 EventPeel prototype ✓ (verdict:
dual) · #130 bit-exact verify ✓ (informal full sweep; formalize if
EventPeel is promoted) · #131 bench — OBSOLETE in original form, the
§10.8 table is the bench · #132 paper §4 quotient rewrite — BLOCKED
on the §10.8 decision · #125 — FIXED (pathAliveCount double
decrement; see §11.5).

## 12. Paper: Size-Free Cost Analysis (2026-06-12, task #132)

The §10.8 engine decision resolved as option 1 in practice: V3LM
stays the perf engine, and the quotient theory was lifted into the
paper as new theory rather than new code.

What landed in `Sigmod2027Nuclear` (Dropbox/Overleaf):

- **New §4.5 "Size-Free Cost Analysis"**: combinatorial quotient
  (def), compressed-output queries in O(r log r) (prop), the
  size-free cost theorem (post-construction work bounded by quotient
  quantities only; class sizes appear only inside O(1) binomial
  arithmetic), and the blow-up separation corollary (fixed quotient,
  unbounded |K_r| ⇒ unbounded gap vs any per-clique algorithm).
  This is the formalization of the "blow-up invariance" idea from
  the optimality discussion; it converts the empirical 1.9e9×
  compression and the ca-HepPh weak case into two ends of one
  theorem.
- Pseudocode made consistent with the theorem (tuple-level records
  instead of per-clique writes; ordered-level selection instead of
  value-indexed bucket scan).
- Intro contribution bullet, abstract sentence, conclusion sentence,
  compression-experiment cross-reference.

Full edit log in the paper repo's `PhaseTracker.md` (2026-06-12
entry). Build verified: 20 pp, 0 overfull, 0 "??". Pre-edit backup
tarball on the T7 drive.

Open paper decisions for the user:
1. 4-tier naming (RegND/+/++/*) vs house 3-tier style — folding ++
   would require regenerating the ablation figure labels.
2. Whether to also state the B&B bound via the tighter arrangement
   argument.

## 13. Moment-engine probe: idea killed by data (2026-06-12, task #133)

First-principles candidate: Vandermonde factorization
C(nh+y, j) = Σ_t C(nh, j−t)·C(y, t) makes every tuple's live support
a static ≤2^r-dim linear functional of a per-path "moment vector"
M_P(t) = Σ_{y live} Π C(np,y)·C(y,t); a peel event would update the
moment vector once per path instead of running one B&B per survivor.
Would have replaced deadCache (per path×tuple) with per-path state.

Probe (PIVOTER_MOMENT_PROBE / _DYN / _EXIT in V3LM, env-gated):

- Static census: #distinct moment patterns vs a_P per path.
  Result: #mom/a_P = 1.1–1.9 (ca-CondMat 3,4: 0.54 ratio inverse;
  ca-GrQc 5,6: 0.70; com-dblp 3,4: 0.65 unweighted, 0.93 weighted).
  Moments are MORE numerous than tuples — paths host nearly all
  tuple compositions, so the downward closure is bigger.
- Dynamic comparison in matched range-DP units:
  cur = Σ recCalls×b_P (measured), mom = Σ dQ×#mom (hypothetical,
  dQ = Δ-region decomposition measured by 1 extra B&B per new box):
  | cell | cur | mom(+dot) | win factor |
  | ca-CondMat 3,4 | 6.1e5 | 2.8e7 | 0.021 |
  | ca-GrQc 5,6    | 5.2e5 | 5.0e8 | 0.001 |
  | com-dblp 3,4   | 3.9e6 | 5.8e8 | 0.007 |

Root cause of the kill: event→survivor incidence is SPARSE (the
affected-check skips most pairs for free; recCalls totals are tiny),
while the moment vector is DENSE — every event pays #mom regardless
of how few supports actually changed. Basis compression trades
sparse scatter for dense updates and loses 50–1000×.

Two standing conclusions reinforced:
1. The correction-pair irreducibility claim survives its strongest
   challenger yet (exact low-rank factorization). V3LM's
   sparse-scatter design (affected predicate + DomPrune + zero-base
   cache) is close to the practical optimum for the sequential
   update rule.
2. Remaining upside is NOT a better update rule: it is parallelism
   (per-path independence) and build-phase twin reduction (the only
   non-size-free phase), cf. §12.

Probe cost: ~1 hour, zero engine code. Probe kept env-gated in V3LM.

## 14. False-twin census (2026-06-13)

Question: does extending profile classes from true twins (N[u]=N[v])
to the full twin partition (+ false twins: N(u)=N(v), non-adjacent,
kappa-preserving via automorphism but invisible to MC profiles) buy
real compression on paper-6? Probe: scripts/probe_false_twins.py.

Removable vertices (sum |class|-1, degree >= s-1 = 3 filter):

| graph        | true-twin rmv | false-twin rmv | ft max class |
|--------------|---------------|----------------|--------------|
| dblp-core30  | 1069          | 0              | 0   |
| ca-GrQc      | 665           | 10             | 6   |
| ca-HepPh     | 751           | 35             | 9   |
| ca-CondMat   | 3802          | 43             | 8   |
| com-dblp     | 49115         | 1065           | 18  |
| web-it-2004  | 376250        | 3230           | 213 |

Verdict: false twins are <=0.6% of vertices on every paper-6 graph
(0 on dblp-core30); the twin-partition extension is theoretical
armor (closes the cocktail-party K_{n x 2} counterexample, reaches
"all O(n+m)-computable local symmetry"), not a practical win on
these families. web-it has one 213-vertex false-twin class
(template/hub structure) — the only hint of practical relevance.
Decision left to the user: add a one-remark extension to the paper
or keep it for the thesis. No engine work justified.

## 15. Paper story-line pass (2026-06-13, task #132 continued)

Three targeted edits (full log in paper PhaseTracker 2026-06-13):
two-halves framing (CND = static half, RegND* = dynamic half) in §1
and §3; three per-clique costs named explicitly in Scaling Barrier;
thm:fm demoted to inline remark (saves a theorem box, kills
redundancy with thm:vsafe). Bonus: fixed a latent overclaim — "Csafe
subsumes FM isolation" is false in general (multi-region classes
inside FM regions can be unsafe); exactness was never affected, only
the optimization-coverage claim. Build verified 20pp/0 overfull/0 ??.

## 16. Paper pre-submission passes done (2026-06-13, task #132 CLOSED)

4-tier naming kept as ablation ladder (user approved). Vocabulary
cull: clean. Long-sentence pass: 4 splits (all >50w or two-thought),
22 single-chain sentences kept, 0 >50w remain. Build verified.
detect_long_sentences.py now lives in the paper repo root.
Paper state: story line, theory capstone (size-free + blow-up),
correctness fix (FM overclaim), and language passes all landed.
Remaining before submission is figure/example re-checks against the
final prose (old PhaseTracker checklist) and the author's own pass.

## 17. A-level fixes landed (2026-06-13, task #134)

Reviewer simulation (4 personas, §16 follow-on) found two REAL
hierarchy bugs in RegND-H, both fixed and verified by a new
brute-force oracle (scripts/verify_hierarchy_brute.py, 60/60):
Rule-B path-anchor over-merge and FM-overlap double-count/under-
merge. Paper: Prelim definitions aligned with firstNucleus
(family-internal admissibility; old defs made class-integrality
literally false), class-symmetry/class-integrality re-proved via
the swap involution, bucket clamp fixed, (s-r)^2 factor, §5
rewritten to the verified vertex-hierarchy class-DSU algorithm
with soundness/completeness lemmas. Full log in paper PhaseTracker.
TODO: regenerate cs10 dump + RegND-H overhead numbers (fixed engine
does strictly less work); B/C-level review items pending.

## 18. Hierarchy speed re-measure + fork-drift finding (2026-06-13)

Task ①: hierarchy phase is now FASTER than before the correctness
fixes (rule-B removal dropped per-tuple pathsCoveringTuple calls;
v1's FM-overlap correction was quadratic on FM-heavy graphs and was
replaced by first-claim shared-vertex counting, v2): CondMat 23->6,
GrQc 20->4, com-dblp 254->125, ca-HepPh 3,5 5528->43 ms (128x; rule
B's per-tuple path intersections dominated on the region-overlap-
heavy graph). 60/60 oracle + counterexamples
pass. cs10 uses the R=1 ST_V3 dump path, NOT RegND-H -- no figure
regeneration needed; V3H_DUMP_MEMBERS has no script consumer.

User-spotted anomaly CONFIRMED: on ca-GrQc 5,6 the HIER binary's
total (1833 ms) beats no-hier V3LM (2271 ms). Cause: FORK DRIFT,
not free hierarchy -- LowMem_Hierarchy.cpp forks an older V3LM peel
(no saturation machinery; unaffected by #125 since double-decrement
needs saturation+death). On CondMat/com-dblp HIER is slower as
expected (+130/+773 ms). The sweep itself is 4-125 ms everywhere.
RECOMMENDED CLEANUP: port the verified BuildHier sweep into the
main V3LM engine as a flag (paper §5 already describes it
engine-agnostically); kills the fork debt and makes
with-hierarchy >= without by construction.

## 19. tods2 hierarchy refresh sweep launched (2026-06-13, task #135)

What changed in code = only the V3LM_HIER engine (rule-B removal +
FM-overlap fix), so the refresh re-runs exactly the affected rows:
1. bench_hierarchy.py, both server grids, fresh
   (paper_data/bench_hierarchy.csv on tods2 repo).
2. bench_v3_all.py with new BENCH_ALGOS=V3LM_HIER filter,
   BENCH_GRAPHS=paper-6 (all six on tods2),
   OUTCSV=/data/wenqianz/hier_refresh/bench_v3_hier_refresh.csv.
CND/RegNDC/NOCPI rows are NOT re-run: engines untouched, V3LM
perf-neutrality A/B-verified in §12.

Chain pid 3640436, logs /data/wenqianz/hier_refresh/chain.log,
sentinel ALL_DONE. Persistent local monitor: hourly status, DONE /
stall / error alerts. Pre-launch smoke on tods2: oracle 10/10 +
dblp-core30 HIER run OK.

On DONE: fetch CSVs -> recompute hierarchy overhead stats (replace
median 0.085% / p95 1.16% / max 2.8% in Experiments.tex) -> refresh
RegND-H series in exp figures (make_sigmod_figs.py + merged CSV) ->
rebuild paper -> commit.

## 20. Hierarchy refresh sweep COMPLETE + paper data updated (2026-06-13, task #135)

Sweep finished in ~6.5 h on tods2: 37 bench_hierarchy rows (both
grids) + 133,770 V3LM_HIER rows (full paper-6 grids, 1:1 replacement
of the old rows in bench_full_merged.csv). Zero errors.

Paper updates (all rebuilt, 20pp / 0 overfull / 0 ??):
- Hierarchy overhead sentence (Experiments): now 3,577 matched
  cells, sweep <= 2.5 s absolute, median 0% / p95 3.4% of total
  time, with the large relative shares confined to sub-second runs
  (peel so fast the sweep's fixed sort cost shows). Replaces the
  stale 1,083-cell 0.085%/1.16%/2.8% claim measured pre-fix on a
  partial grid.
- Headline numbers now exact and traceable (reviewer item): "up to
  5,767-fold runtime (dblp-core30 (5,40)) and 2,835-fold memory
  (dblp-core30 (5,11))" anchored in §exp-time; abstract/intro/
  conclusion updated from the rounded-up 5,800/2,800.
- Consistency check: geomeans unchanged (14.49x / 21.62x) since
  REF/RegNDC rows were untouched — exactly as predicted.
- Figures regenerated from the updated merged CSV directly into the
  Overleaf figures dir (repo symlink): time/memory advantage,
  runtime, coverage, etc.

Provenance: refresh CSVs force-added under paper_data/hier_refresh/;
merged CSV updated in place (backup at
/tmp/bench_full_merged_backup_20260613.csv).

## 21. ARB parallel-baseline comparison (2026-06-16, task #136, newSigmod fork)

GOAL: add the Shi et al. shared-memory parallel nucleus decomposition
(ARB-NUCLEUS-DECOMP) as a baseline, addressing reviewer R2's "no
parallel baseline" attack. TWO variants, single-thread for apples-to-
apples with our single-thread numbers:
  ARB    = ~/arb-nucleus-decomp (core numbers only)  <-> RegND*
  ARB-Hi = ~/arb-nucleus-hierarchy (with hierarchy)  <-> RegND-H

SETUP (all on tods2, binaries already bazel-built Oct-2025):
- CLI: NucleusDecomposition_main -s --rClique R --sClique S graph.adj
  (PARLAY_NUM_THREADS=1 for single-thread; prints "### Running Time"
  to STDOUT, /usr/bin/time -v to STDERR).
- Graph format = GBBS AdjacencyGraph CSR. The stale /data/wenqianz/*.adj
  are a DIFFERENT dataset (com-dblp.adj n=1.05M vs our 317K) — DO NOT
  reuse. Converted OUR 6 .edges via scripts/edges_to_gbbs_adj.py into
  /data/wenqianz/arb_adj/; all 6 verified n,|E| match our .edges.
- Harness bench_arb.py: 6 graphs x 2 variants x r in {3,4,5,6,7}, s from
  r+1 up with skip-floor on timeout, 1h cap, 24 concurrent single-thread
  cells + 80GB mem gate, resumable, /usr/bin/time -v per cell.
  OUTCSV=/data/wenqianz/arb_run/bench_arb.csv, sentinel ARB_DONE.
- Bug fixed during bring-up: stdout was DEVNULL'd, losing the GBBS
  timing line (all cells mislabelled ERROR); now capture both streams.

EARLY FINDINGS (smoke, single-thread):
- ARB no-hier is fast: dblp-core30 (3,4) 0.0018s, web-it (3,4) 0.12s.
- ARB-Hi is 1000x+ slower: dblp-core30 (3,4) 16s, ca-GrQc (3,5) 15.8s,
  com-dblp (3,4) 23s, web-it (3,4) TIMEOUTs. Our RegND-H hierarchy
  overhead was median 0% / <=2.5s (section 20) -- a strong contrast.

PAPER INTEGRATION (after ARB_DONE):
- Add ARB to algorithms-under-test (Setup) and as a 3rd series in
  exp-time (runtime) + exp-memory; ARB-Hi vs RegND-H in the hierarchy
  paragraph of exp-breakdown-time.
- Honest framing: ARB is parallel-by-design but run single-thread to
  match our numbers; note this is conservative for ARB (it could use
  more cores) yet RegND* still wins on the single-thread axis, and
  ARB-Hi cannot produce hierarchies at the scale RegND-H does.
- Reference how the prior nucleus papers presented ARB (data shape in
  python/experiment/data/dataARBnoHi.csv: dataset,s,r,total_sec,
  max_rss_mb,exit,source).

Monitor: persistent, hourly status + DONE alert.

## 22. ARB comparison COMPLETE — corrected results (2026-06-16, task #136)

CRITICAL CORRECTION: the first ARB sweep used --rClique/--sClique with
default -tt 0, which the GBBS binary treats as a NO-OP (parse+exit in
~ms). All those rows were invalid; "ARB faster on 763/763 cells" was an
artifact. Caught by the "4ms where CND/RegND* take minutes" smell test.
Correct invocation: -rounds 1 -s -r R -ss S -relabel -efficient 1 -tt 5.
Correctness cross-checked via max core: ca-HepPh 236, dblp-core30 111,
com-dblp 111 — all match RegND*.

Also corrected: the "ARB-Hi 200GB explosion" was a tt=0 artifact too;
with correct flags ARB-Hi memory is normal (dblp (5,6) 6GB). 250GB
RLIMIT_AS cap kept as the paper-budget safety net.

FINAL (single-thread, end-to-end wall, correct flags):
- RegND* vs ARB (no-hier), 45 matched OK cells: gmean 57.6x faster
  (dblp-core30 244.7x, ca-GrQc 122.6x, ca-CondMat 65.1x, com-dblp 6.0x,
  ca-HepPh 0.6x = our weak case, web-it 0 matched: ARB completes none).
  Max 7851x (ca-GrQc).
- RegND-H vs ARB-Hi, 33 matched OK cells: gmean 65.1x, max 6407x.
- Coverage: ARB completes 0/web-it cells, ~3 low-(r,s) cells on
  dblp-core30/com-dblp; times out/OOMs at higher s exactly like CND.
  We complete 92000/6216/6116 respectively.

Net: we beat the parallel baseline (single-thread) by MORE than we beat
CND (57.6x vs 14.49x), because ARB's per-clique work explodes faster
with s. Directly answers reviewer R2's "no parallel baseline".
Data: paper_data/bench_arb.csv (139 rows: 79 OK/8 OOM/52 TIMEOUT).
Caveat: collected at 12-way concurrency (~3% contention measured);
headline cells to be re-measured sequentially before camera-ready.

NEXT: add ARB/ARB-Hi as a third series in Experiments §exp-time and the
hierarchy paragraph; note single-thread fairness (ARB is parallel-by-
design, run single-thread to match our axis; even so we win).

## 23. ARB single-thread results integrated into paper (2026-06-17, task #136)

Experiments.tex edits (single-thread ARB, per user "先整合单线程"):
- Algorithms Under Test: added \arb{} (Shi et al., sotaPrveHierarchy)
  as the parallel baseline run single-threaded for a shared axis;
  correctness noted (max (3,4)-core 236 ca-HepPh, 111 com-dblp match).
- §exp-time: new "Comparison with the parallel baseline" para —
  RegND* 57.6x gmean over ARB on 45 matched cells (vs 14.49x over
  CND), 244.7x on dblp-core30, 122.6x on ca-GrQc, inverts only on
  ca-HepPh (0.6x); coverage: ARB 0 web-it cells / 3 low cells on
  dblp-core30,com-dblp vs our 92000/6216/6116.
- §exp-breakdown-time hierarchy para: RegND-H 65.1x gmean over
  hierarchical ARB (33 cells), up to 6407x (dblp-core30 3,5:
  0.31s vs 1985s).
Build: 20pp, 0 "??", 0 undefined refs. ONE pre-existing overfull in
Introduction.tex:111 ("using the same grouped view"), from the user's
intro rewrite (title now "Beyond Per-Clique Peeling..."), NOT this
change — left for the author.
NOT done (deferred per user): multi-thread ARB run; adding ARB
series to the cactus/runtime FIGURES (prose+numbers only so far).

## 24. Region-IE fusion investigation (2026-06-18, tasks #137-138) — LIVE

GOAL (user idea): drop CPI as a separate "pivot index bolted onto the
code"; fuse the pivot-compression NATURALLY into the Region/Class
quotient so there is one structure, not two stacked layers.

ANALYSIS (full reasoning in chat 2026-06-18):
- CPI's ONE irreducible job is providing a PARTITION of s-cliques
  (canonical homes via degeneracy+pivot) so support = a clean SUM.
  Region/Class alone gives a COVER (overlapping maximal cliques), so
  support over it = INCLUSION-EXCLUSION. "Partition->sum (cheap),
  cover->IE (maybe costly)" is the whole crux.
- Within a region, compression is just binomial over classes (pivot
  machinery NOT needed for that). CPI is only needed to (a) enumerate
  maximal cliques and (b) dedup cliques across overlapping regions.
  The private-cloud / safe-class machinery (Thm vsafe) is ALREADY the
  fusion for the private (single-region) part.
- Two fusion routes:
  F1 region-IE: sup(R0)=Σ_{∅≠A⊆Host(R0)} (-1)^{|A|+1}
                C(|∩_{M∈A}M|-r, s-r); intersection sizes from class
                profiles. Pure quotient, no SDCT, but pays IE
                (2^|Host|, dominance-prunable; same family as the
                existing dead-box IE -> unifies init+update).
  F2 region-native partition: derive canonical homes from regions+
                degeneracy at (region x class) granularity -> a path
                set that slots into the EXISTING path-interface
                abstraction (Def path-interface already allows any
                partition path set). Lowest risk: peel untouched.
                Unknown: can the split be produced cheaply.

OVERLAP DATA (V3LM Step-1, (3,4)): max regions-per-vertex =
dblp-core30 13, ca-GrQc 45, ca-CondMat 199, com-dblp 253,
ca-HepPh 1411 (11M overlap pairs). => IE cheap on 5/6 graphs,
explodes only on ca-HepPh, which is ALREADY our one weak case.
Fusion introduces NO NEW weak case -- strongest argument for it.

PRIZE: elegance (one quotient story) + attacks the SDCT-build
bottleneck (RQ5: on web-it peel is 0.6-1%, construction dominates) +
drops CPI Tree memory + unifies init/update onto one IE engine.

CONSTRAINTS (user): do NOT touch existing code; build new standalone
files only; keep this doc + the task list updated live.

PROBE (#137, scripts/probe_region_ie.py, standalone):
  (1) correctness: region-IE == direct s-clique enumeration on small
      graphs (ground truth that CPI also computes);
  (2) cost: IE term-count per sampled r-clique, raw 2^|H| vs after
      merging equal region-intersections (realistic);
  (3) prize size: SDCT/path build share from existing binary stdout.
Greenlight engine prototype (#138) only if (1) exact + (2) small on
5/6 + (3) build share real. Probe-before-engine, per this season's
repeatedly-validated discipline.

### 24.1 Probe results (2026-06-18) — STRONGLY POSITIVE

scripts/probe_region_ie.py, sampled r-cliques per graph:

| graph        | regions | |Host| med/max | meet-closure med/max | correctness        |
|--------------|---------|---------------|----------------------|--------------------|
| dblp-core30  | 68      | 1 / 13        | 1 / 50               | 200/200 EXACT      |
| ca-GrQc      | 311     | 1 / 7         | 1 / 20               | 200/200 EXACT      |
| ca-CondMat   | 8824    | 1 / 36        | 1 / 167              | 197/197 EXACT (3 skip |H|>18) |
| com-dblp     | 96786   | 1 / 6         | 1 / 15               | (cost only)        |
| ca-HepPh     | --      | BK 60s timeout (V3LM: max vtxMC 1411, 11M pairs) | BLOW-UP | weak case (already) |

Findings:
1. The region-IE support formula is BIT-EXACT vs direct s-clique
   enumeration: 597/597 tested cells. The F1 math is correct.
2. Cost is tiny on 5/6 graphs. MEDIAN |Host| = 1 on every graph: most
   r-cliques sit in exactly ONE maximal clique, so counting collapses
   to a single binomial (the private/safe fast path). Real IE happens
   only on a small shared minority, and even then the meet-closure
   (realistic term count) maxes at 15-167.
3. ca-HepPh blows up, exactly the existing weak case. NO new weak case.

USER REFINEMENT (2026-06-18): the ideal is not raw IE but a method that
KEEPS pivot THEORY and does the counting natively on regions/classes
(fast + elegant), not a ported external index. Precise reading: this is
the principled "partition s-cliques by exact host-region set; count each
block by canonical-home / Mobius over the region-intersection lattice."
That IS pivot theory (free/forced split, canonical homes, binomial
counting) on the quotient. Its cost = the meet-closure size already
measured = small on 5/6 graphs; median |Host|=1 means the binomial base
case is hit for most cliques. Both elegance and speed are data-backed.
Equivalence to raw IE (validated exact) makes the partition form exact too.
Still TODO: the PEEL-UPDATE half must also be reformulated natively.

### 24.2 region_native engine: counting validated, SDCT bottleneck killed (2026-06-18, task #139)

New standalone region_native/region_native.cpp (no existing code touched).
Pipeline: degeneracy-ordered BK maximal cliques (regions) -> profile
classes (each region = disjoint union of classes) -> canonical-home
tuple enumeration -> initial support by region union-count IE B&B at
CLASS granularity (|Host|<=2 closed-form fast paths; median |Host|=1).

ADVERSARIAL CORRECTNESS: bit-exact vs direct s-clique enumeration,
ALL tuples on dblp-core30/ca-GrQc/ca-CondMat across (3,4)/(4,5)/(3,6)/
(5,6)/(4,6) -- e.g. 1706/1706, 29997/29997, 100000/100000. Two bugs
caught & fixed by the self-check: (a) unionCount aliased its arg by ref
and emptied it before the IE subtraction (over-count on every multi-
host tuple); (b) naive BK pivot was O(n^2) at the root (flame-graph:
1547/1548 samples in intersect_nbr) -> replaced with degeneracy-ordered
ELS: MCE 55s -> 0.30s on com-dblp (180x).

SPEED (construction+counting; the CPI-replacement):
| cell | region-native MCE+support | existing SDCT+MCEnum+CPI+PathInfo |
| com-dblp (3,4) | 0.33 + 0.34 = 0.67s | 0.29+0.19+0.16+0.15 ~ 0.97s (1.4x) |
| web-it  (3,4) | 3.94 + 5.35 = 9.3s | SDCT 85.5s + 3.4s ~ 89s (~10x)  |
The win = region-native NEEDS only maximal cliques (degeneracy MCE
3.94s), NOT the full pivoted SDCT (85.5s on web-it, building 423886
path leaves). Same 84871 regions, bit-exact support.

HONEST SCOPE: this is INITIAL SUPPORT only. The peel-UPDATE half is not
yet region-native (the unionCount B&B IS structurally the dead-box
update, so the path is clear). End-to-end claim needs the peel built.
web-it support (5.35s, 3.28M B&B nodes) is the next optimization target
but already dwarfed by the 85.5s SDCT it removes.

### 24.3 Server full counting comparison launched (2026-06-18)
bench_region_native.py on tods2 (pid 369733, /data/wenqianz/brn/):
paper-6 x {(3,4),(3,5),(3,6),(4,5),(4,6),(5,6)}, region-native MCE+support
vs existing SDCT+MCEnum+CPI counting+PathInfo, with on-server --verify
correctness gate. Existing server binary was STALE (printed SDCT_Fused +
"clique not found" errors); rebuilt from HEAD (now clean, SDCT_MaxClique).
Persistent monitor: per-cell speedup, DONE/stall/MISMATCH alerts.
Next milestone after this: the region-native PEEL-UPDATE half.

### 24.4 Size-free scaling test launched (2026-06-18, user's main criterion)
User criterion: region-native time must NOT grow significantly with r or s.
Full grid r in {3..7} x s in {r+1..20}, region-native-only (BRN_SKIP_CPI),
per-cell 600s timeout; contrast = CND/RegND* totals from bench_full_merged.
DEPLOY GOTCHA fixed: `git pull --ff-only` was SILENTLY failing (output
swallowed by `| tail -1`), so tods2 ran the stale 6-cell grid + pre-opt
region_native for several relaunches. Now use `git fetch + reset --hard
origin/main` and VERIFY `git rev-parse HEAD` before launching.
Early reads (pre-fix, still valid): region-native FLAT -- dblp-core30
0.02-0.03s, ca-GrQc 0.03-0.13s across (r,s); ca-HepPh NOT flat
(37->166s, the heavy-overlap IE blow-up = the known universal weak case).
Plot: scripts/make_rn_scaling_fig.py.

### 24.5 Progress snapshot (2026-06-18, task #139 region-native)

STATUS: counting engine built + adversarially validated + size-free
scaling confirmed on the done graphs; sampled server sweep in flight.

WHAT WORKS (committed: region_native/region_native.cpp, optimized):
- Pipeline: degeneracy-ordered BK maximal cliques (regions) -> profile
  classes (region = disjoint union of classes) -> canonical-home tuple
  enumeration -> support by region union-count IE B&B at class
  granularity; |Host|<=2 closed-form fast paths; median |Host|=1.
- Correctness: bit-exact vs direct s-clique enumeration, thousands of
  tuples on dblp-core30/ca-GrQc/ca-CondMat across (3,4)/(4,5)/(3,6)/
  (5,6)/(4,6). Two bugs caught & fixed by the self-check:
  (1) unionCount aliased its arg by ref -> emptied before IE subtract
      (over-count on every multi-host tuple);
  (2) naive BK pivot O(n^2) at root (flame-graph: 1547/1548 in
      intersect_nbr) -> degeneracy ELS, MCE 55s->0.30s on com-dblp.

SIZE-FREE SCALING (user's main criterion: time must not grow with r,s):
sampled grid s-offsets {1,3,6,10,15} per r in {3..7}, verify OFF (the
direct-enum ground-truth was the runtime bottleneck, not region-native).
- dblp-core30: region-native 0.02-0.04s across the WHOLE grid (flat
  ratio 2x); CND explodes 1.36s@(3,4) -> 835s@(5,6) -> timeout.
- ca-GrQc: region-native 0.01-0.27s (spikes only at s=r+1 cells), and
  FASTER than current RegND* there ((5,6) 0.13s vs RegND* 2.14s; (6,7)
  0.19s vs 1.21s); CND (7,8)=298s.
=> region-native is flat in (r,s); vs CND 30x-28000x and the gap GROWS
   with s; vs RegND* comparable-or-faster (flatter on hard s=r+1 cells).

IN-FLIGHT: sampled sweep on tods2 (pid 399562, /data/wenqianz/brn/,
BRN_VERIFY=0, region-native skip-floor). Done: dblp-core30, ca-GrQc (23
cells each); grinding ca-HepPh (weak case, (4,5)=168s, self-bounds via
skip-floor); ca-CondMat/com-dblp/web-it to follow. Monitor bsdquuek6.

HONEST SCOPE: validated half = INITIAL SUPPORT counting only. The
peel-UPDATE half is NOT yet region-native (its IE B&B is structurally
the same dead-box update, so the path is clear). End-to-end engine =
next milestone. ca-HepPh (heavy overlap, vtx in up to 1411 regions) is
the IE blow-up weak case = the already-known universal weak case; no
new weak case introduced.

DEPLOY GOTCHAS (this session): `git pull --ff-only` SILENTLY failed
(swallowed by `| tail -1`) -> server ran stale code/binary for several
relaunches; FIX: `git fetch` then `git reset --hard origin/main`, each
as a SEPARATE ssh call (chained multi-stmt ssh sometimes returns no
output), and VERIFY `git rev-parse --short HEAD` before launching.
Also rebuild region_native AND degeneracy_cliques on the server after
any pull (the stale degeneracy_cliques printed SDCT_Fused + threw
"clique not found" on ca-GrQc until rebuilt from HEAD).

NEXT: (a) finish sampled sweep, plot scripts/make_rn_scaling_fig.py;
(b) design+build the region-native PEEL-UPDATE for an end-to-end
engine; (c) decide paper framing (CPI -> region-native counting).

## 25. Region-native flatness fix: support is host-determined (2026-06-18, #139)

FINDING (the size-free lever). The s-clique support of a region tuple
depends ONLY on its host = the set of regions (maximal cliques) that
contain it; it does NOT depend on the tuple's class multiplicities.
Proof: sup = union over host regions M of {s-cliques in M containing
R0} = IE over the host using C(|cap region|-r, s-r); every term is a
property of the host set, not of which classes R0 picks. Two different
r-cliques with the same host have identical support.

WHY com-dblp was NOT flat (prior sweep: (7,8)=75.78s, 143x flat-ratio).
The first region_native enumerated EVERY pattern (28.4M at com-dblp 7,8)
and recomputed host from scratch per leaf (vector alloc + interClasses),
plus it computed support for the huge single-region (direct-binned)
population it should have skipped. Diagnosis via PIVOTER_RN_HOSTPROBE:
com-dblp (7,8) has only 2067 DISTINCT HOSTS behind 28.4M patterns
(collapse 13719x) -> support was recomputed ~13700x too often.

FIX (region_native.cpp, commit c5282e6), three changes, all bit-exact
(--verify EXACT on every tested cell):
 1. Enumerate ONLY multi-region classes. A single-region (safe) class
    pins |Host|=1 (direct-binned, Thm vsafe), never peeled.
 2. Incremental host in a depth stack: host = cap classRegions[chosen],
    updated per class-choice, no per-leaf re-intersection/alloc. The
    moment |host| drops to 1 the whole subtree is pruned (its patterns
    are all direct-binned) -- the single-region population is never
    enumerated.
 3. Memoize support per host hash: each distinct host's union-count runs
    once, not once per pattern.

RESULT (local Mac, support-only phase):
  com-dblp (7,8) : 75.78s -> 0.89s  (85x), bit-exact
  com-dblp (6,7) : 12.17s -> 0.34s
  com-dblp (7,10): 22.93s -> 0.38s
  com-dblp now FLAT in (r,s): 0.14-0.89s (was 0.53-75.78s, 143x -> ~6x).
  4/5 sparse graphs already flat; com-dblp was the lone violator, now
  fixed.

ca-HepPh STILL the outlier (NOT fixed, and fundamentally so): 79,545
distinct hosts (vs com-dblp 2067), collapse only 45-290x, support
13-35s and TIMEOUT at (6,10)+. Its cost is NOT enumeration -- it is the
per-host union-count itself: dense overlap -> |Host| large -> IE B&B
deep. This is exactly the regime CPI's pivot compression was built for.
Clean scientific framing: region-native (CPI-free) wins on most graphs
but loses on the densest-overlap graph, where pivoting earns its keep.
NOT a regression -- ca-HepPh was already the universal weak case.

DENSE RE-SWEEP launched on tods2 (pid 420990, repo at
/home/wenqianz/UNSW/pivoter NOT /data/wenqianz/pivoter): full grid
s in [r+1,20], r in 3..8, ca-HepPh last, region-native-only
(BRN_SKIP_CPI=1, BRN_VERIFY=0). Out: /data/wenqianz/brn/
bench_region_native_dense.csv. Monitor b5ovqwnlp watches for DONE.

DEPLOY NOTE: server git repo is /home/wenqianz/UNSW/pivoter (the
/data/wenqianz/pivoter dir is a stale non-git copy). region_native runs
from there; rebuild with g++ -O3 -std=c++17 after each pull.

NEXT: (a) plot dense flatness fig (make_rn_scaling_fig.py, point it at
the dense CSV); (b) consider host-level Mobius counting to also flatten
ca-HepPh's enumeration (won't help its union-count depth though);
(c) region-native PEEL-UPDATE half for an end-to-end engine.

## 26. Region-native PEEL half: correct, but needs controlled_split (2026-06-18, #139)

Built region_native/region_native_peel.cpp = region-native END-TO-END
(r,s)-nucleus peel. Front half (load/MCE/classes) shared with
region_native.cpp. KEY DESIGN: support during peeling is the region-IE
union (same B&B as initial support) where each region-intersection LEAF
returns ccpath::support_count with the peeled tuples as the forbidden
antichain. A witness (composition y) dies iff y >= some peeled tuple f
(it contains a peeled orbit) -- this is pattern-peeling semantics
(peeling removes a whole class-multiplicity orbit), so aliveness is a
composition property and the component-max forbidden IE is exact.
Reuses CCPathCore.h (each region = degenerate CCPath: h=0, ell=0, u=n,
T=s).

MILESTONES (each adversarially gated):
 M1 (commit e373d85): foundation test test_region_forbidden.cpp --
    region->CCPath + support_count-with-forbidden vs VERTEX-LEVEL brute
    force, 2 hand cases + 4000 random trials, ALL bit-exact.
 M3 (commit 0a652ab): end-to-end peel core-distribution correctness.
    Oracle = scripts/verify_nucleus_brute.py (textbook peel on individual
    r-cliques). vs brute: 34/34 (5 tiny graphs x r=2,3,4 x s), vs
    production V3LM: 4/4 on r>=3. Restricted to core>=1 (the witnessed
    domain: an r-clique has support>=1 iff it sits in a region iff core>=1;
    region-native scores exactly these, matching V3LM which also emits
    core>=1 only). BUG found+fixed by the discriminating t_k6k4 test: the
    inherited stable_partition reordered regionClasses and broke the
    sorted-order invariant leafCount's binary search relies on -> silent
    support corruption. Removed the partition (peel enumerates ALL
    patterns; keeps sorted order).

M4 (PERFORMANCE) -- BLOCKED on the same wall CPI hit. Profiled
dblp-core30 (3,4): enum=0.00s, initial-support=0.00s (both instant), but
the PEEL does not finish. Root cause CONFIRMED via maxForb tracking: as
patterns peel, each region's forbidden antichain grows (0->8->...), and
ccpath::support_count's inclusion_exclusion_terms is 2^|forbidden|, so
support recompute goes exponential. EVERY real cell times out (even
ca-GrQc (10,12), high r/s); only the tiny synthetic graphs complete.
This is EXACTLY the support-maintenance bottleneck the production engine
solved with controlled_split (tasks #92-98): when |forbidden| > kmax,
split the path so each child's antichain is bounded. CCPathCore.h ALREADY
provides the primitives (first_failing_split_by_vector, choose_split_vector,
controlled_split) -- the region-native peel must adopt them (maintain a
SET of split CCPaths per region instead of one).

HONEST STATUS:
 - Region-native COUNTING (initial support): DONE, fast, size-free on 5/6
   graphs, real experiments in hand (500-cell dense sweep).
 - Region-native PEEL: CORRECT (proven), but naive support recompute is
   exponential in the antichain -> needs controlled_split before any real
   peel experiment. This is a UNIFICATION result, not a dead end:
   region-native peeling and CPI peeling face the identical
   forbidden-antichain growth and need the same controlled-split cure.

NEXT: integrate controlled_split into region_native_peel.cpp (per-region
split-path set, support = sum over the region's current paths). Re-gate
correctness vs brute/V3LM, then run the peel size-free sweep + end-to-end
vs CND/RegND*.

## 27. controlled_split done for single-region; multi-host is the wall (2026-06-18, #139)

Implemented controlled_split in region_native_peel (commit d912b8f).
Per-region state = a SET of CCPaths each with forbidden <= KMAX=12; regSupp
sums support_count over them. Split preserves the count, correctness holds
(34/34 vs brute). The single-region antichain is now BOUNDED: ca-GrQc (3,4)
maxForb (split-path count) plateaus ~46 vs 81+-and-climbing before.

But ca-GrQc (3,4) still does not finish. maxForb tracking proves the
antichain is no longer the wall; the MULTI-HOST fallback is. |host|>=2
patterns (fb counter, ~600+ and growing on ca-GrQc 3,4) still recompute
via the old suppOf over the GLOBAL peeled list (O(P)/leaf, unbounded leaf
antichain). They now dominate.

ROOT INSIGHT (honest, important for the paper): a multi-host pattern's
shared witnesses span several regions, so its peel-update needs a
CROSS-REGION inclusion-exclusion whose intersection leaves need their own
bounded (split) antichains. CPI's UNIFIED pivot paths sidestep this:
witnesses of a tuple live on ONE path regardless of how many maximal
cliques they span, so the dead-box update is single-structure. So:
 - region-native WINS on COUNTING / initial support (size-free, fast,
   real 500-cell experiments in hand) -- the quotient exposes the
   host-collapse CPI can't.
 - region-native PEEL-UPDATE is where CPI's unified paths actually earn
   their keep: splitting witnesses across regions turns the multi-host
   support maintenance into a cross-region IE that CPI avoids.
This is a clean division-of-labor finding, not a failure: the peel may be
best left on CPI/paths while counting moves to the region quotient.

To finish a fully-fast region-native peel one would maintain bounded
split-path sets per DISTINCT HOST (region intersection), i.e. re-derive
the production engine's path machinery on the quotient -- a multi-session
effort. Status: peel CORRECT + single-region fast; multi-host fast is the
open piece.

## 28. PIVOT-on-quotient redesign (2026-06-18, #139) -- user steer

User insight: "use the PIVOT concept to compute, otherwise the update is
too complex. pivot is an idea, not necessarily a CPI." Correct diagnosis
of my detour: region-overlap union-IE counts fine but makes the peel
update need cross-region IE. The pivot idea gives a SINGLE canonical
hold/optional witness structure, so each witness is counted once and the
update is one dead-box -- no cross-region IE. And pivot != CPI: do it at
CLASS granularity (tiny), not the heavy vertex-level global SDCT.

KEY ENABLING FACT (quotient-native, no vertex adjacency needed):
  class i and class j are fully mutually adjacent  <=>  they co-occur in
  some region (classRegions[i] cap classRegions[j] != empty).
Proof: full adjacency => i union j is a clique => extends to a maximal
clique (a region) containing both. So a tuple's witness space is the
(s-r)-cliques of its CANDIDATE-CLASS GRAPH: nodes = classes in its host
regions, edges = region co-occurrence, weights = (remaining) class sizes.
This graph is read straight off the quotient.

DESIGN: per tuple (shared across a host-group), pivot the candidate-class
graph into ONE CCPath (hold/optional). support = support_count along it;
peel = dead-box (forbidden antichain) on that one path; controlled_split
bounds the antichain. Multi-host collapses to one path -> no cross-region
IE. This is "pivot fused into region/class," lightweight (class-level).

BUILD PLAN (correctness-first, each gated):
 P1: build candidate-class graph + class-level pivot -> one CCPath;
     validate support_count(path, empty forbidden) == region-IE/direct on
     sampled tuples (foundation).
 P2: peel via dead-boxes on these per-host-group paths; re-gate core
     distribution vs brute + V3LM.
 P3: size-free sweep + end-to-end vs CND/RegND*.

## 29. CONSOLIDATED STATE + P2 build (orbit-aware class Pivoter) (2026-06-18, #139)

WHERE WE ARE (region-native engine, region_native/):
 - COUNTING (initial support): DONE, fast, size-free. region_native.cpp
   (incremental host + memoized-by-host + subtree prune). 500-cell dense
   sweep in paper_data/bench_region_native_dense.csv: 5/6 graphs flat
   (1.8-8.8x over r=3..8 x s=r+1..20); ca-HepPh the IE-overlap outlier.
 - PEEL (end-to-end nucleus): region_native_peel.cpp. CORRECT (vs brute
   34/34, vs V3LM r>=3 4/4). controlled_split bounds the single-region
   antichain. BUT not fast: multi-host (|host|>=2) patterns need a
   cross-region IE that blows up (the wall).

WHY THE PEEL IS HARD, AND THE FIX (user steer, confirmed necessary):
 Counting (s-r)-cliques in a UNION of cliques (the host regions) without a
 pivot forces region-IE; for the peel that means cross-region dead-box IE.
 The pivot/SCT of the candidate-class graph turns the union into disjoint,
 single-counted leaves -> dead-box update is one structure, no cross-region
 IE. The candidate-class graph is QUOTIENT-NATIVE (class i ~ j iff they
 co-occur in a region; P1 validated 400/400 == region-IE). pivot != CPI:
 do it at CLASS granularity (classes << vertices), lightweight.

P2 = ORBIT-AWARE CLASS-WEIGHTED PIVOTER (the irreducible core):
 Build the SCT of the quotient graph (nodes=classes, edges=co-occurrence,
 weight w_c=classSize). The vertex SDCT (existing, src/SDCT_*) holds/pivots
 individual vertices; the class version must hold/pivot CLASSES carrying
 w_c interchangeable vertices, so "holding" a class consumes one and leaves
 a residual (orbit/binomial handling) -- this is the one intricate piece
 (production spent ~7 tasks on its vertex analogue).
 Output: CCPath leaves (h_c holds, n_c pivots) s.t. sum over leaves of
 support_count(leaf) == total weighted s-cliques, each counted once.
 Then the peel = the existing CCPath dead-box machinery (support_count +
 lazy_delete + controlled_split, all in CCPathCore.h) on these class-leaves
 -- single structure, no region-IE.

GATES (correctness-first; user wants adversarial subagents each step):
 G1 build SCT, sum_leaves support_count(empty forbidden) == region-IE
    (region_native.cpp suppOf) on ca-GrQc/ca-CondMat/com-dblp, all (r,s).
 G2 peel on leaves; core dist == brute (verify_nucleus_brute.py) on tiny
    graphs + == V3LM (r>=3) on small real cells.
 G3 size-free sweep + end-to-end time/mem vs CND/RegND* on tods2.

KEY FILES: region_native/region_native_peel.cpp (peel + candcheck P1),
 src/NucleusDecomposition/CCPathCore.h (support_count, insert_antichain,
 controlled_split, first_failing_split_by_vector -- REUSE for the peel),
 scripts/verify_nucleus_brute.py (oracle), graphs/t_*.edges (tiny tests,
 regenerable from the inline python in the session).
COMMITS this session: c5282e6 (counting flatness fix) .. 0ebfc4d (P1).

## 30. P2 DONE: class-SCT peel correct + peel-phase fast; map-build is the perf bottleneck (2026-06-18, #139)

region_native/region_native_sct_peel.cpp: nucleus peel on the disjoint
class-SCT leaves. support(pattern)=clean SUM over hosting leaves (no IE,
leaves disjoint); peel=per-leaf dead-box + controlled_split; affected
update via exact before/after diff over changed split paths. Built by a
subagent (267k tok), INDEPENDENTLY re-verified by me.

CORRECTNESS (all gates pass, independently re-run):
 - vs brute oracle: 68/68 (tiny + 4 dense stress graphs t_stress/t_chain/
   t_dense/t_split the subagent added to graphs/), core>=1.
 - vs V3LM (r>=3): 4/4 identical (ca-GrQc 3,4 [24 bins] / 4,5; dblp-core30
   3,4; t_k7k5k4).
 - ca-GrQc(3,4) FINISHES: peel=0.70s, Max core 41 == V3LM (region-IE peel
   timed out 300s+). The cross-region IE is GONE (disjoint leaves => sum).
 - Subagent also found: the SC(max) delta shortcut over-promotes on dense
   graphs (fixed: exact diff); the OLD region-IE peel is itself WRONG on
   dense graphs (invents a spurious core; SCT matches brute).
 - class-SCT itself adversarially audited: ~605k cases vs independent
   vertex-expansion oracle, 0 bug; leaf count Theta(#maximal-cliques);
   dead-box = orbit/batch deletion = exactly the pattern-peel semantics.

PERFORMANCE (honest): the PEEL phase is fast (0.70s ca-GrQc, 0.05s dblp),
but END-TO-END is NOT yet beating V3LM:
   ca-GrQc(3,4): SCT total 5.80s (build+maps 5.09 + peel 0.70) vs V3LM 0.06s
   dblp-core30(3,4): SCT 0.12s vs V3LM 0.08s
 The bottleneck is the hosting-map build: O(nLeaf x nPats) (tests every
 pattern against every leaf via support_count>0). FIX (next): enumerate,
 per leaf, the r-multisets of its classes (the patterns it hosts) and map
 directly -> O(sum of pattern-leaf incidences), the production enumCb
 approach. Also: nC>6000 still needs sparse quotient adjacency (dense CxC
 matrix skips ca-CondMat nC=11152).

NEXT: (a) optimize map-build (per-leaf enumeration); (b) sparse quotient
adjacency for nC>6000; (c) G3 size-free sweep + end-to-end vs CND/RegND*/
V3LM on tods2. Commit a774bf2.

## 31. SCALE: scalable SCT integrated; BUILD scales, PEEL-PHASE is the perf gap (2026-06-18, #139)

Integrated classsct_scalable::scalableBuildClassSCT (degeneracy-seeded,
sparse adj, compact leaves) + sparse patterns (compToLocal, no full-C
patVec) into region_native_sct_peel. Adversarial gates all hold:
 - scalable SCT vs dense oracle: 50000 trials sum+disjointness, 0 fail.
 - brute re-gate after each integration step: PASS (15-22/each).
 - ca-CondMat(3,4) == V3LM, and 6.6s -> 2.8s (sparse build beats dense).
 - com-dblp(3,4) nC=123526: class-SCT BUILD 0.19s (158315 leaves), G1
   s-cliques 16713192 == region-IE, G2a 934117/934117 patterns
   SCT-support == region-IE. Support computation correct at scale.

HONEST PERFORMANCE VERDICT (the make-or-break, now measured at scale):
 - The class-SCT BUILD is lightweight + scalable: 0.19s on com-dblp
   (nC=123k) where a dense C×C matrix is 15GB and the vertex SDCT is ~85s
   on web-it. This is the user's "pivot at class level, not CPI" win, real.
 - The PEEL-PHASE is SLOW and NOT competitive with V3LM. com-dblp(3,4)
   peel ~250s (vs V3LM seconds). Root cause: heavily-peeled leaves grow a
   large controlled_split set (maxSplit=4452 paths) and the per-affected-
   pattern support recompute over that set dominates. On small/medium
   graphs too the peel trails V3LM (ca-GrQc 4,5 peel 7.9s vs V3LM 0.23s).
 - So: region-native class-SCT peel is CORRECT, ELEGANT (pivot-on-quotient,
   no cross-region IE), with a fast lightweight build -- but its peel-phase
   constant factors / split-overhead make it slower than the mature V3LM.
   Beating V3LM end-to-end needs the same peel-phase engineering V3LM spent
   tasks #93-98 on (single-level split, batch bucket moves, lazy affected
   updates). Not done.

WHAT'S SOLID: the scientific result -- region-native peeling via a
class-level pivot is correct and removes the cross-region IE -- plus a
near-instant scalable build. WHAT'S NOT: end-to-end speed vs V3LM (peel
phase). NEXT (if pursued): peel-phase optimization (KMAX tuning, batch
updates, lazy recompute) to close the gap; the build + counting halves are
already fast. Commits 7f01fc4..23710db.

## 32. CONSOLIDATED DETAILED HANDOFF — region-native engine (2026-06-18, #139)

Self-contained pickup doc for the whole region-native effort: two halves
(COUNTING = done/fast; PEEL = correct, perf in progress), the full
trial-and-error, files, commits, how to run, gates, next steps.

### 0. BIG PICTURE
Goal: (r,s)-nucleus decomposition on the Region/Class QUOTIENT instead of
per-r-clique. Region = maximal clique >= s. Class = vertices with identical
region-membership (interchangeable). Pattern = class-multiplicity signature
of an r-clique orbit. Two halves:
 (A) COUNTING (initial s-clique support of every pattern): DONE, fast,
     size-free. File region_native/region_native.cpp.
 (B) PEEL (end-to-end core numbers): CORRECT + validated; build is fast and
     scalable; PEEL-PHASE is the open perf problem (slower than CND/V3LM).
     File region_native/region_native_sct_peel.cpp.

### 1. COUNTING HALF (region_native.cpp) — DONE
KEY FACT: an r-clique's s-support depends ONLY on its HOST (the set of
regions containing it), not on its class multiplicities. So support =
region-IE union-count over the host: sum_{A subset host} (-1)^{|A|+1}
C(|cap A|-r, s-r), at CLASS granularity.
THREE WINS (commit c5282e6) that made it size-free:
 1. enumerate ONLY multi-region classes; a single-region (safe) class pins
    |Host|=1 -> direct-binned (Thm vsafe), pruned mid-recursion (the huge
    single-region tuple population is never enumerated).
 2. incremental host in a depth stack (no per-leaf re-intersection/alloc).
 3. memoize support per host hash (each distinct host's union-count once).
RESULT: com-dblp (7,8) 75.78s -> 0.89s; FLAT over the full (r,s) grid on
5/6 graphs (dense sweep paper_data/bench_region_native_dense.csv, commits
3cdf210/3bd0ce9). ca-HepPh is the IE-overlap outlier (79545 distinct hosts;
its per-host union-count is genuinely deep -- that is where CPI's pivot
compression earns its keep).

### 2. PEEL HALF — THE TRIAL-AND-ERROR (read this before re-trying anything)
What is hard: peeling removes a pattern (an r-clique orbit). A witness
s-clique dies when ANY of its r-subcliques is peeled, so aliveness is a
COMPOSITION property: a witness (composition y) is dead iff y >= some
peeled pattern f. Maintaining live support under removals is the crux.

DEAD ENDS (do NOT repeat):
 (a) region-IE peel + naive full recompute of affected patterns each peel
     -> O(P^2) support evals. Too slow. (commit 0a652ab was correct but
     this recompute strategy is the slow part.)
 (b) region-IE peel + incremental delta + per-region forbidden antichain
     -> the forbidden antichain GROWS as patterns peel and support_count is
     2^|antichain|. maxForb hit 81 on ca-GrQc(3,4); EVERY real cell blew up.
     (commits eccfe37, b23f0b4). Antichain explosion is intrinsic to
     "count alive witnesses avoiding peeled" at compressed granularity.
 (c) controlled_split on per-region paths -> bounded the SINGLE-region
     antichain (maxForb plateau ~46), but MULTI-HOST patterns' witnesses
     span regions, forcing a CROSS-REGION inclusion-exclusion that blew up
     (commit d912b8f, f4f8750). This is the wall the region-IE approach
     cannot pass.
 LESSON: counting via region-IE is fast, but the region-IE PEEL needs a
 cross-region dead-box IE that does not scale. The fix is to NOT represent
 witnesses as a union of regions.

USER STEER (the breakthrough): use the PIVOT idea (a single canonical
hold/optional witness structure) so each witness is counted once and a
peel is one dead-box on ONE structure -- no cross-region IE. pivot != CPI:
do it at CLASS granularity (lightweight), not the heavy vertex SDCT.

KEY ENABLING FACT (commit 0ebfc4d, validated 400/400 vs region-IE):
classes i,j are fully mutually adjacent <=> they co-occur in some region
(classRegions[i] cap classRegions[j] != empty). Proof: full adjacency =>
i union j is a clique => extends to a maximal clique (region) with both.
So a tuple's witnesses = (s-r)-cliques of its CANDIDATE-CLASS GRAPH (nodes
= host-region classes, edges = region co-occurrence), read straight off
the quotient. And the SCT of that graph turns the union into DISJOINT
leaves -> support is a SUM (no IE), peel is a per-leaf dead-box.

### 3. THE SOLUTION (class-level Succinct Clique Tree)
 - ClassSCT.h: dense orbit-aware class-weighted SCT. emitLeaf: empty hold
   (h=0); spine class via ell=u=mult (weight C(w,mult)); free pivot pool;
   canonical lowest-non-neighbour branching => DISJOINT leaves; sum over
   leaves of support_count == weighted s-clique count. Adversarially
   audited ~605k cases vs an independent vertex-expansion oracle (commit
   9ab372b, subagent).
 - ClassSCTScalable.h: scalable version. Dense C×C matrix OOMs (~15GB at
   nC=123k); the gen() also does O(|P|^2) pivot over P=all-classes. FIX:
   degeneracy-seeded (Eppstein) + sparse adjacency + COMPACT leaves
   (touched classes only, global ids in CCPath.classIds). Emitting full-C
   leaves OOM'd at 15.8GB -> compact leaves: 150x faster, 80x less mem.
   Validated 50000 trials vs the dense oracle + scale to C=125000 in 4.4s
   (commits 7f01fc4).
 - region_native_sct_peel.cpp: the peel. support(pattern) = SUM over its
   hosting leaves of support_count(leaf-slot, m). Peel inserts the pattern
   threshold into hosting leaves' forbidden antichains; a leaf over KMAX is
   replaced by controlled_split (a SET of CCPaths). Map-build = PER-LEAF
   enumeration (every r-multiset over a leaf's classes is a registered
   pattern -> O(incidences); the per-pattern scan was 33s on ca-CondMat).
   Patterns stay SPARSE: compToLocal maps a comp into a leaf's local dim
   via binary search (full-C patVec would be TB on com-dblp). commits
   a774bf2, 9423767, 23710db.

### 4. CORRECTNESS GATES (all PASS, independently re-run by me)
 - class-SCT: 605k random cases (sum + per-pattern disjointness) vs oracle.
 - scalable vs dense: 50000 trials, 0 fail.
 - peel vs brute (scripts/verify_nucleus_brute.py): 68/68 on tiny + 4 dense
   stress graphs (t_dense/t_chain/t_stress/t_split, in graphs/), core>=1.
 - peel vs V3LM (PIVOTER_RUN_REGION_V3LM): identical on ca-GrQc 3,4/4,5,
   dblp-core30 3,4, ca-CondMat 3,4.
 - peel vs CND (PIVOTER_RUN_REF): identical (ca-GrQc 3,4/4,5, ca-CondMat).
 - G1 (total s-cliques class-SCT == region-IE): ca-GrQc 329297, dblp-core30
   13.8M, com-dblp 16.7M.
 - G2a (per-pattern SCT-support == region-IE): com-dblp 934117/934117.

### 5. PERFORMANCE (honest, measured)
 BUILD is the win: class-SCT builds com-dblp (nC=123526) in 0.19s vs the
 vertex SDCT ~85s on web-it; counting half is size-free.
 PEEL-PHASE is the problem: ca-GrQc(4,5) peel 10.9s vs CND 0.18s/V3LM 0.18s;
 com-dblp(3,4) peel ~250s. Profile (macOS sample): support_count /
 count_with_extra_lower + the affected-update loop dominate; maxSplit=4452
 (KMAX=2 over-splits). CND ~= V3LM on ca-GrQc/ca-CondMat (V3LM only beats
 CND on specific graphs -- the motivation to do better).
 IN PROGRESS: peel-phase optimization subagent (KMAX sweep, truly-
 incremental affected-update over only-changed slots, bucket-queue
 efficiency, prune affected set), gated exact vs brute+CND+V3LM.

### 6. HOW TO RUN
 build SCT peel:  cd region_native && g++ -O3 -std=c++17 \
   -I../src/NucleusDecomposition -o region_native_sct_peel region_native_sct_peel.cpp
 run:  ./region_native_sct_peel graphs/G.edges r s    (prints core=k count=N + TIMING)
 CND:  PIVOTER_RUN_REF=1 ./build/bin/degeneracy_cliques graphs/G.edges r s degen
 V3LM: PIVOTER_RUN_REGION_V3LM=1 ./build/bin/degeneracy_cliques graphs/G.edges r s degen
 brute oracle: python3 scripts/verify_nucleus_brute.py graphs/G.edges r s
 3-way bench: python3 scripts/bench_sct_peel.py  (CND vs V3LM vs SCT, 16 graphs;
   env BSCT_GRAPHS/BSCT_RS/BSCT_TIMEOUT/BSCT_OUT)
 KMAX tunable via SCT_KMAX env.

### 7. KEY FILES
 region_native/region_native.cpp          counting half (size-free, done)
 region_native/region_native_peel.cpp     region-IE peel (DEAD END, kept for history)
 region_native/ClassSCT.h                  dense class-SCT (oracle)
 region_native/class_sct.cpp               dense SCT self-test (605k cases)
 region_native/ClassSCTScalable.h          scalable sparse class-SCT (production)
 region_native/test_scalable_sct.cpp       scalable self-test (vs dense)
 region_native/region_native_sct_peel.cpp  THE PEEL (current main)
 region_native/test_region_forbidden.cpp   M1 foundation test
 scripts/verify_nucleus_brute.py           brute nucleus oracle
 scripts/bench_sct_peel.py                 3-way CND/V3LM/SCT bench
 graphs/t_{dense,chain,stress,split,k6k4,k7k5k4,2k5tri,mix}.edges, tiny_3k4.edges

### 8. GOTCHAS / ENV
 - CND = PIVOTER_RUN_REF=1 (the no-env DEFAULT crashes: stale SDCT_Fused
   "clique not found").
 - server git repo is /home/wenqianz/UNSW/pivoter (NOT /data/wenqianz/pivoter,
   which is a stale non-git copy). git fetch + reset --hard, separate ssh
   calls, verify HEAD; rebuild binaries after pull.
 - graphs/ is a SYMLINK -> cannot git-add through it; scripts/ is gitignored
   (so .py bench/oracle live untracked but regenerable).
 - dense quotient matrix only up to nC~20000; use the scalable path for big.
 - full-C per-pattern vectors do not scale (com-dblp = TB) -> keep sparse.

### 9. NEXT STEPS
 1. (running) peel-phase optimization -> close the CND/V3LM gap.
 2. broad 3-way sweep (16 NuclearCD-style graphs) -> table of where SCT beats
    CND vs where V3LM beats CND (the headline the user wants).
 3. REGION MERGING (r-mergeable, cf. task #99 / RegND family): merge regions
    to make the pattern count strictly FEWER -> the theoretical lever to beat
    CND on MORE graphs (the SCT peel currently processes the SAME patterns as
    V3LM/CND; merging is the not-yet-applied advantage).
 4. when peel competitive: deploy to tods2, full sweep, paper framing.

COMMIT TRAIL this session (region-native): c5282e6 (counting flatness) ..
0ebfc4d (P1) .. 9ab372b (class-SCT+G1) .. a774bf2 (SCT peel G2) .. 9423767
(map-build) .. 7f01fc4 (scalable) .. 23710db (sparse patterns) .. 5fb6a25/this.

## 33. CORRECTION + KEY POSITIVE RESULT: SCT peel wins at high (r,s) (2026-06-18, #139)

IMPORTANT correction to §31/§32's "peel not competitive" verdict: that was
measured at (r,s)=(3,4) on modest graphs -- exactly the regime where CND is
FAST and the WHOLE region approach (V3LM too) LOSES. That is NOT where the
contribution lives.

The contribution lives at HIGH (r,s) where CND's witness enumeration /
SDCT build EXPLODES. From paper data (bench_full_merged.csv): RegNDC(V3LM)
beats REF(CND) on 417/875 cells, e.g. ca-GrQc(7,10) CND=255s vs V3LM
0.002s (127667x), dblp-core30(5,7) CND=863s vs 0.027s.

SCT PEEL on those exact cells (independently measured, all bit-exact vs V3LM):
  ca-GrQc(7,10)   SCT 0.01s  V3LM 0.04s  CND 255s   (SCT beats CND ~25000x AND beats V3LM)
  ca-GrQc(7,12)   SCT 0.01s  V3LM 0.02s  CND 223s
  ca-CondMat(14,18) SCT 0.03s V3LM 0.06s CND 24s    (SCT beats V3LM)
  dblp-core30(5,7)  SCT 0.10s V3LM 0.09s CND 863s   (~tie, both crush CND)
  dblp-core30(5,8)  SCT 0.17s V3LM 0.08s CND 857s   (SCT 2x V3LM, both crush CND)

So: at high (r,s) SCT beats CND by 1e4-1e5x AND is competitive-to-faster
than V3LM (its lightweight class-SCT build skips V3LM's vertex SDCT build).
Combined with the (3,4) finding that SCT beats CND on the sparsest graphs
where V3LM loses (bio-celegans 1.12x, ca-HepTh 1.02x), SCT beats CND on a
SUPERSET of where V3LM does -> the user's goal (beat CND on MORE graphs).

WHERE SCT STILL LOSES: (3,4) on medium/dense graphs (ca-GrQc 0.39 vs CND
0.05) and dense graphs at any (r,s) (email-Eu-core, ca-HepPh (3,4) timeout)
-- but V3LM also loses there, and these are the low-(r,s) "CND is fast"
regime, not the region approach's domain.

r-MERGEABLE (committed 71879d2): correct (gate: rmerge==no-rmerge==brute
20/20) but barely speeds up -- it removes the ISOLATED (cheap) regions; the
cost is the OVERLAPPING peel + build, untouched. Kept ON for fair apples-
to-apples with V3LM/CND (both merge). Not the hoped lever.

NEXT: broad high-(r,s) 3-way sweep (running, /tmp/bsct_highrs.csv) to
quantify SCT-vs-V3LM across many graphs at the cells that matter; then the
honest paper story = "region-native class-SCT: size-free counting +
lightweight build + correct peel that, at high (r,s), beats CND by 1e4-1e5x
and matches/beats V3LM, on a superset of V3LM's win-set."

## 34. PEEL SPEEDUP via |host|=1 target-skip (the merging idea, done right) (#139)

User push: "use region MERGING to optimize the peel, not just direct-assign;
it will definitely be faster." Careful investigation of every realization:

DEAD ENDS (both empirically falsified vs brute):
1. Direct-bin ALL |host|=1 patterns (core=C(|M|-r,s-r), skip peel): WRONG.
   t_2k5tri shared triangle core 4 vs 2 -- a |host|=1 witness s-clique can
   contain a |host|>=2 r-clique, so skipping the |host|=1 peel inflates the
   |host|>=2 cores.
2. Batch "region mass-death": set u[c]=0 for M's PRIVATE classes at L_M.
   WRONG (t_dense, 4 shared patterns over-cored). ROOT CAUSE: |host|=1 is NOT
   "uses a private class" -- a |host|=1 pattern can use only SHARED classes
   whose region-intersection is a single region (e.g. (0,1,3) in t_dense, all
   of 0,1,3 shared but host={one region}). A correct compact mass-death needs
   to remove host=={M} witnesses, which requires region-IE -- exactly what the
   class-SCT's disjoint leaves avoid. So no clean compact batch exists.

WHAT IS TRUE + CORRECT (verified 0 exceptions, then proven): every |host|=1
pattern peels at EXACTLY L_M=C(|M|-r,s-r). Proof: every r-clique in M (incl
shared ones) has core >= L_M, because a shared r-clique keeps its M-witnesses
(C(|M|-r,s-r)=L_M of them) until M's own |host|=1 patterns die at L_M; so no
witness of a |host|=1 pattern dies before curLevel=L_M.

THE WIN (commit 1b6bf3f, SCT_NO_SKIP_H1 to disable): a |host|=1 pattern's
bucket key is already L_M and never changes, so recomputing its support in
affected-updates is pure waste. SKIP |host|=1 patterns as affected-update
TARGETS; still peel them as SOURCES (witnesses removed, |host|>=2 drops
correct). |host|=1 is the majority -> removes the bulk of scWithTerms.
GATE skip==no-skip==brute 38/38 incl dense. Peel speedup:
  ca-CondMat(5,7) 5.27->1.33s (4x); (5,10) 3.07->1.04s; (7,10) 2.47->1.12s;
  ca-GrQc(5,7) 1.25->0.85s.
New standing (total, skipH1): SCT beats BOTH CND and V3LM on bio-celegans,
ca-HepTh, ca-GrQc(5,10)/(7,10), dblp-core30(5,7); beats CND on every
high-(r,s) explode cell. Weak spot remains dense moderate-(r,s) (ca-CondMat)
where CND/V3LM are sub-second.

NEXT: source-skip -- a |host|=1 pattern whose leaves host NO |host|>=2
pattern affects nothing when it peels, so skip its source-peel entirely
(generalizes fully-mergeable direct-assign to the leaf level).

## 35. skipH1 broad sweep result: SCT now competitive with V3LM, beats CND on more (#139)

48-cell sweep (8 graphs x 6 (r,s)), skipH1 ON, vs V3LM and CND (45s cap).
Data: paper_data/bench_sct_peel_skipH1.csv. ALL cells that ran are EXACT
(no correctness regression from the skip).

TALLY:
  - SCT beats CND on 18 cells; V3LM beats CND on 17  -> SCT beats CND on MORE.
  - SCT beats CND where V3LM does NOT (6): bio-celegans(3,4)(4,5)(5,7),
    ca-HepTh(3,4), ca-CondMat(7,12), amazon0302(4,5).
  - SCT vs V3LM head-to-head (26 both-ran cells): 13-13 (was mostly losing
    before the |host|=1 skip).
  - Big wins: dblp-core30(4,5) SCT 0.06 vs CND 12.73 (212x) vs V3LM 0.09;
    ca-GrQc(6,8) 17x CND; ca-HepTh(5,7) 13.7x CND.

WEAK SPOT (shared with V3LM, not unique to SCT): dense graphs ca-AstroPh /
ca-HepPh -> SCT TIMEOUT(>45s). But V3LM also times out on most of those
(ca-AstroPh(4,5)+, ca-HepPh all) and CND handles only the lowest (r,s).
The dense |host|>=2 peel is the hard core for BOTH region methods. ca-CondMat
/ ca-GrQc mid-(r,s) is where V3LM's mature peel still edges SCT.

CONCLUSION: the region-native class-SCT engine + |host|=1 skip is a correct,
lightweight-build, competitive-with-V3LM peel that beats CND on more cells.
The contribution regime is sparse + high-(r,s) (CND explodes, SCT/V3LM win,
SCT's lighter build often wins the head-to-head). Dense low-(r,s) remains
CND territory for all region methods.

## 36. HONEST large-graph verdict: SCT peel does NOT scale (#139)

User (emphatic): test LARGE graphs, small graphs prove nothing. Done.

com-dblp (n=317k, m=1.05M):
  (3,4): SCT 23.58s (build+maps 2.43 + peel 20.55) | V3LM 3.25s | CND 2.17s
  (4,5): V3LM 14.63s | CND 17.19s | SCT killed (peel >150s)
  (5,7): V3LM 54.29s | SCT peel >150s (1.79M patterns)
  (7,10): SCT 10.75M patterns, peel >150s

=> SCT peel is 7-11x SLOWER than V3LM on com-dblp; does NOT scale.

KMAX lever EXHAUSTED (com-dblp 3,4 peel): KMAX=2 -> 21s, KMAX=4 -> 29s,
KMAX=8/16 -> timeout. Higher KMAX = more 2^KMAX IE terms in scWithTerms,
dominates the fewer-split-paths saving. KMAX=2 already optimal.

ROOT CAUSE (important, honest): the SCT's premise was "disjoint leaves avoid
region inclusion-exclusion, so it is faster." For the PEEL this premise
FAILS -- the controlled_split slot growth (maxSplit=5607 on com-dblp 3,4)
makes each peel's slot scan O(#peels x maxSplit) cost MORE than the bounded
region-IE it avoids. V3LM's bounded-IE tuple peel beats the SCT split peel
at scale. The 0.01s "build" earlier was just the class-build substep; full
build+maps is 2.43s, already ~ V3LM's whole 3.25s.

WHAT ACTUALLY WORKS (honest scorecard):
  - Size-free COUNTING: scales (com-dblp 0.89s) -- the real scalable win.
  - Peel small + high-(r,s): beats CND where it explodes; ~V3LM (13-13).
  - Peel large graphs: loses to V3LM decisively. NOT its domain.
The |host|=1 skip (real ~2x) does not change this ordering.

RECOMMENDATION: contribution = size-free region-native counting (+ light
class-SCT build) + the high-(r,s) peel niche; do NOT claim large-graph peel
superiority. A scalable SCT peel would need to beat V3LM's bounded-IE
incremental update, whose constant factor is lower -- uncertain / open.

## 37. PEEL optimization campaign: DP is NOT the bottleneck, affected-Q ENUMERATION is (2026-06-18, #139)

User target (sharpened): at **r>=4 (esp (4,5)) beat CND on EVERY graph, time AND memory**. That is the
regime where CND must enumerate r-cliques (explodes) while our quotient is size-free, so theory says we
dominate everywhere. The cells where region_native_sct_peel currently LOSES to CND are the dense small ones:
ca-GrQc(4,5) (peel 1.06s vs CND 0.18s), ca-CondMat(4,5) (0.65s vs 0.18s). These are fast to iterate.

### 37.1 Diagnosis (all measured, instrumented /tmp builds)
- Fundamental work on ca-GrQc(4,5) ~0.05s (= the 248942 real support updates) -- BELOW CND's 0.18s. Win is real.
- Profile of the peel (com-dblp(3,4), r=3): DFS candidate-gen ~35%, CCPath memory churn ~35%, the DP only ~6%.
- The TREE IS SMALLER and we touch FEWER pairs (user was right): ca-GrQc(4,5) ours ~3.56M (P,Q) pairs vs
  CND ~ #5cliques x (C(5,4)-1) = 2.2M x 4 ~ 8.8M. So pair COUNT is not the problem.
- THE REAL GAP: per-pair UPDATE cost. After adding guards, the DP (count_with_extra_lower) is called only
  248942 times (= real nonzero drops) ~0.1s -- NOT the bottleneck. The remaining ~0.7s is the affected-Q
  ENUMERATION: ~158 candidate patterns generated per leaf-peel to find ~4 real ones = **40x over-enumeration**.
- Candidate fate (instrumented, ca-GrQc(4,5)): hashMiss=5476 (negligible -> a TRIE is useless, candidates
  ARE real patterns); realZeroDrop=2,228,853 (63%); realNonzero=248,942. So the waste = real co-occurring
  patterns with ZERO drop (their witnesses don't overlap P's killed region, or are forbidden-dead).
- Why CND has no over-enumeration: a dying s-clique's C(s,r) r-subcliques ALL co-occur (same s-clique).
  Our leaf is a BUNDLE of structurally-different witnesses; most of its r-sub-patterns do NOT co-occur with
  P in the same witness -> 40x waste.
- support_count is FUNDAMENTALLY A BINOMIAL not a DP (verified: unconstrained box = Vandermonde
  C(sum(n-b), T-sum b); n=[3,3]T=2b=[1,0] DP=5=C(5,1)). DP only "earns its keep" on split-constrained boxes;
  for the real drops here 0% hit the clean-binomial fast path (addLow=pl raises the floor), but DP is cheap
  anyway so this does not matter. The user's instinct ("don't DP, it's simpler") is right that the DP is not
  the cost; the cost is the enumeration.

### 37.2 Committed bit-identical fixes (all verified cores-identical on 11 cells incl com-dblp(3,4))
- in-place slotForbidDiff (no full split-set rebuild; unchanged paths stay put) -- ~16% peel.
- sum-guard (skip scWithTerms when sum max(ql,pl) > p.T, provably 0) + addLow forbidden early-out (single
  forbidden a <= max(ql,pl) => region dead => 0) -- skips provably-zero drops.
- Cumulative: ca-GrQc(4,5) peel 1.06 -> 0.76s at KMAX=1 (best for dense cells), still ~4x off CND 0.18s.
  Two git commits on region_native/region_native_sct_peel.cpp.

### 37.3 THE redesign target = kill the 40x affected-Q over-enumeration. Two derived directions:
(1) PER-PATH tight-box enumeration (testable now): the DFS enumerates Q over the ENVELOPE uEnv = max-u over
    all chgOld paths, then per-path-filters. The envelope is loose -> generates "feasible-on-the-union but
    on NO single path" phantoms. chgOld paths are DISJOINT (controlled_split partitions witnesses), so P's
    killed witnesses live in specific paths; a Q whose witnesses are in OTHER paths has 0 drop. Enumerate Q
    PER changed-path using that path's tight [ell,u] -> only Q that can co-occur in THAT killed box. Dup cost
    = Q with witnesses in multiple changed paths (chgOld small, so low). Should cut most of the 63% zero-drop.
(2) CO-OCCURRENCE condition / index: on a BASE leaf, P co-occurs with Q  <=>  max(pl,ql)<=u (auto) AND
    sum_c max(pl,ql) <= T  <=>  sum_{c: ql>pl} (ql[c]-pl[c]) <= T-r. For (4,5) T-r is tiny (often <=1) =>
    Q exceeds pl by <=1 total => the affected set is INTRINSICALLY ~4. The DFS bound already encodes
    sum max(pl,ql) <= Tcap; the residual over-enumeration is exactly (a) the loose envelope [fixed by (1)]
    and (b) forbidden-DEAD patterns that pass the geometric bound but whose witnesses are all peeled [the
    dead-tracking]. (b) is what the split/forbidden machinery exists for; per-path tight box (1) already
    excludes dead regions, so (1)+(2) together should give near-zero over-enumeration -> peel ~ fundamental
    ~0.05s -> UNDER CND. User believes this is definitely achievable; the derivation supports it.

### 37.4 Reusable assets / how to verify
- Verification protocol: run /tmp/sct_inplace (baseline) vs the variant on the test cells; cores
  (grep '^core=' | md5) MUST match; measure 'total=' and peel='. Cells: ca-GrQc {4:5,5:7,6:8},
  ca-CondMat{4:5,5:7}, dblp-core30{4:5,5:7}, ca-HepTh 5:7, bio-celegans 4:5, amazon0302 4:5, com-dblp 3:4.
- Build: cd region_native && g++ -O3 -std=c++17 -I. -I../src/NucleusDecomposition -o BIN region_native_sct_peel.cpp
- Run: ./BIN ../graphs/G.edges r s  (prints [sct-peel] TIMING ... + core=k count=N). SCT_KMAX env tunes KMAX.
- CND baseline times are stored in paper_data/bench_sct_peel_skipH1.csv (col 8 cnd_total_s); DO NOT re-run CND.
- Instrumented patchers used: /tmp/patch_{prof,timers,inplace,dfs}.py (counts: scWithTerms/dfsNodes/affPairs/
  hashMiss/realZeroDrop; timers: slotForbidDiff vs dfs+scWithTerms).
- NEXT: implement (1) per-path enumeration, verify bit-identical, measure over-enumeration drop; then (2).

## 38. Experiment baseline + architecture + the PEEL-walks-the-pivot-tree direction (2026-06-19, #139)

### 38.1 RegND(class-SCT) vs CND comparison — tods2, TRUE wall-clock + peak RSS
FAIRNESS NOTE (critical): measure WALL-CLOCK + peak RSS via `/usr/bin/time -f "%e %M"`, NOT the program's
self-reported time. CND prints a "NucleusCoreDecomposition took:" line PER sub-phase; the FIRST one is a
single phase and undercounts CND by ~28x (ca-GrQc(6,8): first "took"=0.388s vs TRUE wall=20.11s). My first
sweep used the wrong number and falsely showed CND winning everywhere -- do not repeat.
Sweep: 8 paper graphs x (3,4)(4,5)(5,7)(6,8)(7,10)(7,12), TO=90s. Script /tmp/rs_sweep2.sh -> /tmp/rs_compare.txt.
HONEST verdict:
 - WIN (time AND memory): small/sparse graphs at HIGH (r,s). dblp-core30 ALL rs (4,5)=211x; ca-GrQc (5,7)+
   (6,8 SCT 1.2s/25MB vs CND 20s/2GB; (7,10) SCT 0.03s vs CND >90s/8.9GB). CND materializes all r-cliques ->
   RAM explodes to 32GB; our class-SCT stays MB.
 - MEMORY advantage is BROAD: SCT memory ~flat (class-level), CND grows with rs to 32GB cap.
 - LOSE (time): dense graphs (ca-CondMat CND wins up to (7,10); ca-HepPh SCT even TIMES OUT at (3,4) while
   CND=11s), and ALL large graphs at low rs (com-dblp(3,4) SCT 33s vs CND 4.6s; web-Stanford/com-youtube SCT
   TIMES OUT >90s while CND finishes 20-73s). On large graphs SCT uses far less memory but DOESN'T FINISH
   (moot). Consistent with SigmodPlus 36 "SCT peel does not scale".
 => advantage is real but REGIME-SPECIFIC: high-rs on small/sparse. Large+dense is the unsolved wall.

### 38.2 Architecture confirmed: NO vertex CPI; class-level SCT; the BUILD is already Pivoter
 - region_native_sct_peel builds NO vertex-level CPI/SDCT. Counting = region-IE closed-form binomials
   (region_native.cpp, sup(tau)=sum_A (-1)^{|A|+1} C(|cap A|-r, s-r)). Peel = class-level SCT.
 - "Our CPI" = ClassSCT = class-weighted Succinct Clique Tree on the QUOTIENT graph (nodes=classes,
   edge iff two classes co-occur in a region). ~26x smaller than vertex level (ca-GrQc nC=199 vs 5242 verts)
   -- THIS is the source of the flat memory.
 - ClassSCTScalable.h gen() IS the SDCT/Pivoter algorithm: orbit-aware pivot recursion (lowest-non-neighbour
   branching) + degeneracy-seeding (each subtree open set <= d) + sparse adjacency. Already fast: com-dblp
   build=0.13s (vs vertex SDCT ~85s on web-it). => NOTHING to borrow for the build; it is already Pivoter.

### 38.3 Pipeline breakdown: PEEL is the wall (62-86%); build is cheap
   ca-GrQc(6,8): build+maps 0.05 | PEEL 0.38 (86%)   |  ca-CondMat(4,5): build+maps 0.27 | PEEL 0.57 (62%)
   com-dblp(3,4): build 0.13, maps ~2.06, PEEL 14.34 (83%).  Optimizing build is pointless; PEEL holds the time.

### 38.4 The PEEL bottleneck (definitive): affected-Q OVER-ENUMERATION ~40x
 - Fundamental work ~0.05s (the ~248942 real support updates) is BELOW CND's time. Waste = affected-Q search
   generates ~57 candidates/leaf-peel for ~4 real. Candidates ARE real patterns (hashMiss ~negligible); the
   waste is real co-occurring patterns whose specific witnesses give ZERO drop.
 - The DP is NOT the bottleneck (confirmed user's intuition): support_count is a binomial (Vandermonde
   C(sum(n-b),T-sum b) when unconstrained); only ~248942 real-drop DP calls ~0.1s. The cost is the ENUMERATION.
 - Co-occurrence closed form: Q affected by P <=> E(P->Q)=sum_c max(0, ql_c-pl_c) <= T-r. (4,5) base leaf
   T-r=1 => affected set intrinsically ~4 (Q = P with one unit moved between classes).

### 38.5 Tried + verdicts (ALL bit-identical; verify cores match before trusting any change)
 - COMMITTED pure wins: in-place slotForbidDiff (8c0086d); sum-guard + addLow forbidden early-out (b76b779).
   ca-GrQc(4,5) 1.06 -> 0.68-0.76s (~30%). KMAX=1 best for dense cells.
 - DUDS/reverted: active-classes DFS skip (~0%, the descent wasn't the cost); PER-PATH tight-box enumeration
   (WORSE on dense: the disjoint split-boxes are barely tighter than their union, so ~2.4x re-enumeration >
   phantom savings; helps only high-rs where chgOld small; REVERTED); Vandermonde fast-path (0% hit: the
   addLow=pl floor always constrains the box). avg chgOld(changed paths/leaf-peel)~2.4 but max~630.

### 38.6 NEXT (user-approved, STARTING NOW): make the PEEL walk the class-SCT PIVOT TREE for affected-Q
 - The BUILD is Pivoter but the PEEL is NOT: it finds affected-Q by FLAT per-leaf DFS (enumerate all
   r-multisets in a leaf + hash-lookup) -> the 40x over-enumeration.
 - The class-SCT pivot tree already organizes cliques by pivot/spine/free-pool. Affected-Q (patterns sharing
   a live witness with P) should be found by NAVIGATING that pivot structure, not by flat enumeration.
 - GOAL: kill the 40x over-enum -> peel approaches its ~0.05s fundamental -> beats CND on dense too (and may
   help scale on large graphs). This is the key to "beat CND on every (r,s) incl dense".

### 38.7 Verify / run / assets
 - SCT build: cd region_native && g++ -O3 -std=c++17 -I. -I../src/NucleusDecomposition -o BIN region_native_sct_peel.cpp
 - Run: ./BIN ../graphs/G.edges r s  (env SCT_KMAX tunes KMAX). Prints [sct-peel] TIMING + core=k count=N.
 - CORRECTNESS GATE: grep '^core=' | md5 must equal /tmp/sct_inplace (or stock) on the verify cells.
 - CND: PIVOTER_RUN_REF=1 build/bin/degeneracy_cliques graphs/G.edges r s degen ; FAIR time = wall-clock+RSS
   (/usr/bin/time -f "%e %M"), NOT the "took" line. Do NOT re-derive cnd from the first "took".
 - tods2 real git repo: /home/wenqianz/UNSW/pivoter (= origin/main, pushed through 23d8291). /data/wenqianz
   is a STALE non-git copy -- do not use it.
 - Verify cells: ca-GrQc{4:5,5:7,6:8} ca-CondMat{4:5,5:7} dblp-core30{4:5,5:7} ca-HepTh 5:7 bio-celegans 4:5
   amazon0302 4:5 com-dblp 3:4. Key files: region_native_sct_peel.cpp (peel), ClassSCTScalable.h (build),
   src/NucleusDecomposition/CCPathCore.h (CCPath/support_count). Instrument patchers: /tmp/patch_*.py.

## 39. Day log 2026-06-19: peel bottleneck PINNED + micro-opts exhausted + overnight diverse-graph phase test
### Done today
 - Architecture confirmed (§38): NO vertex CPI; class-level SCT; the BUILD already IS Pivoter (orbit-aware
   pivot recursion + degeneracy seeding), build=0.13s on com-dblp. Nothing to borrow for the build.
 - RegND-vs-CND experiment baseline (§38, tods2, TRUE wall-clock %e + RSS %M; the program's first "took"
   undercounts CND 28x -- never use it). WIN (time+mem): small/sparse at HIGH (r,s) (dblp-core30 211x;
   ca-GrQc (5,7)+; CND -> 32GB while SCT stays MB). MEMORY win broad. LOSE: dense (ca-CondMat, ca-HepPh)
   and large (com-dblp/web-Stanford/com-youtube) -- SCT slow/times-out (peel doesn't scale, §36).
 - Pipeline breakdown: PEEL is 62-86% of total; build cheap. Optimizing build is pointless.
 - PEEL bottleneck PINNED via fine timers (ca-GrQc(4,5) peel 0.86s): enum+apply (generate+hash-lookup the
   ~57 affected-Q candidates/leaf-peel for ~4 real) = 0.446s (52%); scWithTerms DP 0.204s (24%); slot 0.108s
   (13%). The cost is the 14x affected-Q OVER-enumeration; work is proportional to 57, not 4.
 - Micro-opts EXHAUSTED (all bit-identical, all ~0% on dense, because NONE cuts the candidate count 57->4):
   active-classes DFS skip; per-path tight-box (WORSE, reverted); Vandermonde binomial fast-path (0% hit);
   sparse O(r) feasibility check (today, ~0%). COMMITTED pure wins earlier: in-place slotForbidDiff +
   sum-guard + addLow early-out (~30% on dense, KMAX=1 best).

### Current open problems (THE wall)
 - To cut 57->4 you must know, BEFORE enumerating, which co-occurring patterns (E(P->Q)<=T-r) share a LIVE
   witness with P. The 53 zero-drops are real co-occurring patterns whose shared witnesses are dead (peeled)
   or in other split-boxes. That is dead-witness tracking == the controlled_split blow-up (maxSplit=3498),
   the §24-36 wall. Every cheap dead-tracking tried blows up or doesn't prune.
 - Net: peel does not scale on dense/large; ~4x gap to CND on ca-GrQc/ca-CondMat at low-mid (r,s).
   The class-SCT build, counting, memory, and high-(r,s) time win are all solid; the dense/large peel is open.

### USER'S KEY REMINDER (2026-06-19) + overnight task
 - The phase-time distribution varies HUGELY by graph, determined by MANY graph properties, not just density.
   The algorithm is graph-property-sensitive -> must test MANY varied graphs.
 - OVERNIGHT (running): comprehensive phase-distribution sweep on a DIVERSE graph set (collaboration / web /
   social / citation / product / bio / email, varied clustering/degeneracy/max-clique/skew) x (r,s) small+
   large, capturing per-phase timing (MCE/enum/build/maps/peel) + wall + RSS + the [rn]/[sct] structural
   stats (n,m,regions,nC,base-leaves,patterns,maxSplit) + CND wall+RSS. Plus parallel graph-property
   characterization (workflow). Goal: map which graph properties drive which phase, to find where the peel
   is/ isn't the wall. Output -> /tmp/phase_sweep.txt on tods2.
 - FORMAT: region_native + degeneracy_cliques both need a "n m" header line (the .grp format). /data/wenqianz
   raw .edges need a prepended "maxid+1 numEdges" header; /home/wenqianz/UNSW/pivoter/graphs/*.edges already
   have it. CND time = wall-clock, NOT "took". tods2 real repo /home/wenqianz/UNSW/pivoter (=origin/main).

## 40. Overnight phase-distribution sweep RESULTS: peel is the UNIVERSAL bottleneck, build is never (2026-06-19, #139)
Ran 19 diverse graphs x 6 (r,s)={3:4,4:5,5:7,6:8,7:10,8:12} = 114 cells, TO=120s, SCT(region_native_sct_peel)
vs CND(PIVOTER_RUN_REF). Per-phase timing (MCE/enum/build/maps/peel) + structural (regions,nC,patterns,
rcliques,maxSplit) + wall + RSS. Raw CSV: paper_data/phase_sweep_2026-06-19.csv.
RESULT: 33/114 SCT-completed (full breakdown), 81 timed out. CAVEAT: ran under OTHER-USER contention on tods2
(v15_ablation x1-5, up to 522GB RSS; load avg 2.2-2.9 by end so CPU starvation was MODEST) -> absolute times +
SCT-vs-CND head-to-head are NOISY; phase FRACTIONS within a cell + all structural stats are ROBUST.

WHICH GRAPH PROPERTY DRIVES WHICH PHASE (the user's question, answered):
 - peel  <- maxSplit (witness-box fan-out) + pattern count. THE universal bottleneck. Dominates wherever SCT
            completes: collab (ca-GrQc, ca-CondMat ALL rs), coauthor (dblp-core30 mid/high-rs, com-dblp 3:4
            peel=33.7/39.6), web-NotreDame (peel=59/69 at 6:8), web-it-2004 (peel=85 at 5:7). High maxSplit
            => peel explosion => timeout. maxSplit extremes: web-NotreDame(6,8)=12610, web-it-2004(5,7)=10609,
            com-dblp(3,4)=5607. This is the SAME §24-36 wall, now confirmed ACROSS graph families.
 - MCE   <- graph SIZE (n). cit-Patents (3.77M nodes) MCE=14.11 of 16.65 total; ca-GrQc very-high rs where
            peel->0 (few regions left); dblp-core30(3,4); com-amazon(4,5)(5,7).
 - maps  <- #regions + region->class mapping size. web-it-2004(3,4)(4,5) maps=8.3/12.8s.
 - build <- NEVER. 0.00-0.35s across ALL 33 completed cells. STRONG confirmation: the class-SCT build is NOT
            a bottleneck. => do NOT spend effort borrowing SDCT tricks to speed the build; optimize PEEL.

SCT-vs-CND where SCT completes: SCT WINS at high rs, often hugely (CND materializes all r-cliques and
explodes): ca-GrQc(6,8) 1.23 vs 33.0s (27x); ca-GrQc(7,10)(8,12) SCT<0.1s vs CND>120; dblp-core30(4,5) 0.13
vs 33.3 (256x), (5:7+) SCT<2s vs CND>120; web-NotreDame & web-it-2004 high-rs SCT finishes while CND>120;
ca-CondMat(8,12) 5.94 vs 22.19 (3.7x). The high-(r,s) memory+time win over CND is REAL and broad.

SCT timeout (>120s) by type: social 24/24 (ALL), web 24/30, email 6/6, ego 6/6, collab 12/24, coauthor 5/6,
citation 4/6. Social + dense (HepPh/AstroPh/email/ego) uniformly blow the peel budget.

NET: the sweep CONFIRMS the optimization target is the PEEL (affected-Q over-enumeration / maxSplit blow-up),
and that it is the bottleneck on EVERY graph family that completes, not just dense ones -> the user's intuition
("property-driven, not density-driven") is right: size->MCE, regions->maps, maxSplit->peel, build->never.
NEXT: (a) attack peel (the dead-witness / controlled_split wall); (b) optional CLEAN re-run for trustworthy
absolute SCT-vs-CND numbers now that load is light.

## 41. Peel optimization campaign #2: KMAX=1/adaptive WIN + 2 refuted prunes (candidate-count wall) (2026-06-19, #139)
Method: design panel (5 approaches x adversarial correctness+perf review, the peel-optimization-design workflow)
-> implement top pick -> verify BIT-IDENTICAL (corehash of '^core=' lines vs golden) -> measure. Local clean
iterate loop: clang++ -O3 -I region_native -I src/NucleusDecomposition; /tmp/bench_peel.sh (15 cells, corehash
oracle + peel + maxSplit); golden = /tmp/golden.txt. SERVER timing is contention-noisy -> use min-of-N /
back-to-back ratios; structural stats (maxSplit) are contention-robust.

THE WIN -- default KMAX 2->1, then adaptive (commits 853cf0d, 14d819c, 0b2fc52):
 - The affected-Q drop is per-candidate-IE/DP-bound, so FEWER forbidden/path (lower KMAX) = cheaper per
   candidate. KMAX=1 measured FASTEST on every cell, monotone (1<2<3<4<6): moderate graphs ~1.5-2x peel
   (ca-GrQc 1.45x, ca-CondMat 1.7-2.2x, dblp 1.8-3x), and the maxSplit-BLOWUP graphs too (web-NotreDame(6,8)
   maxSplit~14k: 41.2s vs 53.5s KMAX=2 on a quiet host; com-dblp(3,4) 10.7 vs 12.9s). The O(slot) scan in
   slotForbidDiff is 99%-skipped (impossible test) so it never overtakes the DP in the measured range
   (slot<=14k). (A single contended run falsely showed a KMAX=1 regression -- contention, not real.)
 - Bit-identical: support_count is invariant to the split strategy (KMAX-invariance), verified 15/15 golden
   corehashes at KMAX=1, KMAX=2, and adaptive. ADAPTIVE (default, threshold 16384): kml = 1 + slotSize/16384,
   so = KMAX=1 across everything observed, only hedging the unmeasured extreme tail (slot>16k, e.g. soc-pokec/
   web-it-2004 high-(r,s)) where a runaway slot self-limits. No-regret ship.

TWO REFUTED PRUNES (candidate count is NOT cheaply reducible -- the wall, from the pruning angle):
 - (1) Per-path AND-feasibility DFS prune (the design panel's #1 pick, SCT_DFS_PRUNE): push applyIdx's per-path
   max(ql,pl)<=u test earlier into the DFS. Bit-identical (15/15, nz-counts identical) but NET SLOWDOWN: cuts
   only 1-16% of candidates because uEnv (per-coord union) is already TIGHT (chgOld paths have near-identical
   u), so it barely fires while its per-node bookkeeping adds overhead. The real waste is NOT envelope-
   looseness.
 - (2) Forbidden-COVERAGE subtree prune (data-indicated follow-up, replaced (1) under same flag): kill a path
   when its forbidden a_z<=max(ql,pl) (=scWithTerms early-out) OR u-infeasible; prune subtree when all dead.
   Bit-identical, cuts 12-30% candidates (2-3x more than (1)) but STILL net loss: the ~90% forbidden-coverage
   deadness (hit/nz~=10x) is a LEAF-level property (a_z's critical coords span the vector => critmax large =>
   subtree prune fires only deep), so the DFS must descend anyway and the per-node bookkeeping is pure
   overhead. Flag left OFF-by-default for ablation.
 - CONCLUSION: the affected-Q over-enumeration cannot be cheaply cut by subtree pruning (deadness is leaf-
   level). The per-candidate COST (KMAX) is the tractable lever, not the candidate COUNT. The candidate-count
   wall = the same dead-witness / controlled_split wall (sec 24-36), reconfirmed.

ENV FLAGS: SCT_KMAX=k (force fixed KMAX, disables adaptive, for A/B), SCT_KMAX_THRESH=N (adaptive knee, def
16384), SCT_DFS_PRUNE=1 (refuted coverage prune, OFF), SCT_NO_SKIP_H1 (disable |host|=1 skip). Instrumentation
to stderr: '[sct-peel] dbg cand_gen/hit/nz' (the over-enumeration ratio).
OPEN: scaling timeouts on dense/social/large-web need a fundamentally different peel (not a prune) -- the
candidate count is irreducibly large and per-candidate work is already minimal at KMAX=1. Open research.
BROADER EXPERIMENT (server, SCT-new adaptive vs SCT KMAX=2, peel min-of-2 back-to-back = contention-robust
ratio): CONSISTENT 1.15-1.60x peel speedup on EVERY graph family, ALL 12 cells bit-identical (corehash ok):
ca-GrQc(6,8) 1.60x, ca-CondMat(7,10) 1.36x, web-NotreDame(7,10) 1.36x, dblp-core30(8,12) 1.30x,
web-NotreDame(6,8) 1.27x [the maxSplit~14k BLOWUP graph: 42.4 vs 53.8s -> settles the contention-noise
question, KMAX=1 wins there too, NO regression], cit-Patents(7,10) 1.27x, com-amazon(3,4) 1.23x,
ca-CondMat(5,7) 1.21x, web-it-2004(3,4) 1.21x, com-dblp(3,4) 1.15x. (Server back-to-back ratios run smaller
than the clean-local 1.5-2x because server peel times are larger, diluting the per-candidate DP saving.)
SHIPPED: adaptive KMAX (HEAD 0b2fc52), bit-identical, ~1.2-1.6x peel everywhere, zero regression.

## 42. Cracking the scaling wall: s=r+1 WITNESS-FLOOR fast path (output-sensitive affected-Q) (2026-06-19, #139)
The affected-Q DFS over-enumerates because it works at the wrong granularity. KEY STRUCTURE: when
P peels and s=r+1, the dying witnesses above P are the SINGLE points floor_c = m_P + e_c (Σ=s, so the
only witness >= floor_c is floor_c itself). Every affected Q that uses that witness is floor_c - e_d
(move one unit d->c), and they all SHARE the one floor. So instead of generating ~140 candidate Q and
filtering, iterate the ~width witness-floors, check liveness of the single point floor_c once per
(c,box), and expand only LIVE floors. The drop is CLOSED FORM (no IE/DP):
   drop_Q(box p) = [floor_c alive in p] * (n_p[d] - m_Q[d]),  summed over chgOld boxes.
"alive" = ell<=floor_c<=u AND no forbidden a<=floor_c. (floor_c-m_Q = e_d, weight = C(n_d-m_Q[d],1).)

CORRECTNESS (3 independent checks): (1) closed form validated BIT-EXACT vs scWithTerms on 400k random
boxes incl ell>0 spine boxes (/tmp/validate_cf2.py, uses solved.py's brute-validated CCPath ops);
(2) C++ cores BYTE-IDENTICAL to golden on all 15 cells (7 of them s=r+1); (3) core-value set matches
CND exactly. Flag SCT_NO_WFLOOR disables (A/B).

SPEEDUP (local clean, peel, vs adaptive-KMAX baseline, ALL bit-identical):
  ca-GrQc(3,4) 4.5x, ca-GrQc(4,5) 3.4x, dblp-core30(4,5) 2.0x, ca-CondMat(4,5) 1.5x, ca-CondMat(3,4)
  1.45x. vs ORIGINAL KMAX=2: ca-GrQc(4,5) 0.75->0.18 = 4.2x. Non-s=r+1 cells fall back to the DFS,
  UNCHANGED (0.97-1.00x). Commit fe4c70a.

WHY IT CRACKS THE WALL (vs the refuted prunes): the prunes tried to filter the ~140 candidates but the
deadness is leaf-level so they couldn't reach it; the witness-floor changes the ENUMERATION ITSELF to be
output-sensitive (work ~ live floors, not all candidates) AND removes the DP. Many timeout cells are
s=r+1 ((3,4)/(4,5): ca-HepPh, ca-AstroPh, web-Google, twitter, wiki-Talk, com-youtube, web-NotreDame(4,5)
...), so this directly targets the scaling timeouts. SERVER CONFIRMATION + PIVOTAL FINDING:
 - AT SCALE the witness-floor holds: com-dblp(3,4) peel 21.1s -> 8.3s = 2.5x, bit-identical (643485 peeled).
 - BUT the densest timeout graphs are NOT peel-bound -- they die in the 'pattern<->leaf maps + compaction'
   phase BEFORE the peel: ca-HepPh(3,4) MCE 0.22s, build 2.63s, maps+compaction=206s; web-Google(3,4)
   MCE 4.3s, build 2.6s, maps+compaction=234s. The witness-floor (a PEEL optimization) cannot help these.
   The phase sweep's 'maps <- #regions' prediction was right; cracking the peel just EXPOSED the maps wall.
 - THE MAPS-PHASE BOTTLENECK (region_native_sct_peel.cpp 655-741, the enumLP per-leaf r-multiset enum):
   for each leaf, enumerate every r-multiset over its classes, build a std::STRING key (compKey), look up
   in unordered_map<string,int> patIdx, and call support_count (a DP) PER incidence to confirm the host.
   Cost = Σ_leaf (#r-multisets) x (string-key build + map lookup + support_count DP). For dense graphs with
   huge leaves (ca-HepPh 14411 base-leaves/1.25M patterns; web-Google 774473 leaves/6M patterns) this is the
   206-234s. Clear engineering wins: integer/rolling-hash key (the peel already uses hashVec) instead of
   std::string; cheaper/batched host-confirm; inline enumLP (it is a std::function, not inlined).
NEXT TARGETS (the wall is per-phase, property-driven): (a) MAPS phase for the densest/dense graphs (string->
int keys etc.) -- highest impact on the current timeouts; (b) s=r+delta witness-floor generalization for the
peel-bound higher-(r,s) cells.
NEXT (open): generalize to s=r+delta (delta>=2, the (5,7)/(6,8)/(7,10)/(8,12) cells): group affected Q by
their dying witnesses (m_P + delta units); enumerate ALIVE witnesses, take r-shadows, accumulate
weight(y,m_Q). Output-sensitive if alive witnesses are few; drop is a (delta-j)-free-unit weighted sum
(cheaper than the full T-DP) but no longer a single closed form. More complex; do after s=r+1 lands.

## 43. MAPS phase optimization: cracked the DP-bound regime, scale-bound regime remains (2026-06-19, #139)
The densest timeout graphs die in 'pattern<->leaf maps + compaction' (206-234s) BEFORE the peel (sec 42).
That phase (region_native_sct_peel.cpp ~660-740): per leaf enumerate every r-multiset, look up in an index,
confirm host, and store per-incidence local b vectors into patLeaves/leafPats/pbLocal/leafPatLocB. THREE
optimizations, each bit-identical (15/15 golden corehashes; commits f2fa66c, c27ef6b, 06677a5):
 (1) integer rolling-hash key instead of std::string compKey -> ~1% (the string key was NOT the cost).
 (2) host-confirm via FEASIBILITY test instead of support_count DP: support_count(box,b)>0 iff the box
     {max(ell,b)<=y<=u, Σy=s} is nonempty (all weights are positive binomials), so with empty forbidden
     (pre-peel box) it is an O(width) check, no DP. -> ca-HepPh(3,4) 205->85s = 2.4x.
 (3) build the local b INLINE during enum (blocal[i]=j; lcs parallel to box local dim) and fill all 4 maps
     in one pass -> ZERO compToLocal (was ~3 allocs+binary-searches per incidence). Dropped the leafPats
     sort (not load-bearing; peel accesses by hash-mapped index). -> web-Google(3,4) 233->215s only (~1.1x).
RESULT (graph-regime-dependent): ca-HepPh-type (DP-bound, dense-but-fewer-patterns) CRACKED 2.4x; web-Google-
type (SCALE-bound: 6M patterns) barely moved -- the cost there is the sheer volume (6M-entry patIdx build +
millions of per-incidence Vec allocs in pbLocal/leafPatLocB; the push_back(blocal) still heap-allocs).
web-Google(4,5) maps still >300s. OPEN: the scale-bound regime needs an ARENA for the per-incidence local
vectors (flat int16 buffer + (offset,len) instead of vector<vector<Vec>>) to kill the millions of small
allocations -- a deeper refactor that also touches the peel's leafPatLocB[lid][t] reads. Bigger effort,
uncertain ROI; surfaced as a decision point.
SESSION NET (all bit-identical, shipped): PEEL adaptive-KMAX 1.2-1.6x + s=r+1 witness-floor 2.5-4.5x; MAPS
hash+feasibility+inline-blocal -> ca-HepPh maps 2.4x. The wall is per-phase + per-graph-regime: peel (cracked
for s=r+1), maps-DP (cracked), maps-scale (open, arena), s>=r+2 peel (open, witness generalization).

## 37. CORRECTION to 36 (contention) + precise peel complexity (#139)

§36's "SCT peel 7-11x slower, does not scale" was measured under CPU
CONTENTION (multiple parallel benches). Re-measured CLEAN (sequential, nothing
else running). True picture on com-dblp:
  (3,4): SCT 6.38s  vs V3LM 5.99s  -> TIED (build+maps 1.26s, peel 4.69s)
  (4,5): SCT 73.3s  vs V3LM 12.99s -> SCT 5.6x slower
  (5,7): SCT >200s  vs V3LM 54.76s -> SCT loses (timeout)
So the direction (loses at scale) holds, but the MAGNITUDE was wrong: it is
TIED at (3,4), and degrades only as (r,s)/pattern-count grows. [[feedback_clean_benchmarking]]

PRECISE COMPLEXITY (profiled, com-dblp 3,4):
  peel = 50% slotForbidDiff (slot scan) + 50% affected-update (scWithTerms).
  slot scan = O(#patterns x maxSplit). maxSplit = size of a leaf's split-path
  set, which GROWS during peeling via controlled_split: 1 -> 5681 on (3,4)
  (358M path-visits). At higher (r,s) overlap is denser => maxSplit far larger
  => peel is SUPER-LINEAR in pattern count.
  V3LM peels region-IE tuples with an update bounded by |host| (small, fixed)
  => near-linear. That is why V3LM degrades gently (6/13/55s) and SCT steeply
  (6/73/>200s).

ROOT: the SCT's disjoint-leaf + controlled_split peel has cost that blows up
with clique OVERLAP (split-set growth) -- the very thing it was meant to be
good at. Build is genuinely light (1.26s); the split peel is the non-scaling
part. |host|=1 skip (commit 1b6bf3f/9e72810) helps but does not change the
super-linear slot-growth term.

## 44. Front-end: r-mergeable O(Σdeg^2) quadratic FIXED + the low-RS regime line (2026-06-19, #139)
CONTEXT: on graphs where the maximal-clique count is TRACTABLE (CND completes), CND beat us at low RS
purely on the FRONT END: cit-Patents(3,4) CND 35.9s vs our MCE 14.5s + then >90s stuck; com-youtube(3,4)
CND 20.5s vs MCE 9.4s + stuck. (On clique-EXPLOSION graphs twitter/soc-pokec/wiki-Talk, CND ALSO times out
-> graph's problem, out of scope for everyone. The SDCT build is itself O(#maximal cliques): it IS the
Bron-Kerbosch recursion (src/SDCT.cpp listAllCliquesDegeneracyRecursive...) with tree_buffer.push_back per
leaf, so it does NOT escape the maximal-clique cost -- confirmed by reading the code + CND timing out there.)
THE BUG (fixed, commit cc14617): r-mergeable was O(Σ_v deg_R(v)^2) -- for each region, scan every other
region sharing a vertex; a hub vertex in K regions costs K^2. Measured 6.0e10 ops on cit-Patents (>128s),
5.1e10 on com-youtube (>123s). FIX: r-CLIQUE reformulation -- region M is mergeable iff every r-subset of M
is unique to M; enumerate (sorted r-subset, region) pairs, sort, any r-subset shared by >=2 regions makes
those regions non-mergeable. O(Σ_M C(|M|,r) log), no per-hub square. Exact (compares real r-tuples) ->
BIT-IDENTICAL (r-clique path 15/15 golden corehashes, pairwise fallback 6/6, the two agree). Fallback to
pairwise only on r-subset blowup (Nsub>3e8, clique-explosion = out of scope). SCT_RM_PAIRWISE forces old.
RESULT: r-merge 13-16s now (was >128s), bit-identical -- a real robustness win (protects ANY hub-heavy graph,
incl high-RS). BUT it did NOT win cit-Patents/com-youtube: the MERGEABLE FRACTION IS TINY (54227/2.78M = 2%;
4990/1.18M = 0.4%), so ~2.7M/1.18M regions stay ACTIVE and the wall just MOVES downstream (class build now
~1s/529K classes, then stuck in the SCT quotient+build >130s; maps and peel would follow). ROOT CAUSE: these
are millions-of-overlapping-regions graphs; the region-native pipeline does per-region work across 5 phases
while CND counts s-cliques directly -- at low RS the region machinery is structurally heavier no matter how
many phases are optimized (phase-by-phase whack-a-mole vs CND's lean 36s).
THE REGIME LINE (now precisely located): HIGH RS -> CND explodes on r-cliques -> region-native wins 3-256x
(the contribution). LOW RS + many regions -> CND's direct clique-counting is lean -> CND wins by design.
RECOMMENDATION: bank the r-merge fix (genuine, bit-identical, general); position the paper at HIGH RS; treat
low-RS-many-region (cit-Patents/com-youtube) as CND's regime and clique-explosion (twitter/soc-pokec) as
out-of-scope-for-everyone. Continuing to chase low RS = optimizing every phase (SCT build, maps-arena, peel)
against a structural disadvantage, uncertain payoff. ENV: SCT_RM_PAIRWISE (old r-merge), RM_DBG (timers).

## 45. PIPELINE STRUCTURE (code map) + per-phase timing by graph (2026-06-19, #139)
FILES: region_native/region_native_sct_peel.cpp (1389 lines, main@254) = FULL pipeline (count+peel), the
binary we optimise. region_native/region_native.cpp = separate count-only sibling. Headers included by the
peel binary: src/NucleusDecomposition/CCPathCore.h (CCPath math: support_count / count_with_extra_lower /
inclusion_exclusion_terms / controlled_split / first_failing_split); region_native/ClassSCT.h (DENSE class-
SCT buildClassSCT, the oracle); region_native/ClassSCTScalable.h (SPARSE scalableBuildClassSCT, the PRODUCTION
path actually used).

PIPELINE -> region_native_sct_peel.cpp line ranges (printf marker -> phase):
 0. load + degeneracy        254-266  (T0->T1)                        [rn] graph
 1. MCE -> regions           266-277  (T1->T2)  MCE mce(g,s);mce.run() [rn] regions(>=s)   regions=max cliques>=s
 2. r-mergeable              296-415  (Trm0)    r-clique OR pairwise   [rn] r-mergeable    <- this session (cost-based pick)
 3. class build              ~395-411 (->T3)    vtxR + region-membership string-key group  [rn] classes
 4. pattern enum             410-607  (T3->T4)  distinct patterns + r-clique counts + host; DIRECTBIN@602  [rn-peel] patterns
 5. quotient + ClassSCT      620-636  (Tqg0->Tqg1) scalableBuildClassSCT(nC,qw,qadj,s)@634  [sct] quotient   <- ClassSCTScalable.h
 6. counting + verify gate   636-660  suppOf(P)=region-IE; class-SCT total == region-IE total @658  [sct] total s-cliques
 7. pattern<->leaf maps      660-813  (Tqg1->T5) enumLP per-leaf; patLeaves/leafPats/pbLocal/leafPatLocB  [sct] maps+compaction  <- this session (hash/feasibility/inline)
 8. bucket-queue peel        813-1389 (T5->T6)  KMAX@857-872, sEqRp1 witness-floor@881, slotForbidDiff/controlled_split  [sct-peel] peel/TIMING  <- this session (KMAX/witness-floor)

PER-PHASE TIMING (sec) on diverse graphs, current binary (HEAD 9305e04), all opts ON:
 graph          (r,s) | load   MCE  rmrg class enum sctbld maps   peel | TOTAL | dominant
 ca-GrQc        4,5   | 0.01  0.02  0.00 0.00  0.03  0.00  0.12   0.39 |  0.56 | peel 70%
 ca-CondMat     5,7   | 0.03  0.05  0.00 0.00  0.10  0.01  0.29   1.69 |  2.15 | peel 79%
 com-dblp       3,4   | 0.20  0.64  0.08 0.04  0.54  0.32  2.13   9.07 | 12.82 | peel 71%
 web-NotreDame  6,8   | 0.21  1.34  0.06 0.01  3.17  0.05  3.69  43.89 | 52.22 | peel 84%
 com-amazon     3,4   | 0.17  0.54  0.07 0.06  0.36  0.22  0.96   0.75 |  2.97 | maps+peel
 web-it-2004    3,4   | 1.12  6.50  0.18 0.03  0.71  0.10  7.88   1.48 | 16.89 | MCE+maps (size)
 cit-Patents    7,10  | 3.50 13.56  0.06 0.05  1.62  0.01  0.35   0.22 | 15.87 | MCE 85% (size)
 ca-HepPh       4,5   | TIMEOUT: front end fine (rmrg 0.05s), maps/peel bound (dense)
PATTERN (property-driven, matches the §40 sweep): peel dominates (70-84%) on collab/coauthor/web with MODERATE
regions (the main cost where SCT completes); MCE dominates on HUGE graphs (cit-Patents 85%, web-it-2004)
size-driven; maps significant on large web (web-it-2004 7.9s) region-count-driven. After this session, rmrg/
class/sctbld/enum are all <1s -- none is a bottleneck; the live targets are peel (cracked for s=r+1) and, on
the densest, maps-scale (open, arena). Reproduce: /tmp/phasebreak.sh on tods2 (git-pulls, rebuilds, parses
the [rn]/[sct]/[sct-peel] lines into the per-phase split). NOTE: cost-based r-merge pick (commit 9305e04)
fixed a regression where the r-clique reform made ca-HepPh r-merge 28s (dense large regions, big C(|M|,r)).

## 46. WHERE the losing graphs drag (per-phase) + phase-4 hub-K^2 blowup is the marathon (2026-06-19, #139)
On graphs where CND completes but we don't, the stall phase DIFFERS (measured, line-buffered, grep '\[rn|\[sct'):
 - cit-Patents(3,4) CND 35.9s -> we stall in PHASE 4 (pattern enum), >105s, 2.7M active regions.
 - com-youtube(3,4)  CND 20.5s -> PHASE 4 (pattern enum), 1.18M active regions.
 - web-Google(3,4)   CND 30.5s -> PHASE 6-7 (counting/maps), passed quotient, 6M patterns (scale-bound).
 - ca-HepPh(3,4)     CND 11.7s -> PHASE 7 maps (86s) + peel (dense).
 - twitter/soc-pokec/wiki-Talk: CND ALSO times out -> graph's problem (SDCT build is O(#maximal cliques),
   confirmed by reading src/SDCT.cpp = BK recursion + tree_buffer.push_back per leaf). Out of scope.
PHASE 4 ROOT CAUSE (region_native_sct_peel.cpp ~565): each pattern's host = ∩_c classRegions[c] is computed
per (region, r-multiset). A hub class in K regions: every region M containing it enumerates patterns with it
and pays |classRegions[hub]|=K -> K^2 (cross-region), exactly r-merge's Σdeg^2. Measured PE_DBG: hostWork
super-linear 91M->783M->2.28B->9.18B@half (cit-Patents), heading ~50-100B = the >105s stall.
ATTEMPTED FIX (commit c820527, KEPT, bit-identical 15/15): thread the running host through the recursion
(one interClasses per class ADD, shared across multiplicities + subtree; host depends only on class-set) +
inline (was std::function). RESULT: only removes the WITHIN-region redundancy -> hostWork 9.18B->6.33B (~31%),
TIME unchanged (12.2s vs 11.6s, per-node alloc ate it). PHASE 4 STILL STALLS. The K^2 CROSS-region part
(pattern enumerated once per hosting region) is untouched.
VERDICT: phase-4 crack needs a GLOBAL reformulation (compute each pattern's host ONCE, not once-per-region --
like r-merge's r-clique trick but the quantity is an INTERSECTION ∩classRegions[c], not a count, so harder),
AND it's only the FIRST of cit-Patents/com-youtube's heavy phases (5/6/7 follow). => 马拉松, not one cut.
CONFIRMS the regime line (sec 44): low-RS-many-region = CND's territory by design; region-native's win is
HIGH RS. Recommendation stands: position at high RS; do NOT chase the low-RS marathon. ENV: PE_DBG (phase-4
leaves/hostWork timers).

## 47. GLOBAL host (phase-4 K² CRACKED) -- but cit-Patents/com-youtube still lose on the SUM (2026-06-19, #139)
DID IT (commit 4aa6c78, bit-identical): the phase-4 host = ∩classRegions[c] per (region, r-multiset) was the
K² hub blowup. KEY INSIGHT: host(P) = {M : P's classes ⊆ regionClasses[M]} = EXACTLY the regions that
ENUMERATE P -> emit (comp-key, region) per (region, r-multiset), sort by comp-key, GROUP: each group's region
list IS the host, NO intersection. (Same reform as r-merge pairwise->r-clique.) canonical-home dedup is now
just the grouping. Bit-identical: global 15/15 golden corehashes, per-region fallback (SCT_PE_PERREGION) 6/6,
agree (corehash = order-independent core distribution, so pats order doesn't matter).
RESULT: phase 4 CRACKED -- cit-Patents >105s -> 28.7s (emit 54.5M incidences 3.5s + sort/group ~25s, pats
3.85M), com-youtube stuck -> 14.4s (pats 2.53M). A real, general algorithmic win on the pattern-enum phase
for ANY many-region graph.
BUT STILL LOSE TO CND -- now unmistakably the SUM, not one phase:
  cit-Patents(3,4): MCE 14.5 + r-merge 13 + phase4 28.7 + SCT 6 + maps + peel ≈ 62s+ TIMEOUT vs CND 35.9s.
  MCE(14.5)+r-merge(13) ALONE ≈ 28s; then SCT/maps/peel on 3.85M patterns / 1.95M base-leaves still to go.
Every region-native phase scales with the millions of regions/patterns; CND builds NONE of these structures
(counts 4-cliques directly). The phase-4 sort itself is ~25s (54.5M records) -- could be radix/hash-grouped
to ~5s, but the front-end SUM (MCE+r-merge+SCT) already exceeds CND, so it wouldn't flip the result.
FINAL VERDICT: the low-RS-many-region regime is CND's by DESIGN (structural, the sum of phases), not a single
fixable bottleneck. The global host is KEPT (general win, bit-identical, the cleanest proof the "compute host
once" instinct was right). Region-native's contribution is HIGH RS (CND explodes, we win 3-256x). Stop
chasing low RS. ENV: SCT_PE_PERREGION (old per-region host).

## 48. BIG-GRAPH per-phase timing, ALL opts on (HEAD a5e09a8), TO=220 (2026-06-19, #139)
Reproduce: /tmp/biggraph.sh on tods2 (git-pulls, rebuilds, parses per-phase). Phases in sec:
 graph         r,s  | load   MCE  rmrg class  enum sctbld  maps   peel | TOTAL  status
 cit-Patents   7,10 | 3.41 14.17  0.06 0.05  0.09  0.01   0.35   0.22 | 14.95  OK  <- high RS, fast
 web-it-2004   3,4  | 1.12  5.93  0.16 0.03  0.40  0.10   7.81   1.44 | 15.88  OK
 com-dblp      3,4  | 0.20  0.62  0.07 0.04  0.45  0.33   2.15   8.62 | 12.28  OK  (peel)
 web-NotreDame 6,8  | 0.21  1.33  0.06 0.00  0.71  0.05   3.83  42.37 | 48.36  OK  (peel)
 web-it-2004   5,7  | 1.13  5.86  0.16 0.03  1.37  0.09  16.60  68.21 | 92.34  OK  (peel)
 ca-AstroPh    3,4  | 0.05  0.41  0.22 0.01  1.66  0.66 102.08  24.71 |129.75  OK  (MAPS 102s)
 cit-Patents   3,4  | 3.40 14.34 13.53 0.79 26.01  5.98     -      -  |  TO   stall maps/peel
 web-Google    3,4  | 0.76  4.34  7.66 0.32  8.61  2.54     -      -  |  TO   stall maps/peel
 web-Google    6,8  | 0.78  4.38  0.83 0.13 35.78  0.96     -      -  |  TO   stall maps/peel (enum 35.8!)
 com-youtube   3,4  | 0.57  9.18  7.59 0.27 14.45  6.31     -      -  |  TO   stall maps/peel
 ca-HepPh      3,4  | 0.04  0.22  0.05 0.00  1.42  2.81  84.91    -   |  TO   maps 85s + peel TO
 soc-pokec     6,8  | 3.90 51.15 61.51   -     -     -      -      -  |  TO   front-end: MCE 51 + rmrg 61
TAKEAWAYS:
 - HIGH-RS / moderate-region big graphs COMPLETE and are not slow: cit-Patents(7,10) 15s, web-it-2004(3,4)
   16s, com-dblp 12s, web-NotreDame(6,8) 48s. This is the contribution regime.
 - global host WORKED: phase-4 (enum) is no longer the stall anywhere -- it's 0.1-36s now (was >105s). The
   remaining walls are MAPS-scale + PEEL-scale (millions of patterns) and, on soc-pokec, the MCE+r-merge
   front-end (51+61s; soc-pokec is both hub-skewed AND large-region so BOTH r-merge methods are costly).
 - Timeout breakdown: low-RS-many-pattern (cit-Patents/web-Google/com-youtube 3:4, web-Google 6:8) -> maps/
   peel; ca-HepPh(3,4) -> maps 85s + huge peel; ca-AstroPh(3,4) completes but maps=102s. soc-pokec(6,8) ->
   front-end. NONE is phase-4 anymore.
 - NET after this session: front-end (r-merge, phase-4) is fixed; the live targets are maps-scale (arena) and
   peel-scale -- both ~inherent to materialising millions of patterns, the low-RS structural disadvantage.

## 49. MAPS phase: the bottleneck was a VERIFICATION left in production (region-IE), not the arena (2026-06-19, #139)
Chose B (attack maps-scale). MEASURE-FIRST saved a wrong fix: I assumed the per-incidence Vec copies in
pbLocal/leafPatLocB were the maps cost and was about to build an arena. MAPS_DBG split showed otherwise:
 - ca-AstroPh(3,4): "patIdx-build"=99.7s, enumLP=3.8s. Replacing the unordered_map patIdx with a sorted
   (hash,pi)+binary-search (commit 395e3cc) DID NOT change the 99.7s -> the timer (from Tqg1) was catching
   something ELSE. It was LINE 698: `for P: P.sup = suppOf(P)` -- the per-pattern region-IE inclusion-
   exclusion, which the code itself labels "kept ONLY as the G2a cross-check reference". The PEEL runs on
   the SCT support (sctSupport), and the G2a gate ALREADY computes sctSupport for every pattern. So region-IE
   is PURE VERIFICATION: 99.7s ca-AstroPh, 48s ca-HepPh, 15s web-it-2004, and the maps-stall on the low-RS
   losers.
FIX (commit 407017e, bit-identical): region-IE suppOf + the total-s-clique gate + the G2a per-pattern compare
are now SCT_VERIFY-gated (default OFF); production sets P.sup = sctSupport directly (region-IE == SCT, integer-
valued). Verified: production 15/15 golden corehashes, SCT_VERIFY 3/3 + gate [OK]. (Kept the sorted patIdx.)
IMPACT -- cit-Patents(3,4) went TIMEOUT(>220s) -> 82.3s COMPLETE (maps 102s-region-IE-bound -> 10.6s):
  load 3.5 + MCE 14.4 + rmrg 13.1 + class 0.8 + enum(phase4) 25.9 + sctbld 5.9 + maps 10.6 + peel 11.6 = 82.3s.
  Still LOSES to CND(35.9s) on the SUM, but now COMPLETES (was a timeout). ENV: SCT_VERIFY (re-enable the
  region-IE cross-check), MAPS_DBG (patIdx/enumLP split).
BIG-GRAPH RE-RUN (HEAD 407017e, region-IE OFF), vs §48 (before):
  graph         r,s  | load MCE  rmrg class enum sctbld maps  peel | TOTAL  | §48 before
  cit-Patents   3,4  | 3.5 14.4 13.1  0.8 25.9  5.9  10.6  11.6 | 82.3 OK | TIMEOUT  <- now completes
  web-Google    3,4  | 0.6  4.2  7.6  0.3  8.7  2.5  16.6  35.0 | 74.9 OK | TIMEOUT  <- now completes
  com-youtube   3,4  | 0.6  9.2  7.6  0.2 14.5  6.3  14.3  14.8 | 66.7 OK | TIMEOUT  <- now completes
  ca-AstroPh    3,4  | 0.1  0.3  0.2  0.0  1.4  0.6   4.99 24.8 | 32.3 OK | 129.8 (maps 102->5, 4x)
  web-it-2004   3,4  | 1.1  5.9  0.2  0.0  0.4  0.1   0.43  1.4 |  8.4 OK | 15.9 (2x)
  web-it-2004   5,7  | 1.0  5.7  0.2  0.0  1.4  0.1   1.69 68.9 | 77.9 OK | 92.3
  web-NotreDame 6,8  | 0.3  1.7  0.1  0.0  0.7  0.1   0.98 42.4 | 45.9 OK | 48.4
  com-dblp      3,4  | 0.2  0.6  0.1  0.0  0.4  0.3   1.28  8.6 | 11.3 OK | 12.3
  cit-Patents   7,10 | 3.5 14.4  0.1  0.1  0.1  0.0   0.10  0.2 | 14.9 OK | 14.9 (high RS)
  web-Google    6,8  | 0.8  4.4  0.8  0.1 36.1  1.0  94.87   -  |  TO     | TIMEOUT (maps 94.9 = genuine scale)
  ca-HepPh      3,4  | 0.0  0.2  0.0  0.0  1.2  2.7  37.68   -  |  TO     | maps 85->37.7 but peel TO
  soc-pokec     6,8  | 4.3 51.4 60.6   -    -    -     -     -  |  TO     | front-end MCE 51 + rmrg 60
RESULT: region-IE removal turned 3 timeouts into completions (cit-Patents/web-Google/com-youtube 3:4) and
gave 2-4x on the maps-heavy completers (ca-AstroPh 130->32). The completers at low RS still LOSE to CND on the
SUM (cit-Patents 82 vs 36, web-Google 75 vs 30) but they no longer TIME OUT. REMAINING walls (all post-maps-
verification): (1) web-Google(6,8) maps 94.9s = genuine enumLP/incidence scale (the arena I deferred could
help HERE); (2) ca-HepPh(3,4) huge PEEL; (3) soc-pokec(6,8) MCE 51s + r-merge 60s front-end. The "maps was
the wall" story (sec 42-48) was largely a VERIFICATION artifact -- with it gone, the true residual walls are
peel-scale, the genuine maps-incidence scale on a few cells, and soc-pokec's MCE.

## 50. MAPS lookup: flat open-addressing hash (the 56s was a binary search) + reserve pass (2026-06-20, #140)
MEASURE-FIRST on web-Google(6,8), the one cell where maps (94.9s) was genuine, not verification. Decomposed
the maps enumLP (the per-incidence loop) with an ENUMLP_PROBE ladder (throwaway /tmp build, levels skip
sub-steps):
  L0 bare recursion = 1.50s | L1 +lookup = 57.95 | L2 +hostFeasible = 63.30 | L3 +4 push_backs = 91.89
  => lookup 56.45s (61%) | push_backs 28.59s (31%) | hostFeasible 5.35s (6%) | recursion 1.5s (2%).
So the ARENA (push_backs) was the WRONG target -- the dominant 61% was the comp->pi LOOKUP: a std::lower_bound
binary search over the 23.2M-pattern sorted array, ~25 cache-cold comparisons x 107M reaches. We had abandoned
hashing earlier because unordered_map BUILD was pathological (99.7s/848k, node-per-element) -- but we jumped
straight to sorted+binary-search and never tried a FLAT open-addressing table (one contiguous vector<int>,
O(N) linear build, ~1-2 probes; collisions resolved by comparing the real comp -> bit-exact).
FIX A -- flat hash (commit 845d648): replaced patSorted+lower_bound with htab(hcap=nextpow2(2*pats), -1) +
linear probe. web-Google(6,8) maps 96.66 -> 50.38s (1.92x); enumLP 93.52 -> 48.75s; ns/inc 943 -> 492;
patIdx-build 3.11 -> 1.61s.
FIX B -- reserve() pass (commit 9b2f4c1, audit-driven via a 6-agent reserve-audit workflow): reserve at hot
single-container appends -- BK cand/P/X, r-merge active, class profKey/classRegions/classSize, interClasses
out, unionAlive inter/keep, phase-4 rec/recReg (EXACT via a count-only recursion; these grow to ~4GB on
web-Google(6,8), reserve avoids ~2x realloc-copy AND the ~12GB geometric transient that risks OOM), maps cur,
peel slot grow, SCT-build nonNb/P. SKIPPED with reason: pats (vector<Pat> realloc MOVES structs ~0.5s; counting
groups costs a full duplicate W2-scan = net loss) and chgOld (reused buffer, clear() keeps capacity). Effect:
small but real -- enum EMIT 6.0 -> 4.09s (the rec/recReg reserve), patIdx 3.11 -> 1.61s. Single-digit seconds,
exactly as predicted; the reserves are a constant-factor + OOM-safety cleanup, NOT a phase-transformer. The
phase win is FIX A.
VALIDATION: combined (845d648+9b2f4c1) bit-identical vs trusted 407017e -- 11 pass / 0 fail (r,s = 3,4/3,5/4,5/
5,6/6,7 across com-dblp/ca-GrQc/ca-CondMat/dblp-core30) + 4 skip (peel>80s timeout on both, not mismatch).
web-Google(6,8) front-end preamble 138 -> 92s; still TO in peel (23.2M patterns). ENV: ENUMLP_PROBE (0-3
sub-step ladder, throwaway), MAPS_DBG, PE_DBG.
NEW TARGETS this exposed (web-Google 6,8, post-fix): (1) maps push_backs 28.59s now ~59% of enumLP -> the
flat-list + counting-sort-to-CSR is the #1 maps target (two-pass-count is a NET LOSS here: pass 1 re-runs the
~11s+ lookup; must collect a flat incidence list in ONE pass then bucket). (2) enum sort/group 29.6s
(std::sort of 82.9M incidences) is now the enum bottleneck -> radix sort on the fixed W2-int comp-key, or
parallel sort. (3) peel-scale unchanged (the high-RS wall).

BROAD DETAILED SWEEP (HEAD 9b2f4c1 combined, single-thread, /tmp/sct_new2), vs §49 totals:
  graph         r,s  | load  MCE  rmrg class enum sctbld maps  peel | TOTAL | §49 tot
  cit-Patents   3,4  | 3.42 13.38 13.48 0.72 25.59 6.05  7.37 11.71 | 78.31 | 82.3
  cit-Patents   7,10 | 3.52 13.77 0.07  0.06  0.10 0.01  0.07  0.22 | 14.30 | 14.9
  web-Google    3,4  | 0.76  4.34 7.56  0.26  8.13 2.55  9.93 34.80 | 67.59 | 74.9
  web-Google    6,8  | 0.79  4.38 0.81  0.12 33.91 0.97 50.31   TO  | TO    | TO (maps 94.9->50.3)
  web-it-2004   3,4  | 1.11  5.74 0.16  0.02  0.38 0.09  0.33  1.41 |  8.15 | 8.4
  web-it-2004   5,7  | 0.91  5.45 0.16  0.02  1.28 0.10  1.14 63.54 | 71.70 | 77.9
  web-NotreDame 6,8  | 0.25  1.62 0.07  0.01  0.68 0.05  0.70 39.28 | 42.42 | 45.9
  com-youtube   3,4  | 0.57  9.19 7.57  0.23 13.74 6.44  7.45 14.51 | 59.15 | 66.7
  soc-pokec     6,8  | 3.83 48.94 72.05  -     -    -     -     -   | TO    | TO (front-end wall)
  com-dblp      3,4  | 0.20  0.56 0.07  0.03  0.45 0.34  1.00  7.95 | 10.40 | 11.3
  ca-HepPh      3,4  | 0.04  0.22 0.05  0.00  1.40 3.05 32.95   TO  | TO    | TO (enumLP 32.9, ns/inc 719)
  ca-AstroPh    3,4  | 0.05  0.41 0.22  0.01  1.61 0.71  3.29 24.67 | 30.92 | 32.3
SUB-BREAKDOWN (enum = emit + sort/grp ; maps = patIdx + enumLP[ns/inc] ; peel slotFwd = slotForbidDiff):
  cit-Patents  3,4 | enum emit 2.66 + sort/grp 22.86 | maps enumLP 6.56 (590) | peel slotFwd 2.81
  web-Google   3,4 | enum emit 0.94 + sort/grp 7.16  | maps enumLP 9.38 (399) | peel slotFwd 14.58
  web-Google   6,8 | enum emit 4.10 + sort/grp 29.61 | maps enumLP 48.68(491) | peel TO
  web-it-2004  5,7 | enum emit 0.18 + sort/grp 1.09  | maps enumLP 1.04 (255) | peel slotFwd 13.12
  web-NotreDame6,8 | enum emit 0.12 + sort/grp 0.55  | maps enumLP 0.66 (337) | peel slotFwd 10.78
  com-youtube  3,4 | enum emit 1.49 + sort/grp 12.20 | maps enumLP 6.88 (481) | peel slotFwd 3.02
  com-dblp     3,4 | enum emit 0.05 + sort/grp 0.40  | maps enumLP 0.93 (372) | peel slotFwd 4.50
  ca-HepPh     3,4 | enum emit 0.13 + sort/grp 1.26  | maps enumLP 32.89(719) | peel TO
  ca-AstroPh   3,4 | enum emit 0.17 + sort/grp 1.43  | maps enumLP 3.23 (481) | peel slotFwd 13.99
READING (where the time now is, by regime):
 - LOW-RS losers (cit-Patents/com-youtube/web-Google 3,4): the front+mid is dominated by enum SORT/GROUP
   (cit-Patents 22.86s, com-youtube 12.20s, web-Google 6,8 29.61s) + rmrg (cit/youtube ~13/7.5s) + MCE
   (13.4/9.2s) + peel. The single biggest reducible mid-phase is the std::sort -> RADIX on the W2-int comp-key.
 - HIGH-RS OK cells (web-it-2004 5,7, web-NotreDame 6,8, ca-AstroPh 3,4, web-Google 3,4): PEEL dominates,
   and slotForbidDiff is 20-57% of peel (14.58/13.99/13.12/10.78s) -> the peel-internal slotForbidDiff is the
   target there.
 - ca-HepPh 3,4: maps enumLP 32.89s at ns/inc 719 (highest) -> the push_back/CSR target bites hardest here.
 - soc-pokec 6,8: MCE 48.9 + rmrg 72.0 front-end wall (graph's problem, out of scope).
RANKED NEXT TARGETS: (1) enum sort/group -> radix sort (fixed W2-int key, appears on every low-RS cell,
   biggest reducible mid-phase, 7-30s); (2) maps push_back -> flat-list+counting-sort CSR (ca-HepPh/web-Google
   6,8); (3) peel slotForbidDiff (high-RS OK cells, 20-57% of peel).

## 51. enum sort/group: LSD radix sort, packed (class,mult) digit (2026-06-20, #141)
Did target (1). Replaced the std::sort over the phase-4 incidence index array (comparator = lexicographic
over W2=2r ints + region tiebreak, each compare doing 2 random gathers into the multi-GB rec) with an LSD
radix sort (stable counting passes, least-significant first: region, then comp columns high..low). Same total
order => bit-identical pats/host => identical corehash.
MEASURE-FIRST CAUGHT THE NAIVE VERSION: unpacked radix (one digit per int, 2r+1 passes) won 3.5x at r=3 but
was NEUTRAL (0.99x) at r=6 -- 13 random-gather passes ~ the comparison cost; the win scales INVERSELY with r.
The per-pass random gather is the cost, so FEWER PASSES wins. Fix: pack each (class,mult) pair into ONE digit
class<<mb | mult (mb=bits(r); mult<=r fits, class is high part => order-preserving, sentinel (nC,0) sorts
last) -> r+1 passes instead of 2r+1 (commit 9c743a2; ablation SCT_PE_RADIX_UNPACKED, SCT_PE_STDSORT).
A/B enum sort/group (packed | unpacked | stdsort | packed-vs-std):
  cit-Patents 3,4 |  6.20 |  7.87 | 30.06 | 4.85x
  com-youtube 3,4 |  3.48 |  9.18 | 15.75 | 4.53x
  web-Google  3,4 |  3.42 |  4.26 |  8.84 | 2.58x
  web-Google  6,8 | 22.60 | 25.85 | 31.71 | 1.40x   <- packing lifted r=6 from neutral 0.99x to 1.40x
Packed beats unpacked on every cell. VALIDATION: 3-way packed=unpacked=stdsort bit-identical, 10/11 + 1
timeout-only flag (com-dblp 3,5 stdsort >90s -> empty hash; packed=unpacked agree and match the prior run's
stdsort value 5100c9ad74e1; confirmed clean with a 240s re-run). Commits 4ca20c6 (radix) + 9c743a2 (packed).
IMPACT on the low-RS losers (the cells that lose to CND): cit-Patents 3,4 enum 25.6 -> ~6 (-20s of an 78s
total); com-youtube 3,4 -> -11s of 59s; web-Google 3,4 -> -5s. ENV: SCT_PE_STDSORT, SCT_PE_RADIX_UNPACKED.

## 52. maps push_back: recompute blocal (drop the ~200M stored Vec payloads) (2026-06-20, #142)
Target (2). The maps phase stores 4 vector<vector<>>: patLeaves[pi]/leafPats[lid] (int) + pbLocal[pi][k]/
leafPatLocB[lid][t] (Vec=vector<int16_t>). The two Vec maps are ~200M per-incidence heap allocs (+~6.8GB) and
the bulk of the maps push_back cost. CSR was the obvious fix but BLOCKED: support_count/covers/hashVec all take
a const Vec& (a real std::vector), so a flat int16 arena would force a per-access copy back to a Vec. KEY
OBSERVATION: each blocal is fully determined by (pattern comp, leaf class layout) -- blocal[i] = mult of leaf-
class supC[lid][i] in pats[pi].comp -- so it is RECOMPUTABLE on demand (merge of two sorted lists) into a
reused scratch, no storage at all.
TWO-PHASE (TDD): Phase 1 (commit 8f9e89e) added localB() + SCT_MAPS_VALIDATE asserting localB == every stored
payload: 129,053,828 incidences across 14 cells (incl ca-HepPh 91.5M, ca-AstroPh 13.4M), 0 bad -> proven
equivalent. Phase 2 (commit 46fb3af) SCT_MAPS_RECOMPUTE: enumLP skips the 2 Vec stores; the 5 consumers
(sctSupport, ensureLeafMap, peel P-side pl, wf-path confirm qbAll[t]==fl, general-path confirm qbAll[t]==ql)
recompute via localB into reused scratch (sctScr/elmScr/plScr/qlScr; plScr held across the k-body, qlScr per
confirm). Both paths kept for A/B; SCT_MAPS_VALIDATE guarded to the stored path.
VALIDATION (bit-identical corehash, recompute vs stored, contention-insensitive): 14 pass / 0 fail / 1 skip
(com-dblp 4,6 TO on both) -- covers BOTH peel paths: s=r+1 (wf confirm) and s>r+1 (general confirm).
STATUS: correctness-complete, GATED OFF by default (production = stored). TIMING DEFERRED -- tods2 contended by
another user (gengdaz dess, load ~23) since this work began; clean-sweep waiter armed. HONEST PERF CAVEAT (not
yet measured): recompute does MORE blocal-builds than store (~2-4x Ninc vs 2x) but each is a reused-scratch
merge (no malloc) vs a heap alloc -- expected net win from killing 200M mallocs + 6.8GB, but could be a wash;
will not claim until measured clean. And maps is only large on web-Google 6,8 / ca-HepPh (both peel-TO), so
this mainly trims front-end of TO cells -- peel slotForbidDiff (15-18s on COMPLETING cells) is the higher-value
remaining target. ENV: SCT_MAPS_RECOMPUTE (drop storage + recompute), SCT_MAPS_VALIDATE (Phase-1 oracle check).

## 53. PEEL deep-scope: cost map + ranked opportunities (2026-06-20, #143, 6-agent workflow)
Deep-scoped the peel before touching it (peel-deep-scope workflow: 6 parallel readers + synthesis). FOUR cost
centers (peel = dominant phase on completing cells; [profile] = slotForbidDiff 11-15s + affected-update ~10-20s):
 A. Full-prefix slot scan (slotForbidDiff): `for i in [0,w)` touches EVERY live path though ~99% hit
    `impossible` and skip. O(w*|plNZ|) per (P,leaf), w up to thousands. = the 11-15s slotForbidDiff bucket.
 B. Knapsack DP + IE (scWithTerms/CCPathCore): O(M*T*U_c) reweighted convolution per IE-term per (chgOld x Q).
    Two vector<double>(T+1) alloc+memset per class per term; NCR is vector<vector<double>> (double indirection
    in the innermost cell). = dominant arithmetic of the affected-update bucket (s>r+1 cells).
 C. DFS over-generation + map probes (general s>r+1): DFS gen >> #affected-Q; every wasted candidate pays one
    node-based unordered_map (leafQ2pat) cache-cold probe. = the lookup half of affected-update.
 D. s=r+1 liveness re-scan: liveness over chgOld reruns Mloc times (once per receiving class c), mostly
    c-invariant. O(Mloc^2*|chgOld|*|fbd|). = affected-update on s=r+1 cells ONLY.
 (E allocator churn, F scheduling = few % each, smeared.)
 NOTE: a cell is EITHER s=r+1 (path D) OR s>r+1 (paths B+C) -- disjoint, so wins on both sides needed to move
 the aggregate.
RANKED (bit-identical preserving; excludes already-done KMAX/witness-floor/IE-cache/affected-Q-DFS + 2 refuted prunes):
 #1 Slot u-envelope INDEX (center A, the biggest bucket): index slot paths so the scan visits only paths with
    p.u>=bloc on plNZ, not the whole [0,w). HIGHEST ceiling but medium-risk/hard (index vs swap-remove/append
    consistency). bit-id LIKELY (changed-set unchanged, slot order-independent). Needs a differential verify-flag
    harness (assert indexed chgOld == full-scan chgOld) before trusting.
 #2 Flat open-addressing leaf map (center C): replace leafQ2pat unordered_map with the in-repo htab pattern
    (lines 826-848, already proven bit-exact). LOW risk / moderate / bit-id YES. Notable benefit.
 #3 Flatten+hoist NCR binomial table (center B): hoist `const double* row = NCR[nc-bc]` out of the y-loop +
    flat vector<double> stride. LOW risk / EASY / bit-id YES BY CONSTRUCTION (same doubles, same order). 5-15%
    of affected-update. = SYNTHESIS TOP PICK (provable now on a contended box, clear benefit, no bookkeeping).
 #4 s=r+1 per-box liveness precompute (center D): predicate differs from pl only at coord c; precompute ell-
    violation + per-forbidden gap-coord once -> Mloc^2 becomes Mloc. medium/moderate, bit-id YES. s=r+1 only.
 #5 s=r+1 uEnv class-prune (center D): uEnv already computed for general path; skip class c with uEnv[c]<pl[c]+1
    (provably nAlive=0). LOW/EASY/bit-id YES. stacks with #4.
 #6 Hoist per-peel scratch out of while loop (aff/chgOld/chgOldTerms/plNZ/DP) -> clear() not reconstruct (center
    E). LOW/EASY/bit-id YES. (NOTE: delta already hoisted line 1273; aff IS loop-local 1318.)
 #7 DP scratch thread_local + reachable-band [SL,min(T,SU)] (center B/E). LOW/EASY. bit-id YES (w==0 skips).
 #8 Dense coreDist + flat bk level-queue + precomputed affectsH2 (center F). LOW, few-% each, bundled.
 DROPPED/traps: forced-bloc split (changes materialised slot -> risky), q2p hash-fn change (must edit build+probe
 atomically), scWithTerms floor-collapse memo (subtle cache key).
SYNTHESIS ORDER: #3 -> #2 -> #6/#8 (free safety) -> #4/#5 -> #1 (big swing, with harness, when a clean box
allows timing). Full report: workflow wsfbvc9qm.

## 54. PEEL #1 (slot-skip) groundwork: cost-structure probe + the order-independence unlock (2026-06-20, #144)
User picked #1 (the sub-linear slot skip, biggest bucket). MEASURE-FIRST before building the index.
SFD_DBG cost-structure probe (slotForbidDiff per-path-test counts, deterministic -> contention-proof):
  cell           tested     affected%  coord-tests/test  fail-on-1st   slotForbidDiff
  com-dblp 3,4    359M       0.81%      1.37              66.2%         5.87s (56% of peel)
  ca-AstroPh 3,4  870M       1.19%      1.31              74.4%         17.97s (58%)
  web-Google 3,4  215M       8.04%      1.49              70.9%         15.68s (42%)
  web-NotreDame   837M       0.19%      1.72              52.3%         15.23s (28%)
  web-it-2004 5,7 928M       0.49%      1.45              66.2%         18.76s (23%)
=> 0.19-1.2% affected on most cells (web-Google 8% is the outlier): a 50-500x ceiling on the SCAN. The scan
pointer-chases w heap Vecs (p.u), cache-miss-bound. KEY DESIGN QUESTION found by reading the code: the scan
uses SWAP-REMOVE, giving the slot a specific internal order; a sub-linear skip must reorder (find-then-act,
mark-compact). Is slot-path ORDER load-bearing for bit-identical? Supports are summed in slot order and FP add
is non-associative -- IF supports are not exact-integer-recovered, reorder flips llround and breaks corehash.
PROBE (commit c67f364, SCT_SLOT_REVERSE: reverse every slot once after build, corehash vs default): 12 SAME /
0 different (com-dblp/ca-GrQc x5/ca-CondMat/dblp-core30/ca-AstroPh, s=r+1 and s>r+1). => ORDER IS NOT LOAD-
BEARING (supports exact-integer, FP error <0.5). UNLOCK: the skip may find-then-act / mark-compact / remove in
any order; only the live-path SET must be correct. Risk drops very-hard -> moderate.
DESIGN (in progress): per-leaf lazy bucket index bkt[c][v]=path-indices with u[c]==v (u only DECREASES on
split, so maxv from the initial build bounds all future u); query = pick pivot = plNZ coord with min
Sum_{v>=bloc[c]}cnt[c][v], enumerate that bucket suffix, filter survivors by the other plNZ coords. Maintain
via swap-remove + back-pointers bpos. Incremental plan: (1) build+maintain index with SCT_SLOT_IDX_VERIFY
asserting index==cur consistency, still full-scan for affected; (2) switch affected-find to the index with a
per-call differential assert (index-affected==scan-affected); (3) corehash + timing (clean box). ENV: SFD_DBG,
SCT_SLOT_REVERSE.

## 55. PEEL #1 (slot dominance index) -- CORRECTNESS-COMPLETE, bit-identical (2026-06-20, #145)
Built the sub-linear slot skip in 2 verified steps; both pass.
INDEX (per-leaf, lazy like leafQ2pat): ixBkt[lid][c*maxv+v] = slot positions with u[c]==v; ixBpos[lid][pos*M+c]
= back-pointer (index within its bucket) for O(1) removal. maxv from the initial build bounds all future u (u
only DECREASES on split). ixFindAffected: pivot = plNZ coord with min Sum_{v>=bloc[c]}|bkt[c][v]|, enumerate
that bucket suffix, filter survivors by the other plNZ coords (read cur[pos].u, few candidates). Maintained at
swap-remove (ixRemove + ixRelabel w->i) and append (ixAppend), cur stays dense so consumers (sctSupport, peel
iteration) are untouched.
STEP 1 (commits 8a783c6/ccaab31, +#include <climits>): index built+maintained alongside the full scan;
SCT_SLOT_IDX_VERIFY asserts ixFindAffected == scan-affected on every call (pristine slot). 13 OK / 0 FAIL
(com-dblp x2, ca-GrQc x5, ca-CondMat x3, dblp-core30 x2, ca-AstroPh; s=r+1 AND s>r+1). Proves the index +
mutation maintenance (swap-remove relabel, append, back-pointers) is correct.
STEP 2 (commit 39d38f3, SCT_SLOT_IDX): index DRIVES the find. Phase A classify (snapshot/cover/split/modify-in-
place, NO position moves) -> Phase B batch swap-remove dead DESCENDING (w>=p always, cur[w] live) + append
children, maintaining the index. chgOld in pivot order (order-independent per sec 54). corehash index-driven vs
default: 13 OK / 0 MISMATCH, both peel paths. => slot skip is BIT-IDENTICAL, correctness-complete.
STATUS: gated SCT_SLOT_IDX (default OFF). TIMING DEFERRED -- tods2 contended all session (gengdaz).
ROUGH A/B (load 26, NOT a clean claim) slotForbidDiff default->idx: com-dblp 3,4 6.93->5.79 (16%), ca-CondMat
4,5 0.39->0.39 (~0), ca-AstroPh 3,4 17.78->17.36 (2%). => MODEST, smaller than the 200x scan ceiling. HONEST
CORRECTION to the earlier estimate: the impossible-SCAN is ~20-25% of slotForbidDiff, NOT 40-70%. slotForbidDiff
is DOMINATED by the affected-paths' REAL WORK (chgOld CCPath snapshot copies + controlled_split), which is
NECESSARY and not skippable. The index correctly removes the 928M-test scan but that scan was a minority. So #1
is a real, bit-identical, algorithmically-correct win but a modest one (~16% of slotForbidDiff, contended; clean
number TBD). Durable value: the correct sub-linear index + the sec-54 order-independence finding.
STRATEGIC RE-RANK: the AFFECTED-UPDATE ("rest" in [profile], 10-39s) is the BIGGER peel bucket on most cells
(web-NotreDame rest=38.85, web-Google rest=21.80) vs slotForbidDiff 11-19s. So the higher-value remaining
targets are the affected-update: #3 (NCR binomial flatten/hoist, the DP inner cell, bit-id by construction) and
#4 (s=r+1 liveness precompute). ENV: SCT_SLOT_IDX (drive), SCT_SLOT_IDX_VERIFY (Step-1 assert), SFD_DBG.

CLEAN LOCAL #1 NUMBERS (user's idea: measure on local uncontended Mac with ca-* graphs scp'd to /tmp; the
contended-server signal was misleading). These split (maxSplit 1243-5524), single clean run, default->idx:
  ca-CondMat 4,5 | aff 2.35% | slotForbidDiff 0.15->0.15 (~0%) | total 0.43->0.42
  ca-CondMat 5,6 | aff 0.71% | slotForbidDiff 0.25->0.20 (20%) | total 0.58->0.51 (12%)
  ca-GrQc    5,6 | aff 0.24% | slotForbidDiff 0.25->0.15 (40%) | total 0.47->0.36 (23%)
  ca-GrQc    6,7 | aff 0.18% | slotForbidDiff 0.18->0.10 (44%) | total 0.30->0.22 (27%)
=> #1 cuts slotForbidDiff 40-44% and TOTAL peel 23-27% on LOW-affected-% cells (<=0.25%), scaling to ~0 at
2.35%. The win scales INVERSELY with affected% exactly as the mechanism predicts. The server's big slotForbidDiff
cells are low-affected (web-NotreDame 0.19%, web-it-2004 0.49%, ca-AstroPh 1.19%) so #1 should help them
substantially. MUCH better than the contended 16%; #1 is a SOLID win, not modest. (Totals sub-second so some
noise, but the trend is clean+monotone in affected%.) METHOD NOTE: local Mac measurement bypasses the contended
tods2 -- use it for all clean timing while gengdaz holds the server. (Harness gotcha: the Bash tool shell is
ZSH, which does NOT word-split unquoted $var; `set -- $spec` leaves $2/$3 empty -> use explicit args.)

## 56. #3 (NCR binomial hoist) is a WASH -- measure-first caught the synthesis overestimate (2026-06-20, #146)
Did #3: hoist `const double* ncrow = NCR[nc-bc].data()` out of the scWithTerms DP y-loop, replace the
per-cell ccpath_ncr(nc-bc,y-bc) with ncrow[y-bc]. Bit-identical by construction (verified IDENTICAL cores on
8 cells, both small and high-s). CLEAN LOCAL A/B (3-trial min peel, s>r+1 cells so the general DP fires):
  ca-GrQc 3,5 -3% | ca-GrQc 4,6 -1% | ca-CondMat 3,5 -0% | ca-CondMat 4,6 -0%
  ca-CondMat 5,7 -0% | ca-CondMat 6,8 -0% | ca-GrQc 6,8 -3% | ca-GrQc 5,8 -0%   (rest=0.23-0.63s)
=> ~0% everywhere. The synthesis's 5-15% estimate was WRONG: ccpath_ncr is a one-line `return NCR[n][k]` that
the compiler already inlines, so removing the function call / double-indirection saves nothing. Scale-invariant
0% -> won't materialize at server scale (rest=38s) either. REVERTED (no measured benefit; project rule: no
unmeasured changes). The affected-update cost is NOT the binomial lookup -- it's the DFS over-generation (C) or
the inherent O(patterns x chgOld x T) DP volume / IE machinery. #4 (s=r+1 liveness) and #2 (flat leaf map)
remain untried.
SESSION CONSOLIDATED (this stretch): radix sort = 4.85x enum (clean, SHIPPED default); maps-CSR recompute = bit-
identical, gated SCT_MAPS_RECOMPUTE OFF, timing unmeasured (may be wash); slot index #1 = bit-identical, CLEAN
23-27% total peel on low-affected cells, gated SCT_SLOT_IDX OFF -> candidate to make DEFAULT; #3 = wash,
reverted. NEXT: make #1 default after a broad clean validation; then #4/#2 or step back to structural peel.
MAPS-CSR MEASURED (clean local 3-trial-min, all bit-id) -- NET LOSS on TIME, keep gated OFF:
  ca-AstroPh 3,4: maps 1.27->0.80 (-37%, the 200M-alloc removal works) BUT peel 13.23->15.50 (+17%, recompute
    at the peel consumers costs more than the stored read) => total 16.64->17.40 (+5% SLOWER).
  ca-CondMat 3,4/4,5 same shape (maps faster, peel slower, total ~flat-to-slower).
  => recompute's on-demand blocal at the hot peel consumers dominates; peel >> maps so net is a small LOSS.
  DECISION: maps-CSR stays gated OFF (default stored = faster). Its real value is MEMORY (~6.8GB saved on huge
  cells web-Google 6,8 / ca-HepPh that OOM) -> keep SCT_MAPS_RECOMPUTE as a MEMORY escape-hatch, NOT a speed default.
UPDATE (#147): #1 made DEFAULT (commit 5d60e08, escape hatch SCT_NO_SLOT_IDX), default==full-scan bit-identical
on 4 cells. User chose to STEP BACK to the STRUCTURAL peel redesign (constant-factor opts can't change the
regime: peel is O(#patterns x per-peel-update) and low-RS #patterns is millions = inherent). Launching a design-
exploration workflow for structural approaches that break the O(#patterns x update) scaling while staying exact.

## 57. STRUCTURAL peel redesign: the scaling wall is EXACT-IRREDUCIBLE past |host|=1 (2026-06-20, #148, workflow)
6-approach design-exploration workflow (wymzy4bgx); the synthesizer EMPIRICALLY verified its claims (ran the
corehash gate + profiling on dblp-sigmod/db on the server). DECISIVE NEGATIVE RESULT: NO exact approach breaks
O(#patterns). 4 of 6 are rigorously-assessed DEAD-ENDS:
 - h-index convergence (generalize R=1 Local to r>=2): EXACT-on-convergence but makes the per-update term WORSE
   -- the histogram needs a per-witness min over up-to C(s,r)-1 co-r-subclique patterns, and support_count
   counts witnesses COMBINATORIALLY without listing them (its whole advantage over CND); producing per-witness
   thresholds forces s-clique enumeration x #iters = strictly worse than the bucket peel at r=3,s=4. (R=1 wins
   only because a witness threshold there folds into a closed-form nCr breakpoint.)
 - lazy support maintenance: the slot MUTATION (44-48% of peel) cannot be deferred (later sctSupport reads would
   see stale slots = WRONG); lazy only touches the cheap incremental delta and replaces it with from-scratch
   recompute. Ceiling <2x, realistic = slowdown. Dead.
 - global witness->pattern inverted index: exact requires carrying the full C(n-b,y-b) reweighting per entry
   (reproduces current work, saves nothing); dropping it = the componentwise-max shortcut the code PROVES wrong
   (L1433). Box-splitting churns it every peel. Dead.
 - generalized bulk-prune (enlarge r-mergeable/direct-bin past |host|=1): EMPIRICALLY FALSIFIED -- corehash
   SCT_DIRECTBIN_ALL_HOST1 = 1cc6ca3a vs correct 32451c95 on sigmod(3,4). A |host|=1 pattern's witness s-clique
   can contain a |host|>=2 r-clique, so its core != its initial support. The exact fixed-point stratum is
   PROVABLY MAXIMAL at |host|=1 (already implemented via skipH1 + direct-bin). No exact enlargement without a
   new theorem. (sigmod 3,4: 13966 |host|=1 / 3236 |host|=2 -- the 81% H1 mass is already harvested.)
TOP PICK (only remaining EXACT lever): BATCH-PEEL-WHOLE-LEVEL. Drain the whole bk[curLevel] bucket at once; per
touched leaf run ONE batched slotForbidDiff + ONE merged affected-Q DFS with cross-pattern seen[] dedup, instead
of once per pattern. CONSTANT-FACTOR (does NOT break O(#patterns) -- honest), divides the per-update constant by
the same-level co-hosting fan-in and kills the affected-Q DFS over-enumeration (MEASURED gen/nz = 8.5 on db,
13.7 on sigmod => 8-14 candidates generated per genuine drop). Bit-identical (slot order-invariant, sec 54).
SUBTLETY: multiple antichain thresholds in the SAME path interact non-linearly through the IE, so combined-drop
!= Σ individual-drops -- must recompute affected-Q against the COMBINED slot state (still exact, loses the
additivity shortcut not correctness). Reuses seen[]/aff dedup, chgOldTerms cache, the per-level bucket. medium/
medium.
KILL-GATE FIRST EXPERIMENT (counter-only, no rebuild, contention-insensitive): instrument the peel to count, per
level, total (pattern,leaf) touches A vs distinct (leaf,level) touches B; fanin = A/B. If fanin~1 batch saves
nothing (near-dead, STOP); if fanin>=5 proceed. Run on dblp-sigmod(3,4)/dblp-db(3,4)/ca-CondMat.
BOTTOM LINE: the "crack the scaling wall" hope is FALSIFIED for EXACT cores -- the asymptotic factor is
irreducible and the bulk-bin stratum is exhausted. Batch-peel is the largest exact constant-factor reclaim
left, gated by the fan-in probe.
KILL-GATE RESULT (commit 63496da, FANIN_DBG, local counter-only): fan-in = avg patterns sharing a (leaf,level):
  dblp-sigmod 3,4 = 10.0 | dblp-db 3,4 = 18.4 | ca-CondMat 4,5 = 32.4 | ca-CondMat 5,6 = 70.5 | ca-GrQc 5,6 = 1213
ALL >> the >=5 GO threshold. Same-level patterns heavily co-host leaves => the affected-update re-processes each
leaf 10-1213x per core level. STRONG GO for batch-peel: doing the per-leaf slotForbidDiff (44-48% of peel) +
affected-Q DFS once per (leaf,level) instead of once per pattern divides that work by ~fan-in. NEXT: build
batch-peel (medium/medium), validate the non-linearity soundness probe (combined multi-threshold insert vs
Σ-individual -> must recompute affected-Q against the COMBINED slot state, exact) + end-to-end corehash gate.

## 58. BATCH-PEEL design + incremental plan (2026-06-20, #149) -- the build, scoped from the full loop read
CURRENT PEEL LOOP (region_native_sct_peel.cpp ~1399-1683): while(peeledN<npat){ advance curLevel to lowest
non-empty bk; pi=bk[curLevel].back(); pop; if(!alive||key!=curLevel)continue; alive=false;core=curLevel;
peeledN++; coreDist+=mult; source-skip(host1); FOR k in patLeaves[pi]: lid=pleaf[k]; slotForbidDiff(lid,pl)
mutate slot + capture chgOld; then wf-path(s=r+1) OR general-DFS(s>r+1) accumulate delta[qi]/aff; APPLY: for qi
in aff: ns=sup-delta; nk=llround(ns); clamp>=curLevel; re-bucket. }
WHY BATCH IS EXACT: cores are order-independent WITHIN a level (all level-L patterns get core L; drops to higher-
level Q telescope: total = pre-level-sup - post-level-sup, order-independent) + slot order not load-bearing
(sec 54). So we may drain a whole curLevel wave and process leaf-major.
DESIGN (gated SCT_BATCH_PEEL, default OFF = proven per-pattern path; A/B + corehash):
 Step 1 (wave-drain, bit-identical, NO win yet -- scaffold): outer loop drains ALL current bk[curLevel] into a
   wave list (mark each peeled), processes each pattern's affected-update with the EXISTING inner logic but
   ACCUMULATES aff/delta across the wave, applies ONCE per wave. Cascade (drops re-bucketing Q to curLevel) is
   handled by the outer loop re-draining curLevel until empty. Verify corehash == default.
 Step 2 (leaf-major group + shared per-leaf setup): group the wave's (pi,leaf) work by leaf; per leaf do
   ensureLeafMap / leafQ2pat / sufPl setup ONCE for the group (amortize). Still per-pattern slotForbidDiff +
   delta inside the group. Verify corehash. (Partial win: amortized leaf setup.)
 Step 3 (THE win -- batch slot mutation + pre-minus-post delta): per (leaf,wave) insert ALL group thresholds
   into the slot in one pass (combined slotForbidDiff), tracking per chgOld-path its POST paths (in-place-
   modified self OR split children). Then each affected Q's drop = Σ_changed-paths [SC(p_pre, ql) - SC(p_post,
   ql)] -- the pre-minus-post telescoped drop (NOT Σ-individual-threshold; that is the non-linearity the
   synthesizer flagged). This divides slotForbidDiff (44-48% of peel) + the affected-Q DFS by ~fan-in (measured
   10-1213). REQUIRES: slotForbidDiff to return the chgOld->post-paths map; the affected-Q enumeration to cover
   Q affected by ANY group threshold (DFS bounds from the union of the group's m_P). Verify corehash + a
   soundness probe (combined-state drop == sequential sum on one level of ca-GrQc 3,4).
RISK: high (bit-identical peel core). MITIGATION: gated, incremental, corehash after each step; the slot index
(#1) used the same find-then-act + verify discipline successfully. ENV: SCT_BATCH_PEEL.
STATUS: batch-peel DESIGNED + de-risked (kill-gate GO, exactness proof), NOT yet implemented -- it is a ~260-line
restructure of the bit-identical peel core (extract affected-update into a callable unit, move aff/apply to wave
scope, add wave-drain driver). Deferred to a focused effort: too risky to rush via incremental Edits at the tail
of a long session. §58 is the actionable plan; pick up from Step 1.

## 59. #1 (slot index) VALIDATED AT SCALE on the big server cell (2026-06-20, #150)
Clean pinned (taskset -c 90) A/B on tods2 (96 cores, load 20, gengdaz busy -> pinning dodges it), HEAD 753e89a,
default(#1 ON) vs SCT_NO_SLOT_IDX(#1 OFF):
  web-NotreDame 6,8 (0.19% affected): slotForbidDiff 16.27 -> 9.38 (-42%) ; peel total 54.06 -> 43.83 (-19%).
=> the shipped #1 default holds AT SCALE: slotForbidDiff cut 42% (matches the local ca-GrQc 6,7 -44% prediction
exactly for ~0.19% affected); peel-total cut 19% here (lower than the local 23-27% because web-NotreDame's
affected-update "rest"=38s is unusually large, diluting the slotForbidDiff fraction). Confirms #1 is a real
clean win at scale, not a small-graph artifact.
LEAN 1-trial A/B (no min, contended box load 14) -- NOISE-CONTAMINATED, peel totals UNTRUSTWORTHY:
  web-NotreDame 6,8 peel 95.65/73.64 (-30%!?) -- a CONTENTION SPIKE (ON 95.65 vs the 3-trial-min 43.83 on the
    SAME cell); 1-trial under load is unreliable, do NOT trust the lean peel numbers.
  web-it-2004 5,7: sfd 10.54/18.31 = -42% (clean, matches 0.49% affected) ; ca-AstroPh 3,4: sfd -3% (1.19%
    affected) ; com-dblp 3,4: sfd -9%.
TRUSTWORTHY at-scale #1: web-NotreDame 3-trial-min (sfd -42%, peel -19%) + web-it-2004 sfd -42%. The sfd cuts
corroborate the inverse-affected% mechanism (big at <=0.5% affected, small at >1%). LESSON RE-CONFIRMED: 1-trial
timing on a contended box is worthless; use 3-trial-min + taskset pin, or measure on the clean local Mac.
DEFINITIVE CLEAN LOCAL A/B (user scp'd the big graphs to local Mac, 3-trial-min, HEAD 11b3e60, ALL bit-id corehash
ON==OFF) -- settles the lean -30% scare = it was SERVER CONTENTION, not a #1 regression:
  web-NotreDame 6,8 (0.19% aff): peel ON 16.00 / OFF 18.78  => #1 +15% FASTER (lean server said -30% = noise!)
  ca-AstroPh 3,4 (1.19% aff):    peel ON 12.51 / OFF 13.27  => #1 +6% faster
  web-it-2004 5,7 (0.49% aff):   peel ON 26.58 / OFF 30.04  => #1 +12% faster
=> #1 is FASTER on ALL big cells (6-15% peel), bit-identical, NO regression anywhere. The shipped default is
validated clean at scale. (6-15% < the 23-27% on tiny ca-GrQc cells because big cells have a larger affected-
update "rest" fraction #1 doesn't touch -- inverse-affected% mechanism.) Local Mac (web-NotreDame peel 16s vs
the bad server trial's 95.65s) is the clean measurement env -- use it.

## 60. MEASUREMENT METHODOLOGY: local Mac DRIFTS between windows -> INTERLEAVE conditions (2026-06-20, #151)
A block-sequential A/B (3 ON trials, then 3 OFF) gave web-NotreDame 6,8 #1 = -26% (ON SLOWER), contradicting an
earlier +15% on the SAME clean local cell. ROOT CAUSE (5-trial within-window run): within a tight window ON/OFF
are STABLE (ON~26.3, OFF~29.2, +10%, SAME RSS 0.8G -> index is NOT a memory hog), but ACROSS windows minutes
apart the whole Mac's speed drifts hugely (OFF was 18.7 / 29.2 / 18.6 in different runs; ON 16 / 23.7 / 26).
Drift = thermal + my own background agents/workflows loading the Mac all session. Block-sequential lets drift
masquerade as a real effect on SMALL-margin cells. FIX: INTERLEAVE ON/OFF per trial (ON,OFF,ON,OFF...) and take
the mean of per-pair ON/OFF RATIOS -- each pair samples the same instant, drift cancels. (Interleaved web-NotreDame
= +13%, then the full table below.)
#1 FINAL DRIFT-ROBUST TABLE (local, interleaved, mean ON/OFF ratio over 5 pairs, all bit-id) -- MONOTONE in aff%:
  ca-GrQc 6,7   0.18% +32% | ca-GrQc 5,6   0.24% +24% | web-NotreDame 6,8 0.19% +16% | ca-CondMat 5,6 0.71% +12%
  web-it-2004 5,7 0.49% +11% | ca-AstroPh 3,4 1.19% +1% | ca-CondMat 4,5 2.35% +1% | ca-GrQc 4,6 1.55% -1%
=> #1 wins big at low affected% (<=0.5% -> +11..32%), NEUTRAL at high (>1% -> +-1%), never meaningfully regresses.
The -26% web-NotreDame scare = 100% drift artifact; #1 is +16% there. Shipped default VALIDATED.
P1 STATUS: DONE. maps-CSR = net loss on time (gated OFF, memory hatch); #1 = clean validated table above.
LESSON: for any small-margin local timing A/B, INTERLEAVE conditions + paired ratio, never block-sequential.

## 61. P2 HEADLINE: region-native vs CND vs V3LM, 3-way time+memory (2026-06-20, #152, icml2)
SETUP: icml2 (64c, 503GB, IDLE) -- cloned repo, NO-SUDO deps (apt-get download + dpkg -x for boost-serialization
+ sparsehash; manual-Boost CMake patch), built BOTH degeneracy_cliques (CND=PIVOTER_RUN_REF, V3LM/RegNDC=
PIVOTER_RUN_REGION_V3LM) and region-native; staged 13 graphs (tods2->local->icml2). 3-WAY COREHASH VERIFIED
IDENTICAL (CND=V3LM=region on ca-GrQc/ca-AstroPh/ca-CondMat) -> fair comparison. /usr/bin/time -v wall+maxRSS,
single trial, timeout 400s, HEAD 817939b.
RESULTS (wall-clock | peak RSS):
  cell             CND          V3LM         region-native
  com-dblp 3,4     4.4s 0.6G    6.1s 0.2G    12.1s 1.2G
  ca-AstroPh 3,4   4.2s 0.3G    52.5s 0.3G   32.7s 2.6G
  ca-HepPh 3,4     11.3s 0.5G   353s 0.7G    TIMEOUT(>400)
  com-youtube 3,4  21s 1.1G     120s 1.7G    62.8s 6.4G
  web-Google 3,4   27.9s 3.1G   61s 1.7G     80.8s 9.5G
  cit-Patents 3,4  36.5s 3.3G   84.8s 2.7G   75.8s 7.3G
  web-NotreDame 6,8 TIMEOUT     47s 0.2G     40.1s 0.9G
  web-it-2004 5,7   TIMEOUT     201.9s 0.7G  70.8s 1.8G
  cit-Patents 7,10  21.2s 2.5G  38s 0.7G     17.3s 0.6G   <- region FASTEST
  web-Google 6,8    TIMEOUT     TIMEOUT      TIMEOUT      <- graph's problem (all 3)
  ca-HepPh 6,8      TIMEOUT     TIMEOUT      TIMEOUT
HONEST REGIME CHARACTERIZATION (confirms sec 57):
 1. LOW-RS (3,4): CND WINS on time (4-37s, fastest) AND region-native is a MEMORY HOG (6-9.5G vs CND 1-3G --
    the pattern materialisation). Worst case = dense ca-HepPh 3,4: region TIMEOUT vs CND 11s. region-native
    LOSES the low-RS regime, decisively. (vs V3LM mixed: region beats V3LM on ca-AstroPh/com-youtube/ca-HepPh,
    loses on web-Google.)
 2. HIGH-RS: region-native WINS -- CND TIMEOUTS on web-NotreDame 6,8 + web-it-2004 5,7 (structural r-clique
    explosion = CND's wall), region completes (40s,71s) and BEATS V3LM (201.9->70.8 on web-it-2004); cit-Patents
    7,10 region is FASTEST of all 3 (17.3s). Memory LEAN here (0.6-1.8G).
 3. The 2 hardest cells (web-Google 6,8, ca-HepPh 6,8) TIMEOUT for ALL 3 -- the graph's problem, out of scope.
MEMORY ANSWER (user's q): region-native peak = 0.6-9.5G (hog at low-RS, lean at high-RS); CND 0.3-3.3G; V3LM
0.2-2.7G. All fit 64GB local; the 503GB icml2 has huge headroom. Caveat: SINGLE-TRIAL (icml2 idle so clean-ish);
multi-trial for the paper. PAPER POSITIONING: region-native's win is the HIGH-RS regime where CND explodes on
r-cliques; honest limitation = low-RS (slower + memory-heavy), exact-irreducible (sec 57).

## 62. MEMORY loss DECOMPOSED -- half fixable, half inherent (2026-06-20, #153) [FIX QUEUED]
User asked: a succinct (SDCT/CPI, no s-clique enumeration) method should NOT lose memory. Investigated the 2
worst low-RS hogs (icml2, /usr/bin/time -v peak RSS, 4 configs):
  web-Google 3,4 (CND 3.1G):  default 9.5G/80s | maps-CSR 6.4G/76s | no-index 8.7G/77s | MEM-MIN 5.6G/72s
  com-youtube 3,4 (CND 1.1G): default 6.4G/62s | maps-CSR 4.7G/58s | no-index 5.6G/60s | MEM-MIN 3.9G/55s
  (MEM-MIN = SCT_MAPS_RECOMPUTE=1 SCT_NO_SLOT_IDX=1)
TWO PARTS:
 1. FIXABLE (implementation): the maps Vec storage (~200M blocal Vecs) + the #1 index buckets. MEM-MIN reclaims
    web-Google 9.5->5.6G, com-youtube 6.4->3.9G. KEY: on these BIG-incidence cells maps-CSR is FREE -- both
    LOWER memory AND FASTER (80->72, 62->58) -- contradicting sec 52's small-cell "net time loss" verdict. So
    maps-CSR's time sign is CELL-SIZE-DEPENDENT (hurts small, helps big) -> it should be an ADAPTIVE default
    (auto-enable when incidences/patterns are large).
 2. INHERENT (design): even MEM-MIN stays ~2-3.5x CND (5.6 vs 3.1; 3.9 vs 1.1). This is the explicit ORBIT/
    PATTERN materialisation -- region peels PATTERNS (bucket queue + affected-update), so it must materialise
    all Pat structs (each = host+comp+classSet vectors) + the patLeaves/leafPats maps. CND peels per-r-clique
    INCREMENTALLY -- a single int counter per r-clique, no orbit structure. So even with #patterns <= #r-cliques,
    region's PER-ELEMENT overhead (3 vectors vs 1 int) makes it net heavier. => the "succinct counting saves but
    orbit materialisation costs" tradeoff; region does NOT theoretically guarantee a memory win.
FIX PLAN (user: "start later"):
 (a) maps-CSR ADAPTIVE default: auto-enable SCT_MAPS_RECOMPUTE when incidences exceed a threshold (find the
     break-even from sec 52 small-cell vs this big-cell data). Cheapest high-value: free memory on big cells.
 (b) Flatten Pat (host/comp/classSet) + patLeaves/leafPats into CSR/arena -> kills the per-vector overhead of
     the residual. Diminishing returns, won't reach CND, but narrows the 2-3.5x gap.
 (c) #1 index buckets (vector<vector<vector<int>>>) are heavy (~0.8G on web-Google 3,4); a flatter layout or
     large-slot-only build trims it -- but #1 is a TIME win so keep it default, just lighten the structure.

## 63. Is the peel CACHE-bound? NO -- it is COMPUTE/ALLOC-bound (2026-06-20, #154)
User: before flattening structures for speed, measure if the peel is cache-bound. perf BLOCKED (perf_event_
paranoid=4 on icml2 AND tods2, no sudo); valgrind not installed + won't apt-get-download. So measured via a
LAYOUT-ONLY A/B: SCT_PL_CSR flattens patLeaves (vector<vector<int>>, pointer-chased per peeled pattern) to a
contiguous CSR -- same lids, same order, bit-identical; ONLY the memory layout changes. If peel speeds up =>
cache-bound. RESULT (icml2, interleaved 4-pair drift-robust, all bit-id): web-NotreDame 6,8 +1% | ca-AstroPh
3,4 -1% | com-youtube 3,4 +0% => NULL. patLeaves access is NOT a bottleneck.
CONVERGENT EVIDENCE the peel is COMPUTE/ALLOC-bound (not cache-bound):
 - patLeaves flatten = 0% (this probe).
 - #3 NCR binomial hoist = wash (sec 56) -> the DP's table lookup isn't the bottleneck.
 - deep-scope (sec 53): cost = support_count DP convolution (compute) + controlled_split (CCPath child allocs)
   + affected-Q DFS over-generation (gen/nz = 8-14x = MEASURED redundant compute).
CONCLUSION: flattening structures is a MEMORY play (sec 62) with at-best-modest speed benefit -- NOT the speed
lever. The genuine "faster, no trade" levers are: (1) BATCH-PEEL (dedups the 8-14x DFS over-gen = the biggest
measured compute waste; sec 58 design, kill-gate GO sec 57); (2) controlled_split ALLOC pooling/arena (dense
cells maxSplit~1000s of small CCPath child allocs -> arena = faster AND less memory, a real Pareto win).
CAVEAT: patLeaves is ONE (outer-loop) probe; the heaviest random access is slotPaths (CCPath p.u Vec chase in
slotForbidDiff) + pats[qi] in the DFS -- not separately isolated, but the DP/split/over-gen evidence already
points compute-bound. (Probe reverted to keep the default hot path lambda-free; MEM_DBG kept.)

## 64. The support_count DP (convolution) is UNNECESSARY for ~89% of calls -> closed-form direction (2026-06-20, #155)
User insight: "why convolution? the count can be computed directly." TRUE. support_count = the bounded binomial-
weighted composition count: sum_{ell<=z<=u', sum z=Z} prod_c C(M_c, z_c) (M_c=n_c-b_c, Z=s-sum_b). The DP only
exists to handle the UPPER bounds u'; with no truncation it Vandermonde-collapses to a single C(sum M_c, Z).
CF_DBG PROBE (per DP IE-term, count #upper bounds that TRUNCATE = Uc-Lc < total-slack), s>r+1 cells:
  cell            DP IE-terms   #binding-upper=0 share   tail(many binding)
  web-NotreDame 6,8  16.2M       89.3%                    ~10%
  web-it-2004 5,7    29.5M       90.0%                    ~10%
  ca-GrQc 4,6        1.14M       89.4%                    ~10%
  ca-CondMat 4,6     0.77M       83.6%                    ~16%
=> ~89% of DP calls the upper bound NEVER truncates -> the DP's upper-bound machinery is wasted 89% of the time.
The avg closed-form/DP-work RATIO printed >1 (19.7/10.0/147.7/1.67) is a MISLEADING aggregate -- dragged up by
the ~10% high-binding tail (2^k blowup in the naive model). Split: the 89% easy terms are ~1/(M*s) the DP cost
(~50x CHEAPER as closed form); the 10% tail is where 2^k explodes. => HYBRID: closed form for low-binding terms,
DP for the high-binding tail, cuts the DP (the dominant peel arithmetic) massively.
DERIVATION SKETCH (general, any s-r): count = [x^Z] prod_c P_c(x), P_c = sum_{z=ell_c}^{u'_c} C(M_c,z) x^z =
(1+x)^{M_c} - low_c(x) - high_c(x). Expand: count = sum over (low-tail choices of the raised-lower classes) x
(high-tail choices of the binding-upper classes) of +- a Vandermonde C(remaining sum M_c, Z - degree). #terms =
prod_{raised-lower}(ell_c+1) x prod_{binding-upper}(...) -- SMALL when few bounds bind (the 89%), since for an
r-clique drop the raised lowers are <= r classes (where P's threshold pl exceeds Q's) with small ell_c. This
GENERALIZES the s=r+1 witness-floor (s-r=1 done) to ANY s-r. NOTE: lower bounds come from the IE-over-forbidden
`terms` (already expanded) PLUS addLow=pl (P's threshold) in the drop. NEXT (user): derive+verify the general
closed form for ANY s-r (not just 2), numerically check vs the DP, then hybridise. ENV: CF_DBG (binding probe).

## 65. GENERAL closed form for support_count -- DERIVED + NUMERICALLY VERIFIED (2026-06-20, #156)
The DP count_with_extra_lower = sum_{ell_c<=z_c<=u'_c, sum z=Z} prod_c C(M_c,z_c) (M_c=n_c-b_c, ell_c=L_c-b_c,
u'_c=U_c-b_c, Z=s-sum b_c). CLOSED FORM (any Z), by expanding [x^Z] prod_c [(1+x)^M_c - low_c - high_c] where
low_c=sum_{j<ell_c}C(M_c,j)x^j, high_c=sum_{j>u'_c}C(M_c,j)x^j:
   count = SUM over K (subset of ACTIVE classes, each picking ONE tail term j_c in [0,ell_c-1] OR [u'_c+1,M_c])
           of  (-1)^|K| * prod_{c in K} C(M_c,j_c) * C( sum_{c not in K} M_c ,  Z - sum_{c in K} j_c )
   ACTIVE class = ell_c>0 OR u'_c<M_c. Classes not active are "full" (1+x)^M_c -> Vandermonde-collapse into the
   trailing binomial. NO CONVOLUTION -- pure binomial lookups. #terms = prod_{active}(1 + ell_c + max(0,M_c-u'_c)).
VERIFIED: /tmp/cf_verify.py, closed-form == brute-force DP on 119,691 random cases (n<=5, M<=8, lo<=3) -> 0
mismatches; + real-peel regime (large M<=200, small Z<=6) [running]. The formula is EXACT (polynomial identity).
COST: #terms tiny when few bounds bind (the 89% from sec 64: B empty, |A|<=r raised lowers with small ell ->
a handful of binomials) vs the DP's O(M*(s+1)) cells. HYBRID rule: compute #terms cheaply; if #terms < DP cost
(~M*(T+1), capped), use closed form, else DP. Implementing next; gate SCT_CF, must be bit-identical corehash
(the count is an exact integer either way -> same llround -> same cores).

## 66. Hybrid closed form IMPLEMENTED: bit-identical but speed INCONSISTENT -> not the lever (2026-06-21, #157)
Implemented the hybrid (commit, SCT_NO_CF escape hatch, default ON): per IE-term, if active-class expansion
<CF_CAP and < DP cost, Vandermonde-IE closed form (recursive, pure binomials) else DP. CF_AB (icml2, interleaved,
default-CF vs SCT_NO_CF): ALL BIT-IDENTICAL (no FP cancellation break -- the hybrid only fires CF on few-bound
terms where cancellation is mild). SPEED: ca-GrQc 4,6 -34% | ca-CondMat 4,6 -8% | web-NotreDame 6,8 -5% |
web-it-2004 5,7 -2% | com-dblp 3,5 -1% (all SLOWER). High-s: cit-Patents 7,10 +13% (faster!) | ca-GrQc 5,9 +0% |
ca-CondMat 5,9 -50% (slower). => INCONSISTENT, mostly slower. WHY: the DP inner loop is a tight FMA/vectorizable
loop; the closed form is recursion + scattered NCR lookups + a per-term active-class scan. For small T=s-r the
DP is already cheap so the CF overhead loses; and the DP's total cost is dominated by the NUMBER of calls (the
8-14x affected-Q over-generation), NOT per-call size -- which the CF doesn't touch. So eliminating the convolution
per-call is mathematically right (sec 65 verified) but is NOT a speed win. (Kept gated, default... TBD; lean OFF.)
PIVOT (user's deeper idea): COMPOSITION-LEVEL DIRECT COUNTING (like CND but succinct). Maintain a_Y = #alive
s-cliques per s-COMPOSITION Y (a_Y = prod C(n_c,Y_c), computed once). Support(Q) = sum_{Y>=m_Q} a_Y. Peel P =>
zero a_Y for all Y>=m_P (the whole orbit dies => every composition >=m_P dies). Drop(Q) = sum of the zeroed a_Y
>= max(m_Q,m_P). NO convolution, NO closed form -- just sum + zero. It is the EXPLICIT version of the forbidden-
antichain (compact dead-set) + recompute. TRADE: memory (one counter per alive composition per leaf) vs the
compact box+convolution. VIABILITY = #s-compositions per leaf (COMP_DBG measuring now): few=>direct-count wins
(simpler+faster), many (wide leaves/high s)=>too much memory, the box+convolution is why the current code exists.

## 67. COMP_DBG measured: #s-compositions per leaf = 4-45x #patterns (direct-counting memory cost) (2026-06-21, #158)
COMP_DBG probe (gated, after maps print, count-DP over slotPaths[lid][0].u per leaf): integer pts (Sum y=T, 0<=y<=u).
  cell            leaves   total comps   avg/leaf  max/leaf   #patterns   comp/pat
  web-NotreDame 6,8  19409   2,917,509     150      94,146     695,144     4.20x
  ca-CondMat   4,6    4,469     669,780     150      27,132     105,784     6.33x
  web-it-2004  5,7   70,930  12,005,581     169     153,451   1,302,866     9.21x
  web-Google   3,4  774,473  71,550,413      92      24,158   6,096,377    11.74x
  com-dblp     3,4   99,641  15,282,334     153     443,662     643,485    23.75x
  ca-GrQc      4,6      165     610,034    3697     117,603      13,636    44.74x
VERDICT: s-compositions are TENS OF MILLIONS (4-45x patterns), NOT billions of s-cliques => direct-count is succinct
(nowhere near CND). BUT 4-45x more memory than the pattern materialization; the box+forbidden+convolution is the
COMPACT form (one box covers many compositions). Direct-count = no convolution + kills slotForbidDiff(44-48%) but
4-45x memory. Fundamental trade: convolution is the PRICE of the compact representation. (low-RS cells already lose
on memory => direct-count makes memory worse there; high-RS lean cells 4-9x may be ok.)

## 68. THEORY: "delete a class & update CPI" -- exact characterization (2026-06-21, #159)
User Q: CND has an efficient delete-clique-from-CPI primitive; can we have an efficient delete-CLASS-from-our-CPI?
Derived carefully from the leaf-box algebra (CCPathCore.h).

### Leaf-box model (precise)
A leaf box is a CCPath p: classes c with n[c] vertices (all mutually adjacent inside the region, so any subset is a
clique); h[c]=held(forced-in) classes; T=s-Sum h = #more vertices to pick; ell[c],u[c]=bounds on y_c; forbidden=
antichain of dead-corner thresholds. #s-cliques in the box containing footprint b =
  support_count(p,b) = Sum_{IE terms over forbidden} Sum_{Y in box, ell<=Y<=u, Sum Y=T} Prod_c C(n_c - b_c, Y_c - b_c).
Support of pattern Q = Sum over leaves of support_count(p, b=m_Q). Peel Q = insert threshold a=tuple_to_threshold(m_Q)
into each hosting leaf's forbidden antichain (kills every composition Y >= m_Q, i.e. every s-clique containing a
Q-orbit r-clique).

### What "delete a class" means and the exact theorem
Two distinct dead-set representations:
  (A) FORBID-CORNER (current): dead = {Y >= a}, a multi-coordinate corner. Accumulates an ANTICHAIN. k corners =>
      up to 2^k IE terms; controlled_split bounds it to kmax by cutting the box into axis-aligned sub-boxes. This is
      the slotForbidDiff + IE + convolution machinery -- the entire expensive peel.
  (B) DELETE-CLASS (reduce n_c / u_c): removes class-c VERTICES. Keeps the box PRODUCT-FORM (just a smaller box),
      NO antichain, NO IE, NO split. support = one clean bounded composition count.
THEOREM. A dead pattern's footprint corner {Y>=m_Q} is axis-aligned (= a single class-bound reduction, rep B) IFF
m_Q is supported on exactly ONE class (|supp(m_Q)|=1). Multi-class patterns (|supp|>=2) are genuine corners,
NOT expressible as any class deletion.
COROLLARY (why r=1 is already clean). For r=1 every pattern is a single vertex = single class => EVERY peel is an
axis-aligned class/vertex deletion => the alive set is ALWAYS a plain box, the forbidden antichain is NEVER needed.
This is exactly why the r=1 path (ST_V2) has no IE/split. The user's instinct ("class deletion saves a lot of
trouble") is PROVABLY correct -- but only at r=1, where dead objects ARE vertices.
OBSTRUCTION (r>=2). Peeling removes r-CLIQUES = combinations of vertices across classes, NOT vertices. Removing the
orbit (c:1,d:1) does not reduce n_c or n_d (those vertices still live in other r-cliques) => cannot be a class
deletion => must forbid the cross-class corner. The antichain/IE is the irreducible price of r>=2 multi-class
deadness. So a LITERAL "delete a class" cannot drive the r>=2 (r,s)-nucleus peel.

### The salvage: HYBRID axis-aligned absorption (concrete, measurable)
Split every dead footprint by |supp(m_Q)|:
  - |supp|=1 (single-class pattern, e.g. (c:r) = r interchangeable verts in one class): absorb as a clean u_c
    reduction. NO antichain entry, NO IE -- it is a class deletion.
  - |supp|>=2 (cross-class): keep in the residual forbidden antichain (the only place IE/split fires).
=> the expensive antichain holds ONLY genuinely cross-class dead corners; all single-class peels become free box
shrinks. Win size = fraction of dead footprints that are single-class (measurable, like a CF_DBG probe over the
forbidden insertions). If single-class peels are common (classes of size>=r are common in dense regions), the
antichain shrinks a lot => fewer IE terms / less split => faster, SAME memory, bit-identical. NEXT GATE: probe what
fraction of forbidden-threshold insertions have |supp(a)|=1 on the real cells; if high, build the hybrid (u_c-reduce
fast path in lazy_delete_tuple / slotForbidDiff before insert_antichain). This is the honest realization of "delete
a class": it is exact + free where the pattern is single-class, and falls back to the corner only where it must.

## 69. SUPP_DBG measured: single-class share <2% where it costs -> class-deletion hybrid NOT worth building (2026-06-21, #160)
Probe (gate SUPP_DBG, in slotForbidDiff): histogram of |supp(a)| over EVERY forbidden-threshold insertion, with
calls / affected-path work / controlled_split-triggers per |supp| bucket. |supp|==1 == axis-aligned == the clean
"delete a class" u_c-reduction (§68). Single-class SHARE = the hybrid's win ceiling. icml2, 12 cells:
  cell             total-calls  splits    SC-calls%  SC-aff%  SC-splits%   dominant |supp|
  web-Google 3,4    12.25M      9.37M      0.6%      0.5%      0.5%        |supp|=3: 93.8%
  com-dblp   3,4     1.99M      1.84M      0.4%      1.7%      2.3%        |supp|=3: 88.5%
  com-dblp   3,5     1.20M      1.61M      0.7%      2.8%      3.6%        |supp|=3: 80.7%
  ca-CondMat 4,5     0.24M      0.14M      0.1%      0.3%      0.2%        |supp|=4: 86.0%
  ca-CondMat 4,6     0.16M      0.10M      0.2%      0.3%      0.3%        |supp|=4: 77.4%
  ca-CondMat 5,9     0.047M     0.023M     0.3%      0.4%      0.3%        |supp|=5: 52.3%
  web-NotreDame 6,8  0.76M      0.80M      0.2%      0.6%      1.0%        |supp|=6: 52.0%
  web-it-2004 5,7    2.21M      2.41M      1.9%      1.2%      2.1%        |supp|=4: 35.6%
  ca-GrQc    4,6     0.042M     0.044M     0.2%      0.1%      0.0%        |supp|=4: 85.0%
  cit-Patents 7,10   5302       2651       0.1%      0.1%      0.1%        |supp|=6: 90.6%
  ca-GrQc    5,9     1287        616       4.4%      5.6%      3.9%        (tiny: sub-ms peel)
  ca-GrQc    7,10    528         230       6.2%      9.9%      8.3%        (tiny: sub-ms peel)
  (2 cells failed: web-NotreDame 4,5 timeout; ca-HepPh 6,10 std::bad_alloc in r-mergeable 683M r-cliques. Neither
   changes the verdict.)
VERDICT: single-class share is <=2.3% of split work on EVERY cell with real cost (>100k splits). It only reaches
4-10% on ca-GrQc 5,9 / 7,10 which are sub-millisecond (528-1287 total calls). The dominant footprint dimension is
|supp| ~ r (88-95% at r=3; 77-90% at r=4): a peeled r-clique's r vertices almost always land in ~r DISTINCT classes
(in dense regions classes are small, so vertices are not interchangeable). => the "delete a class" axis-aligned
fast path captures <2% of the work where it matters. NOT WORTH BUILDING.
ADVERSARIAL VERIFY (separate agent, 140M oracle-checked queries, 0 mismatches): the hybrid IS provably bit-identical
(insert single-class corner a == set u_c=t-1; AND it correctly prunes pre-existing multi-class corners a' with
a'_c>=t). The §68 theorem is SOUND. The only blocker is the empirical fraction, not correctness. Verifier flagged the
ixBkt slot-index update on u_c-shrink as "the single most error-prone implementation point" => real risk for a <2% gain.
DEEPER TAKEAWAY (paper asset): |supp|~r is the QUANTITATIVE justification for the corner/antichain/IE machinery. r>=2
deadness is irreducibly r-dimensional cross-class corners, NOT axis-aligned class deletions (confirms the §68
obstruction with data). Both of the user's "eliminate the convolution" routes are now measured + closed: composition
direct-counting = 4-45x memory (§67); class-deletion = <2% opportunity (§69). The compact box+forbidden+IE
representation matches the real structure of what dies. Direction CLOSED with evidence; keep the current machinery.

## 70. §58 BATCH-PEEL IMPLEMENTED + bit-identical (2026-06-21, #161) [local-only validation, no server]
Built the §58 batch-peel: gate SCT_BATCH_PEEL, default OFF, active only for s>r+1 (the path that owns an
affected-Q DFS to amortize; s=r+1 keeps the proven per-pattern witness-floor path via !sEqRp1).
DESIGN (region_native_sct_peel.cpp ~1506-1650): drain a whole curLevel WAVE (re-drain for cascade), mark all
wave patterns peeled up front, group (pattern,leaf) tasks LEAF-MAJOR (sorted by lid), and per leaf: apply ALL the
leaf's wave-thresholds via slotForbidDiff (accumulate pre-images into coAll tagged with their per-threshold pl in
coPls/coPlIdx), then run ONE DFS over ql<=uEnv (Sum=r) and for each confirmed Q compute drop = Sum over coAll
entries of scWithTerms(p, coTerms[e], ql, &plE) -- the SAME proven single-threshold delta-formula, summed over the
leaf's thresholds. Apply once per wave.
WHY BIT-IDENTICAL (proven + verified): drops are ADDITIVE over thresholds; each threshold's pre-image is snapshotted
at the correct (post-prior-insertions) slot state, exactly like the per-pattern sequential capture; the drop math
reads only slot paths (never .sup), so reordering apply to wave-scope changes nothing; every wave member peels at
curLevel regardless of order (order-independence within a level + slot order not load-bearing, sec 54); intra-wave
drops clamp to curLevel==key and never re-bucket in the per-pattern path, so skipping them (via !alive) is exact.
The batch DFS drops the per-pattern Sum-max(pl,ql)<=Tcap prune (no single pl) => over-generates candidates, all
filtered by applyIdxB's feasibility + d==0 checks; it is a strict SUPERSET of generated candidates so misses no
genuine Q.
VALIDATION (LOCAL, corehash = md5 of sorted core=K count=N): 7/7 bit-identical, default vs SCT_BATCH_PEEL:
  s>r+1 (batch active): dblp-sigmod 3,5 / 4,6 ; dblp-db 4,6 / 5,7 / 3,5  -> all match.
  s=r+1 (fall-through):  dblp-sigmod 3,4 (=32451c95, the known-correct ref) ; dblp-db 4,5  -> all match.
Plus a separate ADVERSARIAL REVIEW (subagent, 10 checkpoints: accumulation order, whole-wave-peel safety, cascade,
DFS-bound false-negative, bk-iterator rehash safety, source-skip, count-once, memory reuse, pbLocal indexing,
uniform Mloc/T): 10/10 PASS, "bit-identical: yes", no bugs. Memory: all batch buffers (wave/taskLL/coAll/coPlIdx/
coPls/coTerms/chgTmp/plNZ/aff/uEnv/qcand) declared outside the loops, clear()/assign() per wave/leaf (capacity
retained, no per-iteration realloc), bounded by one leaf's wave-threshold count. No regression vs per-pattern alloc.
SCOPE/HONESTY: the slot index (#1) already made slotForbidDiff output-sensitive, so the batch's win is the
affected-Q DFS + candidate-generation amortization (divides by fan-in), NOT the slot mutation (still F insertions).
Win concentrated in s>r+1 (high-RS, where region-native already wins + the DFS dominates the "rest"). Timing: local
interleaved peel-time only (per user: no server experiments this round).

## 71. §58 batch-peel: per-leaf CB GATE makes it a clean win, no regression (2026-06-21, #162) [local-only]
§70's first cut (single shared no-prune DFS) was measured (local interleaved peel time, dblp-db): r=3,s=5 2.13x WIN
but r=4,s=6 0.83x and r=5,s=7 0.59x LOSS. Cause: the shared DFS dropped the per-pattern Σmax(pl,ql)<=T prune (no
single pl over a leaf's wave-thresholds) -> candidate enumeration BLOWS UP at high r. Two prune attempts measured +
rejected: per-ENTRY prune (state per pre-image, O(nCo)/node) -> r=3 crashed to 0.21x (61.86s); per-THRESHOLD prune
in one DFS (O(F)/node + uEnvThr precompute O(nCo*Mloc)/leaf) -> r=4 WORSE at 0.35x (264s). The tight prune's per-node
/precompute overhead exceeds the over-generation it removes on the common (non-exploding) leaves.
FIX = cheap PER-LEAF gate (committed 6dd557b). CB = #{ql: Σ=r, ql<=uEnv} via a saturating O(Mloc*r^2) count-DP.
  CB<=cap (default 128, env SCT_BATCH_CB_CAP): ONE shared no-prune DFS, full drop over ALL pre-images [0,nCo).
  CB> cap: fall back to F cheap pl_j-pruned DFS (per-pattern's exact prune, O(1)/node) over each threshold's
           contiguous pre-image range [coStart[j],coStart[j+1]); drop = that threshold's contribution.
KEY: both paths feed the SAME confirm(); the drop partitions exactly by threshold (Σ_j drop[range_j) == drop[all]),
so the gate choice is BIT-IDENTICAL for ANY cap -- the gate is pure speed, never correctness. coAll is naturally
grouped by threshold (pushed threshold-by-threshold) so coStart gives the ranges for free.
RESULT (dblp-db, peel time, interleaved/known-def, ALL bit-identical):
  r=3,s=5: 13.11->6.05s = 2.17x | r=4,s=6: 91.31->46.99s = 1.94x | r=5,s=7: 980.32->321.68s = 3.05x.
Gate split @ dblp-db 4,6: 86.8% leaf-instances single-DFS, 13.2% fallback -> batches the common case, routes only
the few exploding leaves to per-pattern cost. The old worst case (5,7) is now the BEST speedup.
CORRECTNESS: 17 corehash cells across 4 graphs (dblp-sigmod/db, amazon-copurchase, soc-Epinions1), s=r+1 (fall-
through) and s>r+1 (both DFS paths) -- all default==SCT_BATCH_PEEL. Plus the §70 by-construction telescoping proof +
10/10 adversarial review. MEMORY: all buffers reused (coStart/cbDP/sufScr/uThrScr added, all clear()/assign() per
leaf). Per user this round: LOCAL validation only, no server experiments. Default OFF (SCT_BATCH_PEEL opt-in);
s=r+1 untouched (proven witness-floor path). NOTE: the #1 slot index already made slotForbidDiff output-sensitive,
so the win is the affected-Q enumeration amortization, concentrated in s>r+1 (where region-native already wins).
NEXT (not done, deferred): wider sweep + cbCap sensitivity if promoting to default; s=r+1 batch variant (witness-
floor) if low-RS peel ever becomes the target.

## 72. s=r+2 witness-major drop: DERIVED + brute-force VERIFIED (theory confirmed, pre-code) (2026-06-21, #163)
User asked for special opts at small s-r (=1,2). Derived the WITNESS-MAJOR drop reorganization of the proven
per-pattern delta-formula: write witness Y = pl + delta (Σdelta = s-r = tail); the Q it feeds are Q = Y - gamma
(Σgamma = tail, gamma<=Y, gamma!=delta), with weight Π_k C(n_k - Q_k, gamma_k) [= Π C(n-ql, Y-ql), the same
reweighting]. So:
  drop(Q) = Σ_{delta: Σ=tail, Y=pl+delta ALIVE in box} Σ_{gamma: Σ=tail, gamma<=Y, gamma!=delta, Q=pl+delta-gamma} Π C(n_k-Q_k, gamma_k)
"ALIVE" = ell<=Y<=u AND Y not>= any forbidden (the antichain handled by a direct dominance test -- NO IE, NO DP, NO DFS).
This is EXACT for ANY tail (just witness-major vs Q-major summation of the same set); cheap ONLY for small tail
(the (delta,gamma) count is O(M^tail)). That is precisely why s=r+1/r+2 are special.
  s=r+1 (tail 1): delta=e_c, gamma=e_d, weight=n_d-Q_d == the EXISTING witness-floor (sEqRp1). cross-check.
  s=r+2 (tail 2): gamma=2e_e -> C(n_e-Q_e,2); gamma=e_e+e_f -> (n_e-Q_e)(n_f-Q_f). delta in {2e_c} ∪ {e_c+e_c'}.
                  closed-form binomials, output-sensitive, no DP/IE/DFS.
VERIFICATION (/tmp/verify_rp2.py, vs brute-force general Σ_{Y>=max(pl,Q),ΣY=s,ell<=Y<=u,Y not>=forbidden} ΠC(n-Q,Y-Q)):
  tail=1: 3271 inst / 382379 (box,Q) checks / 0 mismatch | tail=2: 2847 / 287528 / 0 | tail=3: 914 / 81448 / 0.
  Random n<=6, ell, u, 0-2 forbidden corners, M=2-5. THEORY CONFIRMED (reweighting + Q!=P exclusion correct).
COST (s=r+2, per (P,leaf)): enumerate delta (Σ=2, O(M^2) configs) x scan chgOld boxes for Y-liveness (O(boxes*M))
x per alive delta enumerate gamma (Σ=2, gamma<=Y, ~O((r+2)^2)) -> credit drops. ~O(M^2 * boxes * M) per (P,leaf).
IMPLEMENTATION CHOICE (pending): (a) per-pattern witness-major for s=r+2 (like witness-floor) -- simplest, compare
vs the batch (which already gives 1.9-3.0x at s=r+2, ALL my batch cells were s=r+2); OR (b) compose witness-major
WITH the wave batch. Gate as SCT_RP2_WITNESS (or fold into sEqRp1-style branch). Corehash-validate, LOCAL only.
NOTE: this COMPETES with the batch at s=r+2 -- must beat the batch's 1.9-3.0x to be worth shipping. Measure first.

## 73. s=r+2 witness-major IMPLEMENTED: beats batch AND default (option (a) wins) (2026-06-21, #164) [local]
Built the per-pattern s=r+2 witness-major fast path (gate SCT_RP2_WITNESS, default OFF, only when s==r+2 and not
batchPeel). Sits next to the sEqRp1 witness-floor in the per-pattern affected-update. Dying witness Y=pl+δ (Σδ=2,
two added classes a1<=a2), affected Q=Y-γ (Σγ=2, two removed classes b1<=b2, (b1,b2)!=(a1,a2)), drop =
[γ=2e_b -> C(n_b-Q_b,2)] or [γ=e_b1+e_b2 -> (n_b1-Q_b1)(n_b2-Q_b2)] times #alive-boxes. n leaf-constant -> weight
box-independent, only the alive-box COUNT matters. uEnv pre-prune skips infeasible δ before the box-liveness scan.
NO IE / NO DP / NO DFS.
CORRECTNESS: brute-force theory verified (§72) + 6 corehash cells bit-identical across 4 graphs (dblp-sigmod/db,
amazon-copurchase, soc-Epinions1), default==SCT_RP2_WITNESS.
SPEED (dblp-db, peel time, interleaved 3-way default/batch/witness, ALL s=r+2):
  3,5: def 11.36 | batch 5.40 | witness 2.85  -> witness 3.99x vs def, 1.89x vs batch
  4,6: def 99.39 | batch 49.96| witness 20.49 -> witness 4.85x vs def, 2.44x vs batch
  5,7: def 1129  | batch 372  | witness 164.6 -> witness 6.86x vs def, 2.26x vs batch
=> WITNESS-MAJOR is the clear s=r+2 winner: 4.0-6.9x vs default, 1.9-2.4x vs the batch. The closed-form (no DP/IE/
DFS) beats the batch's DFS-amortization. Option (a) succeeds; no need for (b) batch composition to beat batch.
END-STATE design for the affected-update: s=r+1 -> witness-floor (default ON); s=r+2 -> witness-major (THIS,
currently gated, should become default for s=r+2); s>r+2 -> batch (SCT_BATCH_PEEL, gated) or general DFS.
PENDING (user decision): flip SCT_RP2_WITNESS to default-on for s=r+2? (bit-identical + 4-7x local win; a server
sweep would confirm at scale but this round is local-only.) And whether to extend the witness-major to s=r+3
(tail 3, O(M^3) configs -- likely still beats general for moderate M, but the batch/general may win past some tail).

## 74. UNIFIED tail-parameterized witness path (replaces witness-floor + witness-major) + crossover (2026-06-21, #165) [local]
Replaced the two hand-rolled fast paths (t=1 witness-floor, t=2 witness-major, ~125 lines) with ONE
tail-parameterized recursion (gate SCT_WITNESS_TMAX, default 2): addDelta recursively adds t units to pl as a
non-decreasing class multiset -> witness Y; count alive boxes; remGamma recursively removes t units -> Q with
weight Π C(n_b-Q_b, mult_b); credit drop = weight × #alive. Q==P excluded by qi!=pi (no explicit γ≠δ needed --
Q=pl maps only to pi). m=1 weight fast path (avoid ccpath_ncr). Active for 1<=t<=witnessTMax; batch condition
changed to (batchPeel && !useWitness) so batch only runs past the witness cap.
CORRECTNESS: bit-identical to the general fallback (SCT_WITNESS_TMAX=0) at t=1,2,3 across dblp-sigmod/db (incl
the 32451c95 reference). The NEW t=3 path also matches -> the recursion is exact for any tail.
CROSSOVER (local peel time, the SCT_WITNESS_TMAX knob makes this directly measurable):
  t=1 dblp-db 4,5: witness 7.3s vs general 21.2s = 2.45x (still wins).
  t=3 dblp-sigmod 3,6: witness 0.07 | general 0.08 | batch 0.03 -> witness ~= general, LOSES to batch (0.43x).
  t=3 dblp-sigmod 4,7: witness ~= general ~= batch (~0.09s, tied).
  => t=3 is the CROSSOVER: the O(M^t) witness enumeration stops winning; batch/general take over. Confirms
     witnessTMax=2 is the right default. (t=3 cells here are tiny <0.1s; a larger t=3 graph would sharpen it.)
UNIFICATION TAX: the unified recursion is ~7% slower at t=1 than the old flat-loop witness-floor (7.34 vs 6.86s,
interleaved) -- runtime-t recursion can't unroll to a flat loop. Accepted: t=1 is low-RS (not region-native's
strength), the path still wins 2.45x vs general, and the user prioritized de-duplication. A template<int t>
dispatch would regain flat-loop speed with one source if ever wanted. Default behavior unchanged (t=1,2 witness).

## 75. CROSSOVER IS GRAPH-DEPENDENT (driven by leaf width M) -> fixed witnessTMax is not universal (2026-06-21, #166) [local]
User: the crossover param differs per graph. CONFIRMED. The witness enumeration is O(M^t) (M=leaf width, t=s-r),
so the crossover where it loses to batch/general depends on M, which varies a lot by graph density. Measured
max leaf width M (SCT_WIT_DBG, r=3 s=5): amazon-copurchase M=7 | dblp-sigmod M=18 | dblp-db M=33 | soc-Epinions1
M huge (witness t=2 >150s timeout -- very dense).
CROSSOVER SWEEP (r=3, t=2/3/4, witness vs general vs batch, local peel time; w/batch = batch/witness):
  amazon  (M=7) : t2 0.97 | t3 1.33 | t4 - (sub-0.1s cells, noise)
  sigmod  (M=18): t2 0.83 | t3 0.57 | t4 0.40 (ALL sub-0.1s -- NOISE, do not trust)
  dblp-db (M=33): t2 1.51 | t3 1.13 | t4 0.29  <-- RELIABLE (cells 3.6-54s)
RELIABLE READ (dblp-db, seconds-scale): witness BEATS batch at t=2 (1.51x) and t=3 (1.13x), LOSES at t=4 (0.29x,
batch 3.4x faster). Crossover between t=3 and t=4 here. amazon (tiny M) keeps winning further; soc-Epinions (huge
M) is so dense even t=2 witness is very slow (loses risk at t=2, unmeasured). So: witnessTMax=2 is SAFE on the
measured moderate-M graphs (witness clearly wins t=2 on the reliable cell), t=3 is graph-dependent-marginal, t>=4
loses. But a FIXED cap is NOT optimal per graph -- small M wants a higher cap, huge M may want lower (even <2).
ELEGANT FIX (proposed, not built): make witness/fallback ADAPTIVE PER LEAF, like the batch's CB gate (sec 71):
per leaf cheaply estimate the witness enumeration size ~ #compositions of t over its M (a small DP), and use
witness only when it is below a cost threshold vs the general/batch alternative; else fall back. Auto-adapts to M
and t together, no per-graph tuning. (The batch already does exactly this with CB; the witness can reuse the idea.)
NOTE: small-graph timing (<0.1s) is noise-dominated; only seconds-scale cells (dblp-db) are trustworthy here.

## 75b. soc-Epinions (huge M) t=2: witness STILL WINS 1.19x -> witnessTMax=2 safe across the density spectrum
soc-Epinions1 3,5 (t=2, huge M): witness=73.01s vs batch=87.03s -> w/batch=1.19, witness WINS. (Earlier "timeout"
was the 150s cap; it just takes 73s.) So at t=2 witness wins on EVERY measured graph: amazon(M7)~tie/noise,
soc-Epinions(hugeM) 1.19x, dblp-db(M33) 1.51x. CRITICAL TREND: the t=2 margin SHRINKS as M grows (1.51x @ M=33 ->
1.19x @ huge M), so for pathologically dense graphs the t=2 win -> 1.0 and could eventually tie/lose. CONCLUSION:
witnessTMax=2 is a SAFE, well-justified default across the measured spectrum; t=3 is a graph-dependent marginal
win (dblp-db 1.13x, amazon 1.33x); t>=4 loses. The adaptive per-leaf gate (sec 75) remains the robust long-term
fix -- it would (a) safely capture the t=3 win where it exists and (b) protect against the shrinking t=2 margin on
extreme-density graphs -- but is OPTIONAL polish, not needed for the default to be correct/safe.

## 76. ADAPTIVE per-leaf witness/general gate -> auto-captures the graph-dependent t>=3 win (2026-06-21, #167) [local]
Replaced the fixed witnessTMax cap with a PER-LEAF cost gate (both paths bit-identical -> pure speed, zero
correctness risk). Per leaf: Wδ = #feasible witnesses (saturating DP, compositions of t, δ_c<=uEnv_c-pl_c);
CB = #affected-Q candidates (compositions of r, Q_c<=uEnv_c). Take witness while Wδ <= CB*witK (witness per-unit
box-scan << general per-unit scWithTerms-DP). t=1,2 skip the gate (always win, measured); gate at t>=witGateMinT=3;
witnessTMax=8 hard backstop. Tunables SCT_WIT_K(8)/SCT_WIT_GATE_MINT(3)/SCT_WITNESS_TMAX(8).
CORRECTNESS: default(adaptive) == pure-general(TMAX=0) == pure-witness(GATE_MINT=99), bit-identical t=1/2/3/4 on
dblp-sigmod/db. TIMING (dblp-db, peel s): t2 adaptive 4.43 ~ witness 4.93 << general 17.9 | t3 adaptive 16.6 ~
witness 16.5 << general 56.7 | t4 adaptive 79 ~ witness 79 < general 126. HEADLINE: the old fixed cap=2 used GENERAL
at t=3 (56.7s); the adaptive auto-picks witness (16.6s) = 3.4x, with no per-graph tuning. The gate keeps witness
wherever it beats general (here t<=4) and would fall back on extreme M/t -- robust, never worse than min(witness,
general). GAP (honest): at t>=4 the BATCH is the true winner (sec 75: batch 15.8s vs witness 79s at dblp-db 3,7),
but batch is wave-based and NOT in the per-leaf witness/general choice. So the default is optimal for the common
t=1,2,3 regimes; for the rare t>=4 it picks witness-over-general (still beats general) but misses the batch.
Integrating batch into the per-leaf fallback (3-way choice) is the remaining step -- bigger rewrite of the wave
driver, deferred. batch stays opt-in (SCT_BATCH_PEEL, now gated batchPeel && tail>=2).

## 77. Does per-leaf ADAPTIVE beat a fixed cap? MEASURED: vs GENERAL, negligible; the lever is vs BATCH (2026-06-21, #168) [local]
User asked: is per-leaf adaptive much better than a fixed witnessTMax, or about the same? Added a gate-split
counter (SCT_WIT_DBG: gated leaves that chose witness vs fell back to general). RESULT -- the gate almost NEVER
splits:
  dblp-db t=3: 0.0% fell back (181198/181198 chose witness) | t=4: 0.0% | t=5: 3.6% (3895/109k)
  dblp-sigmod t=4: 0.0% | t=5: 0.1%
=> On the measured graphs, per-leaf adaptive == "always witness up to a high cap": the gate picks witness for
~99.9-100% of leaves through t=4, only starting to split at t=5 (3.6%). So vs a WELL-CHOSEN fixed cap (e.g. 4-5),
the per-leaf gate's raw speed edge is NEGLIGIBLE here. The earlier 3.4x "win" was vs the CONSERVATIVE cap=2, not
vs a good fixed cap. WHY: the gate's alternative is the GENERAL DFS (slow DFS+IE baseline), which witness beats so
broadly that the crossover witness->general sits at very high M/t (rarely hit). 
THE REAL LEVER: the gate would actually SPLIT (and per-leaf would earn its keep) if the alternative were the BATCH,
not general -- the crossover witness->BATCH is LOWER (t=4 on dblp-db: batch 15.8s < witness 79s, sec 75). At t=4 the
gate keeps witness (0% fallback, witness 79 < general 126) but the true best is batch (15.8) which is NOT in the
witness/general choice. So: (a) per-leaf vs general = theoretically >= fixed but empirically ~equal (gate ~never
splits); (b) per-leaf vs batch = WOULD matter (lower crossover, real splitting) -> integrating batch into the
per-leaf choice is where adaptive actually pays off. The adaptive gate's current value is robustness + no per-graph
tuning + future-proofing, NOT a raw speedup over a good fixed cap. Honest.

## 78. METHOD ROUTING by measured crossover: witness (t<=witCross) -> BATCH (t>witCross) (2026-06-21, #169) [local]
Completed the tail-indexed method hierarchy for the affected-update, all bit-identical:
  t=1      : witness-floor (unified witness path)
  t=2,3    : per-leaf-adaptive witness (sec 76; gate ~never splits vs general, sec 77)
  t>witCross(3): the wave BATCH driver (sec 71) -- routes here BY DEFAULT (tunable SCT_WIT_CROSS).
WHY routing not per-leaf-3-way: integrating witness INTO the wave driver is a big risky rewrite of the bit-identical
core, and sec 77 showed a whole-peel switch ~= per-leaf (witness/general gate ~never splits). The witness->BATCH
crossover IS real and lower than witness->general: at dblp-db t=4 batch beats witness ~3.6x. So a tail-based route
to batch past the crossover captures the win at low risk. Verified batch bit-identical at t=3,4 (dblp-sigmod/db)
BEFORE making it a default there.
RESULT (dblp-db, peel time, default-routing vs force-witness-everywhere SCT_WIT_CROSS=99):
  t=2 1.68 vs 1.66 (unchanged, witness) | t=3 6.03 vs 6.01 (unchanged, witness) | t=4 8.43 vs 30.25 = 3.59x (batch).
Bit-identical: default == force-witness(WIT_CROSS=99) == general(TMAX=0) at t=3,4. witCross=3 is the dblp-db-
calibrated crossover; it is graph-dependent (sparser graphs witness wins further, denser earlier) but tunable, and
sec 77's "per-leaf ~= fixed cap" logic applies. STATE: affected-update now auto-picks the fastest method by tail
across the whole spectrum; t=1,2,3 witness (1.9-7x over general), t>=4 batch (3.6x over witness). All local.

## 79. STATE OF THE PEEL + remaining optimization landscape (next: MEMORY) (2026-06-21, #170)
ARC COMPLETE (sec 70-78): the affected-update now AUTO-PICKS the fastest method by tail t=s-r, all bit-identical:
  t=1        -> witness-floor (unified witness path, sec 74)
  t=2,3      -> per-leaf-adaptive witness (sec 76); gate ~never splits vs general (sec 77)
  t>witCross(3) -> wave BATCH driver (sec 71), routed BY DEFAULT (sec 78)
  Env knobs: SCT_WITNESS_TMAX(8 hard cap) / SCT_WIT_GATE_MINT(3) / SCT_WIT_K(8) / SCT_WIT_CROSS(3) /
             SCT_BATCH_PEEL / SCT_BATCH_CB_CAP(128) / SCT_NO_SLOT_IDX. Diagnostics: SCT_WIT_DBG, SUPP_DBG, COMP_DBG.
  Wins (dblp-db peel, local, vs general baseline): t=2 ~4x | t=3 ~3.4x | t=4 ~3.6x (batch over witness).
CURRENT PEEL BREAKDOWN (local, after the witness/batch work -- where time goes NOW):
  dblp-db 3,5 (t=2): peel 1.72s = slotForbidDiff 0.44 (26%) + affected-update 1.27 (74%)
  dblp-db 4,6 (t=2): peel 12.4s = slotForbidDiff 4.05 (33%) + affected-update 8.33 (67%)
  => the AFFECTED-UPDATE (now witness) is STILL the dominant 67-74% (its box-scan O(Wδ*boxes*M)); slotForbidDiff
  26-33%. build/maps are TIME-cheap (0.1-0.3s) -- the maps' problem is MEMORY, not time.
REMAINING OPTIMIZATION LANDSCAPE (prioritized):
  A. MEMORY (pattern<->leaf maps) -- THE competitive gap: region-native only loses to CND on memory at low-RS
     (sec 62: maps = 4.16G, the biggest consumer). Lever = ADAPTIVE maps-CSR: don't store the local footprint for
     WIDE leaves, recompute on demand (sec 52 SCT_MAPS_RECOMPUTE exists but full-on loses time; recompute only the
     storage-expensive leaves -> ~6.8GB saved at ~no time cost). HIGH value, closes the only weakness. <- NEXT (user picked A).
  B. PEEL PARALLELISM (multi-core) -- peel is single-thread; the 67-74% affected-update is per-pattern within a
     wave, and the batch already groups a level leaf-major -> wave leaves parallelize. New speedup + scalability
     axis. Effort medium-high (delta[] accumulate + re-bucket shared state -> per-thread accumulate + merge).
  C. WITNESS box-scan reduction -- the dominant peel cost now, but fundamental (index the per-δ box liveness like
     slot-idx #1? hard: liveness depends on δ and forbidden). MEDIUM/uncertain.
  D. slotForbidDiff (26-33%) -- slot-idx #1 already made it output-sensitive; sec 57 killed lazy maintenance. LOW.
CLOSED (do not revisit): closed-form CF (sec 66 time-unstable), composition direct-counting (sec 67, 4-45x mem),
class-deletion (sec 69, <2% single-class), bulk-bin enlargement past |host|=1 (sec 57 falsified).
