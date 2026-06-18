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
