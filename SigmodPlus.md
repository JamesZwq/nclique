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

| Resource | Default cap     | Where enforced                                  |
|----------|-----------------|-------------------------------------------------|
| Time     | 1800 s (30 min) | `subprocess.run(..., timeout=...)` per cell     |
| Memory   | 200–250 GB      | `--mem-cap-gb` flag; checked via `time -v` line |
| Disk     | unbounded       | logs land under `/data/wenqianz/...`            |

When a cell is killed by timeout the script also calls
`os.killpg(pid, SIGKILL)` on the child group so the binary cannot
survive `/usr/bin/time -v` and become an orphan.

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
  within 30-min budget.  These graphs are an honest infeasibility
  regime for RegND and are not in Fig 3.
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
