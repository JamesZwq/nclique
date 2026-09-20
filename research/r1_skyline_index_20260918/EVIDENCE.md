# Evidence inventory (chain index line)

Every table in RESULTS_FINAL.md, in the paper (`Sigmod2027ChainIndex/make_tables.py`)
and on the artifact page is generated from the files below; nothing is typed by
hand. All runs are one thread. `*.json` hold one record per graph (`result` =
the tool's JSON line; `error` when a graph did not build) plus input sha256,
index-file sha256, host, peak RSS and wall time; `*-logs/` hold the raw
stdout/stderr of every process under `/usr/bin/time`.

## Current module (values as double, per-size widths, 64-bit solver ids)

| File | Machine | Content | Used for |
|---|---|---|---|
| `final.json`, `final-logs/` | laptop | 13 graphs: bytes, build phases, query latencies (three forms), selftest counts | RESULTS_FINAL.md Sections 9, 10; paper Tables size/queries/build (laptop rows); artifact page |
| `more.json`, `more-logs/` | laptop | 4 more graphs (ca-HepTh, email-Eu-core, com-amazon, dblp-coauthor) | same |
| `buildonly.json`, `buildonly-logs/` | laptop | 13 graphs, build only: phases, resident memory per phase, peak RSS | Section 14 (build memory); paper Table prior (build peak) |
| `tods1.json`, `tods1-logs/` | tods1 | 20 graphs; com-lj, hollywood, orkut failed (member ID overflow, before the 64-bit fix) | Section 15; paper Tables (server rows) |
| `tods2.json`, `tods2-logs/` | tods2 | 7 graphs incl. cit-Patents, web-BerkStan, web-NotreDame | same |
| `tods1_big.json`, `tods1_big2.json` + logs | tods1 | the three largest graphs after the 64-bit fixes: com-lj bad_alloc at 350 GB, hollywood OOM at 435 GB, orkut stopped | Section 15 (limits) |
| `prior/prior_<graph>.json` | laptop | the production single-size pipeline (src/, ST_V3 + BuildHier) run once per size on five graphs: per-size build/peel/hierarchy times, wall, peak RSS, branch counts | Section 16; paper Exp "against the existing single-size tool", Table prior |

Merged view: `report_tables.py` (`merged()`) dedups the three machines by
(n, m): 29 distinct graphs.

## Design study (stage 2, laptop, five graphs, integer values)

| File | Content |
|---|---|
| `stages/index_vertices.json` | one S tree per size over vertices: the baseline's measured bytes and memcpy listing latency |
| `stages/index.json` | S trees over twin classes, and the SGL-style skyline dedup |
| `stages/index_chains.json` | S trees over chains (explicit vertex-to-chain map) |
| `stages/index_aligned.json` | chains with aligned labels (first run-array-free layout) |
| `stages/counts.json` | stage 1 counts (canonical nodes, residue cells, F2 checks) |

Paper Table layouts and RESULTS_FINAL.md Section 9.8 come from these.

## Superseded evidence (kept for the record)

| File | What changed afterwards |
|---|---|
| `archive/final_v1.json` | first five-graph run of the module (fixed-width values, scalar fill) |
| `archive/final_v2.json` | 13 graphs, scalar-tail fill |
| `archive/final_v3.json` | 13 graphs, branchless fill, fixed-width values |
| `archive/final_v4_int.json`, `archive/more_v4_int.json` | 17 graphs, integer values with per-size widths (before double) |
| `archive/index_aligned_v1.json` | aligned layout with the slow push_back expansion |

## Additional experiments (2026-09-21, RESULTS_FINAL.md Section 17)

- `compress.json` (laptop, bytes only): xz/zstd on the S trees materialized by `strees_dump` (input labels and aligned labels) and on the index file, 17 graphs; `archive/compress_aligned_only.json` is the first run (aligned labels only).
- `stages/index_vertices_tods2.json` + `stages/index-logs_vertices_tods2/`: S trees latency baseline on the seven tods2 graphs (`run_strees.py`, launcher `tods2_strees.sh`); tods1 (20 graphs) queued behind the CND sweep (`tods1_after_prior.sh`).
- `scale_tods2.json` + `scale_tods2-logs/`: vertex-induced samples at 20/40/60/80 percent of cit-Patents and web-BerkStan (`sample_vertices.py`, seed 20260921), full bench (`tods2_scale.sh`); the 100 percent points are the tods2.json rows.
- `profile_tods2.json` + `profile_tods2-logs/`: query latency profiles (`query_profile`: the fixed workload by level, own-level deciles by answer size, stratified own-level queries per clique size, seed 20260921 for the strata) on the seven tods2 graphs and the eight samples (`tods2_profile.sh`).
- `prior/tods1/`: CND sweep on the twenty tods1 graphs (in progress).

## Not measured (known gaps)

- The S-tree baseline's query latency existed only for the five laptop graphs (`stages/index_vertices.json`) until 2026-09-21; the server runs are listed above.
- Laptop records in `final.json`/`more.json` have no `peak_rss_bytes` field; the laptop peaks are in `buildonly.json` and in the `final-logs/*.log` time output.
- The CND sweep covers the five laptop graphs, the seven tods2 graphs and (in progress) the twenty tods1 graphs.
- Query latencies mix machines across tables; each row names its machine. Server single-core speed is 1.5 to 3 times lower than the laptop's.
