# Experiment data archive — VLDB R=1 paper

Snapshot date: 2026-04-25

This directory consolidates every experimental CSV cited in the paper
(`vldbNuclearR1/main.tex`).  It is committed to the repository so the
data behind every figure and every table is reproducible from origin.

Numbering reflects approximate order of appearance in the paper.

---

## Top-level performance

### `01_main_benchmark_762.csv`
- **Source:** `tods2:nclique/paper_data/01_main_benchmark_all_graphs.csv` (produced by `bench_r1_main.py`)
- **Used in:** §7.2 end-to-end speedup, §7.3 memory reduction, §7.5 per-graph breakdown, §8.2 CS-1
- **Schema:** `graph,r,s,algorithm,time_ms,memory_kB,status`
- **Algorithms:** `Ours_ST` (PIVOTER_RUN_ST=1) and `REF_R1` (no env flag); r=1 only
- **Rows:** ~1530 data rows = 762 (graph, s) configs × 2 algos
- **Graphs (10):** com-amazon, com-dblp, twitter, web-Stanford, com-youtube, web-Google, wiki-Talk, web-it-2004, soc-pokec, com-orkut
- **s ranges per graph:** see Table 2 of the paper (`tab:datasets`)

### `02_breakdown_summary.csv`
- **Source:** `tods2:nclique/results/breakdown/breakdown_summary.csv` (ours-side via `PIVOTER_RUN_ST`)
- **Used in:** §7.4 time-and-memory breakdown by phase
- **Schema:** `graph,s,algo,total_time_ms,peak_rss_kb,load_ms,build_ms,peel_ms`
- **Rows:** 34 (1 header + 33 data) — 3 graphs × 17 (graph,s) × 2 algos
- **Configs covered:**
    - com-youtube: s ∈ {3, 5, 8, 12, 16}
    - web-Stanford: s ∈ {3, 5, 10, 20, 40, 60}
    - web-it-2004:  s ∈ {3, 10, 30, 100, 200, 400}

### `03_breakdown_median.csv`
- **Source:** same experiment as `02_breakdown_summary.csv`, finer granularity.
- **Used in:** detailed phase-level analysis (e.g. per-component memory)
- **Schema:** `graph,s,algo,phase,n_runs,time_ms,rss_kb,delta_rss_kb,component_bytes`
- **Rows:** 222 — one row per (graph, s, algo, phase). Median over 3 runs.
- **Phase names (Ours, PIVOTER_RUN_ST):** `loadAndSort`, `buildSDCT`, `preMutation`, `prepareGraph`, `dispatch_total`
- **Phase names (Ref, SOTA):** `loadAndSort`, `buildSDCT`, `preMutation`, `prepareGraph`, `REF_initSupports`, `REF_heapBuild`, `REF_peel_loop`, `dispatch_total`

### `04_breakdown_raw.tsv`
- **Source:** raw event log of every run × every phase
- **Used in:** rebuilds 02 and 03; provenance for variance / outlier checks
- **Schema (TAB-separated):** `meta\tphase\tduration_ms\trss_kb\tdelta_rss_kb\tcomponent_bytes`
- **`meta` column:** `graph=<g>,r=1,s=<s>,algo=<ours|ref>,run=<0|1|2>`
- **Rows:** 663 (1 header + 102 runs × ~6.5 phases avg)

---

## Case studies

### `05_cs1_speed_keytable.csv`
- **Source:** `case_study/cs1_speed/cs1_key_table.csv`
- **Used in:** §8 motivation; subset of 01 highlighting the headline configs
- **Schema:** `graph,s,REF_ms,ST_ms,Speedup,REF_MB,ST_MB,MemRedux`

### `06_cs2_dblp_granularity.csv`
- **Source:** `case_study/cs2_dblp/cs2_summary.csv`
- **Used in:** §8.2 case study 1 — granularity is real (active-set shrinkage, top-K stability)
- **Schema:** `s,nonzero_count,pct_of_graph,max_core,top100_min,top1K_min`
- **Graph:** com-dblp, r=1 only

### `07_cs3_groundtruth_metrics.csv`
- **Source:** `case_study/cs3_groundtruth/cs3_metrics_K10000.csv`
- **Used in:** §8.3 case study 2 — ground-truth quality vs k-core (Tab 4)
- **Schema:** `method,s,top_size,prec@K,best_F1,best_recall,rec@50%,rec@70%`
- **Graph:** com-dblp; ground truth: SNAP top-5000 communities

### `08_cs3_density.csv`
- **Source:** `case_study/cs3_groundtruth/cs3_density.csv`
- **Used in:** §8.3 induced triangle density column
- **Schema:** `method,K,top_size,edges,density,triangles,tri_per_v`

### `09_cs4_crossr_dblp_s10.csv`
- **Source:** `case_study/cs4_crossr/cs4_crossr.csv`
- **Used in:** §8.4 case study 3 — Pareto cross-r at s=10 (Tab 5)
- **Schema:** `method,K,top_size,prec,rec50,density,tri_per_v,time_ms`

### `10_cs5_bestf1.csv`
- **Source:** `case_study/cs5_bestf1/cs5_bestf1_summary.csv`
- **Used in:** §8 supporting CS-5 (best-F1 per community; weak signal as expected)
- **Schema:** `method,median,mean,top10,top100,n_F1>=0.3,n_F1>=0.5,n_F1>=0.7`

### `11_cs6_grid_dblp.csv`
- **Source:** `case_study/cs6_grid/cs6_grid.csv`
- **Used in:** §8.4 Fig fig:cs6 (Pareto plot over r×s grid on com-dblp)
- **Schema:** `r,s,top,prec,rec50,density,tri_per_v,bf1_top10,runtime_ms`
- **Grid:** r ∈ {1,2,3}, s ∈ {5,10,15,20}; 12 cells

### `12_cs7_youtube_crossr.csv`
- **Source:** `case_study/cs7_youtube/cs7_youtube.csv`
- **Used in:** §8.5 case study 4 — cross-domain validation (Tab 6)
- **Schema:** `method,top,prec,rec50,rec70,density,tri_per_v,runtime_ms`
- **Graph:** com-youtube; methods include k-core + (1,s)/(2,s)/(3,s) for s ∈ {5,10,15}

### `13_cs8_survival_dblp.csv`
- **Source:** `case_study/cs8_survival/cs8_survival.csv`
- **Used in:** §8 supporting CS-8 — survival profile per s (Nuclear-CD Case Study I analog)
- **Schema:** `method,k,verts,ccs`
- **Methods:** s ∈ {2, 3, 5, 10, 15, 20, 25}; sweeps k-rank percentile

---

## Reproduce instructions

| File | Command to reproduce |
|---|---|
| 01 | `python3 benchmark_all.sh` (pre-existing harness; output `benchmark_all_results.csv`) |
| 02–04 | `python3 run_breakdown.py --bin ./build/bin/degeneracy_cliques --graph-dir ./graphs --graphs com-youtube web-Stanford web-it-2004 --s-list <per-graph> --runs 3 --resume --out results/breakdown` |
| 05–13 | `python3 case_study/cs<N>_<name>/analyze_cs<N>.py` (per-CS scripts) |

For 02–04 specifically, the driver script is `scripts_run_breakdown_all.sh`
on tods2 (`/home/wenqianz/nclique/`).

---

## Notes / caveats

- **`Ours_ST` in 01 is `PIVOTER_RUN_ST`**, the immutable-tree variant. The
  tree-free variant (`PIVOTER_RUN_ST_V2`, §5 of the paper) has a known
  peel-phase segfault on `web-Stanford` at $s \ge 8$; running the breakdown
  with ST_V2 instead of ST is a follow-up after that bug is fixed.
- **`02–04` use the same `Ours_ST` algorithm** as `01`, so the per-phase
  breakdown is consistent with the headline 6.16× geometric-mean speedup.
- **All times in milliseconds**, all memory in kilobytes (kB) for raw CSVs;
  the paper reports MB for memory.  Conversion: MB = kB / 1024.
