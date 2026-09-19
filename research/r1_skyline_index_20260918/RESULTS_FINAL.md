# The Chain Index: Final Module, Measured

Date: 2026-09-19. Module `chain_index.hpp` (header-only, `namespace
chainindex`), tool `chain_index_tool.cpp` (`--selftest`, `--bench`),
driver `run_final.py`. Evidence: `final.json`, `final-logs/` (one
`/usr/bin/time -l` log per graph, build and selftest logs), index files
`cx/<graph>.cx` (not committed). Tables are printed by `report_tables.py`
from `final.json` and the stage-2 files `stages/index*.json`. Theory:
[THEORY.md](THEORY.md) (hierarchy, F1-F8, canonical nodes, certified
tail), [CHAINS.md](CHAINS.md) (chains C1-C3, layout C4-C6).

## 1. Problem Summary

Store, for one graph and r = 1, every (1, s)-nucleus hierarchy for every
clique size s at once: the core value kappa_s(v) of every vertex at every
size, and every (s, k)-nucleus (the s-clique-connected components of the
vertices with kappa_s >= k), so that a query returns
- Q1 the value kappa_s(v),
- Q2 the community of v at (s, k) as a vertex set.
(Membership of u in v's community and the ladder of nested communities
above v follow from Q2 and are not separate targets; the module answers
them and the selftest checks them, the report does not measure them.)
Exactness is mandatory; the targets are bytes and query latency; the
build must be one pass of the existing all-size peel.

## 2. Current Baseline And Bottlenecks

The strongest simple exact representation is one S tree per size: the
canonical merge tree of size s with a DFS array of vertex ids, a
community being one contiguous slice (stage 2, mode `vertices`,
[RESULTS_INDEX.md](stages/RESULTS_INDEX.md), `stages/index_vertices.json`). It lists a
community by one memcpy (0.09-0.15 ns per output vertex) and costs
30-60 bytes per vertex (272 KB GrQc to 24.9 MB Stanford, values
included). Its redundancy is the vertex axis: a vertex appears in
omega(v) - 1 arrays and carries omega(v) - 1 values, although vertices
with the same hierarchy position are indistinguishable. The SGL index of
Zong et al. (SIGMOD 2026) removes that redundancy for bipartite cores by
atomic units that are neither nested nor overlapping; at r = 1 the
analogue is stronger, a partition (chains, CHAINS.md C1-C3).

## 3. Candidate Algorithm Ideas

1. Skyline dedup (SGL transfer, stage 2): one entry per vertex per
   non-dominated size. Measured: 0.99-2.28x fewer bytes than S trees over
   the same classes and 2-3x slower queries; a space/time trade only.
2. S trees over twin classes: 1.03-1.28 vertices per class at r = 1,
   the class map costs more than it saves.
3. S trees over chains (stage 2, mode `chains`): 2.5-4.5x smaller than
   over twins and 1.7-8.6x faster listing; the vertex-to-chain map is
   45-58 percent of the bytes on the large graphs.
4. Aligned labels (stage 2, mode `aligned`): relabel vertices so that a
   chain is one id range; the map becomes a bitmap with rank; community
   answers become lists of id ranges. 3.7-8.3x smaller than per-vertex S
   trees; the explicit-id expansion was 2x slower than memcpy.
5. Run arrays with node entry points (this module): order chains
   lexicographically by their own-node preorder tuple, walk each tree by
   smallest rank first, store the DFS array as maximal runs of
   consecutive ids, store per node where its segment enters the run
   array. A community is located in O(1) and reported as its minimal
   list of ranges; at s = 2 one range (CHAINS.md C4-C6).

## 4. Chosen Approach And Rationale

Candidate 5 on top of 3 and 4, with the chain-id DFS array (candidate 3
with aligned labels) kept as the build form and measured as the
ablation. Rationale: chains remove the vertex-axis redundancy exactly
(proved), aligned labels remove the map, runs turn the remaining
per-(chain, size) array into per-run pairs and make the answer a pointer.
The skyline dedup is dropped: dominated once chains exist.

Index blocks (all flat arrays, one file, magic `CHAINX02`):
- Block 1 (map): bitmap of chain starts over the n labels with a rank
  directory (one word per 64 labels), and `start_pos` per chain.
- Block 2 (per chain): omega, sigma (first certified size), trajectory
  (own node id per size, omega - 1 entries), residue values for
  2 <= s < sigma; values at s >= sigma are the binomial C(omega - 1, s - 1)
  (THEORY.md F8).
- Block 3 (per size s): node tops (W bytes), parent, subtree size,
  entry (run index, label) per node plus a sentinel, and the run array
  (lo, hi label pairs). Jump pointers (Myers' skip pointer, one per node)
  are rebuilt at load time and are not stored.

Queries: Q1 = rank on the bitmap, then residue or binomial. Q2 = own
node from the trajectory, climb to the highest ancestor with top >= k
(jump pointers, tops decrease upward), then Lemma C6: head range, whole
runs, tail range.

## 5. Complexity Discussion

Let C be the number of chains, N_T the canonical nodes over all sizes, P
the number of (chain, size) pairs (sum over chains of omega - 1), R the
number of runs over all sizes, and W the count width.
- Bytes: n/8 + 4 n/64 + 4 C (map) + 2 C + 4 P + W x residue cells
  (chains) + (W + 16) N_T + 8 R (layers). Per-vertex S trees cost
  4 x (sum over v of omega(v) - 1) + W x active pairs + (W + 12) N_T.
- Build: one all-size peel (the terminal solver of
  `research/r1_terminal_20260918/`), one union-find pass per size for the
  trees, one map insertion per vertex for the chains (key length
  omega(v) - 1), one sort of items per node for the DFS arrays.
- Queries: Q1 O(1); Q2 O(depth of the climb) with jump pointers plus
  O(1) to locate, O(ranges) to copy, O(vertices) to expand.
The count of chains is bounded by n and by the twin classes; it is below
N_T on every input measured but not by a theorem (CHAINS.md Section 7).

## 6. Implementation Summary

- `chain_index.hpp` (about 180 lines): `ChainIndex<T>` with the blocks
  above, `chain_of`, `own_node`, `climb`, `value`, `community_runs`
  (pointer form), `community_ranges` (vector form, both forms),
  `expand` (branchless fill: every block of eight ids is stored
  unconditionally with 128-bit vector stores and the pointer advances by
  the true length, so the caller's buffer carries eight spare slots),
  `compact_runs`, `save`, `load`, byte counters; `member` and `ladder`
  as derived queries used by the selftest.
- `chain_index_tool.cpp`: `build_chain_index` (solve, trees with initial
  preorder ids, chains keyed by preorder tuples in lexicographic order
  with the inactive chain last, aligned labels, per-chain blocks, per-size
  smallest-rank DFS producing final node ids and the chain-id array),
  `selftest_graph` (brute-force nuclei from clique masks; checks the build
  form, the compact form, and the loaded file), `bench`.
- `run_final.py`: Release and ASan/UBSan builds, both selftests, five
  graphs one at a time under `/usr/bin/time -l`, `final.json`.
- Compiler: Apple clang, `-O3 -std=gnu++23`, one thread
  (`OMP_NUM_THREADS=1`), macOS arm64 (Apple M-series).

## 7. Correctness Validation Summary

Brute force: for every labelled graph on at most 6 vertices, 200 random
graphs on 7-10 vertices, split graphs and K8 (34,075 graphs), nuclei are
computed from clique masks and union-find, and every (v, s, k) with
1 <= k <= kappa_s(v) is queried. Each graph is checked three times: build
form (chain-id arrays), compact form (runs), and the compact form loaded
back from disk. Checks per run: values 1,921,263 (including s = omega + 1
returning 0), community queries 1,816,701 (set equality after relabeling;
ranges well formed and fully merged; pointer form expands to the same
set), membership 10,955,679, ladders 1,039,647 (levels strictly
decreasing, counts equal to the community sizes), and every
(2, k)-community is exactly one range. Release and ASan/UBSan builds both
pass (`final-logs/build-2.log`, `build-asan-2.log`).

Consistency across pipelines on the real graphs: canonical node counts
(1,517 / 3,538 / 33,979 / 53,239 / 65,134) equal those of stage 1 and
stage 2, which used the frozen reference kernel instead of the terminal
solver; chain counts equal those of `chains.cpp` (13,459 dblp, 17,963
Stanford, 40,867 amazon) up to the one chain of vertices without edges
that the module keeps (GrQc 716, HepPh 1,136).

## 8. Experimental Setup

Machine: the local Apple M-series laptop of the earlier stages, one
thread (`OMP_NUM_THREADS=1`), no other timing run concurrent. Seventeen
inputs in five families: collaboration (`data/ca-GrQc.edges`,
`data/ca-HepPh.edges`, `graphs/ca-HepTh.edges`, `graphs/ca-AstroPh.edges`,
`graphs/ca-CondMat.edges`, `data/com-dblp.edges`, and the dense
4.05 M-vertex `graphs/dblp-coauthor.edges` with s_max 450), citation
(`graphs/cit-HepPh.edges`), web (`graphs/web-Stanford.edges`), product
(`graphs/amazon0302.edges`, `graphs/amazon-copurchase.edges`), social and
communication (`graphs/email-Eu-core.edges`, `graphs/loc-Brightkite.edges`,
`graphs/soc-Epinions1.edges`, `graphs/soc-Slashdot0902.edges`,
`graphs/com-youtube.edges`, `graphs/soc-pokec.edges`). The first thirteen
are in `final.json`, the last four (HepTh, email, com-amazon,
dblp-coauthor) in `more.json` (sha256 of every input recorded); both
files are from the double-valued module (Section 14), the integer-valued
run of the same day is `archive/final_v4_int.json` and
`archive/more_v4_int.json`. Solver count width chosen from its rows: 64
bits except com-dblp (128), ca-HepPh (256, s_max 239) and dblp-coauthor
(512, s_max 450); the index stores doubles regardless. Queries (seed 20260918, drawn from active vertices
and sizes 2 <= s <= omega(v)): community regimes own (k = kappa_s(v),
20,000), half (k = max(1, kappa/2), 20,000), root (k = 1, 1,000);
values 200,000. One warm-up pass plus five timed passes, median reported.
(Runs before 2026-09-19 evening also timed membership and ladder queries;
those fields remain in the archived JSON files and are not reported.)
Explicit ids are written into a preallocated caller buffer. The
per-vertex S-tree bytes are computed by the tool with the stage-2
`vertices` accounting (verified equal to `stages/index_vertices.json` on the
five stage-2 graphs, 16,710,168 bytes on dblp); the memcpy listing
baseline exists only for those five graphs (stage-2 run, same protocol,
its own query draw). Evidence lineage: `archive/final_v1.json` / `archive/final-logs_v1/` is the
first five-graph run (before the baseline field and the extra graphs);
`archive/final_v2.json` / `archive/final-logs_v2/` is the 13-graph run with the earlier
fill (eight ids per step, scalar tail); `final.json` / `final-logs/` is
the 13-graph run with the branchless fill that the module now uses.
Bytes, build times and every non-listing latency agree across the runs
within run-to-run noise; the explicit-listing columns differ by the fill.

Commands:
```
cmake -S research/r1_skyline_index_20260918 -B research/r1_skyline_index_20260918/build -DCMAKE_BUILD_TYPE=Release -DSANITIZE=OFF
cmake --build research/r1_skyline_index_20260918/build -j 12 --target chain_index_tool
research/r1_skyline_index_20260918/build/chain_index_tool --selftest
/usr/bin/time -l research/r1_skyline_index_20260918/build/chain_index_tool --bench data/com-dblp.edges research/r1_skyline_index_20260918/cx/com-dblp.cx
python3 research/r1_skyline_index_20260918/run_final.py graphs/ca-AstroPh.edges graphs/ca-CondMat.edges graphs/cit-HepPh.edges graphs/loc-Brightkite.edges graphs/soc-Epinions1.edges graphs/soc-Slashdot0902.edges graphs/com-youtube.edges graphs/soc-pokec.edges
python3 research/r1_skyline_index_20260918/run_final.py --tag more --only graphs/ca-HepTh.edges graphs/email-Eu-core.edges graphs/amazon-copurchase.edges graphs/dblp-coauthor.edges
python3 research/r1_skyline_index_20260918/report_tables.py
```

## 9. Experimental Results

### 9.1 Size
| Graph | n | s_max | W bits | chains | n / chains | canonical nodes | chains / nodes | (chain,s) pairs | runs | map B | chains B | layers B | total B | file B | perm B | build form total B | per-vertex S trees B | ratio | ratio with perm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 5,242 | 44 | 64 | 716 | 7.3 | 1,517 | 0.47 | 2,079 | 624 | 3,856 | 17,321 | 40,151 | 61,328 | 57,136 | 20,968 | 66,904 | 272,564 | 4.44x | 3.31x |
| ca-HepPh | 12,008 | 239 | 256 | 1,136 | 10.6 | 3,538 | 0.32 | 8,468 | 4,969 | 6,808 | 52,117 | 134,157 | 193,082 | 188,801 | 48,032 | 186,952 | 1,899,628 | 9.84x | 7.88x |
| com-dblp | 317,080 | 114 | 128 | 13,459 | 23.6 | 33,979 | 0.40 | 70,208 | 28,784 | 113,304 | 460,177 | 1,068,332 | 1,641,813 | 1,510,643 | 1,268,320 | 1,729,876 | 16,710,168 | 10.18x | 5.74x |
| web-Stanford | 281,903 | 72 | 64 | 17,963 | 15.7 | 53,239 | 0.34 | 151,364 | 97,074 | 124,720 | 1,315,734 | 2,112,252 | 3,552,706 | 3,342,774 | 1,127,612 | 3,737,596 | 24,853,426 | 7.00x | 5.31x |
| amazon0302 | 262,111 | 7 | 64 | 40,867 | 6.4 | 65,134 | 0.63 | 143,567 | 47,380 | 212,628 | 1,122,338 | 1,746,902 | 3,081,868 | 2,821,691 | 1,048,444 | 3,875,704 | 14,414,214 | 4.68x | 3.49x |
| ca-AstroPh | 18,772 | 57 | 64 | 3,331 | 5.6 | 6,462 | 0.52 | 36,381 | 25,742 | 16,860 | 232,019 | 368,353 | 617,232 | 593,793 | 75,088 | 616,488 | 2,216,476 | 3.59x | 3.20x |
| ca-CondMat | 23,133 | 26 | 64 | 1,917 | 12.1 | 3,357 | 0.57 | 8,693 | 3,658 | 12,020 | 61,272 | 104,370 | 177,662 | 165,372 | 92,532 | 204,208 | 1,299,206 | 7.31x | 4.81x |
| cit-HepPh | 34,546 | 31 | 64 | 16,623 | 2.1 | 9,369 | 1.77 | 104,308 | 61,217 | 72,980 | 831,455 | 706,524 | 1,610,959 | 1,574,826 | 138,184 | 2,093,928 | 3,199,324 | 1.99x | 1.83x |
| loc-Brightkite | 58,228 | 53 | 64 | 6,326 | 9.2 | 11,547 | 0.55 | 35,833 | 21,427 | 36,232 | 319,597 | 455,239 | 811,068 | 767,125 | 232,912 | 915,304 | 2,628,872 | 3.24x | 2.52x |
| soc-Epinions1 | 75,879 | 68 | 64 | 8,821 | 8.6 | 12,886 | 0.68 | 55,344 | 38,876 | 49,524 | 456,481 | 611,619 | 1,117,624 | 1,068,940 | 303,516 | 1,304,016 | 3,221,518 | 2.88x | 2.27x |
| soc-Slashdot0902 | 82,168 | 56 | 64 | 6,798 | 12.1 | 6,309 | 1.08 | 30,496 | 17,551 | 42,608 | 269,319 | 289,200 | 601,127 | 578,259 | 328,672 | 734,136 | 2,995,392 | 4.98x | 3.22x |
| com-youtube | 1,134,890 | 52 | 64 | 42,815 | 26.5 | 22,403 | 1.91 | 179,601 | 93,541 | 384,064 | 1,477,821 | 1,241,440 | 3,103,325 | 3,015,917 | 4,539,560 | 3,987,484 | 33,847,412 | 10.91x | 4.43x |
| soc-pokec | 1,632,803 | 48 | 64 | 383,206 | 4.3 | 67,544 | 5.67 | 2,375,352 | 1,212,885 | 1,838,988 | 19,247,860 | 11,304,285 | 32,391,133 | 32,122,997 | 6,531,212 | 44,018,876 | 101,770,646 | 3.14x | 2.61x |
| ca-HepTh | 9,877 | 32 | 64 | 1,104 | 8.9 | 1,583 | 0.70 | 3,112 | 1,044 | 6,284 | 26,475 | 43,388 | 76,147 | 71,199 | 39,508 | 87,244 | 422,086 | 5.54x | 3.65x |
| email-Eu-core | 1,005 | 35 | 64 | 667 | 1.5 | 2,602 | 0.26 | 6,080 | 4,423 | 2,868 | 47,387 | 95,335 | 145,590 | 136,689 | 4,020 | 166,432 | 186,738 | 1.28x | 1.25x |
| amazon-copurchase | 548,552 | 7 | 64 | 56,315 | 9.7 | 86,798 | 0.65 | 163,658 | 54,803 | 328,132 | 1,385,981 | 2,261,230 | 3,975,343 | 3,628,510 | 2,194,208 | 4,840,472 | 19,697,700 | 4.95x | 3.19x |
| dblp-coauthor | 4,049,537 | 450 | 512 | 363,201 | 11.1 | 274,215 | 1.32 | 2,958,470 | 1,761,549 | 2,212,112 | 21,089,986 | 20,982,755 | 44,284,853 | 43,206,515 | 16,198,148 | 47,832,116 | 530,275,990 | 11.97x | 8.77x |

"total B" is in memory and includes the derived jump pointers (4 bytes
per node); "file B" is the disk image. "perm B" is the 4 n-byte
permutation from input labels to aligned labels, needed only if the graph
is not stored in the aligned order. "build form total B" is the same index
before `compact_runs` (chain-id DFS arrays).

### 9.2 Build, save, load
| Graph | solve (all-size peel) | trees | chains + labels | layout | build total | compact | save | load | process wall s | peak RSS MB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 2 | 5 | 0 | 0 | 7 | 0.06 | 1.0 | 0.2 | 0.50 | 10 |
| ca-HepPh | 259 | 213 | 0 | 2 | 477 | 0.15 | 2.6 | 0.2 | 1.90 | 26 |
| com-dblp | 243 | 387 | 2 | 8 | 641 | 0.52 | 5.1 | 0.5 | 19.58 | 155 |
| web-Stanford | 1,089 | 1,542 | 3 | 21 | 2,657 | 1.65 | 6.4 | 1.2 | 31.68 | 334 |
| amazon0302 | 101 | 209 | 1 | 9 | 321 | 1.07 | 3.9 | 0.8 | 16.65 | 123 |
| ca-AstroPh | 56 | 87 | 0 | 2 | 146 | 0.26 | 1.2 | 0.2 | 6.17 | 24 |
| ca-CondMat | 7 | 15 | 0 | 0 | 23 | 0.11 | 2.5 | 0.1 | 1.38 | 15 |
| cit-HepPh | 176 | 228 | 1 | 5 | 411 | 0.74 | 3.2 | 0.4 | 24.73 | 95 |
| loc-Brightkite | 295 | 482 | 0 | 3 | 781 | 0.33 | 4.9 | 0.3 | 4.89 | 87 |
| soc-Epinions1 | 679 | 1,391 | 0 | 4 | 2,075 | 0.44 | 0.9 | 0.3 | 14.00 | 492 |
| soc-Slashdot0902 | 445 | 874 | 0 | 2 | 1,322 | 0.32 | 4.5 | 0.2 | 7.52 | 261 |
| com-youtube | 518 | 977 | 4 | 10 | 1,510 | 1.05 | 6.2 | 0.7 | 44.71 | 555 |
| soc-pokec | 7,073 | 9,009 | 13 | 118 | 16,244 | 19.65 | 33.5 | 8.5 | 632.57 | 3,226 |
| ca-HepTh | 3 | 6 | 0 | 0 | 9 | 0.04 | 2.3 | 0.1 | 0.49 | 11 |
| email-Eu-core | 26 | 38 | 0 | 1 | 65 | 0.10 | 0.7 | 0.1 | 2.80 | 17 |
| amazon-copurchase | 88 | 271 | 3 | 13 | 376 | 1.71 | 5.9 | 1.4 | 18.63 | 158 |
| dblp-coauthor | 122,171 | 64,102 | 58 | 430 | 186,927 | 20.27 | 45.0 | 12.3 | 1067.96 | 8,785 |

Process wall time and peak RSS are for the whole `--bench` run (graph
load, the solver's row index, build, three forms, all query passes);
"build total" excludes the row index, whose time is the `ti_ms` field of
`final.json` (Section 14 itemises both). The peak is dominated by the
solver's row index plus the graph (Section 14).

### 9.3 Community queries (ns per query; compact form loaded from disk)
| Graph | regime | output vertices | ranges | locate (pointer) | ranges copied | explicit ids | per-vertex S trees memcpy | build form ranges | build form explicit |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | own | 1,658 | 14.1 | 24.0 | 46 | 176 | 225 | 178 | 386 |
| ca-GrQc | half | 2,108 | 15.8 | 27.3 | 49 | 200 | 305 | 286 | 462 |
| ca-GrQc | root | 2,208 | 18.6 | 8.5 | 37 | 215 | 206 | 276 | 449 |
| ca-HepPh | own | 5,327 | 129.8 | 14.4 | 49 | 681 | 625 | 819 | 679 |
| ca-HepPh | half | 6,257 | 137.0 | 25.4 | 102 | 668 | 832 | 882 | 1,677 |
| ca-HepPh | root | 7,439 | 169.4 | 7.9 | 100 | 774 | 794 | 948 | 1,958 |
| com-dblp | own | 164,832 | 871.7 | 27.6 | 295 | 12,455 | 23,493 | 9,706 | 20,962 |
| com-dblp | half | 199,796 | 951.0 | 37.9 | 311 | 17,634 | 26,144 | 16,724 | 31,005 |
| com-dblp | root | 229,952 | 1236.8 | 12.6 | 441 | 26,862 | 30,259 | 9,940 | 23,777 |
| web-Stanford | own | 83,227 | 1217.7 | 44.6 | 491 | 6,315 | 11,983 | 7,912 | 20,309 |
| web-Stanford | half | 102,997 | 1339.6 | 76.3 | 522 | 9,595 | 12,157 | 8,315 | 18,571 |
| web-Stanford | root | 115,934 | 1772.8 | 17.5 | 619 | 17,269 | 19,274 | 8,925 | 18,526 |
| amazon0302 | own | 105,997 | 734.6 | 15.2 | 194 | 5,321 | 14,277 | 14,215 | 18,013 |
| amazon0302 | half | 151,391 | 1283.2 | 13.9 | 305 | 7,627 | 24,757 | 23,171 | 27,452 |
| amazon0302 | root | 160,229 | 1529.7 | 5.1 | 338 | 8,792 | 23,022 | 18,515 | 25,197 |
| ca-AstroPh | own | 8,854 | 960.8 | 5.6 | 204 | 1,145 | - | 1,807 | 2,926 |
| ca-AstroPh | half | 9,701 | 984.4 | 21.2 | 247 | 1,240 | - | 2,003 | 4,147 |
| ca-AstroPh | root | 11,775 | 1140.8 | 6.2 | 267 | 1,275 | - | 2,000 | 5,499 |
| ca-CondMat | own | 9,999 | 165.0 | 8.2 | 74 | 860 | - | 660 | 1,676 |
| ca-CondMat | half | 12,008 | 181.1 | 20.7 | 85 | 949 | - | 829 | 1,861 |
| ca-CondMat | root | 14,389 | 214.9 | 7.6 | 81 | 1,220 | - | 719 | 2,271 |
| cit-HepPh | own | 18,252 | 3693.7 | 5.1 | 1,013 | 3,721 | - | 9,594 | 17,585 |
| cit-HepPh | half | 21,348 | 4107.8 | 28.8 | 1,150 | 3,925 | - | 11,032 | 14,123 |
| cit-HepPh | root | 26,461 | 4941.0 | 6.0 | 1,368 | 5,147 | - | 12,015 | 15,934 |
| loc-Brightkite | own | 32,323 | 400.7 | 11.5 | 108 | 1,812 | - | 3,788 | 5,030 |
| loc-Brightkite | half | 39,732 | 415.3 | 11.3 | 115 | 2,210 | - | 3,623 | 5,301 |
| loc-Brightkite | root | 45,279 | 463.9 | 5.0 | 115 | 2,479 | - | 3,700 | 5,699 |
| soc-Epinions1 | own | 49,098 | 922.5 | 7.7 | 228 | 2,577 | - | 6,105 | 9,227 |
| soc-Epinions1 | half | 57,388 | 973.7 | 7.3 | 249 | 2,950 | - | 7,101 | 8,534 |
| soc-Epinions1 | root | 61,009 | 1149.8 | 4.4 | 288 | 2,918 | - | 7,002 | 8,985 |
| soc-Slashdot0902 | own | 47,604 | 469.5 | 8.0 | 139 | 2,704 | - | 5,211 | 6,646 |
| soc-Slashdot0902 | half | 59,057 | 478.6 | 8.4 | 156 | 2,717 | - | 4,725 | 7,174 |
| soc-Slashdot0902 | root | 70,421 | 459.1 | 6.7 | 139 | 2,854 | - | 4,522 | 7,638 |
| com-youtube | own | 774,495 | 959.7 | 12.3 | 292 | 32,014 | - | 26,084 | 58,157 |
| com-youtube | half | 936,661 | 993.2 | 17.6 | 302 | 38,678 | - | 31,353 | 68,615 |
| com-youtube | root | 1,005,937 | 985.7 | 4.5 | 277 | 43,574 | - | 27,305 | 65,871 |
| soc-pokec | own | 827,163 | 40378.2 | 8.2 | 12,447 | 88,853 | - | 344,700 | 505,276 |
| soc-pokec | half | 992,602 | 44650.3 | 37.1 | 13,897 | 102,828 | - | 374,935 | 751,648 |
| soc-pokec | root | 1,240,689 | 51695.1 | 7.2 | 13,230 | 129,016 | - | 423,324 | 811,181 |
| ca-HepTh | own | 3,957 | 38.2 | 14.8 | 24 | 286 | - | 356 | 605 |
| ca-HepTh | half | 5,149 | 44.2 | 13.1 | 40 | 357 | - | 400 | 649 |
| ca-HepTh | root | 5,654 | 52.1 | 5.9 | 33 | 378 | - | 423 | 760 |
| email-Eu-core | own | 647 | 285.1 | 6.7 | 92 | 286 | - | 640 | 759 |
| email-Eu-core | half | 679 | 287.6 | 40.4 | 121 | 338 | - | 689 | 809 |
| email-Eu-core | root | 744 | 276.8 | 7.2 | 96 | 312 | - | 676 | 656 |
| amazon-copurchase | own | 123,947 | 627.0 | 21.7 | 241 | 6,949 | - | 16,745 | 24,462 |
| amazon-copurchase | half | 179,284 | 949.6 | 15.9 | 376 | 11,506 | - | 23,367 | 32,476 |
| amazon-copurchase | root | 187,527 | 1154.4 | 8.5 | 403 | 12,421 | - | 23,446 | 32,605 |
| dblp-coauthor | own | 1,907,463 | 64066.9 | 10.6 | 20,800 | 211,961 | - | 570,612 | 463,097 |
| dblp-coauthor | half | 2,188,825 | 68647.1 | 25.8 | 21,831 | 231,834 | - | 512,946 | 494,870 |
| dblp-coauthor | root | 2,750,068 | 78730.2 | 8.2 | 23,975 | 273,154 | - | 511,610 | 591,312 |

"locate" returns the head range, a pointer to the whole runs and the
tail range (no copy). "ranges copied" materialises the range list.
"explicit ids" writes every vertex id. The build-form columns are the
same index before `compact_runs` (chain-id arrays, tops and residues as
T). The memcpy column is the stage-2 per-vertex S-tree listing, available
for the five stage-2 graphs. All compact-form numbers are measured on the
index loaded back from its file (per-size widths, `CHAINX03`).

### 9.4 Value queries (ns per query)
| Graph | value | max tree depth |
|---|---:|---:|
| ca-GrQc | 23.9 | 26 |
| ca-HepPh | 11.9 | 131 |
| com-dblp | 16.6 | 157 |
| web-Stanford | 19.0 | 2138 |
| amazon0302 | 10.7 | 10 |
| ca-AstroPh | 6.6 | 571 |
| ca-CondMat | 10.1 | 105 |
| cit-HepPh | 7.1 | 1062 |
| loc-Brightkite | 9.7 | 730 |
| soc-Epinions1 | 7.4 | 1614 |
| soc-Slashdot0902 | 7.2 | 587 |
| com-youtube | 8.7 | 1604 |
| soc-pokec | 14.6 | 2413 |
| ca-HepTh | 9.5 | 18 |
| email-Eu-core | 10.5 | 325 |
| amazon-copurchase | 6.1 | 10 |
| dblp-coauthor | 7.7 | 5543 |

### 9.5 Explicit listing cost per output vertex (ns)
| Graph | final module (own) | per-vertex S trees (own) | final module (root) | per-vertex S trees (root) |
|---|---:|---:|---:|---:|
| ca-GrQc | 0.106 | 0.136 | 0.098 | 0.090 |
| ca-HepPh | 0.128 | 0.117 | 0.104 | 0.112 |
| com-dblp | 0.076 | 0.143 | 0.117 | 0.132 |
| web-Stanford | 0.076 | 0.145 | 0.149 | 0.153 |
| amazon0302 | 0.050 | 0.133 | 0.055 | 0.139 |
| ca-AstroPh | 0.129 | - | 0.108 | - |
| ca-CondMat | 0.086 | - | 0.085 | - |
| cit-HepPh | 0.204 | - | 0.195 | - |
| loc-Brightkite | 0.056 | - | 0.055 | - |
| soc-Epinions1 | 0.052 | - | 0.048 | - |
| soc-Slashdot0902 | 0.057 | - | 0.041 | - |
| com-youtube | 0.041 | - | 0.043 | - |
| soc-pokec | 0.107 | - | 0.104 | - |
| ca-HepTh | 0.072 | - | 0.067 | - |
| email-Eu-core | 0.443 | - | 0.419 | - |
| amazon-copurchase | 0.056 | - | 0.066 | - |
| dblp-coauthor | 0.111 | - | 0.099 | - |

Selftest of the module that produced this evidence (both builds, Release
and ASan/UBSan, identical counts): 34,075 graphs, 4,844,536 community
queries, 5,123,368 value checks, 29,215,144 membership checks and
2,772,392 ladder checks across the exact and the double form. Byte ratio
over the 17 graphs: min 1.28x, median 4.95x, max 11.97x.

### 9.6 Top encoding ablation (same process)
| Graph | bytes T tops | bytes packed | climb own T / packed | climb half T / packed | climb root T / packed | locate own T / packed | value T / packed |
|---|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 68,989 | 61,328 | 13.6 / 16.6 | 17.6 / 21.8 | 6.8 / 5.8 | 21.8 / 24.0 | 20.5 / 23.9 |
| ca-HepPh | 199,645 | 193,082 | 13.2 / 16.3 | 27.8 / 17.1 | 7.2 / 7.5 | 22.4 / 14.4 | 6.2 / 11.9 |
| com-dblp | 1,756,069 | 1,641,813 | 15.2 / 18.2 | 16.7 / 31.8 | 8.6 / 8.9 | 20.1 / 27.6 | 12.0 / 16.6 |
| web-Stanford | 3,708,306 | 3,552,706 | 14.2 / 18.0 | 49.0 / 55.1 | 12.3 / 12.3 | 22.6 / 44.6 | 20.1 / 19.0 |
| amazon0302 | 3,537,806 | 3,081,868 | 6.9 / 7.7 | 9.4 / 7.7 | 4.7 / 5.2 | 10.8 / 15.2 | 11.6 / 10.7 |
| ca-AstroPh | 636,199 | 617,232 | 5.9 / 8.1 | 15.3 / 17.7 | 5.8 / 4.7 | 8.1 / 5.6 | 6.2 / 6.6 |
| ca-CondMat | 196,752 | 177,662 | 11.4 / 8.2 | 15.9 / 20.6 | 8.3 / 6.6 | 14.8 / 8.2 | 9.5 / 10.1 |
| cit-HepPh | 1,656,743 | 1,610,959 | 4.9 / 4.8 | 32.8 / 24.9 | 5.7 / 6.3 | 6.0 / 5.1 | 8.4 / 7.1 |
| loc-Brightkite | 850,977 | 811,068 | 8.2 / 7.2 | 12.5 / 9.9 | 2.8 / 5.2 | 9.1 / 11.5 | 9.2 / 9.7 |
| soc-Epinions1 | 1,178,357 | 1,117,624 | 10.3 / 7.7 | 11.0 / 11.1 | 5.7 / 3.2 | 11.0 / 7.7 | 10.3 / 7.4 |
| soc-Slashdot0902 | 629,427 | 601,127 | 9.8 / 7.8 | 8.2 / 8.4 | 5.4 / 4.3 | 10.4 / 8.0 | 8.5 / 7.2 |
| com-youtube | 3,237,905 | 3,103,325 | 6.8 / 9.2 | 5.9 / 6.6 | 4.4 / 4.8 | 8.0 / 12.3 | 8.8 / 8.7 |
| soc-pokec | 32,681,536 | 32,391,133 | 5.3 / 9.8 | 37.8 / 23.6 | 9.4 / 15.8 | 15.3 / 8.2 | 12.0 / 14.6 |
| ca-HepTh | 85,683 | 76,147 | 11.8 / 9.4 | 9.8 / 11.1 | 5.5 / 5.2 | 12.9 / 14.8 | 7.6 / 9.5 |
| email-Eu-core | 158,767 | 145,590 | 6.7 / 7.1 | 31.6 / 29.8 | 7.3 / 7.2 | 8.0 / 6.7 | 7.2 / 10.5 |
| amazon-copurchase | 4,582,929 | 3,975,343 | 13.0 / 13.1 | 14.3 / 10.3 | 7.2 / 6.3 | 20.4 / 21.7 | 9.2 / 6.1 |
| dblp-coauthor | 45,076,102 | 44,284,853 | 8.0 / 9.9 | 16.0 / 37.4 | 5.4 / 6.3 | 10.7 / 10.6 | 16.7 / 7.7 |

### 9.7 Climb only
| Graph | own build / packed | half build / packed | root build / packed |
|---|---:|---:|---:|
| ca-GrQc | 14.8 / 16.6 | 13.1 / 21.8 | 6.3 / 5.8 |
| ca-HepPh | 13.0 / 16.3 | 21.9 / 17.1 | 7.7 / 7.5 |
| com-dblp | 13.3 / 18.2 | 24.2 / 31.8 | 6.0 / 8.9 |
| web-Stanford | 10.0 / 18.0 | 31.9 / 55.1 | 6.0 / 12.3 |
| amazon0302 | 10.0 / 7.7 | 8.8 / 7.7 | 6.0 / 5.2 |
| ca-AstroPh | 5.8 / 8.1 | 20.0 / 17.7 | 4.6 / 4.7 |
| ca-CondMat | 12.0 / 8.2 | 26.5 / 20.6 | 5.2 / 6.6 |
| cit-HepPh | 5.2 / 4.8 | 30.8 / 24.9 | 4.1 / 6.3 |
| loc-Brightkite | 9.6 / 7.2 | 6.1 / 9.9 | 4.2 / 5.2 |
| soc-Epinions1 | 10.0 / 7.7 | 7.1 / 11.1 | 4.8 / 3.2 |
| soc-Slashdot0902 | 9.7 / 7.8 | 7.1 / 8.4 | 4.6 / 4.3 |
| com-youtube | 10.1 / 9.2 | 10.4 / 6.6 | 4.5 / 4.8 |
| soc-pokec | 7.8 / 9.8 | 24.7 / 23.6 | 5.0 / 15.8 |
| ca-HepTh | 8.3 / 9.4 | 5.9 / 11.1 | 3.3 / 5.2 |
| email-Eu-core | 7.1 / 7.1 | 42.9 / 29.8 | 7.8 / 7.2 |
| amazon-copurchase | 11.7 / 13.1 | 20.9 / 10.3 | 4.7 / 6.3 |
| dblp-coauthor | 13.6 / 9.9 | 34.2 / 37.4 | 7.6 / 6.3 |

### 9.8 The four stage-2 layouts for reference (five graphs)

#### Bytes with values (Block D), stage-2 layouts
| Graph | per-vertex S trees | over twins | over chains | aligned chains (index.cpp) | aligned vs per-vertex |
|---|---:|---:|---:|---:|---:|
| ca-GrQc | 272,564 | 256,238 | 103,220 | 65,128 | 4.19x |
| ca-HepPh | 1,899,628 | 1,532,406 | 395,336 | 306,065 | 6.21x |
| com-dblp | 16,710,168 | 16,880,402 | 4,447,786 | 2,024,433 | 8.25x |
| web-Stanford | 24,853,426 | 26,086,394 | 5,762,930 | 3,632,412 | 6.84x |
| amazon0302 | 14,414,214 | 17,079,748 | 5,744,638 | 3,860,362 | 3.73x |

#### Community listing, own level (k = kappa_s(v)), ns per query, explicit ids unless noted
| Graph | output vertices | per-vertex S trees (memcpy) | over twins | over chains | aligned explicit | aligned ranges only |
|---|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 1,656 | 225 | 2,903 | 574 | 822 | 231 |
| ca-HepPh | 5,332 | 625 | 9,668 | 2,116 | 1,844 | 867 |
| com-dblp | 164,595 | 23,493 | 370,861 | 42,952 | 47,205 | 8,694 |
| web-Stanford | 82,833 | 11,983 | 163,655 | 27,140 | 30,132 | 7,507 |
| amazon0302 | 107,061 | 14,277 | 223,981 | 133,725 | 129,013 | 20,811 |

## 10. Analysis Of Runtime And Memory

Measured facts.
- Bytes: 1.28x (email-Eu-core) to 11.97x (dblp-coauthor) below
  per-vertex S trees, median 4.95x over the seventeen graphs; 1.25x to
  8.76x if the label permutation is charged to the index. The ratio
  follows the vertex collapse n / chains: 1.5 vertices per chain on
  email-Eu-core (1,005 vertices, nearly every vertex its own chain), 2.1
  on cit-HepPh, 4.3 on pokec (ratios 1.3x, 2.0x, 3.2x) against 11.1 on
  dblp-coauthor, 23.6 on dblp and 26.5 on youtube (11.7x, 10.3x, 11.2x).
  The largest input, dblp-coauthor (4.05 M vertices, cliques up to 450,
  512-bit solver counts), is also the best case: 44 MB against 530 MB. Where the
  bytes sit changes with the graph: the per-size layers are 58-74 percent
  of the index on the collaboration, web and product graphs (HepPh 74,
  GrQc 67, dblp 66), the per-chain block (trajectories and residue
  values) 44-58 percent on the social and citation graphs (pokec 58,
  cit-HepPh 51, youtube 46, Slashdot 44); the map is 3-13 percent.
- Chains against canonical nodes: 0.26-0.70 on twelve graphs, but 1.08
  (Slashdot), 1.32 (dblp-coauthor), 1.77 (cit-HepPh), 1.91 (youtube) and
  5.67 (pokec). The recombination of CHAINS.md Section 7 is real on
  social, citation and dense collaboration graphs; the index stays
  smaller than per-vertex S trees there because the per-vertex trees pay
  per vertex and per size, and chains still collapse 2-27 vertices each.
- Runs against (chain, size) pairs: 1.4x-3.3x fewer entries. With
  per-size widths the compact form is now 0-31 percent smaller than the
  chain-id build form (the 8-byte entry per node is paid back by the
  narrower tops and residues) and answers faster: range answers 2.4x-109x,
  explicit answers 1.55x-5.5x.
- Locating a community takes 3.5-59 ns on every graph and regime: bitmap
  rank, trajectory lookup, jump-pointer climb, two entry reads. Copying the
  range list takes 30 ns to 23 us (dblp-coauthor root, 78,730 ranges); the
  range list is 2.3x (email-Eu-core) to 807x (youtube) shorter than the
  vertex list.
- Explicit listing runs at 0.04-0.11 ns per vertex on fourteen graphs
  (dblp-coauthor included: 1.9 M vertices per own-level community in
  167 us) and 0.17-0.77 ns on the three most fragmented (AstroPh 9
  vertices per range, cit-HepPh 5, email-Eu-core 2.3). Against the memcpy baseline (0.09-0.15 ns per vertex on
  the five graphs where it was measured) the fill is faster on 13 of 15
  (graph, regime) points: own level 1.24x (GrQc), 1.17x (HepPh), 1.75x
  (dblp), 1.54x (Stanford), 1.68x (amazon); root 1.37x-1.80x on the three
  large graphs, 0.86x and 0.95x on GrQc and HepPh. Per output vertex the
  fill reads 8 bytes per range instead of 4 bytes per vertex, and the
  branchless eight-wide stores remove the per-range branch that made the
  scalar-tail fill of `archive/final_v2.json` 1.3x-2.4x slower on the fragmented
  graphs.
- Top encoding (Section 9.6, same process): with tops packed to per-size
  widths the climb is 0.91x (own) and 0.89x (root) of the T-tops climb in
  the median over the thirteen graphs, the spread (0.4x-1.5x) being the
  run-to-run noise of identical code on this machine; residue packing has
  no visible cost. Against the fixed-width run of `archive/final_v3.json` the
  median ratios of every latency class lie between 0.92 and 1.12 with
  both signs, so the widths change bytes, not time.
- Values 6.1-24 ns: one bitmap rank plus a residue or binomial lookup.
- Build: 15 ms (GrQc) to 34 s (pokec, 1.63 M vertices) and 189 s
  (dblp-coauthor, 4.05 M vertices, 512-bit counts) including the solver's
  row index (pokec: row index 14 s, all-size peel 9 s, trees and trie 11 s;
  dblp-coauthor: 50 s, 91 s, 48 s); chains and layout add at most 0.4 s.
  Peak build memory is 157 MB on dblp, 2.9 GB on pokec and 8.8 GB on
  dblp-coauthor (Section 14), against 1.6, 32 and 45 MB indexes.

Hypotheses (not measured): a prefetch of the next run could shave the
remaining per-range constant on AstroPh and cit-HepPh; a per-node vertex
count (4 bytes per node) would give community sizes in O(1) without
summing runs.

The permutation, computed (not a hypothesis): perm[v] = start_pos[chain(v)]
plus the number of earlier input vertices in the same chain, so it is a
rank query over the chain-id sequence in input order and compresses to
n ceil(log2 C) bits (a wavelet tree): dblp 0.55 MB instead of 1.27 MB,
youtube 2.27 MB instead of 4.54 MB, pokec 3.9 MB instead of 6.5 MB. The
byte ratios against per-vertex S trees would move from 1.39x-4.98x
(4 n permutation) to 1.43x-6.20x (compressed) against 1.47x-8.42x
(aligned storage, no permutation). Each lookup would then cost about
log2 C bit-vector ranks (50-100 ns), more than the 3-51 ns community
locate itself, so the compressed form is a fallback for label-bound
deployments, not the default; storing the graph in the index's order
(this codebase already relabels by degeneracy order at load) removes the
cost entirely.

Algorithmic versus engineering gains: the chain partition and the
lexicographic rank order (one range at s = 2, few ranges elsewhere) are
algorithmic; the run array with entry points is a layout theorem
(CHAINS.md C6) implemented as a constant number of array reads; the
vectorised fill and the Myers skip pointers are engineering.

## 11. Failure Cases / Limitations

- The label permutation is not free: 4 n bytes unless the graph is stored
  in aligned order (60 percent of the dblp index, 113 percent on youtube
  where the index is tiny against n). Reported both ways.
- When almost every vertex is its own chain the index degenerates toward
  per-vertex S trees: email-Eu-core (1,005 vertices, 667 chains, 2.3
  vertices per range) is 1.28x smaller, and its explicit listing costs
  0.77 ns per vertex. The saving is the vertex collapse; without it only
  the run encoding and the certified tail remain.
- Explicit listing costs 0.20-0.24 ns per vertex on the two most
  fragmented graphs (AstroPh, cit-HepPh), about twice the cost elsewhere;
  no memcpy baseline was measured on them. On HepPh (own and half) the
  fill is 5-10 percent slower than the memcpy listing.
- Peak build memory is the solver's row index plus the graph (2.9 GB on
  pokec, of which the row index is 1.46 GB), not the index; the dense
  core matrix and the per-size own arrays of the first version are gone
  (Section 14). With 512-bit counts (dblp-coauthor) the solver's own
  per-size weight arrays (two counts per clique-tree row) add about 4 GB
  on top of the 1.7 GB row index, for 8.8 GB in total; that is inside
  `r1_terminal_20260918`, not in this module.
- The number of chains exceeds the number of canonical nodes on four of
  thirteen graphs (up to 5.7x on pokec); no bound in terms of the
  hierarchy exists (CHAINS.md Section 7).
- The memcpy baseline was measured only on the five stage-2 graphs, in
  the stage-2 binary with its own query draw (same seed, protocol and
  machine); output sizes differ by 0.1-0.3 percent between the draws.
- One machine, seventeen graphs up to 4.05 M vertices, in-memory
  single-process timing; no cold-cache or multi-process protocol; no
  comparison with a reimplemented SGL (its inputs are bipartite).
- com-lj (4.0 M vertices, 34.7 M edges) does not build here: the
  solver's clique-tree row index of `r1_terminal_20260918` stores member
  ids in 32 bits and its `append` stops with "member ID overflow" after
  208 s and 10.8 GB. com-orkut (117 M edges) was not attempted. Both need
  either 64-bit member ids in that solver or a server with more memory;
  neither is a limit of the index itself (dblp-coauthor, larger in
  vertices and far denser, builds in 189 s).

## 12. Final Conclusion

The all-size r = 1 nucleus hierarchy of a graph is stored (exactly below
2^53, to the nearest double above) in 1.3x-12.0x fewer bytes than one S
tree per size (median 4.95x over seventeen graphs of five families; dblp
10.2x, the 4.05 M-vertex dblp-coauthor 12.0x, the 1,005-vertex
email-Eu-core 1.3x), communities
are located in constant time (3.5-59 ns) and listed faster than a memcpy
of the vertex list on 13 of the 15 points where that baseline exists
(0.86x-1.80x), from a partition of the vertices
(chains) that the theory proves exact, aligned labels that make chains id
ranges, per-size run arrays with node entry points that make every
community a head range, whole runs and a tail range, and per-size byte
widths for the stored values. The build streams the solver's rows and
needs O(n + chains x sizes) working memory beyond the solver's own row
index. Correctness is brute-force verified in four forms on 34,075 graphs
plus K_300 (sizes beyond 255) under Release and sanitizers, and the canonical node counts agree with the independent
stage-1/2 pipeline on the five shared graphs. The SGL skyline dedup is
dominated in this setting and dropped. Framing for a paper: the
contribution starts from the chain partition (hierarchy equivalence),
not from storing S trees; S trees over vertices are the baseline it
beats on both axes; for r >= 2 the paper's size forest already plays the
role of the chain structure (CHAINS.md Section 10).

## 13. Next Recommended Improvements

1. Per-node vertex counts (4 bytes per node) for O(1) community sizes
   without summing runs.
2. The branchless fill is in (it cut the per-vertex listing cost 1.3x
   to 1.9x on every graph tried against the scalar-tail fill of
   `archive/final_v2.json`); a prefetch of the next run could shave the remaining
   per-range constant on the fragmented graphs.
3. Store the graph in the aligned order and drop the permutation; the
   compressed wavelet-tree form (numbers in Section 10) is only for
   deployments that must keep input labels and can pay 50-100 ns per
   lookup.
4. Stream the build size by size to cut the 5 GB peak on pokec-sized
   inputs; then com-lj and com-orkut on the server.
5. A bound on the number of chains for graphs with monotone core
   trajectories, or a construction showing none exists beyond
   CHAINS.md Section 7; the four graphs with chains > nodes are the
   test cases.
6. r >= 2 is closed (CHAINS.md Section 10 verdict): no separate chain
   index; the unifying remark belongs in the paper's discussion.
7. com-lj and com-orkut: widen the solver's member ids to 64 bits (or
   split the row index) and run on a server; with 64-bit counts their
   row indexes will be several GB.

## 14. Build Memory And Per-Size Widths (added later on 2026-09-19)

Two changes after the 13-graph run above, both in `chain_index_tool.cpp`
and `chain_index.hpp`; the solver gained an optional row sink
(`research/r1_terminal_20260918/terminal.hpp`, default path unchanged,
its own selftest still passes). Evidence: `buildonly.json`,
`buildonly-logs/` (one `--build-only` process per graph, no query passes).

**Build memory.** The all-size solver now streams one core row at a time
(a two-row window instead of the s_max x n matrix), each row becomes its
canonical tree immediately, and chains are refined size by size in a
prefix trie (one node per distinct own-node prefix) instead of keeping
one own-node array per size; the count width is derived from the
solver's own rows (with an overflow retry) instead of the joint clique
path index that the earlier tool built only for that purpose (1.7 GB and
about 40 s on pokec). The produced index is byte-identical (GrQc, HepPh
sha256 against `final.json`; every selftest form passes). What remains
is the solver's shared row index (`ti`), the graph, and O(n + pairs)
working arrays.

| Graph | row index ms | solve ms | trees + trie ms | chains ms | layout ms | compact ms | total ms | total before ms | peak RSS MB after | before | row index MB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 2 | 2 | 4 | 0 | 0 | 0.0 | 8 | 14 | 5 | 14 | 0 |
| ca-HepPh | 71 | 168 | 138 | 0 | 2 | 0.5 | 379 | 503 | 29 | 178 | 3 |
| com-dblp | 334 | 253 | 505 | 2 | 5 | 0.4 | 1,100 | 1,371 | 157 | 899 | 28 |
| web-Stanford | 2,024 | 1,126 | 1,927 | 4 | 15 | 2.2 | 5,098 | 4,866 | 299 | 869 | 115 |
| amazon0302 | 237 | 155 | 391 | 2 | 17 | 1.6 | 804 | 801 | 134 | 205 | 28 |
| ca-AstroPh | 205 | 159 | 250 | 1 | 14 | 0.5 | 630 | 595 | 28 | 54 | 8 |
| ca-CondMat | 52 | 26 | 44 | 0 | 1 | 0.1 | 124 | 106 | 12 | 26 | 2 |
| cit-HepPh | 879 | 285 | 359 | 1 | 8 | 0.9 | 1,532 | 1,650 | 95 | 220 | 42 |
| loc-Brightkite | 350 | 423 | 690 | 1 | 11 | 3.9 | 1,479 | 1,753 | 95 | 238 | 41 |
| soc-Epinions1 | 4,305 | 1,519 | 3,128 | 1 | 7 | 0.7 | 8,960 | 8,167 | 448 | 1,124 | 213 |
| soc-Slashdot0902 | 1,426 | 1,054 | 2,105 | 1 | 4 | 0.4 | 4,590 | 3,753 | 257 | 634 | 107 |
| com-youtube | 4,652 | 1,037 | 2,469 | 5 | 16 | 1.8 | 8,180 | 7,283 | 557 | 1,663 | 210 |
| soc-pokec | 37,275 | 12,907 | 16,545 | 14 | 211 | 28.0 | 66,980 | 63,386 | 2,889 | 6,100 | 1,455 |

"total before" is `build_ms` of `final.json`, which included the joint
Layout build; times are on the loaded laptop (1-minute load 10-13 while
this sweep ran) and agree with the earlier run within noise, so the build
is not slower; the memory columns are load-independent.

**Per-size widths.** Node tops and residue values are stored at the width
of the widest value of their size, rounded to 1, 2, 4, 8, 16 or 32 bytes
(file format `CHAINX03`); a field of w bytes at index y sits at offset
y w, aligned, so every decode is one constant-size load, dispatched once
per query. Residue packing has no hot-path cost (value queries equal or
faster in every in-process comparison). Top packing puts one width
dispatch (a 5-way jump) in front of each climb: when consecutive queries
alternate between sizes of different widths that jump mispredicts, which
costs about 1-3 ns per query in the worst case and nothing when the size
repeats. In-process two-round comparisons of "compact, tops as T" against
"compact, packed tops" on GrQc, HepPh, amazon and dblp (fields `full_*`
against the unprefixed ones in the bench JSON) are in Section 9.6 for
all thirteen graphs: packed tops climb in 0.91x (own) and 0.89x (root) of
the T-tops time in the median, with a 0.4x-1.5x spread that matches the
run-to-run noise of identical code on this machine (1-minute load 4-7
during the run, after the browser was closed). Packing is on by default
and switchable (`compact_runs(false)` keeps tops as T; the file records
the choice).

| Graph | index B before | index B after | change | vs per-vertex S trees | map B | chains B | layers B |
|---|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 68,560 | 59,896 | -12.6% | 4.55x | 3,856 | 15,889 | 40,151 |
| ca-HepPh | 333,072 | 215,578 | -35.3% | 8.81x | 6,808 | 49,845 | 158,925 |
| com-dblp | 2,138,946 | 1,626,879 | -23.9% | 10.27x | 113,304 | 433,259 | 1,080,316 |
| web-Stanford | 4,086,330 | 3,516,780 | -13.9% | 7.07x | 124,720 | 1,279,808 | 2,112,252 |
| amazon0302 | 3,859,326 | 3,000,134 | -22.3% | 4.80x | 212,628 | 1,040,604 | 1,746,902 |
| ca-AstroPh | 696,534 | 610,570 | -12.3% | 3.63x | 16,860 | 225,357 | 368,353 |
| ca-CondMat | 208,494 | 173,828 | -16.6% | 7.47x | 12,020 | 57,438 | 104,370 |
| cit-HepPh | 2,170,902 | 1,577,713 | -27.3% | 2.03x | 72,980 | 798,209 | 706,524 |
| loc-Brightkite | 977,340 | 798,416 | -18.3% | 3.29x | 36,232 | 306,945 | 455,239 |
| soc-Epinions1 | 1,428,086 | 1,099,982 | -23.0% | 2.93x | 49,524 | 438,839 | 611,619 |
| soc-Slashdot0902 | 764,640 | 587,531 | -23.2% | 5.10x | 42,608 | 255,723 | 289,200 |
| com-youtube | 4,021,798 | 3,017,695 | -25.0% | 11.22x | 384,064 | 1,392,191 | 1,241,440 |
| soc-pokec | 43,724,688 | 31,624,721 | -27.7% | 3.22x | 1,838,988 | 18,481,448 | 11,304,285 |

The ratio against per-vertex S trees is 1.28x (email-Eu-core) to 11.71x
(dblp-coauthor), median 4.95x over seventeen graphs; dblp 10.27x. Section
9 holds the latencies of this version (`final.json`, `more.json`); the
fixed-width run is `archive/final_v3.json`.

**Sizes beyond 255.** The first version kept omega, sigma and the trie
level in one byte and failed on dblp-coauthor (degeneracy 449, s_max 450)
with "trie path length". They are 16-bit now (2 more bytes per chain,
under 0.3 percent of any index); the selftest gained K_300 (one chain,
every size certified, one range per community, values C(299, s-1) up to
512 bits).

**Values as double (user decision, precision not required).** The stored
value type is `double` (format `CHAINX05`; the header records 0 for
double, else the integer width, so a file cannot be read with the wrong
type). The solver still peels with exact integers (the bucket order needs
exact comparisons); the conversion happens once when a value is written
into the index. Per-size widths keep integers below 2^32 at 1, 2 or 4
bytes, so on the fourteen graphs whose counts fit 64 bits nothing changes
(byte-identical value blocks); the 128/256/512-bit graphs lose the wide
fields: HepPh 215,578 -> 193,082 B (-10 percent), dblp-coauthor
45.3 -> 44.3 MB (-2 percent, its wide values were already confined to a
few sizes). Below 2^53 every stored value is exact. Above it a value is
the nearest double; two consecutive levels of a merge tree would have to
agree to 53 bits to be confused, which does not happen on any input
measured (K_300, values up to C(299, 149) ~ 10^88, reproduces every level
and community; the double binomial table agrees with the exact one to
1e-12 relative). Query latencies are unchanged within noise (interleaved
A/B on GrQc, HepPh, dblp, amazon: climbs within +-5 ns, both signs). The
exact integer form remains available as `ChainIndex<T>` with the solver's
count type and is what the selftest runs alongside the double form
(4.84 M community queries in total). The
earlier latency evidence was measured on the index compacted in place
rather than on the loaded copy (a reference bound once before the target
pointer changed); the two hold identical arrays, so those numbers stand,
and the bench now measures the loaded copy.

## 15. Server Runs (tods1, tods2; added 2026-09-19 evening)

Both UNSW servers (Ubuntu 22.04, GCC 11.4, 96 cores, 503 GB; one thread
per run, other users' load 1-8 on 96 cores) built the same sources
(commit 0c7e43a and later), passed the brute-force selftest with the
same counts as the laptop, and ran `run_final.py` on every graph present
there: tods1 20 graphs (`tods1.json`, `tods1-logs/`), tods2 7 graphs
(`tods2.json`, `tods2-logs/`); a second tods1 run covers the three
largest graphs after the 64-bit fix below (`tods1_big.json`). Launchers:
`tods1_run.sh`, `tods1_big.sh`, `tods2_run.sh`.

Portability findings. (1) The index does not depend on the input's edge
order: com-dblp, web-Stanford, com-amazon and web-Google exist on both
servers as files with different sha256 (different edge orders) and give
byte-identical indexes. (2) Laptop and server indexes of the same graph
differ by at most 0.2 percent in bytes (GrQc 61,328 against 61,440): the
per-size trees are built from unstable sorts whose tie order differs
between libc++ and libstdc++, which changes node preorder ids, hence
chain ranks, hence how many runs merge; every form is checked against
brute force on both platforms. (3) Single-thread query latencies on the
servers are 1.5-3x those of the Apple M-series laptop (locate 14-35 ns,
listing 0.10-0.31 ns per vertex, values 15-35 ns), consistent with the
CPUs; ratios between configurations are the same.

The 32-bit member limit. com-lj (4.0 M vertices, 34.7 M edges),
ca-hollywood-2009 (1.07 M vertices, 56.3 M edges, cliques up to 2,209)
and com-orkut (3.07 M vertices, 117 M edges) stopped in the solver's
clique-tree row index with "member ID overflow": the number of (row,
vertex) incidences exceeds 2^32. The offsets into the member array are
64-bit as of commit 0c7e43a (`Row::begin/hold_end/pivot_end/end`,
choice offsets; solver and index selftests pass, dblp index
byte-identical); the second tods1 run uses that build. Results are
appended to the merged table when they finish. hollywood additionally
needs core values up to about 10^660 (C(2208, 1104)), beyond 512-bit
integers and beyond double (10^308); the exact solver will report the
count bound as exceeded, and a long-double or log-domain count type is
the only way to cover such cliques.

New graphs and what they say. web-uk-2005 (cliques up to 500, 155
vertices per chain) and ca-coauthors-dblp (cliques up to 337, 58 per
chain) give the largest ratios so far, 48.7x and 45.8x: large cliques
make many vertices hierarchy-equivalent. web-NotreDame 24.7x and
web-BerkStan 19.2x show the same on web graphs with cliques of 150-200.
cit-Patents (3.8 M vertices) is 12.1x, so the poor 2.0x of cit-HepPh is
not a property of citation graphs but of that small dense one.
wiki-Talk (2.4 M vertices, 81 per chain) is 15.9x. The two slowest
builds, tech-as-skitter (928 s, 13.9 GB) and wiki-Talk (734 s, 20 GB),
spend that time and memory in the solver's row index, not in the index
(28 MB and 4 MB).

### All inputs, one row per distinct graph (29 graphs; laptop = Apple M-series, tods1/tods2 = Ubuntu 22.04, GCC 11, 96 cores, 503 GB; one thread everywhere)
| Graph | machine | n | s_max | W bits | chains | n / chains | index B | per-vertex S trees B | ratio | locate own ns | list ns per vertex | value ns | build s | peak RSS MB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| email-Eu-core | tods1 | 1,005 | 35 | 64 | 667 | 1.5 | 145,590 | 186,738 | 1.28x | 17.0 | 0.823 | 18.5 | 0.2 | 16 |
| cit-HepPh | laptop | 34,546 | 31 | 64 | 16,623 | 2.1 | 1,610,959 | 3,199,324 | 1.99x | 5.1 | 0.204 | 7.1 | 0.9 | 0 |
| soc-Epinions1 | tods1 | 75,879 | 68 | 64 | 8,821 | 8.6 | 1,117,672 | 3,221,518 | 2.88x | 17.1 | 0.190 | 18.0 | 8.9 | 311 |
| soc-pokec-relationships | tods1 | 1,632,803 | 48 | 64 | 383,206 | 4.3 | 32,391,325 | 101,770,646 | 3.14x | 27.0 | 0.313 | 24.8 | 73.3 | 2,807 |
| loc-Brightkite | laptop | 58,228 | 53 | 64 | 6,326 | 9.2 | 811,068 | 2,628,872 | 3.24x | 11.5 | 0.056 | 9.7 | 1.1 | 0 |
| ca-AstroPh | tods1 | 18,772 | 57 | 64 | 3,331 | 5.6 | 617,240 | 2,216,476 | 3.59x | 22.6 | 0.295 | 14.6 | 0.6 | 24 |
| tech-as-skitter | tods1 | 1,694,616 | 112 | 128 | 163,053 | 10.4 | 28,075,132 | 123,683,600 | 4.41x | 22.9 | 0.190 | 22.4 | 928.4 | 13,874 |
| ca-GrQc | tods1 | 5,242 | 44 | 64 | 716 | 7.3 | 61,440 | 272,564 | 4.44x | 27.8 | 0.141 | 18.4 | 0.0 | 10 |
| amazon0302 | laptop | 262,111 | 7 | 64 | 40,867 | 6.4 | 3,081,868 | 14,414,214 | 4.68x | 15.2 | 0.050 | 10.7 | 0.5 | 0 |
| com-amazon.ungraph | tods1 | 334,863 | 7 | 64 | 46,023 | 7.3 | 3,263,670 | 15,545,122 | 4.76x | 31.1 | 0.112 | 19.9 | 0.9 | 115 |
| amazon-copurchase | laptop | 548,552 | 7 | 64 | 56,315 | 9.7 | 3,975,343 | 19,697,700 | 4.95x | 21.7 | 0.056 | 6.1 | 0.5 | 0 |
| soc-Slashdot0902 | laptop | 82,168 | 56 | 64 | 6,798 | 12.1 | 601,127 | 2,995,392 | 4.98x | 8.0 | 0.057 | 7.2 | 2.2 | 0 |
| ca-HepTh | laptop | 9,877 | 32 | 64 | 1,104 | 8.9 | 76,147 | 422,086 | 5.54x | 14.8 | 0.072 | 9.5 | 0.0 | 0 |
| web-Google | tods1 | 875,713 | 45 | 64 | 73,836 | 11.9 | 11,156,521 | 63,649,906 | 5.71x | 31.0 | 0.161 | 21.5 | 6.3 | 355 |
| web-Stanford | tods1 | 281,903 | 72 | 64 | 17,963 | 15.7 | 3,552,754 | 24,853,426 | 7.00x | 24.3 | 0.138 | 21.8 | 5.1 | 228 |
| ca-CondMat | tods1 | 23,133 | 26 | 64 | 1,917 | 12.1 | 177,574 | 1,299,206 | 7.32x | 20.6 | 0.169 | 14.4 | 0.1 | 14 |
| ca-HepPh | tods1 | 12,008 | 239 | 256 | 1,136 | 10.6 | 193,074 | 1,899,628 | 9.84x | 22.1 | 0.128 | 20.1 | 0.7 | 18 |
| ca-MathSciNet | tods1 | 332,689 | 25 | 64 | 15,869 | 21.0 | 1,234,250 | 12,339,606 | 10.00x | 23.2 | 0.099 | 16.9 | 0.8 | 101 |
| com-dblp | tods1 | 317,080 | 114 | 128 | 13,459 | 23.6 | 1,641,981 | 16,710,168 | 10.18x | 24.3 | 0.188 | 15.8 | 2.1 | 113 |
| com-youtube | tods1 | 1,134,890 | 52 | 64 | 42,815 | 26.5 | 3,104,309 | 33,847,412 | 10.90x | 18.8 | 0.200 | 17.1 | 7.5 | 494 |
| dblp-core30 | tods1 | 1,206 | 114 | 128 | 35 | 34.5 | 50,114 | 549,992 | 10.97x | 25.1 | 0.772 | 28.9 | 0.0 | 10 |
| dblp-coauthor | laptop | 4,049,537 | 450 | 512 | 363,201 | 11.1 | 44,284,853 | 530,275,990 | 11.97x | 10.6 | 0.111 | 7.7 | 221.0 | 0 |
| cit-Patents | tods2 | 3,774,768 | 65 | 64 | 203,819 | 18.5 | 12,079,564 | 146,173,956 | 12.10x | 23.8 | 0.157 | 22.0 | 19.5 | 2,005 |
| web-it-2004 | tods1 | 509,338 | 432 | 512 | 35,402 | 14.4 | 9,866,036 | 135,833,300 | 13.77x | 34.9 | 0.202 | 20.0 | 66.5 | 524 |
| wiki-Talk | tods1 | 2,394,385 | 132 | 64 | 29,699 | 80.6 | 3,964,777 | 63,192,530 | 15.94x | 14.0 | 0.195 | 15.1 | 733.6 | 20,141 |
| web-BerkStan | tods2 | 685,230 | 202 | 256 | 39,559 | 17.3 | 8,607,977 | 165,295,380 | 19.20x | 27.0 | 0.213 | 24.0 | 65.6 | 1,069 |
| web-NotreDame | tods2 | 325,729 | 156 | 256 | 9,122 | 35.7 | 1,007,187 | 24,840,466 | 24.66x | 19.0 | 0.175 | 18.4 | 4.2 | 189 |
| ca-coauthors-dblp | tods1 | 540,486 | 337 | 512 | 9,326 | 58.0 | 4,619,939 | 211,523,396 | 45.78x | 23.1 | 0.195 | 34.6 | 106.8 | 876 |
| web-uk-2005 | tods1 | 129,632 | 500 | 512 | 839 | 154.5 | 4,093,454 | 199,316,948 | 48.69x | 21.5 | 0.193 | 144.5 | 69.5 | 296 |

byte ratio over 29 distinct graphs: min 1.28x, median 7.00x, max 48.69x
failed on tods1: com-lj: member ID overflow
failed on tods1: ca-hollywood-2009: member ID overflow
failed on tods1: com-orkut: failed
