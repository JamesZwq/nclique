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
- Q2 the community of v at (s, k) as a vertex set,
- Q3 membership: is u in the community of v at (s, k),
- Q4 the ladder of v at size s: every nucleus containing v with its size.
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
runs, tail range. Q3 = climb from v's own node, then a subtree interval
test on u's own node id. Q4 = walk parents, summing run lengths.

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
  O(1) to locate, O(ranges) to copy, O(vertices) to expand; Q3 O(climb);
  Q4 O(depth x runs per level).
The count of chains is bounded by n and by the twin classes; it is below
N_T on every input measured but not by a theorem (CHAINS.md Section 7).

## 6. Implementation Summary

- `chain_index.hpp` (about 180 lines): `ChainIndex<T>` with the blocks
  above, `chain_of`, `own_node`, `climb`, `value`, `community_runs`
  (pointer form), `community_ranges` (vector form, both forms),
  `expand` (branchless fill: every block of eight ids is stored
  unconditionally with 128-bit vector stores and the pointer advances by
  the true length, so the caller's buffer carries eight spare slots),
  `member`, `ladder`, `compact_runs`, `save`, `load`, byte counters.
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
dblp-coauthor) in `more.json` (sha256 of every input recorded). Count
width chosen from the solver's rows: 64 bits except com-dblp (128),
ca-HepPh (256, s_max 239) and dblp-coauthor (512, s_max 450). Queries (seed 20260918, drawn from active vertices
and sizes 2 <= s <= omega(v)): community regimes own (k = kappa_s(v),
20,000), half (k = max(1, kappa/2), 20,000), root (k = 1, 1,000);
membership 20,001 (random u, mixed regimes); values 200,000; ladders on
the own set. One warm-up pass plus five timed passes, median reported.
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
| ca-GrQc | 5,242 | 44 | 64 | 716 | 7.3 | 1,517 | 0.47 | 2,079 | 624 | 3,856 | 15,889 | 40,151 | 59,896 | 55,704 | 20,968 | 65,472 | 272,564 | 4.55x | 3.37x |
| ca-HepPh | 12,008 | 239 | 256 | 1,136 | 10.6 | 3,538 | 0.32 | 8,468 | 4,969 | 6,808 | 49,845 | 158,925 | 215,578 | 211,297 | 48,032 | 311,136 | 1,899,628 | 8.81x | 7.21x |
| com-dblp | 317,080 | 114 | 128 | 13,459 | 23.6 | 33,979 | 0.40 | 70,208 | 28,784 | 113,304 | 433,259 | 1,080,316 | 1,626,879 | 1,495,709 | 1,268,320 | 2,052,686 | 16,710,168 | 10.27x | 5.77x |
| web-Stanford | 281,903 | 72 | 64 | 17,963 | 15.7 | 53,239 | 0.34 | 151,364 | 97,074 | 124,720 | 1,279,808 | 2,112,252 | 3,516,780 | 3,306,848 | 1,127,612 | 3,701,670 | 24,853,426 | 7.07x | 5.35x |
| amazon0302 | 262,111 | 7 | 64 | 40,867 | 6.4 | 65,134 | 0.63 | 143,567 | 47,380 | 212,628 | 1,040,604 | 1,746,902 | 3,000,134 | 2,739,957 | 1,048,444 | 3,793,970 | 14,414,214 | 4.80x | 3.56x |
| ca-AstroPh | 18,772 | 57 | 64 | 3,331 | 5.6 | 6,462 | 0.52 | 36,381 | 25,742 | 16,860 | 225,357 | 368,353 | 610,570 | 587,131 | 75,088 | 609,826 | 2,216,476 | 3.63x | 3.23x |
| ca-CondMat | 23,133 | 26 | 64 | 1,917 | 12.1 | 3,357 | 0.57 | 8,693 | 3,658 | 12,020 | 57,438 | 104,370 | 173,828 | 161,538 | 92,532 | 200,374 | 1,299,206 | 7.47x | 4.88x |
| cit-HepPh | 34,546 | 31 | 64 | 16,623 | 2.1 | 9,369 | 1.77 | 104,308 | 61,217 | 72,980 | 798,209 | 706,524 | 1,577,713 | 1,541,580 | 138,184 | 2,060,682 | 3,199,324 | 2.03x | 1.86x |
| loc-Brightkite | 58,228 | 53 | 64 | 6,326 | 9.2 | 11,547 | 0.55 | 35,833 | 21,427 | 36,232 | 306,945 | 455,239 | 798,416 | 754,473 | 232,912 | 902,652 | 2,628,872 | 3.29x | 2.55x |
| soc-Epinions1 | 75,879 | 68 | 64 | 8,821 | 8.6 | 12,886 | 0.68 | 55,344 | 38,876 | 49,524 | 438,839 | 611,619 | 1,099,982 | 1,051,298 | 303,516 | 1,286,374 | 3,221,518 | 2.93x | 2.30x |
| soc-Slashdot0902 | 82,168 | 56 | 64 | 6,798 | 12.1 | 6,309 | 1.08 | 30,496 | 17,551 | 42,608 | 255,723 | 289,200 | 587,531 | 564,663 | 328,672 | 720,540 | 2,995,392 | 5.10x | 3.27x |
| com-youtube | 1,134,890 | 52 | 64 | 42,815 | 26.5 | 22,403 | 1.91 | 179,601 | 93,541 | 384,064 | 1,392,191 | 1,241,440 | 3,017,695 | 2,930,287 | 4,539,560 | 3,901,854 | 33,847,412 | 11.22x | 4.48x |
| soc-pokec | 1,632,803 | 48 | 64 | 383,206 | 4.3 | 67,544 | 5.67 | 2,375,352 | 1,212,885 | 1,838,988 | 18,481,448 | 11,304,285 | 31,624,721 | 31,356,585 | 6,531,212 | 43,252,464 | 101,770,646 | 3.22x | 2.67x |
| ca-HepTh | 9,877 | 32 | 64 | 1,104 | 8.9 | 1,583 | 0.70 | 3,112 | 1,044 | 6,284 | 26,475 | 43,388 | 76,147 | 71,199 | 39,508 | 87,244 | 422,086 | 5.54x | 3.65x |
| email-Eu-core | 1,005 | 35 | 64 | 667 | 1.5 | 2,602 | 0.26 | 6,080 | 4,423 | 2,868 | 47,387 | 95,335 | 145,590 | 136,689 | 4,020 | 166,432 | 186,738 | 1.28x | 1.25x |
| amazon-copurchase | 548,552 | 7 | 64 | 56,315 | 9.7 | 86,798 | 0.65 | 163,658 | 54,803 | 328,132 | 1,385,981 | 2,261,230 | 3,975,343 | 3,628,510 | 2,194,208 | 4,840,472 | 19,697,700 | 4.95x | 3.19x |
| dblp-coauthor | 4,049,537 | 450 | 512 | 363,201 | 11.1 | 274,215 | 1.32 | 2,958,470 | 1,761,549 | 2,212,112 | 21,189,882 | 21,883,355 | 45,285,349 | 44,207,011 | 16,198,148 | 140,276,804 | 530,275,990 | 11.71x | 8.62x |

"total B" is in memory and includes the derived jump pointers (4 bytes
per node); "file B" is the disk image. "perm B" is the 4 n-byte
permutation from input labels to aligned labels, needed only if the graph
is not stored in the aligned order. "build form total B" is the same index
before `compact_runs` (chain-id DFS arrays).

### 9.2 Build, save, load
| Graph | solve (all-size peel) | trees | chains + labels | layout | build total | compact | save | load | process wall s | peak RSS MB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 3 | 8 | 0 | 1 | 11 | 0.08 | 1.6 | 0.1 | 0.38 | 13 |
| ca-HepPh | 174 | 131 | 0 | 1 | 312 | 0.69 | 0.9 | 0.5 | 1.71 | 37 |
| com-dblp | 292 | 472 | 2 | 9 | 779 | 0.79 | 9.1 | 0.9 | 18.74 | 140 |
| web-Stanford | 937 | 1,385 | 3 | 21 | 2,349 | 3.63 | 6.7 | 2.2 | 30.26 | 319 |
| amazon0302 | 121 | 326 | 2 | 15 | 466 | 1.84 | 2.3 | 0.9 | 25.27 | 131 |
| ca-AstroPh | 105 | 179 | 0 | 4 | 289 | 0.22 | 3.0 | 0.3 | 9.18 | 32 |
| ca-CondMat | 14 | 30 | 0 | 1 | 45 | 0.15 | 2.6 | 0.2 | 1.74 | 14 |
| cit-HepPh | 309 | 381 | 1 | 7 | 698 | 0.60 | 4.0 | 0.6 | 37.39 | 99 |
| loc-Brightkite | 359 | 605 | 1 | 4 | 969 | 0.27 | 5.1 | 0.4 | 6.93 | 89 |
| soc-Epinions1 | 958 | 1,860 | 1 | 6 | 2,826 | 0.62 | 3.1 | 0.5 | 19.44 | 495 |
| soc-Slashdot0902 | 527 | 1,001 | 1 | 2 | 1,531 | 0.22 | 2.7 | 0.2 | 8.95 | 258 |
| com-youtube | 485 | 910 | 4 | 10 | 1,410 | 1.23 | 4.3 | 0.7 | 44.20 | 553 |
| soc-pokec | 6,687 | 8,698 | 12 | 113 | 15,541 | 16.25 | 27.2 | 7.0 | 495.31 | 3,229 |
| ca-HepTh | 4 | 8 | 0 | 0 | 13 | 0.06 | 0.4 | 0.2 | 0.75 | 11 |
| email-Eu-core | 41 | 48 | 0 | 1 | 91 | 0.16 | 0.6 | 0.2 | 4.88 | 19 |
| amazon-copurchase | 213 | 640 | 5 | 29 | 888 | 2.84 | 11.1 | 1.9 | 32.05 | 160 |
| dblp-coauthor | 91,055 | 47,728 | 50 | 359 | 139,340 | 52.64 | 62.9 | 10.6 | 989.67 | 8,816 |

Process wall time and peak RSS are for the whole `--bench` run (graph
load, the solver's row index, build, three forms, all query passes);
"build total" excludes the row index, whose time is the `ti_ms` field of
`final.json` (Section 14 itemises both). The peak is dominated by the
solver's row index plus the graph (Section 14).

### 9.3 Community queries (ns per query; compact form loaded from disk)
| Graph | regime | output vertices | ranges | locate (pointer) | ranges copied | explicit ids | per-vertex S trees memcpy | build form ranges | build form explicit |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | own | 1,658 | 14.1 | 18.7 | 43 | 181 | 225 | 228 | 288 |
| ca-GrQc | half | 2,108 | 15.8 | 25.3 | 61 | 230 | 305 | 145 | 480 |
| ca-GrQc | root | 2,208 | 18.6 | 8.7 | 30 | 238 | 206 | 118 | 542 |
| ca-HepPh | own | 5,327 | 129.8 | 13.8 | 94 | 533 | 625 | 383 | 1,101 |
| ca-HepPh | half | 6,257 | 137.0 | 30.9 | 102 | 710 | 832 | 1,111 | 1,980 |
| ca-HepPh | root | 7,439 | 169.4 | 12.6 | 95 | 839 | 794 | 1,157 | 1,988 |
| com-dblp | own | 164,832 | 871.7 | 24.2 | 399 | 13,421 | 23,493 | 11,068 | 24,068 |
| com-dblp | half | 199,796 | 951.0 | 29.4 | 429 | 17,343 | 26,144 | 11,373 | 27,780 |
| com-dblp | root | 229,952 | 1236.8 | 8.4 | 557 | 20,913 | 30,259 | 14,061 | 32,508 |
| web-Stanford | own | 83,227 | 1217.7 | 6.5 | 859 | 7,800 | 11,983 | 9,125 | 16,595 |
| web-Stanford | half | 102,997 | 1339.6 | 20.9 | 891 | 10,256 | 12,157 | 8,821 | 20,578 |
| web-Stanford | root | 115,934 | 1772.8 | 6.0 | 710 | 10,705 | 19,274 | 9,819 | 18,977 |
| amazon0302 | own | 105,997 | 734.6 | 28.9 | 335 | 8,520 | 14,277 | 19,294 | 29,894 |
| amazon0302 | half | 151,391 | 1283.2 | 38.4 | 562 | 13,921 | 24,757 | 29,015 | 41,943 |
| amazon0302 | root | 160,229 | 1529.7 | 9.1 | 613 | 16,791 | 23,022 | 30,443 | 48,376 |
| ca-AstroPh | own | 8,854 | 960.8 | 9.0 | 375 | 1,503 | - | 2,954 | 4,614 |
| ca-AstroPh | half | 9,701 | 984.4 | 44.6 | 405 | 1,521 | - | 3,028 | 4,633 |
| ca-AstroPh | root | 11,775 | 1140.8 | 10.4 | 446 | 2,491 | - | 3,710 | 5,075 |
| ca-CondMat | own | 9,999 | 165.0 | 16.9 | 97 | 1,110 | - | 1,097 | 1,908 |
| ca-CondMat | half | 12,008 | 181.1 | 26.5 | 113 | 1,330 | - | 1,204 | 2,270 |
| ca-CondMat | root | 14,389 | 214.9 | 9.9 | 106 | 1,570 | - | 722 | 2,597 |
| cit-HepPh | own | 18,252 | 3693.7 | 10.4 | 1,851 | 5,079 | - | 16,839 | 21,592 |
| cit-HepPh | half | 21,348 | 4107.8 | 40.0 | 2,008 | 6,698 | - | 19,426 | 23,891 |
| cit-HepPh | root | 26,461 | 4941.0 | 10.7 | 2,200 | 9,043 | - | 21,200 | 28,590 |
| loc-Brightkite | own | 32,323 | 400.7 | 7.3 | 168 | 2,491 | - | 5,317 | 7,621 |
| loc-Brightkite | half | 39,732 | 415.3 | 9.3 | 182 | 2,847 | - | 5,538 | 8,547 |
| loc-Brightkite | root | 45,279 | 463.9 | 4.3 | 185 | 2,078 | - | 4,642 | 9,412 |
| soc-Epinions1 | own | 49,098 | 922.5 | 7.1 | 382 | 3,751 | - | 8,835 | 12,791 |
| soc-Epinions1 | half | 57,388 | 973.7 | 8.2 | 408 | 4,589 | - | 9,187 | 13,786 |
| soc-Epinions1 | root | 61,009 | 1149.8 | 3.5 | 483 | 2,796 | - | 10,189 | 13,659 |
| soc-Slashdot0902 | own | 47,604 | 469.5 | 8.2 | 102 | 2,054 | - | 4,068 | 8,816 |
| soc-Slashdot0902 | half | 59,057 | 478.6 | 6.8 | 110 | 3,508 | - | 7,074 | 9,073 |
| soc-Slashdot0902 | root | 70,421 | 459.1 | 4.4 | 101 | 2,680 | - | 4,780 | 10,085 |
| com-youtube | own | 774,495 | 959.7 | 8.1 | 281 | 33,848 | - | 27,297 | 59,264 |
| com-youtube | half | 936,661 | 993.2 | 6.6 | 288 | 39,096 | - | 28,422 | 67,235 |
| com-youtube | root | 1,005,937 | 985.7 | 5.1 | 264 | 43,134 | - | 28,789 | 69,581 |
| soc-pokec | own | 827,163 | 40378.2 | 9.6 | 9,731 | 76,321 | - | 346,782 | 416,385 |
| soc-pokec | half | 992,602 | 44650.3 | 27.0 | 13,033 | 87,064 | - | 374,750 | 449,683 |
| soc-pokec | root | 1,240,689 | 51695.1 | 8.5 | 13,302 | 105,102 | - | 405,789 | 503,011 |
| ca-HepTh | own | 3,957 | 38.2 | 17.9 | 49 | 416 | - | 532 | 889 |
| ca-HepTh | half | 5,149 | 44.2 | 27.4 | 54 | 457 | - | 622 | 1,246 |
| ca-HepTh | root | 5,654 | 52.1 | 7.0 | 43 | 493 | - | 633 | 1,399 |
| email-Eu-core | own | 647 | 285.1 | 12.5 | 179 | 496 | - | 980 | 1,539 |
| email-Eu-core | half | 679 | 287.6 | 46.1 | 186 | 527 | - | 1,188 | 1,365 |
| email-Eu-core | root | 744 | 276.8 | 12.1 | 154 | 495 | - | 993 | 1,342 |
| amazon-copurchase | own | 123,947 | 627.0 | 32.4 | 767 | 10,605 | - | 24,358 | 36,395 |
| amazon-copurchase | half | 179,284 | 949.6 | 59.1 | 829 | 14,281 | - | 35,093 | 56,646 |
| amazon-copurchase | root | 187,527 | 1154.4 | 16.6 | 1,053 | 17,359 | - | 39,934 | 60,956 |
| dblp-coauthor | own | 1,907,463 | 64066.9 | 14.2 | 14,987 | 166,958 | - | 299,530 | 743,364 |
| dblp-coauthor | half | 2,188,825 | 68647.1 | 24.0 | 21,907 | 195,011 | - | 331,545 | 924,366 |
| dblp-coauthor | root | 2,750,068 | 78730.2 | 6.7 | 23,250 | 222,717 | - | 429,993 | 864,781 |

"locate" returns the head range, a pointer to the whole runs and the
tail range (no copy). "ranges copied" materialises the range list.
"explicit ids" writes every vertex id. The build-form columns are the
same index before `compact_runs` (chain-id arrays, tops and residues as
T). The memcpy column is the stage-2 per-vertex S-tree listing, available
for the five stage-2 graphs. All compact-form numbers are measured on the
index loaded back from its file (per-size widths, `CHAINX03`).

### 9.4 Membership, values, ladders (ns per query)
| Graph | member | value | ladder (compact) | ladder steps | ladder (build form) | max depth |
|---|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 20.2 | 13.2 | 27 | 2.08 | 162 | 26 |
| ca-HepPh | 27.0 | 12.1 | 135 | 7.28 | 1,401 | 131 |
| com-dblp | 15.3 | 12.7 | 848 | 5.05 | 13,044 | 157 |
| web-Stanford | 20.1 | 20.2 | 20,851 | 61.20 | 87,207 | 2138 |
| amazon0302 | 27.5 | 16.5 | 778 | 3.22 | 40,651 | 10 |
| ca-AstroPh | 28.3 | 11.5 | 7,523 | 42.85 | 32,041 | 571 |
| ca-CondMat | 14.9 | 12.5 | 328 | 7.10 | 2,249 | 105 |
| cit-HepPh | 28.1 | 15.2 | 24,349 | 48.10 | 140,178 | 1062 |
| loc-Brightkite | 4.9 | 7.8 | 502 | 7.23 | 7,177 | 730 |
| soc-Epinions1 | 6.9 | 7.1 | 9,475 | 17.89 | 35,034 | 1614 |
| soc-Slashdot0902 | 9.5 | 6.3 | 606 | 6.46 | 9,299 | 587 |
| com-youtube | 8.5 | 8.2 | 2,343 | 4.12 | 29,040 | 1604 |
| soc-pokec | 11.0 | 11.2 | 113,261 | 23.04 | 1,637,301 | 2413 |
| ca-HepTh | 17.6 | 13.0 | 63 | 2.35 | 666 | 18 |
| email-Eu-core | 24.7 | 12.8 | 4,328 | 66.31 | 24,216 | 325 |
| amazon-copurchase | 30.9 | 17.6 | 941 | 2.50 | 35,952 | 10 |
| dblp-coauthor | 18.7 | 10.2 | 437,715 | 90.89 | 2,510,380 | 5543 |

### 9.5 Explicit listing cost per output vertex (ns)
| Graph | final module (own) | per-vertex S trees (own) | final module (root) | per-vertex S trees (root) |
|---|---:|---:|---:|---:|
| ca-GrQc | 0.109 | 0.136 | 0.108 | 0.090 |
| ca-HepPh | 0.100 | 0.117 | 0.113 | 0.112 |
| com-dblp | 0.081 | 0.143 | 0.091 | 0.132 |
| web-Stanford | 0.094 | 0.145 | 0.092 | 0.153 |
| amazon0302 | 0.080 | 0.133 | 0.105 | 0.139 |
| ca-AstroPh | 0.170 | - | 0.212 | - |
| ca-CondMat | 0.111 | - | 0.109 | - |
| cit-HepPh | 0.278 | - | 0.342 | - |
| loc-Brightkite | 0.077 | - | 0.046 | - |
| soc-Epinions1 | 0.076 | - | 0.046 | - |
| soc-Slashdot0902 | 0.043 | - | 0.038 | - |
| com-youtube | 0.044 | - | 0.043 | - |
| soc-pokec | 0.092 | - | 0.085 | - |
| ca-HepTh | 0.105 | - | 0.087 | - |
| email-Eu-core | 0.768 | - | 0.665 | - |
| amazon-copurchase | 0.086 | - | 0.093 | - |
| dblp-coauthor | 0.088 | - | 0.081 | - |

selftests (final.json): {"build": {"passed": true, "graphs": 34075, "community_queries": 2422268, "membership_checks": 14607572, "value_checks": 2561684, "ladder_checks": 1386196}, "build-asan": {"passed": true, "graphs": 34075, "community_queries": 2422268, "membership_checks": 14607572, "value_checks": 2561684, "ladder_checks": 1386196}}
selftests (more.json): {"build": {"passed": true, "graphs": 34075, "community_queries": 2422268, "membership_checks": 14607572, "value_checks": 2561684, "ladder_checks": 1386196}, "build-asan": {"passed": true, "graphs": 34075, "community_queries": 2422268, "membership_checks": 14607572, "value_checks": 2561684, "ladder_checks": 1386196}}
byte ratio over 17 graphs: min 1.28x, median 4.95x, max 11.71x

selftests: {"build": {"passed": true, "graphs": 34075, "community_queries": 2422268, "membership_checks": 14607572, "value_checks": 2561684, "ladder_checks": 1386196}, "build-asan": {"passed": true, "graphs": 34075, "community_queries": 2422268, "membership_checks": 14607572, "value_checks": 2561684, "ladder_checks": 1386196}}

selftests: {"build": {"passed": true, "graphs": 34075, "community_queries": 1816701, "membership_checks": 10955679, "value_checks": 1921263, "ladder_checks": 1039647}, "build-asan": {"passed": true, "graphs": 34075, "community_queries": 1816701, "membership_checks": 10955679, "value_checks": 1921263, "ladder_checks": 1039647}}

selftests: {"build": {"passed": true, "graphs": 34075, "community_queries": 1816701, "membership_checks": 10955679, "value_checks": 1921263, "ladder_checks": 1039647}, "build-asan": {"passed": true, "graphs": 34075, "community_queries": 1816701, "membership_checks": 10955679, "value_checks": 1921263, "ladder_checks": 1039647}}

### 9.6 Top encoding ablation (same process)
| Graph | bytes T tops | bytes packed | climb own T / packed | climb half T / packed | climb root T / packed | locate own T / packed | member T / packed | value T / packed | ladder T / packed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 67,557 | 59,896 | 13.1 / 14.6 | 12.3 / 16.2 | 6.9 / 7.5 | 22.2 / 18.7 | 17.0 / 20.2 | 12.5 / 13.2 | 51 / 27 |
| ca-HepPh | 282,285 | 215,578 | 16.1 / 13.9 | 34.7 / 27.2 | 12.6 / 10.5 | 22.4 / 13.8 | 26.9 / 27.0 | 12.2 / 12.1 | 241 / 135 |
| com-dblp | 2,000,983 | 1,626,879 | 17.4 / 17.3 | 12.7 / 21.8 | 8.2 / 9.2 | 23.5 / 24.2 | 18.7 / 15.3 | 12.2 / 12.7 | 878 / 848 |
| web-Stanford | 3,672,380 | 3,516,780 | 11.3 / 5.6 | 19.0 / 20.1 | 9.6 / 10.9 | 17.2 / 6.5 | 11.6 / 20.1 | 17.7 / 20.2 | 20,058 / 20,851 |
| amazon0302 | 3,456,072 | 3,000,134 | 14.7 / 15.0 | 13.5 / 37.1 | 9.0 / 10.7 | 27.5 / 28.9 | 17.1 / 27.5 | 14.9 / 16.5 | 451 / 778 |
| ca-AstroPh | 629,537 | 610,570 | 11.3 / 9.8 | 16.8 / 32.1 | 11.7 / 10.4 | 14.5 / 9.0 | 12.3 / 28.3 | 11.6 / 11.5 | 9,211 / 7,523 |
| ca-CondMat | 192,918 | 173,828 | 14.9 / 13.5 | 10.8 / 18.1 | 8.5 / 9.1 | 7.8 / 16.9 | 11.8 / 14.9 | 11.9 / 12.5 | 324 / 328 |
| cit-HepPh | 1,623,497 | 1,577,713 | 8.3 / 9.9 | 30.0 / 42.3 | 10.8 / 9.5 | 10.4 / 10.4 | 21.1 / 28.1 | 15.9 / 15.2 | 24,006 / 24,349 |
| loc-Brightkite | 838,325 | 798,416 | 17.7 / 7.0 | 11.1 / 6.5 | 7.2 / 3.8 | 20.1 / 7.3 | 10.8 / 4.9 | 17.1 / 7.8 | 928 / 502 |
| soc-Epinions1 | 1,160,715 | 1,099,982 | 13.3 / 6.3 | 9.3 / 5.7 | 6.2 / 3.9 | 15.1 / 7.1 | 11.1 / 6.9 | 12.2 / 7.1 | 9,963 / 9,475 |
| soc-Slashdot0902 | 615,831 | 587,531 | 10.3 / 7.9 | 6.7 / 6.6 | 4.5 / 3.9 | 11.1 / 8.2 | 10.1 / 9.5 | 9.2 / 6.3 | 606 / 606 |
| com-youtube | 3,152,275 | 3,017,695 | 8.6 / 9.9 | 9.5 / 4.5 | 3.7 / 3.2 | 9.5 / 8.1 | 9.9 / 8.5 | 8.0 / 8.2 | 2,304 / 2,343 |
| soc-pokec | 31,915,124 | 31,624,721 | 5.2 / 7.9 | 31.7 / 20.3 | 6.3 / 7.8 | 13.6 / 9.6 | 9.7 / 11.0 | 12.3 / 11.2 | 112,989 / 113,261 |
| ca-HepTh | 85,683 | 76,147 | 17.5 / 14.1 | 10.5 / 12.7 | 7.2 / 6.6 | 19.6 / 17.9 | 10.3 / 17.6 | 13.2 / 13.0 | 70 / 63 |
| email-Eu-core | 158,767 | 145,590 | 4.0 / 11.0 | 23.6 / 42.5 | 10.2 / 10.2 | 8.9 / 12.5 | 18.6 / 24.7 | 13.7 / 12.8 | 4,107 / 4,328 |
| amazon-copurchase | 4,582,929 | 3,975,343 | 18.1 / 20.4 | 37.0 / 41.1 | 11.3 / 10.2 | 37.8 / 32.4 | 57.2 / 30.9 | 29.0 / 17.6 | 647 / 941 |
| dblp-coauthor | 60,532,038 | 45,285,349 | 11.1 / 6.3 | 42.1 / 20.6 | 9.2 / 5.8 | 15.7 / 14.2 | 17.7 / 18.7 | 10.0 / 10.2 | 421,504 / 437,715 |

### 9.7 Climb only
| Graph | own build / packed | half build / packed | root build / packed |
|---|---:|---:|---:|
| ca-GrQc | 16.5 / 14.6 | 8.0 / 16.2 | 6.0 / 7.5 |
| ca-HepPh | 16.7 / 13.9 | 19.4 / 27.2 | 10.8 / 10.5 |
| com-dblp | 17.6 / 17.3 | 13.3 / 21.8 | 10.2 / 9.2 |
| web-Stanford | 13.0 / 5.6 | 20.9 / 20.1 | 10.5 / 10.9 |
| amazon0302 | 8.2 / 15.0 | 32.3 / 37.1 | 10.9 / 10.7 |
| ca-AstroPh | 11.3 / 9.8 | 16.9 / 32.1 | 10.4 / 10.4 |
| ca-CondMat | 14.7 / 13.5 | 10.2 / 18.1 | 7.5 / 9.1 |
| cit-HepPh | 7.1 / 9.9 | 25.4 / 42.3 | 9.1 / 9.5 |
| loc-Brightkite | 15.3 / 7.0 | 8.4 / 6.5 | 5.8 / 3.8 |
| soc-Epinions1 | 12.6 / 6.3 | 8.2 / 5.7 | 5.7 / 3.9 |
| soc-Slashdot0902 | 10.5 / 7.9 | 10.1 / 6.6 | 4.2 / 3.9 |
| com-youtube | 8.6 / 9.9 | 5.4 / 4.5 | 3.8 / 3.2 |
| soc-pokec | 6.5 / 7.9 | 14.7 / 20.3 | 6.6 / 7.8 |
| ca-HepTh | 14.2 / 14.1 | 14.2 / 12.7 | 6.8 / 6.6 |
| email-Eu-core | 8.5 / 11.0 | 38.0 / 42.5 | 9.2 / 10.2 |
| amazon-copurchase | 16.6 / 20.4 | 14.9 / 41.1 | 5.9 / 10.2 |
| dblp-coauthor | 10.8 / 6.3 | 27.0 / 20.6 | 9.2 / 5.8 |

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
- Bytes: 1.28x (email-Eu-core) to 11.71x (dblp-coauthor) below
  per-vertex S trees, median 4.95x over the seventeen graphs; 1.25x to
  8.62x if the label permutation is charged to the index. The ratio
  follows the vertex collapse n / chains: 1.5 vertices per chain on
  email-Eu-core (1,005 vertices, nearly every vertex its own chain), 2.1
  on cit-HepPh, 4.3 on pokec (ratios 1.3x, 2.0x, 3.2x) against 11.1 on
  dblp-coauthor, 23.6 on dblp and 26.5 on youtube (11.7x, 10.3x, 11.2x).
  The largest input, dblp-coauthor (4.05 M vertices, cliques up to 450,
  512-bit counts), is also the best case: 45 MB against 530 MB. Where the
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
  explicit answers 1.55x-5.5x, ladders 3.7x-52x.
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
- Membership 4.9-31 ns and values 6.3-20 ns: two bitmap ranks plus a
  climb, or one rank plus a residue or binomial lookup.
- Build: 15 ms (GrQc) to 34 s (pokec, 1.63 M vertices) and 189 s
  (dblp-coauthor, 4.05 M vertices, 512-bit counts) including the solver's
  row index (pokec: row index 14 s, all-size peel 9 s, trees and trie 11 s;
  dblp-coauthor: 50 s, 91 s, 48 s); chains and layout add at most 0.4 s.
  Peak build memory is 157 MB on dblp, 2.9 GB on pokec and 8.8 GB on
  dblp-coauthor (Section 14), against 1.6, 32 and 45 MB indexes.

Hypotheses (not measured): a per-node vertex count (4 bytes per node)
would make ladders O(depth) (pokec 132 us, cit-HepPh 23 us today); the
build's dense core matrix and per-size own arrays can be streamed size by
size; a prefetch of the next run could shave the remaining per-range
constant on AstroPh and cit-HepPh.

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
- Ladder queries on deep trees (dblp-coauthor depth 5,543, pokec 2,413,
  cit-HepPh 1,062) cost 24-438 us because each level sums its runs.
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

The all-size r = 1 nucleus hierarchy of a graph is stored exactly in
1.3x-11.7x fewer bytes than one S tree per size (median 4.95x over
seventeen graphs of five families; dblp 10.3x, the 4.05 M-vertex
dblp-coauthor 11.7x, the 1,005-vertex email-Eu-core 1.3x), communities
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

1. Per-node vertex counts (4 bytes per node) for O(depth) ladders and
   O(1) community sizes.
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
