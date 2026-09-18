# The Chain Index: Final Module, Measured

Date: 2026-09-19. Module `chain_index.hpp` (header-only, `namespace
chainindex`), tool `chain_index_tool.cpp` (`--selftest`, `--bench`),
driver `run_final.py`. Evidence: `final.json`, `final-logs/` (one
`/usr/bin/time -l` log per graph, build and selftest logs), index files
`cx/<graph>.cx` (not committed). Tables are printed by `report_tables.py`
from `final.json` and the stage-2 files `index*.json`. Theory:
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
[RESULTS_INDEX.md](RESULTS_INDEX.md), `index_vertices.json`). It lists a
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
  `expand` (eight ids per step with 128-bit vector stores, scalar tail),
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
thread (`OMP_NUM_THREADS=1`), no other timing run concurrent. Thirteen
inputs: the five of stage 2 (`data/ca-GrQc.edges`, `data/ca-HepPh.edges`,
`data/com-dblp.edges`, `graphs/web-Stanford.edges`,
`graphs/amazon0302.edges`) plus `graphs/ca-AstroPh.edges`,
`graphs/ca-CondMat.edges`, `graphs/cit-HepPh.edges`,
`graphs/loc-Brightkite.edges`, `graphs/soc-Epinions1.edges`,
`graphs/soc-Slashdot0902.edges`, `graphs/com-youtube.edges`,
`graphs/soc-pokec.edges` (sha256 of every input in `final.json`). Count
width chosen by `count_bound`: 64 bits except ca-HepPh (256, s_max 239)
and com-dblp (128). Queries (seed 20260918, drawn from active vertices
and sizes 2 <= s <= omega(v)): community regimes own (k = kappa_s(v),
20,000), half (k = max(1, kappa/2), 20,000), root (k = 1, 1,000);
membership 20,001 (random u, mixed regimes); values 200,000; ladders on
the own set. One warm-up pass plus five timed passes, median reported.
Explicit ids are written into a preallocated caller buffer. The
per-vertex S-tree bytes are computed by the tool with the stage-2
`vertices` accounting (verified equal to `index_vertices.json` on the
five stage-2 graphs, 16,710,168 bytes on dblp); the memcpy listing
baseline exists only for those five graphs (stage-2 run, same protocol,
its own query draw). The first five-graph run of this module, before the
baseline field and the eight extra graphs were added, is kept as
`final_v1.json` / `final-logs_v1/`; its numbers agree with the tables
below within run-to-run noise.

Commands:
```
cmake -S research/r1_skyline_index_20260918 -B research/r1_skyline_index_20260918/build -DCMAKE_BUILD_TYPE=Release -DSANITIZE=OFF
cmake --build research/r1_skyline_index_20260918/build -j 12 --target chain_index_tool
research/r1_skyline_index_20260918/build/chain_index_tool --selftest
/usr/bin/time -l research/r1_skyline_index_20260918/build/chain_index_tool --bench data/com-dblp.edges research/r1_skyline_index_20260918/cx/com-dblp.cx
python3 research/r1_skyline_index_20260918/run_final.py graphs/ca-AstroPh.edges graphs/ca-CondMat.edges graphs/cit-HepPh.edges graphs/loc-Brightkite.edges graphs/soc-Epinions1.edges graphs/soc-Slashdot0902.edges graphs/com-youtube.edges graphs/soc-pokec.edges
python3 research/r1_skyline_index_20260918/report_tables.py
```

## 9. Experimental Results

### 9.1 Size
| Graph | n | s_max | W bits | chains | n / chains | canonical nodes | chains / nodes | (chain,s) pairs | runs | map B | chains B | layers B | total B | file B | perm B | build form total B | per-vertex S trees B | ratio | ratio with perm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 5,242 | 44 | 64 | 716 | 7.3 | 1,517 | 0.47 | 2,079 | 624 | 3,856 | 16,892 | 47,812 | 68,560 | 64,308 | 20,968 | 65,472 | 272,564 | 3.98x | 3.04x |
| ca-HepPh | 12,008 | 239 | 256 | 1,136 | 10.6 | 3,538 | 0.32 | 8,468 | 4,969 | 6,808 | 100,632 | 225,632 | 333,072 | 328,536 | 48,032 | 311,136 | 1,899,628 | 5.70x | 4.98x |
| com-dblp | 317,080 | 114 | 128 | 13,459 | 23.6 | 33,979 | 0.40 | 70,208 | 28,784 | 113,304 | 571,222 | 1,454,420 | 2,138,946 | 2,007,646 | 1,268,320 | 2,052,686 | 16,710,168 | 7.81x | 4.90x |
| web-Stanford | 281,903 | 72 | 64 | 17,963 | 15.7 | 53,239 | 0.34 | 151,364 | 97,074 | 124,720 | 1,693,758 | 2,267,852 | 4,086,330 | 3,876,310 | 1,127,612 | 3,701,670 | 24,853,426 | 6.08x | 4.77x |
| amazon0302 | 262,111 | 7 | 64 | 40,867 | 6.4 | 65,134 | 0.63 | 143,567 | 47,380 | 212,628 | 1,443,858 | 2,202,840 | 3,859,326 | 3,599,126 | 1,048,444 | 3,793,970 | 14,414,214 | 3.73x | 2.94x |
| ca-AstroPh | 18,772 | 57 | 64 | 3,331 | 5.6 | 6,462 | 0.52 | 36,381 | 25,742 | 16,860 | 292,354 | 387,320 | 696,534 | 673,022 | 75,088 | 609,826 | 2,216,476 | 3.18x | 2.87x |
| ca-CondMat | 23,133 | 26 | 64 | 1,917 | 12.1 | 3,357 | 0.57 | 8,693 | 3,658 | 12,020 | 73,014 | 123,460 | 208,494 | 196,162 | 92,532 | 200,374 | 1,299,206 | 6.23x | 4.32x |
| cit-HepPh | 34,546 | 31 | 64 | 16,623 | 2.1 | 9,369 | 1.77 | 104,308 | 61,217 | 72,980 | 1,345,614 | 752,308 | 2,170,902 | 2,134,722 | 138,184 | 2,060,682 | 3,199,324 | 1.47x | 1.39x |
| loc-Brightkite | 58,228 | 53 | 64 | 6,326 | 9.2 | 11,547 | 0.55 | 35,833 | 21,427 | 36,232 | 445,960 | 495,148 | 977,340 | 933,328 | 232,912 | 902,652 | 2,628,872 | 2.69x | 2.17x |
| soc-Epinions1 | 75,879 | 68 | 64 | 8,821 | 8.6 | 12,886 | 0.68 | 55,344 | 38,876 | 49,524 | 706,210 | 672,352 | 1,428,086 | 1,379,318 | 303,516 | 1,286,374 | 3,221,518 | 2.26x | 1.86x |
| soc-Slashdot0902 | 82,168 | 56 | 64 | 6,798 | 12.1 | 6,309 | 1.08 | 30,496 | 17,551 | 42,608 | 404,532 | 317,500 | 764,640 | 741,700 | 328,672 | 720,540 | 2,995,392 | 3.92x | 2.74x |
| com-youtube | 1,134,890 | 52 | 64 | 42,815 | 26.5 | 22,403 | 1.91 | 179,601 | 93,541 | 384,064 | 2,261,714 | 1,376,020 | 4,021,798 | 3,934,322 | 4,539,560 | 3,901,854 | 33,847,412 | 8.42x | 3.95x |
| soc-pokec | 1,632,803 | 48 | 64 | 383,206 | 4.3 | 67,544 | 5.67 | 2,375,352 | 1,212,885 | 1,838,988 | 30,291,012 | 11,594,688 | 43,724,688 | 43,456,488 | 6,531,212 | 43,252,464 | 101,770,646 | 2.33x | 2.03x |

"total B" is in memory and includes the derived jump pointers (4 bytes
per node); "file B" is the disk image. "perm B" is the 4 n-byte
permutation from input labels to aligned labels, needed only if the graph
is not stored in the aligned order. "build form total B" is the same index
before `compact_runs` (chain-id DFS arrays).

### 9.2 Build, save, load
| Graph | solve (all-size peel) | trees | chains + labels | layout | build total | compact | save | load | process wall s | peak RSS MB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 9 | 5 | 0 | 0 | 14 | 0.05 | 4.3 | 0.1 | 0.27 | 11 |
| ca-HepPh | 362 | 124 | 6 | 1 | 494 | 0.12 | 1.8 | 0.9 | 1.71 | 178 |
| com-dblp | 909 | 433 | 71 | 9 | 1,422 | 0.46 | 4.2 | 1.0 | 17.01 | 884 |
| web-Stanford | 3,625 | 1,316 | 47 | 18 | 5,007 | 1.11 | 3.5 | 1.1 | 29.12 | 872 |
| amazon0302 | 466 | 241 | 52 | 15 | 777 | 0.77 | 3.5 | 1.5 | 23.50 | 201 |
| ca-AstroPh | 443 | 132 | 2 | 2 | 580 | 0.16 | 2.9 | 0.3 | 8.23 | 53 |
| ca-CondMat | 49 | 23 | 4 | 1 | 77 | 0.06 | 0.4 | 0.2 | 1.65 | 25 |
| cit-HepPh | 1,424 | 246 | 6 | 6 | 1,682 | 0.42 | 6.7 | 0.8 | 35.09 | 223 |
| loc-Brightkite | 1,071 | 580 | 6 | 4 | 1,661 | 0.14 | 1.4 | 0.4 | 6.93 | 237 |
| soc-Epinions1 | 6,215 | 1,734 | 9 | 5 | 7,964 | 0.32 | 13.4 | 0.7 | 23.15 | 1,122 |
| soc-Slashdot0902 | 2,812 | 1,069 | 8 | 3 | 3,892 | 0.19 | 3.5 | 0.4 | 12.07 | 629 |
| com-youtube | 6,395 | 1,242 | 86 | 13 | 7,743 | 0.59 | 6.0 | 1.6 | 63.24 | 1,663 |
| soc-pokec | 55,925 | 10,574 | 372 | 243 | 67,162 | 8.18 | 64.9 | 14.2 | 754.30 | 5,125 |

Process wall time and peak RSS are for the whole `--bench` run (graph
load, build, both forms, all query passes); the peak is the build's dense
core matrix (s_max x n words), the per-size own-node arrays and the chain
grouping, not the index.

### 9.3 Community queries (ns per query; compact form loaded from disk)
| Graph | regime | output vertices | ranges | locate (pointer) | ranges copied | explicit ids | per-vertex S trees memcpy | build form ranges | build form explicit |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | own | 1,658 | 14.1 | 12.7 | 52 | 163 | 225 | 146 | 223 |
| ca-GrQc | half | 2,108 | 15.8 | 13.1 | 24 | 261 | 305 | 264 | 370 |
| ca-GrQc | root | 2,208 | 18.6 | 5.2 | 11 | 277 | 206 | 272 | 236 |
| ca-HepPh | own | 5,327 | 129.8 | 9.1 | 72 | 749 | 625 | 385 | 1,630 |
| ca-HepPh | half | 6,257 | 137.0 | 30.3 | 100 | 1,058 | 832 | 977 | 1,290 |
| ca-HepPh | root | 7,439 | 169.4 | 11.0 | 104 | 1,363 | 794 | 966 | 2,710 |
| com-dblp | own | 164,832 | 871.7 | 21.1 | 404 | 15,601 | 23,493 | 9,319 | 24,661 |
| com-dblp | half | 199,796 | 951.0 | 37.6 | 376 | 17,375 | 26,144 | 10,779 | 28,258 |
| com-dblp | root | 229,952 | 1236.8 | 7.5 | 542 | 18,788 | 30,259 | 10,006 | 28,524 |
| web-Stanford | own | 83,227 | 1217.7 | 9.0 | 486 | 9,936 | 11,983 | 8,515 | 19,024 |
| web-Stanford | half | 102,997 | 1339.6 | 31.9 | 617 | 13,097 | 12,157 | 9,807 | 21,188 |
| web-Stanford | root | 115,934 | 1772.8 | 15.6 | 1,161 | 13,381 | 19,274 | 13,216 | 29,418 |
| amazon0302 | own | 105,997 | 734.6 | 15.1 | 286 | 9,014 | 14,277 | 19,444 | 29,314 |
| amazon0302 | half | 151,391 | 1283.2 | 37.2 | 481 | 15,524 | 24,757 | 27,099 | 40,592 |
| amazon0302 | root | 160,229 | 1529.7 | 10.2 | 601 | 20,929 | 23,022 | 34,318 | 35,900 |
| ca-AstroPh | own | 8,854 | 960.8 | 4.8 | 352 | 3,333 | - | 2,367 | 5,164 |
| ca-AstroPh | half | 9,701 | 984.4 | 18.5 | 408 | 2,328 | - | 2,652 | 5,294 |
| ca-AstroPh | root | 11,775 | 1140.8 | 4.6 | 238 | 2,530 | - | 3,369 | 6,856 |
| ca-CondMat | own | 9,999 | 165.0 | 3.7 | 48 | 1,368 | - | 733 | 2,227 |
| ca-CondMat | half | 12,008 | 181.1 | 15.7 | 89 | 1,515 | - | 884 | 2,800 |
| ca-CondMat | root | 14,389 | 214.9 | 3.2 | 94 | 1,546 | - | 1,717 | 3,251 |
| cit-HepPh | own | 18,252 | 3693.7 | 9.3 | 1,436 | 10,402 | - | 15,709 | 24,742 |
| cit-HepPh | half | 21,348 | 4107.8 | 50.5 | 2,038 | 13,025 | - | 17,842 | 29,077 |
| cit-HepPh | root | 26,461 | 4941.0 | 9.9 | 1,227 | 16,355 | - | 20,959 | 32,657 |
| loc-Brightkite | own | 32,323 | 400.7 | 4.9 | 101 | 3,142 | - | 4,600 | 7,250 |
| loc-Brightkite | half | 39,732 | 415.3 | 10.7 | 156 | 3,145 | - | 5,124 | 8,030 |
| loc-Brightkite | root | 45,279 | 463.9 | 4.1 | 184 | 2,271 | - | 4,291 | 9,727 |
| soc-Epinions1 | own | 49,098 | 922.5 | 7.3 | 407 | 4,612 | - | 8,363 | 10,892 |
| soc-Epinions1 | half | 57,388 | 973.7 | 5.4 | 389 | 5,909 | - | 8,202 | 13,887 |
| soc-Epinions1 | root | 61,009 | 1149.8 | 3.4 | 556 | 3,781 | - | 8,926 | 15,753 |
| soc-Slashdot0902 | own | 47,604 | 469.5 | 14.7 | 192 | 3,965 | - | 6,352 | 10,596 |
| soc-Slashdot0902 | half | 59,057 | 478.6 | 22.3 | 210 | 3,279 | - | 6,351 | 10,768 |
| soc-Slashdot0902 | root | 70,421 | 459.1 | 9.3 | 188 | 4,019 | - | 5,442 | 9,153 |
| com-youtube | own | 774,495 | 959.7 | 11.7 | 343 | 48,238 | - | 39,874 | 86,640 |
| com-youtube | half | 936,661 | 993.2 | 14.5 | 316 | 53,390 | - | 42,417 | 98,510 |
| com-youtube | root | 1,005,937 | 985.7 | 5.2 | 207 | 75,538 | - | 58,831 | 96,997 |
| soc-pokec | own | 827,163 | 40378.2 | 6.0 | 9,831 | 217,377 | - | 504,287 | 826,475 |
| soc-pokec | half | 992,602 | 44650.3 | 30.6 | 15,727 | 230,606 | - | 577,158 | 819,594 |
| soc-pokec | root | 1,240,689 | 51695.1 | 8.6 | 16,244 | 271,898 | - | 634,052 | 853,397 |

"locate" returns the head range, a pointer to the whole runs and the
tail range (no copy). "ranges copied" materialises the range list.
"explicit ids" writes every vertex id. The build-form columns are the
same index before `compact_runs`. The memcpy column is the stage-2
per-vertex S-tree listing, available for the five stage-2 graphs.

### 9.4 Membership, values, ladders (ns per query)
| Graph | member | value | ladder (compact) | ladder steps | ladder (build form) | max depth |
|---|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 9.8 | 13.5 | 58 | 2.08 | 200 | 26 |
| ca-HepPh | 21.1 | 11.3 | 194 | 7.28 | 1,207 | 131 |
| com-dblp | 11.9 | 5.3 | 886 | 5.05 | 12,609 | 157 |
| web-Stanford | 14.7 | 15.8 | 18,303 | 61.20 | 79,821 | 2138 |
| amazon0302 | 12.0 | 18.6 | 639 | 3.22 | 35,283 | 10 |
| ca-AstroPh | 12.6 | 5.0 | 8,430 | 42.85 | 31,152 | 571 |
| ca-CondMat | 5.8 | 9.5 | 319 | 7.10 | 2,046 | 105 |
| cit-HepPh | 30.7 | 6.2 | 23,472 | 48.10 | 129,372 | 1062 |
| loc-Brightkite | 3.5 | 6.9 | 501 | 7.23 | 6,867 | 730 |
| soc-Epinions1 | 7.5 | 6.8 | 8,790 | 17.89 | 36,572 | 1614 |
| soc-Slashdot0902 | 20.0 | 11.4 | 631 | 6.46 | 11,980 | 587 |
| com-youtube | 9.5 | 11.8 | 2,805 | 4.12 | 43,353 | 1604 |
| soc-pokec | 12.4 | 10.3 | 132,043 | 23.04 | 1,994,663 | 2413 |

### 9.5 Explicit listing cost per output vertex (ns)
| Graph | final module (own) | per-vertex S trees (own) | final module (root) | per-vertex S trees (root) |
|---|---:|---:|---:|---:|
| ca-GrQc | 0.098 | 0.136 | 0.126 | 0.090 |
| ca-HepPh | 0.141 | 0.117 | 0.183 | 0.112 |
| com-dblp | 0.095 | 0.143 | 0.082 | 0.132 |
| web-Stanford | 0.119 | 0.145 | 0.115 | 0.153 |
| amazon0302 | 0.085 | 0.133 | 0.131 | 0.139 |
| ca-AstroPh | 0.376 | - | 0.215 | - |
| ca-CondMat | 0.137 | - | 0.107 | - |
| cit-HepPh | 0.570 | - | 0.618 | - |
| loc-Brightkite | 0.097 | - | 0.050 | - |
| soc-Epinions1 | 0.094 | - | 0.062 | - |
| soc-Slashdot0902 | 0.083 | - | 0.057 | - |
| com-youtube | 0.062 | - | 0.075 | - |
| soc-pokec | 0.263 | - | 0.219 | - |

selftests: {"build": {"passed": true, "graphs": 34075, "community_queries": 1816701, "membership_checks": 10955679, "value_checks": 1921263, "ladder_checks": 1039647}, "build-asan": {"passed": true, "graphs": 34075, "community_queries": 1816701, "membership_checks": 10955679, "value_checks": 1921263, "ladder_checks": 1039647}}

### 9.6 The four stage-2 layouts for reference (five graphs)

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
- Bytes: 1.47x (cit-HepPh) to 8.42x (com-youtube) below per-vertex S
  trees, median 3.9x over the thirteen graphs; 1.39x to 4.98x if the
  label permutation is charged to the index. The ratio follows the
  vertex collapse n / chains: 2.1 vertices per chain on cit-HepPh, 4.3 on
  pokec, 5.6 on AstroPh (ratios 1.5x, 2.3x, 3.2x) against 23.6 on dblp and
  26.5 on youtube (7.8x, 8.4x). Where the bytes sit also changes with the
  graph: the per-size layers are 55-70 percent of the index on the
  collaboration, web and product graphs (dblp 68, GrQc 70), the per-chain
  block (trajectories and residue values) 49-69 percent on the social and
  citation graphs (pokec 69, cit-HepPh 62, youtube 56, Slashdot 53,
  Epinions 49); the map is 2-10 percent everywhere.
- Chains against canonical nodes: 0.32-0.68 on nine graphs, but 1.08
  (Slashdot), 1.77 (cit-HepPh), 1.91 (youtube) and 5.67 (pokec). The
  recombination of CHAINS.md Section 7 is real on social and citation
  graphs; the index stays smaller than per-vertex S trees there because
  the per-vertex trees pay per vertex and per size, and chains still
  collapse 2-27 vertices each.
- Runs against (chain, size) pairs: 1.4x-3.3x fewer entries. The compact
  form costs 1.1-14.2 percent more bytes than the chain-id form (the
  8-byte entry per node) and pays back in queries: range answers 2.8x-284x
  faster, explicit answers 0.85x-4.3x (one of 39 points slower, GrQc root),
  ladders 3.5x-55x.
- Locating a community takes 3-51 ns on every graph and regime: bitmap
  rank, trajectory lookup, jump-pointer climb, two entry reads. Copying the
  range list takes 11 ns to 16 us (pokec root, 51,695 ranges); the range
  list is 5x (cit-HepPh) to 807x (youtube) shorter than the vertex list.
- Explicit listing runs at 0.06-0.14 ns per vertex on ten graphs, which
  is at or below the memcpy baseline (0.09-0.15 ns per vertex on the five
  graphs where it was measured: own level 0.84x-1.58x of memcpy, root
  0.58x-1.61x, faster on dblp, Stanford, amazon, slower on HepPh), and
  0.26-0.57 ns per vertex on the three fragmented graphs (pokec 20
  vertices per range, AstroPh 9, cit-HepPh 5), where the per-range
  constant (about 3-4 ns) dominates.
- Membership 3.5-31 ns and values 5-19 ns: two bitmap ranks plus a climb,
  or one rank plus a residue or binomial lookup.
- Build: 14 ms (GrQc) to 67 s (pokec, 1.63 M vertices), 60-85 percent in
  the all-size peel (the terminal solver), the rest in the per-size
  union-find trees; chains and layout add at most 0.6 s. Peak build memory
  reaches 5.1 GB on pokec against a 44 MB index.

Hypotheses (not measured): a branchless or prefetching fill would cut the
per-range constant on the fragmented graphs; a per-node vertex count
(4 bytes per node) would make ladders O(depth) (pokec 132 us, cit-HepPh
23 us today); the permutation compresses to about n log2(C) bits (dblp
0.55 MB instead of 1.27 MB) with a wavelet tree over the chain ids in
input order; the build's dense core matrix and per-size own arrays can be
streamed size by size.

Algorithmic versus engineering gains: the chain partition and the
lexicographic rank order (one range at s = 2, few ranges elsewhere) are
algorithmic; the run array with entry points is a layout theorem
(CHAINS.md C6) implemented as a constant number of array reads; the
vectorised fill and the Myers skip pointers are engineering.

## 11. Failure Cases / Limitations

- The label permutation is not free: 4 n bytes unless the graph is stored
  in aligned order (60 percent of the dblp index, 113 percent on youtube
  where the index is tiny against n). Reported both ways.
- Explicit listing is slower than a memcpy when chains are small and
  communities fragmented (cit-HepPh, AstroPh, pokec); the range form is
  always far faster, but it is a different output.
- Ladder queries on deep trees (pokec depth 2,413, cit-HepPh 1,062) cost
  23-132 us because each level sums its runs.
- Peak build memory is the solver's dense s_max x n core matrix plus the
  per-size own arrays (5.1 GB on pokec), not the index.
- The number of chains exceeds the number of canonical nodes on four of
  thirteen graphs (up to 5.7x on pokec); no bound in terms of the
  hierarchy exists (CHAINS.md Section 7).
- The memcpy baseline was measured only on the five stage-2 graphs, in
  the stage-2 binary with its own query draw (same seed, protocol and
  machine); output sizes differ by 0.1-0.3 percent between the draws.
- One machine, thirteen graphs up to 1.6 M vertices, in-memory
  single-process timing; no cold-cache or multi-process protocol; no
  comparison with a reimplemented SGL (its inputs are bipartite).

## 12. Final Conclusion

The all-size r = 1 nucleus hierarchy of a graph is stored exactly in
1.5x-8.4x fewer bytes than one S tree per size (median 3.9x over thirteen
graphs, 3.7x-7.8x on the five stage-2 graphs), communities are located in
constant time (3-51 ns) and listed at memcpy speed or better wherever a
chain holds a few dozen vertices, from a partition of the vertices
(chains) that the theory proves exact, aligned labels that make chains id
ranges, and per-size run arrays with node entry points that make every
community a head range, whole runs and a tail range. Correctness is
brute-force verified in three forms on 34,075 graphs under Release and
sanitizers, and the canonical node counts agree with the independent
stage-1/2 pipeline on the five shared graphs. The SGL skyline dedup is
dominated in this setting and dropped. Framing for a paper: the
contribution starts from the chain partition (hierarchy equivalence),
not from storing S trees; S trees over vertices are the baseline it
beats on both axes; for r >= 2 the paper's size forest already plays the
role of the chain structure (CHAINS.md Section 10).

## 13. Next Recommended Improvements

1. Per-node vertex counts (4 bytes per node) for O(depth) ladders and
   O(1) community sizes.
2. Branchless or prefetching range fill for fragmented graphs
   (cit-HepPh, AstroPh, pokec); measure the per-range constant.
3. Compress the permutation to n log2(C) bits (wavelet tree over chain
   ids), or store the graph in aligned order and drop it.
4. Stream the build size by size to cut the 5 GB peak on pokec-sized
   inputs; then com-lj and com-orkut on the server.
5. A bound on the number of chains for graphs with monotone core
   trajectories, or a construction showing none exists beyond
   CHAINS.md Section 7; the four graphs with chains > nodes are the
   test cases.
6. r >= 2 is closed (CHAINS.md Section 10 verdict): no separate chain
   index; the unifying remark belongs in the paper's discussion.
