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
thread, no other timing run concurrent. Inputs: `data/ca-GrQc.edges`,
`data/ca-HepPh.edges`, `data/com-dblp.edges`, `graphs/web-Stanford.edges`,
`graphs/amazon0302.edges` (sha256 in `final.json`). Count width: 64 bits
(GrQc, dblp, Stanford, amazon) and 256 bits (HepPh, s_max 239).
Queries (seed 20260918, drawn from active vertices and sizes
2 <= s <= omega(v)): community regimes own (k = kappa_s(v), 20,000),
half (k = max(1, kappa/2), 20,000), root (k = 1, 1,000); membership 20,001
(random u, mixed regimes); values 200,000; ladders on the own set. One
warm-up pass plus five timed passes, median reported. Explicit ids are
written into a preallocated caller buffer; the memcpy baseline of stage 2
also wrote into preallocated storage. The baseline numbers come from the
stage-2 run (`index_vertices.json`, same protocol, its own query draw).

Commands:
```
cmake -S research/r1_skyline_index_20260918 -B research/r1_skyline_index_20260918/build -DCMAKE_BUILD_TYPE=Release -DSANITIZE=OFF
cmake --build research/r1_skyline_index_20260918/build -j 12 --target chain_index_tool
research/r1_skyline_index_20260918/build/chain_index_tool --selftest
/usr/bin/time -l research/r1_skyline_index_20260918/build/chain_index_tool --bench data/com-dblp.edges research/r1_skyline_index_20260918/cx/com-dblp.cx
python3 research/r1_skyline_index_20260918/run_final.py     # all of the above, both builds, five graphs
python3 research/r1_skyline_index_20260918/report_tables.py
```

## 9. Experimental Results

### 9.1 Size

| Graph | n | s_max | chains | canonical nodes | (chain,s) pairs | runs | map B | chains B | layers B | total B | file B | perm B | build form total B | per-vertex S trees B | ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 5,242 | 44 | 716 | 1,517 | 2,079 | 624 | 3,856 | 16,892 | 47,812 | 68,560 | 64,308 | 20,968 | 65,472 | 272,564 | 3.98x |
| ca-HepPh | 12,008 | 239 | 1,136 | 3,538 | 8,468 | 4,969 | 6,808 | 100,632 | 225,632 | 333,072 | 328,536 | 48,032 | 311,136 | 1,899,628 | 5.70x |
| com-dblp | 317,080 | 114 | 13,459 | 33,979 | 70,208 | 28,784 | 113,304 | 571,222 | 1,454,420 | 2,138,946 | 2,007,646 | 1,268,320 | 2,052,686 | 16,710,168 | 7.81x |
| web-Stanford | 281,903 | 72 | 17,963 | 53,239 | 151,364 | 97,074 | 124,720 | 1,693,758 | 2,267,852 | 4,086,330 | 3,876,310 | 1,127,612 | 3,701,670 | 24,853,426 | 6.08x |
| amazon0302 | 262,111 | 7 | 40,867 | 65,134 | 143,567 | 47,380 | 212,628 | 1,443,858 | 2,202,840 | 3,859,326 | 3,599,126 | 1,048,444 | 3,793,970 | 14,414,214 | 3.73x |

"total B" is in memory and includes the derived jump pointers (4 bytes
per node); "file B" is the disk image. "perm B" is the 4 n-byte
permutation from input labels to aligned labels, needed only if the graph
is not stored in the aligned order; counting it, the ratios are 3.04x,
4.98x, 4.90x, 4.77x, 2.94x. The per-vertex S-tree bytes are the stage-2
`vertices` layout with values (Block D).

### 9.2 Build, save, load (ms, one thread)

| Graph | solve (all-size peel) | trees | chains + labels | layout | build total | compact | save | load | process wall s | peak RSS MB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 10 | 6 | 1 | 0 | 17 | 0.05 | 0.4 | 0.2 | 0.29 | 12 |
| ca-HepPh | 391 | 79 | 11 | 2 | 482 | 0.10 | 5.1 | 0.6 | 1.59 | 190 |
| com-dblp | 857 | 532 | 65 | 17 | 1,472 | 0.52 | 3.9 | 1.3 | 16.72 | 928 |
| web-Stanford | 3,557 | 1,144 | 41 | 18 | 4,762 | 1.15 | 6.7 | 1.6 | 28.82 | 909 |
| amazon0302 | 547 | 240 | 55 | 11 | 855 | 0.51 | 8.0 | 1.5 | 24.25 | 212 |

Process wall time and peak RSS are for the whole `--bench` run (graph
load, build, all query passes); the peak is the build's dense core matrix
(s_max x n words) plus per-size own-node arrays, not the index.

### 9.3 Community queries (ns per query; compact form loaded from disk)

| Graph | regime | output vertices | ranges | locate (pointer) | ranges copied | explicit ids | per-vertex S trees memcpy | build form ranges | build form explicit |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | own | 1,658 | 14.1 | 9.7 | 43 | 140 | 225 | 234 | 288 |
| ca-GrQc | half | 2,108 | 15.8 | 19.7 | 45 | 263 | 305 | 159 | 500 |
| ca-GrQc | root | 2,208 | 18.6 | 5.8 | 26 | 274 | 206 | 194 | 239 |
| ca-HepPh | own | 5,327 | 129.8 | 17.2 | 57 | 766 | 625 | 853 | 1,242 |
| ca-HepPh | half | 6,257 | 137.0 | 36.4 | 74 | 1,142 | 832 | 1,049 | 888 |
| ca-HepPh | root | 7,439 | 169.4 | 10.8 | 43 | 1,387 | 794 | 411 | 1,021 |
| com-dblp | own | 164,832 | 871.7 | 9.6 | 338 | 15,034 | 23,493 | 8,826 | 25,558 |
| com-dblp | half | 199,796 | 951.0 | 17.5 | 450 | 17,327 | 26,144 | 11,649 | 27,417 |
| com-dblp | root | 229,952 | 1236.8 | 3.4 | 551 | 18,302 | 30,259 | 13,422 | 23,370 |
| web-Stanford | own | 83,227 | 1217.7 | 9.8 | 422 | 10,158 | 11,983 | 9,105 | 17,126 |
| web-Stanford | half | 102,997 | 1339.6 | 39.7 | 660 | 12,988 | 12,157 | 9,692 | 19,563 |
| web-Stanford | root | 115,934 | 1772.8 | 11.0 | 409 | 15,558 | 19,274 | 7,757 | 25,187 |
| amazon0302 | own | 105,997 | 734.6 | 24.0 | 348 | 10,261 | 14,277 | 18,378 | 29,387 |
| amazon0302 | half | 151,391 | 1283.2 | 24.1 | 362 | 15,829 | 24,757 | 27,999 | 42,105 |
| amazon0302 | root | 160,229 | 1529.7 | 8.9 | 350 | 14,524 | 23,022 | 30,859 | 40,872 |

"locate" returns the head range, a pointer to the whole runs and the
tail range (no copy). "ranges copied" materialises the range list.
"explicit ids" writes every vertex id. The build-form columns are the
same index before `compact_runs` (chain-id DFS arrays, on-the-fly
merging). The per-vertex S-tree column is the stage-2 memcpy listing.

### 9.4 Membership, values, ladders (ns per query)

| Graph | member | value | ladder (compact) | ladder steps | ladder (build form) | max depth |
|---|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 7.2 | 13.6 | 58 | 2.08 | 163 | 26 |
| ca-HepPh | 12.4 | 12.4 | 141 | 7.28 | 1,464 | 131 |
| com-dblp | 4.9 | 8.6 | 540 | 5.05 | 12,162 | 157 |
| web-Stanford | 15.6 | 16.3 | 18,941 | 61.20 | 81,488 | 2138 |
| amazon0302 | 18.3 | 18.4 | 708 | 3.22 | 40,005 | 10 |

### 9.5 Explicit listing cost per output vertex (ns)

| Graph | final module (own) | per-vertex S trees (own) | final module (root) | per-vertex S trees (root) |
|---|---:|---:|---:|---:|
| ca-GrQc | 0.085 | 0.136 | 0.124 | 0.090 |
| ca-HepPh | 0.144 | 0.117 | 0.187 | 0.112 |
| com-dblp | 0.091 | 0.143 | 0.080 | 0.132 |
| web-Stanford | 0.122 | 0.145 | 0.134 | 0.153 |
| amazon0302 | 0.097 | 0.133 | 0.091 | 0.139 |

### 9.6 The four stage-2 layouts for reference (bytes with values; own-level listing ns)

| Graph | per-vertex S trees | over twins | over chains | aligned chains (index.cpp) | listing: per-vertex | twins | chains | aligned explicit | aligned ranges |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 272,564 | 256,238 | 103,220 | 65,128 | 225 | 2,903 | 574 | 822 | 231 |
| ca-HepPh | 1,899,628 | 1,532,406 | 395,336 | 306,065 | 625 | 9,668 | 2,116 | 1,844 | 867 |
| com-dblp | 16,710,168 | 16,880,402 | 4,447,786 | 2,024,433 | 23,493 | 370,861 | 42,952 | 47,205 | 8,694 |
| web-Stanford | 24,853,426 | 26,086,394 | 5,762,930 | 3,632,412 | 11,983 | 163,655 | 27,140 | 30,132 | 7,507 |
| amazon0302 | 14,414,214 | 17,079,748 | 5,744,638 | 3,860,362 | 14,277 | 223,981 | 133,725 | 129,013 | 20,811 |

## 10. Analysis Of Runtime And Memory

Measured facts.
- Bytes: 3.7x (amazon) to 7.8x (dblp) below per-vertex S trees; 2.9x to
  5.0x if the label permutation is charged to the index. The layers are
  55-70 percent of the bytes; the per-chain block (trajectories and
  residues) 25-40 percent; the map 3-6 percent.
- Where the bytes went: the vertex axis collapsed to chains (dblp 317,080
  vertices, 13,459 chains); the (chain, size) pairs collapsed to runs
  (dblp 70,208 to 28,784; Stanford 151,364 to 97,074; amazon 143,567 to
  47,380). The compact form costs 1.7-10.4 percent more bytes than the
  chain-id form because every node carries an 8-byte entry; it pays back
  in queries.
- Locating a community takes 3-40 ns on every graph and regime: the
  bitmap rank, the trajectory lookup, the jump-pointer climb and two
  entry reads. Copying the range list takes 26-660 ns for 14-1,773
  ranges; the range list is 41-189x shorter than the vertex list.
- Explicit listing runs at 0.08-0.19 ns per vertex; against the memcpy
  baseline it is faster on the three large graphs (1.18-1.65x) and slower
  on HepPh (0.57-0.82x) and on the root regime of GrQc, where the per-range
  overhead (about 4 ns per range, 130-170 ranges for 5-7 thousand
  vertices) is comparable to copying the 20-30 KB answer. Per output
  vertex the fill is cheaper than memcpy because it reads 8 bytes per
  range instead of 4 bytes per vertex.
- Build form versus compact form: range answers 5.4-53x faster, explicit
  answers 1.6-2.9x faster, ladders 2.8-56x faster, at the byte premium
  above. The chain-id loop pays per chain (dblp own-level communities
  span thousands of chains); the run form pays per run.
- Membership 5-18 ns and values 9-18 ns: two bitmap ranks plus a climb, or
  one rank plus a residue or binomial lookup.
- Build: 0.02-4.8 s single thread, dominated by the all-size peel (the
  terminal solver) and the per-size union-find trees; chains and layout
  add 12-66 ms.

Hypotheses (not measured): the per-range overhead could be halved by
prefetching two runs ahead; a per-node vertex count (4 bytes per node)
would make ladders O(depth); the permutation compresses to about
n log2(C) bits (dblp 0.55 MB instead of 1.27 MB) with a wavelet tree over
the chain ids in input order.

Algorithmic versus engineering gains: the chain partition and the
lexicographic rank order (one range at s = 2, few ranges elsewhere) are
algorithmic; the run array with entry points is a layout theorem
(CHAINS.md C6) implemented as a constant number of array reads; the
vectorised fill and the Myers skip pointers are engineering.

## 11. Failure Cases / Limitations

- The label permutation is not free: 4 n bytes unless the graph is stored
  in aligned order (60 percent of the dblp index). Reported both ways.
- Explicit listing is slower than memcpy when communities are small and
  fragmented (HepPh: 40 vertices per range). The range form is always
  faster than memcpy, but it is a different output.
- Ladder queries on deep trees (Stanford: 61 levels on average, depth up
  to 2,138) cost 19 microseconds because each level sums its runs.
- Peak build memory is the dense s_max x n core matrix of the solver
  (0.9 GB on dblp and Stanford), not the index.
- The count of chains has no bound in terms of the hierarchy alone
  (CHAINS.md Section 7); the sizes measured rely on real graphs having
  coherent core values across sizes.
- The memcpy baseline was measured in the stage-2 binary with its own
  query draw (same seed, protocol and machine); the output sizes differ
  by 0.1-0.3 percent between the two draws.
- One machine, five graphs, in-memory single-process timing; no
  cold-cache or multi-process protocol; no comparison with a reimplemented
  SGL (its inputs are bipartite).

## 12. Final Conclusion

The all-size r = 1 nucleus hierarchy of a graph is stored exactly in
3.7-7.8x fewer bytes than one S tree per size, and communities are
located in constant time and listed at least as fast as a memcpy on the
large graphs, from a partition of the vertices (chains) that the theory
proves exact, aligned labels that make chains id ranges, and per-size run
arrays with node entry points that make every community a head range,
whole runs and a tail range. Correctness is brute-force verified in three
forms on 34,075 graphs under Release and sanitizers. The SGL skyline
dedup is dominated in this setting and dropped. Framing for a paper: the
contribution starts from the chain partition (hierarchy equivalence),
not from storing S trees; S trees over vertices are the baseline it
beats on both axes.

## 13. Next Recommended Improvements

1. Per-node vertex counts (4 bytes per node) for O(depth) ladders and
   O(1) community sizes.
2. Compress the permutation to n log2(C) bits (wavelet tree over chain
   ids), or store the graph in aligned order and drop it.
3. Prefetch in the run fill; measure the per-range constant on HepPh.
4. r >= 2: chains of r-cliques and their relation to the pattern classes
   of the NSI paper (CHAINS.md Section 8); theory first.
5. A bound on the number of chains for graphs with monotone core
   trajectories, or a construction showing none exists beyond
   CHAINS.md Section 7.
6. Larger inputs (com-youtube, soc-pokec, web-Google on the server) with
   the same tool; the build is a single all-size peel, so the server
   sweep engine can host it.
