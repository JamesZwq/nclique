# Skyline Index Counts: Results

Date: 2026-09-18. Counting experiment only. No index was implemented, no
query was timed, no production, paper or default changed. The words
formulas below are a storage model, not a measured file size.

## 1. Problem Summary

For r = 1 and every size s = 2..d+1, decide whether the skyline community
index of [THEORY.md](THEORY.md) would be smaller than S explicit merge
trees, and verify the lemmas the index relies on, on the five local
full-range inputs. The quantities are those of THEORY.md Section 11.

## 2. Current Baseline And Bottlenecks

Baseline representation: one canonical merge tree per size plus one
position per active (class, size) pair ("S trees"). Values come from
`terminal::Solver<T>::solve` on the plain (mode 0) terminal index and are
checked against the frozen `Kernel<T>::fixed_sparse` control on every
graph. Twin classes are closed-neighbourhood classes.

## 3. Candidate Algorithm Ideas

Only the counting program is at stake here; the design candidates are
in THEORY.md Sections 4 and 7 (skyline entries, canonical nodes with one
cross-size pointer, certified-tail collapse, twin quotient).

## 4. Chosen Approach And Rationale

`count.cpp` implements THEORY.md Section 6 Step 1 (per-size union-find
over live rows, canonical nodes created at the level where a component
changes), the Kruskal-Katona shadow with `cpp_int`, the skyline, the
certification point, and the checks F1, F2, F5, F8, L1, L2. The selftest
compares every level's partition and every canonical interval with a
brute-force nucleus computation from explicit cliques.

The first version, written by the Codex subagent, stopped at the shadow
unit test because IMPLEMENTATION.md carried wrong sigma_3 examples
(sigma_3 acts on the 3-cascade of kappa_4; sigma_3(35) = 21, not
sigma_3(21) = 21). The spec was corrected (commit 378e1aa). Review of the
Codex code then fixed: the recorded creator of a node (it took the first
vertex of the level, not a vertex of the component, which invalidates
the L2 check), the partition comparison (it compared label numbers
instead of partitions), the canonical interval comparison (it compared
a top level with an interval and swept levels absent from the data), the
F5 check (first and last class member only), the residue count (missing
the top cell of an uncertified vertex), and the ASan environment
(`detect_leaks` is unsupported on macOS). The colex brute force required
by the spec was added.

## 5. Complexity Discussion

The count run costs one all-size decomposition plus, per size, one pass
over the valid rows with union-find, plus one shadow per (class, size)
pair outside the certified tail. Wall times were 0.02 s (GrQc) to 6.4 s
(Stanford); peak RSS up to 1.58 GB on dblp, dominated by the explicit
128-bit core matrix and the selftest-free trees. These are recorded for
the record, not as a claim.

## 6. Implementation Summary

`count.cpp` (program), `CMakeLists.txt`, `run.py` (serial driver:
Release and ASan/UBSan builds, both selftests, then the five graphs one
at a time under `/usr/bin/time -l`; raw JSON in `counts.json`, logs in
`counts-logs/`), `verify_theory.py` (independent Python brute-force check
of F2, F6, F7, F8, L1, L2, L3 and Theorem Q2, including the full
retrieval procedure, on random graphs up to 11 vertices).

## 7. Correctness Validation Summary

| Check | Release | ASan/UBSan |
|---|---:|---:|
| Selftest graphs (all labelled graphs on <= 6 vertices, 200 random 7..10, split and complete graphs) | 34,075 | 34,075 |
| Core cells compared with the explicit-clique definition | 844,230 | 844,230 |
| Level partitions compared with brute-force nuclei | 93,057 | 93,057 |
| Canonical intervals compared with brute force | every node | every node |

Independent Python verifier: 1,501 random graphs up to 11 vertices,
18,976 F2 checks, 62,100 dominance containments (F6), 17,228 L2 chains,
58,924 retrieval queries (Theorem Q2), no violation.

On the five real graphs (exact wide integers, no clamp):

| Graph | F2 checks | F2 violations | L1 checks | L2 chain checks | control equality |
|---|---:|---:|---:|---:|---|
| ca-GrQc | 14,804 | 0 | 20,045 | 14,126 | equal |
| ca-HepPh | 160,957 | 0 | 172,963 | 156,585 | equal |
| com-dblp | 929,853 | 0 | 1,246,933 | 869,613 | equal |
| web-Stanford | 1,326,086 | 0 | 1,607,989 | 825,077 | equal |
| amazon0302 | 562,880 | 0 | 824,991 | 356,277 | equal |

F1, F5 and F8 held on every vertex and class of every graph.

## 8. Experimental Setup

Local macOS arm64, Clang from /opt/homebrew, `-O3 -DNDEBUG`, one thread,
run from the repository root. Exact commands (also in `counts-logs/*.log`
first lines and in `run.py`):

```
cmake -S research/r1_skyline_index_20260918 -B research/r1_skyline_index_20260918/build -DCMAKE_BUILD_TYPE=Release -DSANITIZE=OFF
cmake --build research/r1_skyline_index_20260918/build -j 12
research/r1_skyline_index_20260918/build/count --selftest
cmake -S research/r1_skyline_index_20260918 -B research/r1_skyline_index_20260918/build-asan -DCMAKE_BUILD_TYPE=Release -DSANITIZE=ON
cmake --build research/r1_skyline_index_20260918/build-asan -j 12
research/r1_skyline_index_20260918/build-asan/count --selftest
/usr/bin/time -l research/r1_skyline_index_20260918/build/count --graph data/ca-GrQc.edges        (and the other four)
python3 research/r1_skyline_index_20260918/verify_theory.py 7 1500
```

Inputs: data/ca-GrQc.edges, data/ca-HepPh.edges, data/com-dblp.edges,
graphs/web-Stanford.edges, graphs/amazon0302.edges. Full size range
s = 2..d+1: 44 / 239 / 114 / 72 / 7. Count widths 64 / 256 / 128 / 64 / 64.

## 9. Experimental Results

Per-class quantities (classes = closed-neighbourhood twins):

| Graph | n | classes | active (class,s) pairs | skyline entries E | avg_traj | mu | avg_traj/mu |
|---|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 5,242 | 4,105 | 12,704 | 4,780 | 3.10 | 1.17 | 2.66 |
| ca-HepPh | 12,008 | 9,559 | 114,694 | 13,909 | 12.00 | 1.46 | 8.25 |
| com-dblp | 317,080 | 256,983 | 928,980 | 316,949 | 3.62 | 1.23 | 2.93 |
| web-Stanford | 281,903 | 270,303 | 1,429,928 | 721,162 | 5.29 | 2.67 | 1.98 |
| amazon0302 | 262,111 | 253,766 | 787,990 | 457,900 | 3.11 | 1.80 | 1.72 |

Canonical trees and certification (per vertex cells):

| Graph | N_T (all sizes) | largest T_s | empty nodes | chain nodes | max depth | residue cells | residue with delta = 0 | certified cells |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 1,517 | 448 | 844 | 353 | 27 | 1,058 (5.3%) | 171 | 18,987 (94.7%) |
| ca-HepPh | 3,538 | 381 | 2,669 | 2,220 | 132 | 6,004 (3.5%) | 895 | 166,959 (96.5%) |
| com-dblp | 33,979 | 6,459 | 23,761 | 7,095 | 158 | 84,010 (6.7%) | 13,524 | 1,162,923 (93.3%) |
| web-Stanford | 53,239 | 5,451 | 14,869 | 44,039 | 2,139 | 872,260 (54.2%) | 238,442 | 735,729 (45.8%) |
| amazon0302 | 65,134 | 27,685 | 31,037 | 21,951 | 11 | 355,255 (43.1%) | 74,957 | 469,736 (56.9%) |

Storage model in 4-byte words (IMPLEMENTATION.md Section 8; S trees =
class positions + 4 words per node; skyline = entries + gamma bytes +
6 words per node including the cross-size pointer and reverse list):

| Graph | words S trees | words skyline | ratio | inequality lhs vs rhs (THEORY.md Sec. 7) | verdict |
|---|---:|---:|---:|---|---|
| ca-GrQc | 18,772 | 15,077 | 1.25x | 7,926 > 4,229 | skyline smaller |
| ca-HepPh | 128,846 | 38,614 | 3.34x | 100,806 > 10,554 | skyline smaller |
| com-dblp | 1,064,896 | 600,060 | 1.77x | 612,031 > 147,195 | skyline smaller |
| web-Stanford | 1,642,884 | 1,220,886 | 1.35x | 708,766 > 286,769 | skyline smaller |
| amazon0302 | 1,048,526 | 963,179 | 1.09x | 330,090 > 244,743 | skyline smaller |

Raw JSON with every field is in `counts.json`.

## 10. Analysis Of Runtime And Memory

Nothing about runtime was measured. On storage the picture is:

- The membership entries shrink by avg_traj/mu: 2.7x to 2.9x on the
  collaboration graphs, 8.3x on HepPh (long trajectories, 96 percent
  certified), 2.0x on Stanford, 1.7x on amazon.
- The canonical-node count N_T is small against n (0.3 to 25 percent of
  the classes) but not negligible against E on amazon (65,134 nodes
  against 457,900 entries), where the model ratio falls to 1.09x. On the
  other four graphs the node term does not decide the outcome.
- Empty canonical nodes are the majority on the four graphs other than
  Stanford (56 to 75 percent): most nodes are containers or connectors
  carrying no skyline entry, which is the certified-region collapse of
  THEORY.md Section 4 in numbers. Removing them would need sibling links
  (THEORY.md Section 8); not done, not modelled.
- Residue values: on the three collaboration graphs 3.5 to 6.7 percent
  of the cells are residue and only 0.5 to 1.1 percent are residue cells
  with delta = 0 (the only cells that Variant B must reconstruct).
  Storing all residue values explicitly (Variant A) keeps point value
  queries at one lookup at a cost of 84,010 wide integers on dblp. On
  Stanford and amazon residue is 43 to 54 percent of the cells, so the
  value block stays comparable to the dense table there.
- Depth: Stanford's merge trees reach depth 2,139, so the walk of Q2
  step 3 and of construction step 3 needs jump pointers there.

The earlier per-vertex estimates from the 2026-05 saturation data are
reproduced (GrQc 3.4x, dblp 3.3x, amazon 1.8x per vertex); the HepPh
estimate from that data (1.9x) was wrong, as flagged, because that
experiment clamped high sizes; the exact value is 10.6x per vertex.

## 11. Failure Cases / Limitations

- Storage is a word model, not a written file; per-node fields could be
  narrower or wider in a real layout.
- amazon0302 is a marginal case (1.09x); with sibling links or any extra
  per-node field the sign could flip.
- Stanford's low certification (46 percent) and depth 2,139 make it the
  adverse input for both the value block and the query walk.
- The S-tree baseline is modelled at one position per active (class,
  size) pair; a per-vertex layout would be avg 1.1 to 1.2 times larger
  for both designs on these graphs (class collapse is weak at r = 1).
- No query was run; the constant-factor loss of skyline retrieval against
  a contiguous interval scan (THEORY.md Section 7) is unmeasured.
- Only mode 0 rows were used for connectivity; factored rows (mode 2)
  need the live condition with choice counters and were not exercised.

## 12. Final Conclusion

All proof gates of THEORY.md hold on the five real graphs in exact
arithmetic (F2: 2,994,580 checks, 0 violations; L1 and L2 on every
active cell), and the implementation matches brute-force nuclei on
34,075 graphs in Release and under sanitizers. Under the storage model,
the skyline index is smaller than S explicit trees on all five inputs,
by 1.09x (amazon) to 3.34x (HepPh); the membership part alone shrinks
1.7x to 8.3x. This passes the build test of THEORY.md Section 7 on every
graph, marginally on amazon. Nothing about query or build speed is
claimed.

## 13. Next Recommended Improvements

1. Decide Variant A against B per graph from the residue share; on the
   collaboration graphs Variant A is cheap and keeps constant-time values.
2. Model and count the sibling-link variant that drops empty nodes (56 to
   75 percent of N_T) before writing any index code.
3. Add level-ancestor jump pointers in the design for deep trees.
4. Only then implement the index and measure bytes and query latency
   against S trees laid out as DFS intervals, the strongest same-content
   baseline.
