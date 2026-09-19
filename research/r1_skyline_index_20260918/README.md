# r = 1 All-Size Community Index (chain index)

Research line, 2026-09-18/19. One index for the whole (1, s)-nucleus
hierarchy of a graph over every clique size s: core values kappa_s(v),
every (s, k)-nucleus as a vertex set, membership, and the ladder of a
vertex. Exact; brute-force verified; measured on 13 graphs.

Headline (17 graphs of five families, one thread, RESULTS_FINAL.md
Sections 9 and 14): 1.3x-11.7x fewer bytes than one S tree per size
(median 4.95x; dblp 10.3x; the 4.05 M-vertex dense dblp-coauthor 11.7x;
the 1,005-vertex email-Eu-core 1.3x is the worst case); a community is
located in 3.5-59 ns and listed at 0.04-0.28 ns per vertex on all but the
most fragmented graphs (faster than a memcpy of the vertex list on 13 of
15 measured points); membership 5-31 ns; values 6-20 ns; build 15 ms
(GrQc) to 34 s (pokec) and 189 s (dblp-coauthor, 512-bit counts) with
157 MB (dblp) to 8.8 GB (dblp-coauthor) peak memory. com-lj and com-orkut
exceed the solver's 32-bit member ids (see RESULTS_FINAL.md Section 11).

## What the index is

1. Chains: vertices with the same own canonical node at every size
   (hierarchy equivalence, proved exact; CHAINS.md C1-C3). Vertices are
   relabelled so that every chain is one id range; the vertex-to-chain map
   is a bitmap with a rank directory.
2. Per size s: the canonical merge tree over chains (tops, parents,
   subtree sizes) and its DFS array stored as maximal runs of consecutive
   labels with one entry point per node, so a community is a head range,
   whole runs and a tail range found in O(1) after the climb (C6); every
   (2, k)-community is exactly one range (C4).
3. Values per chain: omega, the certification point sigma, and residues
   for s < sigma; for s >= sigma the value is C(omega - 1, s - 1).
4. Node tops and residues at per-size byte widths; file format CHAINX04
   (omega and sigma 16-bit, so clique sizes above 255 are fine).

## Layout

| Path | What |
|---|---|
| `chain_index.hpp` | the index: header-only `chainindex::ChainIndex<T>` (build form, compact form, queries, save/load) |
| `chain_index_tool.cpp` | `--selftest`, `--build-only <graph>`, `--bench <graph> <out.cx>` |
| `count.cpp`, `common.hpp` | shared pieces: `make_tree_row` (per-size canonical trees), brute-force helpers, input preparation; also the stage-1 counting program |
| `CMakeLists.txt` | targets `chain_index_tool`, `count`, `index`, `chains`; `-DSANITIZE=ON` for ASan/UBSan |
| `run_final.py [extra graphs]` | Release + sanitizer builds, both selftests, `--bench` on the five stage-2 graphs plus any extra ones; writes `final.json`, `final-logs/` |
| `run_buildonly.py` | `--build-only` on the 13 graphs; writes `buildonly.json`, `buildonly-logs/` |
| `report_tables.py` | prints the tables of RESULTS_FINAL.md from `final.json` and the stage-2 files |
| `THEORY.md` | the hierarchy theory: dominance through Kruskal-Katona shadows, canonical nodes, certified tail, retrieval; Section 14 says what was built |
| `CHAINS.md` | chains: Lemmas C1-C8, chain counts, the final layout (Section 9), the r >= 2 verdict (Section 10) |
| `RESULTS_FINAL.md` | the 13-section report of the final module; Section 14 = build memory and per-size widths |
| `final.json`, `final-logs/` | evidence of the current module (13 graphs, latencies, bytes, build) |
| `more.json`, `more-logs/` | the same for four more graphs (ca-HepTh, email-Eu-core, com-amazon, dblp-coauthor) |
| `buildonly.json`, `buildonly-logs/` | build phases and resident memory of the current module |
| `stages/` | stage 1 (counting gate) and stage 2 (four layouts: per-vertex S trees, twins, chains, aligned; skyline dedup): `IMPLEMENTATION*.md`, `RESULTS.md`, `RESULTS_INDEX.md`, `RESULTS_CHAINS.md`, `index.cpp`, `chains.cpp`, `verify_theory.py`, `run.py`, `run_index.py`, `counts.json`, `index*.json` and their logs |
| `archive/` | superseded evidence of earlier module versions: `final_v1` (five graphs), `final_v2` (13, scalar-tail fill), `final_v3` (13, branchless fill, fixed-width values), `index_aligned_v1` |
| `build/`, `build-asan/`, `cx/` | build directories and index files, not tracked |

The all-size peel is the solver of `../r1_terminal_20260918/terminal.hpp`
(`solve(..., sink)` streams one core row at a time; the default path is
unchanged); `count.cpp` includes `../r1_orderreplay_20260917/shared_harness.inc`
for input preparation.

## Build and run

```
cmake -S research/r1_skyline_index_20260918 -B research/r1_skyline_index_20260918/build -DCMAKE_BUILD_TYPE=Release -DSANITIZE=OFF
cmake --build research/r1_skyline_index_20260918/build -j 12
research/r1_skyline_index_20260918/build/chain_index_tool --selftest
research/r1_skyline_index_20260918/build/chain_index_tool --build-only data/com-dblp.edges
research/r1_skyline_index_20260918/build/chain_index_tool --bench data/com-dblp.edges /tmp/dblp.cx
```

Reproduce the evidence (refuses to overwrite existing files; move the old
ones to `archive/` first):

```
python3 research/r1_skyline_index_20260918/run_final.py graphs/ca-AstroPh.edges graphs/ca-CondMat.edges graphs/cit-HepPh.edges graphs/loc-Brightkite.edges graphs/soc-Epinions1.edges graphs/soc-Slashdot0902.edges graphs/com-youtube.edges graphs/soc-pokec.edges
python3 research/r1_skyline_index_20260918/run_final.py --tag more --only graphs/ca-HepTh.edges graphs/email-Eu-core.edges graphs/amazon-copurchase.edges graphs/dblp-coauthor.edges
python3 research/r1_skyline_index_20260918/run_buildonly.py
python3 research/r1_skyline_index_20260918/report_tables.py
```

Timing runs are single-threaded and must not run concurrently with other
load; the laptop numbers vary 2-3x between rounds of identical code when
the machine is busy.

## Using the index from code

```
#include "chain_index.hpp"
auto ix = chainindex::ChainIndex<uint64_t>::load("dblp.cx");   // T = the count width recorded in the file
uint64_t k = ix.value(v, s);                                    // kappa_s(v), internal label v
chainindex::ChainIndex<uint64_t>::Runs r; uint32_t node;
ix.community_runs(v, s, k, r, node);                            // O(1): head range, whole runs, tail range
std::vector<uint32_t> ids(r_total + 8); ix.expand(r, ids.data()); // explicit labels (8 spare slots)
bool same = ix.member(u, v, s, k);
```

Labels are the aligned internal labels; the permutation from input labels
is written next to the index file by the tool (`<out.cx>.perm`, 4 bytes
per vertex) and is unnecessary if the graph is stored in that order.

## Status

Complete locally. Open: runs on the servers (tods1/tods2) when they are
back; optional per-node vertex counts for O(depth) ladders (bytes versus
latency, a user decision); r >= 2 is closed (CHAINS.md Section 10: the
paper's size forest already is the chain structure). No production code
(`src/`) or paper text was changed.
