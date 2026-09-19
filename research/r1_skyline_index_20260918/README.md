# r = 1 All-Size Community Index (chain index)

The self-contained code lives in `src-r1index/` at the repository root
(index header, hierarchy, tool, flattened solver, drivers); this directory
is the research record: theory, design stages, evidence and the report.

Research line, 2026-09-18/19. One index for the whole (1, s)-nucleus
hierarchy of a graph over every clique size s: core values kappa_s(v) and
every (s, k)-nucleus as a vertex set. Exact; brute-force verified;
measured on 17 graphs of five families.

Headline (17 graphs of five families, one thread, RESULTS_FINAL.md
Sections 9 and 14): 1.3x-11.7x fewer bytes than one S tree per size
(median 4.95x; dblp 10.3x; the 4.05 M-vertex dense dblp-coauthor 11.7x;
the 1,005-vertex email-Eu-core 1.3x is the worst case); a community is
located in 3.5-59 ns and listed at 0.04-0.28 ns per vertex on all but the
most fragmented graphs (faster than a memcpy of the vertex list on 13 of
15 measured points); values 6-24 ns; build 15 ms
(GrQc) to 34 s (pokec) and 189 s (dblp-coauthor, 512-bit counts) with
157 MB (dblp) to 8.8 GB (dblp-coauthor) peak memory. On the servers 29
distinct graphs in total: 1.3x-48.7x, median 7.0x (web-uk-2005 48.7x,
ca-coauthors-dblp 45.8x, web-NotreDame 24.7x). com-lj, hollywood and orkut
exceed the solver's memory (clique trees of more than 350-435 GB; see
RESULTS_FINAL.md Section 15).

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
4. Values are stored as `double` (exact below 2^53; the solver still
   peels with exact integers and converts once at build time); node tops
   and residues at per-size byte widths (1/2/4-byte integers below 2^32,
   else the 8-byte double). File format CHAINX05; omega and sigma are
   16-bit, so clique sizes above 255 are fine.

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
| `tods1.json`, `tods2.json`, `tods1_big*.json` and their `-logs/` | server runs (20 + 7 graphs; the big-graph attempts), launchers `tods1_run.sh`, `tods2_run.sh`, `tods1_big*.sh` |
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
using Index = chainindex::ChainIndex<double>;                    // files written by the tool hold double values
auto ix = Index::load("dblp.cx");
double k = ix.value(v, s);                                       // kappa_s(v), internal label v (exact below 2^53)
Index::Runs r; uint32_t node;
ix.community_runs(v, s, k, r, node);                             // O(1): head range, whole runs, tail range
std::vector<uint32_t> ids(Index::total(r) + Index::kSlack); ix.expand(r, ids.data());   // explicit labels
```

`member(u, v, s, k)` and `ladder(v, s, out)` exist as derived queries
(the selftest uses them as cross-checks); they are not part of the
measured interface.

Labels are the aligned internal labels; the permutation from input labels
is written next to the index file by the tool (`<out.cx>.perm`, 4 bytes
per vertex) and is unnecessary if the graph is stored in that order.

## Status

Complete; measured on the laptop and both servers. r >= 2 is closed (CHAINS.md Section 10: the
paper's size forest already is the chain structure). No production code
(`src/`) or paper text was changed.
