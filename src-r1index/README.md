# src-r1index: the chain index for r = 1 clique-core hierarchies

Self-contained C++23 code of the chain index: one index for a graph's
(1, s)-nucleus hierarchy over every clique size s, answering the value
query (core value of a vertex at a size) and the community query (the
(s, k)-nucleus containing a vertex, as label ranges) from a partition of
the vertices into hierarchy-equivalence chains.  Theory, design record
and the measured evidence are in `research/r1_skyline_index_20260918/`
(THEORY.md, CHAINS.md, RESULTS_FINAL.md, README.md); the paper draft is in
`Sigmod2027ChainIndex/`.  This directory contains everything needed to
build and run, with no include outside it.

## Layout

| Path | What |
|---|---|
| `chain_index.hpp` | the index: `chainindex::ChainIndex<V>` (V = double in the files the tool writes), build form and compact form, queries, save/load; header-only |
| `hierarchy.hpp` | per-size canonical trees from one row of core values (`make_tree_row`), brute-force helpers for the selftest, input preparation |
| `chain_index_tool.cpp` | `--selftest`, `--build-only <graph>`, `--bench <graph> <out.cx>` |
| `solver/` | the all-size vertex peel (pivoting clique tree, per-size bucket peel with cross-size bounds, replay and certificates), flattened from the research lines r1_theory .. r1_terminal; `terminal.hpp` is the entry, `solve(..., sink)` streams one row of core values per size |
| `scripts/run_final.py` | Release and sanitizer builds, both selftests, `--bench` on graphs given on the command line (`--tag NAME --only <graphs>`); writes `NAME.json` and `NAME-logs/` next to the script |
| `scripts/run_buildonly.py` | `--build-only` sweep (bytes, build phases, resident memory) |
| `scripts/report_tables.py` | prints the result tables from a directory of evidence JSON files |
| `CMakeLists.txt` | target `chain_index_tool`; `-DSANITIZE=ON` for ASan/UBSan |

Requirements: a C++23 compiler (Apple clang 21 and GCC 11.4 are tested),
CMake 3.20+, Boost multiprecision headers (`/opt/homebrew/include` on
macOS is on the include path; `/usr/include` on Linux), Python 3.10+ for
the scripts.  Graph files: first line `n m`, then one edge `u v` per line,
vertices 0..n-1.

## Build and check

```
cmake -S src-r1index -B src-r1index/build -DCMAKE_BUILD_TYPE=Release -DSANITIZE=OFF
cmake --build src-r1index/build -j 12
src-r1index/build/chain_index_tool --selftest
src-r1index/build/chain_index_tool --build-only data/com-dblp.edges
src-r1index/build/chain_index_tool --bench data/com-dblp.edges /tmp/dblp.cx
```

The selftest compares every value and community query against brute
force on every labelled graph with at most 6 vertices, 200 random graphs
on 7 to 10 vertices, split graphs, K8 and K300, in four storage forms and
with exact-integer and double values (4,844,536 community queries).  The
index this tree produces is byte-identical to the one recorded in the
research evidence (`final.json`, sha256 of `com-dblp.cx`).

## Query workload (fixed rule)

`--bench` draws its queries with a fixed generator (seed 20260918); this
rule is part of the protocol and is not changed between runs or machines.

- Community queries (v, s, k): v uniform among the active vertices
  (omega(v) >= 2); s uniform in [2, omega(v)] for that v; k in three
  regimes, each with its own draw: own (k = kappa_s(v); 20,000 queries),
  half (k = max(1, floor(kappa_s(v) / 2)); 20,000), root (k = 1; 1,000).
- Value queries (v, s): v uniform among all vertices, s uniform in
  [2, omega(v) + 1] (the last value returns 0); 200,000 queries.
- Timing: one warm-up pass, then five passes; the median is reported.
  Community time is split into locate (climb and entry reads), copying
  the range list, and expanding to explicit labels; the output size in
  vertices and ranges is recorded with it.

## Using the index

```
#include "chain_index.hpp"
using Index = chainindex::ChainIndex<double>;
auto ix = Index::load("dblp.cx");                       // one file, magic CHAINX05
double k = ix.value(v, s);                               // kappa_s(v) for an internal label v (exact below 2^53)
Index::Runs r; uint32_t node;
ix.community_runs(v, s, k, r, node);                     // O(1) after the climb: head range, whole runs, tail range
std::vector<uint32_t> ids(Index::total(r) + Index::kSlack);
ix.expand(r, ids.data());                                // explicit labels, eight per vector store
```

Labels are the index's aligned labels; the tool writes the permutation
from the labels of the input file next to the index file (`<out.cx>.perm`,
4 bytes per vertex: `perm[file label] = internal label`; before 2026-09-21
the file held the permutation from the tool's internal degeneracy-order
labels instead).  It is unnecessary if the graph is stored in the index's
order.

## Limits

Values are doubles: exact below 2^53, nearest double above; the
decomposition itself is computed with exact integers of 64 to 512 bits,
chosen from the graph.  The solver materialises every row of its clique
tree; graphs whose tree exceeds memory (com-lj, ca-hollywood-2009 and
com-orkut on a 503 GB machine) cannot be built.  Membership and ladder
queries exist as derived functions (`member`, `ladder`) and are checked by
the selftest but are not part of the measured interface.
