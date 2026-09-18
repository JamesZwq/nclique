# Stage 2: Build Both Indexes, Measure Bytes And Query Latency

Date: 2026-09-18. Authorized by the user after the counting gate passed
("可以的 继续"). Scope: one C++ program that builds, in memory, the S-tree
baseline and the skyline index of [THEORY.md](THEORY.md) for r = 1 over all
sizes, checks both against brute force and against each other, and
measures exact bytes and query latency on the five local inputs. No
production, paper or default change. Everything stays in this directory.

Reuse: `count.cpp` already contains the correct construction of the
canonical trees (`make_tree`), the shadow, the twin classes, the skyline
and certification computation, and the selftest. Move those pieces into
`common.hpp` without changing behaviour, make `count.cpp` include it,
and prove equality by rerunning `count --selftest` (same JSON summary:
34,075 graphs, 844,230 cells, 93,057 partitions) and `count --graph` on
ca-GrQc (same JSON line as in `counts.json`). Then write `index.cpp`.

## Shared layout (both designs)

Per size s, the canonical nodes of T_s numbered in DFS preorder of the
forest. Within a node, its own attached classes come first, then the
children's ranges, so that the classes of a node and all its descendants
are one contiguous slice of a per-size class array. Node record:
`top` (W bytes, k_hi), `parent` (uint32, absent for roots), `size`
(uint32, number of nodes in the subtree including itself), `bucket`
(uint32, offset of the node's own slice in the per-size class array;
equal to the next node's offset when the node has no own classes).
"Own classes of N at size s" means, for the baseline, all classes c with
X_s(c) = N and kappa_s(c) > 0; for the skyline index, only those with
s in Sky(c). W is the count width chosen by `width()` (64/128/256/512).

## Baseline: S trees

- Per size s: node records as above.
- Per size s: the class array in DFS order (uint32 class ids), one entry
  per active (class, s) pair.
- Per active (class, s) pair: `leaf[s][c]` = id of X_s(c) (uint32),
  stored as one CSR over classes (offset per class, then one id per
  active size, sizes 2..omega(c)).
- Values: see Block D below, identical in both designs.

Queries:
- Community (v, s, k): c = class(v); if s > omega(c) return empty;
  N = leaf[s][c]; while parent(N) exists and top(parent(N)) >= k:
  N = parent(N). Output the class-array slice [bucket(N),
  bucket(N + size(N))) (the offset of the node just after the subtree, or
  the array end), expanded to vertices through the class-to-vertex CSR.
- Membership (u, v, s, k): find N as above from v; d = class(u); if
  s > omega(d) return false; return id(leaf[s][d]) in [id(N), id(N)+size(N)).
- Value (v, s): Block D.

## Skyline index

- Per size s: node records as above, plus `ahi` (uint32, id in T_{s-1}
  of the container A_hi, absent at s = 2) and the reverse cross-size CSR:
  for every node Z of T_s the ids of nodes M of T_{s+1} with ahi(M) = Z,
  ordered by id(Z) so that the lists of an id range are one slice.
- Per size s: the entry array in DFS order of nodes, one entry per class c
  with s in Sky(c), placed in the slice of X_s(c) and sorted within the
  slice by gamma = prev(c, s) ascending; an entry is a uint32 class id
  plus a uint8 gamma stored in a parallel byte array (gamma < 256 holds
  because s_max <= 240 here; assert it).
- Location: one CSR over classes; for class c the ids of X_s(c) for s in
  Sky(c) ascending, and, in a parallel uint8 array, the size s of each
  entry (needed to find nxt(c, s) without per-size data).
- Values: Block D.

Construction of ahi: THEORY.md Section 6 Step 3, using the creator vertex
of each node (already recorded by `make_tree`): l = sigma_s(top(M)),
start at leaf_s[creator(M)], walk up while top(parent) >= l. Cache
sigma by (s, value).

Queries (THEORY.md Section 5; implement literally):
- Community (v, s, k): c = class(v). Scan Location(c) for the first entry
  with size s' >= s; if none, return empty. X = that node id (in T_{s'}).
  Repeat s' - s times: X = ahi(X) (each step moves to the previous size).
  N = X (now in T_s). Walk up while parent exists and top(parent) >= k.
  Then enumerate: at size s, output every entry of the slice
  [bucket(N), bucket(N + size(N))) (no gamma filter is needed at size s
  itself, but applying gamma < s is harmless). Let Adm = [id(N), id(N) +
  size(N)) as an id range of T_s. For t = s, s+1, ...: take the reverse
  lists of the nodes of Adm (for t = s one contiguous slice; afterwards
  one slice per node), that is the list of T_{t+1} nodes whose ahi lies
  in Adm; for each such node M output the entries of M's own slice with
  gamma < s (stop the scan at the first gamma >= s: the slice is sorted);
  the next Adm is the list of those M. Stop when Adm is empty or t+1 is
  the last size. Expand classes to vertices as in the baseline.
- Membership (u, v, s, k): find N for v as above; d = class(u); locate
  X_s(d) by the same Location scan and ahi chain; return id(X_s(d)) in
  [id(N), id(N)+size(N)).
- Value (v, s): Block D.

## Block D (values, both designs)

Per class: omega (uint8 is enough here, assert), sigma (uint8), and the
residue values kappa_s(c) for 2 <= s < sigma(c) stored explicitly as W
bytes each in one CSR over classes (Variant A). Value query: if s >
omega: 0; if s >= sigma: C(omega-1, s-1) computed with the existing
`Combinations` table; else the stored value. Also report the bytes
Variant B would need (only skyline sizes s < sigma) as a number, without
implementing its reconstruction.

## Measurement protocol

Per graph, after building both indexes from the same trees:
1. Bytes: for each design list every array with element size times
   length, and the totals: node records, per-pair or per-entry data,
   location, cross-size CSR, class maps (shared; report once), Block D
   (shared; report once). Print the two totals with and without Block D.
2. Correctness: on every timed query below, run both designs and require
   identical vertex sets (sort and compare) and identical membership and
   value answers. Do this in a separate pass from the timing.
3. Queries: seed 20260918; 20,000 community queries: v uniform over
   vertices with omega >= 2, s uniform in [2, omega(v)], k drawn as one
   of kappa_s(v), max(1, kappa_s(v)/2) (integer halving of the W-byte
   value), 1, cycling through the three; 20,000 membership queries with
   an extra uniform u; 200,000 value queries (v uniform, s uniform in
   [2, omega(v)+2] so that zeros and closed forms are exercised).
   Time each design separately: one warm-up pass, then 5 repetitions of
   the whole batch with std::chrono::steady_clock; report the median of
   the 5 batch times divided by the number of queries (ns/query), and
   for community queries also the total output vertices per batch and
   ns per output vertex. Community queries must write the output into a
   preallocated buffer in both designs (same materialisation cost).
   Membership and value queries must accumulate their answers into a
   checksum that is printed, so that the compiler cannot drop them.
4. Print one JSON line per graph with all of the above; `run_index.py`
   runs the five graphs serially from the repo root with
   `/usr/bin/time -l`, writes `index.json` and `index-logs/`, records
   sha256 of sources and inputs, and refuses to overwrite evidence.
5. Selftest `index --selftest`: on every labelled graph with <= 6
   vertices and 200 random 7..10-vertex graphs, for every v, every s in
   [2, omega(v)], and every k in {1..kappa_s(v)}, both community answers
   must equal the brute-force (s, k)-nucleus containing v (compute
   nuclei as in count.cpp's selftest); every membership answer must
   equal brute force for every u; every value must equal the core
   matrix. Also run under ASan/UBSan (`ASAN_OPTIONS=halt_on_error=1`,
   no detect_leaks on macOS).

## Deliverables

`common.hpp`, `count.cpp` (now including common.hpp, behaviour
unchanged), `index.cpp`, `CMakeLists.txt` (two targets), `run_index.py`,
`index.json`, `index-logs/`, `RESULTS_INDEX.md` in the 13-section
AGENTS.md format with: bytes table per graph and design (with and without
Block D), latency table (ns/query and ns/output vertex, both designs,
three k regimes separated), selftest summary, exact commands, and honest
limits (single machine, warm cache, in-memory index, no disk format).
Update the row of this directory in `../SUMMARY.md`. Commit with the
usual trailer.

## Rules

- Build with at most `-j 12`. Serial runs. One thread.
- No claim beyond what the JSON shows. Report every regression: if the
  skyline index is slower on some query class, say so with the numbers.
- If any correctness check fails, stop, keep the failing case, report.
- Do not touch files outside this directory except the SUMMARY.md row.
