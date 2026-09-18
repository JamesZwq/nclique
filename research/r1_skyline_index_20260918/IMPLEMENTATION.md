# Measurement Gate: Counts Before Any Index Prototype

Date: 2026-09-18. Authorized by the user ("开始吧") after the theory note.
Scope: a C++ counting program plus its selftest, run on the five local
full-range inputs. NO index implementation, NO query benchmark, NO change
to production, paper, or any other research directory.

What to count is Section 11 of [THEORY.md](THEORY.md). The program lives
in this directory as `count.cpp`, built by its own `CMakeLists.txt` in
the style of `../r1_terminal_20260918/CMakeLists.txt` (C++23, Boost from
/opt/homebrew/include, optional SANITIZE). It reuses the terminal harness:
`#include "../r1_terminal_20260918/terminal.hpp"` (which pulls the whole
allsize/fullrange/orderreplay stack), `using namespace orderreplay;` and
`#include "../r1_orderreplay_20260917/shared_harness.inc"` exactly as
`../r1_terminal_20260918/main.cpp` does, so `prepare(path)` gives
`Input{graph, ordinary, d, load_ms, order_ms}` with the graph already
degeneracy-relabelled, `Layout` gives `count_bound` and `width` selects
the count type T among uint64_t, unsigned __int128, uint256_t, uint512_t.

## Inputs

`data/ca-GrQc.edges`, `data/ca-HepPh.edges`, `data/com-dblp.edges`,
`graphs/web-Stanford.edges`, `graphs/amazon0302.edges`, run from the repo
root exactly like `../r1_terminal_20260918/run.py` (which passes the
graph path as argv and uses cwd = repo root). Full range: s_max = d + 1.

## Pipeline per graph

1. `Input input = prepare(path)`; `maximum = max(2, d + 1)`;
   `bits = width(count_bound(graph, Layout(graph, maximum), d))` computed
   the same way as `--paired` in the terminal main; dispatch on T.
2. Build the plain terminal index: `terminal::Index index(maximum);
   terminal::build(graph, index, 0); index.prepare(graph.n)`. Mode 0
   only (rows are (H, Q), no choice sets). Row p has members
   `index.members[row.begin, row.hold_end)` = holds and
   `[row.hold_end, row.pivot_end)` = optional Q; `row.valid(s)` is the
   size filter; `index.touching(v)` yields codes `(p << 2) | role`,
   role 0 = hold, 1 = Q.
3. Core matrix: `auto out = terminal::Solver<T>::solve(graph, index,
   choose, input.ordinary)` with `choose = Kernel<T>::Combinations(d+1,
   maximum)`; `out.common.data.core[s * n + v]` is kappa_s(v) for
   s = 2..maximum (row 0 and 1 unused). Verify against
   `Kernel<T>{}.fixed_sparse<true>(layout, n, choose, ordinary).core`
   once per graph (they must be equal; this is the frozen control).
4. Twin classes: two vertices are twins iff their closed neighbourhoods
   are equal. Rows of `graph.row(v)` are sorted. Hash (degree, row with v
   inserted) then compare exactly within a hash bucket. Output n_cls,
   the class-size histogram summary (max, mean), and the class label
   array. Assert that twins have identical core rows for every s (F5).
5. omega(v) = max{s : kappa_s(v) > 0} (0 or 1 if none); assert support
   nesting: kappa_s(v) > 0 for all 2 <= s <= omega(v) (F1).
6. Shadow and skyline. Implement sigma_s(k) with cpp_int: s-cascade by
   greedy binary search on a with C(a, s) <= remainder, then the shadow
   sum. Unit-test it: sigma_s(C(a, s)) == C(a, s-1) for a in [s, 40];
   sigma_2(3) = 3 (3 = C(3,2), shadow C(3,1)); sigma_2(4) = 4
   (4 = C(3,2) + C(1,1), shadow C(3,1) + C(1,0) = 3 + 1);
   sigma_3(21) = 21 (21 = C(7,3), shadow C(7,2) = 21);
   sigma_3(22) = 23 (22 = C(7,3) + C(2,2), shadow C(7,2) + C(2,1) = 21 + 2);
   check sigma_s(k) nondecreasing on a sweep k = 0..200 for s = 2..6, and
   compare against a brute-force shadow: for k <= 60 and s <= 4 take the
   first k s-subsets of {1..12} in colex order and count their distinct
   (s-1)-subsets; that count must equal sigma_s(k) (Kruskal-Katona is
   tight on colex initial segments).
   Then per vertex, for 2 <= s < omega(v):
   delta_s(v) = kappa_s(v) - sigma_s(kappa_{s+1}(v)); assert >= 0 (F2,
   record the count of checked cells). Cache sigma by (s, value).
   Sky(v) = {s : s = omega(v) or delta_s(v) > 0}.
   Certification: sigma(v) = least s with kappa_s(v) == C(omega(v)-1, s-1)
   (compare with cpp_int), or omega(v)+1 if none; assert that once equal
   it stays equal for all larger s <= omega(v) (F8).
7. Connectivity and canonical trees, per s = 2..maximum, exactly the
   construction of THEORY.md Section 6, Step 1, mode 0:
   - active vertices: kappa_s(v) >= 1, processed by level k descending
     (bucket by value; T may be 128/256/512-bit, so sort the active
     vertices by value descending with std::sort and sweep groups).
   - per row counters activeHolds[p], activeQ[p] (uint32), live flag,
     representative vertex rep[p]. Activating v: for each code in
     `index.touching(v)`: skip rows with !row.valid(s); bump the counter
     for its role; if not yet live and activeHolds == row.holds() and
     row.holds() + activeQ >= s: set live, rep = v, union v with every
     active member of the row (scan members, test active[]); else if
     live: union(v, rep).
   - DSU with path compression + union by size.
   - canonical nodes: keep curNode[root] (absent for fresh singletons)
     and a per-level touched list. On activating v: root(v) is touched.
     On union(ra, rb) with survivor r: touched; collect into
     pendingChildren[r] the current nodes of ra and rb that exist and
     are not already collected (use a stamp per node). At the end of the
     level: for each touched root create a node {k_hi = k (store as T),
     parent = absent, children = pendingChildren}, set the parent of each
     child to the new node, curNode[root] = new node; clear pending.
     leaf_s[v] = curNode[root(v)] read at the end of level kappa_s(v)
     (that is, right after the level that activated v).
   - after the layer: N_T(s) = number of nodes, max depth (parent chain),
     nodes with exactly one child (chain nodes), and assert L1:
     k_hi(leaf_s[v]) == kappa_s(v) for all active v, and that every node
     with children has k_hi strictly below the children's k_hi.
   - empty nodes: after the skyline is known, count nodes N of T_s such
     that no v with s in Sky(v) has leaf_s[v] == N.
   - optional check of L2 (do it; it is cheap): for each active v with
     s < omega(v) and delta_s(v) == 0, compute l = sigma_s(kappa_{s+1}(v))
     (= kappa_s(v)) and verify that walking up from leaf_s[u], where u is
     the vertex that created node leaf_{s+1}(v) (record creator[node] =
     the vertex at whose level the node was made; any vertex activated in
     that level inside that component works), while parent.k_hi >= l,
     ends at leaf_s[v]. Keep leaf arrays for two consecutive layers.
8. Counts to report per graph (JSON line + Markdown table):
   n, m, d, s_max, count_bits, n_cls, max class size,
   active_pairs_vertex = sum_v (omega(v)-1) over omega >= 2,
   active_pairs_class = the same summed once per class,
   E_vertex = sum_v |Sky(v)|, E_class = once per class,
   avg_traj = active_pairs_vertex / n_active, mu = E_vertex / n_active,
   and the class versions,
   N_T total and per-s max, empty_nodes total, chain_nodes total,
   max_depth over s,
   residue_cells = sum_v (sigma(v) - 2), residue_zero_delta = residue
   cells with s < omega(v) and delta_s(v) == 0, certified_cells =
   active_pairs_vertex - residue_cells, F2_checked, F2_violations (must
   be 0), L1 and L2 checks passed (counts),
   words_S_trees = active_pairs_class + 4 * N_T (positions + node words),
   words_skyline = E_class + E_class/4 + 6 * N_T (entries + gamma bytes
   as quarter words + node words including A_hi and reverse list),
   verdict = "skyline smaller" iff words_skyline < words_S_trees,
   plus the raw inequality n_cls*(avg_traj_class - mu_class) versus
   2*N_T + mu_class*n_cls/4 from THEORY.md Section 7.
   Also per-graph time and peak RSS of the count run (/usr/bin/time -l
   on macOS or getrusage), for the record only.
9. `--selftest`: exhaustive labelled graphs on up to 6 vertices, 200
   seeded random graphs on 7..10 vertices, the split and complete graphs
   used by the terminal selftest. For each graph and each s in
   [2, maximum]: compute the core matrix as above; compute ground-truth
   nuclei for every k in the set of positive values by brute force:
   W = {v : kappa_s(v) >= k}, s-cliques of G[W] by `bottomup::clique_masks`
   restricted to W, union-find over cliques; compare with the DSU
   components after level k of Step 7 (snapshot components at the end of
   each level and compare as sorted vertex partitions). Also compare the
   canonical intervals: for every node N, I(N) = {k : N is a nucleus at
   k} from brute force must equal [k_lo, k_hi] from the tree. Verify F2,
   F5, F8, L1, L2 there too. Exit 0 with a JSON summary on success.
   Run the selftest in Release and in the sanitizer build.

## Deliverables in this directory

- `count.cpp`, `CMakeLists.txt`, `run.py` (serial, one graph at a time,
  records commands, source hashes, stdout JSON lines to `counts.json`
  and logs to `counts-logs/`).
- `RESULTS.md` in the 13-section format of AGENTS.md, with the counts
  table, the build-test verdict per graph, the selftest summary, exact
  commands, and an honest statement that nothing about speed or index
  size has been measured.
- One row update for `../SUMMARY.md` (the row for this directory) with
  the outcome.

## Rules

- No timing claims. These are counts.
- Do not modify files outside this directory except the SUMMARY.md row.
- Build with `cmake --build <dir> -j 12` at most.
- If any assertion (F1, F2, F5, F8, L1, L2, control equality, selftest)
  fails, stop, keep the failing case, and report it; do not weaken the
  check.
