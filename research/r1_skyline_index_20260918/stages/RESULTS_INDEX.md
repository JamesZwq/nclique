# Skyline Index, Stage 2: Bytes And Query Latency Against S Trees

Date: 2026-09-18. In-memory implementation of both designs on the five
local full-range inputs, one machine, one thread, warm cache. No disk
format, no production, paper or default change.

## 1. Problem Summary

Build, from the same canonical trees, (a) the S-tree baseline with the
DFS-interval layout and (b) the skyline index of [THEORY.md](../THEORY.md),
for r = 1 over all sizes s = 2..d+1; verify both against brute force and
against each other; measure exact bytes per block and the latency of
community, membership and value queries.

## 2. Current Baseline And Bottlenecks

Baseline: per size, canonical nodes in DFS preorder with own classes
placed first, so a community is one contiguous slice of a per-size
class array; per active (class, size) pair one node id (CSR over
classes). Both designs share the class maps and Block D (values).

## 3. Candidate Algorithm Ideas

Only the two designs of THEORY.md Sections 4 and 7 are at stake here.
Section 13 of THEORY.md settled that empty canonical nodes stay.

## 4. Chosen Approach And Rationale

`index.cpp` includes the validated stage-1 source (`count.cpp`, main
renamed) so that `make_tree`, `shadow`, `twins` and the selftest graph
set are shared byte for byte. The Codex draft was reviewed and rewritten:
its baseline stored one node id per (class, size) in a dense
(s_max+1) x classes array (739 KB of its 863 KB on GrQc), it recovered
node creators by scanning all vertices per node and skyline bucket
offsets by scanning all classes per node (both quadratic, which is why
com-dblp never finished), its selftest compared the two designs only,
not brute force, and its baseline membership called the skyline locate.
All of that is replaced; the k = 1 (root) regime is timed on 1,000
queries instead of 6,667 because those communities hold up to 230,000
vertices.

## 5. Complexity Discussion

Build: one all-size decomposition, one union-find pass per size over the
valid rows, one shadow per uncertified (class, size), one tree walk per
canonical node for the cross-size pointer. Queries: THEORY.md Section 5.
Whole runs took 4 s (GrQc) to 559 s (dblp), almost all of it in the
community batches, whose outputs are 10^3 to 2.3 x 10^5 vertices per query.

## 6. Implementation Summary

`index.cpp` (both designs, selftest, measurement), `common.hpp`
(includes only; the construction lives in `count.cpp`), `run_index.py`
(Release and ASan/UBSan builds, `index --selftest` and `count --selftest`
in both, five graphs serially under `/usr/bin/time -l`, hashes, refusal to
overwrite), `index.json` (all fields), `index-logs/`.

## 7. Correctness Validation Summary

| Check | Release | ASan/UBSan |
|---|---:|---:|
| Selftest graphs (all labelled <= 6 vertices, 200 random 7..10, split, complete) | 34,075 | 34,075 |
| Community queries (both designs) against brute-force nuclei, every v, s, k | 605,476 | 605,476 |
| Membership answers against brute force | 3,651,590 | 3,651,590 |
| Value answers against the core matrix (including s > omega) | 1,439,834 | 1,439,834 |
| `count --selftest` unchanged after the include refactor | 34,075 / 844,230 / 93,057 | same |

On the five real graphs every timed community query (41,000 per graph)
and membership query (20,001 per graph) was answered identically by both
designs in a separate pass before timing; the core matrix equals the
frozen `fixed_sparse` control; F2 and F5 hold on every class.

## 8. Experimental Setup

Local macOS arm64, Clang from /opt/homebrew, `-O3 -DNDEBUG`, one thread,
seed 20260918. Queries: v uniform over vertices with omega >= 2, s uniform
in [2, omega(v)]; k = kappa_s(v) ("own", 20,000 queries), k = max(1,
kappa_s(v)/2) ("half", 20,000), k = 1 ("root", 1,000); membership 20,001
queries with a uniform u over all vertices, k cycling the three regimes;
values 200,000 queries with s uniform in [2, omega(v)+2]. One warm-up
pass, then 5 passes; the median pass time divided by the batch size.
Community outputs are written into preallocated buffers in both designs
and expanded to vertices. Commands: `python3
research/r1_skyline_index_20260918/run_index.py` from the repository
root; per-command lines are in `index-logs/*.log` and `index.json`.

## 9. Experimental Results

Bytes of the class-based layout as implemented (W = count width in bytes;
node records W+12 baseline, W+16 skyline; ids 4 B; gamma and size bytes;
offsets 8 B in the location and Block D CSRs):

| Graph | W | classes | active pairs | skyline entries | nodes | base, no D | skyline, no D | ratio | base, with D | skyline, with D | ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 8 | 4,105 | 12,704 | 4,780 | 1,517 | 206,756 | 185,924 | 1.11 | 256,238 | 235,406 | 1.09 |
| ca-HepPh | 32 | 9,559 | 114,694 | 13,909 | 3,538 | 1,245,768 | 547,422 | 2.28 | 1,532,406 | 834,060 | 1.84 |
| com-dblp | 16 | 256,983 | 928,980 | 316,949 | 33,979 | 12,975,764 | 10,146,478 | 1.28 | 16,880,402 | 14,051,116 | 1.20 |
| web-Stanford | 8 | 270,303 | 1,429,928 | 721,162 | 53,239 | 16,921,860 | 14,410,400 | 1.17 | 26,086,394 | 23,574,934 | 1.11 |
| amazon0302 | 8 | 253,766 | 787,990 | 457,900 | 65,134 | 11,733,624 | 11,803,352 | 0.99 | 17,079,748 | 17,149,476 | 1.00 |

Skyline block breakdown (bytes): shared class maps / nodes / entries /
location / cross-size lists:

| Graph | shared | nodes | entries | location | cross |
|---|---:|---:|---:|---:|---:|
| ca-GrQc | 58,360 | 36,408 | 23,900 | 56,748 | 10,508 |
| ca-HepPh | 134,304 | 169,824 | 69,545 | 146,025 | 27,724 |
| com-dblp | 3,564,576 | 1,087,328 | 1,584,745 | 3,640,617 | 269,212 |
| web-Stanford | 3,336,440 | 1,277,736 | 3,605,810 | 5,768,242 | 422,172 |
| amazon0302 | 3,111,956 | 1,563,216 | 2,289,500 | 4,319,636 | 519,044 |

Latency, ns per query, baseline / skyline, with the mean output size:

| Graph | own | half | root | membership | value (shared) |
|---|---|---|---|---|---:|
| ca-GrQc | 2,903 / 6,589 (1,657 vertices) | 3,449 / 7,962 (2,113) | 4,200 / 9,664 (2,296) | 4.3 / 19.0 | 4.3 |
| ca-HepPh | 9,668 / 22,344 (5,332) | 10,014 / 26,045 (6,227) | 11,129 / 29,247 (7,102) | 13.0 / 63.3 | 5.4 |
| com-dblp | 370,861 / 824,832 (164,595) | 442,375 / 983,270 (199,753) | 498,255 / 1,109,001 (229,884) | 14.1 / 46.3 | 4.6 |
| web-Stanford | 163,655 / 478,185 (82,833) | 209,661 / 578,905 (102,620) | 256,140 / 693,421 (126,320) | 89.3 / 93.2 | 8.2 |
| amazon0302 | 223,981 / 486,636 (107,061) | 326,479 / 738,355 (151,318) | 352,531 / 808,025 (165,616) | 9.7 / 26.8 | 8.5 |

Per output vertex the baseline costs 1.6 to 2.3 ns and the skyline 3.8
to 5.8 ns in every regime on every graph.

## 10. Analysis Of Runtime And Memory

Memory. The implemented layout shrinks the index by 1.11x to 2.28x on
four graphs and not at all on amazon0302 (0.99x). Three effects explain
the gap to the word model of stage 1 (1.25x to 3.34x):

- The class maps are shared and large: 3.1 to 3.6 MB on the three big
  graphs, 30 to 35 percent of the skyline total, because at r = 1 the
  classes are nearly singletons (1.03 to 1.28 vertices per class). The
  maps cost about 12 bytes per vertex and save less than that in
  per-pair or per-entry data on every graph. A per-vertex layout without
  classes is therefore smaller for BOTH designs; computed from the same
  counts (E and active pairs per vertex, identical node and cross-size
  bytes, 4-byte per-vertex offsets):

| Graph | base, no D | skyline, no D | ratio | base, with D | skyline, with D | ratio |
|---|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 211,672 | 127,078 | 1.67 | 251,592 | 166,998 | 1.51 |
| ca-HepPh | 1,587,412 | 409,364 | 3.88 | 1,851,592 | 673,544 | 2.75 |
| com-dblp | 12,195,200 | 6,398,064 | 1.91 | 15,441,844 | 9,644,708 | 1.60 |
| web-Stanford | 15,056,308 | 10,656,644 | 1.41 | 23,725,810 | 19,326,146 | 1.23 |
| amazon0302 | 8,951,056 | 7,817,848 | 1.14 | 13,365,766 | 12,232,558 | 1.09 |

  These per-vertex numbers are arithmetic on measured counts, not a
  second implementation.
- The location block duplicates the entry ids (one copy in the node
  buckets, one per class) and carries 8-byte offsets; with 4-byte
  offsets it loses 1 MB on dblp. The duplicate is what makes locate()
  a scan of a short per-class list instead of a search.
- Block D is common to both designs and is 23 to 39 percent of the
  totals on Stanford and amazon, where 43 to 54 percent of the cells
  are residue.

Speed. The skyline index is slower on every query class except values,
which are shared: community listing 2.1x to 2.9x slower per query and
2.3x to 2.9x per output vertex (scattered bucket prefixes across the
admissible nodes of every larger size against one contiguous slice),
membership 1.0x (Stanford, where the deep climb dominates both) to 4.4x
slower (the location scan and the cross-size chain against one array
read). Community queries are output-bound in both designs: even at a
vertex's own level the mean community has 1.7 x 10^3 to 1.6 x 10^5
vertices, because low sizes have giant nuclei.

## 11. Failure Cases / Limitations

- amazon0302: no memory gain with classes, 1.09x to 1.14x without, and
  2x slower queries; a loss on both axes.
- The gain is bounded by avg_traj/mu and by the node and location
  overheads; it is largest where trajectories are long and certified
  (HepPh 2.3x to 3.9x) and small where trajectories are short.
- Single machine, warm cache, in-memory arrays, no disk format, one
  seed, medians of five passes; latency numbers are indicative, not a
  benchmark protocol with repeated processes.
- Stanford's tree depth (2,139) makes membership 89 ns in the baseline
  too; jump pointers were not implemented.

## 12. Final Conclusion

Both designs are correct against brute force and each other. The
skyline index trades speed for space: with the implemented class layout
it saves 9 to 55 percent of the bytes on four graphs and nothing on
amazon0302, and with a per-vertex layout it would save 9 to 64 percent
(computed); community listing costs 2 to 3x and membership 1 to 4x the
baseline's time. The S-tree baseline with DFS intervals remains the
fastest exact representation, at 30 to 60 bytes per vertex on these
graphs. There is no configuration in which the skyline index is both
smaller and faster.

## 13. Next Recommended Improvements

1. If the memory saving is wanted, drop the twin classes at r = 1 and use
   4-byte offsets; measure the per-vertex layout rather than computing it.
2. For latency, the baseline is the reference; the only skyline cost that
   could be removed is the location scan (store the own-size node id at
   the vertex's first skyline size directly, one extra word per vertex).
3. The decision between the two is the user's: bytes against nanoseconds.
