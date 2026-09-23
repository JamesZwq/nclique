# Merge plan: the VLDB r=1 paper + the chain index (2026-09-24)

User decisions (2026-09-24): the VLDB paper ("Efficient (1,s)-Nucleus Decomposition on Billion-Scale Graphs",
source ~/UNSW/pivoter-quotient-research/vldbNuclearR1) will never be submitted.  Everything it says must be said here,
following its storyline, with the chain index added on top.  Standing rules kept: the model is the clique core (never
"(1,s)-nucleus" or "the r=1 special case"; nucleus only named as the origin), theory-only paper, no em-dashes, SVG
figures of the stored index, every figure carries its baseline, 12 body pages (see open question Q1).

## One-sentence thesis
One static clique tree serves every clique size: it gives the core values of all sizes by change-only peeling that
skips the vertices already on their floor, the hierarchy of every size without re-enumerating a clique, and an index
that stores each vertex once for all sizes and returns every community as a few label ranges.

## Storyline (VLDB's structure: challenges <-> measured limitations <-> contributions)
Problem: answer value and community queries of the clique core at every size s and level k.
The one primitive: set support c_s(v, R) (s-cliques through v inside R).  The one invariant: the clique tree is static
under vertex removal (a path loses a hold -> dead, loses a pivot -> counter - 1).

Challenges (each with a measured number on the baseline):
1. Residual state during peeling (VLDB): CND rewrites paths and a reverse map per pop (17.9x gmean peel gap).
2. Memory of the residual index (VLDB): 1.55x gmean, CND fails at s >= 3 on com-friendster.
3. Construction as the new bottleneck (VLDB): web-it build 137 s vs peel 59 ms at s = 3; and across sizes one
   decomposition per size (CND 15.7 h for the 499 sizes of web-uk).
4. Hierarchy emission from a mutated index (VLDB): level-DSU up to 43x memory, OOM on wiki-Talk s >= 9.
5. Repetition across sizes (ours): S trees 190 MB on web-uk, a vertex stored once per size.
6. Contiguity of communities (ours): a community is not one block at every size.

Contributions (final deliverables only):
- Static clique tree + change-only peeling SPIN* (VLDB) extended across sizes with settled vertices (Theorem: peel
  with settled vertices): the construction of every size from one tree.
- ParaBuild: parallel clique-tree construction (VLDB).
- BuildHier: the tree of every size from the static clique tree by merge events + union-find (VLDB).
- Chains (own nodes decide hierarchy equivalence, coarsest lossless partition) and the chain index (aligned labels,
  chain order, run arrays, floor and tail) with O(1) community location (ours).
- Experiments: peel, memory, phases, parallel build, billion-edge graph, hierarchy (VLDB) + construction of all
  sizes, index size, query latency, scalability (ours) + case studies.

## Section map
| # | Section | From VLDB | From the chain paper | New / to write |
|---|---|---|---|---|
| 1 | Introduction | challenges 1-4, contributions, running example | challenges 5-6, index idea, applications | merged storyline; one running example (ours, Figure 1) |
| 2 | Preliminaries | clique path (H, Q, eta), BuildCPI pseudocode, Sigma | clique core, own node, problem statement | one tree for all sizes: path size interval [lo, hi] |
| 3 | Baselines | mutable-CPI baseline (CND) section | S trees section | both as "what exists and what it costs" |
| 4 | Static clique tree | vertex-removal counter theorem | | stated for every size at once |
| 5 | Core values | set support lemma, SPIN (fixed point), SPIN* (change-only), closed-form loss lemma, correctness, complexity, trace table | settled vertices (Theorem 7.1), example | SPIN* across sizes with settled keys; complexity in terms of the incidences of unsettled paths |
| 6 | Construction | ParaBuild; BuildHier (events + elder rule, correctness, complexity) | Build (refine chains, layout) | Build as the pipeline calling SPIN*, BuildHier per size |
| 7 | Chains | | Section 4 (equivalence, coarsest partition, twins) | |
| 8 | The chain index | | Section 5 (labels, order, runs, values, size) | |
| 9 | Queries | | Section 6 | |
| 10 | Experiments | peel/memory vs CND, phases, ParaBuild scaling, dense synthetic, friendster, input scaling, hierarchy vs level-DSU | construction vs CND per size + settled effect, index size, size by part, design study, latency, latency by size, scalability | one setup block; decide which VLDB experiments are rerun on the chain pipeline |
| 11 | Case studies | granularity knob, ranking vs k-core, ego network (cross-r case dropped: conflicts with the clique-core framing) | best size per query (SNAP), Amazon zoom | pick 2-3 |
| 12 | Related work, conclusion | both | both | merge |

## Experiments: what exists, what must be rerun
- VLDB experiments ran SPIN* per single s on ten graphs (VLDB code, not chain_index_tool).  Numbers can be reused only
  if the paper says they are per-size runs of the same peel; the chain pipeline's peel is the same algorithm with
  settled vertices added.
- ParaBuild: exists in the VLDB code (par_src); chain_index_tool builds the tree single-threaded.  Either port it or
  report ParaBuild on the per-size tree as in VLDB.
- com-friendster: VLDB ran SPIN* at s = 2,3,4; the chain index on friendster was never built (feasibility unknown).
- Hierarchy vs level-DSU: VLDB bench on four graphs; our tree pass is a different (equivalent) algorithm in code.

## Open questions for the user
Q1. Page budget: VLDB alone was 15 pages, the chain paper is 12.  (a) keep 12 body pages, compress the VLDB material
    to definitions, lemmas, pseudocode and the key plots; (b) write the full merged version first, cut afterwards.
Q2. Which VLDB experiments to keep, and whether to rerun them on the chain pipeline (ParaBuild port, friendster).
Q3. Title and framing: the index paper's title, or a title covering decomposition + index.
