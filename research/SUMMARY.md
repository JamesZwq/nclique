# Research Summary And Navigation

Last updated: 2026-09-18. This is the entry point for ALL directories under
`research/`, including unsuccessful attempts. Read [project memory](../PROJECT_MEMORY.md)
and [research rules](../AGENTS.md) before starting another attempt.

Mandatory user workflow: [write the full proof before implementation](PROOF_PROTOCOL.md).
Every new rule needs explicit hypotheses, a forward derivation, both
directions of any claimed equivalence/invariant, boundary cases, safe
transitions and a proof gate. Tests check the implementation, not the theorem.

## Current Decisions

2026-09-18 direction change (user-driven, later the same day): the peel
line is closed for now; the target is the STORED index for r=1 across all
sizes, including community (nucleus) storage, which must be as compact as
the values. Measured verdict on the paper's maximal-clique forest as r=1
community storage: 2-27x LARGER than S explicit merge trees on seven
collaboration/citation/product graphs (cost sum_v |M(v)| against
sum_v top(v)); it wins only on web/mesh inputs. The user pointed to the
SIGMOD 2026 SGL bi-component index (~/Downloads/3802129.pdf); its
zero-redundancy retrieval transfers once dominance is defined through
iterated Kruskal-Katona shadows. [Theory](r1_skyline_index_20260918/THEORY.md)
proves the order, the chain/level-exactness lemmas that reduce SGL's
summary graph to per-size merge trees with one cross-size pointer per
node, the certified-tail collapse and the twin quotient, plus query and
construction correctness and an exact break-even against S trees. An
independent brute-force checker (1,501 graphs up to 11 vertices, 58,924
retrieval queries) found no violation. The user authorized the counting
gate of its Section 11 ("开始吧"); [spec](r1_skyline_index_20260918/IMPLEMENTATION.md).
2026-09-19 outcome: the skyline dedup measured as a space/time trade and was
dropped; the lever is the partition into hierarchy-equivalence CHAINS
([CHAINS.md](r1_skyline_index_20260918/CHAINS.md), Lemmas C1-C6). The final
module ([RESULTS_FINAL.md](r1_skyline_index_20260918/RESULTS_FINAL.md),
`chain_index.hpp`) stores the all-size hierarchy in 3.7x-7.8x fewer bytes than
one S tree per size, locates any community in O(1) (3-40 ns) and lists it at
memcpy speed or better on the large graphs; brute-force verified. Research
code only: no production, no paper edit; paper framing is the user's call.

2026-09-18 implementation authorization: the user explicitly requested C++
and measurements of the latest branch-factoring theory. The isolated
`r1_terminal_20260918` implementation includes full/partial grouping, a
packed unmerged control, and factored compiled replay. Release and sanitizer
V2 correctness suites pass. The [completed report](r1_terminal_20260918/RESULTS.md)
records 320 paired timings and 120 fresh processes on five full-range
inputs, all single-threaded. Old-active/new-partial-replay solver ratios
are 1.134 paired / 1.171 fresh; complete compute ratios are 1.127/1.187.
Against the best observed old control separately per graph, the paired
solver aggregate is almost neutral. Stanford has substantial work and
memory savings, but other gains are small, variable, or negative. Retain
this as a conditional option, not a Base or manuscript update. The pause
below remains the rule for unrelated new candidates.

Current workflow clarification after ordered replay (2026-09-17): pause
new prototypes and benchmarks. The next step is a theory-only assessment
covering correctness, dominant work removed, all additional costs,
break-even conditions, applicability and adverse examples. A correct
local shortcut is not automatically worth implementing. Follow the
separate benefit gate in [the proof protocol](PROOF_PROTOCOL.md) and
discuss the derivation before moving back to implementation.

Latest direction decision (2026-09-17): park complete-join DP as a validated
special case, including the proposed private-region detector. The user wants
more general theory, not more recognition of rare graph structures. Seek
relations valid on arbitrary graphs and general s that remove commonly
repeated computation; include recognition, state and maintenance costs.
General correctness alone does not establish broadly useful performance.

After the full d+1 experiment, the user judged the1.167x stream improvement
insufficient. Keep stream as a benchmark reference, not a major contribution.
Do not tune another queue/filter variant without a new argument for removing
substantial repeated counting or nonzero loss work; preserve all old evidence.

Previous implementation: [certified order replay](r1_orderreplay_20260917/RESULTS.md)
does remove actual nonzero support updates by certifying whole new size
layers from static forward counts. Its active-set implementation improves
on stream1.175x overall, but Stanford/Amazon regress and memory increases.
Only3.32% of all positive support decrements disappear across the five
full-range inputs. Retain the proof/code as a conditional research option,
not a Base replacement or the requested broadly large breakthrough.

Earlier theoretical follow-up: [local order repair](r1_orderreplay_20260917/LOCAL_REPAIR_THEORY.md)
proves a safe descending repair and precise count-reuse conditions, without
new code or experiments. Whole-block relabeling can preserve certificate
counts even when every value changes. Conversely, few failed vertices can
lead to a global cascade, and the specified uncapped full-rescan repair has
a quadratic-work family for every fixed s. Existing cross-size bounds can
settle that family; this is not an impossibility result for general DP.
Incremental repair's useful total-work bound is unresolved;
do not prototype it just because the transition is correct. Direct local
label relaxation overlaps established h-index nucleus theory already tried.

The [split-only follow-up](r1_orderreplay_20260917/SPLIT_CERTIFICATE_THEORY.md)
now proves capped forward replay and maintains certificate lower bounds by
processing only losses from split old ties. Uncertain vertices are refreshed;
only exact failures drive repair. Existing reverse occurrences suffice, with
O((n+I) log(n+1)) comparison/member work per complete round, plus wide
arithmetic. Total repeated-round cost and useful coverage remain unresolved.
No new implementation or measurements; the benefit gate has not passed.

Earlier assessment (2026-09-18): the
[round-barrier proof](r1_orderreplay_20260917/ROUND_BARRIER_THEORY.md) constructs
an arbitrary-s family that survives the existing caps and size expiration.
A legitimate previous peeling order still needs linearly many repairs,
making the specified global scans quadratic; initial-degree clipping does
not remove this obstruction. Reject that full-round realization before
code, while retaining the local lemmas and the unresolved incremental route.
A second paired construction has identical ALL-smaller-size core profiles
and largest-clique labels but different target cores. Numeric core outputs
alone cannot give an exact recurrence; this does not exclude DP retaining
structural/counting state. No new experiment or general impossibility claim.

Earlier continuation (2026-09-18):
[forward-count transfers](r1_orderreplay_20260917/OWNER_TRANSFER_THEORY.md)
derive exact adjacent-swap and bulk-left-move identities. Only cliques whose
first vertex changes transfer their forward contributions; path formulas
avoid explicit enumeration. This is not a two-target update rule for raw
support or core certificates. Size recurrence, order maintenance, discovery,
prefix labels and certificate effects are separately charged. Full-solver
benefit remains unresolved; no code, tests, timings or Base promotion.

Earlier continuation (2026-09-18):
[partial certificate pruning](r1_orderreplay_20260917/PARTIAL_CERTIFICATE_THEORY.md)
holds the candidate vector fixed and removes failed witness vertices. The
greatest feasible survivor subset has exact candidate values, even when
the full vector fails. Each vertex is a source once, but dirty paths can
be processed repeatedly. Known vertices must still supply support events
to unresolved vertices in the restored graph. This specializes existing
generalized-core pruning, not a new basic pruning theorem. Useful total
savings remain unresolved; no implementation, tests or measurements.

Earlier continuation (2026-09-18):
[cross-size extension DP](r1_extensiondp_20260918/THEORY.md) derives the
exact residual-dependent weighted transition from s-cliques to (s+1)-clique
degrees. Deletions must remove faces AND update extensions of surviving
faces. Freezing weights gives wrong final cores on an arbitrary-s family;
initially equal family weights can also diverge after external deletions.
Reject explicit per-small-clique storage before implementation: it expands
the representation and establishes no saving against compressed path
losses. The identity itself was already used in the region-bound work.
A compressed dynamic transition remains unresolved, not ruled impossible.
No new code, test, timing, paper change or promotion.

Earlier continuation (2026-09-18):
[stable extension blocks](r1_extensiondp_20260918/BLOCK_THEORY.md) groups
possible new smallest vertices by their exact neighborhood profiles. It
proves unique one-step coverage and closed-form batch losses with no
deletion-time group splitting. A restricted merge of existing paths also
reduces memberships/source work and never increases nonzero degree writes
under a matched trace. Discovery, coverage and full-solver benefit remain
unresolved; simple repeated extension can grow the representation even
on complete graphs unless additional compression is applied. Static reuse
is not cross-s residual DP; general factorization is known. No code, test,
probe, benchmark, paper edit or Base promotion.

Earlier continuation (2026-09-18):
[factored-loss work bounds](r1_extensiondp_20260918/FACTOR_WORK_THEORY.md)
close the target-read gap left by nonzero-write dominance. A static merged
choice list can revisit dead choices; one lazily cleaned scratch list,
without inverse positions, gives cumulative target-read dominance under
a matched stream trace. Its proof charges each discarded ID to its own
original-path deletion. Two unit counts suffice for the exact three-role
loss formula. Recognition, scratch writes, arithmetic, conversion peak
and active-replay integration remain charged; coverage and full benefit
are unresolved. No C++, test, probe, benchmark or paper edit. This is a
same-size factoring bound, not the requested all-size residual DP.

Pre-implementation assessment (2026-09-18):
[construction-time branch factoring](r1_extensiondp_20260918/TERMINAL_FACTOR_THEORY.md)
avoids post-hoc path matching for locally isolated candidate branches.
The old independent-tail test/direct completion is not new. The proposed
replacement retains one zero-or-one group for a full tail, or exactly one
isolated choice plus the unchanged continuing search for a partial group.
The partial case works with an arbitrary remaining candidate graph but
charges a minimum-degree list in the required degree pass. Residual losses,
size bounds, flat-membership savings and compact joint q_s,z_s state are
derived; a fixed-size lazy list cannot be shared blindly across sizes.
The nonempty-child extension also needs its residual clique structure:
an arbitrary-s construction refutes using only the two counters, child
size and a prior coefficient. Richer shared child state is not ruled out,
but overlaps the existing shared-DAG route and keeps its costs.
Old logs lacked group coverage, leaving total benefit unresolved at that
stage. The authorized [implementation](r1_terminal_20260918/RESULTS.md)
now measures this coverage and shows conditional gains; it does not
implement the joint cross-size counters. This links to shared peeling,
not a new ascending-s residual DP. Production and manuscript remain unchanged.

| Priority | Result | Stage | Status |
|---|---|---|---|
| Conditional measured result, not promoted | Construction-time factoring implementation | Build, all-size stream and compiled replay | Release/sanitizers pass; uniform old-active/new replay ratios 1.134 paired / 1.171 fresh, but small gains reverse and best-old per-graph comparison is almost neutral. Clear Stanford work/memory savings, no production promotion |
| Theory with scoped implementation | Construction-time zero-degree branch factoring | Index output and shared-peeling state | Full/partial local groups avoid post-hoc matching; fixed-size losses implemented with conditional measured gains. Joint cross-size counter scheme remains theoretical |
| Theory only, full benefit unresolved | Amortized factored-loss scan bound | Same-size source/state factoring | One lazy choice list gives cumulative target-read dominance against matched stream; full cost/coverage unknown, no implementation |
| STARRED | Bulk universal pivots | Index construction | Accepted default for the three scoped callback builders; use production V2 evidence |
| STARRED | Root-suffix scans | Index construction | Useful experimental option; NOT promoted to default |
| Retained | Certified clique-block elimination | Before index construction | Earlier accepted research method; NOT the production construction Base |
| Experimental | Whole-class weighted R1 quotient | Before construction and during peeling | Validated structural reduction; conditional wins, no stable real-input win over strongest controls |
| Retained implementation | Hash grouping and bitset search | Grouping and index construction | Same R1 theory/index as V1; measured phase sum 3.184x faster than V1, 2.573x over original controls; not installed in production |
| Theory first | All-s raw-k query foundations | Cross-parameter values and connectivity | Raw values are not globally monotone; proved fixed-k suffix and transformed-threshold nesting; exact checks pass, no full index yet |
| Retained theory, not promoted | Joint all-s frontiers and ordinary-ceiling completion | Shared peeling across sizes | Exact all-size prototypes; fewer events, but final fixed/batch peeling ratio 0.678 geometric mean; no unconditional once-per-vertex result |
| Retained theory, not promoted | Clique-bound-aware all-s completion | Shared peeling across sizes | One completion per vertex on balanced complete multipartite inputs; only 63 extra certified real vertices; bound and phase variants lose overall to independent peeling |
| Retained theory, not promoted | Adjacent-size floor completion | Shared peeling across sizes | Proof-first extension removes 195,352 extra real events (13.03%); same-kernel time ratio 1.163, but fixed/floor 0.962 overall; no universal n-event result |
| Retained theory, not promoted | Sharp shadow inference and permanent incidence pruning | Shared peeling across sizes | Proof-first variants pass independent oracles; shadow skips only 71 real count queries. Pruning removes 33.19% of count reads, but neither it nor a separate scan-order follow-up wins overall |
| Parked by user | Cross-vertex all-size batches | Shared peeling across sizes | Not adopted; preserve evidence and stop this variant. Same-batch target reads fall only 1.26%; snapshot scheduling adds 1.81% vertex events. No overall win |
| Retained theory, not promoted | Shared residual sets with size splitting | Shared peeling across sizes | Two proof-first variants pass independent oracles. V1 fixed/shared ratio 0.843; complete-minimum refinement 0.676. More working state, no overall win; preserve the 14-vertex global-order conflict |
| Retained theory, not promoted | Bottom-up search reuse and positive-label transfer | Small-s to large-s computation | Two proof-first variants. Frontier reuse is 3.158x repeated construction but loses to shared-index controls. Positive transfer cuts local corrections 28.87%, yet its local kernel is 3.598x independent-peel time. No overall win |
| Retained, not promoted | Triangle-component shared peeling | Shared peeling across sizes | Exact independent components, but copies/state erase the event savings; fixed/component ratio 0.680 in the parent sweep. This failed variant does not rule out the whole road |
| Promising research candidate | Bottom-up capped peeling and sorted-order reuse | Small-s to large-s computation | Paired solver fixed/DP 1.307 wall, 1.304 CPU; ordinary-only stream/DP 1.062 wall, 1.076 CPU. Heap comparisons fall 42.72% beyond ordinary bounds; raw count reads do not fall. Not production, not a complete count-state DP |
| Retained, not promoted | Count reuse along inherited thresholds | Small-s to large-s computation | Two proof-first variants pass independent oracles. Counter reuse removes member rescans, but suffix coverage is only 2.65%. Paired stream/reuse is 0.827; guarded refinement 0.823. No overall improvement over the preceding stream |
| Conditional, not promoted | Witnessed output-curve scheduling | Small-s to large-s computation | Known absorption theorem implemented with support-preserving events. Determines 61.11% of positive output pairs; HepPh improves, but paired/fresh-process overall signs differ. Stream/pruned variant is 1.126 paired versus 0.953 fresh; no stable overall or memory win |
| Theory/certificate prototype, not promoted | Cross-(s,k) component intervals and two-anchor inference | Before construction / between solved sizes | Density plus boundary certifies entire query intervals without clique search, but adds only 1,653 covered vertex-size pairs on five real graphs. Two-anchor component sandwich proved and independently checked; not yet a faster all-size solver |
| Retained theory, not promoted | Two-anchor bounded peeling and linear path certification | Cross-size core-value computation | Two proof-first schedules pass independent oracles and sanitizers. Paths can be omitted with unknown members when all clique minima are fixed. Main linear-prune stream/candidate ratio is 0.824 paired and 0.700 fresh; extra checks/state erase savings. Full shared index, production and manuscript unchanged |
| Retained, not promoted | Lower-curve DP and initial-degree floor stream | Cross-size core-value computation | No path-certification pass. Main floor-stream ratio is 0.951 paired / 0.891 fresh. Additional equality test fixes only1,418 positive pairs and saves exactly that many heap pops, no support updates. Proof explains the benefit ceiling; unchanged index, more scratch |
| Retained, not promoted | Initial count-curve aggregation | Cross-size initialization reuse | Equal binomial curves reduce nonzero additions by 84.69%, but not subsequent peeling. Packed follow-up gives stream/candidate ratios 0.838 paired / 0.902 fresh and higher RSS. Both versions pass independent oracles; unchanged index and Base |
| Parked by user | Complete-join all-size core DP | Before construction, across all sizes | Conditional theorem retained as a special case. Only 15 extra real vertices/pairs found; overall ratios 0.958 paired / 0.932 fresh. Do not pursue the proposed private-region detector or promote the structured-family result as a general contribution |
| Retained theory, not promoted | One-sided contributions and permanent target omission | Cross-size bounds and peeling | General-graph omission persists at larger valid s; prefix boundaries move only once. Degree writes fall38.70% versus matched floor control. Beyond basic size expiration, inference saves9.45% more source reads and41.34% boundary probes. Latest total ratios1.032 paired/1.053 fresh are small/noisy, with regressions and more RSS; no Base promotion |
| Retained theory, not promoted | Mandatory-witness source omission | Cross-size source notifications | Complete-batch and strict witness-replacement proofs permit permanent source omission. Source reads fall31.25%, but live pivot decrements fall only0.495% and actual support writes do not fall. Frozen-monotone/new ratios0.950 paired/1.035 fresh; extra4P bytes. No stable overall win; this rule cannot remove nonzero batch losses |
| Completed measurement, not promoted | Full ordinary-core size range | Existing stream versus direct peeling | Every s through d+1:44/239/114/72/7 on five real graphs, exact64/256/128/64/64-bit counts. Twelve paired trials give direct/stream1.167 and empty-tail-stop/stream1.119 geometric means, with noise and regressions. Not an end-to-end ratio or a new algorithm |
| Conditional research option, not promoted | Certified order replay and active-size expiration | Cross-size residual counting | General whole-vector certificate eliminates complete layer peels. Active child gives stream/new1.175 and fixed-stop/new1.315; actual positive decrements fall3.32% overall. Stanford/Amazon regress, scratch increases; preserve both implementations and all evidence |
| Theory retained; full-round prototype rejected | Local order repair, capped prefix counts and split-only certificates | Failed cross-size replay layers | Correct local lemmas retained. New padded family survives existing caps/filtering and a legitimate inherited order yet gives quadratic full-round scans. Fully incremental benefit unresolved. Identical all-smaller-size profiles do not determine the next core row; richer-state DP is not excluded. No new implementation, measurement or novelty claim |
| Theory retained; benefit unresolved | First-vertex transfers for local order changes | Forward counting inside failed replay repair | Adjacent swaps and bulk left moves have exact compressed transfers. They do not remove core-certificate maintenance; dynamic order adds state and work. No proved total gain over stream, implementation or new measurements |
| Theory retained; benefit unresolved | Fixed-label partial certificate pruning | Salvaging a subset after full-vector failure | Greatest feasible survivors have exact values; each verification vertex is removed once, but paths can be revisited. Remainder must restore support and schedule known vertices correctly. Known generalized-core principle; no proved total gain or implementation |
| Theory assessment; explicit prototype rejected | Cross-size extension DP | Small-s to large-s residual counting | Dynamic weighted recurrence proved; frozen weights give wrong cores. Explicit face state and reverse extensions expand the index and recreate target work. Compressed transitions need external-dependency handling and a total benefit proof; no implementation or measurements |
| Theory retained; full benefit unresolved | Stable extension blocks and restricted path factoring | Compressed one-step extension / shared state inside a size | Exact ownership and static profile groups avoid face expansion and deletion-time splitting. Qualifying existing-path groups reduce occurrences and repeated nonzero writes. Construction growth and representative total gain remain open; no implementation or measurements |
| Parked | Child-degree reuse | Index construction | Correct but no useful overall measured gain |

STARRED means the user wants the result preserved and prioritized. It does
not mean every dataset improves, literature novelty is established, or the
code is enabled in production. The two starred results do not change peeling.

## The Two Starred Results

### 1. Bulk Universal Pivots

The required pivot-degree scan identifies candidates adjacent to all other
candidates. Contract their consecutive single-branch search steps together;
do not recount the remaining degrees and repeat intersections at each step.
The ordered output is preserved, with no new persistent index fields.

- Initial standalone theory and measurements: [r1_bulkpivot_20260915/RESULTS.md](r1_bulkpivot_20260915/RESULTS.md).
- Production integration: [r1_bulkbase_20260915/RESULTS.md](r1_bulkbase_20260915/RESULTS.md).
- Actual default: `src/SDCT_Augmented.inl`, `PIVOTER_BULK_PIVOTS=1` or unset.
  `=0` selects the old recursion. Scope is NoTree, FlatEvents and interleaved
  callbacks, not every parallel/universal/retained-tree builder.
- Use `benchmark_v2.json`, [V2 table](r1_bulkbase_20260915/benchmark_v2-table.md),
  [audit_v2.json](r1_bulkbase_20260915/audit_v2.json), and
  [commands](r1_bulkbase_20260915/REPRODUCE.md). V1 is archived, not authoritative.
- Production construction old/new ratios across s=5,8 and two layouts:
  GrQc 1.04-2.18; HepPh 9.92-13.34; DBLP 1.04-1.68;
  Stanford 1.10-1.21; Amazon0302 0.96-2.16.
- Neither the standalone nor production comparison is a CND speedup claim.
  Final index size is unchanged; peak RSS is not guaranteed to decrease.

### 2. Root-Suffix Scans

Root preparation previously scanned a candidate's whole adjacency row and
immediately rejected all earlier neighbors. Sorted rows let it find the
boundary and read only the later-neighbor suffix. It removes exactly
`sum_u a(u)^2` inner-loop reads, where `a(u)` is u's earlier-neighbor count.
With later degree at most degeneracy d, retained scans are at most m*d;
binary searches and one sortedness check remain additional work.

- [Theory](r1_forwardroot_20260916/THEORY.md),
  [implementation](r1_forwardroot_20260916/include/root_suffix.inl),
  [report](r1_forwardroot_20260916/RESULTS.md),
  [commands](r1_forwardroot_20260916/REPRODUCE.md).
- Same-binary controls: `PIVOTER_ROOT_SCAN=base|scan|suffix` with bulk pivots
  enabled. `scan` isolates local-helper placement from prefix skipping.
  Unsorted rows fall back to the original helper without changing the graph.
- [Full table](r1_forwardroot_20260916/benchmark-table.md),
  [counted work](r1_forwardroot_20260916/work-table.md),
  [audit](r1_forwardroot_20260916/audit.json). Raw timings are in `benchmark.json`
  and `benchmark-logs/`; correctness logs and manifests are in the same directory.
- Construction Base/suffix ratios across s=5,8 and two layouts:
  GrQc 0.926-1.040; HepPh 0.799-1.083; DBLP 0.830-0.989;
  Stanford 6.582-8.278; Amazon0302 0.887-1.098.
- Stanford construction-plus-peel improves 5.661-7.198x. This excludes
  hierarchy, input and ordering. The peeler itself is unchanged.
- Stanford root reads fall from 7,831,982,712 to 27,345,326, plus 16,506,688
  boundary comparisons. This read ratio is NOT the runtime ratio.
- Overall construction geometric mean is 1.446x; excluding Stanford it is
  0.960x. Peak-RSS ratios range from 0.932 to 1.108. Keep all regressions.
- Release and sanitizer suites pass; ordered traces cover 34,777 graphs /
  208,377 configurations / three interfaces. All 180 timed runs match full
  per-vertex output. This is still an isolated experimental overlay.

## Complete Inventory

Rows below cover every top-level research directory. Related experiments
outside `research/` are linked separately at the end. The linked report owns
the exact configuration, correctness scope, timing boundary and final dataset.
Historical ratios must NOT be pooled across different controls or versions.

### Output Certificates, Index Avoidance And R1 Reductions

| Directory / report | Idea and stage | Outcome / reuse decision |
|---|---|---|
| [r1_theory_20260905](r1_theory_20260905/THEORY_REVIEW.md) | Initial theory: partial-output certificates; one-time hierarchy activation; loss-triggered peeling | Sound tested rules. Notification prototype is much slower; certificates need useful evidence early. Contains separate C++ iterations listed below |
| [r1_indexfree_20260914](r1_indexfree_20260914/RESULTS.md) | Recount through graph search without retaining a full index; also minimum-region-size pruning | No nonempty real-case total-runtime win over the strongest tested control. Keep the sound MCE size guard; index-free solver remains experimental |
| [r1_simplicial_20260914](r1_simplicial_20260914/RESULTS.md) | Remove certified private clique blocks before construction; preserve boundary lower bounds and hierarchy stars | Measured benefits on GrQc/DBLP, not universal. Retained earlier research option, not production integration |
| [r1_fused_20260914](r1_fused_20260914/RESULTS.md) | Suspend clique search and use certified removals to shrink unfinished search | Full output checked, but no stable improvement over strong preceding controls; parked |
| [r1_multipartite_20260914](r1_multipartite_20260914/RESULTS.md) | Solve private complete-multipartite regions with coefficient formulas before construction | Proven structural separation on a constructed family; few additional regions on real data. Do not repeat broad timing without new applicability evidence |
| [r1_weightedquotient_20260916](r1_weightedquotient_20260916/RESULTS.md) | Collapse true-twin classes before weighted construction and whole-class R1 peeling | Correct and complete. Conditional synthetic wins; real matched-phase geometric mean 0.772 against strongest production controls / 1.036 against singleton ablation. Not promoted |
| [r1_quotientimpl_20260916](r1_quotientimpl_20260916/RESULTS.md) | Same weighted R1 theory, hash grouping and root-local bitset search | Complete: 3.184x V1 / 2.573x original-production phase-sum geometric means. Final index bytes unchanged; RSS can rise. Fast singleton/grouped ratio is only 1.013 overall: separate implementation and theory |

Important sub-attempts within `r1_theory_20260905`:

| Entry | What changed | Outcome |
|---|---|---|
| [PARTIAL_INDEX_RESEARCH.md](r1_theory_20260905/PARTIAL_INDEX_RESEARCH.md) | Minimum-core certificate for omitted cliques, bounded peeling and witness connectivity | Python prototype: conditional benefit, not production evidence |
| [cpp/RESULTS.md](r1_theory_20260905/cpp/RESULTS.md) | C++ partial construction, unrestricted and budgeted witness passes | HepPh improvement over matched full control; costly witness generation causes losses elsewhere |
| [cpp-online/RESULTS.md](r1_theory_20260905/cpp-online/RESULTS.md) | Obtain certificates during search rather than a separate greedy pass | Weaker early evidence; not a replacement for the strongest offline control |
| [cpp-replay/RESULTS.md](r1_theory_20260905/cpp-replay/RESULTS.md) | Replay witnesses, reuse true-twin pivot degrees, simplify certificate tests | Equivalent linear minimum test; no general runtime improvement. Twin grouping here preserves the original labeled index and peeler |

### Peeling And Shared Computation

| Directory / report | Idea and stage | Outcome / reuse decision |
|---|---|---|
| [r1_dominance_20260913](r1_dominance_20260913/RESULTS.md) | Persistent degree order permits deferring updates while a smaller-degree witness survives | Fewer decrements but detection/restoration/target-list costs erase net gains; parked |
| [r1_circuit_20260913](r1_circuit_20260913/RESULTS.md) | Share loss propagation through an arithmetic/set circuit | Synthetic work saving; no nonempty real-case win over the fastest tested control including construction; parked |
| [r1_jointpeel_20260916](r1_jointpeel_20260916/RESULTS.md) | One size frontier per vertex; finish a whole prefix at its ordinary-core ceiling | Proven/tested joint process, not n events in general. Final all-s outputs agree; events fall 2,329,788 to 1,499,270, but no net runtime win. Separate initial, ceiling V1 and final evidence |
| [r1_envelopecore_20260916](r1_envelopecore_20260916/RESULTS.md) | Replace binomial ceilings by colored-shadow/Turan bounds; test threshold phases instead of cross-size rank search | Wider conditional theorem, exhaustive checks pass. Real events fall by only 441 beyond the previous joint rule. Fixed/integer time ratio 0.759; separate noisy phase trial also loses overall. Do not promote |
| [r1_floorbatch_20260916](r1_floorbatch_20260916/RESULTS.md) | Complete consecutive sizes at valid floors; separately test initial-minimum seeding | Detailed proofs precede C++; 34,294-graph release/sanitizer suites pass. Real events 1,499,270 to 1,303,918; entry-count reads unchanged. Improved joint control, but not an overall independent-control win. Retain, no production promotion |
| [r1_countskip_20260916](r1_countskip_20260916/RESULTS.md) | Skip lower-size counts using sharp shadow certificates; remove permanently irrelevant incidences | Two proof gates precede the original and scan-order kernels; both pass 34,289-case Release/sanitizer suites. Only 71 real queries avoided by shadow. Pruning cuts count/source reads but not target work; both matched sweeps show no overall strong-control win. Retain, not promoted |
| [r1_wavebatch_20260917](r1_wavebatch_20260917/RESULTS.md) | Share target updates across justified vertex/size deletion batches | PARKED BY USER; not adopted. Preserve proof, code, 34,306-case suites, 168 timings and nine-vertex example. Fixed/wave ratio 0.719; same-batch sequential/wave 0.785. Do not resume this variant without an explicit request |
| [r1_sizesplit_20260917](r1_sizesplit_20260917/RESULTS.md) | Shared physical residuals; a refinement maintains the complete common-minimum set with monotone eligibility masks | Correct but not promoted. V1 saves 16.35% of events but only 1.31% of target reads; the refinement saves 0.46% more events yet adds 7.58% target reads. Preserve both proofs, 120 V1 timings, 90 authoritative refinement timings, and separate exploratory data |
| [r1_bottomup_20260917](r1_bottomup_20260917/RESULTS.md) | Ascend s with saved search states; separately transfer previous positive core labels to a local solver | Both Release/ASan and full real outputs agree. Preserve 210 frontier and 84 label timings. Extra exclusion saves almost no new search beyond memoization; label transfer helps its own kernel 1.152x but not strong peeling. Not promoted |
| [r1_tworoads_20260917](r1_tworoads_20260917/RESULTS.md) | Compare component-local shared peeling with bottom-up integer-capped peeling; [stream refinement](r1_tworoads_20260917/stream/RESULTS.md) reuses the previous core order | Parent and stream each have 270 timings plus nine full verifications and separate Release/ASan suites. Additional paired check has 180 measured calls and 30 verified warmups. Stream is promising, but the DP-only increment is modest, counts/index size remain unchanged, and no production promotion is made |
| [r1_countreuse_20260917](r1_countreuse_20260917/RESULTS.md) | Reuse threshold membership counters to certify an inherited high-core suffix; [guarded refinement](r1_countreuse_20260917/guarded/RESULTS.md) avoids group insertion on first-member failure | Both Release/ASan suites pass 34,296 graphs. Preserve separate 162+180 and 189+245 broad/paired timings. Guard cuts extra reads 11.497M to 5.223M but not certificate coverage; remains about 1.215x preceding-stream time. No new index or production edit |
| [r1_taildp_20260917](r1_taildp_20260917/RESULTS.md) | Schedule known all-size output curves from clique witnesses; ordinary-only, cross-s and all-known-path controls | Existing absorption theorem, not a new discovery. Release/ASan each pass 36,344 graphs and large-integer boundaries. Preserve 189 broad timings, nine verifications and 245 paired calls. HepPh benefits, Stanford regresses, aggregate signs differ by measurement design. Conditional only; unchanged index and production |
| [r1_skregions_20260917](r1_skregions_20260917/RESULTS.md) | Stronger relations between size, threshold and entire components; [two-anchor consequence](r1_skregions_20260917/SANDWICH_THEORY.md) | Latest scalar shadow bound dominates older ones. Density/boundary C++ passes 33,991 graphs per build and real full-set checks; 72 overhead timings, no all-s speedup. Five real S10 inputs gain 213 non-clique query intervals / 1,653 covered pairs, no extra exact positive core values. Separate two-anchor checker passes 1,261 graphs; actual two-anchor solver/costs are now recorded in the following attempt |
| [r1_anchorpeel_20260917](r1_anchorpeel_20260917/RESULTS.md) | Compute both size anchors, omit paths with fixed clique minima, then bounded-peel the remainder; balanced and consecutive-anchor schedules | Release/ASan each pass 36,337 configurations per version, 37,173 symbolic path tests and 35,738 inverse tests. 771 timed calls; all paired full arrays agree. Parent fixes 786,773 positive pairs but adds 79.7M certificate reads. Linear child reduces order state but loses overall: solver ratios 0.824 paired/0.700 fresh, compute totals 0.948/0.881. Clique synthetic wins are conditional; no Base promotion or hierarchy claim |
| [r1_floorstream_20260917](r1_floorstream_20260917/RESULTS.md) | Reuse high-anchor lower curves by DP; test full initial degree against lower and schedule fixed events using a static stream; zero-source guard as separate control | Release/ASan each pass36,345 configurations; a separate [benefit-ceiling proof](r1_floorstream_20260917/postmortem/THEORY.md) and checker cover34,277 graphs each. 666 timings, 101 main processes plus2 diagnostics. Floor-stream scans5.624M extra positions to save1,418 heap pops; no support decreases saved beyond floor-DP. Main ratios0.951/0.891; zero-guard overall signs differ. No production or paper edits |
| [r1_countcurve_20260917](r1_countcurve_20260917/RESULTS.md) | Sum identical initial count curves once per vertex/type; [packed follow-up](r1_countcurve_20260917/packed/THEORY.md) isolates metadata lookup costs | Each version passes 36,337 Release/ASan configurations; 973 timings and 93 processes. Initial nonzero additions fall 62.063M to 9.502M, but source/target/queue work is identical. Parent ratio 0.835; packed 0.838 paired / 0.902 fresh. Temporary 4I code bytes can increase RSS. No new (s,k) theorem, permanent index reduction or Base promotion |
| [r1_joindp_20260917](r1_joindp_20260917/RESULTS.md) | Complete joins of internally uniform factors yield final core curves; [proof](r1_joindp_20260917/THEORY.md) includes private-boundary extension and a no-gain neighbor-refinement theorem | PARKED BY USER as too special-case. Release/ASan each pass 34,332 configurations. Whole-component recognition adds only 15 real vertices; ratios0.958/0.932, RSS increases. Preserve synthetic savings and negative real evidence, but do not continue the proposed private-region detector without a new request |
| [r1_targetprune_20260917](r1_targetprune_20260917/RESULTS.md) | One-sided clique omission from valid bounds; [permanent cross-size consequence](r1_targetprune_20260917/monotone/THEORY.md) permits shrinking target prefixes and safe source compaction; [range control](r1_targetprune_20260917/range_control/TABLES.md) isolates basic size expiration | Four proof-first variants,34,300 Release/ASan configurations each,540 repeated timing calls plus10 work-only calls. Omission saves degree writes; cross-size inference saves6.409M source reads beyond size expiration. Last aggregate ratios1.032/1.053 are not a large stable win; preserve earlier negative variants. Index bytes unchanged, RSS higher, production/paper unchanged |
| [r1_sourceskip_20260917](r1_sourceskip_20260917/RESULTS.md) | Strengthen path omission using a retained mandatory witness; [proof](r1_sourceskip_20260917/THEORY.md) permits permanent source omission even when the whole path is still needed | Complete: Release/ASan and independent schedules pass;140 repeated timed calls and25 work-only calls. Source reads61.416M->42.223M, but support writes unchanged. Frozen-monotone/new ratios0.950 paired/1.035 fresh, extra4P bytes, no Base promotion. [Benefit boundary](r1_sourceskip_20260917/BENEFIT_BOUNDARY.md) proves that equal batch coefficients leave actual nonzero losses unchanged; stop tuning this as a count-reuse method |
| [r1_fullrange_20260917](r1_fullrange_20260917/RESULTS.md) | Measure every size through maximum ordinary core + 1, using identical exact integer widths and a direct empty-tail-stop control | Complete. Release/ASan each34,048 cases; large analytic cliques validate counts beyond uint64. All180 timed full outputs agree. Direct/stream per-graph ratios1.018/1.311/1.333/0.953/1.274, geometric mean1.167; stronger tail-stop control1.119. Original sources unchanged, no Base/paper promotion |
| [r1_orderreplay_20260917](r1_orderreplay_20260917/RESULTS.md) | Reuse per-path residual membership statistics along an order; certify whole output vectors before omitting peeling; [active refinement](r1_orderreplay_20260917/ACTIVE_THEORY.md) expires irrelevant sizes; [local repair](r1_orderreplay_20260917/LOCAL_REPAIR_THEORY.md), [split-only certificates](r1_orderreplay_20260917/SPLIT_CERTIFICATE_THEORY.md), [round barrier](r1_orderreplay_20260917/ROUND_BARRIER_THEORY.md) and [forward transfers](r1_orderreplay_20260917/OWNER_TRANSFER_THEORY.md) are theory only | Proof-first parent/child Release+ASan each34,048 configurations/150,056 arbitrary-order size checks. All240 repeated timings and5 work-only calls agree with full references. Stream/new1.175, fixed-stop/new1.315; decrements106.015M->102.499M. More memory, clear Stanford/Amazon regressions; conditional only, no Base/paper promotion. Capped full-round repair rejected on the padded family; useful incremental repair unresolved. Exact local forward transfers do not remove certificate costs. Numeric previous outputs alone are insufficient for exact DP. No new measurements |
| [r1_terminal_20260918](r1_terminal_20260918/RESULTS.md) | Construction-time full/partial isolated-branch factoring, packed plain ablation, and factored compiled replay; separate proofs precede implementation | Complete V2 Release/ASan each34,054 configurations,5,160,192 core cells,843,864 structural/residual/order checks. All320 paired full outputs and120 fresh hashes agree. Uniform old-active/new ratios1.134 paired/1.171 fresh; near-neutral paired comparison with best old per-graph controls. Stanford memberships fall44.25%, positive writes30.52%, replay RSS25.94%; HepPh/DBLP RSS slightly increases. Conditional only; preserve V1 evidence, no Base/paper promotion |

### Query Index Design (r=1, All Sizes)

| Directory / report | Idea and stage | Outcome / reuse decision |
|---|---|---|
| [r1_skyline_index_20260918](r1_skyline_index_20260918/RESULTS_INDEX.md) | All-size community/value index for r=1 adapted from the SIGMOD 2026 SGL bi-component index: dominance through iterated Kruskal-Katona shadows, per-size canonical merge trees with one cross-size pointer per node, one entry per class per skyline size, certified tail = one entry, twin quotient; [theory](r1_skyline_index_20260918/THEORY.md), [counts](r1_skyline_index_20260918/RESULTS.md), then both designs implemented in memory and measured | Proofs plus counts plus implementation. Stage 1: Release/ASan selftests match brute-force nuclei on 34,075 graphs; F2 2,994,580 checks, 0 violations on five real graphs. Stage 2: both designs match brute force (605,476 community queries, 3.65M membership, 1.44M values) and each other on every timed query. Measured bytes (class layout): skyline smaller by 1.11x GrQc, 2.28x HepPh, 1.28x dblp, 1.17x Stanford, 0.99x amazon (with values 1.09/1.84/1.20/1.11/1.00); per-vertex layout computed from the same counts 1.14x-3.88x. Latency: skyline 2.1-2.9x slower on community listing (per output vertex 3.8-5.8 ns vs 1.6-2.3 ns), 1.0-4.4x slower on membership, equal on values. A space/time trade, never both; S trees with DFS intervals remain the fastest exact representation. No production, paper or promotion. 2026-09-19 [chains](r1_skyline_index_20260918/CHAINS.md): hierarchy-equivalence classes (same own node at every size) are a partition refined by twins and coarser by 6.4x-23.6x on the five inputs (dblp 317,080 vertices -> 13,459 chains); every nucleus is a union of chains (proved); projected S trees over chains 2.4x-4.7x smaller than per-vertex S trees, 3.8x-8.1x with aligned labels, at S-tree query speed; counts only, chain index not yet built. MEASURED 2026-09-19 ([results](r1_skyline_index_20260918/RESULTS_CHAINS.md)): S trees over chains vs S trees over twins, with values: GrQc 2.48x, HepPh 3.88x, dblp 3.80x, Stanford 4.53x, amazon 2.97x smaller (aligned-label projection 3.9x-8.3x; vs naive dense 2.7x-18x), community listing 1.7x-8.6x FASTER (0.26-0.40 ns per output vertex), membership up to 3.1x faster, brute-force checked (1,210,952 queries, both modes, Release+ASan); skyline over chains is dominated. Both smaller and faster: the SGL effect, obtained from a partition instead of mu-fold units. FINAL MODULE 2026-09-19 ([RESULTS_FINAL.md](r1_skyline_index_20260918/RESULTS_FINAL.md), `chain_index.hpp`, CHAINS.md Section 9 Lemmas C4-C6): aligned labels + lexicographic chain ranks + per-size run arrays with node entry points; brute-force checked in build, compact and loaded forms (1,816,701 community queries, 10.96M membership checks, Release+ASan); bytes vs per-vertex S trees 3.98x GrQc, 5.70x HepPh, 7.81x dblp, 6.08x Stanford, 3.73x amazon (2.9x-5.0x if the 4n-byte label permutation is charged); community located in O(1) (3-40 ns), range answers 26-660 ns, explicit ids 0.6x-1.65x memcpy speed (faster on the three large graphs, 0.08-0.19 ns per vertex), one range per (2,k)-community by construction; membership 5-18 ns, values 9-18 ns; build 0.02-4.8 s one thread including the all-size peel. Compact run form costs 2-10 percent more bytes than chain-id arrays and answers 5x-53x faster |

### Alternative Counting And Construction Representations

| Directory / report | Idea | Outcome / reuse decision |
|---|---|---|
| [r1_searchreuse_20260914](r1_searchreuse_20260914/RESULTS.md) | Reuse equivalent candidate searches; derive prefix/suffix count identities | Stanford has measurable reuse, but probe does not establish full-pipeline improvement |
| [r1_shareddag_20260914](r1_shareddag_20260914/RESULTS.md) | Shared candidate-state DAG constructor | Stanford construction improves locally, with substantial memory overhead and losses on low-reuse graphs; not promoted |
| [r1_conflict_20260915](r1_conflict_20260915/RESULTS.md) | Count using sparse nonedge constraints | Correctness gate passes; target states cover at most 2.65% of measured pivot queries. Joint report is in the next row |
| [r1_overlap_20260915](r1_overlap_20260915/RESULTS.md) | Signed overlap counting instead of disjoint recursive expansion | Different representation, but no established speed/memory win in this first version |
| [r1_overlapcert_20260915](r1_overlapcert_20260915/RESULTS.md) | Reuse structural certificates for signed overlap | Correct, but extra recognition work and representation cost remain; V1 not promoted |
| [r1_overlapcert_v2_20260915](r1_overlapcert_v2_20260915/RESULTS.md) | Reuse required pivot scans for those certificates | Conditional limited full-pipeline wins; larger signed representation fails the user's desired size tradeoff. Not production |

### Construction Without Changing The Final Index

| Directory / report | Idea | Outcome / reuse decision |
|---|---|---|
| [r1_bulkpivot_20260915](r1_bulkpivot_20260915/RESULTS.md) | STARRED: batch consecutive universal pivots | Standalone origin; use the next row for current production status |
| [r1_bulkbase_20260915](r1_bulkbase_20260915/RESULTS.md) | Port universal-pivot batching into production callback builders | Accepted default; authoritative source/evidence is V2, not archived V1 |
| [r1_emptychildren_20260915](r1_emptychildren_20260915/RESULTS.md) | Finish independent candidate sets without child recursion | Correct with constructed-family benefits; mixed real performance, disabled by default |
| [r1_childcounts_20260916](r1_childcounts_20260916/RESULTS.md) | Preserve degrees already computed during child preparation | Removes actual reads, adds O(d) scratch; build geometric mean 1.0019x, isolated control 0.9432x. Not promoted |
| [r1_forwardroot_20260916](r1_forwardroot_20260916/RESULTS.md) | STARRED: skip rejected earlier-neighbor prefixes | Large Stanford win, regressions elsewhere; experimental, not default |

### Correctness Audits And Repairs

| Directory / report | Scope | Outcome |
|---|---|---|
| [code_audit_20260914](code_audit_20260914/REVIEW.md) | Loader, counting limits, validation, fixed-s pruning, hierarchy and storage audit | Findings and reproductions, not a performance variant; follow repair status below |
| [code_fixes_20260914](code_fixes_20260914/FIXES.md) | Implement confirmed repairs and focused guards | Integrated fixes and full validation; [slowdown diagnosis](code_fixes_20260914/SLOWDOWN_DIAGNOSIS.md) qualifies initial timing interpretations |

## Evidence And Version Rules

1. Start with `RESULTS.md` (or the linked audit report), then `THEORY.md` and
   `REPRODUCE.md`. Some older attempts put commands directly in `RESULTS.md`.
2. Follow that report's authoritative JSON filenames. `probe`, `initial`,
   `before`, V1 and intermediate validation files are historical evidence,
   not replacements for explicitly designated final measurements.
3. Raw outputs, source/binary/input hashes, all per-vertex comparisons,
   trial ranges and unsuccessful cells remain with their experiment.
4. A later production Base does not retroactively change an old comparison.
   New candidates must rerun against current strong correct controls.
5. Construction, peeling, hierarchy, input/ordering and process wall time
   are different timing scopes. Index bytes and peak process memory are
   different quantities. Never silently substitute one for another.
6. No performance gain in this inventory establishes literature-wide novelty.
   In particular, the two starred construction rules are not r=1-only rules.

## Lessons For The Next R1-Specific Attempt

- Do not repeat a full witness pass or global certificate refresh without
  first showing it can avoid more work than it adds.
- Fewer stored paths are insufficient if discovering them already did all
  the expensive work. Construction avoidance must happen before that search.
- Preserving every clique count is stronger than preserving all core labels.
  Any reduction exploiting this difference must specify the replacement
  input contract and prove it, not run ordinary peeling on incomplete counts.
- Equal core labels alone do not certify hierarchy connectivity.
- Reusing neighbor degrees and preserving output paths has already been tried.
  A new class-based approach must actually reduce construction/peeling objects,
  not just repeat the old true-twin degree cache under another name. Whole-class
  weighted R1 now does this, but closed-neighborhood grouping and coefficients
  erase much of the real-input benefit. Its [next directions](r1_weightedquotient_20260916/RESULTS.md)
  prioritize broader equivalence coverage before another full solver rewrite.
  The [implementation follow-up](r1_quotientimpl_20260916/RESULTS.md) subsequently
  showed large same-theory speedups: do not treat V1's poor implementation as
  a definitive rejection of the reduction. Fast singleton remains essential
  to identify the part attributable to grouping, rather than bitset search.
  The [ordinary-core census](r1_quotientimpl_20260916/core_filter/RESULTS.md)
  now measures exact classes after (s-1)-core preprocessing on all five
  graphs at s=3..20. At s=8, grouped objects fall by 33.60% on DBLP and
  37.18% on HepPh; Stanford saves 10.68%, Amazon is empty. These are
  structural counts, not new runtime or memory measurements. Earlier
  timings did not apply this prefilter; future controls must both apply it.
  The user's [all-s grouping-index idea](r1_quotientimpl_20260916/core_filter/ALL_S_IDEA.md)
  has a no-splitting/whole-class-exit proof and a linear-size merge-forest
  representation. All 34,268 property tests pass. No shared clique-index
  implementation or speed claim yet; grouping history and clique/core
  results across s are different storage contracts.
  The user then clarified the [full query goal](r1_quotientimpl_20260916/cross_s/THEORY.md):
  maximum s given (k,v), kappa_s(v) given (s,v), and component vertices
  given (k,s,v). They requested monotonicity research before implementation.
  Raw kappa_s is not even unimodal, and same-k components can cross.
  However, fixed-k feasibility and components are nested for all s>=s0(k),
  where s0(k)=min{s>=2: k<=binom(2s-1,s)}. This gives a tested maximum-s
  lookup rule. Binomially transformed thresholds nest across all sizes;
  normalized values decrease, and witnessed equality determines a tail.
  The shadow/absorption foundation already appears in local NSI work; do
  not claim it as a new invention here. The two 34,271-graph suites pass,
  including 620,676 maximum-s queries. Archived s=5,8 outputs also pass
  878,344 vertex-pair and 96,734 hierarchy-edge transfer checks. No new
  performance result, all-s production index, or paper edit in this attempt.
- The [one-peeling follow-up](r1_jointpeel_20260916/RESULTS.md) now implements
  a single joint queue. The shadow bound proves that each vertex's unfinished
  sizes form a prefix. At ordinary ceiling c(v)+1, matching lower/upper
  bounds finish that entire prefix with binom(c(v),s-1), without a clique
  witness. This is conditional: unsaturated higher sizes still need events.
  The initial full-target/bucket variants were slower. The final prefix
  version is 1.629x its single-size-event ablation, but only 0.678x the
  stronger same-index independent control in geometric mean. All modes
  reuse ordinary-core preprocessing for s=2; final queue timings cover
  s=3..10, one thread, three trials, five real graphs. Keep the failed
  performance evidence and both proofs; do not promote or claim a global
  novelty result. Source snapshots preserve the version that re-peeled s=2.
- Keep memory overhead visible. No global speed guarantee follows from r=1
  being a special case.
- The [tighter-bound follow-up](r1_envelopecore_20260916/RESULTS.md) uses a
  certified coloring bound and existing graph shadow theorems. It lifts
  their bounds to core values and completes non-clique extremal regions.
  With a tight coloring bound, every balanced complete multipartite graph
  with at least three parts uses one joint event per vertex after ordinary
  preprocessing. Unbalanced parts 1,2,3,4 still need two events per vertex.
  Only Stanford gains extra real coverage: 63 vertices, 441 fewer events.
  Four other graphs have b=d, so the new bound does not tighten their range.
  Global-frontier and separate phase experiments retain all regressions
  and noisy timings; neither establishes an overall independent-control
  win. Read THEORY.md, phase/THEORY.md, TABLES.md and REPRODUCE.md before
  revisiting this direction. No final all-s query index or production edit.
- The [floor-completion follow-up](r1_floorbatch_20260916/RESULTS.md) expands
  both the foundations and the new rule before implementation. A current
  count at a valid common layer floor proves the answer AND a safe
  minimum-label deletion. Consecutive sizes can share one source update;
  stop at the first failed size and cache that entry. Initial per-size
  minima are a separate optional seed, not free information. Five real
  graphs gain 195,352 extra completions, reducing joint events by 13.03%.
  Source reads fall 13.19%, but entry reads remain 76,085,751. The matched
  control/floor preparation-plus-peel ratio is 1.163; fixed/floor is 0.962.
  Stanford and Amazon remain slower than independent peeling. Seeding
  doubles batch coverage on three small unbalanced multipartite fixtures
  from two events to one per vertex, but adds no real event reduction
  beyond floor and adds substantial initialization traffic. Keep both
  optional variants; do not pool times with earlier experiments or promote
  them. A known label alone does not authorize removal; prove that too.
- The [count-skipping follow-up](r1_countskip_20260916/RESULTS.md) derives
  the sharp unrestricted-shadow threshold before C++ and lifts it to core
  witnesses. A higher answer below its own ceiling can already force a
  lower ceiling; retain the known answer but wait for safe deletion. This
  is an application of existing extremal mathematics, not a new classical
  shadow theorem. On five real inputs, 1,758,227 threshold tests establish
  130 certificates but skip only 71 queries and 423 incidence reads.
  Permanent pruning separately removes 33.19% of count reads and 47.94%
  of source reads; target reads fall only 0.12%. It adds tests, compaction
  writes and 4n bytes, consumes the reverse traversal view, and does not
  free index capacity. Control/prune time ratio is 0.931; fixed/prune is
  0.818. A separately proved scan-order implementation cuts pruning tests
  from 50,830,648 to 26,260,337 but still loses overall: old-prune/new-prune
  is 0.956 in that matched sweep. Do not pool the two sweeps or promote
  either variant. Both pass Release and ASan/UBSan, with 34,289 cases per
  suite, 1,042,090 compared values and 658,368 state audits. All 219 timed
  executions match verified outputs. Next investigate shared target-loss
  updates for already justified removals, with a new proof first; the
  current evidence does not establish that target loops dominate time.
- The [cross-vertex batching follow-up](r1_wavebatch_20260917/RESULTS.md)
  is [parked by explicit user decision](r1_wavebatch_20260917/DECISION.md).
  Keep all evidence, leave production unchanged, and do not continue this
  variant without a new request. The broader all-size goal is separate.
  The completed attempt
  implements equal-key snapshots and whole ordinary-core groups at a
  proved minimum band. Each vertex can finish a different size interval;
  selected survivors must receive batch losses at their cached next size.
  The proof precedes C++, with Release/ASan suites covering 34,306 cases,
  1,044,155 values and 480,094 state audits. All 168 timing runs match
  verified outputs. The same-batch sequential ablation also matches
  batch hashes, events, source/count reads and channel decrements.
  Aggregation saves only 804,767 of 63,855,947 target reads (1.26%),
  yet copies 146,066,522 snapshot channels. Relative to the preceding
  floor solver, events rise 1,303,918 to 1,327,468 and targets rise 4.30%.
  The [nine-vertex example](r1_wavebatch_20260917/SNAPSHOT_TRADEOFF.md)
  shows why: testing all lower layers before applying a batch misses
  completions that the preceding sequential deletion would enable.
  Base uses 9 vertex events, wave 11, with identical answers. Pure
  ordinary-ceiling grouping saves only 4,397 real target reads. Fixed/wave
  time ratio is 0.719; same-batch sequential/wave 0.785. Stanford state
  plus curve rises 37.109 to 84.111 MiB versus base. Do not promote, pool
  times, or mistake 36,493 batches for only 36,493 vertex computations.
  Before further batching, prove how newly available lower-size deletions
  can propagate without losing the intended traversal saving.

Earlier user-directed change of research direction (2026-09-17): use smaller-s
decomposition results in a bottom-up DP, not merely compatible deletions.
[Bottom-up results](r1_bottomup_20260917/RESULTS.md) preserve two attempts:
the first memoizes unresolved pivot-search states and updates future clique-
size bounds from finished core labels; the second initializes size t from
floor(kappa_{t-1}(v)*(c(v)-t+2)/(t-1)) and performs local threshold correction.
Its shadow and h-index ingredients are established prior work. The second
is closer to positive-value reuse but still repeats threshold counts.
Neither variant is a completed fast count-reusing DP or a Base replacement.
Read their proofs and complete failed timing evidence before another attempt.
Next priority: share the COUNT state, with a proof for changing threshold
sets, rather than assume that tighter upper labels alone remove the scans.

### Current Comparison Of The Two Roads

The user clarified that BOTH roads should continue to be considered:
**shared peeling across sizes** versus **bottom-up cross-s recurrence**.
The new [two-road report](r1_tworoads_20260917/RESULTS.md) keeps those
contracts separate. Triangle components permit independent higher-clique
peeling, but materializing local index copies is not a performance win.
Integer bounds from finished lower-s cores avoid many key decreases in
batch peeling, but the parent full heap retains expensive queue work.

The [sorted-order follow-up](r1_tworoads_20260917/stream/RESULTS.md) uses
a monotone transformed bound: the previous exact core order already sorts
the next upper bounds. Only vertices whose raw degree drops BELOW their
bound need a heap. Exceptions may leave out of the old order, so this is
not a common exact peeling order. At s=3 the bound equals the ordinary-core
bound; additional higher-layer information begins at s=4. Raw degrees are
still separately computed for each s. Classical shadow theory and generic
bounded peeling are prior ingredients, not new global novelty claims.

The parent and stream proofs preceded their C++. Release/ASan suites cover
34,296 and 34,292 graphs respectively, with independent definition oracles
and intermediate-state checks. Each sweep has 279 executions and its own
hash manifest; do not pool their timings. Their interactive-host timing
variation remains visible (88/90 and 55/90 noisy cells, max/min >1.2).

The additional [paired experiment](r1_tworoads_20260917/stream/PAIRED_TABLES.md)
uses one immutable index per graph, six balanced rounds, full output checks,
and both elapsed wall and process CPU time. Five real graphs, all s=2..10,
one thread. Graph-median paired ratios give fixed/DP 1.307 wall and 1.304
CPU; ordinary-stream/DP 1.062 wall and 1.076 CPU. GrQc loses to the ordinary
stream, Stanford is close/noisy, and shared floor remains faster on HepPh.
No universal winner or stable end-to-end gain has been established.

Over those five graphs, heap removals fall 2,329,788 -> 671,853 relative
to the matched full heap, or 1,127,614 -> 671,853 beyond the ordinary-only
stream. Comparisons fall 72.13% and 42.72%, respectively. Initial count
reads remain 64,829,935. Stream DP and full-heap DP also have identical
source/target reads. This is saved PRIORITY work, not yet saved cross-size
counting. The index/output sizes are unchanged; Stanford RSS in the broad
sweep is 282.484 MiB versus fixed's 270.500 MiB, so no memory-win claim.

Retain the new stream and its ordinary-only control as research baselines;
do not replace production or paper data. Next DP work must prove count-state
reuse for changing residuals. Next shared-peeling work must avoid copies
and show meaningful target savings. The user-parked wave variant is still
parked. [Commands](r1_tworoads_20260917/stream/REPRODUCE.md) reproduce each
evidence set separately; all successes, regressions and raw samples remain.

### Follow-Up: Count Reuse Did Not Yet Beat The Stream

[Inherited-threshold count reuse](r1_countreuse_20260917/RESULTS.md) maintains
present hold/pivot counts while descending upper-bound groups. An entire
group must support its threshold before its values become final. Accepted
vertices remain as support for unresolved vertices; only their own updates
and wholly accepted paths may be omitted. This is reuse within nested
tests inherited from the previous size, NOT raw-degree reuse between
arbitrary cross-size residuals.

The parent replaces 62,168,700 direct member scans with counter queries,
but finalizes only 61,783/2,329,788 positive pairs (2.65%). It saves
4,985,456 original count/source/target reads and adds 11,497,394 inspection
reads. Paired direct/reuse is 1.331, but stream/reuse is 0.827. The weaker
direct control is not sufficient evidence of a better overall algorithm.

The [guarded version](r1_countreuse_20260917/guarded/RESULTS.md) first checks
one vertex directly. Failure avoids whole-group counter insertion/rollback;
success still requires checking every other group member. Accepted sets,
fallback work and outputs are unchanged. Extra reads fall to 5,222,864,
yet paired old-reuse/guard is only 1.014 and stream/guard 0.823. Keep the
preceding stream; no production or paper replacement. Stanford extra state
is 42.583 versus stream's 37.137 MiB; index/output bytes are unchanged.

Both proofs precede C++. Each Release/ASan suite covers 34,296 graphs,
1,040,619 vertex-size values, 410,134 audits and 1,112 definition cases.
Parent evidence: 162 broad timings + nine verifications, 180 paired calls
and 30 warmups. Guarded evidence: 189 broad timings + nine verifications,
245 paired calls + 35 warmups. Every paired call checks full arrays.
No sweeps or pilots are pooled; noisy broad cells are 35/54 and 42/63.
[Notes](r1_countreuse_20260917/RESEARCH_NOTES.md) preserve failed directions,
the ALREADY known witnessed-binomial-tail theorem and the important fact
that a size-truncated path's whole vertex union need not be a clique.

## Related Work Outside This Directory

- [experiments/st_interval/REPORT.md](../experiments/st_interval/REPORT.md):
  interval/adaptive/compact layout comparisons preceding these research attempts.
- [tests/regression/CMakeLists.txt](../tests/regression/CMakeLists.txt):
  current production regression and explicit core/hierarchy oracles.
- [DO_NOT_REPEAT.md](../DO_NOT_REPEAT.md): broader project traps and NSI work;
  distinguish its general-r experiments from the r=1 research indexed here.
- [PROJECT_MEMORY.md](../PROJECT_MEMORY.md): user decisions, production defaults,
  paper constraints and older design history.

Maintenance rule: add one inventory row when a research directory is created;
update its conclusion after validation, including failure or inconclusive
results. Add a separate starred entry only when the user requests it. Do not
erase superseded runs or replace measured failures with a newer hypothesis.
