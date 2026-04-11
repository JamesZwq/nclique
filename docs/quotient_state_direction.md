# Quotient State Direction for General R

## Why this direction

The weak part of the current `r>=3` line is the state size.

The current algorithms still build and maintain many explicit `r`-cliques.
This is the main reason why the hardest cases become slow or time out.

The new idea is simple.
For a clean leaf, many `r`-cliques play the same role.
We should not keep them as separate states at the start.
We should first group them into quotient classes.

For a clean leaf, the first quotient class is the number of pivot vertices inside the `r`-clique.
If two `r`-cliques have the same number of pivot vertices, then they have the same clean-leaf contribution.

This does not solve the full problem yet.
But it gives a new basic state object.

## What was implemented

A new lab entry was added:

- `PIVOTER_RUN_ST_QUOTIENT_LAB=1`

File:

- [NucleusCoreDecompositionRCliqueST_QuotientLab.cpp](/Users/zhangwenqian/UNSW/pivoter/src/NucleusDecomposition/NucleusCoreDecompositionRCliqueST_QuotientLab.cpp)

This lab does two things.

1. It measures the compression from explicit clean-leaf `r`-cliques to quotient classes.
2. It then falls back to the exact `V20` path.

There is also a fast survey mode:

- `PIVOTER_QUOTIENT_LAB_ONLY=1`

This mode prints the quotient statistics and then stops.
It is useful for very hard graphs.

## Current evidence

All numbers below are from local runs.

### email-Eu-core, r=3, s=4

- Total clean quotient compression: `14.07x`
- Total one-step refined compression: `3.80x`
- Total one-removed delta compression: `1.85x`
- Largest clean class multiplicity: median `9`, p90 `42`, max `680`
- Largest delta class multiplicity: median `2`, p90 `3`, max `4`
- Median per-leaf compression: `5.00x`
- P90 per-leaf compression: `23.00x`
- Max per-leaf compression: `408.00x`

### com-youtube, r=3, s=4

- Total clean quotient compression: `5.77x`
- Total one-step refined compression: `2.04x`
- Total one-removed delta compression: `1.72x`
- Median per-leaf compression: `4.00x`
- P90 per-leaf compression: `9.50x`
- Max per-leaf compression: `340.00x`

### com-youtube, r=4, s=5

- Total clean quotient compression: `11.09x`
- Median per-leaf compression: `6.00x`
- P90 per-leaf compression: `25.50x`
- Max per-leaf compression: `1190.00x`

### web-Stanford, r=4, s=5

- Total clean quotient compression: `455.38x`
- Median per-leaf compression: `11.00x`
- P90 per-leaf compression: `325.00x`
- Max per-leaf compression: `260927.50x`

### web-it-2004, r=3, s=4

- Total clean quotient compression: `34055.01x`
- Total one-step refined compression: `11495.15x`
- Total one-removed delta compression: `2.64x`
- Total one-removed touch fraction: `0.004657%`
- Largest clean class multiplicity: median `120`, p90 `816`, max `13251095`
- Largest delta class multiplicity: median `3`, p90 `3`, max `4`
- Median per-leaf compression: `82.50x`
- P90 per-leaf compression: `484.50x`
- Max per-leaf compression: `6671880.00x`
- Median one-step refined compression: `27.50x`
- P90 one-step refined compression: `161.50x`
- Max one-step refined compression: `2223960.00x`
- Median one-removed delta compression: `3.00x`
- P90 one-removed delta compression: `3.00x`
- Max one-removed delta compression: `3.00x`
- Prototype quotient state memory: `80 MB`
- Prototype clean entries per leaf: avg `2.00`, max `2`
- Prototype delta entries per leaf: avg `5.00`, max `5`

### web-it-2004, r=4, s=5

- Total clean quotient compression: `2720020.37x`
- Total one-step refined compression: `703196.74x`
- Total one-removed delta compression: `3.50x`
- Total one-removed touch fraction: `0.00007795%`
- Largest clean class multiplicity: median `330`, p90 `3060`, max `1417867165`
- Largest delta class multiplicity: median `4`, p90 `4`, max `6`
- Median per-leaf compression: `247.50x`
- P90 per-leaf compression: `1938.00x`
- Max per-leaf compression: `715559130.00x`
- Median one-step refined compression: `26.25x`
- P90 one-step refined compression: `3425.62x`
- Max one-step refined compression: `178889782.50x`
- Median one-removed delta compression: `4.00x`
- P90 one-removed delta compression: `4.00x`
- Max one-removed delta compression: `4.00x`
- Prototype quotient state memory: `76 MB`
- Prototype clean entries per leaf: avg `2.00`, max `2`
- Prototype delta entries per leaf: avg `5.00`, max `5`

## What these numbers mean

These numbers are too large to ignore.

This is not a small constant-factor trick.
The state space is much larger than it needs to be.

The effect also gets stronger when `r` grows.
This is important.
The weakest part of the current system is exactly the `r>=3` regime.

The effect is also much stronger on dense leaves.
This is also important.
The hardest graphs create exactly these dense leaves.

The `web-it-2004` numbers are the strongest signal.
This is one of the hardest graphs in the current benchmark.
It is also where the current general-`r` line becomes weak.
The quotient compression there is enormous.
Even after one representative refinement step, the compression is still enormous.
This is the key sign that a refinement-based quotient algorithm may work.

There is also a second signal.
The one-removed delta classes are very small.
This means the exact update rule after the first damage can stay at the class level.

There is now a third signal.
The clean classes can be huge.
On `web-it-2004`, one clean class can represent millions or even billions of explicit `r`-cliques.
But the corresponding first-damage delta class is still tiny.
This is exactly the shape we want.

There is now a fourth signal.
The total touched part is tiny.
On `web-it-2004`, the first damage touches only `0.004657%` of the explicit state for `r=3,s=4`.
For `r=4,s=5`, it drops to `0.00007795%`.
This is the strongest evidence so far for lazy expansion.

The current lab also verifies a one-removed exact overlap formula on small leaves.
Enable:

- `PIVOTER_QUOTIENT_VERIFY=1`

Current sample checks return:

- `clean=OK`
- `delta=OK`
- `expand=OK`
- `update=OK`

The last two checks are important.
They show that a class can be expanded on demand, and can drive one exact damage update.

The lab can also build the prototype quotient state directly.
Enable:

- `PIVOTER_QUOTIENT_BUILD_STATE=1`

On `web-it-2004`, this prototype stays small.
This is strong evidence that a full quotient-state algorithm is implementable.

The lab now also has a first-round update prototype.
Enable:

- `PIVOTER_QUOTIENT_ONE_ROUND=1`

This prototype does one real peeling layer.
It finds the first removed `r`-cliques.
It collects the changed leaves.
It then emits exact per-clique deltas only for single-removed clean leaves.

Current local results:

- `email-Eu-core, r=3, s=4`
  - removed `1391` r-cliques in the first frontier
  - changed `949` leaves
  - single-removed leaves `577`
  - multi-removed leaves `372`
  - quotient-handled leaves `577`
  - handled explicit entries `2308`
  - handled touched cliques `1731`
  - handled touch fraction `75.0%`
  - counting lookup misses `0`
  - delta lookup misses `0`
- `com-youtube, r=3, s=4`
  - removed `428986` r-cliques in the first frontier
  - changed `253102` leaves
  - single-removed leaves `121545`
  - multi-removed leaves `131557`
  - quotient-handled leaves `121545`
  - handled explicit entries `486180`
  - handled touched cliques `364635`
  - handled touch fraction `75.0%`
  - counting lookup misses `0`
  - delta lookup misses `0`

This gives a more careful picture.
The first-round path already has real coverage.
But on these two graphs, the handled leaves still touch a large part of their explicit state.
So the full gain will not come from the easy first frontier alone.

There is now one more important result.
The old state-to-index alignment issue is gone.
The first-round prototype now uses the same fused full-enum path for counting and local leaf keys.
On both `email-Eu-core` and `com-youtube`, the lookup misses are now zero.

This changes the diagnosis.
The next blocker is no longer local alignment.
The next blocker is the global clique-id layer itself.

There is now another new result.
The one-round prototype also has a sparse-support version.
Enable:

- `PIVOTER_QUOTIENT_ONE_ROUND_SPARSE_SUPPORT=1`

This version does not use `StaticCliqueIndex` to maintain the first-round support map.
It uses a sparse global tuple key instead.

Current local checks:

- `email-Eu-core, r=3, s=4`
  - sparse support states `104730`
  - removed sparse cliques `1391`
  - changed leaves `949`
  - single-removed leaves `577`
  - quotient-handled leaves `577`
  - unique sparse deltas `1583`
- `com-youtube, r=3, s=4`
  - sparse support states `2564310`
  - removed sparse cliques `428986`
  - changed leaves `253102`
  - single-removed leaves `121545`
  - quotient-handled leaves `121545`
  - unique sparse deltas `277866`

These numbers match the clique-id prototype on the same graphs.
So the first-round support layer can also be written without `StaticCliqueIndex`, at least on manageable graphs.

There is now a stronger prototype as well.
Enable:

- `PIVOTER_QUOTIENT_ONE_ROUND_SPARSE_EXACT=1`

This version does an exact first round with a sparse global support map.
It does not stop at the first frontier.
It also applies an exact BK-based local delta on every changed leaf.
After that, it computes the next frontier on the updated sparse support map.

Current local results:

- `email-Eu-core, r=3, s=4`
  - sparse support states after round 1: `103131`
  - removed sparse cliques in round 1: `1391`
  - changed leaves: `949`
  - exact-handled leaves: `949`
  - next frontier size: `192`
- `com-youtube, r=3, s=4`
  - sparse support states after round 1: `2060336`
  - removed sparse cliques in round 1: `428986`
  - changed leaves: `253102`
  - exact-handled leaves: `253102`
  - next frontier size: `37491`

This is the first time the quotient line has a real end-to-end peeling prototype.
It is still only one round.
But it now has:

- a global sparse support state
- exact changed-leaf routing
- exact BK-based local delta
- an updated next frontier

So the work is now past the old “local formula only” stage.

There is now a multi-round version for the no-split regime.
Enable:

- `PIVOTER_QUOTIENT_MULTI_ROUND_SPARSE=1`

This version keeps running as long as the changed leaves all die after each round.
It stops if a true split appears.

Current local results:

- `email-Eu-core, r=3, s=4`
  - completed rounds: `5`
  - split encountered: `NO`
  - dead leaves so far: `1113`
  - remaining sparse states: `102884`
- `com-youtube, r=3, s=4`
  - completed rounds: `3`
  - split encountered: `NO`
  - dead leaves so far: `281067`
  - remaining sparse states: `2013096`

This is a stronger milestone.
The quotient line now has a real peeling prototype that runs for multiple rounds.
It is not full yet.
But it is no longer only a collection of local formulas.

There is now a stronger structural signal.
Many active leaves have `keep > r`.
Current local counts:

- `email-Eu-core, r=3, s=4`
  - leaves with `keep > r`: `15669`
  - explicit entries from these leaves: `62676`
- `com-youtube, r=3, s=4`
  - leaves with `keep > r`: `317076`
  - explicit entries from these leaves: `1268304`

This matters.
It means the general-`r` support layer cannot rely too heavily on a model that fixes the whole keep set and only chooses extra pivots.
The quotient view is better aligned with these leaves.
It counts by state classes, not by one rigid keep/pivot split.

There is now a direct coverage check as well.
Enable:

- `PIVOTER_QUOTIENT_COMPARE_INDEX=1`

This check builds `StaticCliqueIndex` and then asks a simple question:
for positive-contribution cliques that come from `keep > r` leaves, how many are covered by that index?

Current local results:

- `email-Eu-core, r=3, s=4`
  - unique `keep > r` positive cliques: `25561`
  - covered by `StaticCliqueIndex`: `25524`
  - missed by `StaticCliqueIndex`: `37`
- `com-youtube, r=3, s=4`
  - unique `keep > r` positive cliques: `539899`
  - covered by `StaticCliqueIndex`: `533040`
  - missed by `StaticCliqueIndex`: `6859`

This is a stronger result than the earlier size-gap signal.
It shows that the old global index path is not perfectly aligned with these leaves.
So the quotient direction is not only about speed.
It is also a better match to the real state structure.

On `web-it-2004, r=3, s=4`, the one-round prototype reaches the fused build stage and prints:

- `Clique Index: 338761454 cliques, 7178413 vertices.`

This is the next hard fact.
The quotient-state idea already shrinks the leaf-local state by many orders of magnitude.
But the current prototype still asks for a full global explicit clique index.
This is too large on the hardest graph.
So the next algorithm must also weaken or delay the global clique-id layer.

## Current conclusion

The quotient idea is worth continuing.

The next version should not be another small optimization on top of `V20`.
It should build a new exact state layer.

The clean target is:

- do not materialize all explicit `r`-cliques at the start
- maintain quotient classes first
- refine the classes only when symmetry is broken
- use exact overlap-based delta after the first damage

## Next exact algorithm step

The next algorithm should have four stages.

1. Build clean-leaf quotient states.

For each leaf, keep only the classes indexed by `subNumPivot`.
Each class stores its multiplicity and its current clean contribution.

2. Delay explicit `r`-clique materialization.

Do not expand a class into explicit `r`-cliques unless the later updates really need it.

3. Refine only affected leaves.

When deletions break clean symmetry, refine only the touched classes.
Do not rebuild the whole leaf.

For the first removed `r`-clique in a clean leaf, use the exact overlap formula.
The update only depends on:

- the number of pivot vertices in the surviving `r`-clique
- the number of shared pivot vertices with the removed `r`-clique

4. Expand only near the peel frontier.

If a quotient class must interact with exact peeling, expand only that class.
Do not expand the whole leaf.

## What this would contribute

If this works, the new contribution is not just a faster implementation.

The new contribution would be:

- a new exact quotient-state representation for `r>=3`
- a refinement-based peeling process
- a shift from explicit clique states to class states

This is the right level for a stronger journal story.

## Important limitation

The current data also shows one limit.

The biggest gain is in the state size.
The gain in one-step delta compression is smaller.

This means the next algorithm should focus on:

- compact clean-leaf state
- exact overlap-based delta
- lazy expansion near the peel frontier

It should not rely on a huge per-step delta compression alone.

There is also one concrete engineering blocker now.
It is not the local lookup any more.
It is the global support layer on the hardest graphs.
The clique-id version is one instance of this problem.
The sparse-support version shows that the same pressure remains even after we remove `StaticCliqueIndex`.
The next step is not more benchmarking.
The next step is to avoid building the full global first-round support map on the hardest graphs.

## New dynamic prototype result

We now have a split-aware multi-round sparse prototype.
It is implemented in `NucleusCoreDecompositionRCliqueST_QuotientLab.cpp`.

This prototype keeps:

- a global sparse support map
- dynamic active leaves
- a dynamic `vertex -> active leaf` index
- BK-based local split for changed leaves

This is not only a one-round study any more.
It can roll multiple rounds.
It can also continue after true split happens.

## What happens before split

For the first few rounds, the changed leaves often die directly.
They do not split.

We observed:

- `email-Eu-core, r=3, s=4`: 5 rounds, 0 spawned subleaves
- `com-youtube, r=3, s=4`: 6 rounds, 0 spawned subleaves
- `web-Stanford, r=3, s=4`: 4 rounds, 0 spawned subleaves
- `web-Stanford, r=4, s=5`: 3 rounds, 0 spawned subleaves

This is important.
It means the quotient dynamic path already covers the low-core frontier on several graphs.

## What happens after split

The more important result is that split does appear when `s-r > 1`.
We observed this on `email-Eu-core, r=3, s=5`.

The key milestones are:

- first round with core value larger than 1: round 7
- first round with true split: round 14
- first split wave: 409 spawned subleaves in round 14

The prototype can continue after that.
With `200` rounds on the same setting, we get:

- `96` rounds with split
- `1690` split leaves
- `1557` single-child leaves
- `133` multi-child leaves
- `1825` spawned leaves
- maximum spawned leaves in one round: `561`
- maximum observed core value: `8`

So the current prototype is already beyond the no-split regime.
It can survive repeated split rounds on a real graph.

## Updated conclusion

This changes the status of the direction.

The main blocker is no longer the local split logic.
The main blocker is now the global support layer on the hardest graphs.

More precisely:

- local quotient state works
- local split works
- dynamic multi-round update works
- repeated split rounds also work on small graphs

The remaining hard problem is:

- how to avoid building or maintaining a huge global sparse support map on the hardest graphs

So the next exact algorithm step should target the global frontier layer.
It should not go back to tuning local formulas.

## New support-sharing result

We also measured a key structural quantity for the global sparse support layer.

For each sparse `r`-clique key, we counted how many leaves contribute to it.
This is the leaf-incidence of the key.

This test answers one important question:

- can we split the global support layer into
  - mostly unique keys kept inside leaves
  - only a small shared-key map

The answer is no.

### `email-Eu-core, r=3, s=5`

- sparse keys: `102747`
- total leaf-key incidence: `3246772`
- mean leaf incidence: `31.60`
- unique-key ratio: `5.88%`
- shared-incidence ratio: `99.81%`
- incidence median / p90 / max: `13 / 80 / 1032`

### `com-youtube, r=3, s=4`

- sparse keys: `2564310`
- total leaf-key incidence: `14335854`
- mean leaf incidence: `5.59`
- unique-key ratio: `28.78%`
- shared-incidence ratio: `94.85%`
- incidence median / p90 / max: `3 / 13 / 1179`

### `web-Stanford, r=3, s=4`

- sparse keys: `11156923`
- total leaf-key incidence: `115534074`
- mean leaf incidence: `10.36`
- unique-key ratio: `28.88%`
- shared-incidence ratio: `97.21%`
- incidence median / p90 / max: `4 / 16 / 18935`

## What this means

This result is very important.

The main blocker is not just the number of keys.
It is also the amount of cross-leaf sharing.

So the next exact algorithm should not assume:

- most keys are unique
- only a small fraction of keys need a global map

That assumption is false on the tested graphs.

The next direction should be stronger.
It should try to avoid a full global support layer in a different way.

The most likely options now are:

- a frontier-local support layer
- a support layer that is built by core level, not all at once
- a delayed global materialization strategy

## New frontier-sharing result

We then checked a more specific question.

The overall global sharing is strong.
But do the *frontier* keys also have strong sharing?

Here the frontier means the keys with the current minimum core value.

This result is different from the global one.
For the tested points, the first frontier is completely leaf-unique.

### `email-Eu-core, r=3, s=5`

- frontier min core: `1`
- frontier keys: `1981`
- frontier mean incidence: `1.00`
- frontier unique-key ratio: `100%`

### `com-youtube, r=3, s=4`

- frontier min core: `1`
- frontier keys: `428986`
- frontier mean incidence: `1.00`
- frontier unique-key ratio: `100%`

### `web-Stanford, r=3, s=4`

- frontier min core: `1`
- frontier keys: `124970`
- frontier mean incidence: `1.00`
- frontier unique-key ratio: `100%`

## Updated interpretation

This changes the picture again.

Globally, the sparse support layer is highly shared.
But the first peel frontier is not.

So a naive global factorization is still wrong.
But a frontier-local exact algorithm becomes plausible again.

The next promising idea is now:

- do not build the whole global support map
- try to build only the current frontier layer
- then expand to the next needed core level

In other words, the next target is not `global sparse support`.
The next target is `frontier-first exact support`.

## New low-support prototype

We then built a stronger prototype for the frontier idea.

This prototype does not keep exact support for every key.
It only keeps keys whose support is at most `tau`.
Once a key exceeds `tau`, it is moved out of the exact low-support map.

This gives an exact filter for the current low-core layer.

## Exact result for `tau = 1`

On the tested graphs, the `tau=1` prototype exactly matches the true first frontier.

### `email-Eu-core, r=3, s=5`

- low-support keys kept: `1981`
- exact frontier size: `1981`
- exact agreement: `YES`

### `com-youtube, r=3, s=4`

- low-support keys kept: `428986`
- exact frontier size: `428986`
- exact agreement: `YES`

### `web-Stanford, r=3, s=4`

- low-support keys kept: `124970`
- exact frontier size: `124970`
- exact agreement: `YES`

So the first frontier can already be extracted exactly without storing all exact supports.

## Exact result for larger `tau`

We also tested larger low-support layers.

### `email-Eu-core, r=3, s=5`

- `tau=2`: exact low-support keys `3013`, agreement `YES`
- `tau=3`: exact low-support keys `4845`, agreement `YES`

### `com-youtube, r=3, s=4`

- `tau=2`: exact low-support keys `770205`, agreement `YES`

This is important.
The low-support layer grows, but it is still much smaller than the full global support layer.

## Updated next step

We now have stronger evidence for a new exact direction.

The next exact algorithm should try to:

- build the `tau=1` frontier exactly without a full support map
- remove the current frontier
- then raise `tau` or rebuild the next low-support layer as needed

This is now more than a guess.
We already have an exact prototype that isolates the low-support layer.

## New first-frontier exact update prototype

We pushed the low-support idea one step further.

We now have a `first-frontier low-support` prototype.
It does three things:

- build the first frontier using only the exact `tau=1` low-support map
- collect changed leaves from that frontier
- run exact local BK-based delta on those changed leaves

This prototype does **not** build the full global support map first.

### `email-Eu-core, r=3, s=5`

- frontier keys: `1981`
- changed leaves: `1081`
- exact-handled leaves: `1081`
- single-removed leaves: `680`
- multi-removed leaves: `401`
- generated subleaves: `0`
- unique sparse deltas: `8976`

### `com-youtube, r=3, s=4`

- frontier keys: `428986`
- changed leaves: `253102`
- exact-handled leaves: `253102`
- single-removed leaves: `121545`
- multi-removed leaves: `131557`
- generated subleaves: `0`
- unique sparse deltas: `790262`

These numbers match the earlier first-round exact prototype.
The difference is only the support layer:

- old prototype: build full global sparse support first
- new prototype: build only the exact low-support frontier first

## Updated status

This is the strongest result so far for the new direction.

We now have:

- exact quotient clean state
- exact dynamic multi-round split handling
- exact frontier isolation without full support
- exact first-frontier update without full support

So the next algorithm step is now very concrete.

The next exact target is:

- keep the current frontier-first support layer
- update to the next support layer after the frontier is removed
- avoid ever materializing the full global support map on the hardest graphs

## New multi-round low-support prototype

We then pushed the idea one step further.

We now have a `tau`-limited multi-round prototype.
It keeps only the exact keys whose support is at most `tau`.
It does not build the full global support map.

The update rule is:

- remove the current frontier from the low-support layer
- find changed leaves
- recompute only the touched candidate keys exactly
- keep untouched low-support keys unchanged
- rebuild the next low-support layer from these two parts

This is now a real peeling prototype.

## `tau = 2` on `email-Eu-core, r=3, s=5`

The prototype reproduces the exact low-core sequence:

- round 1: `1981 -> 385`
- round 2: `385 -> 103`
- round 3: `103 -> 16`
- round 4: `16 -> 8`
- round 5: `8 -> 1`
- round 6: `1 -> 318` at core value `2`
- round 7: `318 -> 23`
- round 8: `23 -> 3`
- round 9: `3 -> 2`
- round 10: `2 -> 2`
- round 11: `2 -> 24`
- round 12: `24 -> 4`

So `tau=2` is already enough to keep peeling exactly through all support-1 and support-2 layers.

## `tau = 2` on `com-youtube, r=3, s=4`

The same prototype also tracks the low-core frontier for several rounds:

- round 1: `428986 -> 37491`
- round 2: `37491 -> 5522`
- round 3: `5522 -> 1193`
- round 4: `1193 -> 329`
- round 5: `329 -> 90`
- round 6: `90 -> 35`

This again matches the earlier exact multi-round prototype.

## `tau = 3` on `email-Eu-core, r=3, s=5`

This is the strongest result so far.

With `tau=3`, the low-support prototype crosses the first true split wave:

- round 14: support-3 frontier of size `2500`
- generated subleaves: `409`
- next frontier: `836` with support `1`

It then continues further:

- round 15: `836 -> 14`
- round 16: `14 -> 1`
- round 17: `1 -> 33` at support `2`
- round 18: `33 -> 4`
- round 19: `4 -> 566` at support `3`
- round 20: `566 -> 120` with `94` new subleaves

So this prototype is no longer only a first-frontier tool.
It can already peel many exact low-support rounds without a full global support map.

## Updated status

This is now the strongest evidence for the new direction.

We no longer only have:

- quotient state compression
- local delta formulas
- split-aware dynamic local updates

We now also have:

- exact frontier extraction without full support
- exact next-frontier update without full support
- exact multi-round peeling for bounded low-support layers

The next step is now even clearer.

The next exact algorithm should extend this bounded low-support peeling outward.
In other words, it should turn the current `tau`-limited prototype into a full level-by-level peeling algorithm.

## New adaptive-`tau` prototype

We then removed one more restriction.

The previous prototype fixed `tau` in advance.
The new prototype raises `tau` only when the current low-support layer is exhausted.

So it works in phases:

1. peel all exact keys with support at most current `tau`
2. when this layer becomes empty, increase `tau`
3. rebuild the low-support layer from the current active leaves
4. continue peeling

This is the first prototype that starts to look like a complete level-by-level exact algorithm.

### `email-Eu-core, r=3, s=5`

With

- `tau_start = 2`
- `tau_max = 3`

the prototype already reproduces the known exact low-core trajectory:

- rounds `1-13`: all support-1 and support-2 layers
- round `14`: first support-3 split wave, `2500` frontier keys, `409` spawned subleaves
- rounds `15-20`: continues after the split
- round `20`: second support-3 split wave, `566` frontier keys, `94` spawned subleaves

The summary is:

- tau phases used: `2`
- completed rounds: `20`
- total exact leaves handled: `4139`
- total candidate keys recomputed: `33273`
- total spawned subleaves: `503`
- maximum observed frontier support: `3`

## Updated interpretation

This is a real step forward.

We no longer need to fix one `tau` and stop there.
The prototype can now:

- finish one low-support band
- raise the band
- rebuild from the current active leaves
- continue exact peeling

So the next goal is no longer conceptual.
It is now implementation-oriented:

- improve the rebuild cost between `tau` phases
- push the adaptive prototype to larger graphs
- check how far it can go before the low-support layer becomes too large

## Stronger result on `email-Eu-core`

We then pushed the adaptive prototype much further.

With

- `tau_start = 1`
- `tau_max = 8`
- `max_rounds = 200`

the adaptive prototype already reproduces the full 200-round dynamic run on `email-Eu-core, r=3, s=5`.

The summary is:

- tau phases used: `8`
- completed rounds: `200`
- total exact leaves handled: `9757`
- total candidate keys recomputed: `96576`
- total spawned subleaves: `1825`
- maximum observed frontier support: `8`
- total rebuild time: `708 ms`
- total prototype time: `1401 ms`

The phase breakdown is also informative:

- `tau=1`: rebuild `102 ms`, `6` rounds
- `tau=2`: rebuild `98 ms`, `7` rounds
- `tau=3`: rebuild `87 ms`, `25` rounds
- `tau=4`: rebuild `82 ms`, `13` rounds
- `tau=5`: rebuild `88 ms`, `25` rounds
- `tau=6`: rebuild `84 ms`, `51` rounds
- `tau=7`: rebuild `84 ms`, `43` rounds
- `tau=8`: rebuild `83 ms`, `30` rounds

So the adaptive path is no longer a short prefix prototype.
It already behaves like a full level-by-level exact peeling engine on this graph.

## First large-graph match

We then checked whether the adaptive path still agrees with the older full dynamic sparse prototype on larger graphs.

### `com-youtube, r=3, s=4`

Adaptive low-support with

- `tau_start = 1`
- `tau_max = 3`
- `max_rounds = 30`

gives:

- completed rounds: `30`
- total exact leaves handled: `580551`
- total spawned subleaves: `61721`
- maximum observed frontier support: `2`
- total rebuild time: `3004 ms`
- total prototype time: `13748 ms`

The old full dynamic sparse prototype with the same `30` rounds gives:

- completed rounds: `30`
- dead leaves so far: `580551`
- spawned leaves so far: `61721`
- maximum observed core: `2`
- prototype time: `17688 ms`

So the two paths agree on the main exact state counts for these `30` rounds.
The adaptive path is also faster here.

### `web-Stanford, r=3, s=4`

Adaptive low-support with

- `tau_start = 1`
- `tau_max = 3`
- `max_rounds = 30`

gives:

- completed rounds: `30`
- total exact leaves handled: `208144`
- total spawned subleaves: `47770`
- maximum observed frontier support: `2`
- total rebuild time: `12200 ms`
- total prototype time: `13627 ms`

The old full dynamic sparse prototype with the same `30` rounds gives:

- completed rounds: `30`
- dead leaves so far: `208144`
- spawned leaves so far: `47770`
- maximum observed core: `2`
- prototype time: `59606 ms`

So the two paths again agree on the main exact state counts.
Here the adaptive path is much faster.

## New interpretation

This changes the status of the project.

We now have more than a local quotient idea.
We now have an exact peeling prototype that:

- avoids a full global support map
- continues across many rounds
- survives true split waves
- matches the older dynamic prototype on the tested large-graph prefixes

The next bottleneck is no longer the local update formula.
The next bottleneck is the rebuild cost of the adaptive low-support layer on larger graphs.

## New buffered adaptive prototype

We then made one more algorithmic change.

The previous adaptive prototype rebuilt the low-support layer at every new `tau`.
The new buffered version keeps a small lookahead window.

With lookahead `1`, one rebuild keeps all keys with support at most `tau+1`.
So one rebuild can cover two support bands.

This changes the rebuild pattern.

## `email-Eu-core, r=3, s=5`

With

- `tau_start = 1`
- `tau_max = 8`
- `lookahead = 1`
- `max_rounds = 200`

the buffered prototype still matches the same 200-round exact trajectory:

- completed rounds: `200`
- total exact leaves handled: `9757`
- total candidate keys recomputed: `96576`
- total spawned subleaves: `1825`
- maximum observed frontier support: `8`

But the rebuild side becomes much smaller:

- rebuild count: `4` instead of `8`
- total rebuild time: `359 ms` instead of `708 ms`
- total prototype time: `1132 ms` instead of `1401 ms`

So the buffered version is clearly better here.

## `com-youtube, r=3, s=4`

With

- `tau_start = 1`
- `tau_max = 3`
- `lookahead = 1`
- `max_rounds = 30`

the buffered prototype still matches the same `30`-round exact trajectory:

- completed rounds: `30`
- total exact leaves handled: `580551`
- total candidate keys recomputed: `1238336`
- total spawned subleaves: `61721`
- maximum observed frontier support: `2`

The rebuild side improves:

- rebuild count: `1` instead of `2`
- total rebuild time: `1907 ms` instead of `3004 ms`

But the total runtime is almost the same:

- buffered: `13986 ms`
- plain adaptive: `13748 ms`

So buffering is not a universal win.
The larger buffered low-support layer adds extra round cost here.

## `web-Stanford, r=3, s=4`

With the same setting,

- completed rounds: `30`
- total exact leaves handled: `208144`
- total candidate keys recomputed: `287114`
- total spawned subleaves: `47770`
- maximum observed frontier support: `2`

Again the exact state counts match.
But now the speedup is large:

- rebuild count: `1` instead of `2`
- total rebuild time: `7638 ms` instead of `12200 ms`
- total prototype time: `10185 ms` instead of `13627 ms`

So buffering helps a lot when rebuild dominates.

## Updated interpretation

We now have three increasingly strong versions:

- full dynamic sparse support
- adaptive low-support
- buffered adaptive low-support

The buffered version shows that rebuild cost is not just a measurement artifact.
It is a real algorithmic bottleneck.

But it also shows something more subtle.
Reducing rebuild count is not enough by itself.
The buffered low-support map becomes larger.
So the best next step is not "always use more buffering".

The next step should be selective buffering.
It should turn buffering on only when rebuild dominates more than the extra buffered round cost.

## One more check: `lookahead = 2`

We also tested a larger buffer.

This gives one more clear signal.

### `email-Eu-core, r=3, s=5`

With `lookahead = 2`:

- rebuild count: `3`
- total rebuild time: `287 ms`
- total prototype time: `1018 ms`

This is even better than `lookahead = 1`.

### `com-youtube, r=3, s=4`

With `lookahead = 2`:

- rebuild count: `1`
- total rebuild time: `1680 ms`
- total prototype time: `15557 ms`

This is worse than both:

- plain adaptive: `13748 ms`
- buffered `lookahead = 1`: `13986 ms`

So aggressive buffering hurts here.

### `web-Stanford, r=3, s=4`

With `lookahead = 2`:

- rebuild count: `1`
- total rebuild time: `8003 ms`
- total prototype time: `12619 ms`

This is still better than plain adaptive (`13627 ms`).
But it is worse than buffered `lookahead = 1` (`10185 ms`).

## Final interpretation of this stage

The new lesson is simple.

There is no single best fixed lookahead.

- small graphs with many phases can benefit from larger buffering
- some large graphs benefit from `lookahead = 1`
- some graphs get slower if the buffered low-support map grows too much

So the next real algorithmic step is now very clear:

we need an adaptive buffering rule, not a fixed buffering rule.

## First adaptive buffering rule

We then implemented a simple automatic rule.

The rule is:

- for the first phase, use the larger buffer only if the raw clique layer is small
- after that, use the larger buffer only if the previous low-support layer is still small

In the current prototype, we used:

- max lookahead `2`
- low-support threshold `300000`
- raw clique threshold `1000000`

This rule is simple.
But it already improves the fixed-buffer story.

### `email-Eu-core, r=3, s=5`

The auto-buffered version gives:

- completed rounds: `200`
- total exact leaves handled: `9757`
- total candidate keys recomputed: `96576`
- total spawned subleaves: `1825`
- total rebuild time: `258 ms`
- total prototype time: `926 ms`

This is better than all previous versions:

- adaptive: `1401 ms`
- buffered `lookahead=1`: `1132 ms`
- buffered `lookahead=2`: `1018 ms`

### `com-youtube, r=3, s=4`

The auto-buffered version gives:

- completed rounds: `30`
- total exact leaves handled: `580551`
- total candidate keys recomputed: `1238336`
- total spawned subleaves: `61721`
- total rebuild time: `1581 ms`
- total prototype time: `13346 ms`

This is also the best result so far on this prefix:

- adaptive: `13748 ms`
- buffered `lookahead=1`: `13986 ms`
- buffered `lookahead=2`: `15557 ms`

### `web-Stanford, r=3, s=4`

The auto-buffered version gives:

- completed rounds: `30`
- total exact leaves handled: `208144`
- total candidate keys recomputed: `287114`
- total spawned subleaves: `47770`
- total rebuild time: `7055 ms`
- total prototype time: `9473 ms`

Again this is the best result so far on this prefix:

- adaptive: `13627 ms`
- buffered `lookahead=1`: `10185 ms`
- buffered `lookahead=2`: `12619 ms`

## Updated conclusion

This is the strongest result so far for the quotient line.

We no longer only know that rebuild policy matters.
We now have a first adaptive rebuild policy that improves all three tested prefixes.

It is still simple.
So it is probably not the final rule.
But it proves that the next layer of progress is real:

- not a new local formula
- not a new split case
- but a better policy for when and how much low-support state to rebuild

## Banded buffered window

We then changed the buffered layer itself.

The old buffered version keeps one exact map for all keys up to the phase cap.
This reduces rebuilds.
But it also makes each round scan a larger map.

The new banded version uses two exact layers:

- an active low-support layer
- a buffered next layer

The active layer is scanned in the current rounds.
The buffered layer is kept exact, but not scanned yet.
When the active layer is exhausted, we promote the buffered layer without rebuilding.

This gives a cleaner tradeoff:

- fewer rebuilds than the plain adaptive version
- less per-round overhead than a full larger window

### `email-Eu-core, r=3, s=5`

With banded auto-buffering, we get:

- completed rounds: `200`
- total exact leaves handled: `9757`
- total candidate keys recomputed: `96576`
- total spawned subleaves: `1825`
- total rebuild time: `242 ms`
- total prototype time: `870 ms`

This improves over the previous auto-buffered result:

- auto-buffered: `926 ms`

### `com-youtube, r=3, s=4`

With banded auto-buffering, we get:

- completed rounds: `30`
- total exact leaves handled: `580551`
- total candidate keys recomputed: `1238336`
- total spawned subleaves: `61721`
- total rebuild time: `1698 ms`
- total prototype time: `12421 ms`

This improves over the previous auto-buffered result:

- auto-buffered: `13346 ms`

### `web-Stanford, r=3, s=4`

With banded auto-buffering, we get:

- completed rounds: `30`
- total exact leaves handled: `208144`
- total candidate keys recomputed: `287114`
- total spawned subleaves: `47770`
- total rebuild time: `6077 ms`
- total prototype time: `7924 ms`

This also improves over the previous auto-buffered result:

- auto-buffered: `9473 ms`

### What this means

This result is important.

The plain buffered prototype already told us that rebuild count matters.
The banded buffered prototype shows something stronger.

The next support layer is useful.
But it should not be mixed into the current active map too early.

This is why the banded design works:

- it keeps the next layer exact
- it delays the cost of scanning that layer

We also tested fixed banded `lookahead=2` on the larger graphs.
It was not better:

- `com-youtube, 30 rounds`: `13230 ms`
- `web-Stanford, 30 rounds`: `9527 ms`

So the main lesson is not "always use a larger band."
The real lesson is:

- keep the extra layer exact
- but only promote it when it becomes active

## Deeper prefix with the banded auto policy

We then pushed the best current version further.

We used:

- banded buffered low-support
- auto lookahead
- `tau_max=5`
- `200` rounds

### `web-Stanford, r=3, s=4`

This run gives:

- completed rounds: `200`
- total exact leaves handled: `855052`
- total candidate keys recomputed: `1598245`
- total spawned subleaves: `161003`
- max observed round min: `4`
- total rebuild time: `12178 ms`
- total prototype time: `30437 ms`

The phase breakdown is:

- `tau=1`, `cap=2`, `79` rounds
- `tau=3`, `cap=4`, `121` rounds

This is important.

The quotient prototype is no longer confined to the first low-support prefix.
On this graph, it already crosses into the `core=4` regime while still using the same exact frontier-first machinery.

### `com-youtube, r=3, s=4`

This run gives:

- completed rounds: `200`
- total exact leaves handled: `994284`
- total candidate keys recomputed: `3026028`
- total spawned subleaves: `170736`
- max observed round min: `4`
- total rebuild time: `2677 ms`
- total prototype time: `31608 ms`

The phase breakdown is:

- `tau=1`, `cap=2`, `54` rounds
- `tau=3`, `cap=4`, `146` rounds

So the same pattern also appears on `com-youtube`.

The current best quotient prototype can now move beyond the shallow low-support prefix on two large graphs.
It already reaches the `core=4` regime without falling back to the old full global clique state.

## Much deeper prefix on `web-Stanford`

We then pushed the same prototype further on `web-Stanford`.

We used:

- banded buffered low-support
- auto lookahead
- `tau_max=7`
- `400` rounds

This run gives:

- completed rounds: `400`
- total exact leaves handled: `956396`
- total candidate keys recomputed: `2382815`
- total spawned subleaves: `197964`
- max observed round min: `7`
- total rebuild time: `21326 ms`
- total prototype time: `42306 ms`

The phase sequence is:

- `tau=1`, `cap=2`, `79` rounds
- `tau=3`, `cap=4`, `158` rounds
- `tau=5`, `cap=6`, `141` rounds
- `tau=7`, `cap=7`, `22` rounds

This is a much stronger result.

The quotient line is no longer just a first-frontier or low-core prototype.
On this graph, it already moves through several exact core levels in one end-to-end run.

## Longer run on `com-youtube`

We also pushed the same setting on `com-youtube`.

We used:

- banded buffered low-support
- auto lookahead
- `tau_max=7`
- `400` rounds

This run gives:

- completed rounds: `400`
- total exact leaves handled: `1103884`
- total candidate keys recomputed: `3657544`
- total spawned subleaves: `203917`
- max observed round min: `4`
- total rebuild time: `2860 ms`
- total prototype time: `39031 ms`

The phase sequence is:

- `tau=1`, `cap=2`, `54` rounds
- `tau=3`, `cap=4`, `346` rounds

This result is different from `web-Stanford`.

The prototype still reaches `core=4`.
But most of the work is concentrated in one long `tau=3` phase.

So the next bottleneck is now clearer:

- on some graphs, the frontier-first quotient line can already move deep
- but a large mid-core phase can still dominate the total work

## Bucketed active layer

We then attacked the next bottleneck directly.

The banded prototype still scans the whole active low-support map in each round:

- once to find the current minimum
- once to collect the frontier

This is expensive on a long mid-core phase.

So we built a new version with exact support buckets for the active layer.

The idea is simple:

- keep the same exact low and buffered layers
- but store the active layer by support bucket
- update only frontier keys and touched keys
- stop rescanning the whole active map every round

### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

The old banded auto-buffered version took:

- `39031 ms`

The new bucketed version takes:

- `37821 ms`

Other totals stay aligned:

- completed rounds: `400`
- max observed round min: `4`
- total exact leaves handled: `1103883`
- total candidate keys recomputed: `3657546`
- total spawned subleaves: `203916`

So this change does not alter the exact computation.
It reduces the cost of the long mid-core phase itself.

### `email-Eu-core, r=3, s=5`, `tau_max=8`, `200` rounds

The earlier banded auto-buffered version took:

- `870 ms`

The bucketed version takes:

- `840 ms`

The main totals stay aligned:

- completed rounds: `200`
- total exact leaves handled: `9757`
- total candidate keys recomputed: `96576`
- total spawned subleaves: `1825`
- max observed round min: `8`

### `web-Stanford, r=3, s=4`, `tau_max=3`, `30` rounds

The earlier banded auto-buffered version took:

- `7924 ms`

The bucketed version takes:

- `7336 ms`

The main totals stay aligned:

- completed rounds: `30`
- total exact leaves handled: `208144`
- total candidate keys recomputed: `287114`
- total spawned subleaves: `47770`
- max observed round min: `2`

### `web-Stanford, r=3, s=4`, `tau_max=7`, `400` rounds

The earlier banded auto-buffered version took:

- `42306 ms`

The bucketed version takes:

- `32974 ms`

This is a much larger improvement.

The main totals remain very close:

- completed rounds: `400`
- max observed round min: `7`
- total candidate keys recomputed: `2382815`

So the bucketed active layer is now the strongest version of the quotient line.

## Phase profile after bucketed support

We then profiled the current strongest bucketed version.

### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

The current default bucketed version gives:

- prototype time: `37452 ms`
- total rebuild time: `2699 ms`

The phase profile is:

- `tau=1`:
  - `phase_ms=11742`
  - `leaf_ms=1023`
  - `recomp_ms=7750`
  - `bucket_ms=50`
- `tau=3`:
  - `phase_ms=25380`
  - `leaf_ms=1105`
  - `recomp_ms=21836`
  - `bucket_ms=94`

So the new bottleneck is very clear.

It is not rebuild.
It is not split.
It is not bucket maintenance.

It is touched-key support recomputation.

### `web-Stanford, r=3, s=4`, `tau_max=7`, `400` rounds

The current default bucketed version gives:

- prototype time: `33138 ms`
- total rebuild time: `21263 ms`

The phase profile is:

- `tau=1`:
  - `phase_ms=7274`
  - `leaf_ms=349`
  - `recomp_ms=669`
  - `bucket_ms=13`
- `tau=3`:
  - `phase_ms=13546`
  - `leaf_ms=2099`
  - `recomp_ms=3904`
  - `bucket_ms=96`
- `tau=5`:
  - `phase_ms=6582`
  - `leaf_ms=245`
  - `recomp_ms=1263`
  - `bucket_ms=22`
- `tau=7`:
  - `phase_ms=5373`
  - `leaf_ms=94`
  - `recomp_ms=768`
  - `bucket_ms=8`

The same pattern appears here.
Recomputation remains the main non-rebuild cost.

## Negative and mixed experiments

We also tried two delta-based shortcuts.

### Delta fast path

This idea was:

- if a key is already inside the exact window
- update its new support by `old support + local delta`
- avoid full recomputation

This did not help.

On `com-youtube, tau_max=7, 400 rounds`, it made the run slower.
So we do not keep it on by default.

### Delta prune

This idea was:

- if a key is outside the exact window
- and its local delta is non-negative
- then it cannot move into the current window
- so we skip its full recomputation

This is a mixed result.

- on `com-youtube`, it helps
- on `web-Stanford`, it does not

So we keep it only as an experimental switch.

The current mainline conclusion is:

- the bucketed active layer is the strongest exact prototype so far
- the next real bottleneck is touched-key support recomputation

## Negative experiment: pivot-aware recompute index

We also tried a pivot-aware recompute index.

The idea was simple.
We added:

- `activePivotLeafByVertex`
- `activeKeepCount`
- `activePivotCount`

Then we replaced `leafContributionForKey(...)` with an indexed version.

This did not help.

On `web-Stanford, r=3, s=4, 30 rounds`, the time went from about `6575 ms` to about `6842 ms`.
So we did not keep this version.

The reason is also simple.
These leaves are not large enough.
Extra hash lookups are more expensive than a short leaf scan.

## New mainline: over-window exact support cache

The next idea was different.

We noticed that many touched keys are outside the current exact window.
We recomputed them from scratch.
Then we threw them away.
Later we touched them again.

So we added a small exact cache for touched keys whose support is still above the current window.

This cache is exact.
It is not an approximation.
If a cached key is touched again, we update it by:

- `new support = old cached support + local delta`

This avoids a full recomputation.

This idea works very well.

### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

Before the cache:

- prototype time: `37452 ms`
- total rebuild time: `2699 ms`
- `tau=1 recomp_ms=7750`
- `tau=3 recomp_ms=21836`

With the cache:

- prototype time: `18595 ms`
- total rebuild time: `2423 ms`
- `tau=1 recomp_ms=4498`
- `tau=3 recomp_ms=7138`
- `tau=1 over_cache_hit=346625`
- `tau=3 over_cache_hit=1329831`

This is about `2.01x` faster than the previous bucketed mainline.

### `web-Stanford, r=3, s=4`, `tau_max=7`, `400` rounds

Before the cache:

- prototype time: `33138 ms`
- total rebuild time: `21263 ms`

With the cache:

- prototype time: `28243 ms`
- total rebuild time: `19044 ms`
- `tau=1 recomp_ms=515`
- `tau=3 recomp_ms=2483`
- `tau=5 recomp_ms=606`
- `tau=7 recomp_ms=452`

This is about `1.17x` faster than the previous bucketed mainline.

### Short sanity checks

`web-Stanford, r=3, s=4, 30 rounds`:

- before: `6575 ms`
- with cache: `6625 ms` to `6670 ms`

This is roughly flat.
It does not break the short prefix.

`email-Eu-core, r=3, s=5, 200 rounds`:

- before: `840 ms`
- with cache: `635 ms`

This is about `1.32x` faster.

## Updated conclusion

The current strongest exact prototype is now:

- bucketed banded buffered adaptive low-support
- plus an exact over-window support cache

This is the first change after bucketed support that gives a large win on both:

- `com-youtube`
- `web-Stanford`

The next bottleneck is still touched-key recomputation.
But now the repeated over-window part is much smaller.

## New strongest version: use `localDelta` as the touched-key set

We then pushed one more exact change.

Before this step, the bucketed prototype still built a separate `candidateKeys` set.
That set came from explicit old/new leaf scans.

This was unnecessary.

For changed leaves, the exact touched set is already in `localDelta`.
So we changed the prototype to:

- use `localDelta` keys as the touched-key set
- skip the extra `candidateKeys` construction
- keep the exact over-window cache

This change is exact.
It does not change the algorithm output.
It only changes how we route touched keys.

### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

Before this step:

- prototype time: `18595 ms`
- total rebuild time: `2423 ms`
- `tau=1 recomp_ms=4498`
- `tau=3 recomp_ms=7138`

After this step:

- prototype time: `16565 ms`
- total rebuild time: `2356 ms`
- `tau=1 recomp_ms=4200`
- `tau=3 recomp_ms=6256`

This is about `1.12x` faster than the previous over-cache version.
It is about `2.26x` faster than the old pre-cache bucketed version (`37452 ms`).

### `web-Stanford, r=3, s=4`, `tau_max=7`, `400` rounds

Before this step:

- prototype time: `28243 ms`
- total rebuild time: `19044 ms`

After this step:

- prototype time: `25607 ms`
- total rebuild time: `19275 ms`
- `tau=1 recomp_ms=408`
- `tau=3 recomp_ms=1607`
- `tau=5 recomp_ms=385`
- `tau=7 recomp_ms=275`

This is about `1.10x` faster than the previous over-cache version.
It is about `1.29x` faster than the old pre-cache bucketed version (`33138 ms`).

### `email-Eu-core, r=3, s=5`, `tau_max=8`, `200` rounds

Before this step:

- prototype time: `635 ms`

After this step:

- prototype time: `592 ms`

So this change also helps on the small graph.

## Updated current status

The current strongest exact prototype is now:

- bucketed banded buffered adaptive low-support
- exact over-window support cache
- `localDelta`-driven touched-key routing

The bottleneck is now split by graph:

- on `com-youtube`, recomputation is still the main cost
- on `web-Stanford`, rebuild is now the main cost again

So the next real direction is no longer the same on every graph.
The next natural target is phase rebuild, not just touched-key recomputation.

## New phase policy: auto full-band when rebuild is the real bottleneck

We then pushed the phase policy one step higher.

The old auto rule only changed the lookahead width.
That was not enough on `web-Stanford`.

The new rule does two things:

- if the raw clique layer is already small, start with one full-band build
- otherwise start narrow, but promote the next phase to full-band if the current phase shows:
  - a very large `init_over / init_low` ratio
  - rebuild time clearly larger than the real update work

This does not change the exact algorithm.
It only changes the phase schedule.

### `email-Eu-core, r=3, s=5`, `tau_max=8`, `200` rounds

The old strongest version took:

- `592 ms`

The new auto full-band rule takes:

- `506 ms`

It now starts with one full-band phase:

- `Phase tau=1 cap=8 window=8 lookahead=8 full_band=yes`

So this is about `1.17x` faster than the previous strongest version.

### `web-Stanford, r=3, s=4`, `tau_max=7`, `400` rounds

The old strongest version took:

- `25607 ms`

The new auto full-band rule takes:

- `24055 ms`

The phase trace is exactly the one we wanted:

- first phase stays narrow:
  - `Phase tau=1 cap=2 window=2 lookahead=1 full_band=no`
- second phase is promoted to full-band:
  - `Phase tau=3 cap=7 window=7 lookahead=5 full_band=yes`

This cuts the number of rebuild phases from four to two.
It also cuts total rebuild time:

- `19275 ms -> 15298 ms`

So this is about `1.06x` faster than the previous strongest version.

### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

Here the new rule does not trigger full-band.
That is the right behavior.

The clean phase trace is:

- `Phase tau=1 cap=2 window=2 lookahead=1 full_band=no`
- `Phase tau=3 cap=4 window=4 lookahead=1 full_band=no`

We ran a clean control with full-band disabled.
That control took:

- `18549 ms`

We then reran the new auto rule.
That clean run took:

- `17573 ms`

So on this graph the new policy does not hurt.
It stays in the same timing band as the old strongest version.

## Updated current status

The current strongest exact prototype is now:

- bucketed banded buffered adaptive low-support
- exact over-window support cache
- `localDelta`-driven touched-key routing
- auto full-band phase promotion

This new phase rule solves a real weakness on `web-Stanford`.
It also improves the small graph.
It does not need to trigger on `com-youtube`.

So the next bottleneck is now even clearer:

- `web-Stanford`: phase rebuild is much less dominant than before
- `com-youtube`: recomputation is still the main cost

So the next real target should move back to `recomp_ms`.

## New exact change: use stored low/buffered support directly

We then returned to the real bottleneck on `com-youtube`.

The key observation is simple.

For keys that are already in the current low or buffered layer, we do know the exact old support.
So after one local delta, we should not recompute them from scratch.

We changed the prototype to do this:

- if a key is in the current low/buffered state, use `old_support + local_delta`
- if the new value stays in-window, insert it directly
- if the new value goes above the window, move it to the over-window cache
- only keys with no exact old support still go to full recompute

This change is exact.

### `email-Eu-core, r=3, s=5`, `tau_max=8`, `200` rounds

Before this step:

- `506 ms`

After this step:

- `446 ms`

This is about `1.13x` faster.

### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

Before this step:

- `17573 ms`

After this step:

- `17542 ms`

So the wall time is almost flat.
But the internal shape is better:

- `tau=1 full_recomp: 906884 -> 797584`
- `tau=3 full_recomp: 1071142 -> 889797`

We also now handle many keys directly from exact state:

- `tau=1 delta_state=109300`
- `tau=3 delta_state=181345`

So this step is real.
It just does not yet move the total time much on `com-youtube`.

### `web-Stanford, r=3, s=4`, `tau_max=7`, `400` rounds

Before this step:

- `24055 ms`

After this step:

- `18617 ms`

This is a much stronger gain.
It is about `1.29x` faster.

The phase profile also improves:

- `tau=1 full_recomp=226450`
- `tau=3 full_recomp=856053`
- `tau=1 delta_state=41191`
- `tau=3 delta_state=654039`

So this exact state-delta change is a clear positive result.

## Negative follow-up: one more exact band above the window

We also tested a stricter idea.

We added one more exact support band above the current window.
The goal was to reduce full recompute even more.

This idea does reduce `full_recomp`.
But it does not yet improve wall time.

### `com-youtube, r=3, s=4`, `over_band=1`

- prototype time: `18522 ms`
- `tau=1 full_recomp=634875`
- `tau=3 full_recomp=765287`

### `com-youtube, r=3, s=4`, `over_band=2`

- prototype time: `18002 ms`
- `tau=1 full_recomp=511108`
- `tau=3 full_recomp=660206`

Both settings reduce recomputation a lot.
But both are still slower than the current strongest run:

- strongest: `17542 ms`

So this new exact band is not the next mainline direction.
It is too expensive in phase rebuild and state maintenance.

## Updated current status

The current strongest exact prototype remains:

- bucketed banded buffered adaptive low-support
- exact over-window support cache
- `localDelta`-driven touched-key routing
- auto full-band phase promotion
- exact state-delta for keys already in low/buffered state

The new exact over-band is only an experiment.
It is not part of the strongest version.

So the next target is now sharper:

- `web-Stanford`: phase rebuild is already much better
- `com-youtube`: the remaining hard part is still first-time full recompute for keys with no exact old support

## New auto policy: medium graphs get a wider first window

We then revisited the auto rule.

The old rule only used `rawCliqueCount`.
That was too coarse.

It could help `com-youtube`.
But it could also widen the first phase on graphs that should stay conservative.

So we changed the rule again.

We now reuse the clean quotient survey.
The first phase widens only when:

- `rawCliqueCount` is not too large
- and clean `total explicit entries` is in a medium range

The default medium range is:

- `autoWideExplicitMin = 10000000`
- `autoWideExplicitMax = 50000000`

This keeps:

- `com-youtube` inside the widened regime
- `web-Stanford` outside it
- `email-Eu-core` outside it

## New single-run results

### `email-Eu-core, r=3, s=5`, `tau_max=8`, `200` rounds

The new run gives:

- prototype time: `417 ms`
- rebuild count: `1`
- total rebuild time: `87 ms`
- survey explicit entries: `3246772`

This graph stays on the small-graph path.

### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

The new run gives:

- prototype time: `14827 ms`
- rebuild count: `1`
- total rebuild time: `1707 ms`
- survey explicit entries: `14335854`
- phase: `tau=1 cap=3 window=4 lookahead=3 full_band=no`

This is better than the previous strongest result:

- old strongest: `17542 ms`

It is also better than the earlier fixed-wide run:

- fixed `lookahead=3`: `15359 ms`

So the new medium-graph rule is a real improvement on `com-youtube`.

### `web-Stanford, r=3, s=4`, `tau_max=7`, `400` rounds

The new run gives:

- prototype time: `18437 ms`
- rebuild count: `2`
- total rebuild time: `11669 ms`
- survey explicit entries: `115534074`
- phase 1: `tau=1 cap=2 window=2 lookahead=1 full_band=no`
- phase 2: `tau=3 cap=7 window=7 lookahead=5 full_band=yes`

This is also better than the previous strongest result:

- old strongest: `18617 ms`

So the new explicit-range rule keeps the `web-Stanford` behavior we wanted.

## Updated strongest prototype

The current strongest exact prototype is now:

- bucketed banded buffered adaptive low-support
- exact over-window support cache
- `localDelta`-driven touched-key routing
- auto full-band phase promotion
- exact state-delta for keys already in low/buffered state
- a medium-graph first-phase widening rule based on clean explicit-entry range

This is the first auto rule that improves all three representative graphs in clean single runs:

- `email-Eu-core`: `446 ms -> 417 ms`
- `com-youtube`: `17542 ms -> 14827 ms`
- `web-Stanford`: `18617 ms -> 18437 ms`

---

*The following results are from a new contributor.*

## Build-time over-cache seeding

We attacked the remaining `full_recomp` bottleneck on `com-youtube`.

### Problem

The previous over-window cache only stored keys that were touched during peeling rounds.
Keys touched for the first time had no cached support and required expensive full recomputation
via `intersect_dense_sets_multi` + `leafContributionForKey`.

On `com-youtube, r=3, s=4, 400 rounds`, the previous strongest version had:

- `tau=1 full_recomp=797584`
- `tau=3 full_recomp=889797`

Each full recompute does an `r`-way vertex intersection over all active leaves.
This dominated `recomp_ms`.

### Key observation

During `buildBandedLowSupportLayerActive`, the code already computes exact support
for every key.
But for keys exceeding the window, it places them in `overWindow` (a set) and SKIPS
further contributions from other leaves (line 2358-2361).
The exact support value is discarded.

Later, when those keys are touched by a delta, the code cannot use
`old_support + delta` because `old_support` was never stored.
So it falls through to a full recompute.

### Fix

We added a `keepExactOverSupport` flag to `buildBandedLowSupportLayerActive`.
When true, the build function keeps accumulating contributions for over-window keys
in `overLowerBound` instead of skipping them.
This makes `overLowerBound` exact (not just a lower bound).

After the build, we seed the runtime `overSupportCache` from these exact values.
All subsequent touches can then use `old_support + delta` instead of full recompute.

The build cost increases slightly because over-window keys are no longer skipped.
But the build already iterates all leaves anyway.
The extra cost is one hash-map update per skipped incidence.

### Results

#### `email-Eu-core, r=3, s=5`, `tau_max=8`, `200` rounds

- prototype time: `174 ms`
- full_recomp: `24`
- over_cache_hit: `86537`

Previous strongest: `417 ms`. This is `2.40x` faster.

#### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

- prototype time: `9919 ms`
- total rebuild time: `3767 ms`
- `tau=1 recomp_ms=1274, full_recomp=11870, over_cache_hit=1896370`
- `tau=4 recomp_ms=515, full_recomp=236, over_cache_hit=962466`

Previous strongest: `14827 ms`. This is `1.49x` faster.

The `full_recomp` dropped from `797584+889797 = 1.69M` to `11870+236 = 12K`.
That is a `140x` reduction in full recomputes.
The remaining 12K full recomputes are from keys generated by new sub-leaves
that did not exist during the build.

#### `com-dblp, r=3, s=4`, `tau_max=7`, `400` rounds (new graph)

- prototype time: `4947 ms`
- total rebuild time: `3944 ms`
- `full_recomp: 1929+148+19 = 2096`
- `over_cache_hit: 70081+36839+21060 = 127980`

### Why this works

The build already does the expensive work of iterating all leaves and computing all
key contributions.
The only additional cost is storing the result instead of discarding it.
This trades `O(overTauKeys)` extra hash-map memory during build for
eliminating `O(touchedKeys)` full recomputes during peeling.

On `com-youtube`, the build takes `~3.8s` (up from `~1.7s`).
But the peeling recompute drops from `~11s` to `~1.8s`.
Net gain: `~7s`.

### Remaining full recomputes

The remaining `~12K` full recomputes on `com-youtube` come from keys in
sub-leaves created during peeling.
These sub-leaves did not exist at build time, so their keys have no cached support.

Two options to reduce further:

1. When a sub-leaf is created, compute and cache support for all its new keys
2. Use `localDelta` as exact update for sub-leaf keys (they are already in the delta)

Option 2 is the more natural next step.

### Updated strongest prototype

The current strongest exact prototype is now:

- bucketed banded buffered adaptive low-support
- **build-time over-cache seeding** (new)
- exact over-window support cache
- `localDelta`-driven touched-key routing
- auto full-band phase promotion
- exact state-delta for keys already in low/buffered state
- medium-graph first-phase widening rule

Performance summary:

- `email-Eu-core, r=3, s=5, 200 rounds`: `417 ms -> 174 ms` (`2.40x`)
- `com-youtube, r=3, s=4, 400 rounds`: `14827 ms -> 9919 ms` (`1.49x`)

## Zero-support delta shortcut

We then eliminated the remaining full recomputes.

### Observation

The remaining `~12K` full recomputes are for keys that had zero support at build time.
These keys were filtered by `hasPositiveContribution` during the build, so no cache entry was created.

Key insight: if a key had zero support at build time and has never been touched since
(no cache/state entry), its support is still `0`.
So when it first appears in `localDelta`:

- `new_support = 0 + delta = delta`

No intersection-based recompute is needed.

### Why this is exact

1. Build processes ALL active leaves. Zero-support keys have no positive contribution from any leaf.
2. Previous rounds' leaf deaths subtract zero contributions → support stays at `0`.
3. The first time a zero-support key appears in `localDelta`, it is from a sub-leaf with different keep/pivot assignments that make the contribution positive.
4. No other alive leaf contributes positively (they were all processed at build time).

### Results

#### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

- prototype time: `9348 ms`
- `tau=1 recomp_ms=1204` (was `1274`)
- `tau=4 recomp_ms=514` (was `515`)

State counts unchanged (1,103,883 exact leaves, 203,916 spawned).

#### `email-Eu-core, r=3, s=5`, `tau_max=8`, `200` rounds

- prototype time: `196 ms`

State counts unchanged (9,757 exact leaves, 1,825 spawned).

### Updated performance summary

- `email-Eu-core, r=3, s=5, 200 rounds`: `417 ms -> 196 ms` (`2.13x`)
- `com-youtube, r=3, s=4, 400 rounds`: `14827 ms -> 9348 ms` (`1.59x`)

### Current bottleneck profile (`com-youtube, 400 rounds`)

| Component | Time (ms) | Fraction |
|---|---|---|
| Rebuild | 3,470 | 37% |
| Leaf update (BK split) | 1,390 | 15% |
| Recompute (delta + cache) | 1,718 | 18% |
| Bucket update | 86 | 1% |
| Framework overhead | 2,684 | 29% |

The rebuild cost is now the single largest component.
It increased from `~1.7s` to `~3.5s` due to exact over-cache accumulation during build.
The recompute cost dropped from `~11s` to `~1.7s`.

The next bottleneck is rebuild cost.
The build iterates all active leaves and accumulates all key contributions.
The `~8M` extra hash-map updates (previously skipped) add `~2s` to each rebuild.

## web-it-2004: build enumeration barrier

We attempted the prototype on `web-it-2004, r=3, s=4`.

It cannot run.

The build phase (`buildBandedLowSupportLayerActive`) enumerates `C(|L|, r)` r-cliques
per leaf. On this graph:

- 423,886 leaves, largest with 432 vertices (keepC=1, pivotC=431)
- `C(432, 3) = 13,343,760` cliques per dense leaf
- Total leaf-level incidences: `28,870,547,161` (28.8B)
- Even at 5ns per combination: `28.8B * 5ns = 144,000s = 40 hours`

The skip optimization does not help here.
Even though most keys get evicted to `overTau` quickly, the build function still
GENERATES each combination (to compute its hash key and check the map).
The combination generation itself is the bottleneck.

Meanwhile, the number of unique r-clique keys is only `338,884,212` (338M).
The quotient compression is `34,055x`.
And the first frontier has 100% unique-key ratio.

So the build enumerates 28.8B incidences to produce 338M unique keys,
of which most are immediately above `tau`.
This is an `85x` overhead from leaf-level enumeration.

### Why this matters

This is the graph where the quotient approach should have the biggest advantage.
V20/REF both timeout on this graph (1h).
The quotient compression is `34,055x`.
But the current prototype cannot even START because the build is leaf-centric.

### Next direction: vertex-centric low-support build

The fix is to enumerate unique r-clique keys from vertices, not from leaves.

Instead of:
```
for each leaf L:
    for each C(|L|, r) combination:
        compute key, look up in map
```

Do:
```
for each unique r-clique key K (338M total):
    find all leaves containing K via vertex-to-leaf index
    sum contributions to get support
```

This reduces the build from `O(28.8B)` to `O(338M)`.

For `tau=1`, this can be further optimized:
- A key with support=1 appears in exactly one leaf
- Short-circuit: if the vertex intersection has size > 1, skip immediately
- Most keys have support >> 1 on dense graphs, so most are quickly skipped

The vertex-to-leaf index (`activeLeafByVertex`) already exists in the prototype.
The remaining challenge is enumerating the 338M unique keys efficiently
without going through the leaf-level `C(|L|, r)` enumeration.

One approach: enumerate unique r-clique keys from the `treeGraphV` bipartite structure.
For `r=3`, iterate over vertex triples that share at least one common leaf.
Use 3-way `activeLeafByVertex` intersection to count shared leaves.

This is the next concrete algorithmic step.

## Two-phase build result

We implemented and tested a two-phase low-support build.

Phase A (leaf-centric): process the N densest leaves to build overTau set.
Phase B (key-centric): verify remaining lowSupport candidates via intersection.

### web-it-2004, r=3, s=4

Phase A completed in `35s` (50 leaves, 623M enumerated out of 28.8B).
This is a `46x` reduction over full enumeration.

However: `lowKeys=0` at `tau=1`. Also `tau=10`, `tau=100`: all `lowKeys=0`.

At `tau=430`: `lowKeys=52,094,867`, `overKeys=0`.

### Why the initial min support is ~429

Every vertex appears as KEEP in at least one leaf.
A triple `(v1,v2,v3)` gets contribution from `v1`'s keep-leaf:
- `subNumPivot=2` (only `v2,v3` are pivots)
- contribution = `nCr[pivotC-2][s-keepC-2]` = `nCr[~428][1]` = `~428`

Even the most peripheral triple gets `~428` support from its keep-leaf alone.
Adding the pivot-leaf contribution `1` gives `~429`.

### What this means

The `low-support frontier first` strategy degrades on web-it-2004.

At `tau ≈ 429`, the low-support layer contains `52M` out of `338M` unique keys.
This is `15%` of all keys — too large for an efficient frontier layer.

So the bottleneck has shifted again.
The issue is no longer build speed (two-phase handles it).
The issue is that the initial min support is too high for the frontier approach to help.

### Structural root cause

On dense clique-rich graphs:
- every vertex is a keep in its own leaf
- this gives every r-clique a large "base support" from the keep-leaf contribution
- the contribution formula `nCr[pivotC - subNumPivot][s - keepC - subNumPivot]` amplifies this

So the support values are clustered near `omega - r + 1`, not near `1`.
The low-support frontier approach works best when supports are spread across small values.

### Next direction

The quotient approach needs to move beyond individual-key support tracking.

Two possible paths:

1. **Class-level peeling with support factorization**:
   decompose support into per-leaf class contributions + inter-leaf overlap corrections.
   Track corrections compactly instead of per-key support.

2. **Direct quotient peeling**:
   peel at the quotient class level. For a clean leaf, all keys in a class die together
   (they have the same per-leaf contribution from that leaf). The question is whether
   we can determine the TOTAL support ordering at the class level.

## Vertex-level support formula (r=3)

We derived and verified an exact support formula that computes from graph structure.

For `r=3, s=4`, a triangle `K=(va,vb,vc)` where `va` is the earliest vertex
in the degeneracy ordering:

```
support(K) = (pivotC_va - 2) + n_commonEarlierNeighbors(K)
```

where `n_commonEarlierNeighbors` = number of vertices `u` with earlier degeneracy
position that are adjacent to all three.

### Verification

| Graph | Smallest pivotC | Estimated min support | Analysis time |
|---|---|---|---|
| email-Eu-core | 3 | 1 | 0 ms |
| web-it-2004 | 3 | 1 | 11 ms |

This correctly identifies `min support = 1` from `pivotC=3` leaves, in milliseconds.
No leaf enumeration or hash map needed.

### Vertex-level peeling prototype

We implemented a peeling prototype that:
1. Computes support via the formula (no 338M enumeration)
2. Uses lazy level-by-level enumeration (only from small leaves first)
3. Does correct s-clique-based support updates after peeling

Results on `web-it-2004, r=3, s=4` (where V20 timeouts at 1h):
- 200 rounds, 572 triangles peeled, **195 ms**
- No CliqueIndex needed (0 extra memory vs V20's 34GB)

Results on `email-Eu-core, r=3, s=4`:
- 537 rounds, 6179 triangles peeled, **56 ms**
- Correct cascade behavior (support updates create new low-support triangles)

### Current limitations

1. The bucket PQ needs better deduplication (core values are not yet monotonic)
2. Only ~6% of total triangles are peeled (lazy enumeration needs to expand further)
3. Support update correctness needs verification against V20's core distribution

### Progress on correctness

We fixed multiple bugs in the peeling prototype:

1. **S-clique ownership** (fixed): s-cliques from leaf L must include L's keep vertex.
   When K is all-pivot in L, the only valid extra vertex is the keep vertex.
   Earlier code allowed all vertices as extra → massive over-decrement.

2. **Already-destroyed s-clique check** (fixed): when peeling K, s-cliques
   already destroyed by earlier peelings should not cause decrements.
   Added check: if any other triangle in {K∪d} is already peeled, skip.

3. **Proper drain** (fixed): standard peeling requires processing ALL items
   with support ≤ current core level (not just one bucket).
   Items that cascade to lower support join the current frontier.

After these fixes:
- `core=1 peeled=1853` — **exact match with V20** ✓
- But `core=2 peeled=102877` — should be 3082. Over-decrement persists.

### Current diagnosis

The remaining issue: within a single core-level round, 1853 peeled triangles
share s-cliques with many neighbors. Each peeled triangle correctly decrements
its neighbors. But the CUMULATIVE effect across 1853 peelings reduces most
triangles' support to ≤ 2.

The bug is likely in how updates interact within a batch: when multiple peeled
triangles share the SAME s-clique, the decrements may double-count.

Example: s-clique S={a,b,c,d} contains K1={a,b,c} and K2={a,b,d}.
If both K1 and K2 are peeled in the same round:
- K1 peeled → N={a,c,d} gets -1 from S
- K2 peeled → N={a,c,d} should NOT get another -1 from S (already destroyed by K1)
- The `isAlreadyDestroyed` check handles this IF K1 is processed before K2

But the check looks for peeled neighbors of K2 in S, not K1 specifically.
For S={a,b,c,d}: K2={a,b,d}, check replaces each of K2's vertices with d=c:
- (c,b,d) → not K1. (a,c,d)=N → not peeled. (a,b,c)=K1 → IS peeled!
So isAlreadyDestroyed returns true. S is correctly skipped for K2. ✓

The check appears correct. The over-decrement must come from a different source.

### Validation on dblp-core30 (686K triangles)

The prototype now matches V20 on several low-core layers:

| Core | V20 | Prototype | Match |
|---|---|---|---|
| 1 | 15 | 15 | ✅ |
| 2 | 20 | 20 | ✅ |
| 3 | 53 | 53 | ✅ |
| 4 | 48 | 48 | ✅ |
| 5 | 42 | 42 | ✅ |
| 7 | 36 | 36 | ✅ |
| 19 | 420 | split across 16-19 | ❌ |
| 24-25 | 196 | split across 24-35+ | ❌ |

Time: 4.9s for the current complete prototype run (686K triangles, 34 rounds).

The low-core prefix is exact on the listed layers.
The full decomposition is still inexact.
The remaining triangles have visible core offsets.
The total enumerated triangle count is still 686,182.

Our current hypothesis is still the same.
Within a large batch, some s-clique decrements may still double-count.
The `isAlreadyDestroyed` check handles many cases.
It does not yet close the full gap.

### Key milestone

The vertex-level peeling approach now has three clear properties:
1. It is exact on a verified low-core prefix.
2. It runs without `StaticCliqueIndex`.
3. It can finish a full prototype run on `dblp-core30` in 4.9s.

It is still not an exact full decomposition algorithm.

### Remaining work

1. Fix the remaining within-batch s-clique double-decrement
2. Implement lazy enumeration for web-it-2004 (full enumeration is infeasible at 28.8B)
3. Combine vertex-level formula with lazy expansion for large graphs

## Class-level peeling breakthrough

We implemented a class-level prototype for the `p=2` phase.
This is not yet a full class-level decomposition.

### Key idea

Instead of tracking 338M individual r-cliques, track ~775K quotient classes.
Each class = (leafId, subNumPivot p). For r=3, s=4 with keepC=1 leaves:
- Class p=2: C(pivotC, 2) triangles, each with contribution pivotC-2
- Class p=3: C(pivotC, 3) triangles, each with contribution 1

### Phase 1: p=2 class peeling only

Sort leaves by contribP2 ascending. Peel each leaf's entire p=2 class at once.
No individual triangle enumeration needed.

Results on **web-it-2004, r=3, s=4** for Phase 1 only:
- 387,650 leaves processed
- 338,476,979 triangles peeled (p=2 class)
- Max p=2 class core: 429
- **Time: 12 ms**

This number is only for the `p=2` phase.
It is not yet comparable to a full V20 decomposition run.

### Cascade rule

Peeling leaf L's p=2 class destroys each s-clique {keep, pi, pj, pk} from L.
This causes each p=3 triple (pi,pj,pk) from L to lose exactly 1 support.
This is computed in O(1) per class, not O(C(pivotC,3)).

### Phase 2: p=3 class (pending)

After all `p=2` classes are peeled, a very large `p=3` state still remains.
The current prototype only reports this next step.
It does not yet peel these `p=3` classes exactly.

Our target formula is:

```
support_p3(T) = Σ_{vi ∈ T} contribP2(vi's keep-leaf) + n_allPivotLeaves(T)
              - n_leaves_with_p2_already_peeled(T)
```

This is the next planned direction.
It is not implemented end-to-end yet.

## Current clean baseline

Current clean single-run results for the strongest default path:

- `email-Eu-core, r=3, s=5, 200 rounds`: `186 ms`
- `com-youtube, r=3, s=4, 400 rounds`: `10274 ms`
- `web-Stanford, r=3, s=4, 400 rounds`: `18704 ms`

These runs use the current default:
- bucketed banded buffered adaptive low-support
- exact over-window support cache
- `localDelta` touched-key routing
- auto full-band phase promotion
- automatic `over-lower` routing when cache-reband is on

For `com-youtube`, the current default still uses a wide first phase:
- `tau=1 cap=4 window=5 lookahead=4`
- `rebuild count = 1`

For `web-Stanford`, the current default still works best with:
- `phase 1: tau=1 cap=2 window=2 lookahead=1`
- `phase 2: tau=3 cap=7 window=7 lookahead=5 full_band=yes`

## Negative result: auto two-phase first build

We also tested an experimental `two-phase` first rebuild on the strongest path.
It is now kept behind an explicit env flag:

- `PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_TWO_PHASE=1`

This experiment is not on by default.

On `web-Stanford, r=3, s=4, 400 rounds`, it gave:
- first rebuild: `764 ms`
- but phase 1 produced `lowKeys=0`
- total time: `40904 ms`

So this idea reduces the first rebuild itself.
But it destroys too much low-support state and hurts the whole run.
It is not suitable as the default path.

## Phase-level cache reband

We then fixed a real bug in the earlier `cache-reband` idea.

The earlier version tried to rebuild the next phase from the old bucketed state.
That could not work well, because the phase loop clears `state` after each phase.
So the old `cache-reband` path was often inert.

We changed the mechanism.
Now the next phase can rebuild directly from the exact `overSupportCache`.
This is only used when `overExactBand = 0`.
In that case, the cache already holds the exact support of all keys above the previous window.

This is enabled by:

- `PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND=1`

### `web-Stanford, r=3, s=4`, `tau_max=7`, `400` rounds

This is the main positive result.

Without phase-level cache reband, a clean run gave about:

- prototype time: `34733 ms`
- total rebuild time: `25882 ms`
- phase 2 rebuild: `12904 ms`

With the fixed cache-reband path, a clean run gave:

- prototype time: `17145 ms`
- total rebuild time: `10903 ms`
- phase 2 rebuild: `1047 ms`
- phase 2: `cache_reband=yes`

So this is not a small constant-factor tweak.
It removes most of the second rebuild cost.

### `com-youtube, r=3, s=4`, `tau_max=7`, `400` rounds

This graph still finishes in one wide first phase:

- prototype time: `11394 ms`
- rebuild count: `1`
- phase 1: `tau=1 cap=4 window=5 lookahead=4 full_band=no cache_reband=no`

So the new cache-reband path does not help here.
But it also does not disturb the structure of the run.

### `email-Eu-core, r=3, s=5`, `tau_max=8`, `200` rounds

This graph still stays on the small-graph full-band path:

- prototype time: `421 ms`
- rebuild count: `1`
- phase 1: `tau=1 cap=8 window=8 lookahead=8 full_band=yes cache_reband=no`

Again, the new path does not trigger here.

## Negative result: heavy first-phase auto widening

After the cache-reband fix, we also tried a new heuristic:

- if clean compression is very high
- and leaves are very large
- and cache-reband is on
- then widen the first phase to `lookahead=3`

This helps on `web-Stanford` when forced by hand.
But as an automatic rule it changes the later phase schedule in a bad way.
So it is not part of the default path.

It is now behind an explicit env flag:

- `PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_CACHE_REBAND_HEAVY=1`

Current conclusion:

- `cache-reband` is a real step forward
- it materially improves the hard multi-phase case
- but the next weak point is still the first rebuild on `web-Stanford`
- and the heavy first-phase widening rule is not stable enough yet

## New default: cache-reband + over-lower

After the phase-level cache-reband fix, we re-tested the old `over-lower` idea.

This time it became useful.
The reason is simple.
Now the hard multi-phase path keeps a much better exact cache across phases.
So the `over-lower` filter can prune many touched keys before they fall back to full recompute.

We therefore changed the strongest path:

- when `PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND=1` is on
- `over-lower` is now enabled by default

If we need the old behavior for comparison, we can turn it off with:

- `PIVOTER_QUOTIENT_BUCKETED_DELTA_OVER_LOWER_OFF=1`

### Clean results with the new default

- `email-Eu-core, r=3, s=5, 200 rounds`: `186 ms`
- `com-youtube, r=3, s=4, 400 rounds`: `10274 ms`
- `web-Stanford, r=3, s=4, 400 rounds`: `18704 ms`

These are better than the previous cache-reband baseline:

- `email`: `421 ms -> 186 ms`
- `com-youtube`: `11394 ms -> 10274 ms`
- `web-Stanford`: `21698 ms -> 18704 ms`

On `web-Stanford`, the phase summary now shows:

- phase 1: `delta_over_lb = 151`
- phase 2: `delta_over_lb = 5538`
- phase 2 `full_recomp = 57`

On `com-youtube`, the phase summary shows:

- phase 1: `delta_over_lb = 27216`
- phase 1 `full_recomp = 1127`

So the current strongest result is not coming from a new local formula.
It comes from a better interaction between:

- phase-level cache reband
- over-window lower bounds
- exact sparse delta updates

## Negative results after this step

We also tested two more ideas.

### Limited exact cache for cache-reband

We tried:

- `PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND_LIMITED=1`

The idea was to keep exact cache only up to `tauMax`.
Keys far above `tauMax` would only keep an `overWindow` mark.

This did not work well.
Without extra help, it caused too many later recomputations.
Even with `over-lower`, it was not better than the new default path.

So this idea is kept only as an experiment.
It is not part of the default path.

### Larger reserve for build maps

We also tested:

- `PIVOTER_QUOTIENT_BUILD_LARGE_RESERVE=1`

This gave only a tiny change on `web-Stanford`.
It did not materially reduce the first rebuild.

So this is also not part of the default path.

## More negative results on the first rebuild

We kept pushing on the remaining weak point.
That is still the first rebuild on `web-Stanford`.

### Auto limited exact cache for medium graphs

We tried to make `limited exact cache` more selective.
The idea was:

- keep the current default on very large graphs
- but on medium graphs, only keep exact over-window support up to `tauMax`
- use `over-lower` for the rest

This was implemented as an experiment:

- `PIVOTER_QUOTIENT_BUFFERED_CACHE_REBAND_LIMITED_AUTO=1`

We also fixed the implementation so that:

- exact seed cache and lower-bound cache are truly split
- this is no longer the earlier fake version where both paths shared one map

The result is still negative.

On `com-youtube, r=3, s=4, 400 rounds`, it gave:

- prototype time: `10656 ms`
- first phase `full_recomp = 1514879`
- second phase appeared as `tau=6 -> 7`

So this variant changed the phase structure in a bad way.
It did not beat the current default.

On `web-Stanford, r=3, s=4, 400 rounds`, it also did not help:

- prototype time: `20337 ms`

So this idea is not the next main path.

### Sorting leaves by size before the build

We also tested a simple structural idea:

- process larger leaves first
- let hot keys enter `overWindow` earlier

This experiment is behind:

- `PIVOTER_QUOTIENT_BUILD_SORT_LEAVES=1`

It was also negative.

On `web-Stanford, r=3, s=4, 400 rounds`, it gave:

- prototype time: `19379 ms`

This is slower than the current default `18704 ms`.

On `com-youtube, r=3, s=4, 400 rounds`, it gave:

- prototype time: `12029 ms`

This is also slower than the current default `10274 ms`.

So leaf ordering is not the right explanation for the first rebuild cost.

## Current research status

At this point, the strongest path is still:

- bucketed banded buffered adaptive low-support
- phase-level cache reband
- automatic `over-lower`

The current clean single-run results are:

- `email-Eu-core, r=3, s=5, 200 rounds`: `186 ms`
- `com-youtube, r=3, s=4, 400 rounds`: `10274 ms`
- `web-Stanford, r=3, s=4, 400 rounds`: `18704 ms`

The main remaining weak point is still the first rebuild on `web-Stanford`.
But the recent experiments show two things clearly:

- it is not just a cache-size problem
- it is not solved by changing leaf order

So the next useful step is likely not another scheduling tweak.
It is more likely a different first-build representation.

## Single-over-map heavy mode

We then tested a new first-build representation for heavy graphs.

The idea is simple:

- do not keep both `overWindow` and an exact over-support map during the first narrow phase
- use one exact over-support map instead
- keep the first phase narrow
- still let phase 2 go to `full-band + cache-reband`

This first worked only as a manual experiment:

- `PIVOTER_QUOTIENT_BUILD_SINGLE_OVER_MAP=1`

By itself, this reduced the first rebuild a lot on `web-Stanford`.
But it also hid the true over-state size from the phase scheduler.
So phase 2 no longer promoted to full-band.
That made the total run worse.

We fixed this by changing the phase policy.
It now uses tracked over-state size, not only `overWindow.size()`.

After that, the result became much better on `web-Stanford`.
But the best version was not `single-over-map` alone.
It was:

- `single-over-map`
- plus `capped-low`
- only on the heavy first narrow phase

We wrapped this as a heavy-mode experiment:

- `PIVOTER_QUOTIENT_BUILD_SINGLE_OVER_MAP_HEAVY=1`

This mode is meant for:

- large leaf
- high compression
- first narrow phase
- cache-reband path

In our local runs, it gave:

- `web-Stanford, r=3, s=4, 400 rounds`: `16688 ms`
- `com-youtube, r=3, s=4, 400 rounds`: `11144 ms`

The important point is not only the total time.
The first rebuild on `web-Stanford` dropped to about `8.2 s`.
That is much lower than the old first rebuild.

So this is the first first-build representation change that looks truly useful.

## Updated status

The main line is now split into two cases:

- normal graphs: keep the old strongest path
- heavy compressed graphs: try `single-over-map heavy`

The current strongest ideas are:

- bucketed banded buffered adaptive low-support
- phase-level cache reband
- automatic `over-lower`
- heavy first-build mode: `single-over-map + capped-low`

The remaining weak point is still the first rebuild.
But now we finally have one representation change that helps on that point.

## Heavy exact-over reserve

After the heavy first-build mode worked, we looked one layer deeper.

The first narrow phase on `web-Stanford` still had:

- `single-over-map = yes`
- `capped-low = yes`
- about `10.8M` tracked over keys

But the reserve policy for that exact over-support map was still too small.
So the map had to grow many times during the first build.

We changed this in a very narrow way:

- only for the heavy `single-over-map` path
- reserve the exact over-support map much more aggressively

This was a real positive result.

On `web-Stanford, r=3, s=4, 400 rounds`, the new repeated runs gave:

- `11.93 s`
- `11.62 s`

The phase summary is also very clear:

- phase 1 still uses `single_over_map=yes` and `capped_low=yes`
- phase 1 rebuild is now about `5.3-5.5 s`
- phase 2 still uses `full-band + cache-reband`

So the current picture is much clearer now.

The useful heavy path is:

- narrow first phase
- `single-over-map`
- `capped-low`
- large reserve for the heavy exact over-support map
- then `full-band + cache-reband` for the next phase

This is the first change that really cuts the hardest first rebuild on `web-Stanford`.

## Current clean default after the reserve fix

With the heavy reserve fix in place, the current default path gives:

- `email-Eu-core, r=3, s=5, 200 rounds`: `439 ms`
- `com-youtube, r=3, s=4, 400 rounds`: `8250 ms`
- `web-Stanford, r=3, s=4, 400 rounds`: `11624 ms`

The phase summaries also show the intended split:

- `email-Eu-core`: `single_over_map=no`, `capped_low=no`
- `com-youtube`: `single_over_map=no`, `capped_low=no`
- `web-Stanford` phase 1: `single_over_map=yes`, `capped_low=yes`

So the heavy first-build path is now really graph-selective.

## `web-it-2004` is still blocked at phase 1 build

We then brought `web-it-2004` back into the analysis.

The survey confirms that it is the most extreme heavy-build case:

- active leaves: `423886`
- max leaf size: `432`
- total explicit entries: `28870547161`
- total clean quotient compression: `34055.01x`
- total one-step refined compression: `11495.15x`
- total one-removed delta compression: `2.64x`

So `web-it-2004` is not a borderline case.
It is exactly the kind of graph that should need a heavy build path.

We then tried two prefix runs on `web-it-2004, r=3, s=4`:

- current strongest path, forced run, `tau_max=4`, `80` rounds
- an even more aggressive `over-set-only` heavy build, forced run, `tau_max=1`, `20` rounds

Neither run produced a phase summary even after a long wait.
So the blocker is now very clear:

- on `web-it-2004`, the current prototype still gets stuck inside phase 1 build itself
- the problem is no longer phase scheduling
- the problem is no longer cache reuse
- the problem is no longer first-round local delta

This means the next step for `web-it-2004` is not another cache tweak.
It needs build-level class aggregation.

## Exact first-frontier ownership

We then added an exact ownership analysis for the first frontier.

Enable:

- `PIVOTER_QUOTIENT_FRONTIER_OWNERSHIP=1`

This analysis builds the exact first frontier on graphs where the sparse-support map is still manageable.
It then asks a simple question:

- for each exact frontier key,
- which leaf class is its witness,
- and how much of that witness class is fully inside the exact frontier.

The results are very clear.

### email-Eu-core, r=3, s=5

- first frontier keys: `1981`
- exact one witness: `1981`
- only `p=2, keepC=1`: `204`
- only `p=3`: `3`
- `keepC>1` witness: `1774`
- frontier witness classes: `1242`
- full-covered classes: `719`
- frontier in full classes: `1283`
- frontier in partial classes: `698`
- class coverage median/p90: `1.0 / 1.0`

### com-youtube, r=3, s=4

- first frontier keys: `428986`
- exact one witness: `428986`
- only `p=2, keepC=1`: `106029`
- only `p=3`: `6349`
- `keepC>1` witness: `316608`
- frontier witness classes: `289568`
- full-covered classes: `176807`
- frontier in full classes: `301324`
- frontier in partial classes: `127662`
- class coverage median/p90: `1.0 / 1.0`

### web-Stanford, r=3, s=4

- first frontier keys: `124970`
- exact one witness: `124970`
- only `p=2, keepC=1`: `50185`
- only `p=3`: `1349`
- `keepC>1` witness: `73436`
- frontier witness classes: `79670`
- full-covered classes: `39035`
- frontier in full classes: `71189`
- frontier in partial classes: `53781`
- class coverage median/p90: `0.6667 / 1.0`

These results change the picture.

First, the exact first frontier is not shared across leaves on these graphs.
Every exact frontier key has exactly one witness leaf.

Second, the current `vertex-level` idea is too narrow.
Most exact first-frontier keys do not come from `p=2, keepC=1`.
Most of them come from `keepC>1`.

The exact frontier category breakdown makes this clearer.

On `email-Eu-core, r=3, s=5`:

- `keepC=1, q=2`: `204`
- `keepC=1, q=3`: `3`
- `keepC=2, q=1`: `926`
- `keepC=2, q=2`: `136`
- `keepC=3, q=0`: `471`
- `keepC=3, q=1`: `189`

On `com-youtube, r=3, s=4`:

- `keepC=1, q=2`: `106029`
- `keepC=1, q=3`: `6349`
- `keepC=2, q=1`: `202585`
- `keepC=2, q=2`: `20865`
- `keepC=3, q=0`: `70739`
- `keepC=3, q=1`: `18481`

On `web-Stanford, r=3, s=4`:

- `keepC=1, q=2`: `50185`
- `keepC=1, q=3`: `1349`
- `keepC=2, q=1`: `56193`
- `keepC=2, q=2`: `3264`
- `keepC=3, q=0`: `11407`
- `keepC=3, q=1`: `2306`

So the exact frontier is dominated by:

- `keepC=2, q=1`
- `keepC=3, q=0`
- then `keepC=1, q=2`

It is not dominated by `keepC=1, q=3`.

Third, many witness classes are already fully frontier.
On these graphs, full-covered classes already contain about:

- `64.8%` of the frontier on `email-Eu-core`
- `70.2%` of the frontier on `com-youtube`
- `57.0%` of the frontier on `web-Stanford`

The `full-covered` and `partial` split is also structured.

On `email-Eu-core, r=3, s=5`:

- full-covered:
  - `keepC=2, q=1`: `666`
  - `keepC=3, q=0`: `471`
  - `keepC=1, q=2`: `144`
- partial:
  - `keepC=2, q=1`: `260`
  - `keepC=2, q=2`: `136`
  - `keepC=3, q=1`: `189`

On `com-youtube, r=3, s=4`:

- full-covered:
  - `keepC=2, q=1`: `147720`
  - `keepC=3, q=0`: `70739`
  - `keepC=1, q=2`: `74256`
- partial:
  - `keepC=2, q=1`: `54865`
  - `keepC=1, q=2`: `31773`
  - `keepC=2, q=2`: `18743`
  - `keepC=3, q=1`: `18343`

On `web-Stanford, r=3, s=4`:

- full-covered:
  - `keepC=2, q=1`: `40690`
  - `keepC=1, q=2`: `17613`
  - `keepC=3, q=0`: `11407`
- partial:
  - `keepC=1, q=2`: `32572`
  - `keepC=2, q=1`: `15503`
  - `keepC=2, q=2`: `3146`
  - `keepC=3, q=1`: `2294`

This is a useful next-step signal.

- `keepC=2, q=1` is the best first target.
- `keepC=3, q=0` is also very good.
- `keepC=1, q=2` is mixed.
- `keepC=2, q=2` and `keepC=3, q=1` are mostly refinement work.

We also compared this with the fast class-level phase-1 lower bound.
That lower bound has a very different shape.

For `web-it-2004, r=3, s=4`, phase-1 class mass is almost entirely:

- `keepC=1, q=3`: `28.53B` incidences

while the other categories are tiny by comparison.

For `web-Stanford`, phase-1 lower-bound mass is also dominated by:

- `keepC=1, q=3`
- `keepC=2, q=2`

For `com-youtube`, the largest phase-1 lower-bound groups are:

- `keepC=3, q=1`
- `keepC=2, q=2`
- `keepC=4, q=0`

So the fast class lower bound is not yet close to the exact frontier structure.
This is especially true on `web-it-2004`.

So the next exact direction is clearer now.

- Phase 1 should be class-first, not key-first.
- It must support `keepC>1`, not only `keepC=1`.
- It should peel full-covered classes directly.
- It should refine only the remaining partial classes.

## Partial-class gap analysis

We then measured the partial classes more carefully.

Enable:

- `PIVOTER_QUOTIENT_PARTIAL_CLASS_ANALYSIS=1`

For each partial class, we enumerate all members of that class.
We then split the non-frontier members by gap:

- `plus1`: support = frontier + 1
- `plus2`: support = frontier + 2
- `plus3p`: support >= frontier + 3

The results are useful.

### email-Eu-core, r=3, s=5

- `keepC=2, q=1`
  - frontier `260`
  - non-frontier `451`
  - `plus1=310`, `plus2=54`, `plus3p=87`
- `keepC=1, q=2`
  - frontier `60`
  - non-frontier `108`
  - `plus1=70`, `plus2=14`, `plus3p=24`

### com-youtube, r=3, s=4

- `keepC=2, q=1`
  - frontier `54865`
  - non-frontier `54865`
  - `plus1=39448`, `plus2=10390`, `plus3p=5027`
- `keepC=1, q=2`
  - frontier `31773`
  - non-frontier `23535`
  - `plus1=17771`, `plus2=3900`, `plus3p=1864`
- `keepC=2, q=2`
  - frontier `18743`
  - non-frontier `18743`
  - `plus1=3634`, `plus2=2651`, `plus3p=12458`
- `keepC=3, q=1`
  - frontier `18343`
  - non-frontier `33143`
  - `plus1=7337`, `plus2=5756`, `plus3p=20050`

### web-Stanford, r=3, s=4

- `keepC=2, q=1`
  - frontier `15503`
  - non-frontier `15503`
  - `plus1=13296`, `plus2=1176`, `plus3p=1031`
- `keepC=1, q=2`
  - frontier `32572`
  - non-frontier `25886`
  - `plus1=22923`, `plus2=1133`, `plus3p=1830`
- `keepC=2, q=2`
  - frontier `3146`
  - non-frontier `3146`
  - `plus1=452`, `plus2=331`, `plus3p=2363`
- `keepC=3, q=1`
  - frontier `2294`
  - non-frontier `4420`
  - `plus1=668`, `plus2=369`, `plus3p=3383`

This points to a cleaner refinement order.

- `keepC=2, q=1` is the best first refinement target.
- `keepC=1, q=2` is also promising.
- `keepC=2, q=2` and `keepC=3, q=1` are much less clean.
- So the next prototype should not try to refine every partial class at once.
- It should start with:
  - full-covered classes
  - then a one-step refinement for partial `keepC=2, q=1`
  - then partial `keepC=1, q=2`
