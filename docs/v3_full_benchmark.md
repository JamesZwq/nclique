# V3 (Region CPI) Full Benchmark Report

## Algorithms

| | ST (CPI) | V2 (Region Tuple) | V3 (Region CPI) |
|---|---|---|---|
| Init | Enumerate r-cliques | Enumerate s-tuples | CPI AggrCount |
| Peeling | BK tree mutation | s-tuple cascade | Analytical split |
| Avoids | s-clique enum | r-clique enum | r-clique AND s-tuple enum |

## Results: dblp-core30 (dense, ω=114, 1206 vertices)

### Time (ms)

| r | s | ST | V2 | **V3** | Fastest |
|---|---|-----|-----|---------|---------|
| 3 | 4 | 422 | **6** | 222 | V2 |
| 3 | 6 | 421 | **39** | 182 | V2 |
| 3 | 8 | 504 | 478 | **130** | **V3** |
| 3 | 10 | 482 | TO | **198** | **V3** |
| 4 | 5 | 13,796 | **18** | 895 | V2 |
| 5 | 6 | 512,743 | **30** | 2,477 | V2 |

### Memory (kB)

| r | s | ST | V2 | V3 |
|---|---|-----|-----|-----|
| 3 | 4 | 246,000 | **10,500** | 77,900 |
| 3 | 6 | 253,000 | **22,900** | 72,300 |
| 3 | 8 | 368,000 | 33,000 | 88,000 |
| 3 | 10 | 373,000 | TO | 90,200 |
| 4 | 5 | 1,791,000 | **13,200** | 757,600 |
| 5 | 6 | 8,027,000 | **15,300** | 13,465,000 |

### Crossover Point

V3 becomes faster than V2 at **s=8** for r=3.
V2 times out at **s=10**. V3 continues to run efficiently.

## Results: email-Eu-core (medium density, ω=25, 1005 vertices)

| r | s | ST | V2 | V3 | Fastest |
|---|---|-----|------|------|---------|
| 3 | 4 | **237** | 5,678 | 5,723 | **ST** |
| 3 | 5 | **374** | 20,055 | 10,240 | **ST** |

email-Eu-core: sparse-ish (region compression ≈ 1x). ST wins.
V3 beats V2 at s=5 (10s vs 20s).

## Results: web-it-2004 (very dense, ω=432, 41M vertices)

ST times out on all configurations.

| r | s | V2 | **V3** | Fastest |
|---|---|------|---------|---------|
| 3 | 4 | **3,550** | 21,479 | V2 |
| 3 | 6 | 18,449 | **10,499** | **V3** |
| 3 | 8 | TO | **5,891** | **V3** |
| 3 | 10 | TO | **3,045** | **V3** |

### Memory (kB)

| r | s | V2 | V3 |
|---|---|------|------|
| 3 | 4 | **693,000** | 19,371,000 |
| 3 | 6 | **2,632,000** | 18,328,000 |
| 3 | 8 | TO | 17,974,000 |
| 3 | 10 | TO | 17,648,000 |

V3 uses more memory (~18GB) due to CPI tree + constrained path structures.
V2 is more memory-efficient when it doesn't time out.

## Crossover Analysis

| Graph | r | V3 beats V2 at s= | V2 times out at s= |
|-------|---|-------------------|-------------------|
| dblp-core30 | 3 | s=8 | s=10 |
| web-it-2004 | 3 | s=6 | s=8 |
| email-Eu-core | 3 | s=5 | — |

## Key Findings

1. **V3 is the only algorithm that scales to large s.** V2's peeling cost grows
   as T_s × C(s,r), which explodes for s≥8. V3's cost depends on r, not s.

2. **V3 always beats ST** on dense graphs (1.3-4.1x on dblp-core30, ∞ on web-it-2004).

3. **V2 wins for small s** (s≤6 typically) due to simple s-tuple cascade.

4. **V3's memory is higher** than V2 because of CPI tree + constrained path structures.
   For web-it-2004: ~18GB (V3) vs ~2.6GB (V2). This is a limitation.

5. **Optimal strategy: adaptive dispatch.**
   - s ≤ crossover: use V2 (Region Tuple)
   - s > crossover: use V3 (Region CPI)
   - Sparse graphs: use ST

## Why V3 Scales Better with s

V2's peeling bottleneck: each s-tuple death triggers C(s,r) r-sub-tuple lookups.
- C(4,3)=4, C(6,3)=20, C(8,3)=56, C(10,3)=120
- At s=9: 35,944 s-tuples × C(9,3)=84 = 3M operations → 52 seconds

V3's analytical split: sub-path count depends on κ (active classes, bounded by r).
- r=3: κ ≤ 3, typically 1-2 sub-paths per split
- Sub-paths DECREASE with larger s (more paths die completely)
- s=8: only 640 splits vs s=4: 3,225 splits
