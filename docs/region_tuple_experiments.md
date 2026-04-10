# Region Tuple Experiments

## Correctness: EXACT for s = r + 1

All tested graphs with s = r + 1: core distribution matches ST baseline
**line-by-line** (every core level, every count).

| Graph | r | s | Core Levels | Match |
|-------|---|---|-------------|-------|
| debug_small (8V) | 3 | 4 | 2 | EXACT |
| dblp-core30 (1.2K V) | 3 | 4 | 29 | EXACT |
| dblp-core30 | 4 | 5 | 28 | EXACT |
| com-dblp (317K V) | 3 | 4 | 45 | EXACT |
| web-Stanford (281K V) | 3 | 4 | 58 | EXACT |

For s > r + 1: MISMATCH (cascade only handles s = r + 1 correctly).
- dblp-core30 r=3 s=5: minor diff at core=171

## Performance

| Graph | r | s | ST time | Tuple time | Tuples | Speedup |
|-------|---|---|---------|-----------|--------|---------|
| dblp-core30 | 3 | 4 | 380ms | 4ms | 1,706 | **95x** |
| dblp-core30 | 4 | 5 | 13.3s | 17ms | 4,434 | **784x** |
| dblp-core30 | 5 | 6 | 566s | 33ms | 7,211 | **17,000x** |
| com-dblp | 3 | 4 | 809ms | 3,613ms | 934K | 0.2x (slower) |
| web-Stanford | 3 | 4 | 17.7s | 211s | 6.3M | 0.08x (slower) |
| web-it-2004 | 3 | 4 | >1h (TO) | 3.2s | 788K | >1000x |

## Key Observation: When Region Tuple Wins

Region Tuple is faster when: **few large regions → few tuples → massive compression**.

| Graph | Max clique | Regions | Avg region size | Tuples/r-cliques | Fast? |
|-------|-----------|---------|----------------|-----------------|-------|
| dblp-core30 | 114 | 104 | ~12 | 1.7K/686K = 0.2% | YES |
| web-it-2004 | 432 | 57K | ~9 | 788K/338M = 0.2% | YES |
| com-dblp | ~35 | ~many | ~1 | 934K/~900K ≈ 100% | NO |
| web-Stanford | ~70 | ~many | ~1 | 6.3M/~6M ≈ 100% | NO |

When most regions have size 1: tuple count ≈ r-clique count → no compression → slower
(due to hash map overhead vs ST's array-based approach).

## The Cascade Bug (Fixed)

Original bug: `currentLevel` could decrease during cascade, causing core values
to be assigned at lower levels than correct. Fix: `coreLevel = max(coreLevel, currentLevel)`.
One-line fix, verified on all test cases.

## Limitation: s > r + 1

Current cascade only handles s = r + 1 (each (r+1)-clique has r+1 r-cliques,
each differing by 1 vertex). For s > r + 1: s-cliques have C(s,r) r-cliques,
some differing by more than 1 vertex. Need generalized cascade.
