# V2 vs V3 Algorithm Comparison: Detailed Walkthrough

## Running Example

Graph with 4 vertices, 2 maximal cliques:
- M1 = {1, 2, 3, 4} (4-clique)
- M2 = {1, 2, 5, 6} (4-clique)

Parameters: r=2 (edges), s=3 (triangles).

### Step 1: Build Regions

| Vertex | prof(v) | Region (Class) |
|--------|---------|----------------|
| 1 | {M1, M2} | A |
| 2 | {M1, M2} | A |
| 3 | {M1} | B |
| 4 | {M1} | B |
| 5 | {M2} | C |
| 6 | {M2} | C |

Classes: A={1,2} (size 2), B={3,4} (size 2), C={5,6} (size 2).

### Step 2: Build r-tuples (r=2)

From M1 (classes A, B): tuples (A,A), (A,B), (B,B).
From M2 (classes A, C): tuples (A,A), (A,C), (C,C).

| r-tuple | mult | edges |
|---------|------|-------|
| τ1=(A,A) | C(2,2)=1 | {1,2} |
| τ2=(A,B) | 2×2=4 | {1,3},{1,4},{2,3},{2,4} |
| τ3=(A,C) | 2×2=4 | {1,5},{1,6},{2,5},{2,6} |
| τ4=(B,B) | C(2,2)=1 | {3,4} |
| τ5=(C,C) | C(2,2)=1 | {5,6} |

Total: 1+4+4+1+1 = 11 edges.

Verify: M1 has C(4,2)=6 edges, M2 has 6, shared: {1,2} → 6+6-1=11. ✓

---

## V2: s-tuple Peeling

### Step 3a: Build s-tuples (s=3)

From M1 (classes A, B):
- σ1=(A,A,B): mult = C(2,2)×2 = 2. Triangles: {1,2,3}, {1,2,4}.
- σ2=(A,B,B): mult = 2×C(2,2) = 2. Triangles: {1,3,4}, {2,3,4}.

From M2 (classes A, C):
- σ3=(A,A,C): mult = C(2,2)×2 = 2. Triangles: {1,2,5}, {1,2,6}.
- σ4=(A,C,C): mult = 2×C(2,2) = 2. Triangles: {1,5,6}, {2,5,6}.

4 s-tuples total, each with mult=2.

### Step 3b: Build bipartite incidence (r-tuple ↔ s-tuple)

For each s-tuple σ: find all r-sub-multisets → incident r-tuples.

**σ1=(A,A,B)**:
- 2-sub-multisets: (A,A), (A,B) → incident to τ1, τ2
- ext(σ1, τ1) = C(2-2, 1-0)×C(2-0, 1-0) = C(0,1)×C(2,1)... 

Actually, ext(σ, τ) = Π_c C(|c| - j_c^τ, σ_c - j_c^τ):
- ext(σ1, τ1): c=A: C(2-2, 2-2)=1, c=B: C(2-0, 1-0)=2 → ext=2.
- ext(σ1, τ2): c=A: C(2-1, 2-1)=1, c=B: C(2-1, 1-1)=1 → ext=1.

**σ2=(A,B,B)**:
- 2-sub-multisets: (A,B), (B,B) → incident to τ2, τ4
- ext(σ2, τ2): c=A: C(2-1, 1-1)=1, c=B: C(2-1, 2-1)=1 → ext=1.
- ext(σ2, τ4): c=A: C(2-0, 1-0)=2, c=B: C(2-2, 2-2)=1 → ext=2.

**σ3=(A,A,C)**:
- 2-sub-multisets: (A,A), (A,C) → incident to τ1, τ3
- ext(σ3, τ1): c=A: C(0,0)=1, c=C: C(2,1)=2 → ext=2.
- ext(σ3, τ3): c=A: C(1,1)=1, c=C: C(1,0)=1 → ext=1.

**σ4=(A,C,C)**:
- 2-sub-multisets: (A,C), (C,C) → incident to τ3, τ5
- ext(σ4, τ3): c=A: C(1,1)=1, c=C: C(1,0)=1 → ext=1.
- ext(σ4, τ5): c=A: C(2,1)=2, c=C: C(0,0)=1 → ext=2.

### Incidence table:

| | τ1(A,A) | τ2(A,B) | τ3(A,C) | τ4(B,B) | τ5(C,C) |
|---|---------|---------|---------|---------|---------|
| σ1(A,A,B) | ext=2 | ext=1 | | | |
| σ2(A,B,B) | | ext=1 | | ext=2 | |
| σ3(A,A,C) | ext=2 | | ext=1 | | |
| σ4(A,C,C) | | | ext=1 | | ext=2 |

### Step 4: Initial support

support(τ) = Σ_{σ incident} ext(σ, τ):
- τ1: ext(σ1)=2 + ext(σ3)=2 = **4**
- τ2: ext(σ1)=1 + ext(σ2)=1 = **2**
- τ3: ext(σ3)=1 + ext(σ4)=1 = **2**
- τ4: ext(σ2)=2 = **2**
- τ5: ext(σ4)=2 = **2**

Verify: edge {1,2} (τ1) is in triangles {1,2,3},{1,2,4},{1,2,5},{1,2,6} → support=4. ✓
Edge {1,3} (τ2) is in {1,2,3},{1,3,4} → support=2. ✓

### Step 5: Peeling

**Level 2**: τ2, τ3, τ4, τ5 all have support=2.

Peel τ2=(A,B), support=2. coreLevel=2.
- σ1 alive → mark dead.
  - τ1: support -= ext(σ1,τ1)=2 → 4-2=**2**.
- σ2 alive → mark dead.
  - τ4: support -= ext(σ2,τ4)=2 → 2-2=**0**.

Peel τ3=(A,C), support=2. coreLevel=2.
- σ3 alive → mark dead.
  - τ1: support -= ext(σ3,τ1)=2 → 2-2=**0**.
- σ4 alive → mark dead.
  - τ5: support -= ext(σ4,τ5)=2 → 2-2=**0**.

Peel τ4=(B,B), support=0. coreLevel=max(2,0)=2.
- σ2 already dead. No cascade.

Peel τ5=(C,C), support=0. coreLevel=2.

Peel τ1=(A,A), support=0. coreLevel=2.

**Result**: all tuples core=2. All 11 edges have core value 2.

### V2 Summary

- **Data**: 4 s-tuples, 8 incidence edges
- **Operations**: 4 s-tuple deaths, 4 ext subtractions
- **No new structures created during peeling**
- **Each s-tuple dies exactly once**

---

## V3: Analytical Split Peeling

### CPI Paths (from SDCT)

Assume the SDCT produces these paths:

| Path | holds | pivots | h | p |
|------|-------|--------|---|---|
| P1 | {1} | {2,3,4} | 1 | 3 |
| P2 | {1} | {2,5,6} | 1 | 3 |
| P3 (from v=2) | {2} | {3,4} | 1 | 2 |
| P4 (from v=2) | {2} | {5,6} | 1 | 2 |

(Simplified; actual SDCT structure depends on degeneracy ordering.)

Class distribution on each path:

| Path | holds | pivots | nh_A | np_A | nh_B | np_B | nh_C | np_C |
|------|-------|--------|------|------|------|------|------|------|
| P1 | {1} | {2,3,4} | 1 | 1 | 0 | 2 | 0 | 0 |
| P2 | {1} | {2,5,6} | 1 | 1 | 0 | 0 | 0 | 2 |
| P3 | {2} | {3,4} | 1 | 0 | 0 | 2 | 0 | 0 |
| P4 | {2} | {5,6} | 1 | 0 | 0 | 0 | 0 | 2 |

### Step 4: CPI Counting (AggrCount)

For each (tuple, path) pair, compute AggrCount / mult using the Vandermonde convolution.

**τ1=(A,A) on P1**: nh_A=1, np_A=1, j_A=2.
- g_A: b_A from max(0,2-1)=1 to min(2,1)=1. b_A=1.
  g_A[1] = C(1,1)×C(1,1) = 1. → g_A(x) = x.
- Classes not in τ: B with np=2.
  g_B(x) = C(2,0) + C(2,1)x + C(2,2)x² = 1 + 2x + x².
- f(x) = x × (1+2x+x²) = x + 2x² + x³.
- target = s-h = 3-1 = 2. f[2] = 2.
- AggrCount/mult = 2/1 = **2**.

**τ1=(A,A) on P2**: same structure → AggrCount/mult = **2**.

**τ1 on P3**: nh_A=1, np_A=0, j_A=2. Need 2 from A but only nh=1+np=0=1. **Infeasible → 0**.

**τ1 on P4**: same → **0**.

Total support(τ1) = 2+2+0+0 = **4**. ✓

**τ2=(A,B) on P1**: j_A=1, j_B=1.
- g_A: b_A from max(0,1-1)=0 to min(1,1)=1.
  g_A[0] = C(1,1)×C(1,0)=1, g_A[1] = C(1,0)×C(1,1)=1. → g_A = 1+x.
- g_B: nh_B=0, np_B=2. b_B from max(0,1)=1 to min(1,2)=1.
  g_B[1] = C(0,0)×C(2,1)=2... wait. 
  
  Using Vandermonde: g_c[t_c] = C(np_c, t_c) × C(nh_c + t_c, j_c).
  g_B[t]: t from max(0, 1-0)=1 to min(2, 2)=2? No: tMax = min(maxPiv, np) = min(2,2) = 2. tMin = max(minPiv=0, max(0, 1-0)=1) = 1.
  g_B[1] = C(2,1)×C(0+1,1) = 2×1 = 2.
  g_B[2] = C(2,2)×C(0+2,1) = 1×2 = 2.
  g_B = {0, 2, 2} → 2x + 2x².

- f = (1+x)(2x+2x²) = 2x + 2x² + 2x² + 2x³ = 2x + 4x² + 2x³.
- target = 2. f[2] = 4.
- AggrCount/mult = 4/4 = **1**.

**τ2=(A,B) on P3**: j_A=1, j_B=1. nh_A=1, np_A=0. nh_B=0, np_B=2.
- g_A: tMin=max(0,1-1)=0, tMax=min(0,0)=0. t_A=0.
  g_A[0] = C(0,0)×C(1+0,1) = 1×1 = 1. → g_A = {1}.
- g_B: same as P1 → {0, 2, 2}.
- f = {1} × {0,2,2} = {0,2,2}. target = s-h = 3-1 = 2. f[2] = 2.
- AggrCount/mult = 2/4 = **0.5**.

Hmm, 0.5 is not integer. Let me recheck.

Actually, target = s - h = 3 - 1 = 2 for P3 (h=1). f[2] = 2. AggrCount = 2. mult = 4. per-C_r' = 2/4 = 0.5.

But support should be integer! Let me verify manually.

Edges from τ2=(A,B): {1,3},{1,4},{2,3},{2,4}. Focus on P3 (holds={2}, pivots={3,4}).

Edge {2,3}: covered by P3 (both vertices in {2,3,4}). Triangles on P3 containing {2,3}: S must contain {2}(hold) and {3}(pivot of edge) plus one more pivot. Available: {4}. Triangle {2,3,4}. Count=1.

Edge {2,4}: similarly, triangle {2,3,4}. Count=1.

Edge {1,3}: vertex 1 is NOT on P3. Not covered. Count=0.

Edge {1,4}: vertex 1 not on P3. Count=0.

Total (C_r, S) pairs from τ2 on P3: 2. Per-C_r' = 2/4 = 0.5. ✓ (2 out of 4 edges contribute 1 each, other 2 contribute 0. Average = 0.5.)

So per-C_r' support is 0.5 from P3. This is correct — it's an average over all C_r' in the tuple.

Total support(τ2) = 1 (P1) + 0.5 (P3) + 0 (P2, no B) + 0 (P4, no B) = 1.5.

But earlier V2 gave support(τ2) = 2. MISMATCH!

This is because we're missing some paths. The actual SDCT would have more paths covering all edges. Let me add the missing ones.

Actually, the issue is my simplified SDCT doesn't cover all triangles. In reality, the SDCT from pivoter covers ALL cliques exactly once. Let me use a complete set:

For the graph M1={1,2,3,4}, M2={1,2,5,6}: the SDCT from vertex v=1 would produce paths covering cliques containing vertex 1. From v=2: paths covering cliques containing 2 but not 1. Etc.

For a correct example, let me use a simpler graph with just ONE maximal clique.

---

Actually, let me redo with a simpler example that avoids the SDCT complexity. Use ONE maximal clique M={1,2,3,4} with classes A={1,2}, B={3,4}.

SDCT from v=1: one path P1 with holds={1}, pivots={2,3,4}.
From v=2 (X={1}): BK finds sub-clique {2,3,4} ⊂ M. Path P2 with holds={2}, pivots={3,4}.
From v=3 (X={1,2}): {3,4}. Path P3 with holds={3}, pivots={4}. Too small for s=3 (h+p=2 < 3). Skipped.
From v=4 (X={1,2,3}): nothing. Skipped.

So: P1 (h=1,p=3) and P2 (h=1,p=2).

---

### Simplified Example (ONE clique, s=3)

M = {1,2,3,4}. Classes: A={1,2}, B={3,4}. r=2, s=3.

**CPI paths:**
- P1: holds={1}, pivots={2,3,4}. nh_A=1, np_A=1, nh_B=0, np_B=2.
- P2: holds={2}, pivots={3,4}. nh_A=1, np_A=0, nh_B=0, np_B=2.

**r-tuples:**
- τ1=(A,A): mult=1. Edge {1,2}.
- τ2=(A,B): mult=4. Edges {1,3},{1,4},{2,3},{2,4}.
- τ3=(B,B): mult=1. Edge {3,4}.

**s-tuples (V2):**
- σ1=(A,A,B): mult=2. Triangles {1,2,3},{1,2,4}.
- σ2=(A,B,B): mult=2. Triangles {1,3,4},{2,3,4}.

**Support via V2:**
- τ1: ext(σ1,τ1)=2 → support=2.
- τ2: ext(σ1,τ2)=1 + ext(σ2,τ2)=1 → support=2.
- τ3: ext(σ2,τ3)=2 → support=2.

All edges have support 2 (each in 2 triangles). All get core=2.

**Support via V3 (CPI counting):**

τ1=(A,A) on P1: g_A(x)=x (b_A=1 only). g_B(x)=1+2x+x² (free pivots).
f = x(1+2x+x²) = x+2x²+x³. target=2. f[2]=2. AggrCount/mult=2/1=2.

τ1 on P2: j_A=2, nh_A=1, np_A=0. Need 2 from A but only 1 available. **0**.

support(τ1) = 2+0 = **2**. ✓

τ2=(A,B) on P1: 
g_A: t_A ∈ [0,1]. g_A[0]=C(1,0)×C(1,1)=1, g_A[1]=C(1,1)×C(2,1)=2. → {1, 2}.
g_B: t_B ∈ [1,2]. g_B[1]=C(2,1)×C(1,1)=2, g_B[2]=C(2,2)×C(2,1)=2. → {0, 2, 2}.
f = {1,2} × {0,2,2} = {0, 2+0, 2+4, 0+4} = {0, 2, 6, 4}. target=2. f[2]=6.
AggrCount/mult = 6/4 = 1.5.

τ2 on P2:
g_A: nh_A=1, np_A=0. t_A ∈ [0,0]. g_A[0]=C(0,0)×C(1,1)=1. → {1}.
g_B: t_B ∈ [1,2]. g_B[1]=C(2,1)×C(1,1)=2, g_B[2]=C(2,2)×C(2,1)=2. → {0, 2, 2}.
f = {1} × {0,2,2} = {0,2,2}. target = 3-1 = 2. f[2]=2.
AggrCount/mult = 2/4 = 0.5.

support(τ2) = 1.5 + 0.5 = **2**. ✓

τ3=(B,B) on P1:
g_B: t_B ∈ [2,2]. g_B[2]=C(2,2)×C(2,2)=1. → {0,0,1}.
g_A: free. t_A ∈ [0,1]. g_A[0]=C(1,0)=1, g_A[1]=C(1,1)=1. → {1,1}.
f = {0,0,1} × {1,1} = {0,0,1,1}. target=2. f[2]=1.
AggrCount/mult = 1/1 = 1.

τ3 on P2:
g_B: t_B ∈ [2,2]. g_B[2]=C(2,2)×C(2,2)=1. → {0,0,1}.
g_A: nh_A=1, np_A=0. t_A=0. g_A[0]=C(0,0)=1. → {1}.
f = {0,0,1} × {1} = {0,0,1}. target=2. f[2]=1.
AggrCount/mult = 1/1 = 1.

support(τ3) = 1 + 1 = **2**. ✓

Great, all match V2.

---

### V2 Peeling (on simplified example)

Initial: all support=2. All peel at level 2.

Peel τ1=(A,A):
- σ1 alive → dead. τ2: support -= ext(σ1,τ2)=1 → 2-1=1.
- σ2: not incident to τ1. Skip.

Peel τ2=(A,B), support=1:
- σ1 already dead.
- σ2 alive → dead. τ3: support -= ext(σ2,τ3)=2 → 2-2=0.

Peel τ3=(B,B), support=0. coreLevel=max(2,0)=2.

Result: all core=2.

**V2 operations: 2 s-tuple deaths, 2 ext subtractions. No new structures.**

---

### V3 Analytical Split Peeling (on simplified example)

**Initial constrained paths:**
- CP1 = P1 with minPiv=0, maxPiv=original for all classes.
- CP2 = P2 with same.

tuple-to-cpath:
- τ1: {CP1} (not on CP2 because np_A=0 on P2, can't place 2 A-vertices)
- τ2: {CP1, CP2}
- τ3: {CP1, CP2}

**Peel τ1=(A,A) on CP1:**

τ's requirement: m_A = max(0, 2-1) = 1. m_B = 0.
Active classes: {A} (m_A=1 > minPiv_A=0). κ=1.

**Step 1: Compute old contributions for affected tuples (those sharing class A).**
- τ2 shares A: oldContrib(τ2, CP1) = 1.5
- τ3 does NOT share A: unaffected (skip aggrCount).

**Step 2: Split CP1 into 1 sub-path (κ=1).**

Part 1: maxPiv_A = m_A - 1 = 0 (cap A's pivots at 0).

Sub-CP1a: same as CP1 but maxPiv_A = 0.
- A can contribute 0 pivots (was 0-1, now 0-0).
- B unchanged (0-2).

**Step 3: Compute new contributions on Sub-CP1a for affected tuples.**

τ2 on Sub-CP1a: g_A: t_A ∈ [0, min(0, 1)] = [0,0]. 
g_A[0] = C(1,0)×C(1,1)=1. → {1}.
g_B: same as before → {0, 2, 2}.
f = {1}×{0,2,2} = {0,2,2}. target=2. f[2]=2.
newContrib(τ2) = 2/4 = 0.5.

**Step 4: Net delta for τ2.**
delta = oldContrib - newContrib = 1.5 - 0.5 = **1.0**.
dSup[τ2] -= 1.0 → 2 - 1 = **1**.

Verify: edge {1,2} (τ1) was in triangles {1,2,3},{1,2,4}. These die. 
τ2's edges: {1,3} was in {1,2,3} (dies) and {1,3,4} (alive). Net: -1 per C_r'. But only 2 of 4 edges lose 1, others lose 0. Average: -0.5? 

Hmm wait, let me recheck. Per-C_r' support decrease = dying triangles containing that specific C_r'.

{1,3}: in {1,2,3} (dies, contains τ1={1,2}) and {1,3,4} (alive). Decrease = 1.
{1,4}: in {1,2,4} (dies) and {1,3,4} (alive). Decrease = 1.
{2,3}: in {1,2,3} (dies) and {2,3,4} (alive). Decrease = 1.
{2,4}: in {1,2,4} (dies) and {2,3,4} (alive). Decrease = 1.

All 4 edges lose 1 triangle each. Per-C_r' decrease = 1. ✓ (delta = 1.0 for all C_r' in τ2.)

**Unaffected tuple τ3=(B,B):** not sharing A with τ1. oldContrib = newContrib → delta = 0. No update needed.

Verify: {3,4} is in {1,3,4} and {2,3,4}. Neither contains τ1={1,2}. No triangles die for τ3. ✓

**Sub-CP1a created.** tuple-to-cpath now:
- τ2: {Sub-CP1a, CP2}
- τ3: {Sub-CP1a, CP2}

(τ1 removed from CP1. CP1 replaced by Sub-CP1a.)

**Peel τ2=(A,B), support=1:**

On Sub-CP1a: m_A = max(0,1-1)=0. m_B = max(0,1-0)=1.
Active: {B} (m_B=1 > minPiv_B=0). κ=1.

Affected: τ3 shares B. Old: aggrCount(τ3, Sub-CP1a).

τ3 on Sub-CP1a: maxPiv_A=0.
g_B: t_B ∈ [2,2]. g_B[2] = C(2,2)×C(2,2) = 1. → {0,0,1}.
g_A: t_A ∈ [0,0]. g_A[0] = C(1,0) = 1 (no j_A for τ3). → {1}.
f = {0,0,1}×{1} = {0,0,1}. target=2. f[2]=1.
oldContrib(τ3) = 1/1 = 1.

Split: Part 1: maxPiv_B = 0.
Sub-CP1a-1: maxPiv_A=0, maxPiv_B=0.
g_B: t_B ∈ [0, min(0,2)] = [0,0]... but j_B=2 for τ3. Need t_B ≥ 2 but maxPiv_B=0. **Infeasible.**
newContrib(τ3) = 0.

delta(τ3) = 1 - 0 = 1. dSup[τ3] -= 1 → 2-1 = 1.

Hmm wait, but V2 gave support(τ3) = 0 after peeling τ1 and τ2. Let me continue.

On CP2: τ2's requirement: m_A = max(0,1-1)=0, m_B doesn't exist on CP2 (no B class).
Active classes on CP2: none (m_A=0, m_B=N/A). → **Path fully dead** for τ2.

Old contributions for τ3 on CP2:
τ3 on CP2: g_B: nh_B=0, np_B=2. t_B ∈ [2,2]. g_B[2]=1. g_A: np_A=0, j_A=0. g_A[0]=1. f={0,0,1}. target=2. f[2]=1. oldContrib=1.

newContrib(τ3, CP2) = 0 (path dead for τ2 → all s-cliques contain τ2's r-clique).

Wait, "path fully dead" means active=∅ → all s-cliques on CP2 contain τ2's r-clique? 
m_A=0 on CP2 means τ2 needs 0 pivots from A (j_A=1, nh_A=1 → already covered by hold). 
m_B=0 (B not on CP2). m_C: j_C=0 for τ2. So all m=0. 
Active = {c: m_c > minPiv_c} = {} (all m=0). Yes, fully dead.

delta(τ3, CP2) = 1 - 0 = 1. dSup[τ3] -= 1 → 1-1 = **0**.

Total: support(τ3) = 0. ✓ (matches V2).

**Peel τ3=(B,B), support=0.** coreLevel=max(2,0)=2. All core=2. ✓

### V3 operations summary:

| Step | Action | aggrCount calls | sub-paths created |
|------|--------|-----------------|-------------------|
| Peel τ1 on CP1 | old(τ2), new(τ2) on Sub-CP1a | 2 | 1 |
| Peel τ2 on Sub-CP1a | old(τ3), new(τ3) | 2 | 0 (infeasible) |
| Peel τ2 on CP2 | old(τ3) | 1 | 0 (dead path) |
| Peel τ3 | no paths left | 0 | 0 |
| **Total** | | **5** | **1** |

Compare V2: 2 s-tuple deaths, 2 subtractions = **4 operations**.

On this small example: V3 (5 ops) ≈ V2 (4 ops). Similar cost.

---

## Why V3 Explodes on Larger Graphs

On the small example: 3 tuples, 2 paths, 1 sub-path. Manageable.

On dblp-core30 (r=3, s=4):
- 1706 tuples, 1230 paths
- Each peeling creates ~2 sub-paths per affected path
- After all 1706 peelings: 28,000+ sub-paths
- Each tuple inherits parent's tuples → tuple-to-path grows from ~8 to ~300
- Total aggrCount calls: **1 million+**

The problem: **sub-paths inherit ALL parent tuples.** When CP splits into Sub-CP-a and Sub-CP-b:
- ALL 30 tuples on CP go to Sub-CP-a AND Sub-CP-b
- Each tuple now on 2 extra sub-paths
- Next peeling checks all sub-paths → 2x more work

This compounds: after k splits, a tuple can be on ~2^k sub-paths.

V2 avoids this because:
- s-tuples don't split — they just die (O(1) mark)
- Total deaths = T_s (fixed, never grows)
- No inheritance, no compounding

---

## The Fundamental Tradeoff

| | V2 | V3 |
|---|---|---|
| **Tracks** | What DIED (s-tuples) | What SURVIVED (sub-paths) |
| **Unit** | s-tuple (dies once) | Sub-path (splits into more) |
| **Growth** | Fixed at T_s | Grows to 28K+ |
| **Per unit cost** | O(C(s,r)) | O(tuples × convolution) |

**Tracking death = O(n). Tracking survival = O(2^n).**

This is why V2's s-tuple approach fundamentally scales better for peeling,
while V3's analytical split is theoretically novel but practically explosive.
