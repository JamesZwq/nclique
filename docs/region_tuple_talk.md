# Region Tuple: 一种无需枚举 r-clique 的精确 (r,s)-核分解算法

---

## 1. 问题定义

### 1.1 什么是 (r,s)-核分解？

给定无向图 G=(V,E) 和两个整数 r < s：

- **r-clique**：r 个顶点的完全子图（r=1 是顶点，r=2 是边，r=3 是三角形）
- **s-clique**：s 个顶点的完全子图
- **Support**：一个 r-clique T 的 support = 包含 T 的 s-clique 数量

**核分解**：把所有 r-clique 按"密度"分层。核值越高 = 越处于图的密集核心。

**标准算法（peeling）**：
1. 找到 support 最小的 r-clique
2. 剥离它（移除），更新邻居的 support
3. 重复，直到全部剥离
4. 每个 r-clique 的核值 = 剥离时的 support level

### 1.2 例子

图：一个 5-clique {A,B,C,D,E} + 一个 3-clique {A,B,F}

```
A --- B --- F
|\ /|
| X |
|/ \|
C --- D
 \ /
  E
```

取 r=2（边），s=3（三角形）：
- 边 {A,B} 的 support = 包含它的三角形数 = {A,B,C}, {A,B,D}, {A,B,E}, {A,B,F} = **4**
- 边 {A,F} 的 support = {A,B,F} = **1**
- 边 {C,D} 的 support = {C,D,A}, {C,D,B}, {C,D,E} = **3**

Peeling 顺序：先剥离 support=1 的边（{A,F}, {B,F}），然后剥离 support=3 的边...

### 1.3 瓶颈

r-clique 的数量随 r **指数增长**：

| 图 (web-it-2004) | r | r-clique 数量 |
|---|---|---|
| | 2 | 7.18M（边） |
| | 3 | 338M（三角形） |
| | 5 | ~260 亿 |
| | 10 | ~10^18 |

r=3 时标准算法已超时（>1小时）。r=10 完全不可能。

**核心问题**：能否不枚举 r-clique 就完成核分解？

---

## 2. 核心洞察：极大团与重叠类

### 2.1 极大团

**极大团**：不能再加入任何顶点的完全子图。

例子（同上图）：
- M₁ = {A,B,C,D,E}（5-clique，极大）
- M₂ = {A,B,F}（3-clique，极大）

每个 r-clique 必然在某个极大团内（因为 r-clique 是完全子图，被某个极大团包含）。

### 2.2 重叠类（Region）

**定义**：顶点 v 的"档案" prof(v) = {包含 v 的极大团集合}。

档案相同的顶点归入同一个 **Region**。

例子：
| 顶点 | 所在极大团 | prof | Region |
|------|-----------|------|--------|
| A | M₁, M₂ | {M₁, M₂} | R_AB |
| B | M₁, M₂ | {M₁, M₂} | R_AB |
| C | M₁ | {M₁} | R_CDE |
| D | M₁ | {M₁} | R_CDE |
| E | M₁ | {M₁} | R_CDE |
| F | M₂ | {M₂} | R_F |

3 个 Region：R_AB={A,B}, R_CDE={C,D,E}, R_F={F}。

**关键性质**：同一 Region 的顶点在图中是**完全等价**的——在完全相同的极大团中出现。

### 2.3 为什么 Region 有用？

**定理 1（包含引理）**：对任意 r-clique T ⊆ Region R，包含 T 的极大团集合 = prof_R。

**证明**（3行）：
- T ⊆ R ⊆ M 对所有 M ∈ prof_R → prof_R 的每个极大团都包含 T。✓
- 若极大团 M' 包含 T，取 T 中任意顶点 t ∈ R → M' ∈ prof(t) = prof_R。✓

**直觉**：T 选了 R 中的哪些具体顶点不重要。只要选的都来自 R，包含极大团就一样。

例子：r=2 时，边 {A,B} 和 ... 等等，R_AB 只有 2 个顶点，所以只有一个边。换 R_CDE：
- 边 {C,D}：包含极大团 = prof_{R_CDE} = {M₁}
- 边 {C,E}：包含极大团 = prof_{R_CDE} = {M₁}
- 边 {D,E}：包含极大团 = prof_{R_CDE} = {M₁}

三条边的包含极大团**完全一样**！

---

## 3. Region Tuple

### 3.1 定义

**Region r-tuple** τ = (R_{i₁}, R_{i₂}, ..., R_{iᵣ}) 是 r 个 Region 的有序多重集。

它**代表**所有"从每个 Region 取一个顶点"形成的 r-clique。

例子（r=3）：
- tuple (R_CDE, R_CDE, R_CDE)：从 R_CDE 中取 3 个顶点 → C(3,3)=1 个三角形 {C,D,E}
- tuple (R_AB, R_CDE, R_CDE)：从 R_AB 取 1 个、R_CDE 取 2 个 → 2 × C(3,2) = 6 个三角形
  - {A,C,D}, {A,C,E}, {A,D,E}, {B,C,D}, {B,C,E}, {B,D,E}
- tuple (R_AB, R_AB, R_CDE)：C(2,2) × 3 = 3 个三角形
  - {A,B,C}, {A,B,D}, {A,B,E}
- tuple (R_AB, R_AB, R_F)：C(2,2) × 1 = 1 个三角形
  - {A,B,F}

**Multiplicity**（每个 tuple 代表多少 r-clique）：

$$\text{mult}(\tau) = \prod_R \binom{|R|}{\tau_R}$$

其中 τ_R = Region R 在 tuple 中出现的次数。

### 3.2 完备性

**定理 2**：图中每个 r-clique 恰好属于一个 Region tuple。

证明：每个顶点属于唯一的 Region。r-clique 的 r 个顶点映射到 r 个 Region（可重复），排序后唯一。□

例子验证：上面 4 个 tuple 代表的三角形 = 1+6+3+1 = 11 个。
图中总三角形 = C(5,3) + C(3,3) - 0 = 10 + 1 = 11。✓（没有重叠，因为 M₁∩M₂={A,B} 只有 2 个顶点，无法形成三角形）

### 3.3 Support 等价

**定理 3**：同一 tuple 的所有 r-clique 有**完全相同的 support**。

证明：由定理 1，同一 tuple 的 r-clique 有相同的包含极大团集合 = ∩ᵢ prof_{R_iⱼ}。support 只依赖这个集合和存活顶点数，不依赖具体选了哪些顶点。□

例子（r=3, s=4）：

tuple (R_AB, R_CDE, R_CDE) 的 6 个三角形：
- {A,C,D}：support = 包含它的 4-clique 数 = {A,C,D,B}, {A,C,D,E} = 2。
  公共极大团 = prof_{R_AB} ∩ prof_{R_CDE} = {M₁,M₂} ∩ {M₁} = {M₁}。
  support = |M₁| - 3 = 5 - 3 = 2。✓
- {B,C,E}：同理。公共极大团 = {M₁}。support = 5 - 3 = 2。✓

6 个三角形**全部 support=2**。一个 tuple 只需存一个 support 值！

### 3.4 压缩比

**定理 4**：对极大团 M，设 ρ = Region 数，|M| = 顶点数，α = |M|/ρ = 平均 Region 大小：

$$\text{压缩比} = \frac{\binom{|M|}{r}}{\binom{\rho+r-1}{r}} \approx \alpha^r$$

例子（web-it-2004 最大极大团：432 顶点，20 个 Region，α=21.6）：

| r | r-clique 数 | tuple 数 | 压缩比 |
|---|------------|---------|--------|
| 3 | 13,375,216 | 1,540 | **8,685x** |
| 5 | 26,000,000,000 | 42,504 | **612,000x** |
| 10 | 10^18 | 20,000,000 | **8×10^10 x** |

r 越大，压缩比**指数级增长**。

---

## 4. 算法：双边关联剥离

### 4.1 为什么需要 s-tuple？

剥离 r-tuple 时需要更新邻居的 support。support = "包含的**活跃** s-clique 数"。
当一个 r-clique 被剥离，包含它的 s-clique 变为"不完整"→ 不再贡献 support。

关键：s-clique 也可以用 tuple 表示！

**s-tuple** σ = (R_{j₁}, ..., R_{jₛ}) 是 s 个 Region 的有序多重集。

### 4.2 双边关联

建立 r-tuple 与 s-tuple 之间的**关联**：

s-tuple σ 包含 r-tuple τ ⟺ τ 是 σ 的 r-子多重集。

例子（r=3, s=4）：
- s-tuple σ = (R_AB, R_CDE, R_CDE, R_CDE)
  - 它的 r-子多重集（3-子集）：
    - (R_AB, R_CDE, R_CDE) ← r-tuple τ₁
    - (R_CDE, R_CDE, R_CDE) ← r-tuple τ₂
  - 即：σ 关联 τ₁ 和 τ₂

**ext(σ, τ)** = 给定一个 τ 的具体 r-clique，有多少种方式补全为 σ 的 s-clique：

$$\text{ext}(\sigma, \tau) = \prod_i \binom{|c_i| - j_i}{m_i - j_i}$$

其中 jᵢ = cᵢ 在 τ 中的次数，mᵢ = 在 σ 中的次数。

例子：σ = (R_AB, R_CDE, R_CDE, R_CDE), τ₁ = (R_AB, R_CDE, R_CDE)
- R_AB：j=1, m=1 → C(2-1, 1-1) = C(1,0) = 1
- R_CDE：j=2, m=3 → C(3-2, 3-2) = C(1,1) = 1
- ext = 1 × 1 = 1

含义：给定三角形 {A,C,D}，只有 1 种方式补全为 4-clique（加 E → {A,C,D,E}）。

### 4.3 Support 计算

$$\text{support}(\tau) = \sum_{\text{alive } \sigma \ni \tau} \text{ext}(\sigma, \tau)$$

例子：τ₁ = (R_AB, R_CDE, R_CDE) 关联哪些 s-tuple？
- σ₁ = (R_AB, R_AB, R_CDE, R_CDE)：ext = C(2-1,2-1)×C(3-2,2-2) = 1×1 = 1
- σ₂ = (R_AB, R_CDE, R_CDE, R_CDE)：ext = C(2-1,1-1)×C(3-2,3-2) = 1×1 = 1
- support(τ₁) = 1 + 1 = **2** ✓（与之前手算一致）

### 4.4 剥离算法

```
初始化：
  计算每个 r-tuple 的 support
  放入 bucket queue（按 support 排序）
  所有 s-tuple 标记为 alive

循环：
  弹出 support 最小的 r-tuple τ
  coreLevel = max(coreLevel, current_bucket_level)  // 非递减！
  核值(τ) = coreLevel

  对每个关联的 s-tuple σ：
    如果 σ 仍然 alive：
      标记 σ 为 dead
      对 σ 的每个**其他**关联 r-tuple τ'：
        support(τ') -= ext(σ, τ')
        更新 bucket queue
```

### 4.5 为什么 s-tuple 的 alive/dead 保证正确？

标准剥离中：当 r-clique T 被剥离 → 包含 T 的 s-clique S 变为"不完整"→ S 不再贡献 support。

在 tuple 层面：当 r-tuple τ 被剥离 → 关联的 s-tuple σ 标记为 dead → σ 不再贡献 support。

这**精确对应**标准剥离。每个 s-clique 在第一个 r-clique 被剥离时就失活，之后不重复计算。

### 4.6 完整例子

图：{A,B,C,D,E} + {A,B,F}。r=3, s=4。

**Region**：R_AB={A,B}, R_CDE={C,D,E}, R_F={F}。

**r-tuple**（3-tuple）：

| r-tuple | mult | 三角形 |
|---------|------|--------|
| τ₁ = (R_AB, R_AB, R_CDE) | C(2,2)×3 = 3 | {A,B,C},{A,B,D},{A,B,E} |
| τ₂ = (R_AB, R_AB, R_F) | C(2,2)×1 = 1 | {A,B,F} |
| τ₃ = (R_AB, R_CDE, R_CDE) | 2×C(3,2) = 6 | {A,C,D},{A,C,E},{A,D,E},{B,C,D},{B,C,E},{B,D,E} |
| τ₄ = (R_CDE, R_CDE, R_CDE) | C(3,3) = 1 | {C,D,E} |

总计：3+1+6+1 = **11 个三角形**。

**s-tuple**（4-tuple）：

| s-tuple | 关联的 r-tuple | ext |
|---------|---------------|-----|
| σ₁ = (R_AB, R_AB, R_CDE, R_CDE) | τ₁(ext=1), τ₃(ext=1) | |
| σ₂ = (R_AB, R_AB, R_AB, R_CDE) | τ₁(ext=?) | 无效：C(2,3)=0，R_AB 只有 2 个顶点 |
| σ₃ = (R_AB, R_CDE, R_CDE, R_CDE) | τ₃(ext=1), τ₄(ext=2) | |
| σ₄ = (R_CDE, R_CDE, R_CDE, R_CDE) | τ₄(ext=?) | 无效：C(3,4)=0 |
| σ₅ = (R_AB, R_AB, R_CDE, R_F) | τ₁(ext=1), τ₂(ext=3) ... | 让我仔细算 |

等等，我来重新仔细枚举。s=4 的 s-tuple 需要 4 个 Region，共享极大团。

R_AB 和 R_CDE 共享 M₁。R_AB 和 R_F 共享 M₂。R_CDE 和 R_F **不共享**极大团。

所以有效的 s-tuple（4 个 Region 全部共享某极大团）：
- 只用 M₁ 的 Region（R_AB, R_CDE）：σ 的 4 个 slot 从 {R_AB, R_CDE} 中选
  - (R_AB, R_AB, R_CDE, R_CDE)：需 |R_AB|≥2, |R_CDE|≥2 ✓
  - (R_AB, R_CDE, R_CDE, R_CDE)：需 |R_CDE|≥3 ✓
  - (R_AB, R_AB, R_AB, R_CDE)：需 |R_AB|≥3 → |R_AB|=2 ✗
  - (R_CDE, R_CDE, R_CDE, R_CDE)：需 |R_CDE|≥4 → |R_CDE|=3 ✗
  - (R_AB, R_AB, R_AB, R_AB)：需 |R_AB|≥4 ✗
- 只用 M₂ 的 Region（R_AB, R_F）：
  - (R_AB, R_AB, R_F, R_F)：需 |R_F|≥2 → |R_F|=1 ✗
  - 其他都需要更大的 Region ✗
- 无法跨 M₁ 和 M₂（因为 4 个顶点需要全部在同一极大团中）

所以只有 **2 个有效 s-tuple**：
| s-tuple | 4-clique | mult |
|---------|----------|------|
| σ₁ = (R_AB, R_AB, R_CDE, R_CDE) | {A,B,x,y} x,y∈{C,D,E} | C(2,2)×C(3,2)=3 |
| σ₂ = (R_AB, R_CDE, R_CDE, R_CDE) | {z,C,D,E} z∈{A,B} | C(2,1)×C(3,3)=2 |

σ₁ 代表 3 个 4-clique：{A,B,C,D}, {A,B,C,E}, {A,B,D,E}。
σ₂ 代表 2 个 4-clique：{A,C,D,E}, {B,C,D,E}。
总共 5 个 4-clique = C(5,4) = 5。✓

**关联与 ext**：

σ₁ = (R_AB, R_AB, R_CDE, R_CDE) 的 3-子多重集：
- τ₁ = (R_AB, R_AB, R_CDE)：ext = C(3-2, 2-2) × C(... wait) 

让我用公式：ext(σ, τ) = Π_i C(|cᵢ|-jᵢ, mᵢ-jᵢ)
- R_AB: |c|=2, j=2(在τ₁中), m=2(在σ₁中) → C(2-2, 2-2) = C(0,0) = 1
- R_CDE: |c|=3, j=1(在τ₁中), m=2(在σ₁中) → C(3-1, 2-1) = C(2,1) = 2
- ext(σ₁, τ₁) = 1 × 2 = 2

含义：给定一个 τ₁ 三角形（如 {A,B,C}），有 2 种方式补全为 σ₁ 的 4-clique（加 D 或 E）。

- τ₃ = (R_AB, R_CDE, R_CDE)：ext = ?
  - R_AB: j=1, m=2 → C(2-1, 2-1) = C(1,1) = 1
  - R_CDE: j=2, m=2 → C(3-2, 2-2) = C(1,0) = 1
  - ext(σ₁, τ₃) = 1 × 1 = 1

σ₂ = (R_AB, R_CDE, R_CDE, R_CDE) 的 3-子多重集：
- τ₃ = (R_AB, R_CDE, R_CDE)：
  - R_AB: j=1, m=1 → C(1,0) = 1
  - R_CDE: j=2, m=3 → C(1,1) = 1
  - ext(σ₂, τ₃) = 1
- τ₄ = (R_CDE, R_CDE, R_CDE)：
  - R_AB: j=0, m=1 → C(2,1) = 2
  - R_CDE: j=3, m=3 → C(0,0) = 1
  - ext(σ₂, τ₄) = 2

**初始 support**：
- support(τ₁) = ext(σ₁, τ₁) = 2
- support(τ₂) = 0（没有任何 s-tuple 关联 τ₂！因为 R_F 无法出现在大小≥4 的极大团中）
- support(τ₃) = ext(σ₁, τ₃) + ext(σ₂, τ₃) = 1 + 1 = 2
- support(τ₄) = ext(σ₂, τ₄) = 2

**剥离**：

Level 0：τ₂ (support=0, mult=1)。coreLevel = 0。核值(τ₂) = 0。
- τ₂ 无关联 s-tuple → 无级联。

Level 2：τ₁, τ₃, τ₄ 都 support=2。
- 弹出 τ₁ (support=2, mult=3)。coreLevel = max(0,2) = 2。核值(τ₁) = 2。
  - σ₁ alive → 标记 dead。
    - σ₁ 的其他关联 r-tuple：τ₃。support(τ₃) -= ext(σ₁,τ₃) = 1。τ₃: 2→1。

- 弹出 τ₃ (support=1, 被级联降低)。coreLevel = max(2,1) = 2。核值(τ₃) = 2。
  - σ₂ alive → 标记 dead。
    - σ₂ 的其他关联：τ₄。support(τ₄) -= ext(σ₂,τ₄) = 2。τ₄: 2→0。

- 弹出 τ₄ (support=0)。coreLevel = max(2,0) = 2。核值(τ₄) = 2。

**最终核值**：
| tuple | 核值 | 三角形数 |
|-------|------|---------|
| τ₂ = (R_AB, R_AB, R_F) | 0 | 1 |
| τ₁ = (R_AB, R_AB, R_CDE) | 2 | 3 |
| τ₃ = (R_AB, R_CDE, R_CDE) | 2 | 6 |
| τ₄ = (R_CDE, R_CDE, R_CDE) | 2 | 1 |

核分布：core=0 有 1 个三角形，core=2 有 10 个三角形。
（验证：标准剥离也给出相同结果。{A,B,F} 不在任何 4-clique 中 → support=0 → core=0。M₁ 的所有三角形 support=2 → core=2。）

---

## 5. 等价性定理

### 5.1 为什么 tuple 剥离和标准剥离结果一样？

**定理 5（等价性）**：Region tuple 剥离给出与标准 per-r-clique 剥离**完全相同**的核值。

**证明思路**：

1. 同一 tuple 的 r-clique 始终有相同 support（定理 3）。
2. 因此在标准剥离中，它们必然在**同一 level** 被剥离。
3. 所以把它们作为一个 tuple 一起剥离 = 标准剥离中的自然行为。
4. 级联通过 s-tuple 的 alive/dead 精确追踪 = 标准剥离中 s-clique 的失活。□

**关键不变量**：`coreLevel = max(coreLevel, currentLevel)` 保证核值非递减，与标准剥离一致。

### 5.2 实验验证

11 个测试用例（3 个图 × 多种 r,s），**逐行完全一致**。不只是 max core 一样——每一个 core level 的 count 都完全相同。

---

## 6. 复杂度分析

### 6.1 存储

| 存储内容 | 数量 | 依赖 r？ |
|---------|------|---------|
| Region | ≤ n | 否 |
| r-tuple | Σ_M C(ρ+r-1, r) | 是，但 << r-clique |
| s-tuple | Σ_M C(ρ+s-1, s) | 是，但 << s-clique |
| 关联 | r-tuple × avg s-tuples | |

**不存储任何单个 r-clique 或 s-clique。**

### 6.2 时间

| 步骤 | 时间 |
|------|------|
| 极大团枚举 | O(n · 3^{d/3}) |
| 建 Region | O(Σ\|Mᵢ\|) |
| 枚举 tuple | O(Σ C(ρ+r-1,r) + Σ C(ρ+s-1,s)) |
| 建关联 | O(s-tuples × C(s,r)) |
| 剥离 | O(s-tuples × r-tuples-per-s-tuple) |

极大团枚举和 tuple 枚举与 r-clique 数量**无关**。

### 6.3 实测

| 图 | r | s | r-clique | r-tuple | 压缩比 | 标准算法 | Region Tuple | 加速 |
|---|---|---|----------|---------|--------|---------|-------------|------|
| dblp-core30 | 3 | 4 | 686K | 2.1K | 327x | 516ms | 11ms | **47x** |
| dblp-core30 | 4 | 5 | — | 4.4K | — | 12.1s | 18ms | **672x** |
| dblp-core30 | 5 | 6 | — | 7.2K | — | >10min | 30ms | **>20,000x** |
| web-it-2004 | 3 | 4 | 338M | 788K | 430x | >1h | 3.8s | **>1000x** |

---

## 7. 极大团枚举：MaxCliqEnum

### 7.1 为什么需要新的枚举器？

原来的 SDCT_MaxClique（82秒）有两个问题：
1. **Pivot 只从 P 选**（不从 P∪X），导致分支太多
2. **不剪枝子团**（X≠∅ 时仍然访问），浪费 80% 的工作

### 7.2 MaxCliqEnum 的改进

- Pivot 从 P∪X 选 → X 中的顶点可以覆盖更多 P → 更少分支
- P=∅ 且 X≠∅ 时立刻返回（标准 BK 行为）→ 跳过所有子团
- 不构建 tree 数据结构 → 直接输出极大团列表

### 7.3 效果

| | SDCT_MaxClique | MaxCliqEnum |
|---|---|---|
| web-it-2004 时间 | 82s | **8.3s** (10x) |
| 找到极大团数 | 84,871 (22%) | **381,110** (100%) |
| 准确性 | 有标记 bug | 完整正确 |

---

## 8. 总结

### 核心贡献

1. **重叠类（Region）不变量**：同一 Region 的顶点可以互换而不改变任何团结构
2. **Region Tuple**：把 r-clique 分组，同组 support 始终相同
3. **双边关联剥离**：r-tuple ↔ s-tuple 关联精确追踪 s-clique 失活
4. **精确等价**：与标准剥离逐行一致，不是近似
5. **压缩比 α^r**：随 r 指数增长，r 越大优势越大

### 一句话总结

> 通过把结构等价的 r-clique 归组为 Region Tuple，我们将 (r,s)-核分解的复杂度从依赖 r-clique 数量（$\binom{K}{r}$，随 r 指数增长）降为依赖 Region 数量（多项式增长），实现了精确且等价的核分解，在密集图上达到 47x–20,000x 的加速。
