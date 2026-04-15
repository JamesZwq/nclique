Case B / pathSplit 的精确剪枝分析

1. 问题定义

设当前一个待处理 leaf 的 pivot 集合为 P。
由已删除的 r-clique 诱导出一族冲突集：

\mathcal{F} = \{F_1, F_2, \dots, F_m\}, \qquad F_i \subseteq P, \; F_i \neq \emptyset.

语义为：若一个冲突集 F_i 没有被“击中”，则当前 leaf 不合法。
“击中”指删除至少一个属于 F_i 的 pivot。

因此，Case B / pathSplit 的核心子问题可形式化为：

在保留至少 needPivot 个 pivot 的前提下，删除尽量少的 pivots，使得每个冲突集都被击中。

记当前活跃 pivot 数为 a = |P|，则允许删除的最大 pivot 数为

B = a - needPivot.

若最小 hitting set 大小大于 B，则该 leaf 必死，可直接剪枝。

⸻

2. 基线算法

定义最优值：

\tau(\mathcal{F}) = \min \{ |X| : X \subseteq P,\; X \cap F_i \neq \emptyset,\; \forall i \}.

基线递归为：
	1.	若 \mathcal{F} = \emptyset，返回 0。
	2.	选一个未满足的冲突集 F\in\mathcal{F}。
	3.	对每个 p\in F，递归求解“删除 p”后的子问题。
	4.	取最小值。

该方法 exact，但分支树很大，遇到高重叠或大实例时容易超时。

⸻

3. 可证明正确的剪枝规则

⸻

定理 1. 重复冲突集删除

若存在 F_i = F_j，删除其中任意一个，不改变最优值 \tau(\mathcal{F})。

证明
任何 hitting set 对两个相同集合的命中条件完全一致。删除重复约束不改变可行解集合，因此不改变最优值。∎

⸻

定理 2. 包含消解

若存在 F_i \subseteq F_j，则删除 F_j 后最优值不变。

证明
任取 hitting set X。若 X 击中 F_i，则必有 X\cap F_i\neq\emptyset。由于 F_i\subseteq F_j，有 X\cap F_j\neq\emptyset。
因此 F_j 的约束被 F_i 蕴含，删除 F_j 不改变可行解集合，也不改变最优值。∎

⸻

定理 3. 单元传播

若存在 singleton 冲突集 F_i=\{p\}，则任意可行解都必须删除 p。

证明
若某可行解 X 不删除 p，则 X\cap F_i=\emptyset，与可行性矛盾。故 p 必须被删除。∎

⸻

定理 4. 不相交打包下界

若存在一组两两不相交的冲突集

F_{i_1}, F_{i_2}, \dots, F_{i_t},
\qquad F_{i_u}\cap F_{i_v}=\emptyset \; (u\neq v),

则

\tau(\mathcal{F}) \ge t.

证明
每个冲突集至少需要删除一个属于它的 pivot 才能被击中。由于这些集合两两不相交，击中它们至少需要删除 t 个不同 pivot。故 \tau(\mathcal{F}) \ge t。∎

⸻

推论 4.1. Leaf-death 剪枝

设允许删除预算为

B = a - needPivot.

若通过不相交打包得到下界 LB，且

LB > B,

则该 leaf 必死，可直接剪枝。

证明
由定理 4，最优删除数至少为 LB。若 LB>B，则任何可行解都需要删除超过预算的 pivot 数，因此 leaf 不可能存活。∎

⸻

定理 5. 连通分量分解

构造 pivot 冲突图 G_c：
	•	顶点集为当前 pivots；
	•	若两个 pivots 同时出现在某个冲突集中，则在二者之间连边。

若 G_c 分解成连通分量

C_1, C_2, \dots, C_k,

则

\tau(\mathcal{F}) = \sum_{j=1}^k \tau(\mathcal{F}[C_j]),

其中 \mathcal{F}[C_j] 为完全包含于分量 C_j 的冲突集子族。

证明
不同连通分量之间没有共享 pivot，因此不存在跨分量的冲突集。每个分量的命中决策互相独立。
任意全局解可分解为各分量子解之并，总代价为分量代价之和；反之，各分量最优解之并构成全局最优解。故结论成立。∎

⸻

推论 5.1. 分量下界加和

若对每个分量 C_j 都能计算安全下界 LB_j，则

\tau(\mathcal{F}) \ge \sum_{j=1}^k LB_j.

若

\sum_{j=1}^k LB_j > B,

则可直接判定 leaf 必死。

∎

⸻

4. 建议的 exact solver 结构

4.1 预处理层

对当前冲突族 \mathcal{F} 先执行：
	1.	重复冲突集删除
	2.	包含消解
	3.	singleton 传播直到稳定

该层不改变正确性。

⸻

4.2 下界层

对化简后的冲突族计算：
	1.	greedy disjoint packing lower bound
	2.	若检测到多连通分量，则对每个分量分别计算下界并求和

若下界超过预算 B，直接判死。

⸻

4.3 分支层

若未被剪掉：
	1.	选一个冲突集 F
	2.	对每个 p\in F 分支：删除 p
	3.	递归求解

推荐分支变量选择策略：
	•	优先选择最小冲突集
	•	若并列，优先选择参与冲突集数最多的 pivot

⸻

5. Python 实验设计

我在本地 Python 环境中实现并测试了以下 exact solver：
	•	baseline：裸递归
	•	simplify：基线 + 去重/包含消解/singleton
	•	pack_lb：simplify + 不相交打包下界
	•	components：pack_lb + 连通分量分解
	•	memo：components + 状态记忆化

输出文件：
	•	benchmark script￼
	•	regular results csv￼
	•	regular summary￼
	•	hard-instance results csv￼
	•	hard-instance summary￼

⸻

6. 普通难度实验结果

测试 family：
	•	uniform_sparse
	•	uniform_dense
	•	subsumption
	•	components
	•	clustered

每类 30 个实例，总计 150 个实例。

6.1 几何平均 speedup（相对 baseline）

Solver	Geomean speedup
simplify	< 1.0x
pack_lb	1.69x
components	条件型强
memo	无稳定额外收益

6.2 最优方案分布

Family	最优方案	结论
uniform_sparse	pack_lb	通用小幅稳定提升
uniform_dense	pack_lb	明显优于 baseline
subsumption	pack_lb	包含消解收益明显
components	components / memo	分量分解极强
clustered	pack_lb	通用最佳

6.3 节点数削减

pack_lb 对递归节点数的削减：
	•	uniform_dense：约 7.29x fewer nodes
	•	subsumption：约 9.12x fewer nodes
	•	clustered：约 9.06x fewer nodes

结论：

pack_lb 是最稳的通用赢家。
components 仅在冲突图真的能裂成多个分量时，收益会极大。

⸻

7. Hard-instance 实验结果

我补做了“baseline 可能原本跑不出来”的 harder instances，并统一使用 3 秒超时。

测试 family：
	•	uniform_hard
	•	subsumption_hard
	•	components_hard
	•	clustered_hard

7.1 结果总表

Family	baseline	pack_lb	components	结论
uniform_hard	timeout	timeout	timeout	无明显结构时剪枝也难救
subsumption_hard	0.021s	0.0041s	0.0071s	pack_lb 最优
components_hard	timeout	0.561s	0.00205s	分量分解压倒性优势
clustered_hard	timeout	1.117s	1.764s	pack_lb 可将 timeout 救活

7.2 关键结论
	1.	确实存在 baseline 原本超时，而剪枝后能解出的情况。
	2.	最典型的两类：
	•	components_hard：baseline 超时，components 0.002 秒解出
	•	clustered_hard：baseline 超时，pack_lb 1.117 秒解出
	3.	uniform_hard 表明：若 conflict hypergraph 极其均匀、缺乏结构，剪枝也未必足够。

⸻

8. 实验结论

8.1 最值得先做的剪枝

\boxed{\text{subsumption + singleton + packing lower bound}}

理由：
	•	完全 exact
	•	通用
	•	工程复杂度低
	•	在普通难度和 hard instances 上都稳定有效

⸻

8.2 最值得作为条件型 fast path 的剪枝

\boxed{\text{component decomposition}}

理由：
	•	当冲突图能明显分裂时，收益极大
	•	在 components_hard 上从 baseline timeout 直接降到 0.002 秒
	•	不适合无脑全开，但适合条件触发

⸻

8.3 暂不建议优先做的

\boxed{\text{memoization}}

原因：
	•	当前实验里没有稳定增益
	•	实现更复杂
	•	只有在真实 Case B 上确认有大量重复状态时才值得

⸻

9. 对当前 C++ 主实现的建议

9.1 第一优先级

在 Case B / pathSplit 中实现以下 exact pruning：
	1.	conflict-set deduplication
	2.	subset subsumption elimination
	3.	singleton / unit propagation
	4.	greedy disjoint-packing lower bound
	5.	若 LB > a - needPivot，直接 leaf-death

⸻

9.2 第二优先级

构造 pivot 冲突图，若发现可明显分解，则启用：
	6.	connected-component decomposition
	7.	分量下界求和
	8.	对小分量单独精确求解

⸻

9.3 暂缓

不要优先继续做：
	•	更多 bitset / hash containment 常数优化
	•	memoization
	•	不带强下界的普通 BK 微调

⸻

10. 最终结论

定理化结论

当前 Case B / pathSplit 可被严格表述为 hitting-set 子问题。
在该表述下：
	•	重复删除
	•	包含消解
	•	singleton 传播
	•	packing 下界
	•	component decomposition

都属于 exact pruning，不会改变最终 core 结果。

实验化结论

在本地 Python 测试中：
	•	pack_lb 是最稳的通用方案
	•	components 在可分解实例上极强
	•	存在 baseline 原本超时、而加剪枝后可解出的 harder cases

因此，当前最值得落地到 C++ 主实现的方案是：

\boxed{
\text{Case B} \leftarrow
\text{subsumption} +
\text{singleton propagation} +
\text{packing lower bound} +
\text{conditional component decomposition}
}

这是一条有理论保证、且经实验验证有效的剪枝主线。
