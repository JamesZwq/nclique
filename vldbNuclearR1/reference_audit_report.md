# `references.bib` reference audit report

**核查日期**：2026-05-20（Australia/Sydney）  
**核查范围**：`references.bib` 中 51 个非 `@String` reference 条目。  
**核查口径**：逐条按 title / authors / year / venue / volume / issue / pages or article number / DOI 进行线上核对；优先采用出版社页面、ACM/IEEE/Springer/ScienceDirect/PNAS/PLOS/Nature/APS 页面、DBLP、作者/机构页面。  

## 结论

原文件中的 reference 基本都能对应到真实论文或正式出版记录；没有发现明显伪造的文献条目。 但原文件不能视为“每个字段完全匹配”：存在若干字段错误、格式错误、缺失 DOI、以及一个不属于任何 BibTeX entry 的杂散文本行。

我另存了一个带高置信度修正的文件：`references_checked_corrected.bib`。该修正版只改动了我能从线上元数据直接确认的内容；对 `sotaPrveHierarchy` 的年份差异和 `NuclearCD` 的 2026 未来/预发表元数据没有强行改写，而是在文件中保留注释，并在本报告中标出。

## 必须修正的问题

| BibTeX key / 位置 | 原问题 | 核查结果 | 建议修正 |
|---|---|---|---|
| stray line after `nuclearTWEB` | `Using Large Cliques for Hierarchical Dense Subgraph Discovery.` 出现在 entry 外部 | 这不是一个 BibTeX 条目，会破坏/污染 `.bib` 文件 | 删除或注释该行。修正版已注释 |
| `sotaHierarchy` | `month = mar` | ACM reference format / PDF 元数据为 February 2024；entry 自身 `issue_date` 也是 February 2024 | 改为 `month = feb` |
| `hindex` | `abstract` 被截断，内容错误；缺少 DOI | 官方元数据为 PNAS 102(46):16569-16572，DOI `10.1073/pnas.0507655102` | 修正 abstract；添加 DOI |
| `nuclearTWEB` | `issue_date = {August 2017}` | ACM/TWEB 元数据为 July 2017；DOI `10.1145/3057742` | 改为 `issue_date = {July 2017}` |
| `akoglu2015graph` | 标题写作 `Graph-based...`；缺少 DOI | 出版元数据标题为 `Graph based anomaly detection and description: A survey`，DOI `10.1007/s10618-014-0365-y` | 标题按官方元数据修正；添加 DOI |
| `communities_dec3` | DOI 字段写成 URL；`number={3}` | ScienceDirect/Physics Reports 元数据为 vol. 486, issues 3-5, pages 75-174, DOI `10.1016/j.physrep.2009.11.002` | DOI 改为裸 DOI；`number={3--5}`；pages 用 `75--174` |
| `Identifying_influential` | DOI 字段写成 URL；pages 用单 hyphen | ScienceDirect 元数据 DOI 为 `10.1016/j.physa.2011.09.017`，pages 1777-1787 | DOI 改为裸 DOI；pages 用 `1777--1787` |
| `apply_biological_networks_1` | `pages = {1-6}`；`abstract = {null}` | PLOS Computational Biology 使用 article ID `e1000385`；`null` 不是有效 abstract | 改为 `pages = {e1000385}`；删除 `abstract = {null}` |
| `kq-core` | DOI 字段写成 URL | DOI 应为 `10.1016/j.chaos.2023.113645` | DOI 改为裸 DOI |
| `e-shipping_recommendation` | DOI 字段写成 URL；pages 用单 hyphen | DOI 应为 `10.1016/j.ins.2018.07.029`；pages 269-287 | DOI 改为裸 DOI；pages 用 `269--287` |

## 需要决策/不能完全确认的问题

| BibTeX key | 状态 | 说明 | 建议 |
|---|---|---|---|
| `sotaPrveHierarchy` | 年份存在来源差异 | DBLP/部分索引按 2021 记录；PVLDB PDF 的 reference format 写作 PVLDB 15(3):583-596, 2022。原文件 `year={2021}` 不是“凭空错误”，但若严格按 PVLDB 文章 PDF 的 reference format，year 应为 2022 | 按目标会议/期刊格式决定：若使用 PVLDB reference format，改为 2022；若沿用 DBLP 元数据，保留 2021 |
| `NuclearCD` | 只能部分确认 | SIGMOD 2026 accepted-papers/作者页面可确认 title 和 authors；但我没有在可检索的 ACM/DOI 页面独立确认 DOI `10.1145/3802092`、PACMMOD vol. 4 no. 3、article 215、26 pages 等完整出版元数据。考虑到 issue date 是 June 2026，而当前核查日期为 2026-05-20，这可能是尚未公开落地的预发布元数据 | 不建议在最终投稿版中声称该 entry 已完全核实；等 ACM DOI 页面上线后再核对。修正版保留原 entry 并加注释 |

## 已补充 DOI 或出版细节的条目

这些不是“论文不存在”的问题，而是原条目不完整。修正版已添加可确认 DOI/细节。

| BibTeX key | 补充内容 |
|---|---|
| `zhang2023efficient` | DOI `10.1109/ICDMW60847.2023.00135` |
| `io-eff` | DOI `10.1109/ICDE.2016.7498235` |
| `kcoreSeidman` | DOI `10.1016/0378-8733(83)90028-X` |
| `nucleusSariyuce14` | DOI `10.1145/2736277.2741640` |
| `tomita` | DOI `10.1016/j.tcs.2006.06.015` |
| `eppstein` | DOI `10.1007/978-3-642-17517-6_36`；补充 LNCS series/volume |
| `lu2016hindex` | DOI `10.1038/ncomms10168` |
| `communityFinder` | DOI `10.1007/s10115-013-0693-z` |
| `edelsbrunner2010topology` | DOI `10.1090/mbk/069` |
| `cohesive_subgraph_protein_networks` | DOI `10.1126/science.1065103` |
| `cohesive_subgraph_complex_networks` | DOI `10.1038/nphys1746` |
| `DBLP:journals/tkde/WenQZLY19` | DOI `10.1109/TKDE.2018.2833070` |
| `hyper_kcore` | DOI `10.1103/PhysRevE.109.014307` |

## 逐条核查结果

| # | BibTeX key | 核查状态 | 主要核查结论 / 修正说明 |
|---:|---|---|---|
| 1 | `LocalNuclear` | 匹配 | title/authors/PVLDB 12(1):43-56/2018/DOI `10.14778/3275536.3275540` 可核实。 |
| 2 | `fastHierarchy` | 匹配 | title/authors/PVLDB 10(3):97-108/2016/DOI `10.14778/3021924.3021927` 可核实。 |
| 3 | `sotaPrveHierarchy` | 需决策 | title/authors/pages/DOI 匹配；年份存在 PVLDB reference format 与 DBLP/索引差异，见上表。 |
| 4 | `sotaHierarchy` | 已修正 | DOI/title/authors/venue/article no. 匹配；`month` 应为 February，不是 March。 |
| 5 | `hindex` | 已修正 | 论文真实；PNAS 102(46):16569-16572，DOI 已添加；abstract 原文被截断。 |
| 6 | `Pivoter` | 匹配 | WSDM 2020, pages 268-276, DOI `10.1145/3336191.3371839` 可核实。 |
| 7 | `BronKerbosch` | 匹配 | CACM 16(9):575-577, DOI `10.1145/362342.362367` 可核实。 |
| 8 | `girvan2002community` | 匹配 | PNAS 99(12):7821-7826, DOI `10.1073/pnas.122653799` 可核实。 |
| 9 | `kempe2003maximizing` | 匹配 | KDD 2003, pages 137-146, DOI `10.1145/956750.956769` 可核实。 |
| 10 | `akoglu2015graph` | 已修正 | 论文真实；官方标题无 hyphen，且首字母为 `A survey`；DOI 已添加。 |
| 11 | `lu2016vital` | 匹配 | Physics Reports 650:1-63, DOI `10.1016/j.physrep.2016.06.007` 可核实。 |
| 12 | `batagelj2011fast` | 匹配 | Advances in Data Analysis and Classification 5(2):129-145, DOI `10.1007/s11634-010-0079-y` 可核实。 |
| 13 | `nuclearTWEB` | 已修正 | TWEB 11(3), Article 16, DOI `10.1145/3057742` 可核实；issue date 修正为 July 2017。 |
| 14 | `zhang2023efficient` | 已补充 | ICDMW 2023 pages 1023-1031 可核实；DOI 已补充。 |
| 15 | `io-eff` | 已补充 | ICDE 2016 pages 133-144 可核实；DOI 已补充。 |
| 16 | `peel_cd` | 匹配 | PVLDB 9(1):13-23, DOI `10.14778/2850469.2850471` 可核实。 |
| 17 | `kcore_4` | 匹配 | PVLDB 6(6):433-444, DOI `10.14778/2536336.2536344` 可核实。 |
| 18 | `huang2014truss` | 匹配 | SIGMOD 2014 pages 1311-1322, DOI `10.1145/2588555.2610495` 可核实。 |
| 19 | `akbas2017truss` | 匹配 | PVLDB 10(11):1298-1309, DOI `10.14778/3137628.3137640` 可核实。 |
| 20 | `chen2020truss` | 匹配 | PVLDB 13(10):1751-1764, DOI `10.14778/3401960.3401971` 可核实。 |
| 21 | `NuclearCD` | 部分确认 | title/authors/SIGMOD 2026 acceptance 可确认；DOI/volume/article/page metadata 暂未能独立确认。 |
| 22 | `coreAlg` | 匹配 | arXiv/CoRR `cs.DS/0310049`, 2003 可核实；没有 DOI。 |
| 23 | `kcoreSeidman` | 已补充 | Social Networks 5(3):269-287, DOI 已补充。 |
| 24 | `ktrussCohen` | 匹配 | Tech report title/author/institution/year 可核实；未发现 DOI。 |
| 25 | `nucleusSariyuce14` | 已补充 | WWW 2015 pages 927-937 可核实；DOI 已补充。 |
| 26 | `tomita` | 已补充 | Theoretical Computer Science 363(1):28-42 可核实；DOI 已补充。 |
| 27 | `eppstein` | 已补充 | ISAAC 2010 LNCS 6506:403-414 可核实；DOI/series/volume 已补充。 |
| 28 | `lu2016hindex` | 已补充 | Nature Communications 7, article 10168 可核实；DOI 已补充。 |
| 29 | `communityFinder` | 已补充 | Knowledge and Information Systems 42(1):181-213 可核实；DOI 已补充。 |
| 30 | `edelsbrunner2010topology` | 已补充 | AMS book title/authors/year 可核实；DOI 已补充。 |
| 31 | `communities_dec1` | 匹配 | SIGMOD 2014 pages 991-1002, DOI `10.1145/2588555.2612179` 可核实。 |
| 32 | `communities_dec2` | 匹配 | ASONAM 2011 pages 87-93, DOI `10.1109/ASONAM.2011.65` 可核实。 |
| 33 | `communities_dec3` | 已修正 | DOI/issue/pages 规范化；issues 应为 3-5。 |
| 34 | `dense_subgraph` | 匹配 | KDD 2015 tutorial pages 2313-2314, DOI `10.1145/2783258.2789987` 可核实。 |
| 35 | `Identifying_influential` | 已修正 | DOI 规范化；pages 改为 BibTeX `--`。 |
| 36 | `network_visualization` | 匹配 | PVLDB/SIGMOD-style record DOI `10.14778/2535568.2448942` 可核实。 |
| 37 | `cohesive_subgraph_Social_Networks` | 匹配 | NIPS 2005/NeurIPS proceedings pages 41-50 可核实；未发现 DOI。 |
| 38 | `cohesive_subgraph_protein_networks` | 已补充 | Science 296(5569):910-913 可核实；DOI 已补充。 |
| 39 | `cohesive_subgraph_complex_networks` | 已补充 | Nature Physics 6(11):888-893 可核实；DOI 已补充。 |
| 40 | `k-core1` | 匹配 | IEEE TKDE 34(7):3126-3138, DOI `10.1109/TKDE.2020.3023925` 可核实。 |
| 41 | `DBLP:journals/tkde/WenQZLY19` | 已补充 | IEEE TKDE 31(1):75-90 可核实；DOI 已补充。 |
| 42 | `kcore_3` | 匹配 | KDD 2014 pages 1316-1325, DOI `10.1145/2623330.2623655` 可核实。 |
| 43 | `distributed_k-core` | 匹配 | IEEE TPDS 24(2):288-300, DOI `10.1109/TPDS.2012.124` 可核实。 |
| 44 | `nbr_cd` | 匹配 | PVLDB 16(9):2061-2074, DOI `10.14778/3598581.3598582` 可核实。 |
| 45 | `hyper_kcore` | 已补充 | Physical Review E 109:014307 可核实；DOI 已补充。 |
| 46 | `kg-core` | 匹配 | CIKM 2023 pages 4013-4017, DOI `10.1145/3583780.3615275` 可核实。 |
| 47 | `kq-core` | 已修正 | DOI 规范化为裸 DOI；volume/article ID 匹配。 |
| 48 | `apply_social_networks` | 匹配 | arXiv preprint physics/0505137, 2005 可核实；未发现 DOI。 |
| 49 | `apply_biological_networks_1` | 已修正 | PLOS article ID 应为 `e1000385`；`abstract={null}` 删除。 |
| 50 | `apply_cyber` | 匹配 | EICC 2024 pages 161-162, DOI `10.1145/3655693.3661300` 可核实。 |
| 51 | `e-shipping_recommendation` | 已修正 | DOI 规范化；pages 改为 BibTeX `--`。 |

## 建议

1. 若用于投稿，请优先使用 `references_checked_corrected.bib`，但在提交前对 `NuclearCD` 再做一次 ACM DOI 页面确认。
2. 若目标 venue 要求严格使用 publisher metadata，建议把 `sotaPrveHierarchy` 的年份按 PVLDB PDF reference format 改为 2022；若沿用 DBLP/索引元数据，保留 2021 也有来源依据。
3. 原文件中有很多 ACM sample `@String` 定义和注释，并不影响 reference 条目真实性，但如果这是论文最终 `.bib`，建议清理掉无关 sample 内容，降低编译/审稿风险。
