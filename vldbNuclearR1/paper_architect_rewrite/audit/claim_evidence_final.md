# Final Claim-Evidence Ledger

| Claim | Evidence in rewrite | Evidence strength |
|---|---|---|
| \((1,s)\)-nucleus is the vertex-centric higher-order core target. | Definitions in Preliminaries; citations `nucleusSariyuce14`, `coreAlg`. | Strong. |
| CPI compactly encodes \(s\)-cliques through hold and pivot sets. | Definition `def:path`; support formula `lem:supp`; citation `NuclearCD`. | Strong. |
| General mutable CPI is necessary outside \(r=1\) because higher-\(r\) peeling can break path structure. | Preliminaries baseline; StaticCPI boundary paragraph; citation `NuclearCD`. | Strong as scoped prior-work contrast. |
| For \(r=1\), CPI path structure remains valid throughout vertex peeling. | Theorem `thm:vertex-removal`. | Strong. |
| Hold removal kills a path contribution and pivot removal decrements an effective pivot counter. | Theorem `thm:vertex-removal`; Example `ex:peel`; Figures `fig:overview`, `fig:peel`. | Strong. |
| Support losses are closed-form binomial deltas. | Lemma `lem:support-delta`. | Strong. |
| LocalH computes the \((1,s)\)-core fixed point. | Algorithm `alg:localH`; Lemma `lem:hindex-fixed-point`; Theorem `thm:localh-complexity`. | Strong at paper-proof level. |
| PeelH computes the same fixed point as LocalH. | Algorithm `alg:peelH`; Lemma `lem:equivalence`. | Strong, but proof should be checked by theory coauthor. |
| PeelH runs in \(O(\Sig+U+\kmax+n)\). | Theorem `thm:peelh-complexity`; \(U\) defined in HIndex. | Strong if \(U\) accounting matches implementation. |
| ParaBuild produces the same CPI semantics with thread-local path emission and deterministic merge. | Algorithm `alg:parabuild`; correctness paragraph; prior BK/CPI citations. | Moderate to strong; full BkEmit inherited from prior work. |
| BuildHier extracts hierarchy in \(O((n+m)\alpha(n))\). | Algorithm `alg:buildhier`; correctness and complexity paragraphs. | Strong. |
| Matched configurations agree with REF. | Experimental RQ1. | Moderate; TODO remains for full 762 verification. |
| CSV-backed subset speedup is \(6.16\times\) geomean and max \(39.41\times\). | Experimental RQ2; `figures/benchmark_all_results.csv`; Figure `fig:exp-endtoend`. | Strong for three-graph 505-cell subset. |
| Ten-graph headline speedup is \(7.42\times\). | Table `tab:per-graph`; existing figures. | Pending verification; TODO retained. |
| Memory reduction follows from flat incidence layout. | Table `tab:mem`; Figures `fig:exp-endtoend`, `fig:exp-mem`; Table `tab:bd-mem`. | Strong qualitatively; quantitative maxima need source check. |
| Construction becomes important after cheap peeling. | Phase table `tab:bd-time`; ParaBuild motivation. | Strong for reported rows. |
| ParaBuild scales on construction workloads. | Table `tab:par`. | Moderate; TODO retained for raw source verification. |
| LocalH agrees with PeelH but is slower on current single-machine implementation. | Experimental RQ5. | Moderate; TODO retained for raw logs. |
| Dense stress test favors V3 where both complete. | Figures `fig:stress-time`, `fig:stress-mem`; RQ6. | Moderate; TODO retained for exact raw values. |
| \((1,s)\) matches higher-\(r\) quality on reported community tasks at lower cost. | CaseStudies Tables `tab:cs4`, `tab:cs7`; Figures `fig:cs6`, `fig:cs7y`. | Strong for scoped tasks. |
| \(s\) controls granularity. | CaseStudies Table `tab:cs1`, Figure `fig:cs1`, ego-network Table `tab:cs9`. | Strong for reported data. |
