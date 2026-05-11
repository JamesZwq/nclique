# Claim-Evidence Ledger

| Claim | Evidence in current repo | Status |
|---|---|---|
| \((1,s)\)-nucleus generalizes \(k\)-core and triangle-support ranking. | Definitions in `sections/Preliminaries.tex`; citations `nucleusSariyuce14`, `coreAlg`. | Supported. |
| CPI encodes \(s\)-cliques by hold/pivot paths and avoids explicit enumeration. | Definition `def:path`; Lemma `lem:supp`; citation `NuclearCD`. | Supported. |
| General CPI peeling needs mutable residual state for \(r\ge2\). | Baseline description in Preliminaries; citation `NuclearCD`. | Supported as prior-work characterization; avoid overclaiming beyond cited implementation. |
| At \(r=1\), vertex peeling preserves CPI path structure. | Theorem `thm:vertex-removal`. | Supported by proof. |
| A hold removal kills a path contribution, while a pivot removal decrements an effective pivot counter. | Theorem `thm:vertex-removal`; Figure `fig:overview`; Example `ex:peel`. | Supported. |
| Support changes have closed-form binomial deltas. | Lemma `lem:support-delta`; Figure `fig:peel`. | Supported. |
| LocalH computes the \((1,s)\)-core fixed point. | Lemma `lem:hindex-fixed-point`; Algorithm `alg:localH`. | Supported, but proof is concise and should be made more explicit. |
| PeelH computes the same result as LocalH. | Lemma `lem:equivalence`; Algorithm `alg:peelH`. | Supported, but proof needs careful invariant wording. |
| PeelH runs in \(O(\Sig+U+\kmax+n)\). | Theorem `thm:peelh-complexity`. | Supported if \(U\) is defined as total incremental refresh events. |
| LocalH costs \(O(T\Sig\log d_{\max})\). | Theorem `thm:localh-complexity`. | Supported for the presented refresh implementation. |
| ParaBuild is correct because seed-owned clique sets are disjoint. | Algorithm `alg:parabuild`; ParaBuild correctness paragraph; citations `BronKerbosch`, `eppstein`, `tomita`, `NuclearCD`. | Supported at high level; full BkEmit details are delegated to prior work. |
| BuildHier extracts the elder-rule hierarchy in \(O((n+m)\alpha(n))\). | Algorithm `alg:buildhier`; BuildHier complexity paragraph; citation `morozov2007persistence`. | Supported. |
| All reported implementations agree with the mutable baseline. | Experimental methodology text says all 762 agree; local `figures/benchmark_all_results.csv` supports 505 R1 cells on three graphs. | Partially verified locally; add TODO for 762 source. |
| End-to-end speedup is \(7.42\times\) geomean over 762 configurations. | Current Experimental table; figures generated from external `paper_data` path. Local CSV supports 505 configs with geomean \(6.16\times\). | Use with TODO verification against source data. |
| Peak memory reduced up to \(3.25\times\). | Current Experimental tables; local CSV supports max \(2.21\times\) for 505 subset; breakdown table includes \(3.33\times\) for com-youtube. | Use carefully with TODO or table-specific scope. |
| ParaBuild reaches \(10.0\times\) geomean speedup at 24 threads. | Table `tab:par`; existing figure scripts mention scalability data outside repo. | Supported by current table; source not local. |
| LocalH agrees with PeelH but is slower in current single-machine implementation. | Experimental subsection `sec:exp-local`. | Supported by current text; source data not local. |
| Stress test keeps V3 below REF where both finish. | Figures `fig_stress_time`, `fig_stress_mem`; Experimental text. | Supported by existing figures; raw source not local. |
| Case studies show \((1,s)\) matches higher-\(r\) quality at lower cost. | Tables `tab:cs4`, `tab:cs7`; figures `cs6_pareto`, `cs7_pareto`; CaseStudies text. | Supported by current tables; avoid claiming universal beyond evaluated tasks. |
