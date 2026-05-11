# Notation Budget

## Core Notation

\(G=(V,E)\):
Keep.
Needed everywhere.

\(s\):
Keep.
Primary clique size parameter.

\(\sdeg(v)\):
Keep.
Use only for \(s\)-clique support.

\(\kappa_s(v)\):
Keep.
Primary output.

\(\mathcal{P}\):
Keep.
Set of CPI paths.

\(\TPath\):
Keep.
Single CPI path.

\(\Vh(\TPath)\), \(\Vp(\TPath)\):
Keep.
Use consistently as hold and pivot sets.

\(\eta_{\TPath}\):
Keep.
Target pivot count.

\(p_{\TPath}\):
Keep.
Dynamic effective pivot counter in static CPI.

\(\Sig\):
Keep.
Vertex--path incidence count and primary index-size parameter.

## Secondary Notation

\(U\):
Keep in HIndex only.
Total number of incremental refresh events.

\(\kmax\):
Keep for bucket complexity.

\(d_{\max}\):
Keep only for LocalH complexity.

\(L\):
Use carefully.
It means witness multiset in HIndex and deletion log in BuildHier.
Prefer \(M_v\) for witness multiset if rewriting reduces overload.

\(V_{\text{fin}}\):
Keep inside PeelH pseudocode only.

\(\Delta_{\hold}, \Delta_{\pivot}\):
Keep in StaticCPI and HIndex.
These are useful because they bind the theorem to the implementation.

## Overloaded Or Risky Notation

\(k\):
It means nucleus level, \(h\)-threshold, and loop cursor.
This is standard but should be locally reintroduced in each context.

\(\mathcal{S}\):
Currently used for sets of admissible cliques.
Keep only in the nucleus definition.
Do not reuse later.

\(P\):
Avoid as an informal path variable because \(\Vp\) already uses \(P\)-style notation.
Use \(\TPath\) in prose.

## Notation To Remove Or Delay

\(d_{\max}\):
Delay until LocalH cost.

CSR arrays:
Mention once in preliminaries or ParaBuild.
Do not make them part of the formal model.

Generation-stamped bitmap:
Mention in ParaBuild implementation notes only.

SCT:
Define only in RelatedWork or ParaBuild if the paper needs to contrast with Pivoter.
