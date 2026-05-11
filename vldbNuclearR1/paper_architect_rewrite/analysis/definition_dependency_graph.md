# Definition Dependency Graph

Graph notation:
\(G=(V,E)\), induced subgraph \(G[X]\), neighbors \(N_G(v)\), CSR storage.

\(k\)-clique:
Requires graph notation.
Used by \(s\)-clique support and nucleus definitions.

\(s\)-clique support \(\sdeg(v)\):
Requires \(k\)-clique.
Counts \(s\)-cliques containing \(v\).
Used by \((1,s)\)-nucleus, core number, baseline peeling, and initial \(h\)-values.

\((1,s)\)-nucleus:
Requires \(s\)-clique support and connectivity through \(s\)-cliques.
Defines the object whose vertex core numbers are computed.

\(\kappa_s(v)\):
Requires \((1,s)\)-nucleus.
Largest \(k\) for which \(v\) belongs to a \(k\)-\((1,s)\)-nucleus.

CPI path:
Requires clique and parameter \(s\).
A path \(\TPath=(\Vh(\TPath),\Vp(\TPath))\) has hold set and pivot set.

Hold set \(\Vh(\TPath)\):
Requires CPI path.
Vertices included in every \(s\)-clique encoded by the path.

Pivot set \(\Vp(\TPath)\):
Requires CPI path.
Vertices from which the remaining clique slots are chosen.

Target \(\eta_P\):
Requires CPI path.
\(\eta_P=s-|\Vh(P)|\), the number of pivot vertices selected per encoded \(s\)-clique.

CPI support formula:
Requires CPI path, hold set, pivot set, target.
Gives support contribution for holds and pivots using binomial coefficients.

Static CPI theorem:
Requires CPI path and vertex peeling.
States that hold removal kills a path contribution, while pivot removal decrements the effective pivot counter.

Support delta:
Requires static CPI theorem and support formula.
Computes the binomial support loss after one or more pivot removals or path death.

\(\mathcal{H}\)-operator:
Requires multiset of witness values.
Returns the largest \(h\) with at least \(h\) witnesses of value at least \(h\).

LocalH:
Requires \(\mathcal{H}\)-operator, CPI witness count, and initial supports.
Computes the fixed point synchronously.

PeelH:
Requires static CPI theorem, support delta, bucket queue, and LocalH fixed-point semantics.
Computes the same core values in event order.

ParaBuild:
Requires CPI construction, degeneracy ordering, seed ownership, static CPI consumer layout.
Builds the static CPI in parallel with thread-local arenas and deterministic merge.

Deletion log:
Requires PeelH.
Sequence of removed vertices and assigned \(\kappa_s\), sorted by nondecreasing core value.

BuildHier:
Requires deletion log, input graph adjacency, and union-find.
Scans the log in reverse to build the elder-rule join tree.
