# One-Sentence Thesis

For \((1,s)\)-nucleus decomposition, vertex peeling preserves the Clique Path Index structure, so exact decomposition can run over a static path index with only counter updates, which turns mutable residual-index maintenance into a cheaper static-index pipeline for peeling, construction, and hierarchy extraction.

# How Each Section Supports It

Introduction:
Introduce the tension between general CPI power and mutable residual state, then state the \(r=1\) invariant and derive the pipeline.

Preliminaries:
Define the target object, CPI paths, support formula, and mutable baseline so the reader can see exactly what state the new method removes.

StaticCPI:
Prove the invariant that CPI paths stay structurally valid under vertex peeling and reduce dynamic state to counters.

HIndex:
Show that once CPI is static, core values can be computed as an \(h\)-index fixed point and then as an event-driven peel using the same counter semantics.

ParaBuild:
Show that cheap peeling shifts the bottleneck to index construction, and that static output lets construction stream paths into thread-local arenas.

BuildHier:
Show that the peel order already contains the sorted filtration needed for hierarchy construction.

Experimental:
Test correctness against the mutable baseline, measure whether the static pipeline is faster and smaller, and isolate where the gains come from.

CaseStudies:
Defend why the vertex-centric \((1,s)\) case is worth special treatment by showing practical quality against higher-\(r\) alternatives.

RelatedWork:
Position the work as neither a generic CPI replacement nor an enumeration method, but as an exact \(r=1\) specialization enabled by a structural invariant.

Conclusion:
Return to the change in state model: CPI becomes static symbolic state instead of mutable residual state.
