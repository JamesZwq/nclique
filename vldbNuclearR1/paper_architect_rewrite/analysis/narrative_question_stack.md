# Narrative Question Stack

## Whole Paper

1. Why do vertex rankings need higher-order witnesses beyond edges?
2. Why is \((1,s)\)-nucleus decomposition the right target object?
3. Why is exact high-\(s\) decomposition expensive under the current CPI framework?
4. What extra structure appears when \(r=1\)?
5. What state remains dynamic if CPI paths are never modified?
6. How can core values be computed from those counters?
7. How do we avoid the redundancy of synchronous local \(h\)-index iteration?
8. Once peeling is cheap, what becomes the bottleneck?
9. Can the same deletion order also produce the hierarchy?
10. Does the method remain exact and faster on real graphs?
11. Is the special \(r=1\) case practically useful?

## Section Questions

Introduction:
What is the central tension, what invariant resolves it, and what pipeline follows?

Preliminaries:
What exactly is a \((1,s)\)-nucleus, what does CPI encode, and what does the mutable baseline maintain?

StaticCPI:
Why does vertex peeling preserve CPI path structure, and why does this fail for \(r\ge2\)?

HIndex:
Given a static CPI, how do LocalH and PeelH compute the same core values?

ParaBuild:
Why does construction dominate after PeelH, and how can construction be parallelized without changing the index semantics?

BuildHier:
Why is the hierarchy already latent in the deletion log?

Experimental:
Which measured claims are supported: exactness, speed, memory, bottleneck shift, parallel scaling, LocalH semantics, stress behavior?

CaseStudies:
Does \((1,s)\) give useful vertex rankings and structures compared with higher-\(r\) nuclei?

RelatedWork:
Which prior lines solve adjacent problems, and what limitation leaves room for static CPI?

Conclusion:
What changed in the reader's model of CPI for \(r=1\)?
