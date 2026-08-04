# FROZEN LEDGER — the NSI plane-index paper (per paper-architect A.2 + section-writing-protocol)

Every section drafts AGAINST this file; edits to shared facts happen HERE first.
Architecture: SigmodPlus.md §237 (spine/story) + §239 (skeleton). Sentences: §238.

## Names
- Index brand: **NSI** (\nsi, "Nucleus Spectrum Index") — the one branded artifact.
- Competitor: **CND** (\cnd, \cite{NuclearCD}) — the group's published PACMMOD'26 counting-based
  algorithm; serial protocol only (§217 quarantine).
- Helpers (never versioned): \peelcell, \buildalg, \queryalg.
- New-term budget (4+1): class, pattern, certificate, residue, **box** (spent deliberately in
  sec:peelcell -- the s-side needs one container noun; SCT stays out of the main line).
  "clique quotient" = the umbrella concept naming the grouping (introduced in §1, defined in
  sec:quotient).
- NOT in the paper's main line: SCT, "region" as a term (the RegND draft uses it; WE say "maximal
  clique" — budget rule), NSI2/3/4 format names.
- DECIDED notation (group conventions inherited from NuclearCD/RegND):
  * kappa_s(R) — core number; r = |R| is implicit, so one glyph serves the whole plane.
  * deg_{(s)}(R) — s-clique support.
  * M(G) — maximal cliques; omega(R) = size of the largest clique of G containing R
    ("local clique number"; standard omega glyph localized — the spine's one number).
  * Plane P: r in [rmin,rmax], r < s <= smax. Query contract: (R, s) with (|R|, s) in P.
  * Inherited vocabulary (group's published usage, not part of the 4-term budget):
    "scored clique" (= r-clique), "witness clique" (= s-clique).
  * sigma(P) — the certification point: first s with kappa at the floor (Absorbing Chain
    locks all later s; always sigma(P) <= omega(P)+1). Storage contract: store P iff
    sigma(P) > r+1, with its residue values. Defined in sec:index.
- Exact object counts, web-it-2004 (archive accounting, integer-exact):
  3-cliques 338,786,461; 4-cliques 28,530,662,583; 5-cliques 2,140,698,950,272.

## The spine sentence (paper = this, unfolded)
For almost every r-clique R, kappa_s(R) = C(omega(R)-r, s-r) where omega(R) is one number: the size
of the largest clique containing R. The paper certifies where this holds, computes the
residue where it fails, and stores only the quotient + residue (P10: unavoidable).

## Running example (verified by two independent scripts; figures/fig1_example.pdf)
- G: 10 vertices v1..v10; maximal cliques M1={v1..v6}(6), M2={v5..v9}(5), M3={v7..v10}(4).
- Classes: A={v1..v4}[M1], B={v5,v6}[M1,M2], C={v7,v8,v9}[M2,M3], D={v10}[M3].
- Objects: 33 triangles, 21 4-cliques, 7 5-cliques; 61 cliques of size 3..5 total.
- Patterns: r=3 → 7; r=4 → 6; r=5 → 3; total 16 (compression 33/7=4.7x at r=3).
- (3,4) values: AAA/AAB/ABB → 3; BBC/BCC/CCC → 2; CCD → 1. Peel has 3 rounds and ONE cascade
  (CCC support 3 → 2 when CCD dies). All values land on the clique floor.
- Cross-cell: (3,5) {v1,v2,v3}=3; (4,5) {v1,v2,v3,v4}=2. Plane r=3..5,s≤8: 12 cells, ONE peel.

## Headline numbers (single-run; refresh after E2 medians — do NOT requote elsewhere without
## checking this table first)
| claim | number | source |
|---|---|---|
| 5-cliques on web-it-2004 (7.2M edges) | >1.5e12 | §233 archive accounting |
| CND serial, web-it (3,4) | 5,581 s / 101 GB | §216 |
| ours vs CND, web-it (3,4) | 742x time / 226x memory | §216 |
| CND on web-it grid (1800s/300GB) | 0/9 cells | §219 |
| grid vs serial CND, 8 graphs/4 domains | no losing cell | §218–219 |
| explicit plane table, web-it | **50.8 TB decimal** (log prints MiB-as-MB: 46.2 TiB); unbuildable r≥4 | §233 |
| NSI size, web-it whole plane | **1.70 MB decimal** (1,700,530 B; earlier "1.66 MB" was a units slip) | §234 |
| UNIT RULE | the engine logs print MiB/GiB as "MB/GB"; the PAPER uses decimal SI everywhere -- convert at transcription | this file |
| point query | 74 ns (web-it r=3 kernel) | §233 |
| archive/index size ratio range | 14x – 3.0e7 (ratios from BYTES, unit-safe) | §233–234 |
| grid protocol split | standalone matched cell: 742x/226x (web-it (3,4)); amortized grid: 42/42 won, up to 5,806x (DB (5,8)), 30/72 CND-infeasible | §216/§218–219 |
| certificates ablation | no-T3 22.6x (cells 67–5406x); no-T5 2.4–3.0x; answers identical | §236 |
| index bytes that are the quotient | 70.5–99.9% | §234 corrected |
| residue bytes | 0.0–0.1% | §221 |
| soc-pokec compression (predicted negative) | 1.001x | §220/D10 |
| compression endpoints quoted in sec:index | nasasrb 8.3x, web-it 515x (do NOT claim a full-roster range without per-graph numbers) | §219/§220 |
| pattern collapse, real graphs | 84.9x – 1,167,073x | §233 |
| example collapse | 33 → 7 | ledger above |
| query kernel false positives | 0 / 165k near-miss + 700k random | §231 (body only; NOT a
  correctness experiment — correctness is by lemma) |

## Section-5 worked numbers (verified by script 2026-08-03)
- Toy boxes: B1={A,B} free; B2={B,C} with ell_C>=1; B3={C,D} with ell_D>=1. Exact partition at
  every size. s=4 counts 15/5/1 (sum 21 ✓); s=5 counts 6/1/0 (sum 7 ✓).
- sup(A^2B)@(3,4)=3 on B1. sup(C^3)=2(B2)+1(B3)=3; after C^2D removed, antichain {(2,1)} on B3
  kills its contribution -> 2 (the cascade). sup(B^2C)=2.
- pkustk13 (3,4): 2.2M patterns vs 45.2M cliques; r=5: 26.8M vs 2.28B.

## Case-study numbers (frozen 2026-08-04; all STRUCTURAL, machine-independent, local runs)
- com-dblp ladder: papers of 114/102/65/63 authors; kappa_4 = 111/99/62/60; kappa_8 =
  128,164,707 / 71,523,144 / 6,471,002 / 5,461,512. Counts C(114,3)=240,464, C(102,3)=171,700,
  2,016 = C(65,3)-C(64,3) (the 64-shared-author discovery), C(63,3)=39,711.
- com-dblp retrieval (3,4): k=111 -> 1 comp (240,464); k=99 -> 2 comps; k=62 -> 65-paper MERGES
  (242,480 = 240,464+2,016); k=60 -> third separate comp 39,711.
- (3,8) nuclei retrieval on mega-cliques is witness-exponential (C(114,8) 8-cliques) -- do NOT
  attempt; connectivity claims at high s stay reasoned, not run.
- han1 ego (built from name-based coauthor.tsv; n=1,358 == cs9 meta cross-check; m=21,699;
  297,899 triangles; anchor id 484 = Jiawei Han 0001): top level (3,4)=93 (96-author paper,
  C(96,3)=142,880 ✓). Residue: 45 @ (3,5) == 45 flagged by the kappa5 != C(kappa4,2) detector
  (exact match with engine); 53 researchers; thins 45/45/38/9 at s=5..8; only survivor above
  floor AT s=8: Aggarwal-Pei-Liu (omega=8, k4=6, k5=14, k8=2 vs floor 1).
- Stored records across the whole han1 plane: **45, all in row 3** (r=4's 196 and r=5's 1 are
  boundary-diagonal residue that lands ON floor -> sigma=r+2 -> zero records; verified via
  (4,6)/(5,7) residue = 0).
- Exp numbering: case studies are Exp-8 and Exp-9.

## Claim → evidence map (contributions)
1. First index over the plane → reconstruction lemma (sec:index) + §233/§234 sizes + 74ns.
2. Quotient framework + certified transfer → interchangeability thm, coefficient lemma, T3, T5
   (sec:quotient..sec:plane) + E7 ablation §236.
3. Applicability test in advance → W definition + prediction table (sec:index) + §220/D10.
4. Experiments → grid (§218-219), plane build (§235), size/query (§233-234), ablations (§236).

## Style gates (every section, before it counts as drafted)
1. sentence detector: target 10–22, cap 28; licensed >28 only for semicolon citation sweeps,
   the research-question hinge, and line-anchored walkthroughs.
2. build: latexmk clean AND `pdftotext` shows 0 "??" AND 0 "[N?]" citations AND
   `grep -c "undefined" main.log` == 0 (a missing bib key renders as [22?], which the
   plain ?? grep does NOT catch -- learned 2026-08-04).
3. prose-cadence pass + teacher_vocab_diff (register gate) at polish time.
4. No symbols before sec:preliminaries; no correctness-as-experiment claims (skill HC-12);
   no self-exposed weaknesses (advisor rule) — applicability is a capability, framed predictive.

## Venue
SIGMOD 2027 rounds: Jan 17 / Apr 17 / Jul 17 / Oct 17. Target: **2026-10-17** (next open round).
Re-verify on the official page before submission (site was behind a bot-check when checked
2026-08-03; dates from the CFP search snippet + chairs' announcement).
