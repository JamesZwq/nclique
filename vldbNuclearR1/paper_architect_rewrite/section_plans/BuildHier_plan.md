# BuildHier Plan

Section purpose:
Show hierarchy as a natural byproduct of deletion order.

Reader knowledge before:
PeelH emits vertices in nondecreasing core value.

Reader knowledge after:
The reader knows how a reverse scan of the sorted deletion log constructs the join tree.

Paragraph contracts:
1. Give intuition for super-level components.
2. Define deletion log.
3. Present BuildHier state machine and algorithm.
4. State correctness.
5. State complexity and boundary.

Key definitions introduced:
Deletion log, super-level filtration, elder-rule join tree.

Claims made:
No sort is needed; reverse scan costs \(O((n+m)\alpha(n))\); CPI is not consulted.

Evidence used:
Algorithm `alg:buildhier`, citation `morozov2007persistence`.

Sentence-level risks:
Do not assume the reader knows elder-rule terminology before intuition.
