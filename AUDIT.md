# Audit ledger

Audit date: September 14, 2026

## Problem record

- Primary source: https://arxiv.org/html/2607.26214v1
- Abstract and version history: https://arxiv.org/abs/2607.26214
- Title: *Minimal Counterexamples of the MacWilliams Extension Theorem for Stabilizer Codes*
- Author: Ali Assem Mahmoud
- Status: arXiv preprint, v1 submitted July 28, 2026
- Exact location: Theorem 5.3 and Question 5.4
- Targeted claim: decide product-unitary equivalence up to qubit permutation for both concrete pairs in the first sentence of Question 5.4

## Evidence ledger

| Claim | Evidence | Independent check | Status |
|---|---|---|---|
| The two pairs are posed as open | Primary source, Question 5.4 and Section 6.2 | Source-auditor reread of version history and quantifiers | Passed |
| Stabilizer elements and supports | Exhaustive binary symplectic enumeration | Exact Pauli matrices including phases | Passed |
| Five-qubit obstruction | Rank-4 triples intersect in 1 versus 2 qubits | Direct exact partial traces of 32 by 32 projectors | Passed |
| Six-qubit obstruction | Four-subset rank histograms differ | Direct exact partial traces of 64 by 64 projectors | Passed |
| All permutations excluded | 5! and 6! exhaustive profile checks | Hand invariant for five qubits, histogram mismatch for six | Passed |
| Checker is not vacuous | Target pairs rejected | Reversed-qubit equivalent copies accepted | Passed |
| Novelty status | Title, identifier, generator, Question 5.4, and invariant searches | Independent source-auditor search | No later resolution found |

## Novelty search scope

Searches used the exact paper title, arXiv identifier, `Question 5.4`, both complete generator triples, combinations of generator strings, `weight-isometric` with local-unitary terminology, `reduced-rank`, `support invariant`, and GitHub-restricted variants. They returned the source preprint, bibliographic mirrors, and broader prior literature, but no later resolution, revision, repository, or discussion deciding either pair.

The honest novelty statement is: no later resolution was found in targeted searches. Indexing latency and unpublished work remain possible.

## Prior art boundary

The partial-trace rank formula is established stabilizer formalism. Relevant sources include:

- Zeng, Cross, and Chuang, *Transversality versus Universality for Additive Quantum Codes*: https://arxiv.org/abs/0706.1382
- Audenaert and Plenio, *Entanglement on mixed stabiliser states, I: Normal Forms and Reduction Procedures*: https://arxiv.org/abs/quant-ph/0505036
- Pllaha, *Symplectic Isometries of Stabilizer Codes*: https://arxiv.org/abs/1807.09107

The apparent new contribution is the permutation-complete application of the known invariant to both exact pairs in Question 5.4.

## Verification environment

- Python 3.14.3
- Darwin 25.6.0 arm64
- Standard library only
- Exact integer and rational arithmetic
- No randomness and no numerical tolerance

See `REPLAY.txt` for the clean run output.
