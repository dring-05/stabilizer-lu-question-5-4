# Reproducible audit of Question 5.4

## Source

- Ali Assem Mahmoud, *Minimal Counterexamples of the MacWilliams Extension Theorem for Stabilizer Codes*
- arXiv:2607.26214v1, submitted July 28, 2026
- Preprint status as checked September 14, 2026. The arXiv history lists only v1. The HTML manuscript carries an internal August 24 date, but there is no v2 in the arXiv history.
- Exact location: Theorem 5.3 gives the generators, Question 5.4 asks whether either displayed pair is related by arbitrary product unitaries and a qubit permutation
- Primary source: https://arxiv.org/html/2607.26214v1

## Outcome

Solved for both named pairs. Their reduced-projector rank profiles are incompatible under every qubit permutation, so neither pair is locally unitarily equivalent up to permutation.

See `PROOF.md` for the mathematical argument.

## Replay

From this directory, run:

```text
python3 verify_support.py
python3 verify_matrices.py
```

Expected final lines:

```text
RESULT: reduced-rank profiles obstruct local-unitary equivalence under every qubit permutation
RESULT: independent exact matrix ranks reproduce both local-unitary obstructions
```

`verify_support.py` uses binary symplectic Pauli labels, enumerates every stabilizer element, computes the subgroup-supported rank profile, and checks every qubit permutation. It also checks a known relabeled equivalent copy as an acceptance control.

`verify_matrices.py` constructs the full stabilizer projector numerator from exact Pauli matrices over Gaussian rational numbers, performs partial traces directly, and computes ranks by exact elimination. It checks Hermiticity, the projector identity, the expected exceptional subsets, and a relabeled-copy acceptance control.

No floating point arithmetic, numerical tolerance, randomness, optimization solver, or external Python package is used.

## Trusted base

- Python 3 standard library
- elementary Pauli matrix multiplication
- exact rational arithmetic supplied by `fractions.Fraction`
- the theorem that unitary conjugation preserves matrix rank

The two programs share the printed generators and expected claim but implement the decisive calculation differently. The short proof is independently inspectable without either program.

## Novelty and limitations

A targeted search on September 14, 2026 using the paper title, arXiv identifier, Question 5.4, exact generator strings, and the reduced-rank invariant found no later resolution. This supports the statement that the result appears new; it is not proof that no unpublished or unindexed solution exists.

The reduced-rank formula is known in stabilizer theory. The apparent new result is the application that separates both pairs in Question 5.4. The broader LU versus LC existence question remains open. These codes have distance 1, as the source notes, so this does not resolve the paper's separate distance and purity question.

Relevant prior art includes Zeng, Cross, and Chuang on support-defined partial projectors for additive quantum codes (https://arxiv.org/abs/0706.1382), and Audenaert and Plenio on partial traces of mixed stabilizer states and code projectors (https://arxiv.org/abs/quant-ph/0505036). The proof here does not claim the rank invariant itself as new.

## Provenance and publication status

This artifact was produced by Codex using independent source-auditor, solver, and skeptic roles, followed by a main-agent derivation and clean replay. AI assistance should be disclosed in any publication or correspondence.

Published with the user's explicit approval at https://github.com/dring-05/stabilizer-lu-question-5-4 on September 14, 2026. It has not been sent to the paper's author or submitted to a journal.
