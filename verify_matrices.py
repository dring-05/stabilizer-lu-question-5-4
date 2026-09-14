#!/usr/bin/env python3
"""Independent exact matrix audit of the reduced-projector rank obstruction."""

from collections import Counter
from fractions import Fraction
from itertools import combinations


ZERO = (Fraction(0), Fraction(0))
ONE = (Fraction(1), Fraction(0))
I = (Fraction(0), Fraction(1))

PAULI = {
    "I": [[ONE, ZERO], [ZERO, ONE]],
    "X": [[ZERO, ONE], [ONE, ZERO]],
    "Y": [[ZERO, (Fraction(0), Fraction(-1))], [I, ZERO]],
    "Z": [[ONE, ZERO], [ZERO, (Fraction(-1), Fraction(0))]],
}

CASES = {
    "five_left": ("ZZZZI", "XXIIX", "IIXXX"),
    "five_right": ("ZZIZX", "XIZXI", "IXXYI"),
    "six_left": ("XXIIXX", "IIXXXX", "ZZZZIX"),
    "six_right": ("XXZXII", "ZZXYII", "IZZZXX"),
}


def qadd(a, b):
    return a[0] + b[0], a[1] + b[1]


def qsub(a, b):
    return a[0] - b[0], a[1] - b[1]


def qmul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def qdiv(a, b):
    denominator = b[0] * b[0] + b[1] * b[1]
    assert denominator
    return (
        (a[0] * b[0] + a[1] * b[1]) / denominator,
        (a[1] * b[0] - a[0] * b[1]) / denominator,
    )


def is_zero(a):
    return a == ZERO


def kron(a, b):
    ar, ac = len(a), len(a[0])
    br, bc = len(b), len(b[0])
    out = [[ZERO for _ in range(ac * bc)] for _ in range(ar * br)]
    for i in range(ar):
        for j in range(ac):
            for p in range(br):
                for q in range(bc):
                    out[i * br + p][j * bc + q] = qmul(a[i][j], b[p][q])
    return out


def pauli_matrix(word):
    matrix = [[ONE]]
    for ch in word:
        matrix = kron(matrix, PAULI[ch])
    return matrix


def matmul(a, b):
    rows, inner, cols = len(a), len(b), len(b[0])
    out = [[ZERO for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for k in range(inner):
            if is_zero(a[i][k]):
                continue
            for j in range(cols):
                if not is_zero(b[k][j]):
                    out[i][j] = qadd(out[i][j], qmul(a[i][k], b[k][j]))
    return out


def matscale(scalar, matrix):
    return [[qmul(scalar, value) for value in row] for row in matrix]


def conjugate_transpose(matrix):
    return [[(matrix[j][i][0], -matrix[j][i][1]) for j in range(len(matrix))] for i in range(len(matrix[0]))]


def identity(size):
    return [[ONE if i == j else ZERO for j in range(size)] for i in range(size)]


def projector_numerator(generators):
    n = len(generators[0])
    generator_matrices = [pauli_matrix(word) for word in generators]
    group_matrices = []
    for mask in range(1 << len(generators)):
        matrix = identity(1 << n)
        for i, generator in enumerate(generator_matrices):
            if mask >> i & 1:
                matrix = matmul(matrix, generator)
        group_matrices.append(matrix)
    out = [[ZERO for _ in range(1 << n)] for _ in range(1 << n)]
    for matrix in group_matrices:
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                out[i][j] = qadd(out[i][j], value)
    return out


def embed_bits(kept_value, env_value, kept, traced, n):
    bits = [0] * n
    for offset, position in enumerate(reversed(kept)):
        bits[position] = kept_value >> offset & 1
    for offset, position in enumerate(reversed(traced)):
        bits[position] = env_value >> offset & 1
    value = 0
    for bit in bits:
        value = (value << 1) | bit
    return value


def partial_trace(matrix, kept, n):
    kept = tuple(kept)
    traced = tuple(i for i in range(n) if i not in kept)
    dimension = 1 << len(kept)
    env_dimension = 1 << len(traced)
    out = [[ZERO for _ in range(dimension)] for _ in range(dimension)]
    for i in range(dimension):
        for j in range(dimension):
            total = ZERO
            for env in range(env_dimension):
                row = embed_bits(i, env, kept, traced, n)
                col = embed_bits(j, env, kept, traced, n)
                total = qadd(total, matrix[row][col])
            out[i][j] = total
    return out


def exact_rank(matrix):
    work = [row[:] for row in matrix]
    rows, cols = len(work), len(work[0])
    rank = 0
    for col in range(cols):
        pivot = next((row for row in range(rank, rows) if not is_zero(work[row][col])), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][col]
        work[rank] = [qdiv(value, pivot_value) for value in work[rank]]
        for row in range(rows):
            if row == rank or is_zero(work[row][col]):
                continue
            factor = work[row][col]
            work[row] = [qsub(value, qmul(factor, base)) for value, base in zip(work[row], work[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def profile(generators, subset_size):
    n = len(generators[0])
    projector = projector_numerator(generators)
    group_size = 1 << len(generators)
    assert projector == conjugate_transpose(projector)
    assert matmul(projector, projector) == matscale((Fraction(group_size), Fraction(0)), projector)
    return {
        subset: exact_rank(partial_trace(projector, subset, n))
        for subset in combinations(range(n), subset_size)
    }


def main():
    profiles = {
        name: profile(generators, 3 if name.startswith("five") else 4)
        for name, generators in CASES.items()
    }
    expected = {
        "five_left": Counter({8: 8, 4: 2}),
        "five_right": Counter({8: 8, 4: 2}),
        "six_left": Counter({16: 12, 8: 3}),
        "six_right": Counter({16: 14, 4: 1}),
    }
    for name, ranks in profiles.items():
        histogram = Counter(ranks.values())
        print(f"{name}: {dict(sorted(histogram.items()))}")
        print("  exceptional subsets:", [(tuple(i + 1 for i in subset), rank) for subset, rank in ranks.items() if rank < max(ranks.values())])
        assert histogram == expected[name]

    # Acceptance control: a known qubit relabeling preserves the rank data.
    reversed_five_left = tuple(word[::-1] for word in CASES["five_left"])
    reversed_profile = profile(reversed_five_left, 3)
    assert Counter(reversed_profile.values()) == Counter(profiles["five_left"].values())
    reversed_exceptional = {
        subset for subset, rank in reversed_profile.items() if rank == 4
    }
    assert reversed_exceptional == {(0, 3, 4), (0, 1, 2)}

    five_left_sets = [set(subset) for subset, rank in profiles["five_left"].items() if rank == 4]
    five_right_sets = [set(subset) for subset, rank in profiles["five_right"].items() if rank == 4]
    assert len(five_left_sets[0] & five_left_sets[1]) == 1
    assert len(five_right_sets[0] & five_right_sets[1]) == 2
    assert Counter(profiles["six_left"].values()) != Counter(profiles["six_right"].values())
    print("acceptance control: exact matrix audit accepted a reversed-qubit copy")
    print("RESULT: independent exact matrix ranks reproduce both local-unitary obstructions")


if __name__ == "__main__":
    main()
