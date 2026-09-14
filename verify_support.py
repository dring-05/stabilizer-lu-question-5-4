#!/usr/bin/env python3
"""Exact support-subgroup audit for the stabilizer pairs in Question 5.4."""

from collections import Counter
from itertools import combinations, permutations


PAULI_BITS = {"I": (0, 0), "X": (1, 0), "Z": (0, 1), "Y": (1, 1)}
BITS_PAULI = {value: key for key, value in PAULI_BITS.items()}

CASES = {
    "five": {
        "left": ("ZZZZI", "XXIIX", "IIXXX"),
        "right": ("ZZIZX", "XIZXI", "IXXYI"),
        "subset_size": 3,
    },
    "six": {
        "left": ("XXIIXX", "IIXXXX", "ZZZZIX"),
        "right": ("XXZXII", "ZZXYII", "IZZZXX"),
        "subset_size": 4,
    },
}


def symplectic_inner(a, b):
    total = 0
    for pa, pb in zip(a, b):
        xa, za = PAULI_BITS[pa]
        xb, zb = PAULI_BITS[pb]
        total ^= (xa & zb) ^ (za & xb)
    return total


def multiply_mod_phase(a, b):
    out = []
    for pa, pb in zip(a, b):
        xa, za = PAULI_BITS[pa]
        xb, zb = PAULI_BITS[pb]
        out.append(BITS_PAULI[(xa ^ xb, za ^ zb)])
    return "".join(out)


def stabilizer_group(generators):
    n = len(generators[0])
    assert all(len(g) == n for g in generators)
    assert all(symplectic_inner(a, b) == 0 for a, b in combinations(generators, 2))
    group = []
    for mask in range(1 << len(generators)):
        word = "I" * n
        for i, generator in enumerate(generators):
            if mask >> i & 1:
                word = multiply_mod_phase(word, generator)
        group.append(word)
    assert len(set(group)) == 1 << len(generators)
    return tuple(group)


def support_mask(word):
    return sum((ch != "I") << i for i, ch in enumerate(word))


def subset_rank_profile(generators, subset_size):
    group = stabilizer_group(generators)
    n = len(generators[0])
    profile = {}
    for subset in combinations(range(n), subset_size):
        mask = sum(1 << i for i in subset)
        subgroup_size = sum((support_mask(word) & ~mask) == 0 for word in group)
        assert subgroup_size > 0 and subgroup_size & (subgroup_size - 1) == 0
        profile[subset] = (1 << subset_size) // subgroup_size
    return profile


def permute_subset(subset, permutation):
    return tuple(sorted(permutation[i] for i in subset))


def profiles_match_under_permutation(left, right, n):
    witnesses = []
    for permutation in permutations(range(n)):
        if all(left[subset] == right[permute_subset(subset, permutation)] for subset in left):
            witnesses.append(permutation)
    return witnesses


def permute_word(word, permutation):
    out = [None] * len(word)
    for old, new in enumerate(permutation):
        out[new] = word[old]
    return "".join(out)


def group_record(generators):
    records = []
    for word in stabilizer_group(generators):
        support = tuple(i + 1 for i, ch in enumerate(word) if ch != "I")
        records.append((len(support), support, word))
    return sorted(records)


def main():
    for name, case in CASES.items():
        left_generators = case["left"]
        right_generators = case["right"]
        n = len(left_generators[0])
        k = case["subset_size"]
        left = subset_rank_profile(left_generators, k)
        right = subset_rank_profile(right_generators, k)
        matches = profiles_match_under_permutation(left, right, n)

        # Acceptance control: a relabeled copy must be recognized as equivalent.
        control_permutation = tuple(reversed(range(n)))
        relabeled = tuple(permute_word(word, control_permutation) for word in left_generators)
        relabeled_profile = subset_rank_profile(relabeled, k)
        control_matches = profiles_match_under_permutation(left, relabeled_profile, n)
        assert control_matches

        print(f"CASE {name}: n={n}, subset_size={k}")
        print("LEFT GROUP")
        for weight, support, word in group_record(left_generators):
            print(f"  {word} weight={weight} support={support}")
        print("RIGHT GROUP")
        for weight, support, word in group_record(right_generators):
            print(f"  {word} weight={weight} support={support}")
        print(f"left reduced-rank histogram: {dict(sorted(Counter(left.values()).items()))}")
        print(f"right reduced-rank histogram: {dict(sorted(Counter(right.values()).items()))}")
        print(f"permutations preserving the cross-profile: {len(matches)}")
        print(f"acceptance-control permutations: {len(control_matches)}")
        assert not matches

        if name == "five":
            left_exceptional = [set(i + 1 for i in subset) for subset, rank in left.items() if rank == 4]
            right_exceptional = [set(i + 1 for i in subset) for subset, rank in right.items() if rank == 4]
            assert len(left_exceptional) == len(right_exceptional) == 2
            left_intersection = len(left_exceptional[0] & left_exceptional[1])
            right_intersection = len(right_exceptional[0] & right_exceptional[1])
            print(f"exceptional-triple intersection sizes: left={left_intersection}, right={right_intersection}")
            assert (left_intersection, right_intersection) == (1, 2)
        else:
            assert Counter(left.values()) == Counter({16: 12, 8: 3})
            assert Counter(right.values()) == Counter({16: 14, 4: 1})
        print("RESULT: reduced-rank profiles obstruct local-unitary equivalence under every qubit permutation")


if __name__ == "__main__":
    main()
