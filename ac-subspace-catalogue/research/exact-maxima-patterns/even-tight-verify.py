#!/usr/bin/env python3
"""Exact verification of the even-order pure AC equality obstruction.

Requires SymPy. No floating-point tolerances are used.
The r=3 proof uses physical rectangular 3j matrices; a separate 6j sum
checks the recoupling formula and selected further equality obstructions.
"""
import json
from pathlib import Path
import sympy as sp
from sympy.physics.wigner import wigner_3j, wigner_6j


def physical_channel_check(r):
    n = (r + 1) ** 2
    j, target = sp.Rational(n, 2), sp.Rational(n - 2, 2)
    rank = 2 * r + 1
    inputs = [-j + i for i in range(n + 1)]
    outputs = [-target + i for i in range(n - 1)]
    input_tensor = {
        m: (-1) ** (j - m) * sp.sqrt(2 * rank + 1)
        * wigner_3j(j, rank, j, -m, 0, m) for m in inputs
    }
    output_tensor = {
        a: (-1) ** (target - a) * sp.sqrt(2 * rank + 1)
        * wigner_3j(target, rank, target, -a, 0, a) for a in outputs
    }
    contributions, vectors = [], []
    total_output = {a: sp.S.Zero for a in outputs}
    for ell in range(1, r + 1):
        channel_output = {a: sp.S.Zero for a in outputs}
        for q in range(-ell, ell + 1):
            entries = {}
            for a in outputs:
                m = a - q
                if -j <= m <= j:
                    e = (-1) ** (target - a) * sp.sqrt(2 * ell + 1)
                    e *= wigner_3j(target, ell, j, -a, q, m)
                    entries[(a, m)] = e
                    channel_output[a] += e * e * input_tensor[m]
            # The displayed rectangular operators have HS norm one.
            assert sp.simplify(sum(e * e for e in entries.values())) == 1
            vectors.append(entries)
        contribution = sp.simplify(sum(
            output_tensor[a] * channel_output[a] for a in outputs
        ))
        contributions.append(contribution)
        for a in outputs:
            total_output[a] += channel_output[a]
    # Full rectangular-basis HS orthogonality, including distinct ranks.
    for i, a in enumerate(vectors):
        for b in vectors[:i]:
            assert sp.simplify(sum(v * b.get(key, 0)
                                   for key, v in a.items())) == 0
    multiplier = sp.simplify(sum(contributions))
    assert multiplier != 0
    assert all(sp.simplify(total_output[a]
                          - multiplier * output_tensor[a]) == 0
               for a in outputs)
    return multiplier, contributions


def recoupling_multiplier(r):
    n = (r + 1) ** 2
    j, target = sp.Rational(n, 2), sp.Rational(n - 2, 2)
    rank = 2 * r + 1
    return sp.simplify(sum(
        (-1) ** (n + rank + ell + 1) * (2 * ell + 1)
        * wigner_6j(j, j, rank, target, target, ell)
        for ell in range(1, r + 1)
    ))


def main():
    results = []
    for r in range(2, 9):
        multiplier = recoupling_multiplier(r)
        assert multiplier != 0
        entry = {"r": r, "n": (r + 1) ** 2, "t": 2 * r,
                 "next_rank_multiplier": str(multiplier),
                 "pure_state_exists": False}
        if r in (2, 3):
            physical, contributions = physical_channel_check(r)
            assert sp.simplify(physical - multiplier) == 0
            entry["physical_channel_independently_checked"] = True
            entry["rank_contributions"] = [str(v) for v in contributions]
        results.append(entry)
    assert recoupling_multiplier(3) == 5 * sp.sqrt(69) / 8398
    report = json.dumps({"exact_checks_passed": True, "results": results}, indent=2) + "\n"
    Path(__file__).with_name("even-tight-verification.json").write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
