#!/usr/bin/env python3
"""Exact necessary-moment exclusion of C_m symmetry, m>=4, at N=27,t=6.

No third-party dependencies. This excludes a symmetry subclass only and does
not prove nonexistence of an arbitrary sixth-order anticoherent state.
"""
from fractions import Fraction
from math import prod
from pathlib import Path
import json

N, T = 27, 6

def exact_mean(values):
    return sum(values, Fraction()) / (N + 1)

def main():
    cyclic4 = []
    for a in range(4):
        support = list(range(a, N + 1, 4))
        assert len(support) == 7
        middle = support[3]
        other = [x for x in support if x != middle]
        denominator = prod(middle - x for x in other)
        assert denominator == -147456
        mean = exact_mean(Fraction(prod(r - x for x in other), denominator)
                          for r in range(N + 1))
        assert mean < 0
        cyclic4.append(dict(rotationOrder=4, residue=a, support=support,
                            cardinalIndex=3, uniformMean=str(mean),
                            reason='The required squared Dicke amplitude is negative.'))
    higher = []
    # m=28 represents every m>27: every residue support is a singleton.
    for m in range(5, N + 2):
        for a in range(min(m, N + 1)):
            support = list(range(a, N + 1, m))
            assert 1 <= len(support) <= T
            mean = exact_mean(Fraction(prod(r - x for x in support))
                              for r in range(N + 1))
            assert mean != 0
            higher.append(dict(rotationOrder=m, residue=a, support=support,
                               annihilatorDegree=len(support), uniformMean=str(mean)))
    assert len(higher) == 396
    report = dict(n=N, order=T, arithmetic='Exact Python integer and Fraction arithmetic',
                  provenClaim='No pure 27-qubit sixth-order anticoherent state has a proper rotational symmetry of order four or higher.',
                  scope='Symmetry subclass only. This argument alone does not exclude unrestricted sixth-order existence.',
                  fourthOrderCases=cyclic4, higherOrderCases=higher,
                  finiteCasesChecked=len(cyclic4)+len(higher), passed=True)
    target = Path(__file__).with_name('symmetry-verification.json')
    target.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['fourthOrderCases','higherOrderCases']}))

if __name__ == '__main__':
    main()
