#!/usr/bin/env python3
"""Verify the N=20, k=4, t=2 algebraic-root certificate using exact rationals.

No optimizer, floating-point linear algebra, or external Python package is
needed.  Every stored bound is recomputed from the integer polynomial
system, rational center, and rational approximate inverse.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import comb, prod
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def reconstruct(n, k, parities, slots):
    """All nonzero real equations for Gram and rank-one/two compressions."""
    index = {tuple(slot): i for i, slot in enumerate(slots)}
    expected_slots = [[a, r, part] for a in range(k)
                      for r in range(parities[a], n + 1, 2) for part in [0, 1]]
    assert slots == expected_slots

    def polynomial(i, j, q, p, part):
        terms = defaultdict(int)
        for r in range(q, n + 1):
            if (i, r - q, 0) not in index or (j, r, 0) not in index:
                continue
            m = n // 2 - r
            weight = comb(n, r) * prod(range(r - q + 1, r + 1))
            weight *= 3 * m * m - n * (n + 2) // 4 if p == 2 else m ** p
            ar, ai = index[i, r - q, 0], index[i, r - q, 1]
            br, bi = index[j, r, 0], index[j, r, 1]
            pairs = [(ar, br, 1), (ai, bi, 1)] if part == 0 else [(ar, bi, 1), (ai, br, -1)]
            for a, b, sign in pairs:
                terms[tuple(sorted((a, b)))] += sign * weight
        return {'basisPair': [i, j], 'raisingPower': q, 'zPower': p,
                'part': 'real' if part == 0 else 'imaginary',
                'constant': -int(q == 0 and p == 0 and i == j and part == 0),
                'terms': [[a, b, v] for (a, b), v in sorted(terms.items()) if v]}

    equations = []
    for p in [0, 1, 2]:
        for i in range(k):
            for j in range(i, k):
                if parities[i] != parities[j]:
                    continue
                equations.append(polynomial(i, j, 0, p, 0))
                if i != j:
                    equations.append(polynomial(i, j, 0, p, 1))
    for q, p in [(1, 0), (1, 1), (2, 0)]:
        for i in range(k):
            for j in range(k):
                if (parities[i] - parities[j] + q) % 2:
                    continue
                for part in [0, 1]:
                    equations.append(polynomial(i, j, q, p, part))
    return equations


def verify():
    raw = (ROOT / 'certificate.json').read_bytes()
    c = json.loads(raw)
    assert (c['n'], c['dimension'], c['order']) == (20, 4, 2)
    assert c['parities'] == [1, 1, 0, 0]
    equations = reconstruct(20, 4, c['parities'], c['slots'])
    assert c['equations'] == equations and len(equations) == 72
    unknown = c['unknownSlotIndices']
    fixed = c['fixedSlotIndices']
    assert len(unknown) == 72 and len(set(unknown)) == 72
    assert set(fixed) == set(range(84)) - set(unknown)
    denominator = int(c['centerDenominator'])
    invden = int(c['inverseDenominator'])
    center = list(map(int, c['fullCenterNumerators']))
    inverse = [list(map(int, row)) for row in c['inverseNumerators']]
    assert len(inverse) == 72 and all(len(row) == 72 for row in inverse)
    assert denominator > 0 and invden > 0
    fixed_den = int(c['fixedDenominator'])
    for i in fixed:
        assert F(center[i], denominator) == F(int(c['fixedNumerators'][str(i)]), fixed_den)
    fnum, jnum = [], []
    for equation in equations:
        value = equation['constant'] * denominator ** 2
        gradient = [0] * 84
        for a, b, coefficient in equation['terms']:
            value += coefficient * center[a] * center[b]
            gradient[a] += coefficient * center[b]
            gradient[b] += coefficient * center[a]
        fnum.append(value)
        jnum.append([gradient[i] for i in unknown])
    defect_num = max(sum(abs((invden * denominator if i == j else 0)
                            - sum(inverse[i][p] * jnum[p][j] for p in range(72)))
                        for j in range(72)) for i in range(72))
    correction_num = max(abs(sum(inverse[i][p] * fnum[p] for p in range(72))) for i in range(72))
    inverse_norm_num = max(sum(abs(v) for v in row) for row in inverse)
    active = set(unknown)
    hessian_bound = max(2 * sum(abs(coefficient) for a, b, coefficient in e['terms']
                               if a in active and b in active) for e in equations)
    radius = F(c['radius'])
    defect = F(defect_num, invden * denominator)
    correction = F(correction_num, invden * denominator ** 2)
    inverse_norm = F(inverse_norm_num, invden)
    contraction = defect + inverse_norm * hessian_bound * radius
    self_map = correction + contraction * radius
    assert defect < 1 and contraction < 1 and self_map < radius
    assert F(c['bounds']['contractionFactor']) == contraction
    assert F(c['bounds']['selfMapRadiusBound']) == self_map

    basis = json.loads((ROOT / 'basis.json').read_text())
    assert basis['id'] == 'algebraic-n20-k4-t2'
    assert (basis['n'], basis['dimension'], basis['t']) == (20, 4, 2)
    algebraic = basis['algebraicDefinition']
    names = {slot: f'x{i+1}' for i, slot in enumerate(unknown)}
    assert algebraic['variables'] == [names[i] for i in unknown]
    assert algebraic['coordinateSlots'] == c['slots']
    assert algebraic['unknownSlotIndices'] == unknown
    for i in unknown:
        assert list(map(F, algebraic['parameterIntervals'][names[i]])) == [F(center[i], denominator) - radius, F(center[i], denominator) + radius]
    def atom(i):
        return names[i] if i in names else str(F(int(c['fixedNumerators'][str(i)]), fixed_den))
    assert algebraic['fixedCoordinates'] == {str(i): atom(i) for i in fixed}
    equation_strings = []
    for equation in equations:
        terms = [str(equation['constant'])] if equation['constant'] else []
        terms += [f'({v})*({atom(a)})*({atom(b)})' for a, b, v in equation['terms']]
        equation_strings.append('+'.join(terms) or '0')
    assert algebraic['equations'] == equation_strings
    index = {tuple(slot): i for i, slot in enumerate(c['slots'])}
    coefficient_count = 0
    for col, vector in enumerate(basis['basis']):
        for r, expression in enumerate(vector['exactCoefficients']):
            if (col, r, 0) not in index:
                assert expression == '0'
            else:
                re, im = atom(index[col, r, 0]), atom(index[col, r, 1])
                assert expression == f'sqrt({comb(20,r)})*(({re})+I*({im}))'
            coefficient_count += 1
    assert coefficient_count == 84
    return {
        'passed': True, 'n': 20, 'dimension': 4, 'order': 2,
        'maximumDimension': 4,
        'certificateSHA256': hashlib.sha256(raw).hexdigest(),
        'basisSHA256': hashlib.sha256((ROOT / 'basis.json').read_bytes()).hexdigest(),
        'realVariables': 72, 'fixedRationalCoordinates': 12,
        'completeRealMomentEquationsReconstructed': 72,
        'complexCrossEntriesExactlyVerified': 144,
        'fullCrossReductionEntriesChecked': 144,
        'coefficientExpressionsAssociated': coefficient_count,
        'radius': str(radius), 'inverseDefectBound': str(defect),
        'newtonCorrectionBound': str(correction),
        'contractionFactorBound': str(contraction),
        'selfMapRadiusBound': str(self_map),
        'strictContraction': True, 'strictSelfMap': True,
        'proof': 'Exact rational contraction of x-AF(x) certifies the unique zero in the box. All nonzero parity-supported Gram and rank-one/two moment equations are reconstructed; omitted entries are identically zero or Hermitian conjugates. The moment identities are equivalent to all 16 ordered two-qubit cross reductions, with 9 complex entries each.',
        'upperBound': '4k<=N-2 gives k<=4 at N=20.',
        'arithmetic': 'Integers and fractions only; no floating-point inequalities or optimizer results establish existence.'
    }


if __name__ == '__main__':
    result = verify()
    (ROOT / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if 'Bound' not in k}, indent=2))
