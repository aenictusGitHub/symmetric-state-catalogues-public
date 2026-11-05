#!/usr/bin/env python3
"""Exact radical witnesses for N=14, t=4, with a source-state comparison.

No fitted decimal is used to establish anticoherence.  The comparison with
ACLiege.jl is a separate high-precision numerical check of the rotation gauge.
"""
from pathlib import Path
import hashlib
import json
import re

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parents[1] / 'references' / 'ACLiege.jl'


def coefficients(sign=1):
    a = 7 / s.sqrt(440) + sign * s.I * s.sqrt(s.Rational(39431, 320760))
    b = s.sqrt(s.Rational(364, 3645))
    c = -s.sqrt(s.Rational(91, 3000)) + sign * s.I * s.sqrt(s.Rational(73229, 243000))
    values = s.zeros(15, 1)
    for r, value in {1: a, 4: b, 7: c, 10: b, 13: a}.items():
        values[r] = value
    return values


def reduction(vector, t):
    amplitudes = s.zeros(t + 1, 15 - t)
    for a in range(t + 1):
        for b in range(15 - t):
            amplitudes[a, b] = vector[a+b] * s.sqrt(
                s.binomial(t, a) * s.binomial(14-t, b) / s.binomial(14, a+b))
    return (amplitudes * amplitudes.H).applyfunc(s.simplify)


def source_comparison():
    mp.mp.dps = 90
    source = SOURCE.read_text()
    block = re.search(r'\[14,\s*\[\s*0,\s*\n(.*?)\n\s*\],\s*\n\s*\],', source, re.S)
    if block is None:
        raise AssertionError('Cannot locate the N=14 source vector')
    triples = re.findall(r'([-+]?\d+\.\d+)\s*([+-])\s*(\d+\.\d+)\s*\*im', block.group(1))
    original = [mp.mpc(real, ('-' if sign == '-' else '') + imag)
                for real, sign, imag in triples]
    if len(original) != 5:
        raise AssertionError('Expected five source amplitudes')
    a = 7/mp.sqrt(440) + 1j*mp.sqrt(mp.mpf(39431)/320760)
    b = mp.sqrt(mp.mpf(364)/3645)
    c = -mp.sqrt(mp.mpf(91)/3000) + 1j*mp.sqrt(mp.mpf(73229)/243000)
    vector, support = [a, b, c, b, a], [1, 4, 7, 10, 13]
    beta = (mp.arg(original[1]/original[0]) + mp.arg(a))/3
    gamma = mp.arg(original[0]/a) - beta
    errors = [abs(value*mp.exp(1j*(gamma+r*beta))-saved)
              for r, value, saved in zip(support, vector, original)]
    if max(errors) >= mp.mpf('1e-17'):
        raise AssertionError('Source state does not match the radical witness')
    return {
        'passed': True,
        'source': 'references/ACLiege.jl',
        'sourceSha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'precisionDecimalDigits': mp.mp.dps,
        'convention': 'source approximately exp(i*chi) exp(-i*beta*J_z) psi_plus; J_z|D_14^r>=(7-r)|D_14^r>',
        'betaRadians': mp.nstr(beta, 75),
        'chiRadians': mp.nstr(gamma+7*beta, 75),
        'maximumDickeCoefficientError': mp.nstr(max(errors), 25),
        'sourceNormalizationError': mp.nstr(abs(sum(abs(v)**2 for v in original)-1), 25),
        'scope': 'Numerical comparison only. Printed source precision does not imply the same number of accurate digits; exact existence is proved by the symbolic reductions.'
    }


def branch(sign):
    vector = coefficients(sign)
    if s.simplify((vector.H*vector)[0]) != 1:
        raise AssertionError('Normalization')
    rho4 = reduction(vector, 4)
    if rho4 != s.eye(5)/5:
        raise AssertionError('Fourth-order anticoherence')
    rho5 = reduction(vector, 5)
    purity5 = s.simplify(s.trace(rho5*rho5))
    if purity5 != s.Rational(15209, 90750):
        raise AssertionError('Fifth-order purity')
    if s.simplify(purity5-s.Rational(1, 6)) != s.Rational(14, 15125):
        raise AssertionError('Strictly nonzero fifth-order defect')
    ratio = s.simplify(vector[1]*vector[13]/(vector[4]*vector[10]))
    if s.simplify(s.im(ratio)) == 0:
        raise AssertionError('The normalizer invariant must distinguish conjugation')
    if any(s.simplify(vector[r]-vector[14-r]) != 0 for r in range(15)):
        raise AssertionError('Palindromic D3 symmetry')
    # The rank-5 representation contains no tetrahedral invariant:
    # chi_5(0)=11, chi_5(2*pi/3)=-1, chi_5(pi)=-1.
    tetrahedral_multiplicity = s.Rational(11+8*(-1)+3*(-1), 12)
    if tetrahedral_multiplicity != 0:
        raise AssertionError('Tetrahedral rank-five character')
    return {
        'sign': sign, 'passed': True, 'normalizationExact': True,
        'rho4Exact': 'I_5/5', 'exactReductionEntriesChecked': 25,
        'rho5Purity': str(purity5), 'rho5PurityExcessAboveOneSixth': '14/15125',
        'exactACOrder': 4, 'normalizerInvariant': str(ratio),
        'properRotationGroup': 'D3', 'tetrahedralRankFiveInvariantMultiplicity': 0,
    }


def main():
    saved = json.loads((ROOT/'basis.json').read_text())
    expected = coefficients(1)
    actual = s.Matrix([s.sympify(v) for v in saved['basis'][0]['exactCoefficients']])
    if any(s.simplify(a-b) != 0 for a, b in zip(actual, expected)):
        raise AssertionError('Saved basis does not match the exact definition')
    results = [branch(1), branch(-1)]
    report = {
        'passed': True, 'n': 14, 't': 4, 'dimension': 1,
        'arithmetic': 'Exact symbolic rational arithmetic and radicals',
        'basisSha256': hashlib.sha256((ROOT/'basis.json').read_bytes()).hexdigest(),
        'branches': results, 'sourceComparison': source_comparison(),
        'optimalDimension': {'lower': 1, 'upper': 1, 'proof': 'Nine orthogonal rank-zero/one/two tensor images per code vector must fit in dimension 15: 9k<=15.'},
        'scope': 'Two exact one-dimensional AC subspaces. Their combined span is not asserted to be AC. No classification of all N=14 fourth-order states is asserted.'
    }
    (ROOT/'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
