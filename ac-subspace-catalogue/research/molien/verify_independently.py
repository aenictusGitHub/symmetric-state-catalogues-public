#!/usr/bin/env python3
"""Independently verify the explicit Molien AC-subspace bases.

This program does not import the generator.  It checks compressed spin
moments rather than using the generator's reduced-density-matrix test.

Write c_r=(-1)^r p_r/sqrt(binomial(N,r)*h), where p_r lies in
Q(sqrt(-3)) and h=sum_r |p_r|^2/binomial(N,r).
After removing the nonzero normalization factors, the matrix element of
J_+^q J_z^p between two rows is the rational number

  (-1)^q sum_r conjugate(p_i[r+q]) p_j[r]
                         (r+q)!/(r! binomial(N,r)) (r-N/2)^p.

For q=0 it must be delta_ij*h_i times the uniform spin-p moment; for
q>0 it must vanish.  These operators and their adjoints span all spin
tensors of ranks through t when p+q<=t, so the checks include every
off-diagonal compressed matrix element and certify the whole subspace.
"""

from fractions import Fraction as F
from pathlib import Path
from math import comb, factorial
import hashlib
import json
import sympy as s
from dataclasses import dataclass

ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Quadratic:
    """A separate exact implementation of Q(sqrt(-3))."""
    a: F = F(0)
    b: F = F(0)

    @staticmethod
    def of(x):
        return x if isinstance(x, Quadratic) else Quadratic(F(x))

    def __bool__(self):
        return bool(self.a or self.b)

    def __add__(self, other):
        other = self.of(other)
        return Quadratic(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Quadratic(-self.a, -self.b)

    def __sub__(self, other):
        return self + -self.of(other)

    def __mul__(self, other):
        other = self.of(other)
        return Quadratic(self.a * other.a - 3 * self.b * other.b,
                         self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def conjugate(self):
        return Quadratic(self.a, -self.b)

    def __truediv__(self, other):
        other = self.of(other)
        den = other.a * other.a + 3 * other.b * other.b
        require(bool(den), 'Division by zero in quadratic field')
        num = self * other.conjugate()
        return Quadratic(num.a / den, num.b / den)

    def symbolic(self):
        return (s.Rational(self.a.numerator, self.a.denominator)
                + s.I * s.sqrt(3) * s.Rational(self.b.numerator, self.b.denominator))

    def __str__(self):
        return str(self.a) if not self.b else f'({self.a}) + ({self.b})*I*sqrt(3)'


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def matrix_rank(rows):
    a = [list(row) for row in rows]
    pivot = 0
    for col in range(len(a[0])):
        i = next((i for i in range(pivot, len(a)) if a[i][col]), None)
        if i is None:
            continue
        a[pivot], a[i] = a[i], a[pivot]
        lead = a[pivot][col]
        a[pivot] = [v / lead for v in a[pivot]]
        for j in range(pivot + 1, len(a)):
            factor = a[j][col]
            if factor:
                a[j] = [v - factor * w for v, w in zip(a[j], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def character_dimension(n, group):
    """Average SU(2) characters over the corresponding rotation group.

    For odd N the binary central element -I acts as -1.  For even N,
    the rotations in the quotient give the class sums below.  This is
    independent of the generator's monomial/Hilbert-series enumeration.
    """
    if n % 2:
        return 0
    m = n // 2
    parity = (-1) ** m
    chi3 = (1, 0, -1)[m % 3]
    chi4 = (1, 1, -1, -1)[m % 4]
    aliases = {
        'q8': 'Q8', 'quaternion': 'Q8', 'quaternion q8': 'Q8',
        'quaternion-character': 'Q8-character',
        '2t': '2T', 'tetrahedral': '2T', 'binary tetrahedral': '2T',
        '2o': '2O', 'octahedral': '2O', 'binary octahedral': '2O',
        '2i': '2I', 'icosahedral': '2I', 'binary icosahedral': '2I',
        'tetrahedral-character': '2T-character',
        'octahedral-sign': '2O-sign',
    }
    kind = aliases.get(str(group).lower(), group)
    if kind == 'Q8':
        numerator, denominator = n + 1 + 3 * parity, 4
    elif kind == 'Q8-character':
        numerator, denominator = n + 1 - parity, 4
    elif kind == '2T':
        numerator, denominator = n + 1 + 3 * parity + 8 * chi3, 12
    elif kind == '2T-character':
        numerator, denominator = n + 1 + 3 * parity - 4 * chi3, 12
    elif kind == '2O':
        numerator, denominator = n + 1 + 9 * parity + 8 * chi3 + 6 * chi4, 24
    elif kind == '2O-sign':
        numerator, denominator = n + 1 - 3 * parity + 8 * chi3 - 6 * chi4, 24
    elif kind == '2I':
        numerator = n + 1 + 15 * parity + 20 * chi3 + 12 * (2 - m % 5)
        denominator = 60
    else:
        raise ValueError('Unrecognized group for independent character check: ' + str(group))
    require(numerator % denominator == 0, 'Character average must be integral')
    return numerator // denominator


def moment_residuals(rows, norms, n, t):
    weights = [F(1, comb(n, r)) for r in range(n + 1)]
    moments = [sum((F(2 * r - n, 2) ** p for r in range(n + 1)), F(0)) / (n + 1)
               for p in range(t + 1)]
    for q in range(t + 1):
        for p in range(t - q + 1):
            factors = [((-1) ** q) * F(factorial(r + q), factorial(r))
                       * weights[r] * F(2 * r - n, 2) ** p
                       for r in range(n - q + 1)]
            for i, row_i in enumerate(rows):
                for j, row_j in enumerate(rows):
                    observed = sum((row_i[r + q].conjugate() * row_j[r] * factors[r]
                                    for r in range(n - q + 1)), Quadratic())
                    target = norms[i] * moments[p] if q == 0 and i == j else F(0)
                    yield (i, j, q, p), observed - target


def verify(record):
    n = record['n']
    t = record.get('verifiedAchievedOrder', record.get('t', record.get('order')))
    rows = [[Quadratic(F(c[0]), F(c[1])) for c in vector['polynomialCoefficientsQuadratic']]
            if 'polynomialCoefficientsQuadratic' in vector else
            [Quadratic(F(c)) for c in vector['polynomialCoefficients']]
            for vector in record['basis']]
    require(all(len(row) == n + 1 for row in rows), record['id'] + ': row length')
    k = record['dimension']
    require(len(rows) == k, record['id'] + ': basis size')
    rank = matrix_rank(rows)
    require(rank == k, record['id'] + ': exact linear independence')
    norms_q = [sum((c.conjugate() * c / comb(n, r) for r, c in enumerate(row)), Quadratic()) for row in rows]
    require(all(not h.b for h in norms_q), record['id'] + ': Hermitian norms are real')
    norms = [h.a for h in norms_q]
    require(all(h > 0 for h in norms), record['id'] + ': nonzero norms')
    coefficient_checks = 0
    for i, (row, h, vector) in enumerate(zip(rows, norms, record['basis'])):
        require(h == F(vector['normSquared']), record['id'] + ': claimed polynomial norm')
        for r, p in enumerate(row):
            require(s.simplify(s.sympify(vector['polynomialCoefficients'][r]) - p.symbolic()) == 0,
                    f"{record['id']}: polynomial field encoding, row {i}, index {r}")
            exact = s.sympify(vector['exactCoefficients'][r])
            expected = (-1) ** r * p.symbolic() / s.sqrt(
                comb(n, r) * s.Rational(h.numerator, h.denominator))
            require(s.simplify(exact - expected) == 0,
                    f"{record['id']}: Dicke conversion, row {i}, index {r}")
            coefficient_checks += 1
    checks = 0
    for key, residual in moment_residuals(rows, norms, n, t):
        require(not residual, f"{record['id']}: compressed moment {key} = {residual}")
        checks += 1
    group = record['group']
    predicted = character_dimension(n, group)
    require(k == predicted, record['id'] + ': character dimension')
    next_failure = None
    if t < n:
        for key, residual in moment_residuals(rows, norms, n, t + 1):
            if residual:
                next_failure = {'basisPair': list(key[:2]), 'raisingPower': key[2],
                                'zPower': key[3], 'unnormalizedResidual': str(residual)}
                break
        require(next_failure is not None, record['id'] + ': claimed highest AC order')
    return {
        'id': record['id'], 'n': n, 'group': group, 'dimension': k,
        'achievedOrder': t, 'passed': True,
        'exactPolynomialRank': rank,
        'characterAverageDimension': predicted,
        'orthonormal': True,
        'exactDickeCoefficientChecks': coefficient_checks,
        'exactCompressedMomentEntries': checks,
        'nextOrderTested': t + 1 if t < n else None,
        'nextOrderFailure': next_failure,
    }


def main():
    source = ROOT / 'bases.json'
    data = json.loads(source.read_text())
    catalogue = json.loads((ROOT.parents[1] / 'data' / 'catalogue.json').read_text())
    require(data['qubitNumbers'] == catalogue['numbers'], 'All catalogue sizes must be covered')
    expected = {(family['id'], n): character_dimension(n, family['id'])
                for family in data['families'] for n in data['qubitNumbers']}
    claimed = {(item['group'], item['n']): item['dimension'] for item in data['dimensionTable']}
    require(expected == claimed, 'Complete dimension table disagrees with character averages')
    supplied = {(record['group'], record['n']): record['dimension'] for record in data['subspaces']}
    require(supplied == {key: k for key, k in expected.items() if k},
            'Every nonzero family component must supply a basis')
    reports = []
    for record in data['subspaces']:
        report = verify(record)
        reports.append(report)
        print(record['id'], 'PASS', flush=True)
    out = {
        'passed': True,
        'inputSha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'method': ('Independent exact rational compressed spin moments J_+^q J_z^p, '
                   'including every ordered pair of basis vectors; character-average '
                   'dimensions, rational matrix rank, norms, and exact Dicke conversion. '
                   'No import from the construction generator.'),
        'subspacesChecked': len(reports),
        'groupAndSizeDimensionChecks': len(expected),
        'allCatalogueSizesCovered': True,
        'allNonzeroFamilyComponentsHaveVerifiedBases': True,
        'exactCompressedMomentEntries': sum(r['exactCompressedMomentEntries'] for r in reports),
        'exactDickeCoefficientChecks': sum(r['exactDickeCoefficientChecks'] for r in reports),
        'witnesses': reports,
    }
    (ROOT / 'independent-verification.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'witnesses'}, indent=2))


if __name__ == '__main__':
    main()
