#!/usr/bin/env python3
"""Reproducible exact AC subspaces from Klein's binary polyhedral forms.

All polynomial algebra, orthogonalization and full cross reductions use the
quadratic field Q(sqrt(-3)), implemented with Python Fraction. SymPy is used
only to format the final normalized Dicke coefficients and their LaTeX.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from math import comb, gcd, lcm, sqrt
from pathlib import Path
import argparse
import json

ROOT = Path(__file__).resolve().parent
DEFAULT_N = list(range(1, 21)) + [22, 24, 35, 42]


@dataclass(frozen=True)
class K:
    """a+b*sqrt(-3), with the positive imaginary square root."""
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    def __add__(self, other):
        q = as_k(other)
        return K(self.a + q.a, self.b + q.b)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b)

    def __sub__(self, other):
        return self + (-as_k(other))

    def __rsub__(self, other):
        return as_k(other) - self

    def __mul__(self, other):
        q = as_k(other)
        return K(self.a*q.a - 3*self.b*q.b, self.a*q.b + self.b*q.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        q = as_k(other)
        d = q.a*q.a + 3*q.b*q.b
        if not d:
            raise ZeroDivisionError
        z = self*q.conjugate()
        return K(z.a/d, z.b/d)

    def conjugate(self):
        return K(self.a, -self.b)

    def __bool__(self):
        return bool(self.a or self.b)

    def pair(self):
        return [str(self.a), str(self.b)]

    def expression(self):
        if not self.b:
            return str(self.a)
        if not self.a:
            return f'({self.b})*sqrt(3)*I'
        return f'({self.a})+({self.b})*sqrt(3)*I'


def as_k(x):
    return x if isinstance(x, K) else K(x)


def polynomial(degree, terms):
    """Coefficient of x^r y^(degree-r), including homogeneous zero endpoints."""
    return [as_k(terms.get(r, 0)) for r in range(degree+1)]


def multiply(p, q):
    out = [K()] * (len(p)+len(q)-1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                if b:
                    out[i+j] = out[i+j] + a*b
    return out


def power(p, exponent):
    out = [K(1)]
    for _ in range(exponent):
        out = multiply(out, p)
    return out


def add(p, q, factor=1):
    assert len(p) == len(q)
    return [a + factor*b for a, b in zip(p, q)]


def inner(p, q, n):
    return sum((a.conjugate()*b / comb(n, r) for r, (a, b) in enumerate(zip(p, q))), K())


def primitive(p):
    denominator = lcm(*(v.denominator for a in p for v in (a.a, a.b)))
    ints = [(int(a.a*denominator), int(a.b*denominator)) for a in p]
    divisor = gcd(*(v for pair in ints for v in pair))
    if not divisor:
        raise ValueError('Dependent polynomial in proposed free-module basis')
    first = next(v for pair in ints for v in pair if v)
    divisor *= 1 if first > 0 else -1
    return [K(a//divisor, b//divisor) for a, b in ints]


def orthogonalize(polynomials, n):
    basis = []
    norms = []
    for source in polynomials:
        p = source[:]
        for q, norm in zip(basis, norms):
            projection = inner(q, p, n) / norm
            p = [a-projection*b for a, b in zip(p, q)]
        p = primitive(p)
        norm = inner(p, p, n)
        assert not norm.b and norm.a > 0
        basis.append(p)
        norms.append(norm.a)
    for a, p in enumerate(basis):
        for b, q in enumerate(basis):
            assert inner(p, q, n) == K(norms[a] if a == b else 0)
    return basis, norms


def definitions():
    a4 = polynomial(4, {2: 1})
    b4 = polynomial(4, {0: 1, 4: 1})
    f6 = polynomial(6, {1: -1, 5: 1})
    g8 = polynomial(8, {0: 1, 4: 14, 8: 1})
    h12 = polynomial(12, {0: 1, 4: -33, 8: -33, 12: 1})
    f4 = polynomial(4, {0: 1, 2: K(0, 2), 4: 1})
    f12 = polynomial(12, {1: -1, 6: 11, 11: 1})
    g20 = polynomial(20, {0: 1, 5: 228, 10: 494, 15: -228, 20: 1})
    h30 = polynomial(30, {0: 1, 5: -522, 10: -10005, 20: -10005, 25: 522, 30: 1})
    one = [K(1)]
    families = [
        ('quaternion', 'Quaternion fixed subspace', 1, (a4, b4), (one, f6), 'A028242', 'trivial'),
        ('quaternion-character', 'Quaternion character subspace', 1, (a4, b4),
         (polynomial(2, {1: 1}), polynomial(4, {0: -1, 4: 1})), None, 'nontrivial'),
        ('tetrahedral', 'Binary tetrahedral fixed subspace', 2, (f6, g8), (one, h12), 'A008612', 'trivial'),
        ('tetrahedral-character', 'Binary tetrahedral character subspace', 2, (f6, g8),
         (f4, power([c.conjugate() for c in f4], 2)), None, 'one nontrivial character; conjugate gives the other'),
        ('octahedral', 'Binary octahedral fixed subspace', 3, (g8, power(f6, 2)),
         (one, multiply(f6, h12)), 'A008647', 'trivial'),
        ('octahedral-sign', 'Binary octahedral sign subspace', 3, (g8, power(f6, 2)),
         (f6, h12), None, 'sign'),
        ('icosahedral', 'Binary icosahedral fixed subspace', 5, (f12, g20), (one, h30), 'A008651', 'trivial'),
    ]
    # Relations ensure the invariant Hilbert numerators use only two secondaries.
    assert power(f6, 2) == multiply(a4, add(power(b4, 2), power(a4, 2), -4))
    assert power(h12, 2) == add(power(g8, 3), power(f6, 4), -108)
    assert power(h30, 2) == add(power(g20, 3), power(f12, 5), 1728)
    assert multiply(f4, [c.conjugate() for c in f4]) == g8
    return families


def generators_at_n(family, n):
    _, _, _, (a, b), secondaries, _, _ = family
    da, db = len(a)-1, len(b)-1
    polys, labels = [], []
    for c, secondary in enumerate(secondaries):
        remaining = n-(len(secondary)-1)
        for u in range(remaining//da+1):
            rest = remaining-u*da
            if rest >= 0 and rest % db == 0:
                v = rest//db
                polys.append(multiply(multiply(power(a, u), power(b, v)), secondary))
                labels.append({'primaryExponents': [u, v], 'secondaryIndex': c})
    return polys, labels


def named_generators(group):
    names = {
        'quaternion': (['Q_4', 'R_4'], ['1', 'F_6']),
        'quaternion-character': (['Q_4', 'R_4'], ['U_2', 'V_4']),
        'tetrahedral': (['F_6', 'G_8'], ['1', 'H_{12}']),
        'tetrahedral-character': (['F_6', 'G_8'], ['A_4', r'\overline{A}_4^2']),
        'octahedral': (['G_8', 'F_6^2'], ['1', 'F_6 H_{12}']),
        'octahedral-sign': (['G_8', 'F_6^2'], ['F_6', 'H_{12}']),
        'icosahedral': (['I_{12}', 'I_{20}'], ['1', 'I_{30}']),
    }
    return names[group]


def source_label(group, label):
    primaries, secondaries = named_generators(group)
    factors = []
    for name, exponent in zip(primaries, label['primaryExponents']):
        if exponent == 1:
            factors.append(name)
        elif exponent > 1:
            factors.append(r'\left(' + name + r'\right)^{' + str(exponent) + '}')
    secondary = secondaries[label['secondaryIndex']]
    if secondary != '1':
        factors.append(secondary)
    return r'\,'.join(factors) if factors else '1'


def polynomial_latex(p):
    import sympy as s
    x, y = s.symbols('x y')
    degree = len(p)-1
    expr = sum((s.Rational(a.a.numerator, a.a.denominator)+s.I*s.sqrt(3)*s.Rational(a.b.numerator, a.b.denominator))*x**r*y**(degree-r)
               for r, a in enumerate(p))
    return s.latex(expr)


def dimension_from_character(group, n):
    """Independent character averages, not monomial counts."""
    if n % 2:
        return 0
    j = n//2
    a, b, c = (-1)**j, [1, 0, -1][j % 3], [1, 1, -1, -1][j % 4]
    if group == 'quaternion':
        return (n+1+3*a)//4
    if group == 'quaternion-character':
        return (n+1-a)//4
    if group == 'tetrahedral':
        return (n+1+3*a+8*b)//12
    if group == 'tetrahedral-character':
        return (n+1+3*a-4*b)//12
    if group == 'octahedral':
        return (n+1+9*a+8*b+6*c)//24
    if group == 'octahedral-sign':
        return (n+1-3*a+8*b-6*c)//24
    if group == 'icosahedral':
        return j//2+j//3+j//5-j+1
    raise ValueError(group)


def check_reductions(basis, norms, n, t):
    """Every ordered cross-reduction entry, with nonzero radicals removed.

    rho_ab = (-1)^(a+b) sqrt(C(t,a)C(t,b)/(norm_i norm_j))
       * sum_l p_i[a+l] conjugate(p_j[b+l]) C(N-t,l)
                   / (C(N,a+l) C(N,b+l)).
    Thus all checks can be made in Q(sqrt(-3)), exactly.
    """
    checked = 0
    first_failure = None
    binomial_n = [comb(n, r) for r in range(n+1)]
    for i, p in enumerate(basis):
        for j, q in enumerate(basis):
            for a in range(t+1):
                for b in range(t+1):
                    value = sum((p[a+r]*q[b+r].conjugate()*F(comb(n-t, r),
                                 binomial_n[a+r]*binomial_n[b+r])
                                 for r in range(n-t+1)), K())
                    target = K(norms[i] / ((t+1)*comb(t, a))) if i == j and a == b else K()
                    checked += 1
                    if value != target and first_failure is None:
                        first_failure = {'basisA': i, 'basisB': j, 'reducedRow': a,
                                         'reducedColumn': b, 'strippedResidual': (value-target).expression()}
    return {'passed': first_failure is None, 'order': t, 'orthonormal': True,
            'fullRank': True, 'matrixElementsChecked': checked,
            'violations': [] if first_failure is None else [first_failure]}


def serialize_basis(p, norm, n, index):
    import sympy as s
    coefficients, exact, support, latex_terms = [], [], [], []
    for r, a in enumerate(p):
        if not a:
            coefficients.append([0.0, 0.0])
            exact.append('0')
            continue
        sign = (-1)**r
        q = s.Rational(a.a.numerator, a.a.denominator) + s.I*s.sqrt(3)*s.Rational(a.b.numerator, a.b.denominator)
        d = norm*comb(n, r)
        expr = sign*q/s.sqrt(s.Rational(d.numerator, d.denominator))
        exact.append(str(expr))
        coefficients.append([sign*float(a.a)/sqrt(float(d)), sign*sqrt(3)*float(a.b)/sqrt(float(d))])
        support.append(r)
        latex_terms.append('(' + s.latex(expr) + r')|D_{' + str(n) + '}^{' + str(r) + r'}\rangle')
    return {'exactCoefficients': exact, 'coefficients': coefficients, 'support': support,
            'stateLatex': '|u_{' + str(index+1) + r'}\rangle=' + '+'.join(latex_terms),
            'polynomialCoefficients': [a.expression() for a in p],
            'polynomialCoefficientsQuadratic': [a.pair() for a in p],
            'normSquared': str(norm)}


def build(ns=DEFAULT_N, max_order=7):
    records, dimensions, family_metadata = [], [], []
    for family in definitions():
        group, name, guarantee, primaries, secondaries, oeis, character = family
        family_metadata.append({'id': group, 'name': name, 'guaranteedOrder': guarantee,
                                'character': character, 'oeis': oeis,
                                'primaryDegrees': [len(p)-1 for p in primaries],
                                'secondaryDegrees': [len(p)-1 for p in secondaries],
                                'primaryNamesLatex': named_generators(group)[0],
                                'secondaryNamesLatex': named_generators(group)[1],
                                'primaryPolynomialsLatex': [polynomial_latex(p) for p in primaries],
                                'secondaryPolynomialsLatex': [polynomial_latex(p) for p in secondaries],
                                'primaryPolynomials': [[a.expression() for a in p] for p in primaries],
                                'secondaryPolynomials': [[a.expression() for a in p] for p in secondaries]})
        for n in ns:
            source, labels = generators_at_n(family, n)
            expected = dimension_from_character(group, n)
            assert len(source) == expected, (group, n, len(source), expected)
            dimensions.append({'group': group, 'n': n, 'dimension': expected})
            if not source:
                continue
            basis, norms = orthogonalize(source, n)
            verification = None
            history = []
            next_order = None
            for t in range(1, min(max_order, n)+1):
                result = check_reductions(basis, norms, n, t)
                history.append(result)
                if not result['passed']:
                    next_order = result
                    break
                verification = result
            assert verification is not None and verification['order'] >= guarantee, (group, n, history)
            achieved = verification['order']
            rid = f'molien-{group}-n{n}-d{len(basis)}-t{achieved}'
            records.append({'id': rid, 'name': name, 'n': n, 't': achieved, 'order': achieved,
                            'dimension': len(basis), 'group': group, 'character': character,
                            'guaranteedOrder': guarantee, 'verifiedAchievedOrder': achieved,
                            'orderIsExact': next_order is not None,
                            'basis': [serialize_basis(p, norm, n, i) for i, (p, norm) in enumerate(zip(basis, norms))],
                            'sourceMonomials': [dict(label, monomialLatex=source_label(group, label)) for label in labels],
                            'sourcePolynomialBasis': [source_label(group, label) for label in labels],
                            'verification': verification,
                            'verificationByOrder': history, 'nextOrder': next_order,
                            'provenance': 'Klein invariant or single-character polynomial construction; exact full cross reductions and weighted Gram-Schmidt. This dimension is a constructive lower bound, not a global optimality claim.'})
            print(f'{group}: N={n}, dimension={len(basis)}, AC order={achieved}', flush=True)
    data = {'schemaVersion': 1, 'convention': 'p(x,y)=sum_r p_r x^r y^(N-r); c_r=(-1)^r p_r/sqrt(binomial(N,r)*normSquared). Quadratic pairs [a,b] mean a+b*i*sqrt(3).',
            'qubitNumbers': list(ns), 'maximumOrderChecked': max_order,
            'sources': ['https://arxiv.org/html/1610.05841v2', 'https://arxiv.org/pdf/math/0503542',
                        'https://doi.org/10.3390/sym12010177'],
            'families': family_metadata, 'dimensionTable': dimensions, 'subspaces': records}
    report = {'passed': True, 'method': 'Exact arithmetic in Q(sqrt(-3)); primitive weighted Gram-Schmidt; independent group-character dimensions; every ordered Dicke cross reduction through the highest verified order and the first failing order (up to 7).',
              'familyCount': len(family_metadata), 'subspaceCount': len(records),
              'dimensionCountsChecked': len(dimensions), 'generatorRelationsPassed': True,
              'exactCrossEntriesChecked': sum(v['matrixElementsChecked'] for r in records for v in r['verificationByOrder']),
              'witnesses': [{k: r[k] for k in ('id', 'n', 'group', 'dimension', 'guaranteedOrder', 'verifiedAchievedOrder', 'orderIsExact')} for r in records]}
    return data, report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, nargs='+', default=DEFAULT_N)
    parser.add_argument('--max-order', type=int, default=7)
    args = parser.parse_args()
    data, report = build(args.n, args.max_order)
    (ROOT/'bases.json').write_text(json.dumps(data, indent=2)+'\n')
    (ROOT/'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'witnesses'}, indent=2))
