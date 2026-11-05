#!/usr/bin/env python3
"""Independent exact verification of 5k <= N+5 for second-order AC codes.

Constructs physical spin matrices without importing the main verifier,
checks magnetization-block characteristic polynomials, and independently
checks the all-N symbolic block spectrum and full auxiliary compression
on known exact AC subspaces.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

ROOT = Path(__file__).resolve().parent
CATALOGUE = ROOT.parents[1] / 'data' / 'catalogue.json'
x = s.symbols('x')


def zero(v, message):
    if s.simplify(s.expand(v)) != 0:
        raise AssertionError(message + ': ' + str(v))


def spin(n):
    z = s.diag(*[s.Rational(n, 2) - r for r in range(n + 1)])
    plus = s.zeros(n + 1)
    for r in range(1, n + 1):
        plus[r - 1, r] = s.sqrt(r * (n - r + 1))
    return z, plus, plus.T


def physical(n):
    jz, jp, jm = spin(n)
    sz, sp, sm = spin(4)
    coupling = s.kronecker_product(jz, sz) + (
        s.kronecker_product(jp, sm) + s.kronecker_product(jm, sp)) / 2
    witness = coupling ** 2 + s.Rational(n + 6, 2) * coupling - s.Rational(n * (n + 2), 2) * s.eye(coupling.rows)
    return coupling, witness


def symbolic():
    n, w = s.symbols('N w', integer=True)
    block = s.zeros(5)
    for a in range(5):
        block[a, a] = (n / 2 - w + a) * (2 - a)
        if a < 4:
            block[a + 1, a] = s.Rational(4 - a, 2) * (w - a)
        if a > 0:
            block[a - 1, a] = s.Rational(a, 2) * (n - w + a)
    values = [n, n / 2 - 2, -s.Integer(3), -n / 2 - 3, -n - 2]
    expected = s.prod(x - v for v in values)
    zero(block.charpoly(x).as_expr() - expected, 'All-N interior block spectrum')
    for size in range(1, 5):
        edge = block.subs(w, size - 1).extract(range(size), range(size))
        zero(edge.charpoly(x).as_expr() - s.prod(x - v for v in values[:size]),
             f'All-N boundary block of size {size}')
    polynomial = lambda v: s.factor(v * v + (n + 6) * v / 2 - n * (n + 2) / 2)
    actual = list(map(polynomial, values))
    predicted = [n * (n + 2), -(5 * n + 4) / 2, -n * (n + 5) / 2,
                 -n * (n + 2) / 2, -n - 2]
    for a, b in zip(actual, predicted):
        zero(a - b, 'All-N quadratic witness eigenvalue')
    # Positivity/negativity on N>=4 can also be checked coefficientwise
    # after N=q+4, q>=0; all signed polynomial coefficients are positive.
    q = s.symbols('q', nonnegative=True)
    for i, v in enumerate(actual):
        signed = v if i == 0 else -v
        require_positive = s.Poly(s.expand(signed.subs(n, q + 4)), q).all_coeffs()
        if not all(c > 0 for c in require_positive):
            raise AssertionError('Strict sign for entire N>=4 domain')
    return {
        'passed': True,
        'method': 'Exact all-N tridiagonal magnetization-block characteristic polynomials, including every boundary size; exact polynomial substitution and coefficientwise strict signs.',
        'XEigenvaluesDescendingTotalSpin': [str(v) for v in values],
        'WEigenvaluesDescendingTotalSpin': [str(v) for v in actual],
        'multiplicitiesDescendingTotalSpin': ['N+5', 'N+3', 'N+1', 'N-1', 'N-3'],
        'inertiaForNAtLeast4': {'positive': 'N+5', 'negative': '4*N', 'zero': 0},
        'conclusion': '5*k <= N+5 for N>=4 and t>=2',
    }


def physical_spectrum(n):
    coupling, witness = physical(n)
    if witness != witness.H or coupling != coupling.H:
        raise AssertionError('Physical operators must be Hermitian')
    values = [s.Integer(n), s.Rational(n, 2) - 2, -s.Integer(3), -s.Rational(n, 2) - 3, -s.Integer(n + 2)]
    multiplicities = [0] * 5
    for total in range(n + 5):
        indices = [5 * r + a for r in range(n + 1) for a in range(5) if r + a == total]
        block = coupling.extract(indices, indices)
        zero(block.charpoly(x).as_expr() - s.prod(x - v for v in values[:len(indices)]),
             f'N={n}: physical coupling block {total}')
        wblock = witness.extract(indices, indices)
        zero((wblock - (block ** 2 + s.Rational(n + 6, 2) * block
                        - s.Rational(n * (n + 2), 2) * s.eye(len(indices)))).norm(),
             f'N={n}: physical witness block {total}')
        for i in range(len(indices)):
            multiplicities[i] += 1
    if multiplicities != [n + 5, n + 3, n + 1, n - 1, n - 3]:
        raise AssertionError('Physical multiplicity count')
    evals = [s.factor(v ** 2 + s.Rational(n + 6, 2) * v - s.Rational(n * (n + 2), 2)) for v in values]
    if not (evals[0] > 0 and all(v < 0 for v in evals[1:])):
        raise AssertionError('Physical witness signs')
    return {'n': n, 'passed': True, 'dimension': witness.rows,
            'eigenvalues': [str(v) for v in evals], 'multiplicities': multiplicities}


def control(record):
    n, k = record['n'], record['dimension']
    basis = s.Matrix.hstack(*[s.Matrix([s.sympify(c) for c in b['exactCoefficients']])
                             for b in record['basis']])
    for v in basis.H * basis - s.eye(k):
        zero(v, 'Control orthonormality')
    coupling, witness = physical(n)
    lift = s.kronecker_product(basis, s.eye(5))
    for v in lift.H * witness * lift:
        zero(v, 'Control full auxiliary witness compression')
    for v in lift.H * coupling * lift:
        zero(v, 'Control first-moment compression')
    scalar = s.Rational(n * (n + 2), 2)
    for v in lift.H * coupling ** 2 * lift - scalar * s.eye(5 * k):
        zero(v, 'Control second-moment compression')
    return {'id': record['id'], 'n': n, 'dimension': k, 'passed': True,
            'exactFullWitnessCompressionEntries': 25 * k * k,
            'exactFirstAndSecondMomentEntries': 50 * k * k}


def main():
    data = json.loads(CATALOGUE.read_text())
    selected = [next(r for r in data['subspaces'] if r['n'] == n and r['t'] >= 2
                     and r['dimension'] == k and r['evidence'] == 'exact')
                for n, k in [(4, 1), (7, 2)]]
    report = {'passed': True,
              'catalogueSha256': hashlib.sha256(CATALOGUE.read_bytes()).hexdigest(),
              'symbolicAllNVerification': symbolic(),
              'physicalSpectra': [physical_spectrum(n) for n in range(4, 13)],
              'exactPositiveControls': [control(record) for record in selected],
              'scope': 'No implementation or output from the primary verifier is imported.'}
    (ROOT / 'independent-verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'physicalSpectra'}, indent=2))


if __name__ == '__main__':
    main()
