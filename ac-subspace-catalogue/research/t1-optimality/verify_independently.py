#!/usr/bin/env python3
"""Independent exact checks of the first-order spin-one inertia bound.

This verifier imports no proof/construction implementation.  It constructs
the physical spin matrices directly in the normalized Dicke basis, forms
W=J.S, and computes the characteristic polynomial of each magnetization
block.  A separate symbolic computation proves the block factors for
arbitrary N.  Known optimal N=4 and N=7 bases check zero compression of W.
"""

from pathlib import Path
import hashlib
import json
import sympy as s

ROOT = Path(__file__).resolve().parent
CATALOGUE = ROOT.parents[1] / 'data' / 'catalogue.json'
x = s.symbols('x')


def assert_zero(value, context):
    if s.simplify(s.expand(value)) != 0:
        raise AssertionError(context + ': ' + str(value))


def spins(n):
    j = s.Rational(n, 2)
    z = s.diag(*[j - r for r in range(n + 1)])
    plus = s.zeros(n + 1)
    for r in range(1, n + 1):
        plus[r - 1, r] = s.sqrt(r * (n - r + 1))
    minus = plus.T
    return z, plus, minus


def physical_w(n):
    jz, jp, jm = spins(n)
    sz, sp, sm = spins(2)
    return s.kronecker_product(jz, sz) + (
        s.kronecker_product(jp, sm) + s.kronecker_product(jm, sp)) / 2


def symbolic_blocks():
    n, w = s.symbols('N w', integer=True)
    # The binomially scaled spin basis removes radicals.  This is a
    # similarity transformation of the physical matrix, not an assumption
    # about allowed state coefficients.
    block = s.zeros(3)
    for a in range(3):
        block[a, a] = (n / 2 - w + a) * (1 - a)
        if a < 2:
            block[a + 1, a] = s.Rational(2 - a, 2) * (w - a)
        if a > 0:
            block[a - 1, a] = s.Rational(a, 2) * (n - w + a)
    cubic = (x - n / 2) * (x + 1) * (x + n / 2 + 1)
    assert_zero(block.charpoly(x).as_expr() - cubic,
                'General interior magnetization block')
    edge = block.subs(w, 1).extract([0, 1], [0, 1])
    quadratic = (x - n / 2) * (x + 1)
    assert_zero(edge.charpoly(x).as_expr() - quadratic,
                'General boundary magnetization block')
    return {
        'passed': True,
        'interiorBlock': [[str(v) for v in row] for row in block.tolist()],
        'interiorCharacteristicPolynomial': str(s.factor(cubic)),
        'edgeCharacteristicPolynomial': str(s.factor(quadratic)),
        'endpointEigenvalue': 'N/2',
        'multiplicitiesForNAtLeast2': {
            'N/2': 'N+3', '-1': 'N+1', '-(N+2)/2': 'N-1'},
        'inertiaForNAtLeast2': {'positive': 'N+3', 'negative': '2*N', 'zero': 0},
        'dimensionConclusion': '3*k <= min(N+3, 2*N)',
    }


def verify_physical_spectrum(n):
    wmat = physical_w(n)
    if wmat != wmat.H:
        raise AssertionError('Physical W is not Hermitian')
    eigenvalues = [s.Rational(n, 2), -s.Integer(1), -s.Rational(n + 2, 2)]
    observed_mult = [0, 0, 0]
    block_reports = []
    for total in range(n + 3):
        indices = [3 * r + a for r in range(n + 1) for a in range(3) if r + a == total]
        block = wmat.extract(indices, indices)
        poly = block.charpoly(x).as_expr()
        expected = s.Integer(1)
        for i in range(len(indices)):
            expected *= x - eigenvalues[i]
            observed_mult[i] += 1
        assert_zero(poly - expected, f'N={n}, excitation block {total}')
        block_reports.append({'totalExcitation': total, 'size': len(indices),
                              'characteristicPolynomial': str(s.factor(poly))})
    expected_mult = [n + 3, n + 1, max(0, n - 1)]
    if observed_mult != expected_mult:
        raise AssertionError('Physical eigenvalue multiplicity mismatch')
    return {'n': n, 'passed': True, 'matrixDimension': wmat.rows,
            'eigenvalues': [str(v) for v in eigenvalues],
            'multiplicities': observed_mult, 'blocks': block_reports}


def verify_code(record):
    n, k = record['n'], record['dimension']
    basis = s.Matrix.hstack(*[
        s.Matrix([s.sympify(c) for c in row['exactCoefficients']])
        for row in record['basis']])
    gram = basis.H * basis - s.eye(k)
    for value in gram:
        assert_zero(value, record['id'] + ': orthonormality')
    lift = s.kronecker_product(basis, s.eye(3))
    compression = lift.H * physical_w(n) * lift
    for value in compression:
        assert_zero(value, record['id'] + ': full ancillary compression')
    return {'id': record['id'], 'n': n, 'dimension': k, 'passed': True,
            'orthonormalityEntries': k * k,
            'exactAncillaryCompressionEntries': 9 * k * k,
            'ancillaryCodeDimension': 3 * k,
            'positiveSpectralDimension': n + 3,
            'negativeSpectralDimension': 2 * n}


def main():
    data = json.loads(CATALOGUE.read_text())
    ids = ['paper-c7', 'cyclic-radical-n7-k3-t1']
    by_id = {record['id']: record for record in data['subspaces']}
    symbolic = symbolic_blocks()
    spectra = [verify_physical_spectrum(n) for n in range(1, 17)]
    controls = [verify_code(by_id[identifier]) for identifier in ids]
    report = {
        'passed': True,
        'method': ('Independent symbolic magnetization-block characteristic polynomials; '
                   'direct physical Hermitian spin matrices; exact full ancillary '
                   'compression for two known whole-subspace bases. No proof-module imports.'),
        'catalogueSha256': hashlib.sha256(CATALOGUE.read_bytes()).hexdigest(),
        'symbolicAllNVerification': symbolic,
        'physicalSpectra': spectra,
        'exactPositiveControls': controls,
        'scope': ('The symbolic characteristic identities establish the spectrum for all N. '
                  'Finite physical matrices and known codes independently check conventions. '
                  'A nonzero vector with zero projection into the smaller-sign spectral '
                  'space has a strictly signed quadratic form, so projection is injective '
                  'on every totally isotropic subspace.')
    }
    (ROOT / 'independent-verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'symbolicAllN': True,
                      'physicalSpectraChecked': len(spectra),
                      'exactPositiveControls': controls}, indent=2))


if __name__ == '__main__':
    main()
