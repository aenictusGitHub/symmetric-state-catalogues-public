#!/usr/bin/env python3
"""Independent physical-matrix and exact-contraction audit of the N=20 code.

No constructor or numerical root-finder is imported. All certificate checks
use Python integers/Fraction; numerical basis entries are checked only as
display approximations, never as evidence of existence.
"""
import ast
from collections import defaultdict
from fractions import Fraction as F
from math import comb, sqrt
from pathlib import Path
import json
import time


def multiply(a, b):
    return [[sum(a[i][p]*b[p][j] for p in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def complex_add(a, b):
    return (a[0]+b[0], a[1]+b[1])


def complex_mul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def poly_expression(node, parameters):
    """Restricted exact parser, returning complex polynomials in 72 variables."""
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return {(): (F(node.value), F(0))} if node.value else {}
    if isinstance(node, ast.Name):
        if node.id == 'I':
            return {(): (F(0), F(1))}
        assert node.id in parameters
        return {(parameters[node.id],): (F(1), F(0))}
    if isinstance(node, ast.UnaryOp):
        p = poly_expression(node.operand, parameters)
        if isinstance(node.op, ast.USub):
            return {m: (-v[0], -v[1]) for m, v in p.items()}
        assert isinstance(node.op, ast.UAdd)
        return p
    assert isinstance(node, ast.BinOp)
    a, b = poly_expression(node.left, parameters), poly_expression(node.right, parameters)
    out = defaultdict(lambda: (F(0), F(0)))
    if isinstance(node.op, (ast.Add, ast.Sub)):
        sign = -1 if isinstance(node.op, ast.Sub) else 1
        for m, v in a.items():
            out[m] = complex_add(out[m], v)
        for m, v in b.items():
            out[m] = complex_add(out[m], (sign*v[0], sign*v[1]))
    elif isinstance(node.op, ast.Mult):
        for u, x in a.items():
            for v, y in b.items():
                m = tuple(sorted(u+v))
                out[m] = complex_add(out[m], complex_mul(x, y))
    else:
        assert isinstance(node.op, ast.Div) and set(b) == {()} and b[()][1] == 0
        denominator = b[()][0]
        assert denominator
        out.update({m: (v[0]/denominator, v[1]/denominator) for m, v in a.items()})
    return {m: v for m, v in out.items() if any(v)}


def verify(directory=None):
    started = time.monotonic()
    root = Path(directory) if directory else Path(__file__).resolve().parent
    c = json.loads((root/'certificate.json').read_text())
    b = json.loads((root/'basis.json').read_text())
    assert (c['n'], c['dimension'], c['order']) == (20, 4, 2)
    n, k, d = 20, 4, 21
    parities = c['parities']
    assert parities == [1, 1, 0, 0]
    expected_slots = [[a, r, part] for a in range(k) for r in range(parities[a], d, 2) for part in range(2)]
    assert c['slots'] == expected_slots
    slots = {tuple(v): i for i, v in enumerate(expected_slots)}
    unknown, fixed = c['unknownSlotIndices'], c['fixedSlotIndices']
    assert len(unknown) == len(set(unknown)) == 72
    assert len(fixed) == len(set(fixed)) == 12
    assert set(unknown).isdisjoint(fixed) and set(unknown+fixed) == set(range(84))

    # Physical spin matrices on b_r=sqrt(C(20,r))*D_r, with positive metric G.
    identity = [[int(r == s) for s in range(d)] for r in range(d)]
    z = [[(10-r)*int(r == s) for s in range(d)] for r in range(d)]
    plus = [[(21-s)*int(r == s-1) for s in range(d)] for r in range(d)]
    z2 = multiply(z, z)
    quad = [[3*z2[r][s]-110*identity[r][s] for s in range(d)] for r in range(d)]
    operators = {(0, 0): identity, (0, 1): z, (0, 2): quad,
                 (1, 0): plus, (1, 1): multiply(plus, z), (2, 0): multiply(plus, plus)}
    metric = [comb(n, r) for r in range(d)]

    def physical_terms(i, j, q, p, part):
        out = defaultdict(int)
        op = operators[q, p]
        for r in range(parities[i], d, 2):
            for s in range(parities[j], d, 2):
                weight = metric[r]*op[r][s]
                if not weight:
                    continue
                # conj(a+ib)*(c+id)=(ac+bd)+i(ad-bc).
                combos = [(0, 0, 1), (1, 1, 1)] if part == 'real' else [(0, 1, 1), (1, 0, -1)]
                for left, right, sign in combos:
                    pair = tuple(sorted((slots[i, r, left], slots[j, s, right])))
                    out[pair] += sign*weight
        return {pair: value for pair, value in out.items() if value}

    expected_keys = set()
    for (q, p) in operators:
        for i in range(k):
            for j in range(k):
                if (parities[i]-parities[j]+q) % 2:
                    assert not physical_terms(i, j, q, p, 'real')
                    assert not physical_terms(i, j, q, p, 'imaginary')
                    continue
                if q == 0:
                    if i > j:
                        continue
                    expected_keys.add((i, j, q, p, 'real'))
                    if i != j:
                        expected_keys.add((i, j, q, p, 'imaginary'))
                    else:
                        assert not physical_terms(i, j, q, p, 'imaginary')
                else:
                    expected_keys.update((i, j, q, p, part) for part in ('real', 'imaginary'))
    equations = c['equations']
    assert len(equations) == len(expected_keys) == 72
    supplied_keys = set()
    reconstructed = []
    for e in equations:
        i, j = e['basisPair']
        q, p, part = e['raisingPower'], e['zPower'], e['part']
        key = i, j, q, p, part
        assert key in expected_keys and key not in supplied_keys
        supplied_keys.add(key)
        terms = physical_terms(i, j, q, p, part)
        supplied = {(a, b): coeff for a, b, coeff in e['terms']}
        assert len(supplied) == len(e['terms']) and terms == supplied
        const = -int((q, p) == (0, 0) and i == j and part == 'real')
        assert e['constant'] == const
        reconstructed.append((const, terms))
    assert supplied_keys == expected_keys

    D, Q = int(c['centerDenominator']), int(c['inverseDenominator'])
    fixed_D = int(c['fixedDenominator'])
    center = list(map(int, c['fullCenterNumerators']))
    inverse = [list(map(int, row)) for row in c['inverseNumerators']]
    assert len(center) == 84 and len(inverse) == 72 and all(len(row) == 72 for row in inverse)
    fixed_values = {int(i): F(int(v), fixed_D) for i, v in c['fixedNumerators'].items()}
    assert set(fixed_values) == set(fixed)
    assert all(F(center[i], D) == fixed_values[i] for i in fixed)
    radius = F(c['radius'])
    assert radius > 0
    # Exact residual and Jacobian from the independently reconstructed equations.
    fnum, jacobian = [], []
    for const, terms in reconstructed:
        value_num = const*D*D
        gradient = [0]*84
        for (i, j), coefficient in terms.items():
            value_num += coefficient*center[i]*center[j]
            gradient[i] += coefficient*center[j]
            gradient[j] += coefficient*center[i]
        fnum.append(value_num)
        jacobian.append([gradient[i] for i in unknown])
    defect_num = max(sum(abs((Q*D if i == j else 0)-sum(inverse[i][p]*jacobian[p][j] for p in range(72)))
                         for j in range(72)) for i in range(72))
    beta_num = max(abs(sum(inverse[i][p]*fnum[p] for p in range(72))) for i in range(72))
    inverse_norm_num = max(sum(abs(v) for v in row) for row in inverse)
    active = set(unknown)
    hessian_lipschitz = max(sum(2*abs(v) for (a, b), v in terms.items() if a in active and b in active)
                           for _, terms in reconstructed)
    defect = F(defect_num, Q*D)
    beta = F(beta_num, Q*D*D)
    gamma = defect+F(inverse_norm_num, Q)*hessian_lipschitz*radius
    self_map = beta+gamma*radius
    assert defect < 1 and gamma < 1 and self_map < radius
    # Compare independently computed exact numbers with every supplied bound.
    bounds = c['bounds']
    assert defect == F(int(bounds['inverseDefectNumerator']), int(bounds['inverseDefectDenominator']))
    assert beta == F(int(bounds['newtonCorrectionNumerator']), int(bounds['newtonCorrectionDenominator']))
    assert F(inverse_norm_num, Q) == F(int(bounds['inverseNormNumerator']), int(bounds['inverseNormDenominator']))
    assert hessian_lipschitz == int(bounds['jacobianLipschitzBound'])
    assert gamma == F(bounds['contractionFactor']) and self_map == F(bounds['selfMapRadiusBound'])

    # Check exact basis association and all 72 rational isolation intervals.
    assert (b['n'], b['t'], b['dimension']) == (20, 2, 4)
    assert len(b['basis']) == 4
    parameter_names = {f'x{i+1}': slot for i, slot in enumerate(unknown)}
    intervals = b['algebraicDefinition']['parameterIntervals']
    assert set(intervals) == set(parameter_names)
    for name, slot in parameter_names.items():
        lo, hi = map(F, intervals[name])
        assert (lo, hi) == (F(center[slot], D)-radius, F(center[slot], D)+radius)
    checked_coefficients = 0
    max_display_error = 0.0
    for i, row in enumerate(b['basis']):
        assert len(row['exactCoefficients']) == len(row['coefficients']) == 21
        for r, expression in enumerate(row['exactCoefficients']):
            tree = ast.parse(expression, mode='eval').body
            if (i, r, 0) not in slots:
                assert poly_expression(tree, parameter_names) == {}
                assert row['coefficients'][r] == [0.0, 0.0]
                continue
            assert isinstance(tree, ast.BinOp) and isinstance(tree.op, ast.Mult)
            assert isinstance(tree.left, ast.Call) and isinstance(tree.left.func, ast.Name)
            assert tree.left.func.id == 'sqrt' and len(tree.left.args) == 1
            assert isinstance(tree.left.args[0], ast.Constant) and tree.left.args[0].value == comb(n, r)
            actual = poly_expression(tree.right, parameter_names)
            expected = defaultdict(lambda: (F(0), F(0)))
            for part in range(2):
                slot = slots[i, r, part]
                monomial = (slot,) if slot in active else ()
                coefficient = F(1) if slot in active else fixed_values[slot]
                v = (coefficient, F(0)) if part == 0 else (F(0), coefficient)
                expected[monomial] = complex_add(expected[monomial], v)
            expected = {m: v for m, v in expected.items() if any(v)}
            assert actual == expected
            for part in range(2):
                center_value = sqrt(comb(n, r))*float(F(center[slots[i, r, part]], D))
                error = abs(row['coefficients'][r][part]-center_value)
                max_display_error = max(max_display_error, error)
                assert error < 1e-14
            checked_coefficients += 1
    return {'passed': True, 'physicalEquationsReconstructed': 72,
            'variables': 72, 'fixedCoordinates': 12,
            'fullMomentCoverage': ['I', 'Jz', '3Jz^2-110I', 'J+', 'J+Jz', 'J+^2'],
            'omittedEntries': 'Parity zeros, diagonal imaginary zeros, and Hermitian conjugates only.',
            'exactBasisCoefficientsMatched': checked_coefficients,
            'parameterIntervalsMatched': 72,
            'inverseDefect': str(defect), 'newtonCorrection': str(beta),
            'jacobianLipschitzBound': hessian_lipschitz,
            'contractionFactor': str(gamma), 'selfMapRadiusBound': str(self_map),
            'radius': str(radius), 'strictContraction': True, 'strictSelfMap': True,
            'inverseNonsingular': True,
            'conclusion': 'There is a unique real algebraic root in the rational cube. Its four normalized Dicke vectors are orthonormal and span a second-order AC subspace. The exact spin-3/2 upper bound gives maximum dimension four at N=20.',
            'maxDisplayCoefficientError': max_display_error,
            'seconds': time.monotonic()-started}


if __name__ == '__main__':
    report = verify()
    (Path(__file__).resolve().parent/'independent-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('inverseDefect','newtonCorrection','contractionFactor','selfMapRadiusBound')}, indent=2))
