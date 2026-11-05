#!/usr/bin/env python3
"""Exact physical spin-2 quadratic witness verification; standard library only."""
from math import comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent


def multiply(a, b):
    d = len(a)
    return [[sum(a[i][r]*b[r][j] for r in range(d)) for j in range(d)] for i in range(d)]


def charpoly(a):
    """Faddeev--LeVerrier, with exact integer divisions."""
    d = len(a)
    b = [[int(i == j) for j in range(d)] for i in range(d)]
    coeffs = [1]
    for r in range(1, d+1):
        b = multiply(a, b)
        trace = sum(b[i][i] for i in range(d))
        assert trace % r == 0
        c = -trace//r
        coeffs.append(c)
        for i in range(d):
            b[i][i] += c
    assert not any(any(row) for row in b)
    return coeffs[::-1]


def value(p, x):
    v = 0
    for c in p[::-1]:
        v = v*x+c
    return v


def divide_root(p, root):
    q = [0]*(len(p)-1)
    q[-1] = p[-1]
    for i in range(len(q)-2, -1, -1):
        q[i] = p[i+1]+root*q[i+1]
    assert p[0]+root*q[0] == 0
    return q


def check(n):
    assert n >= 4
    branches = [(-4*(n+2), n-3), (-2*n*(n+2), n-1),
                (-2*n*(n+5), n+1), (-10*n-8, n+3), (4*n*(n+2), n+5)]
    expected = {}
    for eigenvalue, multiplicity in branches:
        expected[eigenvalue] = expected.get(eigenvalue, 0)+multiplicity
    counts = dict.fromkeys(expected, 0)
    for total in range(n+5):
        basis = [(r, total-r) for r in range(n+1) if 0 <= total-r <= 4]
        where = {v: i for i, v in enumerate(basis)}
        d = len(basis)
        a = [[0]*d for _ in basis]  # A=2 J.S
        for col, (r, b) in enumerate(basis):
            a[col][col] = (2*r-n)*(b-2)
            if r < n and b > 0:
                a[where[r+1, b-1]][col] = (r+1)*(5-b)
            if r > 0 and b < 4:
                a[where[r-1, b+1]][col] = (n-r+1)*(b+1)
        a2 = multiply(a, a)
        # Four times the proposed W is integer in this rescaled Dicke basis.
        w = [[a2[i][j]+(n+6)*a[i][j]-2*n*(n+2)*int(i == j)
              for j in range(d)] for i in range(d)]
        metric = [comb(n, r)*comb(4, b) for r, b in basis]
        assert all(metric[i]*w[i][j] == metric[j]*w[j][i]
                   for i in range(d) for j in range(d))
        p = charpoly(w)
        for eigenvalue in counts:
            while len(p) > 1 and value(p, eigenvalue) == 0:
                p = divide_root(p, eigenvalue)
                counts[eigenvalue] += 1
        assert p == [1], (n, total, p)
    assert counts == expected, (n, counts, expected)
    positive = sum(m for e, m in counts.items() if e > 0)
    negative = sum(m for e, m in counts.items() if e < 0)
    assert positive == n+5 and negative == 4*n and 0 not in counts
    return {'n': n, 'positive': positive, 'negative': negative, 'zero': 0,
            'upper': (n+5)//5, 'branches4W': [{'eigenvalue': e, 'multiplicity': m} for e, m in branches]}


def verify(max_n=300):
    rows = [check(n) for n in range(4, max_n+1)]
    return {'passed': True,
            'method': 'Physical spin matrices in the rescaled Dicke basis; exact positive-metric self-adjointness and integer characteristic polynomials of every conserved-total-m block. No numerical or symbolic eigensolver.',
            'witness': 'W=X^2+(j+3)X-2j(j+1)I, X=J.S, auxiliary spin 2, j=N/2',
            'universalBound': '5k <= N+5 for every N>=4 and every second-order AC subspace',
            'spectralProof': 'proof.tex',
            'finiteChecksAreNotTheUniversalProof': True,
            'matrixChecks': rows}


if __name__ == '__main__':
    report = verify()
    (ROOT/'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'matrixChecks'}, indent=2))
    print(f'Passed {len(report["matrixChecks"])} exact physical spectra.')
