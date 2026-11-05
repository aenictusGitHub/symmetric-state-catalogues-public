#!/usr/bin/env python3
"""Exact spin-1-ancilla inertia bound, using only the standard library.

In the nonorthonormal rescaled Dicke basis b_r=sqrt(C(N,r))D_r,
J_+ b_r=(r+1)b_(r+1), J_- b_r=(N-r+1)b_(r-1), J_z b_r=(r-N/2)b_r.
We build each conserved-total-m block of 2 X = 2 J.S exactly. A positive
diagonal change of metric makes these matrices Hermitian. The three exact
eigenvalues are -N-2,-2,N and the blockwise characteristic polynomials are
verified without numerical diagonalization.
"""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import json


def perm_sign(p):
    return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))


def multiply(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def charpoly(matrix):
    d = len(matrix)
    result = [F(0)]*(d+1)
    for perm in permutations(range(d)):
        p = [F(perm_sign(perm))]
        for i, j in enumerate(perm):
            p = multiply(p, [-matrix[i][j], F(i == j)])
        for a, v in enumerate(p):
            result[a] += v
    return result


def divide_root(p, root):
    q = [F(0)]*(len(p)-1)
    q[-1] = p[-1]
    for a in range(len(q)-2, -1, -1):
        q[a] = p[a+1]+root*q[a+1]
    assert p[0]+root*q[0] == 0
    return q


def value(p, x):
    ans = F(0)
    for c in p[::-1]:
        ans = ans*x+c
    return ans


def check(n):
    counts = {-n-2: 0, -2: 0, n: 0}
    for total in range(n+3):
        basis = [(r, total-r) for r in range(n+1) if 0 <= total-r <= 2]
        where = {v: i for i, v in enumerate(basis)}
        block = [[F(0) for _ in basis] for _ in basis]
        for col, (r, a) in enumerate(basis):
            block[col][col] = F((2*r-n)*(a-1))
            if r < n and a > 0:
                block[where[r+1, a-1]][col] = (r+1)*(3-a)
            if r > 0 and a < 2:
                block[where[r-1, a+1]][col] = (n-r+1)*(a+1)
        p = charpoly(block)
        for root in counts:
            while len(p) > 1 and value(p, root) == 0:
                p = divide_root(p, root)
                counts[root] += 1
        assert p == [1], (n, total, p)
    assert counts == {-n-2: n-1, -2: n+1, n: n+3}
    inertia_upper = min(2*n, n+3)//3
    upper = min(inertia_upper, n//2)
    return {'n': n, 'dimension': 3*(n+1), 'positive': n+3,
            'negative': 2*n, 'zero': 0, 'inertiaUpper': inertia_upper,
            'oneQubitSplittingUpper': n//2, 'upper': upper,
            'targetMatches': upper == n//3+1 if n >= 4 else upper == [0, 1, 1][n-1]}


def verify(max_n=300):
    rows = [check(n) for n in range(1, max_n+1)]
    assert all(row['targetMatches'] for row in rows)
    report = {'passed': True, 'method': 'Exact block characteristic polynomials in the rescaled Dicke basis, using Fraction; no symbolic or numerical eigensolver.',
              'universalSpectrum2X': [{'eigenvalue': '-N-2', 'multiplicity': 'N-1'},
                                     {'eigenvalue': '-2', 'multiplicity': 'N+1'},
                                     {'eigenvalue': 'N', 'multiplicity': 'N+3'}],
              'universalBound': '3k <= min(2N,N+3) for N>=2',
              'oneQubitSplittingBound': '2k <= N for every N>=1',
              'maximumDimension': '0 at N=1; 1 at N=2,3; floor(N/3)+1 at every N>=4',
              'lowerBoundSource': '../analytical-constructions/construct_t1.py: reflection_pairs, central_odd, verify_general_blocks; the latter verifies the central-block identities and positivity for their complete parameter ranges.',
              'verificationScope': 'This script verifies the upper-bound witness matrices. The universal upper-bound argument is in proof.tex, and general lower-bound validity is established by the cited construction, not inferred from finite matrix checks.',
              'matrixChecks': rows}
    return report


if __name__ == '__main__':
    report = verify()
    (Path(__file__).resolve().parent/'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'matrixChecks'}, indent=2))
    print(f'Passed {len(report["matrixChecks"])} exact matrix spectra, N=1..300.')
