"""Exact bridges among existing optimal AC catalogue representatives.

Run from any directory. Writes bridge-verification.json beside this script.
No floating-point tolerances are used in the assertions.
"""
import json
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
CATALOGUE = json.loads((HERE.parent.parent / "data/catalogue.json").read_text())
SUBSPACES = {item["id"]: item for item in CATALOGUE["subspaces"]}


def zero(matrix):
    assert matrix.applyfunc(s.simplify) == s.zeros(*matrix.shape)


def split(c, n, t):
    return s.Matrix(t + 1, n - t + 1, lambda a, b:
                    c[a + b] * s.sqrt(s.binomial(t, a) * s.binomial(n - t, b)
                                      / s.binomial(n, a + b)))


def source(code_id, n, t):
    basis = [s.Matrix([s.sympify(x) for x in vector["exactCoefficients"]])
             for vector in SUBSPACES[code_id]["basis"]]
    for a, x in enumerate(basis):
        for b, y in enumerate(basis):
            target = s.eye(t + 1) / (t + 1) if a == b else s.zeros(t + 1)
            zero(split(x, n, t) * split(y, n, t).H - target)
            assert s.simplify((x.H * y)[0] - int(a == b)) == 0
    projector = sum((x * x.H for x in basis), s.zeros(n + 1))
    zero(projector * projector - projector)
    return projector


def check_top_rank(y, n):
    j = s.Rational(n, 2)
    jz = s.diag(*[j - r for r in range(n + 1)])
    jp = s.zeros(n + 1)
    for r in range(n):
        jp[r, r + 1] = s.sqrt((r + 1) * (n - r))
    jm = jp.H
    generators = [(jp + jm) / 2, (jp - jm) / (2 * s.I), jz]
    casimir = sum((a * a * y - 2 * a * y * a + y * a * a
                   for a in generators), s.zeros(n + 1))
    zero(casimir - n * (n + 1) * y)


def vectorize_centered(projector, n):
    k, d = s.trace(projector), n + 1
    y = projector - k * s.eye(d) / d
    check_top_rank(y, n)
    norm2 = k * (d - k) / d
    assert s.simplify(s.trace(y * y) - norm2) == 0
    zero(y * y - (1 - 2 * k / d) * y - k * (d - k) * s.eye(d) / d**2)
    f = s.zeros(d)
    for r in range(d):
        f[r, n - r] = (-1)**r
    a = y * f / s.sqrt(norm2)
    coefficients = [s.simplify(sum(
        a[r, ell - r] * s.sqrt(s.binomial(n, r) * s.binomial(n, ell - r)
                               / s.binomial(2 * n, ell))
        for r in range(max(0, ell - n), min(n, ell) + 1)))
        for ell in range(2 * n + 1)]
    zero(a - split(coefficients, 2 * n, n))
    zero(a * a.H - s.eye(d) / d - (d - 2 * k) * y / (k * (d - k)))
    reduction = split(coefficients, 2 * n, n - 1)
    zero(reduction * reduction.H - s.eye(n) / n)
    return coefficients


p4 = source("acliege-Π_4_2_1", 4, 1)
p6 = source("acliege-Π_6_3_1", 6, 1)
p7 = source("paper-e1", 7, 2)
r6 = s.Matrix(7, 7, lambda r, u:
               (s.sqrt((7 - r) * (7 - u)) * p7[r, u]
                + s.sqrt((r + 1) * (u + 1)) * p7[r + 1, u + 1]) / 7)
q = 2 * r6
zero(q * q - q)
assert s.trace(q) == 4
zero(p6 + q - s.eye(7))
zero(r6 * p6)

c8 = vectorize_centered(p4, 4)
c12 = [-x for x in vectorize_centered(p6, 6)]
c12_from_reduction = vectorize_centered(q, 6)
zero(s.Matrix(c12) - s.Matrix(c12_from_reduction))
expected8 = s.zeros(9, 1)
expected8[0] = expected8[8] = s.sqrt(30) / 12
expected8[4] = s.sqrt(21) / 6
expected12 = s.zeros(13, 1)
expected12[1] = -s.sqrt(7) / 5
expected12[6] = s.sqrt(11) / 5
expected12[11] = s.sqrt(7) / 5
zero(s.Matrix(c8) - expected8)
zero(s.Matrix(c12) - expected12)

report = {
    "passed": True,
    "arithmetic": "exact symbolic radicals; no numerical tolerance",
    "sourceSubspaces": ["acliege-Π_4_2_1", "acliege-Π_6_3_1", "paper-e1"],
    "complementIdentity": "P6 = I7 - 2 Tr_1(P7)",
    "normalizedReducedCode": "Tr_1(P7/2) = (I7-P6)/4",
    "kernelDimension": 3,
    "outputStates": [
        {"n": 8, "t": 3, "nonzeroDickeCoefficients":
         {i: str(x) for i, x in enumerate(c8) if x}},
        {"n": 12, "t": 5, "nonzeroDickeCoefficients":
         {i: str(x) for i, x in enumerate(c12) if x}},
    ],
    "checks": [
        "all source orthonormalities and ordered AC cross reductions",
        "all source and reduced projector idempotencies",
        "exact complement and kernel identities",
        "rank-four and rank-six tensor Casimir equations",
        "centered projector quadratic identities",
        "highest-spin vectorization entry by entry",
        "both explicit output coefficient formulas",
        "all entries of rho_3=I4/4 and rho_5=I6/6",
    ],
}
(HERE / "bridge-verification.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
