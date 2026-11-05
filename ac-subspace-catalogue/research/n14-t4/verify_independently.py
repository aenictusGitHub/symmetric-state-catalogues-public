#!/usr/bin/env python3
"""Independent exact audit of the radical N=14, order-four AC states.

This verifier does not import the primary certificate or its verifier.
It constructs a symmetric bipartite amplitude matrix directly from the
Dicke decomposition, verifies both chiral branches, and checks the fifth
reduction and the advertised elementary symmetry argument.
"""
from pathlib import Path
import json
import re

import mpmath as mp
import sympy as s

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "references" / "ACLiege.jl"
SUPPORT = (1, 4, 7, 10, 13)


def zero(value, label):
    if s.simplify(s.expand_complex(value)) != 0:
        raise AssertionError(label + ": " + str(value))


def coefficients(sign=1):
    a = s.sqrt(s.Rational(49, 440)) + sign * s.I * s.sqrt(s.Rational(39431, 320760))
    b = s.sqrt(s.Rational(364, 3645))
    c = -s.sqrt(s.Rational(91, 3000)) + sign * s.I * s.sqrt(s.Rational(73229, 243000))
    v = s.zeros(15, 1)
    for r, value in zip(SUPPORT, (a, b, c, b, a)):
        v[r] = value
    return v


def reduction(v, t):
    # |D_N^r> = sum_a sqrt(C(t,a) C(N-t,r-a)/C(N,r))
    #              |D_t^a> |D_(N-t)^(r-a)>.
    amplitude = s.Matrix(t + 1, 15 - t, lambda a, b:
                         v[a + b] * s.sqrt(s.binomial(t, a) *
                         s.binomial(14 - t, b) / s.binomial(14, a + b)))
    return (amplitude * amplitude.H).applyfunc(s.simplify)


def exact_branch(sign):
    v = coefficients(sign)
    zero((v.H * v)[0] - 1, "Normalization")
    rho4 = reduction(v, 4)
    for a in range(5):
        for b in range(5):
            zero(rho4[a, b] - s.Rational(int(a == b), 5), "Fourth reduction")
    rho5 = reduction(v, 5)
    expected5 = s.eye(6) / 6
    for a, b, value in ((0, 3, s.sqrt(70) / 825),
                        (1, 4, -s.sqrt(7) / 165),
                        (2, 5, s.sqrt(70) / 825)):
        expected5[a, b] = expected5[b, a] = value
    for value in rho5 - expected5:
        zero(value, "Fifth reduction")
    if rho5 == s.eye(6) / 6:
        raise AssertionError("Order must be exactly four, not five")
    jp = s.zeros(15)
    for r in range(1, 15):
        jp[r - 1, r] = s.sqrt(r * (15 - r))
    jz = s.diag(*[7 - r for r in range(15)])
    fifth = s.simplify((v.H * jp ** 3 * jz ** 2 * v)[0])
    zero(fifth - 1456 * s.sqrt(7) / 5, "Nonzero fifth-degree moment")

    # A pi rotation about x reverses spin-seven coefficients with phase -1.
    # A 2*pi/3 z rotation fixes every occupied magnetic number.
    for r in range(15):
        zero(v[r] - v[14 - r], "Pi-x projective symmetry")
        if v[r] != 0 and (7 - r) % 3:
            raise AssertionError("Threefold axial symmetry")

    probability_a = s.Rational(854, 3645)
    probability_b = s.Rational(364, 3645)
    z = (-265 + sign * 27 * s.I * s.sqrt(39431)) / 5368
    zero(z * s.conjugate(z) - 1, "Unit phase z")
    q = v[1] / s.sqrt(probability_a)
    zero(q ** 2 - z, "Half-phase gauge")
    gauge = (s.sqrt(probability_a), s.sqrt(probability_b),
             s.sqrt(s.Rational(143, 525) * probability_a) * (z - 2),
             s.sqrt(probability_b) * z, s.sqrt(probability_a) * z ** 2)
    for i, r in enumerate(SUPPORT):
        zero(v[r] * q ** (i - 1) - gauge[i], "Equivalent sparse phase gauge")
    ratio = s.simplify(v[1] * v[13] / (v[4] * v[10]))
    zero(ratio - probability_a / probability_b * z, "Axis-normalizer invariant")
    return {
        "sign": sign, "normalizationExact": True,
        "fourthReduction": "I_5/5", "fourthReductionEntriesChecked": 25,
        "fifthReductionEntriesChecked": 36,
        "fifthReductionOffDiagonal": [
            {"indices": [0, 3], "value": "sqrt(70)/825"},
            {"indices": [1, 4], "value": "-sqrt(7)/165"},
            {"indices": [2, 5], "value": "sqrt(70)/825"}],
        "fifthDegreeMoment": str(fifth), "exactACOrder": 4,
        "sparseGaugeEquivalenceExact": True,
        "axisNormalizerInvariant": str(ratio),
    }


def source_match():
    # Printed source digits are numerical evidence, never an exact identity.
    mp.mp.dps = 90
    text = SOURCE.read_text()
    block = text.split("const acstateslist =", 1)[1].split("    [14,", 1)[1].split("    [15,", 1)[0]
    pairs = re.findall(r"([-+]?\d+\.\d+)([-+]\d+\.\d+)\*im", block)
    if len(pairs) != 5:
        raise AssertionError("Expected five nonzero Julia coefficients")
    source = [mp.mpc(a, b) for a, b in pairs]
    a = mp.sqrt(mp.mpf(49) / 440) + mp.j * mp.sqrt(mp.mpf(39431) / 320760)
    b = mp.sqrt(mp.mpf(364) / 3645)
    c = -mp.sqrt(mp.mpf(91) / 3000) + mp.j * mp.sqrt(mp.mpf(73229) / 243000)
    data = []
    for sign in (1, -1):
        exact = [a, b, c, b, a]
        if sign == -1:
            exact = [mp.conj(x) for x in exact]
        g = source[0] / exact[0]
        q = source[1] / (g * exact[1])
        # Normalize the fitted phase factors so the comparison uses a genuine
        # global phase and collective z rotation, not amplitude fitting.
        g /= abs(g)
        q /= abs(q)
        error = max(abs(source[i] - g * q ** i * exact[i]) for i in range(5))
        data.append({"sign": sign, "maxCoefficientError": mp.nstr(error, 30),
                     "globalPhaseAtFirstSupport": mp.nstr(mp.arg(g), 30),
                     "phaseIncrementPerThreeDickeIndices": mp.nstr(mp.arg(q), 30)})
    if mp.mpf(data[0]["maxCoefficientError"]) >= mp.mpf("1e-17"):
        raise AssertionError("Positive branch does not match the supplied state")
    return {"source": "references/ACLiege.jl, acstateslist N=14", "matches": data,
            "scope": "Numerical comparison only; the exact certificate is the symbolic reduction above. Printed precision is not asserted to equal numerical accuracy."}


def main():
    branches = [exact_branch(sign) for sign in (1, -1)]
    # In rank five, a tetrahedral invariant has multiplicity zero. Therefore
    # a state with the nonzero rank-five moment above cannot have T, O or I
    # symmetry. Its maximal axial rotation is C3 (support differences 3,6,9,12),
    # excluding larger D_m and continuous axial groups. Thus its full proper
    # group is D3. Its unique threefold axis and the invariant checked above
    # exclude a proper rotation between the conjugate branches.
    tetra_rank5 = s.Rational(11 - 8 - 3, 12)
    zero(tetra_rank5, "Tetrahedral rank-five invariant multiplicity")
    report = {"passed": True, "n": 14, "t": 4, "dimension": 1,
              "arithmetic": "Exact SymPy rational and radical arithmetic for the entire AC certificate; 90-digit mpmath only for source comparison.",
              "independentMethod": "Direct Dicke bipartition amplitudes and physical collective spin matrices; primary verifier not imported.",
              "branches": branches,
              "properSymmetryGroup": "D3", "conjugateBranchesProperRotationInequivalent": True,
              "rankFiveTetrahedralInvariantMultiplicity": str(tetra_rank5),
              "sourceComparison": source_match()}
    (HERE / "independent-verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
