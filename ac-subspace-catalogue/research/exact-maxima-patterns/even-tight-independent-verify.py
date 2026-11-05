"""Independent N=16, t=6 channel audit using discrete polynomials and CG.

This does not import the main even-tight verifier or its 3j/6j tensor formula.
The rank-seven diagonal tensors are reconstructed by exact orthogonality to
lower-degree polynomials on the spin magnetic-number grids.
"""
import json
from pathlib import Path
import sympy as s
from sympy.physics.wigner import clebsch_gordan


def verify():
    x = s.symbols("x")

    def rank_seven_polynomial(j):
        coefficients = s.symbols("a:3")
        polynomial = x**7 + sum(a * x**power for a, power in
                                 zip(coefficients, (5, 3, 1)))
        equations = [sum((polynomial * x**power).subs(x, m)
                         for m in range(-j, j + 1)) for power in (1, 3, 5)]
        polynomial = s.expand(polynomial.subs(s.solve(equations, coefficients)))
        assert all(sum((polynomial * x**power).subs(x, m)
                       for m in range(-j, j + 1)) == 0 for power in range(7))
        return polynomial

    input_polynomial = rank_seven_polynomial(8)
    output_polynomial = rank_seven_polynomial(7)
    operators = []
    for ell in range(1, 4):
        for q in range(-ell, ell + 1):
            operator = s.zeros(15, 17)
            for a in range(-7, 8):
                m = a - q
                if -8 <= m <= 8:
                    operator[a + 7, m + 8] = s.sqrt(s.Rational(2 * ell + 1, 15)) \
                        * clebsch_gordan(8, ell, 7, m, q, a)
            assert s.simplify(s.trace(operator.H * operator)) == 1
            operators.append(operator)
    assert len(operators) == 15
    for i, operator in enumerate(operators):
        for previous in operators[:i]:
            assert s.simplify(s.trace(operator.H * previous)) == 0
    input_tensor = s.diag(*[input_polynomial.subs(x, m) for m in range(-8, 9)])
    output_tensor = s.diag(*[output_polynomial.subs(x, a) for a in range(-7, 8)])
    channel_output = sum((e * input_tensor * e.H for e in operators), s.zeros(15))
    monic_multiplier = s.simplify(channel_output[-1, -1] / output_tensor[-1, -1])
    assert monic_multiplier == s.Rational(115, 8398)
    assert (channel_output - monic_multiplier * output_tensor).applyfunc(s.simplify) == s.zeros(15)
    normalized_multiplier = s.simplify(monic_multiplier * s.sqrt(
        s.trace(output_tensor**2) / s.trace(input_tensor**2)))
    assert normalized_multiplier == 5 * s.sqrt(69) / 8398
    scalar_image = sum((e * e.H for e in operators), s.zeros(15))
    assert scalar_image.applyfunc(s.simplify) == s.eye(15)
    return {
        "passed": True,
        "case": {"n": 16, "t": 6, "maximumDimension": 0},
        "independentMethod": "Clebsch-Gordan operators and monic discrete orthogonal polynomials",
        "inputRankSevenPolynomial": str(input_polynomial),
        "outputRankSevenPolynomial": str(output_polynomial),
        "monicMultiplier": str(monic_multiplier),
        "normalizedMultiplier": str(normalized_multiplier),
        "checks": ["all lower polynomial moments vanish", "all 15 rectangular operators are HS orthonormal",
                   "full channel matrix proportionality", "channel image of identity", "nonzero rank-seven multiplier"],
        "proofImplication": "AC6 would force AC7; 20 orthogonal rectangular-tensor images cannot fit in dimension 16",
    }


if __name__ == "__main__":
    report = verify()
    Path(__file__).with_name("even-tight-independent-verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
