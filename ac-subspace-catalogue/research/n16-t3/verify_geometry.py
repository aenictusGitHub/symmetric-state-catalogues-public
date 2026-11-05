#!/usr/bin/env python3
"""Exact Majorana multiplicities and checked numerical coordinates for N=16.

The sparse polynomials live in QQ(sqrt(385), I*sqrt(10)); declaring that
field explicitly avoids primitive-element discovery on unsimplified radicals.
"""
from pathlib import Path
import cmath
import hashlib
import importlib.util
import json
import math

import sympy as s

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "majorana_helpers", HERE.parents[1] / "scripts" / "majorana_coordinates.py")
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
y = s.Symbol("y")


def zero(value, label):
    if s.simplify(s.expand_complex(value)) != 0:
        raise AssertionError(label + ": " + str(value))


def representative(original, index):
    exact_strings = original["exactCoefficients"]
    coefficients = list(map(s.sympify, exact_strings))
    saved = [complex(a, b) for a, b in original["coefficients"]]
    numeric = [complex(s.N(c, 50)) for c in coefficients]
    formula_error = helpers.ray_error(saved, numeric)
    assert formula_error < 1e-14
    a, b = s.sqrt(385), s.I * s.sqrt(10)
    if index == 0:
        L = (2180*a-768*a*b)/1463
        M = (16350-5760*b)/133
        Z = (-s.Integer(14323)-41856*b)/133133
        values = [s.Integer(1), L, M, L, Z]
        low, high = 0, 16
    else:
        L = (-55*a-4*a*b)/175
        M = (19+8*b)/5
        Z = (-725*a-516*a*b)/35035
        values = [s.Integer(1), L, M, Z]
        low, high = 2, 14
    support = list(range(low, high+1, 4))
    common = (-1)**low * s.sqrt(s.binomial(16, low)) * coefficients[low]
    for r, value in zip(support, values):
        zero((-1)**r*s.sqrt(s.binomial(16, r))*coefficients[r]-common*value,
             "Exact reduced-polynomial coefficient")
    for r in set(range(17))-set(support):
        zero(coefficients[r], "Exact sparse support")
    field = s.QQ.algebraic_field(a, b)
    poly = s.Poly(sum(c*y**r for r, c in enumerate(values)), y, domain=field)
    gcd = poly.gcd(poly.diff())
    assert gcd.degree() == 0 and poly.eval(0) == 1
    _, factors = s.sqf_list(poly)
    assert len(factors) == 1 and factors[0][1] == 1
    # P(z)=common*z^low*Q(z^4), with Q(0)!=0. Since Q is square free,
    # all nonzero finite roots are simple, and the projective poles have
    # the exact multiplicities low and 16-high.
    points = []
    if low:
        points.append(helpers.point(0j, low))
    if high < 16:
        points.append(helpers.point(None, 16-high))
    reduced_roots = s.nroots(poly, n=50, maxsteps=300)
    for root in reduced_roots:
        root = complex(root)
        radius = abs(root)**0.25
        angle = cmath.phase(root)
        for k in range(4):
            points.append(helpers.point(cmath.rect(radius, (angle+2*math.pi*k)/4), 1))
    assert sum(p["multiplicity"] for p in points) == 16
    expected_partition = [1]*16 if index == 0 else [2, 2]+[1]*12
    assert sorted((p["multiplicity"] for p in points), reverse=True) == expected_partition
    restored = helpers.reconstruct(points, 16)
    error = helpers.ray_error(saved, restored)
    assert error < 1e-12
    unit_error = max(abs(sum(p[key]**2 for key in ("x", "y", "z"))-1) for p in points)
    separation = min(math.sqrt(sum((p[key]-q[key])**2 for key in ("x", "y", "z")))
                     for i, p in enumerate(points) for q in points[i+1:])
    assert unit_error < 1e-14 and separation > 0.1
    diagnostics = {
        "coefficientSource": "saved exact radical Dicke coefficients",
        "squareFreeFactorizationExact": True,
        "factorizationField": "QQ(sqrt(385), I*sqrt(10))",
        "reducedPolynomialGcdWithDerivative": str(gcd.as_expr()),
        "zeroMultiplicity": low, "infinityMultiplicity": 16-high,
        "nonzeroRootStep": 4,
        "reducedPolynomialFactors": [{"polynomial": str(poly.as_expr()), "multiplicity": 1}],
        "projectiveDickeReconstructionError": error,
        "exactFormulaVsSavedCoefficientsError": formula_error,
        "maxUnitSphereError": unit_error,
        "minimumDistinctPointDistance": separation,
        "rootPrecisionDecimalDigits": 50,
        "noNumericalRootClustering": True,
    }
    return {"id": original["id"], "n": 16, "label": original["label"],
            "points": points, "diagnostics": diagnostics,
            "sourceExactCoefficients": exact_strings,
            "partition": expected_partition,
            "coordinateMethod": "Explicit quartic algebraic field, exact square-free factorization and 50-digit numerical roots.",
            "multiplicityEvidence": "Exact gcd and sparse homogeneous polynomial support."}


def main():
    raw = (HERE/"basis.json").read_bytes()
    source = json.loads(raw)
    records = {record["id"]: record for record in
               (representative(original, i) for i, original in enumerate(source["basis"]))}
    report = {"passed": True, "sourceBasisSha256": hashlib.sha256(raw).hexdigest(),
              "recordsChecked": len(records), "exactCoefficientIdentitiesChecked": 34,
              "exactSquareFreeReducedPolynomialsChecked": 2,
              "partitions": {key: record["partition"] for key, record in records.items()},
              "maximumProjectiveDickeReconstructionError": max(r["diagnostics"]["projectiveDickeReconstructionError"] for r in records.values()),
              "scope": "Multiplicities are exact. Sphere coordinates are numerical renderings, independently checked by reconstruction; roots are never clustered."}
    output = {"schemaVersion": 1,
              "convention": {"polynomial": "P(z)=sum_r (-1)^r sqrt(binomial(n,r))*c_r*z^r",
                             "sphere": "(2Re(z), -2Im(z), |z|^2-1)/(1+|z|^2); infinity is north."},
              "records": records, "verification": report}
    (HERE/"geometry.json").write_text(json.dumps(output, indent=2)+"\n")
    (HERE/"geometry-verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
