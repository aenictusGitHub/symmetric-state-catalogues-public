#!/usr/bin/env python3
"""Exact multiplicities and numerical Majorana rendering for the N=17 radicals.

The explicit quartic number field avoids expensive automatic primitive-element
discovery from unsimplified quotients of Dicke coefficients.
"""
from pathlib import Path
import cmath
import hashlib
import importlib.util
import json
import math

import numpy as np
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


def branch(source, sign):
    original = source["basis"][0]
    exact_strings = original["exactCoefficients"] if sign == 1 else [
        str(s.conjugate(s.sympify(v))) for v in original["exactCoefficients"]]
    coefficients = list(map(s.sympify, exact_strings))
    saved = [complex(a, sign * b) for a, b in original["coefficients"]]
    exact_numeric = [complex(s.N(v, 50)) for v in coefficients]
    formula_error = helpers.ray_error(saved, exact_numeric)
    assert formula_error < 1e-14
    q = s.sqrt(20995)
    u = 710-9*q
    w = sign*s.I*s.sqrt(12780*q-575935)
    L = -s.Rational(2,277)*(u+w)
    M = -s.Rational(16,1628760)*(2130+u+2*w)*(u+w)
    values = [1,L,M,-M,-L,-1]
    common = -s.sqrt(17)*coefficients[1]
    for r, value in zip((1,4,7,10,13,16), values):
        zero((-1) ** r * s.sqrt(s.binomial(17, r)) * coefficients[r]
             - common * value, "Exact reduced-polynomial coefficient")
    for r in set(range(18)) - {1, 4, 7, 10, 13,16}:
        zero(coefficients[r], "Exact sparse support")
    field = s.QQ.algebraic_field(w)
    poly = s.Poly(sum(c * y ** r for r, c in enumerate(values)), y, domain=field)
    gcd = poly.gcd(poly.diff())
    if gcd.degree() != 0:
        raise AssertionError("The quintic must be square free")
    _, factors = s.sqf_list(poly)
    assert len(factors) == 1 and factors[0][1] == 1
    assert poly.eval(0) == 1
    # P(z) = common*z*Q(z^3) has one simple zero, fifteen simple nonzero
    # roots, and one simple projective root at infinity. Q(0)!=0 and
    # gcd(Q,Q')=1 ensure Q(z^3) remains square free away from zero.
    points = [helpers.point(0j, 1), helpers.point(None, 1)]
    reduced_roots = s.nroots(poly, n=50, maxsteps=300)
    for root in reduced_roots:
        root = complex(root)
        radius = abs(root) ** (1 / 3)
        angle = cmath.phase(root)
        for k in range(3):
            points.append(helpers.point(cmath.rect(radius, (angle + 2 * math.pi * k) / 3), 1))
    assert len(points) == 17 and sum(p["multiplicity"] for p in points) == 17
    restored = helpers.reconstruct(points, 17)
    error = helpers.ray_error(saved, restored)
    assert error < 1e-12
    unit_error = max(abs(sum(p[key] ** 2 for key in ("x", "y", "z")) - 1) for p in points)
    separation = min(math.sqrt(sum((p[key] - q[key]) ** 2 for key in ("x", "y", "z")))
                     for i, p in enumerate(points) for q in points[i + 1:])
    assert unit_error < 1e-14 and separation > 0.1
    suffix = "" if sign == 1 else "-conjugate"
    record_id = "radical-n17-k1-t4" + suffix + "-basis-1"
    diagnostics = {
        "coefficientSource": "saved exact radical Dicke coefficients",
        "squareFreeFactorizationExact": True,
        "factorizationField": "QQ(I*sqrt(12780*sqrt(20995)-575935))",
        "reducedPolynomialGcdWithDerivative": str(gcd.as_expr()),
        "zeroMultiplicity": 1, "infinityMultiplicity": 1, "nonzeroRootStep": 3,
        "reducedPolynomialFactors": [{"polynomial": str(poly.as_expr()), "multiplicity": 1}],
        "projectiveDickeReconstructionError": error,
        "exactFormulaVsSavedCoefficientsError": formula_error,
        "maxUnitSphereError": unit_error,
        "minimumDistinctPointDistance": separation,
        "rootPrecisionDecimalDigits": 50,
        "noNumericalRootClustering": True,
    }
    return {"id": record_id, "n": 17,
            "label": "Exact fourth-order AC state · " + ("positive" if sign == 1 else "conjugate") + " branch",
            "points": points, "diagnostics": diagnostics,
            "sourceExactCoefficients": exact_strings,
            "partition": [1] * 17,
            "coordinateMethod": "Explicit quartic-number-field algebraic field, exact square-free factorization and 50-digit numerical roots.",
            "multiplicityEvidence": "Exact gcd and sparse homogeneous polynomial support."}


def main():
    raw = (HERE / "basis.json").read_bytes()
    source = json.loads(raw)
    records = {record["id"]: record for record in (branch(source, 1), branch(source, -1))}
    report = {"passed": True, "sourceBasisSha256": hashlib.sha256(raw).hexdigest(),
              "recordsChecked": len(records),
              "exactCoefficientIdentitiesChecked": 36,
              "exactSquareFreeQuinticsChecked": 2,
              "allSeventeenRootsSimpleExactly": True,
              "maximumProjectiveDickeReconstructionError": max(r["diagnostics"]["projectiveDickeReconstructionError"] for r in records.values()),
              "scope": "Multiplicities are exact. Sphere coordinates are numerical renderings, independently checked by reconstruction; roots are never clustered."}
    output = {"schemaVersion": 1,
              "convention": {"polynomial": "P(z)=sum_r (-1)^r sqrt(binomial(n,r))*c_r*z^r",
                             "sphere": "(2Re(z), -2Im(z), |z|^2-1)/(1+|z|^2); infinity is north."},
              "records": records, "verification": report}
    (HERE / "geometry.json").write_text(json.dumps(output, indent=2) + "\n")
    (HERE / "geometry-verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
