#!/usr/bin/env python3
"""Independent portable audit of serialized Liège AC refinements.

Inputs are the imported observations.json and the refiner's refined.json.
No refiner module is imported. AC residuals use a direct binomial sum. A small
residual is numerical evidence only, never an exact anticoherence certificate.

For the polynomial of the *serialized* vector, optional interval Rouché checks
certify one simple finite root in each source-centred projective disk. Structural
pole roots are handled as exact zero coefficients. This statement concerns that
nearby numerical vector and does not make its anticoherence exact.

Only Python's standard library and mpmath are required.
"""
import argparse
import hashlib
import json
import math
import re
import time
from collections import Counter
from decimal import Decimal, localcontext, ROUND_CEILING, ROUND_FLOOR
from pathlib import Path

import mpmath as mp


def pairs(values):
    return [mp.mpc(str(a), str(b)) for a, b in values]


def sf(value, digits=30):
    return mp.nstr(value, digits)


def bound_string(value, direction="upper", digits=30):
    """Outward decimal rounding of an exact binary mpmath endpoint."""
    sign, mantissa, exponent, _ = value._mpf_
    numerator = (-1 if sign else 1)*mantissa*(1 << max(exponent,0))
    denominator = 1 << max(-exponent,0)
    with localcontext() as context:
        context.prec = digits
        context.rounding = ROUND_CEILING if direction == "upper" else ROUND_FLOOR
        return str(Decimal(numerator)/Decimal(denominator))


def lower(value):
    return mp.mpf(value._mpi_[0])


def upper(value):
    return mp.mpf(value._mpi_[1])


def direct_reduction(c, t):
    """rho_ab = sqrt(C(t,a)C(t,b)) sum_k C(n-t,k) z_(a+k) z*_(b+k)."""
    n = len(c) - 1
    z = [c[r] / mp.sqrt(math.comb(n, r)) for r in range(n + 1)]
    factors = [mp.sqrt(math.comb(t, a)) for a in range(t + 1)]
    tail = [math.comb(n - t, k) for k in range(n - t + 1)]
    max_error = mp.mpf(0)
    square_error = mp.mpf(0)
    trace = mp.mpf(0)
    for a in range(t + 1):
        for b in range(a + 1):
            value = factors[a] * factors[b] * mp.fsum(
                tail[k] * z[a + k] * mp.conj(z[b + k])
                for k in range(n - t + 1)
            )
            if a == b:
                trace += value.real
                value -= mp.mpf(1) / (t + 1)
            error = abs(value)
            max_error = max(max_error, error)
            square_error += error ** 2 * (1 if a == b else 2)
    norm = mp.fsum(abs(x) ** 2 for x in c)
    return {
        "normSquared": sf(norm, 100),
        "normError": sf(abs(norm - 1)),
        "traceError": sf(abs(trace - 1)),
        "traceMinusNorm": sf(abs(trace - norm)),
        "maximumEntryResidual": sf(max_error),
        "HSResidual": sf(mp.sqrt(square_error)),
        "squaredHSResidual": sf(square_error),
    }


def changes(before, after):
    bnorm = mp.fsum(abs(x) ** 2 for x in before)
    anorm = mp.fsum(abs(x) ** 2 for x in after)
    overlap = mp.fsum(mp.conj(a) * b for a, b in zip(before, after))
    phase = mp.conj(overlap) / abs(overlap) if overlap else mp.mpc(1)
    differences = [b * phase - a for a, b in zip(before, after)]
    return {
        "maximumRawComponentChange": sf(max(abs(b - a) for a, b in zip(before, after))),
        "maximumPhaseAlignedComponentChange": sf(max(map(abs, differences))),
        "phaseAlignedL2Change": sf(mp.sqrt(mp.fsum(abs(x) ** 2 for x in differences))),
        "normalizedInfidelity": sf(max(mp.mpf(0), 1 - abs(overlap) ** 2 / (anorm * bnorm))),
        "alignedPhase": [sf(phase.real), sf(phase.imag)],
    }


def exact_zero(pair):
    return all(Decimal(str(x)) == 0 for x in pair)


def significant_digits(value):
    d = Decimal(str(value))
    return 0 if not d else len(d.as_tuple().digits)


def pole_counts(points):
    south = north = 0
    for p in points:
        x, y, z = (Decimal(str(a)) for a in p)
        if x == 0 and y == 0:
            south += int(z == -1)
            north += int(z == 1)
    return south, north


def source_centres(points, removed_south, removed_north):
    """The polynomial is product_i(a_i+b_i z); use |z|<=1 charts."""
    centres = []
    skip_south, skip_north = removed_south, removed_north
    for i, p in enumerate(points):
        x, y, z = (mp.mpf(str(a)) for a in p)
        if x == y == 0 and z == -1 and skip_south:
            skip_south -= 1
            continue
        if x == y == 0 and z == 1 and skip_north:
            skip_north -= 1
            continue
        scale = mp.sqrt(x*x + y*y + z*z)
        x, y, z = x/scale, y/scale, z/scale
        if z < 0:
            centre = -(x - mp.j*y)/(1-z)
            chart = "z"
        else:
            centre = -(x + mp.j*y)/(1+z)
            chart = "reciprocal"
        # Rational decimal centres, so interval evaluation includes no assumed
        # transcendental root coordinates. 55 digits exceed source precision.
        centre_strings = [sf(centre.real, 55), sf(centre.imag, 55)]
        ix, iy, iz = (mp.iv.mpf(str(a)) for a in p)
        iscale = mp.iv.sqrt(ix*ix+iy*iy+iz*iz)
        ix, iy, iz = ix/iscale, iy/iscale, iz/iscale
        ichart = -(mp.iv.mpc(ix,-iy))/(1-iz) if chart == "z" else -(mp.iv.mpc(ix,iy))/(1+iz)
        centre_error = upper(abs(ichart-interval_pair(centre_strings)))
        centres.append({"sourcePointIndex": i, "chart": chart, "centre": centre_strings,
                        "sourceCentreApproximationErrorUpper":bound_string(centre_error)})
    assert skip_south == skip_north == 0
    return centres


def interval_pair(pair):
    return mp.iv.mpc(mp.iv.mpf(str(pair[0])), mp.iv.mpf(str(pair[1])))


def taylor_coefficients(coefficients, centre):
    """Horner shift, retaining interval enclosures of every Taylor coefficient."""
    shifted = [coefficients[-1]]
    for coefficient in reversed(coefficients[:-1]):
        shifted = [centre*shifted[0] + coefficient] + [
            centre*shifted[k] + shifted[k-1] for k in range(1, len(shifted))
        ] + [shifted[-1]]
    return shifted


def root_disks(serialized, points, removed_south, removed_north):
    """Interval Rouché checks and interval-verified disjoint spherical caps.

    If P(c+w)=sum B_k w^k and |B0|+sum_(k>=2)|Bk|R^k < |B1|R,
    exactly one root, counted with multiplicity, is in |w|<R. The inverse
    stereographic map sends this disk into a spherical cap of chord radius2R.
    Thus centre separation >2(Ri+Rj) suffices even across opposite charts.
    """
    n = len(serialized)-1
    stop = n+1-removed_north
    q = [interval_pair(serialized[k])*mp.iv.sqrt(math.comb(n,k))
         for k in range(removed_south, stop)]
    centres = source_centres(points, removed_south, removed_north)
    degree = n-removed_south-removed_north
    assert len(centres) == degree
    assert not exact_zero(serialized[removed_south])
    assert not exact_zero(serialized[n-removed_north])
    if degree == 0:
        return {"passed": True, "finiteSimpleRoots": 0,
                "intervalPrecisionDigits": mp.iv.dps, "disks": [],
                "maximumSphericalChordDisplacementBound": "0",
                "minimumRoucheRelativeMargin": None,
                "minimumDisjointnessRatio": None,
                "meaning": "Only exact pole roots occur in this serialized polynomial."}
    disks = []
    for item in centres:
        centre = interval_pair(item["centre"])
        shifted = taylor_coefficients(q if item["chart"] == "z" else q[::-1], centre)
        absolute = [abs(a) for a in shifted]
        first_lower = lower(absolute[1])
        if first_lower <= 0:
            raise AssertionError(f"Derivative interval contains zero: {item}")
        eta = upper(absolute[0]) / first_lower
        exponent = max(-45, int(mp.ceil(mp.log10(4*eta))) if eta else -45)
        success = False
        for attempt in range(8):
            radius_string = "1e"+str(exponent+attempt)
            radius = mp.iv.mpf(radius_string)
            higher = mp.iv.mpf(0)
            power = radius*radius
            for a in absolute[2:]:
                higher += a*power
                power *= radius
            lhs = absolute[0]+higher
            rhs = absolute[1]*radius
            ratio = upper(lhs/rhs)
            if upper(lhs) < lower(rhs):
                success = True
                break
        if not success:
            raise AssertionError(f"No interval Rouche disk: {item}, ratio={ratio}")
        disks.append({**item, "radius": radius_string,
                      "roucheLeftUpper": bound_string(upper(lhs)),
                      "roucheRightLower": bound_string(lower(rhs),"lower"),
                      "roucheRatioUpper": bound_string(ratio), "simpleRootCount": 1})
    minimum_ratio = mp.inf
    for i, a in enumerate(disks):
        ac = interval_pair(a["centre"])
        au, av = (ac, mp.iv.mpc(1)) if a["chart"] == "z" else (mp.iv.mpc(1), ac)
        for b in disks[:i]:
            bc = interval_pair(b["centre"])
            bu, bv = (bc, mp.iv.mpc(1)) if b["chart"] == "z" else (mp.iv.mpc(1), bc)
            determinant = au*bv-bu*av
            separation = 2*abs(determinant)/mp.iv.sqrt(
                (abs(au)**2+abs(av)**2)*(abs(bu)**2+abs(bv)**2))
            cap_sum = 2*(mp.iv.mpf(a["radius"])+mp.iv.mpf(b["radius"]))
            ratio = lower(separation/cap_sum)
            minimum_ratio = min(minimum_ratio, ratio)
            if lower(separation) <= upper(cap_sum):
                raise AssertionError(f"Root disks not proven disjoint: {a},{b}")
    # Separate finite-root caps from the structural pole roots as well.
    for a in disks:
        ac = interval_pair(a["centre"])
        au, av = (ac, mp.iv.mpc(1)) if a["chart"] == "z" else (mp.iv.mpc(1), ac)
        for pu, pv, multiplicity in ((0,1,removed_south),(1,0,removed_north)):
            if not multiplicity:
                continue
            separation = 2*abs(au*pv-pu*av)/mp.iv.sqrt(abs(au)**2+abs(av)**2)
            cap = 2*mp.iv.mpf(a["radius"])
            ratio = lower(separation/cap)
            minimum_ratio = min(minimum_ratio, ratio)
            assert lower(separation) > upper(cap), "Finite root cap overlaps a structural pole"
    return {
        "passed": True, "finiteSimpleRoots": degree, "intervalPrecisionDigits": mp.iv.dps,
        "maximumSphericalChordDisplacementBound": bound_string(max(upper(2*(mp.iv.mpf(a["radius"])+
              mp.iv.mpf(a["sourceCentreApproximationErrorUpper"]))) for a in disks)),
        "minimumRoucheRelativeMargin": bound_string(min(lower(1-mp.iv.mpf(a["roucheRatioUpper"])) for a in disks),"lower"),
        "minimumDisjointnessRatio": bound_string(minimum_ratio,"lower") if mp.isfinite(minimum_ratio) else None,
        "meaning": "One simple root per disjoint source-centred disk for the serialized numerical vector; this does not certify exact AC.",
        "disks": disks,
    }


def audit_record(source, refined, target, do_roots):
    started = time.time()
    n = int(refined["n"])
    t = int(refined.get("reportedOrder", source["claimedOrder"]))
    assert n == source["n"] and t == source["claimedOrder"]
    original_serialized = refined["original90DigitCoefficients"]
    final_serialized = refined["refinedCoefficients"]
    assert len(original_serialized) == len(final_serialized) == n+1
    assert all(isinstance(x, str) for p in final_serialized for x in p), "Serialized decimals must be strings"
    before, after = pairs(original_serialized), pairs(final_serialized)
    source_stored = pairs(source["highPrecisionCoefficients"])
    source_match = changes(source_stored, before)
    assert mp.mpf(source_match["maximumPhaseAlignedComponentChange"]) < mp.mpf("1e-33"), "Original vector does not match the imported source"
    original_check = direct_reduction(before, t)
    final_check = direct_reduction(after, t)
    assert mp.mpf(final_check["HSResidual"]) <= target, "Reported-order residual exceeds target"
    assert mp.mpf(final_check["normError"]) <= mp.sqrt(t+1)*target, "Norm residual exceeds trace bound"
    asserted_zero_indices = set(refined.get("structuralZeroIndices", []))
    assert all(exact_zero(final_serialized[k]) for k in asserted_zero_indices), "A structural zero changed"
    south, north = pole_counts(source["rawUnitPoints"])
    required_south = south if south >= 2 else 0
    required_north = north if north >= 2 else 0
    expected_zeros = set(range(required_south)) | set(range(n-required_north+1, n+1))
    assert expected_zeros <= asserted_zero_indices, "Repeated source pole was not constrained"
    assert all(exact_zero(original_serialized[k]) for k in expected_zeros), "Original pole convention disagrees"
    assert all(exact_zero(final_serialized[k]) for k in expected_zeros), "Repeated source pole changed"
    low = 0
    while low <= n and exact_zero(final_serialized[low]): low += 1
    high = 0
    while high <= n-low and exact_zero(final_serialized[n-high]): high += 1
    if required_south: assert low == required_south, "South multiplicity changed"
    if required_north: assert high == required_north, "North multiplicity changed"
    source_partition = sorted(source["partition"], reverse=True)
    inferred_partition = sorted([a for a in (required_south, required_north) if a] +
                                [1]*(n-required_south-required_north), reverse=True)
    assert inferred_partition == source_partition, "Source has a non-pole repeated cluster not covered by this audit"
    precision = Counter(significant_digits(x) for p in final_serialized for x in p if Decimal(x))
    roots = root_disks(final_serialized, refined.get("originalPoints",source["rawUnitPoints"]), required_south, required_north) if do_roots else {"passed": None, "skipped": True}
    return {
        "id": source["id"], "n": n, "reportedOrder": t, "passed": True,
        "before": original_check, "after": final_check, "change": changes(before, after),
        "sourceStoredResidual": source["ACCheck"]["frobeniusError"],
        "originalVectorComparisonToImportedSource":source_match,
        "singularFamilyCorrection":refined.get("singularFamilyCorrection"),
        "serializedNonzeroComponentDigitCounts": dict(sorted(precision.items())),
        "sourcePartition": source_partition, "refinedPartition": inferred_partition,
        "structuralPoleZeros": {"southMultiplicity": required_south,
                                "northMultiplicity": required_north,
                                "indices": sorted(expected_zeros), "preservedExactly": True},
        "finiteRootVerification": roots,
        "runtimeSeconds": round(time.time()-started,3),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--refined", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path(__file__).with_name("verification.json"))
    parser.add_argument("--dps", type=int, default=110)
    parser.add_argument("--target", default="1e-60")
    parser.add_argument("--skip-roots", action="store_true")
    parser.add_argument("--only", nargs="*")
    args = parser.parse_args()
    assert args.dps >= 90
    mp.mp.dps = args.dps
    mp.iv.dps = args.dps
    original_bytes = args.original.read_bytes()
    refined_bytes = args.refined.read_bytes()
    original = json.loads(original_bytes)
    refined = json.loads(refined_bytes)
    original_records = original if isinstance(original, list) else original["records"]
    refined_records = refined if isinstance(refined, list) else refined["records"]
    source_by_id = {r["id"]: r for r in original_records}
    assert len(source_by_id) == len(original_records), "Duplicate source IDs"
    ids = [r["id"] for r in refined_records]
    assert len(set(ids)) == len(ids), "Duplicate refined IDs"
    if not args.only:
        assert set(ids) == set(source_by_id), "Source/refined inventory mismatch"
    selected = [r for r in refined_records if not args.only or r["id"] in args.only]
    results = []
    failures = []
    start = time.time()
    for r in selected:
        try:
            result = audit_record(source_by_id[r["id"]], r, mp.mpf(args.target), not args.skip_roots)
            results.append(result)
            print(r["id"], "PASS", "HS", result["after"]["HSResidual"],
                  "roots", result["finiteRootVerification"]["passed"], flush=True)
        except Exception as error:
            failures.append({"id":r["id"], "error":repr(error)})
            print(r["id"], "FAIL", repr(error), flush=True)
    summary = {
        "expectedRecords": len(original_records), "selectedRecords":len(selected),
        "passedRecords":len(results), "failedRecords":len(failures),
        "targetHSResidual":args.target, "workingPrecisionDigits":args.dps,
        "allRootsChecked":not args.skip_roots,
        "allSelectedPassed":not failures and len(results)==len(selected),
        "completeInventory":set(ids)==set(source_by_id) and not args.only,
        "maximumBeforeHSResidual":sf(max((mp.mpf(r["before"]["HSResidual"]) for r in results),default=mp.mpf(0))),
        "maximumAfterHSResidual":sf(max((mp.mpf(r["after"]["HSResidual"]) for r in results),default=mp.mpf(0))),
        "maximumAfterNormError":sf(max((mp.mpf(r["after"]["normError"]) for r in results),default=mp.mpf(0))),
        "maximumPhaseAlignedComponentChange":sf(max((mp.mpf(r["change"]["maximumPhaseAlignedComponentChange"]) for r in results),default=mp.mpf(0))),
        "maximumPhaseAlignedL2Change":sf(max((mp.mpf(r["change"]["phaseAlignedL2Change"]) for r in results),default=mp.mpf(0))),
        "maximumRootSphericalChordDisplacementBound":str(max((Decimal(r["finiteRootVerification"].get("maximumSphericalChordDisplacementBound","0")) for r in results),default=Decimal(0))),
        "runtimeSeconds":round(time.time()-start,3),
    }
    output = {
        "schemaVersion":1, "scientificStatus":"Numerical anticoherence verification, with interval root-preservation checks for serialized coefficients. No exact AC existence or maximality proof is inferred.",
        "method":"Independent direct binomial reduced-state sum; exact serialized pole-zero checks; interval Taylor/Rouche root disks in two stereographic charts.",
        "inputs": {"originalFilename":args.original.name,"originalSHA256":hashlib.sha256(original_bytes).hexdigest(),
                   "refinedFilename":args.refined.name,"refinedSHA256":hashlib.sha256(refined_bytes).hexdigest()},
        "summary":summary,"records":results,"failures":failures,
    }
    args.out.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)
    if failures: raise SystemExit(1)


if __name__ == "__main__":
    main()
