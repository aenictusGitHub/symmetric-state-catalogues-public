#!/usr/bin/env python3
"""Check sector lengths from normalized Dicke coefficients.

This checks that a representative state attains the recorded value. It does
not prove that the value is a global minimum or maximum.

Requires Python 3.10+ and NumPy. Examples:
    python3 scripts/verify_states.py --self-test
    python3 scripts/verify_states.py data/catalogue.json

Expected JSON: {"entries": [{"id": "coherent-example", "n": 4, "k": 2,
  "numericValue": 6, "stateKey": "coherent"}]}

Representation alternatives: stateKey = "coherent", "ghz", or "dicke".
For "dicke", also provide the integer excitation. Representation fields may
be nested in "state". Explicit coefficients must be normalized and have n+1
entries; each entry is a real number or [real, imaginary] pair. Numeric
targets may be in numericValue or value. N is accepted as an alias for n.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np


def symmetric_purities(coefficients: np.ndarray) -> np.ndarray:
    """Return P_r = tr(rho_r^2) for r=0,...,n in the Dicke basis."""
    c = np.asarray(coefficients, dtype=complex)
    if c.ndim != 1 or len(c) < 2:
        raise ValueError("Expected a one-dimensional vector of n+1 amplitudes, n >= 1")
    if not np.all(np.isfinite(c)):
        raise ValueError("Coefficients must be finite")
    norm = float(np.vdot(c, c).real)
    if not math.isclose(norm, 1.0, rel_tol=0.0, abs_tol=1e-7):
        raise ValueError(f"Dicke coefficients are not normalized (norm squared = {norm:.12g})")
    n = len(c) - 1
    purities = []
    for r in range(n + 1):
        matrix = np.empty((r + 1, n - r + 1), dtype=complex)
        for a in range(r + 1):
            for b in range(n - r + 1):
                matrix[a, b] = c[a + b] * math.sqrt(
                    math.comb(r, a) * math.comb(n - r, b) / math.comb(n, a + b)
                )
        reduced = matrix @ matrix.conj().T
        purities.append(float(np.sum(np.abs(reduced) ** 2)))
    return np.asarray(purities)


def sector_lengths(coefficients: np.ndarray) -> np.ndarray:
    """Return unnormalized squared Pauli sector lengths S_0,...,S_n.

    S_k = C(n,k) sum_(r=0)^k (-1)^(k-r) C(k,r) 2^r P_r.
    No negative-value clipping or other correction is applied.
    """
    purities = symmetric_purities(coefficients)
    n = len(purities) - 1
    return np.asarray([
        math.comb(n, k) * math.fsum(
            (-1) ** (k - r) * math.comb(k, r) * 2**r * purities[r]
            for r in range(k + 1)
        )
        for k in range(n + 1)
    ])


def _integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{label} must be an integer")
    return value


def coefficients_for(entry: dict[str, Any], n: int) -> np.ndarray | None:
    """Read an explicit representative or construct a supported family."""
    state = entry.get("state")
    representation = state if isinstance(state, dict) else entry
    coefficients = representation.get("coefficients", entry.get("coefficients"))
    if coefficients is not None:
        if not isinstance(coefficients, list) or len(coefficients) != n + 1:
            raise ValueError(f"Explicit coefficients must contain {n + 1} amplitudes")
        amplitudes = []
        for x in coefficients:
            if isinstance(x, (int, float)) and not isinstance(x, bool):
                amplitudes.append(complex(x))
            elif isinstance(x, list) and len(x) == 2 and all(
                isinstance(y, (int, float)) and not isinstance(y, bool) for y in x
            ):
                amplitudes.append(complex(*x))
            else:
                raise ValueError("Each coefficient must be a real number or [real, imaginary]")
        return np.asarray(amplitudes, dtype=complex)

    state_key = representation.get("stateKey", entry.get("stateKey"))
    if state_key is None:
        return None
    c = np.zeros(n + 1, dtype=complex)
    if state_key == "coherent":
        c[0] = 1
    elif state_key == "ghz":
        c[0] = c[-1] = 1 / math.sqrt(2)
    elif state_key == "dicke":
        excitation = _integer(representation.get("excitation", entry.get("excitation")), "excitation")
        if not 0 <= excitation <= n:
            raise ValueError("Dicke excitation must lie between 0 and n")
        c[excitation] = 1
    else:
        raise ValueError(f"Unsupported stateKey: {state_key!r}")
    return c


def direct_pauli_sectors(coefficients: np.ndarray) -> np.ndarray:
    """Independent small-n reference: enumerate every Pauli tensor product."""
    n = len(coefficients) - 1
    state = np.asarray([
        coefficients[i.bit_count()] / math.sqrt(math.comb(n, i.bit_count()))
        for i in range(2**n)
    ])
    paulis = [
        np.eye(2), np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]),
    ]
    sectors = np.zeros(n + 1)
    for indices in itertools.product(range(4), repeat=n):
        operator = np.ones((1, 1), dtype=complex)
        for index in indices:
            operator = np.kron(operator, paulis[index])
        expectation = np.vdot(state, operator @ state)
        if abs(expectation.imag) > 1e-12:
            raise AssertionError("Hermitian Pauli expectation unexpectedly complex")
        sectors[sum(index != 0 for index in indices)] += expectation.real**2
    return sectors


def self_test() -> None:
    for n in range(2, 13):
        coherent = coefficients_for({"stateKey": "coherent"}, n)
        expected = np.asarray([math.comb(n, k) for k in range(n + 1)])
        np.testing.assert_allclose(sector_lengths(coherent), expected, rtol=1e-9, atol=1e-8)
        ghz = coefficients_for({"stateKey": "ghz"}, n)
        expected = np.asarray([
            math.comb(n, k) if k % 2 == 0 else 0 for k in range(n + 1)
        ], dtype=float)
        expected[-1] = 2 ** (n - 1) + int(n % 2 == 0)
        np.testing.assert_allclose(sector_lengths(ghz), expected, rtol=1e-9, atol=1e-8)
    print("PASS coherent and GHZ sector identities for n = 2,...,12")

    rng = np.random.default_rng(271828)
    for n in (2, 3, 4):
        coefficients = rng.normal(size=n + 1) + 1j * rng.normal(size=n + 1)
        coefficients /= np.linalg.norm(coefficients)
        sectors = sector_lengths(coefficients)
        np.testing.assert_allclose(sectors, direct_pauli_sectors(coefficients), rtol=1e-9, atol=1e-10)
        np.testing.assert_allclose(symmetric_purities(coefficients), symmetric_purities(coefficients)[::-1], atol=1e-12)
        np.testing.assert_allclose(sectors.sum(), 2**n, atol=1e-10)
    print("PASS direct Pauli enumeration for complex symmetric states, n = 2,3,4")

    t = (83 + 9 * math.sqrt(669)) / 430
    u = (37052 + 171 * math.sqrt(669)) / 210700
    c = np.zeros(10, dtype=complex)
    c[[1, 4, 7]] = [1, math.sqrt(14 * u), 2 * t]
    c /= np.linalg.norm(c)
    s8 = sector_lengths(c)[8]
    np.testing.assert_allclose(s8, 159.441403981044, rtol=1e-10, atol=1e-9)
    print(f"PASS n = 9 explicit representative: S_8 = {s8:.12f}")

    try:
        sector_lengths(np.ones(3))
    except ValueError:
        pass
    else:
        raise AssertionError("Unnormalized coefficient input was not rejected")
    print("PASS malformed normalization is rejected")


def verify_catalogue(path: Path, atol: float, rtol: float) -> int:
    with path.open(encoding="utf-8") as stream:
        catalogue = json.load(stream)
    entries = catalogue.get("entries") if isinstance(catalogue, dict) else catalogue
    if not isinstance(entries, list):
        raise ValueError("Catalogue must be a list or an object with an entries list")
    passed = skipped = failed = 0
    for position, entry in enumerate(entries, start=1):
        if not isinstance(entry, dict):
            print(f"FAIL entry {position}: expected a JSON object")
            failed += 1
            continue
        label = entry.get("id", f"entry {position}")
        try:
            n = _integer(entry.get("n", entry.get("N")), "n")
            k = _integer(entry.get("k"), "k")
            if not 1 <= n <= 50 or not 0 <= k <= n:
                raise ValueError("Require 1 <= n <= 50 and 0 <= k <= n")
            c = coefficients_for(entry, n)
            target = entry.get("numericValue", entry.get("value"))
            if c is None or not isinstance(target, (int, float)) or isinstance(target, bool):
                print(f"SKIP {label}: no supported representative and numeric target pair")
                skipped += 1
                continue
            if not math.isfinite(target):
                raise ValueError("Target value must be finite")
            computed = float(sector_lengths(c)[k])
            # Rounded source values carry a declared, per-entry precision.
            entry_atol = entry.get("absoluteTolerance", atol)
            if not isinstance(entry_atol, (int, float)) or not math.isfinite(entry_atol) or entry_atol < 0:
                raise ValueError("absoluteTolerance must be a finite nonnegative number")
            if not math.isclose(computed, target, abs_tol=max(atol, entry_atol), rel_tol=rtol):
                raise ValueError(f"S_{k} = {computed:.12g}; recorded {target:.12g}; difference {computed-target:.3g}")
            precision = f" (source rounding tolerance {entry_atol:g})" if entry_atol > atol else ""
            print(f"PASS {label}: S_{k} = {computed:.12g}{precision}")
            passed += 1
        except (ValueError, TypeError, OverflowError) as error:
            print(f"FAIL {label}: {error}")
            failed += 1
    print(f"\nAttainment checks: {passed} passed, {skipped} skipped, {failed} failed.")
    print("These checks verify representative-state values, not global optimality.")
    return 1 if failed else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("catalogue", nargs="?", type=Path, help="Catalogue JSON (defaults to data/catalogue.json)")
    parser.add_argument("--self-test", action="store_true", help="Run analytic and independent Pauli checks; no catalogue required")
    parser.add_argument("--atol", type=float, default=1e-7, help="Absolute attainment tolerance (default 1e-7)")
    parser.add_argument("--rtol", type=float, default=1e-9, help="Relative attainment tolerance (default 1e-9)")
    args = parser.parse_args()
    if args.atol < 0 or args.rtol < 0 or not math.isfinite(args.atol + args.rtol):
        parser.error("Tolerances must be finite and nonnegative")
    if args.self_test:
        self_test()
        if args.catalogue is None:
            return 0
    path = args.catalogue or Path(__file__).resolve().parents[1] / "data" / "catalogue.json"
    return verify_catalogue(path, args.atol, args.rtol)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, AssertionError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
