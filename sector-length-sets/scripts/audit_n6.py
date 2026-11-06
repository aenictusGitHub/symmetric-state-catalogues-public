#!/usr/bin/env python3
"""Independently audit exported N=6 states with spin-multipole identities."""
import json
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[1]
data = json.loads((root / "data/n6.json").read_text())
jp = np.zeros((7, 7), complex)
for k in range(1, 7):
    jp[k - 1, k] = np.sqrt(k * (7 - k))
js = [(jp + jp.conj().T) / 2, (jp - jp.conj().T) / (2j), np.diag(np.arange(3, -4, -1))]
worst_norm = worst_spin = worst_linear = 0.0
min_refinement_gap = np.inf
for record in data["witnesses"]:
    c = np.array([complex(*v) for v in record["coefficients"]])
    sectors = np.array(record["sectors"])
    worst_norm = max(worst_norm, abs(np.vdot(c, c) - 1))
    dipole = np.array([np.vdot(c, j @ c).real for j in js])
    q = np.array([[np.vdot(c, ((a @ b + b @ a) / 2) @ c).real - 4 * (i == k)
                   for k, b in enumerate(js)] for i, a in enumerate(js)])
    theta = sum((-1) ** k * c[k] * c[6 - k] for k in range(7))
    x = 2 * np.dot(dipole, dipole) / 3
    y = 5 + 4 * np.sum(q * q) / 15
    z = 8 + 2 * x - 8 * abs(theta) ** 2
    worst_spin = max(worst_spin, np.max(np.abs(sectors[:3] - [x, y, z])))
    expected_tail = [45 + 10 * x - 2 * y - 3 * z, 3 * z - 9 * x, 18 - 2 * x + y - z]
    worst_linear = max(worst_linear, np.max(np.abs(sectors[3:] - expected_tail)))
    gap = 5 / 21 * (y + z - 5 * x - 5) + ((162 + 8 * x - 10 * y - 3 * z) / 231) ** 2 - (8 + 2 * x - z) * (401 + 137 * x + 35 * y + 5 * z) / 51744
    min_refinement_gap = min(min_refinement_gap, gap)

assert worst_norm < 1e-12
assert worst_spin < 1e-10
assert worst_linear < 1e-9
assert min_refinement_gap > -1e-11
for patch in data["outer"]["patches"]:
    assert len(patch["attained"]) == len(patch["triangles"])
    count = len(patch["positions"])
    assert all(0 <= index < count for triangle in patch["triangles"] for index in triangle)
    if patch["kind"] == "fourcopy_lower":
        assert not any(patch["attained"])
result = dict(status="PASS", named_states=len(data["states"]), pure_state_witnesses=len(data["witnesses"]),
              triangles=sum(len(p["triangles"]) for p in data["outer"]["patches"]),
              edge_runs=len(data["outer"]["edges"]), norm_error=float(worst_norm),
              independent_spin_error=float(worst_spin), complementary_sector_error=float(worst_linear),
              minimum_refinement_gap=float(min_refinement_gap))
(root / "data/n6_audit.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
