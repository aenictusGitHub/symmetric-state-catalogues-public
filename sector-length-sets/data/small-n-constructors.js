/* Exact state witnesses for continuous inspection of N=3,4 sets and N=5 boundary.
 * Dicke order k=0,...,N. See SMALL_N_SCIENTIFIC_NOTES.txt for provenance.
 * No N=5 interior inverse is asserted: its coverage proof uses connected fibers.
 */
(function (global) {
  'use strict';
  const root5 = Math.sqrt(5);
  const pc = 2 - 3 / root5;
  const p0 = root5 / 3;
  const xc = (22805 - 9840 * root5) / 11881;
  const clamp = (value, low, high) => Math.min(high, Math.max(low, value));
  const nearRange = (value, low, high) => Number.isFinite(value) && value >= low - 1e-10 && value <= high + 1e-10;
  function record(coefficients, sectors, label, parameter, source) {
    return { coefficients, sectors, label, parameter, status: 'attained', source };
  }
  function n3(x) {
    if (!nearRange(x, 0, 3)) return null;
    x = clamp(x, 0, 3);
    const spin = Math.sqrt(x / 3);
    const c = [[Math.sqrt((1 + spin) / 2), 0], [0, 0], [0, 0], [Math.sqrt((1 - spin) / 2), 0]];
    return record(c, [x, 3, 4 - x], '|GHZ(3, φ)⟩', { S1: x }, 'Exact entire three-qubit interval');
  }
  function n4(x, y) {
    if (!nearRange(x, 0, 4) || !nearRange(y, x + 2, 6)) return null;
    x = clamp(x, 0, 4);
    y = clamp(y, x + 2, 6);
    const Y = 3 - y / 2;
    // Rationalizing the inverse removes 0/0 at GHZ and cancellation near y=6.
    const t = Math.sqrt(Math.max(0, y - x - 2));
    const B2 = Y / (2 * (2 + t));
    const A = Math.sqrt(1 - B2), B = Math.sqrt(B2);
    const q = clamp(x / (4 + 2 * t - Y), 0, 1), v = 1 - q;
    const cosPhi = Math.sqrt(q * (2 - q));
    const cosHalfPhi = Math.sqrt((1 + cosPhi) / 2);
    const sinHalfPhi = v / (2 * cosHalfPhi);
    const c = [[A * cosHalfPhi, 0], [0, 0], [0, B], [0, 0], [A * sinHalfPhi, 0]];
    return record(c, [x, y, 12 + x - 2 * y, 3 - 2 * x + y], 'Four-qubit exact-set state',
      { S1: x, S2: y, theta: 2 * Math.atan2(B, A), phi: Math.atan2(v, cosPhi) }, 'PRA 113, 012415 (2026), Eqs.(47)–(50)');
  }
  function n5Sectors(x, y) { return [x, y, 10 + 2 * x - y, 15 - y, 6 - 3 * x + y]; }
  function n5U(p) {
    const p2 = p * p, p3 = p2 * p, p4 = p2 * p2;
    return -(15 * root5 * p4 + 20 * p3 + 18 * root5 * p2 - 28 * p - 9 * root5) /
      (3 * root5 * p4 + 12 * p3 + 18 * root5 * p2 - 20 * p + 3 * root5);
  }
  function n5Rho(p) { return Math.max(0, 2 * p * (5 * p * p + 6 * root5 * p - 11) / (3 * root5 * p * p + 2 * p + 3 * root5)); }
  function n5Phi(p, r) {
    const p2 = p * p;
    return 2 * (25 * p2 * p2 + 10 * p2 * r + 10 * p2 - 24 * root5 * p * r + 19 * r * r + 58 * r + 9) /
      (5 * (1 + p2 + r) ** 2);
  }
  function n5Lower(x) {
    if (!nearRange(x, 0, 5)) return null;
    x = clamp(x, 0, 5);
    let c, y, parameter;
    if (x < xc) {
      const target = Math.sqrt(5 * x);
      let lo = pc, hi = p0;
      for (let j = 0; j < 64; j++) {
        const mid = (lo + hi) / 2;
        if (n5U(mid) > target) lo = mid; else hi = mid;
      }
      const p = (lo + hi) / 2, r = n5Rho(p), norm = Math.sqrt(1 + p * p + r);
      c = [[-p / norm, 0], [0, 0], [Math.sqrt(r) / norm, 0], [0, 0], [1 / norm, 0], [0, 0]];
      y = n5Phi(p, r);
      parameter = { p, rho: r };
    } else if (x <= 9 / 5) {
      const U = Math.sqrt(5 * x), lambda = Math.sqrt(Math.max(0, (3 - U) / (5 + U)));
      const norm = Math.sqrt(1 + lambda * lambda);
      c = [[0, 0], [1 / norm, 0], [0, 0], [0, 0], [0, 0], [lambda / norm, 0]];
      y = (3 * x - 2 * U + 15) / 4;
      parameter = { lambda };
    } else {
      const sin2Theta = Math.sqrt(Math.max(0, 5 * (5 - x) / 16));
      c = [[Math.sqrt(Math.max(0, 1 - sin2Theta)), 0], [Math.sqrt(sin2Theta), 0], [0, 0], [0, 0], [0, 0], [0, 0]];
      y = 2 * x;
      parameter = { theta: Math.asin(Math.sqrt(sin2Theta)) };
    }
    return record(c, n5Sectors(x, y), 'Five-qubit lower-boundary state', parameter,
      'Dedicated completed five-qubit range proof, explicit equality families');
  }
  function n5Upper(x) {
    if (!nearRange(x, 0, 5)) return null;
    x = clamp(x, 0, 5);
    const spin = Math.sqrt(x / 5);
    const c = [[Math.sqrt((1 + spin) / 2), 0], [0, 0], [0, 0], [0, 0], [0, 0], [Math.sqrt((1 - spin) / 2), 0]];
    return record(c, n5Sectors(x, 10), '|GHZ(5, φ)⟩', { S1: x }, 'Generalized GHZ equality family');
  }
  function n5ZeroY(t) {
    const t2 = t * t;
    return (t * (t * (177 * t2 + 36 * root5 * t - 410) - 60 * root5) + 225) / (20 * (t2 - 1) ** 2);
  }
  function n5Zero(y) {
    if (!nearRange(y, 18 / 5, 10)) return null;
    y = clamp(y, 18 / 5, 10);
    const a = y <= 15 / 4 ? -3 / root5 : Math.sqrt(5 / 3);
    const b = y <= 15 / 4 ? -Math.sqrt((345 + 60 * Math.sqrt(15)) / 289) : root5;
    let lo = 0, hi = 1;
    for (let j = 0; j < 64; j++) {
      const mid = (lo + hi) / 2, t = a + mid * (b - a);
      if (n5ZeroY(t) < y) lo = mid; else hi = mid;
    }
    const t = a + (lo + hi) * (b - a) / 2;
    const norm = 2 * Math.sqrt(t * t - 1);
    const c = [[1 / norm, 0], [0, 0], [Math.sqrt(Math.max(0, 3 * t * t - 5)) / norm, 0], [0, 0], [t / norm, 0], [0, 0]];
    return record(c, n5Sectors(0, y), 'Five-qubit zero-dipole state', { t }, 'Explicit C₂ zero-dipole boundary family');
  }
  function n5Boundary(x, y, tolerance = 1e-8) {
    if (!nearRange(x, 0, 5) || !Number.isFinite(y)) return null;
    if (Math.abs(x) <= tolerance) return n5Zero(y);
    if (Math.abs(y - 10) <= tolerance) return n5Upper(x);
    const lower = n5Lower(x);
    return lower && Math.abs(y - lower.sectors[1]) <= tolerance ? lower : null;
  }
  global.SectorConstructors = Object.freeze({ n3, n4, n5Lower, n5Upper, n5Zero, n5Boundary,
    n5JunctionX: xc, supportsInterior: Object.freeze({ 3: true, 4: true, 5: false }) });
})(globalThis);
