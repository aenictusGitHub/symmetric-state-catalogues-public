#!/usr/bin/env node
'use strict';

// Check the constructed pure state, independently of its returned sector tuple.
// This does not test a numerical optimizer or infer exclusions from failed fits.
const assert = require('node:assert/strict');
const path = require('node:path');
const root = path.resolve(__dirname, '..');
const M = require(path.join(root, 'sector-length-sets/majorana-math.js'));
require(path.join(root, 'sector-length-sets/data/small-n-constructors.js'));
const C = globalThis.SectorConstructors;
assert.equal(typeof C?.n5, 'function', 'SectorConstructors.n5 must be implemented');
assert.equal(C.supportsInterior[5], true, 'N=5 interior inspection must be enabled');

const sectorTolerance = 2e-9;
let cases = 0, majoranaCases = 0;
let maxSectorError = 0, maxRawSectorError = 0, maxSpinError = 0;
let maxNormError = 0, minimumFidelity = 1;
const branches = new Set();

function rawSectorLengths(c) {
  const sectors = [1];
  for (let r = 1; r <= 5; r++) {
    let value = 2 ** r * M.reducedPurity(c, r);
    for (let k = 0; k < r; k++)
      value -= M.binomial(r, k) * sectors[k] / M.binomial(5, k);
    sectors.push(value * M.binomial(5, r));
  }
  return sectors.slice(1); // No small-sector rounding.
}

function spinSector(c) {
  let mz = 0, plusReal = 0, plusImag = 0, norm = 0;
  for (let k = 0; k <= 5; k++) {
    const [a, b] = c[k], p = a * a + b * b;
    norm += p;
    mz += (2.5 - k) * p;
    if (k) {
      const [u, v] = c[k - 1], factor = Math.sqrt(k * (6 - k));
      plusReal += factor * (u * a + v * b);
      plusImag += factor * (u * b - v * a);
    }
  }
  return .8 * (mz * mz + plusReal * plusReal + plusImag * plusImag) / norm ** 2;
}

function check(x, y, label, majorana = false) {
  const record = C.n5(x, y);
  assert(record, `${label}: rejected an attainable coordinate (${x}, ${y})`);
  const c = record.coefficients;
  assert.equal(c?.length, 6, `${label}: six Dicke amplitudes are required`);
  for (const z of c) {
    assert(Array.isArray(z) && z.length === 2 && z.every(Number.isFinite),
      `${label}: invalid complex amplitude`);
  }
  const norm = c.reduce((sum, z) => sum + z[0] ** 2 + z[1] ** 2, 0);
  const normError = Math.abs(norm - 1);
  assert(normError < 3e-13, `${label}: normalization error ${normError}`);
  maxNormError = Math.max(maxNormError, normError);
  const expected = [x, y, 10 + 2 * x - y, 15 - y, 6 - 3 * x + y];
  const sectors = M.sectorLengths(c), raw = rawSectorLengths(c);
  assert.equal(record.sectors?.length, 5, `${label}: missing displayed sector tuple`);
  for (let k = 0; k < 5; k++) {
    const error = Math.abs(sectors[k] - expected[k]);
    const rawError = Math.abs(raw[k] - expected[k]);
    assert(error < sectorTolerance,
      `${label}: S${k + 1} error ${error}; target (${x}, ${y})`);
    assert(rawError < sectorTolerance,
      `${label}: raw S${k + 1} error ${rawError}; target (${x}, ${y})`);
    assert(Math.abs(record.sectors[k] - expected[k]) < sectorTolerance,
      `${label}: displayed S${k + 1} disagrees with the requested tuple`);
    maxSectorError = Math.max(maxSectorError, error);
    maxRawSectorError = Math.max(maxRawSectorError, rawError);
  }
  const spinError = Math.abs(spinSector(c) - x);
  assert(spinError < sectorTolerance, `${label}: independent dipole error ${spinError}`);
  maxSpinError = Math.max(maxSpinError, spinError);
  assert(['parity', 'kramers', 'cat'].includes(record.parameter?.branch),
    `${label}: missing constructive-branch metadata`);
  branches.add(record.parameter.branch);
  if (majorana) {
    const data = M.majoranaStars(c);
    assert.equal(data.stars.length, 5, `${label}: incomplete Majorana constellation`);
    assert(data.reliable, `${label}: unreliable Majorana reconstruction`);
    const fidelity = M.fidelity(c, M.stateFromStars(data.stars));
    assert(fidelity > 1 - 1e-8, `${label}: Majorana fidelity ${fidelity}`);
    minimumFidelity = Math.min(minimumFidelity, fidelity, data.fidelity);
    majoranaCases++;
  }
  cases++;
}

// Entire closed region: both edges and intermediate vertical fractions.
for (let i = 0; i <= 70; i++) {
  const x = 5 * i / 70, lower = C.n5Lower(x).sectors[1];
  for (let j = 0; j <= 24; j++)
    check(x, lower + (10 - lower) * j / 24, `grid ${i}/${j}`);
}

let seed = 20261108;
function random() {
  seed = (1664525 * seed + 1013904223) >>> 0;
  return seed / 2 ** 32;
}
for (let i = 0; i < 2500; i++) {
  const x = 5 * random(), lower = C.n5Lower(x).sectors[1];
  check(x, lower + (10 - lower) * random(), `random ${i}`, i % 125 === 0);
}

// Numerical seams, the true lower-bound junctions and both near-corner limits.
const xs = [0, 1e-14, 1e-10, 1e-6,
  C.n5JunctionX - 1e-12, C.n5JunctionX, C.n5JunctionX + 1e-12,
  .4 - 1e-12, .4, .4 + 1e-12,
  1.8 - 1e-12, 1.8, 1.8 + 1e-12,
  5 - 1e-10, 5 - 1e-14, 5];
for (const x of xs) {
  const lower = C.n5Lower(x).sectors[1];
  const ys = [lower, lower + 1e-10, (lower + 10) / 2,
    3.6 - 1e-12, 3.6, 3.6 + 1e-12,
    8.85 - 1e-12, 8.85, 8.85 + 1e-12, 10 - 1e-10, 10];
  for (const y of ys)
    if (y >= lower && y <= 10) check(x, y, `seam (${x},${y})`);
}
for (const [x, y] of [[0, 3.6], [0, 8.85], [0, 10], [.02, 3.6],
  [.4, 3.6], [.4, 8.85], [.4, 10], [1.8, 3.6], [1.8, 7],
  [4, 9], [5, 10]]) check(x, y, `Majorana (${x},${y})`, true);

// Rejection must follow the proved domain, including nonnumeric inputs.
for (const bad of [NaN, Infinity, -Infinity, undefined, null, '1', [], {}]) {
  assert.equal(C.n5(bad, 8), null, `invalid x ${String(bad)}`);
  assert.equal(C.n5(1, bad), null, `invalid y ${String(bad)}`);
}
for (const [x, y] of [[-1e-8, 8], [5 + 1e-8, 10], [1, 10 + 1e-8],
  [0, 3.6 - 1e-8], [1.8, 3.6 - 1e-8], [5, 10 - 1e-8]])
  assert.equal(C.n5(x, y), null, `outside region (${x},${y})`);
for (const x of [.01, .2, .4, 1, 2, 4])
  assert.equal(C.n5(x, C.n5Lower(x).sectors[1] - 1e-8), null,
    `below the curved boundary at x=${x}`);

// The declared 1e-10 domain tolerance admits and clamps tiny roundoff excursions.
for (const [x, y, targetX, targetY] of [
  [-5e-11, 3.6, 0, 3.6], [5 + 5e-11, 10, 5, 10],
  [1, 10 + 5e-11, 1, 10], [0, 3.6 - 5e-11, 0, 3.6]]) {
  const record = C.n5(x, y);
  assert(record, `roundoff coordinate (${x},${y}) was rejected`);
  const actual = M.sectorLengths(record.coefficients);
  assert(Math.abs(actual[0] - targetX) < sectorTolerance);
  assert(Math.abs(actual[1] - targetY) < sectorTolerance);
}
assert.deepEqual([...branches].sort(), ['cat', 'kramers', 'parity']);
console.log(JSON.stringify({status: 'PASS', cases, majoranaCases,
  branches: [...branches].sort(), maxNormError, maxSectorError,
  maxRawSectorError, maxSpinError, minimumFidelity}, null, 2));
