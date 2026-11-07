/* Regression checks for stable rigid axis guides. Run with Node.js. */
'use strict';
const assert = require('node:assert/strict');
const {axisGuide, surfaceDepth, visibleIntervals} = require('../subspace-render-geometry.js');
const near = (a, b, tolerance = 1e-12) => assert.ok(Math.abs(a - b) < tolerance, `${a} versus ${b}`);
const dot = (a, b) => a.reduce((sum, x, i) => sum + x * b[i], 0);
const matvec = (M, p) => M.map(row => dot(row, p));
const matmul = (A, B) => A.map(row => B[0].map((_, j) => row.reduce((sum, x, i) => sum + x * B[i][j], 0)));
const rz = theta => [[Math.cos(theta), -Math.sin(theta), 0], [Math.sin(theta), Math.cos(theta), 0], [0, 0, 1]];
const rx = theta => [[1, 0, 0], [0, Math.cos(theta), -Math.sin(theta)], [0, Math.sin(theta), Math.cos(theta)]];
const projection = M => p => {
  const q = matvec(M, p);
  return {x: 320 + 180 * q[0], y: 240 - 180 * q[2], z: q[1]};
};

// A horizon crossing previously moved the label from +v to -v, by a full
// sphere diameter. Its material endpoint must remain fixed under any camera.
let previous = null;
for (let i = -200; i <= 200; i++) {
  const angle = i * 1e-4, project = projection(rz(angle));
  const guide = axisGuide([1, 0, 0], project);
  assert.deepEqual(guide.label.point, [1.17, 0, 0]);
  assert.deepEqual(guide.halves[0].b, [1.13, 0, 0]);
  assert.deepEqual(guide.halves[1].b, [-1.13, 0, 0]);
  const anchor = project(guide.label.point);
  if (previous) {
    assert.ok(Math.hypot(anchor.x - previous.anchor.x, anchor.y - previous.anchor.y) < 0.03);
    assert.ok(Math.abs(guide.label.opacity - previous.guide.label.opacity) < 0.001);
    assert.ok(Math.abs(guide.halves[0].solidOpacity - previous.guide.halves[0].solidOpacity) < 0.001);
  }
  for (const half of guide.halves) {
    assert.ok(half.solidOpacity >= 0 && half.solidOpacity <= 0.95);
    assert.ok(half.dashedOpacity >= 0 && half.dashedOpacity <= 0.45);
  }
  previous = {anchor, guide};
}

// The physical unoriented axis and its label have one fixed canonical sign,
// regardless of an equivalent negative or rescaled input representation.
const identityProject = projection(rz(0));
assert.deepEqual(axisGuide([-3, 0, 0], identityProject).label.point, [1.17, 0, 0]);
assert.deepEqual(axisGuide({direction: [2, 0, 0]}, identityProject).label.point, [1.17, 0, 0]);

// Face-on axes lose projected length: hide their labels before they pile up
// at the sphere center. Rear labels retain a faint, rather than reversed, end.
near(axisGuide([0, 1, 0], identityProject).label.opacity, 0);
near(axisGuide([1, 0, 0], projection(rz(-Math.PI / 6))).label.opacity, 0.5);
near(axisGuide([1, 0, 0], projection(rz(Math.PI / 6))).label.opacity, 1);
const left = axisGuide([1, 0, 0], projection(rz(Math.PI)));
assert.ok(left.label.radialOffset[0] < 0);
assert.equal(left.label.textAnchor, 'middle');

// A guide and a star on its axis remain collinear after an independent
// alignment plus camera rotation; label and star direction cannot drift.
const axis = [0.6, 0, 0.8];
for (let i = 0; i < 71; i++) {
  const alignment = rz(0.37), camera = matmul(rx(i * 0.13), rz(-i * 0.19));
  const M = matmul(camera, alignment), project = projection(M);
  const guide = axisGuide(axis, project), center = project([0, 0, 0]);
  const star = project(axis), end = project(guide.halves[0].b), label = project(guide.label.point);
  near(end.x - center.x, 1.13 * (star.x - center.x), 1e-10);
  near(end.y - center.y, 1.13 * (star.y - center.y), 1e-10);
  near(label.x - center.x, 1.17 * (star.x - center.x), 1e-10);
  near(label.y - center.y, 1.17 * (star.y - center.y), 1e-10);
  near(dot(matvec(M, axis), matvec(M, [0, 1, 0])), dot(axis, [0, 1, 0]));
}
assert.throws(() => axisGuide([0, 0, 0], identityProject), RangeError);
assert.throws(() => axisGuide([NaN, 1, 0], identityProject), TypeError);
// Opaque radial lobes may occlude one another even on the front hemisphere.
// A screen depth buffer must select the nearest surface, not simply z>=0.
const triangle = [{x:1,y:1,z:0.4},{x:7,y:1,z:0.4},{x:1,y:7,z:0.4}];
const triangleDepth = surfaceDepth(triangle, [[0,1,2]], 9, 9);
near(triangleDepth.data[2*9+2], 0.4, 1e-7);
assert.equal(triangleDepth.data[7*9+7], -Infinity);
const behind = visibleIntervals({x:0,y:2.5,z:0.2},{x:8,y:2.5,z:0.2},triangleDepth,{tolerance:0});
const contains = (intervals,t) => intervals.some(([lo,hi]) => lo <= t && t <= hi);
assert.equal(contains(behind,0.375),false); // Hidden by the front lobe.
assert.equal(contains(behind,0.02),true); // Visible outside the lobe.
assert.equal(contains(behind,0.95),true);
assert.deepEqual(visibleIntervals({x:0,y:2.5,z:0.6},{x:8,y:2.5,z:0.6},triangleDepth,{tolerance:0}),[[0,1]]);

// Interpolated depth, overlapping triangles, and shared quad diagonals.
const slope = surfaceDepth([{x:1,y:1,z:0},{x:7,y:1,z:1},{x:1,y:7,z:0}],[[0,1,2]],9,9);
near(slope.data[2*9+2],0.25,1e-7);
const square = [{x:0,y:0,z:-0.4},{x:8,y:0,z:-0.4},{x:8,y:8,z:-0.4},{x:0,y:8,z:-0.4}];
const squareDepth = surfaceDepth(square,[[0,1,2,3]],8,8);
assert.ok([...squareDepth.data].every(z => Math.abs(z+0.4)<1e-7));
assert.deepEqual(visibleIntervals({x:0,y:3,z:-0.3},{x:8,y:3,z:-0.3},squareDepth,{tolerance:0}),[[0,1]]);
assert.deepEqual(visibleIntervals({x:0,y:3,z:-0.5},{x:8,y:3,z:-0.5},squareDepth,{tolerance:0}),[]);
const nearSquare = square.map(p=>({...p,z:0.6}));
const overlap = surfaceDepth([...square,...nearSquare],[[0,1,2,3],[4,5,6,7]],8,8);
assert.ok([...overlap.data].every(z=>Math.abs(z-0.6)<1e-7));
assert.deepEqual(visibleIntervals({x:-4,y:4,z:0.8},{x:12,y:4,z:0.8},overlap),[[0.25,0.75]]);
assert.deepEqual(visibleIntervals({x:-4,y:-2,z:0.8},{x:12,y:-2,z:0.8},overlap),[]);
// A line along the viewing direction can cross the nearest surface as well.
const depthCross = visibleIntervals({x:4,y:4,z:0},{x:4,y:4,z:1},overlap,{tolerance:0});
assert.equal(depthCross.length,1);
near(depthCross[0][0],0.6,0.002);
near(depthCross[0][1],1);
console.log('PASS: stable axes, horizon/face-on continuity, rigid alignment, opaque-lobe visibility, depth interpolation and viewport clipping.');
