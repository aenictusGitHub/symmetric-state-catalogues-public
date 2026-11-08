/* Real N=6 geometry/BVH gestures, verified state callbacks and app wiring. */
'use strict';
const assert = require('node:assert/strict'), fs = require('node:fs');
const path = require('node:path'), vm = require('node:vm');
const root = path.resolve(__dirname, '../sector-length-sets');
const T = require(path.join(root, 'vendor/three.min.js'));
const M = require(path.join(root, 'majorana-math.js'));
const Occluder = require(path.join(root, 'static-occlusion.js'));
for (const file of ['n6-surface-bounds', 'n6-cyclic-constructors', 'n6-filtered-paired-constructors',
  'n6-roof-constructors', 'n6-surface-state']) require(path.join(root, file + '.js'));
const G = JSON.parse(fs.readFileSync(path.join(root, 'data/n6-smooth-geometry.json')));
const map = p => [2 * (p[0] + 1) / 9 - 1, 2 * (p[2] + 1) / 23 - 1, 2 * (p[1] + 1.5) / 17.5 - 1];
const surface = new T.Group();
for (const patch of G.patches) {
  const values = [];
  for (const face of patch.triangles) for (const i of face) values.push(...map(patch.positions[i]));
  const geometry = new T.BufferGeometry(); geometry.setAttribute('position', new T.Float32BufferAttribute(values, 3));
  const mesh = new T.Mesh(geometry, new T.MeshBasicMaterial({side: T.DoubleSide}));
  mesh.userData.patch = patch; surface.add(mesh);
}
surface.updateMatrixWorld(true);
const bvh = new Occluder(surface.children.map(m => m.geometry.getAttribute('position').array));
const frames = new Map(); let nextFrame = 0;
const requestFrame = fn => { const id = ++nextFrame; frames.set(id, fn); return id; };
const cancelFrame = id => frames.delete(id);
function flushFrames() { const pending = [...frames.values()]; frames.clear(); pending.forEach(fn => fn()); }
const context = {window: {THREE: T}, requestAnimationFrame: requestFrame, cancelAnimationFrame: cancelFrame};
vm.runInNewContext(fs.readFileSync(path.join(root, 'viewer.js'), 'utf8'), context);
const View = context.window.SectorView;
assert.equal(typeof View.prototype.traceSurfaceDrag, 'function');
const close = (a, b, tolerance = 1e-9) => assert.ok(Number.isFinite(a) && Number.isFinite(b) &&
  Math.abs(a - b) <= tolerance, `${a} != ${b} (tolerance ${tolerance})`);
const same = (a, b, tolerance = 1e-9) => { assert.equal(a.length, b.length); a.forEach((v, i) => close(v, b[i], tolerance)); };
const distance = (a, b) => Math.hypot(...a.map((v, i) => v - b[i]));
let physicalChecks = 0, maximumResidual = 0, gestureCount = 0, grayBridgeChecks = 0, rearChecks = 0;
function verify(state) {
  assert(state); assert.equal(state.coefficients.length, 7);
  close(state.coefficients.reduce((sum, z) => sum + z[0] ** 2 + z[1] ** 2, 0), 1, 2e-12);
  const actual = M.sectorLengths(state.coefficients);
  actual.forEach((v, i) => close(v, state.sectors[i], 1e-6));
  actual.slice(0, 3).forEach((v, i) => close(v, state.targetPoint[i], 1e-7));
  maximumResidual = Math.max(maximumResidual, state.reconstructionResidual);
  if (state.surfaceAdjusted) {
    assert(state.clickedPoint && state.surfaceAdjustment);
    state.surfaceAdjustment.forEach((v, i) => close(v, state.targetPoint[i] - state.clickedPoint[i], 1e-12));
  }
  physicalChecks++; return actual;
}
function makeView(yaw, pitch, zoom = 1) {
  const view = Object.create(View.prototype), captures = new Set();
  Object.assign(view, {n: 6, dimension: 3, data: {}, width: 640, height: 640,
    surface, occlusionIndex: bvh, raycaster: new T.Raycaster(),
    camera: new T.OrthographicCamera(-1.8, 1.8, 1.8, -1.8, .01, 100),
    target: new T.Vector3(), yaw, pitch, zoom, rotatable: true, panAllowed: true,
    selectedWorld: null, drag: null, moved: false, callbacks: [], rejected: 0,
    namedWorld: [], namedProjections: [], projected: [], grayPoints: [],
    canvas: {style: {}, getBoundingClientRect: () => ({left: 11, top: 19, width: 640, height: 640}),
      setPointerCapture(id) { captures.add(id); }, hasPointerCapture(id) { return captures.has(id); },
      releasePointerCapture(id) { captures.delete(id); }},
    captures, render() { this.renderCount = (this.renderCount || 0) + 1; },
    clickCount: 0, onClick() { this.clickCount++; }});
  view.camera.zoom = zoom; view.camera.updateProjectionMatrix(); view.updateCamera();
  view.surfaceDragCallback = function (point, meta) {
    this.callbacks.push({point: Array.from(point), meta: {...meta}});
    if (meta.phase !== 'move') return true;
    if (this.rejectNext) { this.rejectNext = false; this.rejected++; return false; }
    const state = N6SurfaceState.construct(point, meta.surface, meta.kind, meta.anchorPoint);
    if (!state) { this.rejected++; return false; }
    verify(state); this.state = state; this.selectedWorld = map(M.sectorLengths(state.coefficients)); return true;
  };
  return view;
}
function event(view, screen, {offset = [0, 0], pointerId = 73, button = 0, shiftKey = false} = {}) {
  const rect = view.canvas.getBoundingClientRect();
  return {clientX: screen[0] + rect.left + offset[0], clientY: screen[1] + rect.top + offset[1], pointerId, button, shiftKey};
}
function setStateAt(view, x, y) {
  const hit = view.attainedSurfaceAt(x, y); if (!hit) return null;
  const state = N6SurfaceState.construct(hit.point, hit.surface, hit.kind, hit.anchorPoint); if (!state) return null;
  view.state = state; view.selectedWorld = map(M.sectorLengths(state.coefficients));
  const p = view.project(view.selectedWorld);
  if (!view.selectedPointAt(p.x, p.y)) return null;
  return {state, hit, screen: [p.x, p.y], world: view.selectedWorld.slice()};
}
function localFixture(view) {
  for (let y = 80; y <= 580; y += 20) for (let x = 80; x <= 580; x += 20) {
    const start = setStateAt(view, x, y); if (!start) continue;
    for (const direction of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
      const goal = start.screen.map((v, i) => v + 18 * direction[i]);
      const trace = view.traceSurfaceDrag({acceptedScreen: start.screen, acceptedWorld: start.world}, ...goal);
      if (trace && distance(trace.screen, goal) < 1e-8 && N6SurfaceState.construct(trace.hit.point, trace.hit.surface, trace.hit.kind, trace.hit.anchorPoint))
        return {...start, direction};
    }
  }
  throw new Error('No visible local violet drag fixture');
}
const movedScreen = (fixture, amount) => fixture.screen.map((v, i) => v + amount * fixture.direction[i]);
const presets = [[-1.04, .22], [-1.95, .50], [2.10, .377]];
const fixtures = [];
for (const [yaw, pitch] of presets) {
  const view = makeView(yaw, pitch), f = localFixture(view); fixtures.push({yaw, pitch, f});
  const orientation = [view.yaw, view.pitch, view.camera.zoom, ...view.target.toArray()], offset = [3, -2];
  const down = event(view, f.screen, {offset});
  assert(view.beginPointerGesture(down)); assert.equal(view.drag.point.n, 6); assert.equal(view.drag.pan, false);
  view.updatePointerGesture({...down, clientX: down.clientX + 3}); assert.equal(view.callbacks.length, 0);
  const a = event(view, movedScreen(f, 8), {offset}), b = event(view, movedScreen(f, 12), {offset});
  assert.equal(view.updatePointerGesture({...a, pointerId: 99}), false); assert.equal(view.callbacks.length, 0);
  view.updatePointerGesture(a); view.updatePointerGesture(b);
  assert.equal(frames.size, 1); assert.deepEqual(view.callbacks.map(c => c.meta.phase), ['start']);
  same(view.callbacks[0].point, f.state.sectors.slice(0, 3), 1e-7);
  flushFrames(); assert.deepEqual(view.callbacks.map(c => c.meta.phase), ['start', 'move']);
  assert.equal(view.rejected, 0); same(view.drag.point.acceptedScreen, movedScreen(f, 12), 1e-8);
  assert(view.callbacks.at(-1).meta.surface && view.callbacks.at(-1).meta.anchorPoint);
  view.updatePointerGesture(event(view, movedScreen(f, 16), {offset})); assert.equal(frames.size, 1);
  assert.equal(view.endPointerGesture({...b, pointerId: 99}), false);
  view.endPointerGesture(event(view, movedScreen(f, 18), {offset}));
  assert.deepEqual(view.callbacks.map(c => c.meta.phase), ['start', 'move', 'move', 'end']);
  assert.equal(frames.size, 0); assert.equal(view.drag, null); assert.equal(view.captures.size, 0);
  same([view.yaw, view.pitch, view.camera.zoom, ...view.target.toArray()], orientation, 1e-13);
  assert.equal(view.clickCount, 0); gestureCount++;

  // Application rejection does not advance accepted geometry or selection.
  const rejection = makeView(yaw, pitch), rf = localFixture(rejection); rejection.beginPointerGesture(event(rejection, rf.screen));
  rejection.updatePointerGesture(event(rejection, movedScreen(rf, 8))); flushFrames();
  const accepted = {...rejection.drag.point, acceptedScreen: rejection.drag.point.acceptedScreen.slice(),
    lastCoordinates: rejection.drag.point.lastCoordinates.slice()}, oldState = rejection.state;
  rejection.rejectNext = true; rejection.updatePointerGesture(event(rejection, movedScreen(rf, 14))); flushFrames();
  assert.equal(rejection.rejected, 1); assert.equal(rejection.state, oldState);
  same(rejection.drag.point.acceptedScreen, accepted.acceptedScreen);
  same(rejection.drag.point.lastCoordinates, accepted.lastCoordinates);
  rejection.updatePointerGesture(event(rejection, movedScreen(rf, 14))); flushFrames(); assert.notEqual(rejection.state, oldState);
  const last = rejection.drag.point.lastCoordinates.slice();
  rejection.updatePointerGesture(event(rejection, movedScreen(rf, 18))); assert.equal(frames.size, 1);
  assert.equal(rejection.cancelPointerGesture({pointerId: 99}), false); assert.equal(frames.size, 1);
  rejection.cancelPointerGesture({pointerId: 73}); assert.equal(frames.size, 0);
  assert.equal(rejection.callbacks.at(-1).meta.phase, 'cancel'); same(rejection.callbacks.at(-1).point, last);
  const count = rejection.callbacks.length; flushFrames(); assert.equal(rejection.callbacks.length, count); gestureCount++;
}

// Find real visible-violet endpoints separated by a broad gray/missing stripe.
// A single fast pointer move must stop before the stripe, not teleport to the
// endpoint. The test uses the nearest-hit BVH on the unchanged manuscript mesh.
for (const [yaw, pitch] of presets) {
  const view = makeView(yaw, pitch), hits = [];
  for (let y = 80; y <= 580; y += 20) for (let x = 80; x <= 580; x += 20)
    if (view.attainedSurfaceAt(x, y)) hits.push([x, y]);
  let fixture = null;
  for (let i = 0; i < hits.length && !fixture; i += Math.max(1, Math.floor(hits.length / 25))) {
    const start = setStateAt(view, ...hits[i]); if (!start) continue;
    for (let j = 0; j < hits.length; j += Math.max(1, Math.floor(hits.length / 25))) {
      const goal = hits[j], length = distance(start.screen, goal); if (length < 60) continue;
      let firstBlocked = null, blockedLength = 0;
      for (let d = 2; d < length; d += 2) {
        const p = start.screen.map((v, k) => v + d * (goal[k] - v) / length);
        if (!view.attainedSurfaceAt(...p)) { if (firstBlocked === null) firstBlocked = d; blockedLength += 2; }
      }
      if (firstBlocked > 6 && blockedLength > 12) { fixture = {...start, goal, firstBlocked, length}; break; }
    }
  }
  assert(fixture, 'A real gray-separated violet fixture exists in every perspective');
  view.state = fixture.state; view.selectedWorld = fixture.world.slice();
  view.beginPointerGesture(event(view, fixture.screen)); view.updatePointerGesture(event(view, fixture.goal)); flushFrames();
  assert(view.drag.point); assert.equal(view.rejected, 0);
  const acceptedDistance = distance(view.drag.point.acceptedScreen, fixture.screen);
  assert(acceptedDistance <= fixture.firstBlocked + 2.01, 'The trace cannot skip its first gray sample');
  assert(distance(view.drag.point.acceptedScreen, fixture.goal) > 8, 'The far violet endpoint is not reached through gray');
  view.updatePointerGesture(event(view, fixture.goal)); flushFrames();
  assert(distance(view.drag.point.acceptedScreen, fixture.goal) > 8);
  view.endPointerGesture(event(view, fixture.goal)); assert.equal(view.callbacks.at(-1).meta.phase, 'end');
  grayBridgeChecks++; gestureCount++;
}

// Rear violet markers cannot drag through a nearer gray face. Decode an
// actual rear triangle hit and reconstruct its pure state before testing.
function decode(hit) {
  const mesh = hit.object, patch = mesh.userData.patch, face = hit.faceIndex, attr = mesh.geometry.getAttribute('position');
  const a = new T.Vector3().fromBufferAttribute(attr, face * 3), b = new T.Vector3().fromBufferAttribute(attr, face * 3 + 1), c = new T.Vector3().fromBufferAttribute(attr, face * 3 + 2);
  const w = T.Triangle.getBarycoord(hit.point, a, b, c, new T.Vector3()); if (!w) return null;
  const weights = [w.x, w.y, w.z], ids = patch.triangles[face];
  return {point: [0, 1, 2].map(k => ids.reduce((sum, id, i) => sum + weights[i] * patch.positions[id][k], 0)),
    anchorPoint: [0, 1, 2].map(k => ids.reduce((sum, id) => sum + patch.positions[id][k], 0) / 3),
    surface: patch.surface, kind: patch.kind};
}
for (const [yaw, pitch] of presets) {
  const view = makeView(yaw, pitch); let found = false;
  outer: for (let y = 70; y < 600; y += 18) for (let x = 70; x < 600; x += 18) {
    view.raycaster.setFromCamera(new T.Vector2(2 * x / 640 - 1, 1 - 2 * y / 640), view.camera);
    const hits = view.raycaster.intersectObjects(surface.children, false);
    if (!hits.length || hits[0].object.userData.patch.attained[hits[0].faceIndex]) continue;
    for (const rear of hits.slice(1)) {
      if (!rear.object.userData.patch.attained[rear.faceIndex]) continue;
      const meta = decode(rear), state = meta && N6SurfaceState.construct(meta.point, meta.surface, meta.kind, meta.anchorPoint);
      if (!state) continue;
      view.selectedWorld = map(M.sectorLengths(state.coefficients)); const p = view.project(view.selectedWorld);
      if (view.attainedSurfaceAt(p.x, p.y)) continue;
      assert.equal(view.selectedPointAt(p.x, p.y), false);
      view.beginPointerGesture(event(view, [p.x, p.y])); assert.equal(view.drag.point, null);
      const before = view.yaw; view.updatePointerGesture(event(view, [p.x + 15, p.y + 8]));
      assert.notEqual(view.yaw, before); view.endPointerGesture(event(view, [p.x + 15, p.y + 8]));
      assert.equal(view.callbacks.length, 0); found = true; rearChecks++; break outer;
    }
  }
  // This particular front/back color combination need not occur in every
  // perspective; the analytic F_(1/3) state below tests both rear views.
}
assert(rearChecks > 0, 'The real mesh includes a front-gray/rear-violet occlusion fixture');
const filteredThird = JSON.parse(fs.readFileSync(path.join(root, 'data/n6-boundary-samples.json')))
  .samples.find(s => s.id === 'n6-fo-c3-3');
assert(filteredThird, 'The exact filtered-octahedral one-third sample is present');
for (let i = 0; i < presets.length; i++) {
  const view = makeView(...presets[i]); view.selectedWorld = map(M.sectorLengths(filteredThird.coefficients));
  const p = view.project(view.selectedWorld);
  assert.equal(view.selectedPointAt(p.x, p.y), i === 2,
    'The actual F_(1/3) state is rear in Main/Alternative and front in Rear');
  if (i < 2) { view.beginPointerGesture(event(view, [p.x, p.y])); assert.equal(view.drag.point, null); view.cancelPointerGesture(); }
  rearChecks++;
}

// The drag-start depth allowance must not grow with zoom or inherit the
// larger drawing-only occlusion epsilon. A distinct point .01 world units
// behind a violet front ray is ineligible even though the marker renderer
// might still paint it as front. This is a guard test, not a state claim.
for (const [yaw, pitch] of presets) {
  const view = makeView(yaw, pitch), f = localFixture(view);
  view.raycaster.setFromCamera(new T.Vector2(2 * f.screen[0] / 640 - 1, 1 - 2 * f.screen[1] / 640), view.camera);
  const front = bvh.closestHit(view.raycaster.ray); assert(front);
  const behind = view.raycaster.ray.at(front.distance + .01, new T.Vector3());
  view.selectedWorld = behind.toArray(); const p = view.project(view.selectedWorld);
  assert(view.attainedSurfaceAt(p.x, p.y), 'The near face is still violet');
  assert.equal(view.selectedPointAt(p.x, p.y), false, 'A thin rear sheet does not pass the strict start-depth test');
}

// Actual app callback + real selectState/showState: live state/constellation
// updates precede release, failed inverse leaves selection unchanged, and
// release/cancel performs exactly one dropdown/URL commit.
const appSource = fs.readFileSync(path.join(root, 'app.js'), 'utf8');
const displayStart = appSource.indexOf('  function showState('), displayEnd = appSource.indexOf('  function updateFamily(', displayStart);
const callbackStart = appSource.indexOf('sector.surfaceDragCallback='), callbackEnd = appSource.indexOf('sector.pointDragCallback=', callbackStart);
assert(displayStart >= 0 && displayEnd > displayStart && callbackStart >= 0 && callbackEnd > callbackStart);
let appMoves = 0, appCommits = 0;
for (const phase of ['end', 'cancel']) {
  const {yaw, pitch} = fixtures[0], view = makeView(yaw, pitch), f = localFixture(view), moves = [];
  for (const amount of [8, 14]) { const hit = view.attainedSurfaceAt(...movedScreen(f, amount)); assert(hit); moves.push(hit); }
  const calls = {select: [], sphere: [], live: 0, full: 0, history: 0, options: 0, notify: 0, sectors: 0}, elements = new Map();
  const el = id => { if (!elements.has(id)) elements.set(id, {dataset: {}, innerHTML: '', textContent: '', hidden: false,
    children: [], setAttribute() {}, querySelectorAll() { return this.children; }, append(o) { this.children.push(o); calls.options++; }}); return elements.get(id); };
  const app = {n: 6, pointDragging: false, previewing: true, previewTimer: 123, selected: {...f.state, id: 'stale'}, current: f.state,
    document: {querySelector: el, createElement: () => ({dataset: {}, remove() {}})}, $: el,
    SectorMath: M, clearTimeout() {}, notify() { calls.notify++; }, fmt: x => String(x),
    title: s => s.label || 'State', complex: z => z.join(','), esc: String, paintMath(e, s) { e.innerHTML = s; },
    AtlasMath: {state: (_id, label) => label, html: s => s}, showExpression() {}, updateFamily() {},
    renderSectorValues() { calls.sectors++; }, history: {replaceState() { calls.history++; }},
    stateData(s) { return {...s, N: 6, coefficients: M.canonicalState(s.coefficients), sectors: M.sectorLengths(s.coefficients)}; },
    sphere: {setState(c, options) { calls.sphere.push(options); return M.majoranaStars(c); }},
    sector: {setSelected() {}}, window: {StatePicker: {updateLive() { calls.live++; }, update() { calls.full++; }},
      N6SurfaceState: {construct(...args) { return app.reject ? null : N6SurfaceState.construct(...args); }}}};
  vm.runInNewContext(appSource.slice(displayStart, displayEnd), app); const realSelect = app.selectState;
  app.selectState = (...args) => { calls.select.push(args); return realSelect(...args); };
  vm.runInNewContext(appSource.slice(callbackStart, callbackEnd), app);
  const visible = app.current; assert.equal(app.sector.surfaceDragCallback(f.state.sectors.slice(0, 3), {...f.hit, phase: 'start'}), true);
  assert.equal(app.selected, visible); assert.equal(app.pointDragging, true);
  for (const hit of moves) {
    assert.equal(app.sector.surfaceDragCallback(hit.point, {...hit, phase: 'move'}), true);
    assert.equal(calls.select.at(-1)[1], false); assert.equal(calls.select.at(-1)[2].live, true);
    assert.equal(calls.sphere.at(-1).preserveOrientation, true); appMoves++;
  }
  assert.equal(calls.live, 2); assert.equal(calls.sphere.length, 2); assert.equal(calls.sectors, 2);
  assert.equal(calls.history, 0); assert.equal(calls.options, 0); assert.equal(calls.notify, 0);
  const accepted = app.selected, selectedCount = calls.select.length; app.reject = true;
  assert.equal(app.sector.surfaceDragCallback(moves[0].point, {...moves[0], phase: 'move'}), false);
  assert.equal(app.selected, accepted); assert.equal(calls.select.length, selectedCount); app.reject = false;
  assert.equal(app.sector.surfaceDragCallback(accepted.sectors.slice(0, 3), {...moves[1], phase}), true);
  assert.equal(app.pointDragging, false); assert.equal(calls.history, 1); assert.equal(calls.options, 1); assert.equal(calls.full, 1);
  assert.equal(calls.select.at(-1)[0], accepted); assert.equal(calls.select.at(-1)[1], true);
  assert.equal(calls.sphere.at(-1).preserveOrientation, true);
  app.sector.surfaceDragCallback(accepted.sectors.slice(0, 3), {...moves[1], phase}); assert.equal(calls.history, 1); appCommits++;
}
assert.equal(frames.size, 0);
console.log(JSON.stringify({status: 'PASS', perspectives: presets.length, gestureCount, physicalChecks,
  maximumResidual, grayBridgeChecks, rearChecks, appMoves, appCommits,
  callbackRejectionPreservesState: 'PASS', frameCoalescingAndCancellation: 'PASS'}, null, 2));
