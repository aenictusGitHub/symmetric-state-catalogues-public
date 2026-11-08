/* Selected-point gestures on the actual N=3,4,5 camera/projection code. */
'use strict';
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '../sector-length-sets');
const T = require(path.join(root, 'vendor/three.min.js'));
const M = require(path.join(root, 'majorana-math.js'));
require(path.join(root, 'data/small-n-constructors.js'));
const C = globalThis.SectorConstructors;
const frames = new Map();
let nextFrame = 0, delivered = 0, physicalChecks = 0, starChecks = 0;
let gestures = 0, maximumSectorError = 0;
const requestFrame = fn => { const id = ++nextFrame; frames.set(id, fn); return id; };
const cancelFrame = id => frames.delete(id);
function flushFrames() {
  const pending = [...frames.values()]; frames.clear();
  for (const fn of pending) fn();
}
const context = {
  window: {THREE: T, SectorConstructors: C},
  requestAnimationFrame: requestFrame, cancelAnimationFrame: cancelFrame,
  AtlasMath: {state: (_id, label) => label || ''}
};
vm.runInNewContext(fs.readFileSync(path.join(root, 'viewer.js'), 'utf8'), context);
const View = context.window.SectorView;
for (const method of ['beginPointerGesture', 'updatePointerGesture', 'endPointerGesture',
  'cancelPointerGesture', 'beginPointDrag', 'regionCoordinates'])
  assert.equal(typeof View.prototype[method], 'function', `Missing production method ${method}`);
const close = (a, b, tol = 2e-10) => assert.ok(Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) <= tol,
  `${a} differs from ${b} by more than ${tol}`);
function equalCoordinates(actual, expected, tolerance = 2e-10) {
  assert.equal(actual.length, expected.length);
  actual.forEach((v, i) => close(v, expected[i], tolerance));
}
function witness(n, coordinates) {
  const [x, y] = coordinates, r = n === 3 ? C.n3(x) : n === 4 ? C.n4(x, y) : C.n5(x, y);
  assert(r, `Delivered unattainable N=${n} coordinates ${coordinates}`);
  close(r.coefficients.reduce((sum, z) => sum + z[0] ** 2 + z[1] ** 2, 0), 1, 2e-12);
  const actual = M.sectorLengths(r.coefficients);
  actual.forEach((v, i) => {
    maximumSectorError = Math.max(maximumSectorError, Math.abs(v - r.sectors[i]));
    close(v, r.sectors[i], 2e-9);
  });
  coordinates.forEach((v, i) => close(v, actual[i], 2e-9));
  if (physicalChecks % 23 === 0) {
    const stars = M.majoranaStars(r.coefficients);
    assert.equal(stars.stars.length, n);
    assert(stars.reliable && stars.fidelity > 1 - 1e-8);
    starChecks++;
  }
  physicalChecks++;
  return r;
}
function makeView(n, coordinates, tx = 0, ty = 0, zoom = 1) {
  const view = Object.create(View.prototype), captures = new Set();
  Object.assign(view, {n, dimension: n === 3 ? 1 : n === 6 ? 3 : 2, data: {},
    ranges: n === 3 ? [[0, 3], [0, 1]] : n === 4 ? [[0, 4], [2, 6]] : [[0, 5], [3, 10]],
    width: 800, height: 600, raycaster: new T.Raycaster(),
    camera: new T.OrthographicCamera(-2, 2, 1.5, -1.5, .01, 100),
    target: new T.Vector3(tx, ty, 0), yaw: 0, pitch: 0, zoom,
    panAllowed: true, rotatable: n === 6, drag: null, moved: false,
    render() { this.renderCount = (this.renderCount || 0) + 1; },
    canvas: {style: {}, getBoundingClientRect: () => ({left: 17, top: 23, width: 800, height: 600}),
      setPointerCapture(id) { captures.add(id); },
      hasPointerCapture(id) { return captures.has(id); },
      releasePointerCapture(id) { captures.delete(id); }},
    namedWorld: [], namedProjections: [], projected: [], grayPoints: [],
    callbacks: [], clickCount: 0, captures});
  view.camera.zoom = zoom; view.camera.updateProjectionMatrix(); view.updateCamera();
  view.selectedWorld = view.map(coordinates);
  view.onClick = () => view.clickCount++;
  view.pointDragCallback = (coords, {phase}) => {
    const values = Array.from(coords); view.callbacks.push({phase, coordinates: values}); delivered++;
    if (n <= 5) witness(n, values);
    if (phase === 'move') view.selectedWorld = view.map(values);
  };
  return view;
}
function eventAt(view, coordinates, {offset = [0, 0], pointerId = 41, button = 0, shiftKey = false} = {}) {
  const p = view.project(view.map(coordinates)), r = view.canvas.getBoundingClientRect();
  return {clientX: p.x + r.left + offset[0], clientY: p.y + r.top + offset[1], pointerId, button, shiftKey};
}
function cameraSnapshot(view) { return [view.target.x, view.target.y, view.target.z, view.yaw, view.pitch, view.camera.zoom]; }
const initial = n => n === 3 ? [1.2] : n === 4 ? [1, 4] : [1, 6];
const middle = n => n === 3 ? [1.8] : n === 4 ? [2, 5] : [2.7, 7];

// Actual camera inversion, pointer-to-center offsets, threshold, ownership,
// animation-frame coalescing and final synchronous release at several zooms.
for (const n of [3, 4, 5]) for (const [tx, ty, zoom] of [[0, 0, .55], [.23, -.37, 1], [-.4, .3, 4], [.01, -.02, 20]]) {
  const view = makeView(n, initial(n), tx, ty, zoom), before = cameraSnapshot(view), offset = [5, -3];
  const down = eventAt(view, initial(n), {offset});
  assert(view.beginPointerGesture(down)); assert(view.drag.point && !view.drag.pan);
  assert(view.captures.has(41)); assert.equal(view.canvas.style.cursor, 'grabbing');
  assert.equal(view.beginPointerGesture({...down, pointerId: 99}), false, 'A second pointer cannot replace the active drag');
  view.updatePointerGesture({...down, clientX: down.clientX + 3});
  assert.equal(view.callbacks.length, 0, 'Sub-threshold motion keeps ordinary click semantics');
  const goal = eventAt(view, middle(n), {offset});
  assert.equal(view.updatePointerGesture({...goal, pointerId: 99}), false);
  assert.equal(view.callbacks.length, 0, 'Another pointer cannot move the selected point');
  view.updatePointerGesture(goal);
  assert.deepEqual(view.callbacks.map(c => c.phase), ['start']);
  equalCoordinates(view.callbacks[0].coordinates, initial(n));
  const last = n === 3 ? [2.1] : n === 4 ? [2.4, 5.4] : [3.2, 8.1];
  view.updatePointerGesture(eventAt(view, middle(n), {offset}));
  view.updatePointerGesture(eventAt(view, last, {offset}));
  assert.equal(frames.size, 1, 'Many pointer moves share one scheduled frame');
  flushFrames();
  assert.deepEqual(view.callbacks.map(c => c.phase), ['start', 'move']);
  equalCoordinates(view.callbacks.at(-1).coordinates, last);
  equalCoordinates(cameraSnapshot(view), before, 1e-13);
  const final = n === 3 ? [.9] : n === 4 ? [1.4, 4.2] : [1.4, 5.2];
  view.updatePointerGesture(eventAt(view, last, {offset}));
  assert.equal(frames.size, 1);
  assert.equal(view.endPointerGesture({...goal, pointerId: 99}), false);
  assert.equal(frames.size, 1, 'Wrong-pointer release leaves the scheduled frame intact');
  view.endPointerGesture(eventAt(view, final, {offset}));
  assert.deepEqual(view.callbacks.map(c => c.phase), ['start', 'move', 'move', 'end']);
  equalCoordinates(view.callbacks.at(-2).coordinates, final);
  equalCoordinates(view.callbacks.at(-1).coordinates, final);
  assert.equal(frames.size, 0); assert.equal(view.drag, null); assert.equal(view.captures.size, 0);
  assert.equal(view.clickCount, 0, 'Dragging does not trigger a second point click');
  const count = view.callbacks.length; flushFrames(); assert.equal(view.callbacks.length, count);
  equalCoordinates(cameraSnapshot(view), before, 1e-13); gestures++;
}

// Dragging a named selected state is supported, including W's shared sector
// coordinate. It does not require a previously created interactive marker.
for (const [n, coordinate] of [[3, [1 / 3]], [4, [1, 3]], [5, [1.8, 3.6]]]) {
  const view = makeView(n, coordinate); view.selectedState = {id: `n${n}-w`, sectors: coordinate};
  const down = eventAt(view, coordinate); view.beginPointerGesture(down);
  assert(view.drag.point); view.updatePointerGesture({...down, clientX: down.clientX + 12}); flushFrames();
  assert.equal(view.callbacks[0].phase, 'start'); equalCoordinates(view.callbacks[0].coordinates, coordinate);
  assert.equal(view.callbacks[1].phase, 'move'); view.endPointerGesture({...down, clientX: down.clientX + 12}); gestures++;
}

// Pointer capture permits motion beyond the canvas; domain clamping still
// supplies a pure-state witness, including N=5's genuinely curved lower edge.
for (const n of [3, 4, 5]) {
  const view = makeView(n, initial(n), -.2, .1, 2), down = eventAt(view, initial(n));
  view.beginPointerGesture(down);
  const targets = n === 3 ? [[[-10], [0]], [[10], [3]]] : n === 4 ?
    [[[-2, -10], [0, 2]], [[7, -10], [4, 6]], [[1, 20], [1, 6]], [[2, -2], [2, 4]]] :
    [[[-2, -10], [0, 3.6]], [[7, -10], [5, 10]], [[.01, 0], [.01, C.n5Lower(.01).sectors[1]]],
      [[1, 20], [1, 10]], [[.5, C.n5Lower(.5).sectors[1] - .2], [.5, C.n5Lower(.5).sectors[1]]]];
  for (const [outside, expected] of targets) {
    const event = eventAt(view, outside); if (n === 3) event.clientY += 10000;
    view.updatePointerGesture(event); flushFrames();
    equalCoordinates(view.callbacks.at(-1).coordinates, expected);
    assert.equal(view.callbacks.at(-1).phase, 'move');
  }
  view.endPointerGesture(eventAt(view, targets.at(-1)[0])); gestures++;
}

// Cancel/lost capture discards pending motion, commits only the last delivered
// point, and never schedules a stale callback after the gesture has ended.
for (const n of [3, 4, 5]) {
  const view = makeView(n, initial(n)), down = eventAt(view, initial(n)); view.beginPointerGesture(down);
  view.updatePointerGesture(eventAt(view, middle(n))); flushFrames();
  const deliveredPoint = view.callbacks.at(-1).coordinates;
  view.updatePointerGesture(eventAt(view, n === 3 ? [2.5] : n === 4 ? [3, 5.5] : [4, 9]));
  assert.equal(frames.size, 1);
  assert.equal(view.cancelPointerGesture({pointerId: 99}), false); assert.equal(frames.size, 1);
  assert(view.cancelPointerGesture({pointerId: 41}));
  assert.equal(view.callbacks.at(-1).phase, 'cancel'); equalCoordinates(view.callbacks.at(-1).coordinates, deliveredPoint);
  assert.equal(frames.size, 0); assert.equal(view.drag, null); assert.equal(view.captures.size, 0);
  const count = view.callbacks.length; flushFrames(); assert.equal(view.callbacks.length, count);
  assert.equal(view.endPointerGesture(down), false); assert.equal(view.cancelPointerGesture(), false); gestures++;
}

// Plain click and canceled unstarted gestures retain their previous behavior.
for (const n of [3, 4, 5]) {
  const view = makeView(n, initial(n)), down = eventAt(view, initial(n)); view.beginPointerGesture(down); view.endPointerGesture(down);
  assert.equal(view.callbacks.length, 0); assert.equal(view.clickCount, 1);
  view.beginPointerGesture(down); view.cancelPointerGesture(); assert.equal(view.callbacks.length, 0); assert.equal(view.clickCount, 1);
}

// Empty-space dragging, Shift and either nonprimary button still pan the view.
for (const n of [3, 4, 5]) for (const variant of ['empty', 'shift', 'middle', 'right']) {
  const view = makeView(n, initial(n)), before = cameraSnapshot(view);
  const options = variant === 'shift' ? {shiftKey: true} : variant === 'middle' ? {button: 1} : variant === 'right' ? {button: 2} : {};
  let down = eventAt(view, initial(n), options);
  if (variant === 'empty') down = {...down, clientX: down.clientX + 60, clientY: down.clientY + 60};
  view.beginPointerGesture(down); assert(view.drag.pan && !view.drag.point);
  const up = {...down, clientX: down.clientX + 40, clientY: down.clientY + 22}; view.updatePointerGesture(up); view.endPointerGesture(up);
  assert.equal(view.callbacks.length, 0); assert.equal(view.clickCount, 0);
  assert.notDeepEqual(cameraSnapshot(view), before, `${variant} must still pan N=${n}`); gestures++;
}

// N=6 remains a surface/rotation viewer: selected-node gestures are not
// silently converted into N=3--5 coordinate dragging.
const six = makeView(6, [0, 5, 0]), down = eventAt(six, [0, 5, 0]), beforeYaw = six.yaw;
six.beginPointerGesture(down); assert.equal(six.drag.point, null); assert.equal(six.drag.pan, false);
six.updatePointerGesture({...down, clientX: down.clientX + 30, clientY: down.clientY + 10});
assert.notEqual(six.yaw, beforeYaw); six.endPointerGesture({...down, clientX: down.clientX + 30, clientY: down.clientY + 10});
assert.equal(six.callbacks.length, 0); assert.equal(frames.size, 0);

// Execute the current application's actual gesture callback and actual
// selectState/showState functions. Only external DOM/rendering boundaries are
// mocked, so a wiring change that updates only on release cannot pass.
const appSource = fs.readFileSync(path.join(root, 'app.js'), 'utf8');
const buildersStart = appSource.indexOf('  function threeQubitState(');
const buildersEnd = appSource.indexOf('  function notify(', buildersStart);
const displayStart = appSource.indexOf('  function showState(');
const displayEnd = appSource.indexOf('  function updateFamily(', displayStart);
const callbackStart = appSource.indexOf('sector.pointDragCallback=');
const callbackEnd = appSource.indexOf('sector.coordinateCallback=', callbackStart);
assert(buildersStart >= 0 && buildersEnd > buildersStart && displayStart >= 0 &&
  displayEnd > displayStart && callbackStart >= 0 && callbackEnd > callbackStart,
  'The production application integration points must be present');
let appLiveMoves = 0, appCommits = 0, appSphereUpdates = 0;
function makeApp(n) {
  const calls = {select: [], sphere: [], sectors: [], expressions: [],
    pickerLive: [], pickerFull: [], history: [], notifications: [], families: [],
    insertedOptions: 0, removedOptions: 0, clearedTimers: []};
  const elements = new Map();
  function element(id) {
    if (!elements.has(id)) elements.set(id, {id, hidden: false, dataset: {}, innerHTML: '', textContent: '', value: '',
      setAttribute(name, value) { this[name] = value; }, children: [],
      querySelectorAll() { return this.children.filter(o => o.dataset?.interactivePoint); },
      append(option) { this.children.push(option); option.parent = this; calls.insertedOptions++; }});
    return elements.get(id);
  }
  const document = {querySelector: element, createElement() { return {dataset: {}, remove() {
    if (this.parent) this.parent.children.splice(this.parent.children.indexOf(this), 1);
    calls.removedOptions++;
  }}; }};
  const app = {n, selected: null, current: null, previewing: true, pointDragging: false,
    previewTimer: 123, interactiveState: null, SectorConstructors: C, SectorMath: M,
    document, $: element, window: {StatePicker: {
      updateLive(state) { calls.pickerLive.push(state); }, update(state) { calls.pickerFull.push(state); }
    }}, AtlasMath: {state: (_id, label) => label || '', html: text => text},
    fmt: (value, digits = 6) => Number(value.toFixed(digits)).toString(),
    title: state => state.label || 'State', complex: z => z.join(','), esc: String,
    paintMath(el, text) { el.innerHTML = text; },
    showExpression(state, coefficients, options) { calls.expressions.push({state, coefficients, options}); },
    renderSectorValues(state, sectors) { calls.sectors.push({state, sectors}); },
    updateFamily(state) { calls.families.push(state); },
    clearTimeout(id) { calls.clearedTimers.push(id); },
    notify(message) { calls.notifications.push(message); },
    history: {replaceState(...args) { calls.history.push(args); }},
    stateData(state) { return {...state, N: n, coefficients: M.canonicalState(state.coefficients), sectors: M.sectorLengths(state.coefficients)}; },
    sphere: {setState(coefficients, options) {
      const data = M.majoranaStars(coefficients); assert(data.reliable);
      calls.sphere.push({coefficients, options}); return data;
    }},
    sector: {setSelected(state) { this.selectedState = state; }}, calls};
  vm.runInNewContext(appSource.slice(buildersStart, buildersEnd) + '\n' +
    appSource.slice(displayStart, displayEnd), app);
  const realSelect = app.selectState;
  app.selectState = function (...args) { calls.select.push(args); return realSelect(...args); };
  vm.runInNewContext(appSource.slice(callbackStart, callbackEnd), app);
  return app;
}
function appState(app, coordinates) {
  return app.n === 3 ? app.threeQubitState(coordinates[0]) :
    app.n === 4 ? app.triangleState(...coordinates) : app.fiveQubitState(...coordinates);
}
for (const n of [3, 4, 5]) for (const ending of ['end', 'cancel']) {
  const app = makeApp(n), visible = appState(app, middle(n)), stale = appState(app, initial(n));
  app.current = visible; app.selected = stale;
  app.sector.pointDragCallback(middle(n), {phase: 'start'});
  assert.equal(app.selected, visible, 'Drag start pins the currently displayed point, not a stale selection');
  assert.equal(app.previewing, false); assert.equal(app.pointDragging, true);
  assert.equal(app.calls.select.length, 0, 'Starting a drag does not unnecessarily rebuild the state');
  const positions = n === 3 ? [[1.6], [2.2]] : n === 4 ? [[1.6, 4.8], [2.2, 5.5]] : [[1.6, 6.8], [2.2, 7.9]];
  for (const coordinates of positions) {
    app.sector.pointDragCallback(coordinates, {phase: 'move'});
    const invocation = app.calls.select.at(-1);
    assert.equal(invocation[1], false); assert.equal(invocation[2].live, true);
    assert.equal(app.selected.id, app.current.id); assert.equal(app.sector.selectedState.id, app.current.id);
    coordinates.forEach((v, i) => close(v, app.current.sectors[i], 2e-9));
    assert.equal(app.calls.sphere.at(-1).options.preserveOrientation, true);
    assert.equal(app.calls.expressions.at(-1).options.live, true);
    assert.equal(app.calls.sectors.at(-1).state.id, app.selected.id);
    appLiveMoves++;
  }
  assert.equal(app.calls.history.length, 0, 'Live motion never writes the shared URL');
  assert.equal(app.calls.insertedOptions, 0, 'Live motion never rebuilds native options');
  assert.equal(app.calls.pickerFull.length, 0); assert.equal(app.calls.pickerLive.length, positions.length);
  assert.equal(app.calls.notifications.length, 0, 'Live motion never creates notification spam');
  assert.equal(app.calls.families.length, 0);
  assert.equal(app.calls.sphere.length, positions.length, 'The Majorana sphere updates before release');
  const finalState = app.selected;
  app.sector.pointDragCallback(positions.at(-1), {phase: ending});
  assert.equal(app.calls.select.length, positions.length + 1);
  assert.equal(app.calls.select.at(-1)[0], finalState); assert.equal(app.calls.select.at(-1)[1], true);
  assert.equal(app.pointDragging, false); assert.equal(app.calls.history.length, 1);
  assert.equal(app.calls.insertedOptions, 1); assert.equal(app.calls.pickerFull.length, 1);
  assert.equal(app.calls.notifications.length, 1); assert.equal(app.calls.families.length, 1);
  assert.equal(app.calls.sphere.at(-1).options.preserveOrientation, true,
    'The final render retains the same sphere orientation as live motion');
  app.sector.pointDragCallback(positions.at(-1), {phase: ending});
  assert.equal(app.calls.history.length, 1, 'A repeated termination cannot commit twice');
  appCommits++; appSphereUpdates += app.calls.sphere.length;
}

// Run the real StatePicker.updateLive method after initialization: it may
// close the menu and repaint its trigger, but native options and the menu's
// mathematical entries must remain untouched throughout high-frequency moves.
const pickerSource = fs.readFileSync(path.join(root, 'state-picker.js'), 'utf8');
function functionSource(source, name) {
  const start = source.indexOf(`  function ${name}(`), end = source.indexOf('\n  function ', start + 1);
  assert(start >= 0 && end > start, `Missing StatePicker.${name}`);
  return source.slice(start, end);
}
const nativeOptions = Object.freeze([{value: 'named-a'}, {value: 'named-b'}]);
const nativeSelect = new Proxy({options: nativeOptions, value: 'named-a'}, {
  set() { throw new Error('Live picker updates must not mutate the native select'); }
});
const menuChildren = Object.freeze([{math: 'state-a'}, {math: 'state-b'}]);
const forbiddenMutation = () => { throw new Error('Live picker updates must not rebuild the menu'); };
const pickerContext = {select: nativeSelect, menu: {hidden: false, children: menuChildren,
  replaceChildren: forbiddenMutation, append: forbiddenMutation}, trigger: {attributes: {},
  setAttribute(name, value) { this.attributes[name] = value; }, removeAttribute(name) { delete this.attributes[name]; }},
  display: {innerHTML: ''}, search: 'old', currentState: null, emptyText: 'Empty',
  init: () => true, typeset: (_id, label) => `rendered:${label}`,
  rebuild: forbiddenMutation, paint: forbiddenMutation};
vm.runInNewContext(functionSource(pickerSource, 'close') + '\n' + functionSource(pickerSource, 'updateLive'), pickerContext);
for (let i = 0; i < 100; i++) {
  const state = {id: `moving-${i}`, mathLabel: `$|\\psi(${i})\\rangle$`};
  pickerContext.updateLive(state);
  assert.equal(pickerContext.currentState, state);
  assert.equal(pickerContext.display.innerHTML, `rendered:${state.mathLabel}`);
  assert.equal(pickerContext.select.options, nativeOptions); assert.equal(pickerContext.select.value, 'named-a');
  assert.equal(pickerContext.menu.children, menuChildren); assert.equal(pickerContext.menu.hidden, true);
}

console.log(JSON.stringify({status: 'PASS', gestures, deliveredCallbacks: delivered,
  physicalChecks, starChecks, maximumSectorError, animationFrameCoalescing: 'PASS',
  pointerOwnershipAndCancellation: 'PASS', panAndN6Preserved: 'PASS',
  appLiveMoves, appCommits, appSphereUpdates, livePickerRepaints: 100,
  continuousAppUpdatesAndSingleCommit: 'PASS'}, null, 2));
