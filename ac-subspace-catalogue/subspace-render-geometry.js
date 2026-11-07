/* Pure display geometry shared by the canvas and SVG whole-subspace views.
 * The camera looks along +y; project(p).z is the signed camera depth of p.
 * Physical axis endpoints are fixed. Camera motion changes only visibility.
 */
(function (root) {
  'use strict';

  const clamp = (x, lo, hi) => Math.max(lo, Math.min(hi, x));
  const scale = (v, factor) => v.map(x => x === 0 ? 0 : x * factor);

  function smoothstep(lo, hi, x) {
    const t = clamp((x - lo) / (hi - lo), 0, 1);
    return t * t * (3 - 2 * t);
  }

  function canonicalDirection(axis) {
    const source = Array.isArray(axis) ? axis : axis?.axis || axis?.direction;
    if (!Array.isArray(source) || source.length !== 3 || !source.every(Number.isFinite)) {
      throw new TypeError('An axis requires three finite coordinates.');
    }
    const length = Math.hypot(...source);
    if (length < 1e-14) throw new RangeError('An axis direction must be nonzero.');
    let direction = scale(source, 1 / length);
    let pivot = 0;
    for (let i = 1; i < 3; i++) {
      if (Math.abs(direction[i]) > Math.abs(direction[pivot])) pivot = i;
    }
    if (direction[pivot] < 0) direction = scale(direction, -1);
    return direction;
  }

  function axisGuide(axis, project, options = {}) {
    const direction = canonicalDirection(axis);
    const origin = [0, 0, 0];
    const extent = options.extent ?? 1.13;
    const labelExtent = options.labelExtent ?? 1.17;
    const horizon = options.horizon ?? 0.10;
    const offset = options.labelOffset ?? 9;
    if (![extent, labelExtent, horizon, offset].every(Number.isFinite) || !(extent > 0 && labelExtent > 0 && horizon > 0)) {
      throw new RangeError('Axis display extents and horizon width must be positive.');
    }
    const center = project(origin);
    const endpoint = project(direction);
    // Depth remains in unit-sphere coordinates; CSS scale only affects x/y.
    const depth = clamp(endpoint.z - center.z, -1, 1);
    const frontness = smoothstep(-horizon, horizon, depth);
    const halves = [1, -1].map(sign => {
      const front = sign > 0 ? frontness : 1 - frontness;
      return {
        a: origin.slice(),
        b: scale(direction, sign * extent),
        solidOpacity: 0.95 * front,
        dashedOpacity: 0.45 * (1 - front)
      };
    });

    // Do not choose the nearest endpoint on each frame: that teleports labels
    // between antipodes when the axis crosses the camera horizon.
    const point = scale(direction, labelExtent);
    const anchor = project(point);
    const dx = anchor.x - center.x, dy = anchor.y - center.y;
    const radialLength = Math.hypot(dx, dy);
    const radialOffset = radialLength > 1e-9 ? [offset * dx / radialLength, offset * dy / radialLength] : [0, 0];
    const projectedFraction = Math.sqrt(Math.max(0, 1 - depth * depth));
    const faceOn = smoothstep(0.08, 0.22, projectedFraction);
    return {
      direction,
      halves,
      label: {
        point,
        opacity: (0.5 + 0.5 * frontness) * faceOn,
        radialOffset,
        textAnchor: 'middle'
      }
    };
  }

  // Rasterize the nearest opaque surface in the same orthographic camera used
  // by the renderer. CSS pixel centers, rather than device pixels, define it.
  function surfaceDepth(projectedVertices, faces, width, height) {
    if (![width, height].every(x => Number.isFinite(x) && x > 0)) {
      throw new RangeError('A depth buffer requires positive finite dimensions.');
    }
    width = Math.max(1, Math.round(width));
    height = Math.max(1, Math.round(height));
    const data = new Float32Array(width * height);
    data.fill(-Infinity);
    const vertex = index => projectedVertices[index]?.p || projectedVertices[index];
    function triangle(a, b, c) {
      if (![a, b, c].every(p => p && [p.x, p.y, p.z].every(Number.isFinite))) return;
      const denominator = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
      if (Math.abs(denominator) < 1e-12) return;
      const x0 = Math.max(0, Math.ceil(Math.min(a.x, b.x, c.x) - 0.5));
      const x1 = Math.min(width - 1, Math.floor(Math.max(a.x, b.x, c.x) - 0.5));
      const y0 = Math.max(0, Math.ceil(Math.min(a.y, b.y, c.y) - 0.5));
      const y1 = Math.min(height - 1, Math.floor(Math.max(a.y, b.y, c.y) - 0.5));
      for (let y = y0; y <= y1; y++) {
        for (let x = x0; x <= x1; x++) {
          const px = x + 0.5, py = y + 0.5;
          const wa = ((b.y - c.y) * (px - c.x) + (c.x - b.x) * (py - c.y)) / denominator;
          const wb = ((c.y - a.y) * (px - c.x) + (a.x - c.x) * (py - c.y)) / denominator;
          const wc = 1 - wa - wb;
          if (wa < -1e-7 || wb < -1e-7 || wc < -1e-7) continue;
          const depth = wa * a.z + wb * b.z + wc * c.z;
          const index = y * width + x;
          if (depth > data[index]) data[index] = depth;
        }
      }
    }
    for (const face of faces) {
      const ids = Array.isArray(face) ? face : face?.ids;
      if (!ids || ids.length < 3) continue;
      const a = vertex(ids[0]);
      for (let i = 1; i + 1 < ids.length; i++) triangle(a, vertex(ids[i]), vertex(ids[i + 1]));
    }
    return {width, height, data};
  }

  // Return ranges in the original segment parameter t. They can be projected
  // into canvas/SVG paths identically; no 3D hemisphere assumption is made.
  function visibleIntervals(a, b, buffer, options = {}) {
    const {width, height, data} = buffer || {};
    if (!(width > 0 && height > 0 && data && data.length === width * height)) {
      throw new TypeError('A valid surface depth buffer is required.');
    }
    if (![a, b].every(p => p && [p.x, p.y, p.z].every(Number.isFinite))) return [];
    const tolerance = options.tolerance ?? 0.003;
    const pixelStep = options.pixelStep ?? 0.75;
    if (!(Number.isFinite(tolerance) && tolerance >= 0 && Number.isFinite(pixelStep) && pixelStep > 0)) {
      throw new RangeError('Visibility tolerance and pixel step must be finite and nonnegative/positive.');
    }
    const dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
    let lower = 0, upper = 1;
    // Liang–Barsky clipping retains t in the input segment, including when
    // its endpoints are far outside the current viewport.
    for (const [p, q] of [[-dx, a.x], [dx, width - a.x], [-dy, a.y], [dy, height - a.y]]) {
      if (Math.abs(p) < 1e-14) { if (q < 0) return []; continue; }
      const r = q / p;
      if (p < 0) lower = Math.max(lower, r); else upper = Math.min(upper, r);
      if (lower > upper) return [];
    }
    function visible(t) {
      const x = clamp(Math.floor(a.x + dx * t), 0, width - 1);
      const y = clamp(Math.floor(a.y + dy * t), 0, height - 1);
      const surface = data[y * width + x];
      return a.z + dz * t + tolerance >= surface;
    }
    const steps = Math.max(1, Math.ceil(Math.hypot(dx, dy) * (upper - lower) / pixelStep));
    const intervals = [];
    let previousT = lower, previousVisible = visible(lower), start = previousVisible ? lower : null;
    for (let i = 1; i <= steps; i++) {
      const t = lower + (upper - lower) * i / steps, isVisible = visible(t);
      if (isVisible !== previousVisible) {
        let lo = previousT, hi = t;
        for (let iteration = 0; iteration < 10; iteration++) {
          const mid = (lo + hi) / 2;
          if (visible(mid) === previousVisible) lo = mid; else hi = mid;
        }
        const boundary = (lo + hi) / 2;
        if (previousVisible) { intervals.push([start, boundary]); start = null; }
        else start = boundary;
      }
      previousT = t;
      previousVisible = isVisible;
    }
    if (previousVisible) intervals.push([start, upper]);
    return intervals.filter(([lo, hi]) => hi - lo > 1e-12);
  }

  const api = { axisGuide, surfaceDepth, visibleIntervals, version: '1.1.0' };
  root.SubspaceRenderGeometry = api;
  if (typeof module === 'object' && module.exports) module.exports = api;
})(typeof window === 'object' ? window : globalThis);
