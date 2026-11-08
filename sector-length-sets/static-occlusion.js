/* Static, two-sided triangle occlusion for an orthographic atlas camera.
 * Source position arrays are read only; only the private triangle IDs partition.
 */
(function (root) {
  'use strict';
  const T = root.THREE || (typeof require === 'function' ? require('./vendor/three.min.js') : null);
  if (!T) throw new Error('StaticTriangleOccluder requires Three.js.');

  function entry(node, ox, oy, oz, dx, dy, dz, limit) {
    let lo = 0, hi = limit;
    const origins = [ox, oy, oz], directions = [dx, dy, dz];
    for (let axis = 0; axis < 3; axis++) {
      const direction = directions[axis], origin = origins[axis];
      if (direction === 0) {
        if (origin < node.min[axis] || origin > node.max[axis]) return Infinity;
      } else {
        let a = (node.min[axis] - origin) / direction;
        let b = (node.max[axis] - origin) / direction;
        if (a > b) { const tmp = a; a = b; b = tmp; }
        lo = Math.max(lo, a); hi = Math.min(hi, b);
        if (hi < lo) return Infinity;
      }
    }
    return lo < limit && hi >= 0 ? lo : Infinity;
  }

  class StaticTriangleOccluder {
    constructor(triangleArrays, {leafSize = 12} = {}) {
      const started = typeof performance === 'object' ? performance.now() : Date.now();
      if (!Array.isArray(triangleArrays)) throw new TypeError('Expected an array of triangle position arrays.');
      this.arrays = triangleArrays.filter(a => a && a.length);
      for (const array of this.arrays) if (array.length % 9) throw new RangeError('Each triangle requires nine coordinates.');
      this.triangleCount = this.arrays.reduce((sum, a) => sum + a.length / 9, 0);
      this.leafSize = Math.max(2, Math.floor(leafSize));
      this.ids = Uint32Array.from({length: this.triangleCount}, (_, i) => i);
      this.sources = new Uint32Array(this.triangleCount);
      this.offsets = new Uint32Array(this.triangleCount);
      const bounds = new Float64Array(this.triangleCount * 6);
      const centers = new Float64Array(this.triangleCount * 3);
      let triangle = 0;
      this.arrays.forEach((array, source) => {
        for (let offset = 0; offset < array.length; offset += 9, triangle++) {
          this.sources[triangle] = source; this.offsets[triangle] = offset;
          for (let axis = 0; axis < 3; axis++) {
            const a = array[offset + axis], b = array[offset + axis + 3], c = array[offset + axis + 6];
            if (![a, b, c].every(Number.isFinite)) throw new RangeError('Triangle coordinates must be finite.');
            bounds[6 * triangle + axis] = Math.min(a, b, c);
            bounds[6 * triangle + axis + 3] = Math.max(a, b, c);
            centers[3 * triangle + axis] = (a + b + c) / 3;
          }
        }
      });
      this.nodes = [];
      const build = (start, end, depth) => {
        const node = {min: [Infinity, Infinity, Infinity], max: [-Infinity, -Infinity, -Infinity], start, end, left: -1, right: -1};
        const centerMin = [Infinity, Infinity, Infinity], centerMax = [-Infinity, -Infinity, -Infinity];
        for (let i = start; i < end; i++) {
          const id = this.ids[i];
          for (let axis = 0; axis < 3; axis++) {
            node.min[axis] = Math.min(node.min[axis], bounds[6 * id + axis]);
            node.max[axis] = Math.max(node.max[axis], bounds[6 * id + axis + 3]);
            centerMin[axis] = Math.min(centerMin[axis], centers[3 * id + axis]);
            centerMax[axis] = Math.max(centerMax[axis], centers[3 * id + axis]);
          }
        }
        const index = this.nodes.push(node) - 1;
        if (end - start <= this.leafSize || depth >= 48) return index;
        let axis = 0;
        for (let j = 1; j < 3; j++) if (centerMax[j] - centerMin[j] > centerMax[axis] - centerMin[axis]) axis = j;
        if (centerMax[axis] === centerMin[axis]) return index;
        const split = (centerMin[axis] + centerMax[axis]) / 2;
        let i = start, j = end - 1;
        while (i <= j) {
          if (centers[3 * this.ids[i] + axis] < split) i++;
          else { const tmp = this.ids[i]; this.ids[i] = this.ids[j]; this.ids[j--] = tmp; }
        }
        const middle = i > start && i < end ? i : start + Math.floor((end - start) / 2);
        node.left = build(start, middle, depth + 1); node.right = build(middle, end, depth + 1);
        return index;
      };
      if (this.triangleCount) build(0, this.triangleCount, 0);
      this._a = new T.Vector3(); this._b = new T.Vector3(); this._c = new T.Vector3(); this._hit = new T.Vector3(); this._stack = [];
      this.stats = {triangles: this.triangleCount, nodes: this.nodes.length, leafSize: this.leafSize,
        buildMilliseconds: (typeof performance === 'object' ? performance.now() : Date.now()) - started};
    }

    /** Closest visible skin triangle, including its original source and face.
     * Gray triangles participate in depth testing exactly like violet ones.
     */
    closestHit(ray, maxDistance = Infinity) {
      if (!this.nodes.length || !(maxDistance > 0)) return null;
      const {x: ox, y: oy, z: oz} = ray.origin, {x: dx, y: dy, z: dz} = ray.direction;
      const length = Math.hypot(dx, dy, dz);
      if (!length || ![ox, oy, oz, dx, dy, dz].every(Number.isFinite)) return null;
      let distance = maxDistance, result = null;
      const stack = this._stack; stack.length = 0; stack.push(0);
      while (stack.length) {
        const node = this.nodes[stack.pop()];
        if (entry(node, ox, oy, oz, dx, dy, dz, distance / length) === Infinity) continue;
        if (node.left >= 0) {
          const left = entry(this.nodes[node.left], ox, oy, oz, dx, dy, dz, distance / length);
          const right = entry(this.nodes[node.right], ox, oy, oz, dx, dy, dz, distance / length);
          if (left <= right) { if (right !== Infinity) stack.push(node.right); if (left !== Infinity) stack.push(node.left); }
          else { if (left !== Infinity) stack.push(node.left); if (right !== Infinity) stack.push(node.right); }
          continue;
        }
        for (let i = node.start; i < node.end; i++) {
          const id = this.ids[i], source = this.sources[id], array = this.arrays[source], offset = this.offsets[id];
          this._a.fromArray(array, offset); this._b.fromArray(array, offset + 3); this._c.fromArray(array, offset + 6);
          if (!ray.intersectTriangle(this._a, this._b, this._c, false, this._hit)) continue;
          const found = this._hit.distanceTo(ray.origin);
          if (found < distance) {
            distance = found;
            result = {source, faceIndex: offset / 9, distance, point: this._hit.clone()};
          }
        }
      }
      return result;
    }

    /** A hit at the ray origin is counted, as in THREE.Ray.intersectTriangle.
     * maxDistance is a world-space distance, even for a non-unit ray direction.
     */
    occludes(ray, maxDistance = Infinity) {
      if (!this.nodes.length || !(maxDistance > 0)) return false;
      const {x: ox, y: oy, z: oz} = ray.origin, {x: dx, y: dy, z: dz} = ray.direction;
      const length = Math.hypot(dx, dy, dz);
      if (!length || ![ox, oy, oz, dx, dy, dz].every(Number.isFinite)) return false;
      const limit = maxDistance / length, maxSquared = maxDistance * maxDistance;
      if (entry(this.nodes[0], ox, oy, oz, dx, dy, dz, limit) === Infinity) return false;
      const stack = this._stack; stack.length = 0; stack.push(0);
      while (stack.length) {
        const node = this.nodes[stack.pop()];
        if (node.left >= 0) {
          const left = entry(this.nodes[node.left], ox, oy, oz, dx, dy, dz, limit);
          const right = entry(this.nodes[node.right], ox, oy, oz, dx, dy, dz, limit);
          if (left <= right) { if (right !== Infinity) stack.push(node.right); if (left !== Infinity) stack.push(node.left); }
          else { if (left !== Infinity) stack.push(node.left); if (right !== Infinity) stack.push(node.right); }
          continue;
        }
        for (let i = node.start; i < node.end; i++) {
          const id = this.ids[i], array = this.arrays[this.sources[id]], offset = this.offsets[id];
          this._a.fromArray(array, offset); this._b.fromArray(array, offset + 3); this._c.fromArray(array, offset + 6);
          if (ray.intersectTriangle(this._a, this._b, this._c, false, this._hit) && this._hit.distanceToSquared(ray.origin) < maxSquared) return true;
        }
      }
      return false;
    }
  }

  root.StaticTriangleOccluder = StaticTriangleOccluder;
  if (typeof module === 'object' && module.exports) module.exports = StaticTriangleOccluder;
})(typeof globalThis === 'object' ? globalThis : window);
