'use strict';
// Run from anywhere with Node. Geometry and witnesses remain unmodified.
const fs = require('fs');
const path = require('path');
const assert = require('assert/strict');
const {performance} = require('perf_hooks');
const root = path.resolve(__dirname, '..');
const THREE = require(path.join(root, 'vendor/three.min.js'));
const StaticTriangleOccluder = require(path.join(root, 'static-occlusion.js'));
const vec = a => new THREE.Vector3(...a);
const ray = (o,d) => new THREE.Ray(vec(o), vec(d));
const triangle = new Float32Array([-1,-1,0, 1,-1,0, 0,1,0]);
const index = new StaticTriangleOccluder([triangle]);
const cases = [
  ['front', ray([0,0,1],[0,0,-1]), 1.01, true],
  ['back', ray([0,0,-1],[0,0,1]), 1.01, true],
  ['strict cutoff', ray([0,0,1],[0,0,-1]), 1, false],
  ['before cutoff', ray([0,0,1],[0,0,-1]), .99, false],
  ['outside box parallel', ray([2,0,1],[0,0,-1]), 2, false],
  ['on box boundary', ray([0,-1,1],[0,0,-1]), 2, true],
  ['away', ray([0,0,1],[0,0,1]), 2, false],
  ['origin on skin', ray([0,0,0],[0,0,1]), 2, true],
  ['coplanar', ray([0,0,0],[1,0,0]), 2, false],
  ['zero direction', ray([0,0,1],[0,0,0]), 2, false],
  ['nonunit inside', ray([0,0,1],[0,0,-3]), 1.01, true],
  ['nonunit outside', ray([0,0,1],[0,0,-3]), .99, false],
  ['infinite distance', ray([0,0,1],[0,0,-1]), Infinity, true],
  ['negative distance', ray([0,0,1],[0,0,-1]), -1, false],
  ['NaN distance', ray([0,0,1],[0,0,-1]), NaN, false]
];
for (const [name, r, distance, expected] of cases)
  assert.equal(index.occludes(r,distance), expected, name);
assert.equal(new StaticTriangleOccluder([]).occludes(cases[0][1],2), false);
assert.equal(new StaticTriangleOccluder([new Float32Array(9)]).occludes(cases[0][1],2),false);
assert.throws(() => new StaticTriangleOccluder([new Float32Array(8)]), RangeError);
assert.throws(() => new StaticTriangleOccluder([[NaN,0,0,0,0,0,0,0,0]]), RangeError);

function map(v) { return [2*(v[0]+1)/9-1, 2*(v[2]+1)/23-1, 2*(v[1]+1.5)/17.5-1]; }
function arraysFrom(patches) {
  return patches.map(p => {
    const output = new Float32Array(p.triangles.length*9);
    p.triangles.forEach((t,i) => t.forEach((id,j) => output.set(map(p.positions[id]),9*i+3*j)));
    return output;
  });
}
const a = new THREE.Vector3(), b = a.clone(), c = a.clone(), hit = a.clone();
function brute(arrays,r,cutoff) {
  const squared = cutoff*cutoff;
  for (const array of arrays) for (let o=0;o<array.length;o+=9) {
    a.fromArray(array,o); b.fromArray(array,o+3); c.fromArray(array,o+6);
    if (r.intersectTriangle(a,b,c,false,hit) && hit.distanceToSquared(r.origin)<squared) return true;
  }
  return false;
}
let seed = 0x12345678;
function random() { seed = (1664525*seed+1013904223)>>>0; return seed/4294967296; }
const randomRays = [];
for (let i=0;i<250;i++) {
  const origin = new THREE.Vector3(...[0,0,0].map(()=>6*random()-3));
  const target = new THREE.Vector3(...[0,0,0].map(()=>2*random()-1));
  const direction = target.sub(origin).normalize().multiplyScalar(i%7===0?2.4:1);
  randomRays.push([new THREE.Ray(origin,direction), .05+random()*6]);
}
const nodes = JSON.parse(fs.readFileSync(path.join(root,'data/n6-manuscript-nodes.json'),'utf8'));
for (const camera of [[3,2,4],[-3,2,-4],[1,-3,3]]) for (const node of nodes) {
  const target = vec(map(node.point));
  const origin = vec(camera);
  const distance = origin.distanceTo(target);
  randomRays.push([new THREE.Ray(origin,target.sub(origin).normalize()),distance-.025]);
}

const results = {passed:true,unitCases:cases.length,cutoff:'strict, in world units',originHit:'included, matching THREE.Ray.intersectTriangle',comparisons:[],benchmarks:[]};
for (const file of ['n6.json','n6-smooth-geometry.json']) {
  const data = JSON.parse(fs.readFileSync(path.join(root,'data',file),'utf8'));
  const arrays = arraysFrom(data.patches || data.outer.patches);
  const before = arrays.map(array => Buffer.from(array.buffer).toString('base64'));
  const bvh = new StaticTriangleOccluder(arrays);
  const checked = randomRays.map(([r,d],i) => {
    const actual = bvh.occludes(r,d), expected = brute(arrays,r,d);
    assert.equal(actual,expected,`${file}, ray ${i}`);
    return actual;
  });
  arrays.forEach((array,i) => assert.equal(Buffer.from(array.buffer).toString('base64'),before[i],`${file} input immutability`));
  results.comparisons.push({file,triangles:bvh.triangleCount,rays:checked.length,hits:checked.filter(Boolean).length,allMatched:true,inputArraysUnchanged:true});
  // Warm before taking query timings, then use deterministic identical workloads.
  for (let i=0;i<1000;i++) { const [r,d]=randomRays[i%randomRays.length]; bvh.occludes(r,d); }
  const trials=5000;
  let t=performance.now();
  for (let i=0;i<trials;i++) { const [r,d]=randomRays[i%randomRays.length]; bvh.occludes(r,d); }
  const bvhMilliseconds=performance.now()-t;
  t=performance.now();
  for (const [r,d] of randomRays) brute(arrays,r,d);
  const bruteMilliseconds=performance.now()-t;
  results.benchmarks.push({...bvh.stats,file,queries:trials,bvhMilliseconds,meanBvhMicroseconds:1000*bvhMilliseconds/trials,bruteQueries:randomRays.length,bruteMilliseconds,meanBruteMicroseconds:1000*bruteMilliseconds/randomRays.length});
}
fs.writeFileSync(path.join(__dirname,'static-occlusion-regression.json'),JSON.stringify(results,null,2)+'\n');
console.log(JSON.stringify(results,null,2));
