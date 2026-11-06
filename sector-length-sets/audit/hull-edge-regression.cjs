/* Regression: draw polyhedral edges, never triangulation diagonals. */
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const assert = require('assert');
const root = path.resolve(__dirname, '..');
const T = require(path.join(root, 'vendor/three.min.js'));
const M = require(path.join(root, 'majorana-math.js'));
const source = fs.readFileSync(path.join(root, 'viewer.js'), 'utf8');
const start = source.indexOf('  function hullTopology(');
const end = source.indexOf('  class MajoranaView', start);
assert(start >= 0 && end > start, 'viewer.js must define hullTopology');
const scope = { T, v3: p => new T.Vector3(...p) };
vm.runInNewContext(source.slice(start, end) + '\nthis.topology = hullTopology;', scope);
const topology = scope.topology;
const ring = n => Array.from({length:n}, (_,i) => [Math.cos(i*2*Math.PI/n),Math.sin(i*2*Math.PI/n),0]);
const base = ring(4).map(p => [p[0]*Math.sqrt(15)/4,p[1]*Math.sqrt(15)/4,-.25]);
const cases = [
  {name:'two distinct stars',points:[[0,0,1],[0,0,-1]],edges:1,triangles:0},
  {name:'coplanar regular pentagon', points:ring(5), edges:5, triangles:3},
  {name:'coplanar regular hexagon', points:ring(6), edges:6, triangles:4},
  {name:'coplanar square', points:ring(4), edges:4, triangles:2},
  {name:'tetrahedron', points:[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(p=>p.map(v=>v/Math.sqrt(3))),edges:6,triangles:4},
  {name:'octahedron',points:[[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],edges:12,triangles:8},
  {name:'square pyramid',points:[...base,[0,0,1]],edges:8,triangles:6},
  {name:'triangular bipyramid',points:[...ring(3),[0,0,1],[0,0,-1]],edges:9,triangles:6}
];
const small = JSON.parse(fs.readFileSync(path.join(root, 'data/small-n.json'),'utf8'));
const actual = small['5'].curves.flatMap(c=>c.points).find(p=>p.id==='n5-top-126');
assert(actual, 'the annotated N5 state must exist');
cases.push({name:'actual n5-top-126',points:M.majoranaStars(actual.coefficients).groups.map(p=>[p.x,p.y,p.z]),edges:5,triangles:3});
const axis = new T.Vector3(.2,.7,-.4).normalize();
const rotations = [0,.371,1.209];
const report = [];
for (const example of cases) for (const angle of rotations) {
  const points = example.points.map(p=>new T.Vector3(...p).applyAxisAngle(axis,angle).toArray());
  const h = topology(points);
  assert.equal(h.edges.length,example.edges, example.name+' edge count');
  assert.equal(h.triangles.length,example.triangles, example.name+' fill triangulation');
  const keys = h.edges.map(e=>[...e].sort((a,b)=>a-b).join(','));
  assert.equal(new Set(keys).size,keys.length,example.name+' duplicate edges');
  for (const [a,b] of h.edges) {
    assert(Number.isInteger(a)&&Number.isInteger(b)&&a!==b&&a>=0&&b>=0&&a<points.length&&b<points.length);
    if (example.name.includes('pentagon')||example.name.includes('hexagon')||example.name.includes('actual n5')) {
      const distance = new T.Vector3(...points[a]).distanceTo(new T.Vector3(...points[b]));
      const lengths=points.flatMap((p,i)=>points.slice(i+1).map(q=>new T.Vector3(...p).distanceTo(new T.Vector3(...q))));
      assert(Math.abs(distance-Math.min(...lengths))<1e-7,'polygon diagonal was drawn');
    }
  }
  report.push({name:example.name,rotation:angle,edges:h.edges.length,triangles:h.triangles.length});
}
const output={pass:true,cases:report.length,results:report};
fs.writeFileSync(path.join(__dirname,'hull-edge-regression.json'),JSON.stringify(output,null,2));
console.log(JSON.stringify(output,null,2));
