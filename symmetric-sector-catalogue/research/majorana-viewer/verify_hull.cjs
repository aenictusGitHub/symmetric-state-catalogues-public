const fs=require('fs'),vm=require('vm'),path=require('path'),assert=require('assert/strict');
const root=path.resolve(__dirname,'../..'),win={};vm.runInNewContext(fs.readFileSync(path.join(root,'majorana-hull.js'),'utf8'),{window:win});
const build=win.MajoranaHull.build,sub=(a,b)=>a.map((x,i)=>x-b[i]),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]],dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const point=v=>{const n=Math.hypot(...v);return {x:v[0]/n,y:v[1]/n,z:v[2]/n,multiplicity:1};};
function check(points,counts){
 const before=JSON.stringify(points),h=build(points);assert.equal(JSON.stringify(points),before);
 const v=h.vertices.map(p=>[p.x,p.y,p.z]);
 if(counts)assert.deepEqual([h.vertices.length,h.edges.length,h.faces.length],counts);
 if(h.dimension===3){
  const used=new Set(h.faces.flat());assert.equal(used.size,h.vertices.length,'all unit-sphere positions are extreme');assert.equal(used.size-h.edges.length+h.faces.length,2,'Euler characteristic');
  const incidence=new Map();
  for(const f of h.faces){
   assert.ok(f.length>=3);const n=cross(sub(v[f[1]],v[f[0]]),sub(v[f[2]],v[f[0]])),length=Math.hypot(...n);assert.ok(length>1e-10);
   for(const p of v)assert.ok(dot(n,sub(p,v[f[0]]))/length<2e-8,'every face supports all vertices');
   for(let i=0;i<f.length;i++){const key=[f[i],f[(i+1)%f.length]].sort((a,b)=>a-b).join(',');incidence.set(key,(incidence.get(key)||0)+1);}
  }
  assert.ok([...incidence.values()].every(n=>n===2),'closed manifold');
 }
 return h;
}
const tetra=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(point);
check(tetra,[4,6,4]);check([...tetra,tetra[0]], [4,6,4]);
check([-1,1].flatMap(x=>[-1,1].flatMap(y=>[-1,1].map(z=>point([x,y,z])))),[8,12,6]);
check([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]].map(point),[6,12,8]);
const phi=(1+Math.sqrt(5))/2,icosa=[];for(const a of [-1,1])for(const b of [-phi,phi])icosa.push(point([0,a,b]),point([a,b,0]),point([b,0,a]));check(icosa,[12,30,20]);
assert.equal(check([], [0,0,0]).dimension,0);assert.equal(check([point([0,0,1])],[1,0,0]).dimension,0);assert.equal(check([point([0,0,1]),point([0,0,-1])],[2,0,0]).dimension,1);
const ring=Array.from({length:9},(_,i)=>point([Math.cos(i*2*Math.PI/9),Math.sin(i*2*Math.PI/9),.3]));assert.equal(check(ring,[9,9,1]).dimension,2);
let seed=131;const random=()=>((seed=(1664525*seed+1013904223)>>>0)/2**32);check(Array.from({length:150},()=>point([random()-.5,random()-.5,random()-.5])));
let records=0;
for(const name of ['symmetric-sector-catalogue','anticoherent-state-catalogue']){
 const data=JSON.parse(fs.readFileSync(path.join(root,'..',name,'data/majorana.json')));for(const r of Object.values(data.records)){try{check(r.points);}catch(e){console.error("Failed record",r.id);throw e;}records++;}
}
const degeneracy=JSON.parse(fs.readFileSync(path.join(root,'../majorana-degeneracy-catalogue/data/catalogue.json')));
for(const r of [...Object.values(degeneracy.states||{}),...Object.values(degeneracy.numericalStates||{})])if(r.majoranaPoints?.length){try{check(r.majoranaPoints);}catch(e){console.error("Failed record",r.id);throw e;}records++;}
const subspaces=JSON.parse(fs.readFileSync(path.join(root,'../ac-subspace-catalogue/data/catalogue.json')));
for(const subspace of Object.values(subspaces.subspaces))for(const r of subspace.basis||[])if(r.points?.length){try{check(r.points);}catch(e){console.error("Failed record",r.id);throw e;}records++;}
console.log(`PASS: tetrahedron, cube, octahedron, icosahedron, planar/degenerate sets, duplicates, random sphere and ${records} catalogue constellations; supporting faces, Euler identity, closed edges, unchanged coordinates.`);
