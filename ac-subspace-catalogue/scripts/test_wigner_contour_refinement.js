/* Scientific geometry regressions for zoomable Wigner isocontours.
 * node ac-subspace-catalogue/scripts/test_wigner_contour_refinement.js
 */
'use strict';
const assert=require('node:assert/strict');
const OperatorMath=require('../operator-math.js'),SpinWigner=require('../spin-wigner.js'),H=require('../husimi-rendering.js');
const catalogue=require('../data/catalogue.json');
const fixture=catalogue.subspaces.find(s=>s.id==='line-radical-n22-t5');
assert.ok(fixture,'The reported N=22 catalogue fixture must remain available.');
const wigner=SpinWigner.prepare(OperatorMath.fromBasis(fixture));
const near=(a,b,tolerance=2e-13)=>assert.ok(Math.abs(a-b)<=tolerance,`${a} versus ${b}`);
const normalize=a=>{const length=Math.hypot(...a);return a.map(x=>x/length);};
const distance=(a,b)=>Math.hypot(...a.map((x,k)=>x-b[k]));
const pairs=path=>{
 const result=[];
 for(let i=1;i<path.points.length;i++)result.push([path.points[i-1],path.points[i]]);
 if(path.closed&&distance(path.points[0],path.points.at(-1))>1e-12)result.push([path.points.at(-1),path.points[0]]);
 return result;
};
const distanceToSegment=(p,a,b)=>{
 const v=b.map((x,k)=>x-a[k]),length2=v.reduce((s,x)=>s+x*x,0);
 const alpha=length2?Math.max(0,Math.min(1,p.reduce((s,x,k)=>s+(x-a[k])*v[k],0)/length2)):0;
 return distance(p,a.map((x,k)=>x+alpha*v[k]));
};
const distanceToPath=(p,path)=>Math.min(...pairs(path).map(([a,b])=>distanceToSegment(p,a,b)));

// The compact grid changes storage, never the sampled Wigner field or seam.
const full=wigner.grid(20,40),compact=wigner.grid(20,40,{geometry:false});
assert.equal(compact.rows,20);assert.equal(compact.cols,40);
assert.ok(compact.values instanceof Float64Array);
assert.equal(compact.vertices,undefined);assert.equal(compact.faces,undefined);
assert.deepEqual(compact.values,full.values);
near(compact.qMin,full.qMin,0);near(compact.qMax,full.qMax,0);
for(let i=0;i<=compact.rows;i++)near(compact.values[i*41],compact.values[i*41+40],0);
const prescribed=[-.17,0,.31];
const fullContours=H.contours(full,prescribed),compactContours=H.contours(compact,prescribed);
assert.deepEqual(compactContours.levels,fullContours.levels);
assert.equal(compactContours.segments.length,fullContours.segments.length);
const canonicalSegment=s=>[s.a,s.b].map(p=>p.map(x=>Math.round(x/1e-10)).join(',')).sort().join('|')+':'+s.level;
assert.deepEqual(compactContours.segments.map(canonicalSegment).sort(),fullContours.segments.map(canonicalSegment).sort(),'Values-only extraction must agree with the equivalent full grid.');

// Joining a closed latitude crosses the azimuth seam without adding a chord.
const circle=Array.from({length:256},(_,i)=>{const phi=2*Math.PI*i/256,z=.23,r=Math.sqrt(1-z*z);return[r*Math.cos(phi),r*Math.sin(phi),z];});
const circleSegments=circle.map((a,i)=>({a,b:circle[(i+1)%circle.length],level:.23,t:.4}));
const unsimplified=H.contourPaths({segments:circleSegments},{tolerance:0});
assert.equal(unsimplified.length,1);assert.equal(unsimplified[0].closed,true);
const tolerance=1e-3,simplified=H.contourPaths({segments:circleSegments},{tolerance});
assert.equal(simplified.length,1);assert.equal(simplified[0].closed,true);
assert.equal(simplified[0].level,.23);assert.equal(simplified[0].t,.4);
assert.ok(simplified[0].points.length<circle.length*.75,'Simplification should remove unnecessary collinear-scale samples.');
for(const p of circle)assert.ok(distanceToPath(p,simplified[0])<=tolerance+1e-12,'Closed-path RDP deviation must stay bounded in physical 3D coordinates.');
for(const p of simplified[0].points)near(Math.hypot(...p),1,2e-13);

// A branching level set terminates its paths at the junction instead of
// inventing a smooth continuation across two arbitrarily selected branches.
const junction=[0,0,1],leaves=[0,2*Math.PI/3,4*Math.PI/3].map(phi=>[Math.sin(.2)*Math.cos(phi),Math.sin(.2)*Math.sin(phi),Math.cos(.2)]);
const branchPaths=H.contourPaths({segments:leaves.map(b=>({a:junction,b,level:0,t:.5}))},{tolerance:1e-4});
assert.equal(branchPaths.length,3);assert.ok(branchPaths.every(p=>!p.closed));
for(const p of branchPaths)assert.ok(distance(p.points[0],junction)<1e-12||distance(p.points.at(-1),junction)<1e-12);
const separatedLevels=H.contourPaths({segments:[...circleSegments,...circleSegments.map(s=>({...s,level:.7,t:.8}))]},{tolerance});
assert.equal(separatedLevels.length,2);assert.deepEqual(separatedLevels.map(p=>p.level).sort(),[.23,.7]);

// Keep the original displayed levels (including the physically special zero
// contour) while refining only their geometry on the denser scalar field.
const coarseMesh=wigner.grid(128,256),levels=H.defaultLevels(coarseMesh);levels.push(0);levels.sort((a,b)=>a-b);
const coarseContours=H.contours(coarseMesh,levels);
const denseMesh=wigner.grid(512,1024,{geometry:false}),denseContours=H.contours(denseMesh,levels);
assert.deepEqual(denseContours.levels,levels,'Dense sampling must not shift the displayed levels.');
const refinedPaths=H.contourPaths(denseContours,{tolerance:2e-5});
assert.ok(refinedPaths.length>0);assert.deepEqual([...new Set(refinedPaths.map(p=>p.level))].sort((a,b)=>a-b),levels);
assert.ok(refinedPaths.every(p=>p.closed),'Every component of this regular N=22 fixture must close without seam fragments.');
for(const path of refinedPaths){
 assert.ok(path.points.length>=2);
 for(const p of path.points)near(Math.hypot(...p),1,3e-12);
}
const coarsePairs=coarseContours.segments.map(s=>({a:s.a,b:s.b,level:s.level}));
const refinedPairs=refinedPaths.flatMap(p=>pairs(p).map(([a,b])=>({a,b,level:p.level})));
function scalarError(segments){
 let maximum=0,squares=0,count=0,zeroMaximum=0;
 const stride=Math.max(1,Math.floor(segments.length/4000));
 for(let i=0;i<segments.length;i+=stride){
  const {a,b,level}=segments[i];
  for(const p of[a,normalize(a.map((x,k)=>x+b[k]))]){
   const residual=Math.abs(wigner.value(p)-level);maximum=Math.max(maximum,residual);squares+=residual*residual;count++;
   if(level===0)zeroMaximum=Math.max(zeroMaximum,residual);
  }
 }
 return{maximum,rms:Math.sqrt(squares/count),zeroMaximum,samples:count};
}
const coarseError=scalarError(coarsePairs),refinedError=scalarError(refinedPairs);
assert.ok(refinedError.rms<coarseError.rms/8,'Refined contours must materially improve true Wigner level accuracy.');
assert.ok(refinedError.maximum<coarseError.maximum/5,'Largest sampled contour error must improve.');
assert.ok(refinedError.zeroMaximum<coarseError.zeroMaximum/5,'The zero contour requires the same geometric accuracy.');
console.log(JSON.stringify({status:'PASS',fixture:fixture.id,levels,coarseSegments:coarsePairs.length,refinedSegments:refinedPairs.length,refinedPaths:refinedPaths.length,coarseError,refinedError,checks:['compact-grid equivalence','sphere norm','azimuth seam and cycle closure','bounded RDP deviation','branch junctions','separate levels','fixed levels including zero','exact-field accuracy improvement']}));
