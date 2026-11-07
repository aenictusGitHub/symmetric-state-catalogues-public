/* node scripts/test_husimi_rendering.js */
'use strict';
const assert=require('node:assert/strict');
const H=require('../husimi-rendering.js');
function mesh(field,rows=64,cols=128){
 const vertices=[],faces=[];let qMin=Infinity,qMax=-Infinity;
 for(let i=0;i<=rows;i++)for(let j=0;j<=cols;j++){
  const t=Math.PI*i/rows,p=2*Math.PI*j/cols,n=[Math.sin(t)*Math.cos(p),Math.sin(t)*Math.sin(p),Math.cos(t)],q=field(n);
  vertices.push({n,q});qMin=Math.min(qMin,q);qMax=Math.max(qMax,q);
 }
 for(let i=0;i<rows;i++)for(let j=0;j<cols;j++){const a=i*(cols+1)+j,b=a+cols+1;faces.push([a,b,b+1,a+1]);}
 return{vertices,faces,qMin,qMax,rows,cols};
}
assert.equal(H.color(0),'#220d47');assert.equal(H.color(1),'#fdda8e');
for(const name of Object.keys(H.palettes)){assert.deepEqual(H.rgb(-3,name),H.rgb(0,name));assert.deepEqual(H.rgb(3,name),H.rgb(1,name));for(let i=0;i<100;i++)assert.ok(H.rgb(i/99,name).every(v=>v>=0&&v<=255));assert.ok(H.gradient(name).includes('100%'));}
assert.deepEqual(H.rgb(-3),H.rgb(0));assert.deepEqual(H.rgb(3),H.rgb(1));
for(let k=0;k<100;k++)assert.ok(H.rgb(k/99).every(v=>v>=0&&v<=255));
assert.equal(H.defaultLevels(mesh(n=>1/9)).length,0);
assert.equal(H.contours(mesh(n=>1/9)).segments.length,0);
const m=mesh(n=>n[2]);const levels=[-.63,.17,.57],result=H.contours(m,levels);
assert.ok(result.segments.length>100);
for(const level of levels){
 const segments=result.segments.filter(s=>s.level===level),degree=new Map();
 let maxAnalyticalError=0;
 for(const s of segments)for(const p of [s.a,s.b]){
  assert.ok(Math.abs(Math.hypot(...p)-1)<1e-12);
  maxAnalyticalError=Math.max(maxAnalyticalError,Math.abs(p[2]-level));
  const key=p.map(v=>Math.round(v/1e-8)).join(',');degree.set(key,(degree.get(key)||0)+1);
 }
 assert.ok(maxAnalyticalError<4e-4,'Sampled interpolation stays close to exact latitude');
 assert.ok([...degree.values()].every(d=>d===2),'Every contour closes, including the longitude seam');
}
// A level exactly on the equator still has no dangling or duplicate segments.
const equator=H.contours(m,[0]).segments,seen=new Set();
for(const s of equator){const id=[s.a,s.b].map(p=>p.map(v=>Math.round(v/1e-8)).join(',')).sort().join('|');assert.ok(!seen.has(id));seen.add(id);}
console.log(JSON.stringify({palette:'plum-rose-gold',alternativePalettes:2,latitudeLevels:levels,segments:result.segments.length,constantField:'PASS',unitSphere:'PASS',seamClosure:'PASS',equatorDuplicates:'PASS'}));
