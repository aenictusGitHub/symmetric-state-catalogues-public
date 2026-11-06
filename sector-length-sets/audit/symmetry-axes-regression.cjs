'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const root=path.resolve(__dirname,'..');
const math=require(path.join(root,'majorana-math.js'));
const symmetry=require(path.join(root,'symmetry-axes.js'));
const small=JSON.parse(fs.readFileSync(path.join(root,'data/small-n.json')));
const six=JSON.parse(fs.readFileSync(path.join(root,'data/n6.json')));
const groups=coords=>coords.map(p=>({x:p[0],y:p[1],z:p[2],multiplicity:p[3]||1}));
const unit=p=>{const norm=Math.hypot(...p);return p.map(x=>x/norm);};
function rotate(p,axis,angle){axis=unit(axis);const c=Math.cos(angle),s=Math.sin(angle),d=p.reduce((sum,v,i)=>sum+v*axis[i],0),cross=[axis[1]*p[2]-axis[2]*p[1],axis[2]*p[0]-axis[0]*p[2],axis[0]*p[1]-axis[1]*p[0]];return p.map((v,i)=>c*v+s*cross[i]+(1-c)*d*axis[i]);}
function rotated(g,axis,angle){return g.map(p=>{const q=rotate([p.x,p.y,p.z],axis,angle);return {x:q[0],y:q[1],z:q[2],multiplicity:p.multiplicity};});}
const polygon=(n,z=0)=>groups(Array.from({length:n},(_,k)=>{const a=2*Math.PI*k/n,r=Math.sqrt(1-z*z);return [r*Math.cos(a),r*Math.sin(a),z];}));
const cases=[
  ['octahedron',groups([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]),{2:6,3:4,4:3},24],
  ['tetrahedron',groups([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(unit)),{2:3,3:4},12],
  ['equatorial pentagon',polygon(5),{2:5,5:1},10],
  ['latitude pentagon',polygon(5,.2),{5:1},5],
  ['equatorial hexagon',polygon(6),{2:6,6:1},12],
  ['coherent',groups([[0,0,1,6]]),{infinity:1},null],
  ['balanced Dicke',groups([[0,0,1,3],[0,0,-1,3]]),{infinity:1},null],
  ['unbalanced Dicke',groups([[0,0,1,2],[0,0,-1,4]]),{infinity:1},null],
  ['unequal non-antipodal stars',groups([[0,0,1,2],[1,0,0,1]]),{},1],
  ['equal non-antipodal stars',groups([[0,0,1,2],[1,0,0,2]]),{2:1},2],
];
const results=[];
for(const[name,input,expected,count] of cases){
  for(const rotation of [null,[[1,2,3],.713],[[.1,-.5,.3],2.1]]){
    const g=rotation?rotated(input,...rotation):input,out=symmetry.analyze(g);
    const histogram={};for(const a of out.axes)histogram[a.order]=(histogram[a.order]||0)+1;
    assert.deepEqual(histogram,expected,name);assert.equal(out.rotationCount,count,name);
    if(name==='balanced Dicke')assert.equal(out.continuous.transverseC2,true);
    if(name==='unbalanced Dicke')assert.equal(out.continuous.transverseC2,false);
    results.push({name,rotated:!!rotation,histogram,rotationCount:out.rotationCount});
  }
}
function preservationError(groups,direction,order){
  const theta=order==='infinity'?.473:2*Math.PI/order,used=new Set();let error=0;
  for(const g of groups){const q=rotate([g.x,g.y,g.z],direction,theta);let best=Infinity,index=-1;
    groups.forEach((h,j)=>{if(used.has(j)||g.multiplicity!==h.multiplicity)return;const e=Math.hypot(q[0]-h.x,q[1]-h.y,q[2]-h.z);if(e<best){best=e;index=j;}});
    assert.ok(index>=0);used.add(index);error=Math.max(error,best);
  }return error;
}
let maxError=0,states=0;
const nodes=JSON.parse(fs.readFileSync(path.join(root,'data/n6-manuscript-nodes.json'))),ids=new Set(nodes.flatMap(n=>n.stateIds));
const actual=[...Object.values(small).flatMap(d=>d.knownStates).filter(s=>s.id!=='n5-dicke2'),...six.states.filter(s=>ids.has(s.id))];
const chosen=[],near=(a,b)=>Math.hypot(a.sectors[0]-b.sectors[0],a.sectors[1]-b.sectors[1])<.02;
for(const curve of small['5'].curves)for(let j=1;j<10;j++){const p=curve.points[Math.round(j*(curve.points.length-1)/10)];if(!chosen.some(q=>near(p,q))&&!small['5'].knownStates.some(q=>near(p,q)))chosen.push(p);}
assert.equal(chosen.length,41);actual.push(...chosen);assert.equal(actual.length,72);assert.ok(actual.some(s=>s.id==='n5-top-126'));
for(const state of actual){const g=math.majoranaStars(state.coefficients).groups,out=symmetry.analyze(g);for(const a of out.axes){const e=preservationError(g,a.direction,a.order);assert.ok(e<1e-5,state.id);maxError=Math.max(maxError,e);}states++;
  if(state.id==='n6-octa')assert.deepEqual(out.axes.reduce((h,a)=>(h[a.order]=(h[a.order]||0)+1,h),{}),{2:6,3:4,4:3});
  if(state.id==='n4-tetra')assert.deepEqual(out.axes.reduce((h,a)=>(h[a.order]=(h[a.order]||0)+1,h),{}),{2:3,3:4});
  if(state.id==='n5-top-126')assert.deepEqual(out.axes.map(a=>a.order),[5]);
}
const report={pass:true,syntheticCases:results.length,actualStates:states,maximumAxisPreservationError:maxError,results};
fs.writeFileSync(path.join(__dirname,'symmetry-axes-regression.json'),JSON.stringify(report,null,2));
console.log(JSON.stringify({...report,results:undefined},null,2));
