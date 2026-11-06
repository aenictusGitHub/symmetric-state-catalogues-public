'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const root=path.resolve(__dirname,'..'),T=require(path.join(root,'vendor/three.min.js'));
const math=require(path.join(root,'majorana-math.js')),symmetry=require(path.join(root,'symmetry-axes.js'));
const small=JSON.parse(fs.readFileSync(path.join(root,'data/small-n.json'))),six=JSON.parse(fs.readFileSync(path.join(root,'data/n6.json')));
const manuscript=JSON.parse(fs.readFileSync(path.join(root,'data/n6-manuscript-nodes.json'))),ids=new Set(manuscript.flatMap(n=>n.stateIds));
const states=[...Object.values(small).flatMap(d=>d.knownStates).filter(s=>s.id!=='n5-dicke2'),...six.states.filter(s=>ids.has(s.id))];
const chosen=[],near=(a,b)=>Math.hypot(a.sectors[0]-b.sectors[0],a.sectors[1]-b.sectors[1])<.02;
for(const curve of small['5'].curves)for(let j=1;j<10;j++){
  const p=curve.points[Math.round(j*(curve.points.length-1)/10)];
  if(!chosen.some(q=>near(p,q))&&!small['5'].knownStates.some(q=>near(p,q)))chosen.push(p);
}
states.push(...chosen);assert.equal(states.length,72);
const toWorld=p=>new T.Vector3(p.x,p.z,-p.y),target=new T.Vector3(0,1,0);
const hist=a=>a.axes.reduce((out,v)=>(out[v.order]=(out[v.order]||0)+1,out),{});
let maxAxisError=0,maxDotError=0,maxInverseError=0,maxDeterminantError=0,aligned=0,continuous=0;
const noAxes=[];
for(const state of states){
  const originalCoefficients=JSON.stringify(state.coefficients),data=math.majoranaStars(state.coefficients),before=symmetry.analyze(data.groups);
  const highest=before.axes.reduce((best,a)=>!best||(a.order==='infinity'?Infinity:a.order)>(best.order==='infinity'?Infinity:best.order)?a:best,null);
  if(!highest){noAxes.push(state.id);continue;}
  const direction=new T.Vector3(highest.direction[0],highest.direction[2],-highest.direction[1]).normalize();
  const q=new T.Quaternion().setFromUnitVectors(direction,target),inverse=q.clone().invert();
  const alignmentError=direction.clone().applyQuaternion(q).distanceTo(target);
  const determinantError=Math.abs(new T.Matrix3().setFromMatrix4(new T.Matrix4().makeRotationFromQuaternion(q)).determinant()-1);
  assert.ok(alignmentError<2e-14,state.id);assert.ok(determinantError<2e-14,state.id);
  maxAxisError=Math.max(maxAxisError,alignmentError);maxDeterminantError=Math.max(maxDeterminantError,determinantError);
  const original=data.groups.map(toWorld),rotated=original.map(p=>p.clone().applyQuaternion(q));
  original.forEach((p,i)=>{
    const error=p.distanceTo(rotated[i].clone().applyQuaternion(inverse));assert.ok(error<2e-14,state.id);maxInverseError=Math.max(maxInverseError,error);
    original.forEach((other,j)=>{const e=Math.abs(p.dot(other)-rotated[i].dot(rotated[j]));assert.ok(e<2e-14,state.id);maxDotError=Math.max(maxDotError,e);});
  });
  const transformed=data.groups.map((g,i)=>({x:rotated[i].x,y:-rotated[i].z,z:rotated[i].y,multiplicity:g.multiplicity}));
  assert.deepEqual(transformed.map(g=>g.multiplicity),data.groups.map(g=>g.multiplicity));
  const after=symmetry.analyze(transformed);
  assert.deepEqual(hist(after),hist(before),state.id);assert.equal(after.rotationCount,before.rotationCount,state.id);
  assert.equal(JSON.stringify(state.coefficients),originalCoefficients,state.id);
  if(highest.order==='infinity')continuous++;aligned++;
}
// setFromUnitVectors must also handle its antipodal and identity limits.
for(const input of [new T.Vector3(0,1,0),new T.Vector3(0,-1,0)]){
  const q=new T.Quaternion().setFromUnitVectors(input,target);
  assert.ok(input.clone().applyQuaternion(q).distanceTo(target)<1e-14);
  assert.ok(Math.abs(q.length()-1)<1e-14);
}
// A generic constellation supplies the disabled/no-axis branch.
const generic=math.majoranaStars([[1,0],[.2,.3],[-.4,.17],[.35,-.27],[.19,.07]]).groups;
assert.equal(symmetry.analyze(generic).axes.length,0);
const report={pass:true,visibleStateRepresentatives:states.length,aligned,continuousAxes:continuous,noAxes,
  maximumAxisAlignmentError:maxAxisError,maximumPairwiseDotError:maxDotError,maximumInverseError:maxInverseError,
  maximumDeterminantError:maxDeterminantError,unchanged:['star multiplicities','finite rotational group orders/counts','stored Dicke amplitudes'],
  displayConvention:'Bloch (x,y,z) becomes viewer world (x,z,-y); viewer world +Y is displayed physical +z.'};
fs.writeFileSync(path.join(__dirname,'alignment-frame-regression.json'),JSON.stringify(report,null,2));
console.log(JSON.stringify(report,null,2));
