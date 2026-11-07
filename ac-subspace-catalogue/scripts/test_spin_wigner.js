'use strict';
const assert=require('node:assert/strict');
const OperatorMath=require('../operator-math.js'),SpinWigner=require('../spin-wigner.js');
const PI4=4*Math.PI;
const near=(a,b,tolerance=3e-12,message='Mismatch')=>assert.ok(Math.abs(a-b)<=tolerance,`${message}: ${a} versus ${b}`);
const complexMultiply=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];
function fromVectors(n,vectors){return OperatorMath.fromBasis({n,basis:vectors.map(coefficients=>({coefficients}))});}
function pure(coefficients){return fromVectors(coefficients.length-1,[coefficients]);}
function unitVector(n,k){return Array.from({length:n+1},(_,i)=>[i===k?1:0,0]);}
function gauss(count){
 const nodes=[],weights=[];
 for(let i=0;i<count;i++){
  let x=Math.cos(Math.PI*(i+.75)/(count+.5)),derivative=0;
  for(let step=0;step<30;step++){
   let a=1,b=x;for(let l=2;l<=count;l++){const c=((2*l-1)*x*b-(l-1)*a)/l;a=b;b=c;}
   const value=count===1?x:b,previous=count===1?1:a;derivative=count*(x*value-previous)/(x*x-1);
   const shift=value/derivative;x-=shift;if(Math.abs(shift)<1e-15)break;
  }
  nodes.push(x);weights.push(2/((1-x*x)*derivative*derivative));
 }
 return {nodes,weights};
}
function integrate(prepared,power=1,observable=()=>1){
 const polar=gauss(prepared.n+3),cols=2*prepared.n+5;let sum=0;
 for(let i=0;i<polar.nodes.length;i++)for(let j=0;j<cols;j++){
  const z=polar.nodes[i],phi=2*Math.PI*j/cols,r=Math.sqrt(1-z*z),v=[r*Math.cos(phi),r*Math.sin(phi),z];
  sum+=polar.weights[i]*Math.pow(prepared.value(v),power)*observable(v)*2*Math.PI/cols;
 }
 return sum;
}
// Independent closed-form qubit Wigner functions fix all spatial signs.
for(const [state,axis]of [
 [[[1,0],[0,0]],[0,0,1]],
 [[[0,0],[1,0]],[0,0,-1]],
 [[[Math.SQRT1_2,0],[Math.SQRT1_2,0]],[1,0,0]],
 [[[Math.SQRT1_2,0],[0,Math.SQRT1_2]],[0,1,0]]
]){
 const w=SpinWigner.prepare(pure(state));
 for(const v of [[1,0,0],[0,1,0],[0,0,1],[0,0,-1],[.36,-.48,.8]])near(w.value(v),(1+Math.sqrt(3)*v.reduce((s,x,i)=>s+x*axis[i],0))/PI4);
 near(integrate(w),1);near(integrate(w,2),2/PI4);
}
assert.ok(SpinWigner.prepare(pure([[1,0],[0,0]])).value([0,0,-1])<0,'Negative values must be retained.');
// Maximally mixed projectors are isotropic for every dimension checked.
for(let n=0;n<=8;n++){
 const model=fromVectors(n,Array.from({length:n+1},(_,k)=>unitVector(n,k))),w=SpinWigner.prepare(model);
 for(const v of [[1,0,0],[0,-1,0],[0,0,1],[.36,-.48,.8]])near(w.value(v),1/PI4);
 near(integrate(w),1);near(integrate(w,2),1/PI4);
}
// Deterministic complex states test normalization and the SW tracial identity.
for(let n=2;n<=12;n++){
 const model=pure(Array.from({length:n+1},(_,k)=>[Math.sin(.37+1.31*k),Math.cos(.19+.73*k)])),w=SpinWigner.prepare(model);
 near(integrate(w),1,2e-11,'Integral normalization');near(integrate(w,2),(n+1)/PI4,2e-11,'Self-dual purity identity');
 const grid=w.grid(18,37);assert.equal(grid.values.length,19*38);assert.equal(grid.vertices.length,19*38);assert.equal(grid.faces.length,18*37);
 for(let k=0;k<grid.values.length;k++){near(grid.vertices[k].q,grid.values[k],0);near(Math.hypot(...grid.vertices[k].n),1,3e-15);}
 for(const face of grid.faces)assert.ok(face.length===4&&face.every(k=>Number.isInteger(k)&&k>=0&&k<grid.vertices.length));
 for(let i=0;i<=18;i++){
  near(grid.values[i*38],grid.values[i*38+37],0,'Grid seam');
  for(const j of [0,7,19,36])near(grid.values[i*38+j],w.at(Math.PI*i/18,2*Math.PI*j/37),2e-13,'Cached polar grid');
 }
 assert.ok(Number.isFinite(grid.min)&&Number.isFinite(grid.max)&&grid.min<=grid.max);
}
// Arbitrary-axis covariance is independently checked using spin rotations.
{
 const n=6,c=Array.from({length:n+1},(_,k)=>[Math.cos(k*.71),Math.sin(.17+k*.39)]),model=pure(c),axis=[2/3,-1/3,2/3],angle=.697,U=OperatorMath.rotation(n,axis,angle);
 const state=model.basis[0],rotated=U.map(row=>row.reduce((s,z,k)=>{const v=complexMultiply(z,state[k]);return [s[0]+v[0],s[1]+v[1]];},[0,0]));
 const before=SpinWigner.prepare(model),after=SpinWigner.prepare(pure(rotated));
 for(const v of [[1,0,0],[0,0,1],[.36,-.48,.8]]){
  const c=Math.cos(-angle),s=Math.sin(-angle),d=axis.reduce((sum,a,i)=>sum+a*v[i],0),cross=[axis[1]*v[2]-axis[2]*v[1],axis[2]*v[0]-axis[0]*v[2],axis[0]*v[1]-axis[1]*v[0]],inverse=v.map((x,i)=>c*x+s*cross[i]+(1-c)*d*axis[i]);
  near(after.value(v),before.value(inverse),4e-12,'Rotational covariance');
 }
}
// A tetrahedral pure state is 2-AC: its Wigner function has no dipole or quadrupole.
{
 const c=unitVector(4,0);c[0]=[Math.sqrt(1/3),0];c[3]=[Math.sqrt(2/3),0];const model=pure(c),w=SpinWigner.prepare(model);
 assert.ok(model.multipoles[1].weight<1e-13&&model.multipoles[2].weight<1e-13);
 for(let i=0;i<3;i++)near(integrate(w,1,v=>v[i]),0,2e-12);
 for(let i=0;i<3;i++)for(let j=0;j<3;j++)near(integrate(w,1,v=>v[i]*v[j]-(i===j?1/3:0)),0,2e-12);
 assert.ok(w.grid(24,48).min<0);
}
// Rank42 exercises the catalogue limit and the efficient regular-grid path.
{
 const model=pure(unitVector(42,0)),w=SpinWigner.prepare(model),grid=w.grid(128,256);
 assert.equal(grid.values.length,129*257);assert.ok(grid.values.every(Number.isFinite));
 near(integrate(w),1,2e-10);near(integrate(w,2),43/PI4,2e-10);
 near(w.value([0,0,1]),grid.values[0]);assert.ok(w.value([0,0,1])>w.value([0,0,-1]));
}
assert.throws(()=>SpinWigner.prepare({n:1,multipoles:[]}),/Complete/);
assert.throws(()=>SpinWigner.prepare(pure([[1,0],[0,0]])).value([0,0,0]),/nonzero/);
const bad={n:1,multipoles:[{rank:0,coefficients:[[1,0]]},{rank:1,coefficients:[[1,0],[0,0],[1,0]]}]};assert.throws(()=>SpinWigner.prepare(bad),/Hermitian/);
console.log('PASS: analytic qubit axes and negativity, isotropic projectors, unit integral, self-dual purity, arbitrary-axis covariance, AC low moments, cached grid/seams and finite N=42.');
