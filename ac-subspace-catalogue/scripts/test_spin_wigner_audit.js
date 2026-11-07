/* Independent Wigner audit: rotate the diagonal Stratonovich--Weyl kernel.
 * The diagonal tensor elements use a factorial Wigner-3j (Racah) sum, not
 * SpinWigner's spherical-harmonic recurrences or OperatorMath's Hahn sum.
 * node ac-subspace-catalogue/scripts/test_spin_wigner_audit.js
 */
'use strict';
const assert=require('node:assert/strict');
const O=require('../operator-math.js'),W=require('../spin-wigner.js');
const catalogue=require('../data/catalogue.json');
const PI4=4*Math.PI;
const fac=[1];for(let k=1;k<=100;k++)fac[k]=fac[k-1]*k;
const near=(a,b,eps=4e-11)=>assert.ok(Math.abs(a-b)<=eps*Math.max(1,Math.abs(a),Math.abs(b)),`${a} versus ${b}`);
const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]],add=(a,b)=>[a[0]+b[0],a[1]+b[1]],conj=z=>[z[0],-z[1]];
function t0(n,l,k){
 // sqrt(2l+1) times the Wigner-3j symbol (j l j;m 0 -m),
 // with its conventional phase cancelled by the tensor-operator phase.
 const delta=fac[l]*fac[n-l]*fac[l]/fac[n+l+1];
 let sum=0;
 for(let z=Math.max(0,k+l-n);z<=Math.min(l,k);z++)sum+=(z%2?-1:1)/(fac[z]**2*fac[l-z]**2*fac[k-z]*fac[n-k-l+z]);
 return Math.sqrt((2*l+1)*delta)*fac[n-k]*fac[k]*fac[l]*sum;
}
function kernelNorth(n){
 return Array.from({length:n+1},(_,k)=>Array.from({length:n+1},(_,l)=>Math.sqrt((n+1)*(2*l+1))*t0(n,l,k)/PI4).reduce((a,b)=>a+b,0));
}
function directKernel(model,direction){
 const n=model.n,len=Math.hypot(...direction),v=direction.map(x=>x/len),theta=Math.acos(Math.max(-1,Math.min(1,v[2]))),sin=Math.hypot(v[0],v[1]),axis=sin>1e-12?[-v[1]/sin,v[0]/sin,0]:[0,1,0];
 const U=O.rotation(n,axis,theta),kernel=kernelNorth(n);let result=0;
 for(const state of model.basis)for(let k=0;k<=n;k++){
  const a=state.reduce((sum,z,i)=>add(sum,mul(conj(U[i][k]),z)),[0,0]);result+=kernel[k]*(a[0]**2+a[1]**2)/model.dimension;
 }
 return result;
}
const directions=[[0,0,1],[0,0,-1],[1,0,0],[0,1,0],[.36,-.48,.8],[-2/3,2/3,1/3]];
function model(n,basis){return O.fromBasis({n,basis:basis.map(coefficients=>({coefficients}))});}
function rotate(model,axis,angle){
 const U=O.rotation(model.n,axis,angle),basis=model.basis.map(v=>U.map(row=>row.reduce((sum,z,i)=>add(sum,mul(z,v[i])),[0,0])));
 return O.fromBasis({n:model.n,basis});
}
function rodrigues(v,a,angle){const cross=[a[1]*v[2]-a[2]*v[1],a[2]*v[0]-a[0]*v[2],a[0]*v[1]-a[1]*v[0]],dot=a.reduce((s,x,i)=>s+x*v[i],0);return v.map((x,i)=>x*Math.cos(angle)+cross[i]*Math.sin(angle)+a[i]*dot*(1-Math.cos(angle)));}
function gauss(n){
 const result=[];
 for(let k=1;k<=n;k++){
  let z=Math.cos(Math.PI*(k-.25)/(n+.5)),p,previous,derivative;
  for(let q=0;q<30;q++){
   let a=1,b=z;for(let l=2;l<=n;l++){const c=((2*l-1)*z*b-(l-1)*a)/l;a=b;b=c;}
   p=n===1?z:b;previous=n===1?1:a;derivative=n*(z*p-previous)/(z*z-1);const step=p/derivative;z-=step;if(Math.abs(step)<1e-15)break;
  }
  result.push([z,2/((1-z*z)*derivative*derivative)]);
 }
 return result;
}
function integrals(w){
 const polar=gauss(w.n+3),cols=2*w.n+5;let integral=0,square=0;
 for(const [z,weight]of polar)for(let k=0;k<cols;k++){
  const phi=2*Math.PI*k/cols,r=Math.sqrt(1-z*z),q=w.value([r*Math.cos(phi),r*Math.sin(phi),z]),area=weight*2*Math.PI/cols;integral+=q*area;square+=q*q*area;
 }
 return {integral,square};
}
let comparisons=0,maximumKernelError=0;
const fixtures=[];
for(let n=1;n<=12;n++){
 const m=model(n,[Array.from({length:n+1},(_,k)=>[Math.cos(.37+k*.87),Math.sin(.21+k*.63)])]);fixtures.push(m);
 // Independently formed north kernel has unit integral and self-dual norm.
 const kernel=kernelNorth(n);near(kernel.reduce((a,b)=>a+b,0),(n+1)/PI4);
 for(let l=0;l<=n;l++)for(let h=0;h<=n;h++)near(Array.from({length:n+1},(_,k)=>t0(n,l,k)*t0(n,h,k)).reduce((a,b)=>a+b,0),l===h?1:0);
}
for(const id of ['acliege-Π_7_2_2','radical-n13-k2-t3','mod4-reflection-n16-d3-t2']){
 const found=catalogue.subspaces.find(s=>s.id===id);assert.ok(found,`Catalogue fixture missing ${id}`);fixtures.push(O.fromBasis(found));
}
for(const m of fixtures){
 const prepared=W.prepare(m),axis=[2/3,-1/3,2/3],angle=1.173,rotated=rotate(m,axis,angle),wr=W.prepare(rotated);
 for(const direction of directions){
  const actual=prepared.value(direction),expected=directKernel(m,direction),error=Math.abs(actual-expected);maximumKernelError=Math.max(maximumKernelError,error);near(actual,expected);comparisons++;
  near(wr.value(rodrigues(direction,axis,angle)),actual,8e-11);
 }
 const {integral,square}=integrals(prepared);near(integral,1);near(square,(m.n+1)*m.diagnostics.purity/PI4);
}
// The +y Bloch direction fixes the often-missed conjugation sign.
{
 const m=model(1,[[[Math.SQRT1_2,0],[0,Math.SQRT1_2]]]),w=W.prepare(m);
 near(w.value([0,1,0]),(1+Math.sqrt(3))/PI4);near(w.value([0,-1,0]),(1-Math.sqrt(3))/PI4);assert.ok(w.value([0,-1,0])<0);
}
console.log(JSON.stringify({status:'PASS',method:'Independent Racah-sum diagonal SW kernel and physical SU(2) rotations',fixtures:fixtures.length,pointComparisons:comparisons,maximumKernelError,checks:['spatial orientation including +y','arbitrary-axis covariance','unit integral','self-dual purity','tensor orthonormality','complex mixed catalogue projectors']}));
