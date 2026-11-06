'use strict';
const assert=require('node:assert/strict');
const M=require('./majorana-math.js');
const fs=require('node:fs');
let seed=0x61602;function rand(){seed=(1664525*seed+1013904223)>>>0;return seed/4294967296;}
const cmul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];
const conj=z=>[z[0],-z[1]];
const counts={states:0,stars:0,pauliStates:0};let worstFidelity=1,worstResidual=0,worstSectorError=0,worstKnownStarError=0;
function pauliSectors(c){
  const n=c.length-1,len=2**n,v=Array.from({length:len},(_,i)=>{let k=0,t=i;while(t){k+=t&1;t>>=1;}return c[k].map(q=>q/Math.sqrt(M.binomial(n,k)));});
  const sectors=Array(n+1).fill(0);
  for(let word=0;word<4**n;word++){
    let t=word,weight=0;const letters=[];for(let bit=0;bit<n;bit++){const l=t%4;t=Math.floor(t/4);letters.push(l);if(l)weight++;}
    let er=0,ei=0;
    for(let i=0;i<len;i++){
      let j=i,phase=[1,0];for(let bit=0;bit<n;bit++){const b=(i>>bit)&1,l=letters[bit];if(l===1)j^=1<<bit;else if(l===2){j^=1<<bit;phase=cmul(phase,[0,b?-1:1]);}else if(l===3&&b)phase=[-phase[0],-phase[1]];}
      const term=cmul(conj(v[j]),cmul(phase,v[i]));er+=term[0];ei+=term[1];
    }
    assert.ok(Math.abs(ei)<1e-12);sectors[weight]+=er*er;
  }
  return sectors.slice(1);
}
function check(c,label,pauli=false,knownStars=null){
  c=M.normalizeState(c);const n=c.length-1,result=M.majoranaStars(c),sectors=M.sectorLengths(c);
  counts.states++;counts.stars+=n;worstFidelity=Math.min(worstFidelity,result.fidelity);worstResidual=Math.max(worstResidual,result.residual);
  assert.equal(result.stars.length,n,label+' star count');
  assert.ok(result.fidelity>1-2e-8,label+' reconstruction fidelity '+result.fidelity);
  assert.ok(result.residual<1e-7,label+' homogeneous residual '+result.residual);
  assert.ok(sectors.every(s=>s>=-1e-10),label+' sector positivity');
  assert.ok(Math.abs(sectors.reduce((s,v)=>s+v,0)-(2**n-1))<1e-9,label+' pure total');
  if(pauli){counts.pauliStates++;const oracle=pauliSectors(c);const error=Math.max(...oracle.map((s,k)=>Math.abs(s-sectors[k])));worstSectorError=Math.max(worstSectorError,error);assert.ok(error<5e-10,label+' independent Pauli oracle '+error);}
  if(knownStars){
    let best=Infinity;function dfs(i,used,worst){if(worst>=best)return;if(i===n){best=worst;return;}for(let j=0;j<n;j++)if(!(used&(1<<j))){const s=result.stars[i],t=knownStars[j],d=Math.hypot(s.x-t.x,s.y-t.y,s.z-t.z);dfs(i+1,used|(1<<j),Math.max(worst,d));}}dfs(0,0,0);worstKnownStarError=Math.max(worstKnownStarError,best);assert.ok(best<2e-4,label+' known star matching '+best);
  }
}
function point(){const z=2*rand()-1,phi=2*Math.PI*rand(),r=Math.sqrt(1-z*z);return {x:r*Math.cos(phi),y:r*Math.sin(phi),z};}
for(let n=3;n<=6;n++){
  for(let k=0;k<=n;k++){const c=Array.from({length:n+1},(_,i)=>[i===k?1:0,0]);check(c,`Dicke ${n},${k}`,true,[...Array.from({length:n-k},()=>({x:0,y:0,z:1})),...Array.from({length:k},()=>({x:0,y:0,z:-1}))]);}
  const ghz=Array.from({length:n+1},(_,k)=>[k===0||k===n?Math.SQRT1_2:0,0]);check(ghz,`GHZ ${n}`,true);
  for(let sample=0;sample<60;sample++)check(Array.from({length:n+1},()=>[rand()-.5,rand()-.5]),`random coefficients ${n}/${sample}`,sample<3);
  for(let sample=0;sample<40;sample++){
    let stars=Array.from({length:n},()=>point());
    if(sample<6)stars=stars.map(()=>({...stars[0]}));
    else if(sample<12){const s={...stars[0]};stars=stars.map((t,k)=>k<n-1?{...s}:t);}
    else if(sample<18){stars[0]={x:0,y:0,z:1};stars[1]={x:0,y:0,z:-1};}
    check(M.stateFromStars(stars),`known stars ${n}/${sample}`,false,stars);
  }
  // Regression: a tiny leading Majorana coefficient creates one very distant
  // affine root. The adaptive SU(2) chart must preserve the physical stars.
  for(let exponent=2;exponent<=14;exponent+=2){
    const c=Array.from({length:n+1},(_,k)=>k===0?[10**(-exponent),0]:[rand()-.5,rand()-.5]);
    check(c,`ill-conditioned affine chart ${n}/1e-${exponent}`);
  }
}
const octa=Array.from({length:7},(_,k)=>[k===1?Math.SQRT1_2:k===5?-Math.SQRT1_2:0,0]);
check(octa,'octahedron',true,[{x:1,y:0,z:0},{x:-1,y:0,z:0},{x:0,y:1,z:0},{x:0,y:-1,z:0},{x:0,y:0,z:1},{x:0,y:0,z:-1}]);
const output={pass:true,counts,worstFidelity,worstResidual,worstSectorError,worstKnownStarError,convention:M.convention};
fs.writeFileSync(__dirname+'/majorana_math_audit.json',JSON.stringify(output,null,2)+'\n');console.log(JSON.stringify(output,null,2));
