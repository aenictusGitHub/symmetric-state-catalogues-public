'use strict';
// Independent browser-math audit of serialized atlas states and boundary samples.
// Run with Node; it does not import either data-generator implementation.
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const M=require('./majorana-math.js');
const base=__dirname,small=JSON.parse(fs.readFileSync(path.join(base,'data/small-n.json'))),n6=JSON.parse(fs.readFileSync(path.join(base,'data/n6.json')));
const stats={pass:true,allStateRecords:0,majoranaRecords:0,byN:{},maxNormalizationError:0,maxStoredSectorError:0,maxStoredPointError:0,minReconstructionFidelity:1,maxPolynomialResidual:0,maxPureTotalError:0,minRefinedFourCopySlack:Infinity,knownN6:{},n5Boundary:{}};
const problems=[];
function checkState(p,n,label,majorana=false){
  if(!p.coefficients)return;
  stats.allStateRecords++;stats.byN[n]=(stats.byN[n]||0)+1;
  if(p.coefficients.length!==n+1)problems.push(label+': wrong coefficient count');
  const norm=p.coefficients.reduce((s,z)=>s+z[0]*z[0]+z[1]*z[1],0);
  stats.maxNormalizationError=Math.max(stats.maxNormalizationError,Math.abs(norm-1));
  const sectors=M.sectorLengths(p.coefficients);
  if(p.sectors){const error=Math.max(...p.sectors.map((s,k)=>Math.abs(s-sectors[k])));stats.maxStoredSectorError=Math.max(stats.maxStoredSectorError,error);if(error>1e-7)problems.push(label+': sector mismatch '+error);}
  if(p.point){const error=Math.max(...p.point.map((s,k)=>Math.abs(s-sectors[k])));stats.maxStoredPointError=Math.max(stats.maxStoredPointError,error);if(error>1e-7)problems.push(label+': point mismatch '+error);}
  stats.maxPureTotalError=Math.max(stats.maxPureTotalError,Math.abs(sectors.reduce((s,x)=>s+x,0)-(2**n-1)));
  if(n===6){const [x,y,z]=sectors;const f=5*(y+z-5*x-5)/21+((162+8*x-10*y-3*z)/231)**2-(8+2*x-z)*(401+137*x+35*y+5*z)/51744;stats.minRefinedFourCopySlack=Math.min(stats.minRefinedFourCopySlack,f);if(f< -2e-9)problems.push(label+': violates proven four-copy bound '+f);}
  if(majorana){const r=M.majoranaStars(p.coefficients);stats.majoranaRecords++;stats.minReconstructionFidelity=Math.min(stats.minReconstructionFidelity,r.fidelity);stats.maxPolynomialResidual=Math.max(stats.maxPolynomialResidual,r.residual);if(r.stars.length!==n||r.fidelity<1-1e-8||r.residual>1e-7)problems.push(label+': Majorana reconstruction failed '+JSON.stringify(r));}
  return sectors;
}
for(const [key,d]of Object.entries(small)){
  const n=Number(key);
  for(const p of d.knownStates)checkState(p,n,p.id,true);
  for(const curve of d.curves){for(let i=0;i<curve.points.length;i++)checkState(curve.points[i],n,curve.points[i].id,i%7===0||i===curve.points.length-1);}
  for(let i=0;i<d.cloud.length;i++)checkState(d.cloud[i],n,d.cloud[i].id,i%37===0);
}
const families=new Map();for(const p of n6.witnesses){if(!families.has(p.family))families.set(p.family,[]);families.get(p.family).push(p);}
for(const p of n6.states){const s=checkState(p,6,p.id,true);stats.knownN6[p.id]={label:p.label,sectors:s,majoranaFidelity:M.majoranaStars(p.coefficients).fidelity};}
for(const [family,records]of families){const stride=Math.max(1,Math.floor(records.length/24));for(let i=0;i<records.length;i++)checkState(records[i],6,records[i].id,i%stride===0||i===records.length-1);}

// Dedicated N=5 boundary audit. This checks that the serialized witnesses
// implement the claimed exact branches; it does not replace their global proof.
const curves=Object.fromEntries(small['5'].curves.map(c=>[c.id,c]));
const pc=2-3/Math.sqrt(5),p0=Math.sqrt(5)/3,xc=(22805-9840*Math.sqrt(5))/11881;
function lowerSmall(p){const r=Math.max(0,2*p*(5*p*p+6*Math.sqrt(5)*p-11)/(3*Math.sqrt(5)*p*p+2*p+3*Math.sqrt(5)));return M.sectorLengths([[-p,0],[0,0],[Math.sqrt(r),0],[0,0],[1,0],[0,0]]).slice(0,2);}
let errSmall=0,errRoot=0,errLine=0,errTop=0,errLeft=0,maxStepX=0;
for(const p of curves['n5-small-arc'].points){const t=lowerSmall(p.parameter.p);errSmall=Math.max(errSmall,...p.sectors.slice(0,2).map((s,k)=>Math.abs(s-t[k])));}
for(const p of curves['n5-root-arc'].points){const [x,y]=p.sectors;errRoot=Math.max(errRoot,Math.abs(y-(3*x-2*Math.sqrt(5*x)+15)/4));}
for(const p of curves['n5-linear-arc'].points)errLine=Math.max(errLine,Math.abs(p.sectors[1]-2*p.sectors[0]));
for(const p of curves['n5-top'].points)errTop=Math.max(errTop,Math.abs(p.sectors[1]-10));
for(const id of ['n5-left-negative','n5-left-positive'])for(const p of curves[id].points)errLeft=Math.max(errLeft,Math.abs(p.sectors[0]));
const points=curves['n5-small-arc'].points;for(let k=1;k<points.length;k++)maxStepX=Math.max(maxStepX,points[k-1].sectors[0]-points[k].sectors[0]);
const smallStart=points[0].sectors.slice(0,2),smallEnd=points.at(-1).sectors.slice(0,2),rootStart=curves['n5-root-arc'].points[0].sectors.slice(0,2),rootEnd=curves['n5-root-arc'].points.at(-1).sectors.slice(0,2);
assert.ok(Math.abs(smallStart[0])<1e-10&&Math.abs(smallStart[1]-3.6)<1e-10);
assert.ok(Math.abs(smallEnd[0]-xc)<1e-10&&Math.max(...smallEnd.map((s,k)=>Math.abs(s-rootStart[k])))<1e-10);
assert.ok(Math.abs(rootEnd[0]-1.8)<1e-10&&Math.abs(rootEnd[1]-3.6)<1e-10);
assert.ok(Math.max(errSmall,errRoot,errLine,errTop,errLeft,maxStepX)<1e-8);
function lowerY(x){if(x>=1.8)return 2*x;if(x>=xc)return (3*x-2*Math.sqrt(5*x)+15)/4;let lo=pc,hi=p0;for(let i=0;i<48;i++){const mid=(lo+hi)/2;if(lowerSmall(mid)[0]>x)lo=mid;else hi=mid;}return lowerSmall((lo+hi)/2)[1];}
let minCloudGap=Infinity;for(const p of small['5'].cloud){const [x,y]=p.sectors;minCloudGap=Math.min(minCloudGap,y-lowerY(x));if(y<lowerY(x)-1e-8||y>10+1e-8)problems.push(p.id+': outside claimed N5 range');}
stats.n5Boundary={source:'Dedicated completed compact proof, not older unfinished draft',junctionX:xc,smallStart,smallEnd,rootStart,rootEnd,errSmall,errRoot,errLine,errTop,errLeft,maxBackwardsStepInSmallArc:maxStepX,minCloudGap};
stats.families=Object.fromEntries([...families].map(([k,v])=>[k,v.length]));
stats.pass=problems.length===0;stats.problems=problems;
fs.writeFileSync(path.join(base,'atlas_scientific_audit.json'),JSON.stringify(stats,null,2)+'\n');
console.log(JSON.stringify({...stats,knownN6:Object.keys(stats.knownN6)},null,2));
if(problems.length)process.exit(1);
