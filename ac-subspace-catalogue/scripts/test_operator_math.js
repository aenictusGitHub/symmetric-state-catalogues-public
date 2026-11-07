/* Reproducible numerical geometry checks: node scripts/test_operator_math.js. */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const O=require('../operator-math.js');global.MajoranaSymmetry=require('../majorana-symmetry.js');
const catalogue=require('../data/catalogue.json');
const cadd=(a,b)=>[a[0]+b[0],a[1]+b[1]],cmul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]],scale=(a,s)=>a.map(x=>x*s);
const near=(a,b,t=2e-12)=>assert.ok(Math.abs(a-b)<t,`${a} versus ${b}`);
const apply=(U,v)=>U.map(row=>row.reduce((s,z,k)=>cadd(s,cmul(z,v[k])),[0,0]));
const fromVectors=(n,basis)=>O.fromBasis({n,basis:basis.map(coefficients=>({coefficients}))});
const vector=(n,k)=>Array.from({length:n+1},(_,i)=>[i===k?1:0,0]);
const report={method:'Normalized projector P/d; exact-integer Hahn sums for tensor operators; numerical roots and signed-subconstellation reconstruction.',fixtures:{}};
// Complete Hilbert space: identity/(N+1), constant coherent Q and no stars.
for(let n=1;n<=8;n++){const m=fromVectors(n,Array.from({length:n+1},(_,k)=>vector(n,k)));near(m.diagnostics.purity,1/(n+1));for(let i=1;i<=n;i++)assert.ok(m.multipoles[i].weight<1e-13);near(O.husimi(m,[.6,0,.8]),1/(n+1));assert.ok(O.symmetry(m).isotropic);}
// Signed rank1: ±Bloch directions have identical full rootsets, different selected classes.
const north=fromVectors(1,[vector(1,0)]),south=fromVectors(1,[vector(1,1)]),cn=O.constellation(north,1),cs=O.constellation(south,1);
near(north.multipoles[1].weight,Math.SQRT1_2);assert.equal(cn.available,true);assert.equal(cs.available,true);assert.ok(cn.selectedPoints[0].z>.9999);assert.ok(cs.selectedPoints[0].z<-.9999);assert.ok(O.rotationResidual(north,[1,0,0],Math.PI)>1);const ns=O.symmetry(north);assert.ok(ns.continuous);assert.equal(ns.perpendicularHalfTurns,false);report.fixtures.signedRankOne='PASS';
// Known current catalogue example: only rank6, icosahedral12stars and31axes.
const s=catalogue.subspaces.find(s=>s.id==='acliege-Π_7_2_2'),m=O.fromBasis(s),b=m.multipoles[6];near(b.weight,Math.sqrt(3/8));near(b.coefficients[1][0],Math.sqrt(42)/20);near(b.coefficients[6][0],Math.sqrt(66)/20);near(b.coefficients[11][0],-Math.sqrt(42)/20);const c=O.constellation(m,6);assert.ok(c.available);assert.equal(c.points.reduce((s,p)=>s+p.multiplicity,0),12);const sym=O.symmetry(m),counts={};for(const a of sym.axes)counts[a.order]=(counts[a.order]||0)+1;assert.deepEqual(counts,{2:15,3:10,5:6});report.fixtures.N7={weights:m.multipoles.map(x=>x.weight),axisCounts:counts,rootResidual:c.rootResidual,signReconstructionResidual:c.reconstructionResidual};
// Complex change of basis leaves projector, weights, coherent Q invariant.
const u=m.basis[0],v=m.basis[1],mixed=fromVectors(7,[u.map((z,k)=>scale(cadd(z,cmul([0,1],v[k])),Math.SQRT1_2)),u.map((z,k)=>scale(cadd(cmul([0,1],z),v[k]),Math.SQRT1_2))]);
for(let a=0;a<8;a++)for(let b=0;b<8;b++)near(Math.hypot(...m.rho[a][b].map((x,k)=>x-mixed.rho[a][b][k])),0);report.fixtures.basisInvariance='PASS';
// Full rotation covariance, using the physical SU2 representation.
const axis=[.3,.4,.5].map(x=>x/Math.sqrt(.5)),angle=1.137,U=O.rotation(7,axis,angle),rot=fromVectors(7,m.basis.map(v=>apply(U,v)));for(let l=0;l<=7;l++)near(m.multipoles[l].weight,rot.multipoles[l].weight,2e-13);const q=[.6,.8,0],cross=[axis[1]*q[2]-axis[2]*q[1],axis[2]*q[0]-axis[0]*q[2],axis[0]*q[1]-axis[1]*q[0]],dot=axis.reduce((s,x,i)=>s+x*q[i],0),rq=q.map((x,i)=>x*Math.cos(angle)+cross[i]*Math.sin(angle)+axis[i]*dot*(1-Math.cos(angle)));near(O.husimi(m,q),O.husimi(rot,rq),2e-13);const rc=O.constellation(rot,6);assert.ok(rc.available);report.fixtures.rotationCovariance='PASS';
// Every catalogue entry, every populated operator rank. Failure is an explicit
// precision outcome, not silently fabricated stars.
let maxPurityResidual=0,maxProjectorResidual=0,maxGramResidual=0,blocks=0,available=0,maxRootResidual=0,maxSignResidual=0;const failures=[];
for(const subspace of catalogue.subspaces){const mm=O.fromBasis(subspace);maxPurityResidual=Math.max(maxPurityResidual,mm.diagnostics.multipolePurityResidual);maxProjectorResidual=Math.max(maxProjectorResidual,mm.diagnostics.projectorIdempotencyResidual);maxGramResidual=Math.max(maxGramResidual,mm.diagnostics.gramResidual);near(mm.diagnostics.purity,1/mm.dimension,1e-12);assert.ok(mm.diagnostics.multipolePurityResidual<1e-12);for(const r of mm.multipoles.slice(1)){if(r.zero)continue;blocks++;const cc=O.constellation(mm,r.rank);if(cc.available){available++;assert.equal(cc.signReliable,true);assert.equal(cc.points.reduce((s,p)=>s+p.multiplicity,0),2*r.rank);maxRootResidual=Math.max(maxRootResidual,cc.rootResidual);maxSignResidual=Math.max(maxSignResidual,cc.reconstructionResidual);}else failures.push({id:subspace.id,rank:r.rank,error:cc.error});}}
report.catalogue={subspaces:catalogue.subspaces.length,maximumN:Math.max(...catalogue.numbers),populatedBlocks:blocks,availableBlocks:available,unresolvedBlocks:failures,maxPurityResidual,maxProjectorResidual,maxGramResidual,maxRootResidual,maxSignResidual};
fs.writeFileSync(path.join(__dirname,'operator-math-audit.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report,null,2));
