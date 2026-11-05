"""Independent validated exclusion verifier for pure N21 AC5.
Requires NumPy for archives and python-flint/Arb for outward-rounded arithmetic.
No optimization or discovery routines are imported. Floating arrays specify
exact binary rational certificate inputs, not unvalidated numerical answers.
"""
from pathlib import Path
import json,time,hashlib,numpy as np
import flint
from validated_bootstrap import A,M,Q,sixj_ball,frob_upper,square_frob_upper,bound,load_families,bootstrap
flint.ctx.prec=128;flint.ctx.threads=1
from angular import direct_coefficient

def prepare(root,cert):
 assert(cert['schema'],cert['n'],cert['t'],cert['copies'])==('n21-ac5-six-copy-validated-separation-v1',21,5,6)
 fam=load_families(root,cert);parent=fam[5]['meta']['source'];dual=np.load(root/cert['dualFile']);dualM=[];pairs=cert['pairs'];seen=set()
 for i,P in enumerate(pairs):
  si,ti=P['si'],P['ti'];assert 0<=si<=ti<len(parent);S,T=parent[si],parent[ti];assert(P['S'],P['T'],P['ps'],P['pt'])==(S['J'],T['J'],S['p'],T['p']);assert(si,ti)not in seen;seen.add((si,ti));v=dual['R'+str(i)];assert v.shape==(5,S['p'],T['p'])and np.all(np.isfinite(v))
  if si==ti:assert np.array_equal(v,v.swapaxes(1,2))
  dualM.append([M(v[ell].tolist())for ell in range(5)])
 return fam,parent,dualM

def source_matrix(B,fam,parent,pairs,dualM):
 J=B['J'];valid=[S for S in parent if abs(21-S['J'])<=J<=21+S['J']];offset={};off=0
 for S in valid:offset[S['J']]=(off,off+S['p']);off+=S['p']
 H=M(off,off);norm2=A(0)
 for P,D in zip(pairs,dualM):
  S,T=P['S'],P['T']
  if S not in offset or T not in offset:continue
  lo,hi=sorted((S,T));block=M(P['ps'],P['pt'])
  for ell in range(1,6):
   c=direct_coefficient(J,S,T,ell)
   if not c.is_zero():block+=D[ell-1]*c
  aa,bb=offset[S];cc,dd=offset[T]
  for a in range(P['ps']):
   for b in range(P['pt']):
    H[aa+a,cc+b]=block[a,b]
    if S!=T:H[cc+b,aa+a]=block[a,b]
  norm2+=(1 if S==T else 2)*square_frob_upper(block)
 V=M(fam[6]['compact'][J].tolist());G=V.transpose()*H*V
 return G,norm2.sqrt().upper()

def verify_lower_state():
 from fractions import Fraction
 from math import prod
 weights={0:915,4:2366,8:2023,13:2023,17:2366,21:915};sign={0:1,4:1,8:1,13:1,17:-1,21:1};den=10608
 moments=[]
 for k in range(5):
  x=sum(Fraction(w*(2*r-21)**k,den)for r,w in weights.items());u=sum(Fraction((2*r-21)**k,22)for r in range(22));assert x==u;moments.append(str(x))
 for q in range(1,4):assert not any(r+q in weights for r in weights)
 terms=[r for r in weights if r+4 in weights];assert terms==[0,4,13,17]
 square=lambda r,q:weights[r]*weights[r+q]*prod((21-r-k)*(r+k+1)for k in range(q))
 for r,t in [(0,17),(4,13)]:assert square(r,4)==square(t,4)and sign[r]*sign[r+4]==-sign[t]*sign[t+4]
 assert[r for r in weights if r+5 in weights]==[8]and square(8,5)>0
 return dict(exactACOrder=4,normalizationDenominatorSquared=den,doubledMagneticMomentsThrough4=moments,rank4CoherenceCancelledByMirrorPairs=True,rank5CoherenceNonzero=True)

def verify(root,regenerate_bootstrap=True):
 start=time.time();raw=(root/'certificate.json').read_bytes();cert=json.loads(raw);lowerWitness=verify_lower_state()
 if regenerate_bootstrap:delta,boot=bootstrap(root,cert)
 else:
  # This shortcut is for producing candidate congruence factors only.
  boot=json.loads((root/'bootstrap-verification.json').read_text());delta=A(Q(boot['maximumExactSixCopyCFPDistance']['upperRational']))
 fam,parent,dualM=prepare(root,cert);factors=np.load(root/'congruences.npz');checks=[];minimum=A(1)
 for B in fam[6]['meta']['source']:
  J=B['J'];G,hnorm=source_matrix(B,fam,parent,cert['pairs'],dualM);v=factors['J'+str(J)];assert v.shape==(B['p'],B['p'])and np.all(np.isfinite(v));T=M(v.tolist());residual=M(np.eye(B['p']).tolist())-T*G*T.transpose();rho=frob_upper(residual);tnorm2=square_frob_upper(T);gramError=((2*delta+delta*delta)*hnorm).upper();total=(rho+tnorm2*gramError).upper();assert total<1,(J,str(total));gap=((1-total)/tnorm2).lower();assert gap>0
  if gap<minimum:minimum=gap
  checks.append(dict(J=J,dimension=B['p'],congruenceResidual=bound(rho),ambientDualFrobeniusNorm=bound(hnorm),CFPInducedGramPerturbation=bound(gramError),inverseFactorFrobeniusNormSquared=bound(tnorm2),totalCongruenceError=bound(total),positiveEigenvalueLowerBound=dict(lowerRational=str(gap.fmpq()),diagnosticApproximation=float(gap))))
  print('validated positive source block',J,B['p'],'error upper',float(total),'seconds',time.time()-start,flush=True)
 report=dict(passed=True,n=21,t=5,copies=6,inputMomentDimension=296010,spinBlocks=len(checks),inputHighestWeightColumns=sum(b['p']for b in fam[6]['meta']['source']),realInvariantVariables=sum(b['p']*(b['p']+1)//2 for b in fam[6]['meta']['source']),allSourceBlocksStrictlyPositive=True,minimumCertifiedSourceEigenvalue=dict(lowerRational=str(minimum.fmpq()),diagnosticApproximation=float(minimum)),CFPBootstrap=boot,positiveSourceBlocks=checks,precisionBits=flint.ctx.prec,arbitraryComplexPureStatesCovered=True,exactFourthOrderWitness=lowerWitness,maximumPureStateACOrder=4,method='Independent exact-factorial Clebsch–Gordan contractions and exact-integer Racah swap coefficients with Arb enclosures; complete recursive CFP-character validation; rigorous projector/polar perturbation bounds; outward-rounded dyadic inverse-Cholesky congruence residuals including CFP error. No naked floating-point eigenvalue is used.',statement='There is no pure fifth-order anticoherent state on21 symmetric qubits. Consequently all orders at least5 are excluded; the explicit fourth-order witness proves the maximum is exactly4.',elapsedSeconds=time.time()-start,certificateMetadataSHA256=hashlib.sha256(raw).hexdigest())
 files=sorted([root/'certificate.json',root/'validated_bootstrap.py',root/'verify_certificate.py',root/'angular.py',root/'dual.npz',root/'congruences.npz']+[root/('cfp'+str(r)+'.npz')for r in range(2,7)]);report['inputSHA256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in files}
 (root/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print('EXCLUSION VERIFIED',float(minimum),report['elapsedSeconds']);return report
if __name__=='__main__':verify(Path(__file__).resolve().parent)
