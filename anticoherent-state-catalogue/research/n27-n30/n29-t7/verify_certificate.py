"""Validated pure-state AC exclusion: exact dyadic inputs and Arb bounds."""
from pathlib import Path
from math import comb
import json,time,hashlib,numpy as np,flint
from validated_bootstrap import A,M,Q,frob_upper,square_frob_upper,bound,load_families,bootstrap
flint.ctx.prec=128;flint.ctx.threads=1
from angular import direct_coefficient
from ppt_term import PPTTerm
PPT_CONTEXT=None

def prepare(root,cert):
 assert cert['schema']=='pure-ac-moment-validated-separation-v1';n,t,r=cert['n'],cert['t'],cert['copies'];assert n>0 and 1<=t<=n and 2<=r<=6
 fam=load_families(root,cert);parent=fam[r-1]['meta']['source'];dual=np.load(root/cert['dualFile']);dualM=[];seen=set()
 for i,P in enumerate(cert['pairs']):
  si,ti=P['si'],P['ti'];assert 0<=si<=ti<len(parent);S,T=parent[si],parent[ti];assert(P['S'],P['T'],P['ps'],P['pt'])==(S['J'],T['J'],S['p'],T['p']);assert(si,ti)not in seen;seen.add((si,ti));v=dual['R'+str(i)];assert v.shape==(t,S['p'],T['p'])and np.all(np.isfinite(v))
  if si==ti:assert np.array_equal(v,v.swapaxes(1,2))
  dualM.append([M(v[ell].tolist())for ell in range(t)])
 return fam,parent,dualM

def source_matrix(B,fam,parent,pairs,dualM,n,t,r):
 J=B['J'];valid=[S for S in parent if abs(n-S['J'])<=J<=n+S['J']];offset={};off=0
 for S in valid:offset[S['J']]=(off,off+S['p']);off+=S['p']
 H=M(off,off);norm2=A(0)
 for P,D in zip(pairs,dualM):
  S,T=P['S'],P['T']
  if S not in offset or T not in offset:continue
  block=M(P['ps'],P['pt'])
  for ell in range(1,t+1):
   c=direct_coefficient(n,J,S,T,ell)
   if not c.is_zero():block+=D[ell-1]*c
  aa,bb=offset[S];cc,dd=offset[T]
  for a in range(P['ps']):
   for b in range(P['pt']):
    H[aa+a,cc+b]=block[a,b]
    if S!=T:H[cc+b,aa+a]=block[a,b]
  norm2+=(1 if S==T else 2)*square_frob_upper(block)
 extra=A(0)
 if PPT_CONTEXT is not None:
  hp,knorm=PPT_CONTEXT.ambient(B,parent);H-=hp;extra=knorm;hnorm=frob_upper(H)
 else:hnorm=norm2.sqrt().upper()
 V=M(fam[r]['compact'][J].tolist());return V.transpose()*H*V,hnorm,extra

def verify(root):
 start=time.time();raw=(root/'certificate.json').read_bytes();cert=json.loads(raw);n,t,r=cert['n'],cert['t'],cert['copies'];delta,boot=bootstrap(root,cert);fam,parent,D=prepare(root,cert);global PPT_CONTEXT;PPT_CONTEXT=PPTTerm(root,cert,fam);d3=A(Q(next(z for z in boot['stages']if z['copies']==3)['maximumIsometryDistance']['upperRational']));factors=np.load(root/'congruences.npz');checks=[];minimum=A(1)
 for B in fam[r]['meta']['source']:
  J=B['J'];G,hnorm,knorm=source_matrix(B,fam,parent,cert['pairs'],D,n,t,r);v=factors['J'+str(J)];assert v.shape==(B['p'],B['p'])and np.all(np.isfinite(v));T=M(v.tolist());rho=frob_upper(M(np.eye(B['p']).tolist())-T*G*T.transpose());tnorm2=square_frob_upper(T);gramError=((2*delta+delta*delta)*hnorm+(2*d3+d3*d3)*knorm).upper();total=(rho+tnorm2*gramError).upper();assert total<1,(J,str(total));gap=((1-total)/tnorm2).lower();assert gap>0
  if gap<minimum:minimum=gap
  checks.append(dict(J=J,dimension=B['p'],congruenceResidual=bound(rho),ambientDualFrobeniusNorm=bound(hnorm),CFPInducedGramPerturbation=bound(gramError),partialTransposeAmbientDualNorm=bound(knorm),thirdCopyCFPDistance=bound(d3),inverseFactorFrobeniusNormSquared=bound(tnorm2),totalCongruenceError=bound(total),positiveEigenvalueLowerBound=dict(lowerRational=str(gap.fmpq()),diagnosticApproximation=float(gap))))
  print('validated source',n,t,J,B['p'],float(total),time.time()-start,flush=True)
 report=dict(passed=True,n=n,t=t,copies=r,inputMomentDimension=comb(n+r,r),spinBlocks=len(checks),inputHighestWeightColumns=sum(b['p']for b in fam[r]['meta']['source']),realInvariantVariables=sum(b['p']*(b['p']+1)//2 for b in fam[r]['meta']['source']),allSourceBlocksStrictlyPositive=True,partialTimeReversalDualPositiveByGramFactors=True,pptDualFactorColumns={str(o['L']):int(np.load(root/cert['ppt']['factorFile'])['L'+str(o['L'])].shape[1])for o in cert['ppt']['outputs']},pptOutputSpins=[o['L']for o in cert['ppt']['outputs']],minimumCertifiedSourceEigenvalue=dict(lowerRational=str(minimum.fmpq()),diagnosticApproximation=float(minimum)),CFPBootstrap=boot,positiveSourceBlocks=checks,precisionBits=flint.ctx.prec,arbitraryComplexPureStatesCovered=True,statement=f'No pure order-{t} anticoherent state on {n} symmetric qubits exists. Therefore the maximum pure-state AC order is at most {t-1}.',elapsedSeconds=time.time()-start,certificateMetadataSHA256=hashlib.sha256(raw).hexdigest())
 files=sorted([root/f for f in ['certificate.json','validated_bootstrap.py','verify_certificate.py','angular.py','ppt_term.py','ppt-factors.npz','dual.npz','congruences.npz']]+[root/('cfp'+str(k)+'.npz')for k in range(2,r+1)]);report['inputSHA256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in files};(root/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print('EXCLUSION VERIFIED',n,t,float(minimum),report['elapsedSeconds']);return report
if __name__=='__main__':verify(Path(__file__).resolve().parent)
