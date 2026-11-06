"""Validated CFP bootstrap. Every certificate float is an exact binary rational.
All nonrational coefficients and every residual are enclosed by Arb balls.
No floating-point residual is accepted as a proof bound.
"""
from pathlib import Path
from math import factorial,comb
from functools import lru_cache
import json,time,numpy as np
import flint
flint.ctx.prec=128;flint.ctx.threads=1
A=flint.arb;M=flint.arb_mat;Q=flint.fmpq
F=[factorial(k)for k in range(512)]
@lru_cache(maxsize=None)
def sixj_ball(a,b,c,d,e,f):
 triangles=[(a,b,c),(a,e,f),(d,b,f),(d,e,c)]
 if any((x+y+z)%2 or z>x+y or z<abs(x-y)for x,y,z in triangles):return A(0)
 aa=[(a+b+c)//2,(a+e+f)//2,(d+b+f)//2,(d+e+c)//2];bb=[(a+b+d+e)//2,(a+c+d+f)//2,(b+c+e+f)//2];low=max(aa);high=min(bb);num=den=1
 for x,y,z in triangles:num*=F[(x+y-z)//2]*F[(x-y+z)//2]*F[(-x+y+z)//2];den*=F[(x+y+z)//2+1]
 total=0
 for z in range(low,high+1):
  dd=1
  for q in aa:dd*=F[z-q]
  for q in bb:dd*=F[q-z]
  assert F[z+1]%dd==0;total+=(-1 if z%2 else 1)*(F[z+1]//dd)
 if not total:return A(0)
 value=A(Q(total*total*num,den)).sqrt();return value if total>0 else-value
@lru_cache(maxsize=None)
def swap_ball(n,L,S,J,K):return (-1 if(n+(S+K)//2)%2 else 1)*A((S+1)*(K+1)).sqrt()*sixj_ball(n,L,S,n,J,K)
def square_frob_upper(X):
 s=A(0)
 for v in X.entries():q=v.abs_upper();s+=q*q
 return s.upper()
def frob_upper(X):return square_frob_upper(X).sqrt().upper()
def bound(x):return dict(upperRational=str(x.upper().fmpq()),diagnosticApproximation=float(x.upper()))
def multiplicities(n,r):
 dp=[[0]*(n*r+1)for _ in range(r+1)];dp[0][0]=1
 for v in range(n+1):
  for k in range(1,r+1):
   for s in range(v,n*r+1):dp[k][s]+=dp[k-1][s-v]
 out={}
 for J in range(n*r,n*r%2-1,-2):
  p=dp[r][(J+n*r)//2]-(dp[r][(J+2+n*r)//2]if J+2<=n*r else 0)
  if p:out[J]=p
 assert sum((J+1)*p for J,p in out.items())==comb(n+r,r)
 return out
def archive_names(spec):
 return [spec] if isinstance(spec,str) else list(spec)
def load_arrays(root,spec):
 names=archive_names(spec);assert names and len(names)==len(set(names));result={}
 for name in names:
  assert isinstance(name,str) and Path(name).name==name
  with np.load(root/name) as z:
   for key in z.files:
    assert key not in result
    result[key]=z[key]
 return result
def load_families(root,cert):
 n=cert['n'];rmax=cert['copies'];fam={}
 for r in range(2,rmax+1):
  d=cert['stages'][str(r)];z=load_arrays(root,d['file']);lower=[dict(J=int(S),a=a,b=b,p=b-a)for S,(a,b)in d['lower'].items()];lower.sort(key=lambda b:b['a']);expected=multiplicities(n,r);assert{b['J']:b['p']for b in d['source']}==expected;assert[b['J']for b in d['source']]==list(expected);off=0
  for b in d['source']:assert(b['a'],b['b'])==(off,off+b['p']);off=b['b']
  if r==2:assert d['lower']=={str(n):[0,1]}
  else:assert d['lower']=={str(b['J']):[b['a'],b['b']]for b in cert['stages'][str(r-1)]['source']}
  arr={};compact={}
  for b in d['source']:
   valid=[p for p in lower if abs(n-p['J'])<=b['J']<=n+p['J']];ids=[i for p in valid for i in range(p['a'],p['b'])];v=z['J'+str(b['J'])];assert v.shape==(len(ids),b['p'])and np.all(np.isfinite(v));full=np.zeros((lower[-1]['b'],b['p']));full[ids]=v;arr[b['J']]=full;compact[b['J']]=v
  fam[r]=dict(meta=d,lower=lower,full=arr,compact=compact)
 for b in fam[2]['meta']['source']:
  q=fam[2]['compact'][b['J']];assert q.shape==(1,1)and q[0,0]in[-1,1]
 return fam
def bootstrap(root,cert=None):
 start=time.time();cert=cert or json.loads((root/'certificate.json').read_text());n=cert['n'];rmax=cert['copies'];fam=load_families(root,cert);delta=A(0);reports=[]
 for r in range(3,rmax+1):
  parents=fam[r-1]['meta']['source'];current=fam[r]['meta']['source'];lower=fam[r-1]['lower'];prev={}
  for B in parents:
   S=B['J'];full=fam[r-1]['full'][S]
   for L in lower:
    if abs(n-L['J'])<=S<=n+L['J']:prev[S,L['J']]=M(full[L['a']:L['b']].tolist())
  offsets={};cur={};res={}
  for B in current:
   J=B['J'];valid=[P for P in parents if abs(n-P['J'])<=J<=n+P['J']];off=0;offsets[J]={};cur[J]={};res[J]={}
   for P in valid:
    S=P['J'];v=M(fam[r]['compact'][J][off:off+P['p']].tolist());cur[J][S]=v;res[J][S]=-v;offsets[J][S]=(off,off+P['p']);off+=P['p']
  count=0
  for si,B in enumerate(parents):
   S=B['J']
   for C in parents[si:]:
    K=C['J'];common=[L['J']for L in lower if(S,L['J'])in prev and(K,L['J'])in prev]
    Js=[b['J']for b in current if S in cur[b['J']]and K in cur[b['J']]]
    if not common or not Js:continue
    products=[(L,prev[S,L].transpose()*prev[K,L])for L in common]
    for J in Js:
     block=M(B['p'],C['p'])
     for L,product in products:
      scalar=swap_ball(n,L,S,J,K)
      if not scalar.is_zero():block+=product*scalar
     res[J][S]+=block*cur[J][K]
     if S!=K:res[J][K]+=block.transpose()*cur[J][S]
    count+=1
   if si%8==0:print('bootstrap stage',r,'parent',si,'seconds',time.time()-start,flush=True)
  beta=(A(r-1)/r*(2*delta+delta*delta)).upper();stage=[];newdelta=A(0)
  for B in current:
   J=B['J'];v=M(fam[r]['compact'][J].tolist());gram=v.transpose()*v-M(np.eye(B['p']).tolist());zeta=frob_upper(gram);ss=A(0)
   for residual in res[J].values():ss+=square_frob_upper(residual)
   raw=((A(r-1)/r)*ss.sqrt()).upper();eps=(raw+beta*(1+zeta).sqrt()).upper();eta=(zeta+eps*eps).upper();assert eta<1,(r,J,str(eta));distance=(eps+eta/(1+(1-eta).sqrt())).upper()
   if distance>newdelta:newdelta=distance
   stage.append(dict(J=J,multiplicity=B['p'],pointProjectorResidual=bound(raw),pointGramResidual=bound(zeta),exactRangeResidual=bound(eps),polarDistanceBound=bound(distance)))
   del res[J],cur[J]
  reports.append(dict(copies=r,previousIsometryDistance=bound(delta),projectorPerturbation=bound(beta),maximumIsometryDistance=bound(newdelta),checks=stage));delta=newdelta;print('STAGE VERIFIED',r,'distance upper',float(delta),'seconds',time.time()-start,flush=True)
  del prev,cur,res
 report=dict(passed=True,n=n,copies=rmax,precisionBits=flint.ctx.prec,maximumExactInputCFPDistance=bound(delta),stages=reports,seconds=time.time()-start,method='Exact dyadic CFP inputs; exact-integer Racah coefficients enclosed by Arb; validated projector and Gram residuals; spectral-gap-one projection and polar-normalization perturbation bounds. No floating residual is used as a proof bound.')
 (root/'bootstrap-verification.json').write_text(json.dumps(report,indent=2)+'\n');return delta,report
if __name__=='__main__':
 root=Path(__file__).resolve().parent;delta,report=bootstrap(root);print('BOOTSTRAP PASS',float(delta),report['seconds'])
