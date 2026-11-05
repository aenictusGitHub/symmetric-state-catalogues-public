"""Exact covariance, Kramers averaging, vanishing, and invariant-dimension checks."""
from exact_field import *
from fractions import Fraction as F
from collections import Counter,defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb,lcm
import numpy as np,time
ROOT=Path(__file__).parent;start=time.time();cert=json.load(open(ROOT/'sos-certificate.json'));defn=json.load(open(ROOT/'covariants.json'));raw=json.load(open(ROOT/'source-basis.json'))['blocks'];holo=json.load(open(ROOT/'holomorphic-basis.json'))['blocks'];kernels=json.load(open(ROOT/'source-kernels.json'))['blocks'];N=11

def deriv(poly,plus):
 out=defaultdict(int)
 for m,c in poly.items():
  for i,count in Counter(m).items():
   r=i%12;conj=i//12
   if plus:factor=11-r if not conj else -r;to=i+1 if not conj else i-1
   else:factor=r if not conj else -(11-r);to=i-1 if not conj else i+1
   if factor:
    mm=list(m);mm.remove(i);mm.append(to);out[tuple(sorted(mm))]+=c*count*factor
 return {m:c for m,c in out.items()if c}

def right(poly):
 out=defaultdict(int)
 for m,c in poly.items():
  for i,count in Counter(m).items():
   r=i%12;other=(23-r)if i<12 else(11-r);mm=list(m);mm.remove(i);mm.append(other);out[tuple(sorted(mm))]+=c*count*(-1)**r
 return {m:c for m,c in out.items()if c}

for kind,blocks,ids,degree in [('source',raw,defn['sourceBlocks'],4),('holo',holo,range(len(holo)),5)]:
 for ib in ids:
  B0=blocks[ib];L=B0['L2']
  for cov in B0['covariants']:
   P=[{tuple(m):int(c)for m,c in p}for p in cov]
   for r,f in enumerate(P):
    assert all(len(m)==degree and all(0<=i<24 for i in m)and sum(i<12 for i in m)==int(B0.get('kind','50')[0])for m in f)
    assert deriv(f,True)==({m:(r+1)*v for m,v in P[r+1].items()}if r<L else {})
    assert deriv(f,False)==({m:(L-r+1)*v for m,v in P[r-1].items()}if r else {})
  print('COVARIANCE',kind,ib,flush=True)
# Character count for the SU(2) x SU(2) invariant part of Sym^10(V x C^2).
mult={}
for degree in [4,5,6]:
 counts=Counter(sum(m)for m in combinations_with_replacement(range(12),degree));mult[degree]={L:counts[(11*degree-L)//2]-counts[(11*degree-L)//2-1]for L in range((11*degree)%2,11*degree+1,2)}
a=sum(m*m for m in mult[5].values());b=sum(m*mult[4].get(L,0)for L,m in mult[6].items());assert(a,b,a-b)==(1794,1504,290)
print('PASS invariant dimension 1794-1504=290',flush=True)
@lru_cache(None)
def monvalue(m):
 ex=[0]*4;nt=nz=nw=sign=0
 for i in m:
  r=i%12
  if r not in (0,3,6,9,11,8,5,2):return None,K()
  other=r in (11,8,5,2);rr=11-r if other else r;ex[int(other)+2*(i//12)]+=1;nt+=rr==6;nz+=rr==9;nw+=rr in (3,9);sign+=(other and r%2)
 return tuple(ex)+(nw%2,),((-1)**sign)*T**nt*Z**nz*U**(nw//2)
def restriction(p):
 out={}
 for m,c in p:
  key,v=monvalue(tuple(m))
  if v:out[key]=out.get(key,K())+int(c)*v
 return {k:v for k,v in out.items()if v}
for ib in defn['sourceBlocks']:
 B0=raw[ib];cols=[[K(v)for v in col]for col in kernels[ib]['kernelColumns']]
 for r in range(B0['L2']+1):
  vals=[restriction(cov[r])for cov in B0['covariants']];keys=set().union(*(set(v)for v in vals))
  for col in cols:
   for key in keys:assert not sum((c*v.get(key,K())for c,v in zip(col,vals)),K())
 print('VANISHING',ib,flush=True)
# Exact holomorphic restrictions and first derivatives in T,U,Z at the attainer.
@lru_cache(None)
def jetmon(m):
 nt=nz=nw=sign=beta=0
 for r in m:
  if r not in (0,3,6,9,11,8,5,2):return None,None
  other=r in (11,8,5,2);rr=11-r if other else r;beta+=other;nt+=rr==6;nz+=rr==9;nw+=rr in(3,9);sign+=(other and r%2)
 v=(-1)**sign*T**nt*Z**nz*U**(nw//2);return(beta,nw%2),[v,v*nt/T,v*K(F(nw,2))/U,v*nz/Z]
restrs={};variables=[(ib,a,b)for ib,B0 in enumerate(holo)for b in range(B0['dimension'])for a in range(b+1)];chosen=[variables[i]for i in cert['invariantBasisColumns']];scales=cert['invariantBasisScales']
for ib in sorted(set(v[0]for v in chosen)):
 B0=holo[ib];out=[]
 for cov in B0['covariants']:
  cc=[]
  for p in cov:
   d={}
   for m,c in p:
    key,jet=jetmon(tuple(m))
    if key is None:continue
    old=d.setdefault(key,[K()for _ in range(4)])
    for j in range(4):old[j]+=int(c)*jet[j]
   cc.append(d)
  out.append(cc)
 restrs[ib]=out
jets=[]
for ib,a,b in chosen:
 L=holo[ib]['L2'];out=[[K()for _ in range(4)]for _ in range(2)]
 for r in range(L+1):
  aa=restrs[ib][a][r];bb=restrs[ib][b][r]
  for (t,p),v in aa.items():
   for pp in range(2):
    if(t,pp)not in bb:continue
    w=bb[t,pp];factor=K(F(1 if a==b else 2,6*comb(5,t)*comb(L,r))/F(scales[ib][a])/F(scales[ib][b]))*(U if p+pp==2 else 1);part=(p+pp)%2;out[part][0]+=factor*v[0]*w[0]
    for j in range(1,4):out[part][j]+=factor*(v[j]*w[0]+v[0]*w[j])
 assert not any(out[1]);jets.append([(v/Q**5).serial()for v in out[0]])
(ROOT/'exact-jets.json').write_text(json.dumps(jets,separators=(',',':'))+'\n');print('PASS exact invariant jets',time.time()-start,flush=True)
# Reconstruct the target directly from reduced density matrices.
def product(p,q):
 out=defaultdict(int)
 for m,c in p.items():
  for n,d in q.items():out[tuple(sorted(m+n))]+=c*d
 return {m:c for m,c in out.items()if c}
norm={(i,12+i):comb(11,i)for i in range(12)};q2=product(norm,norm);p6=defaultdict(int)
for t in range(7):
 wt=comb(11,6)*(-1)**(6-t)*comb(6,t)*2**t
 for a in range(t+1):
  for b in range(t+1):
   for i in range(12-t):
    for j in range(12-t):p6[tuple(sorted((a+i,12+b+i,b+j,12+a+j)))]+=wt*comb(t,a)*comb(t,b)*comb(11-t,i)*comb(11-t,j)
p6={m:c for m,c in p6.items()if c};assert not deriv(p6,True)and not deriv(p6,False)and not right(p6);assert not right(norm)
# Exact candidate value and stationary constraints, holding BOUND constant.
def candidate_jets(poly):
 out=[[K()for _ in range(4)]for _ in range(2)]
 for mon,c in poly.items():
  if any(i%12 not in(0,3,6,9)for i in mon):continue
  rr=[i%12 for i in mon];nt=rr.count(6);nz=rr.count(9);nw=rr.count(3)+nz;v=K(c)*T**nt*Z**nz*U**(nw//2);jet=[v,v*nt/T,v*K(F(nw,2))/U,v*nz/Z]
  for j in range(4):out[nw%2][j]+=jet[j]
 return out
Pj=candidate_jets(p6);Qj=candidate_jets(q2);assert all(Pj[p][j]==BOUND*Qj[p][j]for p in range(2)for j in range(4));assert K(cert['bound'])==BOUND
(ROOT/'target-polynomial.json').write_text(json.dumps([[list(m),str(c)]for m,c in p6.items()],separators=(',',':'))+'\n')
(ROOT/'exact-checks.json').write_text(json.dumps(dict(covarianceVerified=True,kernelVanishingVerified=True,rightInvarianceVerified=True,exactAttainmentVerified=True,exactStationarityVerified=True,invariantDimension=290,weightZeroMultiplicity=1794,weightTwoMultiplicity=1504),indent=2)+'\n');print('PASS exact target, symmetry, attainment and stationarity',time.time()-start,flush=True)
