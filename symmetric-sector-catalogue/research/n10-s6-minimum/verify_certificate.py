"""Global ten-qubit S6 minimum: independent algebraic/rigorous-interval verifier.
Requires python-flint (Arb/FLINT). No optimization, NumPy or SciPy is used.
All interval operations have outward rounding; the approximate solver mode is forbidden.
"""
import os;os.environ['OPENBLAS_NUM_THREADS']='1'
from exact_field import *
from flint import arb,arb_mat,ctx,fmpq_mat
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict
from math import comb,factorial
import time,hashlib
ctx.prec=768
BASE=Path(__file__).parent;start=time.time()
certificate=json.loads((BASE/'sos-certificate.json').read_text());raw=json.loads((BASE/'basis.json').read_text())['blocks'];kernels=json.loads((BASE/'kernels.json').read_text())['blocks']
assert certificate['n']==10 and certificate['sector']==6 and certificate['degree']==8
assert certificate['fieldPolynomial']==DATA['rootPolynomialCoefficients']
assert certificate['fieldRootInterval']==DATA['rootInterval']
assert K(certificate['bound'])==BOUND
scales=certificate['rawBasisScales'];variables=[(ib,a,b)for ib,B in enumerate(raw)for b in range(B['dimension'])for a in range(b+1)]
newvars=[tuple(v)for v in certificate['variables']];piv=certificate['pivotVariables'];x=[F(v)for v in certificate['freeCoefficients']];keys=[int(v)for v in certificate['coefficientMonomials']];rsel=list(range(126));layout={'blocks':certificate['blocks'],'norms':certificate['kernelScales']}
assert len(keys)==len(set(keys))==len(piv)==len(set(piv))==126
assert all(0<=p<len(x)and x[p]==0 for p in piv)
assert len(x)==len(newvars) and len(raw)==len(scales)==len(kernels)
assert layout['blocks']==[ib for ib,B in enumerate(kernels)if B['dimension']]
assert newvars==[(ib,a,b)for ib in layout['blocks']for b in range(kernels[ib]['dimension'])for a in range(b+1)]
# Verify SU(2) covariance, bidegree and weights directly from the polynomials.
from collections import Counter
from itertools import combinations_with_replacement
from math import lcm

def deriv(p,plus):
 out=defaultdict(int)
 for m,c in p.items():
  for i,count in Counter(m).items():
   r=i%11;conj=i//11
   if plus:factor=10-r if not conj else -r;to=i+1 if not conj else i-1
   else:factor=r if not conj else -(10-r);to=i-1 if not conj else i+1
   if factor:
    mm=list(m);mm.remove(i);mm.append(to);out[tuple(sorted(mm))]+=c*count*factor
 return {m:c for m,c in out.items()if c}
polyblocks=[]
counts=Counter(sum(m)for m in combinations_with_replacement(range(11),4));multiplicities={L:counts[20-L]-counts[19-L]for L in range(21)}
assert sum(m*(m+1)//2 for m in multiplicities.values())==130
seen={}
for ib,B in enumerate(raw):
 L=B['L2'];assert L%2==0 and B['weights']==[str(F(1,comb(L,r)))for r in range(L+1)]
 assert len(B['covariants'])==B['dimension']==len(scales[ib])
 pp=[]
 for cov in B['covariants']:
  assert len(cov)==L+1
  P=[{tuple(m):int(c)for m,c in comp}for comp in cov]
  for r,f in enumerate(P):
   assert all(len(m)==4 and all(0<=i<22 for i in m)and sum(i<11 for i in m)==int(B['kind'][0])for m in f)
   assert deriv(f,True)==({m:(r+1)*v for m,v in P[r+1].items()}if r<L else {})
   assert deriv(f,False)==({m:(L-r+1)*v for m,v in P[r-1].items()}if r else {})
  pp.append(P)
 polyblocks.append(pp)
 if B['kind']=='40':
  mons=sorted(set().union(*(set(c[0])for c in pp)));matrix=fmpq_mat([[p[0].get(m,0)for m in mons]for p in pp]);assert matrix.rank()==B['dimension']==multiplicities[L//2];seen[L//2]=B['dimension']
assert seen=={L:m for L,m in multiplicities.items()if m}
print('PASS exact covariance and complete 130-dimensional real invariant basis',flush=True)
# Compute the requested coefficient functionals exactly, directly from the basis.
split=9**11;encode=lambda m:sum(9**i for i in m);conjugate=lambda m:m//split+(m%split)*split
selected=set(keys);selected|={conjugate(m)for m in keys}
ip=[[[{encode(m):v for m,v in p.items()}for p in c]for c in block]for block in polyblocks]
features=[]
for iv,(ib,a,b)in enumerate(variables):
 out=defaultdict(int)
 if kernels[ib]['dimension'] or raw[ib]['kind']=='40':
  factors=[F(w)/(F(scales[ib][a])*F(scales[ib][b]))for w in raw[ib]['weights']];den=lcm(*(c.denominator for c in factors))
  for pa,pb,w in zip(ip[ib][a],ip[ib][b],factors):
   factor=int(w*den);conj=[(conjugate(m),c)for m,c in pb.items()]
   for m,c in pa.items():
    for mm,v in conj:
     key=m+mm
     if key not in selected:continue
     value=factor*c*v;out[key]+=value
     if a!=b:out[conjugate(key)]+=value
  features.append({m:F(v,den)for m,v in out.items()if v})
 else:features.append({})
 if iv%200==0:print('Exact coefficient feature',iv,'of',len(variables),flush=True)
lo,hi=[fmpq(x)for x in DATA['rootInterval']];assert MOD(lo)*MOD(hi)<0
for i in range(850):
 mid=(lo+hi)/2
 if not MOD(mid):lo=hi=mid;break
 if MOD(lo)*MOD(mid)<0:hi=mid
 else:lo=mid
ball=arb(lo).union(arb(hi));print('root',ball,flush=True)
def ev(k):
 v=arb(0)
 for c in reversed(K(k).p.coeffs()):v=v*ball+arb(c)
 return v
ba=ev(BOUND);assert ev(U)>0 and ev(Q)>0;print('bound',ba,flush=True)
def af(x):return arb(fmpq(str(x)))
# Uniqueness in the original isolating interval via a derivative of fixed sign.
orig=arb(fmpq(DATA['rootInterval'][0])).union(arb(fmpq(DATA['rootInterval'][1])));v=arb(0)
for c in reversed(MOD.derivative().coeffs()):v=v*orig+arb(c)
assert not v.contains(0)
# Exact value/derivative jets in the real support ansatz (t,u,z independent).
@lru_cache(None)
def jet(m):
 if any(i not in (0,3,6,9)for i in m):return (0,[K()for _ in range(4)])
 nw=m.count(3)+m.count(9);nt=m.count(6);nz=m.count(9);v=T**nt*Z**nz*U**(nw//2)
 return nw%2,[v,v*F(nw,2)/U,v*nt/T,v*nz/Z]
def restriction(comp):
 out=[[K()for _ in range(4)]for _ in range(2)]
 for m,c in comp:
  p,v=jet(tuple(sorted(i%11 for i in m)))
  if v[0]:
   for j in range(4):out[p][j]+=int(c)*v[j]
 return out
# Check every kernel column vanishes exactly, in Q(t), at the candidate.
restr={}
for ib,B in enumerate(raw):
 if not kernels[ib]['dimension'] and B['kind']!='40':continue
 restr[ib]=[[restriction(comp)for comp in cov]for cov in B['covariants']]
 for col in kernels[ib]['kernelColumns']:
  col=[K(v)for v in col]
  for r in range(B['L2']+1):
   for p in range(2):assert not sum((c*restr[ib][j][r][p][0]for j,c in enumerate(col)),K())
print('PASS all exact vanishing kernels',time.time()-start,flush=True)
# Evaluate the value and three derivatives of all 130 real holomorphic invariants.
invcols=[i for i,(ib,a,b)in enumerate(variables)if raw[ib]['kind']=='40'];assert len(invcols)==130 and invcols==certificate['invariantBasisColumns']
for ib,B in enumerate(raw):
 if B['kind']=='40':
  for r in range(B['L2']+1):
   parities={p for a in range(B['dimension'])for p in range(2)if any(restr[ib][a][r][p])}
   assert len(parities)<=1,'Mixed sqrt(u) parities in a holomorphic weight component'
jets=[]
for iv in invcols:
 ib,a,b=variables[iv];J=[K()for _ in range(4)]
 for r in range(raw[ib]['L2']+1):
  for p in range(2):
   A=restr[ib][a][r][p];B=restr[ib][b][r][p];factor=F(1 if a==b else 2,comb(raw[ib]['L2'],r))/F(scales[ib][a])/F(scales[ib][b]);fac=K(factor)*(U if p else 1)
   J[0]+=fac*A[0]*B[0]
   for j in range(1,4):J[j]+=fac*(A[j]*B[0]+A[0]*B[j])
 jets.append(J)
# The derivative jets already include derivative of sqrt(u); no extra d(U) term.
C=arb_mat([[af(features[j].get(keys[i],0))for j in invcols]for i in rsel]+[[ev(jet[j])for jet in jets]for j in range(4)])
det=C.det();print('injectivity determinant',det,flush=True);assert not det.contains(0)
# Full target jets vanish exactly; s_star is held constant when differentiating.
def product(p,q):
 out=defaultdict(F)
 for m,a in p.items():
  for n,b in q.items():out[tuple(sorted(m+n))]+=a*b
 return {m:c for m,c in out.items()if c}

def purity(t):
 out=defaultdict(F)
 for a in range(t+1):
  for b in range(t+1):
   for i in range(11-t):
    for j in range(11-t):out[tuple(sorted((a+i,11+b+i,b+j,11+a+j)))]+=comb(t,a)*comb(t,b)*comb(10-t,i)*comb(10-t,j)
 return out
norm={(i,11+i):F(comb(10,i))for i in range(11)};q2=product(norm,norm);q4=product(q2,q2);p6=defaultdict(F)
for t in range(7):
 for m,c in purity(t).items():p6[m]+=210*(-1)**(6-t)*comb(6,t)*2**t*c
rawtarget=product(q2,p6)
# Rational part of the target: q^2 P6 - L q^4, L=58437/800.
for m,c in q4.items():rawtarget[m]=rawtarget.get(m,F(0))-F('58437/800')*c
rawtarget={m:c for m,c in rawtarget.items()if c}
assert not deriv(rawtarget,True) and not deriv(rawtarget,False)
assert not deriv(q4,True) and not deriv(q4,False)
for p in [rawtarget,q4]:
 assert all(len(m)==8 and sum(i<11 for i in m)==4 for m in p)
 assert {tuple(sorted(i+11 if i<11 else i-11 for i in m)):c for m,c in p.items()}==p
J=[K()for _ in range(4)]
for m,c in rawtarget.items():
 p,v=jet(tuple(sorted(i%11 for i in m)))
 if any(v):assert p==0
 for j in range(4):J[j]+=K(c)*v[j]
targetdict={encode(m):c for m,c in rawtarget.items()}
qjets=[Q**4,4*Q**3*(120+10*Z*Z),4*Q**3*420*T,4*Q**3*20*U*Z]
assert all(not (J[j]-(BOUND-K(F('58437/800')))*qjets[j])for j in range(4))
print('PASS exact target vanishing and stationary constraints',flush=True)
def q4coef(key):
 if key%split!=key//split:return 0
 key=key%split;ex=[]
 for i in range(11):ex.append(key%9);key//=9
 from functools import reduce
 return F(factorial(4)*reduce(lambda a,b:a*b,[comb(10,i)**e for i,e in enumerate(ex)],1),reduce(lambda a,b:a*b,[factorial(e)for e in ex],1))
rhs=arb_mat([[af(targetdict.get(keys[i],0))-(ba-af(F('58437/800')))*af(q4coef(keys[i]))]for i in rsel]);system=arb_mat(126,126)
for bn,ib in enumerate(layout['blocks']):
 d=raw[ib]['dimension'];e=kernels[ib]['dimension'];R=arb_mat([[ev(K(kernels[ib]['kernelColumns'][b][a]))*af(F(scales[ib][a]))/af(F(layout['norms'][bn][b]))for b in range(e)]for a in range(d)])
 H=arb_mat([[af(x[newvars.index((ib,min(a,b),max(a,b)))])for b in range(e)]for a in range(e)]);G=R*H*R.transpose();rv=[i for i,(block,a,b)in enumerate(variables)if block==ib];pv=[(j,newvars[p])for j,p in enumerate(piv)if newvars[p][0]==ib]
 lift=arb_mat(len(rv),len(pv)+1)
 for i,v in enumerate(rv):
  _,a,b=variables[v];lift[i,0]=G[a,b]
  for j,(_,(_,c,f))in enumerate(pv):lift[i,j+1]=R[a,c]*R[b,f]+(R[a,f]*R[b,c]if c!=f else 0)
 A=arb_mat([[af(features[v].get(keys[i],0))for v in rv]for i in rsel]);result=A*lift
 for i in range(126):
  rhs[i,0]-=result[i,0]
  for j,(p,_)in enumerate(pv):system[i,p]=result[i,j+1]
print('built exact interval system',time.time()-start,flush=True)
solution=system.solve(rhs,algorithm='precond');xx=[af(v)for v in x]
for i,p in enumerate(piv):xx[p]=solution[i,0]
mins=[]
for ib in layout['blocks']:
 e=kernels[ib]['dimension'];H=[[xx[newvars.index((ib,min(a,b),max(a,b)))]for b in range(e)]for a in range(e)]
 for j in range(e):
  v=H[j][j];assert v>arb(fmpq(1,25)),(ib,j,v);mins.append(float(v.lower()))
  for a in range(j+1,e):
   for b in range(j+1,e):H[a][b]-=H[a][j]*H[j][b]/v
 print('PASS positive block',ib,'size',e,flush=True)
print('SUCCESS positivity',min(mins),'seconds',time.time()-start,flush=True)
report={'n':10,'sector':6,'direction':'min','degree':8,'globalOptimalityProven':True,'fullPolynomialIdentityVerified':True,'identityVerification':'130-dimensional invariant-space completeness; 126 exact coefficient equations and four exact vanishing constraints; certified injectivity','attainmentVerifiedExact':True,'arithmetic':'Exact rational and degree-14 algebraic arithmetic; directed-rounding Arb real-ball nonsingularity and positivity certificates','floatingPointToleranceUsed':False,'invariantDimension':130,'coefficientConstraints':126,'vanishingConstraints':4,'positiveDefiniteGramBlocks':len(layout['blocks']),'squaredForms':sum(kernels[ib]['dimension']*(raw[ib]['L2']+1)for ib in layout['blocks']),'LDLPivotStrictLowerBound':'1/25','precisionBits':ctx.prec,'boundEnclosure':str(ba),'inputSha256':{name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()for name in ['basis.json','kernels.json','sos-certificate.json','state-parameters.json','exact_field.py','verify_certificate.py']}}
(BASE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
