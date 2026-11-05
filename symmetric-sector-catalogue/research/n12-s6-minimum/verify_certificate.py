"""Independent exact-algebraic and directed-rounding proof of the N=12 sector minimum.
Requires python-flint. No optimizer or floating-point tolerance is used.
"""
from exact_field import *
from flint import arb,arb_mat,ctx,fmpq_mat
from fractions import Fraction as F
from functools import lru_cache,reduce
from collections import defaultdict
from math import comb,factorial,lcm
import time,hashlib
ctx.prec=1024;BASE=Path(__file__).resolve().parent;start=time.time();k=DATA['k'];rank=246
certificate=json.loads((BASE/'sos-certificate.json').read_text());raw=json.loads((BASE/'basis.json').read_text())['blocks'];ker=json.loads((BASE/'kernels.json').read_text())['blocks'];kernels=ker
assert certificate['n']==12 and certificate['sector']==k and k in [6,8] and certificate['degree']==8
assert certificate['fieldPolynomial']==DATA['rootPolynomialCoefficients'] and certificate['fieldRootInterval']==DATA['rootInterval'] and K(certificate['bound'])==BOUND
scales=certificate['rawBasisScales'];variables=[(ib,a,b)for ib,B in enumerate(raw)for b in range(B['dimension'])for a in range(b+1)];newvars=[tuple(v)for v in certificate['variables']];piv=certificate['pivotVariables'];xx=[F(v)for v in certificate['freeCoefficients']];keys=list(map(int,certificate['coefficientMonomials']));rsel=list(range(rank));layout={'blocks':certificate['blocks']};whites=certificate['kernelMixingMatrices']
assert len(keys)==len(set(keys))==len(piv)==len(set(piv))==rank
assert all(0<=p<len(xx)and xx[p]==0 for p in piv)
assert len(xx)==len(newvars)and len(raw)==len(scales)==len(ker)
assert layout['blocks']==[ib for ib,B in enumerate(ker)if B['dimension']]
assert newvars==[(ib,a,b)for ib in layout['blocks']for b in range(ker[ib]['dimension'])for a in range(b+1)]
# Verify SU(2) covariance, bidegree and weights directly from the polynomials.
from collections import Counter
from itertools import combinations_with_replacement
from math import lcm

def deriv(p,plus):
 out=defaultdict(int)
 for m,c in p.items():
  for i,count in Counter(m).items():
   r=i%13;conj=i//13
   if plus:factor=12-r if not conj else -r;to=i+1 if not conj else i-1
   else:factor=r if not conj else -(12-r);to=i-1 if not conj else i+1
   if factor:
    mm=list(m);mm.remove(i);mm.append(to);out[tuple(sorted(mm))]+=c*count*factor
 return {m:c for m,c in out.items()if c}
polyblocks=[]
counts=Counter(sum(m)for m in combinations_with_replacement(range(13),4));multiplicities={L:counts[24-L]-counts[23-L]for L in range(25)}
assert sum(m*(m+1)//2 for m in multiplicities.values())==250
seen={}
for ib,B in enumerate(raw):
 L=B['L2'];assert L%2==0 and B['weights']==[str(F(1,comb(L,r)))for r in range(L+1)]
 assert len(B['covariants'])==B['dimension']==len(scales[ib])
 pp=[]
 for cov in B['covariants']:
  assert len(cov)==L+1
  P=[{tuple(m):int(c)for m,c in comp}for comp in cov]
  for r,f in enumerate(P):
   assert all(len(m)==4 and all(0<=i<26 for i in m)and sum(i<13 for i in m)==int(B['kind'][0])for m in f)
   assert deriv(f,True)==({m:(r+1)*v for m,v in P[r+1].items()}if r<L else {})
   assert deriv(f,False)==({m:(L-r+1)*v for m,v in P[r-1].items()}if r else {})
  pp.append(P)
 polyblocks.append(pp)
 if B['kind']=='40':
  mons=sorted(set().union(*(set(c[0])for c in pp)));matrix=fmpq_mat([[p[0].get(m,0)for m in mons]for p in pp]);assert matrix.rank()==B['dimension']==multiplicities[L//2];seen[L//2]=B['dimension']
assert seen=={L:m for L,m in multiplicities.items()if m}
print('PASS exact covariance and complete 250-dimensional real invariant basis',flush=True)
# Compute the requested coefficient functionals exactly, directly from the basis.
split=9**13;encode=lambda m:sum(9**i for i in m);conjugate=lambda m:m//split+(m%split)*split
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
for _ in range(1150):
 mid=(lo+hi)/2
 if MOD(lo)*MOD(mid)<0:hi=mid
 else:lo=mid
ball=arb(lo).union(arb(hi))
def ev(z):
 v=arb(0)
 for c in reversed(K(z).p.coeffs()):v=v*ball+arb(c)
 return v
def af(x):return arb(fmpq(str(x)))
print('bound',ev(BOUND),flush=True)
orig=arb(fmpq(DATA['rootInterval'][0])).union(arb(fmpq(DATA['rootInterval'][1])));v=arb(0)
for c in reversed(MOD.derivative().coeffs()):v=v*orig+arb(c)
assert not v.contains(0) and ev(U)>0 and ev(Q)>0
# Exact candidate value/derivative jets in the three free ansatz parameters.
if k==6:JM={0:(165,1,0,0),2:(1,0,1,1),4:(-1,0,0,0),6:(F(1,14),0,0,1),8:(1,0,0,0),10:(1,0,1,1),12:(-165,1,0,0)}
else:JM={1:(1,0,1,1),4:(1,1,0,0),7:(1,0,0,1),10:(1,0,0,0)}
@lru_cache(None)
def jet(m):
 c=F(1);nx=ny=nu=0
 for i in m:
  if i not in JM:return 0,[K()for _ in range(4)]
  co,a,b,u=JM[i];c*=co;nx+=a;ny+=b;nu+=u
 v=K(c)*X**nx*Y**ny*U**(nu//2);return nu%2,[v,v*nx/X,v*ny/Y,v*F(nu,2)/U]
def restr(poly):
 out=[[K()for _ in range(4)]for _ in range(2)]
 for mon,c in poly:
  p,v=jet(tuple(sorted(i%13 for i in mon)))
  for j in range(4):out[p][j]+=int(c)*v[j]
 return out
restrs={}
for ib,B in enumerate(raw):
 if not ker[ib]['dimension']and B['kind']!='40':continue
 restrs[ib]=[[restr(p)for p in cov]for cov in B['covariants']]
 for col in ker[ib]['kernelColumns']:
  col=[K(v)for v in col]
  for r in range(B['L2']+1):
   for p in range(2):assert not sum((c*restrs[ib][j][r][p][0]for j,c in enumerate(col)),K())
print('PASS every exact vanishing kernel',flush=True)
invcols=[i for i,(ib,a,b)in enumerate(variables)if raw[ib]['kind']=='40'];assert len(invcols)==250 and invcols==certificate['invariantBasisColumns'];jets=[]
for iv in invcols:
 ib,a,b=variables[iv];J=[K()for _ in range(4)]
 for r in range(raw[ib]['L2']+1):
  assert len({p for cc in range(raw[ib]['dimension'])for p in range(2)if any(restrs[ib][cc][r][p])})<=1
  for p in range(2):
   A=restrs[ib][a][r][p];B=restrs[ib][b][r][p];fac=K(F(1 if a==b else 2,comb(raw[ib]['L2'],r))/F(scales[ib][a])/F(scales[ib][b]))*(U if p else 1);J[0]+=fac*A[0]*B[0]
   for j in range(1,4):J[j]+=fac*(A[j]*B[0]+A[0]*B[j])
 jets.append(J)
C=arb_mat([[af(features[j].get(keys[i],0))for j in invcols]for i in rsel]+[[ev(jet[j])for jet in jets]for j in range(4)]);det=C.det();assert not det.contains(0);print('PASS injectivity',str(det),time.time()-start,flush=True)
def product(p,q):
 out=defaultdict(int)
 for m,c in p.items():
  for n,d in q.items():out[tuple(sorted(m+n))]+=c*d
 return {m:c for m,c in out.items()if c}
norm={(i,13+i):comb(12,i)for i in range(13)};q2=product(norm,norm);q4=product(q2,q2);pk=defaultdict(int)
for t in range(k+1):
 wt=comb(12,k)*(-1)**(k-t)*comb(k,t)*2**t
 for a in range(t+1):
  for b in range(t+1):
   for i in range(13-t):
    for j in range(13-t):pk[tuple(sorted((a+i,13+b+i,b+j,13+a+j)))]+=wt*comb(t,a)*comb(t,b)*comb(12-t,i)*comb(12-t,j)
target=product(q2,pk);assert not deriv(target,True)and not deriv(target,False)and not deriv(q4,True)and not deriv(q4,False)
for poly in [target,q4]:
 assert all(len(m)==8 and sum(i<13 for i in m)==4 for m in poly)
 assert {tuple(sorted(i+13 if i<13 else i-13 for i in m)):c for m,c in poly.items()}==poly
J=[K()for _ in range(4)]
for mon in set(target)|set(q4):
 p,v=jet(tuple(sorted(i%13 for i in mon)));co=K(target.get(mon,0))-BOUND*q4.get(mon,0)
 if any(v):assert p==0
 for j in range(4):J[j]+=co*v[j]
assert not any(J);print('PASS exact target jets',time.time()-start,flush=True)
encode=lambda mon:sum(9**i for i in mon);target={encode(m):c for m,c in target.items()};q4={encode(m):c for m,c in q4.items()};rhs=arb_mat([[af(target.get(keys[i],0))-ev(BOUND)*af(q4.get(keys[i],0))]for i in rsel]);system=arb_mat(rank,rank)
for bn,ib in enumerate(layout['blocks']):
 d=raw[ib]['dimension'];e=ker[ib]['dimension'];R=arb_mat([[ev(K(ker[ib]['kernelColumns'][b][a]))*af(F(scales[ib][a]))for b in range(e)]for a in range(d)])*arb_mat([[af(F(v))for v in row]for row in whites[bn]]);H=arb_mat([[af(xx[newvars.index((ib,min(a,b),max(a,b)))])for b in range(e)]for a in range(e)]);G=R*H*R.transpose();rv=[i for i,v in enumerate(variables)if v[0]==ib];pv=[(j,newvars[p])for j,p in enumerate(piv)if newvars[p][0]==ib];lift=arb_mat(len(rv),len(pv)+1)
 for i,iv in enumerate(rv):
  _,a,b=variables[iv];lift[i,0]=G[a,b]
  for j,(_,(_,c,f))in enumerate(pv):lift[i,j+1]=R[a,c]*R[b,f]+(R[a,f]*R[b,c]if c!=f else 0)
 A=arb_mat([[af(features[v].get(keys[i],0))for v in rv]for i in rsel]);result=A*lift
 for i in range(rank):
  rhs[i,0]-=result[i,0]
  for j,(p,_)in enumerate(pv):system[i,p]=result[i,j+1]
 print('system block',ib,flush=True)
solution=system.solve(rhs,algorithm='precond');zz=[af(v)for v in xx]
for i,p in enumerate(piv):zz[p]=solution[i,0]
mins=[]
for ib in layout['blocks']:
 e=ker[ib]['dimension'];H=[[zz[newvars.index((ib,min(a,b),max(a,b)))]for b in range(e)]for a in range(e)]
 for j in range(e):
  v=H[j][j];assert v>0,(ib,j,v);mins.append(float(v.lower()))
  for a in range(j+1,e):
   for b in range(j+1,e):H[a][b]-=H[a][j]*H[j][b]/v
 print('positive',ib,flush=True)
print('SUCCESS',k,min(mins),time.time()-start,flush=True)
report={'n':12,'sector':k,'direction':'min','degree':8,'globalOptimalityProven':True,'fullPolynomialIdentityVerified':True,'identityVerification':'250-dimensional invariant-space completeness; 246 exact coefficient equations and four exact vanishing constraints; certified injectivity','attainmentVerifiedExact':True,'arithmetic':'Exact rational/algebraic arithmetic and directed-rounding Arb real-ball certificates','floatingPointToleranceUsed':False,'invariantDimension':250,'coefficientConstraints':246,'vanishingConstraints':4,'positiveDefiniteGramBlocks':len(layout['blocks']),'squaredForms':sum(ker[ib]['dimension']*(raw[ib]['L2']+1)for ib in layout['blocks']),'positiveLDLPivotLowerBound':min(mins),'precisionBits':ctx.prec,'boundEnclosure':str(ev(BOUND)),'inputSha256':{name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()for name in ['basis.json','kernels.json','sos-certificate.json','state-parameters.json','exact_field.py','verify_certificate.py']}}
(BASE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
