"""Independent rational proof of the refined eleven-qubit S6 lower bound.
Exact coefficient verification on the full SU(2)-invariant polynomial space.
Requires python-flint; no floating-point arithmetic or optimizer is used.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
from math import comb,factorial,lcm
import json,hashlib,sys
from flint import fmpq,fmpq_mat
sys.set_int_max_str_digits(0)
BASE=Path(__file__).resolve().parent
basis_file=BASE/'basis.json';certificate_file=BASE/'sos-certificate.json'
raw=json.loads(basis_file.read_text())['blocks'];certificate=json.loads(certificate_file.read_text())
assert certificate['n']==11 and certificate['sector']==6 and certificate['target']=='min'
scales=certificate['basisScales'];keys=list(map(int,certificate['coefficientMonomials']));assert len(keys)==len(set(keys))==177
variables=[(ib,a,b)for ib,B in enumerate(raw)for b in range(B['dimension'])for a in range(b+1)]
kernels=[{'dimension':B['dimension']}for B in raw]
# Verify SU(2) covariance, bidegree and weights directly from the polynomials.
from collections import Counter
from itertools import combinations_with_replacement
from math import lcm

def deriv(p,plus):
 out=defaultdict(int)
 for m,c in p.items():
  for i,count in Counter(m).items():
   r=i%12;conj=i//12
   if plus:factor=11-r if not conj else -r;to=i+1 if not conj else i-1
   else:factor=r if not conj else -(11-r);to=i-1 if not conj else i+1
   if factor:
    mm=list(m);mm.remove(i);mm.append(to);out[tuple(sorted(mm))]+=c*count*factor
 return {m:c for m,c in out.items()if c}
polyblocks=[]
counts=Counter(sum(m)for m in combinations_with_replacement(range(12),4));multiplicities={L:counts[22-L]-counts[21-L]for L in range(23)}
assert sum(m*(m+1)//2 for m in multiplicities.values())==177
seen={}
for ib,B in enumerate(raw):
 L=B['L2'];assert L%2==0 and B['weights']==[str(F(1,comb(L,r)))for r in range(L+1)]
 assert len(B['covariants'])==B['dimension']==len(scales[ib])
 pp=[]
 for cov in B['covariants']:
  assert len(cov)==L+1
  P=[{tuple(m):int(c)for m,c in comp}for comp in cov]
  for r,f in enumerate(P):
   assert all(len(m)==4 and all(0<=i<24 for i in m)and sum(i<12 for i in m)==int(B['kind'][0])for m in f)
   assert deriv(f,True)==({m:(r+1)*v for m,v in P[r+1].items()}if r<L else {})
   assert deriv(f,False)==({m:(L-r+1)*v for m,v in P[r-1].items()}if r else {})
  pp.append(P)
 polyblocks.append(pp)
 if B['kind']=='40':
  mons=sorted(set().union(*(set(c[0])for c in pp)));matrix=fmpq_mat([[p[0].get(m,0)for m in mons]for p in pp]);assert matrix.rank()==B['dimension']==multiplicities[L//2];seen[L//2]=B['dimension']
assert seen=={L:m for L,m in multiplicities.items()if m}
print('PASS exact covariance and complete 177-dimensional real invariant basis',flush=True)
# Compute the requested coefficient functionals exactly, directly from the basis.
split=9**12;encode=lambda m:sum(9**i for i in m);conjugate=lambda m:m//split+(m%split)*split
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
invcols=[i for i,(ib,a,b)in enumerate(variables)if raw[ib]['kind']=='40'];assert len(invcols)==177
M=fmpq_mat([[fmpq(str(features[j].get(key,0)))for j in invcols]for key in keys]);assert M.det()!=0
print('PASS exact coefficient functionals form a basis on the complete invariant space',flush=True)
result=fmpq_mat(177,1);positivity=[]
assert len(raw)==len(certificate['blocks'])==len(certificate['basisTransforms'])
for ib,B in enumerate(raw):
 d=B['dimension'];block=certificate['blocks'][ib];assert block['kind']==B['kind']and block['L2']==B['L2'];H=fmpq_mat(block['matrix']);W=fmpq_mat(certificate['basisTransforms'][ib]);assert H.nrows()==H.ncols()==W.nrows()==W.ncols()==d;assert H==H.transpose()
 for j in range(1,d+1):assert fmpq_mat([[H[a,b]for b in range(j)]for a in range(j)]).det()>0
 G=W*H*W.transpose();rv=[i for i,v in enumerate(variables)if v[0]==ib];A=fmpq_mat([[fmpq(str(features[v].get(key,0)))for v in rv]for key in keys]);V=fmpq_mat([[G[variables[i][1],variables[i][2]]]for i in rv]);result+=A*V;positivity.append({'kind':B['kind'],'L2':B['L2'],'dimension':d});print('PASS positive block and coefficient contribution',ib,flush=True)
def product(p,q):
 out=defaultdict(int)
 for m,c in p.items():
  for n,d in q.items():out[tuple(sorted(m+n))]+=c*d
 return {m:c for m,c in out.items()if c}
norm={(i,12+i):comb(11,i)for i in range(12)};q2=product(norm,norm);q4=product(q2,q2);p6=defaultdict(int)
for t in range(7):
 wt=comb(11,6)*(-1)**(6-t)*comb(6,t)*2**t
 for a in range(t+1):
  for b in range(t+1):
   for i in range(12-t):
    for j in range(12-t):p6[tuple(sorted((a+i,12+b+i,b+j,12+a+j)))]+=wt*comb(t,a)*comb(t,b)*comb(11-t,i)*comb(11-t,j)
target=product(q2,p6);bound=F(certificate['bound'])
for m,c in q4.items():target[m]=target.get(m,0)-bound*c
target={m:c for m,c in target.items()if c}
assert not deriv(target,True)and not deriv(target,False)
assert all(len(m)==8 and sum(i<12 for i in m)==4 for m in target)
assert {tuple(sorted(i+12 if i<12 else i-12 for i in m)):c for m,c in target.items()}==target
target={encode(m):c for m,c in target.items()};assert result==fmpq_mat([[fmpq(str(target.get(key,0)))]for key in keys])
print('PASS full polynomial identity through exact invariant-space injectivity',flush=True)
report={'n':11,'sector':6,'direction':'min','degree':8,'bound':str(bound),'fullPolynomialIdentityVerified':True,'identityVerification':'177 exact coefficient equations; exact injectivity on the complete 177-dimensional invariant polynomial space; directly checked covariance and target membership','floatingPointUsed':False,'sharpnessClaimed':False,'attainmentVerified':False,'invariantDimension':177,'coefficientConstraints':177,'positiveDefiniteGramBlocks':len(raw),'maxMatrixDimension':max(B['dimension']for B in raw),'positivity':positivity,'inputSha256':{name:hashlib.sha256((BASE/name).read_bytes()).hexdigest()for name in ['basis.json','sos-certificate.json','verify_certificate.py']}}
(BASE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
