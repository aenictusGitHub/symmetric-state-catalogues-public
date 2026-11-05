"""Independent exact sector-bound verifier for complex symmetric ten-qubit states.
Standard library only; reconstructs the polynomial from Dicke reduced purities.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import comb,lcm
import json,sys,hashlib
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
basis_file=BASE/'basis.json';certificate_file=BASE/'sos-certificate.json'
basis=json.loads(basis_file.read_text());certificate=json.loads(certificate_file.read_text())
config=json.loads((BASE/'config.json').read_text())
assert certificate['n']==10 and certificate['sector']==6 and F(certificate['bound'])==F('58437/800')

def product(p,q):
 out=defaultdict(int)
 for m,x in p.items():
  for n,y in q.items():out[tuple(sorted(m+n))]+=x*y
 return {m:x for m,x in out.items() if x}

def purity(t):
 p=defaultdict(F)
 for i in range(t+1):
  for j in range(t+1):
   for a in range(11-t):
    for b in range(11-t):
     p[tuple(sorted((i+a,11+j+a,j+b,11+i+b)))]+=comb(t,i)*comb(t,j)*comb(10-t,a)*comb(10-t,b)
 return dict(p)
Q={(i,11+i):F(comb(10,i)) for i in range(11)}
Q2=product(Q,Q);Q3=product(Q2,Q)
def sector(k):
 out=defaultdict(F)
 for t in range(k+1):
  for mon,value in purity(t).items():out[mon]+=comb(10,k)*(-1)**(k-t)*comb(k,t)*2**t*value
 return {m:c for m,c in out.items() if c}
SK=sector(6)
target=defaultdict(F,{m:F('58437/800')*v for m,v in Q3.items()})
for m,v in product(Q,SK).items():target[m]-=v
target={m:v for m,v in target.items()if v}
target={m:-v for m,v in target.items()}
target=product(target,Q)
# Independent Sylvester test by exact determinants of every leading minor.
def determinant(M):
 A=[row[:] for row in M];det=F(1)
 for j in range(len(A)):
  pivot=next((i for i in range(j,len(A)) if A[i][j]),None)
  if pivot is None:return F(0)
  if pivot!=j:A[j],A[pivot]=A[pivot],A[j];det=-det
  v=A[j][j];det*=v
  for i in range(j+1,len(A)):
   factor=A[i][j]/v
   for k in range(j+1,len(A)):A[i][k]-=factor*A[j][k]
 return det
assert len(basis['blocks'])==len(certificate['blocks'])==len(certificate['basisScales'])
# Encode exponent vectors in base 9: degree eight guarantees no digit carries.
# This is an exact monomial bijection, used only to speed polynomial expansion.
def encode(mon):return sum(9**i for i in mon)
CONJ_SPLIT=9**11
def conjugate(mon):return mon//CONJ_SPLIT+(mon%CONJ_SPLIT)*CONJ_SPLIT
target={encode(m):c for m,c in target.items()}
positivity=[];prepared=[];globalden=lcm(*(v.denominator for v in target.values()));allkeys=set(target)
for B,G,scale in zip(basis['blocks'],certificate['blocks'],certificate['basisScales']):
 assert B['kind']==G['kind'] and B['L2']==G['L2']
 n=B['dimension'];H=[[F(v) for v in row] for row in G['matrix']]
 assert len(H)==n and all(len(row)==n for row in H)
 assert all(H[i][j]==H[j][i] for i in range(n) for j in range(n))
 minors=[determinant([row[:j] for row in H[:j]]) for j in range(1,n+1)]
 assert all(m>0 for m in minors)
 positivity.append({'kind':B['kind'],'L2':B['L2'],'dimension':n,'leadingPrincipalMinorSigns':[1 for m in minors]})
 polys=[]
 for cov in B['covariants']:
  components=[]
  for component in cov:
   P={}
   for mon,coef in component:
    assert len(mon)==4 and all(0<=i<22 for i in mon)
    value=F(coef);assert value.denominator==1
    P[encode(mon)]=value.numerator
   components.append(P)
  assert len(components)==B['L2']+1
  polys.append(components)
 assert len(polys)==n
 weights=[F(w) for w in B['weights']]
 assert weights==[F(1,comb(B['L2'],r)) for r in range(B['L2']+1)]
 factors=[[[w*H[a][b]/(F(scale[a])*F(scale[b])) for b in range(n)] for a in range(n)] for w in weights]
 denominator=lcm(*(c.denominator for plane in factors for row in plane for c in row))
 globalden=lcm(globalden,denominator)
 prepared.append((B,n,polys,factors))
 print('PASS exact positive matrix',B['kind'],B['L2'],'dimension',n,flush=True)
result=defaultdict(int)
for B,n,polys,factors in prepared:
 for r in range(B['L2']+1):
  for a in range(n):
   mixed=defaultdict(int)
   for b in range(n):
    factor=factors[r][a][b]*globalden;assert factor.denominator==1
    if not factor:continue
    for m,c in polys[b][r].items():mixed[conjugate(m)]+=factor.numerator*c
   for m,c in polys[a][r].items():
    for mon,v in mixed.items():result[m+mon]+=c*v
 print('PASS full block expansion',B['kind'],B['L2'],flush=True)
allkeys.update(result);result={m:c for m,c in result.items()if c}
assert all((v*globalden).denominator==1 for v in target.values())
assert result=={m:int(v*globalden)for m,v in target.items()},'Full polynomial coefficient identity failed'
print('PASS complete complex-polynomial identity:',len(result),'nonzero monomials;',len(allkeys),'including cancellations.',flush=True)
checks={'n':10,'sector':6,'bound':'58437/800','direction':'min','degree':8,'targetDerivation':'Direct reduced Dicke purities','independentComplexVariables':True,'floatingPointUsed':False,'fullPolynomialIdentityVerified':True,'nonzeroMonomials':len(result),'allMonomialsIncludingCancellations':len(allkeys),'attainmentVerified':False,'sharpnessClaimed':False,'positiveDefiniteGramBlocks':len(positivity),'squaredForms':sum(B['dimension']*(B['L2']+1) for B in basis['blocks']),'maxMatrixDimension':max(B['dimension'] for B in basis['blocks']),'positivity':positivity,'basisSha256':hashlib.sha256(basis_file.read_bytes()).hexdigest(),'certificateSha256':hashlib.sha256(certificate_file.read_bytes()).hexdigest()}
(BASE/'independent-verification.json').write_text(json.dumps(checks,indent=2)+'\n')
