"""Independent exact sector-bound verifier for complex symmetric eight-qubit states.
Standard library only; reconstructs the polynomial from Dicke reduced purities.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import comb
import json,sys,hashlib
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
basis_file=BASE/'basis.json';certificate_file=BASE/'sos-certificate.json'
basis=json.loads(basis_file.read_text());certificate=json.loads(certificate_file.read_text())
config=json.loads((BASE/'config.json').read_text())
assert certificate['n']==8 and certificate['sector']==5 and F(certificate['bound'])==F('147/2')

def product(p,q):
 out=defaultdict(F)
 for m,x in p.items():
  for n,y in q.items():out[tuple(sorted(m+n))]+=x*y
 return {m:x for m,x in out.items() if x}

def purity(t):
 p=defaultdict(F)
 for i in range(t+1):
  for j in range(t+1):
   for a in range(9-t):
    for b in range(9-t):
     p[tuple(sorted((i+a,9+j+a,j+b,9+i+b)))]+=comb(t,i)*comb(t,j)*comb(8-t,a)*comb(8-t,b)
 return dict(p)
Q={(i,9+i):F(comb(8,i)) for i in range(9)}
Q2=product(Q,Q);Q3=product(Q2,Q)
def sector(k):
 out=defaultdict(F)
 for t in range(k+1):
  for mon,value in purity(t).items():out[mon]+=comb(8,k)*(-1)**(k-t)*comb(k,t)*2**t*value
 return {m:c for m,c in out.items() if c}
SK=sector(5)
target=defaultdict(F,{m:F('147/2')*v for m,v in Q3.items()})
for m,v in product(Q,SK).items():target[m]-=v
target={m:v for m,v in target.items() if v}
excitation=1
def evaluate_w(poly):return sum(c for m,c in poly.items() if all(i%9==excitation for i in m))
norm=F(comb(8,excitation))
assert evaluate_w(Q)==norm and evaluate_w(target)==0
sectors=[evaluate_w(sector(k))/norm**2 for k in range(1,9)]
assert sectors[4]==F('147/2')
print('PASS exact Dicke attaining state and complete sector vector.',flush=True)

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
result=defaultdict(F);positivity=[];allkeys=set(target)
for ib,(B,G,scale) in enumerate(zip(basis['blocks'],certificate['blocks'],certificate['basisScales'])):
 assert B['kind']==G['kind'] and B['L2']==G['L2']
 n=B['dimension'];H=[[F(v) for v in row] for row in G['matrix']]
 assert len(H)==n and all(len(row)==n for row in H)
 assert all(H[i][j]==H[j][i] for i in range(n) for j in range(n))
 minors=[determinant([row[:j] for row in H[:j]]) for j in range(1,n+1)]
 assert all(m>0 for m in minors)
 positivity.append({'kind':B['kind'],'L2':B['L2'],'dimension':n,'leadingPrincipalMinors':[str(m) for m in minors]})
 polys=[[{tuple(mon):F(coef)/F(scale[a]) for mon,coef in component} for component in cov] for a,cov in enumerate(B['covariants'])]
 assert len(polys)==n
 weights=[F(w) for w in B['weights']]
 assert weights==[F(1,comb(B['L2'],r)) for r in range(B['L2']+1)]
 for cov in polys:
  for p in cov:assert evaluate_w(p)==0
 # Sum all ordered matrix entries directly; formal conjugation swaps x_i,y_i.
 for r,w in enumerate(weights):
  for a in range(n):
   for b in range(n):
    coef=w*H[a][b]
    if not coef:continue
    pa=polys[a][r]
    pb={tuple(sorted((i+9)%18 for i in mon)):c for mon,c in polys[b][r].items()}
    for mon,c in product(pa,pb).items():result[mon]+=coef*c;allkeys.add(mon)
 print('PASS positive Gram block',B['kind'],B['L2'],'dimension',n,flush=True)
result={mon:value for mon,value in result.items() if value}
assert result==target,'Full polynomial coefficient identity failed'
print('PASS complete complex-polynomial identity:',len(result),'nonzero monomials;',len(allkeys),'including cancellations.',flush=True)
checks={'n':8,'sector':5,'bound':'147/2','direction':'max','degree':6,'targetDerivation':'Direct reduced Dicke purities','independentComplexVariables':True,'floatingPointUsed':False,'fullPolynomialIdentityVerified':True,'nonzeroMonomials':len(result),'allMonomialsIncludingCancellations':len(allkeys),'attainmentVerified':True,'attainingDickeExcitation':excitation,'attainingStateSectorVector':[str(v) for v in sectors],'positiveDefiniteGramBlocks':len(positivity),'squaredForms':sum(B['dimension']*(B['L2']+1) for B in basis['blocks']),'maxMatrixDimension':max(B['dimension'] for B in basis['blocks']),'positivity':positivity,'basisSha256':hashlib.sha256(basis_file.read_bytes()).hexdigest(),'certificateSha256':hashlib.sha256(certificate_file.read_bytes()).hexdigest()}
(BASE/'independent-verification.json').write_text(json.dumps(checks,indent=2)+'\n')
