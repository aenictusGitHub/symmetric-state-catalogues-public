"""Independently check the N=6 S4 >=45/7 exact SOS certificate.
Standard library only. Put this file beside basis.json and sos-certificate.json,
or pass their containing folder as the first command-line argument.
The target is rebuilt directly from Dicke reduced density matrices, including
all purities t=0,...,4. No optimization, floating point or CG routines are used.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import comb
import json,sys,hashlib
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
basis_file=BASE/'basis.json';certificate_file=BASE/'sos-certificate.json'
basis=json.loads(basis_file.read_text());certificate=json.loads(certificate_file.read_text())
assert certificate['n']==6 and certificate['sector']==4 and F(certificate['bound'])==F(45,7)

def product(p,q):
 out=defaultdict(F)
 for m,x in p.items():
  for n,y in q.items():out[tuple(sorted(m+n))]+=x*y
 return {m:x for m,x in out.items() if x}

def purity(t):
 p=defaultdict(F)
 for i in range(t+1):
  for j in range(t+1):
   for a in range(7-t):
    for b in range(7-t):
     p[tuple(sorted((i+a,7+j+a,j+b,7+i+b)))]+=comb(t,i)*comb(t,j)*comb(6-t,a)*comb(6-t,b)
 return dict(p)
Q={(i,7+i):F(comb(6,i)) for i in range(7)}
Q2=product(Q,Q);Q4=product(Q2,Q2)
S4=defaultdict(F)
for t in range(5):
 for mon,value in purity(t).items():S4[mon]+=15*(-1)**(4-t)*comb(4,t)*2**t*value
S4={mon:value for mon,value in S4.items() if value}
target=defaultdict(F,product(Q2,S4))
for mon,value in Q4.items():target[mon]-=F(45,7)*value
target={mon:value for mon,value in target.items() if value}

# Sylvester's criterion, with exact fraction-preserving elimination, is
# independent of the producer's Schur-complement LDL implementation.
def determinant(M):
 M=[r[:] for r in M];sign=1
 for k in range(len(M)-1):
  if not M[k][k]:
   j=next(j for j in range(k+1,len(M)) if M[j][k])
   M[k],M[j]=M[j],M[k];sign=-sign
  pivot=M[k][k]
  for i in range(k+1,len(M)):
   ratio=M[i][k]/pivot
   for j in range(k+1,len(M)):M[i][j]-=ratio*M[k][j]
   M[i][k]=F(0)
 result=F(sign)
 for i in range(len(M)):result*=M[i][i]
 return result

assert len(basis['blocks'])==len(certificate['blocks'])==len(certificate['basisScales'])
result=defaultdict(F);positivity=[]
for ib,(B,G,scale) in enumerate(zip(basis['blocks'],certificate['blocks'],certificate['basisScales'])):
 assert B['kind']==G['kind'] and B['L']==G['L']
 n=B['dimension'];H=[[F(v) for v in row] for row in G['matrix']]
 assert len(H)==n and all(len(row)==n for row in H)
 assert all(H[i][j]==H[j][i] for i in range(n) for j in range(n))
 minors=[determinant([row[:j] for row in H[:j]]) for j in range(1,n+1)]
 assert all(m>0 for m in minors)
 positivity.append({'kind':B['kind'],'L':B['L'],'dimension':n,'leadingPrincipalMinors':[str(m) for m in minors]})
 polys=[[{tuple(mon):F(coef)/F(scale[a]) for mon,coef in component} for component in cov] for a,cov in enumerate(B['covariants'])]
 assert len(polys)==n
 weights=[F(w) for w in B['weights']]
 assert weights==[F(1,comb(2*B['L'],r)) for r in range(2*B['L']+1)]
 # Sum all ordered matrix entries directly instead of using preassembled Gram
 # features. Formal conjugation exchanges independent x_i and y_i variables.
 for r,w in enumerate(weights):
  for a in range(n):
   for b in range(n):
    coef=w*H[a][b]
    if not coef:continue
    pa=polys[a][r]
    pb={tuple(sorted((i+7)%14 for i in mon)):c for mon,c in polys[b][r].items()}
    for mon,c in product(pa,pb).items():result[mon]+=coef*c
 print('PASS positive Gram block',B['kind'],B['L'],'dimension',n,flush=True)
result={mon:value for mon,value in result.items() if value}
assert result==target,'Full polynomial coefficient identity failed'
print('PASS complete complex-polynomial identity:',len(result),'nonzero monomials.',flush=True)
checks={'n':6,'sector':4,'bound':'45/7','targetDerivation':'Direct reduced Dicke purities t=0,...,4','independentComplexVariables':True,'floatingPointUsed':False,'fullPolynomialIdentityVerified':True,'nonzeroDegreeEightMonomials':len(result),'quarticSectorMonomials':len(S4),'positiveDefiniteGramBlocks':len(positivity),'squaredForms':sum(B['dimension']*(2*B['L']+1) for B in basis['blocks']),'positivity':positivity,'basisSha256':hashlib.sha256(basis_file.read_bytes()).hexdigest(),'certificateSha256':hashlib.sha256(certificate_file.read_bytes()).hexdigest()}
(BASE/'independent-verification.json').write_text(json.dumps(checks,indent=2)+'\n')
