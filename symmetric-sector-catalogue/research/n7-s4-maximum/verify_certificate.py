"""Exact, independent certificate for S4<=265/7 in pure symmetric seven-qubit states.
Standard library only. Run beside basis.json and sos-certificate.json, or pass
that containing folder. The target is rebuilt directly from Dicke reduced
purities t=0,...,4. All complex coefficients are independent formal variables;
no optimization, numerical tolerance or Clebsch-Gordan routine is involved.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import comb
import json,sys,hashlib
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
basis_file=BASE/'basis.json';certificate_file=BASE/'sos-certificate.json'
basis=json.loads(basis_file.read_text());certificate=json.loads(certificate_file.read_text())
assert certificate['n']==7 and certificate['sector']==4 and F(certificate['bound'])==F(265,7)

def product(p,q):
 out=defaultdict(F)
 for m,x in p.items():
  for n,y in q.items():out[tuple(sorted(m+n))]+=x*y
 return {m:x for m,x in out.items() if x}

def purity(t):
 p=defaultdict(F)
 for i in range(t+1):
  for j in range(t+1):
   for a in range(8-t):
    for b in range(8-t):
     p[tuple(sorted((i+a,8+j+a,j+b,8+i+b)))]+=comb(t,i)*comb(t,j)*comb(7-t,a)*comb(7-t,b)
 return dict(p)
Q={(i,8+i):F(comb(7,i)) for i in range(8)}
Q2=product(Q,Q);Q3=product(Q2,Q)
def sector(k):
 out=defaultdict(F)
 for t in range(k+1):
  for mon,value in purity(t).items():out[mon]+=comb(7,k)*(-1)**(k-t)*comb(k,t)*2**t*value
 return {m:c for m,c in out.items() if c}
S4=sector(4)
target=defaultdict(F,{m:F(265,7)*v for m,v in Q3.items()})
for mon,value in product(Q,S4).items():target[mon]-=value
target={mon:value for mon,value in target.items() if value}
# Exact restriction to x = alpha*(0,-1,0,0,sqrt(2),0,0,7)
#                       + beta*(7,0,0,-sqrt(2),0,0,-1,0).
# Keep alpha,beta and their conjugates formally independent, and use Q(sqrt(2)).
mapping={0:(1,7,0),1:(0,-1,0),3:(1,-1,1),4:(0,1,1),6:(1,-1,0),7:(0,7,0)}
def restrict_family(poly):
 out=defaultdict(F)
 for mon,coef in poly.items():
  exponents=[0]*4;radicals=0
  for i in mon:
   if i%8 not in mapping:break
   variable,factor,radical=mapping[i%8]
   exponents[variable+(2 if i>=8 else 0)]+=1;coef*=factor;radicals+=radical
  else:out[tuple(exponents)+(radicals%2,)]+=coef*2**(radicals//2)
 return {m:c for m,c in out.items() if c}
assert not restrict_family(target)
assert restrict_family(Q)=={(1,0,1,0,0):F(126),(0,1,0,1,0):F(126)}
def evaluate_attainer(poly):
 out=[F(0),F(0)]
 for (a,b,c,d,rad),value in restrict_family(poly).items():
  if b==d==0:out[rad]+=value
 assert out[1]==0
 return out[0]
sectors=[evaluate_attainer(sector(k))/126**2 for k in range(1,8)]
assert sectors==list(map(F,['9/7','7','5','265/7','223/7','127/7','181/7']))
print('PASS exact normalized attaining state and whole complex Kramers family.',flush=True)

# Independent Sylvester test by exact determinants of every leading minor.
def determinant(M):
 if len(M)==1:return M[0][0]
 return sum((-1)**j*M[0][j]*determinant([[M[r][k] for k in range(len(M)) if k!=j] for r in range(1,len(M))]) for j in range(len(M)))
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
  for p in cov:assert not restrict_family(p)
 # Sum all ordered matrix entries directly; formal conjugation swaps x_i,y_i.
 for r,w in enumerate(weights):
  for a in range(n):
   for b in range(n):
    coef=w*H[a][b]
    if not coef:continue
    pa=polys[a][r]
    pb={tuple(sorted((i+8)%16 for i in mon)):c for mon,c in polys[b][r].items()}
    for mon,c in product(pa,pb).items():result[mon]+=coef*c;allkeys.add(mon)
 print('PASS positive Gram block',B['kind'],B['L2'],'dimension',n,flush=True)
result={mon:value for mon,value in result.items() if value}
assert result==target,'Full polynomial coefficient identity failed'
print('PASS complete complex-polynomial identity:',len(result),'nonzero monomials;',len(allkeys),'including cancellations.',flush=True)
checks={'n':7,'sector':4,'bound':'265/7','targetDerivation':'Direct reduced Dicke purities t=0,...,4','independentComplexVariables':True,'floatingPointUsed':False,'fullPolynomialIdentityVerified':True,'nonzeroDegreeSixMonomials':len(result),'allMonomialsIncludingCancellations':len(allkeys),'quarticSectorMonomials':len(S4),'wholeComplexKramersFamilyVerified':True,'attainingStateNormSquared':1,'attainingStateSectorVector':[str(v) for v in sectors],'positiveDefiniteGramBlocks':len(positivity),'squaredForms':sum(B['dimension']*(B['L2']+1) for B in basis['blocks']),'positivity':positivity,'basisSha256':hashlib.sha256(basis_file.read_bytes()).hexdigest(),'certificateSha256':hashlib.sha256(certificate_file.read_bytes()).hexdigest()}
(BASE/'independent-verification.json').write_text(json.dumps(checks,indent=2)+'\n')
