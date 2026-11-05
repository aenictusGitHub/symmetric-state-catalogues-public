"""Exact, independent certificate for S6>=7 in pure symmetric seven-qubit states.
Standard library only. Run beside basis.json and sos-certificate.json, or pass
that containing folder. The target is rebuilt directly from Dicke reduced
purities t=0,...,6. All complex coefficients are independent formal variables;
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
assert certificate['n']==7 and certificate['sector']==6 and F(certificate['bound'])==7

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
S6=defaultdict(F)
for t in range(7):
 for mon,value in purity(t).items():S6[mon]+=7*(-1)**(6-t)*comb(6,t)*2**t*value
S6={mon:value for mon,value in S6.items() if value}
# A short purity formulation of the same bound. Complementary reduced
# purities agree coefficientwise for pure (possibly unnormalized) states.
P1,P2,P3,P4,P5,P6=[purity(t) for t in range(1,7)]
assert P1==P6 and P2==P5 and P3==P4
short=defaultdict(F)
for factor,p in [(7,Q2),(364,P1),(-924,P2),(560,P3)]:
 for mon,value in p.items():short[mon]+=factor*value
short={mon:value for mon,value in short.items() if value}
assert short==S6
target=defaultdict(F,product(Q,S6))
for mon,value in Q3.items():target[mon]-=7*value
target={mon:value for mon,value in target.items() if value}
# The whole normalized generalized-GHZ family is attained, including phases.
assert not {m:c for m,c in target.items() if all(i%8 in (0,7) for i in m)}

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
  for p in cov:assert not {m:c for m,c in p.items() if all(i%8 in (0,7) for i in m)}
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
checks={'n':7,'sector':6,'bound':'7','targetDerivation':'Direct reduced Dicke purities t=0,...,6','independentComplexVariables':True,'floatingPointUsed':False,'fullPolynomialIdentityVerified':True,'nonzeroDegreeSixMonomials':len(result),'allMonomialsIncludingCancellations':len(allkeys),'quarticSectorMonomials':len(S6),'wholeComplexGeneralizedGHZFamilyVerified':True,'complementaryPurityIdentitiesVerified':True,'shortPurityIdentityVerified':True,'normalizedSectorIdentity':'S6=7+28(13 tr(rho1^2)+20 tr(rho3^2)-33 tr(rho2^2))','equivalentPurityInequality':'13 tr(rho1^2)+20 tr(rho3^2)>=33 tr(rho2^2)','positiveDefiniteGramBlocks':len(positivity),'squaredForms':sum(B['dimension']*(B['L2']+1) for B in basis['blocks']),'positivity':positivity,'basisSha256':hashlib.sha256(basis_file.read_bytes()).hexdigest(),'certificateSha256':hashlib.sha256(certificate_file.read_bytes()).hexdigest()}
(BASE/'independent-verification.json').write_text(json.dumps(checks,indent=2)+'\n')
