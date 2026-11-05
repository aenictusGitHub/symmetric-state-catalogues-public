"""Independent exact sector-bound verifier for complex symmetric twelve-qubit states.
Standard-library implementation with optional exact FLINT acceleration.
Reconstructs the polynomial from Dicke reduced purities.
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
assert certificate['n']==12 and certificate['sector']==6 and F(certificate['bound'])==F('713439361/5000000')

def product(p,q):
 out=defaultdict(int)
 for m,x in p.items():
  for n,y in q.items():out[tuple(sorted(m+n))]+=x*y
 return {m:x for m,x in out.items() if x}

def purity(t):
 p=defaultdict(F)
 for i in range(t+1):
  for j in range(t+1):
   for a in range(13-t):
    for b in range(13-t):
     p[tuple(sorted((i+a,13+j+a,j+b,13+i+b)))]+=comb(t,i)*comb(t,j)*comb(12-t,a)*comb(12-t,b)
 return dict(p)
Q={(i,13+i):F(comb(12,i)) for i in range(13)}
Q2=product(Q,Q);Q3=product(Q2,Q)
def sector(k):
 out=defaultdict(F)
 for t in range(k+1):
  for mon,value in purity(t).items():out[mon]+=comb(12,k)*(-1)**(k-t)*comb(k,t)*2**t*value
 return {m:c for m,c in out.items() if c}
SK=sector(6)
target=defaultdict(F,{m:F('713439361/5000000')*v for m,v in Q3.items()})
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
def positive_minor_signs(H):
 # Optional exact-arithmetic acceleration; the standard-library path below
 # proves the same leading-minor signs when python-flint is unavailable.
 try:
  from flint import fmpq_mat
 except ImportError:
  fmpq_mat=None
 if fmpq_mat is not None:
  for j in range(1,len(H)+1):assert fmpq_mat([[str(v)for v in row[:j]]for row in H[:j]]).det()>0
  return [1]*len(H)
 # Fraction-free Bareiss elimination. Every pivot is a leading determinant
 # of the matrix scaled by one positive common denominator.
 den=lcm(*(v.denominator for row in H for v in row));M=[[int(v*den)for v in row]for row in H];previous=1;signs=[]
 for k in range(len(M)):
  pivot=M[k][k];assert pivot>0;signs.append(1)
  for i in range(k+1,len(M)):
   for j in range(k+1,len(M)):
    numerator=pivot*M[i][j]-M[i][k]*M[k][j];assert numerator%previous==0;M[i][j]=numerator//previous
  previous=pivot
 return signs
assert len(basis['blocks'])==len(certificate['blocks'])==len(certificate['basisScales'])
# Encode exponent vectors in base 9: degree eight guarantees no digit carries.
# This is an exact monomial bijection, used only to speed polynomial expansion.
def encode(mon):return sum(9**i for i in mon)
CONJ_SPLIT=9**13
def conjugate(mon):return mon//CONJ_SPLIT+(mon%CONJ_SPLIT)*CONJ_SPLIT
target={encode(m):c for m,c in target.items()}
positivity=[];prepared=[];globalden=lcm(*(v.denominator for v in target.values()));allkeys=set(target)
for B,G,scale in zip(basis['blocks'],certificate['blocks'],certificate['basisScales']):
 assert B['kind']==G['kind'] and B['L2']==G['L2']
 n=B['dimension'];H=[[F(v) for v in row] for row in G['matrix']]
 assert len(H)==n and all(len(row)==n for row in H)
 assert all(H[i][j]==H[j][i] for i in range(n) for j in range(n))
 minors=positive_minor_signs(H)
 assert all(m>0 for m in minors)
 positivity.append({'kind':B['kind'],'L2':B['L2'],'dimension':n,'leadingPrincipalMinorSigns':[1 for m in minors]})
 polys=[]
 for cov in B['covariants']:
  components=[]
  for component in cov:
   P={}
   for mon,coef in component:
    assert len(mon)==4 and all(0<=i<26 for i in mon)
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
try:
 from flint import fmpz_mpoly_ctx
except ImportError:
 fmpz_mpoly_ctx=None
if fmpz_mpoly_ctx is not None:
 # Exact sparse-polynomial acceleration. It expands the same identity as the
 # standard-library fallback; no numerical coefficients or tolerances enter.
 ctx=fmpz_mpoly_ctx.get(tuple('v'+str(i)for i in range(26)),ordering='lex')
 def decode(m):
  ex=[]
  for i in range(26):ex.append(m%9);m//=9
  assert m==0
  return tuple(ex)
 poly_result=ctx.constant(0)
 for B,n,polys,factors in prepared:
  for r in range(B['L2']+1):
   pp=[ctx.from_dict({decode(m):c for m,c in P[r].items()})for P in polys]
   pc=[ctx.from_dict({decode(conjugate(m)):c for m,c in P[r].items()})for P in polys]
   for a in range(n):
    mixed=ctx.constant(0)
    for b in range(n):
     factor=factors[r][a][b]*globalden;assert factor.denominator==1
     if factor:mixed+=factor.numerator*pc[b]
    poly_result+=pp[a]*mixed
  print('PASS full exact block expansion',B['kind'],B['L2'],flush=True)
 result={sum(e*9**i for i,e in enumerate(ex)):int(c)for ex,c in poly_result.to_dict().items()}
 # Include all weight-preserving monomials of the target bidegree, including
 # those whose coefficients cancel to zero on both sides.
 from itertools import combinations_with_replacement
 byweight=defaultdict(list)
 for mon in combinations_with_replacement(range(13),4):byweight[sum(mon)].append(encode(mon))
 allkeys.update(x+CONJ_SPLIT*y for group in byweight.values()for x in group for y in group)
else:
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
checks={'n':12,'sector':6,'bound':'713439361/5000000','direction':'min','degree':8,'targetDerivation':'Direct reduced Dicke purities','independentComplexVariables':True,'floatingPointUsed':False,'fullPolynomialIdentityVerified':True,'nonzeroMonomials':len(result),'allMonomialsIncludingCancellations':len(allkeys),'attainmentVerified':False,'sharpnessClaimed':False,'positiveDefiniteGramBlocks':len(positivity),'squaredForms':sum(B['dimension']*(B['L2']+1) for B in basis['blocks']),'maxMatrixDimension':max(B['dimension'] for B in basis['blocks']),'positivity':positivity,'basisSha256':hashlib.sha256(basis_file.read_bytes()).hexdigest(),'certificateSha256':hashlib.sha256(certificate_file.read_bytes()).hexdigest()}
(BASE/'independent-verification.json').write_text(json.dumps(checks,indent=2)+'\n')
