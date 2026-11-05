"""Independent exact checks of packing-equality classification for t=1..4.
No bound implementation is imported. Analytic all-N cutoffs are in the proof.
"""
from pathlib import Path
from itertools import product
import json
import sympy as s
from sympy.physics.wigner import wigner_6j

# Equality plus the auxiliary-spin inequalities gives finite cutoffs.
n=s.symbols('n',integer=True,positive=True)
assert s.expand(3*n/2-(n+3))==(n-6)/2
assert s.expand(5*(n+1)/4-(n+5))==(n-15)/4
assert s.expand(7*n/6-(n+7))==(n-42)/6
assert s.expand(10*(n+1)/9-(n+10))==(n-80)/9

third=[]
for N in range(6,43,6):
 j=s.Rational(N,2);jp=j-s.Rational(1,2);roots=[]
 for L in range(4,N):
  multiplier=s.simplify(sum((-1)**(N+L+q)*(2*q+2)*
   wigner_6j(j,j,L,jp,jp,s.Rational(2*q+1,2))for q in range(2)))
  pell=2*L*(L+1)==(N-1)*(N+2)
  assert (multiplier==0)==pell
  if multiplier==0:roots.append(L)
 assert roots==([4]if N==6 else [])
 third.append({'n':N,'admissibleChannelRoots':roots,'direct6jAndPellAgree':True})

fourth=[]
for N in range(8,81,9):
 j=s.Rational(N,2);Q=N*(N+2);roots=[]
 for L in range(5,N+1):
  multiplier=s.simplify(sum((-1)**(N+q+L)*(2*q+1)*
   wigner_6j(j,j,q,j,j,L)for q in range(3)))
  x=L*(L+1)
  closed=s.Rational(3*(10*x*x+(16-12*Q)*x+3*Q*(Q-3)),(N+1)*Q*(Q-3))
  assert multiplier==closed
  if multiplier==0:roots.append(L)
 assert roots==([16]if N==17 else [])
 fourth.append({'n':N,'admissibleChannelRoots':roots,'direct6jAndRationalPolynomialAgree':True})
# Exact algebra in the cited N17 rank-sixteen lifting exclusion.
alpha=s.Rational(4,17)
assert 1-2*alpha==s.Rational(9,17)
assert alpha-alpha**2==s.Rational(52,289)
assert 17*(alpha-alpha**2)==s.Rational(52,17)
assert (1-2*alpha)/(17*(alpha-alpha**2))==s.Rational(9,52)
assert sum(2*l+1 for l in range(8))==64>33

# t2 N3 fails the pure-state N>=4 bound; N11/15 fail 4k<=N-2.
second=[]
for N in [3,7,11,15]:
 k=(N+1)//4
 allowed=N>=4 and (N<8 or 4*k<=N-2)
 assert allowed==(N==7)
 second.append({'n':N,'k':k,'notExcludedByCitedExactTheorems':allowed})

def spin(N):
 J=s.Rational(N,2);d=N+1;z=s.diag(*[J-a for a in range(d)]);p=s.zeros(d)
 for a in range(1,d):p[a-1,a]=s.sqrt(a*(N-a+1))
 return [(p+p.T)/2,(p-p.T)/(2*s.I),z]

def vector(N,terms):
 v=s.zeros(N+1,1)
 for a,c in terms.items():v[a]=c
 return v

rt=lambda a,b:s.sqrt(s.Rational(a,b))
witnesses=[
 (2,1,[{0:rt(1,2),2:rt(1,2)}]),
 (4,1,[{0:rt(1,2),4:rt(1,2)},{2:s.Integer(1)}]),
 (6,1,[{1:rt(3,5),6:rt(2,5)},{0:rt(2,5),5:-rt(3,5)},{3:s.Integer(1)}]),
 (7,2,[{0:rt(3,10),5:rt(7,10)},{2:rt(7,10),7:-rt(3,10)}]),
 (6,3,[{1:rt(1,2),5:rt(1,2)}])]
checks=[]
for N,t,terms in witnesses:
 U=s.Matrix.hstack(*(vector(N,v)for v in terms));k=U.cols;J=spin(N)
 assert (U.H*U-s.eye(k)).applyfunc(s.simplify)==s.zeros(k)
 words=0
 for q in range(1,t+1):
  for indices in product(range(3),repeat=q):
   A=s.eye(N+1)
   for a in indices:A=A*J[a]
   assert (U.H*A*U-s.trace(A)*s.eye(k)/(N+1)).applyfunc(s.simplify)==s.zeros(k)
   words+=1
 checks.append({'n':N,'t':t,'k':k,'allOrderedSpinWordCompressionsExact':True,'wordCount':words})

result={'passed':True,'classification':{'t1':[[2,1],[4,2],[6,3]],'t2':[[7,2]],'t3':[[6,1]],'t4':[]},'meaning':'Equality in unrounded rectangular-tensor bounds, not in rounded integer upper bounds.','cutoffs':{'t1':6,'t2':15,'t3':42,'t4':80},'t2CandidateChecks':second,'t3CandidateChecks':third,'t4CandidateChecks':fourth,'n17Exclusion':{'proof':'../n17-t4/proof.tex','rankSixteenProjectorLiftArithmeticChecked':True,'pure32QubitOrder15Contradiction':'64>33'},'attainingWitnessChecks':checks,'strictCorollaries':['6K(N,3)<=N-1 for N>=7','9K(N,4)<=N for N>=4'],'importsCatalogueBounds':False}
Path(__file__).with_name('packing-tightness-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
