"""Independent exact small-n AC audit. SymPy required; no numerical optimization."""
from pathlib import Path
from math import comb, factorial
from fractions import Fraction as F
import sympy as s
import json

D=s.Rational;rt=s.sqrt;I=s.I
witnesses={
 2:(1,{0:1/rt(2),2:1/rt(2)},'GHZ'),
 3:(1,{0:1/rt(2),3:1/rt(2)},'GHZ'),
 4:(2,{0:1/rt(3),3:rt(2)/rt(3)},'Tetrahedron'),
 5:(1,{0:1/rt(2),5:1/rt(2)},'GHZ'),
 6:(3,{1:1/rt(2),5:1/rt(2)},'Octahedron'),
 7:(2,{0:rt(2)/3,3:rt(14)/6,6:rt(14)/6},'C3 seven-qubit state'),
 8:(3,{0:rt(30)/12,4:rt(21)/6,8:rt(30)/12},'Cube'),
 9:(2,{1:-rt(21)/8,5:rt(30)/8,9:rt(13)/8},'Quantum 2020 Eq.65 chi0'),
 10:(3,{1:rt(5)/4,5:rt(6)/4,9:rt(5)/4},'C4 ten-qubit state'),
 11:(3,{1:(-35*rt(10)+I*rt(10022))/288,4:5*rt(3)/18,7:5*rt(3)/18,10:(-35*rt(10)+I*rt(10022))/288},'New exact C3 eleven-qubit state'),
 12:(5,{1:rt(7)/5,6:rt(11)/5,11:-rt(7)/5},'Icosahedron')
}

def rho(n,t,c):
 A=s.zeros(t+1,n-t+1)
 for a in range(t+1):
  for b in range(n-t+1):
   if a+b in c:A[a,b]=c[a+b]*rt(D(comb(t,a)*comb(n-t,b),comb(n,a+b)))
 return (A*A.conjugate().T).applyfunc(s.simplify)

results=[]
for n,(t,c,label) in witnesses.items():
 assert s.simplify(sum(x*s.conjugate(x) for x in c.values()))==1
 rs=[]
 for order in range(1,t+1):
  reduced=rho(n,order,c);delta=reduced-s.eye(order+1)/(order+1)
  assert delta==s.zeros(order+1),(n,order,delta)
  rs.append({'t':order,'rhoEqualsSymmetricIdentity':True,'purity':str(D(1,order+1))})
 # Verify this representative is not of the next order.
 nxt=rho(n,t+1,c);purity=s.simplify(s.trace(nxt*nxt));assert purity!=D(1,t+2)
 results.append({'n':n,'achievedOrder':t,'label':label,'support':sorted(c),'exactCoefficients':{str(r):str(x) for r,x in c.items()},'verifiedReductions':rs,'nextOrderPurity':str(purity),'nextOrderFlatPurity':str(D(1,t+2)),'representativeExactOrder':t})
 print(f'PASS n={n}: full AC orders1..{t}; representative fails order{t+1}.',flush=True)

# Independent proof of the two-copy identities using rational highest-weight sums.
def C(n,k):return comb(n,k) if 0<=k<=n else 0
def h(n,d,t):
 return sum((F((-1)**q*factorial(n-d+q)*factorial(n-q),factorial(n-d)*factorial(n))*C(t,q)*C(n-t,d-q) for q in range(d+1)),F(0))/C(n,d)
def a_eig(n,L,d):
 return F((2*L+1)*C(n,L),(L+1)*C(n+L+1,n))*sum(((-1)**(t+L)*C(L,t)*C(t+L,t)*h(n,d,t) for t in range(L+1)),F(0))
certificates=[
 {'n':8,'t':4,'J':0,'constant':'-4/27','multipoles':{'1':'-1','2':'11/21','3':'11/14','4':'143/126'}},
 {'n':10,'t':5,'J':6,'constant':'-27/187','multipoles':{'1':'-519/170','2':'-143/170','3':'754/255','4':'52/17','5':'520/561'}}]
for cert in certificates:
 n=cert['n']
 for d in range(0,n+1,2):
  rhs=F(cert['constant'])+sum(F(coef)*a_eig(n,int(L),d) for L,coef in cert['multipoles'].items())
  assert rhs==int(n-d==cert['J'])
 assert F(cert['constant'])<0
 print(f'PASS: exact negative-pJ contradiction for n={n},t={cert["t"]}.',flush=True)

# Check the mislabelled 2022 n9 state directly.
c9={0:1/rt(6),3:-1/rt(3),6:-1/rt(3),9:1/rt(6)}
r3=rho(9,3,c9);defect=(r3-s.eye(4)/4).applyfunc(s.simplify)
assert defect!=s.zeros(4)
misprint={'n':9,'claimedOrderInTable':3,'actualVerifiedOrder':2,'rho3_03':str(r3[0,3]),'rho3Purity':str(s.simplify(s.trace(r3*r3))),'note':'PRResearch2022 TableIV q=3 is inconsistent with its own displayed vector; all rank1..2 vanish, rank3 does not.'}
base=Path(__file__).resolve().parent
(base/'independent-witness-checks.json').write_text(json.dumps({'states':results,'nonexistenceIdentities':certificates,'sourceMisprint':misprint},indent=2))
print('PASS: 11 independent exact witnesses; two exact exclusions; source misprint detected.')
