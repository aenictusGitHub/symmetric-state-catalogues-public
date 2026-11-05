from pathlib import Path
import sympy as s
from sympy.physics.wigner import wigner_6j
import json
from math import isqrt
checks=[]
for n in range(3,25):
 j=s.Rational(n,2);d=n+1
 c=[s.Rational(n*n*(n-1)**2*(n+2)**2,24),s.Rational(n*(n-1)**2*(n+2)*(3*n*n+8*n+12),40),s.Rational(n*(n-1)*(n+2)*(n+3)*(n*n-2*n+4),12),s.Rational(n*(n-2)*(n-1)*(n+2)*(n+3)*(n+4),20)]
 assert all(x>0 for x in c)
 rr=s.Rational((n-1)*(n+2),2)
 for L in range(n+1):
  z=L*(L+1)
  mus=[(-1)**(n+ell+L)*(n+1)*wigner_6j(j,j,ell,j,j,L)for ell in range(4)]
  f=(n*(n+1)-z)*(z-rr)**2
  assert sum(x*y for x,y in zip(c,mus))==f
  assert f>=0
 assert s.simplify((n+1)*c[0]/(n*(n+1)*rr**2))==s.Rational(n,6)
 checks.append({'n':n,'allRanksVerified':True,'possibleAdditionalRank':[L for L in range(1,n+1)if 2*L*(L+1)==(n-1)*(n+2)]})
exceptions=[]
for n in range(6,100001,6):
 for L in [(isqrt(1+2*(n-1)*(n+2))-1)//2]:
  if 2*L*(L+1)==(n-1)*(n+2):exceptions.append([n,L])
report={'passed':True,'exactChecks':checks,'generalExpansionCoefficients':['N²(N−1)²(N+2)²/24','N(N−1)²(N+2)(3N²+8N+12)/40','N(N−1)(N+2)(N+3)(N²−2N+4)/12','N(N−2)(N−1)(N+2)(N+3)(N+4)/20'],'bound':'6k<=N','saturationNecessaryConditionForNdiv6andN>=12':'2L(L+1)=(N−1)(N+2) for some1<=L<=N','enumeratedPossibleSaturationExceptionsThroughN100000':exceptions}
Path(__file__).with_name('t3-saturation-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
