#!/usr/bin/env python3
"""Exact distinct-root sequence approaching a 3-AC state in another partition."""
from pathlib import Path
import json, math
import sympy as s
z,e=s.symbols('z epsilon',real=True)
Q=z**4-6*z**2+1
R=(z*z+1)**2+e
P=s.Poly(Q**2*R,z)
p=[P.nth(r) for r in range(13)]
norm=s.factor(sum(p[r]**2/s.binomial(12,r) for r in range(13)))
rho=s.zeros(4)
for a in range(4):
 for b in range(4):
  rho[a,b]=s.sqrt(s.binomial(3,a)*s.binomial(3,b))*sum(s.binomial(9,k)*p[a+k]*p[b+k]/(s.binomial(12,a+k)*s.binomial(12,b+k)) for k in range(10))/norm
assert (rho.subs(e,0)-s.eye(4)/4).applyfunc(s.simplify)==s.zeros(4)
hs2=s.factor(s.trace((rho-s.eye(4)/4)**2))
qdisc=s.discriminant(Q,z);rdisc=s.factor(s.discriminant(R,z));resultant=s.factor(s.resultant(Q,R,z))
assert qdisc!=0 and rdisc==256*e**2*(e+1)
assert resultant==(e**2+48*e+64)**2
assert s.limit(hs2,e,0)==0
out={'claim':'The infimum of the third-order Hilbert–Schmidt defect over the exact partition(2^4,1^4) is zero. This does not imply the infimum is attained.',
 'parameter':'epsilon>0','polynomial':'(z^4-6*z^2+1)^2*((z^2+1)^2+epsilon)',
 'doubleFactorDiscriminant':str(qdisc),'simpleFactorDiscriminant':str(rdisc),'factorResultant':str(resultant),
 'finiteDegree':12,'partitionForEveryPositiveParameter':[2]*4+[1]*4,
 'normalizationSquared':str(norm),'thirdOrderHSDefectSquared':str(hs2),
 'limitPolynomial':str(P.as_expr().subs(e,0)),'limitPartition':[2]*6,'limitExactly3AC':True,
 'limitations':'The construction neither supplies a third-order state in the target partition nor proves that none exists.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
