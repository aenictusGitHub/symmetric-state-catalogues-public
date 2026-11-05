"""Independently derive bidegree-(4,4) six-qubit S4 target via RDMs.

Write c_r = sqrt(binomial(6,r))*x_r, and y_r=conjugate(x_r).
Q=sum binomial(6,r)*x_r*y_r is the squared norm.
The unnormalized reduced t-qubit matrix is
 rho_t[i,k]=sqrt(C(t,i)C(t,k))*sum_j C(6-t,j)*x_(i+j)*y_(k+j).
Thus Tr(rho_t^2) is a rational bidegree-(2,2) polynomial.
S4_un=15*sum_t=0..4 (-1)^(4-t)*C(4,t)*2^t*Tr(rho_t^2).
Target Q^2*S4_un-(45/7)*Q^4 equals Q^4*(S4_normalized-45/7).
No spin-coupling/SOS conventions are used in this derivation.
"""
from fractions import Fraction as F
from math import comb
from collections import defaultdict
from pathlib import Path
import json,hashlib
N=6;D=N+1
zero=(0,)*D
def exp(*indices):
 a=[0]*D
 for i in indices:a[i]+=1
 return tuple(a)
def mul(A,B):
 C=defaultdict(F)
 for (a,b),v in A.items():
  for (c,d),w in B.items():
   C[(tuple(x+y for x,y in zip(a,c)),tuple(x+y for x,y in zip(b,d)))]+=v*w
 return {k:v for k,v in C.items() if v}
def serialize(P):
 return [[list(a),list(b),str(v)] for (a,b),v in sorted(P.items())]
Q={(exp(i),exp(i)):F(comb(N,i)) for i in range(D)}
S=defaultdict(F)
purities=[]
for t in range(5):
 P=defaultdict(F)
 for i in range(t+1):
  for k in range(t+1):
   for j in range(N-t+1):
    for l in range(N-t+1):
     P[(exp(i+j,k+l),exp(k+j,i+l))]+=comb(t,i)*comb(t,k)*comb(N-t,j)*comb(N-t,l)
 purities.append(dict(P))
 weight=comb(N,4)*(-1)**(4-t)*comb(4,t)*2**t
 for ab,v in P.items():S[ab]+=weight*v
S={k:v for k,v in S.items() if v}
Q2=mul(Q,Q);Q4=mul(Q2,Q2)
T=defaultdict(F,mul(Q2,S))
for ab,v in Q4.items():T[ab]-=F(45,7)*v
T={k:v for k,v in T.items() if v}
for (a,b),v in T.items():
 assert sum(a)==sum(b)==4
 assert sum(i*e for i,e in enumerate(a))==sum(i*e for i,e in enumerate(b))
 assert T.get((b,a))==v
# Exact checks of sparse families use rational squared amplitudes. Equal degree
# and z-rotation weight force a=b on each two-index support.
for support,xsq in [((2,5),{2:F(4,63),5:F(1,126)}),((0,5),{0:F(4,7),5:F(1,14)})]:
 value=F(0);q=F(0)
 for i,prob in xsq.items():q+=comb(N,i)*prob
 assert q==1
 for (a,b),v in T.items():
  if any((a[i] or b[i]) and i not in support for i in range(D)):continue
  assert a==b
  term=v
  for i in support:term*=xsq[i]**a[i]
  value+=term
 assert value==0
# Coherent pole checks: S4=15, hence target=60/7 for Q=1.
assert T[(exp(0,0,0,0),exp(0,0,0,0))]==F(60,7)
assert T[(exp(6,6,6,6),exp(6,6,6,6))]==F(60,7)
obj={'n':N,'sector':4,'bound':'45/7','bidegree':[4,4],
 'coordinateConvention':'c_r=sqrt(binomial(6,r))*x_r; y_r=conjugate(x_r); r=0,...,6',
 'normConvention':'Q=sum_r binomial(6,r)*x_r*y_r',
 'targetConvention':'Q^2*S4_un-(45/7)*Q^4 = Q^4*(S4_normalized-45/7)',
 'coefficientEncoding':'Each row [alpha,beta,rational] means rational*x^alpha*y^beta.',
 'termCount':len(T),'coefficients':serialize(T),
 'unnormalizedS4':serialize(S),
 'checks':{'bidegrees':True,'zRotationWeights':True,'conjugateSymmetry':True,'bothExactAttainersVanish':True,'coherentPoles':True}}
root=Path(__file__).resolve().parent
(root/'target-exact.json').write_text(json.dumps(obj,indent=2)+'\n')
print('Wrote independent exact target:',len(T),'terms; unnormalized S4:',len(S),'terms.',flush=True)
print('All algebraic sanity checks passed.',flush=True)
