"""Independent bidegree-(3,3) target for seven-qubit S6>=7.

c_r=sqrt(binomial(7,r))*x_r, y_r=conjugate(x_r),
q=sum_r binomial(7,r)*x_r*y_r=<psi|psi>.
Unnormalized reduced matrix:
 rho_t[a,b]=sqrt(C(t,a)C(t,b))*sum_j C(7-t,j)*x_(a+j)*y_(b+j).
P6=7*sum_t=0..6 (-1)^(6-t)*C(6,t)*2^t*tr(rho_t^2).
The target q*P6-7*q^3 equals q^3*(S6(normalized)-7).
Everything is derived in exact rational arithmetic, without SOS/CG input.
"""
from fractions import Fraction as F
from math import comb
from collections import defaultdict
from pathlib import Path
import json
N=7;D=N+1
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
def serialize(P):return [[list(a),list(b),str(v)] for (a,b),v in sorted(P.items())]
q={(exp(i),exp(i)):F(comb(N,i)) for i in range(D)}
P6=defaultdict(F)
for t in range(7):
 weight=comb(N,6)*(-1)**(6-t)*comb(6,t)*2**t
 for i in range(t+1):
  for k in range(t+1):
   for j in range(N-t+1):
    for l in range(N-t+1):
     P6[(exp(i+j,k+l),exp(k+j,i+l))]+=weight*comb(t,i)*comb(t,k)*comb(N-t,j)*comb(N-t,l)
P6={k:v for k,v in P6.items() if v}
q2=mul(q,q);q3=mul(q2,q)
T=defaultdict(F,mul(q,P6))
for ab,v in q3.items():T[ab]-=7*v
T={k:v for k,v in T.items() if v}
for (a,b),v in T.items():
 assert sum(a)==sum(b)==3
 assert sum(i*e for i,e in enumerate(a))==sum(i*e for i,e in enumerate(b))
 assert T.get((b,a))==v
# The entire complex endpoint family alpha D7^0+beta D7^7 attains equality:
# all target monomials surviving restriction to this support must vanish.
assert not any(all((a[i]==b[i]==0) for i in range(1,7)) for a,b in T)
# Independent single-excitation Dicke value S6=145/7, with q=7 when x1=1.
assert T[(exp(1,1,1),exp(1,1,1))]==F(4704)
obj={'n':N,'sector':6,'bound':'7','bidegree':[3,3],
 'coordinateConvention':'c_r=sqrt(binomial(7,r))*x_r; y_r=conjugate(x_r); r=0,...,7',
 'normConvention':'q=sum_r binomial(7,r)*x_r*y_r',
 'targetConvention':'q*P6-7*q^3=q^3*(S6_normalized-7)',
 'coefficientEncoding':'Each row [alpha,beta,rational] means rational*x^alpha*y^beta.',
 'termCount':len(T),'coefficients':serialize(T),'unnormalizedS6':serialize(P6),
 'checks':{'bidegrees':True,'zRotationWeights':True,'conjugateSymmetry':True,'complexEndpointFamilyVanishes':True,'singleExcitationDickeValue':True}}
Path(__file__).with_name('target-exact.json').write_text(json.dumps(obj,indent=2)+'\n')
print('Independent exact target:',len(T),'terms; P6:',len(P6),'terms. All algebraic sanity checks passed.',flush=True)
