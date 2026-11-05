"""Exact odd-seven-qubit even-sector identity and complex GHZ-family values.

Coordinates c_r=sqrt(C(7,r))*x_r; y_r=conjugate(x_r);
q=sum C(7,r)x_ry_r. Each P_k=q^2*S_k(normalized).
The direct RDM check proves P2+P4+P6=63q^2.
On support{0,7}, Pk=C(7,k)q^2 for k=0,2,4,6.
Thus the global bound S6>=7 is equivalent to S2+S4<=56;
this script checks the identity/equality family, not the global bound itself.
"""
from fractions import Fraction as F
from math import comb
from collections import defaultdict
from pathlib import Path
import json
N=7;D=8
def exp(*indices):
 a=[0]*D
 for i in indices:a[i]+=1
 return tuple(a)
def mul(A,B):
 C=defaultdict(F)
 for (a,b),v in A.items():
  for (c,d),w in B.items():C[(tuple(x+y for x,y in zip(a,c)),tuple(x+y for x,y in zip(b,d)))]+=v*w
 return {k:v for k,v in C.items() if v}
def add(*pairs):
 C=defaultdict(F)
 for coeff,P in pairs:
  for ab,v in P.items():C[ab]+=coeff*v
 return {k:v for k,v in C.items() if v}
q={(exp(i),exp(i)):F(comb(N,i)) for i in range(D)};q2=mul(q,q)
purities=[]
for t in range(7):
 P=defaultdict(F)
 for i in range(t+1):
  for k in range(t+1):
   for j in range(N-t+1):
    for l in range(N-t+1):P[(exp(i+j,k+l),exp(k+j,i+l))]+=comb(t,i)*comb(t,k)*comb(N-t,j)*comb(N-t,l)
 purities.append(dict(P))
sectors={k:add(*[(comb(N,k)*(-1)**(k-t)*comb(k,t)*2**t,purities[t]) for t in range(k+1)]) for k in [0,2,4,6]}
assert sectors[0]==q2
assert not add((1,sectors[2]),(1,sectors[4]),(1,sectors[6]),(-63,q2))
def restrict(P):return {ab:v for ab,v in P.items() if all(ab[0][i]==ab[1][i]==0 for i in range(1,7))}
for k in [0,2,4,6]:assert restrict(sectors[k])==add((comb(N,k),restrict(q2)))
obj={'n':7,'allComplexAmplitudes':True,'floatingPointUsed':False,
 'coordinateConvention':'c_r=sqrt(binomial(7,r))*x_r; y_r=conjugate(x_r)',
 'normConvention':'q=sum_r binomial(7,r)*x_r*y_r',
 'identity':'P2+P4+P6=63*q^2','normalizedIdentity':'S2+S4+S6=63',
 'identityVerifiedExactly':True,'nonzeroResidualCoefficients':0,
 'equalityFamily':'alpha D7^0+beta D7^7, |alpha|^2+|beta|^2=1, arbitrary complex alpha,beta',
 'evenSectorVector':{'S0':'1','S2':'21','S4':'35','S6':'7'},
 'familyValuesVerifiedExactly':True,
 'equivalentSharpTradeoff':'S2+S4<=56 iff S6>=7',
 'globalBoundProvedByThisScript':False,
 'note':'The separate global SOS certificate establishes the bound; these identities make the tradeoff a sharp corollary. No full equality-state classification is asserted.'}
Path(__file__).with_name('even-sector-verification.json').write_text(json.dumps(obj,indent=2)+'\n')
print(json.dumps(obj,indent=2))
