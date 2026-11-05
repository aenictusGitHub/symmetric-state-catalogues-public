"""Independent exact check of the spin-moment identity for six qubits.

Dicke convention: Jz |D_6^(r)>=(3-r)|D_6^(r)>.
Write c_r=sqrt(binomial(6,r))*x_r and y_r=conjugate(x_r).
q=<psi|psi>; F_i=<psi|J_i|psi>;
Q_ij=<psi|{J_i,J_j}/2|psi>-4*q*delta_ij;
theta=sum_r (-1)^r c_r*c_(6-r).
Then the unnormalized fourth sector is exactly
11*q^2+(8/3)*sum_i F_i^2-(8/15)*sum_ij Q_ij*Q_ji+24*theta*conjugate(theta).
The comparison below derives S4 independently from reduced purities.
"""
from pathlib import Path
from math import comb
from collections import defaultdict
import json
import sympy as s
N=6;D=7
zero=(0,)*D
def exp(*indices):
 a=[0]*D
 for i in indices:a[i]+=1
 return tuple(a)
def add(A,B,scale=1):
 C=dict(A)
 for key,v in B.items():C[key]=s.expand(C.get(key,0)+scale*v)
 return {k:v for k,v in C.items() if v!=0}
def mul(A,B):
 C=defaultdict(lambda:s.S(0))
 for (a,b),v in A.items():
  for (c,d),w in B.items():
   C[(tuple(x+y for x,y in zip(a,c)),tuple(x+y for x,y in zip(b,d)))]+=v*w
 return {k:s.expand(v) for k,v in C.items() if s.expand(v)!=0}
def bilinear(M):
 return {(exp(j),exp(i)):M[i,j] for i in range(D) for j in range(D) if M[i,j]!=0}
# These are D^-1 J D with D_rr=sqrt(C(6,r)); they have rational entries.
W=s.diag(*[comb(N,r) for r in range(D)])
Ap=s.zeros(D)
Am=s.zeros(D)
for r in range(1,D):Ap[r-1,r]=N+1-r
for r in range(N):Am[r+1,r]=r+1
Az=s.diag(*[3-r for r in range(D)])
As=[(Ap+Am)/2,(Ap-Am)/(2*s.I),Az]
assert As[0]*As[1]-As[1]*As[0]==s.I*As[2]
assert sum((A*A for A in As),s.zeros(D))==12*s.eye(D)
q=bilinear(W)
q2=mul(q,q)
f=[bilinear(W*A) for A in As]
F2={}
for Fi in f:F2=add(F2,mul(Fi,Fi))
quad=[[bilinear(W*((Ai*Aj+Aj*Ai)/2-4*s.eye(D)*int(i==j))) for j,Aj in enumerate(As)] for i,Ai in enumerate(As)]
assert not add(add(quad[0][0],quad[1][1]),quad[2][2])
trQ2={}
for i in range(3):
 for j in range(3):trQ2=add(trQ2,mul(quad[i][j],quad[j][i]))
theta=defaultdict(lambda:s.S(0))
for r in range(D):theta[(exp(r,N-r),zero)]+=(-1)**r*comb(N,r)
thetabar={(b,a):v for (a,b),v in theta.items()}
theta2=mul(theta,thetabar)
moment=add(add(add({},q2,11),F2,s.Rational(8,3)),trQ2,-s.Rational(8,15))
moment=add(moment,theta2,24)
# Independent RDM purity expansion, directly in binomial-scaled coordinates.
rdm=defaultdict(lambda:s.S(0))
for t in range(5):
 weight=comb(N,4)*(-1)**(4-t)*comb(4,t)*2**t
 for i in range(t+1):
  for k in range(t+1):
   for j in range(N-t+1):
    for l in range(N-t+1):
     rdm[(exp(i+j,k+l),exp(k+j,i+l))]+=weight*comb(t,i)*comb(t,k)*comb(N-t,j)*comb(N-t,l)
rdm={k:v for k,v in rdm.items() if v}
residual=add(moment,rdm,-1)
assert not residual
assert all(v.is_Rational for v in moment.values())
result={'n':6,'sector':4,'checkedExactly':True,'rdmTermCount':len(rdm),'momentTermCount':len(moment),'nonzeroResidualCoefficients':len(residual),
 'JzConvention':'J_z |D_6^(r)>=(3-r)|D_6^(r)>',
 'q':'sum_r |c_r|^2', 'F_i':'<psi|J_i|psi>',
 'Q_ij':'<psi|(J_i J_j+J_j J_i)/2|psi>-4*q*delta_ij',
 'theta':'sum_r (-1)^r c_r*c_(6-r)',
 'identity':'S4_un=11*q^2+(8/3)*|F|^2-(8/15)*tr(Q^2)+24*|theta|^2',
 'equivalentSharpInequality':'tr(Q^2)<=(60/7)*q^2+5*|F|^2+45*|theta|^2',
 'identityProved':True,'sharpInequalityProvedByThisScript':False,
 'note':'The identity alone does not prove the sharp inequality; a separate global certificate is required.'}
Path(__file__).with_name('spin-moment-check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
