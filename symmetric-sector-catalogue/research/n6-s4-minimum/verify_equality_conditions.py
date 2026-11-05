"""Exact equality conditions for the N=6, S4>=45/7 SOS certificate.
Standard library only. Use beside basis.json and sos-certificate.json, or pass
that folder as first argument. Also independently verifies the spin-moment
identity by comparison with directly expanded reduced Dicke purities.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import comb
import json,sys
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
basis=json.loads((BASE/'basis.json').read_text());certificate=json.loads((BASE/'sos-certificate.json').read_text())
def mul(p,q):
 out=defaultdict(F)
 for a,x in p.items():
  for b,y in q.items():out[tuple(sorted(a+b))]+=x*y
 return {m:x for m,x in out.items() if x}
def add(*terms):
 out=defaultdict(F)
 for coef,p in terms:
  for m,x in p.items():out[m]+=coef*x
 return {m:x for m,x in out.items() if x}
def conj(p):return {tuple(sorted((i+7)%14 for i in m)):v for m,v in p.items()}
def purity(t):
 p=defaultdict(F)
 for i in range(t+1):
  for j in range(t+1):
   for a in range(7-t):
    for b in range(7-t):p[tuple(sorted((i+a,7+j+a,j+b,7+i+b)))]+=comb(t,i)*comb(t,j)*comb(6-t,a)*comb(6-t,b)
 return dict(p)
q={(r,7+r):F(comb(6,r)) for r in range(7)};q2=mul(q,q)
fz={(r,7+r):F((3-r)*comb(6,r)) for r in range(7) if r!=3}
fp={(r+1,7+r):F((6-r)*comb(6,r)) for r in range(6)}
f2=add((1,mul(fz,fz)),(1,mul(fp,conj(fp))))
theta=defaultdict(F)
for r in range(7):theta[tuple(sorted((r,6-r)))]+=(-1)**r*comb(6,r)
theta2=mul(theta,conj(theta))
qzz={(r,7+r):F(((3-r)**2-4)*comb(6,r)) for r in range(7)}
jpp={(r+2,7+r):F((6-r)*(5-r)*comb(6,r)) for r in range(5)}
zp={(r+1,7+r):F((5-2*r)*(6-r)*comb(6,r)) for r in range(6)}
quad=add((F(3,2),mul(qzz,qzz)),(F(1,2),mul(jpp,conj(jpp))),(F(1,2),mul(zp,conj(zp))))
s4=add(*[(F(15*(-1)**(4-t)*comb(4,t)*2**t),purity(t)) for t in range(5)])
assert s4==add((11,q2),(F(8,3),f2),(F(-8,15),quad),(24,theta2))
# Each row expresses an unscaled scalar covariant in the specified four invariants.
rows=[[F(0),F(0),F(1,4),F(0)], [F(-3),F(-11,12),F(-5),F(3,10)], [F(49,5),F(-1309,20),F(-769,5),F(469,150)]]
block=next(b for b in basis['blocks'] if b['kind']=='22' and b['L']==0)
assert block['dimension']==3
invariants=[q2,f2,theta2,quad]
for row,cov in zip(rows,block['covariants']):
 actual={tuple(m):F(c) for m,c in cov[0]}
 assert actual==add(*zip(row,invariants))

def det(M):
 if len(M)==1:return M[0][0]
 return sum((-1)**j*M[0][j]*det([[M[r][k] for k in range(len(M)) if k!=j] for r in range(1,len(M))]) for j in range(len(M)))
minor=det([row[1:] for row in rows]);assert minor!=0
kernel=[F(1),F(36,49),F(0),F(600,49)]
assert all(sum(a*b for a,b in zip(row,kernel))==0 for row in rows)
gram=next(b for b in certificate['blocks'] if b['kind']=='22' and b['L']==0)
G=[[F(c) for c in row] for row in gram['matrix']]
assert all(det([row[:i] for row in G[:i]])>0 for i in range(1,4))
assert F(11)+F(8,3)*kernel[1]-F(8,15)*kernel[3]+24*kernel[2]==F(45,7)
output={'n':6,'sector':4,'minimum':'45/7','invariantOrder':['q^2','|F|^2','|theta|^2','tr(Q^2)'],'scalarCovariantRows':[[str(v) for v in row] for row in rows],'coefficientMinorDeterminant':str(minor),'rank':3,'normalizedNullVector':[str(v) for v in kernel],'spinMomentIdentityVerifiedAgainstRDM':True,'scalarCovariantIdentitiesVerified':True,'scalarGramPositiveDefinite':True,'equalityNecessaryUsingGlobalSOSCertificate':True,'equalitySufficientUsingMomentIdentity':True,'equalityCondition':'For normalized states: theta=0, |<J>|^2=36/49, tr(Q^2)=600/49, where Qij=<{Ji,Jj}/2>-4 deltaij. This does not classify states up to rotations.'}
(BASE/'equality-verification.json').write_text(json.dumps(output,indent=2)+'\n')
print('PASS exact moment identity, three scalar covariant identities, rank3 and equality iff conditions.',flush=True)
