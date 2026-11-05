"""Independent exact algebra audit of the complete N=7 S6=7 equality family.
Standard library only. Input basis.json in this folder or command-line folder.
The global SOS certificate is checked separately. This audit checks:
(1) its20holomorphiccubic components span all20Hankel3x3 minors;
(2) the decisive mixed cubic's two-power factorization and tangent exclusion;
(3) the coefficient identity for the quadratic apolar differential operator.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from itertools import combinations,permutations
from math import comb
import json,sys,hashlib
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
basis=json.loads((BASE/'basis.json').read_text())
def clean(p):return {m:c for m,c in p.items() if c}
def plus(*ps):
 q=defaultdict(F)
 for p in ps:
  for m,c in p.items():q[m]+=c
 return clean(q)
def mul(p,q):
 out=defaultdict(F)
 for m,c in p.items():
  for n,d in q.items():out[tuple(sorted(m+n))]+=c*d
 return clean(out)
def rank(rows):
 a=[list(map(F,r)) for r in rows];i=0
 for j in range(len(a[0])):
  p=next((p for p in range(i,len(a)) if a[p][j]),None)
  if p is None:continue
  a[i],a[p]=a[p],a[i];v=a[i][j];a[i]=[x/v for x in a[i]]
  for k in range(i+1,len(a)):
   if a[k][j]:
    v=a[k][j];a[k]=[x-v*y for x,y in zip(a[k],a[i])]
  i+=1
  if i==len(a):break
 return i
holo=[]
for b in basis['blocks']:
 if b['kind']=='30':
  assert b['dimension']==1
  holo.extend({tuple(mon):F(c) for mon,c in p} for p in b['covariants'][0])
assert len(holo)==20
minors=[]
for cols in combinations(range(6),3):
 p=defaultdict(F)
 for perm in permutations(range(3)):
  inversions=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
  mon=tuple(sorted(row+cols[perm[row]] for row in range(3)))
  p[mon]+=(-1)**inversions
 minors.append(clean(p))
mons=sorted(set().union(*map(set,holo+minors)))
rh=rank([[p.get(m,0) for m in mons] for p in holo])
rm=rank([[p.get(m,0) for m in mons] for p in minors])
rall=rank([[p.get(m,0) for m in mons] for p in holo+minors])
assert (rh,rm,rall)==(20,20,20)
print('PASS20cubiccomponents span20Hankel3x3minors.',flush=True)
B=next(b for b in basis['blocks'] if b['kind']=='21' and b['L2']==13)
assert B['dimension']==1
mixed=[{tuple(m):F(c) for m,c in p} for p in B['covariants'][0]]
expected={(0,2,13):F(1),(1,1,13):F(-1),(0,3,14):F(1),(1,2,14):F(-1),(1,3,15):F(1),(2,2,15):F(-1)}
assert mixed[0]==expected
# Variables A,B,c,d,Abar,Bbar,cbar,dbar indexed0,...,7.
# xr=A delta(r,0)+B c^(7-r)d^r; yr=formal conjugate.
sub=[]
for r in range(8):
 sub.append(plus({tuple([1]+[2]*(7-r)+[3]*r):F(1)}, {(0,):F(1)} if r==0 else {}))
for r in range(8):
 sub.append(plus({tuple([5]+[6]*(7-r)+[7]*r):F(1)}, {(4,):F(1)} if r==0 else {}))
value=defaultdict(F)
for mon,coef in mixed[0].items():
 term={():coef}
 for i in mon:term=mul(term,sub[i])
 for m,c in term.items():value[m]+=c
# A B Bbar c^4 cbar d^2 dbar^5 (c cbar+d dbar).
pref=tuple(sorted([0,1,5]+[2]*4+[6]+[3]*2+[7]*5))
expected_factor=mul({pref:F(1)},{(2,6):F(1),(3,7):F(1)})
assert clean(value)==expected_factor
restricted={m:c for m,c in mixed[4].items() if all(i%8 in (0,1) for i in m)}
assert restricted=={(1,1,9):F(-5)}
print('PASSexacttwo-powerfactorization and tangent−5x1²conj(x1).',flush=True)
# Coefficients of q0*dXX+q1*dXY+q2*dYY applied to f.
for j in range(6):
 expected=42*comb(5,j)
 assert comb(7,j)*(7-j)*(6-j)==expected
 assert comb(7,j+1)*(6-j)*(j+1)==expected
 assert comb(7,j+2)*(j+2)*(j+1)==expected
print('PASScatalecticant/apolar differential identity.',flush=True)
result={'n':7,'sector':6,'minimum':'7','holomorphicCubicCount':20,'hankelMinorCount':20,'holomorphicRank':rh,'hankelMinorRank':rm,'combinedRank':rall,'hankelConvention':'H[i,j]=x_(i+j), i=0,1,2 and j=0,...,5','quadraticApolarIdentityVerified':True,'apolarIdentity':'(q0 dXX+q1 dXY+q2 dYY)f=42 sum_j C(5,j)(q0 x_j+q1 x_(j+1)+q2 x_(j+2))X^(5-j)Y^j','mixedHighestIdentityVerified':True,'mixedHighest':'(x0*x2-x1^2)y5+(x0*x3-x1*x2)y6+(x1*x3-x2^2)y7','twoPowerSubstitution':'x_r=A delta(r,0)+B c^(7-r)d^r, y_r=conj(x_r)','twoPowerFactorization':'A B conj(B) c^4 conj(c) d^2 conj(d)^5 (|c|^2+|d|^2)','twoPowerFactorizationVerified':True,'tangentComponentIndex':4,'tangentComponent':'-5*x1^2*conj(x1)','tangentExclusionIdentityVerified':True,'classification':'Every normalized equality state is U^(tensor7)(alpha|0>^(tensor7)+beta|1>^(tensor7)), with U in SU(2) and |alpha|^2+|beta|^2=1. The proof additionally uses factorization of a nonzero quadratic over C and the homogeneous PDE normal forms described in the report.','basisSha256':hashlib.sha256((BASE/'basis.json').read_bytes()).hexdigest()}
(BASE/'classification-independent.json').write_text(json.dumps(result,indent=2)+'\n')
