"""Exact algebra behind complete equality classification for N7 S6>=7.

Uses positive-definite blocks in the separately verified global SOS certificate.
The holomorphic cubics span all3x3minors of the3x6 Hankel catalecticant.
Its row kernel is a quadratic apolar differential annihilator, giving two powers
or a tangent form. One mixed cubic covariant excludes nonorthogonal powers
and noncoherent tangents. See equality-classification.txt for the full proof.
"""
from pathlib import Path
from itertools import combinations,permutations,combinations_with_replacement
import sympy as s
import json,sys
ROOT=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT
basis=json.loads((BASE/'basis.json').read_text())
certificate=json.loads((BASE/'sos-certificate.json').read_text())
x=s.symbols('x0:8');y=s.symbols('y0:8')
vars=x+y
# Matrix rows are coefficient vectors of cubic forms in x only.
monomials=list(combinations_with_replacement(range(8),3))
index={m:i for i,m in enumerate(monomials)}
covariants=[]
for b in basis['blocks']:
 if b['kind']!='30':continue
 for channel in b['covariants']:
  for poly in channel:
   row=[s.S(0)]*len(monomials)
   for m,c in poly:row[index[tuple(m)]]+=s.Rational(c)
   covariants.append(row)
columns=list(combinations(range(6),3));minors=[]
for cols in columns:
 row=[s.S(0)]*len(monomials)
 for perm in permutations(range(3)):
  inversions=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
  m=tuple(sorted(i+cols[perm[i]] for i in range(3)))
  row[index[m]]+=(-1)**inversions
 minors.append(row)
C=s.Matrix(covariants);M=s.Matrix(minors)
assert C.shape==M.shape==(20,120)
pivots=M.rref()[1];assert len(pivots)==20
T=C[:,pivots]*M[:,pivots].inv()
assert C==T*M
assert T.det()!=0
# Inspect the decisive mixed cubic block, and independently expand its two-power substitution.
b=next(b for b in basis['blocks'] if b['kind']=='21' and b['L2']==13)
assert b['dimension']==1
polys=[sum(s.Rational(c)*s.prod(vars[i] for i in m) for m,c in poly) for poly in b['covariants'][0]]
h=(x[0]*x[2]-x[1]**2)*y[5]+(x[0]*x[3]-x[1]*x[2])*y[6]+(x[1]*x[3]-x[2]**2)*y[7]
assert s.expand(polys[0]-h)==0
A,B,c,d,Ab,Bb,cb,db=s.symbols('A B c d Ab Bb cb db')
sub={x[r]:A*int(r==0)+B*c**(7-r)*d**r for r in range(8)}
sub.update({y[r]:Ab*int(r==0)+Bb*cb**(7-r)*db**r for r in range(8)})
fact=A*B*Bb*c**4*cb*d**2*db**5*(c*cb+d*db)
assert s.expand(h.subs(sub)-fact)==0
restrict={x[i]:0 for i in range(2,8)}|{y[i]:0 for i in range(2,8)}
assert s.expand(polys[4].subs(restrict)+5*x[1]**2*y[1])==0
# Every block used by the argument has a strictly positive Gram matrix.
positivity=[]
for block in certificate['blocks']:
 if block['kind']!='30' and not(block['kind']=='21' and block['L2']==13):continue
 G=s.Matrix([[s.Rational(v) for v in row] for row in block['matrix']])
 principal=[G[:i,:i].det() for i in range(1,G.rows+1)]
 assert all(v>0 for v in principal)
 positivity.append({'kind':block['kind'],'L2':block['L2'],'leadingPrincipalMinors':[str(v) for v in principal]})
assert len(positivity)==4
result={'n':7,'sector':6,'allChecksExact':True,'holomorphicCubics':20,'hankelMatrixShape':[3,6],
 'hankelEntries':'H_ij=x_(i+j), i=0,1,2; j=0,...,5',
 'holomorphicCubicSpanEqualsHankelMinorSpan':True,'spanDimension':20,
 'basisChangeDeterminant':str(T.det()),'basisChangeRows':[[str(v) for v in row] for row in T.tolist()],
 'minorColumnTriples':[list(c) for c in columns],
 'decisiveMixedBlock':{'kind':'21','L2':13,'component':0,'polynomial':str(h)},
 'twoPowerSubstitution':'x_r=A*delta_(r,0)+B*c^(7-r)*d^r; y=conjugate(x)',
 'twoPowerFactorization':'A*|B|^2*c^4*conjugate(c)*d^2*conjugate(d)^5*(|c|^2+|d|^2)',
 'tangentRestrictionComponent':4,'tangentRestriction':'-5*x_1^2*conjugate(x_1)',
 'neededGramBlocksPositiveDefinite':positivity,
 'classification':'Every normalized equality state is a collective SU(2) rotation of alpha D7^0+beta D7^7, with arbitrary complex alpha,beta satisfying |alpha|^2+|beta|^2=1.',
 'scope':'The necessary zero-covariant conditions use the separately verified global SOS identity. Sufficiency is independently checked in verify_even_sectors.py.'}
(ROOT/'equality-classification-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS rank20 Hankel-minor span, exact mixed-cubic factorization, tangent exclusion, positive Gram blocks.',flush=True)
