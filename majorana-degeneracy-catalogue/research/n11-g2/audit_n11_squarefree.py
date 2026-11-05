#!/usr/bin/env python3
"""Independent polynomial-map construction of the Sylvester certificate."""
from pathlib import Path
import argparse,json,math
import sympy as s
from sympy.polys.matrices import DomainMatrix

parser=argparse.ArgumentParser();parser.add_argument('--witness-dir',type=Path,default=Path(__file__).resolve().parent);args=parser.parse_args()
base=args.witness_dir
d=json.loads((base/'witness.json').read_text());cert=json.loads((base/'majorana-certificate.json').read_text())
v=list(map(s.Rational,d['allCoordinateCenter']));r=s.Rational(d['radius']);z=s.symbols('z')
p=s.Poly(sum((-1)**k*s.binomial(11,k)*(v[k]+s.I*v[k+10])*z**k for k in range(10)),z)
assert p.degree()==9 and abs(v[9])>r
dp=p.diff();images=[p.mul(s.Poly(z**a,z)) for a in range(7,-1,-1)]+[dp.mul(s.Poly(z**a,z)) for a in range(8,-1,-1)]
S=s.Matrix([[pol.nth(k) for k in range(16,-1,-1)] for pol in images])
R=s.Matrix([[s.Rational(a)+s.I*s.Rational(b) for a,b in row] for row in cert['inverseApproximation']])
RD=DomainMatrix.from_Matrix(R).convert_to(s.QQ_I);SD=DomainMatrix.from_Matrix(S).convert_to(s.QQ_I)
E=s.eye(17)-(RD*SD).to_Matrix()
norm=lambda A:max(sum(abs(s.re(A[i,j]))+abs(s.im(A[i,j])) for j in range(A.cols)) for i in range(A.rows))
res=norm(E);invnorm=norm(R)
# Each free real or imaginary coordinate changes by at most radius. Fixed
# coordinates are conservatively allowed the same variation here.
polybound=2*r*sum(s.binomial(11,k) for k in range(10))
derivativebound=2*r*sum(k*s.binomial(11,k) for k in range(10))
variation=max(polybound,derivativebound);bound=res+invnorm*variation
assert bound<s.Rational(6,10**22)<1
result={'independentPolynomialMapConstruction':True,'sylvesterSize':17,'degree9LeadingCoefficientNonzero':True,'exactRationalNeumannBoundLessThan':'6/10^22','finitePolynomialSquarefreeThroughoutBox':True,'infinityMultiplicity':2,'multiplicityPartition':[2]+[1]*9,'floatingDiagnosticOnly':float(bound)}
Path(__file__).with_name('n11-independent-squarefree-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
