#!/usr/bin/env python3
"""Independent reduced-matrix and Banach audit of the N11 double-root state."""
from pathlib import Path
import argparse,json,math
import sympy as s
from sympy.polys.matrices import DomainMatrix

parser=argparse.ArgumentParser();parser.add_argument('--witness-dir',type=Path,default=Path(__file__).resolve().parent);args=parser.parse_args()
SOURCE=args.witness_dir/'witness.json'
d=json.loads(SOURCE.read_text());n=11;t=3
x=s.symbols('x0:10',real=True);y=s.symbols('y0:10',real=True);variables=list(x+y)
coeff=[s.sqrt(s.binomial(n,r))*(x[r]+s.I*y[r]) if r<10 else s.Integer(0) for r in range(n+1)]
# Build the actual 3|8 Schmidt matrix in normalized Dicke bases.
M=s.zeros(4,9)
for a in range(4):
 for b in range(9):M[a,b]=s.simplify(coeff[a+b]*s.sqrt(s.binomial(3,a)*s.binomial(8,b)/s.binomial(11,a+b)))
rho=M*M.conjugate().T;equations=[];labels=[]
for a in range(4):
 for b in range(a,4):
  entry=s.expand((rho[a,b]-(s.Rational(1,4) if a==b else 0))/s.sqrt(s.binomial(3,a)*s.binomial(3,b)))
  equations.append(s.expand(s.re(entry)));labels.append([a,b,'real'])
  if a!=b:equations.append(s.expand(s.im(entry)));labels.append([a,b,'imag'])
assert labels==d['equationLabels']
fixed={variables[int(i)]:s.Rational(z) for i,z in d['fixedCoordinates'].items()}
free=[variables[i] for i in d['variableIndices']];assert sorted(d['variableIndices']+list(map(int,d['fixedCoordinates'])))==list(range(20))
assert [str(v) for v in free]==d['variables']
F=s.Matrix(equations).subs(fixed);c=list(map(s.Rational,d['center']));subs=dict(zip(free,c));r=s.Rational(d['radius']);assert len(c)==16
assert all(s.Poly(f,*free).total_degree()==2 for f in F)
allcenter=[fixed[v] if v in fixed else subs[v] for v in variables];assert allcenter==list(map(s.Rational,d['allCoordinateCenter']))
J=F.jacobian(free).subs(subs);R=DomainMatrix.from_Matrix(J).convert_to(s.QQ).inv().to_Matrix();assert R*J==s.eye(16)
e=max(abs(z) for z in R*F.subs(subs))
H=[s.hessian(f,free) for f in F]
bound=max(sum(abs(R[a,i])*sum(abs(z) for z in H[i]) for i in range(16)) for a in range(16))
assert e<s.Rational(1,10**70);assert bound<12000000
assert e+bound*r*r/2<r;assert bound*r<1
assert abs(allcenter[9])>r
result={'independentSchmidtMatrixReconstruction':True,'rho3Equations':16,'variables':16,'fixedRationalCoordinates':4,'jacobianInverseVerifiedExactly':True,'residualUpperBound':'1/10^70','quadraticHessianInverseUpperBound':'12000000','selfMapVerifiedExactly':True,'contractionVerifiedExactly':True,'uniqueRootInSpecified16VariableBox':True,'coefficient9NonzeroThroughoutBox':True,'coefficient10And11IdenticallyZero':True,'maxMultiplicityExactly2ByUniversalOrder3Bound':True}
Path(__file__).with_name('n11-independent-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
