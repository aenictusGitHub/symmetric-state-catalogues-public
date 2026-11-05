#!/usr/bin/env python3
"""Exact existence certificate for an algebraic n11,t3 state with g=2.
Requires SymPy only. It proves a root exists using rational Banach bounds;
printed decimal amplitudes are not treated as an exact solution.
"""
import json,math
from pathlib import Path
from fractions import Fraction
import sympy as s
from sympy.polys.matrices import DomainMatrix
HERE=Path(__file__).resolve().parent

def equations():
 xs=s.symbols('x0:10',real=True);ys=s.symbols('y0:10',real=True);x=list(xs)+[s.Integer(0)]*2;y=list(ys)+[s.Integer(0)]*2
 fs=[];labels=[]
 for a in range(4):
  for b in range(a,4):
   re=sum(s.binomial(8,q)*(x[a+q]*x[b+q]+y[a+q]*y[b+q]) for q in range(9))
   fs.append(s.expand(re-(s.Rational(1,4*s.binomial(3,a)) if a==b else 0)));labels.append([a,b,'real'])
   if a!=b:
    fs.append(s.expand(sum(s.binomial(8,q)*(y[a+q]*x[b+q]-x[a+q]*y[b+q]) for q in range(9))));labels.append([a,b,'imag'])
 return list(xs+ys),s.Matrix(fs),labels

def verify():
 d=json.loads((HERE/'witness.json').read_text());allvars,F,labels=equations();assert d['n']==11 and d['t']==3 and d['equationLabels']==labels
 selected=d['variableIndices'];fixed={int(i):s.Rational(v) for i,v in d['fixedCoordinates'].items()};assert sorted(selected+list(fixed))==list(range(20)) and len(selected)==16 and len(fixed)==4
 variables=[allvars[i] for i in selected];assert list(map(str,variables))==d['variables'];eq=F.subs({allvars[i]:v for i,v in fixed.items()});center=[s.Rational(x) for x in d['center']];assert len(center)==16
 subs=dict(zip(variables,center));radius=s.Rational(d['radius']);assert radius==s.Rational(1,10**30)
 full=[fixed.get(i,subs.get(allvars[i])) for i in range(20)];assert full==[s.Rational(x) for x in d['allCoordinateCenter']]
 for f in eq:assert s.Poly(f,*variables).total_degree()<=2
 J=eq.jacobian(variables).subs(subs);R=DomainMatrix.from_Matrix(J).convert_to(s.QQ).inv().to_Matrix();assert R*J==s.eye(16)
 e=max(abs(a) for a in R*eq.subs(subs));H=[s.hessian(f,variables) for f in eq];hs=[sum(abs(v) for v in h) for h in H];M=max(sum(abs(R[a,i])*hs[i] for i in range(16)) for a in range(16))
 assert e<s.Rational(1,10**70) and M<12000000
 assert e+M*radius**2<radius and M*radius<1
 # c9 is nonzero throughout the certified box; c10=c11 are identically zero.
 assert abs(full[9])>radius
 assert all(max(abs(full[k]),abs(full[k+10]))>radius for k in range(10))
 result={'verified':True,'n':11,'t':3,'method':'Exact rational Banach fixed-point certificate for the sixteen real reduced-density-matrix equations','uniqueRealAlgebraicRootInBox':True,'radius':'1/10^30','residualBound':'1/10^70','HessianInverseBound':'12000000','contractionFactorBound':'12/10^24','exactDoublePointAtInfinity':True,'normalization':'rho3=I4/4, hence norm1 exactly'}
 print(json.dumps(result,indent=2));return result
if __name__=='__main__':verify()
