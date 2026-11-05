"""Exact existence proof for an N=15, dimension-three, second-order AC subspace.
Run: python3 verify.py. Verification uses rational arithmetic and symbolic identities
exclusively; floating-point arithmetic is used only to check displayed coefficients.
"""
import json,sys,hashlib
from pathlib import Path
import sympy as s
import importlib.util
_spec=importlib.util.spec_from_file_location("n15_algebraic_definition",Path(__file__).with_name("definition.py"))
_def=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_def)
X,E,A,B,C=_def.X,_def.E,_def.A,_def.B,_def.C
BASE=Path(__file__).resolve().parent
OUT=BASE/'certificate.json'
F=s.Matrix(E);J=F.jacobian(X);D=len(X)
H=[s.hessian(e,X) for e in E]
def frac(x):return s.Rational(str(x))
def maxnorm(M):return max(sum(abs(v) for v in M.row(i)) for i in range(M.rows))
def records(e):return [{'powers':list(p),'coefficient':str(c)} for p,c in s.Poly(e,*X).terms()]
def full_cross_verification():
 # Split D_15^r into D_2^i tensor D_13^(r-i). The sqrt(binomial(15,r))
 # in our chosen coordinates cancels the Dicke-splitting denominator.
 # The following polynomials are each cross-reduced entry divided by
 # sqrt(binomial(2,i)*binomial(2,j)).
 T=[]
 for h,V in enumerate([A,B,C]):
  for k,W in enumerate([A,B,C]):
   for i in range(3):
    for j in range(3):
     q=sum(s.binomial(13,l)*V.get(i+l,0)*s.conjugate(W.get(j+l,0)) for l in range(14))
     if h==k and i==j:q-=s.Rational(1,3*s.binomial(2,i))
     T.extend([s.re(q).expand(),s.im(q).expand()])
 monomials=sorted(set(p for e in E+T for p in s.Poly(e,*X).monoms()))
 def row(e):
  d=s.Poly(e,*X).as_dict();return [d.get(p,0) for p in monomials]
 M=s.Matrix([row(e) for e in E]);Q=s.Matrix([row(e) for e in T]);piv=M.rref()[1]
 assert len(piv)==D
 multipliers=Q[:,piv]*M[:,piv].inv();assert multipliers*M==Q
 return {'crossReductionMatrices':9,'complexEntries':81,'realPolynomialIdentities':162,'independentEquations':14,'allCrossReductionsExactlyDeltaIOver3':True,'identityMethod':'Each real and imaginary Dicke-splitting polynomial is an exact rational linear combination of the 14 defining equations.'}
def verify_certificate(data):
 assert data['variables']==[str(v) for v in X]
 assert data['equations']==[records(e) for e in E]
 assert data['n']==15 and data['dimension']==3 and data['t']==2
 for row,V in zip(data['basis'],[A,B,C]):assert row['rootCoordinatePolynomialCoefficients']==[str(V.get(r,0)) for r in range(16)]
 center=s.Matrix(list(map(frac,data['center'])));r=frac(data['radius']);R=s.Matrix([[frac(v) for v in row] for row in data['inverseApproximation']]);subs=dict(zip(X,center));J0=J.subs(subs);F0=F.subs(subs)
 defect=s.eye(D)-R*J0;defectnorm=maxnorm(defect);eta=max(abs(v) for v in R*F0)
 # Hessians are constant. For every x in the radius-r box,
 # ||D[x-RF(x)]||_infinity <= ||I-RJ(center)||_infinity + r*M.
 transformedH=[sum((R[i,l]*H[l] for l in range(D)),s.zeros(D)) for i in range(D)]
 M=max(sum(abs(v) for v in h) for h in transformedH)
 q=defectnorm+r*M
 assert eta+q*r<r and q<1
 assert abs(R.det())>0
 # The full Jacobian is nonsingular throughout the box as ||I-RJ||<1.
 # Thus the isolated real root of the rational polynomial system is algebraic.
 cross=full_cross_verification()
 proof={'centerResidualBound':str(eta),'jacobianApproximationDefect':str(defectnorm),'hessianBound':str(M),'lipschitzBound':str(q),'selfMapSlack':str(r-eta-q*r),'radius':str(r),'contraction':True,'mapsClosedBoxStrictlyIntoItself':True,'uniqueRootInBox':True,'nonsingularJacobianInBox':True,'rootIsAlgebraic':True,'fullCrossVerification':cross,'summaryBounds':{'centerResidualLessThan':'1e-59','jacobianDefectLessThan':'1e-59','hessianBoundLessThan':'1e4','lipschitzLessThan':'1e-40'}}
 assert eta<s.Rational(1,10**59);assert defectnorm<s.Rational(1,10**59);assert M<10**4;assert q<s.Rational(1,10**40)
 return proof
def verify():
 data=json.loads(OUT.read_text());proof=verify_certificate(data)
 assert proof==data['certificate']
 basis=json.loads((BASE/'basis.json').read_text())
 assert basis['id']=='algebraic-n15-k3-t2' and basis['n']==15 and basis['t']==2 and basis['dimension']==3
 assert len(basis['basis'])==3 and len(data['basis'])==3
 assert basis['algebraicDefinition']['variables']==[str(v) for v in X]
 assert basis['algebraicDefinition']['equations']==[str(e) for e in E]
 assert basis['algebraicDefinition']['equationsLatex']==[s.latex(e)+'=0' for e in E]
 center=s.Matrix(list(map(frac,data['center'])));radius=frac(data['radius']);sub=dict(zip(X,center))
 for i,x in enumerate(X):assert basis['algebraicDefinition']['parameterIntervals'][str(x)]==[str(center[i]-radius),str(center[i]+radius)]
 errors=[]
 for row,V in zip(basis['basis'],[A,B,C]):
  expressions=[s.sqrt(s.binomial(15,r))*V.get(r,0) for r in range(16)]
  assert row['exactCoefficients']==[str(e) for e in expressions]
  assert row['coefficientLatex']==[s.latex(e) for e in expressions]
  for r,e in enumerate(expressions):
   pair=[float(s.re(e.subs(sub))),float(s.im(e.subs(sub)))];errors.extend(abs(pair[k]-row['coefficients'][r][k]) for k in range(2))
 assert max(errors)<1e-15
 report={'passed':True,'n':15,'t':2,'dimension':3,'id':basis['id'],'exactnessKind':basis['exactnessKind'],'coefficientRenderingError':max(errors),'certificate':proof}
 return report
if __name__=='__main__':
 report=verify();(BASE/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'passed':report['passed'],'id':report['id'],'coefficientRenderingError':report['coefficientRenderingError'],'boundSummary':report['certificate']['summaryBounds']},indent=2))
