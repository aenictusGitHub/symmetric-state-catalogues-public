#!/usr/bin/env python3
"""Exact phase-family and nonextension checks; no global k=3 exclusion."""
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parent

def verify():
    z=s.Symbol('z', nonzero=True)
    def adj(A): return A.T.applyfunc(lambda x:s.conjugate(x).subs(s.conjugate(z),1/z))
    def clean(A): return A.applyfunc(s.simplify)
    p=s.Rational(5,18);q=s.Rational(13,18)
    u=s.zeros(14,1);v=s.zeros(14,1)
    u[0]=s.sqrt(p);u[9]=s.sqrt(q)
    v[4]=s.sqrt(q);v[13]=z*s.sqrt(p)
    V=s.Matrix.hstack(u,v)
    assert clean(adj(V)*V)==s.eye(2)
    def split(w):return s.Matrix(3,12,lambda a,b:w[a+b]*s.sqrt(s.binomial(2,a)*s.binomial(11,b)/s.binomial(13,a+b)))
    splits=list(map(split,[u,v]))
    for a in range(2):
        for b in range(2):
            assert clean(splits[a]*adj(splits[b]))==s.eye(3)*s.Rational(int(a==b),3)
    Jz=s.diag(*[s.Rational(13,2)-r for r in range(14)])
    Jm=s.zeros(14)
    for r in range(13):Jm[r+1,r]=s.sqrt((r+1)*(13-r))
    Jp=Jm.T
    # These fourteen columns lie in the rank <=2 tensor images of V.
    M=s.Matrix.hstack(u,Jz*u,v,Jz*v,Jp*u,Jp**2*u,Jm*u,
          (Jz*Jm+Jm*Jz)*u,Jm*v,Jm**2*v,Jp*v,
          (Jz*Jp+Jp*Jz)*v,Jm**2*u,Jp**2*v)
    determinant=s.factor(M.det())
    constant=s.simplify(determinant/(z**2*(z-11)))
    assert constant!=0 and not constant.has(z)
    F=s.zeros(14)
    for r in range(14):F[13-r,r]=(-1)**r
    P=V*adj(V)
    time_reversed=F*P.applyfunc(lambda x:s.conjugate(x).subs(s.conjugate(z),1/z))*F.T
    overlap=s.simplify(s.trace(P*time_reversed))
    assert s.simplify(overlap-(2-s.Rational(65,162)*(2+z+1/z)))==0
    third_moment=s.simplify((u.T*Jz**3*u)[0])
    assert third_moment==65
    result={'passed':True,'n':13,'order':2,'dimension':2,
      'family':'u=sqrt(5/18) D0+sqrt(13/18) D9; v=sqrt(13/18) D4+z sqrt(5/18) D13; |z|=1',
      'fullCrossReductionEntriesCheckedSymbolically':36,'parameterAssumption':'conjugate(z)=1/z',
      'spinImageDeterminant':str(determinant),'determinantConstant':str(constant),
      'fullTensorImageRank':14,'nonextension':'Every member is inclusion-maximal: no third vector can be added preserving second-order AC.',
      'timeReversalOverlap':str(overlap),'inequivalentParameterInterval':'phi in [0,pi], z=exp(i phi); overlap increases strictly from 32/81 to 2.',
      'thirdMomentOfU':str(third_moment),'exactOrder':2,
      'globalMaximum':'k_max(13,2)=2, proved separately by the auxiliary-spin inertia theorem. Every phase-family member is optimal.'}
    (ROOT/'family-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':print(json.dumps(verify(),indent=2))
