#!/usr/bin/env python3
"""Exact, self-contained verification of the N17,k3,t2 algebraic witness.
Requires Python 3 and SymPy; no numerical root solver is used.
"""
import json,time,math
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parent
def verify(write_report=True):
    start=time.time();cert=json.loads((ROOT/'certificate.json').read_text())
    assert cert['case']=={'n':17,'dimension':3,'t':2}
    N=17;supports=[[2,5,8,11,14,17],[2,5,8,11,14,17],[0,3,6,9,12,15]]
    assert cert['supports']==supports
    x=s.symbols('x0:36',real=True);index=0;vs=[];maps=[]
    for b,support in enumerate(supports):
     v={}
     for r in support:
      scale=round(math.sqrt(math.comb(N,r)))
      for part in ['re','im']:maps.append({'basis':b,'r':r,'part':part,'scale':scale})
      v[r]=(x[index]+s.I*x[index+1])/scale;index+=2
     vs.append(v)
    assert maps==cert['coordinateMap']
    fixed={int(i):s.Rational(value) for i,value in cert['fixedCoordinates'].items()}
    assert len(fixed)==9
    assert fixed=={1:s.S(0),3:s.S(0),9:-s.Rational(191,1000),12:s.S(0),13:s.S(0),15:s.S(0),22:s.Rational(41,500),25:s.S(0),32:-s.Rational(34,125)}
    unknown=cert['unknownIndices'];assert sorted(unknown)==sorted(set(range(36))-set(fixed))
    u=[x[i] for i in unknown];fixed={x[i]:v for i,v in fixed.items()}
    def reduced_polynomial(a,b,i,j):
     out=sum(s.binomial(15,l)*vs[a].get(i+l,0)*s.conjugate(vs[b].get(j+l,0)) for l in range(16))
     if a==b and i==j:out-=s.Rational(1,3*s.binomial(2,i))
     return s.expand(out)
    polys=[];meta=[]
    def append(a,b,i,j,complexeq):
     expr=reduced_polynomial(a,b,i,j)
     for part in ['re','im'] if complexeq else ['re']:
      p=s.re(expr) if part=='re' else s.im(expr)
      polys.append(s.expand(p.subs(fixed)));meta.append({'a':a,'b':b,'row':i,'column':j,'part':part})
    for a in range(3):
     for i in range(3):append(a,a,i,i,False)
    for i in range(3):append(0,1,i,i,True)
    for a in [0,1]:
     for i,j in [(0,1),(1,2),(2,0)]:append(a,2,i,j,True)
    assert meta==cert['equationMetadata']
    assert len(polys)==27
    assert all(str(p)==q for p,q in zip(polys,cert['equations']))
    # Every exact full cross reduction is either identically zero or one selected
    # real polynomial, possibly with a minus sign; support and conjugation suffice.
    allowed={s.S(0)}|set(polys)|{-p for p in polys};complete=0
    for a in range(3):
     for b in range(3):
      for i in range(3):
       for j in range(3):
        p=reduced_polynomial(a,b,i,j)
        for component in [s.re(p),s.im(p)]:
         assert s.expand(component.subs(fixed)) in allowed,(a,b,i,j,component)
         complete+=1
    center=list(map(s.Rational,cert['center']));assert len(center)==27
    radius=s.Rational(cert['radius']);assert radius==s.Rational(1,10**60)
    A=s.Matrix([[s.Rational(q) for q in row] for row in cert['inverseJacobianApproximation']]);assert A.shape==(27,27)
    sub=dict(zip(u,center));F=s.Matrix(polys);J=F.jacobian(u);Fx=F.subs(sub);Jx=J.subs(sub)
    inf=lambda M:max(sum(abs(M[i,j]) for j in range(M.cols)) for i in range(M.rows))
    eta=inf(A*Fx);delta=inf(s.eye(27)-A*Jx);aNorm=inf(A)
    H=max(sum(abs(s.diff(p,ui,uj)) for ui in u for uj in u) for p in polys)
    q=delta+radius*aNorm*H
    assert delta<1 # A and J(x0) are nonsingular.
    assert q<1
    assert eta+q*radius<radius
    for key,value in [('newtonDisplacement',eta),('inverseDefect',delta),('inverseInfinityNorm',aNorm),('hessianBound',H),('contractionBound',q),('selfMappingMargin',radius-eta-q*radius)]:assert value==s.Rational(cert['bounds'][key])
    # The center's numerical amplitudes only preview the algebraic state. Rebuild them.
    full={**fixed,**sub};coeffs=[]
    for b in range(3):
     row=[]
     for r in range(18):
      amplitude=s.sqrt(s.binomial(17,r))*vs[b].get(r,s.S(0)).subs(full)
      row.append([float(s.re(amplitude)),float(s.im(amplitude))])
     coeffs.append(row)
    assert max(abs(v-cert['coefficients'][b][r][part]) for b,row in enumerate(coeffs) for r,z in enumerate(row) for part,v in enumerate(z))<1e-15
    basis=json.loads((ROOT/'basis.json').read_text())
    assert basis['id']=='algebraic-n17-k3-t2' and basis['exactnessKind']=='isolated-algebraic-root'
    assert (basis['n'],basis['t'],basis['dimension'])==(17,2,3)
    assert basis['algebraicDefinition']['variables']==[str(x[i]) for i in unknown]
    assert basis['algebraicDefinition']['unknownIndices']==unknown
    assert basis['algebraicDefinition']['fixedCoordinates']==cert['fixedCoordinates']
    assert basis['algebraicDefinition']['isolatingBox']=={'center':cert['center'],'radius':cert['radius'],'norm':'infinity'}
    assert basis['algebraicDefinition']['coordinateMap']==maps
    assert basis['algebraicDefinition']['scales']==[round(math.sqrt(math.comb(17,r))) for r in range(18)]
    intervals=basis['algebraicDefinition']['parameterIntervals']
    for xi,value in fixed.items():assert intervals[str(xi)]==[str(value),str(value)]
    for xi,value in sub.items():assert list(map(s.Rational,intervals[str(xi)]))==[value-radius,value+radius]
    for count in range(18):
        scale=round(math.sqrt(math.comb(17,count)))
        assert (2*scale-1)**2<4*math.comb(17,count)<(2*scale+1)**2
    for b,record in enumerate(basis['basis']):
     for r in range(18):
      expected=s.simplify((s.sqrt(s.binomial(17,r))*vs[b].get(r,s.S(0))).subs(fixed))
      supplied=s.sympify(record['exactCoefficients'][r],locals={str(v):v for v in x})
      assert s.simplify(supplied-expected)==0,(b,r)
      assert max(abs(record['coefficients'][r][i]-coeffs[b][r][i]) for i in range(2))<1e-15
    report={'passed':True,'case':cert['case'],'arithmetic':'exact rational, with explicit rational bounds','unknowns':27,'quadraticEquations':27,'fixedRationalCoordinates':9,'completeRealImaginaryReducedMatrixComponents':complete,'radius':str(radius),'newtonDisplacementUpperBound':str(s.Rational(1,10**90)),'contractionUpperBound':str(s.Rational(1,10**56)),'exactDisplacementBoundVerified':bool(eta<s.Rational(1,10**90)),'exactContractionBoundVerified':bool(q<s.Rational(1,10**56)),'selfMapping':True,'contraction':True,'uniqueRootInBox':True,'basisOrthonormality':'Follows exactly by tracing every reduced cross identity.','wholeSubspaceOrderAtLeast':2,'basisDimension':3,'globalMaximalityProved':False,'elapsedSeconds':round(time.time()-start,3)}
    assert report['exactDisplacementBoundVerified'] and report['exactContractionBoundVerified']
    if write_report:(ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
