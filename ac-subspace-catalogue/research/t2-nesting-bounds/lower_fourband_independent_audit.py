#!/usr/bin/env python3
"""Independent all-N audit of the strengthened four-band AC2 construction.

The weights are reconstructed by solving the three moment equations from
the independently specified first weight. All subsequent certificates are
recomputed, rather than read from the construction's verification report.
"""
from pathlib import Path
import sympy as s
import json
HERE=Path(__file__).resolve().parent

def audit():
    x,v=s.symbols('x v',real=True);w0=s.symbols('w0',real=True)
    b,c,d=s.symbols('b c d',real=True)
    weights=[w0,b,c,d]
    solve=s.solve([sum(weights)-1,sum((x-l)*weights[l]for l in range(4)),sum((x-l)**2*weights[l]for l in range(4))-v],[b,c,d])
    q,i=s.symbols('q i',integer=True,nonnegative=True);u,z=s.symbols('u z',nonnegative=True)
    records=[]
    for rho in range(8):
        n=8*q+rho;j=n/2;h=2*(q+int(rho>=6));offset=2+rho//2 if rho<6 else (rho-2)//2
        xx=(j-offset-2*i)/h;vv=j*(j+1)/(3*h*h)
        first=v/4+x*x/4-s.Rational(35,36)*x+s.Rational(5,6)
        if rho==5:first-=(x-s.Rational(3,2))/(6*q)
        W=[s.cancel(first.subs({x:xx,v:vv}))]
        W.extend(s.cancel(solve[a].subs(w0,first).subs({x:xx,v:vv}))for a in [b,c,d])
        assert s.cancel(sum(W)-1)==0
        assert s.cancel(sum((j-offset-2*i-l*h)*W[l]for l in range(4)))==0
        assert s.cancel(sum((j-offset-2*i-l*h)**2*W[l]for l in range(4))-j*(j+1)/3)==0
        denominator=144*h*h
        scaled=[s.cancel(denominator*w)for w in W]
        assert all(s.denom(w)==1 for w in scaled)
        rowcharts=[{i:a,q:z+5}for a in range(3)]+[{i:u+3,q:u+z+5}]
        edgecharts=[{i:a,q:z+5}for a in range(2)]+[{i:u+2,q:u+z+5}]
        for chart in rowcharts:
            for w in scaled:
                p=s.Poly(s.expand(w.subs(chart,simultaneous=True)),u,z)
                assert all(a>=0 for a in p.coeffs()) and p.coeff_monomial(1)>0
        # Bounds covering all q>=5, row 0<=i<=q-2 and adjacent edge 0<=i<=q-3.
        first_support=offset;last_support=s.expand(offset+2*(q-2)+3*h)
        assert first_support>=0 and s.simplify(n-last_support)>=0
        bandgap=s.expand(h-2*(q-2));assert bandgap>=4
        assert h%2==0 if h.is_number else s.expand(h/2).is_integer
        # x decreases with row index, so the endpoint inequalities put
        # exactly two positive and two negative magnetic levels in each row.
        assert s.simplify((2-xx.subs(i,0))*h)>0
        assert s.simplify((xx.subs(i,q-2)-1)*h)>0
        edges=[]
        for chart in edgecharts:
            Ep=[]
            for l in range(4):
                r=offset+2*i+l*h
                square=scaled[l]*scaled[l].subs(i,i+1)*(r+1)*(r+2)*(n-r)*(n-r-1)
                poly=s.Poly(s.expand(square.subs(chart,simultaneous=True)),u,z)
                assert all(a>=0 for a in poly.coeffs())
                Ep.append(poly)
            F=sum(Ep)**2-2*sum(p*p for p in Ep)
            correction=s.Poly.from_dict({monomial:-value for monomial,value in F.terms()if value<0},(u,z))
            assert all(a>=0 for a in (F+correction).coeffs())
            certificate=64*s.prod(Ep)-correction**2
            assert all(a>=0 for a in certificate.coeffs()) and certificate.coeff_monomial(1)>0
            edges.append({'correctionMonomials':len(correction.terms()),'strictProductConstant':str(certificate.coeff_monomial(1))})
        records.append({'rho':rho,'momentEquationsReconstructedAndVerified':True,'rowCharts':len(rowcharts),'edgeCharts':len(edgecharts),
                        'minimumBandGap':str(bandgap),'lastSupportDistanceFromN':str(s.expand(n-last_support)),
                        'twoPositiveAndTwoNegativeLevelsPerRow':True,
                        'weightsPositiveForAllRows':True,'strictQuadrilateralInequalityForAllEdges':True,'edgeChecks':edges})
        print('independent rho',rho,'PASS',flush=True)
    a,b,c,d=s.symbols('a b c d',positive=True);ss=a+b+c+d
    assert s.expand(s.prod(ss-2*v for v in [a,b,c,d])-((a*a+b*b+c*c+d*d)**2-2*(a**4+b**4+c**4+d**4)+8*a*b*c*d))==0
    return {'passed':True,'theorem':'K(N,2)>=floor(N/8)-1 for every integer N>=40',
      'arithmetic':'Exact symbolic moment identities and integer polynomial coefficient inequalities, independently reconstructed',
      'domainCoverage':'Row charts cover i=0,1,2 with q>=5, and i>=3 with q-i>=2. Edge charts cover i=0,1 with q>=5, and i>=2 with q-i>=3.',
      'physicalSupportProof':'All four bands have one parity, supports are disjoint, every inter-band gap is at least four, and only neighboring rows couple through J+^2. Positive moment weights fix Gram, Jz and Jz^2. The strict quadrilateral inequality admits explicit radical phases cancelling each J+^2 cross entry.',
      'firstOrderParentProof':'Each row has Jz signature(2,2). If p=Dplus^(1/2)uplus and n=Dminus^(1/2)uminus have equal norms, their Hermitian-orthogonal partners pperp=(-conj(p2),conj(p1)), nperp=(-conj(n2),conj(n1)) produce an independent neutral vector. Its span with u is Jz-neutral. Row planes are disjoint and parity kills Jx,Jy, giving an exact AC1 parent of dimension2k.',
      'quadrilateralProof':'F+A has nonnegative coefficients, and 64*product(E)-A^2 has nonnegative coefficients and positive constant. Thus F+8sqrt(product(E))>0, the product of the four quadrilateral side deficits. With four positive lengths at most one deficit can be nonpositive, so every deficit is positive.',
      'records':records}

if __name__=='__main__':
    report=audit();(HERE/'lower-fourband-independent-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Independent all-N proof audit passed.')
