#!/usr/bin/env python3
"""Independent exact audit of the prescribed-moment four-band family.
Does not import either constructor. Reconstructs three weights by solving
the physical moment equations after specifying the fourth free weight.
"""
from pathlib import Path
import json
import sympy as s
q,i,U,V=s.symbols('q i U V',integer=True,nonnegative=True)

def nonnegative_poly(expr,strict=True):
    p=s.Poly(s.expand(expr),U,V)
    assert all(a>=0 for a in p.coeffs())
    if strict:assert p.coeff_monomial(1)>0
    return p

def case(rho):
    n=8*q+rho
    h=2*q+(2 if rho==8 else 1)
    offset=2 if rho==8 else (rho-1)//2
    m=s.Matrix([s.Rational(1,2)*n-offset-2*i-ell*h for ell in range(4)])
    variance=n*(n+2)/12
    x=m[0]/h
    # Only the free weight is prescribed. Other weights are solved anew.
    w3=(variance/h**2+x**2)/4-s.Rational(19,36)*x+s.Rational(1,6)
    if rho<8:w3+=(x-s.Rational(3,2))/(3*h)
    matrix=s.Matrix([[1,1,1],[m[0],m[1],m[2]],[m[0]**2,m[1]**2,m[2]**2]])
    rhs=s.Matrix([1-w3,-m[3]*w3,variance-m[3]**2*w3])
    w=list(matrix.inv()*rhs)+[w3]
    w=[s.factor(a) for a in w]
    moments=[s.factor(sum(w)-1),s.factor(sum(a*b for a,b in zip(w,m))),s.factor(sum(a*b*b for a,b in zip(w,m))-variance)]
    assert moments==[0,0,0]
    scale=144*h*h
    W=[s.cancel(scale*a) for a in w]
    assert all(s.denom(a)==1 for a in W)
    # Domain coverage: i=0; or U=i-1,V=q-i-1. Both are nonnegative.
    rowcharts=[{i:0,q:V+2},{i:U+1,q:U+V+2}]
    rowdetails=[]
    for chart in rowcharts:
        polys=[nonnegative_poly(a.subs(chart,simultaneous=True)) for a in W]
        rowdetails.append({'substitution':str(chart),'positiveConstants':[str(a.coeff_monomial(1)) for a in polys]})
    r=[offset+2*i+ell*h for ell in range(4)]
    squared=[]
    for a,ri in zip(W,r):
        # This is the physical Jplus^2 matrix element squared.
        squared.append(a*a.subs(i,i+1)*(ri+1)*(ri+2)*(n-ri)*(n-ri-1))
    chart={i:U,q:U+V+2} # U=i,V=q-i-2 covers every adjacent row.
    E=[nonnegative_poly(a.subs(chart,simultaneous=True)) for a in squared]
    # Independently form the four-side area expression from pair products.
    area=2*sum(E[a]*E[b] for a in range(4) for b in range(a+1,4))-sum(a*a for a in E)
    negative=s.Poly.from_dict({exponent:-coefficient for exponent,coefficient in area.terms() if coefficient<0},(U,V))
    positive=nonnegative_poly((area+negative).as_expr(),strict=False)
    lower=nonnegative_poly((64*s.prod(E)-negative*negative).as_expr())
    # Support endpoints and interband gaps, in q>=2.
    endpoint_slack=s.factor(n-r[-1].subs(i,q-1))
    gap=s.factor(h-2*(q-1))
    assert endpoint_slack>=0 and gap>=3
    return {'rho':rho,'momentResiduals':list(map(str,moments)),'rowCharts':rowdetails,
       'edgeChart':str(chart),'edgePositiveConstants':[str(a.coeff_monomial(1)) for a in E],
       'areaPlusNegativePartNonnegative':True,'strictProductCertificate':True,
       'productPolynomialTerms':len(lower.terms()),'productPositiveConstant':str(lower.coeff_monomial(1)),
       'supportEndpointSlack':str(endpoint_slack),'interbandGap':str(gap)}

def main():
    a,b,c,d=s.symbols('a b c d',positive=True)
    sides=[a,b,c,d];squares=[x*x for x in sides]
    area=2*sum(squares[i]*squares[j] for i in range(4) for j in range(i+1,4))-sum(x*x for x in squares)
    assert s.expand(s.prod(sum(sides)-2*x for x in sides)-area-8*a*b*c*d)==0
    D=s.symbols('D',positive=True)
    ca=(D*D+a*a-b*b)/(2*D*a);cb=(D*D+b*b-a*a)/(2*D*b)
    assert s.factor(a*ca+b*cb-D)==0
    assert s.factor(a*a*(1-ca*ca)-b*b*(1-cb*cb))==0
    cases=[]
    for rho in range(1,9):
        cases.append(case(rho));print(f'remainder {rho}: exact independent audit passed',flush=True)
    report={'passed':True,'theorem':'K(N,2)>=floor((N-1)/8) for every integer N>=17',
        'independence':'No construction module imported. Three physical moment weights reconstructed by a 3x3 Vandermonde solve from one prescribed free weight; area identity reconstructed as pair products.',
        'arithmetic':'Exact symbolic rational and integer polynomials, covering all q>=2 and all eight remainders.',
        'domainCoverage':{'rows':'i=0,q=V+2; or i=U+1,q=U+V+2 with U=i-1,V=q-i-1>=0','edges':'i=U,q=U+V+2 with U=i,V=q-i-2>=0'},
        'quadrilateralLogic':'For positive sides, the four factors sum(side)-2*side have at most one nonpositive member. Their product equals F+8sqrt(product E)>0 by the coefficient certificate, so all are strictly positive. The two diagonal intervals overlap strictly; their intersection midpoint gives radical phases.',
        'fullCrossMomentLogic':'Distinct supports give normalization and diagonal-moment orthogonality. The union has no distance-one pair. Distance-two pairs occur only between consecutive rows in each band. The quadrilateral phases cancel every such cross entry. No other rank-one/two entries survive.',
        'records':cases}
    Path(__file__).with_name('search-independent-family-audit.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
