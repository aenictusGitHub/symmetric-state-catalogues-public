#!/usr/bin/env python3
"""Independent symbolic audit of the all-N four-band lower-bound theorem.
No import from the construction implementation. All polynomial checks use
rational coefficients and quantify over nonnegative formal variables.
"""
from pathlib import Path
import json
import sympy as S

base=Path(__file__).resolve().parent
x,v=S.symbols('x v',real=True)
w=[v/4+x*x/4-S.Rational(35,36)*x+S.Rational(5,6),
   -v/4-x*x/4+S.Rational(5,12)*x+S.Rational(1,2),
   -v/4-x*x/4+S.Rational(13,12)*x-S.Rational(1,2),
   v/4+x*x/4-S.Rational(19,36)*x+S.Rational(1,6)]
moments=[S.expand(sum(w)-1),S.expand(sum((x-l)*w[l] for l in range(4))),
         S.expand(sum((x-l)**2*w[l] for l in range(4))-v)]
assert moments==[0,0,0]
a,b,c,d=S.symbols('a b c d',positive=True)
L=[a,b,c,d]
F=sum(u*u for u in L)**2-2*sum(u**4 for u in L)
identity=S.expand(F+8*a*b*c*d-S.prod(sum(L)-2*u for u in L))
assert identity==0
q,i,z=S.symbols('q i z',nonnegative=True)
A=S.Poly(S.Integer(7421703487488)*(i**14*z**2+i**2*z**14),i,z)
residues=[]
for rem in range(8):
    n=8*q+rem; spacing=2*q;offset=4+rem//2
    wx=[S.cancel(u.subs({x:(n/2-offset-2*i)/spacing,
                        v:n*(n+2)/(12*spacing**2)})) for u in w]
    scaled=[S.cancel(576*q*q*u) for u in wx]
    assert all(S.denom(u)==1 for u in scaled)
    rowpolys=[S.Poly(u.subs(q,i+z+4).expand(),i,z) for u in scaled]
    assert all(all(c>=0 for c in p.coeffs()) and p.coeff_monomial(1)>0 for p in rowpolys)
    edges=[]
    for ell in range(4):
        rr=offset+2*i+spacing*ell
        # 576^2 q^4 times the square of the Jplus^2 edge magnitude.
        edge=scaled[ell]*scaled[ell].subs(i,i+1)*(rr+1)*(rr+2)*(n-rr)*(n-rr-1)
        edges.append(S.Poly(edge.subs(q,i+z+5).expand(),i,z))
    E=sum(edges)**2-2*sum(u*u for u in edges)
    nonnegative_part=E+A
    product_margin=64*S.prod(edges)-A*A
    assert all(c>=0 for c in nonnegative_part.coeffs())
    assert all(c>=0 for c in product_margin.coeffs())
    residues.append({'residue':rem,'strictWeightPositivity':True,
                     'rowWeightCoefficientCounts':[len(p.terms()) for p in rowpolys],
                     'FplusA_nonnegativeCoefficients':True,
                     '64productMinusA2_nonnegativeCoefficients':True,
                     'FplusA_termCount':len(nonnegative_part.terms()),
                     'productMargin_termCount':len(product_margin.terms())})
report={'claim':'K(N,2) >= floor(N/8)-3 for every integer N>=40',
        'exact':True,'momentIdentities':True,'quadrilateralIdentity':True,
        'rowDomain':'q=i+z+4, i,z>=0',
        'edgeDomain':'q=i+z+5, i,z>=0',
        'supportGap':8,
        'argument':'F>=-A and 64 product(E)>=A^2, A>=0 imply F+8 sqrt(product(E))>=0; equivalently the four positive edge lengths close. Forward bands and independent adjacent-edge phase increments give vanishing directed Jplus^2 compression. Other required compressions follow from parity, disjoint support, and the three moment identities.',
        'residues':residues}
(base/'fourband-independent-audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'exact':True,'allResiduesVerified':len(residues),'momentIdentities':True,'quadrilateralIdentity':True}))
