#!/usr/bin/env python3
"""Derive and verify an exact algebraic stationary candidate for min S6, N=12.

Requires SymPy and mpmath. Run from any directory. This reproduces the
quartic from exact reduced purities, eliminates variables, and verifies
stationarity and the exact value. It does NOT prove a global minimum.
"""
from pathlib import Path
from math import comb
import json
import time
import sympy as s
import mpmath as mp

OUT = Path(__file__).resolve().parent

def primitive(expr, variable):
    poly = s.Poly(expr, variable).clear_denoms()[1].primitive()[1]
    return poly.as_expr() if poly.LC() > 0 else -poly.as_expr()

def interval_add(a,b): return (a[0]+b[0],a[1]+b[1])
def interval_mul(a,b):
    products=[x*y for x in a for y in b]
    return min(products),max(products)
def interval_div(a,b):
    assert not b[0] <= 0 <= b[1], 'Interval denominator contains zero'
    return interval_mul(a,(1/b[1],1/b[0]))
def interval_polynomial(expr, variable, interval):
    out=(s.S.Zero,s.S.Zero)
    for coefficient in s.Poly(expr,variable).all_coeffs():
        out=interval_add(interval_mul(out,interval),(coefficient,coefficient))
    return out
def interval_bivariate(expr,x,y,xi,yi):
    out=(s.S.Zero,s.S.Zero)
    for coefficient in s.Poly(expr,y).all_coeffs():
        out=interval_add(interval_mul(out,yi),interval_polynomial(coefficient,x,xi))
    return out

def main():
    started=time.perf_counter()
    a,b,c,d=s.symbols('a b c d', real=True)
    amplitudes=[a,0,b,0,-c,0,d,0,c,0,b,0,-a]
    norm=2*a*a+2*b*b+2*c*c+d*d
    Q=s.S.Zero
    full_gradient=s.zeros(13,1)
    for r in range(7):
        weight=comb(12,6)*(-1)**(6-r)*comb(6,r)*2**r
        factors=s.Matrix(r+1,13-r,lambda i,j:s.sqrt(s.Rational(
            comb(r,i)*comb(12-r,j),comb(12,i+j))))
        matrix=s.Matrix(r+1,13-r,lambda i,j:amplitudes[i+j]*factors[i,j])
        reduced=matrix*matrix.T
        Q+=weight*s.trace(reduced*reduced)
        derivative=4*reduced*matrix
        for i in range(r+1):
            for j in range(13-r):
                full_gradient[i+j]+=weight*derivative[i,j]*factors[i,j]
    Q=s.expand(Q)
    # Exact full-real-space gradient agrees with the restricted gradient.
    expected=[s.diff(Q,a)/2,0,s.diff(Q,b)/2,0,-s.diff(Q,c)/2,0,
              s.diff(Q,d),0,s.diff(Q,c)/2,0,s.diff(Q,b)/2,0,-s.diff(Q,a)/2]
    assert all(s.expand(g-e)==0 for g,e in zip(full_gradient,expected))
    # At a real state, the imaginary gradient is zero by conjugation symmetry.

    x,y,T=s.symbols('x y T', real=True)
    # x=a/(sqrt(55)c), y=b/(sqrt(14)d), T=d^2/c^2.
    polynomial=s.Poly(s.expand(Q.subs({a:s.sqrt(55)*x,b:s.sqrt(14)*y*d,c:1})),d)
    F=polynomial.coeff_monomial(1)
    G=polynomial.coeff_monomial(d**2)
    H=polynomial.coeff_monomial(d**4)
    assert polynomial.as_expr()==F+G*d*d+H*d**4
    P=110*x*x+2
    R=28*y*y+1
    value=s.cancel((4*F*H-G*G)/(4*(F*R*R-G*P*R+H*P*P)))
    # d/dT [(F+TG+T^2H)/(P+TR)^2]=0 is linear in T.
    t_numerator=s.expand(2*R*F-P*G)
    t_denominator=s.expand(2*P*H-R*G)
    original=(F+T*G+T*T*H)/(P+T*R)**2
    assert s.cancel(original.subs(T,t_numerator/t_denominator)-value)==0

    x_factors=s.factor_list(s.fraction(s.cancel(s.diff(value,x)))[0])[1]
    y_factors=s.factor_list(s.fraction(s.cancel(s.diff(value,y)))[0])[1]
    f=max((z for z,m in x_factors),key=lambda z:s.degree(z,x)+s.degree(z,y))
    g=max((z for z,m in y_factors),key=lambda z:s.degree(z,x)+s.degree(z,y))
    resultant=s.resultant(f,g,y)
    factors=s.factor_list(resultant)[1]
    root_interval=(s.Rational('0.09396888997957705'),s.Rational('0.09396888997957706'))
    choices=[primitive(z,x) for z,m in factors if s.Poly(z,x).count_roots(*root_interval)]
    assert len(choices)==1
    root_polynomial=choices[0]
    assert s.degree(root_polynomial,x)==12 and s.Poly(root_polynomial,x).is_irreducible
    assert s.Poly(root_polynomial,x).count_roots(*root_interval)==1

    # Recover y rationally in Q[x]/(root_polynomial), using a compact 5/6 form.
    basis=s.groebner([f,g,root_polynomial],y,x,order='lex')
    linear=next(z.as_expr() for z in basis.polys if s.degree(z.as_expr(),y)==1)
    y_polynomial=-linear.subs(y,0)/linear.coeff(y)
    remainders=[s.rem(x**j*y_polynomial,root_polynomial,x) for j in range(7)]
    matrix=s.Matrix([[z.coeff(x,k) for z in remainders] for k in range(6,12)])
    kernel=matrix.nullspace()
    assert len(kernel)==1
    B=sum(kernel[0][j]*x**j for j in range(7))
    A=s.rem(B*y_polynomial,root_polynomial,x)
    coefficients=s.Poly(A,x).all_coeffs()+s.Poly(B,x).all_coeffs()
    common=s.ilcm(*[z.q for z in coefficients])
    A,B=s.expand(A*common),s.expand(B*common)
    common=s.igcd(*[z for z in s.Poly(A,x).all_coeffs()+s.Poly(B,x).all_coeffs()])
    A,B=A/common,B/common
    if s.Poly(B,x).LC()<0: A,B=-A,-B
    for eq in (f,g):
        degree=s.degree(eq,y)
        numerator=sum(s.expand(eq).coeff(y,j)*A**j*B**(degree-j) for j in range(degree+1))
        assert s.rem(numerator,root_polynomial,x)==0
    assert s.gcd(B,root_polynomial)==1

    yi=interval_div(interval_polynomial(A,x,root_interval),interval_polynomial(B,x,root_interval))
    ti=interval_div(interval_bivariate(t_numerator,x,y,root_interval,yi),
                    interval_bivariate(t_denominator,x,y,root_interval,yi))
    assert root_interval[0]>0 and yi[0]>0 and ti[0]>0
    norm_interval=interval_add(interval_polynomial(P,x,root_interval),
        interval_mul(interval_polynomial(R,y,yi),ti))
    assert norm_interval[0]>0

    # Compute the exact minimal polynomial of the physical sector value.
    def red(z): return s.rem(z,root_polynomial,x)
    def mul(z,w): return red(z*w)
    def inverse(z): return s.invert(z,root_polynomial,x)
    y_field=mul(A,inverse(B))
    def evaluate_field(expr):
        return red(sum(coefficient*mul(red(x**i),red(y_field**j))
            for (i,j),coefficient in s.Poly(expr,x,y).terms()))
    ff,gg,hh,pp,rr=map(evaluate_field,(F,G,H,P,R))
    numerator=red(4*mul(ff,hh)-mul(gg,gg))
    denominator=red(4*(mul(ff,mul(rr,rr))-mul(gg,mul(pp,rr))+mul(hh,mul(pp,pp))))
    value_field=mul(numerator,inverse(denominator))
    multiplication=s.Matrix([[red(value_field*x**j).coeff(x,i) for j in range(12)] for i in range(12)])
    V=s.symbols('V')
    value_polynomial=primitive(multiplication.charpoly(V).as_expr(),V)
    assert s.Poly(value_polynomial,V).is_irreducible
    value_interval=(s.Rational('142.68788213589408'),s.Rational('142.68788213589410'))
    assert s.Poly(value_polynomial,V).count_roots(*value_interval)==1

    mp.mp.dps=130
    pfun=s.lambdify(x,root_polynomial,'mpmath')
    root=mp.findroot(pfun,(mp.mpf(str(root_interval[0].evalf(40))),mp.mpf(str(root_interval[1].evalf(40)))),tol=mp.mpf('1e-120'))
    yy=s.lambdify(x,A/B,'mpmath')(root)
    tt=s.lambdify((x,y),t_numerator/t_denominator,'mpmath')(root,yy)
    cc=1/mp.sqrt(110*root*root+2+(28*yy*yy+1)*tt)
    aa,bb,dd=mp.sqrt(55)*root*cc,mp.sqrt(14)*yy*mp.sqrt(tt)*cc,mp.sqrt(tt)*cc
    sector=s.lambdify((a,b,c,d),Q,'mpmath')(aa,bb,cc,dd)
    stationarity=[s.diff(Q,z)-2*V*s.diff(norm,z) for z in (a,b,c,d)]
    residual=max(map(abs,s.lambdify((a,b,c,d,V),stationarity,'mpmath')(aa,bb,cc,dd,sector)))
    assert residual<mp.mpf('1e-115')
    coeffs=[aa,0,bb,0,-cc,0,dd,0,cc,0,bb,0,-aa]
    symbolic={'n':12,'k':6,'status':'exact algebraic stationary candidate; global optimality unproved',
      'stateLatex':r'a(|D_{12}^{(0)}\rangle-|D_{12}^{(12)}\rangle)+b(|D_{12}^{(2)}\rangle+|D_{12}^{(10)}\rangle)+c(-|D_{12}^{(4)}\rangle+|D_{12}^{(8)}\rangle)+d|D_{12}^{(6)}\rangle',
      'normalization':'2a^2+2b^2+2c^2+d^2=1','quartic':str(Q),
      'rootVariable':'x=a/(sqrt(55)*c)','rootPolynomial':str(root_polynomial),
      'rootPolynomialCoefficientsDescending':[str(z) for z in s.Poly(root_polynomial,x).all_coeffs()],
      'rootInterval':[str(z) for z in root_interval],
      'yDefinition':'y=b/(sqrt(14)*d)=A(x)/B(x)',
      'yNumerator':str(A),'yDenominator':str(B),
      'TDefinition':'T=d^2/c^2=(2*R*F-P*G)/(2*P*H-R*G)',
      'F':str(F),'G':str(G),'H':str(H),'P':str(P),'R':str(R),
      'coefficientRecovery':{'c':'1/sqrt(P+T*R)','a':'sqrt(55)*x*c','b':'sqrt(14)*y*sqrt(T)*c','d':'sqrt(T)*c'},
      'reducedStationaryPolynomials':[str(f),str(g)],
      'sectorValueFormula':'(4*F*H-G^2)/(4*(F*R^2-G*P*R+H*P^2))',
      'sectorValueMinimalPolynomial':str(value_polynomial),
      'sectorValueRootInterval':[str(z) for z in value_interval]}
    verification={'allPassed':True,'n':12,'k':6,'precisionDigits':130,
      'root100digits':mp.nstr(root,105),'y100digits':mp.nstr(yy,105),'T100digits':mp.nstr(tt,105),
      'parameters100digits':{z:mp.nstr(v,105) for z,v in zip(('a','b','c','d'),(aa,bb,cc,dd))},
      'coefficients100digits':[mp.nstr(v,105) for v in coeffs],
      'sectorValue100digits':mp.nstr(sector,105),'stationarityResidual':mp.nstr(residual,10),
      'exactTests':['Quartic independently derived from reduced purities',
        'Full real gradient has precisely the four-parameter ansatz pattern; imaginary gradient vanishes by conjugation symmetry',
        'Eliminating T reproduces the stated reduced rational objective',
        'Degree-12 root polynomial is an exact irreducible resultant factor',
        'Sturm count: exactly one root in the specified rational interval',
        'Both stationary equations vanish exactly after rational recovery modulo the root polynomial',
        'Rational interval arithmetic certifies positive x,y,T,normalization and nonzero recovery denominators',
        'Degree-12 value polynomial is irreducible and has exactly one root in its specified interval'],
      'positiveIntervals':{'y':[str(z) for z in yi],'T':[str(z) for z in ti],'norm':[str(z) for z in norm_interval]},
      'globalOptimalityProven':False,'elapsedSeconds':time.perf_counter()-started,'sympy':s.__version__}
    (OUT/'exact-algebraic.json').write_text(json.dumps(symbolic,indent=2)+'\n')
    (OUT/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print('PASS: exact algebraic construction, normalization and full-space stationarity.')
    print('S6 =',mp.nstr(sector,105))
    print('Global optimality remains unproved.')

if __name__=='__main__': main()
