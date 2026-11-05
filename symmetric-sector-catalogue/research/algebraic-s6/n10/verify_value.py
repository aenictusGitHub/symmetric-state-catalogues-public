"""Exact value and Galois certificate for the N=10, S6 stationary candidate.
Requires SymPy. Run beside algebraic.py, polynomials.py, energy-polynomial.json.
This establishes exact attainment and absence of a rational-radical formula,
not global minimality. Output: value-verification.json.
"""
from pathlib import Path
import json
import sympy as s
from algebraic import a,b,c,d,z,R,Z,U,P,Q,F

BASE=Path(__file__).resolve().parent
# This identity proves the shorter value formula from Euler stationarity.
short_identity=(P-Q*s.diff(P,a)/4).subs(a,1)-10080*b*F[0]-10080*c*F[1]-2520*d*F[2]
assert s.expand(short_identity)==0
u,v=s.symbols('u v')
rpoly=s.Poly(R,c,domain=s.QQ)
def red(expr):
    return s.Poly(expr,c,domain=s.QQ).rem(rpoly).as_expr()
def rat(expr):
    num,den=s.fraction(s.cancel(expr))
    assert s.gcd(s.Poly(den,c,domain=s.QQ),rpoly).degree()==0
    return red(red(num)*s.invert(s.Poly(den,c,domain=s.QQ),rpoly).as_expr())
zz=rat(Z)
uu=rat(U.subs(z,zz))
qq=red(1+120*uu+210*c*c+10*uu*red(zz*zz))
value=rat(red(210*(1+168*uu+66*c*c-2*uu*red(zz*zz)-480*uu*c+96*uu*c*zz))/qq)
vp=s.Poly(s.sympify(json.loads((BASE/'energy-polynomial.json').read_text())['polynomial']),v).clear_denoms()[1].primitive()[1]
assert vp.degree()==14
res=0
for coef in vp.all_coeffs():
    res=red(res*value+coef)
assert res==0

# Compare with the independently derived homogeneous sector polynomial P/Q^2.
quartic=P.subs(a,1).subs(d,b*z).expand().subs(b**4,u**2).subs(b**2,u).subs({u:uu,z:zz})
assert red(quartic-red(value*red(qq*qq)))==0
for f in (F[0]/b,F[1],F[2]/b):
    terms=s.Poly(s.expand(f.subs(d,z*b)),b).as_dict()
    h=terms.get((2,),0)*u+terms.get((0,),0)
    assert red(h.subs({u:uu,z:zz}))==0
print('PASS: value polynomial, direct quartic and three stationary identities.',flush=True)

# Rigorous branch selection and rational interval mapping.
lo=s.Rational('0.34031561845452539230789192544421338373007474986524')
hi=s.Rational('0.34031561845452539230789192544421338373007474986525')
assert rpoly.count_roots(lo,hi)==1
lo,hi=rpoly.refine_root(lo,hi,eps=s.Rational(1,10**150))
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):
    products=[a*b for a in x for b in y]
    return (min(products),max(products))
def powi(x,p):
    result=(s.Integer(1),s.Integer(1))
    for _ in range(p):result=mul(result,x)
    return result
def div(x,y):
    assert y[0]*y[1]>0,'Interval denominator contains zero'
    return mul(x,(1/y[1],1/y[0]))
def ipoly(expr,gens,ranges):
    result=(s.Integer(0),s.Integer(0))
    for powers,coef in s.Poly(expr,*gens).terms():
        term=(coef,coef)
        for power,ran in zip(powers,ranges):term=mul(term,powi(ran,power))
        result=add(result,term)
    return result
def rounded_interval(interval,digits=80):
    scale=10**digits
    return [str(s.floor(interval[0]*scale)/scale),str(s.ceiling(interval[1]*scale)/scale)]
zn,zd=s.fraction(Z)
zrange=div(ipoly(zn,[c],[(lo,hi)]),ipoly(zd,[c],[(lo,hi)]))
un,ud=s.fraction(U)
urange=div(ipoly(un,[c,z],[(lo,hi),zrange]),ipoly(ud,[c,z],[(lo,hi),zrange]))
assert zrange[0]>0 and urange[0]>0
qexpr=1+120*u+210*c*c+10*u*z*z
qr=ipoly(qexpr,[c,u,z],[(lo,hi),urange,zrange])
assert qr[0]>0
vexpr=210*(1+168*u+66*c*c-2*u*z*z-480*u*c+96*u*c*z)
vrange=div(ipoly(vexpr,[c,u,z],[(lo,hi),urange,zrange]),qr)
value_interval=[s.Rational('73.04625181143282862655926527676680932516076615331620794697'),s.Rational('73.04625181143282862655926527676680932516076615331620794698')]
assert value_interval[0]<vrange[0]<=vrange[1]<value_interval[1]
assert vp.count_roots(*value_interval)==1
print('PASS: positive real branch maps strictly inside the unique-root interval.',flush=True)

# Squarefree reductions certify irreducibility and Galois-cycle types.
mods={}
expected={59:[14],83:[1,2,11],151:[1,13]}
for prime in expected:
    mod=s.Poly(vp.as_expr(),v,modulus=prime)
    assert mod.degree()==14 and mod.gcd(mod.diff()).degree()==0
    factors=s.factor_list(mod)[1]
    degrees=sorted([f.degree() for f,e in factors for _ in range(e)])
    assert degrees==expected[prime]
    assert all(e==1 and f.is_irreducible for f,e in factors)
    product=s.Poly(1,v,modulus=prime)
    for f,e in factors:product*=f**e
    assert product.monic()==mod.monic()
    mods[str(prime)]={'factorDegrees':degrees,'monicFactors':[str(f.as_expr()) for f,e in factors], 'squarefree':True}
print('PASS: mod 59 irreducible; mod 83 has (1,2,11); mod 151 has (1,13).',flush=True)

proof=(
'Irreducibility modulo 59 makes the rational polynomial irreducible and its Galois group transitive. '
'The squarefree modulo-151 factorization gives a 13-cycle. A nontrivial block system in degree 14 '
'would have 2 or 7 blocks, each of size 7 or 2. The 13-cycle fixes each block setwise, which cannot '
'accommodate its orbit of length 13. Hence the group is primitive. The modulo-83 cycle type '
'(1,2,11), raised to the 11th power, is a transposition. A primitive permutation group containing '
'a transposition is the full symmetric group: the graph of all conjugate transpositions is a '
'nonempty invariant graph, its connected components are blocks, and connected edge transpositions '
'generate the full symmetric group. Thus the Galois group is S14, which is not solvable. For an '
'irreducible polynomial in characteristic zero, a root expressible by radicals over Q would lie in '
'a solvable normal extension after adjoining the needed roots of unity; its conjugates and splitting '
'field would then also lie there. This contradicts the S14 Galois group. No root, including the '
'specified exact sector value, is expressible by radicals over Q.'
)
out={
 'n':10,'sector':6,'globalOptimalityProven':False,
 'valuePolynomialDegree':14,'valuePolynomialVerifiedInNumberField':True,
 'directQuarticIdentityVerified':True,'stationarityIdentitiesVerified':True,
 'shortValueStationarityIdentityVerified':True,
 'positiveRealRecoveryBranchVerified':True,
 'rootParameterInterval':[str(lo),str(hi)],
 'zInterval':rounded_interval(zrange),'uInterval':rounded_interval(urange),
 'normalizationInterval':rounded_interval(qr),
 'mappedValueInterval':rounded_interval(vrange),
 'isolatingValueInterval':[str(x) for x in value_interval],
 'valueRootCount':1,'modularFactorizations':mods,
 'irreducibleOverRationals':True,'galoisGroup':'S14',
 'expressibleByRadicalsOverRationals':False,'galoisProof':proof
}
(BASE/'value-verification.json').write_text(json.dumps(out,indent=2)+'\n')
print('Wrote value-verification.json',flush=True)
