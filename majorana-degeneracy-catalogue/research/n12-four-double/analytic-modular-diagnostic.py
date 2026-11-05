#!/usr/bin/env python3
"""Bounded exploratory modular Groebner diagnostic; not a real nonexistence proof.

Run with SageMath's Python. A caller should impose a time/resource limit.
The monic chart has one doubled root at0, no root at infinity and a real
nonzero constant term for the remaining cubic Q/z after a z-axis rotation.
Other charts and real positivity conditions remain necessary for a global proof.
"""
from sage.all import GF, PolynomialRing, binomial
import json,time
start=time.monotonic()
field=GF(101);imag=field(-1).sqrt()
names='a b c d e f g h j k l m v'.split()
ring=PolynomialRing(field,names=names,order='degrevlex')
a,b,c,d,e,f,g,h,j,k,l,m,v=ring.gens()
poly=PolynomialRing(ring,'z');z=poly.gen()
Q=z*(a+(b+imag*c)*z+(d+imag*e)*z*z+z**3)
Qbar=z*(a+(b-imag*c)*z+(d-imag*e)*z*z+z**3)
R=f+imag*g+(h+imag*j)*z+(k+imag*l)*z*z+(m+imag*v)*z**3+z**4
Rbar=f-imag*g+(h-imag*j)*z+(k-imag*l)*z*z+(m-imag*v)*z**3+z**4
P=Q*Q*R;Pbar=Qbar*Qbar*Rbar
coeff=[P[i]for i in range(13)];bar=[Pbar[i]for i in range(13)]
rho={}
for A in range(4):
 for C in range(4):
  rho[A,C]=sum(field(binomial(9,t))/field(binomial(12,A+t)*binomial(12,C+t))*coeff[A+t]*bar[C+t]for t in range(10))
equations=[ring(field(binomial(3,A))*rho[A,A]-rho[0,0])for A in range(1,4)]
for A in range(4):
 for C in range(A+1,4):
  equations += [ring((rho[A,C]+rho[C,A])/2),ring((rho[A,C]-rho[C,A])/(2*imag))]
print(json.dumps({'stage':'equations_constructed','fieldCharacteristic':101,'variables':len(names),'equations':len(equations),'degrees':[int(p.total_degree())for p in equations],'terms':[len(p.monomials())for p in equations],'seconds':time.monotonic()-start}),flush=True)
ideal=ring.ideal(equations)
basis=ideal.groebner_basis(algorithm='libsingular:slimgb')
print(json.dumps({'stage':'groebner_completed','basisLength':len(basis),'basisDegrees':[int(p.total_degree())for p in basis],'containsOne':any(p==1 for p in basis),'seconds':time.monotonic()-start,'scope':'Exploratory finite-field unsaturated chart; not a real global proof.'}),flush=True)
