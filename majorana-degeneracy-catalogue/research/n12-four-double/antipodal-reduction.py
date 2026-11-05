from pathlib import Path
import json
import sympy as s
z,h,b,c,d,e,f=s.symbols('z h b c d e f', real=True)
Q=z**4+2*h*z**2+1
R=(b+s.I*d)*z**4+(e+s.I*f)*z**3+c*z**2+(-e+s.I*f)*z+b-s.I*d
P=s.Poly(s.expand(Q**2*R),z);p=[P.nth(i) for i in range(13)]
w=[s.expand(v*s.conjugate(v))/s.binomial(12,r) for r,v in enumerate(p)]
zz=sum(((6-r)**2-14)*w[r] for r in range(13))
jplus=sum(s.conjugate(p[r-1])*p[r]*s.Rational(r,s.binomial(12,r-1))*(13-2*r) for r in range(1,13))
jplus2=sum(s.conjugate(p[r-2])*p[r]*s.Rational(r*(r-1),s.binomial(12,r-2)) for r in range(2,13))
eqs=[zz,s.re(jplus),s.im(jplus),s.re(jplus2),s.im(jplus2)]
primitive=[]
for E in eqs:
 E=s.Poly(s.expand(E),h,b,c,d,e,f).clear_denoms()[1].primitive()[1].as_expr();primitive.append(s.factor(E))
for i,E in enumerate(primitive):print(i,E)
Path(__file__).with_suffix('.json').write_text(json.dumps({'Q':str(Q),'R':str(R),'variables':['h','b','c','d','e','f'],'conditions':list(map(str,primitive))},indent=2)+'\n')
