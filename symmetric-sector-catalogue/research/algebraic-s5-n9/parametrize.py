from polynomials import *
from pathlib import Path
import json
u,z=s.symbols('u z');H=[]
for f in [F[0]/b,F[1],F[2]/b]:
 q=s.Poly(s.expand(f.subs(d,z*b)),b).as_dict();H.append(s.expand(q.get((2,),0)*u+q.get((0,),0)))
A=s.diff(H[0],u);B=H[0].subs(u,0)
r12=s.resultant(H[0],H[1],u);r13=s.resultant(H[0],H[2],u)
R=20500*c**4+39640*c**3-10178*c**2+252*c+5
G=s.groebner([r12,r13,R],z,c,order='lex')
print('GB=',G,flush=True)
linear=[p.as_expr() for p in G.polys if p.degree(z)==1][0]
Z=s.factor(-linear.subs(z,0)/s.diff(linear,z));U=s.factor(-B/A)
print('Z=',Z,'U=',U,flush=True)
# Reduce u and z as field elements for concise expressions.
Us=s.cancel(U.subs(z,Z));num,den=s.fraction(Us)
Uf=s.factor(s.rem(s.rem(num,R,c)*s.invert(s.rem(den,R,c),R,c),R,c))
print('u(c)=',Uf,flush=True)
Path(__file__).with_name('algebraic.py').write_text('from polynomials import *\nz=s.symbols("z")\nR='+repr(R)+'\nZ='+repr(Z)+'\nU='+repr(U)+'\nUF='+repr(Uf)+'\n')
Path(__file__).with_name('algebraic.json').write_text(json.dumps({'rootPolynomial':str(R),'zRecovery':str(Z),'uRecovery':str(U),'uPolynomialRecovery':str(Uf)},indent=2))
