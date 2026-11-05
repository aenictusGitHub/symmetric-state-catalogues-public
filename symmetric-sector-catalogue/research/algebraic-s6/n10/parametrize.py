from polynomials import *
from pathlib import Path
import json
u,z=s.symbols('u z')
H=[]
for f in [F[0]/b,F[1],F[2]/b]:
 q=s.Poly(s.expand(f.subs(d,z*b)),b);dd=q.as_dict();H.append(s.expand(dd.get((2,),0)*u+dd.get((0,),0)))
A=s.diff(H[0],u);B=H[0].subs(u,0)
r12=s.Poly(s.resultant(H[0],H[1],u),c,z).primitive()[1].as_expr();r13=s.Poly(s.resultant(H[0],H[2],u),c,z).primitive()[1].as_expr()
subs=s.subresultants(r12,r13,z)
print([(s.degree(f,z),s.degree(f,c)) for f in subs],flush=True)
for f in subs:
 if s.degree(f,z)==1:
  ZZ=s.cancel(-f.subs(z,0)/s.diff(f,z));print('z(c)=',s.factor(ZZ),flush=True)
  Znum,Zden=s.fraction(ZZ)
  coeff=s.Poly(s.factor(s.resultant(r12,r13,z)),c).factor_list()[1]
  R=[pol.as_expr() for pol,pow in coeff if pol.degree()==14][0]
  print('R=',R,flush=True)
  print('degreesZ',s.degree(Znum,c),s.degree(Zden,c),flush=True)
  U=s.factor(-B/A);print('u(c,z)=',U,flush=True)
  data={'rootVariable':'c','rootPolynomial':str(R),'rootInterval':['34031561845452539230789192544421338373007474986524/10^50','34031561845452539230789192544421338373007474986525/10^50'],'zNumerator':str(Znum),'zDenominator':str(Zden),'uFormula':str(U)}
  (Path(__file__).resolve().parent / 'algebraic.json').write_text(json.dumps(data,indent=2))
  (Path(__file__).resolve().parent / 'algebraic.py').write_text('from polynomials import *\nz=s.symbols("z")\nR='+repr(R)+'\nZ='+repr(ZZ)+'\nU='+repr(U)+'\n')
