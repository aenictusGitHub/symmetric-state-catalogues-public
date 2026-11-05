from algebraic import *
from pathlib import Path
import time,json
u,v=s.symbols('u v');start=time.time()
Ur=s.cancel(U.subs(z,Z))
NQ=1+120*u+210*c*c+10*u*z*z
nu=210*(1+168*u+66*c*c-2*u*z*z-480*u*c+96*u*c*z)/NQ
vr=s.cancel(nu.subs(z,Z).subs(u,Ur));num,den=s.fraction(vr)
num=s.rem(num,R,c);den=s.rem(den,R,c)
print('reduced rational secs',time.time()-start,flush=True)
# The resultant vanishes at the exact sector length, then primitive square-free part.
res=s.Poly(s.resultant(R,v*den-num,c),v).primitive()[1].sqf_part()
print('degree',res.degree(),'secs',time.time()-start,flush=True)
factors=s.factor_list(res)[1];print('factors',[(f.degree(),pow) for f,pow in factors],flush=True)
for f,pow in factors:
 print(f.as_expr(),flush=True)
Path(__file__).resolve().with_name('energy-polynomial.json').write_text(json.dumps({'variable':'v','polynomial':str(res.as_expr()),'coefficients':[str(x) for x in res.all_coeffs()],'degree':res.degree(),'factorDegrees':[f.degree() for f,pow in factors]},indent=2))
