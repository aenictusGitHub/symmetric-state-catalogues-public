from algebraic import *
from pathlib import Path
import json
v=s.symbols('v')
Qx=1+84*UF+84*c*c+UF*Z*Z
V=s.cancel(126*(1+84*UF-4*c*c-UF*Z*Z-320*UF*c+80*UF*c*Z)/Qx)
num,den=s.fraction(V);num=s.rem(num,R,c);den=s.rem(den,R,c)
Vf=s.factor(s.rem(num*s.invert(den,R,c),R,c))
res=s.factor(s.resultant(R,v*den-num,c))
print('Vf=',Vf,flush=True);print('polynomial=',res,flush=True)
print('radicals=',s.solve(s.Poly(res,v).sqf_part().as_expr(),v),flush=True)
Path(__file__).with_name('energy-polynomial.json').write_text(json.dumps({'sectorLengthRecovery':str(Vf),'resultant':str(res)},indent=2))
