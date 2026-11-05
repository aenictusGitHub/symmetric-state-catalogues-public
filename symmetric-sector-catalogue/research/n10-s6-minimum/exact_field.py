"""Exact arithmetic in Q[t]/(P14); no numerical approximation is used here."""
import json
from pathlib import Path
from flint import fmpq_poly,fmpq
DATA=json.loads((Path(__file__).parent/'state-parameters.json').read_text())
MOD=fmpq_poly([int(c)for c in DATA['rootPolynomialCoefficients'][::-1]])
class K:
 __slots__=('p',)
 def __init__(self,p=0):
  if isinstance(p,K):self.p=p.p
  elif isinstance(p,fmpq_poly):self.p=p%MOD
  elif isinstance(p,(list,tuple)):self.p=fmpq_poly([fmpq(str(x))for x in p])%MOD
  else:self.p=fmpq_poly([fmpq(str(p))])
 def __add__(self,b):b=K(b);return K(self.p+b.p)
 __radd__=__add__
 def __neg__(self):return K(-self.p)
 def __sub__(self,b):return self+-K(b)
 def __rsub__(self,b):return K(b)+-self
 def __mul__(self,b):b=K(b);return K(self.p*b.p)
 __rmul__=__mul__
 def inv(self):
  g,s,t=self.p.xgcd(MOD);assert g.is_one();return K(s)
 def __truediv__(self,b):return self*K(b).inv()
 def __rtruediv__(self,b):return K(b)*self.inv()
 def __pow__(self,n):
  if n<0:return self.inv()**(-n)
  a=K(1);x=self
  while n:
   if n%2:a=a*x
   x=x*x;n//=2
  return a
 def __bool__(self):return not self.p.is_zero()
 def __eq__(self,b):return self.p==K(b).p
 def serial(self):return [str(c)for c in self.p.coeffs()]
T=K([0,1]);Z=K(DATA['recoveryNumeratorCoefficients'][::-1])/K(DATA['recoveryDenominatorCoefficients'][::-1]);U=(30*T*T*Z+8*T*T-2*T*Z+20*T-1)/(8*(-30*T*Z+150*T+2*Z*Z-29));Q=1+120*U+210*T*T+10*U*Z*Z
BOUND=210*(1+168*U+66*T*T-2*U*Z*Z-480*U*T+96*U*T*Z)/Q
