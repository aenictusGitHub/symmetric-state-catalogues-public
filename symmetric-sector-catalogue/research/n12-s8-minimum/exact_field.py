"""Exact arithmetic in the specified algebraic number field."""
import json
from pathlib import Path
from flint import fmpq_poly,fmpq
DATA=json.loads((Path(__file__).parent/"state-parameters.json").read_text())
MOD=fmpq_poly([fmpq(c)for c in DATA["rootPolynomialCoefficients"][::-1]])
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
X=K([0,1]);Y=K(DATA['recoveryNumerator'][::-1])/K(DATA['recoveryDenominator'][::-1])
B=X;A=Y;U=(308*A*B*B+64*A*B-1435*B*B+112*B+570)/(158*A*A+768*A*B+672*B-1220);Q=12*A*A*U+495*B*B+792*U+66;BOUND=(118800*A**4*U*U+11155320*A*A*U*B*B+9599040*A*A*U*U-451440*A*A*U-58544640*A*B*B*U-12165120*A*B*U+140295375*B**4+184638960*B*B*U+4021380*B*B-10644480*B*U+160142400*U*U+36352800*U+3771900)/(Q*Q);MAPPING={1:(A,1),4:(B,0),7:(K(1),1),10:(K(1),0)}
