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
F0=11180400*X**4+209440*X*X-K('35840/3')*X+K('7504/3');G=-517440*X*X*Y*Y-4400*X*X-71680*X*Y+K('392000/33')*Y*Y-K('272384/33')*Y+K('22448/33');H=K('3578176/11')*Y**4+K('223328/11')*Y*Y+K('8020/11');P=110*X*X+2;R=28*Y*Y+1;T=(2*R*F0-P*G)/(2*P*H-R*G);U=105*T;Q=54450*X*X+132*Y*Y*U+990+K('33/7')*U;BOUND=(4*F0*H-G*G)/(4*(F0*R*R-G*P*R+H*P*P));MAPPING={0:(165*X,0),2:(Y,1),4:(K(-1),0),6:(K('1/14'),1),8:(K(1),0),10:(Y,1),12:(-165*X,0)}
