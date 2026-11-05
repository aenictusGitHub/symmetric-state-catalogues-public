"""Exact arithmetic in the isolated degree-eleven real field."""
import sys,json
from flint import fmpq_poly,fmpq
from pathlib import Path
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

B=K([0,1]);A=K(DATA['Anumerator'][::-1])/K(DATA['Adenominator'][::-1]);V=(10*A*B*B-10*A*B+706*B*B+440*B+39)/(A*A-84*A*B+1848*B+917)
T=1/A;U=B*B/(A*A*V);Z=1/B;Q=1+165*U+462*T*T+55*U*Z*Z
BOUND=(462+152460*U+49896*T*T-4620*U*Z*Z-1774080*U*T+1774080*U*Z*T+7477470*U*U+16659720*U*T*T+9438660*U*U*Z*Z-39029760*U*Z*T*T+77810040*T**4+8842680*U*T*T*Z*Z+732270*U*U*Z**4)/(Q*Q)
