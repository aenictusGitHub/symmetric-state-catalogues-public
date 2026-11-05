import sys,math,json
import sympy as s
from optimize import ROOT

def exact_sector(n,k,coeff):
 w=[s.S(0)]*(n//2+1)
 for r in range(k+1):w[min(r,n-r)]+=s.binomial(n,k)*(-1)**(k-r)*s.binomial(k,r)*2**r
 value=0
 for r,wr in enumerate(w):
  if not wr:continue
  m=s.Matrix(r+1,n-r+1,lambda a,b:coeff.get(a+b,0)*s.sqrt(s.binomial(r,a)*s.binomial(n-r,b)/s.binomial(n,a+b)))
  rho=m*m.conjugate().T
  purity=sum(s.expand_complex(abs(x)**2) for x in rho)
  value+=wr*purity
 return s.simplify(value)

def two(n,k,a,b):
 va=exact_sector(n,k,{a:s.S(1)})
 vb=exact_sector(n,k,{b:s.S(1)})
 mid=exact_sector(n,k,{a:s.sqrt(s.Rational(1,2)),b:s.sqrt(s.Rational(1,2))})
 p=s.Symbol('p',real=True)
 # p is probability on first support a.
 poly=s.expand(va*p+vb*(1-p)+(4*mid-2*va-2*vb)*p*(1-p))
 optimum=s.solve(s.diff(poly,p),p)[0]
 value=s.factor(poly.subs(p,optimum))
 print((n,k,a,b),'p=',optimum,'1-p=',1-optimum,'value=',value,'polynomial=',poly,flush=True)
 return {'n':n,'k':k,'support':[a,b],'weights':[str(optimum),str(1-optimum)],'value':str(value),'polynomial':str(poly)}

if __name__=='__main__':
 results=[two(*x) for x in [(9,6,0,8),(10,8,0,9),(11,7,1,8),(11,9,0,8),(11,10,0,9),(11,8,0,10),(12,11,0,10)]]
 (ROOT/'exact-two-support.json').write_text(json.dumps(results,indent=2)+'\n')
 for n,k,co in [
  (10,8,{1:s.sqrt(5)/4,5:s.sqrt(6)/4,9:s.sqrt(5)/4}),
  (10,9,{2:s.Rational(1,2),5:1/s.sqrt(2),8:s.Rational(1,2)}),
  (12,8,{2:s.sqrt(7)/4,6:1/s.sqrt(8),10:s.sqrt(7)/4}),
  (12,9,{0:s.Rational(1,2),6:1/s.sqrt(2),12:-s.Rational(1,2)}),
  (12,10,{1:s.sqrt(7)/5,6:s.sqrt(11)/5,11:-s.sqrt(7)/5}),
 ]:print('three support',n,k,exact_sector(n,k,co),flush=True)
