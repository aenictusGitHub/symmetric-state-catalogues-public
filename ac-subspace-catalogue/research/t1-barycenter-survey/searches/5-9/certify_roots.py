from pathlib import Path
from fractions import Fraction as F
from math import isqrt,comb
import sympy as s,json
R=Path(__file__).resolve().parent

def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):
 v=[x*y for x in a for y in b];return min(v),max(v)
def neg(a):return -a[1],-a[0]
def abshi(a):return max(abs(a[0]),abs(a[1]))
def abslo(a):return 0 if a[0]<=0<=a[1]else min(abs(a[0]),abs(a[1]))
def interval(e):
 e=s.sympify(e)
 if e.is_Rational:q=F(int(e.p),int(e.q));return q,q
 if e.is_Add:
  a=(F(0),F(0))
  for x in e.args:a=add(a,interval(x))
  return a
 if e.is_Mul:
  a=(F(1),F(1))
  for x in e.args:a=mul(a,interval(x))
  return a
 if e.is_Pow and e.exp==s.Rational(1,2)and e.base.is_Rational:
  q=F(int(e.base.p),int(e.base.q));scale=10**120;k=isqrt(q.numerator*scale*scale//q.denominator)
  assert F(k*k,scale*scale)<=q<F((k+1)**2,scale*scale)
  return F(k,scale),F(k+1,scale)
 raise ValueError(e)
def cmul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
rows=json.load(open(R/'exact-wronskians.json'));records=[]
for row in rows:
 n=row['n'];state=json.load(open(R/f'n{n}-result.json'));z0=state['roots'][0]['root'];c=tuple(F(x)for x in z0);rad=F(1,10**60)
 coefficients={int(e):(interval(s.re(s.sympify(v))),interval(s.im(s.sympify(v))))for e,v in row['coefficients'].items()};powers=[(F(1),F(0))]
 for _ in range(row['degree']):powers.append(cmul(powers[-1],c))
 def derivative(order):
  out=[(F(0),F(0)),(F(0),F(0))]
  for e,(ar,ai)in coefficients.items():
   if e<order:continue
   fac=1
   for j in range(order):fac*=e-j
   br,bi=powers[e-order];br*=fac;bi*=fac
   out[0]=add(out[0],add(mul(ar,(br,br)),neg(mul(ai,(bi,bi)))))
   out[1]=add(out[1],add(mul(ar,(bi,bi)),mul(ai,(br,br))))
  return out
 val=derivative(0);der=derivative(1);res=sum(abshi(v)for v in val);slope=max(abslo(v)for v in der);cbound=abs(c[0])+abs(c[1]);tail=F(0)
 for e,a in coefficients.items():
  for j in range(2,e+1):tail+=sum(abshi(v)for v in a)*comb(e,j)*cbound**(e-j)*rad**j
 assert res+tail<slope*rad
 report={'n':n,'passed':True,'polynomialCoefficients':row['coefficients'],'squareFree':row['squareFree'],'center':[str(q)for q in c],'radius':str(rad),'exactRootCountInDisk':1,'method':'Exact rational Rouché inequality: |W(c)| + higher Taylor tail < |Wprime(c)| radius.','residualUpper':str(res),'tailUpper':str(tail),'linearBoundaryLower':str(slope*rad),'ratioApproximation':float((res+tail)/(slope*rad)),'exactMultiplicity':state['dimension'],'multiplicityArgument':'W has a simple nonzero root. The derivative-row Wronskian has rank k−1 and the first k−1 rows are independent. Its cofactor vector defines a nonzero polynomial with exactly k vanishing root multiplicity; an extra derivative zero would force Wprime=0.','optimalityProved':False}
 (R/f'n{n}-root-certificate.json').write_text(json.dumps(report,indent=2)+'\n');records.append(report);print(n,'certified',report['ratioApproximation'])
(R/'root-certificates.json').write_text(json.dumps({'passed':True,'records':records},indent=2)+'\n')
