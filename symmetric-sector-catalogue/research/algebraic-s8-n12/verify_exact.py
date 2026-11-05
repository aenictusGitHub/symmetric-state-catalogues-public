from pathlib import Path
from fractions import Fraction as F
from math import comb
import sympy as s,json,hashlib
R=Path(__file__).parent;raw=json.loads((R/'exact-algebraic.json').read_text());A,B,T=s.symbols('A B T',real=True);a,b,c,d=s.symbols('a b c d',real=True);loc=dict(A=A,B=B,T=T,a=a,b=b,c=c,d=d)
p=s.Poly(s.sympify(raw['Bpolynomial'],locals=loc),B);mod=p.as_expr();lo,hi=[s.sympify(z.replace('^','**'))for z in raw['Binterval']];assert p.count_roots(lo,hi)==1 and p.is_irreducible;print('PASS isolated irreducible root',p.degree(),flush=True)
def rem(v):return s.rem(v,mod,B)
def power(v,n):
 r=s.Integer(1)
 for i in range(n):r=rem(r*v)
 return r
def evaluate(poly,vals):
 out=0
 for powers,co in poly.terms():
  z=co
  for v,e in zip(vals,powers):z=rem(z*power(v,e))
  out=rem(out+z)
 return out
an=s.sympify(raw['Anumerator'],locals=loc);ad=s.sympify(raw['Adenominator'],locals=loc);av=rem(an*s.invert(ad,mod,B))
tn=evaluate(s.Poly(s.sympify(raw['Tnumerator'],locals=loc),A,B),[av,B]);td=evaluate(s.Poly(s.sympify(raw['Tdenominator'],locals=loc),A,B),[av,B]);tv=rem(tn*s.invert(td,mod,B))
for text in raw['reducedEquations']:assert evaluate(s.Poly(s.sympify(text,locals=loc),A,B,T),[av,B,tv])==0
print('PASS exact stationary equations',flush=True)
# Derive the full quartic directly from reduced density matrices, not the saved formula.
n=12;k=8;xs={1:a,4:b,7:c,10:d};sector=0
for t in range(k+1):
 purity=0
 for i in range(t+1):
  for j in range(t+1):purity+=comb(t,i)*comb(t,j)*sum(comb(n-t,l)*xs.get(i+l,0)*xs.get(j+l,0)for l in range(n-t+1))**2
 sector+=comb(n,k)*(-1)**(k-t)*comb(k,t)*2**t*purity
sector=s.expand(sector);assert sector==s.sympify(raw['sectorPolynomial'],locals=loc)
Q=12*a*a+495*b*b+792*c*c+66*d*d;Fd=s.diff(sector,d).subs(d,1)
for v,w in [(a,12),(b,495),(c,792)]:
 expr=s.Poly(s.expand((66*s.diff(sector,v).subs(d,1)-w*v*Fd).subs(a,A*c).subs(b,B)),c);low=min(m[0]for m,co in expr.terms());expr=s.Poly(s.cancel(expr.as_expr()/c**low),c);assert all(m[0]%2==0 for m,co in expr.terms());z=sum(co*T**(m[0]//2)for m,co in expr.terms());assert evaluate(s.Poly(z,A,B,T),[av,B,tv])==0
# Exact interval arithmetic fixes all real branches and bounds the attained value.
def add(x,y):return x[0]+y[0],x[1]+y[1]
def mul(x,y):
 v=[a*b for a in x for b in y];return min(v),max(v)
def div(x,y):
 assert y[0]*y[1]>0
 return mul(x,(1/y[1],1/y[0]))
def ivpoly(poly,x):
 r=(F(0),F(0))
 for co in s.Poly(poly,B).all_coeffs():r=add(mul(r,x),(F(co),F(co)))
 return r
bi=(F(lo),F(hi));ai=div(ivpoly(an,bi),ivpoly(ad,bi));ti=div(ivpoly(tn,bi),ivpoly(td,bi));assert ti[0]>0
qv=rem(12*power(av,2)*tv+495*B*B+792*tv+66)
f=s.Poly(s.expand(sector.subs(d,1).subs(a,A*c).subs(b,B)),c);assert all(m[0]%2==0 for m,co in f.terms());f=sum(co*T**(m[0]//2)for m,co in f.terms());fv=evaluate(s.Poly(f,A,B,T),[av,B,tv]);ev=rem(fv*s.invert(power(qv,2),mod,B));qi=ivpoly(qv,bi);assert qi[0]>0;ei=ivpoly(ev,bi)
print('PASS real branches and exact value enclosure',float(ei[0]),float(ei[1]),flush=True)
# The multiplication matrix yields an exact defining polynomial for the attained value.
M=s.Matrix.hstack(*[s.Matrix([s.Poly(rem(ev*B**j),B).nth(i)for i in range(14)])for j in range(14)])
import sys
from flint import fmpq,fmpq_mat
cp=fmpq_mat([[fmpq(str(x))for x in row]for row in M.tolist()]).charpoly();V=s.Symbol('V');vp=s.Poly(sum(s.Rational(str(co))*V**i for i,co in enumerate(cp)),V).clear_denoms()[1].primitive()[1];assert vp.is_irreducible
assert evaluate(s.Poly(vp.as_expr().subs(V,A),A),[ev])==0
(R/'value-polynomial.json').write_text(json.dumps({'polynomial':str(vp.as_expr()),'degree':14,'rootInterval':[str(z)for z in ei]},indent=2)+'\n')
assert vp.count_roots(s.Rational(ei[0]),s.Rational(ei[1]))==1
report={'n':12,'sector':8,'rootPolynomialDegree':14,'rootCount':1,'irreducible':True,'exactReducedPurityDerivation':True,'exactStationarityVerified':True,'fullComplexStationarityExplanation':'Every off-support real derivative vanishes by the C3 residue support. Every imaginary derivative vanishes for this real state by conjugation invariance. The on-support equations hold exactly in the isolated real number field.','positiveSquareRootBranchVerified':True,'positiveNormalizationVerified':True,'denominatorsNonzero':True,'attainedValueDegree':14,'attainedValueInterval':[str(z)for z in ei],'globalOptimalityProven':False,'inputSha256':hashlib.sha256((R/'exact-algebraic.json').read_bytes()).hexdigest()}
(R/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS analytical stationary state and degree-14 attained value; optimality not claimed.',flush=True)
