"""Exact rational interval branch checks, no floating-point arithmetic."""
from fractions import Fraction as F
from pathlib import Path
import json,sympy as s
R=Path(__file__).resolve().parent
A,B,T=s.symbols('A B T',real=True);loc={'A':A,'B':B,'T':T}
d=json.loads((R/'exact-algebraic.json').read_text())
lo=F('0.2498461150300332175000977146562933123120')
hi=F('0.2498461150300332175000977146562933123121')

def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):
 p=[a*b for a in x for b in y];return min(p),max(p)
def inv(x):
 assert not x[0]<=0<=x[1], 'Interval denominator contains zero'
 return (1/x[1],1/x[0])
def div(x,y):return mul(x,inv(y))
def power(x,n):
 r=(F(1),F(1))
 for _ in range(n):r=mul(r,x)
 return r
def poly(expr,variables,intervals):
 p=s.Poly(expr,*variables);out=(F(0),F(0))
 for powers,coef in p.terms():
  v=(F(int(coef.p),int(coef.q)),)*2
  for exponent,iv in zip(powers,intervals):v=mul(v,power(iv,exponent))
  out=add(out,v)
 return out
def rational_bounds(iv,places=12):
 scale=10**places
 lower=(iv[0]*scale).__floor__();upper=(iv[1]*scale).__ceil__()
 return [str(F(lower,scale)),str(F(upper,scale))]

biv=(lo,hi)
niv=poly(s.sympify(d['AnumeratorReduced'],locals=loc),(B,),(biv,))
divisor=poly(s.sympify(d['AdenominatorReduced'],locals=loc),(B,),(biv,))
aiv=div(niv,divisor)
assert aiv[1]<0
tn=poly(s.sympify(d['Tnumerator'],locals=loc),(A,B),(aiv,biv))
td=poly(s.sympify(d['Tdenominator'],locals=loc),(A,B),(aiv,biv))
tiv=div(tn,td)
assert tiv[0]>0
qiv=poly(A*A*T+165*B*B+462*T+55,(A,B,T),(aiv,biv,tiv))
assert qiv[0]>0
checked={'B':biv,'Anumerator':niv,'Adenominator':divisor,'A':aiv,'Tnumerator':tn,'Tdenominator':td,'T':tiv,'Q':qiv}
report={'method':'Exact Fraction interval arithmetic. Reported bounds are rounded outwards to rational multiples of 10^-12.','intervals':{key:rational_bounds(iv) for key,iv in checked.items()},'Anegative':True,'Tpositive':True,'Qpositive':True,'denominatorsExcludeZero':True}
(R/'branch-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
