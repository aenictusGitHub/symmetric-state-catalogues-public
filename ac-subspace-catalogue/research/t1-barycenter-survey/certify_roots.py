#!/usr/bin/env python3
"""Certify isolated simple Wronskian roots with exact rational interval arithmetic.
The polynomial is rebuilt from the saved exact Dicke basis, independently of
all numerical searches. The sparse Wronskian identity uses Vandermonde factors.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,comb,prod
from itertools import product
from functools import lru_cache
import json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parent
PREC=120

def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def mul(a,b):
 p=[x*y for x in a for y in b];return min(p),max(p)
def sqrt_interval(a):
 assert a[0]>=0
 scale=10**PREC;lo=isqrt(a[0].numerator*scale*scale//a[0].denominator);hi=isqrt(a[1].numerator*scale*scale//a[1].denominator)+1
 assert F(lo,scale)**2<=a[0] and F(hi,scale)**2>a[1]
 return F(lo,scale),F(hi,scale)
def mag_hi(a):return max(abs(a[0]),abs(a[1]))
def mag_lo(a):return F(0)if a[0]<=0<=a[1]else min(abs(a[0]),abs(a[1]))
ZERO=(F(0),F(0));ONE=(F(1),F(1))
def cadd(a,b):return add(a[0],b[0]),add(a[1],b[1])
def cmul(a,b):return add(mul(a[0],b[0]),neg(mul(a[1],b[1]))),add(mul(a[0],b[1]),mul(a[1],b[0]))
def cscale(a,x):return mul(a[0],(x,x)),mul(a[1],(x,x))
def cexact(x,y=0):return (F(x),F(x)),(F(y),F(y))
@lru_cache(None)
def real_expr(x):
 if x.is_Rational:return (F(int(x.p),int(x.q)),)*2
 if x.func==s.Add:
  v=ZERO
  for a in x.args:v=add(v,real_expr(a))
  return v
 if x.func==s.Mul:
  v=ONE
  for a in x.args:v=mul(v,real_expr(a))
  return v
 if x.func==s.Pow and x.exp==s.Rational(1,2):return sqrt_interval(real_expr(x.base))
 raise ValueError(str(x))
def complex_expr(text):
 x=s.sympify(text,rational=True).expand(complex=True);return real_expr(s.re(x)),real_expr(s.im(x))

def wronskian(source):
 n=source['n'];k=source['dimension'];cols=[]
 for v in source['basis']:
  col=[]
  for r,text in enumerate(v['exactCoefficients']):
   if text=='0':continue
   x=complex_expr(text);w=sqrt_interval((F(comb(n,r)),)*2);x=cmul(x,(w,ZERO));x=cscale(x,F((-1)**r));col.append((r,x))
  cols.append(col)
 coef={}
 for choice in product(*cols):
  powers=[r for r,_ in choice];vand=prod(powers[j]-powers[i]for i in range(k)for j in range(i+1,k))
  if not vand:continue
  power=sum(powers)-k*(k-1)//2;assert power>=0
  val=cexact(vand)
  for _,x in choice:val=cmul(val,x)
  coef[power]=cadd(coef.get(power,cexact(0)),val)
 return coef

def dyadic(x,upper):
 if not x:return F(0),'0'
 e=x.numerator.bit_length()-x.denominator.bit_length()
 q=F(2)**e
 if upper and q<x:e+=1;q*=2
 if not upper and q>x:e-=1;q/=2
 return q,'2^'+str(e)

def certify(source,root):
 coeff=wronskian(source);degree=max(coeff)
 # Rational centers truncate to45 decimal places; disk radius1e-30.
 scale=10**45;c=tuple(F(int(F(str(v))*scale),scale)for v in root);radius=F(1,10**30)
 C=cexact(*c);powers=[cexact(1)]
 for i in range(degree):powers.append(cmul(powers[-1],C))
 value=cexact(0);deriv=cexact(0)
 for power,a in coeff.items():
  value=cadd(value,cmul(a,powers[power]))
  if power:deriv=cadd(deriv,cscale(cmul(a,powers[power-1]),F(power)))
 residual=mag_hi(value[0])+mag_hi(value[1]);slope=max(mag_lo(deriv[0]),mag_lo(deriv[1]));tail=F(0);c_bound=abs(c[0])+abs(c[1])
 for power,a in coeff.items():
  mag=mag_hi(a[0])+mag_hi(a[1])
  for j in range(2,power+1):tail+=mag*comb(power,j)*c_bound**(power-j)*radius**j
 assert residual+tail<slope*radius,(source['n'],float(residual),float(slope*radius))
 ru,rs=dyadic(residual,True);tu,ts=dyadic(tail,True);ll,ls=dyadic(slope*radius,False)
 assert ru+tu<ll
 return dict(n=source['n'],subspaceId=source['id'],passed=True,rootCountWithMultiplicity=1,rootDiskCenter=[str(x)for x in c],rootDiskRadius=str(radius),residualUpper=rs,tailUpper=ts,linearBoundaryLower=ls,roucheRatio=float((residual+tail)/(slope*radius)),wronskianDegree=degree,method='Exact rational complex-interval coefficient evaluation and Rouche comparison with Wprime(c)(z-c).',consequence='A unique simple Wronskian root exists in the disk. Its derivative-kernel state has a root of exactly multiplicity k, where k is the dimension of the chosen subspace.',globalOptimalityProved=False)

if __name__=='__main__':
 data=json.loads((ROOT/'survey.json').read_text());sources=json.loads((ROOT/'subspaces.json').read_text())
 sources={s['n']:s for s in sources};out=[]
 for row in data['results']:
  if not row.get('wronskianRoot'):continue
  result=certify(sources[row['n']],row['wronskianRoot']);out.append(result);print(row['n'],'certified',result['roucheRatio'],flush=True)
 (ROOT/'root-certificates.json').write_text(json.dumps({'passed':True,'precisionDigitsOfRadicalEnclosures':PREC,'surveySHA256':hashlib.sha256((ROOT/'survey.json').read_bytes()).hexdigest(),'certificates':out},indent=2)+'\n')
