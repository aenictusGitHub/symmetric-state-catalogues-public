"""Direct exact-factorial CG contractions, enclosed by Arb.
All spins and magnetic labels passed to cg_ball are doubled integers.
"""
from functools import lru_cache
from fractions import Fraction
from math import factorial
import flint
A=flint.arb;Q=flint.fmpq
F=[factorial(i)for i in range(512)]
@lru_cache(maxsize=None)
def cg_ball(a,b,c,m,k):
 M=m+k
 if any((x+y)%2 for x,y in [(a,m),(b,k),(c,M)])or abs(m)>a or abs(k)>b or abs(M)>c or(a+b+c)%2 or not abs(a-b)<=c<=a+b:return A(0)
 aa=(a+b-c)//2;u=(a-m)//2;v=(b+k)//2;w=(c-b+m)//2;x=(c-a-k)//2;lo=max(0,-w,-x);hi=min(aa,u,v)
 if lo>hi:return A(0)
 total=Fraction(0)
 for z in range(lo,hi+1):
  den=F[z]*F[aa-z]*F[u-z]*F[v-z]*F[w+z]*F[x+z];total+=Fraction(-1 if z%2 else 1,den)
 if not total:return A(0)
 num=(c+1)*F[aa]*F[(a-b+c)//2]*F[(-a+b+c)//2]
 den=F[(a+b+c)//2+1]
 for z in [(a+m)//2,(a-m)//2,(b+k)//2,(b-k)//2,(c+M)//2,(c-M)//2]:num*=F[z]
 value=A(Q(num*total.numerator**2,den*total.denominator**2)).sqrt();return value if total>0 else-value
p=[]
for ell in range(6):
 v=[Fraction((2*i-21)**ell)for i in range(22)]
 for u in p:
  scale=sum(x*y for x,y in zip(v,u))/sum(y*y for y in u);v=[x-scale*y for x,y in zip(v,u)]
 p.append(v)
f=[[A(Q(x.numerator,x.denominator))/A(Q(sum(y*y for y in v).numerator,sum(y*y for y in v).denominator)).sqrt()for x in v]for v in p]
@lru_cache(maxsize=None)
def direct_coefficient(J,S,T,ell):
 lo,hi=sorted((S,T))
 if not hi-lo<=2*ell<=hi+lo:return A(0)
 q=cg_ball(lo,hi,2*ell,lo,-lo).abs_lower();assert q>0
 value=A(0)
 for i in range(22):
  mi=2*i-21;value+=f[ell][i]*cg_ball(21,lo,J,mi,lo)*cg_ball(21,hi,J,mi,lo)
 # q is an interval lower endpoint above; use the full rigorous absolute
 # ball in the actual denominator, so no endpoint substitution is made.
 qball=abs(cg_ball(lo,hi,2*ell,lo,-lo))
 return value/(A(J+1).sqrt()*qball)
