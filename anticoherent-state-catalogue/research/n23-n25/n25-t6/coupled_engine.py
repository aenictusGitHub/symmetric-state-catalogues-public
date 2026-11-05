"""Rational CG-coupled averaged-word bases and raw Wigner recoupling factors.
All angular momenta are doubled integers. Standard-library arithmetic only.
"""
from fractions import Fraction as Q
from math import comb,factorial,prod,gcd,lcm,isqrt
from functools import lru_cache
F=[factorial(i) for i in range(512)]

@lru_cache(None)
def cg_highest(a,b,c):
 assert (a+b-c)%2==0 and abs(a-b)<=c<=a+b
 k=(a+b-c)//2;lo=max(0,k-b);hi=min(k,a);vals=[Q(1)]
 for q in range(lo,hi):vals.append(-vals[-1]*Q(b-k+q+1,a-q))
 den=lcm(*(v.denominator for v in vals));ints=[int(v*den) for v in vals];gg=gcd(*ints);ints=[v//gg for v in ints]
 assert ints[0]>0
 norm=sum(v*v*comb(a,q)*comb(b,k-q) for q,v in zip(range(lo,hi+1),ints))
 return k,lo,tuple(ints),norm

@lru_cache(None)
def sixj_raw(a,b,c,d,e,f):
 tris=[(a,b,c),(a,e,f),(d,b,f),(d,e,c)]
 if any((x+y+z)%2 or not abs(x-y)<=z<=x+y for x,y,z in tris):return Q(0)
 aa=[sum(t)//2 for t in tris];bb=[(a+b+d+e)//2,(a+c+d+f)//2,(b+c+e+f)//2];lo=max(aa);hi=min(bb)
 num=1;den=1
 for x,y,z in tris:num*=F[(x+y-z)//2]*F[(x-y+z)//2]*F[(-x+y+z)//2];den*=F[(x+y+z)//2+1]
 td=prod(F[lo-v] for v in aa)*prod(F[v-lo] for v in bb);term=F[lo+1]//td;total=(-1 if lo%2 else 1)*term
 for z in range(lo,hi):
  nn=(z+2)*prod(v-z for v in bb);dd=prod(z+1-v for v in aa);assert term*nn%dd==0;term=term*nn//dd;total+=(-1 if(z+1)%2 else 1)*term
 if total==0:return Q(0)
 sq=Q(total*total*num,den)
 for x,y,z in tris:sq*=cg_highest(x,y,z)[3]
 p=isqrt(sq.numerator);q=isqrt(sq.denominator)
 assert p*p==sq.numerator and q*q==sq.denominator,(a,b,c,d,e,f)
 return Q(-p if total<0 else p,q)

def physical_metric(n,word):
 from collections import Counter
 co=Counter(word);return Q(prod(comb(n,i)**v*factorial(v) for i,v in co.items()),factorial(len(word)))

def lower_vector(n,vec,step):
 from collections import Counter
 out={}
 for w,v in vec.items():
  for i,c in Counter(w).items():
   if i==0:continue
   z=list(w);z.remove(i);z.append(i-1);z=tuple(sorted(z));out[z]=out.get(z,0)+v*c*(n-i+1)
 assert all(v%step==0 for v in out.values())
 return {w:v//step for w,v in out.items() if v}

def coupled_vector(A,VA,B,VB,L):
 # VA[q] and VB[q] are divided-power descendants of highest columns.
 k,lo,cs,norm=cg_highest(A,B,L);out={}
 for q,v in zip(range(lo,lo+len(cs)),cs):
  for a,aa in VA[q].items():
   for b,bb in VB[k-q].items():out[a,b]=v*aa*bb
 return out,norm
