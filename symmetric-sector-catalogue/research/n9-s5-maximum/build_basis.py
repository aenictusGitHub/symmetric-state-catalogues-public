"""Rational SU(2) degree-four covariants; Dicke order r=0,...,9.
f(z,w)=sum C(9,r)*x_r*z^(9-r)*w^r, T f=sum (-1)^r C(9,r)*y_(9-r)*z^(9-r)*w^r.
y denotes independent conjugates. Physical c_r=sqrt(C(9,r))*x_r.
Covariant component polynomials are RAW binary-form coefficients: squarednorm=sum |F_r|²/C(2L,r).
"""
from math import comb,factorial,gcd,lcm
from functools import reduce
from collections import defaultdict
from pathlib import Path
import json,time
import sympy as s
BASE=Path(__file__).parent


def padd(out,p,mult=1):
 for mon,v in p.items():out[mon]+=mult*v

def pmul(p,q):
 out=defaultdict(int)
 for m,a in p.items():
  for n,b in q.items():out[tuple(sorted(m+n))]+=a*b
 return {m:a for m,a in out.items() if a}

def normalize(F):
 vals=[a for p in F for a in p.values()]
 if not vals:return F
 den=reduce(lcm,[s.denom(a) for a in vals],1)
 g=reduce(gcd,[int(a*den) for a in vals])
 first=next(iter(next(p for p in F if p).values()))
 if first<0:g=-g
 return [{m:int(a*den)//g for m,a in p.items() if a} for p in F]

def ff(n,k):return factorial(n)//factorial(n-k) if 0<=k<=n else 0

def trans(F,G,r):
 n=len(F)-1;m=len(G)-1
 out=[defaultdict(int) for _ in range(n+m-2*r+1)]
 for k in range(r+1):
  for a,p in enumerate(F):
   ca=(-1)**k*comb(r,k)*ff(n-a,r-k)*ff(a,k)
   if not ca:continue
   for b,q in enumerate(G):
    cb=ff(m-b,k)*ff(b,r-k)
    if cb:
     padd(out[a+b-r],pmul(p,q),ca*cb)
 return normalize([{m:a for m,a in p.items() if a} for p in out])

f=[{(r,):comb(9,r)} for r in range(10)]
g=[{(19-r,):(-1)**r*comb(9,r)} for r in range(10)]
pp={J:trans(f,f,9-J) for J in (1,3,5,7,9)}
dd={J:trans(g,g,9-J) for J in (1,3,5,7,9)}
pd={J:trans(f,g,9-J) for J in range(10)}

def lincomb(forms,coeffs):
 out=[]
 for r in range(len(forms[0])):
  p=defaultdict(int)
  for form,coef in zip(forms,coeffs):padd(p,form[r],coef)
  out.append({m:v for m,v in p.items() if v})
 return normalize(out)

# Exact field Q(sqrt(61)), with t satisfying its selected quadratic factor.
K=s.QQ.algebraic_field(s.sqrt(61));zero=K.zero;one=K.one;rad=K.from_sympy(s.sqrt(61))
A=(991+159*rad)/1025;B=(15+rad)/820

def tadd(x,y):return (x[0]+y[0],x[1]+y[1])
def tmul(x,y):
 c=x[1]*y[1];return (x[0]*y[0]-B*c,x[0]*y[1]+x[1]*y[0]-A*c)
def tscale(x,c):return (c*x[0],c*x[1])
def tp(x,n):
 a=(one,zero)
 for _ in range(n):a=tmul(a,x)
 return a
T=(zero,one);Z=tscale(T,5*(15-rad));U=(zero,zero)
for power,c in enumerate([-257552,-31773836,146668933,73578600]):U=tadd(U,tscale(tp(T,power),K(c)/45433661))
from functools import lru_cache
@lru_cache(None)
def monvalue(mon):
 counts=[mon.count(i)for i in [0,3,6,9]]
 if sum(counts)!=len(mon):return (zero,)*4
 nw=counts[1]+counts[3];v=tmul(tp(T,counts[2]),tmul(tp(Z,counts[3]),tp(U,nw//2)))
 return v+(zero,zero) if nw%2==0 else (zero,zero)+v

def restriction(poly):
 v=[zero]*4
 for mon,c in poly.items():
  a=monvalue(tuple(sorted(i%10 for i in mon)))
  v=[x+K(c)*y for x,y in zip(v,a)]
 return v

def kernel(rows,n):
 M=[r[:]for r in rows if any(r)];piv=[];r=0
 for col in range(n):
  pivot=next((i for i in range(r,len(M))if M[i][col]),None)
  if pivot is None:continue
  M[r],M[pivot]=M[pivot],M[r];v=M[r][col];M[r]=[c/v for c in M[r]]
  for i in range(len(M)):
   if i!=r and M[i][col]:
    a=M[i][col];M[i]=[x-a*y for x,y in zip(M[i],M[r])]
  piv.append(col);r+=1
  if r==len(M):break
 out=[]
 for j in range(n):
  if j in piv:continue
  v=[zero]*n;v[j]=one
  for i,p in enumerate(piv):v[p]=-M[i][j]
  out.append(v)
 return out

def serialize(x):
 # ANP coefficients are high-to-low in sqrt(61).
 coeff=list(x.to_list());coeff=[K.dom.zero]*(2-len(coeff))+coeff
 return [str(coeff[1]),str(coeff[0])]
blocks=[]
for kind,source in [('22',dd),('31',pd),('40',pp)]:
 for L2 in range(0,37,2):
  channels=[];forms=[]
  for J in pp:
   for JJ in source:
    if not abs(J-JJ)<=L2//2<=J+JJ:continue
    FF=trans(pp[J],source[JJ],J+JJ-L2//2)
    if not FF[0]:continue
    channels.append([J,JJ]);forms.append(FF)
  if not forms:continue
  mons=sorted(set.union(*(set(FF[0])for FF in forms)))
  _,piv=s.Matrix([[FF[0].get(m,0)for FF in forms]for m in mons]).rref();forms=[forms[i]for i in piv]
  W=[]
  for r in range(L2+1):
   values=[restriction(FF[r])for FF in forms]
   W.extend([[v[i]for v in values]for i in range(4)])
  ker=kernel(W,len(forms));reduced=[]
  for v in ker:
   components=[]
   for r in range(L2+1):
    p={}
    for FF,c in zip(forms,v):
     for mon,coef in FF[r].items():p[mon]=p.get(mon,zero)+c*K(coef)
    components.append([[list(m),serialize(c)]for m,c in sorted(p.items())if c])
   reduced.append(components)
  if reduced:blocks.append({'kind':kind,'L2':L2,'dimension':len(reduced),'weights':[str(s.Rational(1,comb(L2,r)))for r in range(L2+1)],'covariants':reduced})
  print(kind,L2,'multiplicity',len(forms),'kernel',len(ker),flush=True)
(BASE/'basis.json').write_text(json.dumps({'n':9,'fieldRadicand':61,'blocks':blocks},separators=(',',':'))+'\n')
