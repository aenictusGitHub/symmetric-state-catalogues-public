"""Rational SU(2) degree-four covariants; Dicke order r=0,...,10.
f(z,w)=sum C(10,r)*x_r*z^(10-r)*w^r, T f=sum (-1)^r C(10,r)*y_(10-r)*z^(10-r)*w^r.
y denotes independent conjugates. Physical c_r=sqrt(C(10,r))*x_r.
Covariant component polynomials are RAW binary-form coefficients: squarednorm=sum |F_r|²/C(2L,r).
"""
from math import comb,factorial,gcd,lcm
from functools import reduce
from collections import defaultdict
from pathlib import Path
import json,time
import sympy as s
BASE=Path(__file__).parent
config=json.loads((BASE/"config.json").read_text())

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

f=[{(r,):comb(10,r)} for r in range(11)]
g=[{(21-r,):(-1)**r*comb(10,r)} for r in range(11)]
pp={J:trans(f,f,10-J) for J in (0,2,4,6,8,10)}
dd={J:trans(g,g,10-J) for J in (0,2,4,6,8,10)}
pd={J:trans(f,g,10-J) for J in range(11)}

def lincomb(forms,coeffs):
 out=[]
 for r in range(len(forms[0])):
  p=defaultdict(int)
  for form,coef in zip(forms,coeffs):padd(p,form[r],coef)
  out.append({m:v for m,v in p.items() if v})
 return normalize(out)

blocks=[];t0=time.time()
for kind,source in [('22',dd),('31',pd),('40',pp)]:
 for L in range(21):
  channels=[];forms=[]
  for J in pp:
   for K in source:
    if not abs(J-K)<=L<=J+K:continue
    F=trans(pp[J],source[K],J+K-L)
    if not F[0]:continue
    channels.append([J,K]);forms.append(F)
  if not forms:continue
  mons=sorted(set.union(*(set(F[0]) for F in forms)))
  A=s.Matrix([[F[0].get(m,0) for F in forms] for m in mons])
  _,piv=A.rref();forms=[forms[i] for i in piv];channels=[channels[i] for i in piv]
  W=[]
  def restriction(poly):
   out=defaultdict(int);mapping={int(r):tuple(v) for r,v in config['mapping'].items()}
   for mon,coef in poly.items():
    ex=[0]*4;rad=0
    for i in mon:
     if i%11 not in mapping:break
     v,c,p=mapping[i%11];ex[v+2*(i//11)]+=1;coef*=c;rad+=p
    else:out[tuple(ex)+(rad%2,)]+=coef*config['radicand']**(rad//2)
   return out
  for r in range(2*L+1):
   rr=[restriction(F[r]) for F in forms];keys=set().union(*map(set,rr))
   W.extend([[R.get(m,0) for R in rr] for m in sorted(keys)])
  K=s.Matrix(W).nullspace()
  reduced=[lincomb(forms,list(col)) for col in K]
  if reduced:
   blocks.append({'kind':kind,'L2':2*L,'originalChannels':channels,'originalMultiplicity':len(forms),'dimension':len(reduced), 'weights':[str(s.Rational(1,comb(2*L,r))) for r in range(2*L+1)],'covariants':[[[[list(m),str(a)] for m,a in sorted(p.items())] for p in F] for F in reduced]})
  print(kind,L,'multiplicity',len(forms),'kernel',len(reduced),'elapsed',round(time.time()-t0,1),flush=True)

out={'convention':__doc__,'variables':['x'+str(i) for i in range(11)]+['y'+str(i) for i in range(11)],'monomialFormat':'sorted list of variable indices, repetitions allowed','blocks':blocks}
(BASE/'basis.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
print('WROTE',BASE/'basis.json','bytes',(BASE/'basis.json').stat().st_size,flush=True)
