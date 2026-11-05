"""Exact cubic SU(2) transvectants for N=7 S6; see basis.json conventions."""
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

f=[{(r,):comb(7,r)} for r in range(8)]
g=[{(15-r,):(-1)**r*comb(7,r)} for r in range(8)]
pp={J:trans(f,f,7-J) for J in (1,3,5,7)}

def lincomb(forms,coeffs):
 out=[]
 for r in range(len(forms[0])):
  p=defaultdict(int)
  for form,coef in zip(forms,coeffs):padd(p,form[r],coef)
  out.append({m:v for m,v in p.items() if v})
 return normalize(out)

blocks=[];t0=time.time()
for kind,source in [('21',g),('30',f)]:
 for L2 in range(1,22,2):
  channels=[];forms=[]
  for J in pp:
   if not abs(2*J-7)<=L2<=2*J+7:continue
   F=trans(pp[J],source,(2*J+7-L2)//2)
   if not F[0]:continue
   channels.append(J);forms.append(F)
  if not forms:continue
  mons=sorted(set.union(*(set(F[0]) for F in forms)))
  A=s.Matrix([[F[0].get(m,0) for F in forms] for m in mons])
  _,piv=A.rref();forms=[forms[i] for i in piv];channels=[channels[i] for i in piv]
  W=[]
  for r in range(L2+1):
   mons=sorted(set.union(*(set(m for m in F[r] if all(i%8 in (0,7) for i in m)) for F in forms)))
   W.extend([[F[r].get(m,0) for F in forms] for m in mons])
  K=s.Matrix(W).nullspace() if W else [s.eye(len(forms))[:,i] for i in range(len(forms))]
  reduced=[lincomb(forms,list(col)) for col in K]
  if reduced:
   blocks.append({'kind':kind,'L2':L2,'originalChannels':channels,'originalMultiplicity':len(forms),'dimension':len(reduced),'weights':[str(s.Rational(1,comb(L2,r))) for r in range(L2+1)],'covariants':[[[[list(m),str(a)] for m,a in sorted(p.items())] for p in F] for F in reduced]})
  print(kind,L2,'multiplicity',len(forms),'kernel',len(reduced),'elapsed',round(time.time()-t0,1),flush=True)

out={'n':7,'convention':'f(z,w)=sum C(7,r)x_r z^(7-r)w^r; Tf=sum (-1)^r C(7,r)y_(7-r)z^(7-r)w^r; c_r=sqrt(C(7,r))x_r. Cubic rawbinaryform coefficients have normweights1/C(L2,r). kind21=(f,f) coupledTf; kind30=(f,f) coupledf.','variables':['x'+str(i) for i in range(8)]+['y'+str(i) for i in range(8)],'monomialFormat':'sorted list of variable indices, repetitions allowed','blocks':blocks}
(BASE/'basis.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
print('WROTE',BASE/'basis.json','bytes',(BASE/'basis.json').stat().st_size,flush=True)
