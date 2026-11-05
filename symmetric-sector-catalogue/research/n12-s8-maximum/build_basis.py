"""Exact cubic SU(2) transvectants for N=12 sector extrema; see basis.json conventions."""
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

f=[{(r,):comb(12,r)} for r in range(13)]
g=[{(25-r,):(-1)**r*comb(12,r)} for r in range(13)]
pp={J:trans(f,f,12-J) for J in (0,2,4,6,8,10,12)}

def lincomb(forms,coeffs):
 out=[]
 for r in range(len(forms[0])):
  p=defaultdict(int)
  for form,coef in zip(forms,coeffs):padd(p,form[r],coef)
  out.append({m:v for m,v in p.items() if v})
 return normalize(out)

blocks=[];t0=time.time()
for kind,source in [('21',g),('30',f)]:
 for L2 in range(0,37,2):
  channels=[];forms=[]
  for J in pp:
   if not abs(2*J-12)<=L2<=2*J+12:continue
   F=trans(pp[J],source,(2*J+12-L2)//2)
   if not F[0]:continue
   channels.append(J);forms.append(F)
  if not forms:continue
  mons=sorted(set.union(*(set(F[0]) for F in forms)))
  A=s.Matrix([[F[0].get(m,0) for F in forms] for m in mons])
  _,piv=A.rref();forms=[forms[i] for i in piv];channels=[channels[i] for i in piv]
  W=[]
  def restriction(poly):
   out=defaultdict(int)
   mapping={int(r):tuple(v) for r,v in config['mapping'].items()}
   for mon,coef in poly.items():
    ex=[0]*4;rad=0
    for i in mon:
     if i%13 not in mapping:break
     v,c,p=mapping[i%13];ex[v+2*(i//13)]+=1;coef*=c;rad+=p
    else:out[tuple(ex)+(rad%2,)]+=coef*config['radicand']**(rad//2)
   return out
  for r in range(L2+1):
   rr=[restriction(F[r]) for F in forms];keys=set().union(*map(set,rr))
   W.extend([[R.get(m,0) for R in rr] for m in sorted(keys)])
  K=s.Matrix(W).nullspace() if W else [s.eye(len(forms))[:,i] for i in range(len(forms))]
  reduced=[lincomb(forms,list(col)) for col in K]
  if reduced:
   blocks.append({'kind':kind,'L2':L2,'originalChannels':channels,'originalMultiplicity':len(forms),'dimension':len(reduced),'weights':[str(s.Rational(1,comb(L2,r))) for r in range(L2+1)],'covariants':[[[[list(m),str(a)] for m,a in sorted(p.items())] for p in F] for F in reduced]})
  print(kind,L2,'multiplicity',len(forms),'kernel',len(reduced),'elapsed',round(time.time()-t0,1),flush=True)

out={'n':12,'convention':'f(z,w)=sum C(12,r)x_r z^(12-r)w^r; Tf=sum (-1)^r C(12,r)y_(12-r)z^(12-r)w^r; c_r=sqrt(C(12,r))x_r. Cubic rawbinaryform coefficients have normweights1/C(L2,r). kind21=(f,f) coupledTf; kind30=(f,f) coupledf.','variables':['x'+str(i) for i in range(13)]+['y'+str(i) for i in range(13)],'monomialFormat':'sorted list of variable indices, repetitions allowed','blocks':blocks}
(BASE/'basis.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
print('WROTE',BASE/'basis.json','bytes',(BASE/'basis.json').stat().st_size,flush=True)
