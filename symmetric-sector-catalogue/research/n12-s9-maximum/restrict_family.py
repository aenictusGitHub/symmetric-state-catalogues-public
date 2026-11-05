from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
from math import gcd,lcm
from functools import reduce
import json,sympy as s
dst=Path(__file__).resolve().parent
D=json.loads((dst/'basis-tangent.json').read_text());mons=set(tuple(m)for B in D['blocks']for cov in B['covariants']for comp in cov for m,co in comp);print('monomials',len(mons),flush=True)
def mul(A,B):
 C=[0]*4
 for r in range(4):
  for t in range(4):
   re=(r%2)+(t%2);im=(r//2)+(t//2);C[(re%2)+2*(im%2)]+=A[r]*B[t]*462**(re//2)*(-1)**(im//2)
 return C
states=[]
for zz in [1,2,3]:
 z=F(zz);w=(27*z*z+32)/(21*(27*z*z+64));v=z*(27*z*z+88)/(14*(27*z*z+64));a=-9*(w+v*z);e=-(6+9*z*z)/14;h=(-6*v+9*w*z)/14
 x={0:[0,a,0,0],2:[1,0,0,v],4:[0,w,z,0],6:[e,0,0,h],8:[0,-w,-z,0],10:[1,0,0,v],12:[0,-a,0,0]};den=lcm(*(F(t).denominator for a in x.values()for t in a));x={i:[int(t*den)for t in a]for i,a in x.items()}
 vals={}
 for mon in mons:
  A=[1,0,0,0]
  for i in mon:
   B=x.get(i%13,[0]*4);B=[B[0],B[1],(-1 if i>=13 else 1)*B[2],(-1 if i>=13 else 1)*B[3]];A=mul(A,B)
  vals[mon]=A
 states.append(vals)
 print('state',zz,x,flush=True)
out=[]
for B in D['blocks']:
 cov=B['covariants'];W=[]
 for vals in states:
  for r in range(B['L2']+1):
   VV=[]
   for v in cov:
    a=[0]*4
    for m,c in v[r]:
     ev=vals[tuple(m)]
     for i in range(4):a[i]+=int(c)*ev[i]
    VV.append(a)
   W.extend([[v[i]for v in VV]for i in range(4)])
 K=s.Matrix(W).nullspace();forms=[]
 for col in K:
  polys=[]
  for r in range(B['L2']+1):
   P=defaultdict(F)
   for coef,v in zip(col,cov):
    if not coef:continue
    for m,c in v[r]:P[tuple(m)]+=F(coef)*int(c)
   polys.append({m:c for m,c in P.items()if c})
  den=lcm(*(c.denominator for P in polys for c in P.values()));g=reduce(gcd,(int(c*den)for P in polys for c in P.values()));forms.append([[[list(m),str(int(c*den)//g)]for m,c in sorted(P.items())]for P in polys])
 if forms:
  b=dict(B);b['dimension']=len(forms);b['covariants']=forms;b['kernelRestriction']='Vanishing on three additional exact family representatives';out.append(b)
 print(B['kind'],B['L2'],len(cov),'->',len(forms),flush=True)
D['blocks']=out;(dst/'basis.json').write_text(json.dumps(D,separators=(',',':'))+'\n')
