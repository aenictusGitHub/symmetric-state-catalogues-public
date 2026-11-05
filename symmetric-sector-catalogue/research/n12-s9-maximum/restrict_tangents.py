from pathlib import Path
from collections import defaultdict
from fractions import Fraction as F
from math import gcd,lcm
from functools import reduce
import json,sympy as s
dst=Path(__file__).resolve().parent
D=json.loads((dst/'basis-initial.json').read_text());out=[];rad=462;base={0:(462,0),6:(1,1),12:(-462,0)};tangents=[{2:(1,1),4:(-14,0),8:(-14,0),10:(-1,1)},{2:(1,1),4:(14,0),8:(-14,0),10:(1,1)}]
def derivative(poly,tangent,imag):
 z=[0,0]
 for mon,c in poly:
  for a,i in enumerate(mon):
   if i%13 not in tangent:continue
   co,r=tangent[i%13];co*=int(c)*(-1 if imag and i>=13 else 1)
   for b,j in enumerate(mon):
    if a==b:continue
    if j%13 not in base:break
    v,p=base[j%13];co*=v;r+=p
   else:z[r%2]+=co*rad**(r//2)
 return z
for block in D['blocks']:
 cov=block['covariants'];W=[]
 for im,tangent in enumerate(tangents):
  for r in range(block['L2']+1):
   vals=[derivative(v[r],tangent,im)for v in cov];W.extend([[v[p]for v in vals]for p in range(2)])
 K=s.Matrix(W).nullspace();forms=[]
 for col in K:
  polys=[]
  for r in range(block['L2']+1):
   P=defaultdict(F)
   for coef,v in zip(col,cov):
    for m,c in v[r]:P[tuple(m)]+=F(coef)*F(c)
   polys.append({m:c for m,c in P.items()if c})
  den=lcm(*(c.denominator for P in polys for c in P.values()));g=reduce(gcd,(int(c*den)for P in polys for c in P.values()));forms.append([[[list(m),str(int(c*den)//g)]for m,c in sorted(P.items())]for P in polys])
 if forms:
  b=dict(block);b['dimension']=len(forms);b['covariants']=forms;b['kernelRestriction']='Value and both exact zero-Hessian tangent directions at the attaining representative';out.append(b)
 print(block['kind'],block['L2'],len(cov),'->',len(forms),flush=True)
D['blocks']=out;(dst/'basis-tangent.json').write_text(json.dumps(D,separators=(',',':'))+'\n')
