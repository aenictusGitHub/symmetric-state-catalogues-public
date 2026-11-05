from fractions import Fraction as F
import json
from pathlib import Path

def candidate(n,k,h,o):
 c=F(n*(n+2),12); w=[];R=[]
 for i in range(k):
  rr=[o+2*i+l*h for l in range(3)];m=[F(n,2)-r for r in rr]
  ww=[(c+m[(l+1)%3]*m[(l+2)%3])/((m[l]-m[(l+1)%3])*(m[l]-m[(l+2)%3])) for l in range(3)]
  if min(ww)<=0:return
  w.append(ww);R.append(rr)
 for i in range(k-1):
  A=[w[i][l]*w[i+1][l]*(n-R[i][l])*(n-R[i][l]-1)*(R[i][l]+1)*(R[i][l]+2) for l in range(3)]
  delta=2*(A[0]*A[1]+A[1]*A[2]+A[2]*A[0])-sum(a*a for a in A)
  if delta<0:return
 return {'n':n,'k':k,'h':h,'offset':o,'weights':[[str(v)for v in row]for row in w]}

if __name__=='__main__':
 rows=[]
 for n in range(6,101):
  best=None
  for k in range((n+2)//6,0,-1):
   for h in range(2*k+2,(n-2*(k-1))//2+1,2):
    for o in range(n-2*(k-1)-2*h+1):
     best=candidate(n,k,h,o)
     if best:break
    if best:break
   if best:break
  if best:
   rows.append(best)
   lower=max((n+3)//9,n//8-3)
   if best['k']>lower:print(n,best['k'],lower,best['h'],best['offset'],flush=True)
 Path(__file__).with_name('lower-triples-scan.json').write_text(json.dumps(rows,indent=2)+'\n')
