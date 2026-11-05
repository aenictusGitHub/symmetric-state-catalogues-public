from fractions import Fraction as F
from pathlib import Path
import json

def candidate(n,k,offset,bands):
 c=F(n*(n+2),12);W=[];R=[];ds=[]
 for i in range(k):
  rr=[offset+2*i+b for b in bands];m=[F(n,2)-r for r in rr]
  ww=[(c+m[(l+1)%3]*m[(l+2)%3])/((m[l]-m[(l+1)%3])*(m[l]-m[(l+2)%3]))for l in range(3)]
  if min(ww)<=0:return
  W.append(ww);R.append(rr)
 for i in range(k-1):
  A=[W[i][l]*W[i+1][l]*(n-R[i][l])*(n-R[i][l]-1)*(R[i][l]+1)*(R[i][l]+2)for l in range(3)]
  d=2*(A[0]*A[1]+A[1]*A[2]+A[2]*A[0])-sum(a*a for a in A)
  if d<0:return
  ds.append(str(d))
 return dict(n=n,k=k,offset=offset,bands=bands,support=R,weights=[[str(x)for x in w]for w in W],triangleDiscriminants=ds)

def find(n,k):
 gap=2*k+2;mx=n-2*(k-1)
 for b in range(gap,mx-gap+1,2):
  for c in range(b+gap,mx+1,2):
   for o in range(mx-c+1):
    r=candidate(n,k,o,[0,b,c])
    if r:return r

if __name__=='__main__':
 rows=[]
 for n in range(6,101):
  for k in range((n-2)//6,0,-1):
   r=find(n,k)
   if r:
    rows.append(r);old=max((n+3)//9,n//8-3)
    if k>old:print(n,k,old,r['bands'],r['offset'],flush=True)
    break
 Path(__file__).with_name('lower-unequal-triples.json').write_text(json.dumps(rows,indent=2)+'\n')
