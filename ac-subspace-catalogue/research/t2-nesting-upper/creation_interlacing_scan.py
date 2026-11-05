from fractions import Fraction as F
from pathlib import Path
import json

def upper(n,S2,kind):
 s=F(S2,2);J=F(n+(1 if kind=='creation' else -1),2)
 alpha=F(n+(3 if kind=='creation' else -1),3)
 L2range=range(abs(int(2*J)-S2),int(2*J)+S2+1,2)
 eig=[( (F(L2,2)*(F(L2,2)+1)-J*(J+1)-s*(s+1))/2,L2+1)for L2 in L2range]
 if kind=='creation':lo=-alpha*(s+1)/2;hi=alpha*s/2;lm=S2;hm=S2+2
 else:lo=-alpha*s/2;hi=alpha*(s+1)/2;lm=S2+2;hm=S2
 total=2*(S2+1)
 candidates=[(sum(m for e,m in eig if e<=lo)//lm,'low'),(sum(m for e,m in eig if e>=hi)//hm,'high'),(sum(m for e,m in eig if e>=lo)//total,'all-above-low'),(sum(m for e,m in eig if e<=hi)//total,'all-below-high')]
 b,which=min(candidates)
 return b,which
if __name__=='__main__':
 improvements=[];records=[]
 for n in range(8,201):
  old=min((n-2)//4,(n+5)//5)
  best=(old,None,None,None)
  for S2 in range(1,2*n+5):
   for kind in ['creation','derivative']:
    b,which=upper(n,S2,kind)
    if b<best[0]:best=(b,kind,S2,which)
  if best[0]<old:improvements.append({'n':n,'old':old,'new':best[0],'kind':best[1],'spin':str(F(best[2],2)),'branch':best[3]})
 print(improvements)
 Path(__file__).with_name('creation-interlacing-scan.json').write_text(json.dumps({'scope':'Exact finite test, not an all-N optimality theorem','nRange':[8,200],'auxiliaryTwiceSpinRange':'1..2*N+4','improvements':improvements},indent=2)+'\n')
