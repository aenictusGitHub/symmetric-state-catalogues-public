from fractions import Fraction as F
from collections import defaultdict
from math import comb
from pathlib import Path
import json
BASE=Path(__file__).parent
config=json.loads((BASE/'config.json').read_text());n=config['n'];k=config['k'];bound=F(config['bound'])
def mul(p,q):
 o=defaultdict(F)
 for a,c in p.items():
  for b,d in q.items():o[tuple(sorted(a+b))]+=c*d
 return {m:v for m,v in o.items() if v}
def purity(t):
 p=defaultdict(F)
 for i in range(t+1):
  for j in range(t+1):
   for a in range(n+1-t):
    for b in range(n+1-t):p[tuple(sorted((i+a,n+1+j+a,j+b,n+1+i+b)))]+=comb(t,i)*comb(t,j)*comb(n-t,a)*comb(n-t,b)
 return p
q={(i,n+1+i):F(comb(n,i)) for i in range(n+1)};q2=mul(q,q);q3=mul(q2,q)
S=defaultdict(F)
for t in range(k+1):
 for m,v in purity(t).items():S[m]+=comb(n,k)*(-1)**(k-t)*comb(k,t)*2**t*v
T=defaultdict(F,{m:bound*v for m,v in q3.items()})
for m,v in mul(q,S).items():T[m]-=v
if config['mode']=='min':T={m:-v for m,v in T.items()}
rows=[]
for m,v in sorted(T.items()):
 if v:rows.append([[m.count(i) for i in range(n+1)],[m.count(i) for i in range(n+1,2*n+2)],str(v)])
(BASE/'target-exact.json').write_text(json.dumps({'n':n,'sector':k,'bound':str(bound),'coefficients':rows},indent=2)+'\n')
print('Target',len(rows))
