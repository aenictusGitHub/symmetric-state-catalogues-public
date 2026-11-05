from fractions import Fraction as F
from collections import defaultdict
from math import comb
from pathlib import Path
import json
BASE=Path(__file__).parent

def mul(p,q):
 o=defaultdict(F)
 for a,c in p.items():
  for b,d in q.items():o[tuple(sorted(a+b))]+=c*d
 return {m:v for m,v in o.items() if v}
def purity(t):
 p=defaultdict(F)
 for i in range(t+1):
  for j in range(t+1):
   for a in range(8-t):
    for b in range(8-t):p[tuple(sorted((i+a,8+j+a,j+b,8+i+b)))]+=comb(t,i)*comb(t,j)*comb(7-t,a)*comb(7-t,b)
 return p
q={(i,8+i):F(comb(7,i)) for i in range(8)};q2=mul(q,q);q3=mul(q2,q)
S=defaultdict(F)
for t in range(6):
 for m,v in purity(t).items():S[m]+=21*(-1)**(5-t)*comb(5,t)*2**t*v
T=defaultdict(F,{m:F(267,7)*v for m,v in q3.items()})
for m,v in mul(q,S).items():T[m]-=v
rows=[]
for m,v in sorted(T.items()):
 if v:rows.append([[m.count(i) for i in range(8)],[m.count(i) for i in range(8,16)],str(v)])
(BASE/'target-exact.json').write_text(json.dumps({'n':7,'sector':5,'bound':'267/7','coefficients':rows},indent=2)+'\n')
print('Target',len(rows))
