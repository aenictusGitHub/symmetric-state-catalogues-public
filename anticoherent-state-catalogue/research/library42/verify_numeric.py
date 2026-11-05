"""Independently check every saved 85-digit library vector through its next order.
Requires mpmath. Does not read precomputed residuals or the source library.
"""
from pathlib import Path
import json
from functools import lru_cache
import mpmath as mp
mp.mp.dps=100
HERE=Path(__file__).resolve().parent
@lru_cache(None)
def factors(n,q):
 return [[mp.sqrt(mp.binomial(q,a)*mp.binomial(n-q,b)/mp.binomial(n,a+b)) for b in range(n-q+1)] for a in range(q+1)]
result=[]
for state in json.loads((HERE/'states.json').read_text())['states']:
 n,t=state['n'],state['achievedOrder'];c=[mp.mpc(x[0],x[1]) for x in state['highPrecisionCoefficients']]
 norm=mp.sqrt(sum(abs(x)**2 for x in c));c=[x/norm for x in c];checks=[]
 for q in range(1,t+2):
  f=factors(n,q);M=mp.matrix([[c[a+b]*f[a][b] for b in range(n-q+1)] for a in range(q+1)]);delta=M*M.H-mp.eye(q+1)/(q+1)
  residual=mp.sqrt(sum(abs(x)**2 for x in delta))
  assert (residual<mp.mpf('1e-60'))==(q<=t),(state['id'],q,mp.nstr(residual,30))
  checks.append({'t':q,'frobeniusDistance':mp.nstr(residual,50),'numericallyConsistent':q<=t})
 result.append({'id':state['id'],'n':n,'reportedOrder':t,'reductions':checks})
 print(state['id'],'pass',flush=True)
(HERE/'independent-numeric-check.json').write_text(json.dumps({'workingPrecisionDigits':100,'records':result,'passed':True},indent=2)+'\n')
print('PASS:',len(result),'reported orders and next-order failures checked.')
