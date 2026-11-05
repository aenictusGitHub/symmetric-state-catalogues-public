"""Exact rational evaluations of the averaged holomorphic invariant basis."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,lcm
import json,numpy as np,time
ROOT=Path(__file__).parent;cert=json.load(open(ROOT/'sos-certificate.json'));raw=json.load(open(ROOT/'holomorphic-basis.json'))['blocks'];variables=[(ib,a,b)for ib,B in enumerate(raw)for b in range(B['dimension'])for a in range(b+1)];chosen=[variables[j]for j in cert['invariantBasisColumns']];scales=cert['invariantBasisScales'];points=cert['evaluationPoints'];norms=[sum(comb(11,i)*(a*a+b*b)for i,(a,b)in enumerate(x))**5 for x in points];out=[[None]*290 for _ in points];start=time.time()
for ib in sorted(set(v[0]for v in chosen)):
 B=raw[ib];L=B['L2'];v=np.load(ROOT/'holo-integer-evals'/f'{ib}.npy');E=v[...,0].astype(object)+v[...,1].astype(object)*(1<<32)+v[...,2].astype(object)*(1<<64);den=lcm(*(6*comb(5,t)*comb(L,r)for r in range(L+1)for t in range(6)))
 for col,(jb,a,b)in enumerate(chosen):
  if ib!=jb:continue
  val=np.zeros(len(points),dtype=object)
  for r in range(L+1):
   for t in range(6):val+=(E[:,a,r,t,0]*E[:,b,r,t,0]+E[:,a,r,t,1]*E[:,b,r,t,1])*(den//(6*comb(5,t)*comb(L,r)))
  fac=F(1 if a==b else 2,den)/F(scales[ib][a])/F(scales[ib][b])
  for i in range(len(points)):out[i][col]=str(F(int(val[i]),norms[i])*fac)
 print('BASIS BLOCK',ib,'seconds',round(time.time()-start),flush=True)
assert all(v is not None for row in out for v in row);(ROOT/'basis-evaluations.json').write_text(json.dumps(out,separators=(',',':'))+'\n');print('DONE exact invariant evaluations',flush=True)
