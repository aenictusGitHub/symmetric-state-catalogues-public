#!/usr/bin/env python3
"""Independent 80-decimal Dicke-splitting check of saved numerical vectors.
No optimization or tensor-generation code is imported. These are numerical
residual checks, not symbolic existence certificates.
"""
from pathlib import Path
import json,hashlib,math
import mpmath as mp
HERE=Path(__file__).resolve().parent

def verify(directory=None):
 folder=Path(directory or HERE);raw=(folder/'summary.json').read_bytes();d=json.loads(raw);records=[];count=0
 with mp.workdps(80):
  for rec in d['candidates']:
   n,t=rec['n'],rec['t'];assert rec['dimension']==1
   c=[mp.mpc(mp.mpf(str(a)),mp.mpf(str(b)))for a,b in rec['basis'][0]['coefficients']];assert len(c)==n+1
   norm=sum(abs(z)**2 for z in c);assert abs(norm-1)<mp.mpf('1e-12')
   checks=[]
   for r in range(1,t+1):
    sq=mp.mpf(0);mx=mp.mpf(0)
    for a in range(r+1):
     for b in range(r+1):
      val=sum(c[a+h]*c[b+h].conjugate()*math.comb(n-r,h)*mp.sqrt(mp.mpf(math.comb(r,a)*math.comb(r,b))/(math.comb(n,a+h)*math.comb(n,b+h)))for h in range(n-r+1))
      if a==b:val-=mp.mpf(1)/(r+1)
      sq+=abs(val)**2;mx=max(mx,abs(val));count+=1
    frob=mp.sqrt(sq);assert frob<mp.mpf('1e-10')
    checks.append({'t':r,'frobeniusResidual':mp.nstr(frob,25),'maximumEntryResidual':mp.nstr(mx,25),'entriesChecked':(r+1)**2})
   records.append({'id':rec['id'],'n':n,'t':t,'dimension':1,'normalizationResidual':mp.nstr(abs(norm-1),25),'reductions':checks,'numericalOnly':True})
 result={'passed':True,'precisionDecimalDigits':80,'independentEntryCount':count,'numericalOnly':True,'scope':'High-precision evaluation of exported decimal coefficients; no exact existence or optimality claim.','summarySHA256':hashlib.sha256(raw).hexdigest(),'records':records}
 (folder/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
 result=verify();print(json.dumps({'passed':result['passed'],'candidates':len(result['records']),'entries':result['independentEntryCount']}))
