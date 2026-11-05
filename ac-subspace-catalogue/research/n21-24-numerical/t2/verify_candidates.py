#!/usr/bin/env python3
"""Independent direct Dicke-splitting checks of exported numerical AC2 spans."""
from pathlib import Path
import json,hashlib,math
import numpy as np
HERE=Path(__file__).resolve().parent

def verify():
 records=[]
 for n in [21,22,23]:
  path=HERE/f'n{n}-k4.json';raw=path.read_bytes();d=json.loads(raw);assert d['found'] and d['dimension']==4 and d['t']==2
  vectors=[np.array([complex(z[0],z[1])for z in row['coefficients']])for row in d['basis']];assert len(vectors)==4 and all(v.shape==(n+1,)for v in vectors)
  B=np.column_stack(vectors);gram=B.conj().T@B;gramerr=float(np.max(np.abs(gram-np.eye(4))));assert gramerr<1e-12
  checks=[]
  for t in [1,2,3]:
   rows=[]
   for c in vectors:
    A=np.zeros((t+1,n-t+1),complex)
    for a in range(t+1):
     for b in range(n-t+1):A[a,b]=c[a+b]*math.sqrt(math.comb(t,a)*math.comb(n-t,b)/math.comb(n,a+b))
    rows.append(A)
   errors=[]
   for i,A in enumerate(rows):
    for j,D in enumerate(rows):
     delta=A@D.conj().T
     if i==j:delta-=np.eye(t+1)/(t+1)
     errors.append(float(np.linalg.norm(delta)))
   checks.append({'order':t,'orderedMatrices':16,'complexEntries':16*(t+1)**2,'maxCrossReductionFrobenius':max(errors)})
  assert checks[0]['maxCrossReductionFrobenius']<1e-12 and checks[1]['maxCrossReductionFrobenius']<1e-12
  assert checks[2]['maxCrossReductionFrobenius']>1e-6
  records.append({'id':d['id'],'n':n,'t':2,'dimension':4,'acceptedNumerically':True,'maxGramEntryError':gramerr,'reductions':checks,'sourceSHA256':hashlib.sha256(raw).hexdigest(),'evidence':'numerical only'})
 report={'passed':True,'acceptanceThreshold':1e-12,'candidates':records,'method':'Independent normalized-Dicke coefficient matrices, every ordered cross reduction, orders1–3; no search objective or optimizer code imported. Floating-point checks only.','nonexistenceClaims':False}
 (HERE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':print(json.dumps(verify(),indent=2))
