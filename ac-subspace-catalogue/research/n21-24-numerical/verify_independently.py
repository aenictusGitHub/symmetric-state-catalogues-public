"""Independent high-precision normalized-Dicke cross-reduction audit.
Only exported decimal coefficients are read. No optimization code is imported.
"""
from pathlib import Path
import json,hashlib,math
import mpmath as mp
HERE=Path(__file__).resolve().parent

def verify(directory=None):
 root=Path(directory or HERE);raw=(root/'bases.json').read_bytes();d=json.loads(raw);reports=[];count=0
 with mp.workdps(80):
  tolerance=mp.mpf('1e-10')
  for c in d['candidates']:
   n,t,k=c['n'],c['t'],c['dimension']
   rows=[[mp.mpc(str(a),str(b))for a,b in v['coefficients']]for v in c['basis']]
   assert len(rows)==k and all(len(v)==n+1 for v in rows)
   gram=mp.sqrt(sum(abs(sum(mp.conj(rows[a][r])*rows[b][r]for r in range(n+1))-(1 if a==b else 0))**2 for a in range(k)for b in range(k)))
   assert gram<tolerance,(c['id'],'Gram',str(gram))
   errors=[];entries=0
   for order in range(1,t+1):
    coeff=[[mp.sqrt(mp.mpf(math.comb(order,a)*math.comb(n-order,h))/math.comb(n,a+h))for h in range(n-order+1)]for a in range(order+1)]
    mats=[[[row[a+h]*coeff[a][h]for h in range(n-order+1)]for a in range(order+1)]for row in rows]
    maximum=mp.mpf(0)
    for i in range(k):
     for j in range(k):
      norm=mp.mpf(0)
      for a in range(order+1):
       for b in range(order+1):
        value=sum(mats[i][a][h]*mp.conj(mats[j][b][h])for h in range(n-order+1))
        if i==j and a==b:value-=mp.mpf(1)/(order+1)
        norm+=abs(value)**2;entries+=1
      maximum=max(maximum,mp.sqrt(norm))
    assert maximum<tolerance,(c['id'],order,str(maximum))
    errors.append({'order':order,'maxOrderedCrossReductionFrobenius':mp.nstr(maximum,20)})
   count+=entries
   reports.append({'id':c['id'],'n':n,'t':t,'dimension':k,'passed':True,'orthonormalityFrobenius':mp.nstr(gram,20),'allOrders':errors,'maxCrossReductionFrobenius':max(float(x['maxOrderedCrossReductionFrobenius'])for x in errors),'scalarEntriesChecked':entries,'exactExistenceCertificate':False})
 report={'passed':True,'basisFileSHA256':hashlib.sha256(raw).hexdigest(),'precisionDigits':80,'acceptanceThreshold':'1e-10','candidates':reports,'candidateCount':len(reports),'orderedScalarEntriesChecked':count,'method':'Direct normalized-Dicke splitting at 80 decimal digits from the saved decimal strings; every ordered reduced matrix at every order up to the claimed order, and full Gram matrix. No search routines or spin tensors imported.','scope':'Independent numerical consistency audit. Neither exact existence nor global maximality is certified.'}
 return report
if __name__=='__main__':
 result=verify();(HERE/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'passed':True,'candidateCount':result['candidateCount'],'scalarEntries':result['orderedScalarEntriesChecked']},indent=2))
