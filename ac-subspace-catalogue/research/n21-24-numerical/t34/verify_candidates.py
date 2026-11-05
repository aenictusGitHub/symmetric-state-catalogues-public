#!/usr/bin/env python3
"""Independent 80-digit physical reduction checks of exported numerical bases.
No optimizer, spin-tensor generator, QR procedure or search routine is imported.
These are high-precision numerical checks, not exact existence certificates.
"""
from pathlib import Path
from math import comb
import json,hashlib
import mpmath as mp
HERE=Path(__file__).resolve().parent

def verify():
 mp.mp.dps=80;raw=(HERE/'candidates.json').read_bytes();data=json.loads(raw);reports=[];total=0
 for rec in data['candidates']:
  n,t,k=rec['n'],rec['t'],rec['dimension'];B=[[mp.mpc(mp.mpf(str(x)),mp.mpf(str(y)))for x,y in row['coefficients']]for row in rec['basis']]
  assert len(B)==k and all(len(v)==n+1 for v in B)
  gram=[[sum(mp.conj(B[a][r])*B[b][r]for r in range(n+1))-(1 if a==b else 0)for b in range(k)]for a in range(k)]
  ge=mp.sqrt(sum(abs(z)**2 for row in gram for z in row));assert ge<mp.mpf('1e-13')
  orders=[]
  for order in range(1,t+1):
   errors=[]
   for i in range(k):
    for j in range(k):
     deviation=[]
     for a in range(order+1):
      for b in range(order+1):
       z=sum(B[i][a+h]*mp.conj(B[j][b+h])*comb(n-order,h)*mp.sqrt(mp.mpf(comb(order,a)*comb(order,b))/(comb(n,a+h)*comb(n,b+h)))for h in range(n-order+1))
       z-=mp.mpf(1)/(order+1)if i==j and a==b else 0
       deviation.append(z);total+=1
     errors.append(mp.sqrt(sum(abs(z)**2 for z in deviation)))
   error=max(errors);assert error<mp.mpf('1e-13')
   orders.append({'order':order,'maxCrossReductionFrobenius':mp.nstr(error,25),'orderedComplexEntries':(order+1)**2*k*k})
  reports.append({'id':rec['id'],'n':n,'t':t,'dimension':k,'gramFrobenius':mp.nstr(ge,25),'reductions':orders})
 report={'passed':True,'evidence':'numerical','decimalPrecision':80,'candidatesSHA256':hashlib.sha256(raw).hexdigest(),'candidateCount':len(reports),'orderedComplexReducedEntriesChecked':total,'candidates':reports,'scope':'All saved decimal coefficients are checked independently through normalized-Dicke reductions at 80-digit arithmetic. This is a numerical validation, not an exact algebraic certificate and not a proof of global optimality.'}
 (HERE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');return report
if __name__=='__main__':print(json.dumps(verify(),indent=2))
