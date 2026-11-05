#!/usr/bin/env python3
"""Independent exact verification of the seven newly imported Julia states.
Run from any directory; writes imported-states-independent-audit.json beside it.
Requires SymPy. This reads the saved catalogue, not the extraction code.
"""
from pathlib import Path
import json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parents[2]

def reduction(c,n,t):
 A=s.Matrix(t+1,n-t+1,lambda a,b:c[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(n-t,b)/s.binomial(n,a+b)))
 return A*A.H

def check_zero(matrix):
 for a in range(matrix.rows):
  for b in range(matrix.cols):
   value=s.simplify(s.expand_complex(matrix[a,b]))
   if value!=0:return {'row':a,'column':b,'value':str(value)}
 return None

raw=(ROOT/'data/catalogue.json').read_bytes();data=json.loads(raw);results=[]
for state in data['states']:
 if not state['id'].startswith('acliege-'):continue
 n=state['n'];t=state['achievedOrder'];c=s.Matrix([s.sympify(x) for x in state['exactCoefficients']])
 assert len(c)==n+1 and s.simplify((c.H*c)[0])==1,state['id']
 violation=check_zero(reduction(c,n,t)-s.eye(t+1)/(t+1));assert violation is None,(state['id'],violation)
 next_violation=check_zero(reduction(c,n,t+1)-s.eye(t+2)/(t+2));assert next_violation is not None,state['id']
 results.append({'id':state['id'],'n':n,'exactAchievedOrder':t,'normalizationVerified':True,'highestClaimedReductionVerifiedExactly':True,'lowerOrdersFollowByPartialTrace':True,'nextOrderCounterexample':next_violation})
 print(state['id'],'verified',flush=True)
report={'passed':True,'catalogueSha256':hashlib.sha256(raw).hexdigest(),'statesChecked':len(results),'method':'Independent exact normalized-Dicke partial trace. The reduction at the claimed order is exactly I/(t+1); every lower order follows by tracing it. A nonzero entry at the next order establishes the exact achieved order of this representative.','states':results,'scope':'This verifies each imported pure state, not global maximality of the order at its qubit number.'}
Path(__file__).with_name('imported-states-independent-audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'passed':True,'statesChecked':len(results)}))
