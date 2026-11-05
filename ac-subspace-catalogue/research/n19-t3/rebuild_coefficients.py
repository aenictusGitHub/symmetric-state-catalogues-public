#!/usr/bin/env python3
"""Regenerate expanded coefficients from the finite radical recurrence."""
from pathlib import Path
import json
import sympy as s
from construction import n19
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
 path=HERE/'basis.json';d=json.loads(path.read_text());n,R,S,a,b,_=n19()
 for support,values,row in zip([R,S],[a,b],d['basis']):
  exact=['0']*(n+1);numeric=[[0.,0.]for _ in range(n+1)]
  for r,c in zip(support,values):
   exact[r]=str(c);z=complex(s.N(c,30));numeric[r]=[z.real,z.imag]
  row['exactCoefficients']=exact;row['coefficients']=numeric
 path.write_text(json.dumps(d,indent=2)+'\n')
 print('Expanded coefficients regenerated; rerun verifiers to update source hashes.')
