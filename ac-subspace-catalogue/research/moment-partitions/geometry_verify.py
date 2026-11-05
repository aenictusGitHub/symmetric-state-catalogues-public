#!/usr/bin/env python3
"""Exact checks of graph constructions, including every two-body cross entry."""
from pathlib import Path
import json, math
import sympy as s
from geometry_graph import construct
HERE=Path(__file__).resolve().parent
report=[];examples=[]
for n in (12,21):
    data,Q=construct(n);Q=Q.applyfunc(s.simplify);k=data['k']
    indices=[int(s.Rational(n,2)-m) for m in data['levels']]
    full=s.zeros(n+1,k)
    for i,r in enumerate(indices):full[r,:]=Q[i,:]
    # Dicke splitting map, independent of the graph construction.
    split=[]
    for col in range(k):
        M=s.zeros(3,n-1)
        for a in range(3):
            for b in range(n-1):
                M[a,b]=full[a+b,col]*s.sqrt(s.Rational(math.comb(2,a)*math.comb(n-2,b),math.comb(n,a+b)))
        split.append(M)
    for u in range(k):
        for v in range(k):
            expected=s.eye(3)/3 if u==v else s.zeros(3)
            assert all(s.simplify(x)==0 for x in split[u]*split[v].T-expected),(n,u,v)
    examples.append(dict(n=n,t=2,k=k,dickeIndices=indices,components=[[str(Q[i,a]) for i in range(Q.rows)] for a in range(k)],status='exact radical AC2 basis; no maximality claim'))
    report.append(dict(n=n,k=k,passed=True,exactOrderedCrossEntries=9*k*k,verification='Dicke splitting, all cross reductions equal delta_ab I3/3'))
(HERE/'geometry_examples.json').write_text(json.dumps(examples,indent=2)+'\n')
(HERE/'geometry_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
