#!/usr/bin/env python3
"""Exact diagonal-moment check for the odd-N49 AC6 positive control."""
from fractions import Fraction
import json
from pathlib import Path
p=Path(__file__).with_name('positive-control.json')
d=json.loads(p.read_text());n=d['n'];support=d['support'];w=list(map(Fraction,d['rationalProbabilities']))
assert n==49 and len(w)==len(support) and all(x>0 for x in w)
assert min(b-a for a,b in zip(support,support[1:]))==7
for k in range(7):
    assert sum(x*r**k for x,r in zip(w,support)) == sum(Fraction(r**k,n+1) for r in range(n+1))
print('Exact positive rational moments through order six verified; support gap seven kills all corresponding off-diagonal moments.')
