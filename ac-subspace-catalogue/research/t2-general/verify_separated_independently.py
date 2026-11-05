"""Reconstruct every ordered Dicke cross reduction from the exported radicals."""
from pathlib import Path
from math import comb
import json
import sympy as s
HERE = Path(__file__).resolve().parent
reports = []
for st in json.loads((HERE/'separated-bases.json').read_text())['subspaces']:
    n, k = st['n'], st['dimension']
    cs = [[s.sympify(z) for z in row['exactCoefficients']] for row in st['basis']]
    for i in range(k):
        for j in range(k):
            assert s.simplify(sum(cs[i][r]*s.conjugate(cs[j][r]) for r in range(n+1))-int(i==j)) == 0
            for a in range(3):
                for b in range(3):
                    v = sum(cs[i][a+r]*s.conjugate(cs[j][b+r])*comb(n-2,r)*
                            s.sqrt(s.Rational(comb(2,a)*comb(2,b),comb(n,a+r)*comb(n,b+r)))
                            for r in range(n-1))
                    assert s.simplify(v-s.Rational(int(i==j and a==b),3)) == 0
    reports.append({'n': n, 'dimension': k, 'passed': True,
                    'orthonormalityEntries': k*k, 'complexCrossEntriesExactlyVerified': 9*k*k})
report = {'passed': True, 'method': 'Direct exact Dicke splitting of every ordered cross operator; no construction function imported.',
          'records': reports, 'complexCrossEntriesExactlyVerified': sum(r['complexCrossEntriesExactlyVerified'] for r in reports)}
(HERE/'separated-independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
