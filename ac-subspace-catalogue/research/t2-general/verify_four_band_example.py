"""Independent Dicke-splitting check of the exact four-band N=40 example."""
from pathlib import Path
from math import comb
import json
import sympy as s
from construct_exact import exact_basis

HERE = Path(__file__).resolve().parent
B = exact_basis(40)
matrices = [s.Matrix(3, 39, lambda a, b: B[a+b, c]*s.sqrt(
    s.Rational(comb(2, a)*comb(38, b), comb(40, a+b)))) for c in range(B.cols)]
count = 0
for i, A in enumerate(matrices):
    for j, D in enumerate(matrices):
        R = A*D.H-(s.eye(3)/3 if i == j else s.zeros(3))
        for value in R:
            assert s.simplify(s.expand_complex(value)) == 0
            count += 1
report = {'passed': True, 'n': 40, 'dimension': B.cols,
          'method': 'Independent direct Dicke splitting of the exact exported radical basis.',
          'crossEntriesChecked': count}
(HERE/'four-band-n40-independent.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
