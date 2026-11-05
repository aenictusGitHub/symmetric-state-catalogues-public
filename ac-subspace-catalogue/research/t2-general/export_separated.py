"""Export the universal radical family in the catalogue's basis format."""
from pathlib import Path
import json
import sympy as s
from separated_triples import construction, verify

HERE = Path(__file__).resolve().parent
spaces = []
for n in (35, 42):
    data = construction(n)
    verify(data)
    basis = []
    for i, row in enumerate(data['rows']):
        cs = [s.Integer(0)] * (n+1)
        for r, w in zip(row['indices'], row['weights']):
            cs[r] = s.sqrt(s.Rational(w.numerator, w.denominator))
        terms = [s.latex(cs[r])+rf'|D_{{{n}}}^{{({r})}}\rangle' for r in row['indices']]
        basis.append({'exactCoefficients': list(map(str, cs)),
                      'coefficients': [[float(c), 0] for c in cs],
                      'stateLatex': rf'|u_{{{i+1}}}\rangle='+'+'.join(terms)})
    spaces.append({'id': f'separated-triples-n{n}-k{data["k"]}-t2',
                   'n': n, 't': 2, 'dimension': data['k'],
                   'label': 'Universal three-component radical construction',
                   'exactnessKind': 'radicals',
                   'provenance': 'Explicit member of the all-N separated-triple construction. This certifies a lower bound; maximality is not asserted.',
                   'derivationLatex': rf'k=\lfloor(N+3)/9\rfloor={data["k"]},\qquad h={data["h"]}',
                   'basis': basis})
(HERE/'separated-bases.json').write_text(json.dumps({'subspaces': spaces}, indent=2)+'\n')
