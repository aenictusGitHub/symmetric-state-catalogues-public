#!/usr/bin/env python3
"""Independent exact attainment checks for recovered symmetric states.

Run with Python 3 and SymPy: python3 verify_exact_witnesses.py
The result file is written beside this script. All coefficients and targets
are exact rational/radical expressions. This script checks attainment; it
does not establish any global minimum or maximum. A bound proof is separate.
"""
from functools import lru_cache
from math import comb
from pathlib import Path
import hashlib
import json
import platform
import sys
import time
import sympy as s

SPECS = [{'id': 'n7-s6-max',
  'n': 7,
  'k': 6,
  'mode': 'max',
  'value': '49',
  'exactCoefficients': {'0': 'sqrt(2)/3', '3': 'sqrt(14)/6', '6': 'sqrt(14)/6'},
  'separateSharpBoundCertificate': True},
 {'id': 'n8-s4-min',
  'n': 8,
  'k': 4,
  'mode': 'min',
  'value': '310/21',
  'exactCoefficients': {'2': '4/(3*sqrt(3))', '5': '2/(3*sqrt(3))', '8': 'sqrt(7)/(3*sqrt(3))'},
  'separateSharpBoundCertificate': True},
 {'id': 'n8-s6-max',
  'n': 8,
  'k': 6,
  'mode': 'max',
  'value': '364/3',
  'exactCoefficients': {'0': 'sqrt(30)/12', '4': 'sqrt(21)/6', '8': 'sqrt(30)/12'},
  'separateSharpBoundCertificate': True},
 {'id': 'n8-s7-max',
  'n': 8,
  'k': 7,
  'mode': 'max',
  'value': '592/7',
  'exactCoefficients': {'0': 'sqrt(7)/(3*sqrt(3))', '3': '2/(3*sqrt(3))', '6': '4/(3*sqrt(3))'},
  'separateSharpBoundCertificate': True},
 {'id': 'n9-s6-min',
  'n': 9,
  'k': 6,
  'mode': 'min',
  'value': '420/11',
  'exactCoefficients': {'0': 'sqrt(286)/22', '8': '3*sqrt(22)/22'},
  'separateSharpBoundCertificate': False},
 {'id': 'n10-s5-max',
  'n': 10,
  'k': 5,
  'mode': 'max',
  'value': '252',
  'exactCoefficients': {'0': '1'},
  'separateSharpBoundCertificate': True},
 {'id': 'n10-s7-max',
  'n': 10,
  'k': 7,
  'mode': 'max',
  'value': '320',
  'exactCoefficients': {'0': '1/2', '5': 'sqrt(2)/2', '10': '1/2'},
  'separateSharpBoundCertificate': True},
 {'id': 'n10-s8-max',
  'n': 10,
  'k': 8,
  'mode': 'max',
  'value': '405',
  'exactCoefficients': {'1': 'sqrt(5)/4', '5': 'sqrt(6)/4', '9': 'sqrt(5)/4'},
  'separateSharpBoundCertificate': True},
 {'id': 'n10-s8-min',
  'n': 10,
  'k': 8,
  'mode': 'min',
  'value': '21',
  'exactCoefficients': {'0': 'sqrt(6)/3', '9': 'sqrt(3)/3'},
  'separateSharpBoundCertificate': False},
 {'id': 'n10-s9-max',
  'n': 10,
  'k': 9,
  'mode': 'max',
  'value': '880/3',
  'exactCoefficients': {'2': '1/2', '5': 'sqrt(2)/2', '8': '1/2'},
  'separateSharpBoundCertificate': False},
 {'id': 'n11-s5-max',
  'n': 11,
  'k': 5,
  'mode': 'max',
  'value': '462',
  'exactCoefficients': {'0': '1'},
  'separateSharpBoundCertificate': True},
 {'id': 'n11-s6-max',
  'n': 11,
  'k': 6,
  'mode': 'max',
  'value': '462',
  'exactCoefficients': {'0': '1'},
  'separateSharpBoundCertificate': False},
 {'id': 'n11-s7-max',
  'n': 11,
  'k': 7,
  'mode': 'max',
  'value': '250060/517',
  'exactCoefficients': {'1': 'sqrt(8366)/94', '8': 'sqrt(470)/94'},
  'separateSharpBoundCertificate': False},
 {'id': 'n11-s8-max',
  'n': 11,
  'k': 8,
  'mode': 'max',
  'value': '5271/11',
  'exactCoefficients': {'2': '1'},
  'separateSharpBoundCertificate': False},
 {'id': 'n11-s8-min',
  'n': 11,
  'k': 8,
  'mode': 'min',
  'value': '77',
  'exactCoefficients': {'0': 'sqrt(570)/30', '10': 'sqrt(330)/30'},
  'separateSharpBoundCertificate': False},
 {'id': 'n11-s9-max',
  'n': 11,
  'k': 9,
  'mode': 'max',
  'value': '134645/239',
  'exactCoefficients': {'0': 'sqrt(67159)/478', '8': '15*sqrt(717)/478'},
  'separateSharpBoundCertificate': False},
 {'id': 'n11-s10-max',
  'n': 11,
  'k': 10,
  'mode': 'max',
  'value': '19459/35',
  'exactCoefficients': {'0': 'sqrt(2282)/70', '9': 'sqrt(2618)/70'},
  'separateSharpBoundCertificate': False},
 {'id': 'n11-s10-min',
  'n': 11,
  'k': 10,
  'mode': 'min',
  'value': '11',
  'exactCoefficients': {'0': '1'},
  'separateSharpBoundCertificate': False},
 {'id': 'n12-s4-min',
  'n': 12,
  'k': 4,
  'mode': 'min',
  'value': '99',
  'exactCoefficients': {'1': 'sqrt(7)/5', '6': 'sqrt(11)/5', '11': '-sqrt(7)/5'},
  'separateSharpBoundCertificate': False},
 {'id': 'n12-s5-max',
  'n': 12,
  'k': 5,
  'mode': 'max',
  'value': '792',
  'exactCoefficients': {'0': '1'},
  'separateSharpBoundCertificate': True},
 {'id': 'n12-s6-max',
  'n': 12,
  'k': 6,
  'mode': 'max',
  'value': '924',
  'exactCoefficients': {'0': '1'},
  'separateSharpBoundCertificate': False},
 {'id': 'n12-s7-max',
  'n': 12,
  'k': 7,
  'mode': 'max',
  'value': '946',
  'exactCoefficients': {'1': '1'},
  'separateSharpBoundCertificate': True},
 {'id': 'n12-s8-max',
  'n': 12,
  'k': 8,
  'mode': 'max',
  'value': '15405/11',
  'exactCoefficients': {'2': 'sqrt(7)/4', '6': 'sqrt(2)/4', '10': 'sqrt(7)/4'},
  'separateSharpBoundCertificate': False},
 {'id': 'n12-s9-max',
  'n': 12,
  'k': 9,
  'mode': 'max',
  'value': '1280',
  'exactCoefficients': {'0': '1/2', '6': 'sqrt(2)/2', '12': '-1/2'},
  'separateSharpBoundCertificate': False},
 {'id': 'n12-s10-max',
  'n': 12,
  'k': 10,
  'mode': 'max',
  'value': '1606',
  'exactCoefficients': {'1': 'sqrt(7)/5', '6': 'sqrt(11)/5', '11': '-sqrt(7)/5'},
  'separateSharpBoundCertificate': True},
 {'id': 'n12-s10-min',
  'n': 12,
  'k': 10,
  'mode': 'min',
  'value': '594/19',
  'exactCoefficients': {'0': 'sqrt(247)/19', '11': 'sqrt(114)/19'},
  'separateSharpBoundCertificate': False},
 {'id': 'n12-s11-max',
  'n': 12,
  'k': 11,
  'mode': 'max',
  'value': '1270976/1191',
  'exactCoefficients': {'0': 'sqrt(680458)/1191', '10': '13*sqrt(4367)/1191'},
  'separateSharpBoundCertificate': False}]
INPUT_METADATA = {'files': [{'path': 'research/search/recovered-witnesses.json',
            'sha256': '837b4156a224f3c37ddeb449319b5dc43c9964f8112abaadab079ce74a299385'},
           {'path': 'research/certificates/sharp-certificates.json',
            'sha256': '05159fea1201604ae941abe233abbf29bb98516141ed854a32b0564a3659b43b'}],
 'exactRecoveredCount': 19,
 'skippedNonExactIds': ['n10-s4-min', 'n11-s4-min', 'n11-s6-min', 'n12-s6-min', 'n12-s8-min'],
 'manualSource': 'Exact displayed draft Appendix A states plus coherent products; separate sharp '
                 'certificate targets. All exact coefficient definitions are embedded in this '
                 'standalone script.'}

def parse_exact(text):
    value = s.sympify(text, locals={'sqrt': s.sqrt, 'I': s.I})
    if value.free_symbols or value.has(s.Float):
        raise ValueError(f'Non-exact coefficient or target: {text}')
    return value

@lru_cache(maxsize=None)
def exact_purities(n, sparse_items):
    amplitudes = [s.S.Zero] * (n + 1)
    for index, expression in sparse_items:
        if not 0 <= index <= n:
            raise ValueError(f'Dicke index {index} is outside 0,...,{n}')
        amplitudes[index] = parse_exact(expression)
    norm = s.simplify(sum(s.conjugate(c) * c for c in amplitudes))
    if norm != 1:
        raise AssertionError(f'Normalization is {norm}, expected exactly 1')
    purities = []
    for size in range(n // 2 + 1):
        # Schmidt coefficient matrix for the size | (n-size) bipartition.
        matrix = s.Matrix(size + 1, n - size + 1, lambda i, j:
            amplitudes[i+j] * s.sqrt(s.Rational(
                comb(size,i) * comb(n-size,j), comb(n,i+j))))
        reduced = matrix * matrix.conjugate().T
        purities.append(s.simplify(s.trace(reduced * reduced)))
    # A normalized pure state has the same nonzero Schmidt spectrum on both sides.
    return norm, tuple(purities[min(size,n-size)] for size in range(n+1))

def verify(spec):
    n, k = spec['n'], spec['k']
    sparse = tuple(sorted((int(i), str(c)) for i,c in spec['exactCoefficients'].items()))
    norm, purities = exact_purities(n, sparse)
    # Invert 2^r Tr(rho_r^2) = sum_{j<=r} C(r,j) S_j/C(n,j).
    actual = s.simplify(comb(n,k) * sum(
        (-1)**(k-r) * comb(k,r) * 2**r * purities[r] for r in range(k+1)))
    expected = parse_exact(spec['value'])
    if s.simplify(actual - expected) != 0:
        raise AssertionError(f"{spec['id']}: S{k}={actual}; expected {expected}")
    return {'id':spec['id'], 'n':n, 'k':k, 'mode':spec['mode'],
            'normalization':str(norm), 'exactValue':str(actual),
            'expectedValue':str(expected), 'attainmentVerifiedExact':True,
            'exactCoefficients':spec['exactCoefficients'],
            'separateSharpBoundCertificate':spec.get('separateSharpBoundCertificate',False),
            'globalOptimalityProvenByThisScript':False,
            'exactReducedPurities':[str(p) for p in purities]}

def main():
    started = time.perf_counter()
    results = []
    for spec in SPECS:
        row = verify(spec)
        results.append(row)
        print(f"PASS {row['id']}: norm=1, S{row['k']}={row['exactValue']}", flush=True)
    report = {'method':'Exact symbolic reduced purities and binomial inversion from normalized sparse Dicke amplitudes',
              'scope':'Attainment only. Global optimality requires the separate bound certificate.',
              'allPassed':True, 'verifiedCount':len(results),
              'separateSharpBoundCount':sum(row['separateSharpBoundCertificate'] for row in results),
              'python':sys.version, 'sympy':s.__version__, 'platform':platform.platform(),
              'elapsedSeconds':time.perf_counter()-started,
              'inputs':INPUT_METADATA, 'results':results,
              'scriptSHA256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    dest = Path(__file__).with_name('results.json')
    dest.write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS all {len(results)} exact witnesses; results: {dest}',flush=True)

if __name__ == '__main__':
    main()
