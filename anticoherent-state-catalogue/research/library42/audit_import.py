"""Independently audit the selected refined Liège vectors through the next order.

The source library and the degeneracy catalogue are read only. Existing checks
are reused only when the complete refinement input has the same SHA-256 hash.
All newly imported vectors are checked directly at 100 decimal digits.
"""
from pathlib import Path
from functools import lru_cache
import hashlib
import json
import sys
import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REFINEMENT = ROOT.parent / 'majorana-degeneracy-catalogue/research/liege-refinement/refined.json'
DEGENERACY = ROOT.parent / 'majorana-degeneracy-catalogue/data/catalogue.json'
mp.mp.dps = 100
THRESHOLD = mp.mpf('1e-60')


@lru_cache(None)
def factors(n, q):
    return [[mp.sqrt(mp.binomial(q, a) * mp.binomial(n-q, b) /
                     mp.binomial(n, a+b)) for b in range(n-q+1)]
            for a in range(q+1)]


def audit(record):
    n, t = record['n'], record['reportedOrder']
    values = [mp.mpc(a, b) for a, b in record['refinedCoefficients']]
    norm_squared = sum(abs(x)**2 for x in values)
    values = [x / mp.sqrt(norm_squared) for x in values]
    reductions = []
    for q in range(1, t+2):
        f = factors(n, q)
        schmidt = mp.matrix([[values[a+b] * f[a][b]
                              for b in range(n-q+1)] for a in range(q+1)])
        rho = schmidt * schmidt.H
        delta = rho - mp.eye(q+1)/(q+1)
        distance_squared = sum(abs(x)**2 for x in delta)
        distance = mp.sqrt(distance_squared)
        supported = distance < THRESHOLD
        assert supported == (q <= t), (record['id'], q, mp.nstr(distance, 30))
        reductions.append({
            't': q, 'maxEntryError': mp.nstr(max(abs(x) for x in delta), 60),
            'frobeniusDistance': mp.nstr(distance, 60),
            'frobeniusDistanceSquared': mp.nstr(distance_squared, 60),
            'purity': mp.nstr(sum(abs(x)**2 for x in rho), 60),
            'numericallyMaximallyMixed': supported,
        })
    return {
        'id': record['id'], 'n': n, 'reportedOrder': t,
        'workingPrecisionDigits': 100,
        'serializedCoefficientDigits': record.get('serializedDigits', 85),
        'inputNormSquared': mp.nstr(norm_squared, 100),
        'inputNormalizedBeforeCheck': True, 'reductions': reductions,
        'reportedOrderNumericallySupported': True,
        'nextOrderFailsNumerically': True,
    }


def main(max_n):
    digest = hashlib.sha256(REFINEMENT.read_bytes()).hexdigest()
    previous = json.loads((HERE / 'verification.json').read_text())
    cached = ({r['id']: r for r in previous['records']}
              if previous.get('sourceSHA256') == digest else {})
    records = json.loads(REFINEMENT.read_text())['records']
    if isinstance(records, list):
        records = {r['id']: r for r in records}
    selected = [s for s in json.loads(DEGENERACY.read_text())['numericalStates']
                if s['n'] <= max_n]
    checks = []
    new_count = 0
    for state in selected:
        source_id = state['numerical']['sourceRecordId']
        record = records[source_id]
        assert record['reportedOrder'] == state['achievedOrder']
        if source_id in cached:
            check = cached[source_id]
            assert (check['n'], check['reportedOrder']) == (state['n'], state['achievedOrder'])
        else:
            check = audit(record)
            new_count += 1
            print(source_id, 'independent 100-digit check passed', flush=True)
        checks.append(check)
    selected_ids = {r['id'] for r in checks}
    extras = previous.get('additionalSourceRecords', [])
    for extra in extras:
        if extra['id'] in cached and cached[extra['id']]['n'] <= max_n:
            assert extra['id'] not in selected_ids
            checks.append(cached[extra['id']])
    through = [mp.mpf(r['frobeniusDistance']) for c in checks for r in c['reductions']
               if r['t'] <= c['reportedOrder']]
    following = [mp.mpf(r['frobeniusDistance']) for c in checks for r in c['reductions']
                 if r['t'] == c['reportedOrder']+1]
    report = {
        'schemaVersion': 1,
        'method': 'Independent normalized Schmidt-matrix reconstruction and reduced-state multiplication at 100 decimal digits; no use of stored residuals.',
        'scientificStatus': 'Numerical evidence only, not exact existence or optimality proof.',
        'sourceFile': '../majorana-degeneracy-catalogue/research/liege-refinement/refined.json',
        'sourceSHA256': digest, 'workingPrecisionDigits': 100,
        'acceptanceFrobeniusDistance': '1e-60', 'records': checks,
        'summary': {
            'records': len(checks), 'selectedRefinedRecords': len(selected),
            'newlyCheckedRecords': new_count, 'maximumQubits': max_n,
            'allReportedOrdersSupported': True, 'allNextOrdersFailNumerically': True,
            'maximumReportedOrderFrobeniusDistance': mp.nstr(max(through), 60),
            'minimumNextOrderFrobeniusDistance': mp.nstr(min(following), 60),
        },
        'additionalSourceRecords': extras,
    }
    temporary = HERE / 'verification.temporary.json'
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(HERE / 'verification.json')
    print('PASS:', len(checks), 'numerical representatives;', new_count, 'new checks.', flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else
         json.loads((HERE / 'config.json').read_text())['maxQubits'])
