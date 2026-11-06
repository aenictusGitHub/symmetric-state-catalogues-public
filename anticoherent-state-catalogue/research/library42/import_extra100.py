"""Audit additional local Liège records omitted from the selected refinement.

Only nearby refinements at the filename's reported order are attempted. Every
accepted record remains numerical and must pass independent next-order checks.
The original files and the companion catalogue are read only.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import numpy as np
import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEGENERACY = ROOT.parent / 'majorana-degeneracy-catalogue'
sys.path.insert(0, str(DEGENERACY / 'research/liege-library'))
from analyze_records import real_tokens, read_complex, compare_dicke
from analyze_core import reconstruct, coefficient_strings
from audit_import import audit


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


refiner = module('extra100_refiner', DEGENERACY / 'research/liege-refinement/refine.py')
point_refiner = module('extra100_points', DEGENERACY / 'research/liege-refinement/refine_points.py')
geometry = module('extra100_geometry', ROOT / 'scripts/build_majorana.py')


def main(source):
    filenames = ['acs.65.10.txt', 'acs.70.11.txt', 'acs.72.11.txt',
                 'acs.75.11.txt', 'acs.84.11.txt', 'acs.96.12.txt']
    existing = json.loads((HERE / 'states.json').read_text())['states']
    states, records, checked = [], [], []
    for filename in filenames:
        _, ns, ts, _ = filename.split('.')
        n, t = int(ns), int(ts)
        source_id = filename[:-4]
        path = source / filename
        observation = {'id': source_id, 'n': n, 'reportedOrder': t,
                       'sourceFilename': filename,
                       'sourceSHA256': hashlib.sha256(path.read_bytes()).hexdigest()}
        try:
            with mp.workdps(100):
                values, _ = real_tokens(path)
                assert len(values) == 3*n, ('coordinate count', len(values), n)
                raw_points = [tuple(values[i:i+3]) for i in range(0, len(values), 3)]
                coefficients, points, _ = reconstruct(raw_points)
                selection = '100-digit spinor reconstruction of the original coordinates'
                dicke = source / ('acsdk.' + ns + '.' + ts + '.txt')
                if dicke.exists():
                    raw = read_complex(dicke)
                    assert len(raw) == n+1
                    norm, best = compare_dicke(raw, coefficients)
                    observation['DickeComparisonInfidelity'] = mp.nstr(best[0], 50)
                    observation['DickeConvention'] = best[2]
                    if best[0] < mp.mpf('1e-20'):
                        coefficients = best[1]
                        selection = 'Normalized original Dicke file, convention-matched to the coordinate state'
                poles = []
                zeros = []
                for sign, name in [(-1, 'south'), (1, 'north')]:
                    count = sum(x == 0 and y == 0 and z == sign for x, y, z in points)
                    if count:
                        indices = list(range(count)) if sign < 0 else list(range(n-count+1, n+1))
                        assert all(abs(coefficients[i]) < mp.mpf('1e-20') for i in indices)
                        for i in indices:
                            coefficients[i] = mp.mpc(0)
                        zeros.extend(indices)
                        poles.append({'pole': name, 'multiplicity': count, 'indices': indices})
                partition = sorted([p['multiplicity'] for p in poles] +
                                   [1] * (n-sum(p['multiplicity'] for p in poles)), reverse=True)
                input_record = {
                    'id': source_id, 'n': n, 'reportedOrder': t, 'partition': partition,
                    'originalCoefficients': coefficient_strings(coefficients, 95),
                    'originalPoints': [[mp.nstr(v, 95) for v in p] for p in points],
                    'structuralZeroIndices': zeros, 'repeatedPoles': poles,
                    'sourceCoordinateFilename': filename, 'sourceCoefficientSelection': selection,
                }
                refined = refiner.refine(input_record)
                observation['refinement'] = refined
                if refined['status'] != 'target reached':
                    observation['status'] = 'not accepted: reported-order refinement did not reach the numerical threshold'
                    records.append(observation)
                    continue
                check = audit(refined)
                browser = [[float(mp.mpf(a)), float(mp.mpf(b))] for a, b in refined['refinedCoefficients']]
                browser_vector = [complex(a, b) for a, b in browser]
                duplicate = next((s['id'] for s in existing if s['n'] == n and
                                  s.get('sourceRecordId') != source_id and
                                  geometry.ray_error([complex(a, b) for a, b in s['coefficients']], browser_vector) < 1e-10), None)
                if duplicate:
                    observation.update(status='not imported: duplicate numerical ray', duplicateStateId=duplicate)
                    records.append(observation)
                    continue
            with mp.workdps(110):
                plotted = point_refiner.compute(refined)
            assert not plotted['numericalRootErrors']
            points = plotted['refinedMajoranaPoints']
            reconstruction_error = geometry.ray_error(browser_vector, geometry.reconstruct(points, n))
            assert reconstruction_error < 1e-10
            reductions = [dict(r, numericallyConsistent=r['t'] <= t, identityVerified=None)
                          for r in check['reductions']]
            maximum = max((r['frobeniusDistance'] for r in reductions if r['t'] <= t), key=mp.mpf)
            state = {
                'id': 'liege-' + source_id, 'n': n, 'achievedOrder': t,
                'witnessesOrders': list(range(1, t+1)),
                'label': f'Liège numerical state · reported order {t} · {source_id}',
                'status': 'numerical',
                'stateLatex': rf'|\psi_{{\mathrm{{num}}}}\rangle\simeq\sum_{{r=0}}^{{{n}}}c_r|D_{{{n}}}^{{(r)}}\rangle',
                'support': [i for i, c in enumerate(browser) if c != [0., 0.]],
                'coefficients': browser, 'highPrecisionCoefficients': refined['refinedCoefficients'],
                'majoranaPoints': points, 'partition': partition,
                'normalization': 'Numerically normalized at 100-digit working precision; exact normalization is not claimed.',
                'coefficientConvention': 'Normalized Dicke basis; r counts excitations. Downloads retain 85-digit numerical components.',
                'references': [{'title': 'Liège AC-state library', 'url': 'http://www.oq.ulg.ac.be/libspinACstates.html'},
                               {'title': 'Independent additional-source audit', 'url': './research/library42/extra-source-audit-100.json'}],
                'provenance': 'Additional original Liège-library record, independently refined and checked at its reported order.',
                'notes': 'The reported order and failure at the next order are independently checked at 100-digit precision. These numerical checks are not exact existence or maximum-order proofs.',
                'verification': {'exactAC': False, 'exactNormalization': False, 'globalMaximumOrderProven': False,
                                 'method': 'Independent Schmidt-matrix partial traces from the saved 85-digit components.',
                                 'workingPrecisionDigits': 100, 'serializedDigits': 85,
                                 'maximumResidual': maximum, 'acceptanceThreshold': '1e-60',
                                 'reductions': reductions},
                'sourceRecordId': source_id, 'sourceRefinementResidual': refined['serializedHSResidual'],
            }
            states.append(state)
            checked.append(check)
            observation.update(status='accepted numerical representative', verification=check,
                               geometry=plotted, projectiveReconstructionError=reconstruction_error)
            records.append(observation)
            print(source_id, 'accepted numerical order', t, flush=True)
        except Exception as error:
            observation.update(status='not imported', error=repr(error))
            records.append(observation)
            print(source_id, 'not imported:', repr(error), flush=True)
    audit_report = {'scientificStatus': 'Numerical source checks only; no exact AC existence or optimality claim.',
                    'records': records, 'acceptedRecords': len(states), 'originalFilesUnchanged': True}
    (HERE / 'extra-states-100.json').write_text(json.dumps({'states': states}, indent=2) + '\n')
    (HERE / 'extra-source-audit-100.json').write_text(json.dumps(audit_report, indent=2) + '\n')
    verification = json.loads((HERE / 'verification.json').read_text())
    by_id = {r['id']: r for r in verification['records']}
    by_id.update({r['id']: r for r in checked})
    verification['records'] = list(by_id.values())
    extras = {r['id']: r for r in verification.get('additionalSourceRecords', [])}
    for r in records:
        if r['status'] == 'accepted numerical representative':
            extras[r['id']] = {k: r[k] for k in ['id', 'sourceFilename', 'sourceSHA256']}
            extras[r['id']]['audit'] = 'extra-states-100.json'
    verification['additionalSourceRecords'] = list(extras.values())
    verification['summary']['records'] = len(by_id)
    verification['summary']['additionalRecordsThrough100'] = len(states)
    with mp.workdps(100):
        verification['summary']['maximumReportedOrderFrobeniusDistance'] = mp.nstr(max(mp.mpf(r['frobeniusDistance']) for c in by_id.values() for r in c['reductions'] if r['t'] <= c['reportedOrder']), 60)
        verification['summary']['minimumNextOrderFrobeniusDistance'] = mp.nstr(min(mp.mpf(r['frobeniusDistance']) for c in by_id.values() for r in c['reductions'] if r['t'] == c['reportedOrder']+1), 60)
    (HERE / 'verification.json').write_text(json.dumps(verification, indent=2) + '\n')
    print('Additional source audit complete:', len(states), 'accepted of', len(filenames), flush=True)


if __name__ == '__main__':
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path.home() / 'Documents/libqiss/ACstates')
