"""Attach audited numerical refinements without changing rigorous AC bounds."""
import hashlib
import json
from decimal import Decimal


def load_refinement(root):
    folder = root / 'research/liege-refinement'
    refined_path = folder / 'refined.json'
    refined = json.loads(refined_path.read_text())
    points = json.loads((folder / 'refined-points.json').read_text())
    audit = json.loads((folder / 'audit/verification.json').read_text())
    alternate_audit = json.loads((folder / 'audit/alternate-verification.json').read_text())
    points_audit = json.loads((folder / 'audit/points-verification.json').read_text())
    assert hashlib.sha256(refined_path.read_bytes()).hexdigest() == audit['inputs']['refinedSHA256']
    source = root / 'research/liege-library/observations.json'
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    if source_hash != audit['inputs']['originalSHA256']:
        equivalence = audit['inputs']['packagedInputEquivalence']
        assert equivalence['allConsumedScientificFieldsIdentical'] and equivalence['records'] == 144
        assert source_hash == audit['inputs']['equivalentPackagedOriginalSHA256'] == equivalence['packagedSourceSHA256']
    summary = audit['summary']
    assert summary['allSelectedPassed'] and summary['completeInventory'] and summary['allRootsChecked']
    assert summary['passedRecords'] == 144 and not audit['failures']
    assert alternate_audit['summary']['passedRecords'] == 16
    assert points_audit['summary']['passedRecords'] == 144
    for filename, report, key in [
            ('alternate-inputs.json', alternate_audit, 'originalSHA256'),
            ('alternate-refined.json', alternate_audit, 'refinedSHA256'),
            ('refined.json', points_audit, 'refinedSHA256'),
            ('refined-points.json', points_audit, 'pointsSHA256')]:
        assert hashlib.sha256((folder / filename).read_bytes()).hexdigest() == report['inputs'][key]
    records = {r['id']: r for r in refined['records']}
    point_records = {r['id']: r for r in points['records']}
    audit_records = {r['id']: r for r in audit['records']}
    assert len(records) == 144 and records.keys() == point_records.keys() == audit_records.keys()
    return records, point_records, audit_records, summary


def attach_refined_state(state, record, points, audit, root):
    assert record['n'] == state['n'] and record['reportedOrder'] == state['achievedOrder']
    assert record['partition'] == state['partition'] == audit['refinedPartition']
    assert audit['passed'] and audit['finiteRootVerification']['passed']
    assert Decimal(audit['after']['HSResidual']) < Decimal('1e-60')
    relative = 'research/liege-refinement/coefficients/' + record['id'] + '.refined.txt'
    # Downloaded decimals, rather than browser floats, are the certified vector.
    exported = [line.split() for line in (root / relative).read_text().splitlines()
                if line.strip() and not line.startswith('#')]
    assert exported == record['refinedCoefficients']
    state['coefficients'] = [[float(a), float(b)] for a, b in record['refinedCoefficients']]
    state['majoranaPoints'] = points['refinedMajoranaPoints']
    state['label'] = 'Liège refined numerical state · ' + record['id']
    state['verification'].update(
        exactAC=False, exactMultiplicity=True,
        method='AC residual checked numerically from the 85-digit saved coefficients. '
               'Exact pole zeros and interval root counts certify their full multiplicity partition. '
               'Browser values and sphere points are rounded displays of that saved vector.')
    state['refinement'] = {
        'beforeResidual': audit['before']['HSResidual'],
        'afterResidual': audit['after']['HSResidual'],
        'target': '1e-60',
        'metric': 'Hilbert–Schmidt residual of the reduced state at the reported AC order.',
        'workingPrecisionDigits': record['precisionDigits'],
        'serializedDigits': record['serializedDigits'],
        'coefficientPath': './' + relative,
        'partitionCertificatePath': './research/liege-refinement/audit/verification.json',
        'partitionPreservation': 'The full partition of the saved coefficients is certified by exact '
                                 'pole zeros and disjoint interval root counts. AC remains numerical.',
        'maxComponentChange': record['maxComponentChange'],
        'maxPointDisplacement': points['maxPointDisplacement'],
        'singularFamilyCorrection': record.get('singularFamilyCorrection'),
    }
    state['notes'] = ('Refined at the source-reported AC order; the original source files and checks are retained. '
                      'The downloaded 85-digit coefficients and independently checked sphere points describe the '
                      'refined numerical state. The full partition of that saved vector is certified; '
                      'anticoherence and optimality are not thereby proved. The live source host was unavailable, '
                      'so current online byte identity is unverified.')

