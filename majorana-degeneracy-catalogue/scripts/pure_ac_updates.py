"""Import the certified N27 all-simple maximum from the companion catalogue.

The numerical Liège record remains numerical. Both exact attainment and the
universal upper bound are checked against their independent reports before
this one existing partition is updated. Other entries and states are untouched.
"""
from collections import Counter
from copy import deepcopy
import hashlib
import json
import sympy as s
from partition_catalogue import math_body


def attach_pure_ac_updates(data, root):
    ac_root = root.parent/'anticoherent-state-catalogue'
    companion = json.loads((ac_root/'data/catalogue.json').read_text())
    pure_proof_id = 'n27-t6-six-copy-ppt'
    proof_id = 'pure-ac-n27-maximum'
    excluded = next(e for e in companion['entries'] if (e['n'], e['t']) == (27, 6))
    assert excluded['status'] == 'impossible' and excluded['proofId'] == pure_proof_id

    certificate_root = ac_root/'research/n27-n30/n27-t6'
    certificate_bytes = (certificate_root/'certificate.json').read_bytes()
    exclusion = json.loads((certificate_root/'independent-verification.json').read_text())
    assert exclusion['passed'] and exclusion['allSourceBlocksStrictlyPositive']
    assert exclusion['arbitraryComplexPureStatesCovered']
    assert (exclusion['n'], exclusion['t'], exclusion['copies']) == (27, 6, 6)
    assert exclusion['spinBlocks'] == 81
    assert exclusion['certificateMetadataSHA256'] == hashlib.sha256(certificate_bytes).hexdigest()

    witness_root = ac_root/'research/n27-n30/radical-witnesses'
    witness_bytes = (witness_root/'states.json').read_bytes()
    witness_input = json.loads(witness_bytes)
    witness_report = json.loads((witness_root/'independent-verification.json').read_text())
    assert witness_report['passed']
    assert witness_report['statesSHA256'] == hashlib.sha256(witness_bytes).hexdigest()
    exact = next(st for st in witness_input['states'] if st['n'] == 27)
    check = next(st for st in witness_report['states'] if st['n'] == 27)
    assert check['passed'] and check['t'] == 5
    assert check['majoranaPartition'] == [1]*27
    assert any(c['t'] == 5 and c['identityVerified'] for c in check['partialTraceChecks'])

    source = next(st for st in companion['states'] if st['id'] == 'dihedral-n27-t5')
    assert source['achievedOrder'] == 5 and source['partition'] == [1]*27
    assert source['verification']['exactNormalization']
    reductions = {c['t']: c for c in source['verification']['reductions']}
    assert all(reductions[t]['identityVerified'] for t in range(1, 6))
    assert not reductions[6]['identityVerified']
    assert len(source['exactCoefficients']) == len(exact['exactCoefficients']) == 28
    assert all(s.simplify(s.sympify(a)-s.sympify(b)) == 0
               for a, b in zip(source['exactCoefficients'], exact['exactCoefficients']))
    assert len(source['majoranaPoints']) == 27
    assert all(p['multiplicity'] == 1 for p in source['majoranaPoints'])
    assert source['id'] not in {st['id'] for st in data['states']}

    base = '../anticoherent-state-catalogue/'
    state = deepcopy(source)
    for reference in state.get('references', []):
        if reference.get('url', '').startswith('./'):
            reference['url'] = base+reference['url'][2:]
    state.update(maxMultiplicity=1, distinctPoints=27, partitionKey='-'.join(['1']*27),
                 constructionProofId=proof_id,
                 coefficientLatex=[s.latex(s.sympify(c)) for c in state['exactCoefficients']])
    state['verification'].update(
        exactAC=True, exactMultiplicity=True, globalMaximumOrderProven=True,
        method='Exact symbolic Dicke partial traces and algebraic root-multiplicity checks, combined with the independently validated global sixth-order exclusion.',
        certificatePath=base+'research/n27-n30/radical-witnesses/states.json',
        acCertificatePath=base+'research/n27-n30/radical-witnesses/independent-verification.json')
    state['notes'] = ('This exact order-5 state has 27 simple Majorana points. '
                      'Exact reduced matrices and polynomial checks certify attainment. '
                      'The global sixth-order exclusion proves maximum order five. '
                      'The numerical library representative is retained separately.')
    data['states'].append(state)

    entry = next(e for e in data['entries'] if e['n'] == 27 and e['partition'] == [1]*27)
    previous = {'lower': entry['lower'], 'upper': entry['upper'], 'status': entry['status']}
    entry.update(lower=5, upper=5, value=5, status='proven', attainment='explicit',
                 proofId=proof_id, boundLatex=r't_{\max}(27,(1^{27}))=5',
                 notes='An exactly verified state with 27 simple points attains order five. A rigorous six-copy certificate excludes order six for every pure symmetric 27-qubit state, proving this partition maximum. The source numerical representative is retained separately.')
    entry['stateIds'] = list(dict.fromkeys(entry['stateIds']+[state['id']]))
    entry['proofIds'] = list(dict.fromkeys(entry.get('proofIds', [])+[proof_id]))
    entry['upperReasons'] = [{'order': 6, 'proofId': proof_id,
        'reason': 'A validated six-copy separation excludes all pure order-6 AC states at N=27, including this full partition.'}]
    entry['references'] = entry['references'] + [{'title': 'Exact order-5 witness and rigorous global sixth-order exclusion',
                                                'url': base+'research/n27-n30/index.html#n27'}]

    data['proofs'].insert(0, {
        'id': proof_id,
        'title': '27 qubits: maximum order five with 27 simple points',
        'summary': 'An exact analytical state attains fifth order. A rigorous six-copy certificate excludes sixth order for arbitrary complex pure symmetric states.',
        'formulaLatex': r't_{\max}(27,(1^{27}))=5',
        'bodyHtml': math_body(
            r'<p>The exact dihedral witness has \(\rho_5=I_6/6\). Independent algebraic root checks prove that all 27 Majorana points are simple. Its next-order defect is nonzero.</p>'
            r'<p>The upper bound is global: a six-copy rotation-averaged moment of a hypothetical sixth-order state would obey the AC contraction equations and positive first-leg partial-time-reversal constraints. Exact binary-rational dual factors, outward-rounded Arb bounds, and a validated recursive physical representation prove that all 81 source blocks of the separating operator are strictly positive. This contradicts trace one. Arbitrary complex pure states are covered; no constellation symmetry is assumed.</p>'
            '<p><a href="'+base+'research/n27-n30/radical-witnesses/states.json">Exact Dicke coefficients</a> · '
            '<a href="'+base+'research/n27-n30/radical-witnesses/independent-verification.json">Exact attainment and simple-point verification</a> · '
            '<a href="'+base+'research/n27-n30/n27-t6/proof.tex">Global upper-bound proof</a> · '
            '<a href="'+base+'research/n27-n30/n27-t6/independent-verification.json">Full rigorous verification</a></p>'),
        'references': []})
    data['entries'].sort(key=lambda e: (e['n'], -e['lower'], -(e['reportedOrder'] or 0),
                                       -e['upper'], -e['distinctPoints'], e['partition']))
    counts = Counter(e['status'] for e in data['entries'])
    data['validation'].update(
        stateCount=len(data['states'])+len(data['numericalStates']),
        explicitStateCount=len(data['states']),
        radicalStateCount=sum(st.get('exactCoefficients') is not None for st in data['states']),
        explicitWitnessPartitions=sum(bool(e['stateIds']) for e in data['entries']),
        existenceTheoremPartitions=sum(e['attainment'] == 'existence-theorem' for e in data['entries']),
        provenPartitions=counts['proven'], excludedPartitions=counts['no_ac'],
        boundedPartitions=counts['bounded'], exactPartitions=counts['proven']+counts['no_ac'])
    data['companionPureACUpdate'] = {
        'n': 27, 'partition': [1]*27, 'previous': previous, 'lower': 5, 'upper': 5,
        'stateId': state['id'], 'pureExclusionProofId': pure_proof_id,
        'numericalStatesPreserved': True,
        'exclusionCertificateMetadataSHA256': exclusion['certificateMetadataSHA256'],
        'attainmentStatesSHA256': witness_report['statesSHA256']}
    data['compiledDate'] = '2026-11-06'
    data['schemaVersion'] = '2.4'
    return data
