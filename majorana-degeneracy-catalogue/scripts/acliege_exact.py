"""Attach independently verified exact Julia witnesses after the library import.

Numerical records and reported orders are left intact. An exact witness can
raise only the lower bound for its certified full Majorana partition.
"""
import json
from collections import Counter
from copy import deepcopy
import sympy as s
from partition_catalogue import math_body


def attach_acliege_exact(data, root):
    audit = json.loads((root/'research/acliege/exact-states.json').read_text())
    proof_id = 'acliege-exact-higher-orders'
    entries = {(row['n'],tuple(row['partition'])):row for row in data['entries']}
    known_ids = {state['id'] for state in data['states']}
    improvements = []
    for source in audit['states']:
        state = deepcopy(source)
        n, order = state['n'], state['achievedOrder']
        assert n in (13,15,16,18,20) and state['id'] not in known_ids
        assert state['partition'] == [1]*n
        assert state['geometryVerification']['squareFreeFactorizationExact']
        assert state['verification']['exactNormalization']
        assert len(state['verification']['reductions']) == order
        assert all(row['identityVerified'] for row in state['verification']['reductions'])
        state.update(maxMultiplicity=1, distinctPoints=n,
                     partitionKey='-'.join(['1']*n), constructionProofId=proof_id,
                     coefficientLatex=[s.latex(s.sympify(c)) for c in state['exactCoefficients']])
        state['verification'].update(
            exactAC=True, exactMultiplicity=True,
            method='Exact symbolic reduced matrices and square-free factorization over the algebraic coefficient field, including zero and infinity.',
            certificatePath='./research/acliege/exact-states.json',
            acCertificatePath='./research/acliege/verification.json')
        state['notes'] = ('Analytic formula from the supplied ACLiege.jl source, credited there to M. Grassl. '
                          'Every reduced matrix through the stated order and the full simple-point partition were independently verified exactly. '
                          'The next-order failure concerns this representative only. It is not a maximality proof for the partition. '
                          'The previously imported numerical representative is retained separately.')
        state['references'] = [{'file':'ACLiege.jl',
            'location':f"acstateslist, line {state['source']['line']}; exact verification in research/acliege/"}]
        data['states'].append(state)
        known_ids.add(state['id'])
        entry = entries[n,tuple(state['partition'])]
        old_lower = entry['lower']
        entry['stateIds'].append(state['id'])
        entry['lower'] = max(old_lower,order)
        assert entry['lower'] <= entry['upper']
        entry['status'] = 'proven' if entry['lower']==entry['upper'] else 'bounded'
        entry['value'] = entry['lower'] if entry['lower']==entry['upper'] else None
        entry['attainment'] = 'explicit'
        entry['proofIds'] = list(dict.fromkeys(entry['proofIds']+[proof_id]))
        entry['boundLatex'] = (rf"t_{{\max}}({n},\lambda)={entry['lower']}" if entry['value'] is not None else
                               rf"{entry['lower']}\le t_{{\max}}({n},\lambda)\le {entry['upper']}")
        entry['notes'] = ('An exactly verified analytic representative with this full partition attains the displayed lower bound. '
                          'The upper bound remains the independently established general bound. '
                          'The numerical library representative and its reported order remain separate. '
                          'Only source-observed partitions are listed at this qubit number; the highest achievable order is not yet settled.')
        entry['references'] = entry['references'] + [ref for ref in state['references'] if ref not in entry['references']]
        improvements.append({'n':n,'partition':state['partition'],'previousLower':old_lower,
                             'lower':entry['lower'],'upper':entry['upper'],'stateId':state['id']})
    data['proofs'].insert(0,{
        'id':proof_id,'title':'Five exact larger-size witnesses from the supplied Julia formulas',
        'summary':'All-simple constellations attain order three at 13 and 15 qubits, and order five at 16, 18 and 20 qubits. These improve rigorous lower bounds without asserting maximality.',
        'formulaLatex':r'\rho_t=I_{t+1}/(t+1),\qquad\lambda=(1^n)',
        'bodyHtml':math_body(r'<p>The supplied <code>ACLiege.jl</code> file contains analytic Dicke coefficients, credited in its source comment to M. Grassl. The formulas were read without executing Julia. Exact symbolic partial traces verify every order through the stated order. Exact square-free factorization of the Majorana polynomial, together with its zero and infinity factors, proves that all points are simple.</p>'
            r'<p>For \(n=13,15\), the rigorous lower bound rises from one to three. For \(n=16,18,20\), it rises from three to five. The resulting ranges are respectively \([3,5]\), \([3,6]\), and \([5,7]\) at each of the three even sizes. The next-order defect of a particular witness does not exclude a different higher-order state. Previously imported numerical states are retained and are not promoted to exact states.</p>'
            r'<p><a href="./research/acliege/exact-states.json">Exact coefficients and polynomial factors</a> · <a href="./research/acliege/verification.json">Exact verification</a> · <a href="./research/acliege/verify_exact_states.py">Reproduce the checks</a> · <a href="../ac-subspace-catalogue/index.html">Anticoherent-subspace catalogue</a></p>'),
        'references':[{'file':'ACLiege.jl','location':'acstateslist, lines 137, 157–179, 203 and 229.'}]})
    data['sources'].append({'filename':'ACLiege.jl','title':'Supplied analytic and numerical spin-state library',
        'authors':'Local source supplied by the user; higher-size analytic formulas credited in the source to M. Grassl.',
        'notes':'Only the five independently verified analytic higher-size witnesses described above enter the rigorous bounds. Decimal-only formulas and source optimality labels are not exact certificates.'})
    data['entries'].sort(key=lambda e:(e['n'],-e['lower'],-(e['reportedOrder'] or 0),-e['upper'],-e['distinctPoints'],e['partition']))
    counts=Counter(e['status'] for e in data['entries'])
    data['validation'].update(
        stateCount=len(data['states'])+len(data['numericalStates']), explicitStateCount=len(data['states']),
        radicalStateCount=sum(st.get('exactCoefficients') is not None for st in data['states']),
        explicitWitnessPartitions=sum(bool(e['stateIds']) for e in data['entries']),
        existenceTheoremPartitions=sum(e['attainment']=='existence-theorem' for e in data['entries']),
        provenPartitions=counts['proven'],excludedPartitions=counts['no_ac'],boundedPartitions=counts['bounded'],
        exactPartitions=counts['proven']+counts['no_ac'])
    data['exactHigherSizeUpdate']={'source':'references/ACLiege.jl','improvements':improvements,
                                  'numericalStatesPreserved':True,'upperBoundsUnchanged':True}
    data['compiledDate']='2026-10-29'
    data['schemaVersion']='2.3'
    return data
