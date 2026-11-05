"""Reindex the verified state/bound record by full Majorana partition.

No numerical optimization results are promoted to existence or optimality.
"""
from collections import Counter
import html
import re


def partitions(n, largest=None):
    if not n:
        yield ()
        return
    for first in range(min(n, largest or n), 0, -1):
        for tail in partitions(n-first, first):
            yield (first,)+tail


def tex_partition(partition):
    counts = Counter(partition)
    return '(' + ','.join(str(k) if v == 1 else str(k)+'^{'+str(v)+'}'
                          for k, v in counts.items()) + ')'


def math_body(text):
    return re.sub(r'\\\((.*?)\\\)', lambda m:
                  '<span class="math" data-tex="'+html.escape(m.group(1), quote=True)+
                  '" data-display="false">'+html.escape(m.group(1))+'</span>', text)


def reindex(data):
    legacy = data['entries']
    legacy_by_pair = {(e['n'], e['t']): e for e in legacy}
    old_proofs = {p['id']: p for p in data['proofs']}
    source2014 = [{'file':'PRA90_032314(2014).pdf', 'page':4,
                   'location':'Majorana-degeneracy exclusion criterion. The finite-filter existence argument is independently supplied here.'}]
    proofs = [
        {'id':'partition-first-order-exclusion',
         'title':'Which full partitions exclude all positive AC orders?',
         'summary':'A point containing at least half the qubits rules out first order, except for the balanced two-point partition.',
         'formulaLatex':r't_{\max}(n,\lambda)=0\quad\text{if }g>n/2\text{ or }[g=n/2,\ \lambda\ne(n/2,n/2)]',
         'bodyHtml':math_body(r'<p>Rotate a point of multiplicity \(g\) so the Dicke distribution has support \(K\ge g\). First-order anticoherence requires \(\mathbb E[K]=n/2\). This is impossible if \(g>n/2\). If \(g=n/2\), every nonzero amplitude must have \(K=n/2\), so the state is a rotated balanced Dicke state, with exactly two distinct points of multiplicity \(n/2\).</p><p>All other partitions at the half-size boundary are therefore excluded, even though their largest multiplicity meets the coarse bound. Conversely, the balanced Dicke state is first-order anticoherent. A value zero means no positive AC order; the partition itself is a valid symmetric-state partition.</p>'),
         'references':source2014},
        {'id':'partition-filtering-existence',
         'title':'Every partition strictly below the half-size boundary admits first order',
         'summary':'A finite invertible collective filter produces a first-order AC representative while preserving every multiplicity.',
         'formulaLatex':r'\max\lambda<n/2\quad\Longrightarrow\quad t_{\max}(n,\lambda)\ge1',
         'bodyHtml':math_body(r'<p>Start with any normalized symmetric state having the desired full partition. Minimize \(\|A^{\otimes n}\psi\|^2\) over \(A\in\mathrm{SL}(2,\mathbb C)\). Write \(A=U\operatorname{diag}(e^{s/2},e^{-s/2})V\), with \(s\ge0\) and \(U,V\in\mathrm{SU}(2)\). The left unitary does not affect the norm.</p><p>For every \(V\), the rotated state has a strictly positive total probability in the positive-\(J_z\) subspace. Otherwise all Dicke coefficients with \(K<n/2\) would vanish, producing a pole of multiplicity at least \(\lceil n/2\rceil\), contrary to the hypothesis. This probability is continuous on compact \(\mathrm{SU}(2)\), so it has a strictly positive uniform lower bound \(\epsilon\). The filtered norm is at least \(\epsilon e^{2s m_*}\), where \(m_*>0\) is the smallest positive spin eigenvalue, and hence diverges as \(s\to\infty\).</p><p>A finite minimizer therefore exists. Varying it on the left along each Hermitian traceless direction gives \(\langle J_x\rangle=\langle J_y\rangle=\langle J_z\rangle=0\) in the normalized filtered state. The filter is invertible, so its action on the projective line preserves distinctness and every point multiplicity. This proves existence within the exact partition, without requiring displayed Dicke coefficients.</p><p>The strict inequality matters: at the half-size boundary, a minimizing sequence can change the partition in its limit. That boundary is treated by the separate exclusion argument.</p><p><a href="./research/partition-proofs.md">Detailed partition proofs</a></p>'),
         'references':source2014},
        {'id':'partition-second-order-saturation',
         'title':'Saturating the second-order bound fixes the full partition',
         'summary':'At n = 3g + 1, only one partition with largest part g can attain second order.',
         'formulaLatex':r'n=3g+1,\quad \lambda_{\max}=g,\quad t\ge2\quad\Longrightarrow\quad\lambda=(g,1^{n-g})',
         'bodyHtml':math_body(r'<p>In the pole frame, \((K-g)(n-K)\ge0\), and second-order anticoherence fixes its expectation to \(n(n-1-3g)/6\). At \(n=3g+1\), this expectation is zero. Hence the support is contained in \(\{g,n\}\). Both endpoint probabilities are nonzero because the mean is \(n/2\).</p><p>The Majorana polynomial is therefore \(z^g(a+bz^{n-g})\) with \(ab\ne0\). It has one root of multiplicity \(g\), exactly \(n-g\) simple nonzero roots, and no root at infinity. Consequently every other full partition with largest part \(g\) has maximal AC order at most one.</p><p>This distinguishes, for example, the seven-qubit patterns \((2,1^5)\) and \((2^2,1^3)\): the former attains order two; the latter has maximum one.</p>'),
         'references':[{'file':'PRA92_052333(2015).pdf','page':4,'location':'Moment equations; saturation and the full-partition consequence are independent catalogue deductions.'}]},
        {'id':'partition-tetrahedral-states',
         'title':'Two repeated-tetrahedron witnesses at twelve qubits',
         'summary':'Exact quartic factors supply order-two states in two different full partitions.',
         'formulaLatex':r'T_\pm(z)=z^4\pm2i\sqrt3\,z^2+1,\qquad P_1=T_+^3,\quad P_2=T_+^2T_-',
         'bodyHtml':math_body(r'<p>The quartics \(T_+\) and \(T_-\) are squarefree and coprime, so the two polynomials have full partitions \((3^4)\) and \((2^4,1^4)\), respectively. Converting their coefficients to the normalized Dicke basis gives exact states with \(\rho_2=I_3/3\).</p><p>The first partition has maximum order two by the third-order multiplicity bound. For the second partition the certified range is \(2\le t_{\max}\le3\). The displayed representative fails third order, but that failure alone does not exclude a different third-order state in the same partition.</p><p><a href="./research/partition-tetrahedra/PROOF.md">Derivation</a> · <a href="./research/partition-tetrahedra/formula.tex">LaTeX</a> · <a href="./research/partition-tetrahedra/verify.py">Exact verifier</a></p>'),
         'references':[{'file':'PRA92_052333(2015).pdf','page':4,'location':'Anticoherence equations; the repeated-tetrahedron constructions are independently verified here.'}]},
        {'id':'partition-order-bounds',
         'title':'How each maximum or rigorous interval is obtained',
         'summary':'Witnesses must match the entire partition. Universal bounds can still use its largest part.',
         'formulaLatex':r't_{\max}(n,\lambda)=\max\{t\ge0:\exists\,\psi\text{ with partition }\lambda,\ \rho_t=I/(t+1)\}',
         'bodyHtml':math_body(r'<p>Every integer partition of \(n\) is realized by some symmetric pure state, so order zero is always available. The first-order theorems decide exactly whether a positive order is possible. An attached verified state improves the lower bound only for its own complete multiplicity partition.</p><p>If a point multiplicity \(g=\max\lambda\) exceeds the proved maximum at order \(t\), no state with that partition can be \(t\)-AC. Global AC nonexistence proofs give additional upper bounds. Anticoherence at a higher order implies every lower order, so excluding the next order bounds all higher ones. The saturation theorem gives a further restriction that depends on the complete partition.</p><p>A single number is a proved maximum; an interval records a proved lower and upper bound. Failure of a particular witness to reach the next order does not exclude other states with the same partition. A partition is not generally a unique state, rotation orbit, or SLOCC class.</p><p><a href="./research/partition-audit.json">Independent partition audit</a> · <a href="./research/states-verification.json">Exact-state verification</a></p>'),
         'references':[]}
    ]
    proof_lookup = {p['id']:p for p in proofs+data['proofs']}
    entries = []
    for n in data['numbers']:
        for part in partitions(n):
            g = part[0]
            key = '-'.join(map(str,part))
            candidates = sorted([st for st in data['states'] if st['n']==n and tuple(st['partition'])==part],
                                key=lambda st:(-st['achievedOrder'],st['id']))
            lower = max([st['achievedOrder'] for st in candidates], default=0)
            balanced = len(part)==2 and part[0]==part[1]
            excluded = 2*g>n or (2*g==n and not balanced)
            reasons = []
            if excluded:
                assert lower==0, (n,part,'excluded but has witness')
                upper = 0
                proof_id = 'partition-first-order-exclusion'
                attainment = 'trivial'
                proof_ids = [proof_id]
                notes = 'No state with this exact partition is first-order anticoherent, hence no positive AC order is achievable. Zero is the maximal order by the order-zero convention.'
                reasons.append({'order':1,'proofId':proof_id,'reason':'The largest part is at or above the half-size boundary, and this is not the balanced two-point partition.'})
            else:
                lower = max(1,lower)
                attainment = 'explicit' if candidates and candidates[0]['achievedOrder']==lower else 'existence-theorem'
                upper = max(e['t'] for e in legacy if e['n']==n and e['status']=='proven' and g<=e['value'])
                next_entry = legacy_by_pair[n,upper+1]
                proof_id = next_entry['proofId']
                reason = ('No pure symmetric state of this qubit number reaches the next order.' if next_entry['status']=='no_ac'
                          else 'The full partition contains a point exceeding the universal multiplicity bound at the next order.')
                reasons.append({'order':upper+1,'proofId':proof_id,'reason':reason})
                saturation_part = (g,)+(1,)*(n-g)
                if n==3*g+1 and part!=saturation_part:
                    upper = min(upper,1)
                    proof_id = 'partition-second-order-saturation'
                    reasons.insert(0,{'order':2,'proofId':proof_id,'reason':'Equality in the second-order moment bound permits only the partition (g,1^(n-g)).'})
                proof_ids = [proof_id]
                if 2*g<n:
                    proof_ids.append('partition-filtering-existence')
                elif balanced:
                    proof_ids.append('order-one')
                if lower>1 and candidates:
                    if candidates[0]['id']=='gmax-n11-t3-certified':proof_ids.append('n11-algebraic')
                    if candidates[0]['id'] in ('partition-n12-3x4-tetra','partition-n12-2x4-1x4-tetra'):proof_ids.append('partition-tetrahedral-states')
                    if candidates[0].get('constructionProofId'):proof_ids.append(candidates[0]['constructionProofId'])
                notes = ('A verified state with this exact full partition attains the displayed lower bound.' if attainment=='explicit' else
                         'First-order existence in this exact partition follows from a finite invertible collective filter. An explicit filtered coefficient vector is not attached.')
                notes += (' The matching upper bound proves maximality.' if lower==upper else
                          ' The highest achievable order within this full partition is not yet determined; the interval is rigorous.')
            assert 0<=lower<=upper, (n,part,lower,upper)
            status = 'no_ac' if upper==0 else 'proven' if lower==upper else 'bounded'
            value = lower if lower==upper else None
            if n==12 and part==(2,2,2,2,1,1,1,1):proof_ids.append('n12-four-double-obstruction')
            proof_ids = list(dict.fromkeys(proof_ids+['partition-order-bounds']))
            refs=[]
            for pid in proof_ids:
                for ref in proof_lookup[pid].get('references',[]):
                    if ref not in refs:refs.append(ref)
            bound = rf't_{{\max}}({n},\lambda)={lower}' if value is not None else rf'{lower}\le t_{{\max}}({n},\lambda)\le {upper}'
            entries.append({'id':f'n{n}-p{key}','n':n,'partition':list(part),'partitionKey':key,
                            'partitionLatex':tex_partition(part),'maxMultiplicity':g,'distinctPoints':len(part),
                            'status':status,'value':value,'lower':lower,'upper':upper,
                            'stateIds':[st['id'] for st in candidates], 'proofId':proof_id,'proofIds':proof_ids,
                            'boundLatex':bound,'notes':notes,'references':refs,'attainment':attainment,'upperReasons':reasons})
    entries.sort(key=lambda e:(e['n'],-e['lower'],-e['upper'],-e['distinctPoints'],e['partition']))
    counts = Counter(e['status'] for e in entries)
    assert len(entries)==270 and len({e['id'] for e in entries})==270
    data.update(schemaVersion='2.0', title='Anticoherence by full Majorana multiplicity partition',
                scope='Pure symmetric states, n=2–12. All 270 full multiplicity partitions; highest achievable AC order, or rigorous lower and upper bounds.',
                legacyEntries=legacy, legacyValidation=data['validation'], entries=entries, proofs=proofs+data['proofs'],
                validation={'partitionCount':len(entries),'provenPartitions':counts['proven'],
                            'excludedPartitions':counts['no_ac'],'boundedPartitions':counts['bounded'],
                            'exactPartitions':counts['proven']+counts['no_ac'], 'stateCount':len(data['states']),
                            'radicalStateCount':sum(st['exactCoefficients'] is not None for st in data['states']),
                            'certifiedImplicitStateCount':sum(st['exactCoefficients'] is None for st in data['states']),
                            'explicitWitnessPartitions':sum(bool(e['stateIds']) for e in entries),
                            'existenceTheoremPartitions':sum(e['attainment']=='existence-theorem' for e in entries)})
    data['conventions']['zeroOrder']='Every partition is realized by a pure symmetric state. A displayed maximum zero means no state in that partition has any positive AC order.'
    data['conventions']['partition']='The entire unordered multiplicity partition is fixed, including all simple points. Geometry remains free. Witnesses only support their exact partition; a matching largest part is insufficient.'
    return data
