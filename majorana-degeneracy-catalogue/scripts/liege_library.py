"""Attach source-library observations without promoting numerical claims."""
import json
from pathlib import Path
from collections import Counter
from partition_catalogue import tex_partition, math_body
from library_bounds import rigorous_order_upper, theorem_lower
from liege_refinement import load_refinement, attach_refined_state


def attach_library(data, root):
    folder=root/'research/liege-library'
    manifest=json.loads((folder/'manifest.json').read_text())
    audit=json.loads((folder/'observations.json').read_text())
    refinements, refined_points, refinement_audit, refinement_summary = load_refinement(root)
    core_validation=dict(data['validation'])
    by_key={(e['n'],tuple(e['partition'])):e for e in data['entries']}
    for entry in data['entries']:
        entry.update(coverage='complete-core', numericalStateIds=[], reportedOrder=None)
    numerical=[]
    for obs in audit['records']:
        src=obs['sourceRecord'];n=obs['n'];part=tuple(obs['partition'])
        state_id='liege-'+obs['id']
        urls=[{'title':'Liège '+('degenerate-state' if src['sourcecategory']=='degenerate' else 'AC-state')+' library',
               'url':src['sourcePageURL']}]
        source_files=[]
        for kind,label in [('coordinate','Original Majorana coordinates'),('Dicke','Original Dicke coefficients')]:
            filename=src.get(kind+'Filename')
            if filename and (folder/'raw'/filename).exists():
                source_files.append({'label':label,'href':'./research/liege-library/raw/'+filename})
        state={'id':state_id,'n':n,'label':'Liège numerical library · '+obs['id'],
               'achievedOrder':obs['claimedOrder'],'partition':list(part),'partitionKey':'-'.join(map(str,part)),
               'maxMultiplicity':max(part),'distinctPoints':len(part),'coefficients':obs['coefficients'],
               'exactCoefficients':None,'support':list(range(n+1)), 'majoranaPoints':obs['majoranaPoints'],
               'stateLatex':rf'|\psi_{{\mathrm{{num}}}}\rangle=\sum_{{r=0}}^{{{n}}}c_r|D_{{{n}}}^{{(r)}}\rangle',
               'verification':{'exactAC':False,'exactMultiplicity':False,
                               'method':'High-precision numerical reduced-state check. Full partition inferred from well-separated repeated coordinate rows; neither is an exact existence certificate.',
                               'certificatePath':'./research/liege-library/observations.json'},
               'numerical':{'ACCheck':obs['ACCheck'],'clustering':obs['clustering'],
                            'criterion':obs['criterion'],'claimNumericallySupported':obs['claimNumericallySupported'],
                            'sourceRecordId':obs['id'],'dickeComparison':obs.get('dickeComparison')},
               'sourceFiles':source_files,'references':urls,
               'notes':'Local counterpart of a file listed on the source website. The live host was unavailable, so online byte identity is unverified. Coefficients and constellation have been checked for consistency; a numerical residual does not prove exact anticoherence.'}
        attach_refined_state(state, refinements[obs['id']], refined_points[obs['id']], refinement_audit[obs['id']], root)
        numerical.append(state)
        if (n,part) not in by_key:
            lower=theorem_lower(n,part);upper=rigorous_order_upper(n,part)
            proof_ids=['library-rigorous-bounds','partition-filtering-existence']
            if n%2==0 and part==(1,)*n and n>=8:
                lower=max(lower,3);proof_ids.append('library-exact-families')
            if n==3*part[0]+1 and part==(part[0],)+(1,)*(n-part[0]):
                lower=max(lower,2);proof_ids.append('library-exact-families')
            if part==(n//2,n//2):proof_ids=['library-rigorous-bounds','order-one']
            assert lower<=upper
            status='no_ac' if upper==0 else 'proven' if lower==upper else 'bounded'
            key='-'.join(map(str,part))
            entry={'id':f'n{n}-p{key}','n':n,'partition':list(part),'partitionKey':key,
                   'partitionLatex':tex_partition(part),'maxMultiplicity':part[0],'distinctPoints':len(part),
                   'status':status,'value':lower if lower==upper else None,'lower':lower,'upper':upper,
                   'stateIds':[],'numericalStateIds':[],'reportedOrder':None,'coverage':'source-library',
                   'proofId':'library-rigorous-bounds','proofIds':proof_ids,
                   'boundLatex':rf't_{{\max}}({n},\lambda)={lower}' if lower==upper else rf'{lower}\le t_{{\max}}({n},\lambda)\le{upper}',
                   'notes':'This size lists partitions observed in the source library, rather than every partition. The rigorous bounds follow from separate analytic arguments. The reported numerical order is displayed independently and does not raise the guaranteed lower bound.',
                   'references':urls,'attainment':'existence-theorem' if lower else 'trivial',
                   'upperReasons':[{'order':upper+1,'proofId':'library-rigorous-bounds','reason':'Established dimension and moment inequalities exclude orders above the displayed upper bound.'}]}
            by_key[n,part]=entry;data['entries'].append(entry)
        entry=by_key[n,part]
        entry['numericalStateIds'].append(state_id)
        if obs['claimNumericallySupported']:
            assert obs['claimedOrder']<=entry['upper'], (obs['id'],'numerical order exceeds rigorous upper bound')
            entry['reportedOrder']=max(entry['reportedOrder'] or 0,obs['claimedOrder'])
    data['proofs'][:0]=[
        {'id':'library-rigorous-bounds','title':'Rigorous bounds for the larger library sizes',
         'summary':'Numerical source orders are separate from bounds valid for every state in the full partition.',
         'formulaLatex':r't\le\min\!\left(\left\lfloor n/2\right\rfloor,\ 2\left\lfloor\sqrt{n+1}\right\rfloor-1\right)',
         'bodyHtml':math_body(r'<p>The Schmidt bound gives \(t\le\lfloor n/2\rfloor\). Put \(s=\lfloor t/2\rfloor\). The spin-tensor Gram argument requires \((s+1)^2\le n+1\), giving the second upper bound above. For odd \(n\ge5\), the complement obstruction strengthens the first bound to \((n-3)/2\).</p><p>The largest-part moment bounds and the full-partition second-order saturation restriction can lower these caps further. First-order existence comes from the finite-filter theorem or the balanced Dicke state. The exact families below improve some lower bounds. A source-reported numerical order is never used as an exact witness in these inequalities.</p><p><a href="./scripts/library_bounds.py">Bound calculation</a></p>'),
         'references':[]},
        {'id':'library-exact-families','title':'Exact families supporting larger-size lower bounds',
         'summary':'Some library partitions admit analytic witnesses independently of the imported numerical data.',
         'formulaLatex':r'|\psi\rangle=\sqrt{\frac{n+2}{6n}}\left(|D_n^{(0)}\rangle+|D_n^{(n)}\rangle\right)+\sqrt{\frac{2(n-1)}{3n}}|D_n^{(n/2)}\rangle',
         'bodyHtml':math_body(r'<p>For even \(n\ge8\), this state is third-order AC. Symmetry sets the first and third centered moments to zero; its variance is \(n(n+2)/12\). The three Dicke indices are separated by at least four, so all off-diagonal reduced-matrix conditions through order three vanish.</p><p>Its Majorana polynomial is a quadratic in \(z^{n/2}\). Both roots are nonzero and distinct: the discriminant is proportional to \((n-1)\binom n{n/2}-(n+2)>0\). All \(n\) Majorana points are therefore simple. This proves \(t_{\max}(n,(1^n))\ge3\).</p><p>For \(n=3g+1\), the two-term state on indices \(g,n\), with probabilities \(n/[2(n-g)]\) and \((n-2g)/[2(n-g)]\), is second-order AC and has full partition \((g,1^{n-g})\). The support gap \(n-g=2g+1\ge3\) removes the off-diagonal conditions.</p>'),
         'references':[{'file':'PRA92_052333(2015).pdf','location':'Eq. (48), even-size third-order family.'},
                       {'file':'PRA96_032304(2017).pdf','page':9,'location':'Eq. (60), second-order family.'}]}
    ]
    data['numericalStates']=numerical
    data['numbers']=sorted({e['n'] for e in data['entries']})
    data['entries'].sort(key=lambda e:(e['n'],-e['lower'],-(e['reportedOrder'] or 0),-e['upper'],-e['distinctPoints'],e['partition']))
    data['scope']='All 270 full partitions for n=2–12, plus partitions observed in the Liège source library through n=100. Rigorous bounds and numerical reported orders are separate.'
    data['schemaVersion']='2.2'
    data['library']={'sourcePages':[{'title':'Library of spin anticoherent states','url':'http://www.oq.ulg.ac.be/libspinACstates.html'},
                                    {'title':'Library of spin anticoherent states with degeneracy','url':'http://www.oq.ulg.ac.be/libspinACstatesdeg.html'}],
                     'counts':{'records':len(manifest['records']),'availableRecords':len(numerical),
                               'importedNumericalFiles':manifest['counts']['importedNumericalFiles'],
                               'missingNumericalFiles':manifest['counts']['missingNumericalFiles']},
                     'provenanceNote':'Imported from Documents/libqiss/ACstates, using the filenames linked in saved copies of the primary pages. The live website timed out; current online byte identity is unverified. Original files are preserved.',
                     'missingFiles':[r['filename'] for r in manifest['missing']],
                     'validationNote':'All 144 catalogue vectors were refined at their reported orders, using 100-digit working precision and 85-digit saved coefficients. Independent checks certify the full partitions of the saved vectors; AC remains numerical. Sixteen coordinate/Dicke pairs did not match under the tested conventions: their coordinate-derived representatives appear here, and the original Dicke vectors were refined separately as downloads.',
                     'manifestPath':'./research/liege-library/manifest.json','auditPath':'./research/liege-library/observations.json'}
    data['library']['refinement'] = {
        'primaryStates':144,'alternateStates':16,'improvedPrimaryStates':135,'unchangedPrimaryStates':9,
        'targetHSResidual':'1e-60','worstBeforeHSResidual':refinement_summary['maximumBeforeHSResidual'],
        'worstAfterHSResidual':refinement_summary['maximumAfterHSResidual'],
        'downloadPath':'./research/liege-refinement/refined-components.zip',
        'auditPath':'./research/liege-refinement/audit/verification.json',
        'methodsPath':'./research/liege-refinement/README.md',
        'exactAC':False,'savedVectorPartitionsCertified':True}
    counts=Counter(e['status'] for e in data['entries'])
    data['validation']={'core':core_validation,**core_validation,
                        'partitionCount':len(data['entries']),'totalPartitionCount':len(data['entries']),
                        'corePartitionCount':270,'libraryPartitionCount':len(data['entries'])-270,
                        'provenPartitions':counts['proven'],'excludedPartitions':counts['no_ac'],
                        'boundedPartitions':counts['bounded'],'exactPartitions':counts['proven']+counts['no_ac'],
                        'numericalStateCount':len(numerical),'explicitStateCount':len(data['states']),
                        'stateCount':len(data['states'])+len(numerical)}
    return data
