"""Import independently checked refined library records through N=42.

Uses the saved, audited refinement rather than mismatched original Dicke files.
The original library remains read-only. Exact lower-order representatives are
constructed independently where separated positive moment quadrature succeeds.
"""
from pathlib import Path
import importlib.util,json,hashlib
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SIBLING=ROOT.parent/'majorana-degeneracy-catalogue'

def main():
    report=json.loads((HERE/'verification.json').read_text())
    checks={x['id']:x for x in report['records']}
    deg=json.loads((SIBLING/'data/catalogue.json').read_text())
    refined=json.loads((SIBLING/'research/liege-refinement/refined.json').read_text())['records']
    if isinstance(refined,list):refined={r['id']:r for r in refined}
    states=[]
    for source in deg['numericalStates']:
        if source['n']>42:continue
        sid=source['numerical']['sourceRecordId'];r=refined[sid];audit=checks[sid]
        assert audit['reportedOrder']==source['achievedOrder']
        reductions=[]
        for c in audit['reductions']:
            q=c['t'];ok=q<=source['achievedOrder']
            assert (S.Float(c['frobeniusDistance'],100)<S.Float('1e-60',100))==ok
            reductions.append(dict(c,numericallyConsistent=ok,identityVerified=None))
        coeffs=r['refinedCoefficients']
        state={'id':source['id'],'n':source['n'],'achievedOrder':source['achievedOrder'],'witnessesOrders':list(range(1,source['achievedOrder']+1)),
            'label':f"Liège numerical state · reported order {source['achievedOrder']} · {sid}",'status':'numerical','stateLatex':source['stateLatex'],'support':source['support'],
            'coefficients':source['coefficients'],'highPrecisionCoefficients':coeffs,'majoranaPoints':source['majoranaPoints'],'partition':source['partition'],
            'normalization':'Numerically normalized using 100-digit working precision; exact normalization is not claimed.',
            'coefficientConvention':'Normalized Dicke basis; r counts excitations. Browser values are rounded; downloads retain the saved 85-digit components.',
            'references':[{'title':'Liège AC-state library','url':'http://www.oq.ulg.ac.be/libspinACstates.html'},{'title':'Independent high-precision checks and import notes','url':'./research/library42/index.html'}],
            'provenance':'Liège AC-state library · independently checked numerical representative.',
            'notes':'Numerical consistency at the reported order is not an exact existence proof or proof of maximum order. The next reduced state fails the same numerical test. Saved components are retained at 85 digits; displayed amplitudes and sphere points are rounded.',
            'verification':{'exactAC':False,'exactNormalization':False,'globalMaximumOrderProven':False,'method':'Independent 100-digit Schmidt-matrix partial traces from the saved 85-digit components, through the next AC order.','workingPrecisionDigits':100,'serializedDigits':85,'maximumResidual':max((c['frobeniusDistance'] for c in reductions if c['t']<=source['achievedOrder']),key=lambda x:S.Float(x,100)),'acceptanceThreshold':'1e-60','reductions':reductions},
            'sourceRecordId':sid,'sourceRefinementResidual':r['serializedHSResidual']}
        states.append(state)
    extra=HERE/'extra-state.json'
    if extra.exists():states.append(json.loads(extra.read_text())['state'])
    (HERE/'states.json').write_text(json.dumps({'states':states,'sourceLibrary':'Documents/libqiss/ACstates','recordCount':len(states)},indent=2,ensure_ascii=False)+'\n')
    # Reuse the established exact construction and verify each next-order failure.
    spec=importlib.util.spec_from_file_location('exact_orders',ROOT/'research/exact-orders/build.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    old=[]
    for name in ['research/states.json','research/cross-catalogue/exact-states.json','research/exact-orders/states.json']:
        old.extend(json.loads((ROOT/name).read_text())['states'])
    pairs={(s['n'],s['achievedOrder']) for s in old}
    output=[];unfilled=[]
    for n in range(2,43):
        top=max([s['achievedOrder'] for s in states+old if s['n']==n],default=0)
        for t in range(1,top):
            if (n,t) in pairs:continue
            try:
                if t==1:c=[1/S.sqrt(2) if r in (0,n) else S.S.Zero for r in range(n+1)];method='GHZ state'
                else:c,method=mod.sparse_weights(n,t)
            except ValueError:
                unfilled.append([n,t]);continue
            support=[i for i,x in enumerate(c) if x!=0];reductions=[];next_witness=None
            for q in range(1,t+2):
                delta=mod.reduction(c,n,q)-S.eye(q+1)/(q+1)
                zero=all(S.simplify(x)==0 for x in delta)
                assert zero==(q<=t),(n,t,q)
                dist=S.simplify(sum(S.expand_complex(x*S.conjugate(x)) for x in delta))
                reductions.append({'t':q,'identityVerified':zero,'frobeniusDistanceSquared':str(dist),'purity':str(S.simplify(dist+S.Rational(1,q+1)))})
                if q==t+1:
                    a,b=next((a,b) for a in range(q+1) for b in range(q+1) if delta[a,b]!=0)
                    next_witness={'order':q,'row':a,'column':b,'differenceFromIdentity':str(delta[a,b]),'differenceLatex':S.latex(delta[a,b])}
            formula=' + '.join(S.latex(c[r])+rf'|D_{{{n}}}^{{({r})}}\rangle' for r in support).replace('+ -','- ')
            output.append({'id':f'exact-order-n{n}-t{t}','n':n,'achievedOrder':t,'witnessesOrders':list(range(1,t+1)),'label':f'Exactly {t}-AC · {len(support)} Dicke components','status':'proven-existence','stateLatex':formula,'support':support,'exactCoefficients':list(map(str,c)),'coefficients':[[float(S.re(x)),float(S.im(x))] for x in c],'normalization':'sum_r |c_r|^2 = 1, verified exactly','coefficientConvention':'Normalized Dicke basis; r counts excitations.','references':[{'title':'Exact-order construction and checks through N=42','url':'./research/library42/index.html'},{'title':'General analytical prescription','url':'./research/exact-orders/general/index.html'}],'provenance':'Independent exact-order construction: '+method,'notes':f'This state is {t}-AC but not {t+1}-AC. All reduced-state identities and a nonzero next-order entry are checked exactly.','verification':{'method':'Exact symbolic normalized Dicke-basis partial traces through the next order.','exactNormalization':True,'exactOrderVerified':True,'globalMaximumOrderProven':False,'reductions':reductions,'nextOrderWitness':next_witness}})
            print('New exact',n,t,flush=True)
    (HERE/'exact-lower-states.json').write_text(json.dumps({'states':output,'unfilledStrictOrderPairs':unfilled},indent=2,ensure_ascii=False)+'\n')
    print('Imported',len(states),'numerical records;',len(output),'new exact lower-order representatives;',len(unfilled),'strict-order gaps.')

if __name__=='__main__':main()
