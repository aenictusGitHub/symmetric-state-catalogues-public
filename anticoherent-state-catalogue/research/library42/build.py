"""Import independently checked refined library records through the configured maximum qubit number.

Uses the saved, audited refinement rather than mismatched original Dicke files.
The original library remains read-only. Exact lower-order representatives are
constructed independently where separated positive moment quadrature succeeds.
"""
from pathlib import Path
import importlib.util,json,hashlib
import sympy as S
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SIBLING=ROOT.parent/'majorana-degeneracy-catalogue'
MAX_N=json.loads((HERE/'config.json').read_text())['maxQubits']

def grid_weights(n,t):
    """A bounded, reproducible part of the separated-moment prescription.

    Uniform separated grids keep the reduced Majorana polynomial small. A
    missing grid is recorded as a construction gap, never as nonexistence.
    """
    m=t+1
    target=[sum(S.Integer(r)**q for r in range(n+1))/S.Integer(n+1)
            for q in range(m)]
    scaled=np.array([float(x/n**q) for q,x in enumerate(target)])
    for step in range(t+1,n//t+1):
        for offset in range(n-step*t+1):
            nodes=tuple(offset+i*step for i in range(m))
            mat=np.array([[(x/n)**q for x in nodes] for q in range(m)])
            if min(np.linalg.solve(mat,scaled))<1e-8:continue
            exact=S.Matrix([[S.Integer(x)**q for x in nodes] for q in range(m)]).inv()*S.Matrix(target)
            if not all(x>0 for x in exact):continue
            next_moment=sum(exact[i]*nodes[i]**(t+1) for i in range(m))
            uniform_next=sum(S.Integer(r)**(t+1) for r in range(n+1))/S.Integer(n+1)
            if step>t+1 and next_moment==uniform_next:continue
            c=[S.sqrt(exact[nodes.index(r)]) if r in nodes else S.S.Zero for r in range(n+1)]
            return c,'Positive rational moment weights on a separated uniform Dicke grid'
    raise ValueError((n,t,'No positive separated uniform grid recorded'))

def verify_grid(c,n,t):
    """Check every nonzero reduced entry using exact hypergeometric sums."""
    support=[r for r,x in enumerate(c) if x!=0]
    weights={r:S.simplify(c[r]**2) for r in support}
    assert all(w>0 for w in weights.values()) and sum(weights.values())==1
    assert min((b-a for a,b in zip(support,support[1:])),default=n+1)>t
    reductions=[];witness=None
    for q in range(1,t+2):
        diagonal=[sum(weights[r]*S.binomial(q,a)*S.binomial(n-q,r-a)/S.binomial(n,r)
                      for r in support if 0<=r-a<=n-q)-S.Rational(1,q+1)
                  for a in range(q+1)]
        off=S.S.Zero
        if q==t+1:
            off=sum(S.sqrt(weights[r]*weights[r+q]*S.binomial(n-q,r)**2 /
                           (S.binomial(n,r)*S.binomial(n,r+q)))
                    for r in support if r+q in weights and r<=n-q)
        # With separation at least q, the only possible next-order
        # off-diagonal entries are (0,q) and (q,0).
        zero=all(x==0 for x in diagonal) and off==0
        assert zero==(q<=t),(n,t,q)
        distance=sum(x*x for x in diagonal)+2*off**2
        reductions.append({'t':q,'identityVerified':zero,
                           'frobeniusDistanceSquared':str(distance),
                           'purity':str(distance+S.Rational(1,q+1))})
        if q==t+1:
            nonzero=next((a for a,x in enumerate(diagonal) if x!=0),None)
            if nonzero is None:
                assert off.is_positive
                a,b,value=0,q,off
            else:a=b=nonzero;value=diagonal[a]
            witness={'order':q,'row':a,'column':b,
                     'differenceFromIdentity':str(value),'differenceLatex':S.latex(value)}
    return reductions,witness

def main():
    report=json.loads((HERE/'verification.json').read_text())
    checks={x['id']:x for x in report['records']}
    deg=json.loads((SIBLING/'data/catalogue.json').read_text())
    refined=json.loads((SIBLING/'research/liege-refinement/refined.json').read_text())['records']
    if isinstance(refined,list):refined={r['id']:r for r in refined}
    states=[]
    for source in deg['numericalStates']:
        if source['n']>MAX_N:continue
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
    for extras in sorted(HERE.glob('extra-states-*.json')):
        states.extend(s for s in json.loads(extras.read_text())['states'] if s['n']<=MAX_N)
    (HERE/'states.json').write_text(json.dumps({'states':states,'sourceLibrary':'Documents/libqiss/ACstates','recordCount':len(states)},indent=2,ensure_ascii=False)+'\n')
    # Reuse the established exact construction and verify each next-order failure.
    spec=importlib.util.spec_from_file_location('exact_orders',ROOT/'research/exact-orders/build.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    old=[]
    for name in ['research/states.json','research/cross-catalogue/exact-states.json','research/exact-orders/states.json']:
        old.extend(json.loads((ROOT/name).read_text())['states'])
    pairs={(s['n'],s['achievedOrder']) for s in old}
    previous_file=HERE/'exact-lower-states.json'
    previous={r['id']:r for r in json.loads(previous_file.read_text())['states']} if previous_file.exists() else {}
    output=[r for r in previous.values() if r['n']<=MAX_N]
    completed={(r['n'],r['achievedOrder']) for r in output}
    old_unfilled=set(map(tuple,json.loads(previous_file.read_text()).get('unfilledStrictOrderPairs',[]))) if previous_file.exists() else set()
    unfilled=[]
    for n in range(2,MAX_N+1):
        top=max([s['achievedOrder'] for s in states+old if s['n']==n],default=0)
        for t in range(1,top):
            if (n,t) in pairs or (n,t) in completed:continue
            if (n,t) in old_unfilled:
                unfilled.append([n,t]);continue
            try:
                if t==1:c=[1/S.sqrt(2) if r in (0,n) else S.S.Zero for r in range(n+1)];method='GHZ state'
                else:c,method=grid_weights(n,t)
            except ValueError:
                unfilled.append([n,t]);continue
            support=[i for i,x in enumerate(c) if x!=0]
            reductions,next_witness=verify_grid(c,n,t)
            formula=' + '.join(S.latex(c[r])+rf'|D_{{{n}}}^{{({r})}}\rangle' for r in support).replace('+ -','- ')
            output.append({'id':f'exact-order-n{n}-t{t}','n':n,'achievedOrder':t,'witnessesOrders':list(range(1,t+1)),'label':f'Exactly {t}-AC · {len(support)} Dicke components','status':'proven-existence','stateLatex':formula,'support':support,'exactCoefficients':list(map(str,c)),'coefficients':[[float(S.re(x)),float(S.im(x))] for x in c],'majoranaCoordinatePolicy':'numeric-roots-retain-explicit-poles','normalization':'sum_r |c_r|^2 = 1, verified exactly','coefficientConvention':'Normalized Dicke basis; r counts excitations.','references':[{'title':'Exact-order construction and checks through the configured maximum qubit number','url':'./research/library42/index.html'},{'title':'General analytical prescription','url':'./research/exact-orders/general/index.html'}],'provenance':'Independent exact-order construction: '+method,'notes':f'This state is {t}-AC but not {t+1}-AC. All reduced-state identities and a nonzero next-order entry are checked exactly.','verification':{'method':'Exact symbolic normalized Dicke-basis partial traces through the next order.','exactNormalization':True,'exactOrderVerified':True,'globalMaximumOrderProven':False,'reductions':reductions,'nextOrderWitness':next_witness}})
            if output[-1]['id'] in previous:
                policy=previous[output[-1]['id']].get('majoranaCoordinatePolicy')
                if policy:output[-1]['majoranaCoordinatePolicy']=policy
                else:output[-1].pop('majoranaCoordinatePolicy',None)
            print('Exact',n,t,flush=True)
    output.sort(key=lambda r:(r['n'],r['achievedOrder']))
    temporary=HERE/'exact-lower-states.temporary.json'
    temporary.write_text(json.dumps({'states':output,'unfilledStrictOrderPairs':unfilled},indent=2,ensure_ascii=False)+'\n')
    temporary.replace(HERE/'exact-lower-states.json')
    print('Imported',len(states),'numerical records;',len(output),'new exact lower-order representatives;',len(unfilled),'strict-order gaps.')

if __name__=='__main__':main()
