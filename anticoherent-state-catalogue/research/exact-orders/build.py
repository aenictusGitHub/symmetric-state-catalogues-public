"""Construct exact-order witnesses by positive discrete moment quadrature.

Support separation > t kills off-diagonal t-body reductions. Rational weights
match the first t moments of the uniform Dicke-index distribution. A separate
exact reduced-matrix calculation verifies every claim, including failure at t+1.
Run before scripts/build_data.py and scripts/build_majorana.py.
"""
from pathlib import Path
from itertools import combinations
import json
import numpy as np
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def reduction(c,n,t):
    support={i:x for i,x in enumerate(c) if x!=0}
    result=S.zeros(t+1)
    for a in range(t+1):
        for b in range(a,t+1):
            terms=[]
            for r,x in support.items():
                l=r-a
                if 0<=l<=n-t and b+l in support:
                    terms.append(x*S.conjugate(support[b+l])*S.binomial(n-t,l)*S.sqrt(S.Rational(S.binomial(t,a)*S.binomial(t,b),S.binomial(n,r)*S.binomial(n,b+l))))
            result[a,b]=S.simplify(sum(terms));result[b,a]=S.conjugate(result[a,b])
    return result

def sparse_weights(n,t):
    m=t+1
    # r_i=v_i+i*t, with distinct v_i, enumerates gaps >=t+1.
    target=[sum(S.Integer(r)**q for r in range(n+1))/S.Integer(n+1) for q in range(m)]
    scaled=np.array([float(x/n**q) for q,x in enumerate(target)])
    grids=[tuple(offset+i*step for i in range(m)) for step in range(t+1,n//t+1) for offset in range(n-step*t+1)]
    from itertools import chain
    separated=(tuple(x+i*t for i,x in enumerate(v)) for v in combinations(range(n-t*(m-1)+1),m))
    for nodes in chain(grids,separated):
        mat=np.array([[(x/n)**q for x in nodes] for q in range(m)])
        w=np.linalg.solve(mat,scaled)
        if min(w)<1e-8:continue
        exact=S.Matrix([[S.Integer(x)**q for x in nodes] for q in range(m)]).inv()*S.Matrix(target)
        if not all(x>0 for x in exact):continue
        c=[S.sqrt(exact[nodes.index(r)]) if r in nodes else S.S.Zero for r in range(n+1)]
        # Positive amplitudes make a gap t+1 a strictly nonzero next off-diagonal.
        if any(b-a==t+1 for a,b in zip(nodes,nodes[1:])) or sum(exact[i]*nodes[i]**(t+1) for i in range(m))!=sum(S.Integer(r)**(t+1) for r in range(n+1))/S.Integer(n+1):
            return c,'Positive rational moment weights on separated Dicke indices'
    raise ValueError((n,t,'No separated positive quadrature found'))

def construct(n,t,original):
    if t==1:
        return [1/S.sqrt(2) if r in (0,n) else S.S.Zero for r in range(n+1)],'GHZ state'
    if (n,t)==(16,4):
        c=[S.S.Zero]*(n+1)
        for r,x in {0:S.sqrt(27)/16,4:S.sqrt(21)/8,8:S.sqrt(34)/16,12:-S.sqrt(21)/8,16:S.sqrt(27)/16}.items():c[r]=x
        return c,'Signed five-component fourth-order state'
    if (n,t)==(18,4):
        c=[S.simplify(S.Abs(S.sympify(x))) for x in original['acliege-n18-t5']['exactCoefficients']]
        return c,'Phase deformation of the fifth-order witness'
    return sparse_weights(n,t)

def main():
    source=json.loads((ROOT/'research/states.json').read_text())['states']+json.loads((ROOT/'research/cross-catalogue/exact-states.json').read_text())['states']
    original={s['id']:s for s in source};new=[]
    for n in sorted({s['n'] for s in source}):
        top=max(s['achievedOrder'] for s in source if s['n']==n)
        for t in range(1,top):
            if any(s['n']==n and s['achievedOrder']==t for s in source):continue
            c,method=construct(n,t,original)
            assert S.simplify(sum(x*S.conjugate(x) for x in c))==1
            checks=[];next_witness=None
            for q in range(1,t+2):
                delta=reduction(c,n,q)-S.eye(q+1)/(q+1)
                zero=all(S.simplify(x)==0 for x in delta)
                assert zero==(q<=t),(n,t,q,delta)
                distance=S.simplify(sum(S.expand_complex(x*S.conjugate(x)) for x in delta))
                checks.append({'t':q,'identityVerified':zero,'frobeniusDistanceSquared':str(distance),'purity':str(S.simplify(distance+S.Rational(1,q+1)))})
                if q==t+1:
                    a,b=next((a,b) for a in range(q+1) for b in range(q+1) if delta[a,b]!=0)
                    next_witness={'order':q,'row':a,'column':b,'differenceFromIdentity':str(delta[a,b]),'differenceLatex':S.latex(delta[a,b])}
            support=[r for r,x in enumerate(c) if x!=0]
            terms=[S.latex(c[r])+rf'|D_{{{n}}}^{{({r})}}\rangle' for r in support]
            formula=' + '.join(terms).replace('+ -','- ')
            record={'id':f'exact-order-n{n}-t{t}','n':n,'achievedOrder':t,'witnessesOrders':list(range(1,t+1)),'label':f'Exactly {t}-AC · {len(support)} Dicke components','status':'proven-existence','stateLatex':formula,'support':support,'exactCoefficients':list(map(str,c)),'coefficients':[[float(S.re(x)),float(S.im(x))] for x in c],'normalization':'sum_r |c_r|^2 = 1, verified exactly','coefficientConvention':'Normalized Dicke basis; r counts excitations.','references':[{'title':'Exact-order constructions and independent verification','url':'./research/exact-orders/index.html'}],'provenance':'Independent exact-order construction: '+method,'notes':f'This representative is {t}-AC but not {t+1}-AC. Every reduction through order {t} is verified exactly; a nonzero entry of the next reduced-state deviation certifies strict order.','verification':{'method':'Exact reduced-state matrices in the normalized Dicke basis, including a nonzero deviation at the next order.','exactNormalization':True,'reductions':checks,'exactOrderVerified':True,'nextOrderWitness':next_witness,'globalMaximumOrderProven':False}}
            new.append(record);print(n,t,support,'verified',flush=True)
    (HERE/'states.json').write_text(json.dumps({'states':new},indent=2,ensure_ascii=False)+'\n')
    print('Exact-order witnesses:',len(new))

if __name__=='__main__':main()
