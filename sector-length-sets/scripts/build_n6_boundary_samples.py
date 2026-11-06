#!/usr/bin/env python3
"""Explicit pure-state samples on two already proved N6 boundary loci.

No optimization or interpolation is used. Existing scientific datasets are
read only. The optional exact-expression helper is independently generated.
"""
from pathlib import Path
from fractions import Fraction
from math import sqrt, comb, isfinite
import json

ROOT=Path(__file__).resolve().parents[1]
FO_SOURCE='Julia/N6_filtered_octa_attainment.jl: n6fo_c3_state / n6fo_c3_point'
RANK_SOURCE='Julia/N6_rankfour_equality_segment.jl: n6rank4_equality_state'

def sector_lengths(c):
    """Independent Dicke-marginal purity conversion, no target-sector fit."""
    N=len(c)-1;S=[1.]
    for r in range(1,N+1):
        A=[[c[k+j]*sqrt(comb(r,k)*comb(N-r,j)/comb(N,k+j))
            for j in range(N-r+1)] for k in range(r+1)]
        rho=[[sum(A[k][j]*complex(A[l][j]).conjugate() for j in range(N-r+1))
              for l in range(r+1)] for k in range(r+1)]
        purity=sum(abs(z)**2 for row in rho for z in row)
        S.append(comb(N,r)*(2**r*purity-sum(comb(r,k)/comb(N,k)*S[k] for k in range(r))))
    return [0. if abs(s)<1e-11 else s for s in S[1:]]

def overlap(c):
    return abs(sum((-1)**k*c[k]*c[6-k] for k in range(7)))**2

def fo_state(t):
    r=1-t;s=1-5*r/9;d=sqrt(s*s-16*r*r/81)
    return [sqrt((s+d)/2),0.,0.,sqrt(5*r/9),0.,0.,-sqrt(max((s-d)/2,0.))]

def fo_point(t):
    return [2*t*(t+8)/3,5+10*t*t,20*t*(4-t)/3]

def rank_state(t):
    g=13-7*t*t;den=427-639*t+287*t*t-35*t**3
    a=sqrt(35*(7-5*t*t)/(16*g))
    pop=25*(23*t*t-84*t+91)/(108*g)
    u=-3*(7-3*t)*g/den;v=-24*a*(t-1)*g/(5*den)
    pair=5*a*(1985-3484*t+1594*t*t+84*t**3-119*t**4)/(54*g*g)
    discriminant=pop*pop-4*pair*pair
    assert discriminant>=-1e-12
    d=sqrt((pop+sqrt(max(discriminant,0.)))/2);e=pair/d
    b=u*d+v*e;c=-v*d-u*e
    p=[a/3-c/4-e/20,b/2-1+3*d/10,a-c/4+3*e/4,b-d,
       a+c/4-3*e/4,b/2+1+3*d/10,a/3+c/4+e/20]
    normalization=64*(5-2*t*t)/(27*g)
    return [p[k]/sqrt(normalization*comb(6,k)) for k in range(7)]

def rank_point(q):
    return [1.5-474*q/169,9.9+3234*q/845,11-2300*q/169]

def main():
    exact_path=ROOT/'data/n6-boundary-sample-exact.json'
    exact=json.loads(exact_path.read_text()) if exact_path.exists() else {}
    nodes=json.loads((ROOT/'data/n6-manuscript-nodes.json').read_text())
    originals=[n['point'] for n in nodes]
    result=[];audits=[]
    families=[('fo-c3',[Fraction(j,9) for j in range(1,9)],fo_state),
              ('rank4-segment',[Fraction(j,40) for j in (41,42,43,44,45,46,47)],rank_state)]
    for family,parameters,state_factory in families:
        for j,t0 in enumerate(parameters,1):
            t=float(t0);c=state_factory(t);S=sector_lengths(c);q=overlap(c)
            target=fo_point(t) if family=='fo-c3' else rank_point(q)
            norm=sum(abs(z)**2 for z in c)
            residual=max(abs(S[k]-target[k]) for k in range(3))
            nearest=min(sqrt(sum((S[k]-p[k])**2 for k in range(3))) for p in originals)
            assert abs(norm-1)<2e-12,(family,t,norm)
            assert residual<3e-10,(family,t,residual)
            assert nearest>1e-7,(family,t,'duplicate original node')
            assert all(isfinite(z) for z in S)
            id=f'n6-{family}-{j}'
            if family=='fo-c3':
                name='Filtered-octa C₃ boundary';symbol='F';origin='Exact filtered-octa C₃ boundary on the p₂ = 0 plane'
                boundary_residual=abs(S[2]-5*S[0]-5+S[1])
            else:
                name='Rank-four equality family';symbol='R';origin='Complete sharp-plane equality classification: attained C₅–ω segment'
                boundary_residual=abs(65*S[1]+63*S[2]-217*S[0]-1011)
            rec={'id':id,'N':6,'status':'attained','family':name,
                 'label':f'{name} (t = {t0})','mathLabel':rf'$|{symbol}_{{{t0}}}\rangle$',
                 'parameter':{'t':t},'exactParameter':str(t0),
                 'coefficients':[[z,0.] for z in c],'sectors':S,
                 'proofOrigin':origin,'source':FO_SOURCE if family=='fo-c3' else RANK_SOURCE,
                 'boundaryRole':'attainment-color seam' if family=='fo-c3' else 'attained segment inside a supporting face',
                 'q':q}
            if id in exact:rec.update(exact[id])
            result.append(rec)
            audits.append({'id':id,'normError':abs(norm-1),'sectorTargetError':residual,
                           'boundaryEquationError':boundary_residual,'nearestOriginalNodeDistance':nearest,
                           'exactExpressionsPresent':id in exact})
    payload={'N':6,'samples':result,'paths':[{'id':'rank4-equality-segment','kind':'attained-state-locus','attained':True,
        'points':[rank_point((169/484)*j/256) for j in range(257)],
        'source':RANK_SOURCE,'note':'Proved attained straight sector image, not a geometric body crease.'},
        {'id':'filtered-octa-C3-seam','kind':'attained-state-locus','attained':True,
         'points':[fo_point(j/720) for j in range(721)],'source':FO_SOURCE,
         'note':'Exact attained filtered-octahedral C3 boundary, sampled densely through every displayed F_t point.'}],'families':[
        {'id':'filtered-octa-C3','symbol':'F','count':8,'role':'attained violet/gray seam on p₂ = 0 lower plane',
         'source':FO_SOURCE,'domain':'0 ≤ t ≤ 1','sampleParameters':'t=j/9, j=1,…,8'},
        {'id':'rank4-equality','symbol':'R','count':7,'role':'attained straight C₅–ω segment inside supporting face; not excluded curved crease',
         'source':RANK_SOURCE,'domain':'1 ≤ t ≤ sqrt(7/5)','sampleParameters':'t=(40+j)/40, j=1,…,7'}],
        'notes':['All coefficients are actual normalized pure-state witnesses, not convex interpolations.',
                 'The nearby curved rank-four/spectral outer crease is not declared attained.',
                 'Existing 27 manuscript nodes and geometry are unchanged. Requested samples add separate violet selectable points.']}
    report={'status':'PASS','sampleCount':len(result),'foCount':8,'rank4Count':7,
            'maximumNormError':max(a['normError'] for a in audits),
            'maximumSectorTargetError':max(a['sectorTargetError'] for a in audits),
            'maximumBoundaryEquationError':max(a['boundaryEquationError'] for a in audits),
            'allExactExpressionsPresent':all(a['exactExpressionsPresent'] for a in audits),'sampleChecks':audits,
            'proofNote':'audit/N6_boundary_sample_proofs.txt',
            'retainedCreaseProofAudit':'audit/n6-retained-crease-proof-audit.json'}
    (ROOT/'data/n6-boundary-samples.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n')
    (ROOT/'data/n6-boundary-samples.js').write_text('globalThis.SECTOR_N6_BOUNDARY_SAMPLES = '+json.dumps(payload,ensure_ascii=False,separators=(',',':'))+';\n')
    (ROOT/'audit/n6-boundary-samples-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sampleChecks'},indent=2))

if __name__=='__main__':main()
