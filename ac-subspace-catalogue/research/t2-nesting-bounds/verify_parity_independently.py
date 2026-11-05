#!/usr/bin/env python3
"""Independent Dicke-cross-reduction and rational root certificate audit.
Does not import the constructor or primary verifier.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
from collections import defaultdict
import argparse,json,hashlib
import sympy as s
HERE=Path(__file__).resolve().parent
CONST=(-1,-1)

def clean(p):return {m:c for m,c in p.items() if c}
def linear(*terms):
    out=defaultdict(F)
    for weight,p in terms:
        for m,c in p.items():out[m]+=weight*c
    return clean(out)
def pair(slots1,slots2,weight,conjugate_first=True):
    a,b=slots1;u,v=slots2
    re=defaultdict(F);im=defaultdict(F)
    for i,j,w in [(a,u,weight),(b,v,weight)]:re[tuple(sorted((i,j)))]+=w
    sign=1 if conjugate_first else -1
    for i,j,w in [(a,v,sign*weight),(b,u,-sign*weight)]:im[tuple(sorted((i,j)))]+=w
    return clean(re),clean(im)
def add_into(dst,src):
    for m,c in src.items():dst[m]+=c

def audit(folder):
    folder=Path(folder);c=json.loads((folder/'certificate.json').read_text());b=json.loads((folder/'basis.json').read_text());sys=json.loads((folder/'system.json').read_text())
    n=c['n'];k=c['dimension'];assert (n,k)in[(35,5),(42,6)]
    assert (b['n'],b['t'],b['dimension'])==(n,2,k)
    parity=c['parity'];rs=list(range(parity,n+1,2));d=len(rs);size=2*d*k
    exps=c['scaleExponents'];assert len(exps)==d and all(type(e)is int and e>=0 for e in exps)
    rowof={r:a for a,r in enumerate(rs)};j=F(n,2);variance=F(n*(n+2),12)
    def slots(r,col):
        a=rowof[r]*k+col;return a,a+d*k
    def scaling(r,q):return F(1,2**(exps[rowof[r]]+exps[rowof[q]]))
    # Build physical moment matrices from the magnetic eigenvalues and ladder action.
    moments={}
    for kind in range(4):
        for i in range(k):
            for l in range(k):
                parts=[defaultdict(F),defaultdict(F)]
                for r in rs:
                    q=r+(2 if kind==3 else 0)
                    if q not in rowof:continue
                    w=F(comb(n,r))*scaling(r,q)
                    if kind==1:w*=F(n-2*r,n)
                    if kind==2:w*=((F(n,2)-r)**2-variance)/(j*j)
                    if kind==3:w*=F((n-r)*(n-r-1))/(j*j)
                    pp=pair(slots(r,i),slots(q,l),w)
                    for z in range(2):add_into(parts[z],pp[z])
                if kind==0 and i==l:parts[0][CONST]-=1
                moments[kind,i,l]=tuple(clean(p) for p in parts)
    # Hermitian completion: no omitted diagonal imaginary part or lower triangle remains.
    for kind in range(3):
        for i in range(k):
            assert moments[kind,i,i][1]=={}
            for l in range(i):
                assert moments[kind,i,l][0]==moments[kind,l,i][0]
                assert moments[kind,i,l][1]==linear((-1,moments[kind,l,i][1]))
    # All saved equations are the complete physical system, including full ordered J+^2.
    physical=[]
    for kind in range(3):
        for i in range(k):
            for l in range(i,k):
                for z in range(1 if i==l else 2):physical.append((kind,i,l,z,moments[kind,i,l][z]))
    for i in range(k):
        for l in range(k):
            for z in range(2):physical.append((3,i,l,z,moments[3,i,l][z]))
    equations=sys['equations'];m=len(physical);L=int(c['systemScale'])
    assert m==5*k*k and len(equations)==m and (sys['n'],sys['dimension'])==(n,k)
    for (kind,i,l,z,p),e in zip(physical,equations):
        assert (e['kind'],e['basisPair'],e['part'])==(kind,[i,l],z)
        saved={tuple((a,q)):F(w,L)for a,q,w in e['terms']};saved[CONST]=F(e['constant'],L)
        assert clean(saved)==p
    # Direct Dicke splitting: dividing by sqrt(C(2,a)C(2,b)) leaves rational polynomials.
    # Each of all 9 k^2 complex reduced entries is explicitly reconstructed here.
    ratio=F(n,4*(n-1));cross_count=0
    for i in range(k):
        for l in range(k):
            for a in range(3):
                for q in range(3):
                    parts=[defaultdict(F),defaultdict(F)]
                    for t in range(n-1):
                        r=a+t;v=q+t
                        if r not in rowof or v not in rowof:continue
                        pp=pair(slots(r,i),slots(v,l),F(comb(n-2,t))*scaling(r,v),conjugate_first=False)
                        for z in range(2):add_into(parts[z],pp[z])
                    if i==l and a==q:parts[0][CONST]-=F(1,3*comb(2,a))
                    parts=tuple(clean(p)for p in parts)
                    if a==q:
                        coeff=[(F(1,3),F(1,2),ratio),(F(1,6),F(0),-ratio),(F(1,3),F(-1,2),ratio)][a]
                        expect=tuple(linear(*[(v,moments[h,i,l][z])for h,v in enumerate(coeff)])for z in range(2))
                        expect=(expect[0],linear((-1,expect[1]))) # physical reduction reverses bra/ket ordering
                    elif (a,q)==(0,2):expect=(linear((ratio,moments[3,i,l][0])),linear((-ratio,moments[3,i,l][1])))
                    elif (a,q)==(2,0):expect=tuple(linear((ratio,moments[3,l,i][z]))for z in range(2))
                    else:expect=({}, {})
                    assert parts==expect,(n,i,l,a,q)
                    cross_count+=1
    assert cross_count==9*k*k
    # Basis formulas and the isolated box are tied to exactly these same coordinates.
    U=c['unknownCoordinateIndices'];assert len(U)==m and len(set(U))==m and all(0<=i<size for i in U)
    D=int(c['centerDenominator']);center=list(map(int,c['fullCenterNumerators']));radius=F(c['radius'])
    assert len(center)==size and D>0 and radius>0
    definition=b['algebraicDefinition'];names=definition['variables'];assert len(names)==m and len(set(names))==m
    assert definition['unknownCoordinateIndices']==U and definition['fullCenterNumerators']==c['fullCenterNumerators']
    assert int(definition['centerDenominator'])==D and F(definition['radius'])==radius
    assert definition['scaleExponents']==exps and definition['parityIndices']==rs
    intervals=definition['parameterIntervals'];assert set(intervals)==set(names)
    for name,idx in zip(names,U):assert list(map(F,intervals[name]))==[F(center[idx],D)-radius,F(center[idx],D)+radius]
    xs=s.symbols(' '.join(names),real=True);local=dict(zip(names,xs));atoms=[s.Rational(v,D)for v in center]
    for idx,x in zip(U,xs):atoms[idx]=x
    assoc=0;error=0
    for col,row in enumerate(b['basis']):
        assert len(row['exactCoefficients'])==n+1 and row['support']==rs
        for r,expr in enumerate(row['exactCoefficients']):
            actual=s.sympify(expr,locals=local)
            if r in rowof:
                a,im=slots(r,col);expected=s.sqrt(comb(n,r))/2**exps[rowof[r]]*(atoms[a]+s.I*atoms[im])
                approx=(complex(center[a],center[im])/D)*(comb(n,r)**.5/2**exps[rowof[r]])
            else:expected=s.S.Zero;approx=0j
            assert s.expand(actual-expected)==0
            error=max(error,abs(approx-complex(*row['coefficients'][r])));assoc+=1
    assert error<5e-14
    # Independent integer evaluation of the Banach contraction inequalities.
    inv=[[int(x)for x in row]for row in c['inverseNumerators']];ad=int(c['inverseDenominator'])
    assert len(inv)==m and all(len(row)==m for row in inv) and ad>0
    f=[];J=[]
    for e in equations:
        val=e['constant']*D*D;grad=[0]*size
        for a,q,w in e['terms']:
            val+=w*center[a]*center[q];grad[a]+=w*center[q];grad[q]+=w*center[a]
        f.append(val);J.append([grad[idx]for idx in U])
    columns=[[(a,J[a][q])for a in range(m)if J[a][q]]for q in range(m)]
    defect=F(max(sum(abs((ad*D if i==q else 0)-sum(inv[i][a]*v for a,v in columns[q]))for q in range(m))for i in range(m)),ad*D)
    correction=F(max(abs(sum(inv[i][a]*f[a]for a in range(m)))for i in range(m)),ad*D*D)
    invnorm=F(max(sum(abs(x)for x in row)for row in inv),ad)
    active=set(U);H=max(2*sum(abs(w)for a,q,w in e['terms']if a in active and q in active)for e in equations)
    contraction=defect+invnorm*H*radius;selfmap=correction+contraction*radius
    assert defect<1 and contraction<1 and selfmap<radius
    computed={'inverseDefect':defect,'newtonCorrection':correction,'inverseInfinityNorm':invnorm,'jacobianLipschitzBound':F(H),'contractionFactor':contraction,'selfMapRadiusBound':selfmap}
    assert all(F(c['bounds'][key])==val for key,val in computed.items())
    report={'passed':True,'n':n,'order':2,'dimension':k,'realVariables':m,'fixedRationalCoordinates':size-m,'allPhysicalMomentEquationsReconstructed':True,'fullCrossReductionEntriesChecked':cross_count,'directDickeCrossPolynomialIdentitiesVerified':True,'coefficientExpressionsAssociated':assoc,'isolatingBoxAssociated':True,'strictContraction':True,'strictSelfMap':True,'radius':str(radius),'contractionFactor':str(contraction),'selfMapRadiusBound':str(selfmap),'maximumDisplayCoefficientError':error,'certificateSHA256':hashlib.sha256((folder/'certificate.json').read_bytes()).hexdigest(),'basisSHA256':hashlib.sha256((folder/'basis.json').read_bytes()).hexdigest(),'method':'Independent full-coordinate polynomial reconstruction of every ordered two-qubit cross reduction, exact comparison to the complete Gram/spin system, independent rational Banach bounds, and exact association of exported algebraic coefficient formulas. No constructor or primary verifier imported.','scope':'Proves the stated algebraic lower witness; no unrestricted optimality or radical-expression claim.'}
    (folder/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('folder');a=p.parse_args();r=audit(a.folder);print(json.dumps({k:v for k,v in r.items()if k not in ['contractionFactor','selfMapRadiusBound']},indent=2))
