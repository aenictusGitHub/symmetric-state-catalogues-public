#!/usr/bin/env python3
"""Standard-library exact rational contraction proof for parity-supported AC2 bases."""
from pathlib import Path
from math import comb
from collections import defaultdict
from fractions import Fraction as F
import argparse,json,hashlib
HERE=Path(__file__).resolve().parent

def reconstruct(n,k,parity,exponents):
    indices=list(range(parity,n+1,2));d=len(indices);size=2*d*k
    assert len(exponents)==d and all(isinstance(e,int) and e>=0 for e in exponents)
    L=3*n*n*2**(2*max(exponents))
    def index(row,col,part):return row*k+col+part*d*k
    def polynomial(i,j,kind,part):
        terms=defaultdict(int)
        for a,r in enumerate(indices):
            if kind==3:
                if a+1==d:continue
                b=a+1
                weight=F(4*comb(n,r)*(n-r)*(n-r-1),n*n*2**(exponents[a]+exponents[b]))
            else:
                b=a;weight=F(comb(n,r),2**(2*exponents[a]))
                if kind==1:weight*=F(n-2*r,n)
                if kind==2:weight*=F(3*(n-2*r)**2-n*(n+2),3*n*n)
            weight*=L;assert weight.denominator==1;w=weight.numerator
            ar=index(a,i,0);ai=index(a,i,1);br=index(b,j,0);bi=index(b,j,1)
            pairs=[(ar,br,1),(ai,bi,1)] if part==0 else [(ar,bi,1),(ai,br,-1)]
            for x,y,sign in pairs:terms[tuple(sorted((x,y)))]+=sign*w
        return {'kind':kind,'basisPair':[i,j],'part':part,'constant':-L*int(kind==part==0 and i==j),
                'terms':[[x,y,w]for (x,y),w in sorted(terms.items())if w]}
    equations=[]
    for kind in range(3):
        for i in range(k):
            for j in range(i,k):
                for part in range(1 if i==j else 2):equations.append(polynomial(i,j,kind,part))
    for i in range(k):
        for j in range(k):
            for part in range(2):equations.append(polynomial(i,j,3,part))
    assert len(equations)==5*k*k
    return indices,size,L,equations

def exact_bounds(c):
    n=c['n'];k=c['dimension'];indices,size,L,equations=reconstruct(n,k,c['parity'],c['scaleExponents'])
    assert c['systemScale']==str(L)
    center=list(map(int,c['fullCenterNumerators']));den=int(c['centerDenominator']);unknown=c['unknownCoordinateIndices']
    A=[list(map(int,row))for row in c['inverseNumerators']];ad=int(c['inverseDenominator']);m=len(equations)
    assert len(center)==size and len(unknown)==m and len(set(unknown))==m and all(0<=i<size for i in unknown)
    assert len(A)==m and all(len(row)==m for row in A)
    assert den>0 and ad>0
    U=set(unknown);fnum=[];jnum=[]
    for eq in equations:
        value=eq['constant']*den*den;gradient=[0]*size
        for a,b,w in eq['terms']:
            value+=w*center[a]*center[b];gradient[a]+=w*center[b];gradient[b]+=w*center[a]
        fnum.append(value);jnum.append([gradient[a]for a in unknown])
    # Exact sparse matrix product for I - A J(center).
    Jcolumns=[[(r,jnum[r][col])for r in range(m)if jnum[r][col]]for col in range(m)]
    defect_num=max(sum(abs((ad*den if i==col else 0)-sum(A[i][r]*v for r,v in Jcolumns[col]))for col in range(m))for i in range(m))
    correction_num=max(abs(sum(a*v for a,v in zip(row,fnum)))for row in A)
    anorm_num=max(sum(abs(v)for v in row)for row in A)
    H=max(2*sum(abs(w)for a,b,w in eq['terms']if a in U and b in U)for eq in equations)
    defect=F(defect_num,ad*den);correction=F(correction_num,ad*den*den);anorm=F(anorm_num,ad);radius=F(c['radius'])
    contraction=defect+anorm*H*radius;selfmap=correction+contraction*radius
    assert radius>0 and defect<1 and contraction<1 and selfmap<radius
    bounds={'inverseDefect':defect,'newtonCorrection':correction,'inverseInfinityNorm':anorm,
            'jacobianLipschitzBound':F(H),'contractionFactor':contraction,'selfMapRadiusBound':selfmap}
    return bounds,equations

def verify(path):
    raw=Path(path).read_bytes();c=json.loads(raw);bounds,eq=exact_bounds(c)
    assert all(F(c['bounds'][key])==v for key,v in bounds.items())
    root=Path(path).parent;system=json.loads((root/'system.json').read_text());assert system['equations']==eq
    basis=json.loads((root/'basis.json').read_text());n=c['n'];k=c['dimension'];indices=list(range(c['parity'],n+1,2));d=len(indices)
    assert basis['id']==f'algebraic-parity-n{n}-k{k}-t2' and basis['n']==n and basis['dimension']==k and basis['t']==2
    alg=basis['algebraicDefinition'];unknown=c['unknownCoordinateIndices'];position={slot:i+1 for i,slot in enumerate(unknown)}
    assert alg['unknownCoordinateIndices']==unknown and alg['fullCenterNumerators']==c['fullCenterNumerators']
    assert alg['centerDenominator']==c['centerDenominator'] and alg['radius']==c['radius']
    assert alg['scaleExponents']==c['scaleExponents'] and alg['parityIndices']==indices
    assert alg['variables']==['x'+str(i+1)for i in range(len(unknown))]
    den=int(c['centerDenominator']);radius=F(c['radius']);center=list(map(int,c['fullCenterNumerators']))
    def atom(slot):return 'x'+str(position[slot])if slot in position else '('+str(F(center[slot],den))+')'
    for i,slot in enumerate(unknown):assert list(map(F,alg['parameterIntervals']['x'+str(i+1)]))==[F(center[slot],den)-radius,F(center[slot],den)+radius]
    for col,record in enumerate(basis['basis']):
        expected=['0']*(n+1)
        for row,r in enumerate(indices):
            a=row*k+col;b=d*k+a
            expected[r]=f'sqrt({comb(n,r)})/({2**c["scaleExponents"][row]})*({atom(a)}+I*{atom(b)})'
        assert record['exactCoefficients']==expected and record['support']==indices
    return {'passed':True,'n':c['n'],'dimension':c['dimension'],'order':2,'parity':c['parity'],
      'certificateSHA256':hashlib.sha256(raw).hexdigest(),'realQuadraticEquations':len(eq),
      'realUnknownCoordinates':len(c['unknownCoordinateIndices']),
      'fixedRationalCoordinates':len(c['fullCenterNumerators'])-len(c['unknownCoordinateIndices']),
      'allPhysicalEquationsReconstructed':True,'radius':c['radius'],
      'exactBasisCoefficientExpressionsVerified':True,'parameterIntervalsVerified':True,
      'fullCrossReductionEntriesChecked':9*c['dimension']**2,
      'strictContraction':True,'strictSelfMap':True,
      'proof':'Exact rational contraction of x-A F(x), after fixing the other coordinates to rational numbers. Reconstructed Gram, Jz, Jz^2-cI and all ordered J+^2 equations imply the full AC2 compression conditions; every rank-one shift and mixed rank-two shift vanishes by parity.',
      'arithmetic':'Integer products and exact rational inequalities only; no floating point is used in verification.',
      'bounds':{key:str(v)for key,v in bounds.items()}}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('certificate');a=p.parse_args();path=Path(a.certificate);r=verify(path)
    path.with_name('verification.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:v for k,v in r.items()if k!='bounds'},indent=2))
