#!/usr/bin/env python3
"""Exact exclusion of a two-dimensional 2-AC subspace in Sym^9(C^2).

Run beside kernel.json and certificate.json. Requires SymPy only; no numerical
optimizer, floating-point certificate, or Clebsch--Gordan library is used.
The source basis consists of physical four-copy amplitudes antisymmetric in
positions (1,2) and (3,4), symmetric under exchanging these pairs, and satisfying
the Pluecker relations. The last condition selects the Young shape (2,2).
"""
import itertools
import json
import math
import time
from fractions import Fraction
from pathlib import Path
import sympy as s

BASE=Path(__file__).resolve().parent

def primitive(row):
    den=math.lcm(*(v.denominator for v in row))
    integers=[int(v*den) for v in row]
    divisor=math.gcd(*integers)
    assert divisor
    return [v//divisor for v in integers]

def canonical(word):
    a,b,c,d=word
    if a==b or c==d:return None,0
    sign=1
    if a>b:a,b=b,a;sign=-sign
    if c>d:c,d=d,c;sign=-sign
    if (a,b)>(c,d):a,b,c,d=c,d,a,b
    return (a,b,c,d),sign

def modular_rank(matrix,prime=1000003):
    A=[[int(v)%prime for v in row] for row in matrix]
    if not A:return 0
    rank=0;cols=len(A[0])
    for col in range(cols):
        pivot=next((i for i in range(rank,len(A)) if A[i][col]),None)
        if pivot is None:continue
        A[rank],A[pivot]=A[pivot],A[rank]
        inv=pow(A[rank][col],-1,prime)
        A[rank]=[(v*inv)%prime for v in A[rank]]
        for i in range(rank+1,len(A)):
            q=A[i][col]
            if q:A[i]=[(x-q*y)%prime for x,y in zip(A[i],A[rank])]
        rank+=1
        if rank==min(len(A),cols):break
    return rank

def verify(directory=None):
    start=time.time();directory=Path(directory) if directory else BASE
    data=json.loads((directory/'kernel.json').read_text())
    cert=json.loads((directory/'certificate.json').read_text())
    n=data['n'];assert n==9
    assert cert['n']==9 and cert['order']==2 and cert['subspaceDimension']==2
    assert cert['copies']==4 and cert['youngShape']==[2,2]
    pairs=list(itertools.combinations(range(n+1),2))
    words=[p+q for p in pairs for q in pairs if p<=q]
    byweight={}
    for w in words:byweight.setdefault(2*sum(w)-4*n,[]).append(w)
    index={M:{w:i for i,w in enumerate(v)} for M,v in byweight.items()}

    def ladder(M,direction):
        src=byweight.get(M,[]);dst=byweight.get(M+2*direction,[])
        A=s.zeros(len(dst),len(src))
        # Output coefficient at a canonical physical word. Each contribution
        # changes one of its four indices; antisymmetry supplies the sign.
        for row,w in enumerate(dst):
            for leg in range(4):
                v=list(w);v[leg]-=direction
                if not 0<=v[leg]<=n:continue
                q,sign=canonical(v)
                if q is not None:
                    A[row,index[M][q]]+=sign*(w[leg] if direction==1 else n-w[leg])
        return A

    plucker_rows={}
    for M in byweight:
        rows=[]
        for w in itertools.combinations(range(n+1),4):
            if 2*sum(w)-4*n!=M:continue
            a,b,c,d=w;row=[0]*len(byweight[M])
            for q,sign in [((a,b,c,d),1),((a,c,b,d),-1),((a,d,b,c),1)]:row[index[M][q]]=sign
            rows.append(row)
        plucker_rows[M]=s.Matrix(rows) if rows else s.zeros(0,len(byweight[M]))

    # Independent weight counting of Sym^2(Lambda^2 V) minus Lambda^4 V.
    young_weights={M:len(byweight[M])-plucker_rows[M].rows for M in byweight}
    expected=[(M,young_weights[M]-young_weights.get(M+2,0)) for M in sorted(byweight) if M>=0 and young_weights[M]>young_weights.get(M+2,0)]
    assert [(b['J2'],b['multiplicity']) for b in data['blocks']]==expected
    blocks=[];highest_count=0
    for block in data['blocks']:
        J=block['J2'];mult=block['multiplicity'];H=s.Matrix(block['highest'])
        assert H.shape==(len(byweight[J]),mult)
        assert modular_rank(block['highest'])==mult
        assert ladder(J,1)*H==s.zeros(len(byweight.get(J+2,[])),mult)
        assert plucker_rows[J]*H==s.zeros(plucker_rows[J].rows,mult)
        U={J:H}
        for step in range(1,J+1):
            M=J-2*step;U[M]=ladder(M+2,-1)*U[M+2]/step
            assert all(v.is_Integer for v in U[M])
            assert plucker_rows[M]*U[M]==s.zeros(plucker_rows[M].rows,mult)
        # Norms use the number (4 or 8) of signed physical words in the orbit.
        norm=[sum((4 if w[:2]==w[2:] else 8)*math.prod(math.comb(n,r) for r in w)*int(H[i,a])**2 for i,w in enumerate(byweight[J])) for a in range(mult)]
        assert norm==block['norms'] and all(x>0 for x in norm)
        blocks.append({'J2':J,'multiplicity':mult,'vectors':U})
        highest_count+=mult
    dimension=sum((b['J2']+1)*b['multiplicity'] for b in blocks)
    assert dimension==825 and highest_count==57
    variables=[tuple(x) for x in data['variables']]
    assert variables==[(ib,a,b) for ib,B in enumerate(blocks) for a in range(B['multiplicity']) for b in range(a,B['multiplicity'])]
    assert len(variables)==169
    lowweight={}
    for a in range(n+1):
        for b,c in pairs:lowweight.setdefault(2*(a+b+c)-3*n,[]).append((a,b,c))

    rows=[]
    for ell,M3,a,b in data['selectedConstraints']:
        assert ell in (1,2)
        alpha=lowweight[M3][a];beta=lowweight[M3][b];row=[]
        for ib,A,B in variables:
            block=blocks[ib];J=block['J2'];total=Fraction(0)
            for i in range(n+1):
                M=M3+2*i-n
                if M not in block['vectors']:continue
                wa,sa=canonical((i,)+alpha);wb,sb=canonical((i,)+beta)
                if wa is None or wb is None:continue
                U=block['vectors'][M];va=U.row(index[M][wa]);vb=U.row(index[M][wb])
                z=2*i-n;f=z if ell==1 else 3*z*z-n*(n+2)
                products=int(va[A]*vb[B]+(va[B]*vb[A] if A!=B else 0))
                total+=Fraction(sa*sb*math.comb(n,i)*f*products,math.comb(J,(J-M)//2))
            row.append(total)
        rows.append(primitive(row))
    assert rows==data['selectedRows'] and len(rows)==156
    # Modular full row rank is a rigorous lower bound on rational rank.
    assert modular_rank(rows)==156
    K=s.Matrix([[s.Integer(v) for v in row] for row in data['kernel']])
    assert K.shape==(169,13)
    assert modular_rank(K.tolist())==13
    assert s.Matrix(rows)*K==s.zeros(156,13)
    # Thus these are exactly all candidates satisfying the selected equations.
    diagonal=[i for i,(ib,a,b) in enumerate(variables) if a==b]
    assert diagonal==cert['diagonalVariableIndices'] and len(diagonal)==57
    coefficients=[s.Integer(v) for v in cert['positiveCoefficients']]
    assert len(coefficients)==57 and all(v>0 for v in coefficients)
    assert s.Matrix([coefficients])*K[diagonal,:]==s.zeros(1,13)
    # All invariant PSD blocks have nonnegative diagonal coefficients. The
    # strict positive identity sets all57 diagonals to zero, so each PSD block
    # is zero, contradicting normalization of any nonzero twirled |w w>.
    return {'passed':True,'n':9,'order':2,'excludedSubspaceDimension':2,
        'appliesToArbitraryComplexSubspaces':True,'copies':4,'youngShape':[2,2],
        'youngSpaceDimension':dimension,'highestWeightVectorsChecked':highest_count,
        'realSymmetricInvariantVariables':169,'selectedContractionEquations':156,
        'exactConstraintRank':156,'exactKernelDimension':13,
        'strictlyPositiveDiagonalTerms':57,'positiveIdentityVerifiedExactly':True,
        'allPositiveSemidefiniteCandidatesVanish':True,'rankVerificationPrime':1000003,
        'method':'Exact physical-word Young basis, rational first-leg contractions, integer nullspace checks and strictly positive diagonal certificate.',
        'runtimeSeconds':time.time()-start}

if __name__=='__main__':
    result=verify();(BASE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
