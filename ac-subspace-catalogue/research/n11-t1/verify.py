#!/usr/bin/env python3
"""Exact exclusion of a five-dimensional 1-AC subspace in Sym^11(C^2).

Run beside kernel.json and certificate.json. Requires SymPy only. This verifier
uses exact integer/rational arithmetic; it does not run a numerical optimizer.
The source basis is the physical fifth exterior power of the spin-11/2 space.
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
    if len(set(word))!=len(word):return None,0
    sign=(-1)**sum(word[a]>word[b] for a in range(len(word)) for b in range(a+1,len(word)))
    return tuple(sorted(word)),sign

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
    n=data['n'];k=data['k'];assert (n,k,data['order'],data['exteriorPower'])==(11,5,1,5)
    assert (cert['n'],cert['order'],cert['excludedDimension'],cert['copies'])==(11,1,5,5)
    words=list(itertools.combinations(range(n+1),k));byweight={}
    for w in words:byweight.setdefault(2*sum(w)-k*n,[]).append(w)
    index={M:{w:i for i,w in enumerate(v)} for M,v in byweight.items()}
    def ladder(M,direction):
        src=byweight.get(M,[]);dst=byweight.get(M+2*direction,[])
        A=s.zeros(len(dst),len(src))
        # b_i=sqrt(binomial(n,i))|i> has J_+ b_i=(i+1)b_{i+1}.
        # The following equivalent formula assembles each output word.
        for row,w in enumerate(dst):
            for leg in range(k):
                v=list(w);v[leg]-=direction
                if not 0<=v[leg]<=n:continue
                q,sign=canonical(v)
                if q is not None:
                    A[row,index[M][q]]+=sign*(w[leg] if direction==1 else n-w[leg])
        return A
    expected=[(M,len(byweight[M])-len(byweight.get(M+2,[]))) for M in sorted(byweight) if M>=0 and len(byweight[M])>len(byweight.get(M+2,[]))]
    assert [(b['J2'],b['multiplicity']) for b in data['blocks']]==expected
    blocks=[];highest_count=0;lowered_norm_checks=0
    for block in data['blocks']:
        J=block['J2'];mult=block['multiplicity'];H=s.Matrix(block['highest'])
        assert H.shape==(len(byweight[J]),mult)
        assert modular_rank(block['highest'])==mult
        assert ladder(J,1)*H==s.zeros(len(byweight.get(J+2,[])),mult)
        U={J:H}
        for step in range(1,J+1):
            M=J-2*step;U[M]=ladder(M+2,-1)*U[M+2]/step
            assert all(v.is_Integer for v in U[M])
        norms=[sum(math.factorial(k)*math.prod(math.comb(n,r)for r in w)*int(H[i,a])**2 for i,w in enumerate(byweight[J]))for a in range(mult)]
        assert norms==block['norms'] and all(x>0 for x in norms)
        for M,V in U.items():
            for a in range(mult):
                norm=sum(math.factorial(k)*math.prod(math.comb(n,r)for r in w)*int(V[i,a])**2 for i,w in enumerate(byweight[M]))
                assert norm==math.comb(J,(J-M)//2)*norms[a]
                lowered_norm_checks+=1
        blocks.append({'J2':J,'multiplicity':mult,'vectors':U})
        highest_count+=mult
    dimension=sum((b['J2']+1)*b['multiplicity']for b in blocks)
    assert dimension==data['dimension']==math.comb(12,5)==792 and highest_count==49
    variables=[tuple(x)for x in data['variables']]
    assert variables==[(ib,a,b)for ib,B in enumerate(blocks)for a in range(B['multiplicity'])for b in range(a,B['multiplicity'])]
    assert len(variables)==111
    lowweight={}
    for w in itertools.combinations(range(n+1),k-1):lowweight.setdefault(2*sum(w)-(k-1)*n,[]).append(w)
    rows=[]
    for M4,a,b in data['selectedConstraints']:
        alpha=lowweight[M4][a];beta=lowweight[M4][b];row=[]
        for ib,A,B in variables:
            block=blocks[ib];J=block['J2'];total=Fraction(0)
            for i in range(n+1):
                M=M4+2*i-n
                if M not in block['vectors']:continue
                wa,sa=canonical((i,)+alpha);wb,sb=canonical((i,)+beta)
                if wa is None or wb is None:continue
                U=block['vectors'][M];va=U.row(index[M][wa]);vb=U.row(index[M][wb])
                products=int(va[A]*vb[B]+(va[B]*vb[A]if A!=B else 0))
                total+=Fraction(sa*sb*math.comb(n,i)*(2*i-n)*products,math.comb(J,(J-M)//2))
            row.append(total)
        rows.append(primitive(row))
    assert rows==data['selectedRows'] and len(rows)==100
    assert modular_rank(rows)==100
    K=s.Matrix(data['kernel'])
    assert K.shape==(111,11) and all(x.is_Integer for x in K)
    assert modular_rank(K.tolist())==11
    assert s.Matrix(rows)*K==s.zeros(100,11)
    diagonal=[i for i,(_,a,b)in enumerate(variables)if a==b]
    assert diagonal==cert['diagonalVariableIndices'] and len(diagonal)==49
    coefficients=[s.Integer(v)for v in cert['positiveCoefficients']]
    assert len(coefficients)==49 and all(v>0 for v in coefficients)
    assert s.Matrix([coefficients])*K[diagonal,:]==s.zeros(1,11)
    # Positive semidefinite invariant blocks have nonnegative diagonal
    # coefficients even in these nonorthogonal highest-weight coordinates:
    # X=U C U^dagger with U injective is positive iff C is positive. The strict
    # positive identity forces every diagonal, hence every entire block, to
    # vanish. A twirl of a normalized nonzero exterior vector has trace one.
    return {'passed':True,'n':11,'order':1,'excludedSubspaceDimension':5,
      'appliesToArbitraryComplexSubspaces':True,'copies':5,'exteriorPower':5,
      'exteriorSpaceDimension':dimension,'highestWeightVectorsChecked':highest_count,
      'loweredNormChecks':lowered_norm_checks,'realSymmetricInvariantVariables':111,
      'selectedContractionEquations':100,'exactConstraintRank':100,
      'exactKernelDimension':11,'strictlyPositiveDiagonalTerms':49,
      'positiveIdentityVerifiedExactly':True,'allPositiveSemidefiniteCandidatesVanish':True,
      'rankVerificationPrime':1000003,
      'method':'Exact physical exterior-power basis, rational first-leg contractions, integer nullspace checks and strictly positive diagonal certificate.',
      'runtimeSeconds':time.time()-start}

if __name__=='__main__':
    result=verify();(BASE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
