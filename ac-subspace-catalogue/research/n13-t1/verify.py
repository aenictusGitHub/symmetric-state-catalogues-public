#!/usr/bin/env python3
"""Exact certificate: k_max(13,1)=5 for complex symmetric AC subspaces.

Requires NumPy and SymPy; no optimizer or python-flint is needed. NumPy is used
only for finite-field elimination with integers below 10^6 (products <10^12).
All physical rows, kernel identities and positive weights are checked exactly.
The compressed JSON is a lossless encoding of integers, not numerical data.
"""
import gzip,itertools,json,math,time
from pathlib import Path
import numpy as np
import sympy as s
BASE=Path(__file__).resolve().parent

def canonical(word):
    if len(set(word))!=len(word):return None,0
    return tuple(sorted(word)),(-1)**sum(word[a]>word[b]for a in range(len(word))for b in range(a+1,len(word)))

def modular_rank(matrix,prime=1000003):
    if not matrix:return 0
    A=np.array([[int(v)%prime for v in row]for row in matrix],dtype=np.int64)
    rank=0
    for col in range(A.shape[1]):
        candidates=np.flatnonzero(A[rank:,col])
        if not len(candidates):continue
        pivot=rank+int(candidates[0]);A[[rank,pivot]]=A[[pivot,rank]]
        A[rank,col:]=(A[rank,col:]*pow(int(A[rank,col]),-1,prime))%prime
        if rank+1<A.shape[0]:
            A[rank+1:,col:]=(A[rank+1:,col:]-A[rank+1:,col,None]*A[rank,None,col:])%prime
        rank+=1
        if rank==A.shape[0]:break
    return rank

def verify_lower(directory=None):
    directory=Path(directory)if directory else BASE
    data=json.loads((directory/'lower-basis.json').read_text());assert(data['n'],data['order'],data['dimension'])==(13,1,5)
    V=s.Matrix.hstack(*[s.Matrix([s.sympify(x)for x in row])for row in data['basis']])
    simplify=lambda A:A.applyfunc(lambda x:s.simplify(s.expand_complex(x)))
    a=s.Rational(65,282);b=s.Rational(35,94);c=s.Rational(56,141)
    z=(-21*s.sqrt(182)+s.I*s.sqrt(92794))/416
    assert a+b+c==1 and s.simplify(z*s.conjugate(z))==1
    displayed=[{3:s.sqrt(a)*z,6:-s.sqrt(b),9:s.sqrt(c)},
      {4:s.sqrt(c),7:-s.sqrt(b),10:s.sqrt(a)*z},
      {5:1/s.sqrt(2),8:1/s.sqrt(2)},
      {1:s.sqrt(s.Rational(13,24)),13:s.sqrt(s.Rational(11,24))},
      {0:s.sqrt(s.Rational(11,24)),12:-s.sqrt(s.Rational(13,24))}]
    W=s.Matrix(14,5,lambda r,i:displayed[i].get(r,0))
    assert simplify(V-W)==s.zeros(14,5)
    assert simplify(V.H*V)==s.eye(5)
    n=13;j=s.Rational(n,2);Jz=s.diag(*[j-r for r in range(n+1)])
    Jplus=s.zeros(n+1)
    for r in range(1,n+1):Jplus[r-1,r]=s.sqrt(r*(n-r+1))
    assert simplify(V.H*Jz*V)==s.zeros(5)
    assert simplify(V.H*Jplus*V)==s.zeros(5)
    # The adjoint gives J_minus, and hence Jx and Jy. Direct one-qubit cross
    # reductions also check the full subspace statement in Dicke convention.
    split=[]
    for a in range(5):
        split.append(s.Matrix(2,n,lambda r,q:V[r+q,a]*s.sqrt(s.binomial(n-1,q)/s.binomial(n,r+q))))
    for a in range(5):
        for b in range(5):assert simplify(split[a]*split[b].H-s.eye(2)*s.Rational(int(a==b),2))==s.zeros(2)
    return {'passed':True,'n':13,'order':1,'dimension':5,'orthonormalityExactlyVerified':True,
      'compressedSpinEntriesExactlyVerified':50,'orderedCrossReductionEntriesExactlyVerified':100,'displayedCoefficientComparisons':70}

def verify(directory=None):
    start=time.time();directory=Path(directory)if directory else BASE
    raw=directory/'kernel.json'
    data=json.loads(raw.read_text())if raw.exists()else json.loads(gzip.decompress((directory/'kernel.json.gz').read_bytes()))
    cert=json.loads((directory/'certificate.json').read_text());n=data['n'];k=data['k']
    assert(n,k,data['order'],data['exteriorPower'])==(13,6,1,6)
    assert(cert['n'],cert['order'],cert['excludedDimension'],cert['copies'])==(13,1,6,6)
    words=list(itertools.combinations(range(n+1),k));by={}
    for w in words:by.setdefault(2*sum(w)-k*n,[]).append(w)
    idx={M:{w:i for i,w in enumerate(ws)}for M,ws in by.items()}
    # Certificate coordinates use ascending magnetic labels i=m+j, and raw
    # one-spin basis b_i=sqrt(binomial(n,i))|j,m=i-j>. This convention is
    # equivalent to the reversed usual Dicke ordering used in lower-basis.json.
    def shift(U,M,direction):
        dst=by.get(M+2*direction,[]);mult=len(U[0]);out=[[0]*mult for _ in dst]
        for a,w in enumerate(by[M]):
            for leg in range(k):
                i=w[leg];v=list(w);v[leg]+=direction
                if not 0<=v[leg]<=n:continue
                q,sign=canonical(v)
                if q is None:continue
                row=idx[M+2*direction][q];factor=sign*(i+1 if direction==1 else n-i+1)
                for c in range(mult):out[row][c]+=factor*U[a][c]
        return out
    expected=[(M,len(by[M])-len(by.get(M+2,[])))for M in sorted(by)if M>=0 and len(by[M])>len(by.get(M+2,[]))]
    assert[(b['J2'],b['multiplicity'])for b in data['blocks']]==expected
    weights={M:[math.factorial(k)*math.prod(math.comb(n,i)for i in w)for w in ws]for M,ws in by.items()}
    blocks=[];highest=0;normchecks=0
    for block in data['blocks']:
        J=block['J2'];mult=block['multiplicity'];H=block['highest']
        assert len(H)==len(by[J])and all(len(row)==mult for row in H)
        assert modular_rank(H)==mult
        assert not any(v for row in shift(H,J,1)for v in row)
        U={J:H}
        for step in range(1,J+1):
            M=J-2*step;rows=shift(U[M+2],M+2,-1)
            assert all(v%step==0 for row in rows for v in row)
            U[M]=[[v//step for v in row]for row in rows]
        norms=[sum(weights[J][i]*H[i][a]**2 for i in range(len(H)))for a in range(mult)]
        assert norms==block['norms']and all(v>0 for v in norms)
        for M,V in U.items():
            for a in range(mult):
                assert sum(weights[M][i]*V[i][a]**2 for i in range(len(V)))==math.comb(J,(J-M)//2)*norms[a]
                normchecks+=1
        blocks.append({'J2':J,'multiplicity':mult,'U':U});highest+=mult
    dim=sum((b['J2']+1)*b['multiplicity']for b in blocks)
    assert dim==data['dimension']==math.comb(14,6)==3003 and highest==151
    variables=[tuple(v)for v in data['variables']]
    assert variables==[(ib,a,b)for ib,B in enumerate(blocks)for a in range(B['multiplicity'])for b in range(a,B['multiplicity'])]
    assert len(variables)==735
    low={}
    for w in itertools.combinations(range(n+1),k-1):low.setdefault(2*sum(w)-(k-1)*n,[]).append(w)
    prepared={}
    for Ml in {row[0]for row in data['selectedConstraints']}:
        denominator=math.lcm(*[math.comb(B['J2'],(B['J2']-(Ml+2*i-n))//2)for B in blocks for i in range(n+1)if Ml+2*i-n in B['U']])
        for ib,B in enumerate(blocks):
            for i in range(n+1):
                M=Ml+2*i-n
                if M not in B['U']:continue
                factor=math.comb(n,i)*(2*i-n)*denominator//math.comb(B['J2'],(B['J2']-M)//2)
                V=[]
                for alpha in low[Ml]:
                    w,sign=canonical((i,)+alpha)
                    V.append([sign*x for x in B['U'][M][idx[M][w]]]if w is not None else[0]*B['multiplicity'])
                prepared[Ml,ib,i]=(factor,V)
    rows=[]
    for Ml,a,b in data['selectedConstraints']:
        row=[]
        for ib,B in enumerate(blocks):
            mult=B['multiplicity'];part=[0]*(mult*(mult+1)//2)
            for i in range(n+1):
                item=prepared.get((Ml,ib,i))
                if item is None:continue
                factor,V=item;va=V[a];vb=V[b];pos=0
                for A in range(mult):
                    for C in range(A,mult):
                        part[pos]+=factor*(va[A]*vb[C]+(va[C]*vb[A]if A!=C else 0));pos+=1
            row.extend(part)
        divisor=math.gcd(*row);assert divisor;rows.append([v//divisor for v in row])
    assert rows==data['selectedRows']and len(rows)==716
    assert modular_rank(rows)==716
    K=data['kernel'];assert len(K)==735 and all(len(row)==19 for row in K)
    assert modular_rank(K)==19
    for row in rows:
        nonzero=[(i,v)for i,v in enumerate(row)if v]
        assert all(sum(v*K[i][col]for i,v in nonzero)==0 for col in range(19))
    diagonal=[i for i,(_,a,b)in enumerate(variables)if a==b]
    assert diagonal==cert['diagonalVariableIndices']and len(diagonal)==151
    coeff=cert['positiveCoefficients'];assert len(coeff)==151 and all(isinstance(v,int)and v>0 for v in coeff)
    assert all(sum(c*K[i][col]for c,i in zip(coeff,diagonal))==0 for col in range(19))
    lower=verify_lower(directory)
    return {'passed':True,'n':13,'order':1,'maximumDimension':5,'excludedSubspaceDimension':6,
      'appliesToArbitraryComplexSubspaces':True,'copies':6,'exteriorPower':6,
      'exteriorSpaceDimension':dim,'highestWeightVectorsChecked':highest,'loweredNormChecks':normchecks,
      'realSymmetricInvariantVariables':735,'selectedContractionEquations':716,'exactConstraintRank':716,
      'exactKernelDimension':19,'strictlyPositiveDiagonalTerms':151,'positiveIdentityVerifiedExactly':True,
      'allPositiveSemidefiniteCandidatesVanish':True,'rankVerificationPrime':1000003,
      'lowerConstruction':lower,'runtimeSeconds':time.time()-start}

if __name__=='__main__':
    result=verify();(BASE/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
