#!/usr/bin/env python3
"""Independent exact audit of the N9,t3 four-copy obstruction.

Reconstructs selected contraction equations from explicit symmetrized-word
multiplicities, without importing the discovery or certificate-builder code.
The selected equations alone already yield the contradiction.
"""
import argparse
import itertools
import json
import math
import time
from fractions import Fraction
from functools import reduce
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix


def primitive(row):
    row=[Fraction(x) for x in row]
    denominator=math.lcm(*(x.denominator for x in row))
    integers=[int(x*denominator) for x in row]
    divisor=reduce(math.gcd,(abs(x) for x in integers))
    assert divisor
    return [x//divisor for x in integers]


def multiplicity(word):
    return math.factorial(len(word))//math.prod(math.factorial(word.count(k)) for k in set(word))


def verify(data,certificate):
    start=time.time()
    n=data['n'];copies=data['copies'];t=data['t']
    assert (n,copies,t)==(9,4,3)
    words=list(itertools.combinations_with_replacement(range(n+1),copies))
    byweight={}
    for word in words:
        byweight.setdefault(2*sum(word)-copies*n,[]).append(word)
    lesser=list(itertools.combinations_with_replacement(range(n+1),copies-1))
    lowweight={}
    for word in lesser:
        lowweight.setdefault(2*sum(word)-(copies-1)*n,[]).append(word)

    # Collective lowering/raising by acting on every slot of each word.
    def collective(matrix,weight,shift):
        src=byweight[weight];dst=byweight.get(weight+2*shift,[])
        index={word:i for i,word in enumerate(dst)}
        result=s.zeros(len(dst),matrix.cols)
        for column,word in enumerate(src):
            for position,k in enumerate(word):
                target=k+shift
                if not 0<=target<=n:
                    continue
                out=list(word);out[position]=target;out=tuple(sorted(out))
                coefficient=k+1 if shift==1 else n-k+1
                row=index[out]
                for a in range(matrix.cols):
                    result[row,a]+=coefficient*matrix[column,a]
        return result

    expected_spins=[J for J in sorted(byweight) if J>=0 and len(byweight[J])>len(byweight.get(J+2,[]))]
    assert [b['J2'] for b in data['blocks']]==expected_spins
    blocks=[]
    for block in data['blocks']:
        J2=block['J2'];H=s.Matrix([[s.Integer(x) for x in row] for row in block['highestWeight']])
        mult=block['multiplicity']
        assert block['denominatorScale']>0
        assert H.rows==len(byweight[J2]) and H.cols==mult
        assert H.rank()==mult
        assert collective(H,J2,1)==s.zeros(len(byweight.get(J2+2,[])),mult)
        assert mult==len(byweight[J2])-len(byweight.get(J2+2,[]))
        U={J2:H}
        for step in range(1,J2+1):
            M=J2-2*step
            U[M]=collective(U[M+2],M+2,-1)/step
            assert all(x.is_Integer for x in U[M])
        lookup={M:{word:tuple(int(v) for v in matrix.row(i)) for i,word in enumerate(byweight[M])} for M,matrix in U.items()}
        blocks.append({'J2':J2,'multiplicity':mult,'denominatorScale':block['denominatorScale'],'vectors':lookup})
    dimension=sum((b['J2']+1)*b['multiplicity'] for b in blocks)
    assert dimension==math.comb(n+copies,copies)==715
    variables=[tuple(v) for v in data['variables']]
    assert variables==[(ib,a,b) for ib,block in enumerate(blocks) for a in range(block['multiplicity']) for b in range(a,block['multiplicity'])]
    assert len(variables)==88

    reconstructed=[]
    for ell,M3,a,b in data['selectedConstraints']:
        assert 1<=ell<=t
        alpha=lowweight[M3][a];beta=lowweight[M3][b]
        row=[Fraction(0) for _ in variables]
        for iv,(ib,A,B) in enumerate(variables):
            block=blocks[ib];J2=block['J2']
            for i in range(n+1):
                M=M3+2*i-n
                if M not in block['vectors']:
                    continue
                wa=tuple(sorted((i,)+alpha));wb=tuple(sorted((i,)+beta))
                va=block['vectors'][M][wa];vb=block['vectors'][M][wb]
                # A first-leg contraction of averaged words is the ratio of
                # their orbit sizes; the square-root binomial factors pair up.
                ratio=Fraction(multiplicity(alpha)*multiplicity(beta),multiplicity(wa)*multiplicity(wb))
                z=2*i-n
                f=z if ell==1 else 3*z*z-99 if ell==2 else z*(5*z*z-293)
                scalar=Fraction(block['denominatorScale'],math.comb(J2,(J2-M)//2))*math.comb(n,i)*f*ratio
                row[iv]+=scalar*(va[A]*vb[B]+(va[B]*vb[A] if A!=B else 0))
        reconstructed.append(primitive(row))
    assert reconstructed==[[int(x) for x in row] for row in data['selectedRows']]
    matrix=s.Matrix(reconstructed)
    rank=DomainMatrix.from_Matrix(matrix).convert_to(s.QQ).rank()
    assert rank==85
    K=s.Matrix([[s.Integer(x) for x in row] for row in data['kernel']])
    assert K.shape==(88,3) and K.rank()==3
    assert matrix*K==s.zeros(85,3)

    indices=certificate['selectedVariableIndices']
    coefficients=s.Matrix([[s.Integer(x) for x in certificate['positiveCoefficients']]])
    assert all(x>0 for x in coefficients)
    assert all(variables[i][1]==variables[i][2] for i in indices)
    selected=K[indices,:]
    assert selected.rank()==3
    assert coefficients*selected==s.zeros(1,3)
    # Every invariant PSD candidate has nonnegative selected diagonal entries.
    # Their strictly positive combination is zero, so each is zero. The rank3
    # selected kernel rows then force all3 nullspace parameters to vanish.
    return {'n':n,'t':t,'copies':copies,'verifiedExactly':True,
            'method':'Independent exact reconstruction of contraction equations using symmetrized-word orbit multiplicities; no discovery code imported.',
            'fullSymmetricPowerDimension':dimension,'invariantVariables':len(variables),
            'reconstructedConstraints':len(reconstructed),'exactRank':rank,'kernelDimension':K.cols,
            'positiveDiagonalTerms':len(indices),'selectedKernelRank':3,
            'runtimeSeconds':time.time()-start}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('directory',type=Path,nargs='?',default=Path(__file__).resolve().parent)
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('independent-four-copy-verification.json'))
    args=parser.parse_args()
    result=verify(json.loads((args.directory/'kernel.json').read_text()),json.loads((args.directory/'certificate.json').read_text()))
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
