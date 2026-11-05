#!/usr/bin/env python3
"""Exact certificate for three reconstructed radical anticoherent subspaces.

Run with Python 3 and SymPy. Does not execute Julia or trust small numerical
residuals. Every matrix element of every cross partial trace is reconstructed
from the exact Dicke amplitudes, and orthonormality is checked independently.
"""
from pathlib import Path
import hashlib,json
import sympy as s
ROOT=Path(__file__).resolve().parent

def simp(x):return s.simplify(s.expand_complex(x))

def split(v,t):
    n=len(v)-1
    return s.Matrix(t+1,n-t+1,lambda i,j:v[i+j]*s.sqrt(s.binomial(t,i)*s.binomial(n-t,j)/s.binomial(n,i+j)))

def main():
    path=ROOT/'bases.json';data=json.loads(path.read_text());report=[]
    for record in data['subspaces']:
        n,t,k=record['n'],record['order'],record['dimension']
        basis=[s.Matrix([s.sympify(x)for x in row['exactCoefficients']])for row in record['basis']]
        assert len(basis)==k and all(len(v)==n+1 for v in basis)
        gram=(s.Matrix.hstack(*basis).H*s.Matrix.hstack(*basis)).applyfunc(simp)
        assert gram==s.eye(k),(record['id'],'orthonormality',gram)
        checks=[]
        for order in range(1,t+1):
            matrices=[split(v,order)for v in basis]
            for i,A in enumerate(matrices):
                for j,B in enumerate(matrices):
                    residual=(A*B.H-s.eye(order+1)*int(i==j)/(order+1)).applyfunc(simp)
                    assert residual==s.zeros(order+1),(record['id'],order,i,j,residual)
                    checks.append({'order':order,'i':i,'j':j,'entriesVerified':(order+1)**2})
        # An actual basis vector fails the next-order condition, so the entire
        # subspace has exactly the stated order, not merely a lower bound.
        nextsplit=split(basis[0],t+1)
        nextres=(nextsplit*nextsplit.H-s.eye(t+2)/(t+2)).applyfunc(simp)
        nextdefect=simp(sum(x*s.conjugate(x)for x in nextres))
        assert nextdefect!=0
        probs=[[simp(x*s.conjugate(x))for x in v if x!=0]for v in basis]
        assert [[str(p)for p in row]for row in probs]==record['derivation']['probabilities']
        report.append({'id':record['id'],'orthonormal':True,'allCrossReductionsVerified':True,'checks':checks,
            'scalarEntriesVerified':sum(row['entriesVerified']for row in checks),
            'firstBasisVectorNextOrderDefect':str(nextdefect),'exactSubspaceOrder':t})
        print('Verified',record['id'],'through order',t,flush=True)
    output={'basesSHA256':hashlib.sha256(path.read_bytes()).hexdigest(),'allChecksPassed':True,'verificationMethod':'Exact symbolic Dicke bipartition matrices and their cross products.','subspaces':report}
    (ROOT/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print('All three exact subspaces verified.')
if __name__=='__main__':main()
