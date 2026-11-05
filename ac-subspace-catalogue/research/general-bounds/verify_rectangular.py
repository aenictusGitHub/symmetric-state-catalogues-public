#!/usr/bin/env python3
"""Exact low-dimensional audit of the rectangular-tensor dimension proof.

The theorem is representation-theoretic for arbitrary N; these checks validate
its concrete Clebsch--Gordan convention, normalization, product-rank support,
and rectangular covariant-channel eigenvalues without floating-point numbers.
Requires SymPy only. Output is an exact-arithmetic audit, not a proof by sampling.
"""
from pathlib import Path
from functools import lru_cache
import json,time
import sympy as s
from sympy.physics.wigner import clebsch_gordan,wigner_6j
BASE=Path(__file__).resolve().parent

@lru_cache(None)
def tensor(n,m,L2,q2):
    """Sparse matrix of HS-normalized Hom(H_n,H_m) tensor of rank L2/2."""
    j=s.Rational(n,2);jp=s.Rational(m,2);L=s.Rational(L2,2);q=s.Rational(q2,2)
    entries={}
    for a in range(n+1):
        ma=j-a
        # Output magnetic number must equal input magnetic number + q.
        bb=jp-ma-q
        if bb.is_Integer and 0<=bb<=m:
            b=int(bb);value=(-1)**a*clebsch_gordan(jp,j,L,jp-b,-ma,q)
            if value:entries[b,a]=value
    return entries

def hs(A,B):
    return s.simplify(sum(s.conjugate(v)*B.get(pos,0)for pos,v in A.items()))

def product_adjoint(A,B):
    out={}
    for (r,a),v in A.items():
        for (rr,b),w in B.items():
            if r==rr:out[a,b]=out.get((a,b),0)+s.conjugate(v)*w
    return {pos:s.simplify(v)for pos,v in out.items()if v!=0}

def channel(E,A):
    out={}
    for B in E:
        for (r,a),v in B.items():
            for (z,b),w in B.items():
                x=A.get((a,b),0)
                if x:out[r,z]=out.get((r,z),0)+v*x*s.conjugate(w)
    return {pos:s.simplify(v)for pos,v in out.items()if v!=0}

def verify():
    start=time.time();hs_checks=0;product_checks=0;scalar_checks=0;channel_checks=0;cases=[]
    for n in range(2,7):
        for h in range(n+1):
            m=n-h
            # This list contains every rank that can be used for any t<=N.
            labels=[(L2,q2)for L2 in range(h,n+1,2)for q2 in range(-L2,L2+1,2)]
            operators={label:tensor(n,m,*label)for label in labels}
            for i,a in enumerate(labels):
                for b in labels[i:]:
                    assert hs(operators[a],operators[b])==int(a==b);hs_checks+=1
            for a in labels:
                for b in labels:
                    A=product_adjoint(operators[a],operators[b])
                    # The scalar coefficient is fixed by the Hilbert--Schmidt Gram.
                    assert s.simplify(sum(A.get((r,r),0)for r in range(n+1)))==int(a==b)
                    scalar_checks+=1
                    # q of E_a^dagger E_b equals q_b-q_a; all other tensor
                    # components vanish by the explicit magnetic selection rule.
                    q2=b[1]-a[1];assert q2%2==0;q=q2//2
                    for L in range(max((a[0]+b[0])//2+1,abs(q)),n+1):
                        assert hs(tensor(n,n,2*L,2*q),A)==0;product_checks+=1
            # Audit individual rectangular rank channels against the 6j formula.
            j=s.Rational(n,2);jp=s.Rational(m,2)
            for L2 in range(h,n+1,2):
                ell=s.Rational(L2,2);E=[operators[L2,q2]for q2 in range(-L2,L2+1,2)]
                for L in range(m+1):
                    source=tensor(n,n,2*L,0);target=tensor(m,m,2*L,0)
                    transformed=channel(E,source);actual=hs(target,transformed)
                    predicted=(-1)**int(j+jp+ell+L)*(L2+1)*wigner_6j(j,j,L,jp,jp,ell)
                    assert s.simplify(actual-predicted)==0,(n,h,L2,L,actual,predicted)
                    assert all(s.simplify(transformed.get(pos,0)-predicted*target.get(pos,0))==0 for pos in set(transformed)|set(target))
                    channel_checks+=1
            cases.append({'N':n,'shift':h,'targetDimension':m+1,'operatorsChecked':len(labels)})
    return {'passed':True,'arithmetic':'exact SymPy algebraic numbers','minimumN':2,'maximumN':6,
      'parameterCases':len(cases),'hilbertSchmidtGramEntries':hs_checks,
      'productScalarTraceChecks':scalar_checks,'forbiddenHighRankProductComponents':product_checks,
      'rectangularSixJEigenvalueChecks':channel_checks,'coversBothShiftParities':True,
      'scope':'Independent convention and normalization audit. The representation-theoretic proof, not these finite checks, establishes the bound for every N,t,k.',
      'runtimeSeconds':time.time()-start,'cases':cases}

if __name__=='__main__':
    result=verify();(BASE/'rectangular-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
