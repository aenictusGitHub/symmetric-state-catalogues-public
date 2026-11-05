#!/usr/bin/env python3
"""Exact audit of the Gram-bound saturation obstruction for t>=2.

For a k-dimensional t-anticoherent subspace at spin j=n/2, t>=2 gives
P J_a P=0 and P J_a J_b P=c delta_ab P, c=j(j+1)/3.
The four k-dimensional spaces C,J_x C,J_y C,J_z C are mutually orthogonal,
so 4k<=n+1. If equality holds, I=P+Phi(P)/c, Phi(X)=sum_a J_a X J_a.
On rank-L irreducible operators, Phi has eigenvalue j(j+1)-L(L+1)/2.
Thus every nonzero nonscalar component of P must satisfy
3 L(L+1)=2 n(n+2). If no integer rank 1<=L<=n solves this, P=I/4,
contradicting idempotence. This is a necessary condition for saturation,
not a claim that every permitted rank produces a valid subspace.

The verifier uses exact SymPy matrices at both integer and half-integer
spin and exact integer rank enumeration. It also checks an explicit
saturating n=7,k=2,t=2 example, without relying on the saturation equation
alone as an existence argument.
"""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parent


def spin_matrices(n):
    j=s.Rational(n,2)
    Jz=s.diag(*[j-r for r in range(n+1)])
    Jplus=s.zeros(n+1)
    for r in range(1,n+1):Jplus[r-1,r]=s.sqrt(r*(n-r+1))
    Jminus=Jplus.T
    return j,Jplus,Jminus,(Jplus+Jminus)/2,(Jplus-Jminus)/(2*s.I),Jz


def casimir_checks(n):
    j,Jplus,Jminus,Jx,Jy,Jz=spin_matrices(n)
    c=j*(j+1)/3
    assert (Jx*Jx+Jy*Jy+Jz*Jz).applyfunc(s.simplify)==j*(j+1)*s.eye(n+1)
    phi=lambda X: Jx*X*Jx+Jy*X*Jy+Jz*X*Jz
    highest=s.eye(n+1);rows=[]
    for rank in range(n+1):
        # J_+^L is a nonzero highest-weight operator of irreducible rank L.
        assert highest!=s.zeros(n+1)
        assert (Jz*highest-highest*Jz-rank*highest).applyfunc(s.simplify)==s.zeros(n+1)
        assert (Jplus*highest-highest*Jplus).applyfunc(s.simplify)==s.zeros(n+1)
        eigenvalue=j*(j+1)-s.Rational(rank*(rank+1),2)
        residual=(phi(highest)-eigenvalue*highest).applyfunc(s.simplify)
        assert residual==s.zeros(n+1),(n,rank)
        factor=s.simplify(1+eigenvalue/c)
        expected=s.Rational(2*n*(n+2)-3*rank*(rank+1),1)/(s.Rational(n*(n+2),2))
        assert s.simplify(factor-expected)==0
        rows.append({'rank':rank,'phiEigenvalue':str(eigenvalue),
                     'saturationMultiplier':str(factor),'highestWeightIdentityVerified':True})
        highest=highest*Jplus
    assert highest==s.zeros(n+1)
    return {'n':n,'spin':str(j),'c':str(c),'allRanksVerified':True,'ranks':rows}


def saturation_ranks(n):
    return [rank for rank in range(1,n+1) if 3*rank*(rank+1)==2*n*(n+2)]


def saturation_bound(n):
    gram=(n+1)//4;ranks=saturation_ranks(n)
    excluded=(n+1)%4==0 and not ranks
    return {'n':n,'gramUpper':gram,'allowedSaturationRanks':ranks,
            'saturationExcludedByRankCondition':excluded,
            'upperForEveryOrderAtLeastTwo':gram-int(excluded)}


def verify_n7_equality():
    n=7;j,Jplus,Jminus,Jx,Jy,Jz=spin_matrices(n);c=j*(j+1)/3
    v=s.zeros(8,1);w=s.zeros(8,1)
    v[2]=s.sqrt(s.Rational(7,10));v[7]=s.sqrt(s.Rational(3,10))
    w[0]=-s.sqrt(s.Rational(3,10));w[5]=s.sqrt(s.Rational(7,10))
    V=v.row_join(w);P=V*V.H
    assert (V.H*V).applyfunc(s.simplify)==s.eye(2)
    operators=[Jx,Jy,Jz]
    for a,A in enumerate(operators):
        assert (V.H*A*V).applyfunc(s.simplify)==s.zeros(2)
        for b,B in enumerate(operators):
            assert (V.H*A*B*V-c*int(a==b)*s.eye(2)).applyfunc(s.simplify)==s.zeros(2)
    phiP=sum([A*P*A for A in operators],s.zeros(8))
    assert (P+phiP/c).applyfunc(s.simplify)==s.eye(8)
    nonscalar=P-s.eye(8)/4
    # Its only possible nonscalar rank is L=6, with eigenvalue -c.
    phiQ=sum([A*nonscalar*A for A in operators],s.zeros(8))
    assert nonscalar!=s.zeros(8) and (phiQ+c*nonscalar).applyfunc(s.simplify)==s.zeros(8)
    return {'n':7,'dimension':2,'order':2,'orthonormal':True,
            'allFirstAndSecondCompressedSpinIdentitiesVerified':True,
            'resolutionOfIdentityVerified':True,'onlyNonscalarRank':6,
            'basisDickeCoefficients':[[str(x)for x in vec]for vec in (v,w)]}


def main():
    checked=[casimir_checks(n)for n in [11,12]]
    requested=[saturation_bound(n)for n in [7,11,15,19,35]]
    assert {row['n']:row['upperForEveryOrderAtLeastTwo']for row in requested}=={7:2,11:2,15:3,19:4,35:8}
    assert saturation_ranks(7)==[6]
    assert all(not saturation_ranks(n)for n in [11,15,19,35])
    assert s.Rational(1,4)**2-s.Rational(1,4)==-s.Rational(3,16)
    n7=verify_n7_equality()
    result={'allChecksPassed':True,'method':'Exact spin matrices, highest-weight operator Casimir identities, and integer rank enumeration.',
       'theoremScope':'Every pure symmetric anticoherent subspace of order t>=2.',
       'gramBound':'4k<=n+1',
       'saturationIdentity':'I=P+sum_a J_a P J_a/c, c=n(n+2)/12',
       'necessaryNonscalarRankEquation':'3 L(L+1)=2 n(n+2), 1<=L<=n',
       'scalarContradictionWhenNoRankExists':'P=I/4 implies P^2-P=-3I/16, so P is not a projector.',
       'casimirChecks':checked,'requestedConsequences':requested,
       'integerEnumerationThroughN100':[saturation_bound(n)for n in range(1,101)],
       'n7SaturationWitness':n7,
       'limitations':'An allowed rank is only a necessary condition in general. Existence at n=7 is separately established here by its explicit basis. The strengthened upper bounds do not by themselves prove attainment at other sizes.'}
    (ROOT/'gram-saturation-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Exact Gram-saturation audit passed.')
    for row in requested:print('n=',row['n'],' t>=2 upper=',row['upperForEveryOrderAtLeastTwo'],' allowed saturation ranks=',row['allowedSaturationRanks'])
if __name__=='__main__':main()
