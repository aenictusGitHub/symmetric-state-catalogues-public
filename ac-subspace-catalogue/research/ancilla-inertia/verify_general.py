#!/usr/bin/env python3
"""Symbolic all-N identities and exact physical N13/N7 controls.

The proof of the dimension bound is analytic, not inferred from finite tests.
"""
from pathlib import Path
import json
import sympy as s
from bounds import upper_bound,eigenvalues
ROOT=Path(__file__).resolve().parent

def spin(n):
    j=s.Rational(n,2);d=n+1
    Z=s.diag(*[j-r for r in range(d)]);Up=s.zeros(d)
    for r in range(1,d):Up[r-1,r]=s.sqrt(r*(n-r+1))
    return [(Up+Up.T)/2,(Up-Up.T)/(2*s.I),Z]

def simplify_matrix(M):return M.applyfunc(s.simplify)

def physical_control(n,B):
    k=B.cols;J=spin(n);S=spin(3)
    assert simplify_matrix(B.H*B)==s.eye(k)
    c=s.Rational(n*(n+2),12)
    for A in J:assert simplify_matrix(B.H*A*B)==s.zeros(k)
    for a,A in enumerate(J):
        for b,C in enumerate(J):
            assert simplify_matrix(B.H*A*C*B-c*int(a==b)*s.eye(k))==s.zeros(k)
    X=sum((s.kronecker_product(A,C)for A,C in zip(J,S)),s.zeros(4*(n+1)))
    W=X*X-s.Rational(n-4,2)*X-s.Rational(5*n*(n+2),16)*s.eye(4*(n+1))
    V=s.kronecker_product(B,s.eye(4))
    compression=simplify_matrix(V.H*W*V)
    assert compression==s.zeros(4*k)
    return {'n':n,'dimension':k,'ancillaDimension':4,'physicalWitnessMatrixDimension':4*(n+1),
            'compressedWitnessEntriesCheckedExactly':(4*k)**2,
            'wholeSubspaceFirstAndSecondMomentsCheckedExactly':True,
            'compressedWitnessIsExactlyZero':True}

def verify():
    N=s.symbols('N',integer=True,positive=True);j=N/2;anc=s.Rational(3,2)
    L=[j-anc,j-anc+1,j-anc+2,j-anc+3]
    actual_x=[s.expand((ell*(ell+1)-j*(j+1)-anc*(anc+1))/2)for ell in L]
    expected_x=[-3*(N+2)/4,-(N+8)/4,(N-6)/4,3*N/4]
    assert all(s.simplify(a-b)==0 for a,b in zip(actual_x,expected_x))
    actual_w=[s.factor(x*x-(N-4)*x/2-5*N*(N+2)/16)for x in actual_x]
    expected_w=[(N+2)*(5*N-3)/8,-N*(N-7)/8,-(3*N*N+N+6)/8,-N*(N-7)/8]
    assert all(s.factor(a-b)==0 for a,b in zip(actual_w,expected_w))
    mult=[s.expand(2*ell+1)for ell in L]
    assert mult==[N-2,N,N+2,N+4]
    assert s.simplify(sum(mult)-4*(N+1))==0
    assert s.simplify(sum(mult[1:])-(3*N+6))==0
    # This is the only scalar used in compressing X^2 for every admissible N.
    S=spin(3)
    assert sum((A*A for A in S),s.zeros(4))==s.Rational(15,4)*s.eye(4)
    c=N*(N+2)/12
    assert s.simplify(c*s.Rational(15,4)-5*N*(N+2)/16)==0
    # Signs for the entire domain follow by substituting N=h+8: all remaining
    # polynomial coefficients of the relevant signed eigenvalues are positive.
    h=s.symbols('h',nonnegative=True)
    shifted=[s.Poly(s.expand((a if i==0 else -a).subs(N,h+8)),h)for i,a in enumerate(expected_w)]
    assert all(all(q>0 for q in p.all_coeffs())for p in shifted)
    B13=s.zeros(14,2)
    B13[0,0]=s.sqrt(10)/6;B13[9,0]=s.sqrt(26)/6
    B13[4,1]=-s.sqrt(26)/6;B13[13,1]=s.sqrt(10)/6
    B7=s.zeros(8,2)
    B7[0,0]=s.sqrt(s.Rational(3,10));B7[5,0]=s.sqrt(s.Rational(7,10))
    B7[2,1]=s.sqrt(s.Rational(7,10));B7[7,1]=-s.sqrt(s.Rational(3,10))
    controls=[physical_control(13,B13),physical_control(7,B7)]
    assert sum(m for v,m in eigenvalues(7)if v==0)==18
    checks=[]
    for n in range(8,201):
        ev=eigenvalues(n)
        assert sum(m for v,m in ev if v>0)==n-2
        assert sum(m for v,m in ev if v<0)==3*n+6
        assert all(v!=0 for v,m in ev)
        assert upper_bound(n,2)==(n-2)//4
        checks.append(n)
    report={'passed':True,'theorem':'For every N >= 8, every second-order AC subspace has 4k <= N-2.',
            'scope':'Arbitrary complex symmetric subspaces; every higher AC order also satisfies the bound.',
            'generalSymbolicIdentities':{'totalSpinEigenvalues':list(map(str,actual_x)),
                'quadraticWitnessEigenvalues':list(map(str,actual_w)),
                'multiplicities':list(map(str,mult)),
                'ancillaCasimir':'15/4','compressionScalar':'5*N*(N+2)/16',
                'allNSignsProvedByPositiveShiftedPolynomialCoefficients':True},
            'physicalControls':controls,
            'n13Result':{'upper':2,'inertia':{'positive':11,'negative':45,'zero':0},
                         'eigenvalues':[{'value':str(v),'multiplicity':m}for v,m in eigenvalues(13)]},
            'n7Exception':{'zeroEigenvalueMultiplicity':18,'knownDimensionTwoCodePreserved':True},
            'additionalIntegerImplementationChecks':{'firstN':checks[0],'lastN':checks[-1],'count':len(checks)},
            'interpretation':'All-size algebra and sign assertions are symbolic. Finite controls validate conventions and implementation, not the general theorem.'}
    (ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'n13':report['n13Result'],'physicalControls':controls}))
    return report
if __name__=='__main__':verify()
